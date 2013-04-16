#!/usr/bin/env python
from __future__ import division
import argparse
import functools
import itertools
import math
import sys

import dendropy
from dendropy.treecalc import euclidean_distance, robinson_foulds_distance
from dendropy import treesum,  treesplit

def parse_weighted_trees(fp, **kwargs):
    max_log_weight = sys.float_info.min
    trees = []
    for line in fp:
        log_weight, nwk_str = line.split('\t', 1)
        log_weight = float(log_weight)
        if log_weight > max_log_weight:
            max_log_weight = log_weight
        tree = dendropy.Tree.get_from_string(nwk_str, 'newick', **kwargs)
        tree.log_weight = log_weight
        trees.append(tree)

    for tree in trees:
        tree.weight = math.exp(tree.log_weight - max_log_weight)

    return trees

def parse_tree_file(fp, schema='nexus', burnin=0, **kwargs):
    trees = dendropy.tree_source_iter(stream=fp, schema=schema, tree_offset=burnin, **kwargs)
    for tree in trees:
        tree.weight = 1.0
        yield tree

def compute_expectation(fn, trees):
    """
    Computes the weighted expectation of applying ``fn`` to each tree in ``trees``.

    ``fn`` should return a double.
    """
    log_weights = []
    results = []

    for tree in trees:
        log_weights.append(tree.log_weight)
        results.append(fn(tree))

    max_log_weight = max(log_weights)
    result = 0.0
    weight_sum = 0.0
    for lw, r in itertools.izip_longest(log_weights, results):
        w = math.exp(lw - max_log_weight)
        weight_sum += w
        result += r * w

    return result / weight_sum



def main():
    p = argparse.ArgumentParser()
    p.add_argument('reference_tree', help="""Compare to this tree""", type=argparse.FileType('r'))
    p.add_argument('posterior', nargs='?', help="""Tree file [default:
            stdin]""", default=sys.stdin, type=argparse.FileType('r'))
    p.add_argument('-w', '--weighted', default=False, action='store_true',
            help="""Trees are from STS: format '<log_weight>\t<newick_string>'
            [default: false]""")
    p.add_argument('-b', '--burnin', type=int, default=0, help="""Number of
            trees to discard as burn-in [default: %(default)d]""")
    p.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout)

    a = p.parse_args()

    # Trees need a common taxon set
    taxa = dendropy.TaxonSet()

    with a.reference_tree as fp:
        ref_tree = dendropy.Tree.get_from_stream(fp, 'newick', taxon_set=taxa)

    if a.posterior == sys.stdin and sys.stdin.isatty():
        p.error('Attempting to read tree from TTY.')

    with a.posterior as fp:
        if a.weighted:
            trees = parse_weighted_trees(fp, taxon_set=taxa)
            trees = itertools.islice(trees, a.burnin, None)
        else:
            trees = parse_tree_file(fp, burnin=a.burnin, taxon_set=taxa)

        tsum = treesum.TreeSummarizer()
        tsum.weighted_splits = True

        split_distribution = treesplit.SplitDistribution(taxon_set=taxa)
        for tree in trees:
            treesplit.encode_splits(tree)
            split_distribution.count_splits_on_tree(tree)

        con = tsum.tree_from_splits(split_distribution)

        trees1, trees2 = itertools.tee(trees, 2)
        f1 = functools.partial(euclidean_distance, tree2=ref_tree)
        euc_dist = compute_expectation(f1, trees1)
        f2 = functools.partial(robinson_foulds_distance, tree2=ref_tree)
        rf_dist = compute_expectation(f2, trees2)

        with a.output:
            a.output.write('euclidean_distance,rf_distance\n{0},{1}\n'.format(euc_dist, rf_dist))


if __name__ == '__main__':
    main()
