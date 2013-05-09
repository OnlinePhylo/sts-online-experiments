#!/usr/bin/env python

from __future__ import division
import argparse
import csv
import functools
import json
import itertools
import math
import os.path
import sys

import dendropy
from dendropy.treecalc import euclidean_distance, symmetric_difference, robinson_foulds_distance

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

tree_parsers = {}

def tree_parser(exts):
    def deco(func):
        for ext in exts:
            tree_parsers[ext] = func
        return func
    return deco

def load_trees_dendropy(fp, schema, **kwargs):
    tl = dendropy.TreeList.get_from_stream(fp, schema, **kwargs)
    for tree in tl:
        tree.log_weight = 0.0
    return tl

load_trees_nexus = functools.partial(load_trees_dendropy, schema='nexus')
tree_parser(['.t', 'nex'])(load_trees_nexus)

@tree_parser(['.json'])
def load_from_json(fp, **kwargs):
    kwargs.pop('tree_offset', None) # No burn-in for STS
    root = json.load(fp)
    def parse_tree(json_tree):
        t = dendropy.Tree.get_from_string(str(json_tree['newickString']), 'newick', **kwargs)
        t.log_weight = json_tree['logWeight']
        return t

    trees = root['trees'][kwargs.pop('tree_offset', 0):]
    return [parse_tree(t) for t in trees]

def main():
    p = argparse.ArgumentParser()
    p.add_argument('reference_tree', help="""Compare to this tree""", type=argparse.FileType('r'))
    p.add_argument('compare_trees', nargs='+', help="""Trees to compare with reference tree""")
    p.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout)
    p.add_argument('-b', '--nexus-burnin', default=0, type=int)
    a = p.parse_args()

    # Trees need a common taxon set
    taxa = dendropy.TaxonSet()

    with a.reference_tree as fp:
        ref_tree = dendropy.Tree.get_from_stream(fp, 'newick', taxon_set=taxa)

    def distances(tree):
        fns = (symmetric_difference, robinson_foulds_distance, euclidean_distance)
        return [fn(ref_tree, tree) for fn in fns]

    with a.output as fp:
        w = csv.writer(fp, lineterminator='\n')
        w.writerow(('file', 'log_weight', 'rf_distance', 'weighted_rf', 'euclidean'))
        for path in a.compare_trees:
            ext = os.path.splitext(path)[1]
            parse = tree_parsers[ext]

            with argparse.FileType('r')(path) as fp:
                trees = parse(fp, tree_offset=a.nexus_burnin, taxon_set=taxa)

                for tree in trees:
                    w.writerow([path, getattr(tree, 'log_weight', 0.0)] + distances(tree))

if __name__ == '__main__':
    main()
