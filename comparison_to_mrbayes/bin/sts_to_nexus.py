#!/usr/bin/env python
import argparse
import json
import logging
import math
import sys

import dendropy

def parse_weighted_trees(fp, **kwargs):
    max_log_weight = -sys.float_info.max
    trees = []
    jtrees = json.load(fp)
    for i, tree in enumerate(jtrees['trees']):
        sys.stderr.write('{0:10d}\r'.format(i + 1))
        nwk_str = str(tree['newickString'])  # This can't be unicode for dendropy
        log_weight = tree['logWeight']
        if log_weight > max_log_weight:
            max_log_weight = log_weight
        tree = dendropy.Tree.get_from_string(nwk_str, 'newick', **kwargs)
        tree.log_weight = log_weight
        trees.append(tree)

    logging.info('Maximum log weight: %f (%f)', max_log_weight, math.exp(max_log_weight - max_log_weight))

    for tree in trees:
        tree.weight = math.exp(tree.log_weight - max_log_weight)

    return trees

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--infile', type=argparse.FileType('r'),
            default=sys.stdin)
    p.add_argument('-o', '--outfile', type=argparse.FileType('w'),
            default=sys.stdout)
    a = p.parse_args()

    logging.basicConfig(level=logging.INFO)

    with a.infile as ifp, a.outfile as ofp:
        trees = parse_weighted_trees(ifp)
        tl = dendropy.TreeList(trees)
        tl.write_to_stream(ofp, 'nexus', store_tree_weights=True)

if __name__ == '__main__':
    main()
