#!/usr/bin/env python
from __future__ import division
import argparse
import csv
import functools
import itertools
import math
import sys

import dendropy
from dendropy.treecalc import euclidean_distance, robinson_foulds_distance
from dendropy import treesum,  treesplit


def main():
    p = argparse.ArgumentParser()
    p.add_argument('reference_tree', help="""Compare to this tree""", type=argparse.FileType('r'))
    p.add_argument('compare_trees', nargs='+', help="""Trees to compare with reference tree""")
    p.add_argument('--schema', choices=('nexus', 'newick'), default='newick',
            help="""Tree format [default: %(default)s]""")
    p.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout)

    a = p.parse_args()

    # Trees need a common taxon set
    taxa = dendropy.TaxonSet()

    with a.reference_tree as fp:
        ref_tree = dendropy.Tree.get_from_stream(fp, 'newick', taxon_set=taxa)

    with a.output as ofp:
        w = csv.writer(ofp, lineterminator='\n')
        w.writerow(['reference', 'query', 'euclidean_distance', 'rf_distance'])
        for tree_path in a.compare_trees:
            tree = dendropy.Tree.get_from_path(tree_path, a.schema, taxon_set=taxa)

            euc_dist = euclidean_distance(tree, ref_tree)
            rf_dist = robinson_foulds_distance(tree, ref_tree)
            w.writerow([a.reference_tree.name, tree_path, euc_dist, rf_dist])

if __name__ == '__main__':
    main()
