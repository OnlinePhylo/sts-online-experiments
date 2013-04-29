#!/usr/bin/env python
import argparse
import csv
import json
import sys

import dendropy

def main():
    p = argparse.ArgumentParser()
    p.add_argument('ref_tree')
    p.add_argument('sts_json', nargs='+')
    p.add_argument('-o', '--output', type=argparse.FileType('w'),
        default=sys.stdout)
    a = p.parse_args()

    ref_tree = dendropy.Tree.get_from_path(a.ref_tree, 'newick')

    with a.output as ofp:
        w = csv.writer(ofp, lineterminator='\n')
        w.writerow(['tree', 'posterior', 'generation', 'pruned_taxon', 'pendant_bl', 'ess'])
        for f in a.sts_json:
            with argparse.FileType('r')(f) as fp:
                sts_json = json.load(fp)

                for i, g in enumerate(sts_json['generations']):
                    sequence = g['sequence']
                    ess = g['ess']
                    node = ref_tree.find_node_with_taxon_label(sequence)
                    w.writerow([a.ref_tree, f, i, sequence, node.edge_length, ess])

if __name__ == '__main__':
    main()
