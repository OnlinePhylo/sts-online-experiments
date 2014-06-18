#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import csv
import json
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-k', '--keys', default=['proposal_method', 'trim_count',
                                            'keep_count', 'trim_taxon',
                                            'particle_factor', 'by_length',
                                            'n_taxa', 'tree'],
                   type=lambda s: s.split(','),
                   help="""Additional keys""")
    p.add_argument('control_file', nargs='+')
    p.add_argument('-o', '--outfile', default=sys.stdout,
                   type=argparse.FileType('w'))
    a = p.parse_args()

    control_handles = (argparse.FileType('r')(fname)
                       for fname in a.control_file)
    with a.outfile as ofp:
        w = csv.DictWriter(ofp, a.keys + ['last_ess', 'likelihood_calls',
                                          'sts_json'],
                           lineterminator='\n')
        w.writeheader()
        for fp in control_handles:
            with fp:
                ctrl = json.load(fp)

            online_result_files = ctrl['sts_online']
            row_base = {k: ctrl[k] for k in a.keys}

            for f in online_result_files:
                row = row_base.copy()
                with argparse.FileType('r')(f) as fp:
                    sts_json = json.load(fp)
                last_gen = sts_json['generations'][-1]
                row['sts_json'] = f
                row['last_ess'] = last_gen['ess']
                row['likelihood_calls'] = last_gen['totalUpdatePartialsCalls']
                w.writerow(row)

if __name__ == '__main__':
    main()
