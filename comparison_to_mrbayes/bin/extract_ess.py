#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import csv
import json
import sys
import os

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-k', '--keys', default=['proposal_method_name',
                                            'trim_count',
                                            'trim_taxon', 'particle_factor',
                                            'n_taxa'],
                   type=lambda s: s.split(','),
                   help="""Additional keys""")
    p.add_argument('control_file', nargs='+')
    p.add_argument('-o', '--outfile', default=sys.stdout,
                   type=argparse.FileType('w'))
    a = p.parse_args()

    with a.outfile as ofp:
        w = csv.DictWriter(ofp, a.keys + ['generation','ess','unique_particles','sequence'], lineterminator='\n')
        w.writeheader()
        for cf in a.control_file:
            with open(cf) as fp:
                ctrl = json.load(fp)

            online_result_files = ctrl['sts_online']
            row_base = {k: ctrl[k] for k in a.keys}

            for f in online_result_files:
                f = os.path.join(os.path.dirname(cf), os.path.basename(f))
                try:
                    with argparse.FileType('r')(f) as fp:
                        sts_json = json.load(fp)
                except:
                     print(f)
                     continue

                for r in sts_json['generations']:
                    row = row_base.copy()
                    row['generation'] = r['T']
                    row['ess'] = r['ess']
                    row['unique_particles'] = r['uniqueParticles']
                    row['sequence'] = r['sequence']
                    w.writerow(row)

if __name__ == '__main__':
    main()
