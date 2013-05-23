#!/usr/bin/env python
import argparse
import collections
import csv
import functools
import json
import logging
import math
import re
import subprocess
import sys
import tempfile

import dendropy

log = logging.getLogger('asdsf')

ASDSFResult = collections.namedtuple('ASDSFResult', ['asdsf', 'msdsf', 'lod_table'])

def parse_lod_table(fp):
    return [[float(i) for i in line.split(None)[:-1]] for line in fp]

def pp_table_of_lod_table(lod_table):
    return [[10.0 ** i / (1 + 10.0 ** i) for i in row] for row in lod_table]

def calculate_asdsf_msdsf(tree_path1, tree_path2, min_support=0.1, skip=0):
    """
    returns (ASDSF, MSDSF) pair
    """
    with tempfile.NamedTemporaryFile(prefix='LOD_', suffix='.txt') as tf:
        cmd = ['trees-bootstrap',
                '--min-support', str(min_support),
                '--skip', str(skip),
                tree_path1, tree_path2,
                '--LOD-table', tf.name]
        log.info('Running: %s', ' '.join(cmd))
        output = subprocess.check_output(cmd)
        regex = re.compile(r'ASDSF\[min=\d\.\d+\]\s*=\s*(\d\.\d+)\s+MSDSF\s+=\s+(\d\.\d+)')
        m = regex.search(output)
        assert m
        asdsf, msdsf = map(float, m.groups())
        return ASDSFResult(asdsf, msdsf, parse_lod_table(tf))

def load_sts_trees(fp):
    j = json.load(fp)
    tl = dendropy.TreeList()

    trees = (dendropy.Tree.get_from_string(str(tj['newickString']), 'newick')
             for tj in j['trees'])
    tl.extend(trees)
    return tl

def main():
    p = argparse.ArgumentParser()
    p.add_argument('reference_tree1')
    p.add_argument('reference_tree2')
    p.add_argument('sts_trees', nargs='*')
    p.add_argument('-b', '--burnin', type=int, default=250, help="""Number of
            trees to discard as burn-in [default: %(default)d]""")
    p.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout)
    p.add_argument('-p', '--pp-table', type=argparse.FileType('w'))
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO)

    pp_writer = None
    if a.pp_table:
        pp_writer = csv.writer(a.pp_table, lineterminator='\n')
        pp_writer.writerow(('type', 'file1', 'file2', 'pp1', 'pp2'))

    with a.output:
        w = csv.writer(a.output, lineterminator='\n')
        w.writerow(('type', 'file1', 'file2', 'asdsf', 'msdsf'))

        # First, compare between MrBayes
        mb_asdsf, mb_msdsf, mb_lod = calculate_asdsf_msdsf(a.reference_tree1,
                a.reference_tree2, skip=a.burnin)

        if pp_writer:
            rows = (['mrbayes-mrbayes', a.reference_tree1, a.reference_tree2] + row
                    for row in pp_table_of_lod_table(mb_lod))
            pp_writer.writerows(rows)

        w.writerow(('mrbayes-mrbayes', a.reference_tree1, a.reference_tree2,
                    mb_asdsf, mb_msdsf))

        ref_trees = dendropy.TreeList()
        #for tree_path in (a.reference_tree1, a.reference_tree2):
        for tree_path in [a.reference_tree1]:
            log.info("Reading from %s", tree_path)
            ref_trees.read_from_path(tree_path, 'nexus', tree_offset=a.burnin)

        ntf = functools.partial(tempfile.NamedTemporaryFile, suffix='.nwk', prefix='trees-')

        write_newick = functools.partial(dendropy.TreeList.write_to_stream,
                schema='newick',
                suppress_rooting=True)

        with ntf() as ref_fp:
            log.info('Writing %d reference trees to %s', len(ref_trees), ref_fp.name)

            write_newick(ref_trees, ref_fp)
            ref_fp.flush()

            for path in a.sts_trees:
                with open(path) as fp:
                    sts_trees = load_sts_trees(fp)
                with ntf() as sts_fp:
                    log.info('Writing %s to %s', path, sts_fp.name)
                    write_newick(sts_trees, sts_fp)
                    sts_fp.flush()
                    sts_fp.seek(0)
                    asdsf, msdsf, lod = calculate_asdsf_msdsf(ref_fp.name, sts_fp.name)
                    w.writerow(('mrbayes-sts', path, asdsf, msdsf))

                if pp_writer:
                    rows = (['mrbayes-sts', ref_fp.name, path] + row for row in pp_table_of_lod_table(lod))
                    pp_writer.writerows(rows)

        if pp_writer:
            a.pp_table.close()


if __name__ == '__main__':
    main()
