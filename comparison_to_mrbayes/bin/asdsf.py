#!/usr/bin/env python
import argparse
import csv
import functools
import json
import logging
import re
import subprocess
import sys
import tempfile

import dendropy

log = logging.getLogger('asdsf')

def calculate_asdsf_msdsf(tree_path1, tree_path2, min_support=0.1, skip=0):
    """
    returns (ASDSF, MSDSF) pair
    """
    cmd = ['trees-bootstrap',
            '--min-support', str(min_support),
            '--skip', str(skip),
            tree_path1, tree_path2]
    log.info('Running: %s', ' '.join(cmd))
    output = subprocess.check_output(cmd)
    regex = re.compile(r'ASDSF\[min=\d\.\d+\]\s*=\s*(\d\.\d+)\s+MSDSF\s+=\s+(\d\.\d+)')
    m = regex.search(output)
    assert m
    return map(float, m.groups())

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
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO)

    with a.output:
        w = csv.writer(a.output, lineterminator='\n')
        w.writerow(('type', 'file', 'asdsf', 'msdsf'))

        # First, compare between MrBayes
        mb_asdsf, mb_msdsf = calculate_asdsf_msdsf(a.reference_tree1, a.reference_tree2, skip=a.burnin)

        w.writerow(('mrbayes-mrbayes', None, mb_asdsf, mb_msdsf))

        ref_trees = dendropy.TreeList()
        for tree_path in (a.reference_tree1, a.reference_tree2):
            log.info("Reading from %s", tree_path)
            ref_trees.read_from_path(tree_path, 'nexus', tree_offset=a.burnin)

        ntf = functools.partial(tempfile.NamedTemporaryFile, suffix='.nwk', prefix='trees-')

        write_newick = functools.partial(dendropy.TreeList.write_to_stream,
                schema='newick',
                suppress_rooting=True)

        with ntf() as ref_fp:
            log.info('Writing reference trees to %s', ref_fp.name)

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
                    asdsf, msdsf = calculate_asdsf_msdsf(ref_fp.name, sts_fp.name)
                w.writerow(('mrbayes-sts', path, asdsf, msdsf))


if __name__ == '__main__':
    main()
