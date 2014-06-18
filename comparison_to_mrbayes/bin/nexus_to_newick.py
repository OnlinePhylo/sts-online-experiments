#!/usr/bin/env python
import argparse
import itertools

from Bio import Phylo

def main():
    p = argparse.ArgumentParser()
    p.add_argument('infile', type=argparse.FileType('r'))
    p.add_argument('outfile', type=argparse.FileType('w'))
    p.add_argument('-b', '--burnin', default=0, type=int, help="""Number of
            trees to discard as burn-in [default: %(default)d]""")

    a = p.parse_args()
    with a.infile, a.outfile:
        trees = Phylo.parse(a.infile, 'nexus')
        Phylo.write(itertools.islice(trees, a.burnin, None),
                    a.outfile, 'newick')

if __name__ == '__main__':
    main()
