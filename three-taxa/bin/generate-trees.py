#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import dendropy
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-c', '--count', type=int, default=1000)
    p.add_argument('-d', '--distance', type=float, default=0.1)
    p.add_argument('nexus', type=argparse.FileType('w'))
    p.add_argument('tree', type=argparse.FileType('w'))
    a = p.parse_args()

    tree_full = '((A:1e-6,B:1e-6):0.0,C:{0});'.format(a.distance)
    tree = '(A:1e-6,B:1e-6);'

    with a.tree as fp:
        fp.write(tree_full + '\n')

    tl = dendropy.TreeList()
    for _ in xrange(a.count):
        tl.append(dendropy.Tree.get_from_string(tree, 'newick'))

    with a.nexus as fp:
        tl.write_to_stream(fp, 'nexus')

if __name__ == '__main__':
    main()
