#!/usr/bin/env python
import argparse
import collections
import csv
import itertools
import subprocess
import re
import sys

PosteriorComparison = collections.namedtuple(
        'PosteriorComparison',
        ['tree1', 'tree2',
         'd11_mean', 'd11_lower' ,'d11_upper',
         'd22_mean', 'd22_lower' ,'d22_upper',
         'd12_mean', 'd12_lower' ,'d12_upper',
         'p_d12_g_d11',
         'p_d12_g_d22'])

def extract_distance_range(comparison, text):
    suffix = r'\s+~ ([.0-9]+)\s+\(([.0-9]+),\s+([.0-9]+)\)'
    return map(float, re.findall(comparison + suffix, text)[0])

def extract_p_greater(n, text):
    r = r'P\(D12\s+>\s+D{0}{0}\) = ([.0-9]+)'.format(n)
    return float(re.findall(r, text)[0])

def compare_pair(tree1, tree2):
    cmd = ['trees-distances', 'compare', tree1, tree2]
    out = subprocess.check_output(cmd)

    d11 = extract_distance_range('D11', out)
    d12 = extract_distance_range('D12', out)
    d22 = extract_distance_range('D22', out)

    pd12_g_d11 = extract_p_greater(1, out)
    pd12_g_d22 = extract_p_greater(2, out)
    return PosteriorComparison(tree1, tree2,
            d11[0], d11[1], d11[2],
            d22[0], d22[1], d22[2],
            d12[0], d12[1], d12[2],
            pd12_g_d11,
            pd12_g_d22)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('files', metavar='file', nargs='+')
    p.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'))
    a = p.parse_args()

    with a.output as fp:
        w = csv.writer(fp, lineterminator='\n')

        # Header
        w.writerow(PosteriorComparison._fields)

        for tree1, tree2 in itertools.combinations(a.files, 2):
            w.writerow(compare_pair(tree1, tree2))

if __name__ == '__main__':
    main()
