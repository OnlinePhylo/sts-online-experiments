#!/usr/bin/env python
import argparse
import csv
import sys

def key_val(s):
    if not s.count('=') == 1:
        raise ValueError("Format: x=y")
    k, v = s.split('=', 1)
    return k.strip(), v.strip()

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input', default=sys.stdin, type=argparse.FileType('r'))
    p.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'))
    p.add_argument('kvs', metavar='k=v', nargs='+', type=key_val)
    a = p.parse_args()

    keys, values = zip(*a.kvs)
    values = list(values)

    r = csv.reader(a.input)
    headers = list(next(r))
    w = csv.writer(a.output, lineterminator='\n')
    w.writerow(headers + list(keys))

    for row in r:
        w.writerow(list(row) + values)

if __name__ == '__main__':
    main()
