#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import json
import sys

import numpy as np
import pandas as pd
from scipy.misc import logsumexp

BINS = 250
RANGE = (0.0, 1.0)

def kl(p, q):
    """Kullback-Leibler divergence D(P || Q) for discrete distributions

    Parameters
    ----------
    p, q : array-like, dtype=float, shape=n
        Discrete probability distributions.
    """
    p = np.asarray(p, dtype=np.float)
    q = np.asarray(q, dtype=np.float)

    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

def hist_of_empirical(empirical_data):
    total_weight = logsumexp(empirical_data['posterior'])
    weight = np.exp(empirical_data['posterior'] - total_weight)
    branch_length = empirical_data['branch_length']
    hist = np.histogram(branch_length,
                        bins=BINS,
                        range=RANGE,
                        weights=weight,
                        density=True)
    return hist


def main():
    p = argparse.ArgumentParser()
    p.add_argument('sts_output', type=argparse.FileType('r'))
    p.add_argument('empirical_output', type=argparse.FileType('r'))
    p.add_argument('-o', '--outfile', default=sys.stdout,
                   type=argparse.FileType('w'))
    a = p.parse_args()

    with a.empirical_output as fp:
        empirical = pd.read_csv(fp)

    empirical_hist, empirical_edges = hist_of_empirical(empirical)

    with a.sts_output as fp:
        doc = json.load(fp)
    sts = pd.DataFrame.from_records(doc['trees'])
    sts_hist, sts_edges = np.histogram(sts['treeLength'],
                                       bins=BINS,
                                       range=RANGE,
                                       density=True)

    output = pd.DataFrame.from_dict({'edge': sts_edges[1:],
                                     'sts_p': sts_hist,
                                     'empirical_p': empirical_hist})
    empirical_hist = empirical_hist / empirical_hist.sum()
    sts_hist = sts_hist / sts_hist.sum()

    with a.outfile as ofp:
        r = {'kl': kl(sts_hist, empirical_hist),
             'ess': doc['generations'][0]['ess'],
             'comparisons': output.to_dict('records')}
        json.dump(r, ofp, indent=2)
        #output.to_csv(ofp, index=False)


if __name__ == '__main__':
    main()
