#!/usr/bin/env python
import argparse
import os.path
import sys

from Bio import SeqIO

TEMPLATE = """
begin mrbayes;
    set autoclose=yes nowarn=yes;
    execute {nexus};
    lset nst=1 rates=equal;
    prset statefreqpr=fixed(equal);
    {extra}
    mcmc nruns=2 nchains=3 ngen=5000000 samplefreq=5000 printfreq=50000 file={out_base};
end;
"""

def main():
    p = argparse.ArgumentParser()
    p.add_argument('nexus_path')
    p.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout)
    a = p.parse_args()

    base = os.path.splitext(a.nexus_path)[0]
    n_sequences = sum(True for i in SeqIO.parse(a.nexus_path, 'nexus'))

    extra = ''
    if n_sequences <= 4:
        extra = """propset ParsSPR(Tau,V)$prob=0;
  propset ExtSPR(Tau,V)$prob=0;
  propset ExtTBR(Tau,V)$prob=0;"""

    t = TEMPLATE.format(nexus=a.nexus_path, out_base=base, extra=extra)
    with a.outfile as fp:
        fp.write(t)

if __name__ == '__main__':
    main()
