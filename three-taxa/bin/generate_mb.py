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
    constraint myclade = A B D;
    mcmcp nruns={nruns} nchains={nchains} ngen={length} samplefreq={samplefreq} printfreq={printfreq} file={out_base} diagnfreq={printfreq};
    {extra}
    mcmc;
end;
"""

def main():
    p = argparse.ArgumentParser()
    p.add_argument('nexus_path')
    mb_group = p.add_argument_group('mrbayes')
    mb_group.add_argument('-l', '--length', type=int, default=1000000)
    mb_group.add_argument('-r', '--runs', type=int, default=2)
    mb_group.add_argument('-c', '--chains', type=int, default=3)
    p.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout)
    a = p.parse_args()

    base = os.path.splitext(a.nexus_path)[0]
    n_sequences = sum(True for i in SeqIO.parse(a.nexus_path, 'nexus'))

    extra = ''
    if n_sequences <= 4:
        extra = """propset ParsSPR(Tau,V)$prob=0;
  propset ExtSPR(Tau,V)$prob=0;
  propset ExtTBR(Tau,V)$prob=0;
  propset NNI(Tau,V)$prob=0;
  """

    t = TEMPLATE.format(nexus=a.nexus_path, out_base=base, extra=extra, length=a.length,
                        samplefreq=a.length // 1000, printfreq = a.length // 100,
                        nruns=a.runs, nchains=a.chains,
                        diagnfreq=a.length / 2)
    with a.outfile as fp:
        fp.write(t)

if __name__ == '__main__':
    main()
