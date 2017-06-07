#!/usr/bin/env python

import itertools
import json
import csv
import os
import sys
import subprocess
from subprocess import check_call
from timeit import default_timer as timer
import argparse
import logging
from shutil import copyfile
import re

log = logging.getLogger('sts')

p = argparse.ArgumentParser()
p.add_argument('-o', '--output', default='output')
arg = p.parse_args()

logging.basicConfig(level=logging.INFO)

# directory where bin and trees folders are located
dir = '.'

siteCount = 1000

# MrBayes
burnin = 250
generations = 300000
runCount = 2

# Simulation
taxonCounts = (10, 50, 100)
treeReplicates = xrange(1, 6)
trimCounts = 1, 2, 5)
trimReplicates = 3
particleFactors = (1, 5, 10, 50, 100)

methods = {
    'lcfit': ['--proposal-method', 'lcfit'],
    'uniform-length': ['--proposal-method', 'uniform-length'],
    'guided-parsimony': ['--proposal-method', 'guided-parsimony'],
    'guided': ['--proposal-method', 'guided', '-e', '1'],
    'lcfit2': ['--proposal-method', 'lcfit', '-e', '1']
}

bin = os.path.join(dir, 'bin')

#sts = os.path.join(bin,'sts-online')

# scripts
cmp_pt = os.path.join(bin, 'compare_posterior_topologies.py')
decorate_csv = os.path.join(bin, 'decorate_csv.py')
asdsf = os.path.join(bin, 'asdsf.py')
extract_ess_like = os.path.join(bin, 'extract_ess_like.py')
extract_ess = os.path.join(bin, 'extract_ess.py')
generate_mb = os.path.join(bin, 'generate_mb.py')

# output files
ess_calls_csv = os.path.join(arg.output, 'ess_calls.csv')
ess_csv = os.path.join(arg.output, 'ess.csv')
posterior_comparison_csv = os.path.join(arg.output, 'posterior_comparison.csv')
pp_comparison_csv = os.path.join(arg.output, 'pp_comparison.csv')
mrbayes_csv = os.path.join(arg.output, 'mrbayes.csv')

if not os.path.lexists(arg.output):
    os.mkdir(arg.output)

asdsfs = ['csvstack']
pps = ['csvstack']
ppcomps = ['csvstack']
probs = ['csvstack']
controls = []

MRBAYES_TEMPLATE = """
begin mrbayes;
    set autoclose=yes nowarn=yes;
    execute {nexus};
    lset nst=1 rates=equal;
    prset statefreqpr=fixed(equal);
    prset brlenspr = unconstrained:exponential(10.0); 
    {extra}
    mcmcp nruns={nruns} nchains={nchains} ngen={length} samplefreq={samplefreq} printfreq={printfreq} file={out_base} diagnfreq={printfreq} stopVal=0.01 Stoprule=yes;
    mcmc;
end;
"""


def main():
    mb_csv = open(mrbayes_csv, "w")

    for n in taxonCounts:
        for i in treeReplicates:

            nwk = os.path.join(dir, 'trees', '{}taxon-0{}.nwk'.format(n, i))

            tree_count_dir = os.path.join(arg.output, '{}taxon-0{}'.format(n, i))

            if not os.path.lexists(tree_count_dir):
                os.mkdir(tree_count_dir)

            stem = os.path.join(tree_count_dir, '{}taxon-0{}'.format(n, i))

            mb = stem + '.mb'
            mb_t1 = stem + '.run1.t'
            mb_t2 = stem + '.run2.t'
            mb_p1 = stem + '.run1.p'
            mb_p2 = stem + '.run2.p'
            phyx = stem + '.phyx'
            phyml_t = phyx + '_phyml_tree'
            fasta = stem + '.fasta'
            nex = stem + '.nex'

            check_call(['bppseqgen', '--seed', '0', 'input.tree.file='+nwk, 'input.tree.format=Newick',
                        'output.sequence.file='+fasta, 'output.sequence.format=Fasta', 'alphabet=DNA',
                        'number_of_sites={}'.format(siteCount), 'rate_distribution=Constant', 'model=JC69'])

            check_call(['seqmagick', 'convert', '--alphabet', 'dna', '--output-format', 'nexus', fasta, nex])

            check_call([generate_mb, '--runs', str(runCount), '--length', str(generations),
                        os.path.join(arg.output, '{}taxon-0{}'.format(n, i), '{}taxon-0{}.nex'.format(n, i)),
                        '-o', mb])

            check_call(['seqmagick', 'convert', fasta, phyx])
            check_call(['phyml', '-i', phyx, '-u', nwk, '-c', '1', '-m', 'JC69', '-o', 'l', '-b', '0'])
            check_call(['mb', os.path.join(arg.output, '{}taxon-0{}'.format(n, i), '{}taxon-0{}.mb'.format(n, i))])
            check_call(['sed', '-i', '-e', 's/e+00//g', mb_t1, mb_t2, mb_p1, mb_p2])

            # compare posterior
            cmd = [decorate_csv, 'tree={}'.format(nwk), 'n_taxa={}'.format(n), 'trim_taxon=""', 'trim_count=""',
                   'particle_factor=""', 'proposal_method_name=""', 'proposal_args=""', 'type=MrBayes']

            for j in xrange(1, runCount + 1):
                run_t = stem + '.run{}.t'.format(j)
                p_cmp = subprocess.Popen([cmp_pt, '-b', str(burnin), phyml_t, run_t], stdout=subprocess.PIPE)

                with open(stem + '.run{}.comp.csv'.format(j), 'w') as out:
                    a = subprocess.Popen(cmd, stdin=p_cmp.stdout, stdout=out)
                    a.wait()
                p_cmp.wait()

            asdsf_csv = os.path.join(tree_count_dir, 'asdsf.csv')
            pp_csv = os.path.join(tree_count_dir, 'pp.csv')
            pp_annot_csv = os.path.join(tree_count_dir, 'pp_annot.csv')

            asdsfs.append(asdsf_csv)
            pps.append(pp_annot_csv)

            p_cmp = subprocess.Popen(['python', asdsf, mb_t1, mb_t2, '--pp-table', pp_csv],
                                     stdout=subprocess.PIPE)

            cmd = ['python', decorate_csv, 'tree={}'.format(nwk), 'n_taxa={}'.format(n), 'trim_taxon=0', 'trim_count=0',
                   'particle_factor=0', 'proposal_method_name=0', 'proposal_args=0']

            # create asdsf.csv from the standard output of p_cmp
            with open(asdsf_csv, 'w') as out:
                a = subprocess.Popen(cmd, stdin=p_cmp.stdout, stdout=out)
                a.wait()
            p_cmp.wait()

            # create pp_annot.csv from pp.csv
            cmd = ['python', decorate_csv, 'tree={}'.format(nwk), 'n_taxa={}'.format(n), 'trim_taxon=0', 'trim_count=0',
                   'particle_factor=0', 'proposal_method_name=0', 'proposal_args=0', '-i', pp_csv, '-o', pp_annot_csv]
            check_call(cmd)

            for trim_count in trimCounts:
                combinations = ('-'.join(c) for c in itertools.combinations(['t{}'.format(nn+1) for nn in range(n)], trim_count))
                for c in list(itertools.islice(combinations, 0, trimReplicates)):

                    trim_dir = os.path.join(tree_count_dir, c)

                    if not os.path.lexists(trim_dir):
                        os.mkdir(trim_dir)

                    stem_trim = os.path.join(trim_dir, '{}tax_trim_{}'.format(n, c))

                    trim_nex = stem_trim + '.nex'
                    trim_mb = stem_trim + '.mb'
                    trim_t = stem_trim + '.t'
                    trim_p = stem_trim + '.p'

                    check_call(['seqmagick', 'convert', '--pattern-exclude', '^({})$'.format(c).replace('-', '|'),
                                '--input-format', 'nexus' ,'--alphabet', 'dna', nex, trim_nex])
                    check_call([generate_mb, '--runs', '1', '--length', str(generations),
                                os.path.join(arg.output, '{}taxon-0{}'.format(n, i), c, '{}tax_trim_{}.nex'.format(n, c)), '-o', trim_mb])
                    check_call(['mb', os.path.join(arg.output, '{}taxon-0{}'.format(n, i), c, '{}tax_trim_{}.mb'.format(n, c))])
                    check_call(['sed', '-i', '-e', 's/e+00//g', trim_t, trim_p])

                    for method, particle_factor in itertools.product(methods.keys(), particleFactors):
                        if not os.path.lexists(os.path.join(trim_dir, str(particle_factor))):
                            os.mkdir(os.path.join(trim_dir, str(particle_factor)))

                        proposal_dir = os.path.join(trim_dir, str(particle_factor), method)

                        if not os.path.lexists(proposal_dir):
                            os.mkdir(proposal_dir)

                        jsonf = os.path.join(proposal_dir, '{}tax_trim_{}.sts.json'.format(n, c))
                        control = os.path.join(proposal_dir, 'control.json')

                        cmp_csv = os.path.join(proposal_dir, '{}tax_trim_{}.sts.comp.csv'.format(n, c))
                        asdsf_csv = os.path.join(proposal_dir, 'asdsf.csv')
                        pp_csv = os.path.join(proposal_dir, 'pp.csv')
                        pp_annot_csv = os.path.join(proposal_dir, 'pp_annot.csv')
                        probs_csv = os.path.join(proposal_dir, 'probs.csv')

                        asdsfs.append(asdsf_csv)
                        pps.append(pp_annot_csv)
                        ppcomps.append(cmp_csv)
                        controls.append(control)
                        probs.append(probs_csv)

                        cmd = [sts] + methods[method]

                        cmd.extend(['-p', str(particle_factor), '-b', str(burnin), fasta, trim_t, jsonf])
                        print(' '.join(cmd))

                        start = timer()
                        try:
                            check_call(cmd)
                        except subprocess.CalledProcessError, e:
                            sys.stderr.write('FAILED\n')
                            continue
                        end = timer()
                        total_time = end - start
                        print('Time: {}'.format(total_time))

                        write_control(control, method, trim_count, c, n - trim_count,
                                      particle_factor,
                                      n, nwk, [jsonf], total_time)

                        loggPs = []
                        ids = []
                        with open(jsonf) as js:
                            root = json.load(js)
                            for proposals in root['proposals']:
                                if proposals['T'] == trim_count:
                                    loggPs.append(proposals['newLogLike'])

                            for trees in root['trees']:
                                ids.append(trees['particleID'])

                        with open(probs_csv, 'w') as ll:
                            keys = ['trim_taxon', 'n_taxa', 'trim_count', 'particle_factor', 'proposal_method_name']
                            header = [str(c), str(n), str(trim_count), str(particle_factor), method]
                            row_base = {keys[i]: header[i] for i in range(len(header))}

                            w = csv.DictWriter(ll, keys + ['logP'])
                            w.writeheader()
                            for idx in ids:
                                row = row_base.copy()
                                row['logP'] = loggPs[idx]
                                w.writerow(row)

                        # compare posterior
                        cmd = ['python', decorate_csv, 'tree={}'.format(nwk), 'n_taxa={}'.format(n),
                               'trim_taxon={}'.format(c), 'trim_count={}'.format(trim_count),
                               'particle_factor={}'.format(particle_factor),
                               'proposal_method_name={}'.format(method),
                               'proposal_args=' + ' '.join(methods[method]), 'type=sts-online']

                        # create .sts.comp.csv
                        p_cmp = subprocess.Popen(['python', cmp_pt, phyml_t, jsonf],
                                                 stdout=subprocess.PIPE)
                        with open(cmp_csv, 'w') as out:
                            a = subprocess.Popen(cmd, stdin=p_cmp.stdout, stdout=out)
                            a.wait()
                        p_cmp.wait()

                        # Create pp.csv
                        p_cmp = subprocess.Popen(['python', asdsf, mb_t1, jsonf, '--pp-table', pp_csv],
                                                 stdout=subprocess.PIPE)

                        cmd = ['python', decorate_csv, 'tree={}'.format(nwk), 'n_taxa={}'.format(n),
                               'trim_taxon={}'.format(c),
                               'trim_count={}'.format(trim_count),
                               'particle_factor={}'.format(particle_factor),
                               'proposal_method_name={}'.format(method),
                               'proposal_args=' + ' '.join(methods[method])]

                        # create asdsf.csv from the standard output of p_cmp
                        with open(asdsf_csv, 'w') as out:
                            a = subprocess.Popen(cmd, stdin=p_cmp.stdout, stdout=out)
                            a.wait()
                        p_cmp.wait()

                        # create pp_annot.csv from pp.csv
                        cmd = ['python', decorate_csv, 'tree={}'.format(nwk), 'n_taxa={}'.format(n),
                               'trim_taxon={}'.format(c),
                               'trim_count={}'.format(trim_count),
                               'particle_factor={}'.format(particle_factor),
                               'proposal_method_name={}'.format(method),
                               'proposal_args=' + ' '.join(methods[method]), '-i', pp_csv, '-o', pp_annot_csv]
                        check_call(cmd)

                    if trim_count == 5 and n > 10:
                        order = []
                        with open(jsonf) as js:
                            root = json.load(js)
                            for proposals in root['generations']:
                                order.append(proposals['sequence'])

                        length = 30000000
                        base = '{}taxon-0{}'.format(n, i)
                        nexus = base + '.nex'

                        for cc in xrange(1,trim_count+1):
                            seqCount = n+cc-trim_count
                            tree_count_dir_part = os.path.join(tree_count_dir, str(seqCount))
                            if not os.path.lexists(tree_count_dir_part):
                                os.mkdir(tree_count_dir_part)
                            trim_nex = os.path.join(tree_count_dir_part, '{}taxon-0{}.nex'.format(n, i))
                            trim_mb = os.path.join(tree_count_dir_part, '{}taxon-0{}.mb'.format(n, i))

                            if cc == trim_count:
                                copyfile(nex, trim_nex)
                            else:
                                check_call(['seqmagick', 'convert', '--pattern-exclude',
                                            '^({})$'.format('|'.join(order[cc:])), '--input-format',
                                            'nexus', '--alphabet', 'dna', nex, trim_nex])

                            t = MRBAYES_TEMPLATE.format(nexus=nexus, out_base=base, extra='', length=length,
                                                samplefreq=length // 1000, printfreq=length // 10000,
                                                nruns=2, nchains=1, diagnfreq=length / 2)
                            print(trim_mb)
                            with open(trim_mb,'w') as fp:
                                fp.write(t)

                            start = timer()
                            check_call(['mb', base+'.mb'], cwd=tree_count_dir_part)
                            end = timer()
                            total_time = end - start
                            print('Time: {}'.format(total_time))

                            asdsfv = -1
                            pattern = re.compile(r'.+(\d+\.\d+)$')
                            with open(os.path.join(tree_count_dir_part, '{}taxon-0{}.mcmc'.format(n, i))) as fp:
                                for line in fp:
                                    line = line.rstrip('\n').rstrip('\r')
                                    m = pattern.match(line)
                                    if m:
                                        asdsfv = m.group(1)
                            mb_csv.write('{},{},{},{},{},{}\n'.format(i, n, c, cc, total_time, asdsfv))

    mb_csv.close()

    log.info('Create control.json')
    check_call(['python', extract_ess_like, '-o', ess_calls_csv]+controls, cwd=dir)
    check_call(['python', extract_ess, '-o', ess_csv]+controls, cwd=dir)

    log.info('Create comp.csv')
    with open(posterior_comparison_csv, 'w') as out:
        a = subprocess.Popen(ppcomps, stdout=out)
        a.wait()

    log.info('Create pp_annot.csv')
    with open(os.path.join(arg.output, 'pp_comparison.csv'), 'w') as out:
        a = subprocess.Popen(pps, stdout=out)
        a.wait()

    log.info('Create asdsf.csv')
    with open(os.path.join(arg.output, 'asdsf.csv'), 'w') as out:
        a = subprocess.Popen(asdsfs, stdout=out)
        a.wait()

    log.info('Create probs.csv')
    with open(os.path.join(arg.output, 'probs.csv'), 'w') as out:
        a = subprocess.Popen(probs, stdout=out)
        a.wait()


def write_control(filename, method, trim_count, trim_taxon, keep_count, particle_factor, n_taxa, tree, files, tt):
    data = {'proposal_method_name': method,
            'trim_count': trim_count,
            'keep_count': keep_count,
            'trim_taxon': trim_taxon,
            'particle_factor': particle_factor,
            'n_taxa': n_taxa,
            'tree': tree,
            'sts_online': files,
            'time': tt}

    with open(filename, "w") as f:
        json.dump(data, f)


if __name__ == '__main__':
    main()
