#!/usr/bin/env scons -Q

from __future__ import division

import glob
import os
import os.path
import re
import functools

from nestly import Nest, stripext
from nestly.scons import SConsWrap
from bioscons.slurm import SlurmEnvironment
from SCons.Script import Builder

def joiner(*args):
    return functools.partial(os.path.join, *args)

trees = sorted(map(os.path.abspath, glob.glob('trees/50*.nwk')))

env = SlurmEnvironment(ENV=os.environ.copy())
env.PrependENVPath('PATH', './bin')
env['outdir'] = 'output'

env['MB_NRUNS'] = 4

# Builders
env['BUILDERS']['ConvertToNexus'] = Builder(action='seqmagick convert --alphabet dna --output-format nexus $SOURCE $TARGET', suffix='.nex', src_suffix='.fasta')
env['BUILDERS']['ConvertToPhyx'] = Builder(action='seqmagick convert $SOURCE $TARGET', suffix='.phyx', src_suffix='.fasta')
env['BUILDERS']['NexusToNewick'] = Builder(action='nexus_to_newick.py $SOURCE $TARGET -b 250', suffix='.nwk', src_suffix='.t')
env['BUILDERS']['MrBayesConf'] = Builder(action='generate_mb.py --runs $MB_NRUNS --chains 1 --length 10000000 $SOURCE -o $TARGET', suffix='.mb', src_suffix='.nex')
env['BUILDERS']['StsTrees'] = Builder(action='sts_to_nexus.py -i $SOURCE -o $TARGET', suffix='.trees', src_suffix='.sts')
env['BUILDERS']['StsLikes'] = Builder(action='cut -f 1 $SOURCE > $TARGET', suffix='.txt', src_suffix='.sts')
# End builders

natural_extension = env.SConscript('src/SConscript')

nest = Nest()
w = SConsWrap(nest, dest_dir='output')

target_with_env = functools.partial(w.add_target_with_env, environment=env)

# Aggregate
posterior_comparisons = []
posterior_comparisons_keys = ('tree', 'n_taxa', 'trim_taxon', 'particle_factor', 'tree_moves')

nest.add('tree', trees, label_func=stripext)


def get_n_taxa(c):
    m = re.search(r'(\d+)taxon-.*', c['tree'])
    assert m
    return [int(m.group(1))]
nest.add('n_taxa', get_n_taxa, create_dir=False)

@w.add_aggregate(list)
def pendant_bl(outdir, c, inputs):
    env.Local(os.path.join(outdir, 'pendant_bl_ess.csv'),
            env.Flatten([c['tree'], inputs]),
            'pendant_bl.py $SOURCES -o $TARGET')

@target_with_env()
def fasta(env, outdir, c):
    j = joiner(outdir)
    target = j(stripext(c['tree']) + '.fasta')
    return env.Local(target, c['tree'],
      'bppseqgen '
      'input.tree.file=$SOURCE '
      'input.tree.format=Newick '
      'output.sequence.file=$TARGET '
      'output.sequence.format=Fasta '
      'alphabet=DNA '
      'number_of_sites=1000 '
      'rate_distribution=Uniform '
      'model=JC69 ')[0]

@target_with_env()
def full_phylip(env, outdir, c):
    return env.ConvertToPhyx(c['fasta'])[0]

@target_with_env()
def phyml_tree(env, outdir, c):
    return env.Command('${full_phylip}_phyml_tree.txt',
            ['$full_phylip', '$tree'],
            'phyml -i ${SOURCES[0]} -u ${SOURCES[1]} -c 1 -m JC69 -o l -b 0')[0]

@target_with_env()
def full_nexus(env, outdir, c):
    j = joiner(outdir)
    fasta = str(c['fasta'])
    return env.ConvertToNexus(j(stripext(fasta) + '.nex'), fasta)[0]

@target_with_env()
def full_mrbayes_config(env, outdir, c):
    return env.MrBayesConf(c['full_nexus'])[0]

@target_with_env()
def full_mrbayes_trees(env, outdir, c):
    j = joiner(outdir)
    targets = [j('{0}.run{1}.{2}'.format(stripext(str(c['full_nexus'])), i, t))
               for t in ('t', 'p')
               for i in xrange(1, env['MB_NRUNS'] + 1)]
    res = env.SAlloc(targets, c['full_mrbayes_config'], 'mpirun mb $SOURCE', 4)
    res = {'trees': res[:env['MB_NRUNS']],
           'params': res[env['MB_NRUNS']:]}
    j = joiner(outdir)
    env['kvs'] = ' '.join('{0}="{1}"'.format(k, c.get(k, '')) for k in posterior_comparisons_keys)

    pc = [env.Local(j(stripext(str(t)) + '.comp.csv'),
                    [c['phyml_tree'], t],
                    'compare_posterior_topologies.py -b 250 $SOURCES | decorate_csv.py $kvs type=MrBayes > $TARGET')
            for t in res['trees']]
    posterior_comparisons.extend(pc)

    return res

@target_with_env()
def full_mrbayes_newick(env, outdir, c):
    return [env.NexusToNewick(t)[0] for t in c['full_mrbayes_trees']['trees']]

@target_with_env()
def full_mrbayes_consensus(env, outdir, c):
    return [env.Command('$OUTDIR/' + stripext(str(t)) + '.sum.tre',
                        t,
                        'sumtrees.py -q -b 250 $SOURCE > $TARGET')[0] for t in c['full_mrbayes_trees']['trees']]

nest.add('trim_replicates', [5], create_dir=False)
nest.add('trim_count', [1], create_dir=False)

def trim_taxon(c):
    assert c['trim_count'] == 1
    all_taxa = ['t{0}'.format(i+1) for i in xrange(c['n_taxa'])]
    return all_taxa[:c['trim_replicates']]

nest.add('keep_count', lambda c: [c['n_taxa'] - c['trim_count']],
         create_dir=False)

nest.add('trim_taxon', trim_taxon)

nest.add('trim_base', lambda c: ['{n_taxa:02d}tax_trim_$trim_taxon'.format(**c)],
         create_dir=False)

@target_with_env()
def trimmed_nexus(env, outdir, c):
    return env.Local('$OUTDIR/${trim_base}.nex',
            '$full_nexus',
            "seqmagick convert --pattern-exclude '^($trim_taxon)$'  --input-format "
            'nexus --alphabet dna $SOURCE $TARGET')[0]

@target_with_env()
def trimmed_mrbayes_config(env, outdir, c):
    return env.MrBayesConf(c['trimmed_nexus'])[0]

@target_with_env()
def trimmed_mrbayes_trees(env, outdir, c):
    targets = ['$OUTDIR/${{trim_base}}.run{0}.{1}'.format(i, j)
               for j in ('t', 'p')
               for i in xrange(1, env['MB_NRUNS'] + 1)]
    res = env.SAlloc(targets,
                      ['$trimmed_mrbayes_config', '$trimmed_nexus'],
                      'mpirun mb $SOURCE', 4)
    return {'trees': res[:env['MB_NRUNS']],
            'params': res[env['MB_NRUNS']:]}

#nest.add('tree_moves', (0, 5, 10))
nest.add('tree_moves', [0])

# nest.add('particle_factor', [1, 5, 10])
nest.add('particle_factor', [1])

@target_with_env()
def sts_online(env, outdir, c):
    j = joiner(outdir)
    result = [env.Command(j(stripext(str(treefile)) + '.sts.json'),
                          ['$fasta', treefile],
                          'sts-online -p $particle_factor -b 250 --tree-moves $tree_moves $SOURCES $TARGET')[0]
            for treefile in c['trimmed_mrbayes_trees']['trees']]
    c['pendant_bl'].extend(result)
    return result

#@target_with_env()
#def sts_online_trees(env, outdir, c):
    #return [env.StsTrees(i) for i in c['sts_online']]

#@target_with_env()
#def sts_consensus(env, outdir, c):
    #return [env.Command('$OUTDIR/' + stripext(str(t)) + '.sum.tre',
                        #t,
                        #'sumtrees.py -q --weighted-trees $SOURCE > $TARGET')[0] for t in c['sts_online_trees']]
    ##return [env.Command('$OUTDIR/' + stripext(str(t)) + '.sum.tre',
                        ##t,
                        ##'sumtrees.py $SOURCE > $TARGET')[0] for t in c['sts_online_trees']]

@target_with_env()
def sts_posterior_comparison(env, outdir, c):
    j = joiner(outdir)
    env['kvs'] = ' '.join('{0}="{1}"'.format(k, c.get(k, '')) for k in posterior_comparisons_keys)

    res = [env.Local(j(stripext(str(t)) + '.comp.csv'),
                     [c['phyml_tree'], t],
                     'compare_posterior_topologies.py $SOURCES | decorate_csv.py $kvs type=sts-online > $TARGET')
            for t in c['sts_online']]
    posterior_comparisons.extend(res)
    return res

#@target_with_env()
#def consensus_comparison(env, outdir, c):
    #sources = [c['phyml_tree']] + list(c['sts_consensus']) + c['full_mrbayes_consensus']
    #return env.Local('$OUTDIR/consensus_to_source.csv',
            #sources,
            #'compare_to_source.py $SOURCES -o $TARGET --schema nexus')

all_posterior_comparison = env.Local('$outdir/posterior_comparison.csv',
        posterior_comparisons, 'csvstack $SOURCES > $TARGET')
all_posterior_plot = env.Local('$outdir/posterior_comparison.svg',
        all_posterior_comparison, 'plot_posterior_rf.R $SOURCE $TARGET')
#env.Precious(compare_to_source)
#env.Local(['$outdir/compare_to_source.svg', '$outdir/compare_to_source_rf.svg'],
        #compare_to_source,
        #'plot_cons.R $SOURCE $TARGETS')
w.finalize_all_aggregates()
