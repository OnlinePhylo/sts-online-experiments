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

trees = map(os.path.abspath, glob.glob('trees/*.nwk'))

env = SlurmEnvironment(ENV=os.environ.copy())
env.PrependENVPath('PATH', './bin')

# Builders
env['BUILDERS']['ConvertToNexus'] = Builder(action='seqmagick convert --alphabet dna --output-format nexus $SOURCE $TARGET', suffix='.nex', src_suffix='.fasta')
env['BUILDERS']['NexusToNewick'] = Builder(action='nexus_to_newick.py $SOURCE $TARGET -b 250', suffix='.nwk', src_suffix='.t')
env['BUILDERS']['MrBayesConf'] = Builder(action='generate_mb.py $SOURCE -o $TARGET', suffix='.mb', src_suffix='.nex')
env['BUILDERS']['StsTrees'] = Builder(action='sts_to_nexus.py < $SOURCE > $TARGET', suffix='.trees', src_suffix='.sts')
env['BUILDERS']['StsLikes'] = Builder(action='cut -f 1 $SOURCE > $TARGET', suffix='.txt', src_suffix='.sts')
# End builders

natural_extension = env.SConscript('src/SConscript')

nest = Nest()
w = SConsWrap(nest, dest_dir='output')

target_with_env = functools.partial(w.add_target_with_env, environment=env)

nest.add('tree', trees, label_func=stripext)
def get_n_taxa(c):
    m = re.search(r'(\d+)taxon-.*', c['tree'])
    assert m
    return [int(m.group(1))]
nest.add('n_taxa', get_n_taxa, create_dir=False)

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
      'number_of_sites=500 '
      'rate_distribution=Uniform '
      'model=JC69 ')[0]

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
               for i in (1, 2)]
    return env.SAlloc(targets, c['full_mrbayes_config'], 'mpirun mb $SOURCE', 6)

@target_with_env()
def full_mrbayes_newick(env, outdir, c):
    return [env.NexusToNewick(t)[0] for t in c['full_mrbayes_trees'][:2]]

@target_with_env()
def full_mrbayes_consensus(env, outdir, c):
    return [env.Command('$OUTDIR/' + stripext(str(t)) + '.sum.tre',
                        t,
                        'sumtrees.py -b 250 $SOURCE > $TARGET')[0] for t in c['full_mrbayes_trees'][:2]]

def trim_counts(c):
    n_taxa = c['n_taxa']
    half_n_taxa = n_taxa / 2

    return [i for i in (1, 2, 5, 10) if i <= half_n_taxa]

nest.add('trim_count', trim_counts)
nest.add('keep_count', lambda c: [c['n_taxa'] - c['trim_count']],
         create_dir=False)

nest.add('trim_base', lambda c: ['{n_taxa}tax_trim{trim_count}'.format(**c)],
         create_dir=False)

@target_with_env()
def trimmed_nexus(env, outdir, c):
    return env.Local('$OUTDIR/${trim_base}.nex',
            '$full_nexus',
            'seqmagick convert --head $keep_count --input-format '
            'nexus --alphabet dna $SOURCE $TARGET')[0]

@target_with_env()
def trimmed_mrbayes_config(env, outdir, c):
    return env.MrBayesConf(c['trimmed_nexus'])[0]

@target_with_env()
def trimmed_mrbayes_trees(env, outdir, c):
    targets = ['$OUTDIR/${{trim_base}}.run{0}.{1}'.format(i, j)
               for j in ('t', 'p')
               for i in (1, 2)]
    return env.SAlloc(targets,
                      ['$trimmed_mrbayes_config', '$trimmed_nexus'],
                      'mpirun mb $SOURCE', 6)

#nest.add('tree_moves', (0, 5, 10))
nest.add('tree_moves', [0])

@target_with_env()
def sts_online(env, outdir, c):
    j = joiner(outdir)
    return [env.Command(j(stripext(str(treefile)) + '.sts'),
                        ['$fasta', treefile],
                        'sts-online -p 5 -b 250 --tree-moves $tree_moves $SOURCES > $TARGET')[0]
            for treefile in c['trimmed_mrbayes_trees'][:2]]

@target_with_env()
def sts_online_trees(env, outdir, c):
    return [env.StsTrees(i) for i in c['sts_online']]

@target_with_env()
def sts_consensus(env, outdir, c):
    return [env.Command('$OUTDIR/' + stripext(str(t)) + '.sum.tre',
                        t,
                        'sumtrees.py --weighted-trees $SOURCE > $TARGET')[0] for t in c['sts_online_trees']]

@target_with_env()
def consensus_comparison(env, outdir, c):
    sources = [c['tree']] + list(c['sts_consensus']) + c['full_mrbayes_consensus']
    return env.Local('$OUTDIR/consensus_to_source.csv',
            sources,
            'compare_to_source.py $SOURCES -o $TARGET --schema nexus')

#@target_with_env()
#def lnl_comparison(env, outdir, c):
    #return env.Local('$OUTDIR/${trim_base}_lnl_comp.pdf',
                     #env.Flatten([c['sts_online'], c['full_mrbayes_trees'][2:]]),
                     #'lnl_compare.R $SOURCES $TARGET')[0]


#@target_with_env()
#def compare_trees(env, outdir, c):
    #j = joiner(outdir)
    #all_trees = list(c['sts_online_trees']) + c['full_mrbayes_newick']
    #env.Command(j('comparisons.csv'),
                #all_trees,
                #'bin/trees_compare.py $SOURCES -o $TARGET')

#nest.add('type', ('full', 'trimmed'))

#@w.add_target()
#def mrbayes_trees(outdir, c):
    #j = joiner(outdir)
    #targets = [j(stripext(str(c['nexus'])) + '.run{0}.t'.format(i)) for i in (1,2)]
    #return env.SAlloc(targets, c['mrbayes_config'],
            #'mpirun mb $SOURCE', 6)

#@w.add_target()
#def newick_trees(outdir, c):
    #if c['type'] == 'full':
        #r = [env.NexusToNewick(tree) for tree in c['mrbayes_trees']]
    #else:
        #j = joiner(outdir)
        #r = []
        #for tree in c['mrbayes_trees']:
            #i, = env.Command(j(stripext(str(tree)) + '_online.nwk'),
                    #[c['fasta'], tree],
                    #'bin/sts-online $SOURCES -b 250 -p 3 -m 6 | cut -f 2 > $TARGET')
            #r.append(i)
    #c['compare_trees'][c['type']] = r
    #return r

#@w.add_target()
#def natural_extension_result(outdir, c):
    #if c['type'] != 'full':
        #return None

    #tree = c['mrbayes_trees'][0]
    #res, = env.Command(joiner(outdir)(stripext(str(tree)) + '_natext.csv'),
        #[natural_extension, c['fasta'], tree],
        #'${SOURCES[0]} '
        #'input.sequence.file=${SOURCES[1]} '
        #'natural_extension.trees_nexus=${SOURCES[2]} '
        #'natural_extension.burnin=250 '
        #'natural_extension.prune_taxon=t1 '
        #'natural_extension.output_path=$TARGET')
    #env.Depends(res, 'src/natural_extension')

w.finalize_all_aggregates()
