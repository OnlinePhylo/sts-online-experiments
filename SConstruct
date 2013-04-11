#!/usr/bin/env scons -Q

import glob
import os
import os.path
import functools

from nestly import Nest, stripext
from nestly.scons import SConsWrap
from bioscons.slurm import SlurmEnvironment
from SCons.Script import Builder

def joiner(*args):
    return functools.partial(os.path.join, *args)

trees = map(os.path.abspath, glob.glob('trees/*.nwk'))


env = SlurmEnvironment(ENV=os.environ.copy())

# Builders
env['BUILDERS']['ConvertToNexus'] = Builder(action='seqmagick convert --alphabet dna --output-format nexus $SOURCE $TARGET', suffix='.nex', src_suffix='.fasta')
env['BUILDERS']['NexusToNewick'] = Builder(action='nexus_to_newick.py $SOURCE $TARGET -b 250', suffix='.nwk', src_suffix='.t')
env['BUILDERS']['MrBayesConf'] = Builder(action='bin/generate_mb.py $SOURCE -o $TARGET', suffix='.mb', src_suffix='.nex')
# End builders

natural_extension = env.SConscript('src/SConscript')

nest = Nest()
w = SConsWrap(nest, 'output')

nest.add('tree', trees, label_func=stripext)

@w.add_target()
def fasta(outdir, c):
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

@w.add_target()
def full_nexus(outdir, c):
    j = joiner(outdir)
    fasta = str(c['fasta'])
    return env.ConvertToNexus(j(stripext(fasta) + '.nex'), fasta)[0]

@w.add_target()
def full_mrbayes_config(outdir, c):
    return env.MrBayesConf(c['full_nexus'])[0]

@w.add_target()
def full_mrbayes_trees(outdir, c):
    j = joiner(outdir)
    targets = [j(stripext(str(c['full_nexus'])) + '.run{0}.t'.format(i)) for i in (1,2)]
    return env.SAlloc(targets, c['full_mrbayes_config'], 'mpirun mb $SOURCE', 6)

#@w.add_aggregate(dict)
#def compare_trees(outdir, c, inputs):
    #j = joiner(outdir)
    #all_trees = [tree for v in inputs.values() for tree in v]
    #env.Command(j('comparisons.csv'),
                #all_trees,
                #'bin/trees_compare.py $SOURCES -o $TARGET')

#nest.add('type', ('full', 'trimmed'))

#@w.add_target()
#def nexus(outdir, c):
    #j = joiner(outdir)
    #if c['type'] == 'full':
        #fasta = str(c['fasta'])
        #return env.ConvertToNexus(j(stripext(fasta) + '.nex'), fasta)[0]
    #else:
        #target = j(stripext(c['tree']) + '_trimmed.nex')
        #return env.Local(target, c['fasta'],
                #'seqmagick convert $SOURCE $TARGET --alphabet dna --output-format nexus --pattern-exclude="^t1\\$"')[0]

#@w.add_target()
#def mrbayes_config(outdir, c):
    #return env.MrBayesConf(c['nexus'])[0]

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

#w.finalize_all_aggregates()
