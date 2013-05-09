#!/bin/zsh

set -e
set -u
set -x

MOVE_TYPES="guided
noguided"

P="1
5
10"

for p in $P; do
  for move_type in $MOVE_TYPES; do
    echo "$p $move_type"
    base=50tax_trim_t1.run1.${move_type}.p${p}
    ../../bin/compare_posterior_topologies.py ../../output/50taxon-01/50taxon-01.phyx_phyml_tree.txt ${base}.json | ../../bin/decorate_csv.py -o ${base}.comp.csv particle_factor=$p move_type=$move_type type=sts-online
  done
done

cp ../../output/50taxon-01/50taxon-01.run1.comp.csv .
