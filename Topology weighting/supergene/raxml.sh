#!/bin/sh

module load anaconda3
. /share/apps/centos7/anaconda3/5.2.0/etc/profile.d/conda.sh
conda activate raxml090

# Window number
i=$1

echo "doing window $i"

# Run raxml
mkdir -p results/concatenated_trees/$i

raxml-ng \
  --msa tmp/concatenated_alignments/window_${i}.fasta \
  --model GTR+G \
  --threads 2 \
  --tree pars{10},rand{10} \
  --prefix results/concatenated_trees/$i/window_$i \
  > results/concatenated_trees/$i/window_${i}.raxml.log 2>&1
