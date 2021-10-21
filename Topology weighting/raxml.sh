#!/bin/sh

module load anaconda3
. /share/apps/centos7/anaconda3/5.2.0/etc/profile.d/conda.sh
conda activate raxml090

# Region
chr_name=$1

# Window number
i=$2

echo "doing window $i"

# Run raxml
mkdir -p results/concatenated_trees/${chr_name}_$i

raxml-ng \
  --msa tmp/window/${chr_name}/alignments/window_${i}.fasta \
  --model GTR+G \
  --threads 2 \
  --tree pars{10},rand{10} \
  --prefix results/concatenated_trees/${chr_name}_$i/${chr_name}_$i \
  > results/concatenated_trees/${chr_name}_$i/${chr_name}_$i.raxml.log 2>&1
