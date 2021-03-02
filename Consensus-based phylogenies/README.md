# Coalescence-based tree reconstruction

## Introduction

Our aim was to construct a tree representing the supergene region and a tree representing the species phylogeny (from chromosomes 1 to 15).

We built the tree with sequences from single-copy genes (including exons, introns and 1000 bp upstream and downstream of each gene). For each single-copy gene we created a fasta file with the aligned consensus nucleotide sequences for all the individuals, and a partition file with the boundary coordinates for coding and non-coding sequences.

All the processes in this script are detailed in `phylogeny.txt`.

## Selecting the best-fit model of evolution for each gene

We used `modeltest-ng` to determine the best-fit model for each partition in each gene.

```sh

 modeltest-ng \
   -i gene.fasta \
   -d nt \
   -p $CPUs \
   -q gene.partitionfile.txt \
   -T raxml \
   --rngseed 2 \
   -t fixed-ml-gtr \
   --force \
   --output gene.modeltest

```

## Phylogeny for each tree

We used RAXML to build a phylogeny for each single-copy gene. This was done with the script `raxml.zsh`, using the best-fit model determined in the previous step.

``` sh

raxml-ng --all --threads $THREADS \
  --data-type DNA \
  --msa gene.fasta \
  --seed 2 \
  --tree pars{50},rand{50} \
  --redo \
  --bs-trees 100 \
  --bs-metric fbp,tbe \
  --model gene.modeltest.part.aic \
  --prefix gene

```

## Consensus tree reconstructions

We used ASTRAL to build a consensus tree for the supergene region and the chromosomes 1 to 15.

```sh

# Supergene tree
INPUTTREES="astral.chr16r/input.10SNPs.tre"
OUTPUTFILE="astral.chr16r/output.10SNPs.tre"
java -Xmx250g \
-D"java.library.path=/usr/local/src/Astral/lib/" \
-jar /usr/local/src/Astral/astral.5.14.3.jar \
--input $INPUTTREES \
--output $OUTPUTFILE \
--cpu-only \
--cpu-threads 50 \
--branch-annotate 3 \
--keep completed

# Chr1-15 tree
INPUTTREES="astral.chr1-15/input.10SNPs.tre"
OUTPUTFILE="astral.chr1-15/output.10SNPs.tre"
java -Xmx250g \
-D"java.library.path=/usr/local/src/Astral/lib/" \
-jar /usr/local/src/Astral/astral.5.14.3.jar \
--input $INPUTTREES \
--output $OUTPUTFILE \
--cpu-only \
--cpu-threads 50 \
--branch-annotate 3 \
--keep completed

```
