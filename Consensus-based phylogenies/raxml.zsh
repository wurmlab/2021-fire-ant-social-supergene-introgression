#!/bin/zsh
############################################################################
# filter SNPs
# Eckart Stolle, April.2020
############################################################################
## script: raxml.zsh
## scripts processes 1 file (1 alignment), outputs full ML search with BS's
## 50 random and 50 parsimony starting trees, 200 Bootstraps
## use parallel to run it simultaneous on several files
## uses 4 threads
## echo $GENE | parallel --no-notice -j 4 "echo {} && $HOME/scripts/raxml.zsh $F/merged/{}.fasta.allsamples.fa $F/modeltest/{}.modeltest.NONCODING-CODING.part.aic $OUTDIR/{}"

if [ $# -ne 3 ]; then
    echo $0: usage: ./SNP-genotype-filtering.zsh INPUT.raw.vcf.gz SAMPLEname OUTPUTFOLDER 
	echo "\nINPUT.fa: alignment file"
	echo "\nPARTITIONfile with models, eg from mdeltest-ng"
	echo "\nOutput Prefix"
    exit 1
fi

#set/get variables
ALIGNMENTINPUT=$1
MODELFILE=$2
PREFIX=$3
THREADS=2

#echo $GENE | parallel --no-notice -j 4 "echo {} && raxml-ng --all --threads $THREADS --data-type DNA --msa $ALIGNMENTINPUT --seed 2 --tree pars{50},rand{50} --redo --bs-trees  --bs-metric fbp,tbe --bs-cutoff 0.03 --model GTR+G --prefix $PREFIX"
#--bs-metric fbp,tbe --bs-cutoff 0.03
#--outgroup gem-1-bigB-m.{}
#--model {}/{}.annotation.partitionfile.CDS.txt

raxml-ng --all --threads $THREADS --data-type DNA --msa $ALIGNMENTINPUT --seed 2 --tree pars{50},rand{50} --redo --bs-trees 100 --bs-metric fbp,tbe --model $MODELFILE --prefix $PREFIX
echo $PREFIX
