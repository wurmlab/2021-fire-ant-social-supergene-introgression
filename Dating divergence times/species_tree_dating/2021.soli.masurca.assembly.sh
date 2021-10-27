#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j n
#$ -M e.stolle@leibniz-zfmk.de
#$ -m be
#$ -N ek_genomix_masurca_assembly

######################################################################
## qsub -pe smp 10 -q large.q,smalle.q,medium.q -v "SAMPLE=$SAMPLE" /home/estolle/scripts/2021.soli.masurca.assembly.sh

## $NSLOTS
## qalter
## qstat -f
## qhost
## qdel
## qalter -q large.q 25895
## qdel 25895
## qsub -q large.q,small.q,medium.q sleep.sh

######################################################################

## modules
module load seqtk/1.3
module load R/3.6.2
module load bamtools/2.5.1
module load bedtools/2.29.2
module load blast+/2.9.0
module load bwa-mem2/2.0pre2
module load htslib/1.10.2
module load sambamba/0.7.1
module load samtools/1.10
module load masurca/3.3.7

######################################################################
echo "Masurca v3.3.7 assembly of Illumina PE data, fireant sample "$SAMPLE
echo "beginning at "$(date)

CPUs=$NSLOTS
echo "running "$SAMPLE" with "$CPUs

cd /home/estolle/solimasurca

#SAMPLE="SRR9008173"
#SAMPLE="SRR9008142"
#SAMPLE="SRR9008228"
#SAMPLE="SRR9008253"
#SAMPLE="SRR9008168"
#SAMPLE="SRR9008232"
#SAMPLE="SRR9008217"
#SAMPLE="SRR9008215"
#SAMPLE="SRR9008150"
#SAMPLE="SRR9008158"
#SAMPLE="SRR9008133"
#SAMPLE="SRR9008200"

OUTPUTDIR="$SAMPLE"
mkdir -p $OUTPUTDIR
cd $OUTPUTDIR

RUN="$SAMPLE.masurca.3.3.7"
INPUTFWD="/home/estolle/solimasurca/$SAMPLE.R1.fq.gz"
INPUTREV="/home/estolle/solimasurca/$SAMPLE.R2.fq.gz"
INPUTNANOPORE="none"


mkdir -p $RUN
cd $RUN
cat $HOME/scripts/masurca.3.3.7.ILLUMINA.Solenopsis.cfg | sed s,xxxx,"$INPUTREV",g | sed s,yyyy,"$INPUTFWD",g | sed s,zzzz,"$INPUTNANOPORE",g | sed s,bbbb,"$CPUs",g > $RUN.cfg     
masurca $RUN.cfg
./assemble.sh

cp CA/final.genome.scf.fasta $RUN.fa
INPUTFASTA="$RUN.fa"
samtools faidx $INPUTFASTA
echo $SAMPLE" Masurca assembly (v3.3.7)" > $INPUTFASTA.stats
assembly-stats -s -t $INPUTFASTA >> $INPUTFASTA.stats

#genomesize
GZ=$(cat $INPUTFASTA.fai | awk '{ sum+=$2} END {print sum}')
GZ2=$(printf "%0.2f" $(bc -l <<< scale=6\;($GZ / 1000000))| sed "s/,/./g")

#number of N
NN=$(seqtk comp $INPUTFASTA | cut -f9 | awk '{ sum+=$1} END {print sum}')
NFRAC=$(printf "%0.2f" $(bc -l <<< scale=6\;(100-$NN*100/$GZ))| sed "s/,/./g")
NFRAC2=$(printf "%0.2f" $(bc -l <<< scale=6\;($NN*100/$GZ))| sed "s/,/./g")

echo "Genomesize $GZ ($GZ2 Mb)" >> $INPUTFASTA.stats
echo "Ns in fasta: "$NN >> $INPUTFASTA.stats
#echo "fraction of known bases (%): "$NFRAC >> $INPUTFASTA.stats
echo "fraction of Ns (%): "$NFRAC2 >> $INPUTFASTA.stats
