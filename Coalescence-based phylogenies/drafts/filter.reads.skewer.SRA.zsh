
#!/bin/zsh
############################################################################
# readfiltering with skewer
# Eckart Stolle, April 2020
############################################################################
# usage: ./skewer.zsh PATH/NAME.sra_1.fastq PATH/NAME.sra_2.fastq OUTPUTFOLDER CPUs
# no optical filtering (SRA downloaded files do not have coordinates stored

#set/get variables
INPUTR1=$1
INPUTR2=$2
OUTPUTFOLDER=$3
CPUs=$4

FWD_ADAPTERS="/scratch/genomes/illumina.fwd.adapters.fa"
REV_ADAPTERS="/scratch/genomes/illumina.rev.adapters.fa"
AVG_QUAL_THRESHOLD=20
END_QUAL_THRESHOLD=15
L=80

mkdir -p $OUTPUTFOLDER

TMPFOLDER=/backup/tmp/xxxxx
TMPFOLDER=/backup/tmp/$NAME
mkdir -p $TMPFOLDER

NAME=$(echo $INPUTR1 | rev  | cut -d"/" -f 1 | cut -d"." -f 3- | rev)
INPATH=$(echo $INPUTR1 | rev  | cut -d"/" -f 2- | rev)

#$INPATH/$NAME.sra_1.fastq
#$INPATH/$NAME.sra_2.fastq

skewer -m pe $INPATH/$NAME.sra_1.fastq $INPATH/$NAME.sra_2.fastq -x $FWD_ADAPTERS -y $REV_ADAPTERS --mean-quality $AVG_QUAL_THRESHOLD --end-quality $END_QUAL_THRESHOLD -l $L -n -r 0.1 -z --format auto -t $CPUs -o $OUTPUTFOLDER/$NAME

mv $OUTPUTFOLDER/$NAME'-trimmed-pair1'.fastq.gz $OUTPUTFOLDER/$NAME.R1.fq.gz

mv $OUTPUTFOLDER/$NAME'-trimmed-pair2'.fastq.gz $OUTPUTFOLDER/$NAME.R2.fq.gz

mv $OUTPUTFOLDER/$NAME'-trimmed'.log $OUTPUTFOLDER/$NAME.skewer.log

rm -rf $TMPFOLDER


