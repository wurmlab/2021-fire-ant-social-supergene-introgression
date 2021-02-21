#!/bin/zsh
############################################################################
# readfiltering with skewer
# Eckart Stolle, Sept 2017
############################################################################
# uses zsh shell (not bash), get it by: sudo apt-get install zsh
# usage: ./filter.skewer.Soli18.zsh INPUTFOLDER OUTPUTFOLDER CPUs

# !!! comment/uncomment parts if sickle is used for further/stricter trimming in addition  (trimming/stats/plot section)
# !!! see below and update location of fwd/rev adapters to filter
# !!! see below to adapt filter_QCplot.R file location
# !!! see below to adapt the format of the fwd/rev reads file formats (.R1 or _1)
# !!! adapt optical distance for optical-duplicate removal (depending on Platform)
# !!! L=80 length and Qual thresholds

FWD_ADAPTERS="$HOME/scripts/illumina.fwd.adapters.fa"
REV_ADAPTERS="$HOME/scripts/illumina.rev.adapters.fa"
AVG_QUAL_THRESHOLD=20
END_QUAL_THRESHOLD=15
L=100
OPTICALDIST=2500
##dupedist instead dist HiSeq2000:40 HiSeq4000:2500
#PLOTSCRIPT="$HOME/scripts/filter_QCplot.R"
#--mean-quality 20 --end-quality 15 -l 100 -n yes -r 0.1 -z

#run:
#~/scripts/filter.skewer.Soli18.zsh $INPUTFOLDER $OUTPUTFOLDER $CPUs


if [ $# -ne 3 ]; then
    echo $0: usage: ./skewer.zsh INPUTFOLDER OUTPUTFOLDER CPUs 
	echo "\nINPUTFOLDER: Folder containing forward and reverse reads (.fastq), format aaaaa_1.fastq.gz"
	echo "\nOUTPUTFOLDER: Folder for filtered reads, will be created new if not existing"
	echo "\nCPUs: threads/cores to be used"
    exit 1
fi

#set/get variables
INPUTFOLDER=$1
OUTPUTFOLDER=$2
CPUs=$3

#if adapterfiles do not exist
if [[ ! -f $FWD_ADAPTERS ]]; then
echo "FWD adapter file does not exist"
fi
if [[ ! -f $REV_ADAPTERS ]]; then
echo "FWD adapter file does not exist"
fi

echo "filtering reads (skewer) and fastqc"
which skewer

#####if .R1 or _1 format of fastq names
#if [[ -f $(ls -1 $INPUTFOLDER/*.R[1,2].fastq.gz | head -n1) ]]; then
#echo ".R1 format"
#ls -1 $INPUTFOLDER/*.R[1,2].fastq.gz | rev | cut -d"." -f 4- | cut -d"/" -f 1 | rev | sort | uniq > $INPUTFOLDER/QC.1.samples.lst
#fi
#if [[ -f $(ls -1 $INPUTFOLDER/*_[1,2].fastq.gz | head -n1) ]]; then
#echo "_1 format"
#ls -1 $INPUTFOLDER/*_[1,2].fastq.gz | rev | cut -d"." -f 2- | cut -d"_" -f 2- | cut -d"/" -f 1 | rev | sort | uniq > $INPUTFOLDER/QC.2.samples.lst
#fi



mkdir -p $OUTPUTFOLDER
ls -1 $INPUTFOLDER/*_[1,2].fastq.gz | rev | cut -d"." -f 2- | cut -d"_" -f 2- | cut -d"/" -f 1 | rev | sort | uniq > $OUTPUTFOLDER/samples.lst
#ls -1 $INPUTFOLDER/*_[1,2].fastq.gz | rev | cut -d"." -f 2- | cut -d"_" -f 2- | cut -d"/" -f 1 | rev | sort | uniq > $OUTPUTFOLDER/samples.lst
LISTtoPROCESS="$OUTPUTFOLDER/samples.lst"
N=$(wc -l $LISTtoPROCESS | cut -d" " -f1)
echo "processing $N samples (read-pair files)"

#filtering
#dryrun
if [[ -f $FWD_ADAPTERS ]]; then
echo "FWD adapter file found"
else echo "FWD adapter file not found"
fi
if [[ -f $REV_ADAPTERS ]]; then
echo "REV adapter file found"
else echo "REV adapter file not found"
fi
sleep 3s
cat $LISTtoPROCESS | parallel --no-notice -j 1 "echo $INPUTFOLDER/{}'_1'.fastq.gz $INPUTFOLDER/{}'_2'.fastq.gz"


#######
TMPFOLDER=$HOME/tmp
mkdir -p $TMPFOLDER

echo "___________________"

echo "skewer parameters: \n--mean-quality $AVG_QUAL_THRESHOLD \n--end-quality $END_QUAL_THRESHOLD \n-l $L bp"
sleep 1s
echo "starting filtering"

#dupedist instead dist HiSeq2000:40 HiSeq4000:2500
#main filtering (clumpify, skewer)
cat $LISTtoPROCESS | parallel --no-notice -j 2 "echo {}; clumpify.sh in1=$INPUTFOLDER/{}'_1'.fastq.gz in2=$INPUTFOLDER/{}'_2'.fastq.gz out1=$TMPFOLDER/{}.temp.R1.fastq.gz out2=$TMPFOLDER/{}.temp.R2.fastq.gz dedupe optical dupedist=$OPTICALDIST -Xmx30g passes=1 subs=1 k=31 spantiles=f 2> >(tee $OUTPUTFOLDER/{}.optical_duplicates.log >&1); echo clumpify_done; \

skewer -m pe $TMPFOLDER/{}.temp.R1.fastq.gz $TMPFOLDER/{}.temp.R2.fastq.gz -x $FWD_ADAPTERS -y $REV_ADAPTERS --mean-quality $AVG_QUAL_THRESHOLD --end-quality $END_QUAL_THRESHOLD -l $L -n -r 0.1 -z --format auto -t $CPUs -o $TMPFOLDER/{} ; echo skewer_done ; \

mv -f $TMPFOLDER/{}'-trimmed-pair1'.fastq.gz $OUTPUTFOLDER/{}.R1.fastq.gz ; mv -f $TMPFOLDER/{}'-trimmed-pair2'.fastq.gz $OUTPUTFOLDER/{}.R2.fastq.gz ; mv -f $TMPFOLDER/{}'-trimmed'.log $OUTPUTFOLDER/{}.skewer.log ; echo {} ; \
rm -f $TMPFOLDER/{}.temp.R1.fastq.gz $TMPFOLDER/{}.temp.R2.fastq.gz"

echo "finished filtering"

#ls -1 $INPUTFOLDER/*.fastq.gz | rev | cut -d"/" -f1 | cut -d"." -f3- | rev | parallel --no-notice -j 1 "echo {}; fastqc --noextract -f fastq --threads 10 --nogroup -o $INPUTFOLDER/fastqcs/ $INPUTFOLDER/{}.fastq.gz; echo {}"


#ls -1 $INPUTFOLDER/*.fastq.gz | rev | cut -d"/" -f1 | cut -d"." -f3- | rev | grep -P '705501|708501|712501|705502|706502|711502|705503|706503|708503|706504|707504|708504|709504|711504|712504|707505|712505|706507|708507|709507|710507|709508|710506' | parallel --no-notice -j 4 "echo {}; fastqc --noextract -f fastq --threads 10 --nogroup -o $INPUTFOLDER/fastqcs/ $INPUTFOLDER/{}.fastq.gz; echo {}"

#ls -1 $OUTPUTFOLDER/*.fastq.gz | rev | cut -d"/" -f1 | cut -d"." -f3- | rev | parallel --no-notice -j 4 "echo {}; fastqc --noextract -f fastq --threads 10 --nogroup -o $OUTPUTFOLDER/fastqcs/ $OUTPUTFOLDER/{}.fastq.gz; echo {}"

### summarize fastqc's with multiqc
#multiqc -o $INPUTFOLDER/multiqc/ $INPUTFOLDER/fastqcs/
#multiqc -o $OUTPUTFOLDER/multiqc/ $OUTPUTFOLDER/fastqcs/

#tar -cv $INPUTFOLDER/fastqcs/ | pigz -p $CPUs > $INPUTFOLDER/fastqcs.tar.gz
#tar -cv $OUTPUTFOLDER/fastqcs/ | pigz -p $CPUs > $OUTPUTFOLDER/fastqcs.tar.gz

#mv -f $INPUTFOLDER/fastqcs.tar.gz $INPUTFOLDER/QC
#mv -f $INPUTFOLDER/multiqc/ $INPUTFOLDER/QC

#mv -f $OUTPUTFOLDER/fastqcs.tar.gz $OUTPUTFOLDER/QC
#mv -f $OUTPUTFOLDER/multiqc/ $OUTPUTFOLDER/QC

#mkdir -p $OUTPUTFOLDER/QC/fastqcs
#mkdir -p $OUTPUTFOLDER/QC/multiqc
#mkdir -p $INPUTFOLDER/QC/fastqcs
#mkdir -p $INPUTFOLDER/QC/multiqc
#~/scripts/fastqc.folder.zsh $INPUTFOLDER $CPUs


