
# optical duplicates
# mapping, markdup

BAMF="/scratch/ek/solenopsis/2018/bam"
cd $BAMF
ls -1 {gdo*bam,F*.bam,SR*.bam} | rev | cut -d "." -f2- | rev | wc -l
LIST=$(ls -1 {gdo*bam,F*.bam,SR*.bam} | rev | cut -d "." -f2- | rev)
echo $LIST

CPUs=20
N=$(echo $LIST | wc -l) && echo $N
for (( i = 1 ; i < $N+1 ; i++))
 do
SAMPLE=$(echo $LIST | sed -n $i'p') && echo $SAMPLE
TMPDIR="/dev/shm/$SAMPLE"
mkdir -p $TMPDIR
sambamba markdup -l 9 --tmpdir="$TMPDIR" -p -t $CPUs $BAMF/$SAMPLE.bam $BAMF/$SAMPLE.markdup.bam
sambamba index -p -t $CPUs $BAMF/$SAMPLE.markdup.bam
echo "$SAMPLE done"
done

N=$(echo $LIST | wc -l) && echo $N
for (( i = 1 ; i < $N+1 ; i++))
 do
SAMPLE=$(echo $LIST | sed -n $i'p') && echo $SAMPLE
ll $BAMF/$SAMPLE.bam
ll $BAMF/$SAMPLE.markdup.bam
done



##############################################################################################
# published data (SRA): NEE paper, Keller lab (PRJNA421367, SRP127040), Privman 2018 (SRX3960521), Wang et al. 2013 (2+6 individuals)

cat sample.info.acc.txt | grep -vP 'SRR9008101|SRR9008102|SRR9008107|SRR9008109|SRR9008114' | parallel -j3 --no-notice 'prefetch -v --progress {}'

CPUs=15
TMPSRAFOLDER="/dev/shm/solenopsis.tmp"
OUTFOLDER="/scratch/solenopsis.filtered.reads"
cat sample.info.acc.txt | parallel -j 1 --no-notice "echo {} && fasterq-dump {}/{}.sra --threads $CPUs --skip-technical --print-read-nr --temp /backup/tmp --details --progress --mem 9000Mb --outdir $TMPSRAFOLDER && echo {} && $HOME/scripts/skewer.SRA.zsh $TMPSRAFOLDER/{}.sra_1.fastq $TMPSRAFOLDER/{}.sra_2.fastq $OUTFOLDER $CPUs && rm -f $TMPSRAFOLDER/{}.sra_1.fastq $TMPSRAFOLDER/{}.sra_2.fastq"

FOLDER="privman2018.raw"
INPUTFOLDER="$BASEFOLDER/$FOLDER"
OUTPUTFOLDER="$BASEFOLDER/$FOLDER.skewer"
scripts/filter.skewer.privman2018.zsh $INPUTFOLDER $OUTPUTFOLDER $CPUs




##############################################################################################
# mapping

ulimit -n 10240
REF="/dev/shm/ref/gng20170922.fa"
TMPDIR2="/dev/shm/bams"
mkdir -p $TMPDIR2
BAMOUT="/scratch/solenopsis.bams"
mkdir -p $BAMOUT

LIST=$(ls -1 *.R1.fq.gz | rev | cut -d "." -f4 | rev | tr '\n' '|' | sed "s,|$,,g")

CPUs=50
N=$(echo $LIST | wc -l) && echo $N

for (( i = 1 ; i < $N+1 ; i++))
 do
SAMPLE=$(echo $LIST | sed -n $i'p') && echo $SAMPLE
TMPDIR="/dev/shm/$SAMPLE"
mkdir -p $TMPDIR
#
# ReadGroup Information to include into Bamfile
SAMPLEID=$SAMPLE
ID=$SAMPLE
CLADE=$SAMPLE
LIBRARY="TruSeqPCR"
PLATFORM="IlluminaHiSeq2500"
TIMEOFSEQ="Yan2020"
READGROUPHEADER="@RG\tID:$ID\tSM:$SAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY\tPU:$CLADE\tDT:$TIMEOFSEQ"
echo $READGROUPHEADER
#
bwa-mem2 mem -t $CPUs -B 6 -E 2 -L25,25 -U 50 -R $READGROUPHEADER -v 1 -T 50 -h 4,200 -a -V -Y -M $REF $SAMPLE.R1.fq.gz $SAMPLE.R2.fq.gz | sambamba view --sam-input --with-header --ref-filename="$REF" --format=bam -l 0 /dev/stdin | sambamba sort -m 50G -l 0 --tmpdir="$TMPDIR" -p -t $CPUs -o /dev/stdout /dev/stdin > $TMPDIR2/$SAMPLE.bam
#
sambamba index -p -t $CPUs $TMPDIR2/$SAMPLE.bam $TMPDIR2/$SAMPLE.bai
#
sambamba markdup -l 8 --tmpdir="$TMPDIR" -p -t $CPUs $TMPDIR2/$SAMPLE.bam $BAMOUT/bam/$SAMPLE.bam
#
sambamba index -p -t $CPUs $BAMOUT/bam/$SAMPLE.bam
rm -rf $TMPDIR $TMPDIR2/$SAMPLE.bam $TMPDIR2/$SAMPLE.bai
echo "$SAMPLE done"
done



################################ mapping of haploid males

# check if indices exist

REF="ref/gng20170922.fa"
FILE="$REF.bwt.2bit.64"
if [ -f "$FILE" ]; then
    echo "bwa-mem2 index exists already, continuing..."
else
    echo "bwa-mem2 index does not yet exist, creating..."
    bwa-mem2 index $REF
fi
FILE="$REF.fai"
if [ -f "$FILE" ]; then
    echo "samtools index exists already, continuing..."
else
    echo "samtools index does not yet exist, creating..."
    samtools faidx $REF
fi


## store index in RAM /dev/shm for faster access

mkdir -p /dev/shm/ref
cp ref/gng20170922.f* /dev/shm/ref/
REF="/dev/shm/ref/gng20170922.fa"

##mapping

mkdir bam

#### mapping loop USA samples
ls -1 skewer.reads/reads_USA/*.R1.fastq.gz | cut -d"/" -f3 | cut -d"." -f1 | sort | uniq > mappingsamplenamesUSA.lst

REF="/dev/shm/ref/gng20170922.fa"
CPUs=50
SAMPLELST="mappingsamplenamesUSA.lst"
N=$(wc -l $SAMPLELST | cut -d' ' -f 1) && echo $N
for (( i = 15 ; i < $N+1 ; i++))
 do
 SAMPLE=$(cat $SAMPLELST | sed -n $i'p') && echo $SAMPLE

# inputs
FWD="$SAMPLE.R1.fastq.gz"
REV="$SAMPLE.R1.fastq.gz"

INDIR="skewer.reads/reads_USA"
OUTDIR="bam"
TMPDIR="/backup/34324234234"
mkdir -p $TMPDIR

# ReadGroup Information to include into Bamfile
SAMPLEID=$SAMPLE
ID=$(echo -n $SAMPLE | cut -d"_" -f2)
CLADE=$(echo -n $SAMPLE | cut -d"_" -f3)
LIBRARY="TruSeq"
PLATFORM="IlluminaHiSeq4000"
TIMEOFSEQ="2015-2018"
READGROUPHEADER="@RG\tID:$ID\tSM:$SAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY\tPU:$CLADE\tDT:$TIMEOFSEQ"
echo $READGROUPHEADER

## mapping to bam files
bwa-mem2 mem -t $CPUs -B 6 -E 2 -L25,25 -U 50 -R $READGROUPHEADER -v 1 -T 30 -h 4,200 -a -V -Y -M  $REF $INDIR/$FWD $INDIR/$REV | sambamba view --sam-input --with-header --ref-filename="$REF" --format=bam -l 0 /dev/stdin | sambamba sort -m 30G -l 8 --tmpdir="$TMPDIR" -p -t $CPUs -o $OUTDIR/$ID.bam /dev/stdin
sambamba index -p -t $CPUs $OUTDIR/$ID.bam $OUTDIR/$ID.bai


############################# mapping of Privman 2018 samples
ls -1 skewer.reads/privman2018.raw.skewer/*.R1.fq.gz | cut -d"/" -f3 | cut -d"." -f1 | sort | uniq > mappingsamplenamesPrivman.lst

REF="/dev/shm/ref/gng20170922.fa"
CPUs=50
SAMPLELST="mappingsamplenamesPrivman.lst"
N=$(wc -l $SAMPLELST | cut -d' ' -f 1) && echo $N
for (( i = 1 ; i < $N+1 ; i++))
 do
 SAMPLE=$(cat $SAMPLELST | sed -n $i'p') && echo $SAMPLE

# inputs
FWD="$SAMPLE.R1.fq.gz"
REV="$SAMPLE.R1.fq.gz"

INDIR="skewer.reads/privman2018.raw.skewer"
OUTDIR="bam"
TMPDIR="/backup/tmp1234"

# ReadGroup Information to include into Bamfile
SAMPLEID=$SAMPLE
ID=$(echo -n $SAMPLE | cut -d"_" -f2)
CLADE=$(echo -n $SAMPLE | cut -d"_" -f3)
LIBRARY="TruSeq"
PLATFORM="IlluminaHiSeq4000"
TIMEOFSEQ="2015-2018"
READGROUPHEADER="@RG\tID:$ID\tSM:$SAMPLEID\tPL:$PLATFORM\tLB:$LIBRARY\tPU:$CLADE\tDT:$TIMEOFSEQ"
echo $READGROUPHEADER

## mapping to bam files
bwa-mem2 mem -t $CPUs -B 6 -E 2 -L25,25 -U 50 -R $READGROUPHEADER -v 1 -T 50 -h 4,200 -a -V -Y -M  $REF $INDIR/$FWD $INDIR/$REV | sambamba view --sam-input --with-header --ref-filename="$REF" --format=bam -l 0 /dev/stdin | sambamba sort -m 60G -l 8 --tmpdir="$TMPDIR" -p -t $CPUs -o $OUTDIR/$ID.bam /dev/stdin
sambamba index -p -t $CPUs $OUTDIR/$ID.bam $OUTDIR/$ID.bai

echo "mapping of $SAMPLE done"
done
