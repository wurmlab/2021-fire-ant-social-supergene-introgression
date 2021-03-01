#!/usr/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=15G
#$ -l h_rt=48:00:00
#$ -t 1-14074
#$ -m n     # don't email me
#$ -j y     # Join stdout and stderr
#$ -tc 2000 # allow more simultaneous jobs
#$ -o tmp/filtered_logs_2/

module load anaconda3
conda activate 2020-05-env

# Freebayes parameters
MINALTN=4

# Input reference fasta
REF="input/reference_genome_dir/gng20170922wFex.fa"
if [[ ! -s $REF ]] ; then
    echo "File $REF is absent or empty, aborting."
    exit
fi

# The previous step (variant calling) produced one VCF file per region
# This was split into direcotories:
# split output into 9(first) * 10 (last) = 90 subdirs
INDIR=tmp/genotyping_raw/${SGE_TASK_ID:0:1}${SGE_TASK_ID: -1}

# Job is divided by region (read line number ${SGE_TASK_ID})
REGION=$(sed -n "${SGE_TASK_ID}p" input/reference_genome_dir/gnG.REGIONS.splitby.FILTER2.500bp.regions)
echo region:$REGION

# Check if VCF is present
if [[ ! -s $INDIR/$REGION.vcf.gz ]] ; then
    echo "File $INDIR/$REGION.vcf.gz is absent or empty, aborting."
    exit
fi
echo "input vcf:$INDIR/$REGION.vcf.gz"

# BED file with exact location of used regions
FILTERBED="input/reference_genome_dir/gnG.REGIONS.splitby.FILTER2.500bp.bed"
if [[ ! -s $FILTERBED ]] ; then
    echo "File $FILTERBED is absent or empty, aborting."
    exit
fi

# Output
OUTDIR=tmp/genotyping_filtered/${SGE_TASK_ID:0:1}${SGE_TASK_ID: -1}
mkdir -p $OUTDIR

echo "output vcf:$OUTDIR/$REGION.vcf.gz"

zcat $INDIR/$REGION.vcf.gz \
  | vcfintersect -l -b $FILTERBED \
  | vcffilter -f 'QUAL > 30' \
  | bcftools view -e "NUMALT=1 & ((INFO/SRF)<=($MINALTN/2) | (INFO/SRR)<=($MINALTN/2))" - \
  | vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED \
  | vcffixup - \
  | vcfstreamsort \
  | vt normalize -n -r $REF -q - \
  | vcfuniqalleles \
  | vcfnoindels \
  | vcfnumalt - \
  | vcfstreamsort \
  | vcfuniq \
  | bgzip -f -@ 1 -c /dev/stdin \
  > $OUTDIR/$REGION.vcf.gz

tabix -fp vcf $OUTDIR/$REGION.vcf.gz

echo "ENDS"
