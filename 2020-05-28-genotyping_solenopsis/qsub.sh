#!/usr/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=15G
#$ -l h_rt=48:00:00
#$ -t 1-14074
#$ -m n     # don't email me
#$ -j y     # Join stdout and stderr
#$ -tc 2000 # allow more simultaneous jobs
#$ -o tmp/logs/

module load anaconda3
conda activate 2020-05-env

# Freebayes parameters
MINALTFRAC=0.40
MINALTN=4
MINCOV=4

# Input reference fasta
REF="input/reference_genome_dir/gng20170922wFex.fa"

# Input bam list (one bam per line)
BAM_LIST="tmp/bam.list"
if [[ ! -s $BAM_LIST ]] ; then
    echo "File $BAM_LIST is absent or empty, aborting."
    exit
fi

# Job is divided by region (read line number ${SGE_TASK_ID})
REGION=$(sed -n "${SGE_TASK_ID}p" input/reference_genome_dir/gnG.REGIONS.splitby.FILTER2.500bp.regions)
echo $REGION

# Output
# split output into 9(first) * 10 (last) = 90 subdirs
OUTDIR=tmp/freebayes_chunks/${SGE_TASK_ID:0:1}${SGE_TASK_ID: -1}
mkdir -p $OUTDIR

# Run freebayes
freebayes \
  --region $REGION \
  --fasta-reference $REF \
  --ploidy 2 \
  --report-genotype-likelihood-max \
  --use-mapping-quality \
  --genotype-qualities \
  --use-best-n-alleles 6 \
  --haplotype-length 3 \
  --min-mapping-quality 30 \
  --min-base-quality 25 \
  --min-alternate-fraction $MINALTFRAC \
  --min-alternate-total $MINALTN \
  --min-coverage $MINCOV \
  --use-reference-allele \
  --bam-list $BAM_LIST \
  | bgzip -f -@ 1 -c /dev/stdin \
  > $OUTDIR/$REGION.vcf.gz

tabix -fp vcf $OUTDIR/$REGION.vcf.gz

echo "ENDS"
