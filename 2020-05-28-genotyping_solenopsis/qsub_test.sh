#!/usr/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=15G
#$ -l h_rt=00:59:00
#$ -t 1-10
#$ -m n     # don't email me
#$ -j y     # Join stdout and stderr
#$ -tc 2000 # allow more simultaneous jobs
#$ -o test/logs/

module load anaconda3
conda activate 2020-05-env

# Freebayes parameters
MINALTFRAC=0.40
MINALTN=4
MINCOV=4

# Input (reference and bam files)
REF="input/reference_genome_dir/gng20170922wFex.fa"
# ls -1 input/bams/*.bam > tmp/bam.list

# Job is divided by region (read line number ${SGE_TASK_ID})
REGION=$(sed -n "${SGE_TASK_ID}p" input/reference_genome_dir/gnG.REGIONS.splitby.FILTER2.500bp.regions)
echo $REGION

# Output
# split output into 9(first) * 10 (last) = 90 subdirs
OUTDIR=test/freebayes_chunks/${SGE_TASK_ID:0:1}${SGE_TASK_ID: -1}
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
  --bam-list tmp/bam.list \
  | bgzip -f -@ 1 -c /dev/stdin \
  > $OUTDIR/$REGION.vcf.gz

tabix -fp vcf $OUTDIR/$REGION.vcf.gz
