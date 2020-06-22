#!/usr/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=00:30:00
#$ -t 1-14074
#$ -m n     # don't email me
#$ -j y     # Join stdout and stderr
#$ -tc 2000 # allow more simultaneous jobs
#$ -o tmp/drop_genotype_logs/

module load anaconda3
conda activate 2020-05-env

# The previous step (variant filtering) produced one VCF file per region
# This was split into directories:
# split output into 9(first) * 10 (last) = 90 subdirs
INDIR=tmp/filtered_chunks/${SGE_TASK_ID:0:1}${SGE_TASK_ID: -1}

# Job is divided by region (read line number ${SGE_TASK_ID})
REGION=$(sed -n "${SGE_TASK_ID}p" input/reference_genome_dir/gnG.REGIONS.splitby.FILTER2.500bp.regions)
echo region:$REGION

# Check if VCF is present
if [[ ! -s $INDIR/$REGION.vcf.gz ]] ; then
    echo "File $INDIR/$REGION.vcf.gz is absent or empty, aborting."
    exit
fi
echo "input vcf:$INDIR/$REGION.vcf.gz"

# Output
OUTDIR=tmp/drop_genotype_chunks/${SGE_TASK_ID:0:1}${SGE_TASK_ID: -1}
mkdir -p $OUTDIR

echo "output vcf:$OUTDIR/$REGION.vcf.gz"

zcat $INDIR/$REGION.vcf.gz \
  | vcfkeepinfo - AB AC AF AN DP MEANALT NUMALT RO SRF SRR TYPE NS ODDS \
  | grep -v '#' | cut -f1-9 > $OUTDIR/$REGION.vcf

echo "ENDS"
