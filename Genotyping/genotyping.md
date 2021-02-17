# Genotyping

## Aim

Genotype a group of _Solenopsis sp_ sample.

## Data

Many of the samples were collected and sequenced by Eckart Stolle (and colleagues) from the native South American species range (Argentina, Brazil and Uruguay). DNA was extracted from each sample and sequenced (Illumina HiSeq).

(((We added to these samples a data set produced in the lab of Laurent Keller () and another in the lab of Eyal Privman ().))) Maybe this in main project README

Sequenced reads were filtered (with XXX) and aligned to the reference genome XXX with the aligner XXX.

The reference genome is available here:

```sh

REF=/data/home/btw586/db/genomic/S_invicta/2017-09-22-Si_gng20170922_eckart/gng20170922wFex.fa

```

The reads are available here:

```sh

ls /data2/archive/archive-SBCS-WurmLab/db/genomic/reads/S_invicta

```

The alignments are available here:

```sh

ls /data/archive/archive-SBCS-WurmLab/db/genomic/reads/S_invicta/2020-05-bams_388_eckart/bams_all_renamed

# @PG     ID:bwa  PN:bwa  CL:bwa-mem2 mem -t 30 -B 6 -E 2 -L25,25 -U 50 -R @RG\tID:AdR11-2-bigB-p\tSM:Sol_AdR11-2-bigB-p_clade9A_708508\tPL:IlluminaHiSeq4000\tLB:TruSeq\tPU:clade9A\tDT:2015-2018 -v 1 -T 50 -h 4,200 -a -V -Y -M /dev/shm/ref/gng20170922wFex.fa mappingsamples/Sol_AdR11-2-bigB-p_clade9A_708508.R1.fq.gz mappingsamples/Sol_AdR11-2-bigB-p_clade9A_708508.R1.fq.gz    VN:2.0pre2
# @PG     ID:sambamba     CL:view --sam-input --with-header --ref-filename=/dev/shm/ref/gng20170922wFex.fa --format=bam -l 0 /dev/stdin   PP:bwa  VN:0.6.7

ls -1 -d /data/archive/archive-SBCS-WurmLab/db/genomic/reads/S_invicta/2020-05-bams_388_eckart/bams_all_renamed/*.bam > tmp/bam.list

BAMLST=tmp/bam.list

```

There are bams for 388 samples.

Most of the samples (including all ingroup samples) consist of haploid male ants. There are XXX exceptions, which consist of pools of diploid workers. These exceptions form some of the outgroup samples.

## Overview of protocol

We use freebayes to call variants among the samples. This is done in two steps:
1. Identification of variant sites
2. Genotyping of the variant sites

Following this, the genotypes are filtered.

### Identification of variant sites

The reference genome was divided into regions, over which the steps were parallelised.

```sh

# Get region
#   And region size
REGIONS=/data2/archive/archive-SBCS-WurmLab/db/genomic/S_invicta/2017-09-22-Si_gng20170922_eckart/gnG.REGIONS.splitby.FILTER2.500bp.regions

# Some regions were basically empty
grep -cv ":1-1$" $REGIONS
# 14059

```

This was done with freebayes, following the parameters below, parallelised over each line of the `$REGIONS` file.

```sh

REGION=$(sed -n "1p" $REGIONS) # The first line of the region file
MINALTFRAC=0.40
MINALTN=4
MINCOV=4

freebayes --region REGION \
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
  --bam-list $BAMLST \
  | bgzip -f -@ 1 -c /dev/stdin \
  > site/${REGION}.vcf.gz \

tabix -fp vcf site/${REGION}.vcf.gz

```

To filter the sites, removed any within the ranges of the following BED file, which represent highly repetitive regions.

```sh

REPEAT_REGIONS=/data/home/btw586/db/genomic/S_invicta/2017-09-22-Si_gng20170922_eckart/gng20170922.fa.FILTER2.bed.gz

```

We also filter by quality and by number of reads in the forward and backward strand.


```sh

# qsub filter_sites_qsub.sh

zcat site/${REGION}.vcf.gz \
  | vcfintersect -v -l -b $FILTERBED \
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
  > filtered_${REGION}.vcf.gz

tabix -fp vcf site_filtered/${REGION}.vcf.gz

```

Then, we removed the genotype information from the VCFs.

```sh

# qsub drop_genotypes_qsub.sh

zcat site_filtered/${REGION}.vcf.gz \
  | vcfkeepinfo - AB AC AF AN DP MEANALT NUMALT RO SRF SRR TYPE NS ODDS \
  | grep -v '#' | cut -f1-9 > drop_genotype/$REGION.vcf

```

We then fuse all the regions.


```sh

VCF1=drop_genotype/$(head -n 1 tmp/regions).vcf

# Get the header from one VCF
grep "##" $VCF1 > tmp/fused_initial.vcf
grep -m1 "#CHROM" $VCF1 | cut -f1-9 >> tmp/fused_initial.vcf
# Concatenate all VCF files
cat tmp/drop_genotype_chunks/*/*.vcf >> tmp/fused_initial.vcf

bcftools sort tmp/fused_initial.vcf | bgzip -c > tmp/variant_guide.vcf.gz
tabix -fp vcf tmp/variant_guide.vcf.gz

rm -f tmp/fused_initial.vcf

```

## Genotyping each site

For each site in the VCF above (`tmp/variant_guide.vcf.gz`), we genotype all the individuals. This is parallelised over the region file, as above.

```sh


REGION=$(sed -n "1p" $REGIONS) # The first line of the region file

MINALTFRAC=0.35
MINALTN=2
MINCOV=2

REFVCF=tmp/variant_guide.vcf.gz

freebayes \
  --region {} \
  --fasta-reference $REF \
  --ploidy 2 \
  --haplotype-basis-alleles $REFVCF \
  --report-genotype-likelihood-max \
  --use-mapping-quality \
  --genotype-qualities \
  --use-best-n-alleles 4 \
  --haplotype-length -1 \
  --min-mapping-quality 30 \
  --min-base-quality 28 \
  --min-alternate-fraction $MINALTFRAC \
  --min-alternate-total $MINALTN \
  --min-coverage $MINCOV \
  --use-reference-allele \
  --bam-list $BAMLST \
  | bgzip -f -@ 1 -c /dev/stdin \
  > genotype/${REGION}.vcf.gz \

tabix -fp vcf genotype/${REGION}.vcf.gz

```

We then filter the resulting VCF file.

```sh

zcat $OUTPUTFOLDER/vcf-genotyping/raw/{}.vcf.gz \
  | vcfintersect -v -l -b $FILTERBED \
  | vcffilter -f 'QUAL > 30' \
  | bcftools view -e "NUMALT=1 & ((INFO/SRF)<=($MINALTN/2) | (INFO/SRR)<=($MINALTN/2))" - \
  | vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED \
  | vcffixup - \
  | vcfstreamsort \
  | vt normalize -n -r $REF -q - \
  | vcfuniqalleles \
  | vcfnoindels \
  | vcfnumalt - \
  | vcfnulldotslashdot \
  | vcfstreamsort \
  | vcfuniq \
  | bgzip -f -@ 1 -c /dev/stdin \
  > genotype_filter/{}.vcf.gz

tabix -fp vcf genotype_filter/{}.vcf.gz

```

We then concatenate all the VCF chunks.


```sh

VCF1=genotype_filter/$(head -n 1 tmp/regions).vcf

# Get the header from one VCF
grep "##" $VCF1 > tmp/fused_initial.vcf
grep -m1 "#CHROM" $VCF1 | cut -f1-9 >> tmp/fused_initial.vcf
# Concatenate all VCF files
cat genotype_filter/*.vcf >> tmp/fused_genotypes.vcf

bcftools sort tmp/fused_genotypes.vcf \
  | bgzip -c > genotypes.vcf.gz
tabix -fp vcf genotypes.vcf.gz

rm -f tmp/fused_genotypes.vcf

```

We have a list of samples to remove. We remove these, then we filter each genotype by coverage, and turn each heterozygous genotype into "missing".

```sh
#
bcftools view \
  --samples-file ^tmp/remove_samples.txt \
  genotypes.vcf.gz \
  | bcftools +setGT - -- -t q -i "FMT/DP < 2" -n . \
  | bcftools +setGT - -- -t q -i "FMT/GQ < 1 & FMT/DP < 10 " -n . \
  | bcftools +setGT - -- -t ./x -n . \
  | bcftools +setGT - -- -t q -i 'GT="het"' -n . \
  | bgzip -c > genotypes.het_to_missing.vcf.gz

tabix -p vcf genotypes.het_to_missing.vcf.gz

```

We also remove sites where more than 25% of individuals have "missing".

```sh

NSAMPLES=$(bcftools query --list-samples genotypes.het_to_missing.vcf.gz | wc -l) && echo $NSAMPLES
NMISSING10=$(echo -n $NSAMPLES | awk '{ $1=sprintf("%.0f",$1*0.1*2)} {print $1;}') && echo $NMISSING10
NMISSING=$((echo $NSAMPLES-$NMISSING10/2 | bc -l) | awk '{ $1=sprintf("%.0f",$1)} {print $1;}') && echo $NMISSING

bcftools view -Oz -i 'F_MISSING<0.1' genotypes.het_to_missing.vcf.gz \
  > genotypes.het_filter_miss_filter.vcf.gz
tabix -p vcf genotypes.het_filter_miss_filter.vcf.gz

```
