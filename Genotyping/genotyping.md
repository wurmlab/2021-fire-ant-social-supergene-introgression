# Genotyping

Here we summarise the  process of genotyping the _Solenopsis sp_ samples.

## Data

Sequenced reads were filtered with Skewer (v0.2.2), with minimum length of 80 for 100 bp reads and of 100 for 150 bp reads (--min $L)

```sh
skewer -m pe  \
  $INPATH/$NAME.sra_1.fastq \
  $INPATH/$NAME.sra_2.fastq \
  -x foward_adapters.txt \
  -y reverse_adapters.txt \
  --mean-quality 20 \
  --end-quality 15 \
  --min $L \
  -n \
  -r 0.1 \
  --compress \
  --format auto \
  -t $CPUs \
  --output $OUTPUTFOLDER/$NAME

```

We used bwa-mem2 (v2.0pre2) to align reads to the reference genome Si_gnGA.

```sh
REF=gng20170922wFex.fa
bwa-mem2 mem -t $CPUs -B 6 -E 2 -L25,25 -U 50 -R $READGROUPHEADER -v 1 -T 50 -h 4,200 -a -V -Y -M $REF $SAMPLE.R1.fq.gz $SAMPLE.R2.fq.gz

```

## Overview of protocol

We use Freebayes (v1.2.0) to call variants among the samples. This is done in two steps:
1. Identification of variant sites
2. Genotyping of the variant sites

Following this, the genotypes are filtered.

### Identification of variant sites

The reference genome was divided into regions, over which the steps were parallelised.

```sh

# Get region
#   And region size
REGIONS=regions.txt

```

This was done with Freebayes, following the parameters below, parallelised over each line of the `$REGIONS` file.

```sh

mkdir site

REGION=$(sed -n "1p" $REGIONS) # The first (or nth) line of the region file
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

To filter the sites, we removed any within the ranges of the following BED file, which represent highly repetitive regions.

```sh

REPEAT_REGIONS=repetitive.bed

```

We also filter by quality and by number of reads in the forward and backward strand.


```sh

mkdir site_filtered

zcat site/${REGION}.vcf.gz \
  | vcfintersect -v -l -b $REPEAT_REGIONS \
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
  > site_filtered/${REGION}.vcf.gz

tabix -fp vcf site_filtered/${REGION}.vcf.gz

```

Then, we removed the genotype information from the VCFs.

```sh

mkdir drop_genotype

zcat site_filtered/${REGION}.vcf.gz \
  | vcfkeepinfo - AB AC AF AN DP MEANALT NUMALT RO SRF SRR TYPE NS ODDS \
  | grep -v '#' | cut -f1-9 > drop_genotype/$REGION.vcf

```

We then fuse all the regions.


```sh

# Get the header from one VCF
VCF1=drop_genotype/$(head -n 1 $REGIONS).vcf
grep "##" $VCF1 > tmp/fused_initial.vcf
grep -m1 "#CHROM" $VCF1 | cut -f1-9 >> tmp/fused_initial.vcf
# Concatenate all VCF files
cat drop_genotype/*.vcf >> fused_initial.vcf

bcftools sort fused_initial.vcf | bgzip -c > variant_guide.vcf.gz
tabix -fp vcf variant_guide.vcf.gz

rm -f fused_initial.vcf

```

## Genotyping each site

For each site in the VCF above (`variant_guide.vcf.gz`), we genotype all the individuals. This is parallelised over the region file, as above.

```sh

mkdir genotype

REGION=$(sed -n "1p" $REGIONS) # The first (or nth) line of the region file

MINALTFRAC=0.35
MINALTN=2
MINCOV=2

REFVCF=variant_guide.vcf.gz

freebayes \
  --region $REGION \
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

mkdir genotype_filter

zcat genotype/${REGION}.vcf.gz \
  | vcfintersect -v -l -b $REPEAT_REGIONS \
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
  > genotype_filter/${REGION}.vcf.gz

tabix -fp vcf genotype_filter/${REGION}.vcf.gz

```

We then concatenate all the VCF chunks.


```sh

VCF1=genotype_filter/$(head -n 1 $REGIONS).vcf

# Get the header from one VCF
grep "##" $VCF1 > fused_initial.vcf
grep -m1 "#CHROM" $VCF1 | cut -f1-9 >> fused_initial.vcf
# Concatenate all VCF files
cat genotype_filter/*.vcf >> fused_genotypes.vcf

bcftools sort fused_genotypes.vcf \
  | bgzip -c > genotypes.vcf.gz
tabix -fp vcf genotypes.vcf.gz

rm -f fused_genotypes.vcf

```

We filter each genotype by coverage, and turn each heterozygous genotype into "missing". We removed 18 samples for which more than 25% of variant sites could not be genotyped, and the nine others had particularly high numbers of heterozygous sites indicating that they are likely diploid.

```sh
#
bcftools view \
  --samples-file ^remove_samples.txt \
  genotypes.vcf.gz \
  | bcftools +setGT - -- -t q -i "FMT/DP < 2" -n . \
  | bcftools +setGT - -- -t q -i "FMT/GQ < 1 & FMT/DP < 10 " -n . \
  | bcftools +setGT - -- -t ./x -n . \
  | bcftools +setGT - -- -t q -i 'GT="het"' -n . \
  | bgzip -c > genotypes.het_to_missing.vcf.gz

tabix -p vcf genotypes.het_to_missing.vcf.gz

```

We also remove sites where more than 25% of individuals have a "missing" genotype.

```sh

bcftools view -Oz -i 'F_MISSING<0.25' genotypes.het_to_missing.vcf.gz \
  > genotypes.het_filter_miss_filter.vcf.gz
tabix -p vcf genotypes.het_filter_miss_filter.vcf.gz

```
