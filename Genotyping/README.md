# Genotyping

Here we summarise the  process of genotyping the _Solenopsis sp_ samples.

## Sequence read filtering and mapping

Sequenced reads were filtered with Skewer (v0.2.2), with minimum length of 80 for 100 bp reads and of 100 for 150 bp reads (--min $L). Clumpify was used to remove optical duplicates. This process is detailed in `read_filter.zsh`.

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
  --output $OUTPUTFOLDER/$NAME

```

We used bwa-mem2 (v2.0pre2) to align reads to the reference genome Si_gnGA.

```sh

REF=gng20170922wFex.fa

bwa-mem2 mem -B 6 -E 2 -L25,25 -U 50 -R $READGROUPHEADER -v 1 -T 50 -h 4,200 -a -V -Y -M $REF $SAMPLE.R1.fq.gz $SAMPLE.R2.fq.gz

```

The process of read mapping is presented in detail in `read_mapping.sh`.

## Single-copy gene identification

We identified the Hymenopteran single-copy genes in the _S. invicta_ genome using the tool BUSCO. We used the database of Hymenopteran single-copy genes.

```sh

LINEAGE="hymenoptera_odb10"
LINEAGENAME=$(echo $LINEAGE | rev | cut -d"/" -f1 | rev)
AUGUSTUSSPECIES="bombus_impatiens1"

busco -i ref.fa --force \
  --mode genome \
  --lineage_dataset $LINEAGE \
  --offline \
  --augustus_species $AUGUSTUSSPECIES \
  --evalue 1e-03 \
  --limit 5 \
  --out solenopsis_invicta.busco4

```

We transformed the output to create a BED file of single-copy genes, extended 1000bp upstream and downstream. This process is detailed in `single_copy_regions.sh`.

## Identification of regions with very high coverage

We filtered out variants in regions with very high coverage because these are likely to represent collapsed repeats. Coverage was measured with the tools `mosdepth`.

```sh

# Ran for each $SAMPLE
mosdepth -t $CPUs \
  --use-median \
  --by solenopsis_invicta.busco4/busco4.genes.bed \
  $SAMPLE.mosdepth.busco4 \
  bams/$SAMPLE.bam

```

This process is detailed in `normal_coverage_regions.sh`.

## Variant calling

We use Freebayes (v1.2.0) to call variants among the samples. This is done in two steps:
1. Identification of variant sites
2. Genotyping of the variant sites

Following this, the genotypes are filtered.

We summarise the process of variant calling below. The full script used is detailed in `variant_calling.sh`.

### Identification of variant sites

The reference genome was divided into regions, each representing a single-copy gene, over which the steps were parallelised.

```sh

# Get region
REGIONS=regions.txt

```

This was done with Freebayes, following the parameters below, parallelised over each line of the `$REGIONS` file.

```sh

mkdir site

REGION=$(sed -n "1p" $REGIONS) # The first (or nth) line of the region file
MINALTFRAC=0.40
MINALTN=4
MINCOV=4

freebayes --region $REGION \
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

To filter the sites, we ran the script presented in `site_filter.zsh`. This script includes the removal of any variants within the ranges of a BED file representing highly repetitive regions of the genome.

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

We then filter the resulting VCF file, which includes the removal of SNPs in regions with very high coverage and repetitive regions.

```sh

mkdir genotype_filter

zcat genotype/${REGION}.vcf.gz \
  | vcfintersect -v -l -b $HIGH_COV_REGIONS \
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

For genotype filtering, we divided the VCF into a VCF for each sample.

```sh

mkdir -p vcf-genotyping/samples

vcfsamplenames genotypes.vcf.gz > vcf-genotyping/name_samples

cat vcf-genotyping/name_samples \
  | parallel -j 10 \
  "vcfkeepsamples genotypes.vcf.gz {} | bgzip -c \
    > vcf-genotyping/samples/{}.vcf.gz \
    && tabix -fp vcf vcf-genotyping/samples/{}.vcf.gz"

```

We filter each genotype with the following script:

```sh

mkdir -p vcf-genotyping/samples-filtered

MINCOV=2

cat vcf-genotyping/name_samples \
parallel genotype_filter.zsh \
  vcf-genotyping/samples/{}.vcf.gz \
  $MINCOV \
  {}
  vcf-genotyping/samples-filtered

```

This script sets genotypes with low coverage (DP<`$MINCOV`) and low genotype quality (GQ<1) to missing. It also identifies sites with a heterozygous genotypes where the proportion of reads supporting the reference allele is within 0.25 and 0.75 (i.e., where the call is ambiguous). For the pooled outgroup samples, these genotypes are set to a random allele, and for the haploid samples, these genetypes are set to missing. Heterozygous sites where the reference allele proportion is outside this range are set to homozygous for the allele with the highest proportion of supported reads. Note that this script requires a file named `shuffle_file`, with random 0 or 1 with many lines as there are variants in the VCF file (we made it with the tool `shuffle`).

We merged the filtered genotype VCFs into a single file `genotype_filter.vcf.gz`. We then removed sites where more than 25% of individuals have a "missing" genotype.

```sh

# Identify site count for which more than 25% samples have a missing genotype
NSAMPLES=$(vcfsamplenames genotype_filter.vcf.gz | wc -l)
NMISSING25=$(echo -n $NSAMPLES | awk '{ $1=sprintf("%.0f",$1*0.25*2)} {print $1;}')
NMISSING=$((echo $NSAMPLES-$NMISSING25/2 | bc -l) | awk '{ $1=sprintf("%.0f",$1)} {print $1;}')

zcat genotype_filter.vcf.gz \
  | vcftools --recode --recode-INFO-all -c --vcf - --max-missing-count $NMISSING25 \
  | bgzip -f -c /dev/stdin > genotypes.max_missing25.vcf.gz
tabix -fp vcf genotypes.max_missing25.vcf.gz

```
