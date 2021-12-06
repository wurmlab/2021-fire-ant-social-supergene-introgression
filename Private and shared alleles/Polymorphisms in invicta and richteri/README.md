# BBBAA-BBABA test: alleles shared between SB and Sb of different _S. invicta_ and _S. richteri_

Rodrigo Pracana, QMUL 2021

## Samples

We get the samples from the twisst analysis

```sh

mkdir tmp
mkdir results
mkdir input

ln -sfr /data/archive/archive-SBCS-WurmLab/2021-solenopsis_introgression/2021-08-04-twisst/results/twisst_samples_low_missingness input/samples

```

Then choose the same number of samples per species:

```sh


# species in analysis
analysis_species <- c(
  "invicta",
  "richteri",
  "saevissima"
)

# Get twisst sample names (invicta, richteri, saevissima)
sample_table <- read.table("input/samples")
colnames(sample_table) <- c("sample", "species")

sample_table$variant <- gsub(".*_", "", sample_table$species)
sample_table$variant <- gsub("SB", "bb", sample_table$variant)
sample_table$variant <- gsub("Sb", "lb", sample_table$variant)
sample_table$variant[!sample_table$variant %in% c("lb", "bb")] <- "outgroup"

sample_table$species <- gsub("invicta/macdonaghi", "invicta", sample_table$species)
sample_table$species <- gsub("_S[Bb]", "", sample_table$species)

sample_table <- sample_table[sample_table$species %in% analysis_species, ]

sample_table$species_variant <- paste(sample_table$species,
                                      sample_table$variant,
                                      sep = "_")

table(sample_table$species_variant)
#         invicta_bb          invicta_lb         richteri_bb         richteri_lb
#                134                  59                  28                  32
# saevissima_outgroup
#                 11

# Remove North American samples
north_american_samples <- grep("^F[1-8]", sample_table$sample, value = TRUE)

sample_table <- sample_table[!sample_table$sample %in% north_american_samples, ]

# Vector of samples
invicta_bb  <- sample_table$sample[sample_table$species_variant == "invicta_bb"]
invicta_lb  <- sample_table$sample[sample_table$species_variant == "invicta_lb"]
richteri_bb <- sample_table$sample[sample_table$species_variant == "richteri_bb"]
richteri_lb <- sample_table$sample[sample_table$species_variant == "richteri_lb"]
saevissima  <- sample_table$sample[sample_table$species_variant == "saevissima_outgroup"]

# The smallest sample size is richteri_bb, with 28 samples
invicta_bb  <- invicta_bb[sample(length(invicta_bb), 28)]
invicta_lb  <- invicta_lb[sample(length(invicta_bb), 28)]
richteri_lb <- richteri_lb[sample(length(invicta_bb), 28)]

all_samples <- c(invicta_bb, invicta_lb, richteri_bb, richteri_lb, saevissima)

sample_table_chose <- sample_table[sample_table$sample %in% all_samples, ]
write.table(sample_table_chose, file="results/samples", quote=FALSE)

```

## VCF input

```sh

ln -sf /data/archive/archive-SBCS-WurmLab/rpracana/projects/2020-introgression/2021-04-introgression/results/gt.vcf.gz input/all.vcf.gz
ln -sf /data/archive/archive-SBCS-WurmLab/rpracana/projects/2020-introgression/2021-04-introgression/results/gt.vcf.gz.tbi input/all.vcf.gz.tbi

```

Include bi-allelic only:

```sh

module load bcftools/1.13
bcftools view --max-alleles 2 --exclude-types indels -O z \
  input/all.vcf.gz > tmp/biallelic_snps.vcf.gz
bcftools index -t tmp/biallelic_snps.vcf.gz

bcftools index --nrecords input/all.vcf.gz
# 2650251
bcftools index --nrecords tmp/biallelic_snps.vcf.gz
# 2600879

```

## Supergene location

```sh

ln -sf /data/archive/archive-SBCS-WurmLab/rpracana/projects/south_social_chrom_div/south_nf/analysis/2018-05-11-linkage-map/results/linkage_map_supergene.txt input

```

## GFF

```sh

ln -sfr /data/archive/archive-SBCS-WurmLab/db/genomic/S_invicta/Si_gnGA/annotation/SolInvGFF.corrected.mRNACDSintron.gff input/gnga.gff

ln -sfr /data/archive/archive-SBCS-WurmLab/db/genomic/S_invicta/Si_gnGA/assembly/gnG_20161213.fa input/ref.fa

```

## R

The following script loads the VCF for _S. invicta_, _S. richteri_, and the outgroup _S. saevissima_. For each site, the script identifies the ancestral and the derived allele and measures the derived allele frequency of each population (the SB and the Sb samples, respectively, of each species). The script is divided into five parts to parallelise the process.

```sh

module load R/4.0.2

Rscript allele_freq1.r 1> tmp/allele_freq1.r.out 2> tmp/allele_freq1.r.err &
Rscript allele_freq2.r 1> tmp/allele_freq2.r.out 2> tmp/allele_freq2.r.err &
Rscript allele_freq3.r 1> tmp/allele_freq3.r.out 2> tmp/allele_freq3.r.err &
Rscript allele_freq4.r 1> tmp/allele_freq4.r.out 2> tmp/allele_freq4.r.err &
Rscript allele_freq5.r 1> tmp/allele_freq5.r.out 2> tmp/allele_freq5.r.err &

```

The results are added to the directory `results/by_scaff`.

```sh

head -n1 results/by_scaff/private_and_shared_per_site_10.csv > results/snp_frequencies.csv
tail -n +2 --quiet results/by_scaff/private_and_shared_per_site_* >> results/snp_frequencies.csv

wc -l results/snp_frequencies.csv
# 807657 results/snp_frequencies.csv

grep -c supergene results/snp_frequencies.csv
# 39210

```


## Analyse

We use the rmarkdown script `shared_private.Rmd` to
