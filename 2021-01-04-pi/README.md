# Nucleotide diversity in three Sb clades

We aim to measure nucleotide diversity (pi) for each of three clades of Sb: one carried by S. invicta/macdonaghi populations and two carried by S. richteri populations.

## Coordinates for BUSCO regions

Our study only considers SNPs located in the exons, introns and 1000 basepairs upstream and downstream of BUSCO genes (hymenopteran single copy genes). The coordinates for study are here:

```sh

wc -l input/gng20170922.busco4.complete.extended1000bp.coordinates-fixed.merged.bed
# 4041

```

## VCF for correct regions

```sh

ls input/

```

```sh

wc -l input/richteri1_littleb.txt
wc -l input/richteri2_littleb.txt
wc -l input/invicta1_littleb.txt
# 34 input/richteri1_littleb.txt
# 10 input/richteri2_littleb.txt
# 62 input/invicta1_littleb.txt

```

The files were not formatted properly:

```sh

cat input/richteri1_littleb.txt \
  |  tr "−" "-" | ruby -pe 'gsub(/---/, "-")' > tmp/richteri1_littleb.txt

cat input/richteri2_littleb.txt \
  | tr "−" "-" | ruby -pe 'gsub(/---/, "-")' > tmp/richteri2_littleb.txt

cat input/invicta1_littleb.txt \
  | tr "−" "-" | ruby -pe 'gsub(/---/, "-")' > tmp/invicta1_littleb.txt

```

## Change VCF so it can be open by PopGenome

PopGenome needs a pseudo-diploid VCF. Any call with a non-missing genotype is already pseudo-diploid ("0|0" or "1|1"), but missing calls are still ".".

```sh

mkdir -p tmp

conda activate 2020-05-env

bcftools view input/busco4.merged.filtered.SNP.genotypes.MaxMissing25.vcf.gz \
  | bcftools annotate -x ^FORMAT/GT,FORMAT/DP \
  | ruby -pe 'gsub(/\.:/, ".|.:")' \
  | bgzip -c \
  > tmp/gt.vcf.gz

tabix -p vcf tmp/gt.vcf.gz

bcftools query -l tmp/gt.vcf.gz > tmp/samplelist


mkdir results
bcftools view --samples-file tmp/richteri1_littleb.txt tmp/gt.vcf.gz \
  | bcftools view --min-ac=1 - \
  | bgzip -c > results/richteri1.vcf.gz

bcftools view --samples-file tmp/richteri2_littleb.txt tmp/gt.vcf.gz \
  | bcftools view --min-ac=1 - \
  | bgzip -c > results/richteri2.vcf.gz

bcftools view --samples-file tmp/invicta1_littleb.txt tmp/gt.vcf.gz \
  | bcftools view --min-ac=1 - \
  | bgzip -c > results/invicta1.vcf.gz

tabix -p vcf results/richteri1.vcf.gz  &
tabix -p vcf results/richteri2.vcf.gz &
tabix -p vcf results/invicta1.vcf.gz

cp tmp/richteri1_littleb.txt results
cp tmp/richteri2_littleb.txt results
cp tmp/invicta1_littleb.txt results

```

## Measure pi

```sh

ln -sfr \
  input/gng20170922.busco4.complete.extended1000bp.coordinates-fixed.merged.bed \
  tmp/regions.bed


conda activate popgenome


for p in richteri1 richteri2 invicta1; do

  bcftools view -H results/${p}.vcf.gz \
   | cut -f 1 | uniq > tmp/${p}.scaffolds

done


for p in richteri1 richteri2 invicta1; do

  Rscript vcf_region.R \
    --input_vcf_index tmp/${p}.scaffolds \
    --bed tmp/regions.bed \
    --output_bed tmp/${p}.bed

done


for p in richteri1 richteri2 invicta1; do

  Rscript diversity.R \
    --input_vcf results/${p}.vcf.gz \
    --bed tmp/${p}.bed \
    --output_table tmp/${p}.csv \
    1> tmp/${p}.csv.out \
    2> tmp/${p}.csv.err

done

```

We then collate all the different tables

```r

# Join all populations
inv <- read.csv("tmp/invicta1.csv", header = TRUE)
r1  <- read.csv("tmp/richteri1.csv", header = TRUE)
r2  <- read.csv("tmp/richteri2.csv", header = TRUE)

# A problem is that some regions are not present in the VCF (pi = 0)
busco_regions <- read.table("tmp/regions.bed")
colnames(busco_regions) <- c("scaffold", "start", "end")

add_pi_zero <- function(df, region_df) {
  regions_in_df <- paste(df$scaffold, df$start, df$end)
  bed_regions   <- paste(region_df$scaffold, region_df$start, region_df$end)

  # Find the regions that are not present in the measurements
  absent_regions <- ! bed_regions %in% regions_in_df

  absent_regions <- cbind(
    region_df[absent_regions, ],
    matrix(NA, nrow = sum(absent_regions), ncol = ncol(df) - 3)
  )

  colnames(absent_regions) <- colnames(df)

  # Make some stats 0
  absent_regions$nuc.diversity.within <- 0
  absent_regions$n.segregating.sites  <- 0
  absent_regions$nuc.diversity        <- 0

  # Return full table
    return(rbind(df, absent_regions))
}

inv <- add_pi_zero(inv, busco_regions)
r1  <- add_pi_zero(r1, busco_regions)
r2  <- add_pi_zero(r2, busco_regions)

#  Join all populations
inv$population <- "invicta_1"
r1$population <- "richteri_1"
r2$population <- "richteri_2"

div_measurements <- rbind(inv, r1, r2)

# Add information on supergene region
regions <- read.table("input/chr16.scf-supergene-assignment.txt")
colnames(regions)[1:3] <- c("chr", "start", "end")

library(GenomicRanges)

regions_gr              <- GRanges(regions)
div_measurements_gr     <- div_measurements
div_measurements_gr$chr <- div_measurements_gr$scaffold
div_measurements_gr     <- GRanges(div_measurements_gr)

reg_div_gr <- findOverlaps(regions_gr, div_measurements_gr)


div_measurements$region <- ""
div_measurements$region[subjectHits(reg_div_gr)] <- regions_gr$V6[queryHits(reg_div_gr)]

div_measurements$region[div_measurements$region == "nonrec"] <- "non_recombining"
div_measurements$region[div_measurements$region == "rec"] <- "chr16_recombining"
div_measurements$region[div_measurements$region == ""] <- "chr1_chr15"

write.csv(div_measurements,
          file  = "results/diversity_per_region.csv",
          quote = FALSE,
          row.names = FALSE)

```

Now, we need to aggregate population, per region.


```r

library(tidyverse)

read.csv("results/diversity_per_region.csv", header = TRUE) -> div_measurements

# Deal with populations
#  richteri_1 is tagged as "1" in the masnucript
#  richteri_2 is tagged as "2" in the masnucript
# 34 input/richteri1_littleb.txt
# 10 input/richteri2_littleb.txt
# 62 input/invicta1_littleb.txt

# Calculate aggregate value per population
div_measurements %>%
  group_by(population, region) %>%
  summarise(number_of_regions = n(),
            total_region_size = sum(end - start),
            nucleotide_diversity = sum(nuc.diversity.within) / sum(end - start)) %>%
  mutate(population = recode(population,
                             invicta_1 = "S. invicta/macdonaghi\n(n = 62 males)",
                             richteri_1 = "S. richteri 1\n(n = 10 males)",
                             richteri_2 = "S. richteri 2\n(n = 34 males)")) ->
  div_per_population

font_size <- 11

# Plot aggregate value as a bar plot
div_per_population %>%
  # mutate(class = paste(population, region)) %>%
  filter(region == "non_recombining") %>%
  mutate(region = recode(region,
                         chr1_chr15 = "Chromosomes 1-15",
                         non_recombining = "Social chromosome supergene region")) %>%
  ggplot(aes(x = population, y = nucleotide_diversity,
            label = sprintf("%0.4f", round(nucleotide_diversity, 4)))) +
    geom_bar(stat = "identity", width = 0.5) +
    geom_text(position = position_stack(vjust = 1.05)) +
    theme_bw() +
    theme(axis.text.x = element_text(#angle = 0,
                                     #hjust = 1, vjust = 1,
                                     size = font_size,
                                   colour = "black"),
          axis.text.y = element_text(size = font_size, colour = "black"),
          axis.title.x = element_text(size = font_size + 2,
                                      margin = margin(b = 5, t = 10)),
          axis.title.y = element_text(size = font_size + 2)) +
    xlab("Population") +
    ylab("Nucleotide diversity") -> pi_plot

ggsave(pi_plot, file="results/nucleotide_diversity_nr.pdf",
width = 6, height = 5)

# Plot aggregate value as a bar plot
# Compares non-recombining region with chromosomes 1 to 15
div_per_population %>%
  # mutate(class = paste(population, region)) %>%
  filter(region != "chr16_recombining") %>%
  mutate(region = recode(region,
                         chr1_chr15 = "Chromosomes 1-15",
                         non_recombining = "Social chromosome supergene region")) %>%
  ggplot(aes(x = population, y = nucleotide_diversity)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab("Population") +
    ylab("Nucleotide diversity") -> pi_plot

ggsave(pi_plot, file="results/nucleotide_diversity.pdf")

# Writes out aggregate values to table
div_per_population %>%
  write.csv(file = "results/diversity_per_region_summary.csv",
            quote = FALSE,
            row.names = FALSE)

div_per_population %>%
  filter(region == "non_recombining") %>%
  write.csv(file = "results/diversity_per_region_summary_nr.csv",
            quote = FALSE,
            row.names = FALSE)

# Calculates percentage differences
div_per_population %>%
  pivot_wider(names_from = population, values_from = nucleotide_diversity) -> div_per_population_diffs

div_per_population_diffs$percent_1 <- round(100 * div_per_population_diffs$`S. richteri 1` / div_per_population_diffs$`S. invicta/macdonaghi`)
div_per_population_diffs$percent_2 <- round(100 * div_per_population_diffs$`S. richteri 2` / div_per_population_diffs$`S. invicta/macdonaghi`)
# 68 for S. richteri 1 and 100 for S. richteri 2


## Histogram
div_measurements %>%
  # mutate(class = paste(population, region)) %>%
  mutate(region = recode(region,
                         chr1_chr15 = "Chromosomes 1-15",
                         non_recombining = "Social chromosome supergene region",
                         chr16_recombining = "Social chromosome recombining region")) %>%
  ggplot(aes(x = nuc.diversity)) +
    geom_histogram(binwidth = 0.0005) +
    facet_grid(rows = vars(population), cols = vars(region)) +
    theme_bw() +
    xlab("Nucleotide diversity") -> pi_hist

ggsave(pi_hist, file="results/nucleotide_diversity_histogram.pdf")


div_measurements %>%
  mutate(region_size = end - start) %>%
  select(1:3, region_size, region, population, nuc.diversity) %>%
  pivot_wider(names_from = population, values_from = nuc.diversity) ->
  div_per_population_wide

div_per_population_wide %>%
  filter(region == "non_recombining") -> div_per_population_wide_nr

wilcox.test(div_per_population_wide_nr$invicta_1,
            div_per_population_wide_nr$richteri_1,
            paired=TRUE)

wilcox.test(div_per_population_wide_nr$invicta_1,
            div_per_population_wide_nr$richteri_2,
            paired=TRUE)

```

## Rename some samples in the VCF

A set of samples (AR164-bigB, SRR7028251, SRR7028249, SRR7028250,  SRR7028257 and SRR7028261) should be relabelled to littleb, and their colony social status should be changed from "-m". Samples AR66 and AR142 need to be removed.

```sh

cat samplesrename
# AR164-6-bigB-p AR164-6-littleb-p
# SRR7028251_AL-149-bigB-m SRR7028251_AL-149-littleb-p
# SRR7028249_AL-145-bigB-m SRR7028249_AL-145-littleb-p
# SRR7028250_AL-141-bigB-m SRR7028250_AL-141-littleb-p
# SRR7028257_AL-158-bigB-m SRR7028257_AL-158-littleb-p
# SRR7028261_AL-139-bigB-m SRR7028261_AL-139-littleb-p
# AR18-1-littleb-p AR187-1-littleb-p
# AR187-1-littleb-p AR18-1-littleb-p

bcftools reheader -s samplesrename tmp/gt.vcf.gz > tmp/gt_renamed.vcf.gz
tabix -fp vcf tmp/gt_renamed.vcf.gz

```

To make sure it worked, make list of samples in old and new VCF:

```sh

bcftools query -l tmp/gt.vcf.gz > tmp/original_samples
bcftools query -l tmp/gt_renamed.vcf.gz > tmp/renamed_samples

```

Then use R to compare those list of samples to the expected in `samplesrename`.

```r

# list of samples in old and new VCF
original_samples <- as.character(read.table("tmp/original_samples")$V1)
renamed_samples  <- as.character(read.table("tmp/renamed_samples")$V1)
# Table used to rename samples
rename_table     <- read.table("samplesrename")

# Are the samples changed to the expected value in the new VCF?
to_change <- match(rename_table$V1, original_samples)
stopifnot(renamed_samples[to_change] == rename_table$V2)

# Is any other sample changed? (it should not be)
stopifnot(original_samples[-to_change] == renamed_samples[-to_change])

```

Then move to results.

```sh

mv tmp/gt_renamed.vcf.gz results/gt.vcf.gz
mv tmp/gt_renamed.vcf.gz.tbi results/gt.vcf.gz.tbi

mv tmp/gt.vcf.gz tmp/gt_original_sample_names.vcf.gz
mv tmp/gt.vcf.gz.tbi tmp/gt_original_sample_names.vcf.gz.tbi

```
