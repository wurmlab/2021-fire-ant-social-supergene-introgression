#!/usr/bin/Rscript

library(tidyverse)

read_gt <- function(path, samples, df_name) {
  gt           <- read.table(path)
  colnames(gt) <- samples
  gt$df_name   <- df_name
  return(gt)
}

paths <- c(
  "tmp/2020-07-22-filter/stats/subsampled_pre_filter_gt",
  "tmp/2020-07-22-filter/stats/subsampled_genotypes.het_to_missing_gt",
  "tmp/2020-07-22-filter/stats/subsampled_genotypes.het_filter_miss_filter_gt"
)

df_names <- c(
  "Before filtering",
  "After turning heterzygous to missing",
  "After removing sites with >10% missing genotypes"
)

samples <- c(
  "tmp/2020-07-22-filter/stats/subsampled_pre_filter.vcf.gz.samples",
  "tmp/2020-07-22-filter/stats/subsampled_genotypes.het_to_missing.vcf.gz.samples",
  "tmp/2020-07-22-filter/stats/subsampled_genotypes.het_filter_miss_filter.vcf.gz.samples"
)

samples <- lapply(samples,
  function(x) as.character(read.table(x)$V1))

gt <- lapply(1:3, function(i) read_gt(path = paths[i], samples = samples[[i]], df_name = df_names[i]))

gt <- bind_rows(gt)

# Get all GT possibilities
gather(gt, sample, gt, -df_name) %>% filter(!is.na(gt)) -> gt_table
gt_table <- table(gt_table$gt)

gt_levels <- names(gt_table)

alleles <- strsplit(gt_levels, split = "")

homozygous <- sapply(alleles, function(x) x[1]) == sapply(alleles, function(x) x[3])
homozygous <- homozygous & !grepl("\\.", gt_levels)
homozygous <- gt_levels[homozygous]

missing      <- gt_levels[grepl("\\.", gt_levels)]
heterozygous <- gt_levels[!gt_levels %in% c(missing, homozygous)]

# ---------------------------------------------------------------------------- #
# Heterozygous positions per site
het_per_site <- apply(gt, 1, function(r) sum(r %in% heterozygous))
het_per_site <- data.frame(het_per_site = het_per_site,
                            df_name     = factor(gt$df_name,
                                                 levels = df_names))

pdf("tmp/2020-07-22-filter/stats/het_per_site.pdf", height = 8, width = 12)
ggplot(het_per_site) +
  geom_histogram(aes(x = het_per_site, fill = df_name), binwidth = 1) +
  facet_grid(rows = vars(df_name)) +
  theme_bw() +
  xlab("Number of heterozygous samples per site") +
  ggtitle("Heterozygous samples per site\n(subsampled to every 100th line)") +
  xlim(0, 30)
dev.off()

# Heterozygous positions per sample
gt %>%
  gather(sample, gt, -df_name) %>%
  filter(!is.na(gt)) %>%
  filter(gt %in% heterozygous) %>%
  mutate(df_name = factor(df_name, levels = df_names)) %>%
  mutate(sample = as.factor(sample)) -> het_gt

het_gt %>%
  filter(df_name == "Before filtering") %>%
  group_by(sample) %>%
  summarise(n = n()) -> het_gt_n

het_gt_max <- max(het_gt_n$n) + 200

het_ordered_samples <- as.character(het_gt_n$sample[order(het_gt_n$n, decreasing=TRUE)])
het_gt$sample       <- factor(het_gt$sample, levels = het_ordered_samples)

pdf("tmp/2020-07-22-filter/stats/het_per_sample.pdf", width = 20, height=11)
for (i in seq(1, length(levels(het_gt$sample)), by = 100)) {
  # The for loop is used to print 100 samples per plot
  end_level <- ifelse((i + 99) < length(levels(het_gt$sample)),
                      i + 99,
                      length(levels(het_gt$sample)))

  het_gt %>%
    filter(sample %in% levels(het_gt$sample)[i:end_level]) %>%
    ggplot(aes(sample)) +
      geom_bar(aes(fill=df_name)) +
      facet_grid(rows = vars(df_name)) +
      theme_bw() +
      xlab("Sample") +
      ylab("Number of heterozygous sites") +
      ylim(0, het_gt_max) +
      ggtitle("Heterozygous genotypes per sample\n(subsampled to every 100th variant site)") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) -> het_plot
  print(het_plot)
}
dev.off()

# ---------------------------------------------------------------------------- #
# Missing positions per site
# Number of missing samples per site
missing_per_site <- apply(gt, 1, function(r) sum(grepl("\\.", r)))
missing_per_site <- data.frame(missing_per_site = missing_per_site,
                               df_name      = factor(gt$df_name,
                                                     levels = df_names))

pdf("tmp/2020-07-22-filter/stats/missing_per_site.pdf", height = 8, width = 12)
ggplot(missing_per_site) +
  geom_histogram(aes(x = missing_per_site, fill = df_name), binwidth = 5) +
  facet_grid(rows = vars(df_name)) +
  theme_bw() +
  xlab("Number of missing samples per site") +
  ggtitle("Missing genotypes per site\n(subsampled to every 100th variant site)") +
  xlim(0, 250)
dev.off()

# Number of missing genotypes per sample
gt %>%
  gather(sample, gt, -df_name) %>%
  filter(grepl("\\.", gt)) %>%
  mutate(df_name = factor(df_name, levels = df_names)) %>%
  mutate(sample = as.factor(sample)) -> missing_gt

missing_gt %>%
  filter(df_name == "Before filtering") %>%
  group_by(sample) %>%
  summarise(n = n()) -> missing_gt_n

missing_gt_max <- max(missing_gt_n$n) + 1000

ordered_samples <- as.character(missing_gt_n$sample[order(missing_gt_n$n, decreasing=TRUE)])
missing_gt$sample <- factor(missing_gt$sample, levels = ordered_samples)

pdf("tmp/2020-07-22-filter/stats/missing_per_sample.pdf", width = 20, height=11)
for (i in seq(1, length(levels(missing_gt$sample)), by = 100)) {
  # The for loop is used to print 100 samples per plot
  end_level <- ifelse((i + 99) < length(levels(missing_gt$sample)),
                      i + 99,
                      length(levels(missing_gt$sample)))

  missing_gt %>%
    filter(sample %in% levels(missing_gt$sample)[i:end_level]) %>%
    ggplot(aes(sample)) +
      geom_bar(aes(fill=df_name)) +
      facet_grid(rows = vars(df_name)) +
      theme_bw() +
      ylab("Number of missing genotypes per sample") +
      xlab("Sample") +
      ylim(0, missing_gt_max) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) -> plot
  print(plot)
}
dev.off()

# ---------------------------------------------------------------------------- #
# Summary table
gather(gt, sample, gt, -df_name) %>%
  mutate(df_name = factor(df_name, levels  = df_names))-> long_gt

long_gt$gt_class <- "gt_class"

long_gt$gt_class[long_gt$gt %in% homozygous]   <- "homozygous"
long_gt$gt_class[long_gt$gt %in% missing]      <- "missing"
long_gt$gt_class[long_gt$gt %in% heterozygous] <- "heterozygous"

long_gt %>%
  group_by(df_name, sample, gt_class) %>%
  summarise(n = n()) %>%
  mutate(per = round(n/sum(n)*100, 2)) -> gt_summary

write.csv(gt_summary,
          file      = "tmp/2020-07-22-filter/stats/filter_subset_summary_per_sample.csv",
          row.names = FALSE)

# ---------------------------------------------------------------------------- #
# ENDS
# ---------------------------------------------------------------------------- #
