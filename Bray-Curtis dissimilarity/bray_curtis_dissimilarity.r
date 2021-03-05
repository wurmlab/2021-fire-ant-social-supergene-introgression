#!/usr/bin/env Rscript

library(tidyverse)
library(GenomicRanges)

# Read genotype file
gt  <- read.table("input/supergene.gt")

samples <- as.character(read.table("input/all_samples")$V1)

colnames(gt) <- c("chr", "start", "REF", "ALT", samples)
gt$end       <- gt$start + 1

# ---------------------------------------------------------------------------- #
# Make sure all genotypes are mapped to the supergene
supergene_region <- read.table("input/chr16.scf-supergene-assignment.txt")
colnames(supergene_region) <- c("chr", "start", "end", "lg", "id", "region")

supergene_region <- GRanges(supergene_region)
supergene_region <- supergene_region[supergene_region$region == "nonrec"]
stopifnot(unique(queryHits(findOverlaps(GRanges(gt), supergene_region))) == 1:nrow(gt))

# ---------------------------------------------------------------------------- #
# Make matrix
gt_ma <- gt[, samples]
gt_ma <- as.matrix(gt_ma)

gt_ma[gt_ma == ".|."] <- NA

# ---------------------------------------------------------------------------- #
# bi-allelic only
bi_allelic <- apply(gt_ma, 1, function(x) all(x[!is.na(x)] %in% c("0|0", "1|1")))
gt_ma <- gt_ma[bi_allelic, ]


# ---------------------------------------------------------------------------- #
# Function to subset to specific group of individuals
position_variant <- function(gt_row) {
  gt_row <- gt_row[!is.na(gt_row)]

  if (length(gt_row) == 0) {
    is_it_variant <- FALSE
  } else {
    is_it_variant <- any(gt_row[1] != gt_row)
  }
  return(is_it_variant)

}

keep_variant_positions <- function(gt) {
  variant_positions <- apply(gt, 1, position_variant)
  return(gt[variant_positions,])
}

subset_to_samples <- function(gt, samples) {
  stopifnot(samples %in% colnames(gt))
  gt <- gt[, samples]
  gt <- keep_variant_positions(gt)
  return(gt)
}

matrix_to_long <- function(ma, group) {
  df <- as.data.frame(ma)
  df$species_a <- rownames(df)
  df <- pivot_longer(df,
                     cols = -species_a,
                     names_to = "species_b",
                     values_to = "distance")
  df <- filter(df, species_a != species_b)
  df <- cbind(Group = group, df)
  return(df)
}

# ---------------------------------------------------------------------------- #
# Only the positions variable among Sb
littleb_samples <- grep("littleb", samples, value = TRUE)

gt_ma <- subset_to_samples(gt_ma, samples = littleb_samples)

# nrow(gt_ma)
# [1] 24435

# ---------------------------------------------------------------------------- #
# Genotypes need to be numeric
gt_ma[gt_ma == "0|0"] <- "0"
gt_ma[gt_ma == "1|1"] <- "1"
gt_ma <- apply(gt_ma, 2, as.numeric)

# ---------------------------------------------------------------------------- #
# Bray-curtis dissimilarity

library(vegan)
gt_bray   <- vegdist(t(gt_ma), method="bray", na.rm = TRUE)
gt_bray_m <- as.matrix(gt_bray)

# ---------------------------------------------------------------------------- #
# Long version (Supplementary Data)

matrix_supp_table <- function(ma) {
  # Only count upper triangle
  ma[lower.tri(ma, diag = TRUE)] <- NA

  df         <- as.data.frame(ma)
  df$Sample_a <- rownames(df)
  df <- pivot_longer(df,
                     cols = -Sample_a,
                     names_to = "Sample_b",
                     values_to = "Bray_Curtis_dissimilarity")

  df <- df[!is.na(df$Dissimilarity), ]
  stopifnot(df$Sample_a != df$Sample_b)
  stopifnot(nrow(df) == ((nrow(ma)^2 - nrow(ma)) / 2))

  return(df)
}

matrix_supp_table(gt_bray_m) %>%
  write.csv("results/supp_table.csv", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------- #

# Mean distances for all introgression samples

## S. richteri 1
sri1_introgressed_pop <- as.character(read.table("input/richteri2_littleb.txt")$V1)
sri1_origin_pop       <- c("U57-1-littleb-p", "AR118-1-littleb-p")
rich_1_distances      <- gt_bray_m[sri1_introgressed_pop, sri1_origin_pop]
rich_1_distances      <- matrix_to_long(rich_1_distances, "Introgression 1 - S. richteri")

## S. richteri 2
sri2_introgressed_pop <- as.character(read.table("input/richteri1_littleb.txt")$V1)
sri2_origin_pop       <- c("AR46-6-littleb-p", "SRR7028257_AL-158-littleb-p",
                           "AR102-2-littleb-p", "AR93-1-littleb-p")
rich_2_distances      <- gt_bray_m[sri2_introgressed_pop, sri2_origin_pop]
rich_2_distances      <- matrix_to_long(rich_2_distances, "Introgression 2 - S. richteri")

## S. megergates
meg_introgressed_pop <- c("GCa3-9-littleb-p")
meg_origin_pop       <- c("GCa5-6-littleb-p", "GCa6-3-littleb-p", "GCa8-9-littleb-p")
meg_distances        <- gt_bray_m[meg_introgressed_pop, meg_origin_pop, drop = FALSE]
meg_distances        <- matrix_to_long(meg_distances, "Introgression 3 - S. megergates")

# S. interrupta
sint_introgressed_pop <- c("SRR9008118_int-124-Car-littleb",
                           "SRR9008166_int-125-Car-littleb")
sint_origin_pop       <- c("GCa5-6-littleb-p", "GCa6-3-littleb-p",
                           "GCa8-9-littleb-p", "SRR9008242_inv-224-Bra-littleb",
                           "U14-1-littleb-p")
sint_distances        <- gt_bray_m[sint_introgressed_pop, sint_origin_pop, drop = FALSE]
sint_distances        <- matrix_to_long(sint_distances, "Introgression 4 - S. interrupta")

# S. megergates 2
meg2_introgressed_pop <- c("SRR9008167_meg-127-BZ-littleb")
meg2_origin_pop       <- c("AR46-6-littleb-p", "SRR7028257_AL-158-littleb-p",
                      "AR102-2-littleb-p", "AR93-1-littleb-p") # Same as S richteri 2
meg2_distances        <- gt_bray_m[meg2_introgressed_pop, meg2_origin_pop, drop = FALSE]
meg2_distances        <- matrix_to_long(meg2_distances, "Introgression 5 - S. megergates")

# SRR9008168_xAdR-135-Pos-littleb
xadr_introgressed_pop <- c("SRR9008168_xAdR-135-Pos-littleb")
xadr_origin_pop       <- c("AR6-1-littleb-p", "AR18-1-littleb-p", "AR3-1-littleb-p")
xadr_distances        <- gt_bray_m[xadr_introgressed_pop, xadr_origin_pop, drop = FALSE]
xadr_distances        <- matrix_to_long(xadr_distances, "Introgression 6 - S. xAdR")

# ---------------------------------------------------------------------------- #
## Siblings
sibling_vec <- c(
  "AR164-1-littleb-p",
  "AR164-6-littleb-p",
  "AR186-1-littleb-p",
  "AR186-3-littleb-p",
  "AR187-1-littleb-p",
  "AR187-5-littleb-p",
  "AR57-1-littleb-p",
  "AR57-4-littleb-p",
  "AR57-1-littleb-p",
  "AR57-O12-littleb-p",
  "Mir6-1-littleb-p",
  "Mir6-6-littleb-p"
)

sibling_ma   <- matrix(sibling_vec, ncol = 2, byrow = TRUE)

sib_distance <- gt_bray_m[sibling_ma]

sibling_df <- cbind(Group = "same colony", as.data.frame(sibling_ma), sib_distance)
colnames(sibling_df) <- colnames(rich_2_distances)

# ---------------------------------------------------------------------------- #

all_distances <- rbind(sibling_df, meg_distances, rich_1_distances, rich_2_distances)

# ---------------------------------------------------------------------------- #
# Plot 1:
# Only the distribution

library("scales")
library("ggthemes")

gt_bray_df <- gt_bray_m
gt_bray_df[lower.tri(gt_bray_df, diag = FALSE)] <- NA

matrix_to_long(gt_bray_df, "all") %>%
  filter(!is.na(distance)) %>%
  ggplot() +
    geom_histogram(aes(x = distance),
                   boundary = 0,
                   binwidth = 0.025,
                   fill     = "grey",
                   colour   = "grey") +
    xlab("Bray-Curtis dissimilarity between Sb samples") +
    ylab("Count") +
    theme_bw() -> bc_plot

# ggsave(bc_plot, file="results/bc_distribution.pdf", width = 9, height = 6)

# ---------------------------------------------------------------------------- #
# Colour blind palette
cbPalette <- c("#D81B60", "#1E88E5", "#FFC107")

# Plot with average distances
all_distances %>%
  group_by(Group, species_b) %>%
  summarize(mean_distance = mean(distance)) -> bc_distance_means

# Plot with mean distances
bc_plot +
  geom_vline(data = filter(bc_distance_means, Group != "same colony"),
             aes(xintercept = mean_distance, colour = Group),
             size = 1) +
  geom_segment(data = filter(bc_distance_means, Group == "same colony"),
               aes(x = mean_distance, y = 90, xend = mean_distance, yend = 20),
                  arrow = arrow(length = unit(0.5, "cm")),
                  lineend = "round",
                  linejoin = "round",
                  colour = "#88357D",
                  size = 1) +
  scale_colour_manual(values=cbPalette) -> my_plot
ggsave(my_plot,file="results/supp_figure.pdf", width = 9, height = 6)

# ---------------------------------------------------------------------------- #
