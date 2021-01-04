#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
  make_option(c("-i", "--input_vcf_index"), type="character", default="in_dir",
              help="Bgzipped pseudo-diploid VCF file.", metavar="character"),
  make_option(c("-b", "--bed"), type="character", default="in_dir",
              help="Bed with regions.", metavar="character"),
  make_option(c("-o", "--output_bed"), type="character", default="out_table",
              help="Output directory", metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

# ---------------------------------------------------------------------------- #
# Inputs / parameters

input_vcf_index  <- opt$input_vcf_index
input_scaff      <- opt$bed
output_table     <- opt$output_bed

# ---------------------------------------------------------------------------- #
# Checks
if (!file.exists(input_vcf_index)) {
  stop(paste("File", input_vcf, "not found."))
}

if (file.exists(output_table)) {
  stop(paste("File", output_table, "already exists. Aborting."))
}

# ---------------------------------------------------------------------------- #


vcf_index <- read.table(input_vcf_index)
regions   <- read.table(input_scaff)

regions   <- regions[regions$V1 %in% vcf_index$V1, ]

write.table(regions, output_table, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

# ---------------------------------------------------------------------------- #
