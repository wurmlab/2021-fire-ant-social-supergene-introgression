#!/usr/bin/env Rscript

library("optparse")

option_list <- list(
  make_option(c("-i", "--input_vcf"), type="character", default="in_dir",
              help="Bgzipped pseudo-diploid VCF file.", metavar="character"),
  make_option(c("-b", "--bed"), type="character", default="in_dir",
              help="Bed with regions.", metavar="character"),
  make_option(c("-o", "--output_table"), type="character", default="out_table",
              help="Output directory", metavar="character"))

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

# ---------------------------------------------------------------------------- #
# Inputs / parameters

input_vcf        <- opt$input_vcf
input_scaff      <- opt$bed
output_table     <- opt$output_table

# ---------------------------------------------------------------------------- #
# Checks
if (!file.exists(input_vcf)) {
  stop(paste("File", input_vcf, "not found."))
}

if (file.exists(output_table)) {
  stop(paste("File", output_table, "already exists. Aborting."))
}

# ---------------------------------------------------------------------------- #

library("PopGenome")

# ---------------------------------------------------------------------------- #
# Table listing all the scaffolds present in the VCF -
#     and the size of each scaffold
regions           <- read.table(input_scaff)
colnames(regions) <- c("scaffold", "start", "end")

# ---------------------------------------------------------------------------- #
# The following for loop runs Popgenome functions separately for each scaffold

# Empty list to store results of for loop
diversity_measurements <- list()

# Load variants for each region
for (i in 1:nrow(regions)) {

  region_scaff <- as.character(regions$scaffold[i])
  region_start <- regions$start[i]
  region_end   <- regions$end[i]

  region_size <- region_end - region_start

  my_message <- paste0("\nMeasuring diversity and neutrality stats for ",
                       paste0(region_scaff, ":", region_start, "-", region_end),
                       "\n")

  write(my_message, stdout())

  # Load VCF file
  #   "approx = TRUE" required for haploid data
  variants <- readVCF(filename    = input_vcf,
                      tid         = region_scaff,
                      frompos     = region_start,
                      topos       = region_end,
                      numcols     = 10000,
                      approx      = TRUE)

  # Measure diversity and neutrality stats for the entirety of the population
  variants <- diversity.stats(object = variants)
  variants <- neutrality.stats(object = variants)

  # Store results in table
  diversity_stats  <- as.data.frame(get.diversity(variants)[[1]])
  neutrality_stats <- as.data.frame(get.neutrality(variants)[[1]])

  # Measure nucleotide diversity
  diversity_stats$nuc.diversity <- diversity_stats$nuc.diversity.within/region_size

  # Location of each window
  stopifnot(row.names(diversity_stats) == row.names(neutrality_stats))

  window_location <- data.frame(
    scaffold = region_scaff,
    start    = region_start,
    end      = region_end
  )

  # Make a table with all measurements
  diversity_measurements[[i]] <- cbind(window_location,
                                 diversity_stats,
                                 neutrality_stats)

}


# ---------------------------------------------------------------------------- #

# collapse list into a single data frame
diversity_measurements <- do.call(rbind, diversity_measurements)

# ---------------------------------------------------------------------------- #
# Save out
write.csv(diversity_measurements,
          file      = output_table,
          quote     = FALSE,
          row.names = FALSE)

# ---------------------------------------------------------------------------- #
# Final message
my_message <- paste("Finished measuring diversity and neutrality stats. ",
                    "Saved to ",
                    output_table,
                    ".", sep = "")

write(my_message, stdout())

# ---------------------------------------------------------------------------- #
