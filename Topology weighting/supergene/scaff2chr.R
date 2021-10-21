#!/usr/bin/env Rscript
library("optparse")

option_list <- list(
  make_option(c("-s", "--scaffold"), type="character", default=NULL,
              help="Name the input table (measurements by scaffold).", metavar="character"),
  make_option(c("-m", "--map"), type="character", default=NULL,
              help="Modified AGP file (only scaffold lines).", metavar="character"),
  make_option(c("-r", "--region"), type="logical", default=FALSE,
                help="Add a region column? [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
                help="output file name [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

# ---------------------------------------------------------------------------- #
# WARNING: any information on the strand of the feature will be deleted.
# ---------------------------------------------------------------------------- #

# Input errors

if (is.null(opt$scaffold)){
  print_help(opt_parser)
  stop("You need to supply the name the input table (measurements by scaffold).\n", call.=FALSE)
}

if (!file.exists(opt$scaffold)){
  print_help(opt_parser)
  stop("You need to supply the name of a VCF file (--vcf).\n", call.=FALSE)
}

if (is.null(opt$map)){
  print_help(opt_parser)
  stop("You need to supply the name of a file (--map).\n", call.=FALSE)
}

if (!file.exists(opt$map)){
  print_help(opt_parser)
  stop("You need to supply the name of a file (--map).\n", call.=FALSE)
}

# ---------------------------------------------------------------------------- #

# Test data

# gr      <- GRanges(Rle(c("scaff1","scaff1","scaff2")), IRanges(c(1,5,10), width=10))
# gr$mea1 <- 1:length(gr)
# gr$mea2 <- 5:(4+length(gr))
#
# map_gr  <- GRanges(Rle(c("scaff1","scaff2")), IRanges(c(5,1), width=500))
# map_gr$chr         <- c("chr2", "chr2")
# map_gr$chr_start   <- c(5000, 6000)
# map_gr$chr_end     <- c(5000, 6000) + 500 - 1
# map_gr$orientation <- c("+", "-")

# ---------------------------------------------------------------------------- #
# Input 1: table with measurements
meas_table <- read.table(opt$scaffold, header=TRUE)
expected_colnames <- c("scaffold", "start", "end")

if (!all(expected_colnames %in% colnames(meas_table))) {
  error_message <- paste("Input file (--scaffold) ",
                         "requires the following columns: ",
                         paste(expected_colnames, collapse=", "),
                         ".\n", sep="")
  stop(error_message)
}

# Input 2: table made from AGP file
agp_map <- read.table(opt$map, header=TRUE)

expected_colnames <- c("scaffold", "scaffold_start", "scaffold_end",
                       "chr", "chr_start", "chr_end", "orientation")

if (!all(expected_colnames %in% colnames(agp_map))) {
  error_message <- paste("Modified AGP file (--map) ",
                         "requires the following columns: ",
                         paste(expected_colnames, collapse=", "),
                         ".\n", sep="")
  stop(error_message)
}


if(!all(agp_map$orientation %in% c("+", "-", "?"))) {
  stop("Orientation of scaffolds in each chromosome must be coded as '+', '-' and '?'")
}

# ---------------------------------------------------------------------------- #
if (opt$region == TRUE) {
  if (is.null(agp_map$region)) {
    stop('If region flag is set to TRUE, the map must include a column called "region"')
  }
}
# ---------------------------------------------------------------------------- #
# Turn to Genomic Ranges

# Input 1: table with measurements
library(GenomicRanges)
gr      <- GRanges(Rle(meas_table$scaffold),
                   IRanges(start = meas_table$start,
                           end   = meas_table$end))

mcols(gr) <- meas_table[,!colnames(meas_table) %in% c("scaffold", "start", "end")]

# Input 2: table made from AGP file
map_gr  <- GRanges(Rle(agp_map$scaffold),
                   IRanges(start = agp_map$scaffold_start,
                           end   = agp_map$scaffold_end))

mcols(map_gr) <- agp_map[,!(colnames(agp_map) %in%
                                c("scaffold", "scaffold_start", "scaffold_end"))]

# ---------------------------------------------------------------------------- #

# Non-mapped
gr_in_map <- findOverlaps(gr, map_gr, type="within", ignore.strand=TRUE)
nonmapped <- gr[-queryHits(gr_in_map)]

# Mapped to scaffolds in the positive or uncertain orientation
# pos_map    <- map_gr[map_gr$orientation == c("+", "?")]
gr_in_map <- findOverlaps(gr, map_gr, type="within", ignore.strand=TRUE)
#gr_in_map  <- findOverlaps(gr, pos_map, type="within", ignore.strand=TRUE)

gr_i       <- queryHits(gr_in_map)
map_i      <- subjectHits(gr_in_map)

frompos    <- start(gr[gr_i])
topos      <- end(gr[gr_i])
widths     <- width(gr[gr_i])

chr        <- as.character(map_gr[map_i]$chr)
starts     <- rep(-9, length(map_i))
ends       <- rep(-9, length(map_i))

# Negative orientation
which_neg         <- map_gr[map_i]$orientation == "-"
starts[which_neg] <- map_gr[map_i[which_neg]]$chr_start + end(map_gr[map_i[which_neg]]) - topos[which_neg]
ends[which_neg]   <- starts[which_neg] + widths[which_neg] - 1

# Positive orientation
which_pos         <- map_gr[map_i]$orientation != "-"
starts[which_pos] <- map_gr[map_i[which_pos]]$chr_start - start(map_gr[map_i[which_pos]]) + frompos[which_pos]
ends[which_pos]   <- starts[which_pos] + widths[which_pos] - 1

# Organise into data frame
chr       <- c(chr, rep('unmapped', length(nonmapped)))
chr_start <- c(starts, rep(NA, length(nonmapped)))
chr_end   <- c(ends, rep(NA, length(nonmapped)))

scaffold       <- as.character(c(seqnames(gr[gr_i]), seqnames(nonmapped)))
scaffold_start <- c(start(gr[gr_i]), start(nonmapped))
scaffold_end   <- c(end(gr[gr_i]), end(nonmapped))

scaffold_orientation <- as.character(map_gr[map_i]$orientation)
scaffold_orientation <- c(scaffold_orientation, rep(NA, length(nonmapped)))

stopifnot(starts !=-9)
stopifnot(ends   !=-9)

if (opt$region == TRUE) {
  region <- as.character(map_gr[map_i]$region)
  region <- c(region, rep(NA, length(nonmapped)))
  new_df <- data.frame(chr, chr_start, chr_end, scaffold, scaffold_start, scaffold_end, scaffold_orientation, region)
} else {
  new_df <- data.frame(chr, chr_start, chr_end, scaffold, scaffold_start, scaffold_end, scaffold_orientation)

}

metadata_gr    <- data.frame(rbind(elementMetadata(gr[gr_i]),
                                   elementMetadata(nonmapped)))

new_df <- cbind(new_df, metadata_gr)

write.table(new_df, file=opt$out, quote=FALSE, row.names=FALSE)

# ---------------------------------------------------------------------------- #
