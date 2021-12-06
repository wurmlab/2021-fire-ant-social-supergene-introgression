library(GenomicRanges)
library(VariantAnnotation)
library(GenomicFeatures)

# VCF file
vcf_path <- "tmp/biallelic_snps.vcf.gz"

sample_table <- read.table("results/samples", header = TRUE)

invicta_bb  <- sample_table$sample[sample_table$species_variant == "invicta_bb"]
invicta_lb  <- sample_table$sample[sample_table$species_variant == "invicta_lb"]
richteri_bb <- sample_table$sample[sample_table$species_variant == "richteri_bb"]
richteri_lb <- sample_table$sample[sample_table$species_variant == "richteri_lb"]
saevissima  <- sample_table$sample[sample_table$species_variant == "saevissima_outgroup"]
all_samples <- c(invicta_bb, invicta_lb, richteri_bb, richteri_lb, saevissima)

dir.create("results/by_scaff", recursive=TRUE)

# ---------------------------------------------------------------------------- #
# Get supergene/non supergene region
# ---------------------------------------------------------------------------- #

regions <- read.table("input/linkage_map_supergene.txt", header = TRUE)
regions_gr <- regions[, c("scaffold", "scaffold_start", "scaffold_end", "orientation", "chr", "region")]
colnames(regions_gr) <- c(
  "chr", "start", "end", "orientation", "chr_name", "region"
)
regions_gr <- GRanges(regions_gr)

# ---------------------------------------------------------------------------- #
# Read genome fasta and gff file, to predict coding changes
# ---------------------------------------------------------------------------- #
txdb <- makeTxDbFromGFF("input/gnga.gff", format = "gff3")
gnga <- readDNAStringSet("input/ref.fa")
names(gnga) <- gsub(" .*", "", names(gnga))

# ---------------------------------------------------------------------------- #
# Read VCF bit by bit
# ---------------------------------------------------------------------------- #
scaffs_in_vcf <- headerTabix("tmp/biallelic_snps.vcf.gz")$seqnames
for (i in 1:100) {

  if (as.character(seqnames(regions_gr[i])) %in% scaffs_in_vcf) {
    # ---------------------------------------------------------------------------- #
    # Read VCF (VariantAnnotation)
    # ---------------------------------------------------------------------------- #

    param <- ScanVcfParam(which   = regions_gr[i],
                          geno    = "GT",
                          samples = all_samples)

    vcf <- readVcf(vcf_path, row.names=TRUE, param=param)

    # ---------------------------------------------------------------------------- #
    # Remove weird sites
    # ---------------------------------------------------------------------------- #
    acceptable_gt <- c("0|0", "1|1", ".|.")
    bi_allelic    <- apply(geno(vcf)$GT, 1, function(gt_row) all(gt_row %in% acceptable_gt))
    vcf           <- vcf[bi_allelic,]

    # keep only sites with a single alternative - too slow
    #  now done with bcftools in advance
    # vcf <- vcf[sapply(rowRanges(vcf)$ALT, length) == 1, ]

    # ---------------------------------------------------------------------------- #
    # Recode GT
    # ---------------------------------------------------------------------------- #

    geno(vcf)$GT[geno(vcf)$GT == "0|0"] <- "0"
    geno(vcf)$GT[geno(vcf)$GT == "1|1"] <- "1"
    geno(vcf)$GT[geno(vcf)$GT == ".|."] <- NA

    # ---------------------------------------------------------------------------- #
    # No missing
    # ---------------------------------------------------------------------------- #
    # Analysis will only work if most individuals in all populations carry the SNP

    missing_per_row <- function(gt_row, max_missing) {
      # max_missing should be a proportion
      stopifnot(max_missing < 1, max_missing >= 0)
      max_missing <- max_missing * length(gt_row)

      # retuns TRUE is the number of missing genopypes is smaller than max_missing
      return(sum(is.na(gt_row)) <  max_missing)
    }

    missing_per_pop <- function(gt, pop, max_missing) {
      gt <- gt[,pop]
      # returns vector of TRUE or FALSE
      not_missing <- apply(gt, 1, function(row) missing_per_row(row, max_missing))
    }

    max_missing <- 0.20 # that's 5.6 in 28

    # filter by missing
    not_missing_pop <- missing_per_pop(geno(vcf)$GT, invicta_bb, max_missing) &
      missing_per_pop(geno(vcf)$GT, invicta_lb, max_missing) &
      missing_per_pop(geno(vcf)$GT, richteri_bb, max_missing) &
      missing_per_pop(geno(vcf)$GT, richteri_lb, max_missing) &
      missing_per_pop(geno(vcf)$GT, saevissima, max_missing)

    # sum(not_missing_pop)
    # 11979 in 12222
    vcf <- vcf[not_missing_pop, ]

    # ---------------------------------------------------------------------------- #
    # Needs to be polymorphic
    # ---------------------------------------------------------------------------- #

    alternative_frequency <- function(gt_row) {

      # Deal with missing data
      if (all(is.na(gt_row))) {
        print(gt_row)
        stop("Error: all samples are NA")
      }
      gt_row    <- gt_row[!is.na(gt_row)]

      # Calculate frquency of alternative (0) allele
      aaf <- sum(gt_row == 1) / length(gt_row)

    }

    af <- apply(geno(vcf)$GT, 1, function(row) alternative_frequency(row))

    vcf <- vcf[(af != 0) & (af != 1)]

    # ---------------------------------------------------------------------------- #
    # Remove sites that are polymorphic in the outgroup species
    # ---------------------------------------------------------------------------- #

    af_by_pop <- function(gt, pop) {
      apply(gt[,pop], 1, function(row) alternative_frequency(row))
    }

    af_outgroup <- af_by_pop(gt = geno(vcf)$GT, pop = saevissima)

    # sum((af_outgroup == 0) & (af_outgroup == 1))
    # [1] 3836

    vcf         <- vcf[(af_outgroup == 0) | (af_outgroup == 1), ]
    af_outgroup <- af_outgroup[(af_outgroup == 0) | (af_outgroup == 1)]

    # ---------------------------------------------------------------------------- #
    # Predict coding
    # ---------------------------------------------------------------------------- #
    if (nrow(vcf) > 0) {
      code_prediction <- predictCoding(vcf, txdb, seqSource = gnga)

      # Remove SNPs that have a coding effect in multiple transcripts
      duplicated_SNPs <- names(code_prediction)[duplicated(names(code_prediction))]

      vcf <- vcf[!rownames(vcf) %in% duplicated_SNPs]

      code_prediction <- code_prediction[!names(code_prediction) %in% duplicated_SNPs]
    }

    if (nrow(vcf) > 0) {

      # ---------------------------------------------------------------------------- #
      # Recode SNPs as numeric
      # ---------------------------------------------------------------------------- #

      gt_data <- rowRanges(vcf)
      gt      <- t(apply(geno(vcf)$GT, 1, as.numeric))
      colnames(gt) <- colnames(geno(vcf)$GT)

      # ---------------------------------------------------------------------------- #
      # Recode SNPs
      # ---------------------------------------------------------------------------- #
      #   0 is the allele carried by the outgroup
      #   1 is the other allele

      recode_by_outgroup <- function(gt_row, outgroup) {

        no_missing <- gt_row[!is.na(gt_row)]

        if (!any(no_missing %in% c(0,1))) {
          stop("The only allowed genotypes are 1 and 0")
        }

        if (all(no_missing == no_missing[1])) {
          stop("Site needs to be polymorphic")
        }

        outgroup_gt <- gt_row[outgroup]
        if (all(is.na(outgroup_gt))) {
          print(outgroup_gt)
          stop("Error: all outgroup samples are NA")
        }

        outgroup_gt <- outgroup_gt[!is.na(outgroup_gt)]

        if (any(outgroup_gt != outgroup_gt[1])) {
          print(outgroup_gt)
          stop("Error: outgroup cannot be polymorphic")
        }

        # The following is necessary in case there is missing data

        # If the first element of outgroup is the alternative allele (1),
            # Then 1 -> 0 and 0 -> 1
        if (outgroup_gt[1] == 1) {

          gt_row <- c(1,0)[match(gt_row,c(0,1))]

        }

        return(gt_row)

      }

      recoded_gt <- t(apply(gt, 1, function(row) recode_by_outgroup(row, outgroup = saevissima)))
      colnames(recoded_gt) <- colnames(gt)

      gt <- recoded_gt

      # ---------------------------------------------------------------------------- #
      # Allele frequency of each population
      # ---------------------------------------------------------------------------- #

      allele_frequencies_ma <- cbind(
        invicta_bb = af_by_pop(gt = gt, pop = invicta_bb),
        invicta_lb = af_by_pop(gt = gt, pop = invicta_lb),
        richteri_bb = af_by_pop(gt = gt, pop = richteri_bb),
        richteri_lb = af_by_pop(gt = gt, pop = richteri_lb)
      )

      # ---------------------------------------------------------------------------- #
      # Private to SB
      # ---------------------------------------------------------------------------- #

      private_alleles <- function(freq_matrix, pop) {
        pop_freq   <- freq_matrix[, pop, drop = FALSE]
        other_freq <- freq_matrix[, !colnames(freq_matrix) %in% pop, drop = FALSE]

        # Present in population 1 and no other population
        #    The sum of frequencies will be greater than 0 if the allele is present
        #  in at least one of the other populations
        (apply(pop_freq, 1, sum) > 0) & (apply(other_freq, 1, sum) == 0)

      }

      private_to_invicta_bb  <- private_alleles(allele_frequencies_ma, pop="invicta_bb")
      private_to_richteri_bb <- private_alleles(allele_frequencies_ma, pop="richteri_bb")

      # ---------------------------------------------------------------------------- #
      # Shared between clades
      # ---------------------------------------------------------------------------- #
      shared_between_clades <- function(freq_matrix, pop1, pop2) {
        pop1_freq  <- freq_matrix[, pop1, drop = FALSE]
        pop2_freqs <- freq_matrix[, pop2, drop = FALSE]

        # Present in population 1 and at least one of the other populations
        #    The sum of frequencies will be greater than 0 if the allele is present
        #  in at least one of the other populations
        (apply(pop1_freq, 1, sum) > 0) & (apply(pop2_freqs, 1, sum) > 0)
      }

      shared_lb_invicta_bb  <- shared_between_clades(
                                    allele_frequencies_ma,
                                    pop1 = "invicta_bb",
                                    pop2 = c("invicta_lb", "richteri_lb")
                                  ) &
            private_alleles(allele_frequencies_ma, c("invicta_bb", "invicta_lb", "richteri_lb"))

      shared_lb_richteri_bb <- shared_between_clades(
                                    allele_frequencies_ma,
                                    pop1 = "richteri_bb",
                                    pop2 = c("invicta_lb", "richteri_lb")
                                  ) &
          private_alleles(allele_frequencies_ma, c("richteri_bb", "invicta_lb", "richteri_lb"))

      # ---------------------------------------------------------------------------- #
      # Prepare end data frame
      # ---------------------------------------------------------------------------- #
      df_to_write <- as.data.frame(rowRanges(vcf))

      df_to_write$paramRangeID <- NULL
      df_to_write$ALT <- sapply(unlist(df_to_write$ALT), as.character)

      code_prediction <- as.data.frame(code_prediction)
      code_prediction <- code_prediction[, -c(1:9)]

      tables_match <- match(rownames(df_to_write), rownames(code_prediction))
      df_to_write  <- cbind(df_to_write, code_prediction[tables_match,])

      # Add frequencies
      # Recalculate for outgroup because VCF may have been filtered
      # But use the genotypes in VCF and not the recoded GT
      df_to_write$af_saevissima_outgroup <- af_by_pop(gt = geno(vcf)$GT,
                                                      pop = saevissima)
      df_to_write$derived_allele <- ifelse(df_to_write$af_saevissima_outgroup == 0,
                                           "REF", "ALT")

      colnames(allele_frequencies_ma) <- paste0("daf_", colnames(allele_frequencies_ma))
      df_to_write <- cbind(df_to_write, allele_frequencies_ma)

      # Add measurement of each private/shared
      df_to_write$private_to_invicta_bb  <- private_to_invicta_bb
      df_to_write$private_to_richteri_bb <- private_to_richteri_bb

      df_to_write$shared_lb_invicta_bb  <- shared_lb_invicta_bb
      df_to_write$shared_lb_richteri_bb <- shared_lb_richteri_bb

      # Add region
      df_to_write$chr    <- regions_gr$chr_name[i]
      df_to_write$region <- regions_gr$region[i]

      # ---------------------------------------------------------------------------- #
      out_path <- paste("results/by_scaff/private_and_shared_per_site_", i, ".csv", sep = "")

      write.csv(df_to_write, file = out_path,
                row.names=FALSE, quote=FALSE)
      # ---------------------------------------------------------------------------- #
    }
  }
}
