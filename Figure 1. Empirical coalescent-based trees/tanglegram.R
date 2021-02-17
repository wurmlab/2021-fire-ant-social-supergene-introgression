libraries <- c("phylogram", "dendextend", "wesanderson", "viridis", "reticulate")
for (lib in libraries) {
  if (require(package = lib, character.only = TRUE)) {
    message("Library loaded: ", lib)
  } else {
    print("Installing")
    BiocManager::install(lib)
    library(lib, character.only = TRUE)
    message("Library loaded: ", lib)
  }
}

colour_palette <- wes_palette("Darjeeling2", 5, type = "discrete")

raw_15_tree <- read.dendrogram(file = "tmp/chr15.dnd")
raw_16_tree <- read.dendrogram(file = "tmp/chr16.dnd")


# Load the trees
raw_chr1_15tree <- read.dendrogram(file = "chr1_15.dnd")
raw_chr16tree   <- read.dendrogram(file = "chr16.dnd")


# Subtrees
invicta_r_16   <- read.dendrogram(file = "input/subtree.chr16r.invicta.nwk")
invicta_nr_16  <- read.dendrogram(file = "input/subtree.chr16nr.invicta.nwk")
invicta_15     <- read.dendrogram(file = "input/subtree.chr1-15.invicta.nwk")

richteri_r_16  <- read.dendrogram(file = "input/subtree.chr16r.richteri.nwk")
richteri_nr_16 <- read.dendrogram(file = "input/subtree.chr16nr.richteri.nwk")
richteri_15    <- read.dendrogram(file = "input/subtree.chr1-15.richteri.nwk")

#Macdonaghi samples and outgroup highlights
all_mac_text <- read.csv(file = "input/mac_samples.txt", header = FALSE)
all_mac_pattern <- highlight_pattern <- paste(all_mac_text$V1, collapse = "|")
mac_samples <- labels(raw_chr16tree)[grep(labels(raw_chr16tree), pattern = all_mac_pattern)]

out_highlight_pattern <- "SRR9008118|SRR9008166|SRR9008167|GCa3-9-littleb-p|SRR9008168"
out_highlights <- labels(raw_chr16tree)[grep(labels(raw_chr16tree), pattern = out_highlight_pattern)]

untangled_invicta_r_16 <- untangle(intersect_trees(invicta_r_16, invicta_15), method = "labels")
untangled_invicta_nr_16 <- untangle(intersect_trees(invicta_nr_16, invicta_15), method = "labels")
untangled_richteri_r_16 <- untangle(intersect_trees(richteri_r_16, richteri_15), method = "labels")
untangled_richteri_nr_16 <- untangle(intersect_trees(richteri_nr_16, richteri_15), method = "labels")

new_order_16 <- c(rev(labels(untangled_richteri_r_16[[1]])), labels(untangled_invicta_r_16[[1]]))

new_order_15 <- c(labels(untangled_invicta_r_16[[2]]), labels(untangled_richteri_r_16[[2]]))


chr16_ordered <- rotate(raw_chr16tree, which(labels(raw_chr16tree) %in% new_order_16))

chr15_ordered <- rotate(raw_chr1_15tree, which(labels(raw_chr1_15tree) %in% new_order_15))

chr16_ordered <- prune(chr16_ordered, "AR142-1-bigB-m-1")
chr15_ordered <- prune(chr15_ordered, "AR142-1-bigB-m-1")


# Create colour vectors for the label colours based on species
label_colours_16_ordered <- labels(chr16_ordered)
label_colours_15_ordered <- labels(chr15_ordered)

label_colours_16_ordered[which(label_colours_16_ordered %in% labels(invicta_r_16))] <- colour_palette[4]
label_colours_15_ordered[which(label_colours_15_ordered %in% labels(invicta_r_16))] <- colour_palette[4]

label_colours_16_ordered[which(label_colours_16_ordered %in% labels(richteri_r_16))] <- colour_palette[3]
label_colours_15_ordered[which(label_colours_15_ordered %in% labels(richteri_r_16))] <- colour_palette[3]

label_colours_16_ordered[which(label_colours_16_ordered %in% mac_samples)] <- colour_palette[4]
label_colours_15_ordered[which(label_colours_15_ordered %in% mac_samples)] <- colour_palette[4]

label_colours_16_ordered[-which(label_colours_16_ordered %in% colour_palette)] <- colour_palette[2]
label_colours_15_ordered[-which(label_colours_15_ordered %in% colour_palette)] <- colour_palette[2]


chr16_ordered_coloured <- chr16_ordered %>% set("labels_col", label_colours_16_ordered)
chr15_ordered_coloured <- chr15_ordered %>% set("labels_col", label_colours_15_ordered)

chr16_ordered_coloured_pruned <- prune(chr16_ordered_coloured, "AR142-1-bigB-m-1")
chr15_ordered_coloured_pruned <- prune(chr15_ordered_coloured, "AR142-1-bigB-m-1")



ordered_trees <- dendlist(chr15_ordered_coloured, chr16_ordered_coloured)


pdf(file = "results/ordered_tanglegram.pdf")
tanglegram(ordered_trees, main_left = "chr1-15", main_right = "chr16",
           lwd = 0.5, columns_width = c(4,4,4), dLeaf = -2, lab.cex = 0.1,
           highlight_distinct_edges = FALSE, highlight_branches_lwd = FALSE,
           axes = TRUE, rank_branches = FALSE, margin_inner = 2,
           common_subtrees_color_branches = FALSE, color_lines = label_colours_15_ordered)
dev.off()