# Load libraries
load_cran_pkgs <- function(pkg) {
  new_pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new_pkg)) {
    install.packages(new_pkg, dependencies = TRUE)
  }
  sapply(pkg, require, character.only = TRUE)
}

# Install packages available on CRAN
cran_pkgs <- c(
  "BiocManager", "tidyverse", "RColorBrewer", "devtools", "here", "styler",
  "ggrepel", "ggnewscale", "ape", "clipr"
)
load_cran_pkgs(cran_pkgs)

# Create function to load Bioconductor packages
load_bioconductor_pkgs <- function(pkg) {
  new_pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new_pkg)) {
    BiocManager::install(new_pkg, version = "3.12")
  }
  sapply(pkg, require, character.only = TRUE)
}

bioconductor_pkgs <- c("ggtree")
load_bioconductor_pkgs(bioconductor_pkgs)

# Read input trees including bootstrap support values
tree_chr1to15 <- treeio::read.newick(
  here::here("chr1-15_raxml_bootstrap_support/busco.chr1-15.supermatrix.1.raxml.bestTree.raxml.support"),
  node.label = "label"
)

tree_chr16nr <- treeio::read.newick(
  here::here("chr16nr_raxml_boostrap_support/busco.chr16nr.supermatrix.3.raxml.bestTree.raxml.support"),
  node.label = "label"
)

# Read metadata
info <- read.csv(here::here("solenopsis_info_tiplabels.csv"))

# Drop tips from trees
# "AR142-1-bigB-m-1"
tree_chr1to15 <- treeio::drop.tip(tree_chr1to15, "AR142-1-bigB-m-1")
tree_chr16nr <- treeio::drop.tip(tree_chr16nr, "AR142-1-bigB-m-1")

# "AR66-5-bigB-p"
tree_chr1to15 <- treeio::drop.tip(tree_chr1to15, "AR66-5-bigB-p")
tree_chr16nr <- treeio::drop.tip(tree_chr16nr, "AR66-5-bigB-p")

# Create data frame for renaming tips
label <- tree_chr1to15$tip.label
labelsdf <- tibble::tibble(label)
labelsdf <- dplyr::mutate(labelsdf, newlabel = label)

# Remove "-m" from "little-b" new labels
labelsdf <- labelsdf %>%
  dplyr::mutate(
    newlabel = dplyr::case_when(
      label == "SRR7028251_AL-149-bigB-m" ~ "SRR7028251_AL-149-littleb",
      label == "SRR7028257_AL-158-bigB-m" ~ "SRR7028257_AL-158-littleb",
      label == "SRR7028261_AL-139-bigB-m" ~ "SRR7028261_AL-139-littleb",
      label == "SRR7028249_AL-145-bigB-m" ~ "SRR7028249_AL-145-littleb",
      label == "SRR7028250_AL-141-bigB-m" ~ "SRR7028250_AL-141-littleb",
      label == "AR164-6-bigB-p" ~ "AR164-6-littleb-p",
      label == "AR18-1-littleb-p" ~ "AR187-1-littleb-p",
      label == "AR187-1-littleb-p" ~ "AR18-1-littleb-p",
      TRUE ~ newlabel
    )
  )

# Sanity check
labelsdf %>%
  dplyr::filter(newlabel %in% c(
    "SRR7028251_AL-149-littleb",
    "SRR7028257_AL-158-littleb",
    "SRR7028261_AL-139-littleb",
    "SRR7028249_AL-145-littleb",
    "SRR7028250_AL-141-littleb",
    "AR164-6-littleb-p",
    "AR187-1-littleb-p",
    "AR18-1-littleb-p"
  )) # %>%
  # clipr::write_clip(object_type = "table")

# clipr::read_clip_tbl()

# Rename tree tips
tree_chr1to15 <- treeio::rename_taxa(
  tree = tree_chr1to15, data = labelsdf, key = label, value = newlabel
)

tree_chr16nr <- treeio::rename_taxa(
  tree = tree_chr16nr, data = labelsdf, key = label, value = newlabel
)

# Root the trees with an outgroup
ape::is.rooted(tree_chr1to15)
tree_chr1to15 <- ape::root(
  phy = tree_chr1to15,
  outgroup = "gem-1-bigB-m-majorityallele"
)

tree_chr16nr <- ape::root(
  phy = tree_chr16nr,
  outgroup = "gem-1-bigB-m-majorityallele"
)

# Compare the tip labels in tree and metadata
setdiff(info$tiplabel, tree_chr1to15$tip.label)
# intersect(tree_chr1to15$tip.label, info$tiplabel)

# Filter tips from metadata
to_remove <- c("AR142-1-bigB-m-1", "AR66-5-bigB-p")

# info <- info %>%
#   dplyr::filter(tiplabel != "AR142-1-bigB-m-1")

info <- info %>%
  dplyr::filter(!tiplabel %in% to_remove)

# Change tip labels in metadata
info <- info %>%
  dplyr::mutate(
    tiplabel = dplyr::case_when(
      tiplabel == "SRR7028251_AL-149-bigB-m" ~ "SRR7028251_AL-149-littleb",
      tiplabel == "SRR7028257_AL-158-bigB-m" ~ "SRR7028257_AL-158-littleb",
      tiplabel == "SRR7028261_AL-139-bigB-m" ~ "SRR7028261_AL-139-littleb",
      tiplabel == "SRR7028249_AL-145-bigB-m" ~ "SRR7028249_AL-145-littleb",
      tiplabel == "SRR7028250_AL-141-bigB-m" ~ "SRR7028250_AL-141-littleb",
      tiplabel == "AR164-6-bigB-p" ~ "AR164-6-littleb-p",
      tiplabel == "AR18-1-littleb-p" ~ "AR187-1-littleb-p",
      tiplabel == "AR187-1-littleb-p" ~ "AR18-1-littleb-p",
      TRUE ~ tiplabel
    )
  )

# Compare the tip labels in tree and metadata
setdiff(info$tiplabel, tree_chr1to15$tip.label)

# Plot trees with bootstrap support values
# chr1-15
p1 <- ggtree(tree_chr1to15,
  layout = "roundrect", size = 1, color = "#838383"
) %<+% info +
  xlim(0, 0.8)

# Convert the bootstrap values to percentages
p1d <- p1$data
p1d <- p1d[!p1d$isTip, ]
p1d$label <- as.numeric(p1d$label)
p1d$label <- round(p1d$label * 100)
p1d <- p1d[p1d$label > 80, ]

p1 +
  ggtree::geom_tiplab(size = 2, offset = .001, color = "#838383") +
  ggtree::geom_tippoint(aes(shape = grouping),
    size = 1.5,
    color = "black", fill = "white"
  ) +
  ggtree::geom_point2(
    data = p1d, aes(
      subset = !is.na(as.numeric(label)),
      fill = cut(label, c(0, 80, 90, 100))
    ),
    shape = 21, size = 1.5
  ) +
  scale_fill_manual(
    name = "Bootstrap support",
    values = c("grey", "black"),
    labels = c("80-90", "90-100")
  ) +
  scale_shape_manual(
    values = c(22, 23, 25),
    breaks = c("invicta", "richteri", "outgroup"),
    labels = c(
      expression(italic("S. invicta")),
      expression(italic("S. richteri")),
      "Outgroups"
    )
  ) +
  labs(shape = "Group") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, color = "#838383"),
    legend.text = element_text(size = 12, color = "#838383"),
    legend.text.align = 0
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 4)),
    fill = guide_legend(override.aes = list(size = 4))
  ) +
  ggtree::geom_treescale()

ggsave(
  filename = here::here("busco.chr1-15.5450loci.concatenated.raxml.bootstrap.pdf"),
  width = 10, height = 26
)

# chr16nr
p2 <- ggtree(tree_chr16nr,
  layout = "roundrect", size = 1, color = "#838383"
) %<+% info +
  xlim(0, 0.6)

p2 + geom_text(aes(label = node))
p2 <- ggtree::rotate(p2, 479)

# Convert the bootstrap values to percentages
p2d <- p2$data
p2d <- p2d[!p2d$isTip, ]
p2d$label <- as.numeric(p2d$label)
p2d$label <- round(p2d$label * 100)
# d$label <- round(p2d$label, digits = 2)
p2d <- p2d[p2d$label > 80, ]

p2 +
  ggtree::geom_tiplab(size = 2, offset = .001, color = "#838383") +  
  ggtree::geom_tippoint(aes(shape = grouping),
    size = 1.5,
    color = "black", fill = "white"
  ) +
  ggtree::geom_point2(
    data = p2d, aes(
      subset = !is.na(as.numeric(label)),
      fill = cut(label, c(0, 80, 90, 100))
    ),
    shape = 21, size = 1.5
  ) +
  scale_fill_manual(
    name = "Bootstrap support",
    values = c("grey", "black"),
    labels = c("80-90", "90-100")
  ) +
  scale_shape_manual(
    values = c(22, 23, 25),
    breaks = c("invicta", "richteri", "outgroup"),
    labels = c(
      expression(italic("S. invicta")),
      expression(italic("S. richteri")),
      "Outgroups"
    )
  ) +
  labs(shape = "Group") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, color = "#838383"),
    legend.text = element_text(size = 12, color = "#838383"),
    legend.text.align = 0
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 4)),
    fill = guide_legend(override.aes = list(size = 4))
  ) +
  ggtree::geom_treescale()

ggsave(
  filename = here::here("busco.chr16nr.210loci.concatenated.raxml.bootstrap.rotated.pdf"),
  width = 10, height = 26
)