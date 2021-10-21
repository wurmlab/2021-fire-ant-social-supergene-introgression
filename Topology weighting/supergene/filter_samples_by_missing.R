
# Read input
miss    <- read.table("tmp/supergene.missing", header = T)
samples <- read.csv("samples_overview.csv", header=T)[, 1:22]

colnames(miss)[1] <- "Sample.name"

miss <- merge(miss, samples, by="Sample.name")
miss$clade <- paste(miss$Species, miss$Supergene.Variant, sep = "_")

# Order dataset by the frequency of missing genotypes
miss <- miss[order(miss$F_MISS, decreasing=F), ]

# Colour by clade
chosen_cols <- c(
"#0075DC", #Blue
"#2BCE48", #Green
"#FFA405", #Orpiment
"#5EF1F2", #Sky
"#FF5005", #Zinnia
"#005C31", #Forest
"#00998F", #Turquoise
"#FF0010", #Red
"#9DCC00", #Lime
"#003380", #Navy
"#F0A3FF", #Amethyst
"#740AFF", #Violet
"#426600", #Quagmire
"#C20088", #Mallow
"#94FFB5") #Jade

cols <- chosen_cols[as.numeric(as.factor(miss$clade))]

N=nrow(miss)

pdf("results/supergene_missing_data_per_sample.pdf", width=20)
#make the barplot
barplot(miss$F_MISS,
        col = cols,
        xaxt  = "n",
        las=1,
        ylab  = "Proportion of genotypes missing in sample",
        xlab  = "Sample",
        space = 0.2,
        xlim = c(0.2, 1.2*N))

legend("left",
       legend = levels(as.factor(miss$clade)),
       fill = chosen_cols[1:length(levels(as.factor(miss$clade)))])

dev.off()

# Coloured by chosen
chosen_samples <- read.table("twisst_samples.txt", header = F)
samples_match <- match(miss$Sample.name, chosen_samples$V1)
miss$twisst_chosen <- chosen_samples$V2[samples_match]

cols <- chosen_cols[as.numeric(as.factor(miss$twisst_chosen))]

pdf("results/supergene_missing_data_per_sample_twisst_choice.pdf", width=20)
#make the barplot
barplot(miss$F_MISS,
        col = cols,
        fill = cols,
        xaxt  = "n",
        las=1,
        ylab  = "Proportion of genotypes missing in sample",
        xlab  = "Sample",
        space = 0.2,
        xlim = c(0.2, 1.2*N),
      main = "Frequency of missing genotypes per sample\ncoloured by chosen representatives of each clade")

legend("left",
       legend = levels(as.factor(miss$twisst_chosen)),
       fill = chosen_cols[1:length(levels(as.factor(miss$twisst_chosen)))])

dev.off()


richteri1 <- c("AR169-1-littleb-p", "AR55-3-littleb-p", "AR57-1-littleb-p", "AR57-O12-littleb-p", "AR56-2-littleb-p", "U93-2-littleb-p", "U95-1-littleb-p",
"U94-1-littleb-p", "U13-1-littleb-p")

miss[miss$INDV %in% richteri1,]

sum(miss$F_MISS > 0.05)
# [1] 85

100*sum(miss$F_MISS > 0.05) / nrow(miss)
# [1] 23

# List of filtered genes
library(tidyverse)

miss$source_simplified <- miss$Source

miss$source_simplified[grep("Pracana", miss$source_simplified)] <- "USA"
miss$source_simplified[grep("Wang", miss$source_simplified)] <- "USA"
miss$source_simplified[grep("Wang", miss$source_simplified)] <- "USA"

miss$source_simplified[grep("Martinez", miss$source_simplified)] <- "Wurmlab (new)"
miss$source_simplified[grep("Stolle", miss$source_simplified)] <- "Wurmlab (new)"
miss$source_simplified[grep("This study", miss$source_simplified)] <- "Wurmlab (new)"


pdf("results/supergene_missing_vs_cov.pdf")
ggplot(miss) +
  geom_point(aes(x = Genome.coverage, y = F_MISS, colour = clade, shape = source_simplified)) +
  theme_bw() +
  geom_vline(xintercept = 15) +
  geom_hline(yintercept = 0.05) +
  ylab("Proportion of genotypes missing in sample")

ggplot(miss) +
  geom_point(aes(x = Genome.coverage, y = F_MISS, colour = twisst_chosen)) +
  theme_bw() +
  geom_vline(xintercept = 15) +
  geom_hline(yintercept = 0.05) +
  ylab("Proportion of genotypes missing in sample")

dev.off()

# All individuals from our study
filter(miss, F_MISS > 0.05) %>%
  filter(!Source %in% c("Wang et al 2013", "Privman et al 2018", "Yan Z et al. 2020")) %>%
  select(-source_simplified, -twisst_chosen, -clade, -Notes) %>%
  write.csv(file="results/samples_high_missingness.csv")

filter(miss, Genome.coverage < 15, F_MISS <= 0.05) %>%
  filter(!Source %in% c("Wang et al 2013", "Privman et al 2018", "Yan Z et al. 2020")) %>%
  select(-source_simplified, -twisst_chosen, -clade, -Notes) %>%
  write.csv(file="results/samples_low_cov_but_ok_missingness.csv")

filter(miss, !Source %in% c("Wang et al 2013", "Privman et al 2018", "Yan Z et al. 2020")) %>%
  select(-source_simplified, -twisst_chosen, -clade, -Notes) %>%
  write.csv(file="results/all_our_samples.csv")

# List of samples to use in twisst
filter(miss, F_MISS < 0.05) %>%
  filter(!Species %in% c("S. interrupta", "xAdR", "S.megergates")) %>%
  mutate(Species = gsub("S\\.", "", Species)) %>%
  mutate(Species = gsub(" ", "", Species)) %>%
  mutate(Species = ifelse(Species %in% c("invicta", "macdonaghi"), "invicta/macdonaghi", Species)) %>%
  mutate(clade = ifelse(Species %in% c("invicta/macdonaghi", "richteri"),
                        paste(Species, Supergene.Variant, sep = "_"),
                        Species)) %>%
  select(Sample.name, clade) %>%
  write.table(file="results/twisst_samples_low_missingness", quote=F, row.names=F, col.names=F)
