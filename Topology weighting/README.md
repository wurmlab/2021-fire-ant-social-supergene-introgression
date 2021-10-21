# Topology weighting with Twisst

We use topology weighting ([Martin and Van Belleghem 2017]( https://doi.org/10.1534/genetics.116.194720)) on windows of BUSCO genes to test wether the evidence for introgression between *S. invicta* and *S. richteri* is shared across the supergene.

We have already created non-overlapping windows of four concatenates genes for the supergene and made a tree with R AxML for each window. This is the directory `supergene`.

Below, we do the same for genes mapped to chromosome 1 and the remaining parts of chromosome 16. We add the trees made for the supergene region (`supergene`), and use Twisst to analyse them.  

## Subset genes mapping to chr1 and chr16

```sh

ln -sf /PATH/TO/busco_genes_linkage_map.csv input

ln -sf /PATH/TO/linkage_map_supergene.txt input/gngs_linkage_map.txt

```

```r

genes <- read.csv("input/busco_genes_linkage_map.csv")
genes_to_remove <- as.character(read.table("genes_to_remove")$V1)

genes_subset1 <- genes$chr == "chr1"
genes_subset1 <- genes[genes_subset1, ]

genes_chr16 <- genes$chr == "chr16" & genes$BUSCO_region == "chr16r"
genes_chr16 <- genes[genes_chr16, ]

genes_subset2 <- genes_chr16[genes_chr16$chr_start < 10e6,]

genes_subset3 <- genes_chr16[genes_chr16$chr_start > 10e6,]

stopifnot(!any(genes_chr16$BUSCO_gene_ID_Hymenoptera_ODB10 %in% genes_to_remove))

stopifnot(order(genes_subset1$chr_start) == 1:nrow(genes_subset1))
stopifnot(order(genes_subset2$chr_start) == 1:nrow(genes_subset2))
stopifnot(order(genes_subset3$chr_start) == 1:nrow(genes_subset3))

genes_subset1 <- genes_subset1[,"BUSCO_gene_ID_Hymenoptera_ODB10", drop=FALSE]
genes_subset2 <- genes_subset2[,"BUSCO_gene_ID_Hymenoptera_ODB10", drop=FALSE]
genes_subset3 <- genes_subset3[,"BUSCO_gene_ID_Hymenoptera_ODB10", drop=FALSE]

write.table(genes_subset1, file="tmp/genes_chr1", col.names=FALSE, quote=FALSE, row.names=FALSE)
write.table(genes_subset2, file="tmp/genes_chr16A", col.names=FALSE, quote=FALSE, row.names=FALSE)
write.table(genes_subset3, file="tmp/genes_chr16B", col.names=FALSE, quote=FALSE, row.names=FALSE)

```

## Concatenate genes

Get the gene alignments:

```sh

ln -sf /PATH/TO/ALIGNMENTS input

```

```sh

cut -f 1 -d ' ' results/twisst_samples_low_missingness | sort > tmp/twisst_samples
wc -l tmp/twisst_samples
# 267 tmp/twisst_samples

# Parse the coordinates file
module load samtools/1.10
module load bcftools/1.10.2
module load ruby
module load seqtk/1.2

mkdir -p tmp/gene_sequences

# In the following, we
#   - read a fasta file (one entry per line)
#   - removes gene name from fasta headers
#   - correct sample that was wrongly named in previous analyses
#   - use seqtk to subset to the relevant samples
#   - make sure we end up with the correct number of genes

# Run again, for all genes
cat tmp/genes_chr1 tmp/genes_chr16A tmp/genes_chr16B >> tmp/all_genes
while read p; do
  cat input/per.gene.phylo.input/$p/$p.fasta.allsamples.fa \
    | ruby -pe 'gsub(/\..*/, "")' \
    | ruby -pe 'gsub(/SRR7028261_AL-139-bigB-m/, "SRR7028261_AL-139-littleb-p")' \
    | ruby -pe 'gsub(/SRR7028251_AL-149-bigB-m/, "SRR7028251_AL-149-littleb-p")' \
    | seqtk subseq - tmp/twisst_samples > tmp/gene_sequences/$p.fa
done < tmp/all_genes

```

Then we concatenate the samples for chr1:

```sh

module load phyx/1.1

wc -l tmp/genes_chr1
# 471
for i in {1..471..4}; do
  sh concatenate_genes.sh chr1 $i
done

wc -l tmp/genes_chr16A
for i in {1..168..4}; do
  sh concatenate_genes.sh chr16A $i
done

wc -l tmp/genes_chr16B
for i in {1..8..4}; do
  sh concatenate_genes.sh chr16B $i
done

```

## Build a tree for each window

For chr16_B.

```sh

for i in {1..8..4}; do
  echo $i >> tmp/chr16_B_list
done

module load parallel
mkdir -p results/concatenated_trees
cat tmp/chr16_B_list \
  | parallel -j 2 bash raxml.sh chr16B {}

```

For chr1:

```sh

for i in {1..471..4}; do
  echo $i >> tmp/chr1_list
done

module load parallel
mkdir -p results/concatenated_trees

cat tmp/chr1_list \
  | parallel -j 35 bash raxml.sh chr1 {}

```

For chr16_A.

```sh

for i in {1..168..4}; do
  echo $i >> tmp/chr16_A_list
done

cat tmp/chr16_A_list \
  | parallel -j 36 bash raxml.sh chr16A {}

```

The last window only has 3 genes, so we remove it:

```sh

mv results/concatenated_trees/chr1_469 results/chr1_469_incomplete_window

```

## Make a table with the location of each tree

```sh

for i in {1..468..4}; do
  echo $i >> tmp/window_i
  for j in {0..3}; do
    sed -n "$((i+j))p" tmp/genes_chr1 >> tmp/gene${j}
  done
  echo "chr1" >> tmp/region
done

paste tmp/window_i tmp/gene0 tmp/gene1 tmp/gene2 tmp/gene3 tmp/region > tmp/genes_per_window
rm -f tmp/window_i tmp/gene0 tmp/gene1 tmp/gene2 tmp/gene3 tmp/region

for i in {1..168..4}; do
  echo $i >> tmp/window_i
  for j in {0..3}; do
    sed -n "$((i+j))p" tmp/genes_chr16A >> tmp/gene${j}
  done
  echo "chr16A" >> tmp/region
done

paste tmp/window_i tmp/gene0 tmp/gene1 tmp/gene2 tmp/gene3 tmp/region >> tmp/genes_per_window
rm -f tmp/window_i tmp/gene0 tmp/gene1 tmp/gene2 tmp/gene3 tmp/region

for i in {1..8..4}; do
  echo $i >> tmp/window_i
  for j in {0..3}; do
    sed -n "$((i+j))p" tmp/genes_chr16B >> tmp/gene${j}
  done
  echo "chr16B" >> tmp/region
done

paste tmp/window_i tmp/gene0 tmp/gene1 tmp/gene2 tmp/gene3 tmp/region >> tmp/genes_per_window
rm -f tmp/window_i tmp/gene0 tmp/gene1 tmp/gene2 tmp/gene3 tmp/region

```

In R, we add the chromosome coordinates:

```sh

windows           <- read.table("tmp/genes_per_window")
colnames(windows) <- c("window", paste("gene", 1:4, sep = ""), "region")

coords <- read.csv("input/busco_genes_linkage_map.csv", header = TRUE)

start_gene_match <- match(windows$gene1, coords$BUSCO_gene_ID_Hymenoptera_ODB10)
windows$chrom    <- coords$chr[start_gene_match]
windows$start    <- coords$chr_start[start_gene_match]

end_gene_match <- match(windows$gene4, coords$BUSCO_gene_ID_Hymenoptera_ODB10)
windows$end    <- coords$chr_end[end_gene_match]

windows$mid    <- floor(windows$start + (windows$end - windows$start) / 2)

stopifnot(order(windows$chrom[windows$chrom == "chr1"]) == order(windows$window[windows$chrom == "chr1"]))
stopifnot(order(windows$chrom[windows$chrom == "chr16A"]) == order(windows$window[windows$chrom == "chr16A"]))
# If this was not true, then it would be wrong

windows <- windows[, c(7:10,6,1:5)]

write.table(windows, row.names=FALSE, sep = "\t", quote=FALSE,
            file = "results/window_coordinates")

```

## Run twisst

```sh

mkdir -p results/2021-08-07-twisst

cat results/window_coordinates | awk '$5 == "chr1" {print}' > tmp/coord
for i in $(cut -f 6 tmp/coord); do
  cat results/concatenated_trees/chr1_${i}/chr1_${i}.raxml.bestTree \
    >> results/2021-08-07-twisst/input.tree
done

cat results/window_coordinates | awk '$5 == "chr16A" {print}' > tmp/coord
for i in $(cut -f 6 tmp/coord); do
  cat results/concatenated_trees/chr16A_${i}/chr16A_${i}.raxml.bestTree \
    >> results/2021-08-07-twisst/input.tree
done

cat results/window_coordinates | awk '$5 == "chr16B" {print}' > tmp/coord
for i in $(cut -f 6 tmp/coord); do
  cat results/concatenated_trees/chr16B_${i}/chr16B_${i}.raxml.bestTree \
    >> results/2021-08-07-twisst/input.tree
done

wc -l results/2021-08-07-twisst/input.tree
# 161 results/2021-08-07-twisst/input.tree

```

We add the trees from chromosome 16.

```sh

ln -sf supergene/results/concatenated_trees input/supergene_trees

ln -sf supergene/results/window_coordinates input/supergene_coords

cat input/supergene_coords | awk '$1 == "chr16" {print}' > tmp/coord
for i in $(cut -f 5 tmp/coord); do
  cat input/supergene_trees/${i}/window_${i}.raxml.bestTree \
    >> results/2021-08-07-twisst/input.tree
done

```

We make a file with all window coordinates:

```sh

cat input/supergene_coords \
  | awk 'BEGIN { OFS = "\t" } ; {if($1 == "chr16") print $1, $2, $3, $4, "chr16_supergene", $5, $6, $7, $8, $9}' > tmp/coord

mv results/window_coordinates tmp/window_coordinates

cat tmp/window_coordinates tmp/coord >> results/window_coordinates

wc -l results/window_coordinates
# 214 results/window_coordinates

```

We run twisst:

```sh

# Load the right environment
module load anaconda3
conda activate twisst

# Run the twisst script
python twisst/twisst.py \
  -t results/2021-08-07-twisst/input.tree \
  -w results/2021-08-07-twisst/weights_output.csv.gz \
  --outputTopos results/2021-08-07-twisst/topologies_output.trees \
  -g geminata \
  -g saevissima \
  -g pusillignis \
  -g invicta/macdonaghi_Sb \
  -g invicta/macdonaghi_SB \
  -g richteri_Sb \
  -g richteri_SB \
  --outgroup geminata \
  --method complete \
  --groupsFile results/twisst_samples_low_missingness

```

## Analyse twisst results

This is done with `twisst_interpretation.rmd`.

## Move results around

```sh

mv tmp/window results/

```
