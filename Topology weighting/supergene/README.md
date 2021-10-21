# Twisst run


We run twisst on windows of BUSCO genes to test wether the evidence for introgression between *S. invicta* and *S. richteri* is shared across the supergene.

Here, we create non-overlapping windows of four concatenates genes. We make a tree with RAxML for each window. We use Twisst to analyses these trees.  

## Remove individuals with lots of missing data

We only want to use individuals that do not have too many missing genotypes.

```sh

mkdir tmp
mkdir input
mkdir results

ln -sf /PATH/TO/gt.vcf.gz input/all.vcf.gz

module load bcftools/1.10.2
module load vcftools/0.1.16

```

We first subset the VCF to positions located in the genes of interest in the supergene.

```sh

awk '$6 == "chr16nr" {print}' gene_coords.txt \
  > tmp/nr_genes

cut -f 7-9 tmp/nr_genes > tmp/nr_genes.regions

wc -l tmp/nr_genes
# 210 tmp/nr_genes
wc -l tmp/nr_genes.regions
# 210 tmp/nr_genes.regions

```

We count the number of missing sites per individual using vcftools.

```sh

# For all the VCF
vcftools --missing-indv --gzvcf input/all.vcf.gz
mv out.log tmp/missing.log
mv out.imiss tmp/vcf.missing

# For the genes in the supergene
vcftools --missing-indv \
  --bed tmp/nr_genes.regions \
  --gzvcf input/all.vcf.gz

# After filtering, kept 93373 out of a possible 2650251 Sites
mv out.log tmp/supergene_missing.log
mv out.imiss tmp/supergene.missing

```

We can plot missing data and filter the VCF this using the R script `filter_samples_by_missing.R`.

## Genomic locations of each BUSCO gene

First we get the linkage map.

```sh

ln -sf /PATH/TO/linkage_map_supergene.txt input/gngs_linkage_map.txt

```

This linkage map is for the gnGA assembly, which is not quite equal to the assembly used in this analysis. Particularly scaffold NW_011795053, which is not present.

The script `scaff2chr.R` takes an annotation and gives it a chromosome location. But it requires some non-standard files, which I prepare here.


```r

gene_colnames <- colnames(read.table("gene_coords.txt", header = TRUE))
genes <- read.table("tmp/nr_genes", header = FALSE)
colnames(genes) <- gene_colnames
colnames(genes)[7:9] <- c("scaffold", "start", "end")
colnames(genes)[6] <- "region"

write.table(genes, file = "tmp/nr_genes.for_mapping")

```

A test run of the programme shows that my linkage map does not include the start coordinates of a few scaffolds in the supergene.

```r

lm         <- read.table("input/gngs_linkage_map.txt", header = TRUE)
bad_scaff <- c("NW_011797558.1", "NW_011795053.1", "NW_011800952.1")
lm$scaffold_start[lm$scaffold %in% bad_scaff] <- 1

write.table(lm, file = "tmp/linkage_map_editted")

```

```sh

Rscript scaff2chr.R \
  -s tmp/nr_genes.for_mapping\
  --map tmp/linkage_map_editted \
  -o tmp/nr_genes_by_chromosome

grep -v unmapped tmp/nr_genes_by_chromosome \
  | sort -k 1,1 -k2,2n  > tmp/nr_genes_by_chromosome_mapped

```

## Concatenate genes

Get the gene alignments:


```sh

ln -sf /PATH/TO/per.gene.phylo.input/ input

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
cat input/per.gene.phylo.input/2474at7399/2474at7399.fasta.allsamples.fa \
  | ruby -pe 'gsub(/\..*/, "")' \
  | ruby -pe 'gsub(/SRR7028261_AL-139-bigB-m/, "SRR7028261_AL-139-littleb-p")' \
  | ruby -pe 'gsub(/SRR7028251_AL-149-bigB-m/, "SRR7028251_AL-149-littleb-p")' \
  | seqtk subseq - tmp/twisst_samples | grep ">" | ruby -pe 'gsub(/\>/, "")' | sort > a  
diff a tmp/twisst_samples

# Run again, for all genes
cut -f 8 -d " " tmp/nr_genes_by_chromosome_mapped | tail -n +2 > tmp/gene_list
while read p; do
  cat input/per.gene.phylo.input/$p/$p.fasta.allsamples.fa \
    | ruby -pe 'gsub(/\..*/, "")' \
    | ruby -pe 'gsub(/SRR7028261_AL-139-bigB-m/, "SRR7028261_AL-139-littleb-p")' \
    | ruby -pe 'gsub(/SRR7028251_AL-149-bigB-m/, "SRR7028251_AL-149-littleb-p")' \
    | seqtk subseq - tmp/twisst_samples > tmp/gene_sequences/$p.fa
done < tmp/gene_list

```

Then we concatenate the samples:

```sh

module load phyx/1.1

wc -l tmp/gene_list
# 209

mkdir -p tmp/window_gene_names
mkdir -p tmp/concatenated_alignments

for i in {1..205..4}; do
  echo $i
  # File listing the gene name
  tail -n +${i} tmp/gene_list | head -n 4 > tmp/window_gene_names/window_${i}

  # File listing the fasta files
  while read p; do
    echo tmp/gene_sequences/$p.fa >> tmp/window_gene_names/file_$i
  done < tmp/window_gene_names/window_${i}

  phyx pxcat -f tmp/window_gene_names/file_$i \
    > tmp/concatenated_alignments/window_${i}.fasta

done


```

## Build a tree for each window

```sh

for i in {1..205..4}; do
  echo $i >> tmp/window_number_list
done

module load parallel
mkdir -p results/concatenated_trees
cat tmp/window_number_list \
  | parallel -j 18 bash raxml.sh {}

```

## Make a table with the location of each tree

```sh
for i in {1..205..4}; do
  echo $i >> tmp/window_i
  for j in {0..3}; do
    sed -n "$((i+j))p" tmp/gene_list >> tmp/gene${j}
  done
done

paste tmp/window_i tmp/gene0 tmp/gene1 tmp/gene2 tmp/gene3 > tmp/genes_per_window

rm tmp/window_i tmp/gene0 tmp/gene1 tmp/gene2 tmp/gene3

```

In R, we add the chromosome coordinates:

```sh

windows <- read.table("tmp/genes_per_window")
colnames(windows) <- c("window", paste("gene", 1:4, sep = ""))

coords <- read.table("tmp/nr_genes_by_chromosome_mapped", header = TRUE)

windows$chrom <- coords$chr[match(windows$gene1,
                                  coords$BUSCO_gene_ID_Hymenoptera_ODB10)]
windows$start <- coords$chr_start[match(windows$gene1,
                                        coords$BUSCO_gene_ID_Hymenoptera_ODB10)]
windows$end   <- coords$chr_end[match(windows$gene4,
                                      coords$BUSCO_gene_ID_Hymenoptera_ODB10)]
windows$mid <- floor(windows$start + (windows$end - windows$start) / 2)

stopifnot(order(windows$start) == order(windows$window))
# If this was not true, then it would be wrong

windows <- windows[, c(6:9,1:5)]

write.table(windows, row.names=FALSE, sep = "\t", quote=FALSE,
            file = "results/window_coordinates")

```
