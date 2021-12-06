# Patterns of allele sharing across the tree

We tested whether the SB branch of _S. invicta_ shares more alleles with Sb than the SB branch of each of the other species:

((richteri_b, invicta_b),invicta_B),richteri_B)),saevissima)
((richteri_b, invicta_b),richteri_B),invicta_B)),saevissima)

In the following, we run the programme DFOIL for each of the four focal species. The analysis of allele sharing is in `alleles.rmd`.

An attempt to interpret DFOIL statistic is in `dfoil.rmd`.

## Process

Following the DFOIL intructions, We can represent the species tree as follows, where we expect the introgression of many alleles from P3 to P2.

```

----------------------------------
                            |    |
         |---- AdRX B       | P1 |
     |---|                  |    |
     |   |---- AdRX b       | P2 |
  |--|                      |    |
  |  |   |---- invicta b    | P3 |
  |  |---|                  |    |
  |      |---- invicta B    | P4 |
  |                         |    |
  |----------- outgroup     | O  |
                            |    |
----------------------------------

```

Pease and Hahn (2015) have developed an expansion of the ABBA-BABA D-statistic, named Dfoil, to look at 5 taxa rather than 4 taxa. This method has the advantage of being able of distinguishing the direction of introgression. We aim to use this method on each introgressed sample, to prove that there is support for introgression from invicta Sb (P3) to each of the species with introgression.  

Dfoil works directly from alignments. We can load the alignments previously made for the QuiBL analysis.

```sh

mkdir -p input tmp results
ln -sf /data/archive/archive-SBCS-WurmLab/2021-solenopsis_introgression/2021-07-27-quibl_trees/results/gene_sequences input

```

### Install dfoil dependencies

Dfoil is a python script.

```sh

git clone git@github.com:jbpease/dfoil.git

module load anaconda3

conda create -n dfoil_env python=3.6 scipy numpy matplotlib

```

### S. AdRX

First, we create a fasta alignment for each gene with the following 5 samples:

```
SRR9008232_xAdR-134-Bra-bigB
SRR9008168_xAdR-135-Pos-littleb
SRR9008163_inv-225-Mis-littleb
SRR9008275_inv-223-Cor-bigB
SRR9008158_sae-47-Bel-bigB
```

This is done with seqtk.

```sh

module load seqtk parallel anaconda3
conda activate dfoil_env

mkdir -p tmp/adrx/fasta

ls input/gene_sequences \
  | parallel "seqtk subseq input/gene_sequences/{} adrx.sample > tmp/adrx/fasta/{}"

```

Run dfoil.

```sh

## Prepare Dfoil input file
python dfoil/fasta2dfoil.py \
  tmp/adrx/fasta/*.fa \
  --names SRR9008232_xAdR-134-Bra-bigB,SRR9008168_xAdR-135-Pos-littleb,SRR9008163_inv-225-Mis-littleb,SRR9008275_inv-223-Cor-bigB,SRR9008158_sae-47-Bel-bigB \
  --out tmp/adrx/dfoil.in

## Run dfoil
python dfoil/dfoil.py --infile tmp/adrx/dfoil.in --out tmp/adrx/dfoil.out

## Move to results
mkdir results/dfoil_output
cp tmp/adrx/dfoil.out results/dfoil_output/adrx.dfoil

```

## Run with _S. richteri_


First, we create a fasta alignment for each gene with the following 5 samples:

```
SRR9008197_ric-96-Ros-bigB
SRR9008136_ric-74-Ros-littleb
SRR9008163_inv-225-Mis-littleb
SRR9008275_inv-223-Cor-bigB
SRR9008158_sae-47-Bel-bigB
```

This is done with seqtk.

```sh

module load seqtk parallel

mkdir -p tmp/richteri/fasta

ls input/gene_sequences \
  | parallel "seqtk subseq input/gene_sequences/{} richteri.sample > tmp/richteri/fasta/{}"

```

Run dfoil.

```sh

## Prepare Dfoil input file
python dfoil/fasta2dfoil.py \
  tmp/richteri/fasta/*.fa \
  --names SRR9008197_ric-96-Ros-bigB,SRR9008136_ric-74-Ros-littleb,SRR9008163_inv-225-Mis-littleb,SRR9008275_inv-223-Cor-bigB,SRR9008158_sae-47-Bel-bigB \
  --out tmp/richteri/dfoil.in

## Run dfoil
python dfoil/dfoil.py --infile tmp/richteri/dfoil.in --out tmp/richteri/dfoil.out

## Move to results
mkdir -p results/dfoil_output
cp tmp/richteri/dfoil.out results/dfoil_output/richteri.dfoil

```

## S. megergates

First, we create a fasta alignment for each gene with the following 5 samples:

```
SRR9008214_meg-130-BZ-bigB
SRR9008167_meg-127-BZ-littleb
SRR9008163_inv-225-Mis-littleb
SRR9008275_inv-223-Cor-bigB
SRR9008158_sae-47-Bel-bigB
```

This is done with seqtk.

```sh

module load seqtk parallel

mkdir -p tmp/megergates/fasta

ls input/gene_sequences \
  | parallel "seqtk subseq input/gene_sequences/{} megergates.sample > tmp/megergates/fasta/{}"

```

Run dfoil.

```sh

## Prepare Dfoil input file
python dfoil/fasta2dfoil.py \
  tmp/megergates/fasta/*.fa \
  --names SRR9008214_meg-130-BZ-bigB,SRR9008167_meg-127-BZ-littleb,SRR9008163_inv-225-Mis-littleb,SRR9008275_inv-223-Cor-bigB,SRR9008158_sae-47-Bel-bigB \
  --out tmp/megergates/dfoil.in

## Run dfoil
python dfoil/dfoil.py --infile tmp/megergates/dfoil.in --out tmp/megergates/dfoil.out

## Move to results
mkdir -p results/dfoil_output
cp tmp/megergates/dfoil.out results/dfoil_output/megergates.dfoil

```

## S. interrupta

First, we create a fasta alignment for each gene with the following 5 samples:

```
SRR9008217_int-122-Car-bigB
SRR9008166_int-125-Car-littleb
SRR9008163_inv-225-Mis-littleb
SRR9008275_inv-223-Cor-bigB
SRR9008158_sae-47-Bel-bigB
```

This is done with seqtk.

```sh

module load seqtk parallel

mkdir -p tmp/interrupta/fasta

ls input/gene_sequences \
  | parallel "seqtk subseq input/gene_sequences/{} interrupta.sample > tmp/interrupta/fasta/{}"

```

Run dfoil.

```sh

## Prepare Dfoil input file
python dfoil/fasta2dfoil.py \
  tmp/interrupta/fasta/*.fa \
  --names SRR9008217_int-122-Car-bigB,SRR9008166_int-125-Car-littleb,SRR9008163_inv-225-Mis-littleb,SRR9008275_inv-223-Cor-bigB,SRR9008158_sae-47-Bel-bigB \
  --out tmp/interrupta/dfoil.in

## Run dfoil
python dfoil/dfoil.py --infile tmp/interrupta/dfoil.in --out tmp/interrupta/dfoil.out

## Move to results
mkdir -p results/dfoil_output
cp tmp/interrupta/dfoil.out results/dfoil_output/interrupta.dfoil

```

in R:

```r

library(tidyverse)
dfoil <- read.table("results/dfoil_output/interrupta.dfoil", comment.char = '&', header = T)

dfoil %>%
  select(X.chrom, coord, ends_with("_stat")) %>%
  pivot_longer(ends_with("_stat"), names_to = "D_type", values_to = "D_value") %>%
  mutate(D_type = str_replace_all(D_type, "_stat","")) %>%
  mutate(D_type = fct_relevel(D_type, "DFO", "DIL", "DFI", "DOL")) -> d_stats

ggplot(d_stats, aes(x = D_type, y = D_value)) +
    geom_violin(fill="grey") + geom_jitter() + geom_boxplot(width=0.1) +
    geom_hline(yintercept = 0) +
    theme_bw()

```

## Move gene counts

```sh

mkdir -p results/allele_counts

cp tmp/richteri/dfoil.in results/allele_counts/richteri.counts
cp tmp/interrupta/dfoil.in results/allele_counts/interrupta.counts
cp tmp/megergates/dfoil.in results/allele_counts/megergates.counts
cp tmp/adrx/dfoil.in results/allele_counts/adrx.counts

```
