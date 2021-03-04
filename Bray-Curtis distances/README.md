# Identity between introgressed individuals and closest samples

## Data

We retrieved the genotype matrix for the SNPs in the non-recombining region.

```sh

mkdir -p tmp

awk '$6 == "nonrec" {print $1 "\t" $2+1 "\t" $3+1}' \
  input/chr16.scf-supergene-assignment.txt \
  > tmp/supergene.bed

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
  --regions-file tmp/supergene.bed results/gt.vcf.gz > tmp/supergene.gt
bcftools query -l results/gt.vcf.gz > tmp/all_samples

```

We measured the Bray-Curtis distances between samples in `bray_curtis_distances.r`.
