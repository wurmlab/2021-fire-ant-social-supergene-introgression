## Workflow

1. Use bedtools to intersect the VCF and the separate GFF files corresponding to BUSCO genes
2. Use the resulting VCF files to create the FASTA files (already aligned) with `vcf2phylip.py`
3. Concatenate FASTA files using `pxcat` (from `Phyx` toolkit) 
4. Run separate `raxml-ng` analyses of concatenated SNPs from BUSCO genes located on chromosomes 1-15 and on the non-recombining region of chromosome 16

Create conda environments

```bash
conda create -n variant_analysis
conda activate variant_analysis
conda install -c bioconda bcftools
conda install -c bioconda bedtools

conda create -n genometools python=2.7
conda activate genometools
conda install -c bioconda genometools-genometools
```

GFF files on Dropbox

```
Soli.Roddy.Dec2020/gffs.tar.gz (5,851 GFFs for BUSCO genes)
Soli.Roddy.Dec2020/SolInvGFF.corrected.only.CDS.gff
Soli.Roddy.Dec2020/SolInvGFF.corrected.only.CDS.gff.gz
```

## GFF files

Check integrity of GFF files.

```bash
gt gff3 -retainids -checkids -sort -tidy input.gff3
```

This CDS file produces an error.

```bash
gt gff3 -retainids -checkids -sort -tidy SolInvGFF.corrected.only.CDS.gff.gz
# gt gff3: error: Parent "mRNA1" on line 3 in file "SolInvGFF.corrected.only.CDS.gff.gz" was not defined (via "ID=")
```

mRNA GFF passed format checks.

```
gt gff3 -retainids -checkids -sort -tidy SolInvGFF.corrected.only.mRNA.gff.gz
```

Compute descriptive statistics.

```bash
gt gff3 -sort -tidy SolInvGFF.corrected.only.CDS.gff | gt stat
```

## VCF files

VCF files on Dropbox:

```
2020-11-introgression_paper/other_info/2020-vcf/gt.vcf.gz
Soli.Roddy.Dec2020/busco4.merged.filtered.SNP.genotypes.MaxMissing25.vcf.gz
```

```
bcftools stats gt.vcf.gz
bcftools stats busco4.merged.filtered.SNP.genotypes.MaxMissing25.vcf.gz > \
busco4.merged.filtered.SNP.genotypes.MaxMissing25.vcf.bcftools.stats.txt
```

BUSCO genes separated according to chromosome of origin:

```
gng20170922.busco4.chr1-15.complete.tsv = 5479
gng20170922.busco4.chr16nr.complete.tsv = 178
gng20170922.busco4.chr16r.complete.tsv = 212

5479 + 178 + 212 = 5869
```

## Reference genome used for calling variants?

The `gt.vcf.gz` file includes the following comment line:

```
reference=/scratch/ek/solenopsis/2018/ref/gng20170922wFex.fa
db/genomic/S_invicta/2017-09-22-Si_gng20170922_eckart/gng20170922wFex.fa
```

## Analysis

```
mkdir sinv_busco_raxml && cd sinv_busco_raxml
```

The "gffs" folder contains 5,851 mRNA GFF files

## Intersect CDS and BUSCO mRNA GFF files

```bash
mkdir intersect_cds_mrna && cd intersect_cds_mrna

# Intersect the CDS and BUSCO mRNA features
# -wa	Write the original entry in A for each overlap
# -u	Write original A entry once if any overlaps found in B... report the fact at least one overlap was found in B.
ls ../gffs/*.gff | parallel -j 40 'bedtools intersect -wa -u -a ../SolInvGFF.corrected.only.CDS.gff.gz -b {} > {/.}.CDS.gff' &!

ls intersect_cds_mrna | wc -l
5851

# Concatenate the BUSCO mRNA GFF files
# These files have only "mRNA" features and lack headers
# cat *.gff > busco5851.mrna.annotation.inv.gff
# bedtools intersect -wa -u -a SolInvGFF.corrected.only.CDS.gff.gz -b busco5851.mrna.annotation.inv.gff > \
# SolInvGFF.corrected.only.CDS.busco5851.gff
```

## Intersect the BUSCO CDS and VCF files

```bash
mkdir intersect_vcf_gff && intersect_vcf_gff

# Intersect using the 5851 separate CDS GFFs
ls ../intersect_cds_mrna/*.gff | parallel -j 40 'bedtools intersect -header -wa -a ../gt.vcf.gz -b {} > {/.}.intersect.gt.vcf' &!
```

## Use `vcf2phylip` to convert SNPs in VCF format to FASTA alignments
Reference: https://doi.org/10.5281/zenodo.2540861, https://github.com/edgardomortiz/vcf2phylip
Paper citing `vcf2phylip`: https://www.nature.com/articles/s41559-019-0982-3

Example: 

```bash
python vcf2phylip.py --input input.vcf --phylip-disable --fasta > vcf2phylip.log 2>&1 &!
```

```bash
ls intersect_vcf_gff/*.vcf | parallel -j 40 'python vcf2phylip.py --input {} --phylip-disable --fasta > {.}.log 2>&1' &!
mkdir vcf2fasta
mv intersect_vcf_gff/*.log vcf2fasta/
mv intersect_vcf_gff/*.fasta vcf2fasta/
```

Three files containing the ortholog identifiers according to chromosome of origin ("Soli.Roddy.Dec2020" folder)

```
gng20170922.busco4.chr1-15.complete.tsv
gng20170922.busco4.chr16nr.complete.tsv
gng20170922.busco4.chr16r.complete.tsv
```

Get ortholog identifiers according to chromosome of origin

```bash
# chr16nr
ls vcf2fasta | grep -w "$(awk '{print $4}' gng20170922.busco4.chr16nr.complete.tsv)" > \
gng20170922.busco4.chr16nr.complete.vcf2fasta.list
mkdir vcf2fasta_chr16nr
cd vcf2fasta
while read file; do cp "$file" ../vcf2fasta_chr16nr/; done < ../gng20170922.busco4.chr16nr.complete.vcf2fasta.list
cd vcf2fasta_chr16nr
ls *.fasta | wc -l
# 210
mkdir log_files
mv *.log log_files/

# chr1-15
ls vcf2fasta | grep -w "$(awk '{print $4}' gng20170922.busco4.chr1-15.complete.tsv)" > \
gng20170922.busco4.chr1-15.complete.vcf2fasta.list
mkdir vcf2fasta_chr1-15
cd vcf2fasta
while read file; do cp "$file" ../vcf2fasta_chr1-15/; done < ../gng20170922.busco4.chr1-15.complete.vcf2fasta.list
cd vcf2fasta_chr1-15
ls *.fasta | wc -l
# 5450
mkdir log_files
mv *.log log_files/

# Total: 210 + 5450 = 5660
```

## Concatenate FASTA files

```bash
conda install -c bioconda phyx
# chr16nr
cd vcf2fasta_chr16nr
pxcat -s *.fasta -p busco.chr16nr.parts.model -o busco.chr16nr.supermatrix.fasta

# chr1-15
cd vcf2fasta_chr1-15
pxcat -s *.fasta -p busco.chr1-15.parts.model -o busco.chr1-15.supermatrix.fasta

# Print summaries of supermatrices
pxlssq -s busco.chr16nr.supermatrix.fasta
# File type: fasta
# Number of sequences: 370
# Is aligned: true
# Sequence length: 20067
# --------Nucl TABLE---------
# Nucl   Total      Proportion
# A      693097     0.093349
# C      2.89986e+06 0.390564
# G      2.89953e+06 0.39052
# T      699543     0.0942172
# -      0          0
# ?      232765     0.0313497
# G+C    5.79938e+06 0.781084
# --------Nucl TABLE---------

pxlssq -s busco.chr1-15.supermatrix.fasta
# File type: fasta
# Number of sequences: 370
# Is aligned: true
# Sequence length: 428888
# --------Nucl TABLE---------
# Nucl   Total      Proportion
# A      1.50076e+07 0.0945724
# C      6.19553e+07 0.390421
# G      6.16637e+07 0.388583
# T      1.5125e+07 0.0953125
# -      0          0
# ?      4.93696e+06 0.031111
# G+C    1.23619e+08 0.779004
# --------Nucl TABLE---------
```

## Run `raxml-ng`

```bash
conda install -c bioconda raxml-ng

# chr16nr
mkdir chr16nr_raxml
cd chr16nr_raxml
mv ../vcf2fasta_chr16nr/{busco.chr16nr.parts.model,busco.chr16nr.supermatrix.fasta,phyx.logfile} .

# Parse the supermatrix
raxml-ng --parse --msa busco.chr16nr.supermatrix.fasta --model GTR+G --prefix busco.chr16nr.supermatrix.msacheck

# chr1-15
mkdir chr1-15_raxml
cd chr1-15_raxml
mv ../vcf2fasta_chr1-15/{busco.chr1-15.parts.model,busco.chr1-15.supermatrix.fasta,phyx.logfile} .

# Parse the supermatrix
raxml-ng --parse --msa busco.chr1-15.supermatrix.fasta --model GTR+G --prefix busco.chr1-15.supermatrix.msacheck

# Submit array jobs for searches and bootstraps
conda activate variant_analysis
# chr16nr
###############################################################################
#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=1G
#$ -l h_rt=6:0:0
#$ -t 1-10
#$ -cwd
#$ -V
#$ -j y
#$ -N busco.chr16nr.supermatrix.raxml.array

cd chr16nr_raxml_array

raxml-ng --msa ../busco.chr16nr.supermatrix.fasta \
--model GTR+G \
--prefix busco.chr16nr.supermatrix.${SGE_TASK_ID} \
--threads 8 --workers auto --tree pars{1},rand{1} > \
busco.chr16nr.supermatrix.${SGE_TASK_ID}.raxml.log 2>&1
###############################################################################

###############################################################################
#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=1G
#$ -l h_rt=6:0:0
#$ -t 1-100
#$ -cwd
#$ -V
#$ -j y
#$ -N busco.chr16nr.supermatrix.bootstrap.raxml.array

cd chr16nr_raxml_bootstrap

raxml-ng --bootstrap --msa ../busco.chr16nr.supermatrix.fasta \
--model GTR+G --prefix busco.chr16nr.supermatrix.bootstrap.${SGE_TASK_ID} \
--threads 8 --workers auto --force perf_threads \
--bs-trees 1 > busco.chr16nr.supermatrix.bootstrap.${SGE_TASK_ID}.raxml.log 2>&1
###############################################################################

# Choose the tree with the best likelihood score
# chr16nr
cd chr16nr_raxml_array
grep "Final" *.log
busco.chr16nr.supermatrix.10.raxml.log:Final LogLikelihood: -349760.711747
busco.chr16nr.supermatrix.1.raxml.log:Final LogLikelihood: -349817.170912
busco.chr16nr.supermatrix.2.raxml.log:Final LogLikelihood: -349866.437269
busco.chr16nr.supermatrix.3.raxml.log:Final LogLikelihood: -349704.764196 
busco.chr16nr.supermatrix.4.raxml.log:Final LogLikelihood: -349897.541516
busco.chr16nr.supermatrix.5.raxml.log:Final LogLikelihood: -349780.597999
busco.chr16nr.supermatrix.6.raxml.log:Final LogLikelihood: -349780.597999
busco.chr16nr.supermatrix.7.raxml.log:Final LogLikelihood: -349740.392708
busco.chr16nr.supermatrix.8.raxml.log:Final LogLikelihood: -349888.403673
busco.chr16nr.supermatrix.9.raxml.log:Final LogLikelihood: -349760.711747

grep "Final" *.log | awk -F" " '{print $3}' | sort -g
Highest likelihood score = busco.chr16nr.supermatrix.3.raxml.log:Final LogLikelihood: -349704.764196

cd chr16nr_raxml_bootstrap
cat *.bootstraps > busco.chr16nr.supermatrix.raxml.100.bootstraps.trees

mkdir chr16nr_raxml_boostrap_support
cd chr16nr_raxml_boostrap_support
cp ../chr16nr_raxml_array/busco.chr16nr.supermatrix.3.raxml.bestTree .
cp ../chr16nr_raxml_bootstrap/busco.chr16nr.supermatrix.raxml.100.bootstraps.trees .

# Calcultate support values
raxml-ng --support \
--tree busco.chr16nr.supermatrix.3.raxml.bestTree \
--bs-trees busco.chr16nr.supermatrix.raxml.100.bootstraps.trees \
--bs-metric TBE

# chr1-15nr
cd chr1-15_raxml_array
grep "Final" *.log
busco.chr1-15.supermatrix.10.raxml.log:Final LogLikelihood: -8446885.990171
busco.chr1-15.supermatrix.1.raxml.log:Final LogLikelihood: -8446563.393314
busco.chr1-15.supermatrix.2.raxml.log:Final LogLikelihood: -8447314.352549
busco.chr1-15.supermatrix.3.raxml.log:Final LogLikelihood: -8447739.411846
busco.chr1-15.supermatrix.4.raxml.log:Final LogLikelihood: -8447161.086874
busco.chr1-15.supermatrix.5.raxml.log:Final LogLikelihood: -8447206.986304
busco.chr1-15.supermatrix.6.raxml.log:Final LogLikelihood: -8446978.155766
busco.chr1-15.supermatrix.7.raxml.log:Final LogLikelihood: -8447546.619656
busco.chr1-15.supermatrix.8.raxml.log:Final LogLikelihood: -8446868.076419
busco.chr1-15.supermatrix.9.raxml.log:Final LogLikelihood: -8447551.577616

grep "Final" *.log | awk -F" " '{print $3}' | sort -g
Highest likelihood score = busco.chr1-15.supermatrix.1.raxml.log:Final LogLikelihood: -8446563.393314

cd chr1-15_raxml_bootstrap
cat *.bootstraps > busco.chr1-15.supermatrix.raxml.100.bootstraps.trees

mkdir chr1-15_raxml_boostrap_support
cd chr1-15_raxml_boostrap_support
cp ../chr1-15_raxml_array/busco.chr1-15.supermatrix.1.raxml.bestTree .
cp ../chr1-15_raxml_bootstrap/busco.chr1-15.supermatrix.raxml.100.bootstraps.trees .

# Calcultate support values
raxml-ng --support \
--tree busco.chr1-15.supermatrix.1.raxml.bestTree \
--bs-trees busco.chr1-15.supermatrix.raxml.100.bootstraps.trees \
--bs-metric TBE
```

### Renamed tip labels

```markdown
| label                    | newlabel                  |
|--------------------------|---------------------------|
| SRR7028251_AL-149-bigB-m | SRR7028251_AL-149-littleb |
| SRR7028257_AL-158-bigB-m | SRR7028257_AL-158-littleb |
| SRR7028261_AL-139-bigB-m | SRR7028261_AL-139-littleb |
| SRR7028249_AL-145-bigB-m | SRR7028249_AL-145-littleb |
| SRR7028250_AL-141-bigB-m | SRR7028250_AL-141-littleb |
| AR164-6-bigB-p           | AR164-6-littleb-p         |
| AR18-1-littleb-p         | AR187-1-littleb-p         |
| AR187-1-littleb-p        | AR18-1-littleb-p          |
```

### Removed tips

Removed tips "AR142-1-bigB-m-1" and "AR66-5-bigB-p" from trees.
