# Recurring adaptive introgression of a supergene variant that determines social organisation

Eckart Stolle, Rodrigo Pracana, Marian K. Priebe, Gabriel Luis Hernández, Anurag Priyam, Federico López-Osorio, and Yannick Wurm

## Introduction
We used genomes of haploid male fire ants to create coalescent-based phylogenies of _Solenopsis_ species and of their supergene region. We show that the supergene variant responsible for multiple-queen colonies and associated phenotypes evolved in one species and repeatedly spread across species through introgressive hybridization events. This highlights how supergene architecture can enable a complex adaptive phenotype to permeate species boundaries.

## Samples

We genotyped a group of _Solenopsis_ samples, including _Solenopsis invicta_, _S. richteri_,  _S. macdonaghi_, _S.interrupta_, _S. megergates_, _S. pusillignis_, and _S. saevissima_. Eckart Stolle (and colleagues) collected 88 colonies from South America. We sequenced on Illumina HiSeq and stored the resulting 107 samples in the BioProject PRJNA685290. In addition, we used 260 samples from Bioprojects PRJNA396161, PRJNA542606, PRJNA421367, PRJNA450756, PRJNA182127, SRR621118, SRX021921.

## Genotyping

All steps relating to read filtering, mapping and genotype calling are presented in the directory [Genotyping](Genotyping). This directory also includes the identification of single-copy genes with BUSCO.

## Phylogenetic inferences

We inferred the phylogenetic tree of the species and of the supergene using two different approaches:
1. A coalescent-based approached detailed in the directory [Coalescence-based phylogenies](<Coalescence-based phylogenies>)
2. A concatenated alignment approach detailed in the directory [Phylogenies inferred with concatenated SNPs](<Phylogenies inferred with concatenated SNPs>)
3. A dated phylogeny approach detailed in the directory [Dating divergence times](<Dating divergence times>)
4. Quantification of introgression support using QuIBL [QuIBL analysis](<QuIBL analysis>)
5. Topology weighting analysis using Twisst [Topology weighting](<Topology weighting>)

A PCA of the supergene genotypes is presented in the directory [PCA supergene SNPs](<PCA supergene SNPs>).

## Nucleotide diversity

See directory [Nucleotide diversity (pi)](<Nucleotide diversity (pi)>).

## Bray-Curtis dissimilarity between haplotypes

See directory [Bray-Curtis dissimilarity](<Bray-Curtis dissimilarity>).
