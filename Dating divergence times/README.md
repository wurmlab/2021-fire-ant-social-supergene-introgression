# Dated phylogenies

First, we created a dated phylogeny for the species tree. We included the node separating the species _S. fugax_ from the other _Solenopsis_ species. The directory `species_tree_dating` includes all scripts used for this, including:
* _De novo_ assembly of the genomes of different fire ants
* The use of BUSCO to identify single-copy orthologs between these species
* Multiple-sequence alignment between these single-copy orthologs
* Dated phylogeny with IQ-Tree

We then created dated tree with all samples in our analysis (for the supergene and the rest of the genome), calibrated with nodes retrieved from the tree above. The scripts used in this analysis are in the directory `supergene_dating`.
