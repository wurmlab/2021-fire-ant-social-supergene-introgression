
### 1. get gene lists

# get list of BUSCO genes to use here for reanalysis
cat genes.models.overview.tsv | grep -v "_AIC" | grep "converged" | grep "chr16nr" | awk '$4 >=10' > revision.chr16nr.BUSCOs.lst
cat genes.models.overview.tsv | grep -v "_AIC" | grep "converged" | grep "chr1-15" | awk '$4 >=10' > revision.chr1-15.BUSCOs.lst

cat revision.chr1-15.BUSCOs.lst | cut -f1 > NatCommRevision/our.samples.chr1-15.lst
cat NatCommRevision/our.samples.chr1-15.lst | wc -l
#1625
cat revision.chr16nr.BUSCOs.lst | cut -f1 > NatCommRevision/our.samples.chr16nr.lst
cat NatCommRevision/our.samples.chr16nr.lst | wc -l
#97

## some BUSCO genes were wrongly asigned to chr1-15/16nr (n=19, 3 of which are >= 10 SNPs)
## exclude both from this analysis: busco.genes.wrongly.assigned.exclude.lst
mkdir run5
cat our.samples.chr1-15.lst | grep -w -v -f busco.addChr16nr.exlcudeFromChr1-15.10SNPs.lst | grep -w -v -f busco.exlcudeFromChr16nr.10SNPs.lst  > run5/our.samples.chr1-15.lst
cat our.samples.chr16nr.lst | grep -w -v -f busco.exlcudeFromChr16nr.10SNPs.lst > run5/our.samples.chr16nr.lst

cat run5/our.samples.chr1-15.lst > run5/common.buscos.chr1-15.lst
#1623
cat run5/our.samples.chr16nr.lst > run5/common.buscos.chr16nr.lst
#96

cat run5/our.samples.chr1-15.lst run5/our.samples.chr16nr.lst > run5/our.samples.combined.lst











### 2. make alignments (direct alignments from per-sample fasta which all have the same length)

## run parallel
mkdir run5/dating.msa.chr16nr
tac run5/common.buscos.chr16nr.lst | parallel -j 100 "echo {}; cat ../data/fa.per.busco/{}.fasta.samples.fa > run5/dating.msa.chr16nr/{}.all.fasta"
mkdir run5/dating.msa.chr1-15
tac run5/common.buscos.chr1-15.lst | parallel -j 100 "echo {}; cat ../data/fa.per.busco/{}.fasta.samples.fa > run5/dating.msa.chr1-15/{}.all.fasta"

cat dating.msa.chr16nr/8673at7399.all.fasta | grep '>' | tr -d '>' | sed -r "s#\.[0-9]+{1,5}at[0-9]{4}##g" > samples2keepInMSA.lst

## Chr16nr
cat run5/common.buscos.chr16nr.lst | grep -v -w -f busco.addChr16nr.exlcudeFromChr1-15.lst > run5/common.buscos.chr16nr.without-added-buscos.lst
mkdir run5/dating.msa-pruned.chr16nr
INPUTLIST="run5/common.buscos.chr16nr.without-added-buscos.lst"
N=$(cat $INPUTLIST | wc -l) && echo $N
for (( i = 1 ; i < $N+1 ; i++))
 do
BUSCOGENE=$(cat $INPUTLIST | sed -n $i'p')
echo $i". "$BUSCOGENE
cat run5/dating.msa.chr16nr/$BUSCOGENE.all.fasta | sed -r "s#\.[0-9]+{1,5}at[0-9]{4}##g" > run5/dating.msa-pruned.chr16nr/$BUSCOGENE.fa
done
mkdir run5/dating.msaconcat.chr16nr
pxcat -s run5/dating.msa-pruned.chr16nr/*.fa -p run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.partitions -o run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.supermatrix.fa


## Chr1-15
mkdir run5/dating.msa-pruned.chr1-15
INPUTLIST="run5/common.buscos.chr1-15.lst"
N=$(cat $INPUTLIST | wc -l) && echo $N
for (( i = 1 ; i < $N+1 ; i++))
 do
BUSCOGENE=$(cat $INPUTLIST | sed -n $i'p')
echo $i". "$BUSCOGENE
cat run5/dating.msa.chr1-15/$BUSCOGENE.all.fasta | sed -r "s#\.[0-9]+{1,5}at[0-9]{4}##g" > run5/dating.msa-pruned.chr1-15/$BUSCOGENE.fa
done

mkdir run5/dating.msaconcat.chr1-15
pxcat -s run5/dating.msa-pruned.chr1-15/*.fa -p run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.partitions -o run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.supermatrix.fa


## stats 

pxlssq -s run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.supermatrix.fa > run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.supermatrix.fa.stats
cat run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.supermatrix.fa.stats

Number of sequences: 368
Sequence length: 3305494
--------- Nucl TABLE ----------
Nucl        Total   Proportion
   A    369129108     0.303455
   C    218919190      0.17997
   G    220993309     0.181675
   T    371157209     0.305122
   -            0            0
   ?     36222976    0.0297783
 G+C    439912499     0.361645

pxlssq -s run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.supermatrix.fa > run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.supermatrix.fa.stats
cat run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.supermatrix.fa.stats

Number of sequences: 368
Sequence length: 55585672
--------- Nucl TABLE ----------
Nucl        Total   Proportion
   A  -2147483648      1.85713
   C  -2147483648        1.099
   G  -2147483648      1.09946
   T  -2147483648      1.86643
   -            0            0
   ?   1057001984     0.322684
 G+C  -2147483648      2.19846














### 3. partitions


### format partitions
# get models for each BUSCO (here we use the model for the coding region)
#INPUTLIST="run5/common.buscos.chr16nr.lst"
INPUTLIST="run5/common.buscos.chr16nr.without-added-buscos.lst"
rm -f run5/busco.dating.chr16nr.models.busco.for.partitions.txt
N=$(cat $INPUTLIST | wc -l) && echo $N
for (( i = 1 ; i < $N+1 ; i++))
 do
BUSCOGENE=$(cat $INPUTLIST | sed -n $i'p')
echo $i". "$BUSCOGENE
MODEL=$(cat ../genes.models.overview.tsv | grep -w $BUSCOGENE | cut -f2)
echo $MODEL":"$BUSCOGENE"," >> run5/busco.dating.chr16nr.models.busco.for.partitions.txt
done


cat run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.partitions | sed "s/run5\/dating.msa-pruned.chr16nr\///g" | sed "s/\.fa//g" > run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.partitions.formatted1.raxml
#
cat run5/busco.dating.chr16nr.models.busco.for.partitions.txt | tr -d ' ' > run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.models.busco.for.partitions.txt
#
echo "#nexus" > run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.partitions.formatted2.nex
echo "begin sets;" >> run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.partitions.formatted2.nex
cat run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.partitions.formatted1.raxml | sed "s/DNA, /    charset /g" | sed -r "s/$/;/g" >> run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.partitions.formatted2.nex
echo -n "    charpartition mine = " >> run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.partitions.formatted2.nex
cat run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.models.busco.for.partitions.txt | tr '\n' ' ' | sed "s/, $/;\nend;/g" >> run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.partitions.formatted2.nex
echo "done models chr16nr" >> run5/done.txt
cat run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.partitions.formatted2.nex





## chr1-15
# get models for each BUSCO (here we use the model for the coding region)
INPUTLIST="run5/common.buscos.chr1-15.lst"
rm -f run5/busco.dating.chr1-15.models.busco.for.partitions.txt
N=$(cat $INPUTLIST | wc -l) && echo $N
for (( i = 1 ; i < $N+1 ; i++))
 do
BUSCOGENE=$(cat $INPUTLIST | sed -n $i'p')
echo $i". "$BUSCOGENE
MODEL=$(cat ../genes.models.overview.tsv | grep -w $BUSCOGENE | cut -f2)
echo $MODEL":"$BUSCOGENE"," >> run5/busco.dating.chr1-15.models.busco.for.partitions.txt
done

## Chr1-15 general
cat run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.partitions | sed "s/^AA/DNA/g" | sed "s/run5\/dating.msa-pruned.chr1-15.pxclsq.filtered\///g"  | sed "s/run5\/dating.msa-pruned.chr1-15.pxclsq\///g" | sed "s/run5\/dating.msa-pruned.pxclsq.chr1-15\///g" | sed "s/run5\/dating.msa-pruned.chr1-15\///g" | sed "s/\.fa//g" > run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.partitions.formatted1.raxml


cat run5/busco.dating.chr1-15.models.busco.for.partitions.txt | tr -d ' ' > run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.models.busco.for.partitions.txt
#

#
echo "#nexus" > run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.partitions.formatted2.nex
echo "begin sets;" >> run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.partitions.formatted2.nex
cat run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.partitions.formatted1.raxml | sed "s/DNA, /    charset /g" | sed -r "s/$/;/g" >> run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.partitions.formatted2.nex
echo -n "    charpartition mine = " >> run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.partitions.formatted2.nex
cat run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.models.busco.for.partitions.txt | tr '\n' ' ' | sed "s/, $/;\nend;/g" >> run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.partitions.formatted2.nex
cat run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.partitions.formatted2.nex








### 4. phylogeny

## nano node.ages.inferred.txt
gem-1-bigB-m-majorityallele,CGI1-1-bigB-m -3.91
Lad4-1-bigB-m-majorityallele,CGI1-1-bigB-m -2.67
Par2-5-bigB-m,CGI1-1-bigB-m -1.76
Par2-5-bigB-m,SRR9008150_sae-46-Reg-bigB -1.26
SRR9008217_int-122-Car-bigB,CGI1-1-bigB-m -1.2
SRR9008215_meg-128-BZ-bigB,CGI1-1-bigB-m -1.09
SRR9008215_meg-128-BZ-bigB,SRR9008232_xAdR-134-Bra-bigB -0.97
SRR9008200_ric-88-Ros-bigB,CGI1-1-bigB-m -1.01
SRR9008142_mac-145-Col-bigB,CGI1-1-bigB-m -0.67
SRR9008253_inv-218-Cor-bigB,CGI1-1-bigB-m -0.51

## nano node.ages.inferred.subset2.txt
## only outgroups to avoid forcing ages over the SB-Sb divergence
gem-1-bigB-m-majorityallele,CGI1-1-bigB-m -3.91
Lad4-1-bigB-m-majorityallele,CGI1-1-bigB-m -2.67
Par2-5-bigB-m,CGI1-1-bigB-m -1.76
SRR9008215_meg-128-BZ-bigB,CGI1-1-bigB-m -1.09


## Chr16nr
TREE="../data/Trees_in_newick_format/Astral.10SNP-genes.supergene.nwk"
ALIGNMENT="run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.supermatrix.fa"
PARTITIONFILE="run5/dating.msaconcat.chr16nr/busco.dating.chr16nr.partitions.formatted2.nex"
DATES="node.ages.inferred.subset2.txt"
mkdir run5/2021.09.09.chr16nr.3.datingsubset2
iqtree2 -T 50 -te $TREE -p $PARTITIONFILE -s $ALIGNMENT --date $DATES --date-ci 100 --date-tip 0 --prefix run5/2021.09.09.chr16nr.3.datingsubset2/2021.09.09.chr16nr.3.datingsubset2


##chr1-15 Fig1 tree with full set of node ages
ALIGNMENT="run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.supermatrix.fa"
PARTITIONFILE="run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.partitions.formatted2.nex"
DATES="node.ages.inferred.txt"
##TREE="../data/Trees_in_newick_format/Astral.10SNP-genes.chr1-15.nwk"
TREE="../data/Trees_in_newick_format/Astral.10SNP-genes.chr1-15.collapsed-branches-0.05-Fig1E.nwk"
mkdir run5/2021.09.09.chr1-15
iqtree2 -T 90 -te $TREE -p $PARTITIONFILE -s $ALIGNMENT --date $DATES --date-ci 100 --date-tip 0 --prefix run5/2021.09.09.chr1-15
echo "redating chr1-15 fig1 redating with inferred node ages and 10SNP collapsed Fig1 with full mapping alignment done" >> ~/run5.rerun.done.txt


##chr1-15 Fig1 tree with full set of node ages
ALIGNMENT="run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.supermatrix.fa"
PARTITIONFILE="run5/dating.msaconcat.chr1-15/busco.dating.chr1-15.partitions.formatted2.nex"
DATES="node.ages.inferred.subset2.txt"
TREE="../data/Trees_in_newick_format/Astral.10SNP-genes.chr1-15.collapsed-branches-0.05-Fig1E.nwk"
mkdir run5/2021.09.14.chr1-15.4.datingsubset2
iqtree2 -T 90 -te $TREE -p $PARTITIONFILE -s $ALIGNMENT --date $DATES --date-ci 100 --date-tip 0 --prefix run5/2021.09.14.chr1-15.4.datingsubset2/2021.09.14.chr1-15









