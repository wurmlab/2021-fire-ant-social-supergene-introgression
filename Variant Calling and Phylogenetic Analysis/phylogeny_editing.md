

cd /scratch/ek/solenopsis/2018
FILTEROUTFOLDER="/scratch/ek/solenopsis/2018/SNPIVbusco4all/vcf-genotyping/samples-filtered"
MERGEFOLDER="$FILTEROUTFOLDER/merge"
F="$MERGEFOLDER/consensus.fasta"
mkdir -p $F
VCFfile="$MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.MaxMissing25.vcf.gz"

vcfsamplenames $VCFfile > $F/samples.lst
cat $F/samples.lst

TMPDIR="/dev/shm/VCFprocessing"
mkdir -p $TMPDIR
cp $VCFfile $VCFfile.tbi $TMPDIR/  #copy to RAMdisk for faster access



#################################################################
# first loop: make folder per gene and prep input files per gene
#################################################################

VCFfile="$TMPDIR/busco4.merged.filtered.SNP.genotypes.MaxMissing25.vcf.gz"
BUSCOGENES="ref/gng20170922.fa.busco4.chrs.completeFragmented"
INVICTAGFF="ref/SolInvGFF.corrected.only.mRNA.gff"
INVICTACDSGFF="ref/SolInvGFF.corrected.only.CDS.gff"
REF="ref/gng20170922.fa"

cp $BUSCOGENES.genelist.tsv $F


######### consensus sequences
N=$(wc -l $BUSCOGENES.genelist.tsv | cut -d' ' -f 1) && echo $N     ##5851
i=1
for (( i = 1 ; i < 34+1 ; i++))
do
echo "$i of $N"

SPECIES="Sinvicta"
GENE=$(cat $BUSCOGENES.genelist.tsv | sed -n $i'p') && echo $GENE
mkdir -p $F/$GENE
REGIONINV="$(cat $BUSCOGENES.tsv | sed -n $i'p' | cut -f1-3)" && echo $REGIONINV

#### get mRNA position based on detected BUSCO gene (via intersect of region, also save the gff entries
RNAIDINV=$(echo $REGIONINV | intersectBed -a $INVICTAGFF -b stdin | grep -w "mRNA" | head -n1 | cut -f9 | cut -d';' -f1 | cut -d'=' -f2) && echo "RNA.identified:invicta."$RNAIDINV

# test if Variable is empty (i.e. BUSCO gene not part of the annotation)
if [ -z "$RNAIDINV" ]
then
      echo "\$RNAIDINV is empty, i.e. BUSCO gene not part of the annotation"
      MRNASTART=$(echo $REGIONINV | cut -f2) #&& echo $MRNASTART
      echo $REGIONINV > $F/$GENE/$GENE.annotation.inv.CDS.gff.bed
      REGIONGFFINV=$(echo $REGIONINV | sed "s/\.1\t/.1:/g" | sed "s/\.1A\t/.1A:/g" | sed "s/\.1B\t/.1B:/g" | sed "s/\.1C\t/.1C:/g" | sed "s/\.1D\t/.1D:/g" | sed "s/\t/-/g")
      echo "invicta."$REGIONINV > $F/$GENE/$GENE.region.txt
      samtools faidx $REF $REGIONGFFINV > $TMPDIR/$REGIONGFFINV.fa
      cat $TMPDIR/$REGIONGFFINV.fa | sed "s#>#>$SPECIES.$GENE #g" > $F/$GENE/$GENE.fasta.consensus.gnGA.fa
      cat $F/samples.lst | parallel --no-notice -j 40 "bcftools consensus --sample {} --fasta-ref $TMPDIR/$REGIONGFFINV.fa $VCFfile 1> >(sed 's#>#>{}.$GENE #g' > $F/$GENE/$GENE.fasta.consensus.{}.fa) 2> >(sed 's#Applied #$GENE\t{}\t#g' | sed 's# variants##g' > $F/$GENE/$GENE.fasta.consensus.{}.SNPnumber)"
      cat $F/$GENE/$GENE.fasta.consensus.*.fa | seqtk seq -C -L 0 > $F/$GENE/$GENE.fasta.allsamples.fa
      cat $F/$GENE/$GENE.fasta.consensus.*.SNPnumber > $F/$GENE/$GENE.fasta.allsamples.SNPnumber

else
cat $INVICTAGFF | grep -v "#" | grep -m 1 "$RNAIDINV;" > $F/$GENE/$GENE.annotation.inv.gff
REGIONGFFINV=$(cat $F/$GENE/$GENE.annotation.inv.gff | head -n1 | cut -f1,2,4,5 | sed "s/\t/-/g" | sed "s/-Gnomon-/:/g") && echo $REGIONGFFINV
echo "invicta."$REGIONINV > $F/$GENE/$GENE.region.txt
echo "extract fasta sequences based on mRNA positions"
#
# extraction of consensus 0,440sec per sample
samtools faidx $REF $REGIONGFFINV > $TMPDIR/$REGIONGFFINV.fa
cat $TMPDIR/$REGIONGFFINV.fa | sed "s#>#>$SPECIES.$GENE #g" > $F/$GENE/$GENE.fasta.consensus.gnGA.fa
## for each sample, get consensus from VCFfile, redirect stderr into separate file (how many variants were applied)
## in both files modify header to include samplename and/or gene name
cat $F/samples.lst | parallel --no-notice -j 40 "bcftools consensus --sample {} --fasta-ref $TMPDIR/$REGIONGFFINV.fa $VCFfile 1> >(sed 's#>#>{}.$GENE #g' > $F/$GENE/$GENE.fasta.consensus.{}.fa) 2> >(sed 's#Applied #$GENE\t{}\t#g' | sed 's# variants##g' > $F/$GENE/$GENE.fasta.consensus.{}.SNPnumber)"
cat $F/$GENE/$GENE.fasta.consensus.*.fa | seqtk seq -C -L 0 > $F/$GENE/$GENE.fasta.allsamples.fa
cat $F/$GENE/$GENE.fasta.consensus.*.SNPnumber > $F/$GENE/$GENE.fasta.allsamples.SNPnumber

fi
echo "samples consensus fasta extracted"
#
done


################# Partition files

#### create partitions file for CDS codons

echo "create partitions file for CDS codons, as well as Raxml partitionfile with codons and introns"
N=$(wc -l $F/modeltest/missingfiles.all.lst | cut -d' ' -f 1) && echo $N
i=1
for (( i = 1 ; i < $N+1 ; i++))
do
echo "$i of $N"
GENE=$(cat $F/modeltest/missingfiles.all.lst | sed -n $i'p') && echo $GENE
#10281at7399
REGIONINV="$(cat $BUSCOGENES.tsv | grep -w $GENE | cut -f1-3)" && echo $REGIONINV
RNAIDINV=$(echo $REGIONINV | intersectBed -a $INVICTAGFF -b stdin | grep -w "mRNA" | head -n1 | cut -f9 | cut -d';' -f1 | cut -d'=' -f2) && echo "RNA.identified:invicta."$RNAIDINV
#RNA.identified:invicta.mRNA9147

# test if Variable is empty (i.e. BUSCO gene not part of the annotation)
if [ -z "$RNAIDINV" ]
then
      echo "\$RNAIDINV is empty, i.e. BUSCO gene not part of the annotation"
      MRNASTART=$(echo $REGIONINV | cut -f2) #&& echo $MRNASTART
      echo $REGIONINV > $F/$GENE/$GENE.annotation.inv.CDS.gff.bed
      cat $F/$GENE/$GENE.annotation.inv.CDS.gff.bed | awk -v MRNASTART=$MRNASTART '{FS="\t";OFS="\t"} {print $1, $2-MRNASTART+1, $3-MRNASTART+1, "CODINGorNONCODING"}' > $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed
     echo -n "DNA, GENE = " > $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODING.txt
     cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | cut -f2,3 | tr '\t' '-' | tr '\n' ',' | sed "s/,$//g" >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODING.txt
     cat $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODING.txt > $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
else

REGIONGFFINV=$(cat $F/$GENE/$GENE.annotation.inv.gff | head -n1 | cut -f1,2,4,5 | sed "s/\t/-/g" | sed "s/-Gnomon-/:/g") && echo $REGIONGFFINV

cat $INVICTACDSGFF | grep -w $RNAIDINV > $F/$GENE/$GENE.annotation.inv.CDS.gff
paste <(cat $F/$GENE/$GENE.annotation.inv.gff $F/$GENE/$GENE.annotation.inv.CDS.gff | cut -f1,4,5) <(cat $F/$GENE/$GENE.annotation.inv.gff $F/$GENE/$GENE.annotation.inv.CDS.gff | cut -f3) > $F/$GENE/$GENE.annotation.inv.CDS.gff.bed
MRNASTART=$(cat $F/$GENE/$GENE.annotation.inv.CDS.gff.bed | head -n1 | cut -f2) #&& echo $MRNASTART
cat $F/$GENE/$GENE.annotation.inv.CDS.gff.bed | awk -v MRNASTART=$MRNASTART '{FS="\t";OFS="\t"} {print $1, $2-MRNASTART+1, $3-MRNASTART+1, $4}' | awk '{FS="\t";OFS="\t"} {print $1, $2, $3, $4, $3-$2}' | awk '!(($5) < 3)' | cut -f1-4 > $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed

echo -n "DNA, CDS = " > $F/$GENE/$GENE.annotation.partitionfile.CDS.txt
cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n+2 | cut -f2,3 | tr '\t' '-' | tr '\n' ',' | sed "s/,$//g" >> $F/$GENE/$GENE.annotation.partitionfile.CDS.txt

#introns:
MRNAEND=$(cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | head -n1 | cut -f3) #&& echo $MRNAEND
CDSEND=$(cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n1 | cut -f3) #&& echo $CDSEND

MRNARELSTART=$(cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | head -n1 | cut -f2) #&& echo $MRNARELSTART
CDSSTART=$(cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | head -n2 | tail -n1 | cut -f2) #&& echo $CDSSTART

if [ $(( $MRNARELSTART-$CDSSTART )) -eq 0 ] ; then
   #echo "start of CDS = start of mRNA"
   cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n+2 | awk '{OFS="\t"; if($4!=id) {if(e!="") print chr,s,e,id; chr=$1;s=$3;id=$4;e="";} else {e=$2;print chr,s,e,id;s=$3;e="";}}END{if(e!="") print chr,s,e,id;}' | awk '{OFS="\t"}{print $1,$2+1,$3-1,"INTRON"}' > $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.introns.bed
else
   #echo "there is a 5UTR"
   cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n+2 | head -n1 | cut -f1,2 | awk '{OFS="\t"}{print $1,1,$2-1,"5UTR"}' > $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.introns.bed
cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n+2 | awk '{OFS="\t"; if($4!=id) {if(e!="") print chr,s,e,id; chr=$1;s=$3;id=$4;e="";} else {e=$2;print chr,s,e,id;s=$3;e="";}}END{if(e!="") print chr,s,e,id;}' | awk '{OFS="\t"}{print $1,$2+1,$3-1,"INTRON"}' >> $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.introns.bed
fi

if [ $(( $MRNAEND-$CDSEND )) -eq 0 ] ; then
   #echo "end of CDS = end of mRNA"
else
   #echo "mRNA extends beyong end of CDS, i.e. UTR"
   cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n1 | awk -v MRNAEND=$MRNAEND -v CDSEND=$CDSEND '{OFS="\t"}{print $1,CDSEND+1,MRNAEND,"3UTR"}' >> $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.introns.bed
fi


if [ -s $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.introns.bed ]
then
    echo -n "DNA, NONCODING = " > $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODING.txt
    cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.introns.bed | cut -f2,3 | tr '\t' '-' | tr '\n' ',' | sed "s/,$//g" >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODING.txt
    echo >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODING.txt
    echo -n "DNA, CODING = " >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODING.txt
    cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n+2 | cut -f2,3 | tr '\t' '-' | tr '\n' ',' | sed "s/,$//g" >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODING.txt
#
    echo -n "DNA, NONCODING = " > $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
    cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.introns.bed | cut -f2,3 | tr '\t' '-' | tr '\n' ',' | sed "s/,$//g" >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
    echo -n "\nDNA, CDS1 = " >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
    cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n+2 | cut -f2,3 | tr '\t' '-' | tr '\n' ',' | sed "s/,/\\\3,/g" | sed "s/,$//g" >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
    echo -n "\nDNA, CDS2 = " >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
    cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n+2 | cut -f2,3 | awk '{OFS="\t"}{print $1+1,$2}' | tr '\t' '-' | tr '\n' ',' | sed "s/,/\\\3,/g" | sed "s/,$//g" >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
    echo -n "\nDNA, CDS3 = " >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
    cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n+2 | cut -f2,3 | awk '{OFS="\t"}{print $1+2,$2}' | tr '\t' '-' | tr '\n' ',' | sed "s/,/\\\3,/g" | sed "s/,$//g" >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt

else
        echo "noncoding is empty"
echo -n "DNA, CODING = " > $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODING.txt
cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n+2 | cut -f2,3 | tr '\t' '-' | tr '\n' ',' | sed "s/,$//g" >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODING.txt
        echo -n "DNA, CDS1 = " > $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n+2 | cut -f2,3 | tr '\t' '-' | tr '\n' ',' | sed "s/,/\\\3,/g" | sed "s/,$//g" >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
echo -n "\nDNA, CDS2 = " >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n+2 | cut -f2,3 | awk '{OFS="\t"}{print $1+1,$2}' | tr '\t' '-' | tr '\n' ',' | sed "s/,/\\\3,/g" | sed "s/,$//g" >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
echo -n "\nDNA, CDS3 = " >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt
cat $F/$GENE/$GENE.annotation.inv.CDS.gff.relative.coordinates.bed | tail -n+2 | cut -f2,3 | awk '{OFS="\t"}{print $1+2,$2}' | tr '\t' '-' | tr '\n' ',' | sed "s/,/\\\3,/g" | sed "s/,$//g" >> $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODONS.txt

fi

fi

done



################################################################################
######################  modeltest per partition (Coding, Noncoding)
```sh

#F="/scratch/ek/solenopsis/2018/SNPIIIbusco/phylo"
mkdir -p modeltest
BUSCOGENES="gng20170922.fa.busco4.chrs.completeFragmented"

#CPUs=4
#cat $BUSCOGENES.genelist.tsv \
 #| parallel --no-notice -j 15 "
 #echo {} &&
   modeltest-ng \
   -i {$BUSCOGENES}/{$BUSCOGENES}.fasta.allsamples.fa \
     -d nt \
     -p $CPUs \
     -q {$BUSCOGENES}/{$BUSCOGENES}.annotation.partitionfile.NONCODING-CODING.txt \
     -T raxml \
     --rngseed 2 \
     -t fixed-ml-gtr \
     --force \
     --output modeltest/{$BUSCOGENES}.modeltest.NONCODING-CODING
     #"
#echo "modeltesting done"

```
#############################################################################


#######################  collect .aic modeltest files per gene to be used in phylogeny
# copy previous modeltest partition files
N=$(cat $BUSCOGENES.genelist.tsv | wc -l) && echo $N     ##5851
i=1
for (( i = 1 ; i < $N+1 ; i++))
do
echo "$i of $N"
GENE=$(cat $BUSCOGENES.genelist.tsv | sed -n $i'p')
mkdir -p $F/$GENE
cp /scratch/ek/solenopsis/2018/SNPIIIbusco/phylo/modeltest/$GENE.modeltest.NONCODING-CODING.part.aic $F/$GENE
done


# noncoding
      5 F81
     28 F81+I
    127 F81+I+G4
     16 GTR
   1311 GTR+G4
   2681 GTR+I
    355 GTR+I+G4
     42 HKY
    239 HKY+G4
    911 HKY+I
    111 HKY+I+G4

# coding
     62 F81
      1 F81+G4
     56 F81+I
    242 F81+I+G4
     10 GTR
    197 GTR+G4
   1362 GTR+I
   2256 GTR+I+G4
     22 HKY
     84 HKY+G4
    931 HKY+I
    611 HKY+I+G4

### 27 genes do not produce Models as outut for the partitionfile, fetch the model from the gene's logfiles

N=$(wc -l $F/modeltest/missingfiles.nc-coding.lst | cut -d' ' -f 1) && echo $N
i=1
for (( i = 1 ; i < $N+1 ; i++))
do
echo "$i of $N"
GENE=$(cat $F/modeltest/missingfiles.nc-coding.lst | sed -n $i'p') && echo $GENE
MODEL=$(cat $F/modeltest/$GENE.modeltest.NONCODING-CODING.log | grep -A2 "Best model according to AIC$" | tail -n1 | sed "s#              ##g" | cut -d":" -f2) #&& echo $MODEL
cat $F/$GENE/$GENE.annotation.partitionfile.NONCODING-CODING.txt | sed -r "s#DNA#$MODEL#g" > $F/modeltest/$GENE.modeltest.NONCODING-CODING.part.aic
cat $F/modeltest/$GENE.modeltest.NONCODING-CODING.part.aic
done


####################


####### prep gene trees for consensus

#### chr16 supergene assignment

NW_011797094.1	0	1802099	16	469	rec
NW_011794933.1B	0	4024115	16	470	rec
NW_011794943.1	0	5603	16	471	rec
NW_011803760.1	0	1482883	16	472	rec
NW_011795620.1	0	65617	16	473	nonrec
NW_011800952.1	0	104328	16	474	nonrec
NW_011803911.1	0	659464	16	475	nonrec
NW_011800721.1	0	500137	16	476	nonrec
NW_011799419.1	0	1115537	16	477	nonrec
NW_011801243.1	0	1297335	16	478	nonrec
NW_011794409.1	0	445462	16	479	nonrec
NW_011794634.1	0	526590	16	480	nonrec
NW_011795711.1	0	493224	16	481	nonrec
NW_011797558.1	0	340320	16	482	nonrec
NW_011794567.1	0	1946445	16	483	nonrec
NW_011795727.1	0	632960	16	484	nonrec
NW_011795053.1	0	1052406	16	485	nonrec
NW_011795861.1	0	4991	16	486	nonrec
NW_011795747.1	0	6167	16	487	nonrec
NW_011794844.1	0	1251831	16	488	nonrec
NW_011794623.1	0	3037	16	489	nonrec
NW_011801067.1	0	337463	16	490	nonrec
NW_011796111.1	0	560000	16	491	nonrec
NW_011796111.1	560000	1594552	16	492	rec
NW_011797002.1	0	54053	16	493	rec


paste <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep -P 'NW_011797094.1|NW_011794933.1B|NW_011794943.1|NW_011803760.1|NW_011797002.1' | cut -f3-5) <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep -P 'NW_011797094.1|NW_011794933.1B|NW_011794943.1|NW_011803760.1|NW_011797002.1' | cut -f1) <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep -P 'NW_011797094.1|NW_011794933.1B|NW_011794943.1|NW_011803760.1|NW_011797002.1' | cut -f6-) > gng20170922.busco4.chr16r.complete.tsv

paste <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep -P 'NW_011795620.1|NW_011800952.1|NW_011803911.1|NW_011800721.1|NW_011799419.1|NW_011801243.1|NW_011794409.1|NW_011794634.1|NW_011795711.1|NW_011797558.1|NW_011794567.1|NW_011795727.1|NW_011795053.1|NW_011795861.1|NW_011795747.1|NW_011794844.1|NW_011794623.1|NW_011794623.1|NW_011801067.1' | cut -f3-5) <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep -P 'NW_011795620.1|NW_011800952.1|NW_011803911.1|NW_011800721.1|NW_011799419.1|NW_011801243.1|NW_011794409.1|NW_011794634.1|NW_011795711.1|NW_011797558.1|NW_011794567.1|NW_011795727.1|NW_011795053.1|NW_011795861.1|NW_011795747.1|NW_011794844.1|NW_011794623.1|NW_011794623.1|NW_011801067.1' | cut -f1) <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep -P 'NW_011795620.1|NW_011800952.1|NW_011803911.1|NW_011800721.1|NW_011799419.1|NW_011801243.1|NW_011794409.1|NW_011794634.1|NW_011795711.1|NW_011797558.1|NW_011794567.1|NW_011795727.1|NW_011795053.1|NW_011795861.1|NW_011795747.1|NW_011794844.1|NW_011794623.1|NW_011794623.1|NW_011801067.1' | cut -f6-) > gng20170922.busco4.chr16nr.complete.tsv

# scf NW_011796111 until 560kb is non rec, rest is rec

paste <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep "NW_011796111.1" | grep -P '3930at7399|6460at7399|9619at7399|14530at7399|16796at7399|25173at7399|27665at7399|31243at7399' | cut -f3-5) <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep "NW_011796111.1" | grep -P '3930at7399|6460at7399|9619at7399|14530at7399|16796at7399|25173at7399|27665at7399|31243at7399' | cut -f1) <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep "NW_011796111.1" | grep -P '3930at7399|6460at7399|9619at7399|14530at7399|16796at7399|25173at7399|27665at7399|31243at7399' | cut -f6-) >> gng20170922.busco4.chr16r.complete.tsv

paste <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep "NW_011796111.1" | grep -P '1355at7399|2095at7399|3170at7399|3567at7399|3677at7399|4819at7399|7238at7399|7882at7399|9167at7399|14848at7399|15542at7399|16295at7399|16999at7399|19757at7399|20918at7399|21170at7399|21661at7399|23474at7399|25772at7399|29725at7399|30802at7399' | cut -f3-5) <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep "NW_011796111.1" | grep -P '1355at7399|2095at7399|3170at7399|3567at7399|3677at7399|4819at7399|7238at7399|7882at7399|9167at7399|14848at7399|15542at7399|16295at7399|16999at7399|19757at7399|20918at7399|21170at7399|21661at7399|23474at7399|25772at7399|29725at7399|30802at7399' | cut -f1) <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep "NW_011796111.1" | grep -P '1355at7399|2095at7399|3170at7399|3567at7399|3677at7399|4819at7399|7238at7399|7882at7399|9167at7399|14848at7399|15542at7399|16295at7399|16999at7399|19757at7399|20918at7399|21170at7399|21661at7399|23474at7399|25772at7399|29725at7399|30802at7399' | cut -f6-) >> gng20170922.busco4.chr16nr.complete.tsv

paste <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep -vP 'NW_011796111.1|NW_011797094.1|NW_011794933.1B|NW_011794943.1|NW_011803760.1|NW_011797002.1|NW_011795620.1|NW_011800952.1|NW_011803911.1|NW_011800721.1|NW_011799419.1|NW_011801243.1|NW_011794409.1|NW_011794634.1|NW_011795711.1|NW_011797558.1|NW_011794567.1|NW_011795727.1|NW_011795053.1|NW_011795861.1|NW_011795747.1|NW_011794844.1|NW_011794623.1|NW_011794623.1|NW_011801067.1' | cut -f3-5) <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep -vP 'NW_011796111.1|NW_011797094.1|NW_011794933.1B|NW_011794943.1|NW_011803760.1|NW_011797002.1|NW_011795620.1|NW_011800952.1|NW_011803911.1|NW_011800721.1|NW_011799419.1|NW_011801243.1|NW_011794409.1|NW_011794634.1|NW_011795711.1|NW_011797558.1|NW_011794567.1|NW_011795727.1|NW_011795053.1|NW_011795861.1|NW_011795747.1|NW_011794844.1|NW_011794623.1|NW_011794623.1|NW_011801067.1' | cut -f1) <(cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | grep -vP 'NW_011796111.1|NW_011797094.1|NW_011794933.1B|NW_011794943.1|NW_011803760.1|NW_011797002.1|NW_011795620.1|NW_011800952.1|NW_011803911.1|NW_011800721.1|NW_011799419.1|NW_011801243.1|NW_011794409.1|NW_011794634.1|NW_011795711.1|NW_011797558.1|NW_011794567.1|NW_011795727.1|NW_011795053.1|NW_011795861.1|NW_011795747.1|NW_011794844.1|NW_011794623.1|NW_011794623.1|NW_011801067.1' | cut -f6-) > gng20170922.busco4.chr1-15.complete.tsv

cat gng20170922.busco4.chr16r.complete.tsv | wc -l
cat gng20170922.busco4.chr16nr.complete.tsv | wc -l
cat gng20170922.busco4.chr1-15.complete.tsv | wc -l
#172 #178 #4
#210 #212 #0
#5229 #5479 #236
#slight differences becasue for rerun busco v 4.1.4 was used instead 4.0.5

#############################################################################

####### raxml
``` sh
#cd /scratch/ek/solenopsis/2018
F="/scratch/ek/solenopsis/2018/SNPIIIbusco/phylo"
BUSCOGENES="ref/gng20170922.fa.busco4.chrs.completeFragmented"
N=$(wc -l $BUSCOGENES.genelist.tsv | cut -d' ' -f 1) && echo $N     ##5851
OUTDIR="$F/genetrees"
mkdir -p $OUTDIR

cat /scratch/ek/solenopsis/2018/ref/gng20170922.fa.busco4.chr1-15.ComplFragm.lst \
| parallel --no-notice -j 20 \
"echo {} && $HOME/scripts/raxml.zsh \
$F/{}/{}.fasta.allsamples.fa \
$F/modeltest/{}.modeltest.NONCODING-CODING.part.aic \
$OUTDIR/{}"

# chr16nr
cat /scratch/ek/solenopsis/2018/ref/gng20170922.fa.busco4.chr16nr.ComplFragm.lst \
| parallel --no-notice -j 30 \
"echo {} && $HOME/scripts/raxml.zsh\
 $F/{}/{}.fasta.allsamples.fa \
 $F/modeltest/{}.modeltest.NONCODING-CODING.part.aic \
 $OUTDIR/{}"

# chr16r
cat /scratch/ek/solenopsis/2018/ref/gng20170922.fa.busco4.chr16r.ComplFragm.lst\
 | parallel --no-notice -j 30\
  "echo {} && $HOME/scripts/raxml.zsh \
  $F/{}/{}.fasta.allsamples.fa \
  $F/modeltest/{}.modeltest.NONCODING-CODING.part.aic \
  $OUTDIR/{}"

```

############################################################################


##
ls -1 besttrees.with.support/*.supportFBP | rev | cut -d "/" -f1 | cut -d "." -f3 | rev > finished.tees.lst
cat finished.tees.lst | wc -l
#5555

## which tree did not converge/ failed the maximum likelihood tree inference?
cd /home/ek/scratch/solenopsis/2018/SNPIVbusco4all/phylo/tmp/besttrees
cat /home/ek/scratch/solenopsis/2018/SNPIVbusco4all/vcf-genotyping/samples-filtered/merge/per.gene.phylo.input/gng20170922.fa.busco4.chrs.completeFragmented.genelist.tsv | parallel -k -j1 "test ! -f {}.raxml.bestTree && echo {} >> ../missing.trees.lst"
cat missing.trees.lst | wc -l
# 296

cat ../../../ref/gng20170922.fa.busco4.chr1-15.ComplFragm.lst
cat ../../../ref/gng20170922.fa.busco4.chr16nr.ComplFragm.lst
cat ../../../ref/gng20170922.fa.busco4.chr16r.ComplFragm.lst

cat ../../../ref/gng20170922.fa.busco4.chr16r.ComplFragm.lst | grep -wf missing.trees.lst > missing.trees.chr16r.lst

3 genes
1767at7399
13474at7399
24355at7399

cat ../../../ref/gng20170922.fa.busco4.chr16nr.ComplFragm.lst | grep -wf missing.trees.lst > missing.trees.chr16nr.lst

#7 genes
8796at7399
9291at7399
14784at7399
15093at7399
17602at7399
18251at7399
22958at7399

cat ../../../ref/gng20170922.fa.busco4.chr1-15.ComplFragm.lst | grep -wf missing.trees.lst > missing.trees.chr1-15.lst
# 286


#### filter genes with less than 10 variants
cd /scratch/ek/solenopsis/2018/SNPIVbusco4all/vcf-genotyping/samples-filtered/merge/consensus.fasta.SNPnumbers

# collect SNPnumbers for each BUSCO genes (this is bcftools consensus output stored in a file), here we exclude outgroup samples as they contain many more SNPs

for filename in ./*.fasta.allsamples.SNPnumber; do
        cat $filename | grep -vP 'Toc5|gem-1|Copa|Lad4|AR223|USP3|xAdR|_int-|SS1|_meg-|Par1|Par2|_sae-' | cut -f1,3 | bioawk -t '{ sum += $2 } END { if (NR > 0) print $1,sum / NR }' | tr ',' '.' >> /scratch/ek/solenopsis/2018/SNPIVbusco4all/SNPnumbers.tsv
done

# from total of 5851 genes (input, single copy BUSCOs)
cat SNPnumbers.tsv | bioawk -t ' $2 > 10 ' > genes.with.more.than.10.SNPs.lst # 1938 genes

#make gene lists of 10-SNP-genes
cat /scratch/ek/solenopsis/2018/SNPIVbusco4all/genes.with.more.than.10.SNPs.lst | cut -f1 | grep -vwf finished.genes.chr16together.lst | grep -vwf missing.trees.lst > finished.genes.chr1-15.10SNPs.lst
cat /scratch/ek/solenopsis/2018/SNPIVbusco4all/genes.with.more.than.10.SNPs.lst | cut -f1 | grep -vwf finished.genes.chr1-15.lst | grep -vwf missing.trees.lst | grep -vwf finished.genes.chr16r.lst> finished.genes.chr16nr.10SNPs.lst
cat /scratch/ek/solenopsis/2018/SNPIVbusco4all/genes.with.more.than.10.SNPs.lst | cut -f1 | grep -vwf finished.genes.chr1-15.lst | grep -vwf missing.trees.lst | grep -vwf finished.genes.chr16nr.lst> finished.genes.chr16r.10SNPs.lst

cat finished.genes.chr1-15.10SNPs.lst | parallel "echo {}; cp trees.chr1-15/{}.* trees.chr1-15.10SNPs"
cat finished.genes.chr16r.10SNPs.lst | parallel "echo {}; cp trees.chr16r/{}.* trees.chr16r.10SNPs"
cat finished.genes.chr16nr.10SNPs.lst | parallel "echo {}; cp trees.chr16nr/{}.* trees.chr16nr.10SNPs"

cat trees.chr16nr.10SNPs/*.bestTree.pruned.namefix.tre > astral.chr16nr/input.10SNPs.pruned.tre
cat trees.chr16r.10SNPs/*.bestTree.pruned.namefix.tre > astral.chr16r/input.10SNPs.pruned.tre
cat trees.chr1-15.10SNPs/*.bestTree.pruned.namefix.tre > astral.chr1-15/input.10SNPs.pruned.tre
cat trees.chr16nr.10SNPs/*.bestTree.namefix.tre > astral.chr16nr/input.10SNPs.tre
cat trees.chr16r.10SNPs/*.bestTree.namefix.tre > astral.chr16r/input.10SNPs.tre
cat trees.chr1-15.10SNPs/*.bestTree.namefix.tre > astral.chr1-15/input.10SNPs.tre

################################################################################
### ASTRAL consensus trees
```sh

cd /home/ek/scratch/solenopsis/2018/SNPIVbusco4all/phylo/tmp/
mkdir -p astral.chr16nr
mkdir -p astral.chr16r
mkdir -p astral.chr1-15

cat trees.chr16nr/*.bestTree.namefix.tre > astral.chr16nr/input.tre
cat trees.chr16r/*.bestTree.namefix.tre > astral.chr16r/input.tre
cat trees.chr1-15/*.bestTree.namefix.tre > astral.chr1-15/input.tre

#Chr16nr 97 gene trees
INPUTTREES="astral.chr16nr/input.10SNPs.tre"
OUTPUTFILE="astral.chr16nr/output.10SNPs.tre"
java -Xmx250g \
-D"java.library.path=/usr/local/src/Astral/lib/" \
-jar /usr/local/src/Astral/astral.5.14.3.jar \
--input $INPUTTREES \
--output $OUTPUTFILE \
--cpu-only \
--cpu-threads 50 \
--branch-annotate 3 \
--keep completed

#Chr16r 53 gene trees
INPUTTREES="astral.chr16r/input.10SNPs.tre"
OUTPUTFILE="astral.chr16r/output.10SNPs.tre"
java -Xmx250g \
-D"java.library.path=/usr/local/src/Astral/lib/" \
-jar /usr/local/src/Astral/astral.5.14.3.jar \
--input $INPUTTREES \
--output $OUTPUTFILE \
--cpu-only \
--cpu-threads 50 \
--branch-annotate 3 \
--keep completed

#Chr1-15 1631 gene trees
INPUTTREES="astral.chr1-15/input.10SNPs.tre"
OUTPUTFILE="astral.chr1-15/output.10SNPs.tre"
java -Xmx250g \
-D"java.library.path=/usr/local/src/Astral/lib/" \
-jar /usr/local/src/Astral/astral.5.14.3.jar \
--input $INPUTTREES \
--output $OUTPUTFILE \
--cpu-only \
--cpu-threads 50 \
--branch-annotate 3 \
--keep completed
```
