

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

######################  modeltest per partition (Coding, Noncoding)
F="/scratch/ek/solenopsis/2018/SNPIIIbusco/phylo"
mkdir -p $F/modeltest
BUSCOGENES="ref/gng20170922.fa.busco4.chrs.completeFragmented"

CPUs=4
cat $BUSCOGENES.genelist.tsv | parallel --no-notice -j 15 "echo {} && modeltest-ng -i $F/{}/{}.fasta.allsamples.fa -d nt -p $CPUs -q $F/{}/{}.annotation.partitionfile.NONCODING-CODING.txt -T raxml --rngseed 2 -t fixed-ml-gtr --force --output $F/modeltest/{}.modeltest.NONCODING-CODING"
echo "modeltesting done"


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
   #    5 F81
   #   28 F81+I
   #  127 F81+I+G4
   #   16 GTR
   # 1311 GTR+G4
   # 2681 GTR+I
   #  355 GTR+I+G4
   #   42 HKY
   #  239 HKY+G4
   #  911 HKY+I
   #  111 HKY+I+G4

# coding
   #   62 F81
   #    1 F81+G4
   #   56 F81+I
   #  242 F81+I+G4
   #   10 GTR
   #  197 GTR+G4
   # 1362 GTR+I
   # 2256 GTR+I+G4
   #   22 HKY
   #   84 HKY+G4
   #  931 HKY+I
   #  611 HKY+I+G4

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

# NW_011797094.1	0	1802099	16	469	rec
# NW_011794933.1B	0	4024115	16	470	rec
# NW_011794943.1	0	5603	16	471	rec
# NW_011803760.1	0	1482883	16	472	rec
# NW_011795620.1	0	65617	16	473	nonrec
# NW_011800952.1	0	104328	16	474	nonrec
# NW_011803911.1	0	659464	16	475	nonrec
# NW_011800721.1	0	500137	16	476	nonrec
# NW_011799419.1	0	1115537	16	477	nonrec
# NW_011801243.1	0	1297335	16	478	nonrec
# NW_011794409.1	0	445462	16	479	nonrec
# NW_011794634.1	0	526590	16	480	nonrec
# NW_011795711.1	0	493224	16	481	nonrec
# NW_011797558.1	0	340320	16	482	nonrec
# NW_011794567.1	0	1946445	16	483	nonrec
# NW_011795727.1	0	632960	16	484	nonrec
# NW_011795053.1	0	1052406	16	485	nonrec
# NW_011795861.1	0	4991	16	486	nonrec
# NW_011795747.1	0	6167	16	487	nonrec
# NW_011794844.1	0	1251831	16	488	nonrec
# NW_011794623.1	0	3037	16	489	nonrec
# NW_011801067.1	0	337463	16	490	nonrec
# NW_011796111.1	0	560000	16	491	nonrec
# NW_011796111.1	560000	1594552	16	492	rec
# NW_011797002.1	0	54053	16	493	rec


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

###########################

####### raxml

cd /scratch/ek/solenopsis/2018
F="/scratch/ek/solenopsis/2018/SNPIIIbusco/phylo"
BUSCOGENES="ref/gng20170922.fa.busco4.chrs.completeFragmented"
N=$(wc -l $BUSCOGENES.genelist.tsv | cut -d' ' -f 1) && echo $N     ##5851
OUTDIR="$F/genetrees"
mkdir -p $OUTDIR

cat /scratch/ek/solenopsis/2018/ref/gng20170922.fa.busco4.chr1-15.ComplFragm.lst | parallel --no-notice -j 20 "echo {} && $HOME/scripts/raxml.zsh $F/{}/{}.fasta.allsamples.fa $F/modeltest/{}.modeltest.NONCODING-CODING.part.aic $OUTDIR/{}"

# chr16nr
cat /scratch/ek/solenopsis/2018/ref/gng20170922.fa.busco4.chr16nr.ComplFragm.lst | parallel --no-notice -j 30 "echo {} && $HOME/scripts/raxml.zsh $F/{}/{}.fasta.allsamples.fa $F/modeltest/{}.modeltest.NONCODING-CODING.part.aic $OUTDIR/{}"

# chr16r
cat /scratch/ek/solenopsis/2018/ref/gng20170922.fa.busco4.chr16r.ComplFragm.lst | parallel --no-notice -j 30 "echo {} && $HOME/scripts/raxml.zsh $F/{}/{}.fasta.allsamples.fa $F/modeltest/{}.modeltest.NONCODING-CODING.part.aic $OUTDIR/{}"






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

# 3 genes
# 1767at7399
# 13474at7399
# 24355at7399

cat ../../../ref/gng20170922.fa.busco4.chr16nr.ComplFragm.lst | grep -wf missing.trees.lst > missing.trees.chr16nr.lst

#7 genes
# 8796at7399
# 9291at7399
# 14784at7399
# 15093at7399
# 17602at7399
# 18251at7399
# 22958at7399

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


### ASTRAL consensus trees
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
java -Xmx250g -D"java.library.path=/usr/local/src/Astral/lib/" -jar /usr/local/src/Astral/astral.5.14.3.jar --input $INPUTTREES --output $OUTPUTFILE --cpu-only --cpu-threads 50 --branch-annotate 3 --keep completed

#Chr16r 53 gene trees
INPUTTREES="astral.chr16r/input.10SNPs.tre"
OUTPUTFILE="astral.chr16r/output.10SNPs.tre"
java -Xmx250g -D"java.library.path=/usr/local/src/Astral/lib/" -jar /usr/local/src/Astral/astral.5.14.3.jar --input $INPUTTREES --output $OUTPUTFILE --cpu-only --cpu-threads 50 --branch-annotate 3 --keep completed

#Chr1-15 1631 gene trees
INPUTTREES="astral.chr1-15/input.10SNPs.tre"
OUTPUTFILE="astral.chr1-15/output.10SNPs.tre"
java -Xmx250g -D"java.library.path=/usr/local/src/Astral/lib/" -jar /usr/local/src/Astral/astral.5.14.3.jar --input $INPUTTREES --output $OUTPUTFILE --cpu-only --cpu-threads 50 --branch-annotate 3 --keep completed




### ASTRAL bootstrap support
#CHR16 nr (supergene)

# make output folder
mkdir -p bootstraps.chr16nr.10SNPs

# show the genes
cat finished.genes.chr16nr.10SNPs.lst | tr '\n' '\|' | sed "s/|$//g" | sed "s/|/|\\\./g" | sed "s/^/\\\./g"

# root each bootstrap to geminata and create file with PATHs to each bootstrap tree file (1 per gene)
cat finished.genes.chr16nr.10SNPs.lst | parallel -k -j100 "cat bootstraps/{}.raxml.bootstraps | sed -E -e 's#\.10177at7399|\.10600at7399|\.10654at7399|\.11817at7399|\.11872at7399|\.1235at7399|\.13232at7399|\.13350at7399|\.133at7399|\.1355at7399|\.13603at7399|\.1363at7399|\.13801at7399|\.14031at7399|\.14657at7399|\.1469at7399|\.14848at7399|\.15120at7399|\.15229at7399|\.15542at7399|\.16895at7399|\.17022at7399|\.17425at7399|\.18373at7399|\.1839at7399|\.19334at7399|\.1934at7399|\.1941at7399|\.1948at7399|\.19843at7399|\.19867at7399|\.19978at7399|\.20409at7399|\.20776at7399|\.21000at7399|\.21430at7399|\.21834at7399|\.21849at7399|\.21911at7399|\.2203at7399|\.22521at7399|\.22579at7399|\.2283at7399|\.2286at7399|\.23493at7399|\.24124at7399|\.2463at7399|\.2474at7399|\.2480at7399|\.27445at7399|\.27941at7399|\.30722at7399|\.31769at7399|\.3461at7399|\.3614at7399|\.3616at7399|\.3682at7399|\.4012at7399|\.4426at7399|\.4634at7399|\.4808at7399|\.4819at7399|\.4824at7399|\.4903at7399|\.4960at7399|\.5075at7399|\.5077at7399|\.5090at7399|\.5172at7399|\.5198at7399|\.5281at7399|\.531at7399|\.5524at7399|\.5654at7399|\.5744at7399|\.5998at7399|\.6046at7399|\.6076at7399|\.6105at7399|\.6122at7399|\.6199at7399|\.6471at7399|\.6483at7399|\.7024at7399|\.7129at7399|\.7238at7399|\.7882at7399|\.8126at7399|\.8673at7399|\.867at7399|\.8738at7399|\.8815at7399|\.8985at7399|\.9453at7399|\.9500at7399|\.9598at7399|\.9641at7399##g' | nw_reroot - gem-1-bigB-m-majorityallele > bootstraps.chr16nr.10SNPs/{}.bs && echo {}"
ls -1 bootstraps.chr16nr.10SNPs/* > bootstraps.chr16nr.10SNPs.paths.lst

#un Astral, Chr16nr, 97 gene trees, took ca 8hrs
INPUTTREES="astral.chr16nr/input.10SNPs.tre"
INPUTBOOTSTRAPS="bootstraps.chr16nr.10SNPs.paths.lst"
OUTPUTFILE="astral.chr16nr/output.10SNPs.BS100.nwk"
java -Xmx150g -D"java.library.path=/usr/local/src/Astral/lib/" -jar /usr/local/src/Astral/astral.5.14.3.jar --input $INPUTTREES --output $OUTPUTFILE --cpu-only --cpu-threads 80 --keep completed --bootstraps $INPUTBOOTSTRAPS --seed 12345 --reps 100

#take out the two last trees (greedy consensus and the main tree from Astral for viz)
cat $OUTPUTFILE.nwk | tail -n2 > $OUTPUTFILE.greedyconsensus.main.nwk

### same for Chr1-15 gene trees
mkdir -p bootstraps.chr1-15.10SNPs
cat finished.genes.chr1-15.10SNPs.lst | parallel -k -j100 "cat bootstraps/{}.raxml.bootstraps | sed -E -e 's#\.0at7399|\.10001at7399|\.10019at7399|\.1004at7399|\.10055at7399|\.10058at7399|\.10068at7399|\.10094at7399|\.10099at7399|\.10127at7399|\.10128at7399|\.1012at7399|\.10148at7399|\.10153at7399|\.10170at7399|\.1018at7399|\.1022at7399|\.10233at7399|\.10248at7399|\.10270at7399|\.102at7399|\.10332at7399|\.1035at7399|\.10381at7399|\.10395at7399|\.103at7399|\.10402at7399|\.10403at7399|\.10405at7399|\.10406at7399|\.10466at7399|\.10475at7399|\.10481at7399|\.1048at7399|\.1052at7399|\.10535at7399|\.10555at7399|\.10559at7399|\.10598at7399|\.105at7399|\.10601at7399|\.10610at7399|\.10633at7399|\.10673at7399|\.10679at7399|\.10694at7399|\.10695at7399|\.10696at7399|\.10712at7399|\.10716at7399|\.1072at7399|\.10766at7399|\.1078at7399|\.10823at7399|\.10828at7399|\.10857at7399|\.10871at7399|\.10875at7399|\.1087at7399|\.10891at7399|\.1089at7399|\.10913at7399|\.10922at7399|\.10923at7399|\.10927at7399|\.10944at7399|\.10948at7399|\.10955at7399|\.10982at7399|\.10994at7399|\.11047at7399|\.11048at7399|\.11114at7399|\.11127at7399|\.11138at7399|\.1115at7399|\.1116at7399|\.1117at7399|\.11185at7399|\.11190at7399|\.1119at7399|\.11217at7399|\.1121at7399|\.11226at7399|\.1124at7399|\.11274at7399|\.11279at7399|\.112at7399|\.11316at7399|\.11363at7399|\.11417at7399|\.11422at7399|\.11426at7399|\.11445at7399|\.1144at7399|\.11468at7399|\.11472at7399|\.1148at7399|\.11491at7399|\.11511at7399|\.1154at7399|\.11562at7399|\.11563at7399|\.1158at7399|\.11590at7399|\.11597at7399|\.1159at7399|\.11602at7399|\.11631at7399|\.1163at7399|\.11645at7399|\.11647at7399|\.1164at7399|\.11662at7399|\.11666at7399|\.11672at7399|\.11677at7399|\.1167at7399|\.11694at7399|\.1172at7399|\.11741at7399|\.11752at7399|\.11770at7399|\.11780at7399|\.11854at7399|\.11865at7399|\.11866at7399|\.11891at7399|\.11898at7399|\.118at7399|\.1190at7399|\.11914at7399|\.1192at7399|\.11941at7399|\.11959at7399|\.12023at7399|\.12077at7399|\.12080at7399|\.12158at7399|\.12166at7399|\.12169at7399|\.12171at7399|\.12215at7399|\.12216at7399|\.12227at7399|\.12249at7399|\.12272at7399|\.1227at7399|\.1232at7399|\.12330at7399|\.12331at7399|\.12347at7399|\.123at7399|\.12411at7399|\.12415at7399|\.12424at7399|\.1242at7399|\.12437at7399|\.12443at7399|\.12453at7399|\.12462at7399|\.1246at7399|\.12470at7399|\.12506at7399|\.12507at7399|\.12528at7399|\.12545at7399|\.12588at7399|\.12619at7399|\.1261at7399|\.12638at7399|\.12643at7399|\.12664at7399|\.12668at7399|\.12676at7399|\.12680at7399|\.12690at7399|\.12713at7399|\.12734at7399|\.12753at7399|\.12769at7399|\.12770at7399|\.12779at7399|\.1279at7399|\.1282at7399|\.1286at7399|\.12873at7399|\.12879at7399|\.1287at7399|\.1288at7399|\.12892at7399|\.12894at7399|\.12899at7399|\.12908at7399|\.12917at7399|\.1294at7399|\.12960at7399|\.12965at7399|\.1297at7399|\.1298at7399|\.12997at7399|\.1299at7399|\.13007at7399|\.13021at7399|\.13022at7399|\.13026at7399|\.13035at7399|\.1303at7399|\.13055at7399|\.1307at7399|\.13082at7399|\.13103at7399|\.13111at7399|\.13114at7399|\.13120at7399|\.1313at7399|\.13156at7399|\.13189at7399|\.1318at7399|\.131at7399|\.13218at7399|\.1322at7399|\.13234at7399|\.13244at7399|\.13252at7399|\.1325at7399|\.13263at7399|\.13273at7399|\.13288at7399|\.13303at7399|\.1330at7399|\.1334at7399|\.13358at7399|\.13378at7399|\.13397at7399|\.1342at7399|\.13437at7399|\.1343at7399|\.13449at7399|\.13461at7399|\.13472at7399|\.13486at7399|\.13495at7399|\.13499at7399|\.13550at7399|\.13598at7399|\.13609at7399|\.13617at7399|\.1362at7399|\.13655at7399|\.13658at7399|\.13661at7399|\.13680at7399|\.1369at7399|\.13708at7399|\.13745at7399|\.13746at7399|\.13761at7399|\.13774at7399|\.13793at7399|\.13808at7399|\.13814at7399|\.13818at7399|\.13852at7399|\.13854at7399|\.13894at7399|\.138at7399|\.13919at7399|\.1393at7399|\.13940at7399|\.13951at7399|\.13963at7399|\.1397at7399|\.13985at7399|\.14008at7399|\.14009at7399|\.1403at7399|\.14043at7399|\.14052at7399|\.14075at7399|\.14087at7399|\.14095at7399|\.14135at7399|\.14183at7399|\.14192at7399|\.1419at7399|\.14204at7399|\.14217at7399|\.14261at7399|\.14273at7399|\.1427at7399|\.14287at7399|\.14293at7399|\.14332at7399|\.14337at7399|\.14395at7399|\.14419at7399|\.14453at7399|\.14454at7399|\.1446at7399|\.14483at7399|\.14485at7399|\.1448at7399|\.14527at7399|\.14546at7399|\.14551at7399|\.14557at7399|\.14574at7399|\.14580at7399|\.1460at7399|\.14650at7399|\.14674at7399|\.1467at7399|\.14685at7399|\.14725at7399|\.14746at7399|\.1475at7399|\.14768at7399|\.14775at7399|\.14793at7399|\.14806at7399|\.14821at7399|\.14822at7399|\.1483at7399|\.14878at7399|\.14884at7399|\.148at7399|\.14920at7399|\.14940at7399|\.1494at7399|\.14961at7399|\.14986at7399|\.15017at7399|\.15032at7399|\.1503at7399|\.15062at7399|\.1507at7399|\.1508at7399|\.15099at7399|\.15138at7399|\.15174at7399|\.15204at7399|\.15210at7399|\.15214at7399|\.15222at7399|\.15230at7399|\.15248at7399|\.15250at7399|\.15256at7399|\.15267at7399|\.15274at7399|\.15295at7399|\.1529at7399|\.152at7399|\.1532at7399|\.15338at7399|\.15340at7399|\.15362at7399|\.15377at7399|\.15429at7399|\.1542at7399|\.15432at7399|\.15444at7399|\.15459at7399|\.15467at7399|\.15475at7399|\.1550at7399|\.15514at7399|\.15516at7399|\.15531at7399|\.15532at7399|\.15556at7399|\.15603at7399|\.1565at7399|\.1567at7399|\.15681at7399|\.15687at7399|\.15700at7399|\.15717at7399|\.1572at7399|\.15732at7399|\.15736at7399|\.15752at7399|\.15763at7399|\.15764at7399|\.1576at7399|\.15801at7399|\.15808at7399|\.15829at7399|\.15855at7399|\.15879at7399|\.15890at7399|\.15911at7399|\.15939at7399|\.15950at7399|\.1597at7399|\.15989at7399|\.159at7399|\.15at7399|\.16012at7399|\.1602at7399|\.16094at7399|\.16167at7399|\.1616at7399|\.16173at7399|\.16187at7399|\.16215at7399|\.16226at7399|\.16230at7399|\.16232at7399|\.16233at7399|\.16272at7399|\.16283at7399|\.16285at7399|\.16306at7399|\.16312at7399|\.16323at7399|\.16337at7399|\.16350at7399|\.16364at7399|\.16379at7399|\.16390at7399|\.16392at7399|\.16393at7399|\.16401at7399|\.1640at7399|\.16415at7399|\.16425at7399|\.16433at7399|\.16435at7399|\.16437at7399|\.16476at7399|\.1648at7399|\.16496at7399|\.16514at7399|\.16534at7399|\.16535at7399|\.16549at7399|\.16569at7399|\.1659at7399|\.16620at7399|\.166at7399|\.16708at7399|\.16761at7399|\.16768at7399|\.1678at7399|\.16816at7399|\.16843at7399|\.1686at7399|\.1688at7399|\.16902at7399|\.1694at7399|\.16962at7399|\.1696at7399|\.16998at7399|\.17020at7399|\.17027at7399|\.1702at7399|\.17046at7399|\.17063at7399|\.17099at7399|\.17128at7399|\.1712at7399|\.17139at7399|\.1713at7399|\.17143at7399|\.17160at7399|\.17161at7399|\.17208at7399|\.1720at7399|\.17234at7399|\.1723at7399|\.17267at7399|\.1728at7399|\.17310at7399|\.1734at7399|\.17396at7399|\.17411at7399|\.17467at7399|\.17492at7399|\.17496at7399|\.17498at7399|\.17503at7399|\.17506at7399|\.17525at7399|\.17530at7399|\.17534at7399|\.17575at7399|\.17581at7399|\.17587at7399|\.175at7399|\.17609at7399|\.17635at7399|\.17638at7399|\.1763at7399|\.17661at7399|\.17677at7399|\.17681at7399|\.17697at7399|\.17732at7399|\.17753at7399|\.17768at7399|\.17775at7399|\.17801at7399|\.17817at7399|\.1783at7399|\.17875at7399|\.17890at7399|\.1789at7399|\.178at7399|\.17946at7399|\.1794at7399|\.17965at7399|\.17987at7399|\.1799at7399|\.18007at7399|\.1802at7399|\.18039at7399|\.1803at7399|\.18059at7399|\.18073at7399|\.1810at7399|\.1811at7399|\.18151at7399|\.18166at7399|\.18175at7399|\.18177at7399|\.18212at7399|\.18214at7399|\.18215at7399|\.18222at7399|\.18247at7399|\.1824at7399|\.1825at7399|\.18275at7399|\.1827at7399|\.1828at7399|\.1831at7399|\.18324at7399|\.18335at7399|\.18339at7399|\.18368at7399|\.1836at7399|\.18403at7399|\.18411at7399|\.1841at7399|\.18427at7399|\.18447at7399|\.1846at7399|\.18482at7399|\.18557at7399|\.18558at7399|\.18572at7399|\.18576at7399|\.18585at7399|\.18597at7399|\.18647at7399|\.1867at7399|\.18710at7399|\.18729at7399|\.18733at7399|\.18773at7399|\.18811at7399|\.18833at7399|\.18867at7399|\.18887at7399|\.18898at7399|\.18899at7399|\.1889at7399|\.18917at7399|\.1891at7399|\.18948at7399|\.1895at7399|\.18988at7399|\.1898at7399|\.18993at7399|\.1899at7399|\.19001at7399|\.1903at7399|\.19068at7399|\.19071at7399|\.19090at7399|\.1911at7399|\.19175at7399|\.19224at7399|\.19226at7399|\.1923at7399|\.19240at7399|\.1926at7399|\.19279at7399|\.19283at7399|\.19316at7399|\.19342at7399|\.19350at7399|\.19356at7399|\.1939at7399|\.19450at7399|\.19457at7399|\.19483at7399|\.1949at7399|\.19518at7399|\.19567at7399|\.19583at7399|\.19588at7399|\.19594at7399|\.19600at7399|\.19624at7399|\.19636at7399|\.19657at7399|\.19658at7399|\.19666at7399|\.1966at7399|\.1969at7399|\.19748at7399|\.19753at7399|\.19810at7399|\.19861at7399|\.19897at7399|\.19954at7399|\.19992at7399|\.19999at7399|\.1999at7399|\.1at7399|\.20007at7399|\.20048at7399|\.2007at7399|\.20128at7399|\.20138at7399|\.2016at7399|\.20311at7399|\.20392at7399|\.20396at7399|\.2040at7399|\.20419at7399|\.20465at7399|\.20475at7399|\.20492at7399|\.20494at7399|\.20502at7399|\.2055at7399|\.20560at7399|\.20579at7399|\.205at7399|\.20607at7399|\.20633at7399|\.2064at7399|\.20656at7399|\.20685at7399|\.20693at7399|\.20735at7399|\.20745at7399|\.2079at7399|\.2083at7399|\.20854at7399|\.20864at7399|\.20892at7399|\.20919at7399|\.20925at7399|\.20926at7399|\.2097at7399|\.20980at7399|\.2098at7399|\.20992at7399|\.21022at7399|\.2102at7399|\.21119at7399|\.21137at7399|\.2113at7399|\.21151at7399|\.21155at7399|\.21171at7399|\.21173at7399|\.21181at7399|\.21184at7399|\.21216at7399|\.21232at7399|\.21283at7399|\.21296at7399|\.21366at7399|\.21380at7399|\.21415at7399|\.2142at7399|\.21448at7399|\.21473at7399|\.21495at7399|\.21530at7399|\.21533at7399|\.21586at7399|\.21617at7399|\.21629at7399|\.2162at7399|\.21641at7399|\.2164at7399|\.21651at7399|\.2171at7399|\.21740at7399|\.21744at7399|\.21747at7399|\.2176at7399|\.21825at7399|\.21883at7399|\.21892at7399|\.21907at7399|\.21925at7399|\.21926at7399|\.219at7399|\.22040at7399|\.22103at7399|\.22109at7399|\.22114at7399|\.22185at7399|\.22265at7399|\.222at7399|\.22305at7399|\.22311at7399|\.2231at7399|\.2233at7399|\.2236at7399|\.223at7399|\.2248at7399|\.22499at7399|\.2252at7399|\.22541at7399|\.22542at7399|\.2256at7399|\.2257at7399|\.22587at7399|\.2262at7399|\.2264at7399|\.2267at7399|\.22694at7399|\.226at7399|\.22727at7399|\.2274at7399|\.22803at7399|\.22806at7399|\.22820at7399|\.22829at7399|\.22915at7399|\.2293at7399|\.22942at7399|\.2294at7399|\.22951at7399|\.22959at7399|\.2298at7399|\.2299at7399|\.2306at7399|\.23074at7399|\.23077at7399|\.23124at7399|\.23236at7399|\.2324at7399|\.23254at7399|\.23276at7399|\.2328at7399|\.232at7399|\.2336at7399|\.23382at7399|\.233at7399|\.23404at7399|\.23410at7399|\.23435at7399|\.23456at7399|\.2345at7399|\.2346at7399|\.23483at7399|\.2353at7399|\.23586at7399|\.23587at7399|\.23601at7399|\.23618at7399|\.2363at7399|\.23646at7399|\.2367at7399|\.23681at7399|\.23691at7399|\.23705at7399|\.23758at7399|\.23843at7399|\.23857at7399|\.23901at7399|\.23905at7399|\.23918at7399|\.23928at7399|\.23934at7399|\.2393at7399|\.23949at7399|\.2397at7399|\.23985at7399|\.24115at7399|\.24126at7399|\.24155at7399|\.24172at7399|\.24193at7399|\.2419at7399|\.241at7399|\.24204at7399|\.24232at7399|\.2428at7399|\.24300at7399|\.2430at7399|\.24348at7399|\.24356at7399|\.2436at7399|\.24378at7399|\.24398at7399|\.2441at7399|\.2444at7399|\.24462at7399|\.24465at7399|\.244at7399|\.24537at7399|\.24563at7399|\.24589at7399|\.24615at7399|\.2461at7399|\.24624at7399|\.2468at7399|\.2476at7399|\.24788at7399|\.24795at7399|\.24836at7399|\.2484at7399|\.24906at7399|\.2499at7399|\.249at7399|\.25004at7399|\.25005at7399|\.2500at7399|\.25012at7399|\.25041at7399|\.2504at7399|\.25052at7399|\.25112at7399|\.25128at7399|\.25189at7399|\.251at7399|\.25200at7399|\.25222at7399|\.25243at7399|\.25282at7399|\.25286at7399|\.2530at7399|\.2536at7399|\.2543at7399|\.2548at7399|\.25514at7399|\.25516at7399|\.25543at7399|\.2554at7399|\.2556at7399|\.255at7399|\.2560at7399|\.25659at7399|\.2569at7399|\.25708at7399|\.25719at7399|\.2571at7399|\.25729at7399|\.2572at7399|\.25730at7399|\.25747at7399|\.25788at7399|\.25807at7399|\.2582at7399|\.2583at7399|\.2588at7399|\.25939at7399|\.25969at7399|\.25987at7399|\.26041at7399|\.2605at7399|\.2609at7399|\.26116at7399|\.26140at7399|\.26142at7399|\.26183at7399|\.26212at7399|\.2626at7399|\.26289at7399|\.2630at7399|\.2638at7399|\.2641at7399|\.2646at7399|\.2660at7399|\.26639at7399|\.26669at7399|\.26725at7399|\.2674at7399|\.26757at7399|\.2676at7399|\.2678at7399|\.267at7399|\.2682at7399|\.26848at7399|\.26885at7399|\.26903at7399|\.2694at7399|\.26951at7399|\.26966at7399|\.27006at7399|\.2702at7399|\.27030at7399|\.27034at7399|\.27039at7399|\.27069at7399|\.2707at7399|\.27119at7399|\.2713at7399|\.2717at7399|\.2720at7399|\.27290at7399|\.27352at7399|\.27380at7399|\.2741at7399|\.27440at7399|\.27482at7399|\.2748at7399|\.27572at7399|\.27587at7399|\.27595at7399|\.27677at7399|\.27713at7399|\.2778at7399|\.2780at7399|\.27813at7399|\.27814at7399|\.27833at7399|\.27852at7399|\.27869at7399|\.2790at7399|\.27954at7399|\.2796at7399|\.28004at7399|\.28077at7399|\.2809at7399|\.28150at7399|\.281at7399|\.28208at7399|\.2821at7399|\.28247at7399|\.282at7399|\.2836at7399|\.2848at7399|\.28495at7399|\.28517at7399|\.28565at7399|\.285at7399|\.2870at7399|\.28744at7399|\.2874at7399|\.28805at7399|\.28820at7399|\.2885at7399|\.2892at7399|\.28932at7399|\.28949at7399|\.29010at7399|\.2908at7399|\.29095at7399|\.290at7399|\.2922at7399|\.2927at7399|\.292at7399|\.29397at7399|\.29435at7399|\.2943at7399|\.29553at7399|\.29579at7399|\.29717at7399|\.29727at7399|\.2976at7399|\.297at7399|\.29817at7399|\.2983at7399|\.29906at7399|\.29925at7399|\.2996at7399|\.2at7399|\.30075at7399|\.30088at7399|\.3010at7399|\.30196at7399|\.30268at7399|\.3032at7399|\.3034at7399|\.3041at7399|\.3043at7399|\.30463at7399|\.3053at7399|\.30699at7399|\.306at7399|\.3089at7399|\.30955at7399|\.31022at7399|\.31067at7399|\.31184at7399|\.311at7399|\.3122at7399|\.3130at7399|\.3132at7399|\.3134at7399|\.3144at7399|\.31526at7399|\.31553at7399|\.3158at7399|\.3163at7399|\.3167at7399|\.31749at7399|\.3174at7399|\.317at7399|\.3185at7399|\.3186at7399|\.3187at7399|\.321at7399|\.32312at7399|\.324at7399|\.3274at7399|\.3286at7399|\.32at7399|\.3311at7399|\.3332at7399|\.3337at7399|\.3338at7399|\.3345at7399|\.3349at7399|\.3354at7399|\.3356at7399|\.335at7399|\.3365at7399|\.3375at7399|\.3384at7399|\.3385at7399|\.3402at7399|\.3411at7399|\.3414at7399|\.341at7399|\.3423at7399|\.343at7399|\.3441at7399|\.3455at7399|\.3467at7399|\.3468at7399|\.346at7399|\.3473at7399|\.3504at7399|\.3509at7399|\.3516at7399|\.3522at7399|\.3523at7399|\.352at7399|\.3538at7399|\.3564at7399|\.3578at7399|\.359at7399|\.3606at7399|\.3612at7399|\.3630at7399|\.3635at7399|\.3641at7399|\.3645at7399|\.3648at7399|\.3651at7399|\.3654at7399|\.3660at7399|\.3668at7399|\.3669at7399|\.3670at7399|\.3678at7399|\.3679at7399|\.3692at7399|\.3707at7399|\.3713at7399|\.3721at7399|\.3733at7399|\.3736at7399|\.3739at7399|\.3754at7399|\.3763at7399|\.3768at7399|\.376at7399|\.3776at7399|\.3777at7399|\.3780at7399|\.3804at7399|\.3817at7399|\.3838at7399|\.3860at7399|\.3867at7399|\.3873at7399|\.3883at7399|\.3887at7399|\.3899at7399|\.3902at7399|\.3904at7399|\.3906at7399|\.3914at7399|\.3921at7399|\.3928at7399|\.393at7399|\.3945at7399|\.3952at7399|\.3957at7399|\.395at7399|\.3966at7399|\.3973at7399|\.3978at7399|\.3995at7399|\.4002at7399|\.4007at7399|\.4008at7399|\.400at7399|\.4011at7399|\.4022at7399|\.4026at7399|\.4039at7399|\.404at7399|\.4052at7399|\.4069at7399|\.4075at7399|\.4080at7399|\.4087at7399|\.4091at7399|\.4094at7399|\.4098at7399|\.4104at7399|\.4105at7399|\.4117at7399|\.412at7399|\.4135at7399|\.415at7399|\.416at7399|\.4171at7399|\.4176at7399|\.4198at7399|\.4199at7399|\.420at7399|\.4223at7399|\.4228at7399|\.4233at7399|\.4258at7399|\.425at7399|\.4278at7399|\.4284at7399|\.4291at7399|\.4296at7399|\.4297at7399|\.4298at7399|\.4302at7399|\.430at7399|\.4318at7399|\.4323at7399|\.4324at7399|\.4329at7399|\.4340at7399|\.4341at7399|\.4347at7399|\.435at7399|\.436at7399|\.4382at7399|\.4383at7399|\.4388at7399|\.438at7399|\.4391at7399|\.4393at7399|\.4417at7399|\.4441at7399|\.4448at7399|\.4450at7399|\.4473at7399|\.4487at7399|\.4490at7399|\.4516at7399|\.4533at7399|\.453at7399|\.4542at7399|\.4555at7399|\.4557at7399|\.4563at7399|\.4565at7399|\.4568at7399|\.4581at7399|\.4582at7399|\.4591at7399|\.4606at7399|\.4608at7399|\.4621at7399|\.462at7399|\.4643at7399|\.4646at7399|\.4659at7399|\.4678at7399|\.4689at7399|\.468at7399|\.4714at7399|\.4743at7399|\.4758at7399|\.4761at7399|\.4768at7399|\.476at7399|\.4803at7399|\.4811at7399|\.4865at7399|\.4867at7399|\.4873at7399|\.4875at7399|\.4890at7399|\.4921at7399|\.4949at7399|\.4955at7399|\.4962at7399|\.4975at7399|\.4977at7399|\.497at7399|\.4988at7399|\.4997at7399|\.49at7399|\.5000at7399|\.501at7399|\.5028at7399|\.5048at7399|\.504at7399|\.5059at7399|\.5086at7399|\.5122at7399|\.5123at7399|\.5146at7399|\.5161at7399|\.5165at7399|\.5175at7399|\.5190at7399|\.5191at7399|\.5203at7399|\.5214at7399|\.5221at7399|\.523at7399|\.5241at7399|\.5243at7399|\.527at7399|\.5290at7399|\.5291at7399|\.5299at7399|\.5307at7399|\.5317at7399|\.5318at7399|\.5323at7399|\.5337at7399|\.534at7399|\.5353at7399|\.5363at7399|\.5406at7399|\.5410at7399|\.5414at7399|\.5421at7399|\.5431at7399|\.5434at7399|\.5469at7399|\.5498at7399|\.5500at7399|\.551at7399|\.5520at7399|\.552at7399|\.5531at7399|\.5536at7399|\.5558at7399|\.5566at7399|\.556at7399|\.5577at7399|\.558at7399|\.5600at7399|\.5608at7399|\.560at7399|\.5618at7399|\.5621at7399|\.5626at7399|\.5632at7399|\.5636at7399|\.5639at7399|\.5656at7399|\.5679at7399|\.569at7399|\.5725at7399|\.572at7399|\.5730at7399|\.5754at7399|\.5766at7399|\.5789at7399|\.5800at7399|\.5827at7399|\.5832at7399|\.5863at7399|\.5880at7399|\.5891at7399|\.590at7399|\.5949at7399|\.5963at7399|\.5981at7399|\.5988at7399|\.5990at7399|\.6038at7399|\.6063at7399|\.6085at7399|\.6088at7399|\.60at7399|\.6101at7399|\.6104at7399|\.6115at7399|\.6156at7399|\.6159at7399|\.6174at7399|\.61at7399|\.6200at7399|\.6202at7399|\.6204at7399|\.6239at7399|\.6254at7399|\.6258at7399|\.6264at7399|\.6271at7399|\.6273at7399|\.6274at7399|\.6305at7399|\.6318at7399|\.6322at7399|\.6326at7399|\.6331at7399|\.6334at7399|\.6351at7399|\.6352at7399|\.6358at7399|\.635at7399|\.6364at7399|\.6384at7399|\.6400at7399|\.6410at7399|\.6413at7399|\.6420at7399|\.6427at7399|\.6437at7399|\.6442at7399|\.6444at7399|\.6477at7399|\.6488at7399|\.6496at7399|\.6508at7399|\.6510at7399|\.6522at7399|\.6549at7399|\.6560at7399|\.6564at7399|\.6571at7399|\.6595at7399|\.6604at7399|\.6605at7399|\.6606at7399|\.6616at7399|\.6617at7399|\.6641at7399|\.6644at7399|\.6645at7399|\.6647at7399|\.6652at7399|\.6659at7399|\.665at7399|\.6680at7399|\.6704at7399|\.670at7399|\.6731at7399|\.6735at7399|\.674at7399|\.6758at7399|\.6768at7399|\.6814at7399|\.6821at7399|\.6829at7399|\.683at7399|\.684at7399|\.685at7399|\.6892at7399|\.6901at7399|\.6910at7399|\.6918at7399|\.691at7399|\.6940at7399|\.6949at7399|\.6955at7399|\.6963at7399|\.6966at7399|\.6997at7399|\.700at7399|\.7017at7399|\.7022at7399|\.7030at7399|\.7033at7399|\.7034at7399|\.7049at7399|\.7057at7399|\.7083at7399|\.7089at7399|\.7104at7399|\.7108at7399|\.711at7399|\.7148at7399|\.7178at7399|\.7184at7399|\.7212at7399|\.7214at7399|\.7216at7399|\.7218at7399|\.7225at7399|\.7229at7399|\.7237at7399|\.7256at7399|\.7263at7399|\.7282at7399|\.7284at7399|\.728at7399|\.7293at7399|\.7296at7399|\.7318at7399|\.7319at7399|\.7324at7399|\.734at7399|\.7356at7399|\.7382at7399|\.7390at7399|\.739at7399|\.7432at7399|\.7444at7399|\.7474at7399|\.7476at7399|\.7561at7399|\.7572at7399|\.7587at7399|\.7593at7399|\.759at7399|\.7600at7399|\.760at7399|\.7613at7399|\.7617at7399|\.7640at7399|\.764at7399|\.7675at7399|\.767at7399|\.7707at7399|\.7728at7399|\.7735at7399|\.7743at7399|\.7750at7399|\.7751at7399|\.7755at7399|\.775at7399|\.7774at7399|\.7778at7399|\.7779at7399|\.7783at7399|\.7784at7399|\.7786at7399|\.7789at7399|\.778at7399|\.7812at7399|\.7815at7399|\.7819at7399|\.7830at7399|\.7835at7399|\.7839at7399|\.7844at7399|\.7859at7399|\.787at7399|\.7887at7399|\.788at7399|\.7901at7399|\.7910at7399|\.7926at7399|\.7950at7399|\.7966at7399|\.7980at7399|\.7998at7399|\.8029at7399|\.8037at7399|\.8055at7399|\.8056at7399|\.8059at7399|\.8060at7399|\.8079at7399|\.8084at7399|\.8097at7399|\.8130at7399|\.8135at7399|\.8152at7399|\.8153at7399|\.8157at7399|\.8168at7399|\.8215at7399|\.8217at7399|\.8231at7399|\.8247at7399|\.8263at7399|\.8268at7399|\.826at7399|\.8282at7399|\.82at7399|\.8304at7399|\.830at7399|\.8326at7399|\.832at7399|\.8335at7399|\.833at7399|\.8344at7399|\.8350at7399|\.8364at7399|\.8373at7399|\.8473at7399|\.8479at7399|\.8482at7399|\.848at7399|\.8504at7399|\.850at7399|\.8517at7399|\.8529at7399|\.8551at7399|\.8555at7399|\.855at7399|\.8561at7399|\.8572at7399|\.8578at7399|\.8586at7399|\.8602at7399|\.8615at7399|\.8631at7399|\.8632at7399|\.863at7399|\.8657at7399|\.8659at7399|\.8669at7399|\.8672at7399|\.8682at7399|\.8706at7399|\.8720at7399|\.872at7399|\.8731at7399|\.8734at7399|\.8744at7399|\.874at7399|\.8773at7399|\.877at7399|\.8780at7399|\.8792at7399|\.8813at7399|\.882at7399|\.8851at7399|\.8867at7399|\.8895at7399|\.8897at7399|\.88at7399|\.892at7399|\.8966at7399|\.9001at7399|\.9032at7399|\.9051at7399|\.908at7399|\.9091at7399|\.9106at7399|\.9134at7399|\.9182at7399|\.9188at7399|\.9194at7399|\.9215at7399|\.921at7399|\.922at7399|\.9235at7399|\.923at7399|\.9247at7399|\.926at7399|\.9279at7399|\.92at7399|\.9300at7399|\.9301at7399|\.9307at7399|\.9309at7399|\.9336at7399|\.934at7399|\.9404at7399|\.9443at7399|\.9479at7399|\.9483at7399|\.948at7399|\.9498at7399|\.9499at7399|\.9511at7399|\.9514at7399|\.9534at7399|\.9545at7399|\.9562at7399|\.9570at7399|\.9585at7399|\.9586at7399|\.958at7399|\.9595at7399|\.9608at7399|\.9613at7399|\.9618at7399|\.9620at7399|\.964at7399|\.9675at7399|\.9677at7399|\.968at7399|\.971at7399|\.9723at7399|\.9736at7399|\.9761at7399|\.9763at7399|\.9780at7399|\.9785at7399|\.9799at7399|\.9813at7399|\.9825at7399|\.9830at7399|\.9855at7399|\.987at7399|\.9888at7399|\.988at7399|\.9904at7399|\.9917at7399|\.9921at7399|\.9932at7399|\.9944at7399|\.9950at7399|\.9980at7399|\.9984at7399|\.9997at7399##g'  | nw_reroot - gem-1-bigB-m-majorityallele > bootstraps.chr1-15.10SNPs/{}.bs"
ls -1 bootstraps.chr1-15.10SNPs/* > bootstraps.chr1-15.10SNPs.paths.lst

cat astral.chr1-15/input.10SNPs.tre | nw_reroot - gem-1-bigB-m-majorityallele > astral.chr1-15/input.10SNPs.tre

#Chr1-15, 1631 gene trees, ca 148 hrs computation time
INPUTTREES="astral.chr1-15/input.10SNPs.AR142-AR66-pruned.tre"
INPUTBOOTSTRAPS="bootstraps.chr1-15.10SNPs.paths.lst"
OUTPUTFILE="astral.chr1-15/output.10SNPs.BS100.nwk"
java -Xmx150g -D"java.library.path=/usr/local/src/Astral/lib/" -jar /usr/local/src/Astral/astral.5.14.3.jar --input $INPUTTREES --output $OUTPUTFILE --cpu-only --cpu-threads 70 --keep completed --bootstraps $INPUTBOOTSTRAPS --seed 12345 --reps 100 

cat $OUTPUTFILE.nwk | tail -n2 > $OUTPUTFILE.greedyconsensus.main.nwk
