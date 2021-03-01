# ---------------------------------------------------------------------------- #
# Define single-copy gene regions using BUSCO

conda activate busco4

LINEAGE="/home/ek/progz/hymenoptera_odb10"  	## modify
LINEAGENAME=$(echo $LINEAGE | rev | cut -d"/" -f1 | rev)
AUGUSTUSSPECIES="bombus_impatiens1"  				  			## modify
echo "this run is uses this Database: $LINEAGE"
echo "this run is used this starting species instead of Dmel: $AUGUSTUSSPECIES"
echo "this run is used this Database: $LINEAGE" > $INPUTFASTA.stats
echo "this run is used this starting species instead of Dmel: $AUGUSTUSSPECIES" >> $INPUTFASTA.stats

# run BUSCO
busco -i $INPUTFOLDER/$INPUTFILE --force --mode genome --lineage_dataset $LINEAGE --offline --augustus_species $AUGUSTUSSPECIES --cpu $CPUs --evalue 1e-03 --limit 5 --out $INPUTFILE.busco4

cat $INPUTFOLDER/$INPUTFILE.busco4/short_summary.specific.$LINEAGENAME.$INPUTFILE.busco4.txt >> $INPUTFASTA.stats

cat $INPUTFOLDER/$INPUTFILE.busco4/'run_'$LINEAGENAME/full_table.tsv > $INPUTFASTA.busco4.full.table.tsv
cat $INPUTFOLDER/$INPUTFILE.busco4/'run_'$LINEAGENAME/full_table.tsv > $INPUTFASTA.busco4.missing.tsv

# extract the BUSCO gene sequences for later use (e.g. in phylogenies)
for GENE in $(ls -1 $INPUTFOLDER/$INPUTFILE.busco4/'run_'$LINEAGENAME/busco_sequences/fragmented_busco_sequences/*.faa | rev | cut -d"/" -f1 | cut -d"." -f2- | rev)
do
   cat $INPUTFOLDER/$INPUTFILE.busco4/'run_'$LINEAGENAME/busco_sequences/fragmented_busco_sequences/$GENE.faa | seqtk seq -l 0 -C | sed -r "s,^>,>`echo -n $GENE`_`echo -n $INPUTFILENAME`_,g" >> $INPUTFASTA.busco4.fragmented.faa
      cat $INPUTFOLDER/$INPUTFILE.busco4/'run_'$LINEAGENAME/busco_sequences/fragmented_busco_sequences/$GENE.fna | seqtk seq -l 0 -C | sed -r "s,^>,>`echo -n $GENE`_`echo -n $INPUTFILENAME`_,g" >> $INPUTFASTA.busco4.fragmented.fna
done

for GENE in $(ls -1 $INPUTFOLDER/$INPUTFILE.busco4/'run_'$LINEAGENAME/busco_sequences/multi_copy_busco_sequences/*.faa | rev | cut -d"/" -f1 | cut -d"." -f2- | rev)
do
   cat $INPUTFOLDER/$INPUTFILE.busco4/'run_'$LINEAGENAME/busco_sequences/multi_copy_busco_sequences/$GENE.faa | seqtk seq -l 0 -C | sed -r "s,^>,>`echo -n $GENE`_`echo -n $INPUTFILENAME`_,g" >> $INPUTFASTA.busco4.duplicated.faa
      cat $INPUTFOLDER/$INPUTFILE.busco4/'run_'$LINEAGENAME/busco_sequences/multi_copy_busco_sequences/$GENE.fna | seqtk seq -l 0 -C | sed -r "s,^>,>`echo -n $GENE`_`echo -n $INPUTFILENAME`_,g" >> $INPUTFASTA.busco4.duplicated.fna
done

for GENE in $(ls -1 $INPUTFOLDER/$INPUTFILE.busco4/'run_'$LINEAGENAME/busco_sequences/single_copy_busco_sequences/*.faa | rev | cut -d"/" -f1 | cut -d"." -f2- | rev)
do
   cat $INPUTFOLDER/$INPUTFILE.busco4/'run_'$LINEAGENAME/busco_sequences/single_copy_busco_sequences/$GENE.faa | seqtk seq -l 0 -C | sed -r "s,^>,>`echo -n $GENE`_`echo -n $INPUTFILENAME`_,g" >> $INPUTFASTA.busco4.singlecopy.faa
      cat $INPUTFOLDER/$INPUTFILE.busco4/'run_'$LINEAGENAME/busco_sequences/single_copy_busco_sequences/$GENE.fna | seqtk seq -l 0 -C | sed -r "s,^>,>`echo -n $GENE`_`echo -n $INPUTFILENAME`_,g" >> $INPUTFASTA.busco4.singlecopy.fna
done


# could also use pigz for faster compression of large files: /share/pool/ek/busco4.databases/pigz -9 -p $CPUs
gzip -9 $INPUTFASTA.busco4.duplicated.faa
gzip -9 $INPUTFASTA.busco4.fragmented.faa
gzip -9 $INPUTFASTA.busco4.singlecopy.faa
gzip -9 $INPUTFASTA.busco4.duplicated.fna
gzip -9 $INPUTFASTA.busco4.fragmented.fna
gzip -9 $INPUTFASTA.busco4.singlecopy.fna

echo "Busco v4 done:"
cat $INPUTFOLDER/$INPUTFILE.busco4/'run_'$LINEAGENAME/short_summary.txt

# compress the BUSCO results to save space
tar -czf $INPUTFILE.busco4.tar.gz $INPUTFILE.busco4

cat gng20170922.busco4.full.table.tsv | grep -P 'Complete|Fragmented' | tabtk cut -r -f3,4,5,1,2,6,7,8 > gng20170922.busco4.complete.tsv
cat gng20170922.busco4.complete.tsv | cut -f4 > gng20170922.busco4.complete.genelist.tsv

cat gng20170922.busco4.complete.tsv | awk 'BEGIN{FS="\t"}{OFS="\t"}{print $1,$2 - 1000,$3 + 1000,$4,$5,$6,$7,$8}' | awk 'BEGIN{FS="\t"}{OFS="\t"}{print $1,$2<0?1:$2,$3,$4,$5,$6,$7,$8}' > gng20170922.busco4.complete.extended1000bp.tsv

## check coordinates > than scaffolds
INPUTFILE="/scratch/ek/solenopsis/2018/ref/gng20170922.busco4/gng20170922.busco4.complete.extended1000bp"
N=$(cat $INPUTFILE.tsv | wc -l) && echo $N
for (( i = 1 ; i < $N+1 ; i++))
 do
#echo $i
SCF=$(cat $INPUTFILE.tsv | sed -n $i'p' | cut -f 1)
GENEEND=$(cat $INPUTFILE.tsv | sed -n $i'p' | cut -f 3)
SCFEND=$(cat /scratch/ek/solenopsis/2018/ref/gng20170922.fa.fai | grep $SCF | cut -f2)
echo $i"\t"$SCF"\t"$GENEEND"\t"$SCFEND
cat $INPUTFILE.tsv | sed -n $i'p' | awk -v SCFEND=$SCFEND 'BEGIN{FS="\t"}{OFS="\t"} { if (SCFEND<$3) {print $1,$2,SCFEND,$4,$5,$6,$7,$8} else {print $0} }' >> $INPUTFILE.coordinates-fixed.tsv
done
#######

###############    make region files for samtools/VCF
cat $INPUTFILE.coordinates-fixed.tsv | cut -f1-3 | sort -k1,1 -k2,2g | sed "s,\.1\t,.1:,g" | sed "s,\.1A\t,.1A:,g" | sed "s,\.1B\t,.1B:,g" | sed "s,\.1C\t,.1C:,g" | sed "s,\t,-,g" > $INPUTFILE.coordinates-fixed.region

############### bed merge regions to avoid overlaps, merge also if <= 1kb distance
cat $INPUTFILE.coordinates-fixed.region | sed "s,-,\t,g" | sed "s,:,\t,g" > $INPUTFILE.coordinates-fixed.bed
bedtools merge -d 1000 -i $INPUTFILE.coordinates-fixed.bed > $INPUTFILE.coordinates-fixed.merged.bed

cat $INPUTFILE.coordinates-fixed.merged.bed | sed "s,\.1\t,.1:,g" | sed "s,\.1A\t,.1A:,g" | sed "s,\.1B\t,.1B:,g" | sed "s,\.1C\t,.1C:,g" | sed "s,\t,-,g" > $INPUTFILE.coordinates-fixed.merged.region


# ---------------------------------------------------------------------------- #
# Define single-copy gene regions using BUSCO

# mosdepth
BAMOUT="/scratch/solenopsis.bams"

F="$BAMOUT"
cd $F

CPUs=4
OPATH="$F/mosdepth"
mkdir -p $OPATH

REF="/scratch/ek/solenopsis/2018/ref/gng20170922.fa"
BUSCOTABLE="$REF.busco4.full.table.tsv"


cat $BUSCOTABLE | grep -E "^#" -v | awk '{ if ($2 == "Fragmented" || $2 == "Complete") print }' | cut -f 3,4,5 > $OPATH/busco4.genes.bed


LIST=$(ls -1 bam/*.bam | rev | cut -d "/" -f 1 | cut -d "." -f2 | rev)
echo $LIST | tr '\n' '|' | sed "s,|$,,g"


CPUs=2
N=$(echo $LIST | wc -l) && echo $N
for (( i = 1 ; i < $N+1 ; i++))
 do
echo "$i of $N"
SAMPLE=$(echo $LIST | sed -n $i'p') && echo $SAMPLE
ll bam/$SAMPLE.bam
mosdepth -t $CPUs --use-median --by $OPATH/busco4.genes.bed $OPATH/$SAMPLE.mosdepth.busco4 $F/bam/$SAMPLE.bam
mosdepth -t $CPUs --by 10 $OPATH/$SAMPLE.mosdepth.10bp $F/bam/$SAMPLE.bam
#
VAR=$(zcat $OPATH/$SAMPLE.mosdepth.busco4.regions.bed.gz | awk '{a[NR]=$4}END{asort(a,b);median=NR%2?b[(NR+1)/2]:(b[NR/2]+b[NR/2+1])/2;for(i=0;i++<NR;)sum1+=b[i];average=sum1/NR;for(i=0;i++<NR;)sum2+=(b[i]-average)**2;SD=sqrt(sum2/NR);printf "\t%f\t%f\t",median,SD}')
MEDIAN=$(echo $VAR | cut -f 2)
STANDARD_DERIVATION=$(echo $VAR | cut -f 3)
echo $SAMPLE"\t"$MEDIAN"\t"$STANDARD_DERIVATION
echo -e $SAMPLE"\t"$MEDIAN"\t"$STANDARD_DERIVATION >> $OPATH/mosdepth.busco4.stats
#
zcat $OPATH/$SAMPLE.mosdepth.10bp.per-base.bed.gz | tee \
>(awk -v med=$MEDIAN -v std=$STANDARD_DERIVATION '{ if (med-3*std<=$4 && $4<=med+3*std) print }' > $OPATH/$SAMPLE.mosdepth.10bp.normal_coverage.bed) \
>(awk -v med=$MEDIAN -v std=$STANDARD_DERIVATION '{ if ($4>med+3*std) print }' > $OPATH/$SAMPLE.mosdepth.10bp.high_coverage.bed) \
>(awk -v med=$MEDIAN -v std=$STANDARD_DERIVATION '{ if ($4<3) print }' > $OPATH/$SAMPLE.mosdepth.10bp.low_coverage.bed) | \
awk -v med=$MEDIAN -v std=$STANDARD_DERIVATION '{ if ($4>2) print }' > $OPATH/$SAMPLE.mosdepth.10bp.callable.bed
#
rm -f $OPATH/$SAMPLE.mosdepth.10bp.per-base.bed
#
pigz -9 -p 4 $OPATH/$SAMPLE.mosdepth.10bp.mosdepth.global.dist.txt
pigz -9 -p 4 $OPATH/$SAMPLE.mosdepth.10bp.mosdepth.region.dist.txt
pigz -9 -p 4 $OPATH/$SAMPLE.mosdepth.10bp.mosdepth.summary.txt
#
pigz -9 -p 4 $OPATH/$SAMPLE.mosdepth.busco4.mosdepth.global.dist.txt
pigz -9 -p 4 $OPATH/$SAMPLE.mosdepth.busco4.mosdepth.region.dist.txt
pigz -9 -p 4 $OPATH/$SAMPLE.mosdepth.busco4.mosdepth.summary.txt
#
pigz -9 -p 4 $OPATH/$SAMPLE.mosdepth.10bp.normal_coverage.bed
pigz -9 -p 4 $OPATH/$SAMPLE.mosdepth.10bp.high_coverage.bed
pigz -9 -p 4 $OPATH/$SAMPLE.mosdepth.10bp.low_coverage.bed
pigz -9 -p 4 $OPATH/$SAMPLE.mosdepth.10bp.callable.bed
#
done


#### merge / union everything low/high
mkdir -p $OPATH.beds
mkdir -p $OPATH.mergedbeds

# remove outgroup coverage files as they have likely different genomic histories such as duplications etc which could distort the coverage results
ls -1 $OPATH.beds/{Lad4*,Toc*,GCa3*,Par1*,Par2*,SS1*,AR142*,AR223*,Ira3*,USP3*,Copa2*,gem*} | sed "s,^/,rm -f /,g"
# also remove previously published samples belonging to outgroups or the resequence samples from Wang et al. and the original SB and Sb males (different sequencing tec with lots of PCR)

bedops --intersect $OPATH.beds/*.bed > $OPATH.mergedbeds/intersect.bed
bedops -u $OPATH.beds/*.bed > $OPATH.mergedbeds/union.bed
bedops -m $OPATH.beds/*.bed > $OPATH.mergedbeds/merge.bed
bedops --symmdiff $OPATH.beds/*.bed > $OPATH.mergedbeds/symmdiff.bed

bedmap --echo --echo-map-id --delim '\t' --multidelim '\t' $OPATH.mergedbeds/merge.bed $OPATH.mergedbeds/union.bed > $OPATH.mergedbeds/all.bed

for fn in `ls $OPATH.beds/*.bed`; do echo $fn; bedops -n 1 $OPATH.mergedbeds/merge.bed $fn | awk '{ print $0"\tNA" }' | bedops -u - $fn | cut -f4 > $fn.map; done

paste $OPATH.mergedbeds/merge.bed $OPATH.beds/*.bed.map > $OPATH.mergedbeds/mergedbeds.mtx

cat $OPATH.mergedbeds/intersect.bed | wc -l
2055
cat $OPATH.mergedbeds/union.bed | wc -l
362097071
cat $OPATH.mergedbeds/merge.bed | wc -l
202373
cat $OPATH.mergedbeds/symmdiff.bed | wc -l
341679
cat $OPATH.mergedbeds/mergedbeds.mtx | wc -l
7545180
cat $OPATH.mergedbeds/mergedbeds.mtx | grep "^NW*" | wc -l
#185531


# number of columns (3 columns are scf from to)
COLUMNNUMBER=$(echo $(cat $OPATH.mergedbeds/mergedbeds.mtx | grep "^NW*" | head -n 1 | awk '{print NF}')-3 | bc -l) && echo $COLUMNNUMBER
#150
## calculate number of NA's and substract from number of columns --> number of samples in which this region had high coverage
#cat $OPATH.mergedbeds/mergedbeds.mtx | grep "^NW*" | head | grep -o -n 'NA' | cut -d : -f 1 | uniq -c | sed -r "s,^\s+,,g" | cut -d " " -f1 | awk -v CNR=$COLUMNNUMBER ' {print CNR-$1} '
## make new table
paste <(cat $OPATH.mergedbeds/mergedbeds.mtx | grep "^NW*" | cut -f1,2,3) <(cat $OPATH.mergedbeds/mergedbeds.mtx | grep "^NW*" | grep -o -n 'NA' | cut -d : -f 1 | uniq -c | sed -r "s,^\s+,,g" | cut -d " " -f1 | awk -v CNR=$COLUMNNUMBER ' {print CNR-$1} ') | awk ' {OFS="\t"} {print $1,$2,$3,$4,$3-$2} ' > $OPATH.mergedbeds/high.cov.occurence.nr.samples.bed


cat $OPATH.mergedbeds/high.cov.occurence.nr.samples.bed | awk '{s+=$5}{print s}'
# 23054025
# 23.05 Mb

cat $OPATH.mergedbeds/high.cov.occurence.nr.samples.bed | awk '$4 > 5 {sum += $5} END {print sum}'
# 23054013

### regions present in less than 20 samples: only 3 small regions, all other regions consistently occurring in more samples
cat $OPATH.mergedbeds/high.cov.occurence.nr.samples.bed | awk ' $4 < 20 '

NW_011794160.1	132	137	4	5
NW_011794163.1	2236	2243	3	7
NW_011794164.1	29	458	13	429

#######

cat $OPATH.mergedbeds/intersect.bed | awk ' {OFS="\t"} {print $1,$2,$3,$3-$2} ' | awk '{sum += $4} END {print sum}'
# 329919
# 329 kb
