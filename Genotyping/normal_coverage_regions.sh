# Calculate read depth across the genome for each sample
# Regions with very high read depth are not sued because they are likely to represent collapsed repeats

BAMOUT="solenopsis.bams"

F="$BAMOUT"
cd $F

CPUs=4
OPATH="$F/mosdepth"
mkdir -p $OPATH

REF="gng20170922.fa"
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
