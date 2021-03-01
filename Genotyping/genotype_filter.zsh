#!/bin/zsh
############################################################################
# filter SNPs
# Eckart Stolle, April 2020
############################################################################
## script: SNP-genotype-filtering.xxx.zsh
## scripts processes 1 vcf file (1 sample), outputs filtered VCF
## use parallel to run it simultaneous on several files
## uses 2-5 compression threads
## makes a temporary folder in the RAM, thus watch the RAM usage

if [ $# -ne 4 ]; then
    echo $0: usage: ./genotype_filter.zsh INPUT.raw.vcf.gz 4 SAMPLEname OUTPUTFOLDER
	echo "\nINPUT.raw.vcf.gz: vcf file to filter (bgzip-ed/tabix-indexed)"
	echo "\nMinCov: minimum reads coverage to consider a site"
	echo "\nSAMPLEname: Sample Name, will be used for output VCF file"
	echo "\OUTPUTFOLDER: output folder for filtered VCF file"
    exit 1
fi

#set/get variables
VCFINPUT=$1
COVINCL=$2
SAMPLE=$3
FILTEROUTFOLDER=$4
CPUs=5



#### set limits for filtering
FRACTION_AO_RO=0.25
FRACTION_AO_DP=$(echo "1-$FRACTION_AO_RO" | bc -l)
FRACTION_total=$(echo "1/(1/$FRACTION_AO_RO+1)" | bc -l) #0.2
FRACTIONuLIMIT=$(echo "1-$FRACTION_total" | bc -l) #0.8



# !! Check your VCFs genotype declaration
echo "genotype field format/abbreviation used here: GT:GQ:DP:AD:RO:QR:AO:QA:GL"

# make a temporary folder in the RAM
TF="/dev/shm/$SAMPLE"
mkdir -p $TF

# we are changing DIR into that folder to prevent bgzip mangeling all the temporary files created in parallel from different vcf's. those temp names are all the same, causing problems
ORIGFOLDER=$PWD && echo $ORIGFOLDER
cd $TF



tabix -H $VCFINPUT | sed "s#ID=DPR,#ID=AD,#g" > $TF/$SAMPLE.filtered.header
LINESBEFORE=$(zcat $VCFINPUT | grep -v "#" | wc -l)



## set low cov & low GQ to missing, separate missing from nonmissing
zcat $VCFINPUT | sed "s#GT:GQ:DP:DPR:RO:QR:AO:QA:GL#GT:GQ:DP:AD:RO:QR:AO:QA:GL#g" |\
bcftools +setGT - -- -t q -i "FMT/DP < $COVINCL" -n . | bcftools +setGT - -- -t q -i "FMT/GQ < 1 & FMT/DP < 10 " -n . |\
bcftools +setGT - -- -t q -i 'GT="het" & AVG(FMT/AO)=0' -n . |\
tee >(grep -v "#" | grep -P 'GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$|\t\.\|\.|\t\.\/\.|\.:\.:\.:\.:\.|GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.:' |\
sed -r "s#\t\.:\.:\.:\.:\.:\.:\.:\.:\.\$#\t.|.#g" | sed -r "s#\t\.\|\.\$#\t.|.#g" |\
sed -r "s#\t\.\|\.#\t.#g" | sed -r "s#:GL\t\.:[0-9].*#:GL\t\.#g" | sed -r "s#:GL\t\.:-[0-9].*#:GL\t\.#g" > $TF/$SAMPLE.filtered.missing) |\
grep -vP 'GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.$|\t\.\|\.|\t\.\/\.|\.:\.:\.:\.:\.|GT:GQ:DP:AD:RO:QR:AO:QA:GL\t\.:' > $TF/$SAMPLE.vcf

NONMISSINGGENOTYPES=$(cat $TF/$SAMPLE.vcf | grep -v "#" | wc -l) && echo "nonmissing: "$NONMISSINGGENOTYPES
MISSINGGENOTYPES=$(cat $TF/$SAMPLE.filtered.missing | grep -v "#" | wc -l) && echo "missing: "$MISSINGGENOTYPES
FIRSTSUM=$(echo "$NONMISSINGGENOTYPES+$MISSINGGENOTYPES" | bc -l) && echo "$FIRSTSUM of $LINESBEFORE"


#########################################################################################
### separate HMZ from HTZ, set HMZ to phased
#cat $TF/$SAMPLE.vcf | bcftools view -i 'GT="hom"' | bcftools +setGT - -- -t a -n p | grep -v "#" | cut -f10- | tr '\t' '\n' | cut -d":" -f1 | sort | uniq

cat $TF/$SAMPLE.vcf | tee >(bcftools view -i 'GT="hom"' /dev/stdin | bcftools +setGT - -- -t a -n p > $TF/$SAMPLE.filtered.homozygous) >(bcftools view -g miss /dev/stdin > $TF/$SAMPLE.tmp.partiallymissing) | bcftools view -i 'GT="het"' | bcftools +setGT - -- -t a -n u > $TF/$SAMPLE.tmp.htz

HTZ=$(cat $TF/$SAMPLE.tmp.htz | grep -v "#" | wc -l) && echo "HTZ: "$HTZ
HMZ=$(cat $TF/$SAMPLE.filtered.homozygous | grep -v "^#" | wc -l) && echo "HMZ: "$HMZ
PARTMISSING=$(cat $TF/$SAMPLE.tmp.partiallymissing | grep -v "^#" | wc -l) && echo "partially missing: "$PARTMISSING
SECONDSUM=$(echo "$HTZ+$HMZ+$PARTMISSING" | bc -l)
echo "$SECONDSUM of $NONMISSINGGENOTYPES nonmissing"



##########################################################################################
###################### set partially missing to truly missing

### this is the quick solution: set all partially missing GTs to missing (but often partially missing are just wrongly set
#cat $TF/$SAMPLE.tmp.partiallymissing | bcftools +setGT - -- -t ./x -n . | grep -v "#" | sed -r "s#\t\.\/\.:#\t.:#g" > $TF/$SAMPLE.filtered.partiallymissing.missing
#touch $TF/$SAMPLE.filtered.partiallymissing.ref-hmz
#touch $TF/$SAMPLE.filtered.partiallymissing.alt-hmz
#touch $TF/$SAMPLE.filtered.partiallymissing.heterozygous


cat $TF/$SAMPLE.tmp.partiallymissing | bcftools view -e "(FMT/RO+FMT/AO)/FMT/DP<=$FRACTIONuLIMIT" |\
tee >(bcftools view -i "(FMT/RO/FMT/DP)<=$FRACTIONuLIMIT & ((FMT/DP-FMT/RO)/FMT/DP)<=$FRACTIONuLIMIT" | grep -v "^#" |\
sed -r 's#:GL\t\.\|0:|:GL\t\.\|1:#:GL\t0|1:#g' |\
sed -r 's#:GL\t\.\|2:#:GL\t0|2:#g' | sed -r 's#:GL\t\.\|3:#:GL\t0|3:#g' > $TF/$SAMPLE.filtered.partiallymissing.heterozygous) \
>(bcftools view -i "((FMT/DP-FMT/RO)/FMT/DP)>$FRACTIONuLIMIT" | grep -v "^#" |\
sed -E -e 's#\.\|0|0\|\.|\.\|1|1\|\.|\.\|2|2\|\.|\.\|3|3\|\.#1|1#g' > $TF/$SAMPLE.filtered.partiallymissing.alt-hmz) |\
bcftools view -i "(FMT/RO/FMT/DP)>$FRACTIONuLIMIT" /dev/stdin | bcftools +setGT - -- -t a -n 0p | grep -v "^#" > $TF/$SAMPLE.filtered.partiallymissing.ref-hmz

cat $TF/$SAMPLE.tmp.partiallymissing | bcftools view -i "(FMT/RO+FMT/AO)/FMT/DP<=$FRACTIONuLIMIT" |\
bcftools +setGT - -- -t ./x -n . | grep -v "#" | sed -r "s#\t\.\/\.:#\t.:#g" > $TF/$SAMPLE.filtered.partiallymissing.missing


MISSINGGENOTYPES_partiallyMISSING=$(cat $TF/$SAMPLE.filtered.partiallymissing.missing | grep -v "#" | wc -l)
echo "additionally missing: "$MISSINGGENOTYPES_partiallyMISSING
PARTMISSREFHMZ=$(cat $TF/$SAMPLE.filtered.partiallymissing.ref-hmz | grep -v "#" | wc -l) && echo "partially missing = REF HMZ: "$PARTMISSREFHMZ
PARTMISSALTHMZ=$(cat $TF/$SAMPLE.filtered.partiallymissing.alt-hmz | grep -v "#" | wc -l) && echo "partially missing = ALT HMZ: "$PARTMISSALTHMZ
PARTMISSHTZ=$(cat $TF/$SAMPLE.filtered.partiallymissing.heterozygous | grep -v "#" | wc -l) && echo "partially missing = HTZ: "$PARTMISSHTZ
PARTMISSINGSUM=$(echo "$PARTMISSREFHMZ+$PARTMISSALTHMZ+$PARTMISSHTZ+$MISSINGGENOTYPES_partiallyMISSING" | bc -l)
echo "processed partially missing variants $PARTMISSINGSUM of $PARTMISSING"



##########################################################################################
####### HTZ with NUMALT=1
#true HTZ #HMZ (ALT or REF)
rm -f $TF/$SAMPLE.filtered.htz.ref-hmz $TF/$SAMPLE.filtered.htz.alt-hmz $TF/$SAMPLE.filtered.htz.heterozygous

#### one version to process: uses AO vs DP/RO counts. Problem: decomposed Variants have complex AO counts and cannot be addressed properly with bcftools
#cat $TF/$SAMPLE.tmp.htz | awk '/#/{print;next}{if($5 !~ /,/ ){print}}' | bcftools view -i 'NUMALT=1' | bcftools view -i "FMT/RO=0" | grep -v "^#" | sed -E -e 's#0\/1|1\/0|0\|1|1\|0#1|1#g' | sed -E -e 's#0\/2|2\/0|0\|2|2\|0#2|2#g' > $TF/$SAMPLE.filtered.htz.alt-hmz
#cat $TF/$SAMPLE.tmp.htz | awk '/#/{print;next}{if($5 !~ /,/ ){print}}' | bcftools view -i 'NUMALT=1' | bcftools view -e "FMT/RO=0" |\
#tee >(bcftools view -i "(FMT/RO/FMT/AO)>=$FRACTION_AO_RO && (FMT/AO/FMT/RO)>=$FRACTION_AO_RO" /dev/stdin > $TF/$SAMPLE.filtered.htz.heterozygous)\
# >(bcftools view -e "(FMT/RO/FMT/AO)>=$FRACTION_AO_RO && (FMT/AO/FMT/RO)>=$FRACTION_AO_RO" /dev/stdin | bcftools view -i "(FMT/RO/FMT/AO)<$FRACTION_AO_RO" | grep -v "^#" | sed -E -e 's#0\/1|1\/0|0\|1|1\|0#1|1#g' | sed -E -e 's#0\/2|2\/0|0\|2|2\|0#2|2#g' >> $TF/$SAMPLE.filtered.htz.alt-hmz) |\
#bcftools view -e "(FMT/RO/FMT/AO)>=$FRACTION_AO_RO && (FMT/AO/FMT/RO)>=$FRACTION_AO_RO" | bcftools view -e "(FMT/RO/FMT/AO)<$FRACTION_AO_RO" | bcftools +setGT - -- -t a -n 0p | grep -v "^#" > $TF/$SAMPLE.filtered.htz.ref-hmz


##### Simple filtering version: based on RO vs DP counts, should in reality
# (with , as here, limiting analysis to biallelic sites, work fine and better than above)
SIMPLEHTZ=$(cat $TF/$SAMPLE.tmp.htz | awk '/#/{print;next}{if($5 !~ /,/ ){print}}' | bcftools view -i 'NUMALT=1' | grep -v "#" | wc -l)
## RO=0 --> inevitably this biallelic SNP is ALT HMZ then
cat $TF/$SAMPLE.tmp.htz | awk '/#/{print;next}{if($5 !~ /,/ ){print}}' | bcftools view -i 'NUMALT=1' | bcftools view -i "FMT/RO=0" | grep -v "^#" | sed -E -e 's#0\/1|1\/0|0\|1|1\|0#1|1#g' | sed -E -e 's#0\/2|2\/0|0\|2|2\|0#2|2#g' > $TF/$SAMPLE.filtered.htz.alt-hmz
## high REF numbers: REF HMZ
cat $TF/$SAMPLE.tmp.htz | awk '/#/{print;next}{if($5 !~ /,/ ){print}}' | bcftools view -i 'NUMALT=1' | bcftools view -e "FMT/RO=0" | bcftools view -i "FMT/RO/FMT/DP>=$FRACTION_AO_DP" | bcftools +setGT - -- -t a -n 0p | grep -v "^#" > $TF/$SAMPLE.filtered.htz.ref-hmz
## low REF numbers: ALT HMZ (REF not more than 75% && REF less than 25% --> must be ALT HMZ, otherwise REF should be around 50%
cat $TF/$SAMPLE.tmp.htz | awk '/#/{print;next}{if($5 !~ /,/ ){print}}' | bcftools view -i 'NUMALT=1' | bcftools view -e "FMT/RO=0" | bcftools view -e "FMT/RO/FMT/DP>=$FRACTION_AO_DP" | bcftools view -i "FMT/RO/FMT/DP<=$FRACTION_AO_RO" | grep -v "^#" | sed -E -e 's#0\/1|1\/0|0\|1|1\|0#1|1#g' | sed -E -e 's#0\/2|2\/0|0\|2|2\|0#2|2#g' >> $TF/$SAMPLE.filtered.htz.alt-hmz
## REF between 25-75% --> HTZ
cat $TF/$SAMPLE.tmp.htz | awk '/#/{print;next}{if($5 !~ /,/ ){print}}' | bcftools view -i 'NUMALT=1' | bcftools view -e "FMT/RO=0" | bcftools view -e "FMT/RO/FMT/DP>=$FRACTION_AO_DP" | bcftools view -e "FMT/RO/FMT/DP<=$FRACTION_AO_RO" | grep -v "^#" > $TF/$SAMPLE.filtered.htz.heterozygous

FALSEHTZALT=$(cat $TF/$SAMPLE.filtered.htz.alt-hmz | grep -v "^#" | wc -l) && echo "false HTZ (ALT): "$FALSEHTZALT
FALSEHTZREF=$(cat $TF/$SAMPLE.filtered.htz.ref-hmz | grep -v "^#" | wc -l) && echo "false HTZ (REF): "$FALSEHTZREF
TRUEHTZ=$(cat $TF/$SAMPLE.filtered.htz.heterozygous | grep -v "^#" | wc -l) && echo "true HTZ: "$TRUEHTZ
HTZSUM=$(echo "$FALSEHTZALT+$FALSEHTZREF+$TRUEHTZ" | bc -l) && echo "simple HTZ processed: $HTZSUM of $SIMPLEHTZ"




########### multialleic sites

echo "taking out multiallelic sites for processing"

cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.tmp.htz3
cat $TF/$SAMPLE.tmp.htz | awk '/#/{print;next}{if($5 ~ /,/ ){print}}' | tee >(grep -vP '1\|2|2\|1|1\|3|3\|1|2\|3|3\|2|1/2|2/1|1/3|3/1|2/3|3/2' > $TF/$SAMPLE.tmp.htz2) |\
grep -P '1\|2|2\|1|1\|3|3\|1|2\|3|3\|2|1/2|2/1|1/3|3/1|2/3|3/2' >> $TF/$SAMPLE.tmp.htz3

HTZ2=$(cat $TF/$SAMPLE.tmp.htz2 | grep -v "#" | wc -l)
HTZ3=$(cat $TF/$SAMPLE.tmp.htz3 | grep -v "#" | wc -l)
HTZsplitSUM=$(echo "$HTZSUM+$HTZ2+$HTZ3" | bc -l)
echo "$HTZsplitSUM of $HTZ total HTZ variants split into $HTZSUM simple HTZs, $HTZ2 multiallelic (REF/ALT) and $HTZ3 multiallelic (ALT-ALT) HTZs (splitted sum: $HTZsplitSUM"








### HTZ2
rm -f $TF/$SAMPLE.filtered.htz2.ref-hmz $TF/$SAMPLE.filtered.htz2.alt-hmz $TF/$SAMPLE.filtered.htz2.heterozygous

cat $TF/$SAMPLE.tmp.htz2 | bcftools view -i "FMT/RO=0" | grep -v "^#" | sed -E -e 's#0\/1|1\/0|0\|1|1\|0#1|1#g' | sed -E -e 's#0\/2|2\/0|0\|2|2\|0#2|2#g' | sed -E -e 's#0\/3|3\/0|0\|3|3\|0#3|3#g' >> $TF/$SAMPLE.filtered.htz2.alt-hmz

cat $TF/$SAMPLE.tmp.htz2 | bcftools view -e "FMT/RO=0" | tee \
>(bcftools view -i "(FMT/RO/FMT/DP)>=$FRACTION_total & (FMT/RO/FMT/DP)<=$FRACTIONuLIMIT" |\
grep -v "^#" >> $TF/$SAMPLE.filtered.htz2.heterozygous) \
>(bcftools view -i "(FMT/RO/FMT/DP)<$FRACTION_total" | grep -v "^#" |\
sed -E -e 's#0\/1|1\/0|0\|1|1\|0#1|1#g' | sed -E -e 's#0\/2|2\/0|0\|2|2\|0#2|2#g' | sed -E -e 's#0\/3|3\/0|0\|3|3\|0#3|3#g' >> $TF/$SAMPLE.filtered.htz2.alt-hmz) |\
bcftools view  -i  "(FMT/RO/FMT/DP)>$FRACTIONuLIMIT" | bcftools +setGT - -- -t a -n 0p |\
grep -v "^#" >> $TF/$SAMPLE.filtered.htz2.ref-hmz

HTZ2REFHMZ=$(cat $TF/$SAMPLE.filtered.htz2.ref-hmz | grep -v "#" | wc -l) && echo "HTZ2 = REF HMZ: "$HTZ2REFHMZ
HTZ2ALTHMZ=$(cat $TF/$SAMPLE.filtered.htz2.alt-hmz | grep -v "#" | wc -l) && echo "HTZ2 = ALT HMZ: "$HTZ2ALTHMZ
HTZ2HTZ=$(cat $TF/$SAMPLE.filtered.htz2.heterozygous | grep -v "#" | wc -l) && echo "HTZ2 = HTZ: "$HTZ2HTZ
HTZ2SUM=$(echo "$HTZ2REFHMZ+$HTZ2ALTHMZ+$HTZ2HTZ" | bc -l)
echo "processed HTZ2 variants $HTZ2SUM of $HTZ2"






### HTZ3
rm -f $TF/$SAMPLE.filtered.htz3.heterozygous $TF/$SAMPLE.filtered.htz3.ref-hmz $TF/$SAMPLE.filtered.htz3.heterozygous

cat $TF/$SAMPLE.tmp.htz3 | bcftools view -i "FMT/RO=0" | grep -v "^#" >> $TF/$SAMPLE.filtered.htz3.heterozygous

cat $TF/$SAMPLE.tmp.htz3 | bcftools view -e "FMT/RO=0" | tee \
>(bcftools view -i "(FMT/RO/FMT/DP)<=$FRACTIONuLIMIT & ((FMT/DP-FMT/RO)/FMT/DP)<=$FRACTIONuLIMIT" /dev/stdin |\
grep -v "^#" >> $TF/$SAMPLE.filtered.htz3.heterozygous) \
>(bcftools view -i "((FMT/DP-FMT/RO)/FMT/DP)>$FRACTIONuLIMIT" /dev/stdin | grep -v "^#" >> $TF/$SAMPLE.filtered.htz3.heterozygous) |\
bcftools view -i "(FMT/RO/FMT/DP)>$FRACTIONuLIMIT" /dev/stdin | bcftools +setGT - -- -t a -n 0p |\
grep -v "^#" >> $TF/$SAMPLE.filtered.htz3.ref-hmz

## create formerly used file here so that its present, I stopped using it because multiallelic sites are too difficult, so I leave them as heterozygous. They would need to require special attention to fix the genotypes in which an ALT allele is HMZ vs true ALT-ALT HTZ.
touch $TF/$SAMPLE.filtered.htz3.alt-hmz

HTZ3REFHMZ=$(cat $TF/$SAMPLE.filtered.htz3.ref-hmz | grep -v "#" | wc -l) && echo "HTZ3 = REF HMZ: "$HTZ3REFHMZ
HTZ3ALTHMZ=$(cat $TF/$SAMPLE.filtered.htz3.alt-hmz | grep -v "#" | wc -l) && echo "HTZ3 = ALT HMZ: "$HTZ3ALTHMZ
HTZ3HTZ=$(cat $TF/$SAMPLE.filtered.htz3.heterozygous | grep -v "#" | wc -l) && echo "HTZ3 = HTZ: "$HTZ3HTZ
HTZ3SUM=$(echo "$HTZ3REFHMZ+$HTZ3ALTHMZ+$HTZ3HTZ" | bc -l)
echo "processed HTZ3 variants $HTZ3SUM of $HTZ3"


### fuse HTZ-derived files to count
cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.filtered.htz.fused.vcf
cat \
$TF/$SAMPLE.filtered.htz.ref-hmz $TF/$SAMPLE.filtered.htz.alt-hmz $TF/$SAMPLE.filtered.htz.heterozygous \
$TF/$SAMPLE.filtered.htz2.ref-hmz $TF/$SAMPLE.filtered.htz2.alt-hmz $TF/$SAMPLE.filtered.htz2.heterozygous \
$TF/$SAMPLE.filtered.htz3.ref-hmz $TF/$SAMPLE.filtered.htz3.alt-hmz $TF/$SAMPLE.filtered.htz3.heterozygous |\
grep -v "#" >> $TF/$SAMPLE.filtered.htz.fused.vcf
vcfstreamsort -a $TF/$SAMPLE.filtered.htz.fused.vcf | bgzip -f -@ 5 -c /dev/stdin > $TF/$SAMPLE.filtered.htz.fused.vcf.gz && tabix -fp vcf $TF/$SAMPLE.filtered.htz.fused.vcf.gz

FILTEREDHTZ=$(cat $TF/$SAMPLE.filtered.htz.fused.vcf | grep -v "#" | wc -l)
echo "processed $FILTEREDHTZ of $HTZ HTZ variants"



##intersect to check for leftovers
rm -f $TF/$SAMPLE.tmp.htz-leftover

cat $TF/$SAMPLE.tmp.htz | bgzip -f -@ 5 -c /dev/stdin > $TF/$SAMPLE.tmp.htz.gz && tabix -fp vcf $TF/$SAMPLE.tmp.htz.gz

bcftools isec -p $TF -Oz $TF/$SAMPLE.tmp.htz.gz $TF/$SAMPLE.filtered.htz.fused.vcf.gz
zcat $TF/0001.vcf.gz | grep -v "^#" > $TF/$SAMPLE.tmp.htz-leftovertmp
if [ -s $TF/$SAMPLE.tmp.htz-leftovertmp ]
then
echo "leftover variants detected, set to missing"
zcat $TF/0001.vcf.gz | bcftools +setGT - -- -t a -n . | grep -v "^#" > $TF/$SAMPLE.tmp.htz-leftover
else
echo "no leftover, weird heterozygous variants"
touch $TF/$SAMPLE.tmp.htz-leftover
fi




##################################################

## combine pieces
cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.filtered.combined
cat \
$TF/$SAMPLE.filtered.homozygous \
$TF/$SAMPLE.filtered.htz.alt-hmz \
$TF/$SAMPLE.filtered.htz.heterozygous \
$TF/$SAMPLE.filtered.htz.ref-hmz \
$TF/$SAMPLE.filtered.htz2.alt-hmz \
$TF/$SAMPLE.filtered.htz2.heterozygous \
$TF/$SAMPLE.filtered.htz2.ref-hmz \
$TF/$SAMPLE.filtered.htz3.alt-hmz \
$TF/$SAMPLE.filtered.htz3.heterozygous \
$TF/$SAMPLE.filtered.htz3.ref-hmz \
$TF/$SAMPLE.filtered.partiallymissing.ref-hmz \
$TF/$SAMPLE.filtered.partiallymissing.alt-hmz \
$TF/$SAMPLE.filtered.partiallymissing.missing \
$TF/$SAMPLE.filtered.partiallymissing.heterozygous \
$TF/$SAMPLE.filtered.missing \
$TF/$SAMPLE.tmp.htz-leftover |\
grep -v "#" | sort -k1,1 -k2,2g | sed -r "s#\t\./\.#\t.#g" | sed -r "s#\t\.\|\.#\t.#g" >> $TF/$SAMPLE.filtered.combined
echo -n "combined varinats: "
cat $TF/$SAMPLE.filtered.combined | grep -v "#" | wc -l

cat $TF/$SAMPLE.filtered.header > $TF/$SAMPLE.filtered.combined.HTZ
cat \
$TF/$SAMPLE.filtered.htz.heterozygous \
$TF/$SAMPLE.filtered.htz2.heterozygous \
$TF/$SAMPLE.filtered.htz3.heterozygous \
$TF/$SAMPLE.filtered.partiallymissing.heterozygous |\
grep -v "#" | sort -k1,1 -k2,2g | sed -r "s#\t\.\|\.#\t.#g" >> $TF/$SAMPLE.filtered.combined.HTZ
echo -n "HTZ varinats: "
cat $TF/$SAMPLE.filtered.combined.HTZ | grep -v "#" | wc -l

echo "combined VCF parts"



##### sort, HMZ phased, HTZ unphased
vcfstreamsort -a $TF/$SAMPLE.filtered.combined |\
vcffixup - | bcftools +setGT - -- -t q -i 'GT="hom" ' -n p | bcftools +setGT - -- -t q -i 'GT="het" ' -n u | bcftools +setGT - -- -t ./. -n . |\
sed -r "s#\t\.\|\.\|\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" | sed -r "s#\t\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" |\
sed -r "s#\t\.\|\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" | sed -r "s#\t\./\.#\t.#g" | sed -r "s#\t\.\|\.#\t.#g" |\
bgzip -f -@ 5 -c /dev/stdin > $FILTEROUTFOLDER/$SAMPLE.filtered.vcf.gz
tabix -fp vcf $FILTEROUTFOLDER/$SAMPLE.filtered.vcf.gz
echo -n "filtered output file: "
zcat $FILTEROUTFOLDER/$SAMPLE.filtered.vcf.gz | grep -v "#" | wc -l


## HTZ to missing
vcfstreamsort -a $TF/$SAMPLE.filtered.combined |\
vcffixup - | bcftools +setGT - -- -t q -i 'GT="hom" ' -n p | bcftools +setGT - -- -t q -i 'GT="het" ' -n . | bcftools +setGT - -- -t ./. -n . |\
sed -r "s#\t\.\|\.\|\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" | sed -r "s#\t\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" |\
sed -r "s#\t\.\|\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" | sed -r "s#\t\./\.#\t.#g"| sed -r "s#\t\.\|\.#\t.#g" |\
bgzip -f -@ 5 -c /dev/stdin > $FILTEROUTFOLDER/$SAMPLE.filtered.heterozygous2missing.vcf.gz
tabix -fp vcf $FILTEROUTFOLDER/$SAMPLE.filtered.heterozygous2missing.vcf.gz
echo -n "filtered output file (HTZ2missing): "
zcat $FILTEROUTFOLDER/$SAMPLE.filtered.heterozygous2missing.vcf.gz | grep -v "#" | wc -l



## HTZ separately
vcfstreamsort -a $TF/$SAMPLE.filtered.combined.HTZ |\
vcffixup - | bcftools +setGT - -- -t q -i 'GT="het" ' -n u | bcftools +setGT - -- -t ./. -n . |\
bgzip -f -@ 5 -c /dev/stdin > $FILTEROUTFOLDER/$SAMPLE.filtered.HTZ.vcf.gz
tabix -fp vcf $FILTEROUTFOLDER/$SAMPLE.filtered.HTZ.vcf.gz
echo -n "HTZ output (separately): "
zcat $FILTEROUTFOLDER/$SAMPLE.filtered.HTZ.vcf.gz | grep -v "#" | wc -l




####  choose randomly one of the two alleles at heterozygous loci (especially for diploids, 10x randomly)
####  alternative: pick most common allele at heterozygous loci (best for pools of individuals)
SHUFFLE="shuffle_file"    # a files generated with random 0 or 1 # 4059290 lines
#cat $SHUFFLE | wc -l

##### heterozygous calls file, choose majority allele to set HMZ Genotype

cat $TF/$SAMPLE.filtered.combined.HTZ | grep -v "#" > $TF/$SAMPLE.filtered.combined.HTZbody
HTZFILE="$TF/$SAMPLE.filtered.combined.HTZbody"
M=$(cat $HTZFILE | wc -l) && echo $M

if [ -s $HTZFILE ]
then
        echo "HTZ file not empty, making HTZ2random vcf files"

# choose a random allele for each diploid site, based on a list of rendomly generated 0s and 1s (via shuffle)

p=1
REPEATS=3
for (( p = 1 ; p < $REPEATS+1 ; p++))
 do
shuf $SHUFFLE | head -$M > $TF/shuffle.$p
cat $HTZFILE | cut -f10 | cut -d":" -f1 | sed "s/\//,/g" | sed "s/|/,/g" | cut -d"," -f1 > $HTZFILE.B1
cat $HTZFILE | cut -f10 | cut -d":" -f1 | sed "s/\//,/g" | sed "s/|/,/g" | cut -d"," -f2 > $HTZFILE.B2
cat $HTZFILE | cut -f10 | cut -d":" -f2- > $HTZFILE.C
paste $HTZFILE $TF/shuffle.$p $HTZFILE.B1 $HTZFILE.B2 $HTZFILE.C > $HTZFILE.tmp.$p
#originalVCF	shuffle$11	allele1$12	allele2$13	last_part$14
cat $HTZFILE.tmp.$p | awk 'BEGIN{OFS=FS="\t"}; {if($11 == 0) {ALLELE=$12;} else if($11 == 1) {ALLELE=$13;}; {print $1,$2,$3,$4,$5,$6,$7,($8";SELECTEDALLELE"ALLELE),$9,(ALLELE"|"ALLELE":"$14);} }' > $TF/$SAMPLE.filtered.heterozygous2randomallele.$p
rm -f $HTZFILE.tmp.$p $HTZFILE.B1 $HTZFILE.B2 $HTZFILE.C $TF/shuffle.$p


cat $TF/$SAMPLE.filtered.header | grep -v "^#CHROM" > $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p
print '##INFO=<ID=SELECTEDALLELE0,Number=A,Type=String,Description="Heterozygous call, chosen allele: 0">' >> $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p
print '##INFO=<ID=SELECTEDALLELE1,Number=A,Type=String,Description="Heterozygous call, chosen allele: 1">' >> $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p
print '##INFO=<ID=SELECTEDALLELE2,Number=A,Type=String,Description="Heterozygous call, chosen allele: 2">' >> $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p
print '##INFO=<ID=SELECTEDALLELE3,Number=A,Type=String,Description="Heterozygous call, chosen allele: 3">' >> $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p
cat $TF/$SAMPLE.filtered.header | grep "^#CHROM" | sed "s/$/&-$p/" >> $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p
cat $TF/$SAMPLE.filtered.heterozygous2randomallele.$p | sort -k1,1 -k2,2g >> $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p


# store these changed HTZ in separate file
vcfstreamsort -a $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p | grep -v "Dummy" |\
sed -r "s#\t\./\.#\t.#g"| sed -r "s#\t\.\|\.#\t.#g" |\
bgzip -f -@ 5 -c /dev/stdin > $FILTEROUTFOLDER/$SAMPLE.filtered.heterozygous2randomallele.$p.vcf.gz
tabix -fp vcf $FILTEROUTFOLDER/$SAMPLE.filtered.heterozygous2randomallele.$p.vcf.gz

#output full vcf with the heterozygous-derived alleles
cat $TF/$SAMPLE.filtered.header | grep -v "^#CHROM" > $TF/$SAMPLE.filtered.combined.$p
print '##INFO=<ID=SELECTEDALLELE0,Number=A,Type=String,Description="Heterozygous call, chosen allele: 0">' >> $TF/$SAMPLE.filtered.combined.$p
print '##INFO=<ID=SELECTEDALLELE1,Number=A,Type=String,Description="Heterozygous call, chosen allele: 1">' >> $TF/$SAMPLE.filtered.combined.$p
print '##INFO=<ID=SELECTEDALLELE2,Number=A,Type=String,Description="Heterozygous call, chosen allele: 2">' >> $TF/$SAMPLE.filtered.combined.$p
print '##INFO=<ID=SELECTEDALLELE3,Number=A,Type=String,Description="Heterozygous call, chosen allele: 3">' >> $TF/$SAMPLE.filtered.combined.$p
cat $TF/$SAMPLE.filtered.header | grep "^#CHROM" | sed "s/$/&-$p/" >> $TF/$SAMPLE.filtered.combined.$p

cat \
$TF/$SAMPLE.filtered.homozygous \
$TF/$SAMPLE.filtered.heterozygous2randomallele.$p \
$TF/$SAMPLE.filtered.htz.alt-hmz \
$TF/$SAMPLE.filtered.htz.ref-hmz \
$TF/$SAMPLE.filtered.htz2.alt-hmz \
$TF/$SAMPLE.filtered.htz2.ref-hmz \
$TF/$SAMPLE.filtered.htz3.alt-hmz \
$TF/$SAMPLE.filtered.htz3.ref-hmz \
$TF/$SAMPLE.filtered.partiallymissing.ref-hmz \
$TF/$SAMPLE.filtered.partiallymissing.alt-hmz \
$TF/$SAMPLE.filtered.partiallymissing.missing \
$TF/$SAMPLE.filtered.missing \
$TF/$SAMPLE.tmp.htz-leftover | grep -v "#" | sort -k1,1 -k2,2g |\
sed -r "s#\t\.\|\.#\t.#g" | sed -r "s#\t\./\.#\t.#g" >> $TF/$SAMPLE.filtered.combined.$p

vcfstreamsort -a $TF/$SAMPLE.filtered.combined.$p |\
vcffixup - | bcftools +setGT - -- -t q -i 'GT="hom" ' -n p | bcftools +setGT - -- -t q -i 'GT="het" ' -n u | bcftools +setGT - -- -t ./. -n . |\
grep -v "Dummy" |\
sed -r "s#\t\.\|\.\|\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" | sed -r "s#\t\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" |\
sed -r "s#\t\.\|\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" | sed -r "s#\t\./\.#\t.#g" | sed -r "s#\t\.\|\.#\t.#g" |\
bgzip -f -@ 5 -c /dev/stdin > $FILTEROUTFOLDER/$SAMPLE.filtered.$p.vcf.gz
tabix -fp vcf $FILTEROUTFOLDER/$SAMPLE.filtered.$p.vcf.gz
#echo -n "outputfile HTZ set to major allele: "
#zcat $FILTEROUTFOLDER/$SAMPLE.filtered.$p.vcf.gz | grep -v "#" | wc -l

done


####   same thing, but choose majority allele instead of random (more applicable for pools of individuals)
p="majorityallele"
shuf $SHUFFLE | head -$M > $TF/shuffle.$p
cat $HTZFILE | cut -f10 | cut -d":" -f1 | sed "s/\//,/g" | sed "s/|/,/g" | cut -d"," -f1 > $HTZFILE.B1
cat $HTZFILE | cut -f10 | cut -d":" -f1 | sed "s/\//,/g" | sed "s/|/,/g" | cut -d"," -f2 > $HTZFILE.B2
cat $HTZFILE | cut -f10 | cut -d":" -f3 > $HTZFILE.C-DP
cat $HTZFILE | cut -f10 | cut -d":" -f5 > $HTZFILE.C-RO
cat $HTZFILE | cut -f10 | cut -d":" -f2- > $HTZFILE.C
paste $HTZFILE $TF/shuffle.$p $HTZFILE.B1 $HTZFILE.B2 $HTZFILE.C $HTZFILE.C-DP $HTZFILE.C-RO > $HTZFILE.tmp.$p
#originalVCF	shuffle$11	allele1$12	allele2$13	last_part$14	DP$15	RO$16
cat $HTZFILE.tmp.$p | awk 'BEGIN{OFS=FS="\t"}; {if($16/$15 > 0.5) {ALLELE=$12;} else if($16/$15 < 0.5) {ALLELE=$13;} else if($16/$15 == 0.5 && $11 == 0) {ALLELE=$12;} else if($16/$15 == 0.5 && $11 == 1) {ALLELE=$13;}; {print $1,$2,$3,$4,$5,$6,$7,($8";SELECTEDALLELE"ALLELE),$9,(ALLELE"|"ALLELE":"$14);} }' > $TF/$SAMPLE.filtered.heterozygous2randomallele.$p
rm -f $HTZFILE.tmp.$p $HTZFILE.B1 $HTZFILE.B2 $HTZFILE.C $HTZFILE.C-DP $HTZFILE.C-RO $TF/shuffle.$p

#output only the heterozygous-derived alleles
cat $TF/$SAMPLE.filtered.header | grep -v "^#CHROM" > $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p
print '##INFO=<ID=SELECTEDALLELE0,Number=A,Type=String,Description="Heterozygous call, chosen allele: 0">' >> $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p
print '##INFO=<ID=SELECTEDALLELE1,Number=A,Type=String,Description="Heterozygous call, chosen allele: 1">' >> $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p
print '##INFO=<ID=SELECTEDALLELE2,Number=A,Type=String,Description="Heterozygous call, chosen allele: 2">' >> $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p
print '##INFO=<ID=SELECTEDALLELE3,Number=A,Type=String,Description="Heterozygous call, chosen allele: 3">' >> $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p
cat $TF/$SAMPLE.filtered.header | grep "^#CHROM" | sed "s/$/&-$p/" >> $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p
cat $TF/$SAMPLE.filtered.heterozygous2randomallele.$p | sort -k1,1 -k2,2g >> $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p

# store these changed HTZ in separate file
vcfstreamsort -a $TF/$SAMPLE.filtered.combined.heterozygous2randomallele.$p | grep -v "Dummy" |\
sed -r "s#\t\./\.#\t.#g"| sed -r "s#\t\.\|\.#\t.#g" |\
bgzip -f -@ 5 -c /dev/stdin > $FILTEROUTFOLDER/$SAMPLE.filtered.heterozygous2randomallele.$p.vcf.gz
tabix -fp vcf $FILTEROUTFOLDER/$SAMPLE.filtered.heterozygous2randomallele.$p.vcf.gz

#output full vcf with the heterozygous-derived alleles
cat $TF/$SAMPLE.filtered.header | grep -v "^#CHROM" > $TF/$SAMPLE.filtered.combined.$p
print '##INFO=<ID=SELECTEDALLELE0,Number=A,Type=String,Description="Heterozygous call, chosen allele: 0">' >> $TF/$SAMPLE.filtered.combined.$p
print '##INFO=<ID=SELECTEDALLELE1,Number=A,Type=String,Description="Heterozygous call, chosen allele: 1">' >> $TF/$SAMPLE.filtered.combined.$p
print '##INFO=<ID=SELECTEDALLELE2,Number=A,Type=String,Description="Heterozygous call, chosen allele: 2">' >> $TF/$SAMPLE.filtered.combined.$p
print '##INFO=<ID=SELECTEDALLELE3,Number=A,Type=String,Description="Heterozygous call, chosen allele: 3">' >> $TF/$SAMPLE.filtered.combined.$p
cat $TF/$SAMPLE.filtered.header | grep "^#CHROM" | sed "s/$/&-$p/" >> $TF/$SAMPLE.filtered.combined.$p

cat \
$TF/$SAMPLE.filtered.homozygous \
$TF/$SAMPLE.filtered.heterozygous2randomallele.$p \
$TF/$SAMPLE.filtered.htz.alt-hmz \
$TF/$SAMPLE.filtered.htz.ref-hmz \
$TF/$SAMPLE.filtered.htz2.alt-hmz \
$TF/$SAMPLE.filtered.htz2.ref-hmz \
$TF/$SAMPLE.filtered.htz3.alt-hmz \
$TF/$SAMPLE.filtered.htz3.ref-hmz \
$TF/$SAMPLE.filtered.partiallymissing.ref-hmz \
$TF/$SAMPLE.filtered.partiallymissing.alt-hmz \
$TF/$SAMPLE.filtered.partiallymissing.missing \
$TF/$SAMPLE.filtered.missing \
$TF/$SAMPLE.tmp.htz-leftover | grep -v "#" | sort -k1,1 -k2,2g |\
sed -r "s#\t\.\|\.#\t.#g" | sed -r "s#\t\./\.#\t.#g" >> $TF/$SAMPLE.filtered.combined.$p

vcfstreamsort -a $TF/$SAMPLE.filtered.combined.$p |\
vcffixup - | bcftools +setGT - -- -t q -i 'GT="hom" ' -n p | bcftools +setGT - -- -t q -i 'GT="het" ' -n u | bcftools +setGT - -- -t ./. -n . |\
grep -v "Dummy" |\
sed -r "s#\t\.\|\.\|\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" | sed -r "s#\t\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" |\
sed -r "s#\t\.\|\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" | sed -r "s#\t\./\.#\t.#g" | sed -r "s#\t\.\|\.#\t.#g" |\
bgzip -f -@ 5 -c /dev/stdin > $FILTEROUTFOLDER/$SAMPLE.filtered.$p.vcf.gz
tabix -fp vcf $FILTEROUTFOLDER/$SAMPLE.filtered.$p.vcf.gz
echo -n "outputfile HTZ set to major allele: "
zcat $FILTEROUTFOLDER/$SAMPLE.filtered.$p.vcf.gz | grep -v "#" | wc -l

else
        echo "HTZ file empty, not making HTZ2random files"

fi





###############################################

FUSED=$(cat $TF/$SAMPLE.filtered.combined | grep -v "#" | wc -l)
FINAL=$(zcat $FILTEROUTFOLDER/$SAMPLE.filtered.vcf.gz | grep -v "#" | wc -l)


echo "\
IN nonmissing missing sum \n\
$LINESBEFORE input\n\
$NONMISSINGGENOTYPES nonmissing\n\
$MISSINGGENOTYPES missing\n\
$FIRSTSUM sum\n\
 \n\
NONMISSINGGENOTYPES HMZ HTZ HTZ2/3(multialelic) SUM \n\
$HMZ homozygous\n\
$HTZ heterozygous (total), of which \n\
$HTZSUM are 0/1 REF/ALT \n\
$HTZ2 are multiallelic REF/ALT\n\
$HTZ3 are multiallelic ALT/ALT\n\
$PARTMISSING partially missing GT\n\
$SECONDSUM sum\n\
 \n\
HTZ REF-HMZ ALT_HMZ HTZ SUM \n\
$FALSEHTZREF heterozygous locus actually homozygous REF\n\
$FALSEHTZALT heterozygous locus actually homozygous ALT\n\
$TRUEHTZ true heterozygous\n\
$HTZSUM sum\n\
 \n\
HTZ2 REF-HMZ ALT_HMZ HTZ SUM \n\
$HTZ2REFHMZ heterozygous locus actually homozygous REF\n\
$HTZ2ALTHMZ heterozygous locus actually homozygous ALT\n\
$HTZ2HTZ true heterozygous\n\
$HTZ2SUM sum\n\
 \n\
HTZ3 REF-HMZ ALT_HMZ HTZ SUM \n\
$HTZ3REFHMZ heterozygous locus actually homozygous REF\n\
$HTZ3ALTHMZ heterozygous locus actually homozygous ALT\n\
$HTZ3HTZ true heterozygous\n\
$HTZ3SUM sum\n\
\n\
partially missing REF-HMZ ALT_HMZ HTZ SUM \n\
$PARTMISSREFHMZ partially missing GT locus actually homozygous REF\n\
$PARTMISSALTHMZ partially missing GT locus actually homozygous ALT\n\
$PARTMISSHTZ partially missing GT locus actually heterozygous\n\
$MISSINGGENOTYPES_partiallyMISSING partially missing GT set to missing for too low coverage\n\
$PARTMISSINGSUM sum\n\
"



echo "\
IN nonmissing missing sum \n\
$LINESBEFORE input\n\
$NONMISSINGGENOTYPES nonmissing\n\
$MISSINGGENOTYPES missing\n\
$FIRSTSUM sum\n\
 \n\
NONMISSINGGENOTYPES HMZ HTZ HTZ2/3(multialelic) SUM \n\
$HMZ homozygous\n\
$HTZ heterozygous (total), of which \n\
$HTZSUM are 0/1 REF/ALT \n\
$HTZ2 are multiallelic REF/ALT\n\
$HTZ3 are multiallelic ALT/ALT\n\
$PARTMISSING partially missing GT\n\
$SECONDSUM sum\n\
 \n\
HTZ REF-HMZ ALT_HMZ HTZ SUM \n\
$FALSEHTZREF heterozygous locus actually homozygous REF\n\
$FALSEHTZALT heterozygous locus actually homozygous ALT\n\
$TRUEHTZ true heterozygous\n\
$HTZSUM sum\n\
 \n\
HTZ2 REF-HMZ ALT_HMZ HTZ SUM \n\
$HTZ2REFHMZ heterozygous locus actually homozygous REF\n\
$HTZ2ALTHMZ heterozygous locus actually homozygous ALT\n\
$HTZ2HTZ true heterozygous\n\
$HTZ2SUM sum\n\
 \n\
HTZ3 REF-HMZ ALT_HMZ HTZ SUM \n\
$HTZ3REFHMZ heterozygous locus actually homozygous REF\n\
$HTZ3ALTHMZ heterozygous locus actually homozygous ALT\n\
$HTZ3HTZ true heterozygous\n\
$HTZ3SUM sum\n\
\n\
partially missing REF-HMZ ALT_HMZ HTZ SUM \n\
$PARTMISSREFHMZ partially missing GT locus actually homozygous REF\n\
$PARTMISSALTHMZ partially missing GT locus actually homozygous ALT\n\
$PARTMISSHTZ partially missing GT locus actually heterozygous\n\
$MISSINGGENOTYPES_partiallyMISSING partially missing GT set to missing for too low coverage\n\
$PARTMISSINGSUM sum\n" > $FILTEROUTFOLDER/$SAMPLE.filtered.txt

TOTALMISSING=$(echo "$MISSINGGENOTYPES+$MISSINGGENOTYPES_partiallyMISSING" | bc -l)
TOTALHTZ=$(echo "$TRUEHTZ+$HTZ2HTZ+$HTZ3HTZ+$PARTMISSHTZ" | bc -l)
TOTALHMZ=$(echo "$HMZ+$HTZ2REFHMZ+$HTZ2ALTHMZ+$HTZ3REFHMZ+$HTZ3ALTHMZ+$FALSEHTZREF+$FALSEHTZALT+$PARTMISSALTHMZ+$PARTMISSREFHMZ" | bc -l)
TOTAL=$(echo "$TOTALHTZ+$TOTALHMZ+$TOTALMISSING" | bc -l)

echo "total missing $TOTALMISSING"
echo "total HTZ $TOTALHTZ"
echo "total HMZ $TOTALHMZ"
echo "total $TOTAL,$FUSED,$FINAL from input $LINESBEFORE"

echo "total missing $TOTALMISSING" >> $FILTEROUTFOLDER/$SAMPLE.filtered.txt
echo "total HTZ $TOTALHTZ" >> $FILTEROUTFOLDER/$SAMPLE.filtered.txt
echo "total HMZ $TOTALHMZ" >> $FILTEROUTFOLDER/$SAMPLE.filtered.txt
echo "total $TOTAL from input $LINESBEFORE" >> $FILTEROUTFOLDER/$SAMPLE.filtered.txt




###############################################
echo "original: $LINESBEFORE ; total: $TOTAL ; fused: $FUSED ; final: $FINAL"
echo "original: $LINESBEFORE ; total: $TOTAL ; fused: $FUSED ; final: $FINAL" >> $FILTEROUTFOLDER/$SAMPLE.filtered.txt

MISSINGBEFORE=$(zcat $VCFINPUT | grep -v "#" | grep -c "\.:\.:\.:\.:\.") && echo $MISSINGBEFORE
MISSINGAFTER=$TOTALMISSING

MISSINGDIFF=$(echo "$MISSINGGENOTYPES-$MISSINGBEFORE" | bc -l)
echo "number of missing genotypes before/after/difference processing/filtering/setting Htz to missing: "$MISSINGDIFF

echo "processed $FILTEREDHTZ of $HTZ HTZ variants"


echo "$SAMPLE : $TOTAL : $FUSED : $FINAL : $MISSINGBEFORE : $MISSINGAFTER : $MISSINGDIFF : $TOTALHTZ : $TOTALHMZ : $HTZ : $TRUEHTZ : $FALSEHTZREF : $FALSEHTZALT : $HTZ2 : $HTZ3"\
>> $FILTEROUTFOLDER/stats.filtered.txt

#cleanup
cd $ORIGFOLDER
rm -rf $TF/

echo "############################################################################"
echo "##########       $SAMPLE done"
echo "############################################################################  \n \n \n"
### END OF GT-SCRIPT
