#!/bin/zsh
############################################################################
# filter SNPs
# Eckart Stolle, May.2020
############################################################################
## script: bcftools2.zsh
## filter 1 VCF
## $HOME/scripts/bcftools2.zsh $OUTPUTFOLDER/vcf-genotyping/raw/{}.vcf.gz $HICOVBED $SSRBED $OUTPUTFOLDER/vcf-genotyping/filtered MINALTN $REF.capitalized.fa

if [ $# -ne 6 ]; then
    echo $0: usage: ./bcftools2.zsh INPUT.raw.vcf.gz FilterOut1.bed FilterOut2.bed OUTPUTFOLDER 
    exit 1
fi

#set/get variables
VCFfile=$1
VCFfileName=$(echo $VCFfile | rev | cut -d"/" -f1 | cut -d"." -f3- | rev)
FILTER1=$2
FILTER2=$3
OUPUTFOLDER=$4
MINALTN=$5
REF=$6

zcat $VCFfile | vcfintersect -v -l -b $FILTER1 | vcfintersect -v -l -b $FILTER2 | vcffilter -f 'QUAL > 30' |\
bcftools view -e "NUMALT=1 & ((INFO/SRF)<=($MINALTN/2) | (INFO/SRR)<=($MINALTN/2))" - |\
vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED |\
vcffixup - | vcfstreamsort | vt normalize -n -r $REF -q - | vcfuniqalleles | vcfnoindels | vcfnumalt - |\
vcfnulldotslashdot | vcfstreamsort | vcfuniq |\
bgzip -f -@ 1 -c /dev/stdin > $OUPUTFOLDER/$VCFfileName.vcf.gz

tabix -fp vcf $OUPUTFOLDER/$VCFfileName.vcf.gz
echo "$VCFfileName filter done"
##################################
