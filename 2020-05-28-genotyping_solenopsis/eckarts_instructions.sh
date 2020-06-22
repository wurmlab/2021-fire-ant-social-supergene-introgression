
#######################

##Freebayes SNPcall



F="/data2/archive/archive-SBCS-WurmLab/db/genomic/reads/S_invicta/2020-05-bams_388_eckart"
cat $F/ref/gnG.REGIONS.splitby.FILTER2.500bp.bed | cut -f1-3 | sed "s/\.1\t/.1:/g" | sed "s/\.1A\t/.1A:/g" | sed "s/\.1B\t/.1B:/g" | sed "s/\.1C\t/.1C:/g" | sed 
"s/\.1D\t/.1D:/g" | tr '\t' '-' | sed "s/:0-/:1-/g" > $F/ref/gnG.REGIONS.splitby.FILTER2.500bp.regions
REGIONS="$F/ref/gnG.REGIONS.splitby.FILTER2.500bp.regions"   #14074 regions
REF="$F/ref/gng20170922wFex.fa"
OUTPUTFOLDER="$F/2020-05-12.SNPcall"
mkdir -p $OUTPUTFOLDER
mkdir -p $OUTPUTFOLDER/vcf-master/raw
mkdir -p $OUTPUTFOLDER/vcf-master/filtered
mkdir -p $OUTPUTFOLDER/vcf-master/fused
ls -1 $F/bam{1,2}/*.bam > $OUTPUTFOLDER/bam.lst
CPUs=20
ulimit -n 10240
FILTERBED="$F/ref/gnG.REGIONS.splitby.FILTER2.500bp.bed"
BAMLST="$OUTPUTFOLDER/bam.lst"
##################################
MINALTFRAC=0.40
MINALTN=4
MINCOV=4
#
#
cat $REGIONS | parallel -k -j $CPUs "echo {} && freebayes --dry-run --region {} \
--fasta-reference $REF \
--ploidy 2 \
--report-genotype-likelihood-max \
--use-mapping-quality \
--genotype-qualities \
--use-best-n-alleles 6 \
--haplotype-length 3 \
--min-mapping-quality 30 \
--min-base-quality 25 \
--min-alternate-fraction $MINALTFRAC \
--min-alternate-total $MINALTN \
--min-coverage $MINCOV \
--use-reference-allele \
--bam-list $BAMLST |\
bgzip -f -@ 1 -c /dev/stdin > $OUTPUTFOLDER/vcf-master/raw/{}.vcf.gz \
&& tabix -fp vcf $OUTPUTFOLDER/vcf-master/raw/{}.vcf.gz"

# filter
cat $REGIONS | parallel --dry-run --no-notice -k -j $CPUs \
"echo {} && zcat $OUTPUTFOLDER/vcf-master/raw/{}.vcf.gz | vcfintersect -v -l -b $FILTERBED | vcffilter -f 'QUAL > 30' |\
bcftools view -e "NUMALT=1 & ((INFO/SRF)<=($MINALTN/2) | (INFO/SRR)<=($MINALTN/2))" - |\
vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED |\
vcffixup - | vcfstreamsort | vt normalize -n -r $REF -q - | vcfuniqalleles | vcfnoindels | vcfnumalt - |\
vcfnulldotslashdot | vcfstreamsort | vcfuniq |\
bgzip -f -@ 1 -c /dev/stdin > $OUTPUTFOLDER/vcf-master/filtered/{}.vcf.gz && \
tabix -fp vcf $OUTPUTFOLDER/vcf-master/filtered/{}.vcf.gz"

### drop genotypes
cat $REGIONS | parallel --dry-run --no-notice -k -j $CPUs \
"echo {}; zcat $OUTPUTFOLDER/vcf-master/filtered/{}.vcf.gz | vcfkeepinfo - AB AC AF AN DP MEANALT NUMALT RO SRF SRR TYPE NS ODDS | grep -v '#' | cut -f1-9 >> 
$OUTPUTFOLDER/vcf-master/fused/{}.tmp.vcf"

#fuse  ## run directly with 20 Treads for compression
tabix -H $OUTPUTFOLDER/vcf-master/filtered/$(cat $REGIONS | head -n1).vcf.gz | grep "##" > $OUTPUTFOLDER/vcf-master/fused/fused.vcf
tabix -H $OUTPUTFOLDER/vcf-master/filtered/$(cat $REGIONS | head -n1).vcf.gz | grep -m1 "#CHROM" | cut -f1-9 >> $OUTPUTFOLDER/vcf-master/fused/fused.vcf
cat $OUTPUTFOLDER/vcf-master/fused/*.tmp.vcf >> $OUTPUTFOLDER/vcf-master/fused/fused.vcf

cat $OUTPUTFOLDER/vcf-master/fused/fused.vcf | vcfstreamsort -a | bgzip -f -@ 20 -c /dev/stdin > $OUTPUTFOLDER/vcf-master/fused/fused.vcf.gz
tabix -fp vcf $OUTPUTFOLDER/vcf-master/fused/fused.vcf.gz

rm -f $OUTPUTFOLDER/vcf-master/fused/fused.vcf
rm -f $OUTPUTFOLDER/vcf-master/fused/*.tmp.vcf

vt peek $OUTPUTFOLDER/vcf-master/fused/fused.vcf.gz 2> $OUTPUTFOLDER/vcf-master/fused/fused.vcf.stats
cat $OUTPUTFOLDER/vcf-master/fused/fused.vcf.stats

#####################################################################################################
## genotyping

## decompress reference SNP list into RAMdisk for faster access   #dunno if Apocrita allows this, then on scratch
## cp $OUTPUTFOLDER/vcf-master/fused/fused.vcf{gz,gz.tbi} /dev/shm/
##############################################
CPUs=20
ulimit -n 10240
FILTERBED="$F/ref/gnG.REGIONS.splitby.FILTER2.500bp.bed"
F="/data2/archive/archive-SBCS-WurmLab/db/genomic/reads/S_invicta/2020-05-bams_388_eckart"
REGIONS="$F/ref/gnG.REGIONS.splitby.FILTER2.500bp.regions"   #14074 regions
REF="$F/ref/gng20170922wFex.fa"
OUTPUTFOLDER="$F/2020-05-12.SNPcall"
BAMLST="$OUTPUTFOLDER/bam.lst"
mkdir -p $OUTPUTFOLDER/vcf-genotyping/raw
mkdir -p $OUTPUTFOLDER/vcf-genotyping/filtered
mkdir -p $OUTPUTFOLDER/vcf-genotyping/fused
##############################################
MINALTFRAC=0.35
MINALTN=2
MINCOV=2

CPUs=25
REFVCF="$OUTPUTFOLDER/vcf-master/fused/fused.vcf"

cat $REGIONS | parallel --dry-run -k -j $CPUs "echo {} && freebayes --region {} \
--fasta-reference $REF \
--ploidy 2 \
--haplotype-basis-alleles $REFVCF \
--report-genotype-likelihood-max \
--use-mapping-quality \
--genotype-qualities \
--use-best-n-alleles 4 \
--haplotype-length -1 \
--min-mapping-quality 30 \
--min-base-quality 28 \
--min-alternate-fraction $MINALTFRAC \
--min-alternate-total $MINALTN \
--min-coverage $MINCOV \
--use-reference-allele \
--bam-list $BAMLST |\
bgzip -f -@ 1 -c /dev/stdin > $OUTPUTFOLDER/vcf-genotyping/raw/{}.vcf.gz \
&& tabix -fp vcf $OUTPUTFOLDER/vcf-genotyping/raw/{}.vcf.gz"

# filter
cat $REGIONS | parallel --dry-run --no-notice -k -j $CPUs \
"echo {} && zcat $OUTPUTFOLDER/vcf-genotyping/raw/{}.vcf.gz | vcfintersect -v -l -b $FILTERBED | vcffilter -f 'QUAL > 30' |\
bcftools view -e "NUMALT=1 & ((INFO/SRF)<=($MINALTN/2) | (INFO/SRR)<=($MINALTN/2))" - |\
vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED |\
vcffixup - | vcfstreamsort | vt normalize -n -r $REF -q - | vcfuniqalleles | vcfnoindels | vcfnumalt - |\
vcfnulldotslashdot | vcfstreamsort | vcfuniq |\
bgzip -f -@ 1 -c /dev/stdin > $OUTPUTFOLDER/vcf-genotyping/filtered/{}.vcf.gz && \
tabix -fp vcf $OUTPUTFOLDER/vcf-genotyping/filtered/{}.vcf.gz"
echo "filter done"

# fuse
zcat $OUTPUTFOLDER/vcf-genotyping/filtered/$(cat $REGIONS | head -n1).vcf.gz | grep "#" > $OUTPUTFOLDER/vcf-genotyping/fused/fused.vcf
cat $REGIONS | parallel -k -j 1 "echo {}; zcat $OUTPUTFOLDER/vcf-genotyping/filtered/{}.vcf.gz | grep -v '#' >> $OUTPUTFOLDER/vcf-genotyping/fused/fused.vcf" 
## replaced by vcfstreamsort -a, roughly same speed # bcftools sort -Ov $OUTPUTFOLDER/vcf-genotyping/fused/fused.vcf -o 
$OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.vcf
cat $OUTPUTFOLDER/vcf-genotyping/fused/fused.vcf | vcfstreamsort -a | bgzip -f -@ 20 -c /dev/stdin > $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.vcf.gz
tabix -fp vcf $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.vcf.gz
rm -f $OUTPUTFOLDER/vcf-genotyping/fused/fused.vcf
echo "fusing done"

#### here we need to prep a new header (samplenames) and exchange, using bash or bcftools reheader
vcfsamplenames $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.vcf.gz > $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.samplenames.txt
cat $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.samplenames.txt  |\
cut -d"_" -f2,4 | sed -r "s,(_[0-9])\w+,,g" |\
sed "s,SRR9008101,SRR9008101_inv-153-For-bigB,g" |\
sed "s,SRR9008102,SRR9008102_inv-158-For-bigB,g" |\
sed "s,SRR9008107,SRR9008107_inv-231-BZ-bigB,g" |\
sed "s,SRR9008109,SRR9008109_inv-240-Cor-littleb,g" |\
sed "s,SRR9008114,SRR9008114_ric-113-Gua-littleb,g" |\
sed "s,SRR9008116,SRR9008116_ric-115-Gua-littleb,g" |\
sed "s,SRR9008120,SRR9008120_inv-182-For-bigB,g" |\
sed "s,SRR9008130,SRR9008130_inv-242-Cor-bigB,g" |\
sed "s,SRR9008131,SRR9008131_ric-117-Sal-littleb,g" |\
sed "s,SRR9008137,SRR9008137_ric-71-Ros-littleb,g" |\
sed "s,SRR9008139,SRR9008139_inv-220-Cor-bigB,g" |\
sed "s,SRR9008140,SRR9008140_inv-211-Cor-bigB,g" |\
sed "s,SRR9008141,SRR9008141_mac-152-Arr-bigB,g" |\
sed "s,SRR9008142,SRR9008142_mac-145-Col-bigB,g" |\
sed "s,SRR9008143,SRR9008143_inv-193-Cor-bigB,g" |\
sed "s,SRR9008145,SRR9008145_inv-222-Cor-bigB,g" |\
sed "s,SRR9008148,SRR9008148_ric-77-Ros-littleb,g" |\
sed "s,SRR9008152,SRR9008152_ric-67-Ros-bigB,g" |\
sed "s,SRR9008154,SRR9008154_ric-69-Ros-bigB,g" |\
sed "s,SRR9008156,SRR9008156_sae-44-Naz-bigB,g" |\
sed "s,SRR9008157,SRR9008157_sae-45-Kou-bigB,g" |\
sed "s,SRR9008159,SRR9008159_ric-118-Sal-bigB,g" |\
sed "s,SRR9008163,SRR9008163_inv-225-Mis-littleb,g" |\
sed "s,SRR9008166,SRR9008166_int-125-Car-littleb,g" |\
sed "s,SRR9008167,SRR9008167_meg-127-BZ-littleb,g" |\
sed "s,SRR9008170,SRR9008170_ric-107-Bol-littleb,g" |\
sed "s,SRR9008172,SRR9008172_inv-233-BZ-bigB,g" |\
sed "s,SRR9008173,SRR9008173_mac-146-Col-littleb,g" |\
sed "s,SRR9008174,SRR9008174_inv-196-Cor-littleb,g" |\
sed "s,SRR9008178,SRR9008178_ric-85-Ros-bigB,g" |\
sed "s,SRR9008188,SRR9008188_ric-116-Bol-littleb,g" |\
sed "s,SRR9008190,SRR9008190_ric-120-Bol-littleb,g" |\
sed "s,SRR9008192,SRR9008192_ric-106-Bol-littleb,g" |\
sed "s,SRR9008193,SRR9008193_ric-102-Bol-littleb,g" |\
sed "s,SRR9008195,SRR9008195_ric-103-Bol-littleb,g" |\
sed "s,SRR9008199,SRR9008199_ric-98-Ros-bigB,g" |\
sed "s,SRR9008203,SRR9008203_ric-82-Ros-bigB,g" |\
sed "s,SRR9008207,SRR9008207_ric-83-Ros-bigB,g" |\
sed "s,SRR9008218,SRR9008218_ric-99-Bol-bigB,g" |\
sed "s,SRR9008220,SRR9008220_inv-157-For-bigB,g" |\
sed "s,SRR9008226,SRR9008226_inv-219-Cor-littleb,g" |\
sed "s,SRR9008232,SRR9008232_xAdR-134-Bra-bigB,g" |\
sed "s,SRR9008238,SRR9008238_inv-185-For-bigB,g" |\
sed "s,SRR9008244,SRR9008244_mac-151-Per-bigB,g" |\
sed "s,SRR9008245,SRR9008245_mac-139-Gua-bigB,g" |\
sed "s,SRR9008250,SRR9008250_mac-137-Gua-bigB,g" |\
sed "s,SRR9008254,SRR9008254_inv-203-Cor-bigB,g" |\
sed "s,SRR9008258,SRR9008258_inv-194-LaP-bigB,g" |\
sed "s,SRR9008260,SRR9008260_inv-200-Cor-bigB,g" |\
sed "s,SRR9008262,SRR9008262_inv-195-Cor-bigB,g" |\
sed "s,SRR9008268,SRR9008268_inv-191-Con-bigB,g" |\
sed "s,SRR9008271,SRR9008271_inv-209-Cor-bigB,g" |\
sed "s,SRR9008272,SRR9008272_inv-207-Cor-bigB,g" |\
sed "s,SRR9008273,SRR9008273_inv-197-Cor-bigB,g" |\
sed "s,SRR9008275,SRR9008275_inv-223-Cor-bigB,g" |\
sed "s,SRR9008230,SRR9008230_inv-198-Cor-littleb,g" |\
sed "s,SRR9008105,SRR9008105_inv-235-BZ-bigB,g" |\
sed "s,SRR9008108,SRR9008108_inv-234-BZ-bigB,g" |\
sed "s,SRR9008110,SRR9008110_ric-66-Per-littleb,g" |\
sed "s,SRR9008119,SRR9008119_ric-79-Ros-littleb,g" |\
sed "s,SRR9008132,SRR9008132_ric-119-Sal-littleb,g" |\
sed "s,SRR9008136,SRR9008136_ric-74-Ros-littleb,g" |\
sed "s,SRR9008146,SRR9008146_inv-226-Cor-bigB,g" |\
sed "s,SRR9008150,SRR9008150_sae-46-Reg-bigB,g" |\
sed "s,SRR9008171,SRR9008171_ric-104-Bol-littleb,g" |\
sed "s,SRR9008240,SRR9008240_inv-239-Cor-bigB,g" |\
sed "s,SRR9008255,SRR9008255_inv-155-Ped-bigB,g" |\
sed "s,SRR9008264,SRR9008264_inv-162-For-bigB,g" |\
sed "s,SRR9008104,SRR9008104_inv-227-BZ-bigB,g" |\
sed "s,SRR9008111,SRR9008111_ric-111-Gua-littleb,g" |\
sed "s,SRR9008117,SRR9008117_ric-112-Gua-littleb,g" |\
sed "s,SRR9008144,SRR9008144_inv-161-For-bigB,g" |\
sed "s,SRR9008155,SRR9008155_ric-73-Ros-bigB,g" |\
sed "s,SRR9008165,SRR9008165_mac-143-Gua-littleb,g" |\
sed "s,SRR9008169,SRR9008169_inv-210-Cor-bigB,g" |\
sed "s,SRR9008176,SRR9008176_ric-110-Gua-bigB,g" |\
sed "s,SRR9008194,SRR9008194_ric-101-Bol-littleb,g" |\
sed "s,SRR9008200,SRR9008200_ric-88-Ros-bigB,g" |\
sed "s,SRR9008204,SRR9008204_ric-93-Ros-bigB,g" |\
sed "s,SRR9008210,SRR9008210_meg-132-BZ-bigB,g" |\
sed "s,SRR9008213,SRR9008213_meg-131-BZ-bigB,g" |\
sed "s,SRR9008221,SRR9008221_ric-108-Bol-littleb,g" |\
sed "s,SRR9008225,SRR9008225_inv-244-Cor-littleb,g" |\
sed "s,SRR9008228,SRR9008228_inv-199-Cor-littleb,g" |\
sed "s,SRR9008237,SRR9008237_inv-186-For-bigB,g" |\
sed "s,SRR9008241,SRR9008241_inv-245-Cor-bigB,g" |\
sed "s,SRR9008247,SRR9008247_mac-142-Gua-bigB,g" |\
sed "s,SRR9008249,SRR9008249_mac-144-Vil-bigB,g" |\
sed "s,SRR9008256,SRR9008256_inv-154-Pon-bigB,g" |\
sed "s,SRR9008266,SRR9008266_inv-221-Cor-bigB,g" |\
sed "s,SRR9008267,SRR9008267_inv-181-For-bigB,g" |\
sed "s,SRR9008270,SRR9008270_inv-213-Cor-bigB,g" |\
sed "s,SRR9008274,SRR9008274_inv-237-Cor-bigB,g" |\
sed "s,SRR9008098,SRR9008098_inv-228-BZ-bigB,g" |\
sed "s,SRR9008099,SRR9008099_inv-229-BZ-bigB,g" |\
sed "s,SRR9008100,SRR9008100_inv-156-For-bigB,g" |\
sed "s,SRR9008106,SRR9008106_inv-236-BZ-bigB,g" |\
sed "s,SRR9008112,SRR9008112_ric-75-Ros-littleb,g" |\
sed "s,SRR9008113,SRR9008113_ric-65-Per-littleb,g" |\
sed "s,SRR9008151,SRR9008151_inv-238-Cor-littleb,g" |\
sed "s,SRR9008161,SRR9008161_inv-205-Ver-littleb,g" |\
sed "s,SRR9008223,SRR9008223_inv-188-For-littleb,g" |\
sed "s,SRR9008231,SRR9008231_inv-215-Cor-littleb,g" |\
sed "s,SRR9008236,SRR9008236_inv-165-For-bigB,g" |\
sed "s,SRR9008263,SRR9008263_inv-180-For-bigB,g" |\
sed "s,SRR9008118,SRR9008118_int-124-Car-littleb,g" |\
sed "s,SRR9008180,SRR9008180_ric-89-Ros-bigB,g" |\
sed "s,SRR9008181,SRR9008181_ric-90-Ros-bigB,g" |\
sed "s,SRR9008183,SRR9008183_mac-150-Per-littleb,g" |\
sed "s,SRR9008186,SRR9008186_mac-140-Gua-littleb,g" |\
sed "s,SRR9008198,SRR9008198_ric-97-Ros-bigB,g" |\
sed "s,SRR9008208,SRR9008208_ric-84-Ros-bigB,g" |\
sed "s,SRR9008212,SRR9008212_meg-129-BZ-bigB,g" |\
sed "s,SRR9008227,SRR9008227_inv-241-Cor-littleb,g" |\
sed "s,SRR9008261,SRR9008261_inv-230-Sao-bigB,g" |\
sed "s,SRR9008129,SRR9008129_inv-212-Cor-bigB,g" |\
sed "s,SRR9008133,SRR9008133_ric-91-Ros-littleb,g" |\
sed "s,SRR9008135,SRR9008135_ric-76-Ros-littleb,g" |\
sed "s,SRR9008138,SRR9008138_ric-72-Ros-littleb,g" |\
sed "s,SRR9008147,SRR9008147_inv-201-Cor-littleb,g" |\
sed "s,SRR9008158,SRR9008158_sae-47-Bel-bigB,g" |\
sed "s,SRR9008168,SRR9008168_xAdR-135-Pos-littleb,g" |\
sed "s,SRR9008175,SRR9008175_ric-86-Ros-bigB,g" |\
sed "s,SRR9008177,SRR9008177_ric-80-Ros-bigB,g" |\
sed "s,SRR9008187,SRR9008187_ric-109-Bol-littleb,g" |\
sed "s,SRR9008189,SRR9008189_ric-105-Bol-littleb,g" |\
sed "s,SRR9008191,SRR9008191_ric-100-Bol-littleb,g" |\
sed "s,SRR9008205,SRR9008205_inv-202-Cor-bigB,g" |\
sed "s,SRR9008206,SRR9008206_inv-208-Cor-bigB,g" |\
sed "s,SRR9008211,SRR9008211_meg-126-BZ-bigB,g" |\
sed "s,SRR9008215,SRR9008215_meg-128-BZ-bigB,g" |\
sed "s,SRR9008222,SRR9008222_inv-189-For-littleb,g" |\
sed "s,SRR9008224,SRR9008224_inv-217-Cor-littleb,g" |\
sed "s,SRR9008239,SRR9008239_inv-187-For-bigB,g" |\
sed "s,SRR9008242,SRR9008242_inv-224-Bra-littleb,g" |\
sed "s,SRR9008253,SRR9008253_inv-218-Cor-bigB,g" |\
sed "s,SRR9008257,SRR9008257_mac-147-Col-littleb,g" |\
sed "s,SRR9008164,SRR9008164_inv-160-For-littleb,g" |\
sed "s,SRR9008153,SRR9008153_ric-68-Ros-bigB,g" |\
sed "s,SRR9008182,SRR9008182_ric-87-Ros-bigB,g" |\
sed "s,SRR9008196,SRR9008196_inv-164-For-bigB,g" |\
sed "s,SRR9008202,SRR9008202_ric-92-Bol-bigB,g" |\
sed "s,SRR9008229,SRR9008229_inv-183-For-littleb,g" |\
sed "s,SRR9008234,SRR9008234_inv-206-Ver-bigB,g" |\
sed "s,SRR9008235,SRR9008235_inv-184-For-bigB,g" |\
sed "s,SRR9008248,SRR9008248_mac-138-Rin-bigB,g" |\
sed "s,SRR9008251,SRR9008251_mac-148-Gua-bigB,g" |\
sed "s,SRR9008252,SRR9008252_mac-149-Gua-bigB,g" |\
sed "s,SRR9008269,SRR9008269_inv-214-Cor-bigB,g" |\
sed "s,SRR9008103,SRR9008103_inv-232-BZ-bigB,g" |\
sed "s,SRR9008115,SRR9008115_ric-114-Gua-littleb,g" |\
sed "s,SRR9008134,SRR9008134_ric-81-Ros-littleb,g" |\
sed "s,SRR9008160,SRR9008160_inv-163-For-bigB,g" |\
sed "s,SRR9008185,SRR9008185_ric-70-Ros-bigB,g" |\
sed "s,SRR9008162,SRR9008162_inv-216-Pay-littleb,g" |\
sed "s,SRR9008179,SRR9008179_ric-95-Ros-bigB,g" |\
sed "s,SRR9008184,SRR9008184_ric-78-Ros-bigB,g" |\
sed "s,SRR9008197,SRR9008197_ric-96-Ros-bigB,g" |\
sed "s,SRR9008201,SRR9008201_ric-94-Ros-bigB,g" |\
sed "s,SRR9008209,SRR9008209_xAdR-136-Pos-bigB,g" |\
sed "s,SRR9008214,SRR9008214_meg-130-BZ-bigB,g" |\
sed "s,SRR9008216,SRR9008216_int-123-Car-bigB,g" |\
sed "s,SRR9008217,SRR9008217_int-122-Car-bigB,g" |\
sed "s,SRR9008219,SRR9008219_inv-159-For-bigB,g" |\
sed "s,SRR9008233,SRR9008233_inv-190-Vil-bigB,g" |\
sed "s,SRR9008243,SRR9008243_inv-192-Con-bigB,g" |\
sed "s,SRR9008246,SRR9008246_mac-141-Gua-bigB,g" |\
sed "s,SRR9008259,SRR9008259_inv-204-Cor-bigB,g" |\
sed "s,SRR9008265,SRR9008265_inv-243-Cor-bigB,g" |\
sed "s,SRR7028245,SRR7028245_AL-151-bigB,g" |\
sed "s,SRR7028246,SRR7028246_AL-150-bigB,g" |\
sed "s,SRR7028247,SRR7028247_AL-154-bigB,g" |\
sed "s,SRR7028248,SRR7028248_AL-153-bigB,g" |\
sed "s,SRR7028249,SRR7028249_AL-145-bigB,g" |\
sed "s,SRR7028250,SRR7028250_AL-141-bigB,g" |\
sed "s,SRR7028251,SRR7028251_AL-149-bigB,g" |\
sed "s,SRR7028252,SRR7028252_AL-148-bigB,g" |\
sed "s,SRR7028253,SRR7028253_HE-117-bigB,g" |\
sed "s,SRR7028254,SRR7028254_HE-116-bigB,g" |\
sed "s,SRR7028255,SRR7028255_HE-112-bigB,g" |\
sed "s,SRR7028256,SRR7028256_HE-109-bigB,g" |\
sed "s,SRR7028257,SRR7028257_AL-158-bigB,g" |\
sed "s,SRR7028258,SRR7028258_AL-156-bigB,g" |\
sed "s,SRR7028259,SRR7028259_HE-103-bigB,g" |\
sed "s,SRR7028260,SRR7028260_HE-99-bigB,g" |\
sed "s,SRR7028261,SRR7028261_AL-139-bigB,g" |\
sed "s,SRR7028262,SRR7028262_AL-134-bigB,g" |\
sed "s,SRR7028263,SRR7028263_HE-120-bigB,g" |\
sed "s,SRR7028264,SRR7028264_AL-132-bigB,g" |\
sed "s,SRR7028265,SRR7028265_AL-131-bigB,g" |\
sed "s,SRR7028266,SRR7028266_AL-130-bigB,g" |\
sed "s,SRR7028267,SRR7028267_AL-129-bigB,g" |\
sed "s,SRR7028268,SRR7028268_AL-127-bigB,g" |\
sed "s,SRR7028269,SRR7028269_AL-126-bigB,g" |\
sed "s,SRR7028270,SRR7028270_AL-124-bigB,g" |\
sed "s,SRR7028271,SRR7028271_HE-118-bigB,g" |\
sed "s,SRR7028272,SRR7028272_AL-133-bigB,g" |\
sed "s,SRR7028273,SRR7028273_HE-106-bigB,g" |\
sed "s,SRR7028274,SRR7028274_HE-94-bigB,g" |\
sed "s,SRR7028275,SRR7028275_HE-98-bigB,g" |\
sed "s,SRR7028276,SRR7028276_HE-104-bigB,g" |\
sed "s,SRR7028277,SRR7028277_HE-85-bigB,g" |\
sed "s,SRR7028278,SRR7028278_HE-86-bigB,g" |\
sed "s,SRR7028279,SRR7028279_HE-88-bigB,g" |\
sed "s,SRR7028280,SRR7028280_HE-89-bigB,g" |\
sed "s,SRR7028281,SRR7028281_HE-90-bigB,g" |\
sed "s,SRR7028282,SRR7028282_HE-91-bigB,g" |\
sed "s,SRR7028283,SRR7028283_HE-92-bigB,g" |\
sed "s,SRR7028284,SRR7028284_HE-93-bigB,g" |\
sed "s,-bigB-bigB,-bigB,g" |\
sed "s,_Privman2018,,g" |\
sed "s,_na-Wurm2011,,g" |\
sed "s,_na-Wang2013,,g" > $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.newsamplenames.txt

bcftools reheader --samples $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.newsamplenames.txt $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.vcf.gz | bgzip -f -@ 
20 -c /dev/stdin > $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.gz
tabix -fp vcf $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.gz

vcfsamplenames $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.gz > $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.samplenames
diff $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.vcf.samplenames $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.samplenames

vt peek $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.gz 2> $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.stats
cat $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.stats

## split into single sample VCF
mkdir -p $OUTPUTFOLDER/vcf-genotyping/samples
cat $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.samplenames | parallel -j 10 \
"echo {}; vcfkeepsamples $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.gz {} | bgzip -f -@ 5 -c /dev/stdin > 
$OUTPUTFOLDER/vcf-genotyping/samples/{}.vcf.gz && tabix -fp vcf $OUTPUTFOLDER/vcf-genotyping/samples/{}.vcf.gz"
echo "sample vcf done"

##### transfer and I do some more filtering .... or I setup whats needed to filter it there

