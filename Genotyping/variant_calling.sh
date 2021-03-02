
# target regions for variant call, from the BUSCO4 results (full.table.tsv)
# each BUSCO gene coordinate was extended by 1kb upstream and downstream
# gng20170922wFex.fa.busco4.chrs.completeANDfragmented.extended1000bp.coordinates-fixed.merged.region


###############    SNPcall genes
ulimit -n 10240
INPUTFOLDER="bam_files"
OUTPUTFOLDER="$INPUTFOLDER/SNPIVbusco4all"
PROJECT="$OUTPUTFOLDER"
mkdir -p $OUTPUTFOLDER/vcf-master/raw
mkdir -p $OUTPUTFOLDER/vcf-master/filtered
mkdir -p $OUTPUTFOLDER/vcf-master/fused

ls -1 $INPUTFOLDER/bam/*.bam | rev | cut -d"." -f 2- | cut -d"/" -f 1 | rev > $OUTPUTFOLDER/allsamples.lst
ls -1 $BAMOUT/bam/*.bam | rev | cut -d"." -f 2- | cut -d"/" -f 1 | rev >> $OUTPUTFOLDER/allsamples.lst
ls -1 $INPUTFOLDER/bam/*.bam > $OUTPUTFOLDER/allvariantsamples.lst
ls -1 $BAMOUT/bam/*.bam >> $OUTPUTFOLDER/allvariantsamples.lst
BAMLST="$OUTPUTFOLDER/allvariantsamples.lst"
cat $BAMLST | wc -l  #388

BAMLST="$OUTPUTFOLDER/allvariantsamples.lst"
REF="$INPUTFOLDER/ref/gng20170922wFex.fa"
REGIONS="$REF.busco4.chrs.completeANDfragmented.extended1000bp.coordinates-fixed.merged.witMTGP9.region"
HICOVBED="$INPUTFOLDER/mosdepth.mergedbeds/high.cov.occurence.nr.samples.bed"
SSRBED="$REF.SSRs.ext5bp.bed"



#### call variants for Master SNP list
## min cov 4 because MANY SNPs have this number of ALT or REF reads only AND are the only individual with this Variant
## Fraction 0.4 to make call more stringent for detection of real variants in at least 1 sample. Other samples may have lower fraction when heterozygous - these will be called alongside if at least ONE other individual provides GOOD evidence that this SNP exists

CPUs=10
ulimit -n 10240
cd $INPUTFOLDER

MINALTFRAC=0.40
MINALTN=4
MINCOV=4

cat $REGIONS | parallel -k -j $CPUs "echo {} && freebayes --region {} \
--fasta-reference $REF \
--ploidy 2 \
--report-genotype-likelihood-max \
--use-mapping-quality \
--genotype-qualities \
--use-best-n-alleles 8 \
--haplotype-length 1 \
--min-mapping-quality 40 \
--min-base-quality 30 \
--min-alternate-fraction $MINALTFRAC \
--min-alternate-total $MINALTN \
--min-coverage $MINCOV \
--use-reference-allele \
--bam-list $BAMLST |\
bgzip -f -@ 4 -c /dev/stdin > $OUTPUTFOLDER/vcf-master/raw/{}.vcf.gz \
&& tabix -fp vcf $OUTPUTFOLDER/vcf-master/raw/{}.vcf.gz"


########################## filter each vcf roughly and decompose
######   filtering (QUAL > 30: to allow SNPs for which only 1 out of 200 indiv with 0 REF vs 4 ALT reads has a 1/1 genotype, while every other indiv has 0/0), decomposing, fixup, noindels, uniq

CPUs=30
cat $REGIONS | parallel --no-notice -j $CPUs \
"echo {} && echo filter; zcat $PROJECT/vcf-master/raw/{}.vcf.gz | vcffilter -f 'QUAL > 30' |\
vcfallelicprimitives --keep-info --keep-geno -t DECOMPOSED |\
vcffixup - | vcfstreamsort | vt normalize -n -r $REF -q - | vcfuniqalleles | vcfnoindels | vcfnumalt - |\
vcfnulldotslashdot | vcfstreamsort -a | vcfuniq |\
bgzip -f -@ 4 -c /dev/stdin > $PROJECT/vcf-master/filtered/{}.vcf.gz \
&& tabix -fp vcf $PROJECT/vcf-master/filtered/{}.vcf.gz"
echo "filtering done"

#merge per-scaffold-vcf's
zcat $PROJECT/vcf-master/filtered/$(cat $REGIONS | head -n1).vcf.gz | grep "#" > $PROJECT/vcf-master/fused/fused.vcf
cat $REGIONS | parallel -k -j 1 "echo {}; zcat $PROJECT/vcf-master/filtered/{}.vcf.gz | grep -v '#' >> $PROJECT/vcf-master/fused/fused.vcf "

#filter based on region
echo "filtering SSRs and Ns"
cat $PROJECT/vcf-master/fused/fused.vcf | awk '!(($4) == "N")' | vcfintersect -v -l -b $SSRBED | vcfintersect -v -l -b $HICOVBED | vcfstreamsort -a | vcfuniq | bgzip -f -@ 15 -c /dev/stdin > $PROJECT/vcf-master/fused/fused.sorted.vcf.gz
tabix -fp vcf $PROJECT/vcf-master/fused/fused.sorted.vcf.gz
echo "SSRs and HighCov eliminated"

# advanced stats filter #filter low coverage sites in biallelic sites when only represented by FWD or REV reads (SRF, SRR)
tabix -H $PROJECT/vcf-master/fused/fused.sorted.vcf.gz | grep "##" > $PROJECT/vcf-master/fused/fused.sorted.light.vcf
tabix -H $PROJECT/vcf-master/fused/fused.sorted.vcf.gz | grep -m1 "#CHROM" | cut -f1-9 >> $PROJECT/vcf-master/fused/fused.sorted.light.vcf
vcfkeepinfo $PROJECT/vcf-master/fused/fused.sorted.vcf.gz AB AC AF AN DP MEANALT NUMALT RO SRF SRR TYPE NS ODDS | grep -v "#" | cut -f1-10 >> $PROJECT/vcf-master/fused/fused.sorted.light.vcf
cat $PROJECT/vcf-master/fused/fused.sorted.light.vcf | bgzip -f -@ 15 -c /dev/stdin > $PROJECT/vcf-master/fused/fused.sorted.light.vcf.gz
tabix -fp vcf $PROJECT/vcf-master/fused/fused.sorted.light.vcf.gz

bcftools view -e "NUMALT=1 & ((INFO/SRF)<=($MINALTN/2) | (INFO/SRR)<=($MINALTN/2))" $PROJECT/vcf-master/fused/fused.sorted.light.vcf.gz | bgzip -f -@ 20 -c /dev/stdin > $PROJECT/vcf-master/fused/fused.sorted.light.filtered.vcf.gz
tabix -fp vcf $PROJECT/vcf-master/fused/fused.sorted.light.filtered.vcf.gz
#6519 sites eliminated eliminated
#3041720 sites kept

echo "vcf peek"
vt peek $PROJECT/vcf-master/fused/fused.sorted.light.filtered.vcf.gz 2> $PROJECT/vcf-master/fused/fused.sorted.light.filtered.vcf.stats
cat $PROJECT/vcf-master/fused/fused.sorted.light.filtered.vcf.stats
       no. of chromosomes                 :        772
       no. of SNP                         :    3041720
           2 alleles                      :         2990208 (4.34) [2430070/560138]
           3 alleles                      :           51174 (0.83) [46448/55900]
           4 alleles                      :             338 (0.50) [338/676



######### use MASTEr SNP list to genotype samples with slightly relxaed thresholds in appreciation of the low coverage of some samples

mkdir -p $OUTPUTFOLDER/vcf-genotyping/raw
mkdir -p $OUTPUTFOLDER/vcf-genotyping/filtered
mkdir -p $OUTPUTFOLDER/vcf-genotyping/fused

## decompress reference SNP list into RAMdisk for faster access
cp $OUTPUTFOLDER/vcf-master/fused/fused.sorted.light.filtered.vcf.{gz,gz.tbi} /dev/shm/

MINALTFRAC=0.35
MINALTN=2
MINCOV=2


CPUs=25
REFVCF="/dev/shm/fused.sorted.light.filtered.vcf.gz"
cat $REGIONS | parallel -k -j $CPUs "echo {} && freebayes --region {} \
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
bgzip -f -@ 4 -c /dev/stdin > $OUTPUTFOLDER/vcf-genotyping/raw/{}.vcf.gz \
&& tabix -fp vcf $OUTPUTFOLDER/vcf-genotyping/raw/{}.vcf.gz"


# filter
CPUs=30
cat $REGIONS | parallel --no-notice -k -j $CPUs \
"echo {} && $HOME/scripts/bcftools2.zsh $OUTPUTFOLDER/vcf-genotyping/raw/{}.vcf.gz $HICOVBED $SSRBED $OUTPUTFOLDER/vcf-genotyping/filtered $MINALTN $REF.capitalized.fa"

# merge-per-scf-VCFs
zcat $OUTPUTFOLDER/vcf-genotyping/filtered/$(cat $REGIONS | head -n1).vcf.gz | grep "#" > $OUTPUTFOLDER/vcf-genotyping/fused/fused.vcf
cat $REGIONS | parallel -k -j 1 "echo {}; zcat $OUTPUTFOLDER/vcf-genotyping/filtered/{}.vcf.gz | grep -v '#' >> $OUTPUTFOLDER/vcf-genotyping/fused/fused.vcf"
cat $OUTPUTFOLDER/vcf-genotyping/fused/fused.vcf | vcfstreamsort -a | bgzip -f -@ 20 -c /dev/stdin > $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.vcf.gz
tabix -fp vcf $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.vcf.gz
echo "merging done"

vt peek $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.gz 2> $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.stats
cat $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.stats

stats: no. of samples                     :        388
       no. of chromosomes                 :        772
       no. of SNP                         :    3044545
           2 alleles                      :         2993576 (4.36) [2434992/558584]
           3 alleles                      :           50635 (0.83) [45963/55307]
           4 alleles                      :             334 (0.50) [334/668]



### further filtering of ambiguous genotypes and heterozygous calls
# split into single sample VCF
mkdir -p $OUTPUTFOLDER/vcf-genotyping/samples
cat $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.samplenames | parallel -j 10 \
"echo {}; vcfkeepsamples $OUTPUTFOLDER/vcf-genotyping/fused/fused.sorted.renamedSamples.vcf.gz {} | bgzip -f -@ 5 -c /dev/stdin > $OUTPUTFOLDER/vcf-genotyping/samples/{}.vcf.gz && tabix -fp vcf $OUTPUTFOLDER/vcf-genotyping/samples/{}.vcf.gz"

### filter with filterscript (genotype filtering)

FILTEROUTFOLDER="$OUTPUTFOLDER/vcf-genotyping/samples-filtered"
mkdir -p $FILTEROUTFOLDER
ls -1 $OUTPUTFOLDER/vcf-genotyping/samples/*.vcf.gz | rev | cut -d "/" -f1 | cut -d"." -f3- | rev > $FILTEROUTFOLDER/vcf.lst
MINCOV=2
CPUs=15
cat $FILTEROUTFOLDER/vcf.lst | parallel --no-notice -j $CPUs "echo {} && genotype_filter.zsh $OUTPUTFOLDER/vcf-genotyping/samples/{}.vcf.gz $MINCOV {} $FILTEROUTFOLDER; echo {} && echo"

#check stats of missing sample (e.g. to ID samples of low coverage or diploidy):
cat $FILTEROUTFOLDER/stats.filtered.txt | tr -d ' ' | tr ':' '\t' | cut -f1,2 | cut -f2 | sort | uniq


### selection of filtered VCFs:
# pooled samples (n=3: gem-1-bigB|Lad4-1-bigB|Toc5-1-bigB): for each heterozygous site pick majority allele (the allele which has more reads within the respective sample)
# outgroup males (USP3-1-bigB-m|AR223-1-bigB-m|Copa2-1-bigB-m|Par1-1-bigB-m|AR112-1-bigB-p|AR142-1-bigB-m|SS1-2-bigB-m) which show elevated numbers of heterozygous sites: for each male and SNP one of the two alleles chosen for a heterozygous site.
# haploid samples: remove heterozygous genotypes (set to missing in the respective sample)
# remove low coverage samples entirely (i.e. if a sample had more than 25% variants missing, i.e. >=761k SNPs missing out of 3044545 total SNP: n=9: AR171-7-littleb|AR54-1-bigB|AR209-1-bigB|AR171-6-littleb|AR136-1-bigB|AR65-1-bigB|AR23-1-bigB|AR95-1-bigB|SRR7028253)
# diploid samples (diploid males, n=9: AR7-1-littleb-p|CGIn9-1-littleb-p|Mir8-6-littleb-p|Mir9-2-littleb-p|AR118-10-bigB-p|CaGr1A-1-bigB-m|U52-1-bigB-m|AdR11-2-bigB-p|CGIn1-1-littleb-p): remove sample entirely, evidenced by high levels of heterozygous sites. similar to pooled samples

cat $FILTEROUTFOLDER/merge-vcf.lst | wc -l  #370
# 18 samples removed (low cov 9x, diploids 9x)

MERGEFOLDER="$FILTEROUTFOLDER/merge"
mkdir -p $MERGEFOLDER
## merge 370 VCFs into file: $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.vcf
## here we did this by extracting the genotype columns per sample which (checked!) have identical SCF/POSITION/VARIANT entries.

cat $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.vcf | vcfuniq | vcffixup - | sed -r "s#\t\.:\.:\.:\.:\.:\.:\.:\.:\.#\t.#g" |\
bgzip -f -@ 35 -c /dev/stdin > $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.vcf.gz
tabix -fp vcf $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.vcf.gz


### get some stats

vt peek $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.vcf.gz 2> $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.vcf.stats
cat $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.vcf.stats

stats: no. of samples                     :        370
       no. of chromosomes                 :        772
       no. of SNP                         :    3044545
           2 alleles                      :         2993576 (4.36) [2434992/558584]
           3 alleles                      :           50635 (0.83) [45963/55307]
           4 alleles                      :             334 (0.50) [334/668]



########## remove sites for which more than 25% samples have a missing genotype
#max 25% missing, = 75% present)
NSAMPLES=$(vcfsamplenames $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.vcf.gz | wc -l) && echo $NSAMPLES
NMISSING25=$(echo -n $NSAMPLES | awk '{ $1=sprintf("%.0f",$1*0.25*2)} {print $1;}') && echo $NMISSING25
NMISSING=$((echo $NSAMPLES-$NMISSING25/2 | bc -l) | awk '{ $1=sprintf("%.0f",$1)} {print $1;}') && echo $NMISSING
zcat $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.vcf.gz | vcftools --recode --recode-INFO-all -c --vcf - --max-missing-count $NMISSING25 |\
bgzip -f -@ $CPUs -c /dev/stdin > $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.MaxMissing25.vcf.gz
tabix -fp vcf $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.MaxMissing25.vcf.gz

vt peek $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.MaxMissing25.vcf.gz 2> $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.MaxMissing25.vcf.stats
cat $MERGEFOLDER/busco4.merged.filtered.SNP.genotypes.MaxMissing25.vcf.stats

stats: no. of samples                     :        370
       no. of chromosomes                 :        772
       no. of SNP                         :    3039672
           2 alleles                      :         2988873 (4.36) [2431435/557438]
           3 alleles                      :           50469 (0.83) [45831/55107]
           4 alleles                      :             330 (0.50) [330/660]