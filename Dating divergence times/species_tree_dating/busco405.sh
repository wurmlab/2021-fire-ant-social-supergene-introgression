#!/bin/bash
set -o pipefail

#conda activate busco_old_405
#INPUTFASTA="$PATH/name.fa"
# bash ~/scripts/busco405.sh -t $CPUs -i $INPUTFASTA

while getopts r:i:o:t:n: flag
do
    case "${flag}" in
        t) threads=${OPTARG};;
        i) input=${OPTARG};;
    esac
done

echo "threads: $threads";
echo "input / query fasta: $input";
echo "this script uses/requires busco4, seqtk"

## assign variables 
INPUTFASTA=$input
CPUs=$threads

CDIR="$PWD"

LINEAGE="/scratch/db/busco_odb10/hymenoptera_odb10"     ## modify
LINEAGENAME=$(echo $LINEAGE | rev | cut -d"/" -f1 | rev)
AUGUSTUSSPECIES="bombus_impatiens1"                             ## modifyecho "this run is uses this Database: $LINEAGE"

## get input's filename without folders and test if these exist

if [ -f "$INPUTFASTA" ]; then
    echo "INPUTFASTA exists"
else
    echo "INPUTFASTA does not exist: $INPUTFASTA" && exit 1
fi   

# check if input is gzipped
if [[ $INPUTFASTA == *.gz ]]
then 
 echo "gzipped input, temporary unpacking"
 pigz -dk $INPUTFASTA
 INPUTFILE=$(echo $INPUTFASTA | rev | cut -d"/" -f1 | cut -d"." -f2- | rev) && echo $INPUTFILE
 INPUTFILENAME=$(echo -n $INPUTFILE | rev | cut -d"." -f2- | rev ) && echo $INPUTFILENAME
 INPUTFOLDER=$(echo $INPUTFASTA | rev | cut -d"/" -f2- | rev) && echo $INPUTFOLDER
else
 echo "not gzipped input"
  INPUTFILE=$(echo $INPUTFASTA | rev | cut -d"/" -f1 | cut -d"." -f1- | rev) && echo $INPUTFILE
  INPUTFILENAME=$(echo -n $INPUTFILE | rev | cut -d"." -f2- | rev ) && echo $INPUTFILENAME
  INPUTFOLDER=$(echo $INPUTFASTA | rev | cut -d"/" -f2- | rev) && echo $INPUTFOLDER
fi

#INPUTFILE=$(echo $INPUTFASTA | rev | cut -d"/" -f1 | rev) && echo $INPUTFILE
#INPUTFILENAME=$(echo -n $INPUTFASTA | rev | cut -d"/" -f1 | cut -d"." -f2- | rev ) && echo $INPUTFILENAME
#INPUTFOLDER=$(echo $INPUTFASTA | rev | cut -d"/" -f2- | rev) && echo $INPUTFOLDER

cd $INPUTFOLDER
OUTPUTFOLDER="$INPUTFILE.busco405"
mkdir -p $OUTPUTFOLDER

## some basic stats
echo "$INPUTFILE" > $OUTPUTFOLDER/$INPUTFILE.stats
stats -s -t $INPUTFILE >> $OUTPUTFOLDER/$INPUTFILE.stats
echo


## making/checking .fai
if [ -f "$INPUTFILE.fai" ]; then
    echo "samtools index exists already, continuing..."
else
    echo "samtools index does not yet exist, creating..."
    samtools faidx $INPUTFILE
fi

## more basic stats
echo "10 longest scaffolds/contigs/chromosomes:" >> $OUTPUTFOLDER/$INPUTFILE.stats
cat $INPUTFILE.fai | cut -f1,2 | sort -k2,2nr | head >> $OUTPUTFOLDER/$INPUTFILE.stats
echo

echo "this run is used this starting species instead of Dmel: $AUGUSTUSSPECIES"
echo "this run is used this Database: $LINEAGE" >> $OUTPUTFOLDER/$INPUTFILE.stats
echo "this run is used this starting species instead of Dmel: $AUGUSTUSSPECIES" >> $OUTPUTFOLDER/$INPUTFILE.stats
echo $(busco -v) >> $OUTPUTFOLDER/$INPUTFILE.stats

#run BUSCO
busco -i $INPUTFOLDER/$INPUTFILE --force --mode genome --lineage_dataset $LINEAGE --offline --augustus_species $AUGUSTUSSPECIES --cpu $CPUs --evalue 1e-03 --limit 5 --out $OUTPUTFOLDER



cat $OUTPUTFOLDER/short_summary.specific.$LINEAGENAME.$INPUTFILE.busco.txt >> $OUTPUTFOLDER/$INPUTFILE.stats

cat $OUTPUTFOLDER/'run_'$LINEAGENAME/full_table.tsv > $OUTPUTFOLDER/$INPUTFILE.busco.full.table.tsv
cat $OUTPUTFOLDER/'run_'$LINEAGENAME/full_table.tsv | grep "Missing" > $OUTPUTFOLDER/$INPUTFILE.busco.missing.tsv

# extract the BUSCO gene sequences for later use (e.g. in phylogenies)
for GENE in $(ls -1 $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences/fragmented_busco_sequences/*.faa | rev | cut -d"/" -f1 | cut -d"." -f2- | rev)
do
cat $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences/fragmented_busco_sequences/$GENE.faa | seqtk seq -l 0 -C | sed -r "s,^>,>`echo -n $GENE`_`echo -n $INPUTFILENAME`_,g" >> $OUTPUTFOLDER/$INPUTFILE.busco.fragmented.faa
#.fna files only generated in busco4, not 5 anymore
cat $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences/fragmented_busco_sequences/$GENE.fna | seqtk seq -l 0 -C | sed -r "s,^>,>`echo -n $GENE`_`echo -n $INPUTFILENAME`_,g" >> $OUTPUTFOLDER/$INPUTFILE.busco.fragmented.fna
done

for GENE in $(ls -1 $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences/multi_copy_busco_sequences/*.faa | rev | cut -d"/" -f1 | cut -d"." -f2- | rev)
do
cat $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences/multi_copy_busco_sequences/$GENE.faa | seqtk seq -l 0 -C | sed -r "s,^>,>`echo -n $GENE`_`echo -n $INPUTFILENAME`_,g" >> $OUTPUTFOLDER/$INPUTFILE.busco.duplicated.faa
cat $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences/multi_copy_busco_sequences/$GENE.fna | seqtk seq -l 0 -C | sed -r "s,^>,>`echo -n $GENE`_`echo -n $INPUTFILENAME`_,g" >> $OUTPUTFOLDER/$INPUTFILE.busco.duplicated.fna
done

for GENE in $(ls -1 $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences/single_copy_busco_sequences/*.faa | rev | cut -d"/" -f1 | cut -d"." -f2- | rev)
do
cat $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences/single_copy_busco_sequences/$GENE.faa | seqtk seq -l 0 -C | sed -r "s,^>,>`echo -n $GENE`_`echo -n $INPUTFILENAME`_,g" >> $OUTPUTFOLDER/$INPUTFILE.busco.singlecopy.faa
cat $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences/single_copy_busco_sequences/$GENE.fna | seqtk seq -l 0 -C | sed -r "s,^>,>`echo -n $GENE`_`echo -n $INPUTFILENAME`_,g" >> $OUTPUTFOLDER/$INPUTFILE.busco.singlecopy.fna
done



#busco5
#tar -czf $OUTPUTFOLDER/$INPUTFILE.busco.fragmented_busco_sequences.tar.gz $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences/fragmented_busco_sequences/
#tar -czf $OUTPUTFOLDER/$INPUTFILE.busco.multi_copy_busco_sequences.tar.gz $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences/multi_copy_busco_sequences/
#tar -czf $OUTPUTFOLDER/$INPUTFILE.busco.single_copy_busco_sequences.tar.gz $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences/single_copy_busco_sequences/

#busco5
#cat $OUTPUTFOLDER/'run_'$LINEAGENAME/metaeuk_output/rerun_results/$INPUTFILE.codon.fas | seqtk seq -l 0 > $OUTPUTFOLDER/$INPUTFILE.busco.sequences.fa
#cat $OUTPUTFOLDER/'run_'$LINEAGENAME/metaeuk_output/rerun_results/$INPUTFILE.modified.fas | seqtk seq -l 0 > $OUTPUTFOLDER/$INPUTFILE.busco.sequences.faa


#rm -rf $OUTPUTFOLDER/'run_'$LINEAGENAME/hmmer_output
#rm -rf $OUTPUTFOLDER/'run_'$LINEAGENAME/busco_sequences
#tar -czf $OUTPUTFOLDER/'run_'$LINEAGENAME.tar.gz $OUTPUTFOLDER/'run_'$LINEAGENAME/
#rm -rf $OUTPUTFOLDER/'run_'$LINEAGENAME

#pigz -9 $OUTPUTFOLDER/$INPUTFILE.busco.sequences.faa
#pigz -9 $OUTPUTFOLDER/$INPUTFILE.busco.sequences.fa
#pigz -9 $OUTPUTFOLDER/$INPUTFILE.busco.duplicated.faa
#pigz -9 $OUTPUTFOLDER/$INPUTFILE.busco.singlecopy.faa
#pigz -9 $OUTPUTFOLDER/$INPUTFILE.busco.fragmented.faa

echo "Busco done:"
#cat $OUTPUTFOLDER/'run_'$LINEAGENAME/short_summary.txt
#### compress the BUSCO results to save space
#tar -czf $OUTPUTFOLDER/$INPUTFILE.busco.tar.gz $OUTPUTFOLDER/'run_'$LINEAGENAME

# show results
cat $OUTPUTFOLDER/$INPUTFILE.stats
#rm -rf $INPUTFOLDER/$INPUTFILE.busco
#rm -rf busco_downloads

cd $CDIR
