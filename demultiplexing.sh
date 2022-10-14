#!/bin/bash

# Demultiplexing Script for Influenza GbBSeq
# Script: demultiplexing.sh
# This file is shell script that downloads and processes raw sequencing reads
# Requires three arguments that should be passed from the command line:
#    run information, cross information, and experiment folders

#This script is part of the following manuscript:
#Influenza A virus reassortment is strain dependent
#Kishana Y. Taylor | Ilechukwu Agu  | Ivy José | Sari Mäntynen | A.J. Campbell | Courtney Mattson |  Tsui-wen Chou | Bin Zhou | David Gresham | Elodie Ghedin |  Samuel L. Díaz Muñoz

#Verify we have 3 arguments

if [ "$#" -ne 3 ]; then
    echo "You must enter exactly 3 command line arguments: run, cross, and experiment"
    echo "run has format 'runA'"
    echo "File name for cross list. Cross list file is just a plain text file one line per cross: cross1 [break] cross9 etc. "
    echo "Experiment can be any folder name you desire, e.g. 'infection_conditions'"
    exit
fi

#Checking command line parameters:
echo "Checking command line arguments: run, cross, and experiment"
echo $1
echo $2
echo $3

#Download FASTQ sequencing files from SRA (Illumnia PE read files)
#wget https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/021880/SRR21880043/GBSPCR2b_S1_L001_R1_001.fastq.gz -O data/runA_R1.fastq.gz
#wget https://sra-download.be-md.ncbi.nlm.nih.gov/vast/sra01/SRZ/021880/SRR21880043/GBSPCR2b_S1_L001_R2_001.fastq.gz -O data/runA_R2.fastq.gz

#Unzip files
#gunzip data/runA_*.fastq.gz

#This script conducts demultiplexing for one run, which is indicated in the first command-line argument: runA,
# the second command line argument includes a cross list (in this case detailing the experimental coinfections),
# and the third command line argument groups results into a directory (in this case we call this script twice to separate control experiments and pairwise experimental coinfections)

#Change into working directory, which in this case is the base directory where the script resides
BASEDIR=$(dirname "$0")
cd $BASEDIR

#Make a directory to hold all the intermediate files etc.
INTERMED_DIR=mkdir "$1""_""$3"
mkdir $INTERMED_DIR

## First Copy Barcodes only for relevant crosses 
cat $2 | grep -f - $BASEDIR/shared/barcodes5.fasta -A1 -w --no-group-separator > $INTERMED_DIR/barcodes5.fasta
#Repeat with reverse barcodes
cat $2 | grep -f - $BASEDIR/shared/barcodes5rev.fasta -A1 -w --no-group-separator > $INTERMED_DIR/barcodes5rev.fasta

#Now, copy entire files for 3' barcodes because we want all samples 
cp $BASEDIR/shared/barcodes3.fasta $INTERMED_DIR/barcodes3.fasta
cp $BASEDIR/shared/barcodes3rev.fasta $INTERMED_DIR/barcodes3rev.fasta

#START HERE OCT 13, 2022: Need to:
#1. Check needed files in shared/
#2. Check output directory, if any?? Because output comes from amplicon curation, right?

#Change into Run directory
cd $INTERMED_DIR

#First, I'll note that the barcode files have been generated.
#Now will demultiplex cross and samples.
cutadapt -g file:barcodes5.fasta -G file:barcodes3rev.fasta -j 0 -O 8 --no-indels --action=none --discard-untrimmed -o trimmed_{name1}_{name2}_$1_R1.fastq -p trimmed_{name1}_{name2}_$1_R2.fastq ~/GbBSeq_human_IAV/data/$1_R1.fastq ~/GbBSeq_human_IAV/data/$1_R2.fastq > cutadapt_demultiplexing_report.txt

#Checking on the success of this step
head -n 25 cutadapt_demultiplexing_report.txt

#Now let's trim those barcodes.

for file in trimmed_cross*_sample*_$1_R1.fastq;
do
    prefix=${file%_$1_R1.fastq} #Get file prefix
    cross_sample=${prefix#trimmed_} #Get cross sample combination
    cross=${cross_sample%_sample*}
    sample=${cross_sample#cross*_}
    echo "$file"
    echo "$prefix"
    echo "$cross_sample"
    echo "$cross"
    echo "$sample"

    barcode5=$(grep ${cross} -A1 -w --no-group-separator barcodes5.fasta | grep -v ${cross})
    barcode5rev=$(grep ${cross} -A1 -w --no-group-separator barcodes5rev.fasta | grep -v ${cross})
    barcode3=$(grep ${sample} -A1 -w --no-group-separator barcodes3.fasta | grep -v ${sample})
    barcode3rev=$(grep ${sample} -A1 -w --no-group-separator barcodes3rev.fasta | grep -v ${sample})

#Trimming barcodes noted above from entire run sequencingfile
    cutadapt -a ${barcode5}...${barcode3} -A ${barcode3rev}...${barcode5rev} -O 8 -o ${cross}_${sample}_R1.fastq -p ${cross}_${sample}_R2.fastq trimmed_${cross}_${sample}_$1_R1.fastq trimmed_${cross}_${sample}_$1_R2.fastq >> cutadapt_amplicon_trimming.txt
done

#Check on status
grep "Total written (filtered):" cutadapt_amplicon_trimming.txt
#High seventies to 90's for the most part

#Now we will trim the uni13 primer on all samples
for file in cross*_sample*_R1.fastq;
do
  prefix=${file%_R1.fastq} #Get file prefix
  cutadapt -a "CCTTGTTTCTACT" -G "^AGTAGAAACAAGG" -o uni13_${prefix}_R1.fastq -p uni13_${prefix}_R2.fastq ${prefix}_R1.fastq ${prefix}_R2.fastq  > cutadapt_${prefix}_uni13_trimming.txt
done

grep "Total written (filtered):" *uni13_trimming.txt
#All in the 90's pretty much


#Now we move to primer demultiplexing by locus:

for file in uni13_cross*_sample*_R1.fastq;
do
  prefix=${file%_R1.fastq} #Get file prefix
  cross_sample=${prefix#uni13_} #Get cross sample combination
  
  cutadapt -g file:$BASEDIR/shared/segment_specific_primers.fasta -A file:$BASEDIR/shared/segment_specific_primers_rev.fasta -e 0.0 -O 10 --no-indels -o ${cross_sample}_{name}_R1.fastq -p ${cross_sample}_{name}_R2.fastq ${prefix}_R1.fastq ${prefix}_R2.fastq > cutadapt_${cross_sample}_primer_demultiplexing.txt

done

#Check
grep "Total written (filtered):" *_primer_demultiplexing.txt
#Most all in the 90's

echo "Demultiplexing done!"