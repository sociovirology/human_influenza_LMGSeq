#!/bin/bash

#Demultiplexing Script

############################### Below are commands to automate and abstract script, but for now hardcoded ##########
#Move into directory
#cd ~/GbBSeq_human_IAV/experiments/demultiplexing/2019_11_13

#Move barcodes file to analysis directory, but only copy the barcodes used in the current run

#First barcodes5.fasta
#Get column 1 from the cross_data_runD file, remove the first line, add cross to the beginning of the number so we can do an exact search with grep, use those strings as input to search barcodes file and output only the crosses in the run
#cut -d "," -f 1 ~/GbBSeq_human_IAV/data/cross_data_runD.csv | tail -n +2 | awk '$0="cross"$0' | grep -f - ~/GbBSeq_human_IAV/shared/barcodes5.fasta -A1 -w --no-group-separator > barcodes5.fasta 

#Now barcodes5rev.fasta
#cut -d "," -f 1 ~/GbBSeq_human_IAV/data/cross_data_runD.csv | tail -n +2 | awk '$0="cross"$0' | grep -f - ~/GbBSeq_human_IAV/shared/barcodes5rev.fasta -A1 -w --no-group-separator > barcodes5rev.fasta 

#Copy 3' barcodes. Always need all 96 of these, so copying as is
#cp ~/GbBSeq_human_IAV/shared/barcodes3.fasta .
#cp ~/GbBSeq_human_IAV/shared/barcodes3rev.fasta .

####################################################################################

#Demultiplexing Script for Influenza GbBSeq

#This script conducts demultiplexing for one run, which is indicated in the command-line
#For now, harcoding in the script with runD to verify it all works and get stuff done

#Change into working directory
cd ~/GbBSeq_human_IAV

#Make a directory to hold all the intermediate files etc.
mkdir runD #will reach and replace runD with the variable 

## First Copy Barcodes only for relevant crosses 
#Only need cross15 parentals plate, select from shared barcodes file 
grep cross15 ~/GbBSeq_human_IAV/shared/barcodes5.fasta -A1 -w --no-group-separator > runD/barcodes5.fasta
#Repeat with 
grep cross15 ~/GbBSeq_human_IAV/shared/barcodes5rev.fasta -A1 -w --no-group-separator > runD/barcodes5rev.fasta

#Now, similarly just need a few samples from 3' barcodes 

#This piped command series does the following: 
# creates a sequence of numbers in each line |
# awk takes numbers from pipe and adds sample |
# sample/number's are fed to grep which picks these samples from the barcode file  >
# barcodes output to FASTA in current directory 

for i in {7..12}  {19..24} {43..48} {55..60}; do echo $i; done | awk '$0="sample"$0' | grep -f - ~/GbBSeq_human_IAV/shared/barcodes3.fasta -A1 -w --no-group-separator > barcodes3.fasta 

#Same for reverse 3 (sample) barcodes
for i in {7..12}  {19..24} {43..48} {55..60}; do echo $i; done | awk '$0="sample"$0' | grep -f - ~/GbBSeq_human_IAV/shared/barcodes3rev.fasta -A1 -w --no-group-separator > barcodes3rev.fasta


#Change into Run directory
cd runD

#First, I'll note that the barcode files have been generated.
#Now will demultiplex cross and samples.
cutadapt -g file:barcodes5.fasta -G file:barcodes3rev.fasta  -O 8 --no-indels --action=none --discard-untrimmed -o trimmed_{name1}_{name2}_runD_R1.fastq -p trimmed_{name1}_{name2}_runD_R2.fastq ~/GbBSeq_human_IAV/data/runD_R1.fastq ~/GbBSeq_human_IAV/data/runD_R2.fastq > cutadapt_demultiplexing_report.txt

#Checking on the success of this step
head -n 25 cutadapt_demultiplexing_report.txt

#Now let's trim those barcodes.

for file in trimmed_cross*_sample*_runD_R1.fastq;
do
    prefix=${file%_runD_R1.fastq} #Get file prefix
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
    cutadapt -a ${barcode5}...${barcode3} -A ${barcode3rev}...${barcode5rev} -O 8 -o ${cross}_${sample}_R1.fastq -p ${cross}_${sample}_R2.fastq trimmed_${cross}_${sample}_runD_R1.fastq trimmed_${cross}_${sample}_runD_R2.fastq >> cutadapt_amplicon_trimming.txt
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
  
  cutadapt -g file:/home/sldmunoz/GbBSeq_human_IAV/shared/segment_specific_primers.fasta -A file:/home/sldmunoz/GbBSeq_human_IAV/shared/segment_specific_primers_rev.fasta -e 0.0 -O 10 --no-indels -o ${cross_sample}_{name}_R1.fastq -p ${cross_sample}_{name}_R2.fastq ${prefix}_R1.fastq ${prefix}_R2.fastq > cutadapt_${cross_sample}_primer_demultiplexing.txt

done

#Check
grep "Total written (filtered):" *_primer_demultiplexing.txt
#Most all in the 90's

echo "Demultiplexing done!"