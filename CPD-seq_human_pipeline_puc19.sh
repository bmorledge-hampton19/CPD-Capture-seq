#!/bin/bash

#Code adapted from code on Taylor lab github site

export NA=$1
mkdir ${NA}
cd ${NA}
ln -s ../${NA}.fastq

# initial alignment and processing of files
OUTPUT="$(bowtie2 -x /Users/johnwyrick/bin/bowtie2-2.2.6/hg19_puc19 -U ${NA}.fastq -S ${NA}.sam)"
echo "${OUTPUT}"
samtools view -b -S ${NA}.sam >${NA}.bam
bamToBed -i ${NA}.bam >${NA}.bed

# get damage coordinate and sequence of damage
perl ../finddamagecoord.pl <${NA}.bed >${NA}_puc19_damage.bed

# filter puc19 sequences
perl ../filter_pUC19.pl <${NA}_puc19_damage.bed >${NA}_damage.bed

# remove previous genome index file
rm ../hg19_puc19.fa.fai

# get nucleotides of damage sites
fastaFromBed -s -name -fi ../hg19_puc19.fa -bed ${NA}_damage.bed -fo ${NA}_damage.fa

# count damage nucleotides
perl ../count_dinucleotides.pl <${NA}_damage.fa >${NA}_nucount.txt

# extract dipy reads
printf "${NA}_damage.bed\n${NA}_damage.fa\n" | perl ../extract_dipyrimidine_reads.pl >${NA}_dipy.bed

# sort dipy reads
cat ${NA}_dipy.bed | sort -k1,1 -k2,2n -k 6 >${NA}_dipy_sorted.bed

# split strands
printf "${NA}_dipy_sorted.bed\n" | perl ../split_strands.pl 

# rest of processing should be done mannually --> count function of IGV tools, withi window size 1, zoom 10, to generate .wig files
# then run perl set_background.pl >[bkgd set wig file]

# extract dipy reads
printf "${NA}_damage.bed\n${NA}_damage.fa\n" | perl ../extract_Cdipy_reads.pl >${NA}_Cdipy.bed

# sort dipy reads
cat ${NA}_Cdipy.bed | sort -k1,1 -k2,2n -k 6 >${NA}_Cdipy_sorted.bed

# split strands
printf "${NA}_Cdipy_sorted.bed\n" | perl ../split_strands.pl 

# rest of processing should be done manually --> count function of IGV tools, withi window size 1, zoom 10, to generate .wig files
# then run perl set_background.pl >[bkgd set wig file]
