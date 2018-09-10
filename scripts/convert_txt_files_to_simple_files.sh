#!/bin/bash

cd ..

test_name=$1;
output_path="references/vcf_files/single/";
outputFile="references/vcf_files/single/"$test_name"_indels.genome";
vcf_path="references/vcf_files/$test_name/";
#mkdir $output_path

rm -f $outputFile;

for fileName in $vcf_path*.txt; do  
#for fileName in $vcf_path*.genome; do 
path=${fileName%/*}
xbase=${fileName##*/}
xfext=${fileName##*.}

cat $fileName | egrep -v "#" | awk -F "\t" -v name=$test_name 'BEGIN{OFS="\t"}{chrV=$6; gsub(/chr/,"",chrV); sample = $2; mutType = $5; refBeg = $7; redEnd = $8; ref = $9; mut = $10;
print name,sample,"NCI","GRCh37",mutType,chrV,refBeg,redEnd,ref,mut,"SOMATIC"}' | egrep -v "," >> $outputFile
done
