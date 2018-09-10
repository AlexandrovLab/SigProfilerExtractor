#!/bin/bash

cd ..

test_name=$1;
output_path="references/vcf_files/single/";
outputFile="references/vcf_files/single/"$test_name"_indels.genome";
vcf_path="references/vcf_files/$test_name/";
mkdir $output_path

rm -f $outputFile;

for fileName in $vcf_path*.vcf; do 
path=${fileName%/*}
xbase=${fileName##*/}
xfext=${fileName##*.}
xpref=$(basename "${fileName%.*}")

cat $fileName | egrep -v "#" | awk -F "\t" -v sample=$xpref -v name=$test_name 'BEGIN{OFS="\t"}{chrV=$1; gsub(/chr/,"",chrV); mutType = "INDEL"; refBeg = $2; redEnd = $2+length($4)-1; ref = $4; mut = $5;
print name,sample,"NCI","GRCh37",mutType,chrV,refBeg,redEnd,ref,mut,"SOMATIC"}' | egrep -v "," >> $outputFile
done
