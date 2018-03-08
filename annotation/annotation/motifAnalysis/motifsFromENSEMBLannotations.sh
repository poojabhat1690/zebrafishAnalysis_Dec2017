#!/bin/bash


ml bedtools/2.27.1-foss-2017a
ml homer

grep -e "\+" /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/ensembl_dr10_Ensembl_Genes_88/proteinCoding_annotatedUTRs.bed > /scratch/pooja/proteinCoding_annotatedUTRs_plus.bed
cat /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/ensembl_dr10_Ensembl_Genes_88/proteinCoding_annotatedUTRs.bed  | awk '$6 == "-" { print $0 }' > /scratch/pooja/proteinCoding_annotatedUTRs_minus.bed

awk -vOFS="\t" '{print $1, $3-40, $3-5, $4, $5, $6}'  /scratch/pooja/proteinCoding_annotatedUTRs_plus.bed  > /scratch/pooja/proteinCoding_annotatedUTRs_plus_modified.bed
awk -vOFS="\t" '{print $1, $2+5, $2+40, $4, $5, $6}'  /scratch/pooja/proteinCoding_annotatedUTRs_minus.bed  > /scratch/pooja/proteinCoding_annotatedUTRs_minus_modified.bed

cat  /scratch/pooja/proteinCoding_annotatedUTRs_plus_modified.bed /scratch/pooja/proteinCoding_annotatedUTRs_minus_modified.bed > /scratch/pooja/proteinCoding_annotatedUTRs_modified.bed
 
bedtools getfasta -fi /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -bed /scratch/pooja/proteinCoding_annotatedUTRs_modified.bed > /scratch/pooja/proteinCoding_annotatedUTRs_modified.fasta
 
findMotifs.pl /scratch/pooja/proteinCoding_annotatedUTRs_modified.fasta fasta /groups/ameres/Pooja/  -fasta /groups/ameres/Pooja//backgroundSample.fasta -b -len 6
 

