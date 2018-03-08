#!/bin/bash
#$ -t 1-87

#$ -pe smp 2-24
#$ -l hostname=!compute-6-3&!compute-6-14&!compute-6-18&!compute-6-20&!compute-6-3*&!compute-6-4*&!compute-6-5*


########$ -wd /home/<your_UCL_id>/Scratch/output

# 8. Parse parameter file to get variables.
number=$SGE_TASK_ID
#paramfile=/clustertmp/pooja/mapping_dr10_backgroundCW/sampleInfo.txt
paramFile=/groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/quantseq/allData/quantseq/sampleInfo.txt
index=`sed -n ${number}p $paramFile | awk '{print $1}'`
variable1=`sed -n ${number}p $paramfile | awk '{print $2}'`
USENAME=`sed -n ${number}p $paramfile | awk '{print $1}'`
####variable2=`sed -n ${number}p $paramfile | awk '{print $3}'`
##variable3=`sed -n ${number}p $paramfile | awk '{print $4}'`

ml cutadapt

OUTDIR=/clustertmp/pooja/mapping/

INDIR=/clustertmp/pooja/mapping_dr10_backgroundCW/



module load python
pip install --user IntervalTree
module load joblib
module load pysam
module load R/3.2.2
module load samtools/1.3.1

#/groups/ameres/bioinformatics/tools/slamdunk-new-1/slamdunk/bin/slamdunk all -r /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -o "$OUTDIR" -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/output/final90percent/allAnnotations.bed  -fb /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/ensembl_dr10_Ensembl_Genes_88/proteinCoding_annotatedUTRs.bed -a 5 -5 12 -n 100 -t 15 -mq 0 -mi 0.95 -m -rl 88 "$INDIR"/"$index" 


#### i also created counting windows from ends we throw out whrn we apply the 90%thresholding.... I want to check if there is expression in these counting windows at any stage

#OUTDIR_1=/clustertmp/pooja/mapping_lowExpression/

#/groups/ameres/bioinformatics/tools/slamdunk-new-1/slamdunk/bin/slamdunk count -r /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa  -o "$OUTDIR_1" -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/countingWindows_lowexpression.bed -l 88 -s /clustertmp/pooja/mapping/snp/ /clustertmp/pooja/mapping/filter/"$index"_adapterTrimmed_slamdunk_mapped_filtered.bam

#/groups/ameres/bioinformatics/tools/slamdunk-new-1/slamdunk/bin/slamdunk all -r /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -o "$OUTDIR" -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/output/final90percent/allAnnotations.bed  -fb /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10/ensembl_dr10_Ensembl_Genes_88/proteinCoding_annotatedUTRs.bed -a 5 -5 12 -n 100 -t 15 -mq 0 -mi 0.95 -m -rl 88 "$INDIR"/"$index" 





