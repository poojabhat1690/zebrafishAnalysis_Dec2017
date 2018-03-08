#!/bin/bash

#$ -t 1-26

#$ -pe smp 1-24

########$ -wd /home/<your_UCL_id>/Scratch/output

# 8. Parse parameter file to get variables.
number=$SGE_TASK_ID
paramfile=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data//summary_ucscBrowser.txt

index=`sed -n ${number}p $paramfile | awk '{print $16}'`
variable1=`sed -n ${number}p $paramfile | awk '{print $2}'`
FAC=`sed -n ${number}p $paramfile | awk '{print $15}'`

USENAME=`sed -n ${number}p $paramfile | awk '{print $1}'`
####variable2=`sed -n ${number}p $paramfile | awk '{print $3}'`
##variable3=`sed -n ${number}p $paramfile | awk '{print $4}'`





ml deeptools/2.5.0.1-python2.7.3
ml pysam/0.10.0


cd  /clustertmp/pooja//mapping_backgrouncountingWindows_final/filter/


OUTDIR=/clustertmp/pooja/mapping_backgrouncountingWindows_final/ucscBrowser/
mkdir "$OUTDIR"

ml samtools



bamCoverage  -b "$index" -o "$OUTDIR"/"$index"_plus.bg --filterRNAstrand reverse --scaleFactor "$FAC" --binSize 1 -of bedgraph

bamCoverage  -b "$index" -o "$OUTDIR"/"$index"_minus.bg --filterRNAstrand forward --scaleFactor "$FAC" --binSize 1 -of bedgraph


awk '{OFS="\t"; print $1,$2,$3,"-"$4}' "$OUTDIR"/"$index"_minus.bg >"$OUTDIR"/"$index"_neg_minus.bg


ml kent-ucsc/3.8f6f5e0a1cb75

bedGraphToBigWig "$OUTDIR"/"$index"_neg_minus.bg /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.chrom.sizes  "$OUTDIR"/"$index"_neg_minus.bg.bigWig


bedGraphToBigWig "$OUTDIR"/"$index"_plus.bg /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.chrom.sizes  "$OUTDIR"/"$index"_plus.bg.bigWig





cp "$OUTDIR"/*bigWig  /groups/ameres/wwwdata/seq-hub/dr_slamSeq/dr10/

