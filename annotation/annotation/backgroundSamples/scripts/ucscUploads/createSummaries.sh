#!/bin/bash



module load python
pip install --user IntervalTree
module load joblib
module load pysam
module load R/3.2.2
module load samtools/1.3.1





DIR=/clustertmp/pooja/mapping_backgrouncountingWindows_final/
mkdir -p "$DIR"/summaries/


cd  "$DIR"/filter/

#for i in *.bam

#do

#/groups/ameres/Veronika/bin/slamdunk/bin/alleyoop summary -o "$DIR"/summaries/"$i".txt -t "$DIR"/count/  "$DIR"/filter/"$i"

#sed '1d' "$DIR"/summaries/"$i".txt > "$DIR"/summaries/"$i"_tmpfile.txt; mv "$DIR"/summaries/"$i"_tmpfile.txt "$DIR"/summaries/"$i".txt

#done

mkdir -p /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/summaries/

cp "$DIR"/summaries/*_slamdunk_mapped_filtered.bam.txt /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/summaries/
