#!/bin/bash



#### want to loop through these and get the ensembl and refSeq overlapping ends and rename them

files=( "1dpf" "2dpf" "4dpf" "256cell" "2cell" "bud" "dome" "oocyte" "sphere" "testis" )
OUTDIR="/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/bedFiles/"
for i in "${files[@]}"
do
   echo "$i"
   cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/2016_12_QuantSeq_Annotaton/raw_data/"$i"/output_includingPAS_ensembl/polyAmapping_allTimepoints/n_100_global_a0/
   pwd
   cp refSeq_overlapping.bed.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/bedFiles/refSeq_overlapping_"$i".bed.gz
   cp ensembl_overlapping.bed.gz /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/bedFiles/ensembl_overlapping_"$i".bed.gz
 
   gunzip /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/bedFiles/refSeq_overlapping_"$i".bed.gz
   gunzip /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/bedFiles/ensembl_overlapping_"$i".bed.gz
  
   cp /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/2016_12_QuantSeq_Annotaton/raw_data/"$i"/output_includingPAS_ensembl/final90percent/ends_greater90percent_intergenic_n100.bed /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/bedFiles/ends_greaterthan90percent_"$i".bed

done



