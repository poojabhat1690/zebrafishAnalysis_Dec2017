#!/bin/bash                                                                                                                                                                                                                          

######### i want to upload the following : 
	### all stages 
	### final annotation using the final sample
	### all priming sites overlapping with UTR annotations
	### ends after 90% filtering

##### just put all these in one folder 

 cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/bedFiles/




 DIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/bedFiles/



module load ucsc-kent-utils/356-foss-2017a
rm *.sorted.bed
rm *_rearranged.bed
rm *.bb

 for i in refSeq_overlapping_*.bed

 do
 
 	sort -k1,1 -k2,2n "$i" > "$i".sorted.bed
	awk -v OFS='\t' '{print $1, $2, $3, $4, 0, $6}' "$i".sorted.bed > "$i"_rearranged.bed
 	 #bedToBigBed -tab  "$i"_rearranged.bed /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.chrom.sizes "$i".bb       
 done




 for i in ensembl_overlapping_*.bed

	  do
		   
		          sort -k1,1 -k2,2n "$i" > "$i".sorted.bed
			          awk -v OFS='\t' '{print $1, $2, $3, $4, 0, $6}' "$i".sorted.bed > "$i"_rearranged.bed
				           #bedToBigBed -tab  "$i"_rearranged.bed /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.chrom.sizes "$i".bb       
	 done



###### merging the refSeq and ENSEMBL files together 

	 files=( "1dpf" "2dpf" "4dpf" "256cell" "2cell" "bud" "dome" "oocyte" "sphere" "testis" )
	
	 for i in "${files[@]}"

	 do
		cat refSeq_overlapping_"$i".bed_rearranged.bed ensembl_overlapping_"$i".bed_rearranged.bed > refSeq_ensembl_"$i".bed
		sort -k1,1 -k2,2n refSeq_ensembl_"$i".bed > refSeq_ensembl_"$i"_sorted.bed
		bedToBigBed -tab  refSeq_ensembl_"$i"_sorted.bed  /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.chrom.sizes refSeq_ensembl_"$i".bb     
	done





	################# making similar bed files with all the ends (>90 percent) are... 



for i in ends_greaterthan90percent_*


do
	sort -k1,1 -k2,2n "$i" > "$i".sorted.bed
	sed -i '/many/d' "$i".sorted.bed
	awk -v OFS='\t' '{print $1, $2, $3, $4, 0, $6}' "$i".sorted.bed > "$i"_rearranged.bed
	bedToBigBed -tab  "$i"_rearranged.bed  /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.chrom.sizes "$i".bb
done


