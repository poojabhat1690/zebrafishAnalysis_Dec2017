#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-8:00:00     # 2 minutes
#SBATCH --output=addingToIntergenicEnds.txt
#SBATCH --job-name=allData_intergenicEnds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at

BOut="//scratch/bioinfo/pooja/SLAMannotation/dr_allData//output/"
BIn="//groups/ameres/Pooja/Projects/zebrafishAnnotation/dr10_data/quantseq/allData/"

module load bedtools/2.27.1-foss-2017a 
module load parallel/20171122-foss-2017a

args1=200
args2=0.1

OUTparallel="$BOut"/parallel/
mkdir -p "$OUTparallel"
 
cd "$BIn"/rnaseq/
 
find *bam | parallel 'bedtools multicov -split -p -bams {} -bed /scratch/bioinfo/pooja/SLAMannotation/dr_allData//output//ExtendingINtergenicRegions/customAnnotation_longestTranscripts_100IntoUTR.bed > /scratch/bioinfo/pooja/SLAMannotation/dr_allData//output/parallel/{}.txt'


cd "$OUTparallel"


for f in *.sortedByCoord.out.bam.txt ; do cut -f9 $f > $f.tmp ; done


paste *.sortedByCoord.out.bam.txt.tmp > tmp_cumulative.txt


FILEname=`ls *.sortedByCoord.out.bam.txt | head -n 1`


awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' "$FILEname" > "$OUTparallel"/constant.txt


paste constant.txt tmp_cumulative.txt > customAnnotation_longestTranscripts_IntoUTR_coverage_"$args1"_"$args2".bed

###### this is the into UTR part... 


awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}' tmp_cumulative.txt #### this is creating a sum of the signal 
awk '{ s = 0; for (i = 1; i <= NF; i++) s += $i; print (NF > 1) ? s / (NF - 1) : 0; }'  tmp_cumulative.txt > mean_counts.txt

paste customAnnotation_longestTranscripts_IntoUTR_coverage_"$args1"_"$args2".bed mean_counts.txt > customAnnotation_longestTranscripts_IntoUTR_coverage_"$args1"_"$args2"_withMean.bed

ncol=`awk -F'\t' '{print NF; exit}' customAnnotation_longestTranscripts_IntoUTR_coverage_"$args1"_"$args2"_withMean.bed`  #### storing number of columns

######## first I want to combine the mean counts with the custom annotaion and then subset for the means. 

#### i want only a subset of the mean count file ... i.e the gene name and the mean counts  (in UTR)

###cut -f4,"$ncol"  customAnnotation_longestTranscripts_IntoUTR_coverage_"$args1"_"$args2"_withMean.bed > subset_mean_name.txt

###"$BOut"/toExtend_longestEnsembl_refSeq_n100_sorted_distances.bed

paste /scratch/bioinfo/pooja/SLAMannotation/dr_allData//output//ExtendingINtergenicRegions/customAnnotation_longestTranscripts_100IntoUTR.bed mean_counts.txt > customAnnotationWithMean.txt


THRESHOLD=10
##### gettign all entries with mean >10 
ncol_customAnnotation=`awk -F'\t' '{ print NF; exit}' customAnnotationWithMean.txt`  #### storing number of columns

awk -vT=$THRESHOLD '{ if ($9 >=T) print  }' customAnnotationWithMean.txt  >  customAnnotationWithMean_"$THRESHOLD".bed 


########## so only these will be extended...


seq 20 20 39800 > starts.txt
seq 220 20 40000 > ends.txt

paste starts.txt ends.txt > thresholds.txt

#################################### now I want to create a loop to do some stuff


#### splitting my custom annotation with mean... 

cat customAnnotationWithMean_"$THRESHOLD".bed   | awk '$6 == "+" { print $0 }' > customAnnotationWithMean_"$THRESHOLD"_plus.bed 
cat customAnnotationWithMean_"$THRESHOLD".bed   | awk '$6 == "-" { print $0 }' > customAnnotationWithMean_"$THRESHOLD"_minus.bed 

#### i want to loop over the number of lines of the file thresholds 

nLines=`wc -l thresholds.txt`

i=1
while [ "$i" -le 1990 ]

do
	
	awk 'BEGIN {OFS="\t"} {print $1,$2+20,$3+20,$4,$5,$6,$7,$8,$9}' customAnnotationWithMean_"$THRESHOLD"_plus.bed > customAnnotationWithMean_"$THRESHOLD"_plus_tmp.bed   && mv customAnnotationWithMean_"$THRESHOLD"_plus_tmp.bed customAnnotationWithMean_"$THRESHOLD"_plus.bed
	awk 'BEGIN {OFS="\t"} {print $1,$2-20,$3-20,$4,$5,$6,$7,$8,$9}' customAnnotationWithMean_"$THRESHOLD"_minus.bed > customAnnotationWithMean_"$THRESHOLD"_minus_tmp.bed  && mv customAnnotationWithMean_"$THRESHOLD"_minus_tmp.bed customAnnotationWithMean_"$THRESHOLD"_minus.bed
	
	
	cat customAnnotationWithMean_"$THRESHOLD"_plus.bed customAnnotationWithMean_"$THRESHOLD"_minus.bed > customAnnotationWithMean_"$THRESHOLD".bed

#### just checking if all coordinated are still positive
	awk '{ if ($2 > 0) print  }' customAnnotationWithMean_"$THRESHOLD".bed > customAnnotationWithMean_"$THRESHOLD"_tmp.bed && mv  customAnnotationWithMean_"$THRESHOLD"_tmp.bed customAnnotationWithMean_"$THRESHOLD".bed
	awk '{ if ($3 > 0) print  }' customAnnotationWithMean_"$THRESHOLD".bed > customAnnotationWithMean_"$THRESHOLD"_tmp.bed && mv  customAnnotationWithMean_"$THRESHOLD"_tmp.bed customAnnotationWithMean_"$THRESHOLD".bed

	nLinesTMP=`wc -l customAnnotationWithMean_"$THRESHOLD".bed`	
	
	if [ "$nLinesTMP" == 0 ] ; then
      echo 'No more genes'
      break
   else
      echo 'condition not met still searching'
   fi
	   
#  awk 'FNR == 5 {print $1}' thresholds.txt

  startTMP=`awk -v i="$i" -v j=1 'FNR == i {print $j}' thresholds.txt`
  endTMP=`awk -v i="$i" -v j=2 'FNR == i {print $j}' thresholds.txt`

   cd "$BIn"/rnaseq/
   

	
   find *bam | parallel 'bedtools multicov -split -p -bams {} -bed /scratch/bioinfo/pooja/SLAMannotation/dr_allData/output/parallel/customAnnotationWithMean_10.bed > /scratch/bioinfo/pooja/SLAMannotation/dr_allData//output/parallel/{}_params.txt' 

	##### manipulaitng this to put all the counts together 
	cd "$OUTparallel"

	mkdir -p OUTnew
	
	for f in *_params.txt ; do cut -f10 $f > OUTnew/$f.tmp; done

	paste OUTnew/*t.bam_params.txt.tmp > OUTnew/tmp_cumulative.txt


	FILEname=`ls *_params.txt | head -n 1`
	awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' "$FILEname" > "$OUTparallel"/OUTnew/constant.txt

	awk '{ s = 0; for (i = 1; i <= NF; i++) s += $i; print (NF > 1) ? s / (NF - 1) : 0; }'  OUTnew/tmp_cumulative.txt > OUTnew/mean_counts.txt
	
	paste "$OUTparallel"/OUTnew/constant.txt OUTnew/mean_counts.txt > OUTnew/counts_offSet.txt
	awk -v OFS='\t' '{$11 = $10 / $9}1' OUTnew/counts_offSet.txt > OUTnew/counts_offSet_tmp.txt && mv OUTnew/counts_offSet_tmp.txt OUTnew/counts_offSet.txt
	
	awk '{ if ($11 > 0.1) print  }' OUTnew/counts_offSet.txt > OUTnew/counts_offSet_tmp.txt && mv  OUTnew/counts_offSet_tmp.txt OUTnew/counts_offSet.txt

   
   awk ' $8 > "$startTMP" ' OUTnew/counts_offSet.txt > OUTnew/counts_offSet_tmp.txt && mv  OUTnew/counts_offSet_tmp.txt OUTnew/counts_offSet.txt
	mkdir -p extensions
	cp OUTnew/counts_offSet.txt  extensions/counts_extension_"$startTMP"_"$endTMP".txt   
   cat OUTnew/counts_offSet.txt   | awk '$6 == "+" { print $0 }' > customAnnotationWithMean_"$THRESHOLD"_plus.bed 
   cat OUTnew/counts_offSet.txt  | awk '$6 == "-" { print $0 }' > customAnnotationWithMean_"$THRESHOLD"_minus.bed 

 i=`expr "$i" + 1`
done
		
