#!/bin/bash

##### getting the number of 'unique primed sites'
##### number of primed sites after thresholding
##### number after merging


#### i want to have these statistics for 'all samples' and the background samples

### for all samples : 

DIR=/clustertmp/bioinfo/pooja/SLAMannotation/dr_allData/output/polyAmapping_allTimepoints/n_100_global_a0/




print_something () {

count_MINUS=$(wc -l <$1) 
count_PLUS=$(wc -l <$2) 
SUM=`expr $count_MINUS + $count_PLUS`

echo "$SUM"


}


### for all data 
touch "$DIR"/allData_primingSiteStats.txt

primingSITES=`print_something "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed` 
echo $primingSITES uniquePrimingSites >"$DIR"/allData_primingSiteStats.txt


Thresholded=`print_something "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan2.bed "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan2.bed`
echo $Thresholded thresholded >> "$DIR"/allData_primingSiteStats.txt




merged=`print_something "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan2.bed_sorted_merged.bed "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan2.bed_sorted_merged.bed` 
echo "$merged" merged >> "$DIR"/allData_primingSiteStats.txt



################ for background samples 


DIR=/clustertmp/bioinfo/pooja/SLAMannotation/dr_backgroundSamples/output/polyAmapping_allTimepoints/n_100_global_a0/

touch "$DIR"/allData_primingSiteStats.txt

primingSITES=`print_something "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed`
echo $primingSITES uniquePrimingSites >"$DIR"/allData_primingSiteStats.txt


Thresholded=`print_something "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan2.bed "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan2.bed`
echo $Thresholded thresholded >> "$DIR"/allData_primingSiteStats.txt




merged=`print_something "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan2.bed_sorted_merged.bed "$DIR"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan2.bed_sorted_merged.bed`
echo "$merged" merged >> "$DIR"/allData_primingSiteStats.txt


















