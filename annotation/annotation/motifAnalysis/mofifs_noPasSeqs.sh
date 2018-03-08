#!/bin/bash
ml homer/4.9-foss-2017a

cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/
rm *_proper.fasta

for i in noPASaccepted_seque*.fasta

do

awk '{print ">\n" $0;}' "$i" > "$i"_proper.fasta
mkdir /scratch/pooja/"$i"
#findMotifs.pl targets.fa fasta motifResults/ -fasta background.fa
findMotifs.pl  "$i"_proper.fasta fasta /scratch/pooja/"$i" -fasta /groups/ameres/Pooja//backgroundSample.fasta -b 

done


#### doing the same for all the sequences of noPAS



#awk '{print ">\n" $0;}' allNoPAS_SLAMdunkExperiment.fasta > allNoPAS_SLAMdunkExperiment.fasta_proper.fasta


#findMotifs.pl  allNoPAS_SLAMdunkExperiment.fasta_proper.fasta fasta /clustertmp/pooja/allNoPAS_SLAMdunkExperiment -fasta /groups/ameres/Pooja//backgroundSample.fasta -b -len 4
