#!/bin/bash


ml R/3.4.0
ml bedtools


Rscript /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/externalData/nucleotideProfiles_ulitskyEtal/scripts/create120ntWindow.R

bedtools getfasta -s -fi /groups/ameres/bioinformatics/references/danio_rerio/dr10/danRer10.fa -bed /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/externalData/nucleotideProfiles_ulitskyEtal/data/peaks_120nts.bed > /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/externalData/nucleotideProfiles_ulitskyEtal/data/peaks_120nts.fa

Rscript /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/externalData/nucleotideProfiles_ulitskyEtal/scripts/sequences_120nts.R

