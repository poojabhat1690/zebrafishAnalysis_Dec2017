library(GenomicRanges)
totalPeaks_bed = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/ends_lessThan10percent.bed",sep="\t",stringsAsFactors = F)
totalPeaks_bed_positive = totalPeaks_bed[which(totalPeaks_bed$strand == "+"),]
annotation_custom_positive = totalPeaks_bed_positive
colnames(annotation_custom_positive) = paste0("V",c(1:6))
annotation_custom_positive$V2 = annotation_custom_positive$V3 -250
annotation_custom_positive$V2 = annotation_custom_positive$V2 + 1

annotation_custom_positive_split = split(annotation_custom_positive,f = annotation_custom_positive$V4,drop = T )

positive_ranges = lapply(annotation_custom_positive_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4)))
allAnnotations_plus_ranges_reduced = lapply(positive_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )

reducedToDf = function(reduced){
  reduced <- data.frame(seqnames=seqnames(reduced),
                        starts=start(reduced),
                        ends=end(reduced),
                        names=c(names(reduced)),
                        scores=0,strand = strand(reduced))
  return(reduced)
}

allAnnotations_plus_ranges_reduced_df = lapply(allAnnotations_plus_ranges_reduced,function(x) reducedToDf(x))

totalPeaks_bed_negative = totalPeaks_bed[which(totalPeaks_bed$strand == "-"),]
annotation_custom_negative = totalPeaks_bed_negative
colnames(annotation_custom_negative) = paste0("V",c(1:6))
annotation_custom_negative$V3 = annotation_custom_negative$V2 + 250

## changing to 1 based

annotation_custom_negative$V2 = annotation_custom_negative$V2 + 1

annotation_custom_negative_split = split(annotation_custom_negative,f = annotation_custom_negative$V4,drop = T )

negative_ranges = lapply(annotation_custom_negative_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4)))


allAnnotations_minus_ranges_reduced = lapply(negative_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )

allAnnotations_minus_ranges_reduced_df = lapply(allAnnotations_minus_ranges_reduced,function(x) reducedToDf(x))




allAnnotations_plus_ranges_reduced_df = do.call(rbind,allAnnotations_plus_ranges_reduced_df)
allAnnotations_minus_ranges_reduced_df = do.call(rbind,allAnnotations_minus_ranges_reduced_df)

allAnnotations = rbind(allAnnotations_plus_ranges_reduced_df,allAnnotations_minus_ranges_reduced_df)

### converting back to 0 based annotations : 

allAnnotations$starts = allAnnotations$starts -1

#valid chromosomes
chromosomes_mm = paste0("chr",1:25)
allAnnotations = allAnnotations[is.element(allAnnotations$seqnames, chromosomes_mm),]
write.table(allAnnotations,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/countingWindows_lowexpression.bed",quote = F,row.names = F,col.names = F,sep="\t")
