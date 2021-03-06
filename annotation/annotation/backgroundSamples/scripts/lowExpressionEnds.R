
library(checkmate)


## readin in peaks, intergenic, refSeq and all other ends

BOut = "//scratch/bioinfo/pooja/SLAMannotation/dr_backgroundSamples/output/"
intergenicPeaks = read.table(paste0(BOut, "/coverage/allIntergenicPeaks_n100_new.txt"),sep="\t",stringsAsFactors = F,header = F)
#assertDataFrame(intergenicPeaks,ncols = 26)

ends_all= read.delim(paste0(BOut, "/polyAmapping_allTimepoints/n_100_global_a0/ends_all_10threshold_n100.txt"),stringsAsFactors = F,header=T)
assertDataFrame(ends_all, ncols = 22)





## rearranging the data frames

intergenicPeaks$V4 = intergenicPeaks$V16
intergenicPeaks_rearranged = cbind.data.frame(intergenicPeaks$V1,intergenicPeaks$V2, intergenicPeaks$V3, intergenicPeaks$V4, intergenicPeaks$V5, intergenicPeaks$V6, intergenicPeaks$V19,V10="intergenic") 
colnames(intergenicPeaks_rearranged) = paste("V",c(1:ncol(intergenicPeaks_rearranged)),sep="")
assertDataFrame(intergenicPeaks_rearranged,nrows = nrow(intergenicPeaks),ncols = 8)

ends_all$V4 = ends_all$V14
ends_all_rearranged = cbind.data.frame(ends_all$V1, ends_all$V2, ends_all$V3, ends_all$V4, ends_all$V5, ends_all$V6, ends_all$overlap, "UTRoverlapping")
colnames(ends_all_rearranged) = paste("V",c(1:ncol(ends_all_rearranged)),sep="")
assertDataFrame(ends_all_rearranged,nrows = nrow(ends_all),ncols = 8)

querySubject = rbind(intergenicPeaks_rearranged, ends_all_rearranged)




## get sum of the polyA reads


### getting sums of the merged positions : 

newSums = c()
for(i in 1:nrow(querySubject)){
  takeThis = sum(as.numeric(strsplit(as.character(querySubject[i,]$V5),split = ",",fixed = T)[[1]]))
  
  newSums = c(newSums,takeThis)
}
querySubject$sumCounts = newSums

cat("this worked")
### splitting data based on the gene name to compare per gene polyA read ends

querySubject_split = split(x = querySubject,f = querySubject$V4,drop = T)

## ordering this based on the number of counts - decreasing number of counts

querySubject_split_order = lapply(querySubject_split,function(x) x[order(x$sumCounts,decreasing=T),])

counts_sum = lapply(querySubject_split_order,function(x) x$sumCounts/sum(x$sumCounts))
a  =Map(cbind,querySubject_split_order,counts_sum)
a = lapply(a,function(x) x[order(x[,10],decreasing=F),])

cat("this worked")

b = lapply(a,function(x) cumsum(x[,10]))
b = Map(cbind,a,b)
b = lapply(b, function(x) x[which(x[,11]<=0.1),])
b  =lapply( b , setNames , nm = c("chr","startPeak","endPeak","external_gene_name","counts","strand","origin","peakKind","sumCounts","fractionCounts","cumSumCounts") )
cat("this worked")
totalPeaks = do.call(rbind,b)
totalPeaks_bed = totalPeaks[,c(1:6)]


onlyIntergenic = totalPeaks[which(totalPeaks$peakKind == "intergenic"),]


#onlyIntergenic_bed$peakName = onlyIntergenic$external_gene_name
onlyIntergenic_bed = onlyIntergenic[,c(1:6)]



## writing the tables



write.table(totalPeaks_bed,"//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/ends_lessThan10percent.bed",sep="\t",quote = F)

library(GenomicRanges)
