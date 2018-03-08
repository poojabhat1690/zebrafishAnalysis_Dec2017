#### this is for creating a combined summary file for quantSeq samples 

### reading in the summary files : 


samples_zebrafish = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/summaries/",pattern = ".txt")
sample_zebrafish_path = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/summaries//", samples_zebrafish)
allSummedUP = lapply(sample_zebrafish_path,function(x) read.delim(x,header=T,stringsAsFactors = F))

referenceFile = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/differentialExpressionAnalysis/info/sampleInfo.txt",header = F,stringsAsFactors = F)
referenceFile$V1 = paste0(referenceFile$V1,"_adapterTrimmed_slamdunk_mapped_filtered.bam")

colnames(referenceFile) = c("getFileNames","V2","V3")
#allSummedUP = lapply(allSummedUP,function(x) x[-1,])
allSummedUP = do.call(rbind,allSummedUP)
#colnames(allSummedUP) = c("FileName","SampleName","SampleType","SampleTime","Sequenced","Mapped","Deduplicated","MQ-Filtered","Identity-Filtered","NM-Filtered","Multimap-Filtered","Retained","Annotation")

allSummedUP$RPMfac = 1000000/allSummedUP$Retained

getFileNames = strsplit(allSummedUP$FileName,"/",T)
getFileNames = unlist(lapply(getFileNames,function(x) x[7]))
allSummedUP= cbind(allSummedUP,getFileNames)
allSummedUP$getFileNames= as.character(allSummedUP$getFileNames)

allSummedUP = plyr::join(allSummedUP,referenceFile)
allSummedUP$SampleType = 1
allSummedUP = allSummedUP[complete.cases(allSummedUP),]
write.table(allSummedUP,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data//summary_ucscBrowser.txt",sep="\t",quote = F,row.names = F,col.names = F)
