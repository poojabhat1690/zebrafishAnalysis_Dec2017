#### comparing pre-processing of data between all samples and background samples



allOutput = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/allData/data/logs/",pattern = "preProcessingNumbers")
allOutFiles = paste("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/allData/data/logs/",allOutput,sep="")

allOutputData = lapply(allOutFiles,function(x) read.table(x,sep=":"))
names(allOutputData) = allOutput

intialFiles = unlist(lapply(allOutputData,function(x) x[1,2]))
adapterTrimmed = unlist(lapply(allOutputData,function(x) x[2,2]))
fivePrimeTrimmig = unlist(lapply(allOutputData,function(x) x[3,2]))
polyAcontaining = unlist(lapply(allOutputData,function(x) x[4,2]))
totalStats = cbind.data.frame(intialFiles,adapterTrimmed,fivePrimeTrimmig,polyAcontaining)
colnames(totalStats) = c("InitialFiles","AdapterTrimmed","FivePrimeTrimmed","PolyAReads")
totalStats_allData = totalStats
sum_reads = colSums(totalStats)
mean.n <- function(x){
  return(c(y = median(x)*0.4, label = round(sum(x),2))) 
  # experiment with the multiplier to find the perfect position
}


library(ggplot2)
library(reshape)
totalStats_melt = melt(totalStats)
pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/allData_preProcessing.pdf",height=5,width=7)
p = ggplot(totalStats_melt,aes(x=variable,y=value)) +geom_jitter(fill="grey") + xlab("Pre processing step") + ylab("Number of reads") + theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "red")
p = p + theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(margin=margin(5,15,10,5,"pt")), axis.ticks.length = unit(-0.25 , "cm")) + ggtitle("Number of reads in the samples in preprocessing steps")
print(p)
dev.off()
totalStats_melt_allData = totalStats_melt

###### 


allOutput = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples//data/logs/",pattern = "preProcessingNumbers")
allOutFiles = paste("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples//data/logs/",allOutput,sep="")

allOutputData = lapply(allOutFiles,function(x) read.table(x,sep=":"))
names(allOutputData) = allOutput

intialFiles = unlist(lapply(allOutputData,function(x) x[1,2]))
adapterTrimmed = unlist(lapply(allOutputData,function(x) x[2,2]))
fivePrimeTrimmig = unlist(lapply(allOutputData,function(x) x[3,2]))
polyAcontaining = unlist(lapply(allOutputData,function(x) x[4,2]))
totalStats_background = cbind.data.frame(intialFiles,adapterTrimmed,fivePrimeTrimmig,polyAcontaining)
colnames(totalStats_background) = c("InitialFiles","AdapterTrimmed","FivePrimeTrimmed","PolyAReads")

sum_reads = colSums(totalStats)
mean.n <- function(x){
  return(c(y = median(x)*4, label = round(sum(x),2))) 
  # experiment with the multiplier to find the perfect position
}


library(ggplot2)
library(reshape)
totalStats_melt = melt(totalStats_background)

p = ggplot(totalStats_melt,aes(x=variable,y=value)) +geom_jitter(fill="grey") + xlab("Pre processing step") + ylab("Number of reads") + theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "red")
p = p + theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(margin=margin(5,15,10,5,"pt")), axis.ticks.length = unit(-0.25 , "cm")) + ggtitle("Number of reads in the samples in preprocessing steps")
print(p)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/backgroundData_preProcessing.pdf",height=5,width=13)
p = ggplot(totalStats_melt,aes(x=variable,y=value)) +geom_jitter(fill="grey") + xlab("Pre processing step") + ylab("Number of reads") + theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "red")
p = p + theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(margin=margin(5,15,10,5,"pt")), axis.ticks.length = unit(-0.25 , "cm")) + ggtitle("Number of reads in the samples in preprocessing steps")
print(p)
dev.off()

totalStats_melt_allData$category = "allData"
totalStats_melt$category = "background"

allStats = rbind(totalStats_melt_allData,totalStats_melt)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/preProcessingStats.pdf",height=5,width=16)
p = ggplot(allStats,aes(x=variable,y=value,col=category,fill=category)) +geom_jitter(fill="grey") + facet_grid(.~category,scales = "free")+ xlab("Pre processing step") + ylab("Number of reads") + theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+stat_summary(fun.data = mean.n, geom = "text", fun.y = mean, colour = "red")
p = p + theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(margin=margin(5,15,10,5,"pt")), axis.ticks.length = unit(-0.25 , "cm")) + ggtitle("Number of reads in the samples in preprocessing steps")
p + scale_color_manual(values = c("allData" = "blue","background" = "red"))   + theme(axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),axis.title  = element_text(size=14),legend.text=element_text(size=0),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length = unit(-0.25 , "cm"))  
dev.off()               

###.293377 polyA reads for all data
### 0.2637299 polyA reads fpr backfrpind samples

