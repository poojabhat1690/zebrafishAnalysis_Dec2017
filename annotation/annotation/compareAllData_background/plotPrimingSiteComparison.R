##### comparing primin sites.. 
library(ggplot2)

allData = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/primingSites_allData.txt",stringsAsFactors = F)
background = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/primingSites_onlyBackground.txt",stringsAsFactors = F)

allData$type = "alldata"
background$type = "background"

allData$V2 = factor(allData$V2,levels = allData$V2)
background$V2 = factor(background$V2,levels = background$V2)
allPrimingSites = rbind(allData,background)
allPrimingSites$V1 = allPrimingSites$V1/1000000

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/primingSitesComparison.pdf",height = 4, width=8)
p = ggplot(allPrimingSites,aes(x=V2,y=V1,group=type,fill=type)) + geom_bar(stat = "identity",position = "dodge") + theme_bw() + xlab(" ") + ylab("Number of priming sites (in million)")
p = p + theme(axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),axis.title  = element_text(size=14),legend.position=c(0.8,.8),legend.text=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length = unit(-0.25 , "cm"))  
p + scale_fill_manual(values = c("alldata" = "blue","background" = "red"))
dev.off()


######  comparing the number of priming sites overlapping with different annotations... 



background_geneAnnotation_allData = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/nucleotideOverlap_summarry_allData.txt")
background_geneAnnotation_background = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/nucleotideOverlap_summarry_background.txt")

### fraction of signal

background_geneAnnotation_allData$type = "alldata"
background_geneAnnotation_background$type = "background"

background_geneAnnotation_allData$category = factor(background_geneAnnotation_allData$category, levels = background_geneAnnotation_allData$category)
background_geneAnnotation_background$category = factor(background_geneAnnotation_background$category, levels = background_geneAnnotation_background$category)
allSamples = rbind(background_geneAnnotation_allData,background_geneAnnotation_background)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/fractionOfendsAndSignal.pdf",height=4,width=9)
p = ggplot(allSamples,aes(x=category,y=fractionOfReads,group=type,fill=type)) + geom_bar(stat="identity",position="dodge") + theme_bw() + xlab(" ") + theme(axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),axis.title  = element_text(size=14),legend.position=c(0.25,.8),legend.text=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length = unit(-0.25 , "cm"))  
p = p + scale_fill_manual(values = c("alldata" = "blue","background" = "red"))
print(p)
q = ggplot(allSamples,aes(x=category,y=fractionOfPeaks,group=type,fill=type)) + geom_bar(stat="identity",position="dodge") + theme_bw() + xlab(" ") + theme(axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),axis.title  = element_text(size=14),legend.position=c(0.25,.8),legend.text=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length = unit(-0.25 , "cm"))  
q  = q + scale_fill_manual(values = c("alldata" = "blue","background" = "red"))
print(q)
dev.off()


######

######## plotting the fraction of accepted and rejected sites for each

differentOverlaps_allData = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/differentCategories_nucleorideProfiles_allData.txt")

PAScontaining = differentOverlaps_allData[which(differentOverlaps_allData$PASsite == "AATAAA"| differentOverlaps_allData$PASsite == "ATTAAA"|differentOverlaps_allData$PASsite == "APA" ),]
noPAScontaining = differentOverlaps_allData[which(differentOverlaps_allData$PASsite == "noPAS"),]
PAScontaining_accepted = PAScontaining[c(1:9),]
noPAScontaining_accepted = noPAScontaining[c(1:2),]
allAccepted = rbind(PAScontaining_accepted,noPAScontaining_accepted)

#### i want to find, out of each category, how many are accepted... i.e how many ends out of the total intergenic ends are accepted.
sums_category  = colSums(differentOverlaps_allData[,c(3:7)])
sums_accepted = colSums(allAccepted[,c(3:7)])
fraction_accepted_allData = sums_accepted/sums_category

fraction_accepted_allData = as.data.frame(fraction_accepted_allData)
fraction_accepted_allData$catogory = row.names(fraction_accepted_allData)
differentOverlaps_background = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/differentCategories_nucleorideProfiles_background.txt")


PAScontaining = differentOverlaps_background[which(differentOverlaps_background$PASsite == "AATAAA"| differentOverlaps_background$PASsite == "ATTAAA"|differentOverlaps_background$PASsite == "APA" ),]
noPAScontaining = differentOverlaps_background[which(differentOverlaps_background$PASsite == "noPAS"),]
PAScontaining_accepted = PAScontaining[c(1:9),]
noPAScontaining_accepted = noPAScontaining[c(1:2),]
allAccepted = rbind(PAScontaining_accepted,noPAScontaining_accepted)

#### i want to find, out of each category, how many are accepted... i.e how many ends out of the total intergenic ends are accepted.
sums_category  = colSums(differentOverlaps_background[,c(3:7)])
sums_accepted = colSums(allAccepted[,c(3:7)])
fraction_accepted_bg = sums_accepted/sums_category
fraction_accepted_bg = as.data.frame(fraction_accepted_bg)
fraction_accepted_bg$catogory = row.names(fraction_accepted_bg)

fraction_accepted_allData$type = "alldata"
fraction_accepted_bg$type = "background"
colnames(fraction_accepted_allData) = c("V1","V2","V3")
colnames(fraction_accepted_bg) = c("V1","V2","V3")
fraction_accepted_allData$V2 = factor(fraction_accepted_allData$V2, levels = fraction_accepted_allData$V2)
fraction_accepted_bg$V2 =  factor(fraction_accepted_bg$V2, levels = fraction_accepted_bg$V2)

totalFraction = rbind(fraction_accepted_allData,fraction_accepted_bg)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/fractionOfEndsSurvivingCutoff.pdf",height=4,width=10)
r = ggplot(totalFraction,aes(x=V2,y=V1,group=V3,fill=V3)) + geom_bar(stat = "identity",position = "dodge")+ theme_bw() + xlab(" ") + theme(axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),axis.title  = element_text(size=14),legend.position=c(0.8,.8),legend.text=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length = unit(-0.25 , "cm"))  
r + scale_fill_manual(values = c("alldata" = "blue","background" = "red")) + ylab("Fraction of ends accepted") + xlab(" ")
dev.off()

####### now i also want to check fraction of surviving signal

######### 
backgroundSignal = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/signal_summary_background.txt")
allSignal = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/signal_summary_allData.txt")


###################### allData
PAScontaining = allSignal[which(allSignal$L2 == "AATAAA"| allSignal$L2 == "ATTAAA"|allSignal$L2 == "APA" ),]
noPAScontaining = allSignal[which(allSignal$L2 == "noPAS"),]
PAScontaining_accepted = PAScontaining[c(1:9),]
noPAScontaining_accepted = noPAScontaining[c(1:2),]
allAccepted = rbind(PAScontaining_accepted,noPAScontaining_accepted)

sums_category  = colSums(allSignal[,c(3:7)])
sums_accepted = colSums(allAccepted[,c(3:7)])
fraction_accepted_all = sums_accepted/sums_category
fraction_accepted_all = as.data.frame(fraction_accepted_all)
fraction_accepted_all$catogory = row.names(fraction_accepted_all)


##### background data

PAScontaining = backgroundSignal[which(backgroundSignal$L2 == "AATAAA"| backgroundSignal$L2 == "ATTAAA"|backgroundSignal$L2 == "APA" ),]
noPAScontaining = backgroundSignal[which(backgroundSignal$L2 == "noPAS"),]
PAScontaining_accepted = PAScontaining[c(1:9),]
noPAScontaining_accepted = noPAScontaining[c(1:2),]
allAccepted = rbind(PAScontaining_accepted,noPAScontaining_accepted)

sums_category  = colSums(backgroundSignal[,c(3:7)])
sums_accepted = colSums(allAccepted[,c(3:7)])
fraction_accepted_bg = sums_accepted/sums_category
fraction_accepted_bg = as.data.frame(fraction_accepted_bg)
fraction_accepted_bg$catogory = row.names(fraction_accepted_bg)

fraction_accepted_all$type = "alldata"
fraction_accepted_bg$type = "background"
colnames(fraction_accepted_all) = c("V1","V2","V3")
colnames(fraction_accepted_bg) = c("V1","V2","V3")
fraction_accepted_all$V2 = factor(fraction_accepted_all$V2, levels = fraction_accepted_all$V2)
fraction_accepted_bg$V2 =  factor(fraction_accepted_bg$V2, levels = fraction_accepted_bg$V2)
totalFraction = rbind(fraction_accepted_all,fraction_accepted_bg)
pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/fractionOfacceptedSignal.pdf",height = 4,width=15)
r = ggplot(totalFraction,aes(x=V2,y=V1,group=V3,fill=V3)) + geom_bar(stat = "identity",position = "dodge")+ theme_bw() + xlab(" ") + theme(axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),axis.title  = element_text(size=14),legend.position=c(0.8,.8),legend.text=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length = unit(-0.25 , "cm"))  
r + scale_fill_manual(values = c("alldata" = "blue","background" = "red")) + ylab("Fraction of signal accepted") + xlab(" ")
dev.off()
###################### just making cumulative plots instead

library(dplyr)

overlappers_annotations_allData = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/differentCategories_nucleorideProfiles_allData.txt",stringsAsFactors = F)
overlappers_annotations_background = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/differentCategories_nucleorideProfiles_background.txt",stringsAsFactors = F)

plotCululative = function(overlappers_annotations){
overlappers_annotations= cbind.data.frame(overlappers_annotations)
colnames(overlappers_annotations) = c("Acontent", "PASsite" ,"refSeqOverlappers", "ensemblOverlappers","exonOverlappers", "intronOverlappers","nonOverlapping")
                                    
overlappers_annotations_PAS = overlappers_annotations  %>% filter(PASsite != "noPAS") %>%group_by(Acontent) %>% mutate(sum_refSeqOverlappers = sum(refSeqOverlappers))  %>% mutate(sum_ensemblOverlappers = sum(ensemblOverlappers)) %>% mutate(sum_exonOverlappers = sum(exonOverlappers)) %>%  mutate(sum_intronOverlappers = sum(intronOverlappers)) %>% mutate(sum_nonOverlapping = sum(nonOverlapping))
overlappers_annotations_PAS = overlappers_annotations_PAS[!duplicated(overlappers_annotations_PAS[,c("Acontent","sum_refSeqOverlappers","sum_ensemblOverlappers","sum_exonOverlappers","sum_intronOverlappers","sum_nonOverlapping")]),]


overlappers_annotations_PAS= overlappers_annotations_PAS[,c("Acontent","sum_refSeqOverlappers","sum_ensemblOverlappers","sum_exonOverlappers","sum_intronOverlappers","sum_nonOverlapping")] 
rowSums_overlappers_PAS = rowSums(overlappers_annotations_PAS[,c("sum_refSeqOverlappers","sum_ensemblOverlappers","sum_exonOverlappers","sum_intronOverlappers","sum_nonOverlapping")])
overlappers_annotations_PAS[,c("sum_refSeqOverlappers","sum_ensemblOverlappers","sum_exonOverlappers","sum_intronOverlappers","sum_nonOverlapping")]/rowSums_overlappers_PAS
overlappers_annotations_PAS_cumSum  = cumsum(overlappers_annotations_PAS[,c("sum_refSeqOverlappers","sum_ensemblOverlappers","sum_exonOverlappers","sum_intronOverlappers","sum_nonOverlapping")])
sum_refSeqEnsemblAccepted = sum(overlappers_annotations_PAS[c(1:3),c(2:3)])

cumSum_PAS = apply(overlappers_annotations_PAS_cumSum,2,function(x) x/max(x))
row.names(cumSum_PAS) = c("0-0.12","0.12-0.24","0.24-0.36","0.36-0.48","0.48-0.60","0.60-0.72","0.72-0.84","0.84-1.00")
colnames(cumSum_PAS) = c("refSeq 3' UTR overlapping","ENSEMBL 3' UTR overlapping","Exon overlapping","Intron overlapping","Intergenic overlapping")
cumSum_PAS_melt = melt(cumSum_PAS)

p = ggplot(cumSum_PAS_melt,aes(x=X1,y=value,group=X2)) +geom_line(aes(col=X2),size=1) +theme_bw() + ylab("Cumulative fraction (number of priming sites)") + xlab("Downstream A content")  
p = p+ theme(legend.title=element_blank(),axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),legend.position=c(0.8,.4),legend.text=element_text(size=14), axis.ticks.length = unit(-0.25 , "cm"))  + theme(axis.text=element_text(size=12),axis.title=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("PAS containing")+ geom_vline(xintercept = 3,linetype="dashed",color="red")
p = p + annotate(geom = "text",x = 6,y = 0.1,label = paste0("Number of accepted UTR overlapping ends=", sum_refSeqEnsemblAccepted))


overlappers_annotations_PAS = overlappers_annotations  %>% filter(PASsite == "noPAS") %>%group_by(Acontent) %>% mutate(sum_refSeqOverlappers = sum(refSeqOverlappers))  %>% mutate(sum_ensemblOverlappers = sum(ensemblOverlappers)) %>% mutate(sum_exonOverlappers = sum(exonOverlappers)) %>%  mutate(sum_intronOverlappers = sum(intronOverlappers)) %>% mutate(sum_nonOverlapping = sum(nonOverlapping))
overlappers_annotations_PAS = overlappers_annotations_PAS[!duplicated(overlappers_annotations_PAS[,c("Acontent","sum_refSeqOverlappers","sum_ensemblOverlappers","sum_exonOverlappers","sum_intronOverlappers","sum_nonOverlapping")]),]


overlappers_annotations_PAS= overlappers_annotations_PAS[,c("Acontent","sum_refSeqOverlappers","sum_ensemblOverlappers","sum_exonOverlappers","sum_intronOverlappers","sum_nonOverlapping")] 
rowSums_overlappers_PAS = rowSums(overlappers_annotations_PAS[,c("sum_refSeqOverlappers","sum_ensemblOverlappers","sum_exonOverlappers","sum_intronOverlappers","sum_nonOverlapping")])
overlappers_annotations_PAS[,c("sum_refSeqOverlappers","sum_ensemblOverlappers","sum_exonOverlappers","sum_intronOverlappers","sum_nonOverlapping")]/rowSums_overlappers_PAS
overlappers_annotations_PAS_cumSum  = cumsum(overlappers_annotations_PAS[,c("sum_refSeqOverlappers","sum_ensemblOverlappers","sum_exonOverlappers","sum_intronOverlappers","sum_nonOverlapping")])
sum_refSeqEnsemblAccepted = sum(overlappers_annotations_PAS[c(1:3),c(2:3)])

cumSum_PAS = apply(overlappers_annotations_PAS_cumSum,2,function(x) x/max(x))
row.names(cumSum_PAS) = c("0-0.12","0.12-0.24","0.24-0.36","0.36-0.48","0.48-0.60","0.60-0.72","0.72-0.84","0.84-1.00")
colnames(cumSum_PAS) = c("refSeq 3' UTR overlapping","ENSEMBL 3' UTR overlapping","Exon overlapping","Intron overlapping","Intergenic overlapping")
cumSum_PAS_melt = melt(cumSum_PAS)

q = ggplot(cumSum_PAS_melt,aes(x=X1,y=value,group=X2)) +geom_line(aes(col=X2),size=1) +theme_bw() + ylab("Cumulative fraction (number of priming sites)") + xlab("Downstream A content")  
q = q+ theme(legend.title=element_blank(),axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),legend.position=c(0.8,.4),legend.text=element_text(size=14), axis.ticks.length = unit(-0.25 , "cm"))  + theme(axis.text=element_text(size=12),axis.title=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
q = q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("noPAS containing") + geom_vline(xintercept = 2,linetype="dashed",color="red")
q = q  + annotate(geom = "text",x = 6,y = 0.1,label = paste0("Number of accepted UTR overlapping ends=", sum_refSeqEnsemblAccepted))
figs = list(p,q)
return(figs)
}
##### 
pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/cumulativePlots_numberOfsites.pdf",height = 5,width=9)
plotCululative(overlappers_annotations = overlappers_annotations_allData)
plotCululative(overlappers_annotations = overlappers_annotations_background)
dev.off()
#### signal overlappers

overlappers_signal_allData = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/signal_summary_allData.txt",stringsAsFactors = F)
overlappers_signal_background = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/signal_summary_background.txt",stringsAsFactors = F)


overlappers_signal_allData = overlappers_signal_allData[,c("L1","L2","refSeqOverlapping_signal", "ensemblOverlapping_signal", "exonOverlapping_signal", "intronOverlapping_signal", "nonOverlapping_signal")]
overlappers_signal_background = overlappers_signal_background[,c("L1","L2","refSeqOverlapping_signal", "ensemblOverlapping_signal", "exonOverlapping_signal", "intronOverlapping_signal", "nonOverlapping_signal")]

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/cumulativePlots_signal.pdf",height = 5,width=9)
plotCululative(overlappers_annotations = overlappers_signal_allData)
plotCululative(overlappers_annotations = overlappers_signal_background)
dev.off()


############################### comparing the number of ends retained after 90% thresholding 


thresholding_90_allData = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/ends_greater90percent_intergenic_n100_allData.bed",sep="\t",stringsAsFactors = F)
thresholding_90_background = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/ends_greater90percent_intergenic_n100_background.bed",sep="\t",stringsAsFactors = F)

ninetyPercentThreshold = list(thresholding_90_allData, thresholding_90_background)
names(ninetyPercentThreshold) = c("allData","background")

nrow_thresholding = melt(lapply(ninetyPercentThreshold,nrow))

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/highConfidenceThresholds.pdf")
p = ggplot(nrow_thresholding,aes(x=L1,y=value,fill=L1)) + geom_bar(stat="identity")+ scale_fill_manual(values = c("allData" = "blue","background" = "red"))+ theme(legend.title=element_blank(),axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),legend.position=c(0.8,.9),legend.text=element_text(size=14), axis.ticks.length = unit(-0.25 , "cm"))  + theme(axis.text=element_text(size=12),axis.title=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p = p + ylab("Number of high confidence ends")
print(p)
dev.off()

##### distribution of number of high confidence ends per gene


table_numberOfends = as.data.frame((table(thresholding_90_allData$V4)))
table_numberOfends$type = "alldata"
table_numberOfends_bg = as.data.frame((table(thresholding_90_background$V4)))
table_numberOfends_bg$type = "background"

totalEnds = rbind(table_numberOfends,table_numberOfends_bg)


## set working directory
getwd()

## call required packages
pkgs <- c("tidyverse", "ggthemes", "skimr")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

#### plotting the half violins

source("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/geom_violin.R")
library("ggplot2")
library("ggthemes")
install.packages("skimr")
library(skimr)

theme_set(theme_few())


pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/numberOfendsPergene.pdf",height = 5,width=6.5)
ggplot(data = totalEnds, 
       mapping = aes(x = type, y = Freq, fill = type)) + 
  geom_flat_violin(scale = "count", trim = FALSE) + 
  geom_dotplot(binaxis = "y", dotsize = 0.1, stackdir = "down", binwidth = 0.1, 
               position = position_nudge(-0.025)) + 
  theme(legend.position = "none") + 
  labs(x = " ", y = "Number of high confidence ends per gene") +  scale_fill_manual(values = c("alldata" = "blue","background" = "red"))
dev.off()


############for how many genes does all data capture ends but background doesnt

length(intersect(thresholding_90_allData$V4, thresholding_90_background$V4))


notCapturedInbackground = (setdiff(thresholding_90_allData$V4,thresholding_90_background$V4))
notCapturedInbackground =  thresholding_90_allData[ thresholding_90_allData$V4 %in% notCapturedInbackground,]
signal  = unlist(lapply(strsplit(notCapturedInbackground$V5,",",T),function(x) sum(as.numeric(x))))
notCapturedInbackground$signal = signal
signal = as.data.frame(signal)
signal$type = "nonIntersectingSignal" 
colnames(signal) = c("V1","V2")


capturedInBackground_alldata = (intersect(thresholding_90_allData$V4,thresholding_90_background$V4))
capturedInBackground_alldata = thresholding_90_allData[ thresholding_90_allData$V4 %in% capturedInBackground_alldata,]
signal_captured_alldata = unlist(lapply(strsplit(capturedInBackground_alldata$V5,",",T),function(x) sum(as.numeric(x))))
signal_captured_alldata = as.data.frame(signal_captured_alldata)
signal_captured_alldata$type = "intersecting_signal_allData"
colnames(signal_captured_alldata) = c("V1","V2")

capturedInBackground_background = (intersect(thresholding_90_allData$V4,thresholding_90_background$V4))
capturedInBackground_background = thresholding_90_background[ thresholding_90_background$V4 %in% capturedInBackground_background,]
signal_captured_background = unlist(lapply(strsplit(capturedInBackground_background$V5,",",T),function(x) sum(as.numeric(x))))
signal_captured_background = as.data.frame(signal_captured_background)
signal_captured_background$type = "intersecting_signal_backgorund"
colnames(signal_captured_background) = c("V1","V2")


allSignal  = rbind(signal,signal_captured_alldata,signal_captured_background)
allSignal$V2 = as.factor(allSignal$V2)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/signal_overlapping_nonOverlapping.pdf")
 p = ggplot(allSignal,aes(log10(V1),group=V2,col=V2,fill=V2,alpha=0.05)) + geom_density() + xlab("log10(polyA signal)")
 p = p  + theme(legend.title=element_blank(),axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),legend.position=c(0.7,.7),legend.text=element_text(size=14), axis.ticks.length = unit(-0.25 , "cm"))  + theme(axis.text=element_text(size=12),axis.title=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
 print(p)
dev.off()


############ now comparing the counting windows

countingWindows_allData = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/allAnnotations_allData.bed",stringsAsFactors = F,header=F)
countingWindows_background = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/allAnnotations_backgrounds.bed",stringsAsFactors = F,header = F)



contributedByensembl_allData = countingWindows_allData[-which((countingWindows_allData$V3 - countingWindows_allData$V2)==250),]
contributedByensembl_background = countingWindows_background[-which((countingWindows_background$V3 - countingWindows_background$V2)==250),]


countingWindows = list(contributedByensembl_allData, contributedByensembl_background)
names(countingWindows) = c("allData","background")

nrow_cw = melt(lapply(countingWindows,nrow))

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/numberOfCountingWindows.pdf")
p = ggplot(nrow_cw,aes(x=L1,y=value,fill=L1)) + geom_bar(stat="identity")+ scale_fill_manual(values = c("allData" = "blue","background" = "red"))+ theme(legend.title=element_blank(),axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),legend.position=c(0.8,.9),legend.text=element_text(size=14), axis.ticks.length = unit(-0.25 , "cm"))  + theme(axis.text=element_text(size=12),axis.title=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p = p + ylab("Number of counting windows") + xlab(' ')
print(p)
dev.off()

####### Overlapping the counting windows to check how many are exclusive to allData 
countingWindows_allData$id = paste0(countingWindows_allData$V1,countingWindows_allData$V2,countingWindows_allData$V3,countingWindows_allData$V4,countingWindows_allData$V6)
countingWindows_background$id = paste0(countingWindows_background$V1,countingWindows_background$V2,countingWindows_background$V3,countingWindows_background$V4,countingWindows_background$V6)

countingWindows_allData_granges = makeGRangesFromDataFrame(df = countingWindows_allData,keep.extra.columns = T,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)
countingWindows_background_granges = makeGRangesFromDataFrame(df = countingWindows_background,keep.extra.columns = T,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)

overlaps_allData_background = findOverlaps(query = countingWindows_allData_granges,subject = countingWindows_background_granges)

overlapping_countingWindows_allData  = countingWindows_allData[queryHits(overlaps_allData_background),]
overlapping_countingWindows_allData = overlapping_countingWindows_allData[!duplicated(overlapping_countingWindows_allData),]
nonOverlappint = setdiff(countingWindows_allData$id,overlapping_countingWindows_allData$id)
nonOverlapping_allData = countingWindows_allData[countingWindows_allData$id %in% nonOverlappint,]

overlappers_countingWindows_background = countingWindows_background[subjectHits(overlaps_allData_background),]
overlappers_countingWindows_background = overlappers_countingWindows_background[!duplicated(overlappers_countingWindows_background),]
nonOverlappint = setdiff(countingWindows_background$id,overlappers_countingWindows_background$id)
nonOverlapping_background = countingWindows_background[countingWindows_background$id %in% nonOverlappint,]



####### what is the singal contributing the non overlapping CWs


### nonoverlapping background CWs
nonOverlapping_background_granges = makeGRangesFromDataFrame(df = nonOverlapping_background,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)
ends_nonOverlapping_granges = makeGRangesFromDataFrame(df = thresholding_90_background,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)

overlap_ends_cw_background = findOverlaps(query = nonOverlapping_background_granges,subject = ends_nonOverlapping_granges)
endsCreatingBgNonOverlapping = thresholding_90_background[subjectHits(overlap_ends_cw_background),]


#### nonOverlapping allData CWs
nonOverlapping_all_granges = makeGRangesFromDataFrame(df = nonOverlapping_allData,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)
ends_nonOverlapping_granges = makeGRangesFromDataFrame(df = thresholding_90_allData,keep.extra.columns = T,ignore.strand = F,seqnames.field = "V1",start.field = "V2",end.field = "V3",strand.field = "V6",starts.in.df.are.0based = T)

overlap_ends_cw_allData = findOverlaps(query = nonOverlapping_all_granges,subject = ends_nonOverlapping_granges)
endsCreatingAllNonOverlapping = thresholding_90_allData[subjectHits(overlap_ends_cw_allData),]

nonOverlappinbg_signal = unlist(lapply(strsplit(endsCreatingBgNonOverlapping$V5,",",T),function(x) sum(as.numeric(x))))
nonOverlappinallData_signal = unlist(lapply(strsplit(endsCreatingAllNonOverlapping$V5,",",T),function(x) sum(as.numeric(x))))

nonOverlappinbg_signal = as.data.frame(nonOverlappinbg_signal)
nonOverlappinbg_signal$type = "background"
colnames(nonOverlappinbg_signal) = c("V1","V2")

nonOverlappinallData_signal = as.data.frame(nonOverlappinallData_signal)
nonOverlappinallData_signal$type = "alldata"
colnames(nonOverlappinallData_signal) = c("V1","V2")

allSignal_nonOverlappingCw = rbind(nonOverlappinbg_signal,nonOverlappinallData_signal)
pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/plots/signal_overlapping_nonOverlappingCountingWindows.pdf")
p = ggplot(allSignal_nonOverlappingCw,aes(log10(V1),group=V2,col=V2,fill=V2,alpha=0.05)) + geom_density() + xlab("log10(polyA signal)")
p = p  + theme(legend.title=element_blank(),axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),legend.position=c(0.7,.7),legend.text=element_text(size=14), axis.ticks.length = unit(-0.25 , "cm"))  + theme(axis.text=element_text(size=12),axis.title=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
print(p)


ggplot(data = allSignal_nonOverlappingCw, 
       mapping = aes(x = V2, y = log10(V1), fill = V2)) + 
  geom_flat_violin(scale = "count", trim = FALSE) + 
  geom_dotplot(binaxis = "y", dotsize = 0.1, stackdir = "down", binwidth = 0.1, 
               position = position_nudge(-0.025)) + 
  theme(legend.position = "none") + 
  labs(x = " ", y = "log10(polyA signal)") +  scale_fill_manual(values = c("alldata" = "blue","background" = "red")) 


dev.off()

##### compaing only the retained intergenic ends 

intergenicEnds_allData = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/onlyIntergenicEnds_alldata.bed",stringsAsFactors = F,sep="\t")
intergenicEnds_background = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/compareAllData_background/data/onlyIntergenicEnds_background.bed",stringsAsFactors = F,sep="\t")
