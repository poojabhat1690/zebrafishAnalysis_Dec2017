######## exonucleolytic deccay... 


#### reading in priming sites which map to the UTR from the different stages of development... 


############# the different stahes are : 

#stages = as.data.frame(c("1dpf", "2dpf", "4dpf" ,"256cell" ,"2cell" ,"bud", "dome", "oocyte" ,"sphere" ,"testis" ))
stages = as.data.frame(c("testis", "oocyte","2cell" ,"256cell" , "sphere","dome"  ,"bud","1dpf", "2dpf", "4dpf"  ))

path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/2016_12_QuantSeq_Annotaton/raw_data/"
stages$path = paste0(path,stages[,1],"/output_includingPAS_ensembl/polyAmapping_allTimepoints/n_100_global_a0/")
colnames(stages) = c("stage","path")

stages$refSeqOverlapping = paste0(stages$path,"/refSeq_overlapping.bed.gz")
stages$ensemblOverlapping = paste0(stages$path,"/ensembl_overlapping.bed.gz")


stages_ensembl = lapply(stages$ensemblOverlapping,function(x) read.table(x,stringsAsFactors = F,sep="\t"))
stages_refSeq = lapply(stages$refSeqOverlapping,function(x) read.table(x,stringsAsFactors = F,sep="\t"))

refSeqEnsembl = vector("list",length(stages_ensembl))
names(refSeqEnsembl) = stages$stage
for(i in 1:length(stages_ensembl)){
  refSeqEnsembl[[i]] = rbind(stages_ensembl[[i]],stages_refSeq[[i]])
}


splitGeneNames = lapply(refSeqEnsembl,function(x) strsplit(x$V19,"_",T))
splitGeneNames = lapply(splitGeneNames,function(x) lapply(x,function(y) y[[1]]))
splitGeneNames = lapply(splitGeneNames,function(x) unlist(x))

###### now I want to look at the number of priming sites per UTR per stage 
gettingAllTheGenes = as.data.frame(unique(do.call(c,lapply(refSeqEnsembl,function(x) x$V19))))
colnames(gettingAllTheGenes) = "x"

tabulatingGenes  = lapply(splitGeneNames,function(x) as.data.frame(table(x)))
tabulatingGenes = lapply(tabulatingGenes,function(x) plyr::join(gettingAllTheGenes,x))
tabulatingGenes_numbers = lapply(tabulatingGenes,function(x) x$Freq)
tabulatingGenes_numbers = do.call(cbind.data.frame,tabulatingGenes_numbers)
row.names(tabulatingGenes_numbers) = tabulatingGenes$`1dpf`$x


##### overall fraction relative to testis 

allGenes_withRespectToTestis = tabulatingGenes_numbers/tabulatingGenes_numbers$testis
allGenes_withRespectToTestis = allGenes_withRespectToTestis[rowSums(is.na(allGenes_withRespectToTestis))!=ncol(allGenes_withRespectToTestis), ]
allGenes_withRespectToTestis_melt = melt(allGenes_withRespectToTestis)

pdf("/Volumes/groups/ameres/Pooja/exo.pdf")
ggplot(allGenes_withRespectToTestis_melt,aes(y=value,x=variable,group=variable)) + geom_boxplot() + ggtitle("allGenes")+ ylim(c(0,10))


###### noq looking only at the genes that have a high u content in oocyte 

oocyteSpecific = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/noPASaccepted_oocyte.txt")
oocyteSpecific_wrtTestis = allGenes_withRespectToTestis[row.names(allGenes_withRespectToTestis) %in% oocyteSpecific$V1,]
oocyteSpecific_wrtTestis_melt = melt(oocyteSpecific_wrtTestis)
ggplot(oocyteSpecific_wrtTestis_melt,aes(x=variable,y=value,group=variable)) + geom_boxplot()+ ggtitle("oocyte_ucontaining")+ ylim(c(0,10))

### 2 cell 


###### noq looking only at the genes that have a high u content in oocyte 

cell2Specific = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/noPASaccepted_2cell.txt")
cell2Specific_wrtTestis = allGenes_withRespectToTestis[row.names(allGenes_withRespectToTestis) %in% cell2Specific$V1,]
cell2Specific_wrtTestis_melt = melt(cell2Specific_wrtTestis)
ggplot(cell2Specific_wrtTestis_melt,aes(x=variable,y=value,group=variable)) + geom_boxplot()+ ggtitle("2cell_ucontaining")+ ylim(c(0,10))


### bud 

budSpecific = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/noPASaccepted_bud.txt")
budSpecific_wrtTestis = allGenes_withRespectToTestis[row.names(allGenes_withRespectToTestis) %in% budSpecific$V1,]
budSpecific_wrtTestis_melt = melt(budSpecific_wrtTestis)
ggplot(budSpecific_wrtTestis_melt,aes(x=variable,y=value,group=variable)) + geom_boxplot() + ggtitle("bud_ucontaining")+ ylim(c(0,10))

### dome
domeSpecific = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/noPASaccepted_dome.txt")
domeSpecific_wrtTestis = allGenes_withRespectToTestis[row.names(allGenes_withRespectToTestis) %in% domeSpecific$V1,]
domeSpecific_wrtTestis_melt = melt(domeSpecific_wrtTestis)
ggplot(domeSpecific_wrtTestis_melt,aes(x=variable,y=value,group=variable)) + geom_boxplot() + ggtitle("dome_ucontaining")+ ylim(c(0,10))

### dome
sphereSpecific = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/noPASaccepted_sphere.txt")
sphereSpecific_wrtTestis = allGenes_withRespectToTestis[row.names(allGenes_withRespectToTestis) %in% sphereSpecific$V1,]
sphereSpecific_wrtTestis_melt = melt(sphereSpecific_wrtTestis)
ggplot(sphereSpecific_wrtTestis_melt,aes(x=variable,y=value,group=variable)) + geom_boxplot() + ggtitle("sphere_ucontaining")+ ylim(c(0,10))


##### 1dpf specific 

dpf1Specific = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/noPASaccepted_1dpf.txt")
dpf1Specific_wrtTestis = allGenes_withRespectToTestis[row.names(allGenes_withRespectToTestis) %in% dpf1Specific$V1,]
dpf1Specific_wrtTestis_melt = melt(dpf1Specific_wrtTestis)
ggplot(dpf1Specific_wrtTestis_melt,aes(x=variable,y=value,group=variable)) + geom_boxplot() + ggtitle("1dpf_ucontaining")+ ylim(c(0,10))

### 2dpf 
dpf2Specific = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/noPASaccepted_2dpf.txt")
dpf2Specific_wrtTestis = allGenes_withRespectToTestis[row.names(allGenes_withRespectToTestis) %in% dpf2Specific$V1,]
dpf2Specific_wrtTestis_melt = melt(dpf2Specific_wrtTestis)
ggplot(dpf2Specific_wrtTestis_melt,aes(x=variable,y=value,group=variable)) + geom_boxplot() + ggtitle("2dpf_ucontaining")+ ylim(c(0,10))


##### 4dpf specific 

dpf4Specific = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/noPASaccepted_4dpf.txt")
dpf4Specific_wrtTestis = allGenes_withRespectToTestis[row.names(allGenes_withRespectToTestis) %in% dpf4Specific$V1,]
dpf4Specific_wrtTestis_melt = melt(dpf4Specific_wrtTestis)
ggplot(dpf4Specific_wrtTestis_melt,aes(x=variable,y=value,group=variable)) + geom_boxplot()+ ggtitle("4dpf_ucontaining") + ylim(c(0,10))


cell256Specific = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/noPASaccepted_256cell.txt")
cell256Specific_wrtTestis = allGenes_withRespectToTestis[row.names(allGenes_withRespectToTestis) %in% cell256Specific$V1,]
cell256Specific_wrtTestis_melt = melt(cell256Specific_wrtTestis)
ggplot(cell256Specific_wrtTestis_melt,aes(x=variable,y=value,group=variable)) + geom_boxplot()+ ggtitle("256_ucontaining") + ylim(c(0,10))


dev.off()


####
cell256Specific_wrtTestis_melt$type = "U genes 256cell"
dpf2Specific_wrtTestis_melt$type = "U genes 2dpf"
budSpecific_wrtTestis_melt$type = "U genes bud"
domeSpecific_wrtTestis_melt$type = "U genes dome"
sphereSpecific_wrtTestis_melt$type = "U genes sphere"
dpf4Specific_wrtTestis_melt$type = "U genes 4dpf"
dpf1Specific_wrtTestis_melt$type = "U genes 1dpf"
cell2Specific_wrtTestis_melt$type = "U genes 2cellStage"
oocyteSpecific_wrtTestis_melt$type = "U genes Oocyte"
allGenes_withRespectToTestis_melt$type = "allGenes"

total = rbind(dpf2Specific_wrtTestis_melt,budSpecific_wrtTestis_melt,domeSpecific_wrtTestis_melt,sphereSpecific_wrtTestis_melt,dpf4Specific_wrtTestis_melt,dpf1Specific_wrtTestis_melt,cell2Specific_wrtTestis_melt,oocyteSpecific_wrtTestis_melt,allGenes_withRespectToTestis_melt)
ggplot(total,aes(x=variable,y=value,fill=type)) + geom_boxplot()

cell256Specific_wrtTestis_melt_median  = cell256Specific_wrtTestis_melt %>% group_by(variable) %>% mutate(median = mean(value,na.rm=T))
median_256cell = cell256Specific_wrtTestis_melt_median[!duplicated(cell256Specific_wrtTestis_melt_median[,c("variable","type","median")]),]

dpf2Specific_wrtTestis_melt_median  = dpf2Specific_wrtTestis_melt %>% group_by(variable) %>% mutate(median = mean(value,na.rm=T))
median_2dpf = dpf2Specific_wrtTestis_melt_median[!duplicated(dpf2Specific_wrtTestis_melt_median[,c("variable","type","median")]),]

budSpecific_wrtTestis_melt_median  = budSpecific_wrtTestis_melt %>% group_by(variable) %>% mutate(median = mean(value,na.rm=T))
median_bud = budSpecific_wrtTestis_melt_median[!duplicated(budSpecific_wrtTestis_melt_median[,c("variable","type","median")]),]

domeSpecific_wrtTestis_melt_median  = domeSpecific_wrtTestis_melt %>% group_by(variable) %>% mutate(median = mean(value,na.rm=T))
median_dome = domeSpecific_wrtTestis_melt_median[!duplicated(domeSpecific_wrtTestis_melt_median[,c("variable","type","median")]),]

sphereSpecific_wrtTestis_melt_median  = sphereSpecific_wrtTestis_melt %>% group_by(variable) %>% mutate(median = mean(value,na.rm=T))
median_sphere = sphereSpecific_wrtTestis_melt_median[!duplicated(sphereSpecific_wrtTestis_melt_median[,c("variable","type","median")]),]

dpf4Specific_wrtTestis_melt_median  = dpf4Specific_wrtTestis_melt %>% group_by(variable) %>% mutate(median = mean(value,na.rm=T))
median_4dpf = dpf4Specific_wrtTestis_melt_median[!duplicated(dpf4Specific_wrtTestis_melt_median[,c("variable","type","median")]),]


dpf1Specific_wrtTestis_melt_median  = dpf1Specific_wrtTestis_melt %>% group_by(variable) %>% mutate(median = mean(value,na.rm=T))
median_1dpf = dpf1Specific_wrtTestis_melt_median[!duplicated(dpf1Specific_wrtTestis_melt_median[,c("variable","type","median")]),]

cell2Specific_wrtTestis_melt_median  = cell2Specific_wrtTestis_melt %>% group_by(variable) %>% mutate(median = mean(value,na.rm=T))
median_2cell = cell2Specific_wrtTestis_melt_median[!duplicated(cell2Specific_wrtTestis_melt_median[,c("variable","type","median")]),]

oocyteSpecific_wrtTestis_melt_median  = oocyteSpecific_wrtTestis_melt %>% group_by(variable) %>% mutate(median = mean(value,na.rm=T))
median_oocyte = oocyteSpecific_wrtTestis_melt_median[!duplicated(oocyteSpecific_wrtTestis_melt_median[,c("variable","type","median")]),]

allGenes_withRespectToTestis_melt_median  = allGenes_withRespectToTestis_melt %>% group_by(variable) %>% mutate(median = mean(value,na.rm=T))
median_allGenes = allGenes_withRespectToTestis_melt_median[!duplicated(allGenes_withRespectToTestis_melt_median[,c("variable","type","median")]),]


### plotting the number of sites being considered in analysis : 
nrow(cell256Specific_wrtTestis)
nrow(dpf2Specific_wrtTestis)
nrow(budSpecific_wrtTestis)
nrow(domeSpecific_wrtTestis)
nrow(sphereSpecific_wrtTestis)
nrow(dpf4Specific_wrtTestis)
nrow(dpf1Specific_wrtTestis)
nrow(cell2Specific_wrtTestis)
nrow(oocyteSpecific_wrtTestis)
nrow(allGenes_withRespectToTestis)



### plottingOnlt Medians
allMedians = rbind(median_allGenes,median_oocyte,median_2cell,median_1dpf,median_4dpf,median_2dpf,median_bud,median_dome,median_sphere,median_256cell)
pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/plot/fractionOfSites_ucontainingGenes.pdf",height=4,width=6)
ggplot(allMedians,aes(x=variable,y=median,group=type,col=type)) + geom_line(size=1) + ylab("number sites/number sites(testis)") + xlab(" ")
dev.off()

#### plotting the number of priming sites per stage

numberOfPrimingSites = melt(lapply(refSeqEnsembl,nrow))

### for a fair comparison I also need to compare the number of ends that actually pass the A thershold 




checkPAS = function(query_threshold)
{
  
  query_threshold$upstreamSequences = substr(query_threshold$V7,start = 20,stop = 55)
  query_threshold$upstreamSequences = toupper(query_threshold$upstreamSequences)
  #toMatch <- c("TATAAA", "AGTAAA", "AATACA","CATAAA","AATATA","GATAAA","AATGAA","AAGAAA","ACTAAA","AATAGA","AATAAT","AACAAA","ATTACA","ATTATA","AACAAG","AATAAG")
  # using the motifs identified in ulitsky et al. 2012 ... with the motifs from the mouse data
  toMatch = c("TATAAA", "AGTAAA", "TTTAAA", "CATAAA", "AATACA", "AATGAA", "AATATA", "GATAAA", "TGTAAA","AATACA","AAGAAA","ACTAAA","AATAGA","AATAAT","AACAAA","ATTATA","AACAAG","AATAAG")
  
  query_AATAAA = query_threshold[grepl("AATAAA",query_threshold$upstreamSequences),]
  query_threshold = query_threshold[!grepl("AATAAA",query_threshold$upstreamSequences),]
  query_ATTAAA = query_threshold[grepl("ATTAAA",query_threshold$upstreamSequences),]
  query_threshold = query_threshold[!grepl("ATTAAA",query_threshold$upstreamSequences),]
  query_APA = query_threshold[grepl(paste(toMatch,collapse="|"),query_threshold$upstreamSequences),]
  query_minusAPA= query_threshold[!grepl(paste(toMatch,collapse="|"),query_threshold$upstreamSequences),]
  queryStats = list(query_AATAAA,query_ATTAAA, query_APA, query_minusAPA)
  names(queryStats) = c("AATAAA","ATTAAA","APA","noPAS")
  return(queryStats)
  
}

noPASaccepted = function(checked){
  checked_noPAS = checked$noPAS
  checked_noPAS_accepted = checked_noPAS[which(checked_noPAS$V10 < 0.24),]
  return(checked_noPAS_accepted)
}



PASaccepted = function(checked){
  checked_PAS = rbind(checked$AATAAA,checked$ATTAAA,checked$APA)
  checked_PAS_accepted = checked_PAS[which(checked_PAS$V10 < 0.36),]
  return(checked_PAS_accepted)
}

PASaccepted_stages = lapply(refSeqEnsembl,function(x) PASaccepted(checkPAS(x)))
noPASaccepted_stages = lapply(refSeqEnsembl,function(x) noPASaccepted(checkPAS(x)))

PASAcceptedSites = melt(lapply(PASaccepted_stages,nrow))
PASAcceptedSites$L2 = "PAS"
noPASacceptedSites = melt(lapply(noPASaccepted_stages,nrow))
noPASacceptedSites$L2 = "noPAS"
allAcceptedPAS_noPAS = rbind(PASAcceptedSites,noPASacceptedSites)
allAcceptedPAS_noPAS$L1 <- factor(allAcceptedPAS_noPAS$L1, levels = c("testis","oocyte","2cell","256cell","sphere","dome","bud","1dpf","2dpf","4dpf"))

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/plot/numberOfAcceptedPrimingSites.pdf",height=6,width=10)
p = ggplot(allAcceptedPAS_noPAS,aes(x=L1,y=value,group=L2,col=L2,fill=L2)) + geom_bar(stat = "identity",position = "dodge")  + ylab("Number of priming sites") + xlab("stages") + theme_bw()
p +  theme(axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),axis.title  = element_text(size=14),legend.text=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length = unit(-0.25 , "cm"),legend.position=c(0.8,0.8))  
dev.off()

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/plot/numberOfPrimingSites.pdf")
ggplot(numberOfPrimingSites,aes(x=L1,y=value)) + geom_bar(stat = "identity")
dev.off()

