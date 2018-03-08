##### U-rich ends : 
  ### are they the only end of the gene?
  ### are they proximal or distal?


###########################################################################

####### Do genes with U-rich ends, show PAS ends in addition? (consider all ends that pass the A threshold )
      ##### doing this for all (new annotation)
      ##### stage specific annotation
###########################################################################

##### readin in all the refSeq and ensembl overlapping ends 

allData_refSeq = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/output/polyAmapping_allTimepoints/n_100_global_a0/refSeq_overlapping.bed.gz",stringsAsFactors = F,sep="\t")
allData_ensembl = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/output/polyAmapping_allTimepoints/n_100_global_a0/ensembl_overlapping.bed.gz",stringsAsFactors = F,sep="\t")



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


noPASaccepted_stages = lapply(refSeqEnsembl,function(x) noPASaccepted(checkPAS(x)))
PASaccepted_stages = lapply(refSeqEnsembl,function(x) PASaccepted(checkPAS(x)))


##### now I want to split the names of each of the things 

splitNames = function(dataframe){
  
  splitGeneNames = strsplit(dataframe$V19,"_",T)
  splitGeneNames = lapply(splitGeneNames,function(x) x[[1]])
  splitGeneNames = unlist(splitGeneNames)
  
  dataframe$name = splitGeneNames
  return(dataframe)
}


noPASaccepted_stages = lapply(noPASaccepted_stages,function(x) splitNames(x))
PASaccepted_stages = lapply(PASaccepted_stages,function(x) splitNames(x))


####### DIVIDING ALL GENES INTO CLASSES 


##### writing a function for this 

getNumberOfPrimingSites = function(PAScontaining,noPAScontaining,name_stage){
  BothPASandNOPAS = intersect(noPAScontaining$name,PAScontaining$name)
  BothPASandNOPAS_list = vector("list",2)
  names(BothPASandNOPAS_list) = c("inPAS","inNoPAS")
  BothPASandNOPAS_list$inPAS =PAScontaining[PAScontaining$name %in% BothPASandNOPAS,]
  BothPASandNOPAS_list$inNoPAS =noPAScontaining[noPAScontaining$name %in% BothPASandNOPAS,]
  OnlyNoPAS = setdiff(noPAScontaining$name,PAScontaining$name)
  OnlyNoPAS_genes = noPAScontaining[noPAScontaining$name %in% OnlyNoPAS,]
  OnlyNoPAS_greaterThan1  = names(which(table(OnlyNoPAS_genes$name)>1))
  OnlyNoPAS_genes = OnlyNoPAS_genes[OnlyNoPAS_genes$name %in% OnlyNoPAS_greaterThan1,]
  
  OnlyPAS = setdiff(PAScontaining$name,noPAScontaining$name)
  OnlyPAS_genes = PAScontaining[PAScontaining$name %in% OnlyPAS,]
  OnlyPAS_greaterThan1  = names(which(table(OnlyPAS_genes$name)>1))
  OnlyPAS_genes = OnlyPAS_genes[OnlyPAS_genes$name %in% OnlyPAS_greaterThan1,]
  
  
  allAccepted_PAS_noPAS = do.call(rbind,BothPASandNOPAS_list)
  numberPrimingSites = as.data.frame(table(allAccepted_PAS_noPAS$name))
  numberPrimingSites$type = "PAS_noPAS"
  OnlynoPAS_table = as.data.frame(table(OnlyNoPAS_genes$name))
  OnlynoPAS_table$type = "onlyNoPAS"
  OnlyPAS_table = as.data.frame(table(OnlyPAS_genes$name))
  OnlyPAS_table$type = "onlyPAS"
  
  
  allPrimingSites = rbind.data.frame(numberPrimingSites,OnlynoPAS_table,OnlyPAS_table)
p =   ggplot(allPrimingSites,aes(x=type,y=Freq,group=type,col=type)) + geom_violin() + ggtitle(name_stage)
returnthis = vector("list",2)
names(returnthis) = c("allSites","plot")
returnthis[[1]] = allPrimingSites
returnthis[[2]] = p
return(returnthis)
}


pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/plot/numberOfPrimingSites_PASnoPAS.pdf")
oocytePlot = getNumberOfPrimingSites(PAScontaining = PASaccepted_stages$oocyte,noPAScontaining = noPASaccepted_stages$oocyte,"oocyte")
print(oocytePlot[[2]])

testisPlot = getNumberOfPrimingSites(PAScontaining = PASaccepted_stages$testis,noPAScontaining = noPASaccepted_stages$testis,"testis")
print(testisPlot[[2]])

cell2Plot = getNumberOfPrimingSites(PAScontaining = PASaccepted_stages$`2cell`,noPAScontaining = noPASaccepted_stages$`2cell`,"2cell")
print(cell2Plot[[2]])

cell256Plot = getNumberOfPrimingSites(PAScontaining = PASaccepted_stages$`256cell`,noPAScontaining = noPASaccepted_stages$`256cell`,"256cell")
print(cell256Plot[[2]])

spherePlot = getNumberOfPrimingSites(PAScontaining = PASaccepted_stages$sphere,noPAScontaining = noPASaccepted_stages$sphere,"sphere")
print(spherePlot[[2]])

domePlot = getNumberOfPrimingSites(PAScontaining = PASaccepted_stages$dome,noPAScontaining = noPASaccepted_stages$dome,"dome")
print(domePlot[[2]])


budPlot = getNumberOfPrimingSites(PAScontaining = PASaccepted_stages$bud,noPAScontaining = noPASaccepted_stages$bud,"bud")
print(budPlot[[2]])


dpf1Plot = getNumberOfPrimingSites(PAScontaining = PASaccepted_stages$`1dpf`,noPAScontaining = noPASaccepted_stages$`1dpf`,"1dpf")
print(dpf1Plot[[2]])

dpf2Plot = getNumberOfPrimingSites(PAScontaining = PASaccepted_stages$`2dpf`,noPAScontaining = noPASaccepted_stages$`2dpf`,"2dpf")
print(dpf2Plot[[2]])

dpf4Plot = getNumberOfPrimingSites(PAScontaining = PASaccepted_stages$`4dpf`,noPAScontaining = noPASaccepted_stages$`4dpf`,"4dpf")
print(dpf4Plot[[2]])

dev.off()


##### comparing only PAS, onlynoPAs and PAS+noPAS from each stage 

allSitesAndStages = list(oocytePlot[[1]],testisPlot[[1]],cell2Plot[[1]],cell256Plot[[1]],spherePlot[[1]],domePlot[[1]],budPlot[[1]],dpf1Plot[[1]],dpf2Plot[[1]],dpf4Plot[[1]])
names(allSitesAndStages) = c("oocyte","testis","2cell","256cell","sphere","dome","bud","dpf1","dpf2","dpf4")
allSitesAndStages_melt = melt(allSitesAndStages)
allSitesAndStages_melt$L1 <- factor(allSitesAndStages_melt$L1, levels = c("testis","oocyte","2cell","256cell","sphere","dome","bud","dpf1","dpf2","dpf4"))

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/plot/numberOfPrimingSites_allStages.pdf",height=5,width=22)
p = ggplot(allSitesAndStages_melt,aes(x=L1,y=value,group=L1,col=L1)) + geom_violin() + facet_wrap(~type) + ylim(c(0,20)) + ylab("Distribution of number of priming sites per gene") + theme_bw()
p = p + theme(axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),axis.title  = element_text(size=14),legend.text=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length = unit(-0.25 , "cm"),legend.position="none")  
p
dev.off()

##### ###### now I have genes where there is a PAS site and a noPAS site, so are these proximal or distal ?



getNOPAS_PASsites = function(PAScontaining,noPAScontaining){
  BothPASandNOPAS = intersect(noPAScontaining$name,PAScontaining$name)
  BothPASandNOPAS_list = vector("list",2)
  names(BothPASandNOPAS_list) = c("inPAS","inNoPAS")
  BothPASandNOPAS_list$inPAS =PAScontaining[PAScontaining$name %in% BothPASandNOPAS,]
  BothPASandNOPAS_list$inNoPAS =noPAScontaining[noPAScontaining$name %in% BothPASandNOPAS,]
  return(BothPASandNOPAS_list)
}


PASnoPAS_oocyte  = getNOPAS_PASsites(PAScontaining =PASaccepted_stages$oocyte,noPAScontaining = noPASaccepted_stages$oocyte )

PASnoPAS_testis  = getNOPAS_PASsites(PAScontaining =PASaccepted_stages$testis,noPAScontaining = noPASaccepted_stages$testis )

PASnoPAS_2cell  = getNOPAS_PASsites(PAScontaining =PASaccepted_stages$`2cell`,noPAScontaining = noPASaccepted_stages$`2cell` )
PASnoPAS_256cell  = getNOPAS_PASsites(PAScontaining =PASaccepted_stages$`256cell`,noPAScontaining = noPASaccepted_stages$`256cell` )

PASnoPAS_sphere  = getNOPAS_PASsites(PAScontaining =PASaccepted_stages$sphere,noPAScontaining = noPASaccepted_stages$sphere )

PASnoPAS_dome  = getNOPAS_PASsites(PAScontaining =PASaccepted_stages$dome,noPAScontaining = noPASaccepted_stages$dome )
PASnoPAS_bud  = getNOPAS_PASsites(PAScontaining =PASaccepted_stages$bud,noPAScontaining = noPASaccepted_stages$bud )
PASnoPAS_1dpf  = getNOPAS_PASsites(PAScontaining =PASaccepted_stages$`1dpf`,noPAScontaining = noPASaccepted_stages$`1dpf` )
PASnoPAS_2dpf  = getNOPAS_PASsites(PAScontaining =PASaccepted_stages$`2dpf`,noPAScontaining = noPASaccepted_stages$`2dpf` )
PASnoPAS_4dpf  = getNOPAS_PASsites(PAScontaining =PASaccepted_stages$`4dpf`,noPAScontaining = noPASaccepted_stages$`4dpf` )



##### PAS containing sites contribute to the increased number of priming sites 

PASnoPAS_allStages = list(PASnoPAS_oocyte,PASnoPAS_testis,PASnoPAS_2cell,PASnoPAS_256cell,PASnoPAS_sphere,PASnoPAS_dome,PASnoPAS_bud,PASnoPAS_1dpf,PASnoPAS_2dpf,PASnoPAS_4dpf)
names(PASnoPAS_allStages) = c("oocyte","testis","2cell","256cell","sphere","dome","bud","1dpf","2dpf","4dpf")

PASnoPAS_allStages_NumberOfPrimingSites =  lapply(PASnoPAS_allStages,function(x) lapply(x,nrow))
PASnoPAS_allStages_NumberOfPrimingSites_melt = melt(PASnoPAS_allStages_NumberOfPrimingSites)
PASnoPAS_allStages_NumberOfPrimingSites_melt$L1 <- factor(PASnoPAS_allStages_NumberOfPrimingSites_melt$L1, levels = c("testis","oocyte","2cell","256cell","sphere","dome","bud","1dpf","2dpf","4dpf"))

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/plot/numberOfprimingSites_containingPAS_noPAS.pdf",height=6,width=10)
p = ggplot(PASnoPAS_allStages_NumberOfPrimingSites_melt,aes(x=L1,y=value,group=L2,col=L2,fill=L2)) + geom_bar(stat = "identity",position = "dodge")  + ylab("Number of priming sites") + xlab("stages") + theme_bw()
p +  theme(axis.text.x = element_text(size=14,margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(size=14,margin=margin(5,15,10,5,"pt")),axis.title  = element_text(size=14),legend.text=element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.length = unit(-0.25 , "cm"),legend.position=c(0.8,0.8))  
dev.off()

########## now there are more priming sites that are contributed by the 'PAS ends', lets look if the noPAS sites is proximal or distal... 

##### i want to check the fraction of PAS containing sites that are proximal to the noPAS site and the fraction of PAS containing sites that are distal to the noPAS site. 


  #### so I need to calculate : 
    ### total number of PAS containing sites for each gene. 
    ### for each noPAS site, what is the fraction of proximal PAS sites
    ### for each noPAS site, what is the fraction of proximal PAS sites

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/plot/distal_proximalBias.pdf")
for(k in 1:length(PASnoPAS_allStages)){
  

noPAS_positive = PASnoPAS_allStages[[k]]$inNoPAS %>% filter(V6 == "+")
noPAS_negative = PASnoPAS_allStages[[k]]$inNoPAS %>% filter(V6 == "-")

distances_PAS_noPAS_positive = vector("list",nrow(noPAS_positive))

for(i in 1:nrow(noPAS_positive)){
  a =noPAS_positive[i,] 
  consider_subset = PASnoPAS_allStages[[k]]$inPAS[which( PASnoPAS_allStages[[k]]$inPAS$name == a$name),]
  distances_a = consider_subset$V2 - a$V2
  distances_PAS_noPAS_positive[[i]] = distances_a
}

distances_PAS_noPAS_positive_proximal = unlist(lapply(distances_PAS_noPAS_positive,function(x) length(which(x<0))))
distances_PAS_noPAS_positive_distal = unlist(lapply(distances_PAS_noPAS_positive,function(x) length(which(x>0))))

totalPASsites_positive = distances_PAS_noPAS_positive_proximal + distances_PAS_noPAS_positive_distal



proximal_distal_positive = cbind.data.frame(distances_PAS_noPAS_positive_proximal,distances_PAS_noPAS_positive_distal)
colnames(proximal_distal_positive) = c("proximal","distal")
#proximal_distal_positive = proximal_distal_positive/(proximal_distal_positive$proximal + proximal_distal_positive$distal)
proximal_distal_positive = melt(proximal_distal_positive)

distances_PAS_noPAS_negative = vector("list",nrow(noPAS_negative))

for(i in 1:nrow(noPAS_negative)){
  a = noPAS_negative[i,] 
  consider_subset = PASnoPAS_allStages[[k]]$inPAS[which( PASnoPAS_allStages[[k]]$inPAS$name == a$name),]
  distances_a = consider_subset$V2 - a$V2
  distances_PAS_noPAS_negative[[i]] = distances_a
}

distances_PAS_noPAS_negative_proximal = unlist(lapply(distances_PAS_noPAS_negative,function(x) length(which(x>0))))
distances_PAS_noPAS_negative_distal = unlist(lapply(distances_PAS_noPAS_negative,function(x) length(which(x<0))))

totalPASsites_negative = distances_PAS_noPAS_negative_proximal + distances_PAS_noPAS_negative_distal
proximal_distal_negative = cbind.data.frame(distances_PAS_noPAS_negative_proximal,distances_PAS_noPAS_negative_distal)
colnames(proximal_distal_negative) = c("proximal","distal")
#proximal_distal_negative = proximal_distal_negative/(proximal_distal_negative$proximal + proximal_distal_negative$distal)
proximal_distal_negative = melt(proximal_distal_negative)



proximal_distal_bothStrands = rbind.data.frame(proximal_distal_positive,proximal_distal_negative)
p = ggplot(proximal_distal_bothStrands,aes(value,group=variable,col=variable))  + geom_density() + ggtitle(names(PASaccepted_stages)[k]) + ylim(c(0,2.5))
print(p)
}
dev.off()