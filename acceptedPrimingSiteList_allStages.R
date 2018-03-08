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


noPASaccepted_stages = lapply(noPASaccepted_stages,function(x) x[,c(1:6)])
PASaccepted_stages = lapply(PASaccepted_stages,function(x) x[,c(1:6)])
 
refSeqEnsembl_checkPAS = lapply(refSeqEnsembl,function(x) checkPAS(x))
refSeqEnsembl_checkPAS_noPAS = lapply(refSeqEnsembl_checkPAS,function(x) x$noPAS)
refSeqEnsembl_checkPAS_PAS = lapply(refSeqEnsembl_checkPAS,function(x) rbind.data.frame(x$AATAAA,x$ATTAAA,x$APA))


for(i in 1:length(noPASaccepted_stages)){
  
  write.table(noPASaccepted_stages[[i]],paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/forIGV/noPASaccepted_",names(noPASaccepted_stages)[i],".bed"),quote = F,col.names = F,row.names = F,sep="\t")
  write.table(PASaccepted_stages[[i]],paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/forIGV/PASaccepted_",names(PASaccepted_stages)[i],".bed"),quote = F,col.names = F,row.names = F,sep="\t")
  write.table(refSeqEnsembl_checkPAS_noPAS[[i]],paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/forIGV/noPASALL_",names(refSeqEnsembl_checkPAS_noPAS)[i],".bed"),quote = F,col.names = F,row.names = F,sep="\t")
  write.table(refSeqEnsembl_checkPAS_PAS[[i]],paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/forIGV/PASALL_",names(refSeqEnsembl_checkPAS_PAS)[i],".bed"),quote = F,col.names = F,row.names = F,sep="\t")
  
}