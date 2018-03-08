##### comparing genes from the different stages, that are accepted with noPAS from different stages... 


checkPAS = function(query_threshold)
{
  
  query_threshold$upstreamSequences = substr(query_threshold$V7,start = 20,stop = 60)
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


differentCategories = c("1dpf","2cell","2dpf","4dpf","256cell","bud","dome","oocyte","sphere","testis")

#### 

for(i in 1:length(differentCategories)){
  

samples_refSeq = read.table(paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/2016_12_QuantSeq_Annotaton/raw_data/",differentCategories[i],"//output_includingPAS_ensembl/polyAmapping_allTimepoints/n_100_global_a0/refSeq_overlapping.bed.gz"),stringsAsFactors = F,sep="\t")
samples_ensembl = read.table(paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/2016_12_QuantSeq_Annotaton/raw_data/",differentCategories[i],"//output_includingPAS_ensembl/polyAmapping_allTimepoints/n_100_global_a0/ensembl_overlapping.bed.gz"),stringsAsFactors = F,sep="\t")
totalSamples = rbind(samples_refSeq,samples_ensembl)




checkPAS_samples = checkPAS(query_threshold =totalSamples )
checkPAS_samples_noPAS = checkPAS_samples$noPAS
accepted_samplesNoPAS = checkPAS_samples_noPAS[which(checkPAS_samples_noPAS$V10 < 0.24),]
accepted_samplesNoPAS_1 = accepted_samplesNoPAS
accepted_samplesNoPAS_sequences = as.data.frame(accepted_samplesNoPAS$upstreamSequences)
accepted_samplesNoPAS = as.data.frame(unique(accepted_samplesNoPAS$V19))
accepted_samplesNoPAS = strsplit(as.character(accepted_samplesNoPAS$`unique(accepted_samplesNoPAS$V19)`),"_",T)

totalSamplePAS =  as.data.frame(unlist(lapply(accepted_samplesNoPAS,function(x) x[1])))
fileName = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/noPASaccepted_",differentCategories[i],".txt")
write.table(totalSamplePAS,fileName,quote=F,col.names = F,row.names = F)
fileName = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/noPASaccepted_sequences",differentCategories[i],".fasta")
write.table(accepted_samplesNoPAS_sequences,fileName,quote=F,col.names = F,row.names = F)

geneNames_accepted = unlist(accepted_samplesNoPAS)


}


#### for getting genes that have an AATAAA motif 

for(i in 1:length(differentCategories)){
  
  
  samples_refSeq = read.table(paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/2016_12_QuantSeq_Annotaton/raw_data/",differentCategories[i],"//output_includingPAS_ensembl/polyAmapping_allTimepoints/n_100_global_a0/refSeq_overlapping.bed.gz"),stringsAsFactors = F,sep="\t")
  samples_ensembl = read.table(paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/2016_12_QuantSeq_Annotaton/raw_data/",differentCategories[i],"//output_includingPAS_ensembl/polyAmapping_allTimepoints/n_100_global_a0/ensembl_overlapping.bed.gz"),stringsAsFactors = F,sep="\t")
  totalSamples = rbind(samples_refSeq,samples_ensembl)
  checkPAS_samples = checkPAS(query_threshold =totalSamples )
  checkPAS_samples_noPAS = checkPAS_samples$AATAAA
  accepted_samplesNoPAS = checkPAS_samples_noPAS[which(checkPAS_samples_noPAS$V10 < 0.36),]
  accepted_samplesNoPAS = as.data.frame(unique(accepted_samplesNoPAS$V19))
  accepted_samplesNoPAS = strsplit(as.character(accepted_samplesNoPAS$`unique(accepted_samplesNoPAS$V19)`),"_",T)
  totalSamplePAS =  as.data.frame(unlist(lapply(accepted_samplesNoPAS,function(x) x[1])))
  fileName = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/AATAAA_",differentCategories[i],".txt")
  write.table(totalSamplePAS,fileName,quote=F,col.names = F,row.names = F)
  
}


################## getting the noPAS sequences from the new SLAMseq datasets... 

refeqOverlapping = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/output/polyAmapping_allTimepoints/n_100_global_a0/refSeq_overlapping.bed.gz",stringsAsFactors = F)
EnsemblOverlapping = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/output/polyAmapping_allTimepoints/n_100_global_a0/ensembl_overlapping.bed.gz",sep="\t",stringsAsFactors = F)
totalSamples = rbind(refeqOverlapping,EnsemblOverlapping)
checkPAS_samples = checkPAS(query_threshold =totalSamples )
checkPAS_samples_noPAS = checkPAS_samples$noPAS
accepted_samplesNoPAS = checkPAS_samples_noPAS[which(checkPAS_samples_noPAS$V10 < 0.24),]
accepted_samplesNoPAS_sequences = as.data.frame(accepted_samplesNoPAS$upstreamSequences)
accepted_samplesNoPAS = as.data.frame(unique(accepted_samplesNoPAS$V19))
accepted_samplesNoPAS = strsplit(as.character(accepted_samplesNoPAS$`unique(accepted_samplesNoPAS$V19)`),"_",T)
write.table(accepted_samplesNoPAS_sequences,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/noPASacceptedSamples/geneLists_noPASaccepted/allNoPAS_SLAMdunkExperiment.fasta",quote = F,row.names = F,col.names = F)
