##### for thw QuantSeq data 


sampleInfo = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/backgroundSamples/data/summary_ucscBrowser.txt",header = F,stringsAsFactors = F)
sampleInfo$V17 = paste0(sampleInfo$V17,c(1:nrow(sampleInfo)))
### main header

for(i in 1:nrow(sampleInfo)){
  trackName = paste("track",sampleInfo[i,17])
  container =  "container multiWig"
  shortLab = paste("shortLabel", paste0("Sample",sampleInfo[i,17],"RPM"),sep=" ")
  longLab = paste("longLabel", paste0("Sample",sampleInfo[i,17],"RPM"),sep=" ")
  aggregate = "aggregate transparentOverlay"
  type = "type bigWig"
  
  autoSc = "autoscale on"
  visibility ="visibility full"
  colorUse = "color 0,0,130"
  parent = paste("parent",sampleInfo[i,17])
  
  plusTrack = paste0("http://amereslab.imba.oeaw.ac.at/seq-hub/dr_slamSeq/dr10/",sampleInfo[i,16],"_plus.bg.bigWig",sep="")
  minusTrack = paste0("http://amereslab.imba.oeaw.ac.at/seq-hub/dr_slamSeq/dr10/",sampleInfo[i,16],"_neg_minus.bg.bigWig")
  
  ### plus strand
  
  track_plus = paste("track",paste(sampleInfo[i,17],"plus",sep="_"))
  bigDataUrl_plus = paste("bigDataUrl",plusTrack)
  
  ##minusStrand
  
  track_minus = paste("track",paste(sampleInfo[i,17],"minus",sep="_"))
  bigDataUrl_minus = paste("bigDataUrl",minusTrack)
  
  
  ### main format
  
  main  = paste(trackName,container,shortLab,longLab,aggregate,type,sep="\n")
  main = paste(main,"\n")
  plusStrand = paste(track_plus,parent,type,bigDataUrl_plus,autoSc,visibility,colorUse,sep="\n")
  plusStrand = paste(plusStrand,"\n")
  minusStrand = paste(track_minus,parent,type,bigDataUrl_minus,autoSc,visibility,colorUse,sep="\n")
  
  main = paste(main,plusStrand,minusStrand,sep = "\n")
  
  write.table(main,paste0("/Volumes/groups/ameres/wwwdata/seq-hub/dr_slamSeq/dr10/",sampleInfo[i,17],".txt"),quote = F,row.names = F,col.names = F)
  
}



includeFiles1 = paste0(paste("include",sampleInfo[,17],sep=" "),".txt")
#write.table(includeFiles,"/Volumes/groups/ameres/wwwdata/seq-hub/dr_slamSeq/dr10/include_quantSeq_backgroundCWs.txt",quote = F,col.names = F,row.names = F)

###### writing similar file for the bed files for high confidence ends 

bbFiles_ends = list.files (path = "/Volumes/groups/ameres/wwwdata/seq-hub/dr_slamSeq/dr10/",pattern = "ends_")
names_timepoints = unlist(lapply(strsplit(bbFiles_ends,"_",T),function(x) x[[3]]))

for(i in 1:length(bbFiles_ends)){
  trackName = paste("track",names_timepoints[i])
  type = "type bigBed"
  url=paste("bigDataUrl", paste0("http://amereslab.imba.oeaw.ac.at/seq-hub//dr_slamSeq/dr10/",bbFiles_ends[i]))
  shortLab = paste("shortLabel", paste0("Sample",names_timepoints[i],"ends"),sep=" ")
  longLab =  paste("longLabel", paste0("Sample",names_timepoints[i],"high confidence ends"),sep=" ")
  visibility= paste("visibility full")
  colorUse="colorByStrand 255,0,0 0,0,0"
  main  = paste(trackName,type,url,shortLab,longLab,visibility,colorUse,sep="\n")
  write.table(main,paste0("/Volumes/groups/ameres/wwwdata/seq-hub/dr_slamSeq/dr10/highConfidenceEnds_",names_timepoints[i],".txt"),quote = F,row.names = F,col.names = F)
  
}
includeFiles2 = paste0("include ","highConfidenceEnds_",names_timepoints,".txt")
###### writing similar file for the bed files for refSeq+ ensembl overlapping ends  

bbFiles_ends = list.files (path = "/Volumes/groups/ameres/wwwdata/seq-hub/dr_slamSeq/dr10/",pattern = "refSeq_ensembl_")
names_timepoints = unlist(lapply(strsplit(bbFiles_ends,"_",T),function(x) x[[3]]))

for(i in 1:length(bbFiles_ends)){
  trackName = paste("track",names_timepoints[i])
  type = "type bigBed"
  url=paste("bigDataUrl", paste0("http://amereslab.imba.oeaw.ac.at/seq-hub//dr_slamSeq/dr10/",bbFiles_ends[i]))
  shortLab = paste("shortLabel", paste0("Sample",names_timepoints[i],"primedSites"),sep=" ")
  longLab =  paste("longLabel", paste0("Sample",names_timepoints[i],"primed sites overlapping with UTRs"),sep=" ")
  visibility= paste("visibility full")
  colorUse="colorByStrand 255,0,0 0,0,0"
  main  = paste(trackName,type,url,shortLab,longLab,visibility,colorUse,sep="\n")
  write.table(main,paste0("/Volumes/groups/ameres/wwwdata/seq-hub/dr_slamSeq/dr10/UTRoverlappingPrimingSites_",names_timepoints[i],".txt"),quote = F,row.names = F,col.names = F)
  
}



includeFiles3 = paste0("include ","UTRoverlappingPrimingSites_",names_timepoints,".txt")

allInclude = c(includeFiles1,includeFiles2, includeFiles3)
write.table(allInclude,"/Volumes/groups/ameres/wwwdata/seq-hub/dr_slamSeq/dr10/include_quantSeq_backgroundCWs.txt",quote = F,col.names = F,row.names = F)
