library(DESeq2)
library(corrplot)

##### "differential gene expression" analysis using only one replicate

options(scipen=999)
sampleInfo = read.delim("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/differentialExpressionAnalysis/info/sampleInfo.txt",sep="\t",header = F,stringsAsFactors = F)
sampleInfo$countData = paste0(sampleInfo$V1,"_adapterTrimmed_slamdunk_mapped_filtered_tcount.tsv")

sampleInfo$path = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/differentialExpressionAnalysis/data/countData/",sampleInfo$countData)
sampleInfo_data = lapply(sampleInfo$path,function(x) read.table(x,stringsAsFactors = F,sep="\t",header = T))
names(sampleInfo_data) = sampleInfo$V2

readsCPM = do.call(cbind,lapply(sampleInfo_data,function(x) x$ReadsCPM))
colnames(readsCPM) = sampleInfo$V2

rawReads = do.call(cbind,lapply(sampleInfo_data,function(x) x$ReadCount))

############### vst and PCA plots for all the samples  -- this has to be done using the raw data
colnames(rawReads) = paste0(sampleInfo$V2,c(1:nrow((sampleInfo))))
colData_samples = as.data.frame(colnames(rawReads))
colData_samples$condition = sampleInfo$V3
colData_samples$type = sampleInfo$V2
row.names(colData_samples) = NULL
colnames(colData_samples) = c("Samples","condition","type")

readCounts_deSeq = DESeqDataSetFromMatrix(countData = rawReads,colData = colData_samples,~1)
rld = vst(object = readCounts_deSeq,blind = FALSE)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/differentialExpressionAnalysis/plot/PCA_basedOntimepoint.pdf",height=8,width=8)
DESeq2::plotPCA(rld,intgroup=c("condition")) 
dev.off()

############################## for some samples I have multiple replicates... i want to average these #######################
categories_samples = colnames(readsCPM)
categories_samples = unique(categories_samples)
##### slpitting into different categories based on the number of samples 

data_differentCategories = vector("list",length(categories_samples))
names(data_differentCategories) = categories_samples

#### splitting the different time points into categories..

for(i in 1:length(categories_samples)){
  
  data_differentCategories[[i]] = as.data.frame(readsCPM[,which(colnames(readsCPM) == categories_samples[i])])
  
}

numberReplicates  = unlist(lapply(data_differentCategories,function(x) ncol(x)))
data_differentCategories_moreReplicates = data_differentCategories[which(numberReplicates >1)]

### briefly checking how these "replicates" correlate amongst each other... 

data_differentCategories_moreReplicates_cpm0 = lapply(data_differentCategories_moreReplicates,function(x) x[which(apply(x,1,mean)>0),])

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/differentialExpressionAnalysis/plot/scatter_multipleReplicates.pdf")
lapply(data_differentCategories_moreReplicates_cpm0,function(x) pairs(log10(x)))
dev.off()

###### 

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/differentialExpressionAnalysis/plot/corr_multipleReplicates.pdf")
for(i in 1:length(data_differentCategories_moreReplicates_cpm0)){
  cor_tmp = cor(data_differentCategories_moreReplicates_cpm0[[i]],method = "spearman")
  corrplot(cor_tmp,method = "number")
}
dev.off()



### getting the mean of the differenti 'replicates'

data_differentCategories_moreReplicates_mean = lapply(data_differentCategories_moreReplicates,function(x) apply(x,1,mean))
data_differentCategories_moreReplicates_mean = lapply(data_differentCategories_moreReplicates_mean,function(x) as.data.frame(x))
data_differentCategories_moreReplicates_mean_df = do.call(cbind.data.frame,data_differentCategories_moreReplicates_mean)
colnames(data_differentCategories_moreReplicates_mean_df) = names(data_differentCategories_moreReplicates_mean)

##### now getting those that have only one replicate 
data_differentCategories_oneReplicate  = data_differentCategories[which(numberReplicates == 1)]
data_differentCategories_oneReplicate_df = do.call(cbind.data.frame,data_differentCategories_oneReplicate)
colnames(data_differentCategories_oneReplicate_df) = names(data_differentCategories_oneReplicate)

##### all the data 
allData = cbind.data.frame(data_differentCategories_oneReplicate_df,data_differentCategories_moreReplicates_mean_df)
allData_accordingToTime = allData[,c("testis","oocyte","2cell","64cell","128cell","256cell","512cell","1000cell","high-oblong","oblong","sphere","lateSphere","cap","dome","bud","1dpf","2dpf","4dpf")]
allData_accordingToTime_timepoints = allData_accordingToTime[,c(3:18)]
checkColumns = cbind.data.frame(c(1:15),c(2:16))
colnames(checkColumns) = c("col1","col2")

consecutiveColumnCheck = matrix(data = "NA",nrow = nrow(allData),ncol = nrow(checkColumns))
for(i in 1:nrow(checkColumns)){
  
  smallpart = allData_accordingToTime[,c(checkColumns[i,1]:checkColumns[i,2])]
  smallPart_greater = apply(smallpart,1,function(x) length(which(x>2)))
  smallPart_greater[which(smallPart_greater == 2)] <- "YES"
  smallPart_greater[which(smallPart_greater < 2)] <- "NO"
  consecutiveColumnCheck[,i] = smallPart_greater
}

numberYES = apply(consecutiveColumnCheck,1,function(x) length(which(x=="YES")))
numberYES[which(numberYES>0)]<-"KEEP"
numberYES[which(numberYES==0)]<-"NOKEEP"

allData_accordingToTime$threshold = numberYES
metaDataTable = sampleInfo_data[[1]][,c(1:6)]
metaDataTable$threshold = numberYES
metaDataTable_keep = metaDataTable[which(metaDataTable$threshold == "KEEP"),]
allData_accordingToTime_keep = allData_accordingToTime[which(allData_accordingToTime$threshold == "KEEP"),]
##################### now I have a threshold for the further calculations... 


allData_accordingToTime_pre = allData_accordingToTime_keep[,c("2cell","64cell","128cell","256cell","512cell","1000cell")]
allData_accordingToTime_post = allData_accordingToTime_keep[,c("high-oblong","oblong","sphere","lateSphere","cap","dome","bud","1dpf","2dpf","4dpf")]
allData_accordingToTime_late = allData_accordingToTime_keep[,c("1dpf","2dpf","4dpf")]


####### putting the three dfs in a list 


stages_zf = vector("list",3)
stages_zf[[1]] = allData_accordingToTime_pre
stages_zf[[2]] = allData_accordingToTime_post
stages_zf[[3]] = allData_accordingToTime_late
names(stages_zf) = c("pre","post","late")

stages_zf_mean  = lapply(stages_zf,function(x) apply(x,1,mean))
stages_zf_mean = do.call(cbind.data.frame,stages_zf_mean)
pairs(log10(stages_zf_mean))
cor(stages_zf_mean,method = "spearman")
names_all = metaDataTable_keep$Name

stages_zf_mean$foldChange = stages_zf_mean$post/stages_zf_mean$pre
stages_zf_mean  = cbind.data.frame(stages_zf_mean,metaDataTable_keep)
#stages_zf_mean$Name =names_all
stages_zf_mean$classification = NA
stages_zf_mean$classification[which(stages_zf_mean$foldChange>1.2)] <- "ZYGOTIC"
stages_zf_mean$classification[which(stages_zf_mean$foldChange<0.8)] <- "MATERNAL"
stages_zf_mean$classification[-which(stages_zf_mean$classification == "ZYGOTIC" | stages_zf_mean$classification == "MATERNAL")] <- "STABLE"

samples_split = split(stages_zf_mean,stages_zf_mean$Name,drop = T)

samples_split_multiple = samples_split[which(lapply(samples_split,nrow)>1)]
names_increaseFoldChange = unlist(lapply(samples_split_multiple,function(x) length(which(x$foldChange>1.2))))
genes_increasedFC = which(names_increaseFoldChange>0)
genes_increasedFC = names(genes_increasedFC)
names_decreaseFoldChange = unlist(lapply(samples_split,function(x) length(which(x$foldChange<0.8))))
genes_decreasedFC = which(names_decreaseFoldChange>0)
names_decreaseFoldChange = names(genes_decreasedFC)


changeInIsoform = intersect(genes_increasedFC,names_decreaseFoldChange)
changeInIsoform = stages_zf_mean[stages_zf_mean$Name %in% changeInIsoform,]

#######  now to classify into proximal and distal transcripts 


###### plus and minus strand separation 

changeInIsoform_plus = changeInIsoform[which(changeInIsoform$Strand == "+"),]
changeInIsoform_minus = changeInIsoform[which(changeInIsoform$Strand == "-"),]

########### plus strand

changeInIsoform_plus_split = split(changeInIsoform_plus,changeInIsoform_plus$Name,T)
changeInIsoform_plus_split = lapply(changeInIsoform_plus_split,function(x) x[order(x$Start,decreasing = F),])

changeInIsoform_plus_split_proximal = lapply(changeInIsoform_plus_split,function(x) x[1,])
changeInIsoform_plus_split_distal = lapply(changeInIsoform_plus_split,function(x) x[2:nrow(x),])

######### minus strand
changeInIsoform_minus_split = split(changeInIsoform_minus,changeInIsoform_minus$Name,T)
changeInIsoform_minus_split = lapply(changeInIsoform_minus_split,function(x) x[order(x$Start,decreasing = T),])

changeInIsoform_minus_split_proximal = lapply(changeInIsoform_minus_split,function(x) x[1,])
changeInIsoform_minus_split_distal = lapply(changeInIsoform_minus_split,function(x) x[2:nrow(x),])


################################### total 

proximalIsoforms = rbind(do.call(rbind,changeInIsoform_plus_split_proximal),do.call(rbind,changeInIsoform_minus_split_proximal))
distalIsoforms = rbind(do.call(rbind,changeInIsoform_plus_split_distal),do.call(rbind,changeInIsoform_minus_split_distal))

#### fishers exact test to check for over-representation of maternal isoforms in profimal dataset
proximalIsoforms_maternal=table(proximalIsoforms$classification)[which(names(table(proximalIsoforms$classification))=="MATERNAL")]
proximalIsoforms_zygotic=table(proximalIsoforms$classification)[which(names(table(proximalIsoforms$classification))=="ZYGOTIC")]
distalIsoforms_maternal=table(distalIsoforms$classification)[which(names(table(distalIsoforms$classification))=="MATERNAL")]
distalIsoforms_zygotic=table(distalIsoforms$classification)[which(names(table(distalIsoforms$classification))=="ZYGOTIC")]

contingency_proximal = c(proximalIsoforms_maternal,proximalIsoforms_zygotic)
contingency_distal = c(distalIsoforms_maternal,distalIsoforms_zygotic)
contingencyTable = rbind.data.frame(contingency_proximal,contingency_distal)
colnames(contingencyTable) = c("Maternal","Zygotic")
rownames(contingencyTable) =  c("Proximal","Distal")
fisher.test(contingencyTable)


######### 


write.table(changeInIsoform,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/differentialExpressionAnalysis/data/changeInIsoform_oneReplicate.txt",sep="\t",quote = F,row.names = F,col.names = T)

#### comparing this with genes that show a difference in 5' end usage from haeble et.al 

changeIsoform = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/differentialExpressionAnalysis/data/changeInIsoform_oneReplicate.txt",sep="\t",header = T)
promoterSwitch = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/externalData/haebleEtal/switchedPromoters_rearranged.txt",sep="\t",header=T,stringsAsFactors = F)

length(intersect(changeIsoform$Name,promoterSwitch$external_gene_id))
