#2717 mouse embryonic stem cells��with 24174 observed unique transcripts
library(edgeR)
library(DESeq)
library(ggplot2)
library(igraph)
library(slingshot)
library(Rtsne)
# library(umap)
library(cluster)
library(seewave)
library(factoextra)
library(mclust)
setwd("C:/Users/Documents/R_projects/mouse_ESC")

data_d0<-read.table("GSM1599494_ES_d0_main.csv",sep=",",header=T,row.names=1)
data_d2<-read.table("GSM1599497_ES_d2_LIFminus.csv",sep=",",header=T,row.names=1)
data_d4<-read.table("GSM1599498_ES_d4_LIFminus.csv",sep=",",header=T,row.names=1)
data_d7<-read.table("GSM1599499_ES_d7_LIFminus.csv",sep=",",header=T,row.names=1)

colnames(data_d0)<-paste("d0_",colnames(data_d0),sep="")
colnames(data_d2)<-paste("d2_",colnames(data_d2),sep="")
colnames(data_d4)<-paste("d4_",colnames(data_d4),sep="")
colnames(data_d7)<-paste("d7_",colnames(data_d7),sep="")

mouse_ESC_data<-cbind(data_d0,data_d2,data_d4,data_d7)
mouse_ESC_data<- mouse_ESC_data[rowMeans(mouse_ESC_data > 1) > 0.3,]

#CPM normalization
group = c(rep("d0",length(data_d0)),rep("d2",length(data_d2)),rep("d4",length(data_d4)),rep("d7",length(data_d7)))
dge <- DGEList(
  counts = mouse_ESC_data, 
 norm.factors = rep(1, length(mouse_ESC_data[1,])), 
  group = group
)
# dge<-calcNormFactors(dge)
group_edgeR <- factor(group)
design <- model.matrix(~ group_edgeR)
dge <- estimateDisp(dge, design = design, trend.method = "none")
data_cpm = cpm(dge)
write.table(data_cpm,"forSRfilter.txt")
log_data_cpm<-log2(data_cpm+1)
write.table(log_data_cpm,"filter_data.txt")

## dimensionality reduction
filter_data<-read.table("filter_data.txt",row.names = 1)
data<-filter_data
tsne_data<-Rtsne(t(as.matrix(data)),pca = T,theta = 0.0)
plot(tsne_data$Y)
tsne_data<-as.data.frame(tsne_data$Y)
rownames(tsne_data)<-colnames(data)
colnames(tsne_data)<-c("PC1","PC2")

#umap
# umap_data<-umap(t(as.matrix(data)))
# umap_data<-as.data.frame(umap_data$layout)
# plot(umap_data)
# rownames(umap_data)<-colnames(data)
# colnames(umap_data)<-c("PC1","PC2")
##clustering
gap_stat<-clusGap(tsne_data,kmeans,K.max = 8,B=200,verbose=interactive())
best_k<-maxSE(f = gap_stat$Tab[, "gap"], SE.f = gap_stat$Tab[, "SE.sim"])
# #fviz_gap_stat(gap_stat)
clustered_data<-kmeans(tsne_data,best_k)
cluster_result<-clustered_data$cluster

## compute potency states
ESC_SR<-read.table("ESC_SR.txt")
rownames(ESC_SR)<-colnames(filter_data)
colnames(ESC_SR)<-c("sample_name","SR")

infer_SR_states<-function(logSR,max_state_num=3)
{
  print("Fit Gaussian Mixture Model to Signaling Entropies")
  model.o=Mclust(logSR,G=seq_len(max_state_num))  
  potS <- model.o$classification
  nPS <- length(levels(as.factor(potS)))
  print(paste("Identified ",nPS," potency states",sep=""))
  #plot(model.o,what = "classification")
  return(model.o)
}
logSR<-log2(ESC_SR$SR/(1-ESC_SR$SR))
model.o<-infer_SR_states(logSR)
table(model.o$classification)

#The probability distribution of heterogeneous states in each cell cluster was calculated
cluster_result<-as.data.frame(cluster_result)
cluster_result$cellid<-rownames(cluster_result)
km_SRstates<-cbind(cluster_result,model.o$classification)
colnames(km_SRstates)<-c("cluster","cellid","states")
print(best_k)
#cluster1
# cons=0.0001
cluster1_freq1=nrow(km_SRstates[which(km_SRstates$cluster==1&km_SRstates$states==1),])/length(km_SRstates[which(km_SRstates$cluster==1),]$cluster)
cluster1_freq2=nrow(km_SRstates[which(km_SRstates$cluster==1&km_SRstates$states==2),])/length(km_SRstates[which(km_SRstates$cluster==1),]$cluster)
cluster1_freq3=nrow(km_SRstates[which(km_SRstates$cluster==1&km_SRstates$states==3),])/length(km_SRstates[which(km_SRstates$cluster==1),]$cluster)
clu1_distrib<-data.frame(states<-c(1:3),freq<-c(cluster1_freq1,cluster1_freq2,cluster1_freq3))
colnames(clu1_distrib)<-c("states","freq")
print(clu1_distrib)
#cluster2
cluster2_freq1=nrow(km_SRstates[which(km_SRstates$cluster==2&km_SRstates$states==1),])/length(km_SRstates[which(km_SRstates$cluster==2),]$cluster)
cluster2_freq2=nrow(km_SRstates[which(km_SRstates$cluster==2&km_SRstates$states==2),])/length(km_SRstates[which(km_SRstates$cluster==2),]$cluster)
cluster2_freq3=nrow(km_SRstates[which(km_SRstates$cluster==2&km_SRstates$states==3),])/length(km_SRstates[which(km_SRstates$cluster==2),]$cluster)
clu2_distrib<-data.frame(states<-c(1:3),freq<-c(cluster2_freq1,cluster2_freq2,cluster2_freq3))
colnames(clu2_distrib)<-c("states","freq")
print(clu2_distrib)
#cluster3
cluster3_freq1=nrow(km_SRstates[which(km_SRstates$cluster==3&km_SRstates$states==1),])/length(km_SRstates[which(km_SRstates$cluster==3),]$cluster)
cluster3_freq2=nrow(km_SRstates[which(km_SRstates$cluster==3&km_SRstates$states==2),])/length(km_SRstates[which(km_SRstates$cluster==3),]$cluster)
cluster3_freq3=nrow(km_SRstates[which(km_SRstates$cluster==3&km_SRstates$states==3),])/length(km_SRstates[which(km_SRstates$cluster==3),]$cluster)
clu3_distrib<-data.frame(states<-c(1:3),freq<-c(cluster3_freq1,cluster3_freq2,cluster3_freq3))
colnames(clu3_distrib)<-c("states","freq")
print(clu3_distrib)

all_distrib<-list(as.matrix(clu1_distrib),as.matrix(clu2_distrib)
                  ,as.matrix(clu3_distrib)
)
#Calculate the JS divergence distance
JS_compute<-function(prob_distrib1,prob_distrib2){
  prob1=prob_distrib1[,2]
  prob2=prob_distrib2[,2]
  M=(prob1+prob2)/2
  return(0.5*(kl.dist(prob1,M,2)$D1)+0.5*(kl.dist(prob2,M,2)$D1))
}
JS_dist=matrix(0,nrow=length(all_distrib),ncol=length(all_distrib))
for (i in seq_len(length(levels(as.factor(km_SRstates$cluster))))) {
  for(j in seq_len(length(levels(as.factor(km_SRstates$cluster))))){
    JS_dist[i,j]=JS_compute(all_distrib[[i]],all_distrib[[j]])
  }
}

merge_data<-cbind(ESC_SR,cluster_result)
merge_data<-merge_data[,c(1,2,3)]
colnames(merge_data)<-c("cellid","SR","cluster")
#write.csv(merge_data,"SR_cluster.csv")
avgSR_cluster1<-mean(merge_data[which(merge_data$cluster=="1"),]$SR)
avgSR_cluster2<-mean(merge_data[which(merge_data$cluster=="2"),]$SR)
avgSR_cluster3<-mean(merge_data[which(merge_data$cluster=="3"),]$SR)

avgSR_cluster<-data.frame(avgSR_cluster1,avgSR_cluster2,avgSR_cluster3)
start_cluster<-colnames(avgSR_cluster)[which(avgSR_cluster==max(avgSR_cluster),arr.ind = T)[,2]]
#write.table(avgSR_cluster,"averageSR_cluster.txt",col.names = F)
##########################################################
#######compute MST
##########################################################
# reduceDimdata_matrix=read.table("tsne_reducedata.txt",row.names = 1)
reduceDimdata_matrix=as.matrix(tsne_data)
JS_dist<-as.matrix(JS_dist)
gp <- graph.adjacency(JS_dist, mode = "undirected", weighted = TRUE)
dp_mst <- minimum.spanning.tree(gp) 
print(dp_mst)
###########################################################
###########Construct the simultaneous curves and compute pseudotime values
###########################################################
clusLabel=cluster_result$cluster_result
clusLabel=t(clusLabel)
clusLabel=as.numeric(clusLabel)
names(clusLabel)<-cluster_result$cellid
lineage<-getLineages(reduceDimdata_matrix,clusLabel)
lineage_list=list(Lineage1=c("3","1","2"))
#Update the parameters of slingparams
params=slingParams(lineage)
slingParams=params
slingParams$start.clus<-"3"
slingParams$end.clus<-"2"
slingParams$end.given<-FALSE
inputofcurves<-newSlingshotDataSet(reducedDim = reduceDimdata_matrix,clusterLabels = clusLabel,lineages=lineage_list,adjacency=JS_dist,curves=list(),slingParams=slingParams)
outofcurves<-getCurves(inputofcurves,extend="n")
pseudotime<-outofcurves@curves$curve1$lambda
pseudotime.sort<-sort(pseudotime,decreasing = F)
write.table(pseudotime,"ESC_pseudotime.txt")
#
