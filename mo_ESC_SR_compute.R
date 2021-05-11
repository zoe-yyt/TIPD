#mouse ESC 2717 samples
#SR compute
library(clusterProfiler)
library(org.Mm.eg.db)#mouse
library(stringr)
library(Matrix)
library(base)
library(AnnotationDbi)
library(igraph)
library(ggplot2)
setwd("C:/Users/yyt/Documents/R_projects/mouse_ESC")
options(stringsAsFactors = FALSE) 
protein_interactions<-read.table("10090.protein.links.v11.0.txt",header = T)
#idmapping<-read.table("HUMAN_9606_idmapping_selected.tab",header=F,as.is=T,sep="\t")
filter_data<-read.table("filter_data.txt",header=T)
geneid<-rownames(filter_data)

gene_symbol_to_ensp<-select(org.Mm.eg.db,keys=geneid,column="ENSEMBLPROT",keytype="SYMBOL")
removeNAindex<-which(is.na(gene_symbol_to_ensp[,2]))
gene_symbol_to_ensp<-gene_symbol_to_ensp[-removeNAindex,]
protein_interactions$protein1<-str_replace(protein_interactions$protein1,"10090.","")
protein_interactions$protein2<-str_replace(protein_interactions$protein2,"10090.","")
protein_interactions<-protein_interactions[which(protein_interactions[,1] %in% gene_symbol_to_ensp[,2] & protein_interactions[,2] %in% gene_symbol_to_ensp[,2]),]
# remove duplicates
protein_interactions$identical<-paste0(protein_interactions[,1],protein_interactions[,2])
protein_interactions<-protein_interactions[!duplicated(protein_interactions$identical),]
protein_interactions$identical<-paste0(protein_interactions[,2],protein_interactions[,1])
protein_interactions<-protein_interactions[!duplicated(protein_interactions$identical),]
ppi_network<-protein_interactions[,c(1,2)]
write.table(ppi_network,"ppi_network.txt")

forSRfilter<-read.table("forSRfilter.txt",header=T,row.names = 1)#Gene expression data that has been gene-filtered but not logarithmic
ppi_network<-read.table("ppi_network.txt",header=T)
forSRfilter<-forSRfilter[which(rownames(forSRfilter) %in% geneid),]
#ensp tranfer to symbol
ppi_network1<-as.character(ppi_network[,1])
ppi_network2<-as.character(ppi_network[,2])
annot1<-select(org.Mm.eg.db,keys=ppi_network1,columns="SYMBOL",keytype="ENSEMBLPROT")
annot2<-select(org.Mm.eg.db,keys=ppi_network2,columns="SYMBOL",keytype="ENSEMBLPROT")
ppi_network<-cbind(annot1,annot2)
ppi_network<-ppi_network[,c(2,4)]
removeNAindex1<-which(is.na(ppi_network[,1]))
removeNAindex2<-which(is.na(ppi_network[,2]))

edgelist = ppi_network[,c(1,2)]
write.csv(edgelist,file="edgelist.csv",row.names = F)
node_union<-union(ppi_network[,1],ppi_network[,2])

genes<-rownames(forSRfilter)
SR_filter_data<-forSRfilter[which(genes %in% node_union),]
setdiff(node_union,rownames(SR_filter_data))
node_union=union(edgelist[,1],edgelist[,2])
setdiff(node_union,rownames(SR_filter_data))
removelist<-c(setdiff(node_union,rownames(SR_filter_data)))
index1<-which(edgelist[,1] %in% removelist)
index2<-which(edgelist[,2] %in% removelist)
index<-c(index1,index2)
index<-unique(index)
edgelist<-edgelist[-index,]
node_union=union(edgelist[,1],edgelist[,2])
setdiff(node_union,rownames(SR_filter_data))
g1<-graph.data.frame(edgelist,directed=F)
g2<-get.adjacency(g1,sparse=FALSE)

SR_filter_data<forSRfilter[which(rownames(forSRfilter) %in% rownames(SR_filter_data)),]

DoIntegPPI <- function(exp.m, 
                       ppiA.m,
                       log_trans = TRUE)
{
  # set input data matrix class
  data.class <- class(exp.m)
  
  # select log-normalization methods based on data class
  if (data.class == "SingleCellExperiment") {
    sizeFactors(exp.m) <- scater::librarySizeFactors(exp.m)
    exp.m <- scater::normalize(exp.m, log_exprs_offset = 1.1)
    data.m <- Matrix::as.matrix(SummarizedExperiment::assay(exp.m, i = "logcounts"))
  }else if (data.class == "CellDataSet") {
    exp.m <- BiocGenerics::estimateSizeFactors(exp.m)
    data.m <- Matrix::as.matrix(t(t(Biobase::exprs(exp.m)) / 
                                    Biobase::pData(exp.m)[, 'Size_Factor']))
    data.m <- log2(data.m + 1.1)
  }else{
    data.m <- as.matrix(exp.m)
    if (log_trans) {
      TRC.v <- colSums(exp.m)
      maxC <- max(TRC.v)
      for (i in seq_len(dim(data.m)[2])) {
        temp <- maxC / TRC.v[i]
        data.m[, i] <- log2(exp.m[, i] * temp + 1.1)
      }
    }
  }
  
  if (min(data.m) == 0) {
    stop("Input matrix must have non-zero minimal value, please set 
         log_trans = TRUE!")
  }
  
  commonEID.v <- intersect(rownames(ppiA.m),rownames(data.m))
  
  if (DelayedArray::isEmpty(commonEID.v) == TRUE) {
    stop("scRNA-seq data should have the same gene identifier with the network!")
  }
  
  match(commonEID.v,rownames(data.m)) -> map1.idx
  expPIN.m <- data.m[map1.idx,]
  
  match(commonEID.v,rownames(ppiA.m)) -> map2.idx
  adj.m <- ppiA.m[map2.idx,map2.idx]
  
  gr.o <- igraph::graph.adjacency(adj.m,mode="undirected")
  comp.l <- igraph::clusters(gr.o)
  cd.v <- summary(factor(comp.l$member))
  mcID <- as.numeric(names(cd.v)[which.max(cd.v)])
  maxc.idx <- which(comp.l$member==mcID)
  adjMC.m <- adj.m[maxc.idx,maxc.idx]
  expMC.m <- expPIN.m[maxc.idx,]
  
  if (data.class == "SingleCellExperiment") {
    return(list(data.sce = exp.m, expMC = expMC.m, adjMC = adjMC.m, data = data.m))
  }else if (data.class == "CellDataSet") {
    return(list(data.cds = exp.m, expMC = expMC.m, adjMC = adjMC.m, data = data.m))
  }else{
    return(list(expMC = expMC.m, adjMC = adjMC.m, data = data.m))
  }
  }
Integration.l <- DoIntegPPI(SR_filter_data,g2)

CompSRana <- function(Integration.l,
                      local = FALSE,
                      mc.cores=1)
{
  ### compute maxSR for SR normalization
  Integration.l <- CompMaxSR(Integration.l)
  maxSR <- Integration.l$maxSR
  
  idx.l <- as.list(seq_len(ncol(Integration.l$expMC)))
  out.l <- mclapply(idx.l, CompSRanaPRL, 
                    exp.m=Integration.l$expMC, 
                    adj.m=Integration.l$adjMC,
                    local=local,
                    maxSR=maxSR,
                    mc.cores=mc.cores)
  SR.v <- sapply(out.l, function(v) return(v[[1]]))
  invP.v <- sapply(out.l, function(v) return(v[[2]]))
  S.v <- sapply(out.l, function(v) return(v[[3]]))
  NS.v <- sapply(out.l, function(v) return(v[[4]]))
  
  Integration.l$SR <- SR.v 
  Integration.l$inv <- invP.v
  Integration.l$s <- S.v
  Integration.l$ns <- NS.v
  
  if (!is.null(Integration.l$data.sce)) {
    colData(Integration.l$data.sce)$SR <- SR.v
  }else if (!is.null(Integration.l$data.cds)) {
    pData(Integration.l$data.cds)$SR <- SR.v
  }
  return(Integration.l)
}

CompSRana <- function(Integration.l,
                      local = FALSE,
                      mc.cores=1)
{
  ### compute maxSR for SR normalization
  Integration.l <- CompMaxSR(Integration.l)
  maxSR <- Integration.l$maxSR
  
  idx.l <- as.list(seq_len(ncol(Integration.l$expMC)))
  out.l <- mclapply(idx.l, CompSRanaPRL, 
                    exp.m=Integration.l$expMC, 
                    adj.m=Integration.l$adjMC,
                    local=local,
                    maxSR=maxSR,
                    mc.cores=mc.cores)
  SR.v <- sapply(out.l, function(v) return(v[[1]]))
  invP.v <- sapply(out.l, function(v) return(v[[2]]))
  S.v <- sapply(out.l, function(v) return(v[[3]]))
  NS.v <- sapply(out.l, function(v) return(v[[4]]))
  
  Integration.l$SR <- SR.v
  Integration.l$inv <- invP.v
  Integration.l$s <- S.v
  Integration.l$ns <- NS.v
  
  if (!is.null(Integration.l$data.sce)) {
    colData(Integration.l$data.sce)$SR <- SR.v
  }else if (!is.null(Integration.l$data.cds)) {
    pData(Integration.l$data.cds)$SR <- SR.v
  }
  return(Integration.l)
}

CompMaxSR <- function(Integration.l){
  
  adj.m <- Integration.l$adjMC
  
  # find right eigenvector of adjacency matrix
  fa <- function(x,extra=NULL) {
    as.vector(adj.m %*% x)
  }
  ap.o <- igraph::arpack(fa,options=list(n=nrow(adj.m),nev=1,which="LM"), sym=TRUE)
  v <- ap.o$vectors
  lambda <- ap.o$values
  
  # maximum entropy
  MaxSR <- log(lambda)
  
  Integration.l$maxSR <- MaxSR
  
  return(Integration.l)
}

CompSRanaPRL <- function(idx,
                         exp.m,
                         adj.m,
                         local=TRUE,
                         maxSR=NULL)
{
  
  # compute outgoing flux around each node
  exp.v <- exp.m[,idx];
  sumexp.v <- as.vector(adj.m %*% matrix(exp.v,ncol=1));
  invP.v <- exp.v*sumexp.v;
  nf <- sum(invP.v);
  invP.v <- invP.v/nf;
  p.m <- t(t(adj.m)*exp.v)/sumexp.v;
  S.v <- apply(p.m,1,CompS);
  SR <- sum(invP.v*S.v);
  # if provided then normalise relative to maxSR
  if(is.null(maxSR)==FALSE){
    SR <- SR/maxSR;
  }
  if(local){
    NS.v <- apply(p.m,1,CompNS);
  }
  else {
    NS.v <- NULL;
  }
  return(list(sr=SR,inv=invP.v,s=S.v,ns=NS.v));
}

CompNS <- function(p.v){
  
  tmp.idx <- which(p.v>0);
  if(length(tmp.idx)>1){
    NLS <- -sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )/log(length(tmp.idx));
  }
  else {
    # one degree nodes have zero entropy, avoid singularity.
    NLS <- 0;
  }
  return(NLS);
}

CompS <- function(p.v){
  
  tmp.idx <- which(p.v>0);
  LS <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )
  return(LS);
}
Integration.l<-CompSRana(Integration.l)
print(Integration.l$SR)
summary(Integration.l$SR)

SR<-as.data.frame(Integration.l$SR)
rownames(SR)<-colnames(filter_data)
write.table(SR,"ESC_SR.txt",quote=F,col.names = F)
#