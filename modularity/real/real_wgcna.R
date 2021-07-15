setwd("~/Dropbox/Phd/R/Modularity/Real/WGCNA")
require(WGCNA)
require(flashClust)
require(igraph)
require(org.Sc.sgd.db)
require(limma)
require(qvalue)

source("~/Dropbox/Phd/R/Modularity/ModularityFunctions.R")



exp.path<-("~/Dropbox/Phd/Datasets/Anaerobic_timeseries/GSE22832_expression.txt")
exprs <-as.matrix(read.table(exp.path, header = TRUE, sep = "\t", row.names = 1, as.is=TRUE))
time0<-c("GSM564463","GSM564470")
time0.2<-c("GSM564464","GSM564471")
time1<-c("GSM564465","GSM564472")
time3<-c("GSM564466","GSM564473")
time8<-c("GSM564467","GSM564474")
time24<-c("GSM564468","GSM564475")
time79<-c("GSM564469","GSM564476")

exprs.time<-exprs[,c(time0,time0.2,time1,time3,time8,time24,time79)]

exprs.unique<-unique(rownames(exprs.time))
exprs.time<-exprs.time[exprs.unique,]


exprsFile2 <-("~/Dropbox/Phd/Datasets/Anaerobic_timeseries/GSE22832_ID_ORF.csv")
ID.ORF <-as.matrix(read.table(exprsFile2, header = TRUE, sep = "\t", as.is=TRUE))
ORF<-ID.ORF[match(rownames(exprs.time),ID.ORF),2]
rownames(exprs.time)<-ORF


#--- Step 3: build design matrix
t<-factor(c(0,0,0.2,0.2,1,1,3,3,8,8,24,24,79,79))
design <- model.matrix(~0+t)

#--- Step 4: fit lm model
fit <- lmFit(exprs.time,design)
fit$genes$RefSeq <- rownames(exprs.time)
fit$genes$Symbol <-rownames(exprs.time)

#--- Step 5: set contrasts


contrasts <- makeContrasts(t0.2-t0, t1-t0.2, t3-t1,t8-t3,t24-t8,t79-t24, levels=design)

#--- Step 6: eBayes correction and decision
ct.fit <- eBayes(contrasts.fit(fit, contrasts))
p.value<-ct.fit$p.value

##calculating q values from p values
q.value<-qvalue(p.value, lambda = 0.05, pi0.method = "bootstrap")

# replacing p values with q values back
ct.fit$p.value<-q.value$qvalues

#res.fit1<-decideTests(ct.fit,method="separate", adjust.method="BH", p.value=0.01)
res.fit<-decideTests(ct.fit,method="global", adjust.method="BH", p.value=0.01)

#--- Step 12: get the final data object
yap1.limma <- data.frame(ENTREZ = fit$genes$RefSeq,Symbol = fit$genes$Symbol,logFC = ct.fit$coef, q.value = ct.fit$p.value, degenes = unclass(res.fit), stringsAsFactors = FALSE)

#--- Step 3: get all differentially expressed genes from the yap1.limma data object (see Suppl. Box 1)
idx <- rowSums(yap1.limma[,c(15,16,17,18,19,20)]!=0)
dat <- yap1.limma[idx>0,]

#--- Step 4: get the gene expression matrix for the same set of genes
exp <- exprs.time[match(dat$Symbol,rownames(exprs.time)),]

#--- Step 5: compute a co-expression gene network for genes differentially expressed at 20min (i.e. early response)
dat.time0.2 <- dat[dat$degenes.t0.2!=0, "ENTREZ"]
dat.time1 <- dat[dat$degenes.t1!=0, "ENTREZ"]
dat.time3 <- dat[dat$degenes.t3!=0, "ENTREZ"]
dat.time8 <- dat[dat$degenes.t8!=0, "ENTREZ"]
dat.time24 <- dat[dat$degenes.t24!=0, "ENTREZ"]
dat.time79 <- dat[dat$degenes.t79!=0, "ENTREZ"]
dat.time<-union(c(dat.time0.2,dat.time1,dat.time3,dat.time8,dat.time24),dat.time79)
length(dat.time)

exp.data<-exp[dat.time,]



yap1.exp<-t(exp.data)
dim(yap1.exp)

##Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious
##outliers. We use the function flashClust that provides faster hierarchical clustering than the standard function
##hclust

sampleTree = flashClust(dist(yap1.exp), method = "average");
print("clustering data..")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.

jpeg(file = "sampleClustering_threshold.jpg", width = 700, height = 700);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
cut.height=65;
# Plot a line to show the cut
abline(h=cut.height, col = "red")
dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut.height, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = yap1.exp[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##In the plot, shown in Fig. 2, white means a low value, red a high value, and grey a missing entry.
##The last step is to save the relevant expression and trait data for use in the next steps of the tutorial

print("Calculate softthreshold...")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=50, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


# Plot the results:

jpeg(file = "soft_theshold.jpg", width = 800, height = 500);
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.65,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

### Actual calculations 
softPower = softPower;

wgcna.clust<-function(dataExpr,softPower,minModuleSize,MEDissThre)
{
  ### can strat from using yap1.exp 
  adjacency = adjacency(datExpr, power = softPower);
  
  # Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency);
  dissTOM = 1-TOM
  
  test<-colnames(adjacency)
  TOM.temp<-TOM
  colnames(TOM.temp)<-test
  rownames(TOM.temp)<-test
  
  ## Threshold
  TOM.temp[TOM.temp>1e-4] <- 0
  g <- graph.adjacency(TOM.temp, weighted=TRUE, mode="directed")
  edge.list<-get.edgelist(g)
  weight<-E(g)$weight
  int.edgelist<-cbind(edge.list,weight)
  dim(int.edgelist)
  
  
  # Call the hierarchical clustering function
  geneTree = flashClust(as.dist(dissTOM), method = "average");
  # Plot the resulting clustering tree (dendrogram)
  
  jpeg(file = "gentree_dendogram.jpg", width = 700, height = 500);
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04);
  dev.off()
  
  
  # We like large modules, so we set the minimum module size relatively high:
  minModuleSize = minModuleSize;
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  table(dynamicMods)
  
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  
  png(file = paste("DynamicTreeCut.",minModuleSize,".png",sep=""), width = 1000, height = 1000, res=140);
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  dev.off()
  
  print("Calculate eigengenes...")
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = flashClust(as.dist(MEDiss), method = "average");
  
  
  # Plot the result
  png(file = paste("Eigengenes_modules.png.",minModuleSize,".png",sep=""), width = 1000, height = 1000, res=140);
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  MEDissThres = MEDissThre
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  dev.off()
  
  
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
  
  
  png(file = paste("sampleClustering.",minModuleSize,".png",sep=""), width = 1000, height = 800, res=100);
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  
  ### Visualizing the gene network
  # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
  # the color palette; setting the diagonal to NA also improves the clarity of the plot
  plotDiss = dissTOM^10;
  diag(plotDiss) = NA;
  png(file = paste("Network_heatplot.",minModuleSize,".png",sep=""), width = 1000, height = 800, res=100);
  TOMplot(plotDiss, geneTree, mergedColors, main = "Network heatmap plot, Modular genes")
  dev.off()
  
  
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("blue", standardColors(50));  
  moduleLabels = match(moduleColors, colorOrder)-1;
  MEs = mergedMEs;
  # Save module colors and labels for use in subsequent parts
  save(MEs, moduleLabels, moduleColors, geneTree, 
       file = paste("test-networkConstruction-stepByStep.",minModuleSize,".RData",sep=""))
  
  
  ## locus link 
  ORF<-colnames(datExpr)
  
  ## extracting genes
  ORF.colors<-cbind(ORF,moduleColors,deparse.level = 1)
  
  clust.wgcna<-vector("list", length(unique(ORF.colors[,2])))
  #names(clust.wgcna)<-unique(ORF.colors[,2])
  names(clust.wgcna)<-sprintf("M%i",1:length(clust.wgcna))
  
  for(i in 1:length(unique(ORF.colors[,2])))
  {
    clust.wgcna[[i]]<-(as.vector(subset(ORF.colors, as.character(moduleColors)==unique(ORF.colors[,2])[i], select=ORF)))
  }
  
  return(clust.wgcna)
  
}

#######
setwd("~/Dropbox/Phd/R/Modularity/Real/WGCNA")

wgcna.path<-function(dataExpr,softPower,minModuleSize,MEDissThre,subDir,mainDir)
{
  subDir<-subDir
  dir.create(file.path(mainDir, subDir), showWarnings = TRUE)
  setwd(file.path(mainDir, subDir))
  
  wgcna.cluster<-wgcna.clust(dataExpr,softPower,minModuleSize,MEDissThre)
  
  return(wgcna.cluster)
  
}

softPower=6
setwd("~/Dropbox/Phd/R/Modularity/Real/WGCNA")
wgcna.clust.4<-wgcna.path(dataExpr=yap1.exp,softPower=softPower,minModuleSize=20,MEDissThre=0.1,
                          mainDir<-getwd(),subDir="clust.4")

setwd("~/Dropbox/Phd/R/Modularity/Real/WGCNA")
wgcna.clust.8<-wgcna.path(dataExpr=yap1.exp,softPower=softPower,minModuleSize=10,MEDissThre=0.1,
                          mainDir<-getwd(),subDir="clust.8")

setwd("~/Dropbox/Phd/R/Modularity/Real/WGCNA")
wgcna.clust.12<-wgcna.path(dataExpr=yap1.exp,softPower=20,minModuleSize=10,MEDissThre=0.1,
                           mainDir<-getwd(),subDir="clust.12")

setwd("~/Dropbox/Phd/R/Modularity/Real/WGCNA")
wgcna.clust.16<-wgcna.path(dataExpr=yap1.exp,softPower=25,minModuleSize=6,MEDissThre=0.07,
                           mainDir<-getwd(),subDir="clust.16")

setwd("~/Dropbox/Phd/R/Modularity/Real/WGCNA")

save(wgcna.clust.4,wgcna.clust.8,wgcna.clust.12,wgcna.clust.16,file="clust_wgcna.RData")

#simulation time

require(clusterProfiler)

## comapare cluster for biological process
GO.BP.4 <- compareCluster(wgcna.clust.4, ont="BP", organism="yeast", pvalueCutoff=0.05)
GO.BP.8 <- compareCluster(wgcna.clust.8, ont="BP", organism="yeast", pvalueCutoff=0.05)
GO.BP.12 <- compareCluster(wgcna.clust.12, ont="BP", organism="yeast", pvalueCutoff=0.05)
GO.BP.16 <- compareCluster(wgcna.clust.16, ont="BP", organism="yeast", pvalueCutoff=0.05)


png(file = "clusters_BP.4.png", width = 1000, height = 700, r=120);
plot(GO.BP.4,title =" GO enrichment Biological Process",font.size=10)
dev.off()

png(file = "clusters_BP.8.png", width = 1000, height = 700, r=120);
plot(GO.BP.8,title =" GO enrichment Biological Process",font.size=10)
dev.off()

png(file = "clusters_BP.12.png", width = 1000, height = 700, r=120);
plot(GO.BP.12,title =" GO enrichment Biological Process",font.size=10)
dev.off()

png(file = "clusters_BP.16.png", width = 1000, height = 700, r=120);
plot(GO.BP.16,title =" GO enrichment Biological Process",font.size=10)
dev.off()



mod.4.5<-unlist(Modular.Score(GO.BP.4,TopProcesses=5))
mod.4.4<-unlist(Modular.Score(GO.BP.4,TopProcesses=4))
mod.4.3<-unlist(Modular.Score(GO.BP.4,TopProcesses=3))
mod.4.2<-unlist(Modular.Score(GO.BP.4,TopProcesses=2))
mod.4.1<-unlist(Modular.Score(GO.BP.4,TopProcesses=1))

mod.8.5<-unlist(Modular.Score(GO.BP.8,TopProcesses=5))
mod.8.4<-unlist(Modular.Score(GO.BP.8,TopProcesses=4))
mod.8.3<-unlist(Modular.Score(GO.BP.8,TopProcesses=3))
mod.8.2<-unlist(Modular.Score(GO.BP.8,TopProcesses=2))
mod.8.1<-unlist(Modular.Score(GO.BP.8,TopProcesses=1))

mod.12.5<-unlist(Modular.Score(GO.BP.12,TopProcesses=5))
mod.12.4<-unlist(Modular.Score(GO.BP.12,TopProcesses=4))
mod.12.3<-unlist(Modular.Score(GO.BP.12,TopProcesses=3))
mod.12.2<-unlist(Modular.Score(GO.BP.12,TopProcesses=2))
mod.12.1<-unlist(Modular.Score(GO.BP.12,TopProcesses=1))

mod.16.5<-unlist(Modular.Score(GO.BP.16,TopProcesses=5))
mod.16.4<-unlist(Modular.Score(GO.BP.16,TopProcesses=4))
mod.16.3<-unlist(Modular.Score(GO.BP.16,TopProcesses=3))
mod.16.2<-unlist(Modular.Score(GO.BP.16,TopProcesses=2))
mod.16.1<-unlist(Modular.Score(GO.BP.16,TopProcesses=1))



####################################
## Top 1 GO process for each module
png(file = "Top.modularscore.png", width = 1500, height = 1000, res=120);
par(mar=c(7.1, 5.1, 4.1, 2.1))
dat.1<-c(mod.4.1,mod.8.1,mod.12.1,mod.16.1)
w.1 = c(rep(2,length(mod.4.1)-1),4,rep(2,length(mod.8.1)-1),4,rep(2,length(mod.12.1)-1),4,rep(2,length(mod.16.1)-1),4)
space.1=c(rep(0.2,length(mod.4.1)),1,rep(0.2,length(mod.8.1)-1),1,rep(0.2,length(mod.12.1)-1),1,rep(0.2,length(mod.16.1)-2))
col.1=c(rep("blue",length(mod.4.1)-1),"red",rep("blue",length(mod.8.1)-1),"red",rep("blue",length(mod.12.1)-1),"red",
        rep("blue",length(mod.16.1)-1),"red")

barplot(dat.1,width=w.1, space=space.1, main="Top 1 GO process- WGCNA",
        col=col.1, las=2,ylab="-log10(p-value)",cex.axis=1.5,cex.lab=1.7,cex.main=2)
legend("topleft", legend=c("Module score","Model score"), fill=c("blue","red"), bty = "n",cex=1.5)

line=-37
cex=1.60
add=14
mtext(expression(bold("4 clusters")),at=7,line=line,cex=cex)
mtext(expression(bold("8 clusters")),at=25,line=line,cex=cex)
mtext(expression(bold("12 clusters")),at=40+add,line=line,cex=cex)
mtext(expression(bold("16 clusters")),at=75+add,line=line,cex=cex)
dev.off()


###################################
####################################
## Top 5 to top 2 GO process for each module

png(file = "Modularity_score.png", width = 2500, height = 1000, res=100);
par(mfrow=c(1,4),mar=c(12.1, 4.1, 4.1, 2.1))

module.4=c(mod.4.5,mod.4.4,mod.4.3, mod.4.2)
w.4 = c(rep(2,length(mod.4.5)-1),4,rep(2,length(mod.4.4)-1),4,rep(2,length(mod.4.3)-1),4,rep(2,length(mod.4.2)-1),4)
space.4=c(rep(0.5,length(mod.4.5)),1,rep(0.5,length(mod.4.4)-1),1,rep(0.5,length(mod.4.3)-1),1)
barplot(module.4,las=2,space=space.4, width=w.4, cex.axis=2.5, cex.names=1.7, cex.lab=1.7,cex.main=3,
        col=c(rep("blue",length(mod.4.5)-1),"red"), main=" 4 clusters", ylim=c(0,0.7))
legend("topleft", legend=c("Module score","Model score"), fill=c("blue","red"),bty = "n",cex=3)


module.8=c(mod.8.5,mod.8.4,mod.8.3, mod.8.2)
w.8 = c(rep(2,length(mod.8.5)-1),4,rep(2,length(mod.8.4)-1),4,rep(2,length(mod.8.3)-1),4,rep(2,length(mod.8.2)-1),4)
space.8=c(rep(0.5,length(mod.8.5)),1,rep(0.5,length(mod.8.4)-1),1,rep(0.5,length(mod.8.3)-1),1)
col.8<-c(rep("blue",length(mod.8.5)-1),"red",c(rep("blue",length(mod.8.4)-1),"red"),c(rep("blue",length(mod.8.3)-1),"red"),c(rep("blue",length(mod.8.2)-1),"red"))
barplot(module.8,las=2,space=space.8, width=w.8, cex.axis=2.5, cex.names=1.7, cex.lab=1.7,cex.main=3,
        col=col.8, main=" 8 clusters", ylim=c(0,0.7))

module.12=c(mod.12.5,mod.12.4,mod.12.3, mod.12.2)
w.12 = c(rep(2,length(mod.12.5)-1),4,rep(2,length(mod.12.4)-1),4,rep(2,length(mod.12.3)-1),4,rep(2,length(mod.12.2)-1),4)
space.12=c(rep(0.5,length(mod.12.5)),1,rep(0.5,length(mod.12.4)-1),1,rep(0.5,length(mod.12.3)-1),1)
col.12<-c(rep("blue",length(mod.12.5)-1),"red",c(rep("blue",length(mod.12.4)-1),"red"),c(rep("blue",length(mod.12.3)-1),"red"),c(rep("blue",length(mod.12.2)-1),"red"))
barplot(module.12,las=2,space=space.12,width=w.12,cex.axis=2.5, cex.names=1.7, cex.lab=1.7,cex.main=3,
        col=col.12,main=" 12 clusters",ylim=c(0,0.7))

module.16=c(mod.16.5,mod.16.4,mod.16.3, mod.16.2)
w.16 = c(rep(2,length(mod.16.5)-1),4,rep(2,length(mod.16.4)-1),4,rep(2,length(mod.16.3)-1),4,rep(2,length(mod.16.2)-1),4)
space.16=c(rep(0.5,length(mod.16.5)),1,rep(0.5,length(mod.16.4)-1),1,rep(0.5,length(mod.16.3)-1),1)
col.16<-c(rep("blue",length(mod.16.5)-1),"red",c(rep("blue",length(mod.16.4)-1),"red"),c(rep("blue",length(mod.16.3)-1),"red"),c(rep("blue",length(mod.16.2)-1),"red"))

barplot(module.16,las=2,space=space.16,width=w.16,cex.axis=2.5, cex.names=1.7, cex.lab=1.7,cex.main=3,
        col=col.16,main="16 clusters",ylim=c(0,0.7))



line=-71.5
cex=1.30
add.1<--150
mtext(expression(bold("Top 5 GO")),at=-280+add.1,line=line,cex=cex)
mtext(expression(bold("Top 4 GO")),at=-258+add.1+10,line=line,cex=cex)
mtext(expression(bold("Top 3 GO")),at=-238+add.1+20,line=line,cex=cex)
mtext(expression(bold("Top 2 GO")),at=-215+add.1+30,line=line,cex=cex)

add<-150
mtext(expression(bold("Top 5 GO")),at=-280+add+add.1,line=line,cex=cex)
mtext(expression(bold("Top 4 GO")),at=-258+add+add.1+10,line=line,cex=cex)
mtext(expression(bold("Top 3 GO")),at=-238+add+add.1+20,line=line,cex=cex)
mtext(expression(bold("Top 2 GO")),at=-215+add+add.1+30,line=line,cex=cex)

add.2<-145
mtext(expression(bold("Top 5 GO")),at=-280+add+add.2+add.1,line=line,cex=cex)
mtext(expression(bold("Top 4 GO")),at=-258+add+add.2+add.1+10,line=line,cex=cex)
mtext(expression(bold("Top 3 GO")),at=-238+add+add.2+add.1+20,line=line,cex=cex)
mtext(expression(bold("Top 2 GO")),at=-215+add+add.2+add.1+30,line=line,cex=cex)

add.3<-145
mtext(expression(bold("Top 5 GO")),at=-280+add+add.2+add.1+add.3,line=line,cex=cex)
mtext(expression(bold("Top 4 GO")),at=-258+add+add.2+add.1+add.3+10,line=line,cex=cex)
mtext(expression(bold("Top 3 GO")),at=-238+add+add.2+add.1+add.3+20,line=line,cex=cex)
mtext(expression(bold("Top 2 GO")),at=-215+add+add.2+add.1+add.3+30,line=line,cex=cex)

dev.off()

###################################################################
# function to calculating number of enriched GO terms

cutoff<-0.05

GO.logp.4<-GO.enrich(GO.BP=GO.BP.4,cutoff=cutoff,points=15)
GO.logp.8<-GO.enrich(GO.BP=GO.BP.8,cutoff=cutoff,points=15)
GO.logp.12<-GO.enrich(GO.BP=GO.BP.12,cutoff=cutoff,points=15)
GO.logp.16<-GO.enrich(GO.BP=GO.BP.16,cutoff=cutoff,points=15)



png(file = "GO enriched.png", width = 1000, height = 1000, res=200);
plot(GO.logp.4,col="blue",pch=6, type="o",ylab="No of enriched GO terms",
     xlab=expression(paste(-log[10],"(p-value)")), ylim=c(0,max(GO.logp.4,GO.logp.8,GO.logp.12,GO.logp.16)))

points(GO.logp.8,col="red",pch=16, type="o",ylab="No of enriched GO terms",
       xlab=expression(paste(-log[10],"(p-value)")))

points(GO.logp.12,col="blue",pch=12,type="o", ylab="No of enriched GO terms",
       xlab=expression(paste(-log[10],"(p-value)")))

points(GO.logp.16,col="black",pch=23, type="o",ylab="No of enriched GO terms",
       xlab=expression(paste(-log[10],"(p-value)")))



legend("topright",c("Clust.4","Clust.8","Clust.12","Clust.16"),pch=c(6,16,12,23),col=c("blue","red","blue","black"))

dev.off()



############################
## percentage of clusters that had atleast one enriched GO term

perc.mod.4<-perc.GO(GO.BP=GO.BP.4,cutoff=cutoff,points=15,clust=wgcna.clust.4)
perc.mod.8<-perc.GO(GO.BP=GO.BP.8,cutoff=cutoff,points=15,clust=wgcna.clust.8)
perc.mod.16<-perc.GO(GO.BP=GO.BP.16,cutoff=cutoff,points=15,clust=wgcna.clust.16)
perc.mod.12<-perc.GO(GO.BP=GO.BP.12,cutoff=cutoff,points=15,clust=wgcna.clust.12)

png(file = "Perc_Modules.png", width = 1000, height = 1000, res=200);

plot(perc.mod.4,col="blue",pch=6, type="o", ylab="% annotated modules",
     xlab=expression(paste(-log[10],"(p-value)")),ylim=c(0,100))

points(perc.mod.8,col="red",pch=16, type="o", ylab="% annotated modules",
       xlab=expression(paste(-log[10],"(p-value)")))

points(perc.mod.12,col="blue",pch=12,type="o", ylab="% annotated modules",
       xlab=expression(paste(-log[10],"(p-value)")))

points(perc.mod.16,col="black",pch=23, type="o",ylab="% annotated modules",
       xlab=expression(paste(-log[10],"(p-value)")))

legend("topright",c("Clust.4","Clust.8","Clust.12","Clust.16"),pch=c(6,16,22,23),col=c("blue","red","blue","black"))

dev.off()


