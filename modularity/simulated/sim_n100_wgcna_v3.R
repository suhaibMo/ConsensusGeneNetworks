setwd("~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/")
require(WGCNA)
require(flashClust)
require(igraph)
require(org.Sc.sgd.db)

source("~/Dropbox/Phd/R/Modularity/ModularityFunctions.R")
## 500 nodes simulated expression data
exp.path<-("~/Dropbox/Phd/R/Consensus/Simulated/syntren/N100/samples/100/Noise1/nn100_nbgr0_hop0.3_bionoise0.1_expnoise0.1_corrnoise0.1_neighAdd_dataset.txt")
exp.data<-(as.matrix(read.table(exp.path, header = TRUE, sep = "\t", as.is=TRUE)))
exp.data<-t(exp.data)
dim(exp.data)

####################################
## changing gene names to ORF's
rownames(exp.data)<-gsub("_.*","",rownames(exp.data))

if ("MFALPHA2" %in% rownames(exp.data))
{
  rownames(exp.data)<-gsub("MFALPHA2","MF(ALPHA)2",rownames(exp.data))
}

orfs<-matrix()
for (i in 1:nrow(exp.data))
{
  orfs[i]<-get(rownames(exp.data)[i],org.Sc.sgdCOMMON2ORF)
}  

rownames(exp.data)<-orfs
head(exp.data)

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
setwd("~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP/")

wgcna.path<-function(dataExpr,softPower,minModuleSize,MEDissThre,subDir,mainDir)
{
  subDir<-subDir
  dir.create(file.path(mainDir, subDir), showWarnings = TRUE)
  setwd(file.path(mainDir, subDir))

  wgcna.cluster<-wgcna.clust(dataExpr,softPower,minModuleSize,MEDissThre)

  return(wgcna.cluster)

 }

softPower=6
                           
setwd("~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP/")
W.100.clust.4<-wgcna.path(dataExpr=yap1.exp,softPower=softPower,minModuleSize=4,MEDissThre=0.9,
                           mainDir<-getwd(),subDir="clust.4")

setwd("~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP")
W.100.clust.8<-wgcna.path(dataExpr=yap1.exp,softPower=softPower,minModuleSize=6,MEDissThre=0.2,
                           mainDir<-getwd(),subDir="clust.8")

setwd("~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP")
W.100.clust.12<-wgcna.path(dataExpr=yap1.exp,softPower=10.5,minModuleSize=0,MEDissThre=0,
                           mainDir<-getwd(),subDir="clust.12")

setwd("~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP")
W.100.clust.16<-wgcna.path(dataExpr=yap1.exp,softPower=25,minModuleSize=0,MEDissThre=0,
                          mainDir<-getwd(),subDir="clust.16")

setwd("~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA")

save(W.100.clust.4,W.100.clust.8,W.100.clust.12,W.100.clust.16,
     file="W.clusters.100.rda")


wgcna.clust<-list(W.100.clust.4,W.100.clust.8,W.100.clust.12,W.100.clust.16)

#simulation time

	
############################################################################################
require(fpc)
require(cluster)
W.100.mean.sil<-c()
W.100.avg.sil1<-c()
W.100.avg.sil2<-c()
W.100.dunn<-c()
W.100.dunn2<-c()
W.100.sep.idx<-c()
k=c(4,8,12,16)
for (i in 1:length(k)) {
    print (i)
    dface<-dist(exp.data,method="euclidean")
    complete16<-clust.vector(wgcna.clust[[i]],exp.data)
    silwidths<-cluster.stats(d=dface,complete16,silhouette = TRUE)
    W.100.mean.sil[i]<-mean(as.numeric(silwidths$clus.avg.silwidths))
    W.100.avg.sil1[i]<-as.numeric(silwidths$avg.silwidth)
    W.100.dunn[i]<-as.numeric(silwidths$dunn)
    W.100.dunn2[i]<-as.numeric(silwidths$dunn2)
    W.100.sep.idx[i]<-as.numeric(silwidths$sindex)
}

par(mfrow=c(3,3))
plot(W.100.mean.sil,ylim=c(min(W.100.mean.sil),max(W.100.mean.sil)),type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules',
     ylab="mean silhouette width")
axis(1, at=1:length(k), labels=c(k))

plot(W.100.avg.sil1,ylim=c(min(W.100.avg.sil1),max(W.100.avg.sil1)),type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules',
     ylab="Avg silhouette width")
axis(1, at=1:length(k), labels=c(k))
plot(W.100.dunn,ylim=c(min(W.100.dunn),max(W.100.dunn)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
     ylab="Dunn index")
axis(1, at=1:length(k), labels=c(k))
plot(W.100.dunn2,ylim=c(min(W.100.dunn2),max(W.100.dunn2)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
     ylab="Dunn2 index")
axis(1, at=1:length(k), labels=c(k))

plot(W.100.sep.idx,ylim=c(min(W.100.sep.idx),max(W.100.sep.idx)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
     ylab="sep index")
axis(1, at=1:length(k), labels=c(k))

save(W.100.mean.sil,W.100.avg.sil1,W.100.dunn,W.100.dunn2,W.100.sep.idx,file="W.cluster.valid.100.rda")


require(clusterProfiler)

## comapare cluster for biological process
W.100.GO.BP.4 <- compareCluster(W.100.clust.4, ont="BP", organism="yeast", pvalueCutoff=0.05)
W.100.GO.BP.8 <- compareCluster(W.100.clust.8, ont="BP", organism="yeast", pvalueCutoff=0.05)
W.100.GO.BP.12 <- compareCluster(W.100.clust.12, ont="BP", organism="yeast", pvalueCutoff=0.05)
W.100.GO.BP.16 <- compareCluster(W.100.clust.16, ont="BP", organism="yeast", pvalueCutoff=0.05)


setwd("~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP/")
png(file = "clusters_BP.4.png", width = 1000, height = 700, r=120);
plot(W.100.GO.BP.4,title =" GO enrichment Biological Process - WGCNA",font.size=10)
dev.off()

png(file = "clusters_BP.8.png", width = 1000, height = 700, r=120);
plot(W.100.GO.BP.8,title =" GO enrichment Biological Process - WGCNA",font.size=10)
dev.off()

png(file = "clusters_BP.12.png", width = 1000, height = 700, r=120);
plot(W.100.GO.BP.12,title =" GO enrichment Biological Process - WGCNA",font.size=10)
dev.off()

png(file = "clusters_BP.16.png", width = 1000, height = 700, r=120);
plot(W.100.GO.BP.16,title =" GO enrichment Biological Process - WGCNA",font.size=10)
dev.off()



####
# biological process
W.100.mod.BP.4.1<-unlist(Modular.Score(W.100.GO.BP.4,TopProcesses=1))
W.100.mod.BP.8.1<-unlist(Modular.Score(W.100.GO.BP.8,TopProcesses=1))
W.100.mod.BP.12.1<-unlist(Modular.Score(W.100.GO.BP.12,TopProcesses=1))
W.100.mod.BP.16.1<-unlist(Modular.Score(W.100.GO.BP.16,TopProcesses=1))

save(W.100.mod.BP.4.1,
     W.100.mod.BP.8.1,W.100.mod.BP.12.1,W.100.mod.BP.16.1,
     file="~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP/WGCNA.100.BP.Modscore.rda")

#######################
W.100.dat.BP<-c(W.100.mod.BP.4.1["Model.score"],W.100.mod.BP.8.1["Model.score"],W.100.mod.BP.12.1["Model.score"],
         W.100.mod.BP.16.1["Model.score"])
         
save(W.100.dat.BP,file="~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP/W.100.dat.BP.rda")

par(mfrow=c(1,2))
plot(W.100.dat.BP,ylim=c(min(W.100.dat.BP),max(W.100.dat.BP)),type="o",xlab="Number of Modules",xaxt = "n",ylab="Model score", main="BP")
axis(1, at=1:length(k), labels=c(k))

####################################
## Top 1 GO process for each module
png(file = "~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP/Top.modularscore.BP.png", width = 1500, height = 1000, res=120);
par(mar=c(7.1, 5.1, 4.1, 2.1))
dat.1<-c(W.100.mod.BP.4.1,W.100.mod.BP.8.1,W.100.mod.BP.12.1,W.100.mod.BP.16.1)
w.1 = c(rep(2,length(W.100.mod.BP.4.1)-1),4,rep(2,length(W.100.mod.BP.8.1)-1),4,rep(2,length(W.100.mod.BP.12.1)-1),4,rep(2,length(W.100.mod.BP.16.1)-1),4)
space.1=c(rep(0.2,length(W.100.mod.BP.4.1)),1,rep(0.2,length(W.100.mod.BP.8.1)-1),1,rep(0.2,length(W.100.mod.BP.12.1)-1),1,rep(0.2,length(W.100.mod.BP.16.1)-2))
col.1=c(rep("blue",length(W.100.mod.BP.4.1)-1),"red",rep("blue",length(W.100.mod.BP.8.1)-1),"red",rep("blue",length(W.100.mod.BP.12.1)-1),"red",
        rep("blue",length(W.100.mod.BP.16.1)-1),"red")

bp<-barplot(dat.1,width=w.1, space=space.1, main="Top 1 GO:Biological Process - WGCNA",
            col=col.1, las=2,ylab=expression(paste(-log[10],"(p-value)")),cex.axis=1.5,cex.lab=1.7,cex.main=2,ylim=c(0,max(as.numeric(dat.1))+0.2))
text(bp, dat.1, labels = round(dat.1, digit=3), pos=3, cex=0.9)

legend("topright", legend=c("Module score","Model score"), fill=c("blue","red"), bty = "n",cex=1.5)

line=-37
cex=1.60
add=3
mtext(expression(bold("4 Modules")),at=7,line=line,cex=cex)
mtext(expression(bold("8 Modules")),at=23,line=line,cex=cex)
mtext(expression(bold("12 Modules")),at=40+add,line=line,cex=cex)
mtext(expression(bold("16 Modules")),at=55+add,line=line,cex=cex)
dev.off()


###################################################################
# function to calculating number of enriched GO terms

cutoff<-0.05

W.100.BP.logp.4<-GO.enrich(GO.BP=W.100.GO.BP.4,cutoff=cutoff,points=15)
W.100.BP.logp.8<-GO.enrich(GO.BP=W.100.GO.BP.8,cutoff=cutoff,points=15)
W.100.BP.logp.12<-GO.enrich(GO.BP=W.100.GO.BP.12,cutoff=cutoff,points=15)
W.100.BP.logp.16<-GO.enrich(GO.BP=W.100.GO.BP.16,cutoff=cutoff,points=15)

save(W.100.BP.logp.4,W.100.BP.logp.8,W.100.BP.logp.12,W.100.BP.logp.16,
     file="~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP/WGCNA.100.GO.annot.BP.rda")

png(file = "~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP/GO.BP_enriched.png", 
    width = 1000, height = 1000, res=200);

plot(W.100.BP.logp.4,col="blue",pch=4, type="o",ylab="No of enriched GO terms",
     xlab=expression(paste(-log[10],"(p-value)")),ylim=c(0,max(W.100.BP.logp.4)),main="GO:BP-WGCNA")

points(W.100.BP.logp.8,col="red",pch=16, type="o",ylab="No of enriched GO terms",
       xlab=expression(paste(-log[10],"(p-value)")))

points(W.100.BP.logp.12,col="blue",pch=12,type="o", ylab="No of enriched GO terms",
       xlab=expression(paste(-log[10],"(p-value)")))

points(W.100.BP.logp.16,col="black",pch=23, type="o",ylab="No of enriched GO terms",
       xlab=expression(paste(-log[10],"(p-value)")))

legend("topright",c("Modules-4","Modules-8","Modules-12","Modules-16"),pch=c(6,16,12,23),col=c("blue","red","blue","black"))

dev.off()

############################
## percentage of clusters that had atleast one enriched GO term
############################
## percentage of clusters that had atleast one enriched GO term
############################
## percentage of clusters that had atleast one enriched GO term BP
perc.W.100.BP.4<-perc.GO(GO.BP=W.100.GO.BP.4,cutoff=cutoff,points=15,clust=W.100.clust.4)
perc.W.100.BP.8<-perc.GO(GO.BP=W.100.GO.BP.8,cutoff=cutoff,points=15,clust=W.100.clust.8)
perc.W.100.BP.16<-perc.GO(GO.BP=W.100.GO.BP.16,cutoff=cutoff,points=15,clust=W.100.clust.16)
perc.W.100.BP.12<-perc.GO(GO.BP=W.100.GO.BP.12,cutoff=cutoff,points=15,clust=W.100.clust.12)

save(perc.W.100.BP.4,perc.W.100.BP.8,perc.W.100.BP.12,perc.W.100.BP.16,
     file="~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP/WGCNA.100.Perc..BP.rda")

png(file = "~/Dropbox/Phd/R/Modularity/Simulated/N100/WGCNA/BP/Perc_modules_BP.png", width = 1000, height = 1000, res=200);

plot(perc.W.100.BP.4,col="blue",pch=6, type="o", ylab="% annotated Modules",
       xlab=expression(paste(-log[10],"(p-value)")),ylim=c(0,100))

points(perc.W.100.BP.8,col="red",pch=16, type="o", ylab="% annotated Modules",
       xlab=expression(paste(-log[10],"(p-value)")))

points(perc.W.100.BP.12,col="blue",pch=12,type="o", ylab="% annotated Modules",
       xlab=expression(paste(-log[10],"(p-value)")))

points(perc.W.100.BP.16,col="black",pch=23, type="o",ylab="% annotated Modules",
       xlab=expression(paste(-log[10],"(p-value)")))

legend("topright",c("Modules-4","Modules-8","Modules-12","Modules-16"),pch=c(6,16,22,23),col=c("blue","red","blue","black"))

dev.off()
