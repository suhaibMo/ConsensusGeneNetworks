require(RedeR)
require(org.Sc.sgd.db)
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

######################################
res <- cea(exp.data, sig=1e-2, nper=1000, plotcea=T, ptype=4)

ceg <- graph.adjacency(res, diag=FALSE, mode="directed", weighted=TRUE)
edge.list<-get.edgelist(ceg) # interaction object
int.wt<-cbind(edge.list,E(ceg)$weight) # intercation object and edge strength
dim(int.wt)

#--- Step 7: generate the same RData object included in the package
yap1.deg <- list(dat=rownames(exp.data), exp=exp.data, ceg=ceg)

#--- Step 2: start the interface
rdp<-RedPort()
calld(rdp)

#--- Step 3: get a data frame, gene expression matrix and an igraph object
dt <- yap1.deg$dat
sg <- yap1.deg$ceg
gx <- yap1.deg$exp

#--- Step 4: map the data frame to the graph

#--- Step 5: set attributes to RedeR (i.e. gene symbols and fold change)

#--- Step 6: add graph to the app
addGraph(rdp,sg)

#...compute clustering on the gene expression matrix

hc<- hclust(dist(get.adjacency(sg,attr="weight")))
#...check results (e.g. significant clusters)
plot(hc)

#--- Step 8: superimpose significant gene modules onto the network as mapped in clustering analysis
nesthc(rdp,hc, cutlevel=0.5, metric="height", nmemb=2, cex=0.3)


#--- Step 9: assign edges t=3o containers
mergeOutEdges(rdp)

#--- Step 10: relax the network (p.s. fine-tune layout and container size interactively!)
relax(rdp)

#--- Step 11: add color and size legends

#--- reset graph

require(clusterProfiler)

## cutting tree to produce i groups of clusters 
tree<-cutree(hc,k=4)

## splitting tree for corresponding ORFs
clust<-split(names(tree), tree)
names(clust)<-sprintf("M%i",1:length(clust))

## comapring cluster for biological process
GO.MF <- compareCluster(clust, ont="MF", organism="yeast", pvalueCutoff=0.05)

png(file = "clusters_MF.png", width = 1000, height = 700, r=100);

#jpeg("test.jpg", width=100, height=100, units="in", res=500)

plot(GO.MF,by="percentage", title =" GO enrichment Molecular Function")
#p + scale_colour_gradient( low="blue", high= "yellow")
dev.off()

## comapring cluster for biological process
GO.BP <- compareCluster(clust, ont="BP", organism="yeast", pvalueCutoff=0.05)

png(file = "clusters_BP_test.png", width = 1000, height = 700, r=100);

plot(GO.BP,title =" GO enrichment Biological Process",font.size=30)
dev.off()

pdf(file = "clusters_BP_test.pdf", width = 10, height = 10)
par(cex = 1)
q<-plot(GO.BP,by="percentage", title =" GO enrichment Biological Process",font.size=30)
q + scale_colour_gradient( low="blue", high= "yellow")
mtext("Modules",side=1,col="black",line=2)
dev.off()


## comapring cluster for biological process
GO.CC <- compareCluster(clust, ont="CC", organism="yeast", pvalueCutoff=0.05)
png(file = "clusters_CC.png", width = 1000, height = 700, r=100);
p<-plot(GO.CC,by="percentage", title =" GO enrichment Cellular component")
p + scale_colour_gradient( low="blue", high= "yellow")
dev.off()

y <- enrichGO(gene=clust[[1]], organism="yeast", ont="MF", pvalueCutoff=0.01, readable=TRUE)
head(summary(y))
png(file = "cluster_1.png", width = 1000, height = 700, r=100);
plot(y)
dev.off()

KEGG <- compareCluster(clust, fun="enrichKEGG", organism="yeast", pvalueCutoff=0.05)
head(summary(KEGG))
png(file = "KEGG_1.png", width = 1000, height = 700, r=100);
p<-plot(KEGG,by="dot", title =" KEGG enrichment")
p + scale_colour_gradient( low="blue", high= "yellow")
dev.off()


#source("RedeR_data1_q_0.005_clustProf.R")
#  R CMD BATCH --no-save --no-verbose RedeR_data1_q_0.01_clustProf.R outputFile.Rout

