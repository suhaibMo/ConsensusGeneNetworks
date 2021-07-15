setwd("~/Dropbox/Phd/R/Modularity/Simulated/N100/RedeR/")
require(RedeR)
require(org.Sc.sgd.db)

source("~/Dropbox/Phd/R/Modularity/ModularityFunctions.R")
## 500 nodes simulated expression data
## 100 nodes simulated expression data
#exp.path<-("~/Documents/Softwares/Network/SyNTReN/data/results/size100nn100_nbgr0_hop0.3_bionoise0.1_expnoise0.1_corrnoise0.1_clustAdd_dataset.txt")
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


# orfs.path<-("~/Documents/Softwares/Network/SyNTReN/data/results/ORfs.100.clusAdd.csv")
# orfs<-(as.matrix(read.table(orfs.path, header = FALSE, sep = "\t")))
# rownames(exp.data)<-orfs[,1]
# head(exp.data)

rownames(exp.data)<-orfs
head(exp.data)

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
nesthc(rdp,hc, cutlevel=0.9, metric="height", nmemb=2, cex=0.3)

#--- Step 9: assign edges t=3o containers
mergeOutEdges(rdp)

#--- Step 10: relax the network (p.s. fine-tune layout and container size interactively!)
relax(rdp)

#--- Step 11: add color and size legends

#--- reset graph

require(clusterProfiler)

## cutting tree to produce i groups of clusters 
#number of clusters
NumClust<-function(hc.data,clusters)
          {
            ## cutting tree to produce i groups of clusters 
            tree<-cutree(hc.data,k=clusters)
            
            ## splitting tree for corresponding ORFs
            clust<-split(names(tree), tree)
            names(clust)<-sprintf("M%i",1:length(clust))
            return(clust)
          }

R.100.clust.4<-NumClust(hc,4)
R.100.clust.8<-NumClust(hc,8)
R.100.clust.12<-NumClust(hc,12)
R.100.clust.16<-NumClust(hc,16)

save(R.100.clust.4,R.100.clust.8,R.100.clust.12,R.100.clust.16,file="R.clusters.100.rda")


############################################################################################
require(fpc)
require(cluster)
R.100.mean.sil<-c()
R.100.avg.sil1<-c()
R.100.avg.sil2<-c()
R.100.dunn<-c()
R.100.dunn2<-c()
R.100.sep.idx<-c()

k=c(4,8,12,16)

for (i in 1:length(k)) {
  n=k
  print (i)
  print (n[i])
    dface<-dist(exp.data,method="euclidean")
    complete16<-cutree(hc,k=n[i])
    silwidths<-cluster.stats(d=dface,complete16,silhouette = TRUE)
    R.100.mean.sil[i]<-mean(as.numeric(silwidths$clus.avg.silwidths))
    R.100.avg.sil1[i]<-as.numeric(silwidths$avg.silwidth)
    R.100.dunn[i]<-as.numeric(silwidths$dunn)
    R.100.dunn2[i]<-as.numeric(silwidths$dunn2)
    R.100.sep.idx[i]<-as.numeric(silwidths$sindex)
  }

par(mfrow=c(3,3))
plot(R.100.mean.sil,ylim=c(min(R.100.mean.sil),max(R.100.mean.sil)),type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules',
     ylab="mean silhouette width")
axis(1, at=1:length(k), labels=c(k))

plot(R.100.avg.sil1,ylim=c(min(R.100.avg.sil1),max(R.100.avg.sil1)),type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules',
      ylab="Avg silhouette width")
 axis(1, at=1:length(k), labels=c(k))
 plot(R.100.dunn,ylim=c(min(R.100.dunn),max(R.100.dunn)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
      ylab="Dunn index")
 axis(1, at=1:length(k), labels=c(k))
plot(R.100.dunn2,ylim=c(min(R.100.dunn2),max(R.100.dunn2)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
     ylab="Dunn2 index")
axis(1, at=1:length(k), labels=c(k))

plot(R.100.sep.idx,ylim=c(min(R.100.sep.idx),max(R.100.sep.idx)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
     ylab="sep index")
axis(1, at=1:length(k), labels=c(k))

save(R.100.mean.sil,R.100.avg.sil1,R.100.dunn,R.100.dunn2,R.100.sep.idx,file="R.cluster.valid.100.rda")


## comapare cluster for biological process
R.100.GO.BP.4 <- compareCluster(R.100.clust.4, ont="BP", organism="yeast", pvalueCutoff=0.05)
R.100.GO.BP.8 <- compareCluster(R.100.clust.8, ont="BP", organism="yeast", pvalueCutoff=0.05)
R.100.GO.BP.12 <- compareCluster(R.100.clust.12, ont="BP", organism="yeast", pvalueCutoff=0.05)
R.100.GO.BP.16 <- compareCluster(R.100.clust.16, ont="BP", organism="yeast", pvalueCutoff=0.05)


setwd("~/Dropbox/Phd/R/Modularity/Simulated/N100/RedeR/BP/")

png(file = "clusters_BP.4.png", width = 1000, height = 700, r=120);
plot(R.100.GO.BP.4,title =" GO enrichment Biological Process - RedeR",font.size=10)
dev.off()

png(file = "clusters_BP.8.png", width = 1000, height = 700, r=120);
plot(R.100.GO.BP.8,title =" GO enrichment Biological Process - RedeR",font.size=10)
dev.off()

png(file = "clusters_BP.12.png", width = 1000, height = 700, r=120);
plot(R.100.GO.BP.12,title =" GO enrichment Biological Process - RedeR",font.size=10)
dev.off()

png(file = "clusters_BP.16.png", width = 1000, height = 700, r=120);
plot(R.100.GO.BP.16,title =" GO enrichment Biological Process - RedeR",font.size=10)
dev.off()




setwd("~/Dropbox/Phd/R/Modularity/Simulated/N100/RedeR")
################################################################

####
# biological process
R.100.mod.BP.4.1<-unlist(Modular.Score(R.100.GO.BP.4,TopProcesses=1))
R.100.mod.BP.8.1<-unlist(Modular.Score(R.100.GO.BP.8,TopProcesses=1))
R.100.mod.BP.12.1<-unlist(Modular.Score(R.100.GO.BP.12,TopProcesses=1))
R.100.mod.BP.16.1<-unlist(Modular.Score(R.100.GO.BP.16,TopProcesses=1))

save(R.100.mod.BP.4.1,
     R.100.mod.BP.8.1,R.100.mod.BP.12.1,R.100.mod.BP.16.1,
     file="~/Dropbox/Phd/R/Modularity/Simulated/N100/RedeR/BP/Reder.100.BP.Modscore.rda")



#######################
R.100.dat.BP<-c(R.100.mod.BP.4.1["Model.score"],R.100.mod.BP.8.1["Model.score"],
         R.100.mod.BP.12.1["Model.score"],R.100.mod.BP.16.1["Model.score"])


save(R.100.dat.BP,file="~/Dropbox/Phd/R/Modularity/Simulated/N100/RedeR/BP/R.100.dat.BP.rda")

par(mfrow=c(1,2))
plot(R.100.dat.BP,ylim=c(min(R.100.dat.BP),max(R.100.dat.BP)),type="o",xlab="Number of Modules",xaxt = "n",ylab="Model score", main="BP")
axis(1, at=1:length(k), labels=c(k))

####################################
## Top 1 GO process for each module
png(file = "~/Dropbox/Phd/R/Modularity/Simulated/N100/RedeR/BP/Top.modularscore.BP.png", width = 1500, height = 1000, res=120);
par(mar=c(7.1, 5.1, 4.1, 2.1))
dat.1<-c(R.100.mod.BP.4.1,R.100.mod.BP.8.1,R.100.mod.BP.12.1,R.100.mod.BP.16.1)
w.1 = c(rep(2,length(R.100.mod.BP.4.1)-1),4,rep(2,length(R.100.mod.BP.8.1)-1),4,rep(2,length(R.100.mod.BP.12.1)-1),4,rep(2,length(R.100.mod.BP.16.1)-1),4)
space.1=c(rep(0.2,length(R.100.mod.BP.4.1)),1,rep(0.2,length(R.100.mod.BP.8.1)-1),1,rep(0.2,length(R.100.mod.BP.12.1)-1),1,rep(0.2,length(R.100.mod.BP.16.1)-2))
col.1=c(rep("blue",length(R.100.mod.BP.4.1)-1),"red",rep("blue",length(R.100.mod.BP.8.1)-1),"red",rep("blue",length(R.100.mod.BP.12.1)-1),"red",
        rep("blue",length(R.100.mod.BP.16.1)-1),"red")

bp<-barplot(dat.1,width=w.1, space=space.1, main="GO:Biological Process - RedeR",
            col=col.1, las=2,ylab=expression(paste(-log[10],"(p-value)")),cex.axis=1.5,cex.lab=1.7,cex.main=2,ylim=c(0,max(as.numeric(dat.1))+0.2))
text(bp, dat.1, labels = round(dat.1, digit=3), pos=3, cex=0.9)

legend("topright", legend=c("Module score","Model score"), fill=c("blue","red"), bty = "n",cex=1.5)

line=-37
cex=1.60
add=0
mtext(expression(bold("4 Modules")),at=8,line=line,cex=cex)
mtext(expression(bold("8 Modules")),at=23+add,line=line,cex=cex)
mtext(expression(bold("12 Modules")),at=40+add,line=line,cex=cex)
mtext(expression(bold("16 Modules")),at=56+add,line=line,cex=cex)
dev.off()

###################################################################
# function to calculating number of enriched GO terms

cutoff<-0.05

R.100.BP.logp.4<-GO.enrich(GO.BP=R.100.GO.BP.4,cutoff=cutoff,points=15)
R.100.BP.logp.8<-GO.enrich(GO.BP=R.100.GO.BP.8,cutoff=cutoff,points=15)
R.100.BP.logp.12<-GO.enrich(GO.BP=R.100.GO.BP.12,cutoff=cutoff,points=15)
R.100.BP.logp.16<-GO.enrich(GO.BP=R.100.GO.BP.16,cutoff=cutoff,points=15)

save(R.100.BP.logp.4,R.100.BP.logp.8,R.100.BP.logp.12,R.100.BP.logp.16,
     file="~/Dropbox/Phd/R/Modularity/Simulated/N100/RedeR/BP/Reder.100.GO.annot.BP.rda")

png(file = "~/Dropbox/Phd/R/Modularity/Simulated/N100/RedeR/BP/GO.BP_enriched.png", width = 1000, height = 1000, res=200);

plot(R.100.BP.logp.4,col="blue",pch=6, type="o",ylab="No of enriched GO terms",
     xlab=expression(paste(-log[10],"(p-value)")),main="GO:BP-Reder")

points(R.100.BP.logp.8,col="red",pch=16, type="o",ylab="No of enriched GO terms",
       xlab=expression(paste(-log[10],"(p-value)")))

points(R.100.BP.logp.12,col="blue",pch=12,type="o", ylab="No of enriched GO terms",
       xlab=expression(paste(-log[10],"(p-value)")))

points(R.100.BP.logp.16,col="black",pch=23, type="o",ylab="No of enriched GO terms",
       xlab=expression(paste(-log[10],"(p-value)")))

legend("topright",c("Modules-4","Modules-8","Modules-12","Modules-16"),pch=c(6,16,12,23),col=c("blue","red","blue","black"))

dev.off()

############################
## percentage of clusters that had atleast one enriched GO term
perc.R.100.mod.4<-perc.GO(GO.BP=R.100.GO.BP.4,cutoff=cutoff,points=15,clust=R.100.clust.4)
perc.R.100.mod.8<-perc.GO(GO.BP=R.100.GO.BP.8,cutoff=cutoff,points=15,clust=R.100.clust.8)
perc.R.100.mod.16<-perc.GO(GO.BP=R.100.GO.BP.16,cutoff=cutoff,points=15,clust=R.100.clust.16)
perc.R.100.mod.12<-perc.GO(GO.BP=R.100.GO.BP.12,cutoff=cutoff,points=15,clust=R.100.clust.12)

save(perc.R.100.mod.4,perc.R.100.mod.8,perc.R.100.mod.12,perc.R.100.mod.16,
     file="Reder.100.Perc.rda")

############################
## percentage of clusters that had atleast one enriched GO term BP
perc.R.100.BP.4<-perc.GO(GO.BP=R.100.GO.BP.4,cutoff=cutoff,points=15,clust=R.100.clust.4)
perc.R.100.BP.8<-perc.GO(GO.BP=R.100.GO.BP.8,cutoff=cutoff,points=15,clust=R.100.clust.8)
perc.R.100.BP.16<-perc.GO(GO.BP=R.100.GO.BP.16,cutoff=cutoff,points=15,clust=R.100.clust.16)
perc.R.100.BP.12<-perc.GO(GO.BP=R.100.GO.BP.12,cutoff=cutoff,points=15,clust=R.100.clust.12)

save(perc.R.100.BP.4,perc.R.100.BP.8,perc.R.100.BP.12,perc.R.100.BP.16,
     file="~/Dropbox/Phd/R/Modularity/Simulated/N100/RedeR/BP/Reder.100.Perc..BP.rda")

png(file = "~/Dropbox/Phd/R/Modularity/Simulated/N100/RedeR/BP/Perc_modules_BP.png", width = 1000, height = 1000, res=200);

plot(perc.R.100.BP.4,col="blue",pch=6, type="o", ylab="% annotated Modules",
     xlab=expression(paste(-log[10],"(p-value)")),ylim=c(0,100))

points(perc.R.100.BP.8,col="red",pch=16, type="o", ylab="% annotated Modules",
       xlab=expression(paste(-log[10],"(p-value)")))

points(perc.R.100.BP.12,col="blue",pch=12,type="o", ylab="% annotated Modules",
       xlab=expression(paste(-log[10],"(p-value)")))

points(perc.R.100.BP.16,col="black",pch=23, type="o",ylab="% annotated Modules",
       xlab=expression(paste(-log[10],"(p-value)")))

legend("topright",c("Clust.4","Clust.8","Clust.12","Clust.16"),pch=c(6,16,22,23),col=c("blue","red","blue","black"))

dev.off()


