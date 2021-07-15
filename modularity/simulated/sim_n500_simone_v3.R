setwd("~/Dropbox/Phd/R/Modularity/Simulated/N500/SIMoNe")

require(simone)
require(org.Sc.sgd.db)

source("~/Dropbox/Phd/R/Modularity/ModularityFunctions.R")

## 100 nodes simulated expression data
## 500 nodes simulated expression data
exp.path<-("~/Dropbox/Phd/R/Consensus/Simulated/syntren/N500/100/500nn500_nbgr0_hop0.3_bionoise0.1_expnoise0.1_corrnoise0.1_neighAdd_dataset.txt")
exp.data<-(as.matrix(read.table(exp.path, header = TRUE, sep = "\t", as.is=TRUE)))
exp.data<-t(exp.data)
dim(exp.data)

####################################
## changing gene names to ORF's

orf.path<-("~/Dropbox/Phd/R/Modularity/Simulated/N500/gene.list.n500.csv")
orfs<-read.csv(orf.path,header=FALSE)
rownames(exp.data)<-orfs[,2]
head(exp.data)

#transpose the matrix
X<-t(exp.data)

simone.clust<-function(X,clusters,type)
    {
    if (type=="steady-state"){
        edges.steady="graphical.lasso";
        edges.sym.rule = "AND";
        }
        else if (type=="time-course"){
        edge.steady="VAR(1) inference"
        edges.sym.rule = "NO";
        }   
    
		## setting the control options for time series
		control <- setOptions(normalize = TRUE,
		verbose = TRUE,
		penalties = NULL,
		penalty.min = NULL,
		penalty.max = NULL,
		n.penalties = 100,
		edges.max = Inf,
		edges.sym.rule = edges.sym.rule,
		edges.steady = edges.steady,
		edges.coupling = "coopLasso",
		clusters.crit = "AIC",
		clusters.meth = "bayesian",
		clusters.qmin = clusters,
		clusters.qmax = clusters)


		## executing simone algorithm for time course data and clustering option
		res.cl <- simone(X, type=type, clustering=TRUE,control=control)

		###extracting network having maximum edges having lowest penalty score
		g.cl<-getNetwork(res.cl ,
		selection = res.cl$n.edges[length(res.cl$networks)],
		nodes = NULL)

		## saving the cluster plot 
		print(paste("Clusters plot saved at ",getwd()))
		png(file = paste("Network_circle",clusters,".png"), width = 1000, height = 800, res=100);
		plot(g.cl, type="circles",main = paste("SIMoNe network","clusters=",clusters))
		dev.off()

		print("cluster table")
		print(table(res.cl$clusters))
		print(paste("Number of cluster found =",length(table(res.cl$clusters))))
		

		ORF<-rownames(exp.data)

		## extracting genes
		ORF.clust<-cbind(ORF,res.cl$clusters,deparse.level = 1)

		### making emply list vector as per number of clusters
		clust.simone<-vector("list", length(unique(ORF.clust[,2])))
		names(clust.simone)<-unique(ORF.clust[,2])
		names(clust.simone)<-sprintf("M%i",1:length(clust.simone))

		## running loop to extract genes belonging to respective cluster number
		for(i in 1:length(unique(ORF.clust[,2])))
				{
					clust.simone[[i]]<-(as.vector(subset(ORF.clust, as.numeric(res.cl$clusters)==unique(ORF.clust[,2])[i], select=ORF)))
				}


         return(clust.simone)
    }

S.500.clust.4<-simone.clust(X,clusters=4,type="steady-state")
S.500.clust.8<-simone.clust(X,clusters=8,type="steady-state")
S.500.clust.12<-simone.clust(X,clusters=12,type="steady-state")
S.500.clust.16<-simone.clust(X,clusters=16,type="steady-state")            

load("/Users/suhaibmohammed/Dropbox/Phd/R/Modularity/Simulated/N500/SIMoNe/simone.clust.n500.rda")
S.500.clust.4<-simone.clust.4
S.500.clust.8<-simone.clust.8
S.500.clust.12<-simone.clust.12
S.500.clust.16<-simone.clust.16

save(S.500.clust.4,S.500.clust.8,S.500.clust.12,S.500.clust.16,file="S.500.clust.n500.rda")

simone.clust<-list(S.500.clust.4,S.500.clust.8,S.500.clust.12,S.500.clust.16)

############################################################################################
require(fpc)
require(cluster)
S.500.mean.sil<-c()
S.500.avg.sil1<-c()
S.500.avg.sil2<-c()
S.500.dunn<-c()
S.500.dunn2<-c()
S.500.sep.idx<-c()

k=c(4,8,12,16)

for (i in 1:length(k)) {
  n=k
  print (i)
  print (n[i])
  dface<-dist(exp.data,method="euclidean")
  complete16<-clust.vector(simone.clust[[i]],exp.data)
  silwidths<-cluster.stats(d=dface,complete16,silhouette = TRUE)
  S.500.mean.sil[i]<-mean(as.numeric(silwidths$clus.avg.silwidths))
  S.500.avg.sil1[i]<-as.numeric(silwidths$avg.silwidth)
  S.500.dunn[i]<-as.numeric(silwidths$dunn)
  S.500.dunn2[i]<-as.numeric(silwidths$dunn2)
  S.500.sep.idx[i]<-as.numeric(silwidths$sindex)
}

par(mfrow=c(3,3))
plot(S.500.mean.sil,ylim=c(min(S.500.mean.sil),max(S.500.mean.sil)),type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules',
     ylab="mean silhouette width")
axis(1, at=1:length(k), labels=c(k))

plot(S.500.avg.sil1,ylim=c(min(S.500.avg.sil1),max(S.500.avg.sil1)),type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules',
     ylab="Avg silhouette width")
axis(1, at=1:length(k), labels=c(k))
plot(S.500.dunn,ylim=c(min(S.500.dunn),max(S.500.dunn)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
     ylab="Dunn index")
axis(1, at=1:length(k), labels=c(k))
plot(S.500.dunn2,ylim=c(min(S.500.dunn2),max(S.500.dunn2)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
     ylab="Dunn2 index")
axis(1, at=1:length(k), labels=c(k))

plot(S.500.sep.idx,ylim=c(min(S.500.sep.idx),max(S.500.sep.idx)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
     ylab="sep index")
axis(1, at=1:length(k), labels=c(k))

save(S.500.mean.sil,S.500.avg.sil1,S.500.dunn,S.500.dunn2,S.500.sep.idx,file="S.cluster.valid.500.rda")


## comapare cluster for biological process
S.500.GO.BP.4 <- compareCluster(S.500.clust.4, ont="BP", organism="yeast", pvalueCutoff=0.05)
S.500.GO.BP.8 <- compareCluster(S.500.clust.8, ont="BP", organism="yeast", pvalueCutoff=0.05)
S.500.GO.BP.12 <- compareCluster(S.500.clust.12, ont="BP", organism="yeast", pvalueCutoff=0.05)
S.500.GO.BP.16 <- compareCluster(S.500.clust.16, ont="BP", organism="yeast", pvalueCutoff=0.05)


setwd("~/Dropbox/Phd/R/Modularity/Simulated/N500/SIMoNe/BP/")

png(file = "clusters_BP.4.png", width = 1000, height = 700, r=120);
plot(S.500.GO.BP.4,title =" GO enrichment Biological Process - SIMoNe",font.size=10)
dev.off()

png(file = "clusters_BP.8.png", width = 1000, height = 700, r=120);
plot(S.500.GO.BP.8,title =" GO enrichment Biological Process - SIMoNe",font.size=10)
dev.off()

png(file = "clusters_BP.12.png", width = 1000, height = 700, r=120);
plot(S.500.GO.BP.12,title =" GO enrichment Biological Process - SIMoNe",font.size=10)
dev.off()

png(file = "clusters_BP.16.png", width = 1000, height = 700, r=120);
plot(S.500.GO.BP.16,title =" GO enrichment Biological Process - SIMoNe",font.size=10)
dev.off()




setwd("~/Dropbox/Phd/R/Modularity/Simulated/N500/SIMoNe/BP/")
################################################################

####
# biological process
S.500.mod.BP.4.1<-unlist(Modular.Score(S.500.GO.BP.4,TopProcesses=1))
S.500.mod.BP.8.1<-unlist(Modular.Score(S.500.GO.BP.8,TopProcesses=1))
S.500.mod.BP.12.1<-unlist(Modular.Score(S.500.GO.BP.12,TopProcesses=1))
S.500.mod.BP.16.1<-unlist(Modular.Score(S.500.GO.BP.16,TopProcesses=1))

save(S.500.mod.BP.4.1,
     S.500.mod.BP.8.1,S.500.mod.BP.12.1,S.500.mod.BP.16.1,
     file="~/Dropbox/Phd/R/Modularity/Simulated/N500/SIMoNe/BP/SIMoNe.500.BP.Modscore.rda")



#######################
S.500.dat.BP<-c(S.500.mod.BP.4.1["Model.score"],S.500.mod.BP.8.1["Model.score"],
                S.500.mod.BP.12.1["Model.score"],S.500.mod.BP.16.1["Model.score"])


save(S.500.dat.BP,file="~/Dropbox/Phd/R/Modularity/Simulated/N500/SIMoNe/BP/S.500.dat.BP.rda")

par(mfrow=c(1,2))
plot(S.500.dat.BP,ylim=c(min(S.500.dat.BP),max(S.500.dat.BP)),type="o",xlab="Number of Modules",xaxt = "n",ylab="Model score", main="BP")
axis(1, at=1:length(k), labels=c(k))

####################################
## Top 1 GO process for each module
png(file = "~/Dropbox/Phd/R/Modularity/Simulated/N500/SIMoNe/BP/Top.modularscore.BP.png", width = 1600, height = 1000, res=120);
par(mar=c(7.1, 5.1, 4.1, 2.1))
dat.1<-c(S.500.mod.BP.4.1,S.500.mod.BP.8.1,S.500.mod.BP.12.1,S.500.mod.BP.16.1)
w.1 = c(rep(2,length(S.500.mod.BP.4.1)-1),4,rep(2,length(S.500.mod.BP.8.1)-1),4,rep(2,length(S.500.mod.BP.12.1)-1),4,rep(2,length(S.500.mod.BP.16.1)-1),4)
space.1=c(rep(0.2,length(S.500.mod.BP.4.1)),1,rep(0.2,length(S.500.mod.BP.8.1)-1),1,rep(0.2,length(S.500.mod.BP.12.1)-1),1,rep(0.2,length(S.500.mod.BP.16.1)-2))
col.1=c(rep("blue",length(S.500.mod.BP.4.1)-1),"red",rep("blue",length(S.500.mod.BP.8.1)-1),"red",rep("blue",length(S.500.mod.BP.12.1)-1),"red",
        rep("blue",length(S.500.mod.BP.16.1)-1),"red")

bp<-barplot(dat.1,width=w.1, space=space.1, main="GO:Biological Process - SIMoNe",
            col=col.1, las=2,ylab=expression(paste(-log[10],"(p-value)")),cex.axis=1.5,cex.lab=1.7,cex.main=2,ylim=c(0,max(as.numeric(dat.1))+0.3))
text(bp, dat.1, labels = round(dat.1, digit=3), pos=3, cex=0.9)

legend("topright", legend=c("Module score","Model score"), fill=c("blue","red"), bty = "n",cex=1.5)

line=-37
cex=1.60
add=17
mtext(expression(bold("4 Modules")),at=7,line=line,cex=cex)
mtext(expression(bold("8 Modules")),at=25,line=line,cex=cex)
mtext(expression(bold("12 Modules")),at=40+add,line=line,cex=cex)
mtext(expression(bold("16 Modules")),at=75+add,line=line,cex=cex)
dev.off()

###################################################################
# function to calculating number of enriched GO terms

cutoff<-0.05

S.500.BP.logp.4<-GO.enrich(GO.BP=S.500.GO.BP.4,cutoff=cutoff,points=15)
S.500.BP.logp.8<-GO.enrich(GO.BP=S.500.GO.BP.8,cutoff=cutoff,points=15)
S.500.BP.logp.12<-GO.enrich(GO.BP=S.500.GO.BP.12,cutoff=cutoff,points=15)
S.500.BP.logp.16<-GO.enrich(GO.BP=S.500.GO.BP.16,cutoff=cutoff,points=15)

save(S.500.BP.logp.4,S.500.BP.logp.8,S.500.BP.logp.12,S.500.BP.logp.16,
     file="~/Dropbox/Phd/R/Modularity/Simulated/N500/SIMoNe/BP/SIMoNe.500.GO.annot.BP.rda")

png(file = "~/Dropbox/Phd/R/Modularity/Simulated/N500/SIMoNe/BP/GO.BP_enriched.png", width = 1000, height = 1000, res=200);

plot(S.500.BP.logp.4,col="blue",pch=6, type="o",ylab="No of enriched GO terms",
     xlab=expression(paste(-log[10],"(p-value)")),main="GO:BP-SIMoNe")

points(S.500.BP.logp.8,col="red",pch=16, type="o",ylab="No of enriched GO terms",
       xlab=expression(paste(-log[10],"(p-value)")))

points(S.500.BP.logp.12,col="blue",pch=12,type="o", ylab="No of enriched GO terms",
       xlab=expression(paste(-log[10],"(p-value)")))

points(S.500.BP.logp.16,col="black",pch=23, type="o",ylab="No of enriched GO terms",
       xlab=expression(paste(-log[10],"(p-value)")))

legend("topright",c("Modules-4","Modules-8","Modules-12","Modules-16"),pch=c(6,16,12,23),col=c("blue","red","blue","black"))

dev.off()

############################
## percentage of clusters that had atleast one enriched GO term BP
perc.S.500.BP.4<-perc.GO(GO.BP=S.500.GO.BP.4,cutoff=cutoff,points=15,clust=S.500.clust.4)
perc.S.500.BP.8<-perc.GO(GO.BP=S.500.GO.BP.8,cutoff=cutoff,points=15,clust=S.500.clust.8)
perc.S.500.BP.16<-perc.GO(GO.BP=S.500.GO.BP.16,cutoff=cutoff,points=15,clust=S.500.clust.16)
perc.S.500.BP.12<-perc.GO(GO.BP=S.500.GO.BP.12,cutoff=cutoff,points=15,clust=S.500.clust.12)

save(perc.S.500.BP.4,perc.S.500.BP.8,perc.S.500.BP.12,perc.S.500.BP.16,
     file="~/Dropbox/Phd/R/Modularity/Simulated/N500/SIMoNe/BP/SIMoNe.500.Perc..BP.rda")

png(file = "~/Dropbox/Phd/R/Modularity/Simulated/N500/SIMoNe/BP/Perc_modules_BP.png", width = 1000, height = 1000, res=200);

plot(perc.S.500.BP.4,col="blue",pch=6, type="o", ylab="% annotated Modules",
     xlab=expression(paste(-log[10],"(p-value)")),ylim=c(0,100))

points(perc.S.500.BP.8,col="red",pch=16, type="o", ylab="% annotated Modules",
       xlab=expression(paste(-log[10],"(p-value)")))

points(perc.S.500.BP.12,col="blue",pch=12,type="o", ylab="% annotated Modules",
       xlab=expression(paste(-log[10],"(p-value)")))

points(perc.S.500.BP.16,col="black",pch=23, type="o",ylab="% annotated Modules",
       xlab=expression(paste(-log[10],"(p-value)")))

legend("topright",c("Modules-4","Modules-8","Modules-12","Modules-16"),pch=c(6,16,22,23),col=c("blue","red","blue","black"))

dev.off()


