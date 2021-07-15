setwd("~/Dropbox/Phd/R/Modularity/Real/SIMoNe")

require(simone)
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
        edges.steady="VAR(1) inference"
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

simone.clust.4<-simone.clust(X,clusters=4,type="time-course")
simone.clust.8<-simone.clust(X,clusters=8,type="time-course")
simone.clust.12<-simone.clust(X,clusters=12,type="time-course")
simone.clust.16<-simone.clust(X,clusters=16,type="time-course")            

save(simone.clust.4,simone.clust.8,simone.clust.12,simone.clust.16,file="simone.clust.n500.rda")
load("simone.clust.n500.rda")

require(clusterProfiler)

## comapare cluster for biological process
GO.BP.4 <- compareCluster(simone.clust.4, ont="BP", organism="yeast", pvalueCutoff=0.05)
GO.BP.8 <- compareCluster(simone.clust.8, ont="BP", organism="yeast", pvalueCutoff=0.05)
GO.BP.12 <- compareCluster(simone.clust.12, ont="BP", organism="yeast", pvalueCutoff=0.05)
GO.BP.16 <- compareCluster(simone.clust.16, ont="BP", organism="yeast", pvalueCutoff=0.05)

png(file = "clusters_BP.4.png", width = 1000, height = 700, r=120);
plot(GO.BP.4,title =" GO enrichment Biological Process-SIMoNe",font.size=10)
dev.off()

png(file = "clusters_BP.8.png", width = 1000, height = 700, r=120);
plot(GO.BP.8,title =" GO enrichment Biological Process-SIMoNe",font.size=10,showCategory=5)
dev.off()

png(file = "clusters_BP.12.png", width = 1000, height = 700, r=120);
plot(GO.BP.12,title =" GO enrichment Biological Process-SIMoNe",font.size=10)
dev.off()

png(file = "clusters_BP.16.png", width = 1000, height = 700, r=120);
plot(GO.BP.16,title =" GO enrichment Biological Process-SIMoNe",font.size=10)
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

barplot(dat.1,width=w.1, space=space.1, main="Top 1 GO process- RedeR",
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

perc.mod.4<-perc.GO(GO.BP=GO.BP.4,cutoff=cutoff,points=15,clust=simone.clust.4)
perc.mod.8<-perc.GO(GO.BP=GO.BP.8,cutoff=cutoff,points=15,clust=simone.clust.8)
perc.mod.16<-perc.GO(GO.BP=GO.BP.16,cutoff=cutoff,points=15,clust=simone.clust.16)
perc.mod.12<-perc.GO(GO.BP=GO.BP.12,cutoff=cutoff,points=15,clust=simone.clust.12)

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


