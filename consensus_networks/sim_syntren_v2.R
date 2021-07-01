### 
### date modified: 03/02/2015

setwd("~/Dropbox/Phd/R/Consensus/sampling_v2")

## this is the version used to generate results for thesis for 1000 perm
# R CMD BATCH --no-save --no-verbose sim_syntren_26Aug.R outputFile.Rout

exp.path<-("~/Dropbox/Phd/R/Consensus/Simulated/syntren/N100/samples/10/Noise1/nn100_nbgr0_hop0.3_bionoise0.1_expnoise0.1_corrnoise0.1_neighAdd_dataset.txt")
gs.path<-("~/Dropbox/Phd/R/Consensus/Simulated/syntren/N100/samples/10/Noise1/nn100_nbgr0_hop0.3_bionoise0.1_expnoise0.1_corrnoise0.1_neighAdd_network.sif")
exp.data<-(as.matrix(read.table(exp.path, header = TRUE, sep = "\t", as.is=TRUE)))
gs.data<-(as.matrix(read.table(gs.path, header = F, sep = "\t", as.is=TRUE)))
exp.data<-t(exp.data)
dim(exp.data)
dim(gs.data)

## function to converting matrix to edgelist
  matrix2edge <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- (!is.na(m))
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             p.val=m[ut])
}

## function to convert correlation matrix to p-values. See http://goo.gl/nahmV for documentation of this function
  cor.pval <- function (X, dfr) {
  R <-X
 # above <- row(R) < col(R)
 # r2 <- R[above]^2
  r2<- R^2
  Fstat <- r2 * dfr/(1 - r2)
  R <- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <- 0
  R
  }
  
######################## testing with correlation based network
 require(RedeR)
res <- cea(exp.data, sig=1, nper=1000, plotcea=F, ptype=4)
res[row(res) == col(res)] <- 1

reder.corr<-matrix2edge(res)

reder.pval<-cor.pval(X=res,dfr=ncol(exp.data)-2)

#################################
require(impute)
require(preprocessorCore)
require(WGCNA)

datExpr<-t(exp.data)


### Actual calculations 
softPower = 6;
wgcna.adjacency = adjacency(datExpr, type="signed", corFnc = "cor", corOptions = "use = 'p',method = 'pearson'",power =softPower);

###### new conversion back to original correlation value
wgcna2corr<-2*wgcna.adjacency^(1/softPower)-1

wgcna.corr <- graph.adjacency(wgcna2corr, diag=TRUE, mode="directed", weighted=TRUE)
wgcna.edge.list<-get.edgelist(wgcna.corr) # interaction object
edge.wgcna<-cbind(wgcna.edge.list,E(wgcna.corr)$weight) # intercation object and edge strength


## conversion of corr to pval
wgcna2corr.pval<-cor.pval(X=wgcna2corr,dfr=ncol(exp.data)-2)

#################

############# method 2  
require(minet)
require(parmigene)

############################################################
	# sampling algorithm for permuting expression vector of each gene with replacement
	# N represents number of radom sampling
	
 rand.sample<-function(exp.data,N)
    {
	colnames(exp.data)<-seq(1,dim(exp.data)[2],1)
	col<-ncol(exp.data)
	M<-nrow(exp.data)
 
	## mapping sampling columns to real expression data 
	 samp.col<-matrix(nrow=M,ncol=col)
	 samp.data<-list(samp.col)
 
 	for (i in 1:N)  
 	  {
  		for (j in 1:M)
  		  {
 			samp.col[j,]<-sample(exp.data[j,],replace=FALSE)
	 		rownames(samp.col)<-rownames(exp.data)
    	  } 
            samp.data[[i]]<-samp.col
 	  }
 	    return(samp.data)
 		
 	}
 		
 ############################################################
		
N=100
samp.data<-rand.sample(exp.data,N=N)


# running ARACNE algorithms with real data
eps=0.5
tau=0.9
k=3
nbins=10

real.mim <- round(knnmi.all(exp.data,k=k), digit=6)
real.rel <- round(build.mim(t(exp.data),estimator="mi.empirical",disc="equalwidth",nbins=nbins),digit=6)
#real.rel <- round(build.mim(t(exp.data),estimator="pearson"),digit=6)
aracne.real <- round(aracne(real.rel, eps=eps),digit=6)
clr.real <- round(clr(real.rel),digit=6)
mrnet.real <- round(mrnetb(real.rel), digit=6)
       
# running ARACNE algorithms N number times with sampled dataset
samp.mim<-list()
rel.net<-list()
aracne.net<-list()
clr.net<-list()
mrnet.net<-list()
reder.net<-list()

for(i in 1:N)
{
samp.mim[[i]] <- round(knnmi.all(samp.data[[i]],k=k),digit=6)
rel.net[[i]] <- round(build.mim(t(samp.data[[i]]),estimator="mi.empirical",disc="equalwidth",nbins=nbins),digit=6)
aracne.net[[i]] <- round(aracne(rel.net[[i]],eps=eps),digit=6)
clr.net[[i]] <- round(clr(rel.net[[i]]),digit=6)
mrnet.net[[i]] <- round(mrnetb(rel.net[[i]]),digit=6)
}


############################################################
# function to calculate p-val calculations using sampling data

rand.pval<-function(samp.freq,real.freq,N)
{
	rand.pval<-matrix(nrow=nrow(real.freq),ncol=ncol(real.freq))
	rownames(rand.pval)<-rownames(real.freq)
	colnames(rand.pval)<-colnames(real.freq)
	sum<-0

 for (i in 1:nrow(real.freq))
 {
   for (j in 1:ncol(real.freq))
   {
     for (k in 1:N)
     {
      	if (samp.freq[[k]][i,j]>=real.freq[i,j]) {
       	sum=sum+1
      }
    }
        rand.pval[i,j]<-sum/N
        sum<-0
   }
  }
  return(rand.pval)
}

#############################################################

aracne.pval<-rand.pval(aracne.net,aracne.real,N=N)
clr.pval<-rand.pval(clr.net,clr.real,N=N)
mrnet.pval<-rand.pval(mrnet.net,mrnet.real,N=N)



table(aracne.pval<0.05)
table(clr.pval<0.05)
table(mrnet.pval<0.05)

 png(file = "MI2pval_V2.png", width = 1000, height = 2000, res=200);
 par(mfrow=c(3,1))
 hist(aracne.pval, xlab="P val", main ="A. ARACNE p-val distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 hist(clr.pval, xlab="P val", main ="B. CLR p-val distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 hist(mrnet.pval, xlab="P val", main ="C. MRNET p-val distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 dev.off()
 

## Validation test
####################
### Actual calculations of correlations values

wgcna.adjacency = adjacency(datExpr, type="signed", corFnc = "cor", 
                           corOptions = "use = 'p',method = 'pearson'",power =softPower);

###### new conversion back to original correlation value
wgcna2corr<-abs(2*wgcna.adjacency^(1/softPower)-1)

### Acid test to verify if the permutation algorithm really works

wgcna.adjacency<-list()
wgcna.net<-list()

for(i in 1:N)
 {
    wgcna.adjacency[[i]] = round(adjacency(t(samp.data[[i]]), type="signed", corFnc = "cor", corOptions = "use = 'p',method = 'pearson'",power =softPower),digit=6)
    wgcna.net[[i]]<-abs(round((2*wgcna.adjacency[[i]]^(1/softPower)-1),digit=6))
 }

wgcna.pval<-rand.pval(wgcna.net,wgcna2corr,N=N)
wgcna.pval[row(wgcna.pval)==col(wgcna.pval)]<- 0

 png(file = "Corr_pval_V2.png", width = 1000, height = 1000, res=100);
 par(mfrow=c(2,2))
 hist(wgcna2corr.pval, xlab="P-val", main ="A. WGCNA p-val by R function",col="cornflowerblue",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)

 hist(wgcna.pval, xlab="P-val", main ="B. WGCNA p-val by New Permutation",col="cornflowerblue",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)

 plot(wgcna2corr.pval,wgcna.pval, main ="C. WGCNA R function vs New Permutation", xlab="p-val R-function",ylab="p-val New Permutation",
                cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 dev.off()


  
  #################### testing above pval computation
 png(file = "MI_pvalV2.png", width = 1000, height = 1000, res=100);
 par(mfrow=c(3,2))
 plot(abs(res),reder.pval, xlab="Correlation Coefficient",ylab="p val", main="RedeR")
 plot(abs(wgcna2corr),wgcna2corr.pval, xlab="Correlation Coefficient",ylab="p val", main="WGCNA")
 plot(aracne.real,aracne.pval, xlab="Mutual information",ylab="p val", main="ARACNE")
 plot(clr.real,clr.pval, xlab="Mutual information",ylab="p val", main="CLR")
 plot(mrnet.real,mrnet.pval, xlab="Mutual information",ylab="p val", main="MRNET")
 dev.off()
  ############################
	
######################################


edge_pval.reder<-matrix2edge(reder.pval)  
edge_pval.wgcna<-matrix2edge(wgcna2corr.pval) 
#edge_pval.rel<-matrix2edge(rel.pval)
edge_pval.aracne<-matrix2edge(aracne.pval)
edge_pval.clr<-matrix2edge(clr.pval)
edge_pval.mrnet<-matrix2edge(mrnet.pval)


 ##### merging all edges and weights in a data.frame
reder.list <- paste(edge_pval.reder[,1],edge_pval.reder[,2], sep = '_') 
reder.edgelist<-data.frame(edge=reder.list,reder.pval=as.numeric(edge_pval.reder[,3]))

wgcna.list<-paste(edge_pval.wgcna[,1],edge_pval.wgcna[,2], sep = '_')
wgcna.edgelist<-data.frame(edge=wgcna.list,wgcna.pval=as.numeric(edge_pval.wgcna[,3]))

#rel.list <- paste(edge_pval.rel[,1],edge_pval.rel[,2], sep = '_') 
#rel.edgelist<-data.frame(edge=rel.list,rel.pval=as.numeric(edge_pval.rel[,3]))

aracne.list <- paste(edge_pval.aracne[,1],edge_pval.aracne[,2], sep = '_') 
aracne.edgelist<-data.frame(edge=aracne.list,aracne.pval=as.numeric(edge_pval.aracne[,3]))

clr.list <- paste(edge_pval.clr[,1],edge_pval.clr[,2], sep = '_') 
clr.edgelist<-data.frame(edge=clr.list,clr.pval=as.numeric(edge_pval.clr[,3]))

mrnet.list <- paste(edge_pval.mrnet[,1],edge_pval.mrnet[,2], sep = '_') 
mrnet.edgelist<-data.frame(edge=mrnet.list,mrnet.pval=as.numeric(edge_pval.mrnet[,3]))

## combining all the edges: Naive method
temp1<-merge(reder.edgelist,wgcna.edgelist)
#temp2<-merge(temp1,rel.edgelist)
temp3<-merge(temp1,aracne.edgelist)
temp4<-merge(temp3,clr.edgelist)
DF<-merge(temp4,mrnet.edgelist)
  
  
  ####test

  
    ###################################################################################################
    ################################################################################################### 
    #                                                                                                 #
    #                                         Consensus network                                       #                 
    #                                                                                                 #
    ################################################################################################### 
 
red_wg_p<-DF[,2:(dim(DF)[2])]

fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)

F.score=matrix()

 for (i in 1:dim(red_wg_p)[1])
  {
   F.score[i]<-fishersMethod(red_wg_p[i,])
  }
  
###########################################
##################
RWQ.fscore<-cbind(DF[,1:dim(DF)[2]],F.score)
consensus.fscore<-RWQ.fscore
#consensus.fscore<-RWQ.fscore[order(RWQ.fscore[,7]),]
head(consensus.fscore)
  
  
 require(qvalue)

###calculating q values from fisher values
FQ.value<-qvalue(consensus.fscore$F.score, lambda = 0.05, fdr.level=0.05, pi0.method = "bootstrap")
Q.score<-cbind(FQ.value$qvalues, FQ.value$significant)

consensus.fscore<-data.frame(consensus.fscore,Q.score=Q.score[,1], sig=Q.score[,2])

####signficant edges
consensus.sig<-consensus.fscore[consensus.fscore$sig>0,]
head(consensus.sig)
dim(consensus.sig)	

  ## Reder WGCNA MI and F.score plots
 
 png(file = "pvals_MI_Fscores_all_prodV2.png", width = 1000, height = 1000, res=100);
 par(mfrow=c(3,2))

 plot(consensus.fscore$reder.pval, consensus.fscore$F.score, xlab="RedeR ", ylab="F score", main ="A. RedeR v F score",col="cornflowerblue",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 plot(consensus.fscore$wgcna.pval, consensus.fscore$F.score, xlab="WGCNA ", ylab="F score", main ="B. WGCNA v F score",col="cornflowerblue",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 plot(consensus.fscore$aracne.pval, consensus.fscore$F.score, xlab="ARACNE ", ylab="F score", main ="C. ARACNE v F score",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 #plot(consensus.fscore$rel.pval, consensus.fscore$F.score, xlab="RELEV ", ylab="F score", main ="C. RELEV v F score",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 plot(consensus.fscore$clr.pval, consensus.fscore$F.score, xlab="CLR ", ylab="F score", main ="D. CLR v F score",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 plot(consensus.fscore$mrnet.pval, consensus.fscore$F.score, xlab="MRNET ", ylab="F score", main ="E. MRNET v F score",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 dev.off()
 
 png(file = "pvals_dist_all_prodV2.png", width = 1000, height = 1000, res=100);
 par(mfrow=c(3,2))
 
 hist(consensus.fscore$reder.pval, xlab="P val", main ="A. RedeR p-val distribution",col="cornflowerblue",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 hist(consensus.fscore$wgcna.pval, xlab="P val", main ="B. WGCNA p-val distribution",col="cornflowerblue",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 hist(consensus.fscore$aracne.pval, xlab="P val", main ="C. ARACNE p-val distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 #hist(consensus.fscore$rel.pval, xlab="F score", main ="C. RELEV p-val distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 hist(consensus.fscore$clr.pval, xlab="P val", main ="D. CLR p-val distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 hist(consensus.fscore$mrnet.pval, xlab="P val", main ="E. MRNET p-val distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 hist(consensus.fscore$F.score, xlab="F score", main ="F. Fisher combined distribution",col="darkolivegreen",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
dev.off()
  
  
  png(file = "pvals_MI_V2.png", width = 1000, height = 1000, res=100);
 par(mfrow=c(3,3))
 
 hist(aracne.real, xlab="Mutual Information", main ="A. ARACNE MI distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 hist(clr.real, xlab="Mutual Information", main ="B. CLR p-val distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 hist(mrnet.real, xlab="Mutual Information", main ="C. MRNET p-val distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 
 plot(aracne.real,aracne.pval, xlab="Mutual information",ylab="p val", main="ARACNE",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 plot(clr.real,clr.pval, xlab="Mutual information",ylab="p val", main="CLR",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 plot(mrnet.real,mrnet.pval, xlab="Mutual information",ylab="p val", main="MRNET",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 
 hist(consensus.fscore$aracne.pval, xlab="P val", main ="G. ARACNE p-val distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 hist(consensus.fscore$clr.pval, xlab="P val", main ="H. CLR p-val distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
 hist(consensus.fscore$mrnet.pval, xlab="P val", main ="I. MRNET p-val distribution",col="brown3",cex.axis=1.6, cex.lab=1.55,cex.main=1.6)
dev.off()
  
 #####################################################
####################################################
#### calculating ROC, AUROC, PR and AUPR for different algorithms

consensus.edges <- data.frame(do.call('rbind', strsplit(as.character(consensus.fscore$edge),'_',fixed=TRUE)))
consensus.df<-data.frame(x1=consensus.edges[,1],x2=consensus.edges[,2],F.score=(consensus.fscore$F.score))
 
 ### #################################################
### Gold standard 50 genes


gs.data<-gsub("re","ac",gs.data)
gs.data <-data.frame(do.call('rbind', strsplit(as.character(gs.data),'ac',fixed=F)))
gs.100<-gs.data
#### extracting gold stadard network
gs.edge <- paste(gs.100[,1],gs.100[,2], sep = '_') 
gs.tp<-gsub(" ","", gs.edge , fixed=TRUE)

### Method 1 converting to binary vector as edgelist
gs.m1<-matrix(nrow=nrow(consensus.df),ncol=2,byrow=TRUE)
gs.m1[,1]<-as.matrix(consensus.fscore$edge)
gs.bin<-(gs.m1[,1] %in% gs.tp)+0

save(consensus.fscore,gs.bin, file="cons.fscore_goldstd.rda")

#################################################
################################################
######### Plotting ROC curve and area under ROC
require(ROCR)


pred.cons <- prediction((1-consensus.fscore$F.score),gs.bin)
pred.reder <- prediction((1-consensus.fscore$reder.pval),gs.bin)
pred.wgcna <- prediction((1-consensus.fscore$wgcna.pval),gs.bin)
pred.aracne <- prediction((1-consensus.fscore$aracne.pval),gs.bin)
#pred.rel <- prediction((1-consensus.fscore$rel.pval),gs.bin)
pred.clr <- prediction((1-consensus.fscore$clr.pval),gs.bin)
pred.mrnet <- prediction((1-consensus.fscore$mrnet.pval),gs.bin)


perf.cons <- performance(pred.cons,"tpr","fpr")
perf.reder <- performance(pred.reder,"tpr","fpr")
perf.wgcna <- performance(pred.wgcna,"tpr","fpr")
perf.aracne <- performance(pred.aracne,"tpr","fpr")
#perf.rel <- performance(pred.rel,"tpr","fpr")
perf.clr <- performance(pred.clr,"tpr","fpr")
perf.mrnet <- performance(pred.mrnet,"tpr","fpr")
#perf <- performance(pred, "prec", "rec")


png(file = "ROC_allV2.png", width = 800, height = 800, res=100);
# changing params for the ROC plot - width, etc
par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
# plotting the ROC curve

plot(perf.cons, lwd=4,col="black",main="Nodes=100, New permutation") 
plot(perf.reder,lwd=1, col="blue",add=TRUE) 
plot(perf.wgcna, lwd=1,col="green",add=TRUE)  
plot(perf.aracne, lwd=1,col="red",add=TRUE)
#plot(perf.rel, lwd=1,col="darkgreen",add=TRUE)
plot(perf.clr, lwd=1,col="brown",add=TRUE)  
plot(perf.mrnet, lwd=1,col="cornflowerblue",add=TRUE)  
abline(0,1,col="grey50",lty="dashed")

# calculating AUROC
auc.cons <- unlist(slot(performance(pred.cons,"auc"), "y.values"))
auc.reder <- unlist(slot(performance(pred.reder,"auc"), "y.values"))
auc.wgcna <- unlist(slot(performance(pred.wgcna,"auc"), "y.values"))
auc.aracne <- unlist(slot(performance(pred.aracne,"auc"), "y.values"))
#auc.rel <- unlist(slot(performance(pred.rel,"auc"), "y.values"))
auc.clr <- unlist(slot(performance(pred.clr,"auc"), "y.values"))
auc.mrnet <- unlist(slot(performance(pred.mrnet,"auc"), "y.values"))

# adding min and max ROC AUC to the center of the plot
maxauc.cons<-max(round(auc.cons, digits = 3))
maxauc.reder<-max(round(auc.reder, digits = 3))
maxauc.wgcna<-max(round(auc.wgcna, digits = 3))
maxauc.aracne<-max(round(auc.aracne, digits = 3))
#maxauc.rel<-max(round(auc.rel, digits = 3))
maxauc.clr<-max(round(auc.clr, digits = 3))
maxauc.mrnet<-max(round(auc.mrnet, digits = 3))
maxauct.cons <- paste(c("Consensus = "),maxauc.cons,sep="")
maxauct.reder <- paste(c("RedeR(Corr) = "),maxauc.reder,sep="")
maxauct.wgcna <- paste(c("WGCNA(Corr) = "),maxauc.wgcna,sep="")
maxauct.aracne <- paste(c("ARACNE(MI) = "),maxauc.aracne,sep="")
#maxauct.rel <- paste(c("RELEV(MI) = "),maxauc.rel,sep="")
maxauct.clr <- paste(c("CLR(MI) = "),maxauc.clr,sep="")
maxauct.mrnet <- paste(c("MRNET(MI) = "),maxauc.mrnet,sep="")
legend(0.7,0.6,c('Cons','RedeR','WGCNA','ARACNE','CLR','MRNET'),col=c('black','blue','green','red','brown','cornflowerblue'),lwd=2) 
legend(0.6,0.35,c(maxauct.cons,maxauct.reder,maxauct.wgcna,maxauct.aracne,maxauct.clr,maxauct.mrnet,"\n"),
,border="white",cex=0.9,box.col = "white")
dev.off()

vect1.roc<-c(maxauc.cons,maxauc.reder,maxauc.wgcna,maxauc.aracne,maxauc.clr,maxauc.mrnet)
size100.10auroc<-c(maxauc.cons,maxauc.reder,maxauc.wgcna,maxauc.aracne,maxauc.clr,maxauc.mrnet)
save(size100.10auroc,file="size100.10auroc.rda")

color<-c('gray15','cornflowerblue','cornflowerblue','brown3','brown3','brown3')
## Gene coverage
png(file = "~/Dropbox/Phd/R/Consensus/Simulated/syntren/N100/samples/10/Noise1/aurocV2.png", width = 800, height = 800, res=100);
barplot(vect1.roc*100, xlab=" ", ylab="AUROC(%)", 
col=color,ylim=c(0,100),cex.axis=1.6, cex.names=1.5, cex.lab=1.5,cex.main=1.5,axes = T)
axis(1, at=c(0.8,2,3.14,4.20,5.5,6.8,8),
labels=c("Consensus","RedeR","WGCNA","ARACNE","RELEV","CLR","MRNET"), las = 2, lwd=2.5)
legend("topright", c("Consensus","Correlation","Mutual Information"),cex=1,fill=color[c(1,2,4)])
abline(h=vect1.roc[1]*100,col="black",lty="dashed",lwd=2.5)
dev.off()

##### plotting Precision recall curve
pr.cons <- performance(pred.cons,"prec", "rec")
pr.reder <- performance(pred.reder,"prec", "rec")
pr.wgcna <- performance(pred.wgcna,"prec", "rec")
pr.aracne <- performance(pred.aracne,"prec", "rec")
#pr.rel <- performance(pred.rel,"prec", "rec")
pr.clr <- performance(pred.clr,"prec", "rec")
pr.mrnet <- performance(pred.mrnet,"prec", "rec")

############## fucntion to construct gold standard adjacency matrix
gold2adj<-function(matrix,target) {
adj<-(!is.na(matrix))+0-1

b<-matrix()
 for (i in 1:nrow(adj))
 {
   for (j in 1:ncol(adj))
   {
      b<-paste(rownames(adj)[i],colnames(adj)[j],sep = '_')
      if (b%in%target)
         adj[i,j]<-1
    }
    
 } 
 return(adj) 
 }
 
###### consensus pval symmetric matix from edgelist
cons.edgelist<-as.matrix(data.frame(consensus.fscore$edge,consensus.fscore$F.score))
cons.mat<-(!is.na(reder.pval))+0-1
b<-matrix()
 for (i in 1:nrow(cons.mat))
 {
   for (j in 1:ncol(cons.mat))
   {
      b<-paste(rownames(cons.mat)[i],colnames(cons.mat)[j],sep = '_')
      cons.mat[i,j]<-as.numeric(cons.edgelist[which(cons.edgelist==b),2])
    }
    
 } 
############## gold standard adjacency matrix
gs.adj<-gold2adj(matrix=reder.pval,target=gs.tp)
table(gs.adj)

###### consensus pval symmetric matix from edgelist
#cons.edge <- data.frame(do.call('rbind', strsplit(as.character(consensus.fscore$edge),'_',fixed=TRUE)))
#cons.edgelist<-cbind(cons.edge,consensus.fscore$F.score)
#cons.g<-graph.data.frame(cons.edge, directed = T)
#E(cons.g)$weight<-consensus.fscore$F.score
#cons.mat<-as.matrix(get.adjacency(cons.g, attr="weight"))


###### computing performance using minet functions
cons.table<-validate(1-cons.mat, gs.adj)
reder.table <- validate(1-reder.pval, gs.adj)
wgcna.table <- validate(1-wgcna2corr.pval, gs.adj)
#rel.table <- validate(1-rel.pval, gs.adj)
aracne.table <- validate(1-aracne.pval, gs.adj)
clr.table <- validate(1-clr.pval, gs.adj)
mrnet.table <- validate(1-mrnet.pval, gs.adj)

color<-c('black','blue','green','red','darkgreen','brown','cornflowerblue')

dev <- show.roc(cons.table, col=color[1], type="b")
dev <- show.roc(reder.table, col=color[2],device=dev, type="l")
dev <- show.roc(wgcna.table, device=dev, col=color[3], type="l")
#dev <- show.roc(rel.table, device=dev,col=color[4], type="l")
dev <- show.roc(aracne.table, device=dev, col=color[5], type="l")
dev <- show.roc(clr.table, device=dev, col=color[6], type="l")
show.roc(mrnet.table, device=dev, col=color[7], type="l")

dev <- show.pr(cons.table, col=color[1], type="b")
dev <- show.pr(reder.table, col=color[2],device=dev, type="l")
dev <- show.pr(wgcna.table, device=dev, col=color[3], type="l")
#dev <- show.pr(rel.table, device=dev,col=color[4], type="l")
dev <- show.pr(aracne.table, device=dev, col=color[5], type="l")
dev <- show.pr(clr.table, device=dev, col=color[6], type="l")
show.pr(mrnet.table, device=dev, col=color[7], type="l")

vect2.roc<-c(auc.roc(cons.table),auc.roc(reder.table),auc.roc(wgcna.table),auc.roc(aracne.table),auc.roc(clr.table)
,auc.roc(mrnet.table))
vect.pr<-c(auc.pr(cons.table),auc.pr(reder.table),auc.pr(wgcna.table),auc.pr(aracne.table),auc.pr(clr.table)
,auc.pr(mrnet.table))
vect1.roc
vect2.roc
vect.pr


color<-c('gray15','cornflowerblue','cornflowerblue','brown3','brown3','brown3')
## Gene coverage
png(file = "~/Dropbox/Phd/R/Consensus/Simulated/syntren/N100/samples/10/Noise1/auroc_prV2.png", width = 800, height = 800, res=100);
par(mfrow=c(2,1))
barplot(vect2.roc*100, xlab=" ", ylab="AUROC(%)", 
col=color,ylim=c(0,100),cex.axis=1.7, cex.names=1.5, cex.lab=1.5,cex.main=1.5)
axis(1, at=c(0.8,2,3.14,4.20,5.5,6.8),
labels=c("Consensus","RedeR","WGCNA","ARACNE","CLR","MRNET"), las = 2, lwd=2.5)
abline(h=vect2.roc[1]*100,col="black",lty="dashed",lwd=2.5)

## edge coverage
barplot(vect.pr*100, xlab=" ", ylab="AUPR(%)", 
col=color,ylim=c(0,15),cex.axis=1.7, cex.names=1.5, cex.lab=1.5,cex.main=1.5)
axis(1, at=c(0.8,2,3.14,4.20,5.5,6.8),
labels=c("Consensus","RedeR","WGCNA","ARACNE","CLR","MRNET"), las = 2, lwd=2.5)
abline(h=vect.pr[1]*100,col="black",lty="dashed",lwd=2.5)
legend("topright", c("Consensus","Correlation","Mutual Information"),cex=1,fill=color[c(1,2,4)])
dev.off()


### calculating AUPR

png(file = "~/Dropbox/Phd/R/Consensus/Simulated/syntren/N100/samples/10/Noise1/PR_allV2.png", width = 800, height = 800, res=100);
# changing params for the ROC plot - width, etc
par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
# plotting the ROC curve

plot(pr.cons, lwd=4,col="black",main="Nodes=50, synthetic network") 
plot(pr.reder,lwd=1, col="blue",add=TRUE) 
plot(pr.wgcna, lwd=1,col="green",add=TRUE)  
plot(pr.aracne, lwd=1,col="red",add=TRUE)
#plot(pr.rel, lwd=1,col="darkgreen",add=TRUE)
plot(pr.clr, lwd=1,col="brown",add=TRUE)  
plot(pr.mrnet, lwd=1,col="cornflowerblue",add=TRUE)  

# calculating AUPR

dev.off()
 
