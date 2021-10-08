setwd("~/Dropbox/Phd/R/Consensus/Qualitative/Syntren/N100/10")
 
exp.path<-("~/Dropbox/Phd/R/Consensus/Simulated/syntren/N100/samples/10/Noise1/nn100_nbgr0_hop0.3_bionoise0.1_expnoise0.1_corrnoise0.1_neighAdd_dataset.txt")
gs.path<-("~/Dropbox/Phd/R/Consensus/Simulated/syntren/N100/samples/10/Noise1/nn100_nbgr0_hop0.3_bionoise0.1_expnoise0.1_corrnoise0.1_neighAdd_network.sif")
exp.data<-(as.matrix(read.table(exp.path, header = TRUE, sep = "\t", as.is=TRUE)))
gs.data<-(as.matrix(read.table(gs.path, header = F, sep = "\t", as.is=TRUE)))
exp.data<-t(exp.data)
dim(exp.data)
dim(gs.data)

 
load("~/Dropbox/Phd/R/Consensus/Simulated/syntren/N100/samples/10/Noise1/cons.fscore_goldstd.rda")
DF<-consensus.fscore[,1:7]
head(DF)


##################################################
overlap.ratio<-function(DF,top.perc) 
   	{
   	### number of edges to extract
   	no.edges<-(top.perc*dim(DF)[1])/100
   	
   	### consensys FCPT number of edges to extract
   	fcpt.edges<-DF[order(DF$F.score),]$edge[1:no.edges]
   	
   	#I ndividual algorithm top edges
   	reder.edges<-DF[order(DF$reder.pval),]$edge[1:no.edges]
   	wgcna.edges<-DF[order(DF$wgcna.pval),]$edge[1:no.edges]
  	aracne.edges<-DF[order(DF$aracne.pval),]$edge[1:no.edges]
   	clr.edges<-DF[order(DF$clr.pval),]$edge[1:no.edges]
    mrnet.edges<-DF[order(DF$mrnet.pval),]$edge[1:no.edges]

    # overlap ratio or similarity index.
   	reder.sim<-round(length(intersect(fcpt.edges,reder.edges))/no.edges,digit=2)
   	wgcna.sim<-round(length(intersect(fcpt.edges,wgcna.edges))/no.edges,digit=2)
   	aracne.sim<-round(length(intersect(fcpt.edges,aracne.edges))/no.edges,digit=2)
   	clr.sim<-round(length(intersect(fcpt.edges,clr.edges))/no.edges,digit=2)
   	mrnet.sim<-round(length(intersect(fcpt.edges,mrnet.edges))/no.edges,digit=2)

  	sim.idx<-data.frame(RedeR=reder.sim,WGCNA=wgcna.sim,ARACNE=aracne.sim,
  						CLR=clr.sim, MRNETB=mrnet.sim) 

	return(sim.idx)
	
	}

##################################################
### overlap ratio or similarity Index at different top percentages
simIndex.1<-overlap.ratio(DF,1)
simIndex.5<-overlap.ratio(DF,5)
simIndex.10<-overlap.ratio(DF,10)
simIndex.20<-overlap.ratio(DF,20)
simIndex.50<-overlap.ratio(DF,50)

Reder.100.10<-c(simIndex.1$RedeR,simIndex.5$RedeR,simIndex.10$RedeR,simIndex.20$RedeR,simIndex.50$RedeR)
WGCNA.100.10<-c(simIndex.1$WGCNA,simIndex.5$WGCNA,simIndex.10$WGCNA,simIndex.20$WGCNA,simIndex.50$WGCNA)
ARACNE.100.10<-c(simIndex.1$ARACNE,simIndex.5$ARACNE,simIndex.10$ARACNE,simIndex.20$ARACNE,simIndex.50$ARACNE)
CLR.100.10<-c(simIndex.1$CLR,simIndex.5$CLR,simIndex.10$CLR,simIndex.20$CLR,simIndex.50$CLR)
MRNETB.100.10<-c(simIndex.1$MRNETB,simIndex.5$MRNETB,simIndex.10$MRNETB,simIndex.20$MRNETB,simIndex.50$MRNETB)

save(Reder.100.10,WGCNA.100.10,ARACNE.100.10,CLR.100.10,MRNETB.100.10,file="OverlapRatio.100.10.rda")


### Unique edges
###############################################################
unique.edges<-function(DF,top.perc) 
   	{
   	### number of edges to extract
   	no.edges<-(top.perc*dim(DF)[1])/100
   	
   	### consensus FCPT number of edges to extract
   	fcpt.edges<-DF[order(DF$F.score),]$edge[1:no.edges]
   	
   	#Individual algorithm top edges
   	reder.edges<-DF[order(DF$reder.pval),]$edge[1:no.edges]
   	wgcna.edges<-DF[order(DF$wgcna.pval),]$edge[1:no.edges]
  	aracne.edges<-DF[order(DF$aracne.pval),]$edge[1:no.edges]
   	clr.edges<-DF[order(DF$clr.pval),]$edge[1:no.edges]
    mrnet.edges<-DF[order(DF$mrnet.pval),]$edge[1:no.edges]
    
    union.edges =(intersect(intersect(intersect(intersect(reder.edges,wgcna.edges),aracne.edges),clr.edges),mrnet.edges))
    cons.unique<-setdiff(fcpt.edges,union.edges)
    
    return(cons.unique)
    }
###############################################################
### Unique edges by consensus at different top percentages
consUnique.100.10.1<-length(unique.edges(DF,1))
consUnique.100.10.5<-length(unique.edges(DF,5))
consUnique.100.10.10<-length(unique.edges(DF,10))
consUnique.100.10.20<-length(unique.edges(DF,20))
consUnique.100.10.50<-length(unique.edges(DF,50))

save(consUnique.100.10.1,consUnique.100.10.5,consUnique.100.10.10,consUnique.100.10.20,consUnique.100.10.50,
    file="ConsUnique.100.10.rda")

###############################################################
## Sensitivity and specificity
#### performance 
### Gold standard - Truth

gs.data<-gsub("re","ac",gs.data)
gs.data <-data.frame(do.call('rbind', strsplit(as.character(gs.data),'ac',fixed=F)))
gs.100<-gs.data
#### extracting gold stadard network
gs.edge <- paste(gs.100[,1],gs.100[,2], sep = '-') 
gs.tp<-gsub(" ","", gs.edge , fixed=TRUE)

##################################################
sens.spec<-function(DF,gs.tp,top.perc,method)
	{
	### number of edges to extract
   	no.edges<-(top.perc*dim(DF)[1])/100
   	
   	if(method=='fcpt') {
   	   	   	print ("fcpt chosen")

   	### consensus FCPT number of edges to extract
   	predicted<-DF[order(DF$F.score),]$edge[1:no.edges]
   	
   	
   	} else if(method=='reder'){
   	   		print ("reder chosen")
   	predicted<-DF[order(DF$reder.pval),]$edge[1:no.edges]
   	
   	
   	}else if(method=='wgcna'){
   	   		print ("wgcna chosen")
	predicted<-DF[order(DF$wgcna.pval),]$edge[1:no.edges];
   	
   	
   	}else if(method=='aracne'){
   	   		print ("aracne chosen")

 	predicted<-DF[order(DF$aracne.pval),]$edge[1:no.edges];
  	
  	}else if (method=='clr'){
  	   		print ("clr chosen");
	predicted<-DF[order(DF$clr.pval),]$edge[1:no.edges];
   	
   	}else if (method=='mrnetb'){
   			print ("mrnetb chosen")
	predicted<-DF[order(DF$mrnet.pval),]$edge[1:no.edges];
    
   		} else 
   			print ('no method chosen');
   	  		
   	  	
	### Method 1 converting to binary vector as edgelist
	gs.m1<-matrix(nrow=length(predicted),ncol=2,byrow=TRUE)
	gs.m1[,1]<-as.matrix(predicted)
	gs.bin<-(gs.m1[,1] %in% gs.tp)+0

	### Sensitivity
	TP<-length(which(gs.bin==1))
	FP<-length(which(gs.bin==0))
	P<-length(gs.tp)
	FN<-(P-TP)
	sensitivity<-TP/(TP+FN)
	TN<-dim(DF)[1]-length(union(predicted,gs.tp))

	specificity<-TN/(FP+TN)

	return(c(sensitivity,specificity))
}

 ####################
 ## consensus
 	consSensSpec.100.10.1<-sens.spec(DF,gs.tp,1,method='fcpt')
	consSensSpec.100.10.5<-sens.spec(DF,gs.tp,5,method='fcpt')
	consSensSpec.100.10.10<-sens.spec(DF,gs.tp,10,method='fcpt')
	consSensSpec.100.10.20<-sens.spec(DF,gs.tp,20,method='fcpt')
	consSensSpec.100.10.50<-sens.spec(DF,gs.tp,50,method='fcpt')

 ## reder
	rederSensSpec.100.10.1<-sens.spec(DF,gs.tp,1,method='reder')
	rederSensSpec.100.10.5<-sens.spec(DF,gs.tp,5,method='reder')
	rederSensSpec.100.10.10<-sens.spec(DF,gs.tp,10,method='reder')
	rederSensSpec.100.10.20<-sens.spec(DF,gs.tp,20,method='reder')
	rederSensSpec.100.10.50<-sens.spec(DF,gs.tp,50,method='reder')

 ## wgcna
	wgcnaSensSpec.100.10.1<-sens.spec(DF,gs.tp,1,method='wgcna')
	wgcnaSensSpec.100.10.5<-sens.spec(DF,gs.tp,5,method='wgcna')
	wgcnaSensSpec.100.10.10<-sens.spec(DF,gs.tp,10,method='wgcna')
	wgcnaSensSpec.100.10.20<-sens.spec(DF,gs.tp,20,method='wgcna')
	wgcnaSensSpec.100.10.50<-sens.spec(DF,gs.tp,50,method='wgcna')


 ## aracne
	aracneSensSpec.100.10.1<-sens.spec(DF,gs.tp,1,method='aracne')
	aracneSensSpec.100.10.5<-sens.spec(DF,gs.tp,5,method='aracne')
	aracneSensSpec.100.10.10<-sens.spec(DF,gs.tp,10,method='aracne')
	aracneSensSpec.100.10.20<-sens.spec(DF,gs.tp,20,method='aracne')
	aracneSensSpec.100.10.50<-sens.spec(DF,gs.tp,50,method='aracne')

 ## clr
	clrSensSpec.100.10.1<-sens.spec(DF,gs.tp,1,method='clr')
	clrSensSpec.100.10.5<-sens.spec(DF,gs.tp,5,method='clr')
	clrSensSpec.100.10.10<-sens.spec(DF,gs.tp,10,method='clr')
	clrSensSpec.100.10.20<-sens.spec(DF,gs.tp,20,method='clr')
	clrSensSpec.100.10.50<-sens.spec(DF,gs.tp,50,method='clr')

   	
 ## mrnetb
	mrnetbSensSpec.100.10.1<-sens.spec(DF,gs.tp,1,method='mrnetb')
	mrnetbSensSpec.100.10.5<-sens.spec(DF,gs.tp,5,method='mrnetb')
	mrnetbSensSpec.100.10.10<-sens.spec(DF,gs.tp,10,method='mrnetb')
	mrnetbSensSpec.100.10.20<-sens.spec(DF,gs.tp,20,method='mrnetb')
	mrnetbSensSpec.100.10.50<-sens.spec(DF,gs.tp,50,method='mrnetb')
	
	
save(consSensSpec.100.10.1,
	consSensSpec.100.10.5,
	consSensSpec.100.10.10,
	consSensSpec.100.10.20,
	consSensSpec.100.10.50,
	rederSensSpec.100.10.1,
	rederSensSpec.100.10.5,
	rederSensSpec.100.10.10,
	rederSensSpec.100.10.20,
	rederSensSpec.100.10.50,
	wgcnaSensSpec.100.10.1,
	wgcnaSensSpec.100.10.5,
	wgcnaSensSpec.100.10.10,
	wgcnaSensSpec.100.10.20,
	wgcnaSensSpec.100.10.50,
	aracneSensSpec.100.10.1,
	aracneSensSpec.100.10.5,
	aracneSensSpec.100.10.10,
	aracneSensSpec.100.10.20,
	aracneSensSpec.100.10.50,
	clrSensSpec.100.10.1,
	clrSensSpec.100.10.5,
	clrSensSpec.100.10.10,
	clrSensSpec.100.10.20,
	clrSensSpec.100.10.50,
	mrnetbSensSpec.100.10.1,
	mrnetbSensSpec.100.10.5,
	mrnetbSensSpec.100.10.10,
	mrnetbSensSpec.100.10.20,
	mrnetbSensSpec.100.10.50, file="TopSensSpec.100.10.rda")
	
	
	   	
