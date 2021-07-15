setwd("~/Dropbox/Phd/R/Network_Inference/RedeR/Dataset1_20perc/q_0.01")
  
#xx<-as.data.frame(org.Sc.sgdCOMMON2ORF)
#x.orf<-xx[match(x,xx$gene_name),]
  
require(clusterProfiler)

load("HC_Anerobe_20perc_q0.005.RData")

## cutting tree to produce i groups of clusters 
#number of clusters
NumClust<-function(hc.data,clusters)
{
  tree<-cutree(hc.data,k=clusters)
  clust<-split(names(tree), tree)
  names(clust)<-sprintf("M%i",1:length(clust))
  return(clust)
}

clust.4<-NumClust(hc,4)
clust.8<-NumClust(hc,8)
clust.16<-NumClust(hc,16)
clust.32<-NumClust(hc,32)



## splitting tree for corresponding ORFs


## comapring cluster for biological process
GO.BP.4 <- compareCluster(clust.4, ont="BP", organism="yeast", pvalueCutoff=0.05)
GO.BP.8 <- compareCluster(clust.8, ont="BP", organism="yeast", pvalueCutoff=0.05)
GO.BP.16 <- compareCluster(clust.16, ont="BP", organism="yeast", pvalueCutoff=0.05)
GO.BP.32 <- compareCluster(clust.32, ont="BP", organism="yeast", pvalueCutoff=0.05)


M<-list()
add<-c()

for(i in 1:length(clust))
  {
      M[i]<-length(clust[[i]])
      add[i]<-(length(clust[[i]]))
  }

print(M)
print((add))
total<-sum(add)

 
###############################################
### Score function is used to extract top GO enriched features and stats 
### This function is part of clusterProfiler package
      
      
          score<-function(x,
                   type="dot",
                   title="",
                   font.size=12,
                   showCategory=5,
                   by="geneRatio",
                   colorBy="p.adjust") {

              clProf.df <- summary(x)

              ## get top 5 (default) categories of each gene cluster.
              if (is.null(showCategory)) {
                  result <- clProf.df
              } else {
                  Cluster <- NULL # to satisfy codetools
                  result <- ddply(.data = clProf.df,
                                  .variables = .(Cluster),
                                  .fun = function(df, N) {
                                      if (length(df$Count) > N) {
                                          idx <- order(df$Count, decreasing=T)[1:N]
                                          return(df[idx,])
                                      } else {
                                          return(df)
                                      }
                                  },
                                  N=showCategory
                                  )
              }

              ## remove zero count
              result$Description <- as.character(result$Description) ## un-factor
              GOlevel <- result[,c(2,3)] ## GO ID and Term
              GOlevel <- unique(GOlevel)

              result <- result[result$Count != 0, ]
              result$Description <- factor(result$Description,
                                           levels=rev(GOlevel[,2]))


              if (by=="rowPercentage") {
                  Description <- Count <- NULL # to satisfy codetools
                  result <- ddply(result,
                                  .(Description),
                                  transform,
                                  Percentage = Count/sum(Count),
                                  Total = sum(Count))

                  ## label GO Description with gene counts.
                  x <- mdply(result[, c("Description", "Total")], paste, sep=" (")
                  y <- sapply(x[,3], paste, ")", sep="")
                  result$Description <- y

                  ## restore the original order of GO Description
                  xx <- result[,c(2,3)]
                  xx <- unique(xx)
                  rownames(xx) <- xx[,1]
                  Termlevel <- xx[as.character(GOlevel[,1]),2]

                  ##drop the *Total* column
                  result <- result[, colnames(result) != "Total"]

                  result$Description <- factor(result$Description,
                                               levels=rev(Termlevel))

              } else if (by == "count") {
                  ## nothing
              } else if (by == "geneRatio") {
                  ## first try
                  ## cls <- as.character(result$Cluster)
                  ## clsu <- unique(cls)
                  ## idx <- sapply(clsu, function(i) which(i == cls)[1])
                  ## gcSize <- result$Count[idx]
                  ## names(gcSize) <- clsu
                  ## result$GeneRatio <- result$Count / gcSize[cls]
                  ## result$Cluster <- paste(cls, "(", gcSize[cls], ")", sep="")

                  ## second try
                  ## result$GeneRatio <- sapply(as.character(result$GeneRatio), function(i) eval(parse(text=i)))

                  ## final way
                  gsize <- as.numeric(sub("/\\d+$", "", as.character(result$GeneRatio)))
                  gcsize <- as.numeric(sub("^\\d+/", "", as.character(result$GeneRatio)))
                  result$GeneRatio = gsize/gcsize
                  result$Cluster <- paste(as.character(result$Cluster),"\n", "(", gcsize, ")", sep="")
              } else {
                  ## nothing
              }
             # p <- plotting.clusterProfile(result,
             #    							type=type,
             #                              by=by,
             #                              colorBy=colorBy,
             #                              title=title,
             #                              font.size=font.size)
              p<-result
              return(p)
          }
          

 require(plyr)

 TopProcesses<-5      
 GO.unlist<-score(GO.BP.4, by="count", showCategory= TopProcesses)  
 names(which(table(GO.unlist$Cluster)==TopProcesses))
 
 
 ###############################################
 ### modular score and model score function ####
 # This function computes modular score from GO enriched biological process 
 # or molecular function for each cluster of genes by extracting p-values associated for each GO term
 
 Modular.Score<-function(GO.BP,TopProcesses)
	{
	   ### Using score function from ClusterProfiler to extract enriched GO modules statistics
  	   GO.unlist<-score(GO.BP, by="count", showCategory= TopProcesses)  
  	   clust.top<-c(names(which(table(GO.unlist$Cluster)==TopProcesses)))
  	   print(paste("Modules that has top processes =",TopProcesses))
  	   print(clust.top)
        
       i.test<- gsub("M","",clust.top)

        ### transforming all p.value associated with GO with -log10(pvalue)
    	GO.unlist$logp<- -log10(GO.unlist$pvalue)
    
        ### Declaring objects
    	max.tmp<-matrix(nrow=0,ncol=length(clust.top));
    	sum.tmp<-matrix(nrow=0,ncol=length(clust.top));
    	module.score<-matrix(nrow=1,ncol=length(clust.top))
    	colnames(module.score)<-c(i.test)
        model.score<-c()
        max.tmp<-0
        sum.tmp<-0
 
   	for (i in 1:length(clust.top))
  	    	{  
          	  max.tmp<-max(GO.unlist[which(GO.unlist$Cluster==clust.top[i]),]$logp)
          	  sum.tmp<-sum(GO.unlist[which(GO.unlist$Cluster==clust.top[i]),]$logp)
 
                if (TopProcesses!=1){
                
           		    ### Module score calculation
           		    module.score[1,i]<-max.tmp/sum.tmp
           		    print(paste("module",clust.top[i],"=",module.score[1,i]))
           		    }
           		    else if (TopProcesses==1) {
           		    module.score[i]<-max.tmp
           		    print(paste("module",clust.top[i],"=",module.score[1,i]))

           		    }
           		    
           		  max.tmp<-0
                  sum.tmp<-0
                }
                
                  
   	       print(order(module.score, decreasing=T))
   	       print(paste("Model score=",sum(module.score)/length(module.score)))
  	       model.score <- sum(module.score)/length(module.score)
   
  	       scores<-data.frame(Module=module.score,Model.score=model.score)
  	 
  	  return(scores)
  	  
  	}
  	
  	
#######################################
  	   
            
            print(paste("module score=",module.score[i]))
 	# using p vlaues
	j=1;
	k=0
	tmp<-0;
	max.tmp<-matrix(nrow=0, ncol=clust.enrich);
	module.score<-matrix(nrow=0, ncol=clust.enrich);
	model.score<-c()
	p.k<-matrix()
	p.k <- (-log10(GO.unlist$pvalue))

 	for (i in i.test)
  	 {  	 
  	  for(j in j:(j+(TopProcesses-1)))
   	    {
   	     tmp <- tmp + p.k[j]
        }
	
    	 k= j-(TopProcesses-1)
     
        max.tmp[i]<-max(p.k[k:(k+(TopProcesses-1))])    

            if (TopProcesses!=1){
     
     	     	module.score[i]<- max.tmp[i]/tmp
     	    	}
     	    else {
     	    	 
     	   print(paste("k=",k))
     	   module.score[i]<-max.tmp[i]
     	   print(paste("module score=",module.score[i]))
  
            }
            
                   	tmp<-0;
     				j=j+1
    	 }
   
  	       print(module.score)
   	       print(order(module.score, decreasing=T))
           print(paste("Sum module scores=",sum(module.score)))
   	       print(paste("Model score=",sum(module.score)/length(module.score)))
  	       model.score <- sum(module.score)/length(module.score)
   
  	       scores<-list(Module=module.score,Model.score=model.score)
  	 
  	  return(scores)
   
   }
   
   #######################################################
   require(plyr)
   
  
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

  mod.16.5<-unlist(Modular.Score(GO.BP.16,TopProcesses=5))
  mod.16.4<-unlist(Modular.Score(GO.BP.16,TopProcesses=4))
  mod.16.3<-unlist(Modular.Score(GO.BP.16,TopProcesses=3))
  mod.16.2<-unlist(Modular.Score(GO.BP.16,TopProcesses=2))
  mod.16.1<-unlist(Modular.Score(GO.BP.16,TopProcesses=1))

  mod.32.5<-unlist(Modular.Score(GO.BP.32,TopProcesses=5))
  mod.32.4<-unlist(Modular.Score(GO.BP.32,TopProcesses=4))
  mod.32.3<-unlist(Modular.Score(GO.BP.32,TopProcesses=3))
  mod.32.2<-unlist(Modular.Score(GO.BP.32,TopProcesses=2))
  mod.32.1<-unlist(Modular.Score(GO.BP.32,TopProcesses=1))

  
   save(Model.score,m.score,mod.score,p.score,max.tmp,BP.score,file="Model_score_20perc_001_BP.RData")




####################################
## Top 1 GO process for each module
png(file = "~/Dropbox/Phd/R/Modularity/Top.modularscore.png", width = 1500, height = 1000, res=120);
par(mar=c(7.1, 5.1, 4.1, 2.1))
dat.1<-c(mod.8.1,mod.16.1,mod.32.1)
w.1 = c(rep(2,length(mod.8.1)-1),4,rep(2,length(mod.16.1)-1),4,rep(2,length(mod.32.1)-1),4)
space.1=c(rep(0.2,length(mod.8.1)),1,rep(0.2,length(mod.16.1)-1),1,rep(0.2,length(mod.32.1)-2))
col.1=c(rep("blue",length(mod.8.1)-1),"red",rep("blue",length(mod.16.1)-1),"red",
        rep("blue",length(mod.32.1)-1),"red")

barplot(dat.1,width=w.1, space=space.1, main="Top 1 GO process",
        col=col.1, las=2,ylab="-log10(p-value)",cex.axis=1.5,cex.lab=1.7,cex.main=2)
legend("topright", legend=c("Module score","Model score"), fill=c("blue","red"), bty = "n",cex=2)

line=-37
cex=1.60
mtext(expression(bold("8 clusters")),at=8,line=line,cex=cex)
mtext(expression(bold("16 clusters")),at=28,line=line,cex=cex)
mtext(expression(bold("32 clusters")),at=50,line=line,cex=cex)
dev.off()

####################################
## Top 5 to top 2 GO process for each module

png(file = "~/Dropbox/Phd/R/Modularity/Modularity_score.png", width = 2000, height = 1000, res=100);
par(mfrow=c(1,3),mar=c(12.1, 4.1, 4.1, 2.1))
module.8=c(mod.8.5,mod.8.4,mod.8.3, mod.8.2)
w.8 = c(rep(2,length(mod.8.5)-1),5)
space.8=c(rep(0.5,length(mod.8.5)),1,rep(0.5,length(mod.8.5)-1),1,rep(0.5,length(mod.8.5)-1),1)
barplot(module.8,las=2,space=space.8, width=w.8, cex.axis=2, cex.names=1.7, cex.lab=1.7,cex.main=3,
        col=c(rep("blue",length(mod.8.5)-1),"red"), main=" 8 clusters", ylim=c(0,0.7))
legend("topleft", legend=c("Module score","Model score"), fill=c("blue","red"),bty = "n",cex=3)


module.16=c(mod.16.5,mod.16.4,mod.16.3, mod.16.2)
w.16 = c(rep(2,length(mod.16.5)-1),5)
space.16=c(rep(0.5,length(mod.16.5)),1,rep(0.5,length(mod.16.5)-1),1,rep(0.5,length(mod.16.5)-1),1)
barplot(module.16,las=2,space=space.16,width=w.16,cex.axis=2, cex.names=1.7, cex.lab=1.7,cex.main=3,
        col=c(rep("blue",length(mod.16.5)-1),"red"),main="16 clusters",ylim=c(0,0.7))

module.32=c(mod.32.5,mod.32.4,mod.32.3, mod.32.2)
w.32 = c(rep(2,length(mod.32.5)-1),5)
space.32=c(rep(0.5,length(mod.32.5)),1,rep(0.5,length(mod.32.5)-1),1,rep(0.5,length(mod.32.5)-1),1)

barplot(module.32,las=2,space=space.32,width=w.32,cex.axis=2, cex.names=1.7, cex.lab=1.7,cex.main=3,
        col=c(rep("blue",length(mod.32.5)-1),"red"),main=" 32 clusters",ylim=c(0,0.7))

line=-71.5
cex=1.60
mtext(expression(bold("Top 5 GO")),at=-275,line=line,cex=cex)
mtext(expression(bold("Top 4 GO")),at=-245,line=line,cex=cex)
mtext(expression(bold("Top 3 GO")),at=-215,line=line,cex=cex)
mtext(expression(bold("Top 2 GO")),at=-185,line=line,cex=cex)

add<-145
mtext(expression(bold("Top 5 GO")),at=-275+add,line=line,cex=cex)
mtext(expression(bold("Top 4 GO")),at=-245+add,line=line,cex=cex)
mtext(expression(bold("Top 3 GO")),at=-215+add,line=line,cex=cex)
mtext(expression(bold("Top 2 GO")),at=-185+add,line=line,cex=cex)

add.2<-145
mtext(expression(bold("Top 5 GO")),at=-275+add+add.2,line=line,cex=cex)
mtext(expression(bold("Top 4 GO")),at=-245+add+add.2,line=line,cex=cex)
mtext(expression(bold("Top 3 GO")),at=-215+add+add.2,line=line,cex=cex)
mtext(expression(bold("Top 2 GO")),at=-185+add+add.2,line=line,cex=cex)

dev.off()

######################
# function to calculating number of enriched GO terms

GO.enrich<-function(GO.BP,cutoff,points)
  {
   sum.GO.BP<-summary(GO.BP)
   sum.GO.BP$logp<- -log10(sum.GO.BP$pvalue)
   p.seq<-seq(-log10(cutoff),max(sum.GO.BP$logp),length.out=points)

   GO.count<-matrix()
   log.p<-matrix()
   
    for (i in 1:length(p.seq)){
           GO.count[i]<-length(sum.GO.BP[sum.GO.BP$logp>p.seq[i],]$ID)
           log.p[i]<- p.seq[i]
      }
          return(cbind(log.p,GO.count))
  }

cutoff<-0.05
GO.logp.8<-GO.enrich(GO.BP=GO.BP.8,cutoff=cutoff,points=15)

GO.logp.16<-GO.enrich(GO.BP=GO.BP.16,cutoff=cutoff,points=15)

GO.logp.32<-GO.enrich(GO.BP=GO.BP.32,cutoff=cutoff,points=15)


plot(GO.logp.8,col="red",pch=16, type="o",ylab="No of enriched GO terms",
     xlab=expression(paste(-log[10],"(p-value)")),xlim=c(1,17))

points(GO.logp.16,col="black",pch=23, type="o",ylab="No of enriched GO terms",
     xlab=expression(paste(-log[10],"(p-value)")))

points(GO.logp.32,col="blue",pch=12,type="o", ylab="No of enriched GO terms",
     xlab=expression(paste(-log[10],"(p-value)")))

legend("topright",c("Clust.8","Clust.16","Clust.32"),pch=c(16,23,12),col=c("red","black","blue"))


############################
## percentage of clusters that had atleast one enriched GO term

sum.GO.BP.8<-summary(GO.BP.8)
perc.GO<-(length(table(sum.GO.BP[sum.GO.BP.8$pvalue<0.000005,]$Cluster))/length(clust.8))*100

perc.GO<-function(GO.BP,cutoff,points,clust)
{
  sum.GO.BP<-summary(GO.BP)
  sum.GO.BP$logp<- -log10(sum.GO.BP$pvalue)
  p.seq<-seq(-log10(cutoff),max(sum.GO.BP$logp),length.out=points)
  
  perc.mod<-matrix()
  log.p<-matrix()
  
  for (i in 1:length(p.seq)){
    perc.mod[i]<-(length(table(sum.GO.BP[sum.GO.BP$logp>p.seq[i],]$Cluster))/length(clust))*100
    log.p[i]<- p.seq[i]
  }
  return(cbind(log.p,perc.mod))
}

############################

  perc.mod.8<-perc.GO(GO.BP=GO.BP.8,cutoff=cutoff,points=15,clust=clust.8)
  perc.mod.16<-perc.GO(GO.BP=GO.BP.16,cutoff=cutoff,points=15,clust=clust.16)
  perc.mod.32<-perc.GO(GO.BP=GO.BP.32,cutoff=cutoff,points=15,clust=clust.32)


plot(perc.mod.8,col="red",pch=16, type="o", ylab="% annotated modules",
     xlab=expression(paste(-log[10],"(p-value)")),xlim=c(1,17),ylim=c(0,100))

points(perc.mod.16,col="black",pch=23, type="o",ylab="% annotated modules",
       xlab=expression(paste(-log[10],"(p-value)")))

points(perc.mod.32,col="blue",pch=12,type="o", ylab="% annotated modules",
       xlab=expression(paste(-log[10],"(p-value)")))

legend("topright",c("Clust.8","Clust.16","Clust.32"),pch=c(16,23,12),col=c("red","black","blue"))

