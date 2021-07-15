## This script contains all the functions used in extarcting GO enriched modules
## computing module score, model score  


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
  #      						type=type,
  #                              by=by,
  #                              colorBy=colorBy,
  #                              title=title,
  #                              font.size=font.size)
  p<-result
  return(p)
}

###############################################
### Normalised value
norm.scale<-function(x){
    #Normalized Data
    normalized = (x-min(x))/(max(x)-min(x))
    
    print("Orginial values")
    print(x)
    
    print("Normalised values")
    normalized
    return(normalized)
}

# This function output vector of cluster number associated to specific genes
clust.vector<-function(clust.list,exp.data){
	genes<-rownames(exp.data)
	vect<-c()
	for (i in 1:length(genes))
	  {
	    for (j in 1:length(clust.list)){
			 if (genes[i] %in% clust.list[[j]])
			     vect[i]<-j
	       }
	  } 
	  return(vect)
}

###############################################
### modular score and model score function ####
# This function computes modular score from GO enriched biological process 
# or molecular function for each cluster of genes by extracting p-values associated for each GO term

require(plyr)

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
    	GO.unlist$logp<- -log10(GO.unlist$p.adjust)
    
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
   	       #module.score<-norm.scale(module.score)
           model.score <- sum(module.score)/length(module.score)
   	       print(paste("Model score=",model.score))
   	
  	       scores<-data.frame(Module=module.score,Model.score=model.score)
  	 
  	  return(scores)
  	  
  	}

#######################################################

######################################################
# function to calculating number of enriched GO terms

GO.enrich<-function(GO.BP,cutoff,points)
    {
      sum.GO.BP<-summary(GO.BP)
      sum.GO.BP$logp<- -log10(sum.GO.BP$p.adjust)
      p.seq<-seq(-log10(cutoff),max(sum.GO.BP$logp),length.out=points)
      
      GO.count<-matrix()
      log.p<-matrix()
      
      for (i in 1:length(p.seq)){
              GO.count[i]<-length(sum.GO.BP[sum.GO.BP$logp>p.seq[i],]$ID)
              log.p[i]<- p.seq[i]
          }
            return(cbind(log.p,GO.count))
    }

######################################################
# function to % of annotated modules

perc.GO<-function(GO.BP,cutoff,points,clust)
          {
            sum.GO.BP<-summary(GO.BP)
            sum.GO.BP$logp<- -log10(sum.GO.BP$p.adjust)
            p.seq<-seq(-log10(cutoff),max(sum.GO.BP$logp),length.out=points)
            
            perc.mod<-matrix()
            log.p<-matrix()
            
            for (i in 1:length(p.seq)){
                  perc.mod[i]<-(length(table(sum.GO.BP[sum.GO.BP$logp>p.seq[i],]$Cluster))/length(clust))*100
                  log.p[i]<- p.seq[i]
            }
                 return(cbind(log.p,perc.mod))
          }
######################################################
# modules GO table for only the top scored modules
  
Modular.table<-function(GO.BP,TopProcesses,file)
  {
    ### Using score function from ClusterProfiler to extract enriched GO modules statistics
    GO.unlist<-score(GO.BP, by="count", showCategory= TopProcesses)  
    clust.top<-c(names(which(table(GO.unlist$Cluster)==TopProcesses)))
    print(paste("Processing GO table top processes =",TopProcesses))
    print(clust.top)
    
    i.test<- gsub("M","",clust.top)
    
    ### transforming all p.value associated with GO with -log10(pvalue)
    GO.unlist$logp<- -log10(GO.unlist$p.adjust)
    
    GO<-data.frame(Module=GO.unlist$Cluster,GO.ID=GO.unlist$ID, Process=GO.unlist$Description,geneID=GO.unlist$geneID, 
                   Count=GO.unlist$Count,P.value=signif(GO.unlist$p.adjust,5),Module.score=round(GO.unlist$logp,digit=3))
    print(GO)
    path=getwd()
    write.table(GO,file=paste(path,"/",file,".","xls",sep=""),sep="\t",col.names=NA)
  }


