#############
cons.network.plot<-function(edges,gs.tp,size,samples)
{
  ## network of predicted true positives
  pred.tp<-edges[which(edges %in% gs.tp)]
  pred.foo <- data.frame(do.call('rbind', strsplit(as.character(pred.tp),'-',fixed=TRUE)))
  g1.tp<- graph.data.frame(pred.foo, directed = T)
  pred.nodes<-length(V(g1.tp))
  print(paste("predicted network nodes = ",pred.nodes))
  pred.edges<-ecount(g1.tp)
  print(paste("predicted network edges = ",pred.edges))
  
  png(file = paste("consensus.network","size",size,"_.png",sep=""), width = 1000, height = 1000, res=100);
  l <- layout.kamada.kawai(g1.tp)
  l <- layout.norm(l,-1,1, -1,1)
  
  vertex.label.dist=0 
  vertex.label.cex=0.5 
  vertex.label.font=3
  edge.arrow.size=0.4
  vertex.size=6
  vertex.label.color='black'
  vertex.frame.color="gray50"
  vertex.color="deepskyblue1"
  
  plot(g1.tp,layout=l,vertex.size=vertex.size, vertex.color=vertex.color,
       edge.color="red", 
       main=paste("Predicted consensus network","-","Size-",size,",","samples-",samples,sep=" "),
       vertex.label=V(g1.tp)$name,
       vertex.size=vertex.size,
       vertex.label.dist=vertex.label.dist,
       vertex.label.cex=vertex.label.cex,
       vertex.label.font=vertex.label.font,
       edge.arrow.size=edge.arrow.size,
       vertex.size= vertex.size,
       vertex.label.color=vertex.label.color,
       vertex.frame.color=vertex.frame.color,
       rescale=F,xlim=range(l[,1]), ylim=range(l[,2]))
  legend("bottom",(paste("True edges = ",pred.edges)),cex=1)
  
  dev.off()
}


#######################################################

#######################################################
#######################################################
### true network and predicted network function
#######################################################
### true network and predicted network function
plot.network<-function(edges,gs.tp,method)
{
  ## network of true network
  gs.foo <- data.frame(do.call('rbind', strsplit(as.character(gs.tp),'-',fixed=TRUE)))
  g.tp<- graph.data.frame(gs.foo, directed = T)
  gs.nodes<-length(V(g.tp))
  print(paste("Gold standard network nodes = ",gs.nodes))
  gs.edges<-ecount(g.tp)
  print(paste("Gold standard network edges = ", gs.edges))
  
  ## network of predicted true positives
  pred.tp<-edges[which(edges %in% gs.tp)]
  pred.foo <- data.frame(do.call('rbind', strsplit(as.character(pred.tp),'-',fixed=TRUE)))
  g1.tp<- graph.data.frame(pred.foo, directed = T)
  pred.nodes<-length(V(g1.tp))
  print(paste("predicted network nodes = ",pred.nodes))
  pred.edges<-ecount(g1.tp)
  print(paste("predicted network edges = ",pred.edges))
  
  
  png(file = paste("network_tp",method,"_.png",sep=""), width = 1500, height = 1000, res=120);
  par(mfrow=c(1,2))
  set.seed(100)
  l <- layout.kamada.kawai(g.tp)
  l <- layout.norm(l,-1,1, -1,1)
  
  vertex.label.dist=0 
  vertex.label.cex=0.7 
  vertex.label.font=3
  edge.arrow.size=0.4
  vertex.size=13
  vertex.label.color='black'
  vertex.frame.color="gray50"
  vertex.color="deepskyblue1"
  
  plot(g.tp,layout=l,vertex.size=vertex.size, vertex.color=vertex.color,
       edge.color="black", 
       main="Gold standard network",
       vertex.label=V(g.tp)$name,
       vertex.size=vertex.size,
       vertex.label.dist=vertex.label.dist,
       vertex.label.cex=vertex.label.cex,
       vertex.label.font=vertex.label.font,
       edge.arrow.size=edge.arrow.size,
       vertex.size= vertex.size,
       vertex.label.color=vertex.label.color,
       vertex.frame.color=vertex.frame.color,
       rescale=F,xlim=range(l[,1]), ylim=range(l[,2]))
  legend("topleft", c("True interactions","Predicted true interactions"), col=c("black", "red"), lwd=c(2,2),cex=0.8,bty="n")
  legend("bottom",(paste("True edges = ",gs.edges)),cex=1)
  
  set.seed(100)
  m <- layout.kamada.kawai(g1.tp)
  m <- layout.norm(m,-1,1, -1,1)
  
  plot(g1.tp,layout=m,vertex.color=vertex.color,
       edge.color="red", 
       edge.arrow.size=edge.arrow.size,
       main=paste("Predicted",method,"network",sep=" "),
       vertex.label=V(g1.tp)$name,
       vertex.size=vertex.size,
       vertex.label.dist=vertex.label.dist,
       vertex.label.cex=vertex.label.cex,
       vertex.label.font=vertex.label.font,
       edge.arrow.size=edge.arrow.size,
       vertex.size= vertex.size,
       vertex.label.color=vertex.label.color,
       vertex.frame.color=vertex.frame.color,
       rescale=F,xlim=range(l[,1]), ylim=range(l[,2]))
  legend("bottom",(paste("True edges = ",pred.edges)),cex=1)
  
  dev.off()
}

