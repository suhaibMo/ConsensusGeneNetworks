### date modified: 03/02/2015

setwd("~/Dropbox/Phd/R/Consensus/Qualitative/Syntren/N100/10")

# R CMD BATCH --no-save --no-verbose size500_syntren_INV2.R outputFileINV2.txt


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
  
load("~/Dropbox/Phd/R/Consensus/Simulated/syntren/N100/samples/10/Noise1/cons.fscore_goldstd.rda")
DF<-consensus.fscore[,1:8]
head(DF)

##### 
## FCPT - Fishers combined probability test 
## p<0.05

p.val<-0.05

Fisher.df<-subset(DF, DF$Q.score<p.val, select=c(edge,F.score))
dim(Fisher.df)

##### 
## Intersection - Common edges predicted
## pval0.05


##### function to intersection all the edges

Intesect<-function(DF,cut.off)
	{
	reder.df<-subset(DF, DF$reder.pval<cut.off, select=c(edge,reder.pval))
	wgcna.df<-subset(DF, DF$wgcna.pval<cut.off, select=c(edge,wgcna.pval))
	aracne.df<-subset(DF, DF$aracne.pval<cut.off, select=c(edge,aracne.pval))
	clr.df<-subset(DF, DF$clr.pval<cut.off, select=c(edge,clr.pval))
	mrnet.df<-subset(DF, DF$mrnet.pval<cut.off, select=c(edge,mrnet.pval))

	common.edges=(intersect(intersect(intersect(intersect(reder.df$edge,wgcna.df$edge),aracne.df$edge),clr.df$edge),mrnet.df$edge))
	
	return(common.edges)
	}

intersect.edges<-Intesect(DF,cut.off=p.val)

length(intersect.edges)

##### Union of edges
Union<-function(DF,cut.off)
	{
	reder.df<-subset(DF, DF$reder.pval<cut.off, select=c(edge,reder.pval))
	wgcna.df<-subset(DF, DF$wgcna.pval<cut.off, select=c(edge,wgcna.pval))
	aracne.df<-subset(DF, DF$aracne.pval<cut.off, select=c(edge,aracne.pval))
	clr.df<-subset(DF, DF$clr.pval<cut.off, select=c(edge,clr.pval))
	mrnet.df<-subset(DF, DF$mrnet.pval<cut.off, select=c(edge,mrnet.pval))

	union.edges =(union(union(union(union(reder.df$edge,wgcna.df$edge),aracne.df$edge),clr.df$edge),mrnet.df$edge))
	
	return(union.edges)
	}
	
union.edges<-Union(DF,cut.off=p.val)

length(union.edges)

############################################################

### random network


#### performance 

### Gold standard - Truth

gs.data<-gsub("re","ac",gs.data)
gs.data <-data.frame(do.call('rbind', strsplit(as.character(gs.data),'ac',fixed=F)))
gs.100<-gs.data
#### extracting gold stadard network
gs.edge <- paste(gs.100[,1],gs.100[,2], sep = '-') 
gs.tp<-gsub(" ","", gs.edge , fixed=TRUE)


sens.spec<-function(predicted,gs.tp)
	{
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
	
	### Specificity
	TN<-dim(DF)[1]-length(union(predicted,gs.tp))
	specificity<-TN/(FP+TN)

		return(c(sensitivity,specificity))
	}

#######################
Inter.100.10<-sens.spec(intersect.edges,gs.tp)
Union.100.10<-sens.spec(union.edges,gs.tp)
Fisher.100.10<-sens.spec(Fisher.df$edge,gs.tp)


#####################
# random network
require(igraph)
consensus.foo <- data.frame(do.call('rbind', strsplit(as.character(Fisher.df$edge),'-',fixed=TRUE)))

g1 <- graph.data.frame(consensus.foo, directed = T)
con.nodes<-length(V(g1))
con.nodes
con.edges<-ecount(g1)
con.edges

res<-rownames(exp.data)

rand.g1.edgelist<-list()
rand.sens.spec<-list()
sens.rand<-matrix()
spec.rand<-matrix()

for (i in 1:100)
{
  rand.g1 <- erdos.renyi.game(con.nodes,con.edges,directed=T, type="gnm")
  rand.edge.g1<-get.edgelist(rand.g1)
  source.rand.g1<-sample(res)[rand.edge.g1[,1]]
  target.rand.g1<-sample(res)[rand.edge.g1[,2]]
  edge.rand.g1<-cbind(source.rand.g1,target.rand.g1)
  rand.g1.edgelist[[i]] <- paste(edge.rand.g1[,1],edge.rand.g1[,2], sep = '-') 
  
  rand.sens.spec[[i]]<-sens.spec(rand.g1.edgelist[[i]], gs.tp)
  sens.rand[i]<-rand.sens.spec[[i]][1]
  spec.rand[i]<-rand.sens.spec[[i]][2]
}

rand.100.10<-c(mean(sens.rand),mean(spec.rand))

save(Inter.100.10,Union.100.10,Fisher.100.10,rand.100.10, file ="SensSpec.100.10.rda")


#######################################################
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
  vertex.label.cex=0.6 
  vertex.label.font=3
  edge.arrow.size=0.4
  vertex.size=10
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

plot.network(Fisher.df$edge,gs.tp,method="consensus")
######################
#Overlap of edges between consensys and individual network inference



### Venn diagran

require(VennDiagram)

### edges 
t.value<-p.val
t.reder<-DF$edge[DF$reder.pval<t.value]
t.wgcna<-DF$edge[DF$wgcna.pval<t.value]
t.aracne<-DF$edge[DF$aracne.pval<t.value]
t.clr<-DF$edge[DF$clr.pval<t.value]
t.mrnet<-DF$edge[DF$mrnet.pval<t.value]


area1 = length(t.reder)
area2 = length(t.wgcna)
area3 = length(t.aracne)
area4 = length(t.clr)
area5 = length(t.mrnet)
n12 = length(intersect(t.reder,t.wgcna))
n13 = length(intersect(t.reder,t.aracne))
n14 = length(intersect(t.reder,t.clr))
n15 = length(intersect(t.reder,t.mrnet))
n23 = length(intersect(t.wgcna,t.aracne))
n24 = length(intersect(t.wgcna,t.clr))
n25 = length(intersect(t.wgcna,t.mrnet))
n34 = length(intersect(t.aracne,t.clr))
n35 = length(intersect(t.aracne,t.mrnet))
n45 = length(intersect(t.clr,t.mrnet))
n123 = length(intersect(intersect(t.reder,t.wgcna),t.aracne))
n124 = length(intersect(intersect(t.reder,t.wgcna),t.clr))
n125 = length(intersect(intersect(t.reder,t.wgcna),t.mrnet))
n134 = length(intersect(intersect(t.reder,t.aracne),t.clr))
n135 = length(intersect(intersect(t.reder,t.aracne),t.mrnet))
n145 = length(intersect(intersect(t.reder,t.clr),t.mrnet))
n234 = length(intersect(intersect(t.wgcna,t.aracne),t.clr))
n235 = length(intersect(intersect(t.wgcna,t.aracne),t.mrnet))
n245 = length(intersect(intersect(t.wgcna,t.clr),t.mrnet))
n345 = length(intersect(intersect(t.aracne,t.clr),t.mrnet))
n1234 = length(intersect(intersect(intersect(t.reder,t.wgcna),t.aracne),t.clr))
n1235 = length(intersect(intersect(intersect(t.reder,t.wgcna),t.aracne),t.mrnet))
n1245 = length(intersect(intersect(intersect(t.reder,t.wgcna),t.clr),t.mrnet))
n1345 = length(intersect(intersect(intersect(t.reder,t.aracne),t.clr),t.mrnet))
n2345 = length(intersect(intersect(intersect(t.wgcna,t.aracne),t.clr),t.mrnet))
n12345 = length(intersect(intersect(intersect(intersect(t.reder,t.wgcna),t.aracne),t.clr),t.mrnet))



# Reference five-set diagram
venn.plot <- draw.quintuple.venn(
area1 = area1,
area2 = area2,
area3 = area3,
area4 = area4,
area5 = area5,
n12 = n12,
n13 = n13,
n14 = n14,
n15 = n15,
n23 = n23,
n24 = n24,
n25 = n25,
n34 = n34,
n35 = n35,
n45 = n45,
n123 = n123,
n124 = n124,
n125 = n125,
n134 = n134,
n135 = n135,
n145 = n145,
n234 = n234,
n235 = n235,
n245 = n245,
n345 = n345,
n1234 = n1234,
n1235 = n1235,
n1245 = n1245,
n1345 = n1345,
n2345 = n2345,
n12345 = n12345,
category = c("RedeR", "WGCNA", "ARACNE", "CLR", "MRNETB"),
fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
cat.cex = 1.5,
margin = 0.15,
cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.5),
ind = TRUE);


png(file = "size100.10_Venn_diagramP0.05.png", width = 1000, height = 1000, res=200);
grid.draw(venn.plot);
dev.off();


