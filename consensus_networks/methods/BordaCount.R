setwd("~/Dropbox/Phd/R/Consensus/Simulated/syntren/N10/10")

# R CMD BATCH --no-save --no-verbose sim_syntren_24March.R outputFile.txt

exp.path<-("~/Dropbox/Phd/R/Consensus/Simulated/syntren/N10/10/N10nn10_nbgr0_hop0.3_bionoise0.1_expnoise0.1_corrnoise0.1_neighAdd_dataset.txt")
gs.path<-("~/Dropbox/Phd/R/Consensus/Simulated/syntren/N10/10/N10nn10_nbgr0_hop0.3_bionoise0.1_expnoise0.1_corrnoise0.1_neighAdd_network.sif")
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
reder.corr[,3]<-abs(reder.corr[,3])

 ##### merging all edges and weights in a data.frame
reder.list <- paste(reder.corr[,1],reder.corr[,2], sep = '_') 
reder.edgelist<-data.frame(edge=reder.list,reder.corr=as.numeric(reder.corr[,3]))

#### sorting and ranking the edgelist based on confidence scores 
reder.sort<-reder.edgelist[order(reder.edgelist[,2], decreasing=T),]
reder.rank<-data.frame(reder.sort,rank=rank(-reder.sort[,2]))	



## borda count election method to average all ranks


#################################
require(WGCNA)
yap1.exp<-t(exp.data)

##Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious
##outliers. We use the function flashClust that provides faster hierarchical clustering than the standard function
##hclust
sampleTree = flashClust(dist(yap1.exp), method = "average");
print("clustering data..")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.

cut.height=65;
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut.height, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = yap1.exp[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##In the plot, shown in Fig. 2, white means a low value, red a high value, and grey a missing entry.
##The last step is to save the relevant expression and trait data for use in the next steps of the tutorial

print("Calculate softthreshold...")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=50, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


### Actual calculations 
softPower = 6;
wgcna.adjacency = adjacency(datExpr, type="signed", corFnc = "cor", corOptions = "use = 'p',method = 'pearson'",power =softPower);


###### new conversion back to original correlation value
wgcna2corr<-2*wgcna.adjacency^(1/softPower)-1


############# method 2  
require(minet)
require(parmigene)

# running ARACNE algorithms with real data
eps=0.5
tau=0.9
k=2
nbins=10

real.mim <- round(knnmi.all(exp.data,k=k), digit=6)
real.rel <- round(build.mim(t(exp.data),estimator="mi.empirical",disc="equalwidth",nbins=nbins),digit=6)
#real.rel <- round(build.mim(t(exp.data),estimator="pearson"),digit=6)
#aracne.real <- round(aracne.m(real.mim, tau=tau),digit=6)
aracne.real <- round(aracne(real.rel, eps=eps),digit=6)
clr.real <- round(clr(real.mim),digit=6)
mrnet.real <- round(mrnetb(real.mim), digit=6)


#### matrix to edgelist

matrix2edgelist<-function(edge.matrix)
{
edge.list<-matrix2edge(edge.matrix)
edge.list[,3]<-abs(edge.list[,3])

return(edge.list)
}

#### Borda count method sorts edges by ranks based on the edge confidence scores 
#### and then aggregates

SortRank<-function(edge.list) 
 {
    #### merging the source and target genes with confidence scores in a dataframe 
    edge.merge <- paste(edge.list[,1],edge.list[,2], sep = '_') 
    net.edgelist<-data.frame(edge=edge.merge,edge.value=as.numeric(edge.list[,3]))

    #### sorting and ranking the edgelist based on confidence scores 
    net.sort<-net.edgelist[order(net.edgelist[,2], decreasing=T),]
    net.rank<-data.frame(net.sort,rank=rank(-net.sort[,2]))
    
    return(net.rank)
 }

reder_edge.list<-edge.list<-matrix2edgelist(res)
reder.rank<-SortRank(reder_edge.list)

wgcna_edge.list<-edge.list<-matrix2edgelist(wgcna2corr)
wgcna.rank<-SortRank(wgcna_edge.list)

aracne_edge.list<-edge.list<-matrix2edgelist(aracne.real)
aracne.rank<-SortRank(aracne_edge.list)

clr_edge.list<-edge.list<-matrix2edgelist(clr.real)
clr.rank<-SortRank(clr_edge.list)

mrnet_edge.list<-edge.list<-matrix2edgelist(mrnet.real)
mrnet.rank<-SortRank(mrnet_edge.list)


#### Function for agregate of edge ranks ######## 
#### from all network inference algorithms  ###############
#### ##################

AverageRank<-function(net1,net2,net3,net4,net5)
 {
 if (dim(net1)[1]!=dim(net2)[1])
 print("the networks are not in same dimensions")

 else
 
    rank.agg<-matrix(nrow=nrow(net1),ncol=ncol(net1)-1)
    rank.agg[,1]<-as.character(net1$edge)
    colnames(rank.agg)<-c("Edge", "AverageRank")

#### sorting tables to have same order of edges from all dataframes

net1<-net1
net2<-net2[match(net1$edge,net2$edge),]
net3<-net3[match(net1$edge,net3$edge),]
net4<-net4[match(net1$edge,net4$edge),]
net5<-net5[match(net1$edge,net5$edge),]

df<-data.frame(
edge=net1$edge,
net1.rank=net1$rank,
net2.rank=net2$rank,
net3.rank=net3$rank,
net4.rank=net4$rank,
net5.rank=net5$rank)


## loop to calculate average edge ranks from all network inference algorithms.

x<-df[,2:(dim(df)[2])]

row.means=c()

for (i in 1:dim(df)[1])
 {
row.means[i]<-rowMeans(x[i,])
 }
 rank.agg[,2]<-row.means
 
 return(rank.agg[order(as.numeric(rank.agg[,2])),])

 }
 
 ##################################################################
 
community.rank<-AverageRank(net1=reder.rank,net2= wgcna.rank, net3=aracne.rank,net4=clr.rank,net5=mrnet.rank)


#### calculating ROC for different algorithms

 
 ### #################################################
### Gold standard network

gs.data<-gsub("re","ac",gs.data)
gs.data <-data.frame(do.call('rbind', strsplit(as.character(gs.data),'ac',fixed=F)))
gs.100<-gs.data

#### extracting true positives from gold standard network
gs.edge <- paste(gs.100[,1],gs.100[,2], sep = '_') 
gs.tp<-gsub(" ","", gs.edge , fixed=TRUE)

### Method 1 converting to binary vector as edgelist
gs.m1<-matrix(nrow=nrow(community.rank),ncol=2,byrow=TRUE)
gs.m1[,1]<-as.matrix(community.rank[,1])
gs.bin<-(gs.m1[,1] %in% gs.tp)+0


#################################################
################################################
######### Plotting ROC curve and area under ROC
require(ROCR)

pred.cons <- prediction((1-as.numeric(community.rank[,2])),gs.bin)


perf.cons <- performance(pred.cons,"tpr","fpr")
auc.cons <- unlist(slot(performance(pred.cons,"auc"), "y.values"))
maxauc.cons<-max(round(auc.cons, digits = 3))
maxauc.cons
size10.10aurocBorda<-maxauc.cons
save(size10.10auroc, file="size10.10aurocBorda.rda")

P=10
N=9
T=P+N
r=1
b=0.5
p.pos<- (1/T)+(b/P)*(1-(2*(r-1)/(T-1)))
p.pos

