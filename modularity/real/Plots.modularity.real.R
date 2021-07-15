
load("~/Dropbox/Phd/R/Modularity/Real/RedeR/R.cluster.valid.real.rda")
load("~/Dropbox/Phd/R/Modularity/Real/WGCNA/W.cluster.valid.real.rda")
load("~/Dropbox/Phd/R/Modularity/Real/SIMoNe/S.cluster.valid.real.rda")

k=c(4,8,12,16)

## cluster validation indices 
# Silhouette width scores


png(file = "~/Dropbox/Phd/R/Modularity/Real/cluster.validations.real.png", width = 1000, height = 1000, res=120);
par(mar=c(7.1, 5.1, 4.1, 2.1))
par(mfrow=c(3,2))
# mean sil width
plot(R.real.mean.sil,ylim=c(min(R.real.mean.sil,W.real.mean.sil,S.real.mean.sil),
     max(R.real.mean.sil,W.real.mean.sil,S.real.mean.sil)),type="o",pch=16,col="red",xaxt = "n",
     xlab='Number of Modules', ylab="mean Silhouette width")
axis(1, at=1:length(k), labels=c(k))
points(W.real.mean.sil,ylim=c(min(W.real.mean.sil),max(W.real.mean.sil)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
     ylab="")
points(S.real.mean.sil,ylim=c(min(S.real.mean.sil),max(S.real.mean.sil)),type="o",pch=16,col="black",xaxt = "n", xlab='Number of Modules',
     ylab="")
legend("topright", legend=c("RedeR","WGCNA","SIMoNe"), fill=c("red","blue","black"), bty = "n",cex=1)

# avg sil width
plot(R.real.avg.sil1,type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules',
     ylab="Avg silhouette width",ylim=c(min(R.real.avg.sil1,W.real.avg.sil1,S.real.avg.sil1),
     max(R.real.avg.sil1,W.real.avg.sil1,S.real.avg.sil1)))
axis(1, at=1:length(k), labels=c(k))
points(W.real.avg.sil1,ylim=c(min(W.real.avg.sil1),max(W.real.avg.sil1)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
     ylab="Avg silhouette width")
axis(1, at=1:length(k), labels=c(k))
points(S.real.avg.sil1,ylim=c(min(S.real.avg.sil1),max(S.real.avg.sil1)),type="o",pch=16,col="black",xaxt = "n", xlab='Number of Modules',
     ylab="Avg silhouette width")
axis(1, at=1:length(k), labels=c(k))

# dunn index
plot(R.real.dunn,ylim=c(min(R.real.dunn,W.real.dunn,S.real.dunn),max(R.real.dunn,W.real.dunn,S.real.dunn)),
     type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules', ylab="Dunn index")
axis(1, at=1:length(k), labels=c(k))
points(W.real.dunn,ylim=c(min(W.real.dunn),max(W.real.dunn)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
     ylab="Dunn index")
points(S.real.dunn,ylim=c(min(S.real.dunn),max(S.real.dunn)),type="o",pch=16,col="black",xaxt = "n", xlab='Number of Modules',
     ylab="Dunn index")

# dunn2 index
plot(R.real.dunn2,ylim=c(min(R.real.dunn2,W.real.dunn2,S.real.dunn2),max(R.real.dunn2,W.real.dunn2,S.real.dunn2)),
     type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules', ylab="Dunn2 index")
axis(1, at=1:length(k), labels=c(k))
points(W.real.dunn2,ylim=c(min(W.real.dunn2),max(W.real.dunn2)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
       ylab="dunn2 index")
points(S.real.dunn2,ylim=c(min(S.real.dunn2),max(S.real.dunn2)),type="o",pch=16,col="black",xaxt = "n", xlab='Number of Modules',
       ylab="dunn2 index")

# sep index
plot(R.real.sep.idx,ylim=c(min(R.real.sep.idx,W.real.sep.idx,S.real.sep.idx),max(R.real.sep.idx,W.real.sep.idx,S.real.sep.idx)),
     type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules', ylab="sep index")
axis(1, at=1:length(k), labels=c(k))
points(W.real.sep.idx,ylim=c(min(W.real.sep.idx),max(W.real.sep.idx)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
       ylab="sep.idx index")
points(S.real.sep.idx,ylim=c(min(S.real.sep.idx),max(S.real.sep.idx)),type="o",pch=16,col="black",xaxt = "n", xlab='Number of Modules',
       ylab="sep.idx index")

# biological model score
################################
## Model scores

load("~/Dropbox/Phd/R/Modularity/Real/RedeR/BP/R.real.dat.BP.rda")
load("~/Dropbox/Phd/R/Modularity/Real/WGCNA/BP/W.real.dat.BP.rda")
load("~/Dropbox/Phd/R/Modularity/Real/SIMoNe/BP/S.real.dat.BP.rda")

plot(R.real.dat.BP,ylim=c(min(R.real.dat.BP,W.real.dat.BP,S.real.dat.BP)-0.2,max(R.real.dat.BP,W.real.dat.BP,S.real.dat.BP)+0.2),
     type="o",xlab="Number of Modules",xaxt = "n",
     ylab="Model score", main="Size 500",pch=16,col="red")
axis(1, at=1:length(k), labels=c(k))

points(W.real.dat.BP,type="o",xaxt = "n",pch=16,col="blue")

points(S.real.dat.BP,type="o",xaxt = "n",pch=16,col="black")

dev.off()


################################
## Model scores

load("~/Dropbox/Phd/R/Modularity/Real/RedeR/BP/R.real.dat.BP.rda")
load("~/Dropbox/Phd/R/Modularity/Real/WGCNA/BP/W.real.dat.BP.rda")
load("~/Dropbox/Phd/R/Modularity/Real/SIMoNe/BP/S.real.dat.BP.rda")
R.real.dat.BP
W.real.dat.BP
S.real.dat.BP
k=c(4,8,12,16)

png(file = "~/Dropbox/Phd/R/Modularity/Real/Model.score.n500.png", 
       width = 700, height = 700, res=180);
plot(R.real.dat.BP,ylim=c(min(R.real.dat.BP,W.real.dat.BP,S.real.dat.BP)-0.2,max(R.real.dat.BP,W.real.dat.BP,S.real.dat.BP)+0.2),
     type="o",xlab="Number of Modules",xaxt = "n",
     ylab="Model score", main="S.cerevisiae: Biological process",pch=16,col="red")
axis(1, at=1:length(k), labels=c(k))

points(W.real.dat.BP,type="o",xaxt = "n",pch=16,col="blue")

points(S.real.dat.BP,type="o",xaxt = "n",pch=16,col="black")

legend("topleft", legend=c("RedeR","WGCNA","SIMoNe"), fill=c("red","blue","black"), bty = "n",cex=0.8)

dev.off()

################################
## GO annotated scores

load("~/Dropbox/Phd/R/Modularity/Real/RedeR/BP/Reder.real.GO.annot.BP.rda")
load("~/Dropbox/Phd/R/Modularity/Real/WGCNA/BP/WGCNA.real.GO.annot.BP.rda")
load("~/Dropbox/Phd/R/Modularity/Real/SIMoNe/BP/SIMoNe.real.GO.annot.BP.rda")

png(file = "~/Dropbox/Phd/R/Modularity/Real/GO.annoted.n500.png", 
    width = 1000, height = 1000, res=180);
par(mfrow=c(2,2))
##   4 clusters
plot(R.real.BP.logp.4,col="red",pch=6, type="o",ylab="No of enriched GO terms",
     xlab=expression(paste(-log[10],"(p-value)")),main="4 Modules",
     ylim=c(0,max(R.real.BP.logp.4,W.real.BP.logp.4,S.real.BP.logp.4)),
     xlim=c(1,max(R.real.BP.logp.4[,1],W.real.BP.logp.4[,1],S.real.BP.logp.4[,1])))
points(W.real.BP.logp.4,col="blue",pch=23, type="o")
points(S.real.BP.logp.4,col="black",pch=12, type="o")

legend("topright",c("RedeR","WGCNA","SIMoNe"),pch=c(6,23,12),col=c("red","blue","black"),cex=0.8)
abline(v=-log10(0.05),lwd=1.5, lty="dashed" ,col="grey50")

# 8 clusters
plot(R.real.BP.logp.8,col="red",pch=6, type="o",ylab="No of enriched GO terms",
     xlab=expression(paste(-log[10],"(p-value)")),main="8 Modules",
     ylim=c(0,max(R.real.BP.logp.8,W.real.BP.logp.8,S.real.BP.logp.8)),
     xlim=c(1,max(R.real.BP.logp.8[,1],W.real.BP.logp.8[,1],S.real.BP.logp.8[,1])))
points(W.real.BP.logp.8,col="blue",pch=23, type="o")
points(S.real.BP.logp.8,col="black",pch=12, type="o")
abline(v=-log10(0.05),lwd=1.5, lty="dashed" ,col="grey50")

# 12 clusters
plot(R.real.BP.logp.12,col="red",pch=6, type="o",ylab="No of enriched GO terms",
     xlab=expression(paste(-log[10],"(p-value)")),main="12 Modules",
     ylim=c(0,max(R.real.BP.logp.12,W.real.BP.logp.12,S.real.BP.logp.12)),
     xlim=c(1,max(R.real.BP.logp.12[,1],W.real.BP.logp.12[,1],S.real.BP.logp.12[,1])))

points(W.real.BP.logp.12,col="blue",pch=23, type="o")
points(S.real.BP.logp.12,col="black",pch=12, type="o")

abline(v=-log10(0.05),lwd=1.5, lty="dashed" ,col="grey50")

# 16 clusters
plot(R.real.BP.logp.16,col="red",pch=6, type="o",ylab="No of enriched GO terms",
     xlab=expression(paste(-log[10],"(p-value)")),main="16 Modules",
     ylim=c(0,max(R.real.BP.logp.16,W.real.BP.logp.16,S.real.BP.logp.16)),
     xlim=c(1,max(R.real.BP.logp.16[,1],W.real.BP.logp.16[,1],S.real.BP.logp.16[,1])))

points(W.real.BP.logp.16,col="blue",pch=23, type="o")
points(S.real.BP.logp.16,col="black",pch=12, type="o")

abline(v=-log10(0.05),lwd=1.5, lty="dashed" ,col="grey50")

dev.off()

#############################################
## GO percentage annotated scores

load("~/Dropbox/Phd/R/Modularity/Real/RedeR/BP/Reder.real.Perc..BP.rda")
load("~/Dropbox/Phd/R/Modularity/Real/WGCNA/BP/WGCNA.real.Perc..BP.rda")
load("~/Dropbox/Phd/R/Modularity/Real/SIMoNe/BP/SIMoNe.real.Perc..BP.rda")

png(file = "~/Dropbox/Phd/R/Modularity/Real/GO.Percentage.n500.png", 
    width = 1000, height = 1000, res=180);
par(mfrow=c(2,2))
##   4 clusters
plot(perc.R.real.BP.4,col="red",pch=6, type="o",ylab="% GO annotated Modules",
     xlab=expression(paste(-log[10],"(p-value)")),main="4 Modules",
     ylim=c(0,max(perc.R.real.BP.4,perc.W.real.BP.4,perc.S.real.BP.4)),
     xlim=c(1,max(perc.R.real.BP.4[,1],perc.W.real.BP.4[,1],perc.S.real.BP.4[,1])))
points(perc.W.real.BP.4,col="blue",pch=23, type="o")
points(perc.S.real.BP.4,col="black",pch=12, type="o")

abline(v=-log10(0.05),lwd=1.5, lty="dashed" ,col="grey50")

# 8 clusters
plot(perc.R.real.BP.8,col="red",pch=6, type="o",ylab="% GO annotated Modules",
     xlab=expression(paste(-log[10],"(p-value)")),main="8 Modules",
     ylim=c(0,max(perc.R.real.BP.8,perc.W.real.BP.8,perc.S.real.BP.8)),
     xlim=c(1,max(perc.R.real.BP.8[,1],perc.W.real.BP.8[,1],perc.S.real.BP.8[,1])))
points(perc.W.real.BP.8,col="blue",pch=23, type="o")
points(perc.S.real.BP.8,col="black",pch=12, type="o")
abline(v=-log10(0.05),lwd=1.5, lty="dashed" ,col="grey50")
legend("topright",c("RedeR","WGCNA","SIMoNe"),pch=c(6,23,12),col=c("red","blue","black"),cex=0.8)

# 12 clusters
plot(perc.R.real.BP.12,col="red",pch=6, type="o",ylab="% GO annotated Modules",
     xlab=expression(paste(-log[10],"(p-value)")),main="12 Modules",
     ylim=c(0,max(perc.R.real.BP.12,perc.W.real.BP.12,perc.S.real.BP.12)),
     xlim=c(1,max(perc.R.real.BP.12[,1],perc.W.real.BP.12[,1],perc.S.real.BP.12[,1])))

points(perc.W.real.BP.12,col="blue",pch=23, type="o")
points(perc.S.real.BP.12,col="black",pch=12, type="o")

abline(v=-log10(0.05),lwd=1.5, lty="dashed" ,col="grey50")

# 16 clusters
plot(perc.R.real.BP.16,col="red",pch=6, type="o",ylab="% GO annotated Modules",
     xlab=expression(paste(-log[10],"(p-value)")),main="16 Modules",
     ylim=c(0,max(perc.R.real.BP.16,perc.W.real.BP.16,perc.S.real.BP.16)),
     xlim=c(1,max(perc.R.real.BP.16[,1],perc.W.real.BP.16[,1],perc.S.real.BP.16[,1])))

points(perc.W.real.BP.16,col="blue",pch=23, type="o")
points(perc.S.real.BP.16,col="black",pch=12, type="o")

abline(v=-log10(0.05),lwd=1.5, lty="dashed" ,col="grey50")

dev.off()

