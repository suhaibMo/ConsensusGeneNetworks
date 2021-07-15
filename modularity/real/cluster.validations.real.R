
load("/Users/suhaibmohammed/Dropbox/Phd/R/Modularity/Real/RedeR/R.cluster.valid.real.rda")
load("/Users/suhaibmohammed/Dropbox/Phd/R/Modularity/Real/WGCNA/W.cluster.valid.real.rda")
load("/Users/suhaibmohammed/Dropbox/Phd/R/Modularity/Real/SIMoNe/S.cluster.valid.real.rda")


## cluster validation indices 
# Silhouette width scores

k=c(4,8,12,16)

png(file = "~/Dropbox/Phd/R/Modularity/Real/cluster.validations.real.png",
    width = 1200, height = 500, res=180);
=par(mar=c(7.1, 5.1, 4.1, 2.1))
par(mfrow=c(1,3))
# # mean sil width
# plot(R.real.mean.sil,ylim=c(min(R.real.mean.sil,W.real.mean.sil,S.real.mean.sil),
#      max(R.real.mean.sil,W.real.mean.sil,S.real.mean.sil)),type="o",pch=16,col="red",xaxt = "n",
#      xlab='Number of Modules', ylab="Avg Silhouette width")
# axis(1, at=1:length(k), labels=c(k))
# points(W.real.mean.sil,ylim=c(min(W.real.mean.sil),max(W.real.mean.sil)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
#      ylab="")
# points(S.real.mean.sil,ylim=c(min(S.real.mean.sil),max(S.real.mean.sil)),type="o",pch=16,col="black",xaxt = "n", xlab='Number of Modules',
#      ylab="")

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
legend("bottomleft", legend=c("RedeR","WGCNA","SIMoNe"), fill=c("red","blue","black"), bty = "n",cex=1)

# dunn index
plot(R.real.dunn,ylim=c(min(R.real.dunn,W.real.dunn,S.real.dunn),max(R.real.dunn,W.real.dunn,S.real.dunn)),
     type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules', ylab="Dunn index")
axis(1, at=1:length(k), labels=c(k))
points(W.real.dunn,ylim=c(min(W.real.dunn),max(W.real.dunn)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
       ylab="Dunn index")
points(S.real.dunn,ylim=c(min(S.real.dunn),max(S.real.dunn)),type="o",pch=16,col="black",xaxt = "n", xlab='Number of Modules',
       ylab="Dunn index")

# # dunn2 index
# plot(R.real.dunn2,ylim=c(min(R.real.dunn2,W.real.dunn2,S.real.dunn2),max(R.real.dunn2,W.real.dunn2,S.real.dunn2)),
#      type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules', ylab="Dunn like index")
# axis(1, at=1:length(k), labels=c(k))
# points(W.real.dunn2,ylim=c(min(W.real.dunn2),max(W.real.dunn2)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
#        ylab="Dunn like index")
# points(S.real.dunn2,ylim=c(min(S.real.dunn2),max(S.real.dunn2)),type="o",pch=16,col="black",xaxt = "n", xlab='Number of Modules',
#        ylab="Dunn like index")

# sep index
plot(R.real.sep.idx,ylim=c(min(R.real.sep.idx,W.real.sep.idx,S.real.sep.idx),max(R.real.sep.idx,W.real.sep.idx,S.real.sep.idx)),
     type="o",pch=16,col="red",xaxt = "n", xlab='Number of Modules', ylab="Sep index")
axis(1, at=1:length(k), labels=c(k))
points(W.real.sep.idx,ylim=c(min(W.real.sep.idx),max(W.real.sep.idx)),type="o",pch=16,col="blue",xaxt = "n", xlab='Number of Modules',
       ylab="sep.idx index")
points(S.real.sep.idx,ylim=c(min(S.real.sep.idx),max(S.real.sep.idx)),type="o",pch=16,col="black",xaxt = "n", xlab='Number of Modules',
       ylab="sep.idx index")


dev.off()