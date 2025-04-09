#--------------------------------------------------------
#########################################################
# plot GONE results
#########################################################
#--------------------------------------------------------

# plot the different GONE reps for the
par(xpd=FALSE)
library(scales)
library(matrixStats)

pdf("~/Dropbox/CCGP/CCGP_bioinformatics_analyses/Pop_genomic_analyses/Callipepla/Demography/DemographyResults_Sept24/GONE/December2024GONE_results/GONE_results_QuailSpecies_9Dec.pdf")
plot(c(0,150),c(0,100000),type="n",ylab=expression(paste("Historical ",italic(""*N*"")[e],sep="")), xlab="",cex.lab=1.5, xaxt="n")
#,xlab="Generations back in time"     
axis1<-seq(0,156,12) # years
axis2<-seq(0,150,25) # generations
axis(1, at=axis1, labels=axis1*2.5)
axis(1, at=axis2, col.axis="gray", line = 1, tick=FALSE)
    

##################################################################################
# Gambel's Quail
##################################################################################
setwd("~/Dropbox/CCGP/CCGP_bioinformatics_analyses/Pop_genomic_analyses/Callipepla/Demography/DemographyResults_Sept24/GONE/December2024GONE_results/GAM_DEC5_TEMPORARY_FILES/outfileLD_TEMP")
files <- paste("outfileLD_",1:100,"_GONE_Nebest",sep="")

NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# get CI
NeCI <- matrix(NA,nrow=200,ncol=2)
NeCI2 <- matrix(NA,nrow=200,ncol=2)

for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.025,0.975))
  NeCI2[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
  }

lines(5:200,rowMedians(NeMat[5:200,]),col="#05BFC4",lwd=4)
#polygon(x=c(1:200,rev(1:200)),y=c(NeCI[1:200,1],rev(NeCI[1:200,2])),col=adjustcolor("#05BFC4",alpha.f=0.2),border=NA)
polygon(x=c(5:200,rev(5:200)),y=c(NeCI2[5:200,1],rev(NeCI2[5:200,2])),col=adjustcolor("#05BFC4",alpha.f=0.15),border=NA)

##################################################################################
# Mountain Quail
##################################################################################
setwd("~/Dropbox/CCGP/CCGP_bioinformatics_analyses/Pop_genomic_analyses/Callipepla/Demography/DemographyResults_Sept24/GONE/December2024GONE_results/MTN_TEMPORARY_FILES_5DEC/outfileLD_TEMP")
files <- paste("outfileLD_",1:100,"_GONE_Nebest",sep="")

NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# get CI
NeCI <- matrix(NA,nrow=200,ncol=2)
NeCI2 <- matrix(NA,nrow=200,ncol=2)

for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.025,0.975))
  NeCI2[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

lines(5:200,rowMedians(NeMat[5:200,]),col="#169501",lwd=4)
#polygon(x=c(1:200,rev(1:200)),y=c(NeCI[1:200,1],rev(NeCI[1:200,2])),col=adjustcolor("#169501",alpha.f=0.2),border=NA)
polygon(x=c(5:200,rev(5:200)),y=c(NeCI2[5:200,1],rev(NeCI2[5:200,2])),col=adjustcolor("#169501",alpha.f=0.15),border=NA)

##################################################################################
# California Quail
##################################################################################
setwd("~/Dropbox/CCGP/CCGP_bioinformatics_analyses/Pop_genomic_analyses/Callipepla/Demography/DemographyResults_Sept24/GONE/December2024GONE_results/CalDec5_TEMPORARY_FILES/outfileLD_TEMP")
files <- paste("outfileLD_",1:100,"_GONE_Nebest",sep="")

NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# get CI
NeCI <- matrix(NA,nrow=200,ncol=2)
NeCI2 <- matrix(NA,nrow=200,ncol=2)

for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.025,0.975))
  NeCI2[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

lines(5:200,rowMedians(NeMat[5:200,]),col="#F7766D",lwd=4)
#polygon(x=c(1:200,rev(1:200)),y=c(NeCI[1:200,1],rev(NeCI[1:200,2])),col=adjustcolor("#F7766D",alpha.f=0.2),border=NA)
polygon(x=c(5:200,rev(5:200)),y=c(NeCI2[5:200,1],rev(NeCI2[5:200,2])),col=adjustcolor("#F7766D",alpha.f=0.15),border=NA)

dev.off()

##################################################################################
# California Quail_all
##################################################################################

pdf("~/Dropbox/CCGP/CCGP_bioinformatics_analyses/Pop_genomic_analyses/Callipepla/Demography/DemographyResults_Sept24/GONE/December2024GONE_results/GONE_results_CalQuail_only_9Dec.pdf")

plot(c(0,150),c(0,125000),type="n",ylab=expression(paste("Historical ",italic(""*N*"")[e],sep="")),xlab="",xaxt='n',cex.lab=1.5)
#,xlab="Generations back in time"

axisb<-seq(0,156,12) # years
axisc<-seq(0,150,25) # generations
axis(1, at=axisb, labels=axisb*2.5)
axis(1, at=axisc, line = 1, col.axis="gray", tick=FALSE)

lines(5:200,rowMedians(NeMat[5:200,]),col="gray",lwd=4)
#polygon(x=c(1:200,rev(1:200)),y=c(NeCI[1:200,1],rev(NeCI[1:200,2])),col=adjustcolor("gray",alpha.f=0.2),border=NA)
polygon(x=c(5:200,rev(5:200)),y=c(NeCI2[5:200,1],rev(NeCI2[5:200,2])),col=adjustcolor("gray",alpha.f=0.15),border=NA)

##################################################################################
# California Quail_NW
##################################################################################
setwd("~/Dropbox/CCGP/CCGP_bioinformatics_analyses/Pop_genomic_analyses/Callipepla/Demography/DemographyResults_Sept24/GONE/December2024GONE_results/NWCal_DEC5_TEMPORARY_FILES/outfileLD_TEMP")
files <- paste("outfileLD_",1:100,"_GONE_Nebest",sep="")

NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# get CI
NeCI <- matrix(NA,nrow=200,ncol=2)
NeCI2 <- matrix(NA,nrow=200,ncol=2)

for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.025,0.975))
  NeCI2[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

lines(5:200,rowMedians(NeMat[5:200,]),col="#FFC024",lwd=4)
#polygon(x=c(1:200,rev(1:200)),y=c(NeCI[1:200,1],rev(NeCI[1:200,2])),col=adjustcolor("#FFC024",alpha.f=0.2),border=NA)
polygon(x=c(5:200,rev(5:200)),y=c(NeCI2[5:200,1],rev(NeCI2[5:200,2])),col=adjustcolor("#FFC024",alpha.f=0.15),border=NA)

##################################################################################
# California Quai_SW
##################################################################################
setwd("~/Dropbox/CCGP/CCGP_bioinformatics_analyses/Pop_genomic_analyses/Callipepla/Demography/DemographyResults_Sept24/GONE/December2024GONE_results/SWQuail_5DEC_TEMPORARY_FILES/outfileLD_TEMP")
files <- paste("outfileLD_",1:100,"_GONE_Nebest",sep="")

NeMat <- NULL
for(i in 1:100){
  dat <- read.table(files[i],skip=2)
  NeMat <- cbind(NeMat,dat[,2])
}

# get CI
NeCI <- matrix(NA,nrow=200,ncol=2)
NeCI2 <- matrix(NA,nrow=200,ncol=2)

for(i in 1:200){
  NeCI[i,] <- quantile(NeMat[i,],probs=c(0.025,0.975))
  NeCI2[i,] <- quantile(NeMat[i,],probs=c(0.05,0.95))
}

lines(5:200,rowMedians(NeMat[5:200,]),col="#F7766D",lwd=4)
#polygon(x=c(1:200,rev(1:200)),y=c(NeCI[1:200,1],rev(NeCI[1:200,2])),col=adjustcolor("#F7766D",alpha.f=0.2),border=NA)
polygon(x=c(5:200,rev(5:200)),y=c(NeCI2[5:200,1],rev(NeCI2[5:200,2])),col=adjustcolor("#F7766D",alpha.f=0.15),border=NA)
dev.off()