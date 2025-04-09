library(conStruct)
library(elevatr)
library(raster)

#import csv file with admixture Qvalues.
QuailData<-read.csv("~/Dropbox/CCGP/CCGP_bioinformatics_analyses/Pop_genomic_analyses/Callipepla/CCGP_Callipepla_SampleData_formap.csv", header=TRUE)
QuailCoords<-data.frame(QuailData$Longitude, QuailData$Latitude)
colnames(QuailCoords)<-c("Longitude", "Latitude")
QuailCoords<-as.matrix(QuailCoords)
#create matrix for K=3 and K=4
K3<-data.frame(QuailData$Q31, QuailData$Q32, QuailData$Q33)
K4<-data.frame(QuailData$Q43, QuailData$Q42, QuailData$Q44, QuailData$Q41)
colnames(K3)<-c("Q31", "Q32", "Q33")
colnames(K4)<-c("Q43", "Q42", "Q44", "Q41")

K3mat<-as.matrix(K3)
K4mat<-as.matrix(K4)

#import polygon of california from the getData function in raster
USA<-getData("GADM",country="United States", level=1)
#west_states = "California"
Cal<-USA[USA$NAME_1 %in% "California",]

#import elevation raster and use mask function to clip raster to just California.
#elev<-get_elev_raster(Cal, z=9)
#elev_mask<-mask(elev,Cal)
#plot admixture pies on top of elevation and california map.
#col<-c("#F8766D", "#05BFC4", "#994000")
plot(elev_mask, col=grey.colors(25, start=0.95,end=0))
lines(Cal)
make.admix.pie.plot(K3mat, QuailCoords, add=TRUE, layer.colors=c("#05BFC4", "#F8766D", "#FFC024"))

