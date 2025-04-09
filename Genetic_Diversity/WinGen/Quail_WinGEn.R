install.packages(c('wingen', 'terra', 'raster', "viridis","SpatialKDE","sf"), repos="https://cloud.r-project.org", dependencies=TRUE)

library(wingen)
library(terra)
library(raster)
library(viridis)
library(SpatialKDE)
library(sf) 

setwd("/media/data2/Quail_reanalysis/WinGen")
QuailVCF<-"Quail_Autosome_NoHybrid_11Sept.recode.vcf"
QuailData<-read.csv("CCGP_Callipepla_SampleData.csv", header=TRUE)

Cal_lyr<- rast(ncols = 100, nrows = 100,nlyrs=1, xmin = -125, xmax = -113, ymin = 32, ymax = 43,  crs = "+proj=longlat +datum=WGS84",vals=c(-1,1)) 
QuailCoords<-data.frame(QuailData$Longitude, QuailData$Latitude)
colnames(QuailCoords)<-c("Longitude", "Latitude")
QuailCoords_longlat<-st_as_sf(QuailCoords,coords=c("Longitude", "Latitude"), crs="+proj=longlat")

projection_mercator <- "+proj=merc +datum=WGS84" 
QuailCoords_eq<-st_transform(QuailCoords_longlat, crs= "+proj=merc +datum=WGS84" )
CalLyr_eq<-project(Cal_lyr, "+proj=merc +datum=WGS84" )

#import polygon of California from raster package
USA<-getData("GADM",country="United States", level=1)
west_states = "California"
Cal<-USA[USA$NAME_1 %in% west_states,]
Cal_map <- spTransform(Cal, CRS(projection_mercator))

wgd_ho <- window_gd(QuailVCF,
QuailCoords_eq,
CalLyr_eq,
stat="Ho",
fact=3,
wdim=3,
rarify=TRUE,
rarify_n=3,
rarify_nit=5)

wgd_pi <- window_gd(QuailVCF,
QuailCoords_eq,
CalLyr_eq,
stat="pi",
fact=3,
wdim=3,
rarify=TRUE,
rarify_n=3,
rarify_nit=5)

#add colored points to map for different quail pops.
quailSp<-as.factor(QuailData$Species)
col<-c("#F8766D", "#05BFC4", "#994000")
#points(QuailCoords_eq, pch=21, cex=1.5, col="black",  bg=col[quailSp])

#add map outline of California
#ines(Cal_map)

#dev.off()


pdf("Quail_wingen_pi_17Dec2024_wdim3_fact3.pdf")
#spatial KDE plot
kgd_pi<-krig_gd(wgd_pi,index=1,grd= CalLyr_eq, disagg_grd = 2)
grid_pi <- raster(kgd_pi)
kde_lyr_pi<-kde(QuailCoords_eq, kernel="quartic", band_width=200000,decay=1, grid= grid_pi)
mgd_pi <- mask_gd(kgd_pi, kde_lyr_pi, minval = 1)
mgd_pi2<-mask_gd(mgd_pi, Cal_map)
plot_gd(mgd_pi2, col=viridis::viridis(100))

#add colored points to map for different quail pops.
points(QuailCoords_eq, pch=21, cex=1.5, col="black",  bg=col[quailSp])

#add map outline of California
lines(Cal_map)

dev.off()

pdf("Quail_wingen_ho_17Dec2024_wdim3_fact3.pdf")
#spatial KDE plot
kgd_pi<-krig_gd(wgd_ho,index=1,grd= CalLyr_eq, disagg_grd = 2)
grid_pi <- raster(kgd_pi)
kde_lyr_pi<-kde(QuailCoords_eq, kernel="quartic", band_width=200000,decay=1, grid= grid_pi)
mgd_pi <- mask_gd(kgd_pi, kde_lyr_pi, minval = 1)
mgd_pi2<-mask_gd(mgd_pi, Cal_map)
plot_gd(mgd_pi2, col=viridis::viridis(100))

#add colored points to map for different quail pops.
points(QuailCoords_eq, pch=21, cex=1.5, col="black",  bg=col[quailSp])

#add map outline of California
lines(Cal_map)

dev.off()
