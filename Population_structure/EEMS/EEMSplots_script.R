library("rworldmap")
library("rworldxtra")
library("rEEMSplots")
library("rgdal")
library("raster")

#specify path to eems data results
eems_results_chain1<-"~/Dropbox/CCGP/CCGP_bioinformatics_analyses/Pop_genomic_analyses/Callipepla/CalQuail_eems_output/CalQuail_chain1"
eems_results_chain2<-"~/Dropbox/CCGP/CCGP_bioinformatics_analyses/Pop_genomic_analyses/Callipepla/CalQuail_eems_output/CalQuail_chain2"
eems_results_chain3<-"~/Dropbox/CCGP/CCGP_bioinformatics_analyses/Pop_genomic_analyses/Callipepla/CalQuail_eems_output/CalQuail_chain3"
name_figures = "~/Dropbox/CCGP/CCGP_bioinformatics_analyses/Pop_genomic_analyses/Callipepla/CalQuail_eems_output/"
#create vector of the results from eems MCMC chains started from different starting parameters.
eems_data<-c(eems_results_chain1,eems_results_chain2,eems_results_chain3)

projection_none <- "+proj=longlat +datum=WGS84" 
projection_mercator <- "+proj=merc +datum=WGS84" 
USA<-getData("GADM",country="United States", level=1)
west_states = c("Arizona","California","Idaho", "Nevada","Oregon","Utah")
West<-USA[USA$NAME_1 %in% west_states,]
West <- spTransform(West, CRS(projection_mercator))

#run eems-plots to generate output figures
eems.plots(mcmcpath = eems_data,
plotpath = paste0(name_figures, "_EEMS_out"), longlat = TRUE, add.grid=FALSE,add.demes=TRUE,min.cex.demes=1, max.cex.demes=3, projection.in = projection_none, projection.out = projection_mercator, add.map=TRUE, col.map="black", lwd.map=1,m.plot.xy = {
    plot(West, col = NA, add = TRUE,lwd=1)
  },
  q.plot.xy = {
    plot(West, col = NA, add = TRUE,lwd=1)
  }
)