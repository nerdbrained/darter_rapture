### making ultrametric tree using calibration dates

library(ape)
library(phytools)

#calibrating with 3 points
#spectabile-caeruleum split (node 107)= 22.5; flabellare-cragini split (node 58)= 21.5; etheostoma root (node 57)=29
#gleaned from Kelly et al. 2015 tree figure, added 5 million year buffer on each side
### tree estimated in IQtree from concatenated data

darter.tree<-read.tree("~/Desktop/rapture/manuscript/secapr/iq_out/boot/populations.var.phylip.treefile")

darter.root2<-reroot(darter.tree,node.number=103,position=0.05)

### use interactive calibration to add all 3 points 
darter.calibration <- makeChronosCalib(darter.root2, interactive=TRUE)

### make time tree using 3 different models and check likelihood - correlated is highest

darter.timetree3.correlated <- chronos(darter.root2, lambda = 1, model = "correlated", calibration = darter.calibration, control = chronos.control() )

darter.timetree3.relaxed <- chronos(darter.root2, lambda = 1, model = "relaxed", calibration = darter.calibration, control = chronos.control() )

darter.timetree3.discrete <- chronos(darter.root2, lambda = 1, model = "discrete", calibration = darter.calibration, control = chronos.control() )

### get support values, color-code

supn=as.integer(darter.timetree3.correlated$node.label)
supn[supn>99]<-"green"
supn[supn>95]<-"yellow"
supn[supn>80]<-"blue"
supn[supn>0]<-"black"
supc = na.replace()

supn=as.integer(darter.timetree3.correlated$node.label)
supn[is.na(supn)]<-"black"
supn[7]<-"white"
supn[as.integer(supn)>99]<-"green"
supn[as.integer(supn)>95]<-"blue"
supn[as.integer(supn)>80]<-"purple"
supn[as.integer(supn)>50]<-"yellow"
supn[as.integer(supn)>0]<-"white"

### plot tree with support values

plot(darter.timetree3.correlated,cex=0.5)
nodelabels(text=NULL,pch=21,bg=supn,cex=0.8)
axisPhylo()
legend(x=0,y=20,legend=c("100%","95-99%","80-95%","50-80%","<50%"),box.lwd=0,pch=c(19,19,19,19,1),col=c("green","blue","purple","yellow","black"))
text(x=3,y=22,labels="Bootstrap Support")

### write time-calibrated tree for use in phydesign

write.tree(darter.timetree3.correlated,"~/Desktop/rapture/manuscript/secapr/iq_out/darter_timetree3_correlated.newick")

### phylogeography: plotting IQtree result on map
### need treefile plus some files connecting sites, individuals, and metapopulations (sitelist.csv, indsites_meta2.csv, metalist2.csv)
### for plotting river networks, also need the "USA_Rivers_and_Streams" shapefile (google it!)

library(ape)
library(phytools)
library(plotrix)
library(stringr)
library(rgdal)
library(raster)
library(rgeos)

tree<-read.tree(file="~/Desktop/rapture/manuscript/secapr/iq_out/darter_timetree3_correlated.newick")

### create and modify list of IDs - remove metapopulation/species IDs

ID=tree$tip.label
ID=str_remove(ID,"Ecr_")
ID=str_remove(ID,"Esp_")
ID=str_remove(ID,"Eca_")
ID=str_remove(ID,"Esp_")
ID=str_remove(ID,"SFN_")
ID=str_remove(ID,"NFN_")
ID=str_remove(ID,"Col_")
ID=str_remove(ID,"Chik_")
ID=str_remove(ID,"Cim_")
ID=str_remove(ID,"Salt_")
ID=str_remove(ID,"Med_")
ID=str_remove(ID,"LArk_")
ID=str_remove(ID,"Wal_")
ID=str_remove(ID,"Spring_")
ID=str_remove(ID,"Ark_")
ID=str_remove(ID,"Ratt_")

### convert taxon IDs to site IDs

sites=gsub('.{2}$', '', ID)

### list of taxon IDs from ultrametric tree

treelen$tip.label=ID

### extract location and metapopulation IDs for E. cragini sites

sitelocs=read.csv("~/Desktop/rapture/ECRapture45/sitelocs.csv")
indsites=read.csv("~/Desktop/rapture/ECRapture45/indsites_meta2.csv",header=T)
indsites$long<-sitelocs$x[match(indsites$site,sitelocs$SITE.NAME)]
indsites$lat<-sitelocs$y[match(indsites$site,sitelocs$SITE.NAME)]
metalist=read.csv("~/Desktop/rapture/ECRapture45/metalist2.csv",header=F)
metalist$V2=rainbow(n=14)
indsites$col<-metalist$V2[match(indsites$meta,metalist$V1)]
indsites$ind=NULL
indsites_uni=unique(indsites)
poplocs=data.frame(sites,row.names=ID)
poplocs$lat=indsites_uni$lat[match(poplocs$sites,indsites_uni$site)]
poplocs$long=indsites_uni$long[match(poplocs$sites,indsites_uni$site)]
poplocs$col=indsites_uni$col[match(poplocs$sites,indsites_uni$site)]

### added in sites for other Etheostoma species

poplocs$long[49:56]=c(-88.762853,-88.762853,-88.031704,-88.031704,-85.3584,-85.35245,-85.35245,-85.3584)
poplocs$lat[49:56]=c(41.440032,41.440032,40.152534,40.152534,42.3241,42.41763,42.41763,42.3241)
poplocs$col[49:56]=c("black","black","black","black","black","black","black","black")

### have to use setNames to get the colors to display correctly

cols=setNames(poplocs$col,treelen$tip.label)

### plot tree on map!

obj<-phylo.to.map(tree=treelen,coords=poplocs[,2:3],type="phylogram",database="state",xlim=c(min(poplocs$long),max(poplocs$long)),ylim=c(min(poplocs$lat),max(poplocs$lat)),colors=cols)
plot(obj,xlim=c(min(poplocs$long),max(poplocs$long)),ylim=c(min(poplocs$lat),max(poplocs$lat)),colors=cols)

###add streams
USstreams=readOGR(dsn="~/Desktop/kansasGIS/USA_Rivers_and_Streams", layer="USA_Rivers_and_Streams")
x=c(-104.83,-94.1)
y=c(36.08,39.52)
xy <- cbind(x,y)
S <- SpatialPoints(xy)
darterbox=bbox(S)

gClip <- function(shp, bb){
  if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else b_poly <- as(extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}


ArkRiv=USstreams[which(USstreams$NAME == "Arkansas River"),]
CimRiv=USstreams[which(USstreams$NAME == "Cimarron River"),]
ChikRiv=USstreams[which(USstreams$NAME == "Chikaskia River"),]
MedRiv=USstreams[which(USstreams$NAME == "Medicine Lodge River"),]
SaltRiv=USstreams[which(USstreams$NAME == "Salt Fork Arkansas River"),]
NFNRiv=USstreams[which(USstreams$NAME == "North Fork Ninnescah River"),]
SFNRiv=USstreams[which(USstreams$NAME == "South Fork Ninnescah River"),]
WalnutCrk=USstreams[which(USstreams$NAME == "Walnut Creek"),]
RattCrk=USstreams[which(USstreams$NAME == "Rattlesnake Creek"),]
IllRiv=USstreams[which(USstreams$NAME == "Illinois River"),]
SpringRiv=USstreams[which(USstreams$NAME == "Spring River"),]
BigSandy=USstreams[which(USstreams$NAME == "Big Sandy Creek"),]
Rush=USstreams[which(USstreams$NAME == "Rush Creek"),]


ArkRiv_c=gClip(ArkRiv,darterbox)
CimRiv_c=gClip(CimRiv,darterbox)
ChikRiv_c=gClip(ChikRiv,darterbox)
MedRiv_c=gClip(MedRiv,darterbox)
SaltRiv_c=gClip(SaltRiv,darterbox)
NFNRiv_c=gClip(NFNRiv,darterbox)
SFNRiv_c=gClip(SFNRiv,darterbox)
WalnutCrk_c=gClip(WalnutCrk,darterbox)
RattCrk_c=gClip(RattCrk,darterbox)
IllRiv_c=gClip(IllRiv,darterbox)
SpringRiv_c=gClip(SpringRiv,darterbox)
BigSandy_c=gClip(BigSandy,darterbox)
Rush_c=gClip(Rush,darterbox)

plot(ArkRiv_c,lwd=0.4,col="black",add=T)
plot(CimRiv_c,lwd=0.4,col="black",add=T)
plot(ChikRiv_c,lwd=0.4,col="black",add=T)
plot(MedRiv_c,lwd=0.4,col="black",add=T)
plot(SaltRiv_c,lwd=0.4,col="black",add=T)
plot(NFNRiv_c,lwd=0.4,col="black",add=T)
plot(SFNRiv_c,lwd=0.4,col="black",add=T)
plot(WalnutCrk_c,lwd=0.4,col="black",add=T)
plot(RattCrk_c,lwd=0.4,col="black",add=T)
plot(IllRiv_c,lwd=0.4,col="black",add=T)
plot(SpringRiv_c,lwd=0.4,col="black",add=T)
plot(BigSandy_c,lwd=0.4,col="black",add=T)
plot(Rush_c,lwd=0.4,col="black",add=T)

### re-plot sites on map

points(x=poplocs$long,y=poplocs$lat,col=cols,pch=19,cex=1.2)