### R code for plotting admixture results
### Read in metadata and create reference tables first

sitelocs=read.csv("~/Desktop/rapture/ECRapture45/sitelocs.csv")
indsites=read.csv("~/Desktop/rapture/ECRapture45/indsites_meta2.csv",header=T)
indsites$long<-sitelocs$x[match(indsites$site,sitelocs$SITE.NAME)]
indsites$lat<-sitelocs$y[match(indsites$site,sitelocs$SITE.NAME)]
#plot(x=NULL, y=NULL, xlim=c(-105,-94),ylim=c(36,40))
#points(indsites$long,indsites$lat)
#text(indsites$long,indsites$lat,labels=indsites$site,cex=0.3)

metalist=read.csv("~/Desktop/rapture/ECRapture45/metalist2.csv",header=F)
metalist$V2=rainbow(n=14)

indsites$col<-metalist$V2[match(indsites$meta,metalist$V1)]


###########Example with Rapture data aligned to E. cragini genome############


library(dplyr)

### read Q matrix output from PCAngsd
q=read.csv("~/Desktop/rapture/manuscript/revision/angsdout_rev/EcrAdalign_indf_admix.csv",header=F)

### Some data cleanup - removing one sample of uncertain provenance, some hatchery fish, and renaming one metapopulation 
q=q[-1676,]
q$metapop=indsites$meta
q$metapop=gsub('Big Sandy / Fountain Creek','Big Sandy / Rush Creeks',q$metapop)
q$ind=indsites$ind
q$site=indsites$site
q$long=indsites$long
q=q[-which(q$metapop=="Hugo Ponds Hatchery Population"),]
cols=rainbow(ncol(q))
Qsortlong=q[order(q$metapop,q$long),]
sitecounts<-as.data.frame(dplyr::count(Qsortlong,metapop))
sitecounts$delinpos<-cumsum(sitecounts$n)
sitecounts$labelpos<-(sitecounts$delinpos-(sitecounts$n/2))
Qsortlong$metapop=NULL
Qsortlong$ind=NULL
Qsortlong$site=NULL
Qsortlong$long=NULL

### make barplot
par(mai=c(1.5,0.8,0.2,0.5),xpd=T)
barplot(t(Qsortlong),col=cols,space=0,border=NA,xlab="",xaxt="n",ylab="Admixture proportions")
axis(side=1, at=sitecounts$labelpos, labels = FALSE,lwd=0,lwd.ticks=1)
axis(side=1, at=c(0,1775), labels = FALSE,lwd=1,lwd.ticks=0)
abline(v=sitecounts$delinpos,col="white",lwd=1.5)
text(x=sitecounts$labelpos,y=-0.05,labels=sitecounts$metapop,cex=0.8,adj=0,col="black",srt=315)


### Aggregate admixture proportions by site. Change # of values summarized (P1-PX) according to number of columns in q file (= # of populations)
aggQ <-ddply(q,~site,summarise,p1=mean(V1),p2=mean(V2),p3=mean(V3),p4=mean(V4),p5=mean(V5),p6=mean(V6),p7=mean(V7),p8=mean(V8),p9=mean(V9),P10=mean(V10),P11=mean(V11),P12=mean(V12),P13=mean(V13),P14=mean(V14),P15=mean(V15),P16=mean(V16))

meltQ <- melt(aggQ)

meltQ$x=sitelocs$x[match(meltQ$site, sitelocs$SITE.NAME)]
meltQ$y=sitelocs$y[match(meltQ$site, sitelocs$SITE.NAME)]
qxyz=make.xyz(x=meltQ$x,y=meltQ$y,z=meltQ$value,group=meltQ$variable)


### for plotting streams on map, need USA_Rivers_and_Streams shapefile (Google it!)

USstreams=readOGR(dsn="~/Desktop/kansasGIS/USA_Rivers_and_Streams", layer="USA_Rivers_and_Streams")
x=c(-104.83,-94.1)
y=c(36.08,39.52)
xy <- cbind(x,y)
S <- SpatialPoints(xy)
bbox(S)

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



darterbox=bbox(S)

gClip <- function(shp, bb){
  if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else b_poly <- as(extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}

USstreams_clipped <- gClip(USstreams,darterbox)
par(mai=c(0.1,0.1,0.1,0.1))
plot(USstreams_clipped,lwd=0.05,col="black")
map('state',fill=FALSE,add=TRUE,lwd=2,col='lightblue')
plot(ArkRiv,lwd=0.25,col="black",add=T)
plot(CimRiv,lwd=0.25,col="black",add=T)
plot(ChikRiv,lwd=0.25,col="black",add=T)
plot(MedRiv,lwd=0.25,col="black",add=T)
plot(SaltRiv,lwd=0.25,col="black",add=T)
plot(NFNRiv,lwd=0.25,col="black",add=T)
plot(SFNRiv,lwd=0.25,col="black",add=T)
plot(WalnutCrk,lwd=0.25,col="black",add=T)
plot(RattCrk,lwd=0.25,col="black",add=T)
plot(IllRiv,lwd=0.25,col="black",add=T)
plot(SpringRiv,lwd=0.25,col="black",add=T)
plot(BigSandy,lwd=0.25,col="black",add=T)
plot(Rush,lwd=0.25,col="black",add=T)

### Drawing pie charts on map

draw.pie(x=qxyz$x,y=qxyz$y,z=qxyz$z,radius=0.05,col=cols)

text(x=mean(subindsitessplit$'Upper Arkansas River'$long),y=c(mean(subindsitessplit$'Upper Arkansas River'$lat)-0.35),labels=expression(italic('Upper Arkansas River')),cex=0.7,srt=-18)

text(x=mean(subindsitessplit$'Middle Arkansas River'$long),y=c(mean(subindsitessplit$'Middle Arkansas River'$lat)-0.15),labels=expression(italic('Middle Arkansas River')),cex=0.7,srt=-5)

text(x=mean(subindsitessplit$'Big Sandy / Rush Creeks'$long),y=c(mean(subindsitessplit$'Big Sandy / Rush Creeks'$lat)-0.2),labels=expression(italic('Big Sandy / Rush Creeks')),cex=0.7,srt=-15)

text(x=mean(subindsitessplit$'Walnut Creek'$long),y=c(mean(subindsitessplit$'Walnut Creek'$lat)-0.15),labels=expression(italic('Walnut Creek')),cex=0.7,srt=-18)

text(x=mean(subindsitessplit$'Rattlesnake Creek'$long)-0.32,y=mean(subindsitessplit$'Rattlesnake Creek'$lat)-0.02,labels=expression(italic('Rattlesnake Creek')),cex=0.7,srt=40)

text(x=mean(subindsitessplit$'Lower Arkansas River'$long),y=mean(subindsitessplit$'Lower Arkansas River'$lat)+0.3,labels=expression(italic('Lower Arkansas River')),cex=0.7,srt=-28)

text(x=mean(subindsitessplit$'South Fork Ninnescah River'$long)+0.9,y=mean(subindsitessplit$'South Fork Ninnescah River'$lat)-0.55,labels=expression(italic('Ninnescah River')),cex=0.7,srt=-55)

text(x=mean(subindsitessplit$'Chikaskia River'$long)+0.3,y=mean(subindsitessplit$'Chikaskia River'$lat)-0.3,labels=expression(italic('Chikaskia River')),cex=0.7,srt=-45)

text(x=mean(subindsitessplit$'Medicine Lodge River'$long)+0.48,y=mean(subindsitessplit$'Medicine Lodge River'$lat)-0.45,labels=expression(italic('Medicine Lodge River')),cex=0.7,srt=-65)

text(x=mean(subindsitessplit$'Salt Fork Arkansas River'$long)+0.2,y=mean(subindsitessplit$'Salt Fork Arkansas River'$lat)-0.4,labels=expression(italic('Salt Fork Arkansas River')),cex=0.7,srt=-35)

text(x=mean(subindsitessplit$'Cimarron River'$long),y=mean(subindsitessplit$'Cimarron River'$lat)-0.2,labels=expression(italic('Cimarron River')),cex=0.7)

text(x=mean(subindsitessplit$'Illinois River'$long)-0.65,y=mean(subindsitessplit$'Illinois River'$lat),labels=expression(italic('Illinois River')),cex=0.7)

text(x=mean(subindsitessplit$'Spring River'$long)+0.3,y=mean(subindsitessplit$'Spring River'$lat)+0.2,labels=expression(italic('Spring River')),cex=0.7)

######## Scripts for plotting PCA to further visualize population structure and check for batch effects ########

### Sample table with site/batch information
infeaux=read.csv("~/Desktop/rapture/manuscript/revision/samptable_final.csv")


#### Covariance table from PCAngsd output
cov=read.table("~/Desktop/rapture/manuscript/revision/angsdout_rev/Ecr_adalign_indf_sel.cov")
cov$V1676=NULL
cov=cov[-1676,]
e=eigen(cov)
PC1 <- e$vectors[,1]
PC2 <- e$vectors[,2]
PC3 <- e$vectors[,3]
PC4 <- e$vectors[,4]
PC5 <- e$vectors[,5]
PC6 <- e$vectors[,6]


PCA<-data.frame(PC1,PC2,PC3,PC4,PC5,PC6,indsites$col,indsites$site)
PCA$batch=infeaux$Rapture_plate[match(PCA$indsites.site,infeaux$SiteID)]

par(mai=c(0.7,0.7,0.1,0.1),mfrow=c(2,2))
plot(PCA$PC1,PCA$PC2,col=as.character(PCA$indsites.col),xlab="PC1",ylab="PC2")
plot(PCA$PC3,PCA$PC4,col=as.character(PCA$indsites.col),xlab="PC3",ylab="PC4")
plot(PCA$PC5,PCA$PC6,col=as.character(PCA$indsites.col),xlab="PC5",ylab="PC6")
plot(x=NULL,y=NULL,xlim=c(0,0.5),ylim=c(0,1),xaxt=F,yaxt=F,xlab="",ylab="")
legend(0.1,0.9,legend=metalist$V1,col=metalist$V2,cex=0.7,pch=1,box.lwd=0)

par(mai=c(0.7,0.7,0.1,0.1),mfrow=c(2,2))
plot(PCA$PC1,PCA$PC2,col=PCA$batch,xlab="PC1",ylab="PC2")
plot(PCA$PC3,PCA$PC4,col=PCA$batch,xlab="PC3",ylab="PC4")
plot(PCA$PC5,PCA$PC6,col=PCA$batch,xlab="PC5",ylab="PC6")
plot(x=NULL,y=NULL,xlim=c(0,0.5),ylim=c(0,1),xaxt=F,yaxt=F,xlab="",ylab="")
legend(0.1,0.9,legend=c("Batch 1","Batch 2","Batch 3","Batch 4","Batch 5"),col=c(1,2,3,4,5),cex=0.7,pch=1,box.lwd=0)