### read in selection statistics and SNP site information from PCAngsd as well as file containing names and lengths of chromosomes 

### doing E cragini subsampled data first

arkRaptsel=read.csv("~/Desktop/rapture/manuscript/revision/sel/EcrRaptSampSBlast.csv",header=F)
arkRaptsites=read.table("~/Desktop/rapture/manuscript/revision/sel/EcrRaptSampSBlast.sites")
chrlen=read.table("~/Desktop/rapture/manuscript/nu_selectah/chrlengths.txt")

### Extract chromosome name and site information. Identify loci on baits putatively under selection
library(stringr)
sites1=str_split_fixed(arkRaptsites$V1, "_", 3)
sites2=str_split_fixed(arkRaptsites$V1, ";", 2)
arkRaptsel$chrom=sites2[,1]
arkRaptsel$site=sites1[,3]
arkRaptsel$type=c(rep("N",8492),rep("S",89))
cumlen=cumsum(chrlen)

### calculate P values from selection statistics and check if any are below significance threshold
arkRaptsig=subset(arkRaptsel, -log(pchisq(arkRaptsel$V1,df=1,lower.tail=FALSE),10) > -log10(0.05/8694))

### for Manhattan plotting, give each site a position = position on chromosome + length of previous chromosomes
chr1=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_2913"),]
chr2=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_476"),]
chr2$site=as.integer(chr2$site)+cumlen[1,]
chr3=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_1855"),]
chr3$site=as.integer(chr3$site)+cumlen[2,]
chr4=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_2954"),]
chr4$site=as.integer(chr4$site)+cumlen[3,]
chr5=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_4655"),]
chr5$site=as.integer(chr5$site)+cumlen[4,]
chr6=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_1670"),]
chr6$site=as.integer(chr6$site)+cumlen[5,]
chr7=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_3598"),]
chr7$site=as.integer(chr7$site)+cumlen[6,]
chr8=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_4661"),]
chr8$site=as.integer(chr8$site)+cumlen[7,]
chr9=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_537"),]
chr9$site=as.integer(chr9$site)+cumlen[8,]
chr10=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_906"),]
chr10$site=as.integer(chr10$site)+cumlen[9,]
chr11=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_1510"),]
chr11$site=as.integer(chr11$site)+cumlen[10,]
chr12=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_1453"),]
chr12$site=as.integer(chr12$site)+cumlen[11,]
chr13=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_3977"),]
chr13$site=as.integer(chr13$site)+cumlen[12,]
chr14=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_3998"),]
chr14$site=as.integer(chr14$site)+cumlen[13,]
chr15=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_2244"),]
chr15$site=as.integer(chr15$site)+cumlen[14,]
chr16=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_1570"),]
chr16$site=as.integer(chr16$site)+cumlen[15,]
chr17=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_4660"),]
chr17$site=as.integer(chr17$site)+cumlen[16,]
chr18=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_3430"),]
chr18$site=as.integer(chr18$site)+cumlen[17,]
chr19=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_978"),]
chr19$site=as.integer(chr19$site)+cumlen[18,]
chr20=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_4472"),]
chr20$site=as.integer(chr20$site)+cumlen[19,]
chr21=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_313"),]
chr21$site=as.integer(chr21$site)+cumlen[20,]
chr22=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_4666"),]
chr22$site=as.integer(chr22$site)+cumlen[21,]
chr23=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_868"),]
chr23$site=as.integer(chr23$site)+cumlen[22,]
chr24=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_4557"),]
chr24$site=as.integer(chr24$site)+cumlen[23,]

### open file and prep for plotting multiple selection graphs 
jpeg("~/Desktop/rapture/manuscript/revision/Ecrselection.jpg", width = 1600, height = 1200,res=300)
par(mfrow=c(3,1),mai=c(0.1,0.6,0.2,0.1))

### create manhattan plot with p-values on log scale
plot(x=NULL,y=NULL,xlim=c(0,cumlen[24,]),ylim=c(0,10),ylab=expression(-log[10]*italic("p")),xaxt="n",xlab="",main=expression('Rapture (subsampled) aligned to '*italic('E. cragini')), cex.main=1,cex.axis=0.7)
palette(c("gold","black"))
points(y=-log(pchisq(chr1[,1],df=1,lower.tail=FALSE),10),x=chr1$site,cex=0.05,pch=19,col=as.factor(chr1$type))
points(y=-log(pchisq(chr2[,1],df=1,lower.tail=FALSE),10),x=chr2$site,cex=0.05,pch=19,col=as.factor(chr2$type))
points(y=-log(pchisq(chr3[,1],df=1,lower.tail=FALSE),10),x=chr3$site,cex=0.05,pch=19,col=as.factor(chr3$type))
points(y=-log(pchisq(chr4[,1],df=1,lower.tail=FALSE),10),x=chr4$site,cex=0.05,pch=19,col=as.factor(chr4$type))
points(y=-log(pchisq(chr5[,1],df=1,lower.tail=FALSE),10),x=chr5$site,cex=0.05,pch=19,col=as.factor(chr5$type))
points(y=-log(pchisq(chr6[,1],df=1,lower.tail=FALSE),10),x=chr6$site,cex=0.05,pch=19,col=as.factor(chr6$type))
points(y=-log(pchisq(chr7[,1],df=1,lower.tail=FALSE),10),x=chr7$site,cex=0.05,pch=19,col=as.factor(chr7$type))
points(y=-log(pchisq(chr8[,1],df=1,lower.tail=FALSE),10),x=chr8$site,cex=0.05,pch=19,col=as.factor(chr8$type))
points(y=-log(pchisq(chr9[,1],df=1,lower.tail=FALSE),10),x=chr9$site,cex=0.05,pch=19,col=as.factor(chr9$type))
points(y=-log(pchisq(chr10[,1],df=1,lower.tail=FALSE),10),x=chr10$site,cex=0.05,pch=19,col=as.factor(chr10$type))
points(y=-log(pchisq(chr11[,1],df=1,lower.tail=FALSE),10),x=chr11$site,cex=0.05,pch=19,col=as.factor(chr11$type))
points(y=-log(pchisq(chr12[,1],df=1,lower.tail=FALSE),10),x=chr12$site,cex=0.05,pch=19,col=as.factor(chr12$type))
points(y=-log(pchisq(chr13[,1],df=1,lower.tail=FALSE),10),x=chr13$site,cex=0.05,pch=19,col=as.factor(chr13$type))
points(y=-log(pchisq(chr14[,1],df=1,lower.tail=FALSE),10),x=chr14$site,cex=0.05,pch=19,col=as.factor(chr14$type))
points(y=-log(pchisq(chr15[,1],df=1,lower.tail=FALSE),10),x=chr15$site,cex=0.05,pch=19,col=as.factor(chr15$type))
points(y=-log(pchisq(chr16[,1],df=1,lower.tail=FALSE),10),x=chr16$site,cex=0.05,pch=19,col=as.factor(chr16$type))
points(y=-log(pchisq(chr17[,1],df=1,lower.tail=FALSE),10),x=chr17$site,cex=0.05,pch=19,col=as.factor(chr17$type))
points(y=-log(pchisq(chr18[,1],df=1,lower.tail=FALSE),10),x=chr18$site,cex=0.05,pch=19,col=as.factor(chr18$type))
points(y=-log(pchisq(chr19[,1],df=1,lower.tail=FALSE),10),x=chr19$site,cex=0.05,pch=19,col=as.factor(chr19$type))
points(y=-log(pchisq(chr20[,1],df=1,lower.tail=FALSE),10),x=chr20$site,cex=0.05,pch=19,col=as.factor(chr20$type))
points(y=-log(pchisq(chr21[,1],df=1,lower.tail=FALSE),10),x=chr21$site,cex=0.05,pch=19,col=as.factor(chr21$type))
points(y=-log(pchisq(chr22[,1],df=1,lower.tail=FALSE),10),x=chr22$site,cex=0.05,pch=19,col=as.factor(chr22$type))
points(y=-log(pchisq(chr23[,1],df=1,lower.tail=FALSE),10),x=chr23$site,cex=0.05,pch=19,col=as.factor(chr23$type))
points(y=-log(pchisq(chr24[,1],df=1,lower.tail=FALSE),10),x=chr24$site,cex=0.05,pch=19,col=as.factor(chr24$type))

### plot chromosome breaks and significance threshold
abline(v=cumlen[1:23,],lty=2,lwd=0.2)
abline(h=-log10(0.05/8581), col="red",lty=2,lwd=0.2)

### Plot E cragini Rapture full dataset 

arkRaptsel=read.csv("~/Desktop/rapture/manuscript/revision/sel/EcrRaptSBlast.csv",header=F)
arkRaptsites=read.table("~/Desktop/rapture/manuscript/revision/sel/EcrRaptSBlast.sites")
library(stringr)
sites1=str_split_fixed(arkRaptsites$V1, "_", 3)
sites2=str_split_fixed(arkRaptsites$V1, ";", 2)
arkRaptsel$chrom=sites2[,1]
arkRaptsel$site=sites1[,3]
arkRaptsel$type=c(rep("N",8605),rep("S",89))
chrlen=read.table("~/Desktop/rapture/manuscript/nu_selectah/chrlengths.txt")
cumlen=cumsum(chrlen)

arkRaptsig=subset(arkRaptsel, -log(pchisq(arkRaptsel$V1,df=1,lower.tail=FALSE),10) > -log10(0.05/16401))

chr1=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_2913"),]
chr2=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_476"),]
chr2$site=as.integer(chr2$site)+cumlen[1,]
chr3=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_1855"),]
chr3$site=as.integer(chr3$site)+cumlen[2,]
chr4=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_2954"),]
chr4$site=as.integer(chr4$site)+cumlen[3,]
chr5=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_4655"),]
chr5$site=as.integer(chr5$site)+cumlen[4,]
chr6=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_1670"),]
chr6$site=as.integer(chr6$site)+cumlen[5,]
chr7=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_3598"),]
chr7$site=as.integer(chr7$site)+cumlen[6,]
chr8=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_4661"),]
chr8$site=as.integer(chr8$site)+cumlen[7,]
chr9=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_537"),]
chr9$site=as.integer(chr9$site)+cumlen[8,]
chr10=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_906"),]
chr10$site=as.integer(chr10$site)+cumlen[9,]
chr11=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_1510"),]
chr11$site=as.integer(chr11$site)+cumlen[10,]
chr12=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_1453"),]
chr12$site=as.integer(chr12$site)+cumlen[11,]
chr13=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_3977"),]
chr13$site=as.integer(chr13$site)+cumlen[12,]
chr14=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_3998"),]
chr14$site=as.integer(chr14$site)+cumlen[13,]
chr15=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_2244"),]
chr15$site=as.integer(chr15$site)+cumlen[14,]
chr16=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_1570"),]
chr16$site=as.integer(chr16$site)+cumlen[15,]
chr17=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_4660"),]
chr17$site=as.integer(chr17$site)+cumlen[16,]
chr18=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_3430"),]
chr18$site=as.integer(chr18$site)+cumlen[17,]
chr19=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_978"),]
chr19$site=as.integer(chr19$site)+cumlen[18,]
chr20=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_4472"),]
chr20$site=as.integer(chr20$site)+cumlen[19,]
chr21=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_313"),]
chr21$site=as.integer(chr21$site)+cumlen[20,]
chr22=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_4666"),]
chr22$site=as.integer(chr22$site)+cumlen[21,]
chr23=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_868"),]
chr23$site=as.integer(chr23$site)+cumlen[22,]
chr24=arkRaptsel[which(arkRaptsel$chrom=="ScbMSFa_4557"),]
chr24$site=as.integer(chr24$site)+cumlen[23,]

plot(x=NULL,y=NULL,xlim=c(0,cumlen[24,]),ylim=c(0,10),ylab=expression(-log[10]*italic("p")),xaxt="n",xlab="",main=expression('Rapture full dataset aligned to '*italic('E. cragini')), cex.main=1,cex.axis=0.7)
palette(c("lightblue","black"))
points(y=-log(pchisq(chr1[,1],df=1,lower.tail=FALSE),10),x=chr1$site,cex=0.05,pch=19,col=as.factor(chr1$type))
points(y=-log(pchisq(chr2[,1],df=1,lower.tail=FALSE),10),x=chr2$site,cex=0.05,pch=19,col=as.factor(chr2$type))
points(y=-log(pchisq(chr3[,1],df=1,lower.tail=FALSE),10),x=chr3$site,cex=0.05,pch=19,col=as.factor(chr3$type))
points(y=-log(pchisq(chr4[,1],df=1,lower.tail=FALSE),10),x=chr4$site,cex=0.05,pch=19,col=as.factor(chr4$type))
points(y=-log(pchisq(chr5[,1],df=1,lower.tail=FALSE),10),x=chr5$site,cex=0.05,pch=19,col=as.factor(chr5$type))
points(y=-log(pchisq(chr6[,1],df=1,lower.tail=FALSE),10),x=chr6$site,cex=0.05,pch=19,col=as.factor(chr6$type))
points(y=-log(pchisq(chr7[,1],df=1,lower.tail=FALSE),10),x=chr7$site,cex=0.05,pch=19,col=as.factor(chr7$type))
points(y=-log(pchisq(chr8[,1],df=1,lower.tail=FALSE),10),x=chr8$site,cex=0.05,pch=19,col=as.factor(chr8$type))
points(y=-log(pchisq(chr9[,1],df=1,lower.tail=FALSE),10),x=chr9$site,cex=0.05,pch=19,col=as.factor(chr9$type))
points(y=-log(pchisq(chr10[,1],df=1,lower.tail=FALSE),10),x=chr10$site,cex=0.05,pch=19,col=as.factor(chr10$type))
points(y=-log(pchisq(chr11[,1],df=1,lower.tail=FALSE),10),x=chr11$site,cex=0.05,pch=19,col=as.factor(chr11$type))
points(y=-log(pchisq(chr12[,1],df=1,lower.tail=FALSE),10),x=chr12$site,cex=0.05,pch=19,col=as.factor(chr12$type))
points(y=-log(pchisq(chr13[,1],df=1,lower.tail=FALSE),10),x=chr13$site,cex=0.05,pch=19,col=as.factor(chr13$type))
points(y=-log(pchisq(chr14[,1],df=1,lower.tail=FALSE),10),x=chr14$site,cex=0.05,pch=19,col=as.factor(chr14$type))
points(y=-log(pchisq(chr15[,1],df=1,lower.tail=FALSE),10),x=chr15$site,cex=0.05,pch=19,col=as.factor(chr15$type))
points(y=-log(pchisq(chr16[,1],df=1,lower.tail=FALSE),10),x=chr16$site,cex=0.05,pch=19,col=as.factor(chr16$type))
points(y=-log(pchisq(chr17[,1],df=1,lower.tail=FALSE),10),x=chr17$site,cex=0.05,pch=19,col=as.factor(chr17$type))
points(y=-log(pchisq(chr18[,1],df=1,lower.tail=FALSE),10),x=chr18$site,cex=0.05,pch=19,col=as.factor(chr18$type))
points(y=-log(pchisq(chr19[,1],df=1,lower.tail=FALSE),10),x=chr19$site,cex=0.05,pch=19,col=as.factor(chr19$type))
points(y=-log(pchisq(chr20[,1],df=1,lower.tail=FALSE),10),x=chr20$site,cex=0.05,pch=19,col=as.factor(chr20$type))
points(y=-log(pchisq(chr21[,1],df=1,lower.tail=FALSE),10),x=chr21$site,cex=0.05,pch=19,col=as.factor(chr21$type))
points(y=-log(pchisq(chr22[,1],df=1,lower.tail=FALSE),10),x=chr22$site,cex=0.05,pch=19,col=as.factor(chr22$type))
points(y=-log(pchisq(chr23[,1],df=1,lower.tail=FALSE),10),x=chr23$site,cex=0.05,pch=19,col=as.factor(chr23$type))
points(y=-log(pchisq(chr24[,1],df=1,lower.tail=FALSE),10),x=chr24$site,cex=0.05,pch=19,col=as.factor(chr24$type))
abline(v=cumlen[1:23,],lty=2,lwd=0.2)
abline(h=-log10(0.05/8694), col="red",lty=2,lwd=0.2)

### plot E cragini WGA dataset (this will take longer because of increased # of sites)

arkWGSsel=read.csv("~/Desktop/rapture/manuscript/revision/sel/EcrWGSSel.csv",header=F)
arkWGSsites=read.table("~/Desktop/rapture/manuscript/revision/sel/EcrWGSSel.sites",header=F)
library(stringr)
sites1=str_split_fixed(arkWGSsites$V1, "_", 3)
sites2=str_split_fixed(arkWGSsites$V1, ";", 2)
arkWGSsel$chrom=sites2[,1]
arkWGSsel$site=sites1[,3]
chrlen=read.table("~/Desktop/rapture/manuscript/nu_selectah/chrlengths.txt")
cumlen=cumsum(chrlen)

chr1=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_2913"),]
chr2=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_476"),]
chr2$site=as.integer(chr2$site)+cumlen[1,]
chr3=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_1855"),]
chr3$site=as.integer(chr3$site)+cumlen[2,]
chr4=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_2954"),]
chr4$site=as.integer(chr4$site)+cumlen[3,]
chr5=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_4655"),]
chr5$site=as.integer(chr5$site)+cumlen[4,]
chr6=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_1670"),]
chr6$site=as.integer(chr6$site)+cumlen[5,]
chr7=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_3598"),]
chr7$site=as.integer(chr7$site)+cumlen[6,]
chr8=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_4661"),]
chr8$site=as.integer(chr8$site)+cumlen[7,]
chr9=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_537"),]
chr9$site=as.integer(chr9$site)+cumlen[8,]
chr10=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_906"),]
chr10$site=as.integer(chr10$site)+cumlen[9,]
chr11=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_1510"),]
chr11$site=as.integer(chr11$site)+cumlen[10,]
chr12=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_1453"),]
chr12$site=as.integer(chr12$site)+cumlen[11,]
chr13=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_3977"),]
chr13$site=as.integer(chr13$site)+cumlen[12,]
chr14=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_3998"),]
chr14$site=as.integer(chr14$site)+cumlen[13,]
chr15=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_2244"),]
chr15$site=as.integer(chr15$site)+cumlen[14,]
chr16=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_1570"),]
chr16$site=as.integer(chr16$site)+cumlen[15,]
chr17=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_4660"),]
chr17$site=as.integer(chr17$site)+cumlen[16,]
chr18=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_3430"),]
chr18$site=as.integer(chr18$site)+cumlen[17,]
chr19=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_978"),]
chr19$site=as.integer(chr19$site)+cumlen[18,]
chr20=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_4472"),]
chr20$site=as.integer(chr20$site)+cumlen[19,]
chr21=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_313"),]
chr21$site=as.integer(chr21$site)+cumlen[20,]
chr22=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_4666"),]
chr22$site=as.integer(chr22$site)+cumlen[21,]
chr23=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_868"),]
chr23$site=as.integer(chr23$site)+cumlen[22,]
chr24=arkWGSsel[which(arkWGSsel$chrom=="ScbMSFa_4557"),]
chr24$site=as.integer(chr24$site)+cumlen[23,]


arkWGSsel$pPC1=(-log10(pchisq(arkWGSsel$V1,df=1)))
arkWGSsel$pPC2=(-log10(pchisq(arkWGSsel$V2,df=1)))
arkWGSsig1=subset(arkWGSsel,arkWGSsel$pPC1 > -log10(0.05/5759437))
arkWGSsig1a=subset(arkWGSsig1,arkWGSsig1$pPC1 < Inf)
arkWGSsig2=subset(arkWGSsel,arkWGSsel$pPC2 > -log10(0.05/5759437))
arkWGSsig2a=subset(arkWGSsig2,arkWGSsig1$pPC2 < Inf)


plot(x=NULL,y=NULL,xlim=c(0,cumlen[24,]),ylim=c(0,10),ylab=expression(-log[10]*italic("p")),xaxt="n",xlab="",main=expression('WGS aligned to '*italic('E. cragini')), cex.main=1,cex.axis=0.7,col="red")
points(y=-log(pchisq(chr1[,1],df=1,lower.tail=FALSE),10),x=chr1$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr2[,1],df=1,lower.tail=FALSE),10),x=chr2$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr3[,1],df=1,lower.tail=FALSE),10),x=chr3$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr4[,1],df=1,lower.tail=FALSE),10),x=chr4$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr5[,1],df=1,lower.tail=FALSE),10),x=chr5$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr6[,1],df=1,lower.tail=FALSE),10),x=chr6$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr7[,1],df=1,lower.tail=FALSE),10),x=chr7$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr8[,1],df=1,lower.tail=FALSE),10),x=chr8$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr9[,1],df=1,lower.tail=FALSE),10),x=chr9$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr10[,1],df=1,lower.tail=FALSE),10),x=chr10$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr11[,1],df=1,lower.tail=FALSE),10),x=chr11$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr12[,1],df=1,lower.tail=FALSE),10),x=chr12$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr13[,1],df=1,lower.tail=FALSE),10),x=chr13$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr14[,1],df=1,lower.tail=FALSE),10),x=chr14$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr15[,1],df=1,lower.tail=FALSE),10),x=chr15$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr16[,1],df=1,lower.tail=FALSE),10),x=chr16$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr17[,1],df=1,lower.tail=FALSE),10),x=chr17$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr18[,1],df=1,lower.tail=FALSE),10),x=chr18$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr19[,1],df=1,lower.tail=FALSE),10),x=chr19$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr20[,1],df=1,lower.tail=FALSE),10),x=chr20$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr21[,1],df=1,lower.tail=FALSE),10),x=chr21$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr22[,1],df=1,lower.tail=FALSE),10),x=chr22$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr23[,1],df=1,lower.tail=FALSE),10),x=chr23$site,cex=0.05,pch=19,col="red")
points(y=-log(pchisq(chr24[,1],df=1,lower.tail=FALSE),10),x=chr24$site,cex=0.05,pch=19,col="red")
abline(v=cumlen[1:23,],lty=2,lwd=0.2)
abline(h=-log10(0.05/5759437), col="red",lty=2,lwd=0.2)

### save plot
dev.off()