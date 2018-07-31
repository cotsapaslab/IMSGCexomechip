#### INFO ####
# This script makes PCAs plots for the any dataset with binary phenotype that is projected on HapMap3.
# USAGE: R CMD BATCH -evec -pops -out ~nikolaos/codes/R/pcaplot1KGbin.R
# LAST UPDATE: 2010-06-09 01:01:21
#############

## READS input options
rm(list=ls())
file=commandArgs()[7]
file=substr(file, 2, nchar(file))

pops=commandArgs()[8]
pops=substr(pops, 2, nchar(pops))

out=commandArgs()[9]
out=substr(out,2,nchar(out))


## Load the data for the distribution of pairs
data1 <- read.table(file, header=F)
data2 <- read.table(pops, header=T)
names(data1)[1]<-"IID"
data <- merge(data1, data2, by="IID") 

# subset CEU and start making the plots for the 10 first PCAs

acb <- subset(data, data$POP=="ASB")
asw <- subset(data, data$POP=="ASW")
beb <- subset(data, data$POP=="BEB")
cdx <- subset(data, data$POP=="CDX")
ceu <- subset(data, data$POP=="CEU")
chb <- subset(data, data$POP=="CHB")
chs <- subset(data, data$POP=="CHS")
clm <- subset(data, data$POP=="CLM")
esn <- subset(data, data$POP=="ESN")
fin <- subset(data, data$POP=="FIN")
gbr <- subset(data, data$POP=="GBR")
gih <- subset(data, data$POP=="GIH")
gwd <- subset(data, data$POP=="GWD")
ibs <- subset(data, data$POP=="IBS")
itu <- subset(data, data$POP=="ITU")
jpt <- subset(data, data$POP=="JPT")
khv <- subset(data, data$POP=="KHV")
lwk <- subset(data, data$POP=="LWK")
msl <- subset(data, data$POP=="MSL")
mxl <- subset(data, data$POP=="MXL")
pel <- subset(data, data$POP=="PEL")
pjl <- subset(data, data$POP=="PJL")
pur <- subset(data, data$POP=="PUR")
stu <- subset(data, data$POP=="STU")
tsi <- subset(data, data$POP=="TSI")
yri <- subset(data, data$POP=="YRI")
CASES <- subset(data, data$POP=="DATA" & TYPE==2 )
CTRLS <- subset(data, data$POP=="DATA" & TYPE==1)


colors <- c("green", "green1", "red", "red1", "blue", "red2", "red3", "yellow", "green2", "blue2", "blue3", "red4", "green3", "blue4", "orangered", "orange", "orange1", "green4", "lawngreen", "yellow2", "yellow3", "orange2", "yellow4", "orange3", "dodgerblue4", "limegreen", "#ff0000", "#ff0000")
groups <- c("ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI", "CASES", "CTRLS"  )
symbols <- c()
symbols[1:26] <- 19
symbols[27:28] <- c(5, 3)
for(i in 1:9){
	pdf(paste(out,"_PCA",i,"_PCA",i+1,".pdf", sep=""))
	nf<- layout(matrix(c(0,0,1,2),2,2,byrow=TRUE), c(3,1),c(0,4), TRUE)
	xrange <- range(eval(parse(text=paste("data$V",i+1,sep=""))))
	yrange <- range(eval(parse(text=paste("data$V",i+2,sep=""))))
	# plot everything in white to set the axes correctly
	par(mar=c(5.1,4.1,4.1,2.1))
	plot(eval(parse(text=paste("data$",i+1,sep="V"))), eval(parse(text=paste("data$",i+2,sep="V"))), col="#ffffff",  xlab=paste("PCA",i,sep=""), ylab=paste("PCA",i+1,sep=""), xlim=xrange, ylim=yrange)
	points(eval(parse(text=paste("acb$",i+1,sep="V"))), eval(parse(text=paste("acb$",i+2,sep="V"))), col=colors[1], pch=1, xlab="", ylab="")
        points(eval(parse(text=paste("asw$",i+1,sep="V"))), eval(parse(text=paste("asw$",i+2,sep="V"))), col=colors[2], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("beb$",i+1,sep="V"))), eval(parse(text=paste("beb$",i+2,sep="V"))), col=colors[3], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("cdx$",i+1,sep="V"))), eval(parse(text=paste("cdx$",i+2,sep="V"))), col=colors[4], pch=1, xlab="", ylab="")
        points(eval(parse(text=paste("ceu$",i+1,sep="V"))), eval(parse(text=paste("ceu$",i+2,sep="V"))), col=colors[5], pch=1, xlab="", ylab="")
        points(eval(parse(text=paste("chb$",i+1,sep="V"))), eval(parse(text=paste("chb$",i+2,sep="V"))), col=colors[6], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("chs$",i+1,sep="V"))), eval(parse(text=paste("chs$",i+2,sep="V"))), col=colors[7], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("clm$",i+1,sep="V"))), eval(parse(text=paste("clm$",i+2,sep="V"))), col=colors[8], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("esn$",i+1,sep="V"))), eval(parse(text=paste("esn$",i+2,sep="V"))), col=colors[9], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("fin$",i+1,sep="V"))), eval(parse(text=paste("fin$",i+2,sep="V"))), col=colors[10], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("gbr$",i+1,sep="V"))), eval(parse(text=paste("gbr$",i+2,sep="V"))), col=colors[11], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("gih$",i+1,sep="V"))), eval(parse(text=paste("gih$",i+2,sep="V"))), col=colors[12], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("gwd$",i+1,sep="V"))), eval(parse(text=paste("gwd$",i+2,sep="V"))), col=colors[13], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("ibs$",i+1,sep="V"))), eval(parse(text=paste("ibs$",i+2,sep="V"))), col=colors[14], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("itu$",i+1,sep="V"))), eval(parse(text=paste("itu$",i+2,sep="V"))), col=colors[15], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("jpt$",i+1,sep="V"))), eval(parse(text=paste("jpt$",i+2,sep="V"))), col=colors[16], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("khv$",i+1,sep="V"))), eval(parse(text=paste("khv$",i+2,sep="V"))), col=colors[17], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("lwk$",i+1,sep="V"))), eval(parse(text=paste("lwk$",i+2,sep="V"))), col=colors[18], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("msl$",i+1,sep="V"))), eval(parse(text=paste("msl$",i+2,sep="V"))), col=colors[19], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("mxl$",i+1,sep="V"))), eval(parse(text=paste("mxl$",i+2,sep="V"))), col=colors[20], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("pel$",i+1,sep="V"))), eval(parse(text=paste("pel$",i+2,sep="V"))), col=colors[21], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("pjl$",i+1,sep="V"))), eval(parse(text=paste("pjl$",i+2,sep="V"))), col=colors[22], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("pur$",i+1,sep="V"))), eval(parse(text=paste("pur$",i+2,sep="V"))), col=colors[23], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("stu$",i+1,sep="V"))), eval(parse(text=paste("stu$",i+2,sep="V"))), col=colors[24], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("tsi$",i+1,sep="V"))), eval(parse(text=paste("tsi$",i+2,sep="V"))), col=colors[25], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("yri$",i+1,sep="V"))), eval(parse(text=paste("yri$",i+2,sep="V"))), col=colors[26], pch=1, xlab="", ylab="")
	points(eval(parse(text=paste("CASES$",i+1,sep="V"))), eval(parse(text=paste("CASES$",i+2,sep="V"))), col="#ff0000", xlab="", ylab="", pch=5)
	points(eval(parse(text=paste("CTRLS$",i+1,sep="V"))), eval(parse(text=paste("CTRLS$",i+2,sep="V"))), col="#ff0000", xlab="", ylab="", pch=3)
	# add a legend
	par(mar=c(0,0,0,0))
	plot(1,1, axes=FALSE, col=rgb(0, 0, 255, 0, maxColorValue=255),  xlab="", ylab="")
	legend("center", groups , cex=0.8, title="Populations", pch=symbols, col=colors)
	dev.off()
}


