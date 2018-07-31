#### INFO ####
# This script makes PCAs plots for datasets with binary phenotypes.
# USAGE: R CMD BATCH -evec -out ~nikolaos/codes/R/pcaplotBin.R
# LAST UPDATE: 2010-06-09 01:01:21
#############

## READS input options
rm(list=ls())
file=commandArgs()[7]
file=substr(file, 2, nchar(file))

out=commandArgs()[8]
out=substr(out,2,nchar(out))


## Load the data for the distribution of pairs
data <- read.table(file, header=F)

# get the column number with the information of case/control status
n<-ncol(data)
# subset cases/controls and make plots for the 10 first PCAs
CASES <- subset(data[,2:(n-1)], data[, n]=="Case")
CTRLS <- subset(data[,2:(n-1)], data[, n]=="Control")

## make a vector with the color scheme
colors <- c("blue", "red")
## vector with names
groups <- c("Cases", "Controls")
#nf<- layout(matrix(c(0,0,1,2),2,2,byrow=TRUE), c(3,1),c(0,4), TRUE)
for(i in 1:9){
	pdf(paste(out,"_PCA",i,"_PCA",i+1,".pdf", sep=""))
	nf<- layout(matrix(c(0,0,1,2),2,2,byrow=TRUE), c(3,1),c(0,4), TRUE)
	xrange <- range(eval(parse(text=paste("data$V",i+1,sep=""))))
	yrange <- range(eval(parse(text=paste("data$V",i+2,sep=""))))
	# plot everything in white to set the axes correctly
	par(mar=c(5.1,4.1,4.1,2.1))
	plot(eval(parse(text=paste("data$",i+1,sep="V"))), eval(parse(text=paste("data$",i+2,sep="V"))), col="#ffffff",  xlab=paste("PCA",i,sep=""), ylab=paste("PCA",i+1,sep=""), xlim=xrange, ylim=yrange)
	points(eval(parse(text=paste("CASES$",i+1,sep="V"))), eval(parse(text=paste("CASES$",i+2,sep="V"))), col=colors[1], xlab="", ylab="", pch=5)
	points(eval(parse(text=paste("CTRLS$",i+1,sep="V"))), eval(parse(text=paste("CTRLS$",i+2,sep="V"))), col=colors[2], xlab="", ylab="", pch=3)
	# add a legend
	par(mar=c(0,0,0,0))
	plot(1,1, axes=FALSE, col=rgb(0, 0, 255, 0, maxColorValue=255),  xlab="", ylab="")
	legend("center", groups , cex=0.8, title="Populations", pch=c(5, 3), col=colors)
	dev.off()
}

## Logistic regression to test statistical significance of the PCAs
data1<-data[,2:n]
data1[,(n-1)]<-ifelse(data1[,(n-1)]=="Case", 1, 0)
data2<-data1[,1:(n-2)]
for(i in 1:ncol(data2)){
 names(data2)[i]<-paste("PCA", i, sep="")
}
l0<-glm(data1[,(n-1)] ~ 1, data=data2, family="binomial")
l<-glm(data1[,(n-1)] ~ ., data=data2, family="binomial")
lf<-step(l0, scope = list( upper=l, lower=~1 ), direction = "forward", trace=FALSE)
write.table(summary(lf)$coefficients, file=paste(out,"regression","text",sep="."), col.names=F, quote=F)
