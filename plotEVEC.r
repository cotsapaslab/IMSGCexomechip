#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

C <- read.table(args[1],header=FALSE) #,check.names=FALSE)

png(args[2],width = 1000,height = 1000)

layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(0.7,0.3), heights=c(1))

plot(C$V2,C$V3,col=C$V12,xlab="PC1",ylab="PC2",type="n")
minX <- min(C$V2)
maxX <- max(C$V2)
minY <- min(C$V3)
maxY <- max(C$V3)
l <- levels(C$V12)

colors <- rep(c("blue","red","darkgreen","purple","black","orange"),10)
symbols <- c(rep(20,6),rep(8,6),rep(3,6),rep(7,6),rep(5,6),rep(11,6),rep(10,6),rep(13,6))
#symbols <- c(rep(12,5), rep(2,5),rep(5,5),rep(1,5),rep(13,5))
for (i in 1:length(l)){
    K <- subset(C, C$V12 == l[i])
    points(K$V2, K$V3, col = colors[i],pch = symbols[i])
}
abline(v= 0.006,lwd=1.2,lty=2,col="red")
abline(v= -0.02,lwd=1.2,lty=2,col="red")
abline(h=-0.027,lwd=1.2,lty=2,col="red")
abline(h=0.02,lwd=1.2,lty=2,col="red")
plot(1:3, rnorm(3), pch = 1, lty = 1, ylim=c(-2,2), type = "n", axes = FALSE, ann = FALSE)
legend("center",levels(C$V12),col=colors[1:length(l)],pch=symbols[1:length(l)],cex=1.5)
dev.off()

png(args[3],width = 1000,height = 1000)

layout(matrix(c(1,2,5,3,4,5), 2, 3, byrow = TRUE),widths=c(1,1,0.5), heights=c(1,1))

plot(C$V2,C$V3,col=C$V12,xlab="PC1",ylab="PC2",type="n")
minX <- min(C$V2)
maxX <- max(C$V2)
minY <- min(C$V3)
maxY <- max(C$V3)
l <- levels(C$V12)

#colors <- rep(c("blue","red","darkgreen","purple","black","orange"),4)
#symbols <- c(rep(20,6),rep(8,6),rep(3,6),rep(7,6),rep(5,6),rep(11,6),rep(10,6))
#symbols <- c(rep(12,5), rep(2,5),rep(5,5),rep(1,5))
for (i in 1:length(l)){
    K <- subset(C, C$V12 == l[i])
    points(K$V2, K$V3, col = colors[i],pch = symbols[i])
}

plot(C$V2,C$V4,col=C$V12,xlab="PC1",ylab="PC3",type="n")
minX <- min(C$V2)
maxX <- max(C$V2)
minY <- min(C$V4)
maxY <- max(C$V4)
l <- levels(C$V12)

#colors <- rep(c("blue","red","darkgreen","purple","black","orange"),4)
#symbols <- c(rep(20,6),rep(8,6),rep(3,6),rep(7,6),rep(5,6),rep(11,6),rep(10,6))
#symbols <- c(rep(12,5), rep(2,5),rep(5,5),rep(1,5))
for (i in 1:length(l)){
    K <- subset(C, C$V12 == l[i])
    points(K$V2, K$V4, col = colors[i],pch = symbols[i])
}

plot(C$V2,C$V5,col=C$V12,xlab="PC1",ylab="PC4",type="n")
minX <- min(C$V2)
maxX <- max(C$V2)
minY <- min(C$V5)
maxY <- max(C$V5)
l <- levels(C$V12)

#colors <- rep(c("blue","red","darkgreen","purple","black","orange"),4)
#symbols <- c(rep(20,6),rep(8,6),rep(3,6),rep(7,6),rep(5,6),rep(11,6),rep(10,6))
#symbols <- c(rep(12,5), rep(2,5),rep(5,5),rep(1,5))
for (i in 1:length(l)){
    K <- subset(C, C$V12 == l[i])
    points(K$V2, K$V5, col = colors[i],pch = symbols[i])
}

plot(C$V3,C$V4,col=C$V12,xlab="PC2",ylab="PC3",type="n")
minX <- min(C$V3)
maxX <- max(C$V3)
minY <- min(C$V4)
maxY <- max(C$V4)
l <- levels(C$V12)

#colors <- rep(c("blue","red","darkgreen","purple","black","orange"),4)
#symbols <- c(rep(20,6),rep(8,6),rep(3,6),rep(7,6),rep(5,6),rep(11,6),rep(10,6))
#symbols <- c(rep(12,5), rep(2,5),rep(5,5),rep(1,5))
for (i in 1:length(l)){
    K <- subset(C, C$V12 == l[i])
    points(K$V3, K$V4, col = colors[i],pch = symbols[i])
}

plot(1:3, rnorm(3), pch = 1, lty = 1, ylim=c(-2,2), type = "n", axes = FALSE, ann = FALSE)
legend("center",levels(C$V12),col=colors[1:length(l)],pch=symbols[1:length(l)],cex=1.5)
dev.off()


n <- length(levels(C$V12))
c <- ceiling(n/4)

png(args[4],width=2000,c*500)

layout(matrix(seq(1,c*4,1), c, 4, byrow = TRUE),widths=c(1,1,1,1), heights=rep(1,c))

minX <- min(C$V2)
maxX <- max(C$V2)
minY <- min(C$V3)
maxY <- max(C$V3)
l <- levels(C$V12)

#colors <- rep(c("blue","red","darkgreen","purple","cyan","orange"),7)
#symbols <- c(rep(20,6),rep(8,6),rep(3,6),rep(7,6),rep(5,6),rep(11,6),rep(10,6))
#symbols <- c(rep(12,6), rep(2,5),rep(5,5),rep(1,5))
for (i in 1:length(l)){
    K <- subset(C, C$V12 == l[i])
    plot(K$V2, K$V3, col = colors[i],pch = symbols[i],xlab="PC1",ylab="PC2",main=l[i],xlim=c(minX,maxX),ylim=c(minY,maxY))
}
dev.off()