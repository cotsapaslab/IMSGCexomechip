#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

X <- read.table(args[1],header=TRUE,sep="\t")
Y <- read.table(args[2],header=TRUE,sep="\t")

png(args[3],width=800,height=400)
layout(matrix(c(1,2),nrow=1))
plot(Y$percHET,X$percMISS,xlab="Rare - % Het Calls", ylab="Common - % No Calls",main="All Points")
plot(Y$percHET,X$percMISS,ylim=c(0,2),xlim=c(0.1,1.2),xlab="Rare - % Het Calls", ylab="Common - % No Calls",main="Zoomed In",cex=0.3,pch=20)
abline(h=args[6],col="red")
abline(v = args[4],col="red")
abline(v=args[5],col="red")
dev.off()