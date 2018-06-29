#!/usr/bin/env Rscript
## Calculate MDS components from identity-by-missingness distance matrix, cluster on the resulting components,
## plot the results and output a list of outliers
## CopyLeft 2013 Mitja Mitrovic, Cotsapas Lab, Yale School of Medicine
## Usage: Rscript cluster.ibm.R infile.fam

## Take the input argument and chop off the file extension
args <- commandArgs(TRUE)
argument1 <- paste(args[1], sep="")
filename <- data.frame(strsplit(argument1, "\\.fam"))[1,]
matrixfilename <- paste(filename,".mdist.missing", sep="")

## Load the required libraries
library(mclust)
library(ggplot2)

## Read-in the IBM similarity matrix and convert similarities to dissimilarities 
similarity.matrix <- read.table(matrixfilename, quote = "", header = F)
distance.matrix <- 1 - similarity.matrix

## Calculate first four MDS solutions of the distance matrix
mds <- cmdscale(distance.matrix, k=4)
colnames(mds) <- c("C1", "C2","C3","C4") 

## Merge FID and SID from the fam file with mds data frame into a new data frame
fam.file <- read.table(args[1], quote = "", header = F)
colnames(fam.file) <- c("FID", "SID","MID","PID", "SEX", "PHENO") 
full.table <- data.frame(FID = fam.file$FID, SID = fam.file$SID, mds)

## Remove any missing data (i.e. NAs/nans etc)
#full.table <- na.omit(full.table)

## Cluster samples based on first two/four MDS components
m <- full.table[,c("C1", "C2")]
#m <- full.table[,c("C1", "C2", "C3","C4")] # if one wants to cluster on four components 
m_Mclust <- Mclust(m, G=2)

## Assign cluster classification number to the samples and output a table
full.table$CLUSTER <- factor(m_Mclust$classification)
filename1 <- paste(filename,".clustering.table",".txt", sep="")
write.table(full.table, file= filename1, quote=FALSE, row.names = F)

## Plot the first four components in pairs
filename2 <- paste(filename,".ibm",".pairs",".png", sep="")
png(filename2,width=1000,height=705)
pairs(full.table[3:6], pch=20,cex=0.7)
dev.off()

## Count occurrences of "1" and "2" in the "CLUSTER" column
ones <- subset(full.table, CLUSTER == 1)
n.ones<- round((nrow(ones)), digits = 0)
twos <- subset(full.table, CLUSTER == 2)
n.twos <- round((nrow(twos)), digits = 0)

## Assuming there is more "OK samples" than outilers in the data,
## plot the 1st two components and show the outliers

filename3 <- paste(filename,".ibm",".png", sep="")

if (n.ones > n.twos) {
  good.samples <- subset(full.table, CLUSTER == 1)
  bad.samples <- subset(full.table, CLUSTER == 2)
  } else {
  good.samples <- subset(full.table, CLUSTER == 2)
  bad.samples <- subset(full.table, CLUSTER == 1)
}

if (n.ones > n.twos) {
  png(filename3,width=1000,height=705)
  ggplot(data=full.table, aes(x=C1, y=C2, color=CLUSTER))+ 
    geom_point()+
    xlab("C1")+
    ylab("C2")+
    scale_color_manual(values=c("limegreen","brown1"),labels=c("OK samples","Outliers"))
  } else {
  png(filename3,width=1000,height=705)
  ggplot(data=full.table, aes(x=C1, y=C2, color=CLUSTER))+ 
    geom_point()+
    xlab("C1")+
    ylab("C2")+
    scale_color_manual(values=c("brown1","limegreen"),labels=c("Outliers","OK samples"))
}
dev.off()

## Output a PLINK-friendly list of outliers  
outliers <- subset(bad.samples, select = c(FID,SID))
filename4 <- paste(filename,".ibm.outliers",".txt",sep="")
write.table(outliers, file=filename4, quote=FALSE, row.names = F, col.names = F)

print("\"Always pass on what you have learned\" - Yoda", justify = "right", quote = F, row.names = F, col.names = F)
