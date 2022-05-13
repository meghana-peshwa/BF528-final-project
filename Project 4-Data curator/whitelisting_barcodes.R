#import libraries
library(dplyr)
library(ggplot2)

#reading required data
data1 <- read.table("SRR3879604_counts")
data2 <- read.table("SRR3879605_counts")
data3 <- read.table("SRR3879606_counts")

combined_counts <- read.table("combined_counts")

#filtering out the barcodes that have length greater than 19
combined_counts <- combined_counts[order(combined_counts$V1), , drop=FALSE]
combined_counts <- combined_counts[nchar(combined_counts$V2) == 19, , drop=FALSE]

data1 <- data1[order(data1$V1), , drop=FALSE]
data1 <- data1[nchar(data1$V2) == 19, , drop=FALSE]

data2 <- data2[order(data2$V1), , drop=FALSE]
data2 <- data2[nchar(data2$V2) == 19, , drop=FALSE]

data3 <- data3[order(data3$V1), , drop=FALSE]
data3 <- data3[nchar(data3$V2) == 19, , drop=FALSE]

#Plotting CDF for all samples to determine the threshold for filtering barcodes
plot(ecdf(data1$V1), main='CDF plot - SRR3879604_1', ylab="Propotion",
     xlab='Barcode counts', xlim=c(0,500))
plot(ecdf(data2$V1), main='CDF plot - SRR3879605_1', ylab="Propotion",
     xlab='Barcode counts', xlim=c(0,500))
plot(ecdf(data3$V1), main='CDF plot - SRR3879606_1', ylab="Propotion",
     xlab='Barcode counts', xlim=c(0,500))
plot(ecdf(combined_counts$V1), main='CDF plot - Combined', ylab="Propotion",
     xlab='Barcode counts', xlim=c(0,500))

#knee seems to be near 30 so filter counts less than 30
filtered_combined_counts <- combined_counts[combined_counts$V1>30,]$V2
#writeLines(filtered_combined_counts,"whitelist.txt")