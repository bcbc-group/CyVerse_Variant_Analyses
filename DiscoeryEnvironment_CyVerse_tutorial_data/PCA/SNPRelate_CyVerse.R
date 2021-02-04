setwd("~/PCA")

##For cluster files, use the below lines to make sure packages are installed and loaded
#if (!require("pacman")) install.packages("pacman")
#pacman::p_load(data.table, ggplot2, ggpmisc)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SNPRelate")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("VariantAnnotation")

#load necessary packages
library(gdsfmt)
library(SNPRelate)
library(VariantAnnotation)

#call VCF file
vcf.fn <- "Washingtonia_test_samples.recode.p.snps.vcf"

#convert to GDS format for SNPRelate analysis
snpgdsVCF2GDS(vcf.fn, "/home/rstudio/test.gds", method="biallelic.only", ignore.chr.prefix="PDK_30s")

#summarize data
snpgdsSummary("/home/rstudio/test.gds")
genofile <- snpgdsOpen("/home/rstudio/test.gds", readonly=FALSE)

#LD prune
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=1, maf=NaN, missing.rate=NaN, slide.max.n=5000,start.pos="random")
names(snpset)
snpset.id <- unlist(snpset)

#run initial PCA
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)
pdf("/home/rstudio/PCA_plot_eigenvectors.pdf")
plot(pca)
dev.off()

# determine the variance proportion (%) explained
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
rounded <- head(round(pc.percent, 2))

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the third eigenvector
                  EV4 = pca$eigenvect[,4],    # the third eigenvector
                  stringsAsFactors = FALSE)
head(tab)

#population code for color coding
pop_code <- add.gdsn(genofile, "popmap.txt")

#test plot
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

# Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- scan("pop.txt", what=character())
table(pop_code)
head(cbind(sample.id, pop_code))

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the third eigenvector
                  EV4 = pca$eigenvect[,4],    # the third eigenvector
stringsAsFactors = FALSE)
head(tab)


#color choice. Make sure the number of colors matches the number of states in the popmap
palette <- c("#FF7C26",  "#04B404", "#8B5600", "#9999FF", "#FFFB23")

#plot initial PC without colors
pdf("/home/rstudio/PCA_by_population.pdf")
plot(tab$EV1,tab$EV2, col=as.integer(tab$pop), xlab="eigenvector 1", ylab="eigenvector 2")
par(xpd=TRUE)
legend("topleft",ncol=2, cex=0.8, legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))
dev.off()

#plot PC 1 vs PC 2
pdf("/home/rstudio/PCA_unique_colors_PC1vsPC2.pdf",10,10)
par(mar=c(6,6,1.5,2))
col.list <- palette
plot(tab$EV1,tab$EV2, col=col.list[as.integer(tab$pop)], cex=4, pch=16,cex.lab=3,cex.axis=2, main="SNP PCA", xlab=paste0("eigenvector 1 (",rounded[1],"%)"), ylab=paste0("eigenvector 2 (",rounded[2],"%)"))
legend("topright", legend=levels(tab$pop), cex=1.6,pch=16, col=col.list[1:nlevels(tab$pop)])
dev.off()

#plot PC 2 vs PC 3
pdf("/home/rstudio/PCA_unique_colors_PC2vsPC3.pdf",10,10)
par(mar=c(6,6,1.5,2))
col.list <- palette
plot(tab$EV2,tab$EV3, col=col.list[as.integer(tab$pop)], cex=4, pch=16,cex.lab=3,cex.axis=2, main="SNP PCA", xlab=paste0("eigenvector 2 (",rounded[2],"%)"), ylab=paste0("eigenvector 3 (",rounded[3],"%)"))
legend("bottomleft", legend=levels(tab$pop), cex=1.6,pch=16, col=col.list[1:nlevels(tab$pop)])
dev.off()

#plot PC 3 vs PC 4
pdf("/home/rstudio/PCA_unique_colors_PC3vsPC4.pdf",10,10)
par(mar=c(6,6,1.5,2))
col.list <- palette
plot(tab$EV3,tab$EV4, col=col.list[as.integer(tab$pop)], cex=4, pch=16,cex.lab=3,cex.axis=2, main="SNP PCA", xlab=paste0("eigenvector 3 (",rounded[3],"%)"), ylab=paste0("eigenvector 4 (",rounded[4],"%)"))
legend("topleft", legend=levels(tab$pop), cex=1.6,pch=16, col=col.list[1:nlevels(tab$pop)])
dev.off()