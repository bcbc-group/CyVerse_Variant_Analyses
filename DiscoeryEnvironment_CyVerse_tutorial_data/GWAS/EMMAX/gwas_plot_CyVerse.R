setwd("~/EMMAX")

#For each plot, we will need to know position of SNPs
pos<-read.table("position.txt",stringsAsFactors = F)
#name variables for writing to a text document with top hits
chr <-pos[,1]
chrpos <- pos[,2]

#Color Chromosomes, colored alternating
col_chr<-pos[,1]
col_chr[col_chr%in%c(1,3,5,7,9)]<-"lightgrey"
col_chr[col_chr%in%c(2,4,6,8,10)]<-"darkgrey"

#for each trait used for GWAS, open .ps file (contains p-values)
#Trait 1
pvals<-read.table("manhattan.plot",stringsAsFactors = F,header=TRUE)
thresh<--log((0.05)/nrow(pvals),10)
phenosig<--log(pvals[,4],10)>thresh
thiscol<-col_chr
thiscol[phenosig]<-"darkblue"
pdf(file="/home/rstudio/Trait1_GWAS.pdf")
plot(-log(pvals[,4],10),pch=16,cex=0.7,col=thiscol,main="Trait 1 hits", xlab="Position along chromosome", ylab="log10 p-values") + abline(h=thresh,col="red")
dev.off()

#write text document with chromosomal position of hits
pv <- pvals[,4]
logPV <- -log10(pv)
output <- cbind(chr, chrpos, logPV)
cleaned <- na.omit(output)
write.table(cleaned, file="/home/rstudio/Trait1_positional_hits.txt", row.names=FALSE, quote=FALSE, sep='\t')


