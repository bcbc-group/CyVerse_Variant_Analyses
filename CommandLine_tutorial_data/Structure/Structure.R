setwd("/Users/jacoblandis/Desktop/Webinar_2021/data/Structure")

#you will need to install this packages if its the first time on a new computer
#if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#BiocManager::install("LEA")

#install.packages(c("fields","mapplots"))
#source("http://bioconductor.org/biocLite.R")
#biocLite("LEA")


library(fields)
library(LEA)

#additional files need from the structure format
source("Conversion.R")
source("POPSutilities.R")

#import input file, if this is coming directly from Stacks you will need to delete the first two rows of the file. Should be left with just the data, no other information
struct2geno(file = "Washingtonia_test_samples.recode.p.structure", TESS = FALSE, diploid = TRUE, FORMAT = 2,extra.row = 0, extra.col = 0, output = "Wash.geno")


#test for best K
obj.snmf = snmf("Wash.geno", K = 1:25,  ploidy = 2, entropy = T, alpha = 100, project = "new")
pdf(file="delta_K.pdf")
plot(obj.snmf, col = "blue4", cex = 2.0, pch = 19,cex.axis=1.4,cex.lab=1.4)
dev.off()


##########
### K = 2
##########
obj.snmf = snmf("Wash.geno", K = 2, alpha = 100, project = "new", iterations=500)
qmatrix = Q(obj.snmf, K = 2)
palette4 <- c("#FF9326","#FFFB23")
names <- c("NACA14","NACA16","NACA18","CACA11","CACA17","CACA18","COMO11","COMO14","COMO18","MULE16","MULE19","MULE20","ROFO12","ROFO16","ROFO17","SJD11","SJD19","SJD20","BERRE12","BERRE15","BERRE19","BOCA14","BOCA19","BOCA20","BOR13","BOR14","BOR18","CATA1","CATA5","CATA9","STAMA12","STAMA19","STAMA20","JOSH14","JOSH18","JOSH19","KOFA12","KOFA13","KOFA18","PALM11","PALM17","PALM18","PALO13","PALO16","PALO20","ALA13","ALA15","ALA18")
#plot classic structure diagram
pdf(file="K_equals_2.pdf",10,5)
barplot(t(qmatrix), col = c(palette4), border = NA, space = 0,names.arg=names,las=2,cex.names=0.7,
        xlab = "Individuals", ylab = "Admixture coefficients")
dev.off()


#######
## K = 3
#######
obj.snmf = snmf("Wash.geno", K = 3, alpha = 100, project = "new", iterations=500)
qmatrix = Q(obj.snmf, K = 3)
palette4 <- c("#FFFB23","#FF9326","#DF0101")
names <- c("NACA14","NACA16","NACA18","CACA11","CACA17","CACA18","COMO11","COMO14","COMO18","MULE16","MULE19","MULE20","ROFO12","ROFO16","ROFO17","SJD11","SJD19","SJD20","BERRE12","BERRE15","BERRE19","BOCA14","BOCA19","BOCA20","BOR13","BOR14","BOR18","CATA1","CATA5","CATA9","STAMA12","STAMA19","STAMA20","JOSH14","JOSH18","JOSH19","KOFA12","KOFA13","KOFA18","PALM11","PALM17","PALM18","PALO13","PALO16","PALO20","ALA13","ALA15","ALA18")
#plot classic structure diagram
pdf(file="K_equals_3.pdf",10,5)
barplot(t(qmatrix), col = c(palette4), border = NA, space = 0,names.arg=names,las=2,cex.names=0.7,
        xlab = "Individuals", ylab = "Admixture coefficients")
dev.off()


orginal_palette4 <- c("#FFFB23","#FF9326","#DF0101","#A945FF")

#run with best K of 4
obj.snmf = snmf("Wash.geno", K = 4, alpha = 100, project = "new", iterations=500)
qmatrix = Q(obj.snmf, K = 4)
palette4 <- c("#FF9326","#DF0101","#FFFB23","#A945FF")
names <- c("NACA14","NACA16","NACA18","CACA11","CACA17","CACA18","COMO11","COMO14","COMO18","MULE16","MULE19","MULE20","ROFO12","ROFO16","ROFO17","SJD11","SJD19","SJD20","BERRE12","BERRE15","BERRE19","BOCA14","BOCA19","BOCA20","BOR13","BOR14","BOR18","CATA1","CATA5","CATA9","STAMA12","STAMA19","STAMA20","JOSH14","JOSH18","JOSH19","KOFA12","KOFA13","KOFA18","PALM11","PALM17","PALM18","PALO13","PALO16","PALO20","ALA13","ALA15","ALA18")
#plot classic structure diagram
pdf(file="K_equals_4.pdf",7,3)
barplot(t(qmatrix), col = c(palette4), border = NA, space = 0,names.arg=names,las=2,cex.names=0.7,
        xlab = "Individuals", ylab = "Admixture coefficients")
dev.off()

#####
# K = 5
#####
obj.snmf = snmf("Wash.geno", K = 5, alpha = 100, project = "new", iterations=500)
qmatrix = Q(obj.snmf, K = 5)

palette5 <- c("#A945FF","#FFFB23","#FF9326","#DF0101","#9999FF")
names <- c("NACA14","NACA16","NACA18","CACA11","CACA17","CACA18","COMO11","COMO14","COMO18","MULE16","MULE19","MULE20","ROFO12","ROFO16","ROFO17","SJD11","SJD19","SJD20","BERRE12","BERRE15","BERRE19","BOCA14","BOCA19","BOCA20","BOR13","BOR14","BOR18","CATA1","CATA5","CATA9","STAMA12","STAMA19","STAMA20","JOSH14","JOSH18","JOSH19","KOFA12","KOFA13","KOFA18","PALM11","PALM17","PALM18","PALO13","PALO16","PALO20","ALA13","ALA15","ALA18")
#plot classic structure diagram
pdf(file="K_equals_5.pdf",10,5)
barplot(t(qmatrix), col = c(palette5), border = NA, space = 0,names.arg=names,las=2,cex.names=0.7,
        xlab = "Individuals", ylab = "Admixture coefficients")
dev.off()

#####
# K = 6
#####
obj.snmf = snmf("Wash.geno", K = 6, alpha = 100, project = "new", iterations=500)
qmatrix = Q(obj.snmf, K = 6)

palette6 <- c("#FF9326","#A945FF","#FFFB23","#9999FF","#04B404","#DF0101")
names <- c("NACA14","NACA16","NACA18","CACA11","CACA17","CACA18","COMO11","COMO14","COMO18","MULE16","MULE19","MULE20","ROFO12","ROFO16","ROFO17","SJD11","SJD19","SJD20","BERRE12","BERRE15","BERRE19","BOCA14","BOCA19","BOCA20","BOR13","BOR14","BOR18","CATA1","CATA5","CATA9","STAMA12","STAMA19","STAMA20","JOSH14","JOSH18","JOSH19","KOFA12","KOFA13","KOFA18","PALM11","PALM17","PALM18","PALO13","PALO16","PALO20","ALA13","ALA15","ALA18")
#plot classic structure diagram
pdf(file="K_equals_6.pdf",10,5)
barplot(t(qmatrix), col = c(palette6), border = NA, space = 0,names.arg=names,las=2,cex.names=0.7,
        xlab = "Individuals", ylab = "Admixture coefficients")
dev.off()

#####
# K = 7
#####
obj.snmf = snmf("Wash.geno", K = 7, alpha = 100, project = "new", iterations=500)
qmatrix = Q(obj.snmf, K = 7)

#original
#palette7 <- c("#FF9326", "#04B404", "#2121D9", "#9999FF", "#FFFB23", "#DF0101", "#A945FF")


palette7 <- c("#04B404","#9999FF","#A945FF","#2121D9","#FFFB23","#FF9326","#DF0101")
names <- c("NACA14","NACA16","NACA18","CACA11","CACA17","CACA18","COMO11","COMO14","COMO18","MULE16","MULE19","MULE20","ROFO12","ROFO16","ROFO17","SJD11","SJD19","SJD20","BERRE12","BERRE15","BERRE19","BOCA14","BOCA19","BOCA20","BOR13","BOR14","BOR18","CATA1","CATA5","CATA9","STAMA12","STAMA19","STAMA20","JOSH14","JOSH18","JOSH19","KOFA12","KOFA13","KOFA18","PALM11","PALM17","PALM18","PALO13","PALO16","PALO20","ALA13","ALA15","ALA18")
#plot classic structure diagram
pdf(file="K_equals_7.pdf",10,5)
barplot(t(qmatrix), col = c(palette7), border = NA, space = 0,names.arg=names,las=2,cex.names=0.7,
        xlab = "Individuals", ylab = "Admixture coefficients")
dev.off()

#####
# K = 8
#####
obj.snmf = snmf("Wash.geno", K = 8, alpha = 100, project = "new", iterations=500)
qmatrix = Q(obj.snmf, K = 8)

#original
#palette7 <- c("#FF9326", "#04B404", "#2121D9", "#9999FF", "#FFFB23", "#DF0101", "#A945FF")

palette8 <- c("#9999FF","#FFFB23","#A945FF","#2121D9","#FF9326","gray","#04B404","#DF0101")
names <- c("NACA14","NACA16","NACA18","CACA11","CACA17","CACA18","COMO11","COMO14","COMO18","MULE16","MULE19","MULE20","ROFO12","ROFO16","ROFO17","SJD11","SJD19","SJD20","BERRE12","BERRE15","BERRE19","BOCA14","BOCA19","BOCA20","BOR13","BOR14","BOR18","CATA1","CATA5","CATA9","STAMA12","STAMA19","STAMA20","JOSH14","JOSH18","JOSH19","KOFA12","KOFA13","KOFA18","PALM11","PALM17","PALM18","PALO13","PALO16","PALO20","ALA13","ALA15","ALA18")
#plot classic structure diagram
pdf(file="K_equals_8.pdf",10,5)
barplot(t(qmatrix), col = c(palette8), border = NA, space = 0,names.arg=names,las=2,cex.names=0.7,
        xlab = "Individuals", ylab = "Admixture coefficients")
dev.off()

