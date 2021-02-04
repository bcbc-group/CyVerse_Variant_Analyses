library(adegenet)
library(ade4)
library(MASS)

#import genepop file from Stacks with each individual
#had to remove the "epop" at the end of the file name to make work
gen <- import2genind(file="Washingtonia_test_samples.recode.p.snps.gen")


#import lat, long in XY format (long,lat)
xylist <- read.table(file="latlong.txt", header=TRUE)
rownames(xylist) <- xylist[,1]
xylist[,1] <- NULL
colnames(xylist) <- c("x","y")

#load lat and long into @other slot
gen@other$xy <- xylist

#view lat/long storeed in genind
xylist

#create necessary variable to run the analyses
toto <- genind2genpop(gen,process.other=TRUE)
toto

#calculate genetic and geographic distance
Dgen <- dist.genpop(toto, method=2)
Dgen
Dgeo <- dist(toto$other$xy,method="euclidean")
Dgeo
dim(as.matrix(Dgen))
dim(as.matrix(Dgeo))

#perform the isolation by distance analyses with 10,000 repetitions
ibd <- mantel.randtest(Dgen,Dgeo,nrepet=10000)
ibd

#plot the histogram of permutations to for significance
pdf(file="/home/rstudio/IBD_plot.pdf")
plot(ibd,main="Isolation by distance",xlab="Permuted simulated correlations")
dev.off()

#plot genetic distance against geographic distance
pdf(file="/home/rstudio/GeoDistance_vs_geneticDistance.pdf")
plot(Dgeo,Dgen, xlab="Geographic distance (km)", ylab="Genetic distance")
dev.off()

#color plot of genetic distance vs geographic distance
pdf(file="/home/rstudio/IBD_color_plot.pdf")
dens <- kde2d(Dgeo,Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen~Dgeo))
title("Isolation by distance plot")
dev.off()

#write a file with the genetic distance between individuals
mydf <- as.data.frame(as.matrix(Dgen))
write.csv(mydf, file="/home/rstudio/genetic_distance_matrix.csv")
