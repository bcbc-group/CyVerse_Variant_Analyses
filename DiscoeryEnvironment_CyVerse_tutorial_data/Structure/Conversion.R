struct2geno = function(file = NULL, output.format = "geno", TESS = FALSE, diploid = TRUE, FORMAT = 1, extra.row = 0, extra.col = 0, output = "genotype.geno"){
  
  ##### Description
  #
  # Export files from the TESS2.3 and STRUCTURE format to the ".geno" and ".lfmm" formats. 
  #
  #####
  
  
  ##### Arguments
  #
  # file = a file in the STRUCTURE or TESS format
  #
  # output.format = a character string. Allowed values are "geno" and "lfmm" 
  #
  # TESS = TRUE if the TESS 2.3 format is used (geographic coordinates (Lon,Lat) are binded to the left of the genotypic matrix).
  #
  # diploid = TRUE for diploids, = FALSE for haploids 
  #
  # FORMAT = 1 for markers encoded using one row of data for each individual
  #
  # FORMAT = 2 for markers encoded using two rows of data for each individual
  #
  # extra.row = an integer indicating the number of extra rows in the input file (e.g., marker ids)
  #
  # extra.col = an integer indicating the number of extra columns in the input file (e.g., individual ids, pop ids, phenotypes, etc)
  # geographic coordinates must be considered as extra columns if the flag TESS = FALSE   
  #
  # Missing data are encoded as "-9" or any negative values
  #
  #####
  
  ##### Value
  #
  # Output files in the ".geno" or ".lfmm" formats. Geographic coordinates in a separate file. 
  #
  #####
  
  if (!diploid & FORMAT == 2) stop("FORMAT = 2 is for diploids only")
  
  dat = read.table(file)
  
  
  if (TESS == FALSE){ 
    if (extra.row > 0) dat = dat[-(1:extra.row),]
    if (extra.col > 0) dat = dat[,-(1:extra.col)]
    n = dim(dat)[1]
    L = dim(dat)[2]
    if (FORMAT == 1 & diploid == FALSE) {n.ind = n; n.loc = L}
    if (FORMAT == 1 & diploid == TRUE) {n.ind = n; n.loc = L/2}
    if (FORMAT == 2 & diploid == TRUE) {n.ind = n/2; n.loc = L}    
    cat("Input file in the STRUCTURE format. The genotypic matrix has", n.ind, "individuals and", n.loc,
        "markers.","\n")
    cat("The number of extra rows is", extra.row,"and the number of extra columns is",extra.col,".\n")
  }
  
  
  if (TESS == TRUE){ 
    if (extra.row > 0) dat = dat[-(1:extra.row),]
    if (extra.col > 0) dat = dat[,-(1:extra.col)]
    
    n = dim(dat)[1]
    L = dim(dat)[2]
    
    if (FORMAT == 1 & diploid == FALSE) {
      n.ind = n; n.loc = L-2;
      coord = dat[,1:2]
      dat = dat[,-(1:2)]
      write.table(file = "coordinates.coord", coord, quote = F, row.names = F, col.names = F)
    }
    if (FORMAT == 1 & diploid == TRUE) {
      n.ind = n; n.loc = L/2 - 1; 
      coord = dat[,1:2]
      dat = dat[,-(1:2)]
      write.table(file = "coordinates.coord", coord, quote = F, row.names = F, col.names = F)      
    }
    if (FORMAT == 2 & diploid == TRUE) {
      n.ind = n/2; n.loc = L - 2;
      coord = dat[seq(1,n,by = 2), 1:2]
      dat = dat[,-(1:2)]
      write.table(file = "coordinates.coord", coord, quote = F, row.names = F, col.names = F)       
    }    
    cat("Input file in the TESS format. The genotypic matrix has", n.ind, "individuals and", n.loc,
        "markers.","\n")
    cat("The number of extra rows is", extra.row, 
        "and the number of extra columns is",extra.col,".\n")
  }
  
  dat = as.matrix(dat)
  
  unique.dat = unique(as.numeric(dat))
  missing.dat = unique.dat[unique.dat < 0]
  
  if (length(missing.dat) == 0)  cat("The input file contains no missing genotypes","\n")
  if (length(missing.dat) == 1)  cat("Missing alleles are encoded as",missing.dat,"\n")
  if (length(missing.dat) > 1) stop("Multiple values for missing data","\n") 
  
  
  
  # Convert allelic data into absence/presence data at each locus
  # Results are stored in the "dat.binary" object
  
  L = dim(dat)[2]
  
  if (FORMAT == 1 & diploid == FALSE) {
    
    dat.binary = NULL
    for (j in 1:L){
      allele = sort(unique(dat[,j]))
      for (i in allele[allele >= 0]) dat.binary=cbind(dat.binary, dat[,j]== i )
      LL = dim(dat.binary)[2]
      ind = which(dat[,j] < 0)
      if (length(ind) != 0){dat.binary[ind, (LL - length(allele) + 2):LL] = -9}
    }}
  
  if (FORMAT == 1 & diploid == TRUE) {
    dat.2 = matrix(NA, ncol = L/2, nrow = 2*n)
    for (ii in 1:n){
      dat.2[2*ii-1,] = dat[ii, seq(1,L,by = 2)]
      dat.2[2*ii,] = dat[ii,seq(2,L,by = 2)]
    }
    L = dim(dat.2)[2]
    
    dat.binary = NULL
    
    for (j in 1:L){
      allele = sort(unique(dat.2[,j]))
      for (i in allele[allele >= 0]) dat.binary=cbind(dat.binary, dat.2[,j]==i)
      LL = dim(dat.binary)[2]
      ind = which(dat.2[,j] < 0)
      if (length(ind) != 0){dat.binary[ind, (LL - length(allele) + 2):LL] = -9}
    }}
  
  
  
  if (FORMAT == 2 & diploid == TRUE) {
    dat.binary = NULL
    for (j in 1:L){
      allele = sort(unique(dat[,j]))
      for (i in allele[allele >= 0]) dat.binary=cbind(dat.binary, dat[,j]==i)
      LL = dim(dat.binary)[2]
      ind = which(dat[,j] < 0)
      if (length(ind) != 0){dat.binary[ind, (LL - length(allele) + 2):LL] = -9}
    }}
  
  # Compute a genotype count for each allele (0,1,2 or 9 for a missing value)
  # The results are stored in 'genotype'  
  
  n = dim(dat.binary)[1]
  
  if (diploid == TRUE){
    n = n/2
    genotype = matrix(NA,nrow=n,ncol=dim(dat.binary)[2])
    for(i in 1:n){
      genotype[i,]= dat.binary[2*i-1,]+dat.binary[2*i,]
      genotype[i, (genotype[i,] < 0)] = 9
    }}
  
  
  if (FORMAT == 1 & diploid == FALSE){
    genotype = dat.binary
    for(i in 1:n){
      genotype[i, (genotype[i,] < 0)] = 9
    }}  
  
  cat("The output file is",output, "\n") 
  
  if (output.format == "geno"){
    # Export to the "geno" format
    write.table(file = output,t(genotype), row.names=F,col.names=F,quote=F, sep = "")
  } else {
    # Export to the "lfmm" format
    write.table(file = output, genotype, row.names=F,col.names=F,quote=F)
  }
  
}

fst = function(project,run = 1, K, ploidy = 2){
  library(LEA)
  ll = dim(G(project, K = K, run = run))[1]

  if (ploidy == 2) {freq = G(project, K = K, run = run)[seq(2,ll,by = 3),]/2 + G(project, K = K, run = run)[seq(3,ll,by = 3),] } 
      else {freq = G(project, K = K, run = run)[seq(2,ll,by = 2),]}
   
   q = apply(Q(project, K = K, run = run), MARGIN = 2, mean)
   H.s = apply(freq*(1-freq), MARGIN = 1, FUN = function(x) sum(q*x) )
   P.t = apply(freq, MARGIN = 1, FUN = function(x) sum(q*x) )
   
   return(1-H.s/P.t/(1-P.t))
}

