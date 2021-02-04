###############################################################################
#EMMAX options can be found here https://genome.sph.umich.edu/wiki/EMMAX
###############################################################################

#make folders necessary later on
mkdir Phenotype
mkdir GWAS_output

plink --vcf Washingtonia_modified_GWAS_genotype_data.vcf --maf 0.05 --double-id --recode 12 --out plink_filter --output-missing-genotype 0 --transpose --allow-extra-chr

#Calculate kinship matrix
./emmax-kin-intel64 -v -d 10 plink_filter

#Make individual phenotype files from transformed phenotypes file containing all phenotypes
#make sure that Phenotype file does not have a header
#Very important, make sure that all sample names are in the same order as produced from the plink run
awk '{print $1,$1,$2}' Phenotype.txt > Phenotype/Trait1.txt
awk '{print $1,$1,$3}' Phenotype.txt > Phenotype/Trait2.txt
awk '{print $1,$1,$4}' Phenotype.txt > Phenotype/Trait3.txt
awk '{print $1,$1,$5}' Phenotype.txt > Phenotype/Trait4.txt


#run GWAS test
./emmax-intel64 -v -d 10 -t plink_filter -p Phenotype/Trait1.txt -k plink_filter.aBN.kinf -o GWAS_output/Trait1
./emmax-intel64 -v -d 10 -t plink_filter -p Phenotype/Trait2.txt -k plink_filter.aBN.kinf -o GWAS_output/Trait2
./emmax-intel64 -v -d 10 -t plink_filter -p Phenotype/Trait3.txt -k plink_filter.aBN.kinf -o GWAS_output/Trait3
./emmax-intel64 -v -d 10 -t plink_filter -p Phenotype/Trait4.txt -k plink_filter.aBN.kinf -o GWAS_output/Trait4

#create the position text file for plotting
awk '{print $1,$2}' plink_filter.map > position.txt
