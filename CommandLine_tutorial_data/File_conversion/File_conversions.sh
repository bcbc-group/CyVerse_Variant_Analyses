############################
######     Stacks    #######
############################

mkdir data_set_individuals
mkdir data_set_SNAPP
mkdir data_set_populations


#after vcftools filters individuals
populations --batch_size 1 -V Washingtonia_test_samples.recode.vcf -O data_set_individuals -M pop_map_individuals.txt -t 2 --ordered-export --vcf --genepop 

#after vcftools filters populations for Structure
populations --batch_size 1 -V Washingtonia_test_samples.recode.vcf -O data_set_populations -M pop_map_populations.txt -t 2 --ordered-export --vcf --structure 

############################
######  VCF2Phylip   #######
############################

#convert to nexus and fasta format, for SplitsTree and SNAPP respectively
python3 vcf2phylip.py -i Washingtonia_test_samples.recode.p.snps.vcf -n -f


python3 vcf2phylip.py -i Washingtonia_test_SNAPP.recode.vcf -n -f
