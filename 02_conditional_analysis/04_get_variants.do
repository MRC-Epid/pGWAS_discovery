* do-file to obtain SNP information

* set the correct folder
cd /home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/01_conditional_analysis/input

* read the file with SNPS to be collected
import delimited Conditional.Variants.GCTA.20200807.1.txt, varn(1)

* get the SNPs from the data set
snpdosage rsid, s(omics-fenland) imp(hrcuk10k) outfile(conditional.GCTA.omics.1)

* clear and load next data
clear

* read the file with SNPS to be collected
import delimited Conditional.Variants.GCTA.20200807.2.txt, varn(1)

* get the SNPs from the data set
snpdosage rsid, s(omics-fenland) imp(hrcuk10k) outfile(conditional.GCTA.omics.2)


