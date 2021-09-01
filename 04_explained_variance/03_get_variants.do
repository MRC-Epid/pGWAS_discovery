* do-file to obtain SNP information

* set the correct folder
cd /home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/03_explained_variance/data

* read the file with SNPS to be collected
import delimited query.snps.explained.variance.1.txt, varn(1)

* get the SNPs from the data set
snpdosage rsid, s(omics-fenland) imp(hrcuk10k) outfile(explained.omics.1)

* second data set
clear
import delimited query.snps.explained.variance.2.txt, varn(1)

* get the SNPs from the data set
snpdosage rsid, s(omics-fenland) imp(hrcuk10k) outfile(explained.omics.2)

* third data set
clear
import delimited query.snps.explained.variance.3.txt, varn(1)

* get the SNPs from the data set
snpdosage rsid, s(omics-fenland) imp(hrcuk10k) outfile(explained.omics.3)


