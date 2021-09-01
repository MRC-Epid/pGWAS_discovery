* do-file to obtain SNP information

* set the correct folder
cd /home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/15_GTEx/data

* read the file with SNPS to be collected
import delimited Look.up.snps.QTLs.txt, varn(1)

* get the SNPs from the data set
snpdosage rsid, s(omics-fenland) imp(hrcuk10k) outfile(snps.QTLs)

* second data set
* clear
* import delimited query.snps.2.txt, varn(1)

* get the SNPs from the data set
* snpdosage rsid, s(omics-fenland) imp(hrcuk10k) outfile(snps.omics.2)

* third data set
* clear
* import delimited query.snps.3.txt, varn(1)

* get the SNPs from the data set
* snpdosage rsid, s(omics-fenland) imp(hrcuk10k) outfile(snps.omics.3)

* fourth data set
* clear
* import delimited query.snps.4.txt, varn(1)

* get the SNPs from the data set
* snpdosage rsid, s(omics-fenland) imp(hrcuk10k) outfile(snps.omics.4)

* fith data set
* clear
* import delimited query.snps.5.txt, varn(1)

* get the SNPs from the data set
* snpdosage rsid, s(omics-fenland) imp(hrcuk10k) outfile(snps.omics.5)
