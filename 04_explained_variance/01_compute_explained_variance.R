############################################################
#### Exlpained variance for pQTLs from SOMAscan         ####
#### Maik Pietzner                           17/07/2020 ####
############################################################

rm(list=ls())
setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/03_explained_variance/data/")
options(stringsAsFactors = F)
load(".RData")

############################################
####  import results from both studies  ####
############################################

## load updated signal list
res.gwas  <- read.table("../../Supplemental.table.1.20210225.txt", sep="\t", header=T)

############################################
####   create a list to query all SNPs  ####
############################################

## create file
query <- data.frame(rsid=unique(res.gwas$rsID))
## 5,442 SNPs to separate files

## export to do in STATA
write.table(query[1:2000, , drop=F], "query.snps.explained.variance.1.txt", sep="\t", row.names=F)
write.table(query[2001:4000, , drop=F], "query.snps.explained.variance.2.txt", sep="\t", row.names=F)
write.table(query[4001:nrow(query), , drop=F], "query.snps.explained.variance.3.txt", sep="\t", row.names=F)

############################################
####     import the phenotype data      ####
############################################

## load phenotype data (residuals, similar as in the GWAS)
require(data.table)
pheno      <- fread("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Ellie/SomaLogic/GWAS_discovery_subsets/Post_GWAS/conditional_analysis/Fenland-OMICS_res_invnX_8350.sample", sep="\t", header=T, data.table=F)
## skip first line and replace missing values
pheno      <- pheno[-1,]
## only those needed
pheno      <- pheno[, c("ID_1", unique(res.gwas$GWAS.id))]
## make numeric
pheno[,-1] <- apply(pheno[,-1], 2, function(x){
  x            <- as.numeric(x)
  x[x == -999] <- NA
  return(x)
})

## load the SNP data
require(readstata13)
snps1       <- read.dta13("explained.omics.1.dta")
snps2       <- read.dta13("explained.omics.2.dta")
snps3       <- read.dta13("explained.omics.3.dta")
## combine both
snps        <- merge(snps1, snps2)
snps        <- merge(snps, snps3)
rm(snps1); rm(snps2); rm(snps3); gc()

## load label
snps.label1 <- read.table("explained.omics.1.dosage_code", sep="\t", header=T)
snps.label2 <- read.table("explained.omics.2.dosage_code", sep="\t", header=T)
snps.label3 <- read.table("explained.omics.3.dosage_code", sep="\t", header=T)
snps.label  <- rbind(snps.label1, snps.label2, snps.label3)
## create markername column
snps.label$markername <- apply(snps.label[, c(1,4:6)], 1, function(x){
  paste0("chr", as.numeric(x[1]), ":", as.numeric(x[2]), "_", paste(sort(x[3:4]), collapse = "_"))
})

## test whether all SNPs have been included
snps.label  <- subset(snps.label, rsid1 %in% query$rsid)

## merge the data
pheno       <- merge(pheno, snps, by.x="ID_1", by.y="id")
save.image()

############################################
####     compute explained variance     ####
############################################

## all proteins with at least one significant hit
soma.id  <- unique(res.gwas$GWAS.id) ## 4,030 SOMAmers

## add column with SNP IDs to the 
res.gwas <- merge(res.gwas, snps.label, by.x="rsID", by.y="rsid1")

## run analysis for each aptamers
require(doMC)
registerDoMC(10)

res.var  <- mclapply(soma.id, function(x){
  
  ## get all relevant SNPs
  ii <- subset(res.gwas, GWAS.id == x)
  
  ## lead cis-pQTL
  if(nrow(subset(ii, cis.trans == "cis")) > 0){
    
    ## get the lead signal
    kk         <- subset(ii, cis.trans == "cis" & signal.locus == 1)$rsid2
    e.lead.cis <- summary(lm(as.formula(paste(x, kk, sep=" ~ ")), data = pheno))$r.squared
    n.lead.cis <- length(kk)
    
  }else{
    e.lead.cis <- 0
    n.lead.cis <- 0
  }
  
  ## secondary cis-pQTL
  if(nrow(subset(ii, cis.trans == "cis" & signal.locus > 1)) > 0){
    
    ## get the lead signal
    kk        <- subset(ii, cis.trans == "cis" & signal.locus > 1)$rsid2
    e.sec.cis <- summary(lm(as.formula(paste(x, paste(kk, collapse = " + "), sep=" ~ ")), data = pheno))$r.squared
    n.sec.cis <- length(kk)
    
  }else{
    e.sec.cis <- 0
    n.sec.cis <- 0
  }
  
  ## all cis-pQTL
  if(nrow(subset(ii, cis.trans == "cis")) > 0){
    
    ## get the lead signal
    kk        <- subset(ii, cis.trans == "cis")$rsid2
    e.all.cis <- summary(lm(as.formula(paste(x, paste(kk, collapse = " + "), sep=" ~ ")), data = pheno))$r.squared
    n.all.cis <- length(kk)
    
  }else{
    e.all.cis <- 0
    n.all.cis <- 0
  }
  
  ## all trans-pQTLs
  if(nrow(subset(ii, cis.trans != "cis")) > 0){
    
    ## get the lead signal
    kk          <- subset(ii, cis.trans != "cis")$rsid2
    e.all.trans <- summary(lm(as.formula(paste(x, paste(kk, collapse = " + "), sep=" ~ ")), data = pheno))$r.squared
    n.all.trans <- length(kk)
    
  }else{
    e.all.trans <- 0
    n.all.trans <- 0
  }
  
  ## specific trans-pQTLs
  if(nrow(subset(ii, cis.trans != "cis" & !(pQTL.classification %in% c("none", "unspecific_effect")))) > 0){
    
    ## get the lead signal
    kk           <- subset(ii, cis.trans != "cis" & !(pQTL.classification %in% c("none", "unspecific_effect")))$rsid2
    e.spec.trans <- summary(lm(as.formula(paste(x, paste(kk, collapse = " + "), sep=" ~ ")), data = pheno))$r.squared
    n.spec.trans <- length(kk)
    
  }else{
    e.spec.trans <- 0
    n.spec.trans <- 0
  }
  
  ## unspecific trans-pQTLs
  if(nrow(subset(ii, cis.trans != "cis" & pQTL.classification %in% c("none", "unspecific_effect"))) > 0){
    
    ## get the lead signal
    kk             <- subset(ii, cis.trans != "cis" & pQTL.classification %in% c("none", "unspecific_effect"))$rsid2
    e.unspec.trans <- summary(lm(as.formula(paste(x, paste(kk, collapse = " + "), sep=" ~ ")), data = pheno))$r.squared
    n.unspec.trans <- length(kk)
    
  }else{
    e.unspec.trans <- 0
    n.unspec.trans <- 0
  }
  
  ## all pQTLs
  if(nrow(ii) > 0){
    
    ## get the lead signal
    kk    <- ii$rsid2
    e.all <- summary(lm(as.formula(paste(x, paste(kk, collapse = " + "), sep=" ~ ")), data = pheno))$r.squared
    n.all <- nrow(ii)
    
  }else{
    e.all <- 0
    n.all <- 0
  }
  
  ## create combined data frame
  return(data.frame(pheno = x, n.lead.cis = n.lead.cis, e.lead.cis = e.lead.cis,
                    n.sec.cis = n.sec.cis, e.sec.cis = e.sec.cis,
                    n.all.cis = n.all.cis, e.all.cis = e.all.cis,
                    n.all.trans = n.all.trans, e.all.trans = e.all.trans,
                    n.spec.trans = n.spec.trans, e.spec.trans = e.spec.trans,
                    n.unspec.trans = n.unspec.trans, e.unspec.trans = e.unspec.trans,
                    n.all = n.all, e.all = e.all))
  
}, mc.cores = 10)
res.var <- do.call(rbind, res.var)

## write to file
write.table(res.var, "Explained.variance.SOMAscan.20210701.txt", sep="\t", row.names=F)

############################################################################################
####                                    END OF SCRIPT                                   ####
############################################################################################

