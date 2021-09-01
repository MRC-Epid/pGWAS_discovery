###############################################
#### Conditional analysis SOMAscan pGWAS   ####
#### Maik Pietzner              01/06/2020 ####
###############################################

rm(list=ls())
setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/01_conditional_analysis/input/")
options(stringsAsFactors = F)
load(".RData")

#####################################
####    load sentinel variants   ####
#####################################

## regionel sentinels from Ellie (update: 06/08/2020), restrict to MAF > 1%
res.soma      <- read.table("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Ellie/SomaLogic/GWAS_discovery_subsets/Post_GWAS/01_regional_sentinel_variants/combined_regional_sentinels.out", sep="\t", header=T)

## exlude regions with MAF < 0.01! & MHC region
res.soma$MAF  <- ifelse(res.soma$Freq1 > .5, 1-res.soma$Freq1, res.soma$Freq1)

## omit HLA region and X-chromosome
write.table(res.soma[which(res.soma$MAF > .01 & !(res.soma$chr == 6 & res.soma$pos >= 28477797 & res.soma$pos <= 33448354) & res.soma$chr != 23), c("phenos", "chr", "region_start", "region_end")],
            "SOMAscan.regions.conditional.analysis.txt", sep="\t", col.names=F, row.names=F, quote=F)
## 7,660 regions taking forward

#####################################
####    check slurm scripts      ####
#####################################

ii <- dir("../")
ii <- grep("slurm", ii, value=T)

## get the slurm job matching
slurm.prot <- lapply(ii, function(x){
  ## open a connection
  con <- file(paste0("../", x), "r")
  fl  <- readLines(con)
  ## grep line with information require - misses header for some
  print(x)
  ## extract region
  ii <- grep("^--out output", fl)
  ll <- fl[ii]
  close(con)
  ## return slurm number and phenotype tested
  ll  <- strsplit(ll, "\\.")[[1]]
  return(data.frame(slurm=x, pheno=ll[2], chr=ll[3], start=ll[4], end=ll[5]))
})
slurm.prot     <- do.call(rbind, slurm.prot)
## generate number
slurm.prot$num <- as.numeric(gsub(".*-([[:digit:]]+).*", "\\1", slurm.prot$slurm))
## create simple identifier

#####################################
####   collect the output files  ####
#####################################

## --> SOMAscan <-- ##
soma.cond <- dir("../output//")
soma.cond <- grep("jma", soma.cond, value=T)

## read files and combine
soma.cond <- lapply(soma.cond, function(x){
  ## N.B. force character entries for some of the columns
  tmp       <- read.table(paste0("../output/", x), sep="\t", header=T, colClasses = c("numeric", "character", "numeric", "character", rep("numeric", 10)))
  ## order by association strength in primary data
  tmp       <- tmp[order(abs(tmp$b/tmp$se), decreasing = T),]
  ## add ordering
  tmp$hit   <- 1:nrow(tmp)
  ## add protein, position
  jj               <- strsplit(x, "\\.")[[1]][c(2,4,5)]
  tmp$pheno        <- jj[1]
  tmp$region_start <- jj[2]
  tmp$region_end   <- jj[3]
  return(tmp)
})
soma.cond <- do.call(rbind, soma.cond)

## create ID by using the lead SNP and phenotype
names(soma.cond) <- gsub("Chr", "chr", names(soma.cond))
names(res.soma)  <- gsub("phenos", "pheno", names(res.soma))
## add MA results to those from the conditional analysis
soma.cond        <- merge(soma.cond, unique(res.soma[, c("pheno", "chr", "region_start", "region_end", "MarkerName")]), by=c("pheno", "chr", "region_start", "region_end"))

## drop proteins with no secondary signal (first create identifier)
soma.cond$p_locus <- paste(soma.cond$MarkerName, soma.cond$pheno, sep="$")
## now count
tmp               <- aggregate(hit ~ p_locus, soma.cond, max)

## subset to those with at least one second signal
soma.cond         <- subset(soma.cond, p_locus %in% tmp$p_locus[which(tmp$hit > 1)])
length(unique(soma.cond$p_locus))


## add rsIDs to the SOMAscan results
require(data.table)
mapping   <- fread("~/rds/rds-rjh234-mrc-epid/Studies/People/Maik/Fenland.OMICs.MarkerName.rsID.chr.pos.txt", sep=" ", data.table=F)

## add to the data
soma.cond <- merge(soma.cond, mapping[, c("MarkerName", "rsid")], by.x="SNP", by.y="MarkerName", all.x=T)

## clean up
rm(mapping); gc()

#############################################################################
####              run final joint model in Fenland-OMICs                 ####
#############################################################################

## export list of variants to be loaded from Fenland-OMICs, split
tmp <- data.frame(rsid=unique(soma.cond$rsid))
write.table(tmp[1:1800,,drop=F], "Conditional.Variants.GCTA.20200807.1.txt", sep="\t", row.names=F)
write.table(tmp[1801:nrow(tmp),,drop=F], "Conditional.Variants.GCTA.20200807.2.txt", sep="\t", row.names=F)

## load phenotype data (residuals, similar as in the GWAS)
require(data.table)
pheno      <- fread("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Ellie/SomaLogic/GWAS_discovery_subsets/Post_GWAS/conditional_analysis/Fenland-OMICS_res_invnX_8350.sample", sep="\t", header=T, data.table=F)
## skip first line and replace missing values
pheno      <- pheno[-1,]
## restrict to what is needed
pheno      <- pheno[, c("ID_1", unique(soma.cond$pheno))] ## 1,137 somamers
## make numeric
pheno[,-1] <- apply(pheno[,-1], 2, function(x){
  x            <- as.numeric(x)
  x[x == -999] <- NA
  return(x)
})

## load the SNP data
require(readstata13)
snps1       <- read.dta13("conditional.GCTA.omics.1.dta")
snps2       <- read.dta13("conditional.GCTA.omics.2.dta")
## combine both
snps        <- merge(snps1, snps2)
rm(snps1); rm(snps2); gc()

## load label
snps.label1 <- read.table("conditional.GCTA.omics.1.dosage_code", sep="\t", header=T)
snps.label2 <- read.table("conditional.GCTA.omics.2.dosage_code", sep="\t", header=T)
snps.label  <- rbind(snps.label1, snps.label2)

## create markername column
snps.label$markername <- apply(snps.label[, c(1,4:6)], 1, function(x){
  paste0("chr", as.numeric(x[1]), ":", as.numeric(x[2]), "_", paste(sort(x[3:4]), collapse = "_"))
})
rm(snps.label1); rm(snps.label2); gc()

## merge the data
pheno       <- merge(pheno, snps, by.x="ID_1", by.y="id")

#--------------------------------------#
##--      compute joint models      --##
#--------------------------------------#

res.joint <- lapply(unique(soma.cond$p_locus), function(x){
  ## obtain all SNPs needed for the analysis
  snps <- subset(soma.cond, p_locus == x)$SNP
  ## map to rsIDs
  snps <- sapply(snps, function(k) snps.label$rsid2[which(snps.label$markername == k)], simplify = T)
  ## some SNPs miss rsIDs, omit those
  if(mode(snps) == "character"){
    ## outcome to be used
    outc <- strsplit(x, "\\$")[[1]][2]
    ## create the formula
    ff   <- paste(outc, "~", paste(snps, collapse = " + "))## run the model
    ff   <- summary(lm(as.formula(ff), data = pheno))$coefficients
    ## prepare for storage (omit Intercept)
    ff          <- data.frame(ff[-1,])
    ff$sentinel <- strsplit(x, "\\$")[[1]][1]
    ff$rsid2    <- rownames(ff)
    ff$pheno    <- outc
    ff$p_locus  <- x
    return(ff)
  }else{
    cat("missing snps for", x, "\n")
  }
  
})
res.joint        <- do.call(rbind, res.joint)
save.image()

## edit names
names(res.joint) <- c("b", "se", "tval", "p", "sentinel", "rsid2", "pheno", "p_locus")
## add SNP identifiers
res.joint        <- merge(res.joint, snps.label, by="rsid2")
## rename some 
names(res.joint) <- gsub("markername", "SNP", names(res.joint))

## create combined data set
res.all          <- merge(soma.cond, res.joint, by=c("p_locus", "SNP", "pheno"), suffixes = c(".ma", ".joint.omics"))

## order and look what didn't remained significant
res.all          <- res.all[order(res.all$p_locus, res.all$hit),]

## indicate loci with difficult observations
res.all$cond_lost <- ifelse(res.all$p_locus %in% unique(subset(res.all, p.joint.omics > 5e-8)$p_locus), 1, 0) 

## skip all variants not passing the threshold
res.all           <- subset(res.all, p.joint.omics < 5e-8)
## reassign hits
res.all           <- res.all[order(res.all$p_locus, -abs(res.all$b.ma/res.all$se.ma)),]
res.all$hit       <- unlist(lapply(unique(res.all$p_locus), function(x) seq(1, nrow(subset(res.all, p_locus == x)))))

## summarize and drop those with no remaining hits
tmp               <- aggregate(hit ~ p_locus, data = res.all, max)
summary(tmp)
## susbet 
res.all          <- subset(res.all, p_locus %in% subset(tmp, hit > 1)$p_locus)
## summarize 
tmp              <- aggregate(hit ~ p_locus, data = res.all, max)
summary(tmp)
## 1285 loci with secondary signals

#############################################################################
####                           write to file                             ####
#############################################################################

write.table(res.all, "Results.conditional.analysis.pGWAS.GCTA.QCed.20200807.txt", sep="\t", row.names=F)

