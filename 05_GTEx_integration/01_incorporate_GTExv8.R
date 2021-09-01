##############################################
#### Incorporate eQTL data from Erin      ####
#### Maik Pietzner             04/11/2020 ####
##############################################

rm(list=ls())
setwd("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/15_GTEx/data/")
options(stringsAsFactors = F)
load(".RData")

##############################################
####    load results tables for pQTLs     ####
##############################################

res.soma  <- read.table("../../06_locus_definition_novelty/data/Regional.sentinels.pGWAS.SomaLogic.annotated.locus.novelty.20200813.txt", sep="\t", header=T)
cond.soma <- read.table("../../06_locus_definition_novelty/data/Conditional.results.pGWAS.SomaLogic.annotated.novelty.20200813.txt", sep="\t", header=T)
## rename to ease merging
names(cond.soma)[11:14] <- c("Freq1", "Effect", "StdErr", "P.value")

## generat common data set to map to
ii            <- c("p_locus", "cis_trans", "rsid", "MarkerName", "chr", "pos", "Allele1", "Allele2", "Freq1","pheno", "Target", "TargetFullName", "Effect", "StdErr", "P.value", "published.pqtls.1",  
                   "Ensemblgeneid", "hgnc_symbol", "consequence.lead", "gene.lead", "symbol.lead", "rsid.ld.proxy", "r2.ld.proxy", "gene.ld.proxy", "symbol.ld.proxy",  "consequence.ld.proxy")
res.pqtls     <- res.soma[, ii]
res.pqtls$hit <- 1
res.pqtls     <- rbind(res.pqtls, cond.soma[, c(ii, "hit")])
res.pqtls     <- res.pqtls[order(res.pqtls$p_locus, res.pqtls$hit),]

##############################################
####      import list from look-up        ####
##############################################

## define header
hd      <- c("rsid", "MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "Pvalue", "Direction", "HetSq", "HetChiSq",
             "HetDf", "HetPVal", "TotalSampleSize", "chr", "pos")

## read in the look_up data
ii      <- dir("../look_up/")
look.up <- lapply(ii, function(x){
  ## read the file (some SNPs have been filtered)
  tmp        <- read.table(paste0("../look_up/", x), sep="\t", header=F)
  ## assign header
  names(tmp) <- hd
  ## add phenotype
  tmp$pheno  <- strsplit(x, "\\.")[[1]][1]
  return(tmp)
}) 
look.up <- do.call(rbind, look.up)

## write to file
write.table(look.up, "Summary.stats.proxy.pQTLs.eQTLs.txt", sep="\t", row.names=F)


##############################################
####       including sQTL look-up         ####
##############################################

## read in the look_up data
ii      <- dir("../look_up/")
look.up <- lapply(ii, function(x){
  ## read the file (some SNPs have been filtered)
  tmp        <- read.table(paste0("../look_up/", x), sep="\t", header=F)
  ## assign header
  names(tmp) <- hd
  ## add phenotype
  tmp$pheno  <- strsplit(x, "\\.")[[1]][1]
  return(tmp)
}) 
look.up <- do.call(rbind, look.up)

## write to file
write.table(look.up, "Summary.stats.proxy.pQTLs.eQTLs.sQTLs.txt", sep="\t", row.names=F)

##############################################
####         prepare eQTL coloc           ####
##############################################

## import list
eqtl.coloc <- read.table("eQTL.pQTL.coloc.list.lead.txt", sep="\t", header=T)

## add positions (250kb flanking region)
eqtl.coloc$region_start <- eqtl.coloc$pos.pGWAS - 25e4
eqtl.coloc$region_end   <- eqtl.coloc$pos.pGWAS + 25e4
## add negative borders
eqtl.coloc$region_start <- ifelse(eqtl.coloc$region_start < 0, 0, eqtl.coloc$region_start)

## create list for look up
write.table(eqtl.coloc[, c("pheno", "gene_id", "chr.pGWAS", "region_start", "region_end")], "eQTL.coloc.list.txt", sep="\t", col.names = F, row.names = F, quote=F)

#---------------------------------#
##--       import results      --##
#---------------------------------#

ii <- dir("../tissue_coloc/")

## create indicator, which pairs are missing
eqtl.coloc$missing <- apply(eqtl.coloc[, c("pheno", "gene_id", "chr.pGWAS", "region_start", "region_end")], 1, function(x){
  ## combine to match file names
  x <- paste(x[1], x[2], as.numeric(x[3]), as.numeric(x[4]), as.numeric(x[5]), "txt", sep=".")
  # print(x)
  ## test whether in the file list
  return(x %in% ii)
})
table(eqtl.coloc$missing)
## missing is due to differential X chromosome naming, that is not really missing at all

## rerun failed jobs
# write.table(eqtl.coloc[which(eqtl.coloc$missing == F), c("pheno", "gene_id", "chr.pGWAS", "region_start", "region_end")], 
#             "eQTL.coloc.list.rerun.txt", sep="\t", col.names = F, row.names = F, quote=F)

## import all results
res.eqtl.coloc <- lapply(ii, function(x){
  ## read the results
  tmp               <- read.table(paste0("../tissue_coloc/", x), header=T, sep="\t")
  ## add phenotype
  tmp$pheno        <- strsplit(x, "\\.")[[1]][1]
  tmp$region_start <- strsplit(x, "\\.")[[1]][4]
  tmp$region_end   <- strsplit(x, "\\.")[[1]][5]
  tmp$gene         <- strsplit(x, "\\.")[[1]][2]
  return(tmp)
})
res.eqtl.coloc <- do.call(rbind, res.eqtl.coloc)

#---------------------------------#
##--      process the data     --##
#---------------------------------#

## create ID to make unique
res.eqtl.coloc$id <- paste(res.eqtl.coloc$pheno, res.eqtl.coloc$gene, sep="$")

## package to perform meta-analysis
require(metafor)

## aggregate the results 
res.eqtl.coloc.top <- lapply(unique(res.eqtl.coloc$id), function(x){
  ## get the coloc results
  tmp <- subset(res.eqtl.coloc, id == x)
  ## get the top coloc signal
  ii  <- which.max(tmp$PP.H4.abf)
  ## get a heterogeneity estimate across tissues
  ht  <- rma(yi=tmp$slope, sei=tmp$slope_se, method="REML")
  ## how many with evidence from coloc
  nc  <- nrow(subset(tmp, PP.H4.abf > .7))
  ## how many with evidence from LD
  nd  <- nrow(subset(tmp, ld.check.top > .8))
  ## return collated information
  return(data.frame(tmp[ii,], n.coloc=nc, n.ld=nd, pval.hetero=ht$QEp, iSq=ht$I2))
})
res.eqtl.coloc.top <- do.call(rbind, res.eqtl.coloc.top)

## flag inconsistent results
res.eqtl.coloc.top$flag <- ifelse(res.eqtl.coloc.top$ld.check.sens < .8, 1, 0)

## write to file
write.table(res.eqtl.coloc.top, "Results.cis.pQTL.eQTL.coloc.top.20210121.txt", sep="\t", row.names=F)

## write more comprehensive results as well for clustering
write.table(res.eqtl.coloc, "Results.cis.pQTL.eQTL.coloc.all.20210121.txt", sep="\t", row.names=F)

##############################################
####       prepare cond. eQTL coloc       ####
##############################################

## import list
eqtl.coloc.cond <- read.table("eQTL.pQTL.coloc.list.secondary.txt", sep="\t", header=T)

## add positions (250kb flanking region)
eqtl.coloc.cond$region_start <- eqtl.coloc.cond$pos.pGWAS - 25e4
eqtl.coloc.cond$region_end   <- eqtl.coloc.cond$pos.pGWAS + 25e4
## add negative borders
eqtl.coloc.cond$region_start <- ifelse(eqtl.coloc.cond$region_start < 0, 0, eqtl.coloc.cond$region_start)

## create list for look up
write.table(eqtl.coloc.cond[, c("pheno", "gene_id", "chr.pGWAS", "region_start", "region_end")], 
            "eQTL.coloc.cond.list.txt", sep="\t", col.names = F, row.names = F, quote=F)

## export lead signals at each locus to create mathching conditional stats
for(j in eqtl.coloc.cond$pheno){
  ## identify the matching lead cis-signal
  ii <- subset(res.soma, pheno == j & cis_trans == "cis")
  ## write to file to be used for adjustment
  write.table(ii$MarkerName, paste0("../input_cond/cond.", j, ".txt"), col.names = F, row.names = F, quote = F, sep="\t")
}

#---------------------------------#
##--       import results      --##
#---------------------------------#

ii <- dir("../tissue_coloc_cond//")

## import all results
res.eqtl.coloc.cond <- lapply(ii, function(x){
  ## read the results
  tmp               <- read.table(paste0("../tissue_coloc_cond//", x), header=T, sep="\t")
  ## add phenotype
  tmp$pheno        <- strsplit(x, "\\.")[[1]][1]
  tmp$region_start <- strsplit(x, "\\.")[[1]][4]
  tmp$region_end   <- strsplit(x, "\\.")[[1]][5]
  tmp$gene         <- strsplit(x, "\\.")[[1]][2]
  return(tmp)
})
res.eqtl.coloc.cond <- do.call(rbind, res.eqtl.coloc.cond)

#---------------------------------#
##--      process the data     --##
#---------------------------------#

## create ID to make unique
res.eqtl.coloc.cond$id <- paste(res.eqtl.coloc.cond$pheno, res.eqtl.coloc.cond$gene, sep="$")

## package to perform meta-analysis
require(metafor)

## aggregate the results 
res.eqtl.coloc.cond.top <- lapply(unique(res.eqtl.coloc.cond$id), function(x){
  ## get the coloc results
  tmp <- subset(res.eqtl.coloc.cond, id == x)
  ## get the top coloc signal
  ii  <- which.max(tmp$PP.H4.abf)
  ## get a heterogeneity estimate across tissues
  ht  <- rma(yi=tmp$slope, sei=tmp$slope_se, method="REML")
  ## how many with evidence from coloc
  nc  <- nrow(subset(tmp, PP.H4.abf > .7))
  ## how many with evidence from LD
  nd  <- nrow(subset(tmp, ld.check.top > .8))
  ## return collated information
  return(data.frame(tmp[ii,], n.coloc=nc, n.ld=nd, pval.hetero=ht$QEp, iSq=ht$I2))
})
res.eqtl.coloc.cond.top <- do.call(rbind, res.eqtl.coloc.cond.top)

## flag inconsistent results
res.eqtl.coloc.cond.top$flag <- ifelse(res.eqtl.coloc.cond.top$ld.check.sens < .8, 1, 0)

## write to file
write.table(res.eqtl.coloc.cond.top, "Results.cis.pQTL.eQTL.coloc.cond.top.20210122.txt", sep="\t", row.names=F)

## write more comprehensive results as well for clustering
write.table(res.eqtl.coloc.cond, "Results.cis.pQTL.eQTL.coloc.cond.all.20210122.txt", sep="\t", row.names=F)

##############################################
####         prepare sQTL coloc           ####
##############################################

## import list
sqtl.coloc        <- read.table("sQTL.pQTL.coloc.list.lead.txt", sep="\t", header=T)

## make unique for introns
sqtl.coloc$intron <- sapply(sqtl.coloc$phenotype_id, function(x){
 ## get coordinates of the intron
 x <- strsplit(x, ":")[[1]][1:3]
 ## re-assemble
 return(paste(x[1], as.numeric(x[2]), as.numeric(x[3]), sep=":"))
})

## make unique
sqtl.coloc              <- unique(sqtl.coloc[,-which(names(sqtl.coloc) == "phenotype_id")]) 

## add positions (250kb flanking region)
sqtl.coloc$region_start <- sqtl.coloc$pos.pGWAS - 25e4
sqtl.coloc$region_end   <- sqtl.coloc$pos.pGWAS + 25e4
## add negative borders
sqtl.coloc$region_start <- ifelse(sqtl.coloc$region_start < 0, 0, sqtl.coloc$region_start)

## create list for look up
write.table(sqtl.coloc[, c("pheno", "intron", "gene_id", "chr.pGWAS", "region_start", "region_end")], "sQTL.coloc.list.txt", sep="\t", col.names = F, row.names = F, quote=F)

#---------------------------------#
##--       import results      --##
#---------------------------------#

ii <- dir("../tissue_coloc_sQTL//")

## import all results
res.sqtl.coloc <- lapply(ii, function(x){
  ## read the results
  tmp               <- read.table(paste0("../tissue_coloc_sQTL//", x), header=T, sep="\t")
  ## add phenotype
  tmp$pheno        <- strsplit(x, "\\.")[[1]][1]
  tmp$region_start <- strsplit(x, "\\.")[[1]][5]
  tmp$region_end   <- strsplit(x, "\\.")[[1]][6]
  tmp$gene         <- strsplit(x, "\\.")[[1]][3]
  tmp$intron       <- strsplit(x, "\\.")[[1]][2]
  return(tmp)
})
res.sqtl.coloc <- do.call(rbind, res.sqtl.coloc)

#---------------------------------#
##--      process the data     --##
#---------------------------------#

## create ID to make unique
res.sqtl.coloc$id <- paste(res.sqtl.coloc$pheno, res.sqtl.coloc$intron, res.sqtl.coloc$gene, sep="$")

## package to perform meta-analysis
require(metafor)

## aggregate the results 
res.sqtl.coloc.top <- lapply(unique(res.sqtl.coloc$id), function(x){
  ## get the coloc results
  tmp <- subset(res.sqtl.coloc, id == x)
  ## get the top coloc signal
  ii  <- which.max(tmp$PP.H4.abf)
  ## get a heterogeneity estimate across tissues
  ht  <- rma(yi=tmp$slope, sei=tmp$slope_se, method="REML")
  ## how many with evidence from coloc
  nc  <- nrow(subset(tmp, PP.H4.abf > .7))
  ## how many with evidence from LD
  nd  <- nrow(subset(tmp, ld.check.top > .8))
  ## return collated information
  return(data.frame(tmp[ii,], n.coloc=nc, n.ld=nd, pval.hetero=ht$QEp, iSq=ht$I2))
})
res.sqtl.coloc.top <- do.call(rbind, res.sqtl.coloc.top)

## flag inconsistent results
res.sqtl.coloc.top$flag <- ifelse(res.sqtl.coloc.top$ld.check.sens < .8, 1, 0)
table(res.sqtl.coloc.top$flag) ## 19 flagged

## how many meet the threshold
nrow(subset(res.sqtl.coloc.top, PP.H4.abf > .8))

## write to file
write.table(res.sqtl.coloc.top, "Results.cis.pQTL.sQTL.coloc.top.20210128.txt", sep="\t", row.names=F)

## write more comprehensive results as well for clustering
write.table(res.sqtl.coloc, "Results.cis.pQTL.sQTL.coloc.all.20210128.txt", sep="\t", row.names=F)

##############################################
####       prepare cond. eQTL coloc       ####
##############################################

## import list
sqtl.coloc.cond        <- read.table("sQTL.pQTL.coloc.list.secondary.txt", sep="\t", header=T)

## make unique for introns
sqtl.coloc.cond$intron <- sapply(sqtl.coloc.cond$phenotype_id, function(x){
  ## get coordinates of the intron
  x <- strsplit(x, ":")[[1]][1:3]
  ## re-assemble
  return(paste(x[1], as.numeric(x[2]), as.numeric(x[3]), sep=":"))
})

## make unique
sqtl.coloc.cond        <- unique(sqtl.coloc.cond[,-which(names(sqtl.coloc.cond) == "phenotype_id")]) 

## add positions (250kb flanking region)
sqtl.coloc.cond$region_start <- sqtl.coloc.cond$pos.pGWAS - 25e4
sqtl.coloc.cond$region_end   <- sqtl.coloc.cond$pos.pGWAS + 25e4
## add negative borders
sqtl.coloc.cond$region_start <- ifelse(sqtl.coloc.cond$region_start < 0, 0, sqtl.coloc.cond$region_start)

## create list for look up
write.table(sqtl.coloc.cond[, c("pheno", "intron", "gene_id", "chr.pGWAS", "region_start", "region_end")], "sQTL.coloc.cond.list.txt", sep="\t", col.names = F, row.names = F, quote=F)

## export lead signals at each locus to create mathching conditional stats
for(j in sqtl.coloc.cond$pheno){
  ## identify the matching lead cis-signal
  ii <- subset(res.soma, pheno == j & cis_trans == "cis")
  ## write to file to be used for adjustment
  write.table(ii$MarkerName, paste0("../input_cond/cond.", j, ".txt"), col.names = F, row.names = F, quote = F, sep="\t")
}

#---------------------------------#
##--       import results      --##
#---------------------------------#

ii <- dir("../tissue_coloc_sQTL_cond///")

## import all results
res.sqtl.coloc.cond <- lapply(ii, function(x){
  ## read the results
  tmp               <- read.table(paste0("../tissue_coloc_sQTL_cond///", x), header=T, sep="\t")
  ## add phenotype
  tmp$pheno        <- strsplit(x, "\\.")[[1]][1]
  tmp$region_start <- strsplit(x, "\\.")[[1]][5]
  tmp$region_end   <- strsplit(x, "\\.")[[1]][6]
  tmp$gene         <- strsplit(x, "\\.")[[1]][3]
  tmp$intron       <- strsplit(x, "\\.")[[1]][2]
  return(tmp)
})
res.sqtl.coloc.cond <- do.call(rbind, res.sqtl.coloc.cond)

#---------------------------------#
##--      process the data     --##
#---------------------------------#

## create ID to make unique
res.sqtl.coloc.cond$id <- paste(res.sqtl.coloc.cond$pheno, res.sqtl.coloc.cond$intron, res.sqtl.coloc.cond$gene, sep="$")

## package to perform meta-analysis
require(metafor)

## aggregate the results 
res.sqtl.coloc.cond.top <- lapply(unique(res.sqtl.coloc.cond$id), function(x){
  ## get the coloc results
  tmp <- subset(res.sqtl.coloc.cond, id == x)
  ## get the top coloc signal
  ii  <- which.max(tmp$PP.H4.abf)
  ## get a heterogeneity estimate across tissues
  ht  <- rma(yi=tmp$slope, sei=tmp$slope_se, method="REML")
  ## how many with evidence from coloc
  nc  <- nrow(subset(tmp, PP.H4.abf > .7))
  ## how many with evidence from LD
  nd  <- nrow(subset(tmp, ld.check.top > .8))
  ## return collated information
  return(data.frame(tmp[ii,], n.coloc=nc, n.ld=nd, pval.hetero=ht$QEp, iSq=ht$I2))
})
res.sqtl.coloc.cond.top <- do.call(rbind, res.sqtl.coloc.cond.top)

## flag inconsistent results
res.sqtl.coloc.cond.top$flag <- ifelse(res.sqtl.coloc.cond.top$ld.check.sens < .8, 1, 0)

## write to file
write.table(res.sqtl.coloc.cond.top, "Results.cis.pQTL.sQTL.coloc.cond.top.20210201.txt", sep="\t", row.names=F)

## write more comprehensive results as well for clustering
write.table(res.sqtl.coloc.cond, "Results.cis.pQTL.sQTL.coloc.cond.all.20210201.txt", sep="\t", row.names=F)