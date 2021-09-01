#!/usr/bin/env Rscript

## script to run cond. coloc for mapping protein gene pairs using GTEx sQTL data
## Maik Pietzner 28/01/2021
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/15_GTEx/")

## --> import parameters <-- ##

## aptamer
soma  <- args[1]
## intron
intr  <- args[2]
## gene
gene  <- args[3]
## chromosome
chr.s <- as.numeric(args[4])
## start position of the region
pos.s <- as.numeric(args[5])
## end position of the region
pos.e <- as.numeric(args[6])

#-----------------------------------------#
##-- 	       import protein data       --##
#-----------------------------------------#

## read the relevant data
res.soma         <- data.table::fread(paste0("input_cond/", soma, ".cma.cojo"), sep = "\t", header = T, data.table = F)
## do some renaming to align with the results (rename bC to actual beta!)
names(res.soma)  <- c("chr", "MarkerName", "pos", "Allele1", "Freq1", "beta_ma", "se_ma", "pval_ma", "TotalSampleSize", "freq_geno", "Effect", "StdErr", "Pvalue")
## add rsID
map              <- data.table::fread("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/Fenland.OMICs.MarkerName.rsID.chr.pos.txt", sep=" ", header=T, data.table=F)
## add to the data
res.soma         <- merge(res.soma, map)
## omit missing values
res.soma         <- subset(res.soma, !is.na(Effect))
rm(map); gc()

## apply MAF filter
res.soma$MAF     <- ifelse(res.soma$Freq1 > .5, 1 - res.soma$Freq1, res.soma$Freq1)

## create Allele2
res.soma$Allele2 <- apply(res.soma[, c("MarkerName", "Allele1")], 1, function(x){
  ## split MarkerName
  ii <- strsplit(x[1], "_")[[1]][2:3]
  return(ii[ii != x[2]])
})

## keep only what is needed
res.soma         <- res.soma[, c("MarkerName", "rsid", "chr", "pos", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "Pvalue", "TotalSampleSize")]

## subset tor region of interest
res.soma         <- subset(res.soma, pos >= pos.s & pos <= pos.e)

#-----------------------------------------#
##--        lift coordinates over      --##
#-----------------------------------------#

## change chromosome if needed
if(chr.s == 23){
  chr.s <- "X"
}

# BiocManager::install("liftOver")
require(liftOver)

## import chain
chain      <- import.chain("hg19ToHg38.over.chain")

## create GRanges object
grObject   <- GRanges(seqnames = paste0("chr", chr.s), ranges=IRanges(start = res.soma$pos, end = res.soma$pos, names=res.soma$MarkerName))

## now do the mapping
tmp        <- as.data.frame(liftOver(grObject, chain))
## rename
names(tmp) <- c("group", "MarkerName", "seqnames", "pos.hg38", "end", "width", "strand")
## add to the data
res.soma   <- merge(res.soma, tmp[, c("MarkerName", "pos.hg38")])

#-----------------------------------------#
##--        import tissue data         --##
#-----------------------------------------#

## get all GTEx files
ii          <- dir("/home/mdp50/rds/rds-rjh234-cc-mrc-epid/People/Ellie/GTEx_v8_download/GTEx_Analysis_v8_EUR_sQTL_all_associations_csv/") 

## get unqiue tissues
tissue      <- unique(sapply(ii, function(x) strsplit(x, "\\.")[[1]][1]))
# tissue      <- tissue[-which(tissue == "slurm_output")] 

## get the Stats
res.gtex    <- lapply(tissue, function(x){
  ## get all information for the respective gene
  tmp        <- paste0("grep ", intr, " ", "/home/mdp50/rds/rds-rjh234-cc-mrc-epid/People/Ellie/GTEx_v8_download/GTEx_Analysis_v8_EUR_sQTL_all_associations_csv/",
                       x, ".v8.EUR.sqtl_allpairs.chr",chr.s,".csv | grep ", gene, " - ")
  # print(tmp)
  tmp        <- data.table::fread(cmd=tmp, sep=",", header=F, data.table = F)
  ## proceed only if gene is available in the respective tissue
  if(nrow(tmp)>1){
    ## add header
    names(tmp) <- c("nr", "phenotype_id", "variant_id", "tss_distance", "maf", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se")
    ## add tissue
    tmp$tissue <- x
    ## compute N
    tmp$N      <- (tmp$ma_count/tmp$maf)/2
    return(tmp[, c("phenotype_id", "variant_id", "N", "slope", "slope_se", "pval_nominal", "tissue")])
  }else{
    cat(gene, "not found in", x, "\n")
  }
  
})
## combine into one data frame
res.gtex <- do.call(rbind, res.gtex)

## add gene and intron instead of phenotype ID from LeafCutter
res.gtex$intron  <- intr
res.gtex$gene_id <- gene

## delete phenotype ID
res.gtex$phenotype_id <- NULL

## reshape
library(tidyr)
res.gtex <- pivot_wider(res.gtex, id_cols = c("variant_id", "intron", "gene_id"),
                        values_from = c("slope", "slope_se", "pval_nominal", "N"),
                        names_from = "tissue", names_sep = ".")
## convert to data frame
res.gtex <- data.frame(res.gtex)

## subset to tissues with available gene expression
tissue   <- grep("slope_se", names(res.gtex), value=T)
tissue   <- gsub("slope_se\\.", "", tissue)

## create position and allele columns
res.gtex$pos.hg38 <- as.numeric(sapply(res.gtex$variant_id, function(x) strsplit(x, "_")[[1]][2]))
# coding allele
res.gtex$Allele1  <- sapply(res.gtex$variant_id, function(x) strsplit(x, "_")[[1]][4])
res.gtex$Allele2  <- sapply(res.gtex$variant_id, function(x) strsplit(x, "_")[[1]][3])

## delete non-biallelic variants
ii <- table(res.gtex$pos.hg38)
ii <- names(ii[which(ii > 1)])
if(length(ii) > 0){
  res.gtex <- res.gtex[-which(res.gtex$pos %in% ii),]
}

#-----------------------------------------#
##--         combine everything        --##
#-----------------------------------------#

res.all              <- merge(res.gtex, res.soma, by="pos.hg38", suffixes = c(".gtex", ".soma"))

## align effect estimates
res.all$Allele1.soma <- toupper(res.all$Allele1.soma)

## recode INDELs
res.all$Allele1.gtex <- apply(res.all[, c("Allele1.gtex", "Allele2.gtex")], 1, function(x){
  if(nchar(x[1]) > 1 | nchar(x[2]) > 1){
    if(nchar(x[1]) > 1){
      return("I")
    }else{
      return("D")
    }
  }else{
    return(x[1])
  }
})

## create aligned effect estimate
res.all$Effect.aligned <- ifelse(res.all$Allele1.gtex == res.all$Allele1.soma, res.all$Effect, -res.all$Effect)

#-----------------------------------------#
##--       import genotype data        --##
#-----------------------------------------#

## read the dosage file
require(data.table)
tmp           <- fread(paste0("tmp_input//tmp.", soma, ".", ifelse(chr.s == "X", 23, chr.s), ".", pos.s, ".", pos.e, ".dosage"), sep=" ", header=T, data.table=F)
## transpose
rownames(tmp) <- tmp$rsid
## store allele information to realign effect estimates afterwards
tmp.info      <- tmp[, 1:6]
tmp           <- t(tmp[,-c(1:6)])
## retransform to data set (keep X in mind for variable names)
tmp           <- data.frame(ID_1=rownames(tmp), tmp)

## create another column to info to map names from the SNP data set
tmp.info$id         <- sapply(tmp.info$rsid, function(x) ifelse(substr(x, 1, 2) == "rs", x, paste0("X", gsub(":", ".", x))))
## edit some IDs (X-chromosome)
tmp.info$id         <- gsub("XX", "X", tmp.info$id)
tmp.info$id         <- gsub("XAffx-", "Affx.", tmp.info$id)
## create MarkerName column as well
tmp.info$MarkerName <- apply(tmp.info, 1, function(x){
  paste0("chr", as.numeric(x[1]), ":", as.numeric(x[4]), "_", paste(sort(x[5:6]), collapse = "_"))
})
## reduce to SNPs in the stats
tmp.info            <- subset(tmp.info, MarkerName %in% res.soma$MarkerName)

## rename
pheno <- tmp
rm(tmp); gc()

## add to ease mapping of LD matrix
res.all     <- merge(res.all, tmp.info[, c("MarkerName", "id")], by="MarkerName")

## omit NAs
res.all     <- na.omit(res.all)

## sort
res.all     <- res.all[order(res.all$pos.hg38),]

## store the top hit for the protein
top.snp     <- which.max(abs(res.soma$Effect/res.soma$StdErr))
top.snp     <- res.soma$MarkerName[top.snp]

## get new top snp in the merged data set
top.merged  <- which.max(abs(res.all$Effect/res.all$StdErr))
top.merged  <- res.all$MarkerName[top.merged]

## compute the LD between both
ld.sens     <- cor(pheno[, tmp.info$id[which(tmp.info$MarkerName == top.snp)]], pheno[, tmp.info$id[which(tmp.info$MarkerName == top.merged)]], m="p", u="p")^2

#-----------------------------------------#
##--            run coloc              --##
#-----------------------------------------#

## compute MAF
res.all$MAF <- ifelse(res.all$Freq1 > .5, 1-res.all$Freq1, res.all$Freq1)

require(coloc)

## run coloc in parallele
res.coloc <- lapply(tissue, function(x){
  
  print(x)
  
  #-----------------------------------------#
  ##-- 	         sanity check            --##
  #-----------------------------------------#
  
  ## get the top SNP for the outcome
  io     <- res.all$MarkerName[which.max(abs(res.all[, paste0("slope.",x)]/res.all[, paste0("slope_se.",x)]))]
  
  ## keep names
  ld.top <- cor(pheno[, tmp.info$id[which(tmp.info$MarkerName == top.merged)]], pheno[, tmp.info$id[which(tmp.info$MarkerName == io)]], m="p", u="p")^2
  
  #-----------------------------------------#
  ##-- 	            run coloc            --##
  #-----------------------------------------#
  
  ## prepare input for protein
  D1          <- list(beta=res.all$Effect.aligned, 
                      varbeta=res.all$StdErr^2, 
                      type="quant", 
                      N=10708, 
                      sdY=1,
                      MAF=res.all$MAF,
                      snp=res.all$id,
                      position=1:nrow(res.all))
  
  ## prepare input for GTEx
  D2          <- list(beta=res.all[, paste0("slope.",x)], 
                      varbeta=res.all[, paste0("slope_se.",x)]^2,
                      # pvalues=res.all[, paste0("pval_nominal.",x)],
                      type="quant",
                      N=max(res.all[, paste0("N.", x)]),
                      MAF=res.all$MAF,
                      snp=res.all$id,
                      sdY=1,
                      position=1:nrow(res.all))
  
  ## do naive coloc as well
  naive.coloc <- coloc.signals(D1, D2, method="single", p12=1e-5)
  
  ## add checks to the data
  naive.coloc$summary$ld.check.top  <- ld.top
  naive.coloc$summary$ld.check.sens <- ld.sens
  
  ## add the trait id
  naive.coloc$summary$tissue        <- x
  
  #-----------------------------------------#
  ##-- 	        draw selected            --##
  #-----------------------------------------#
  
  if(naive.coloc$summary$PP.H4.abf > .7 | naive.coloc$summary$ld.check.top > .8){
    source("scripts/plot_coloc_results.R")
    png(paste0("graphics/", soma,".", gene, ".", x, ".", chr.s, ".", pos.s, ".", pos.e, ".png"), width=16, height=8, units="cm", res=100)
    par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=.01, cex.main=.6, font.main=2)
    layout(matrix(c(1,1,2,3),2,2))
    plot.regional.coloc(res.all, naive.coloc$summary, pheno, "Effect", "StdErr", 
                        paste0("slope.", x), paste0("slope_se.", x),
                        gene, soma, tmp.info)
    dev.off()
  }
  
  
  ## write results to file
  naive.coloc <- naive.coloc$summary
  # print(naive.coloc)
  
  ## give back data for top pQTL
  ii          <- which(res.all$MarkerName == top.merged) 
  
  ## add top SNP from phenotype with estimates
  naive.coloc <- data.frame(naive.coloc, 
                            res.all[ii, c("intron", "gene_id", "chr", "pos.hg38", "pos", "variant_id", "MarkerName", "rsid",
                                          "Allele1.gtex", "Allele2.gtex", "MAF", "Effect.aligned", "StdErr", "Pvalue",
                                          paste(c("slope", "slope_se", "pval_nominal"), x, sep="."))])
  
  ## edit names
  names(naive.coloc) <- gsub(paste0("\\.", x),  "", names(naive.coloc))
  return(naive.coloc)
  
})
res.coloc <- do.call(rbind, res.coloc)

## store results
write.table(res.coloc, paste0("tissue_coloc_sQTL_cond/", soma,".", intr, ".", gene, ".", chr.s, ".", pos.s, ".", pos.e, ".txt"), sep="\t", row.names=F)
