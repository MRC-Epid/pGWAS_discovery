################################################
#### naive cis-pQTL coloc PheWAS            ####
#### Maik Pietzner               29/10/2020 ####
################################################

rm(list=ls())
setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/14_PheWAS_2/input/")
options(stringsAsFactors = F)
load(".RData")

##############################################
####        load SomaLogic targets        ####
##############################################

## import prep from COVID work
soma.targets <- read.table("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/COVID19_v2/02_naive_cis_coloc/input/SomaLogic.Targets.cis.regions.txt", sep="\t", header=T)

##############################################
####        create input for coloc        ####
##############################################

## subet to promising regions
coloc.regions        <- soma.targets
coloc.regions$Pvalue <- as.numeric(soma.targets$Pvalue) 
coloc.regions        <- subset(coloc.regions, Pvalue < 1e-5)

## create input file
write.table(coloc.regions[, c("pheno", "chr", "region_start", "region_end")], "coloc.regions.txt", sep="\t", row.names=F, col.names=F, quote=F)

##############################################
####            import results            ####
##############################################

ii        <- dir("../output/")
res.coloc <- lapply(ii, function(x){
  ## read results
  tmp              <- read.table(paste0("../output/", x), sep="\t", header=T)
  ## add some identifier
  tmp$region_start <- as.numeric(strsplit(x, "\\.")[[1]][5]) 
  tmp$region_end   <- as.numeric(strsplit(x, "\\.")[[1]][6])
  tmp$chr          <- as.numeric(strsplit(x, "\\.")[[1]][4])
  ## indicate no coloc
  if(is.na(tmp$id.ieu[1])){
    tmp$id.ieu <- "no evidence"
  }
  return(tmp)
})
require(plyr)
res.coloc <- do.call(rbind.fill, res.coloc)

## combine to see results
res.coloc <- merge(coloc.regions, res.coloc, by=c("pheno", "chr", "region_start", "region_end"), all = T, suffixes = c(".soma", ".covid"))

##############################################
####            filter results            ####
##############################################

## drop missing runs
res.coloc           <- subset(res.coloc, !is.na(id.ieu) & id.ieu != "no evidence")
## drop bad stats (missing top protein hit)
res.coloc           <- subset(res.coloc, ld.check.sens > .9 & ld.lead.soma > .8)

## flag interesting results
res.coloc$PP.cand   <- ifelse(res.coloc$PP.H4.abf > .8, 1, 0)
res.coloc$LD.cand   <- ifelse(res.coloc$ld.check.top > .8, 1, 0)
table(res.coloc$PP.cand, res.coloc$LD.cand) 

## write results to file
write.table(res.coloc, "Results.naive.coloc.cis.regions.SomaLogic.txt", sep="\t", row.names = F)

##############################################
####     prepare conditional analysis     ####
##############################################

## second attempt but now conditioning on top signal
## reduce to genome-wide primary signals!
coloc.cond.regions <- subset(coloc.regions, Pvalue < 5e-8)
## N = 1,906

## create a file for each SNP to condition on
for(j in 1:nrow(coloc.cond.regions)){
  write.table(coloc.cond.regions$MarkerName[j], 
              paste("../input_cond/cond", coloc.cond.regions$pheno[j], coloc.cond.regions$chr[j], coloc.cond.regions$region_start[j], coloc.cond.regions$region_end[j], "txt", sep="."), 
              sep="\t", row.names=F, col.names = F, quote = F)
}

## now write regions to file
write.table(coloc.cond.regions[, c("pheno", "chr", "region_start", "region_end")], "coloc.cond.regions.txt", sep="\t", row.names=F, col.names=F, quote=F)

## does include some rare variants

##############################################
####           import results             ####
##############################################

ii             <- dir("../output_cond//") ## 190 regions missing
res.coloc.cond <- lapply(ii, function(x){
  ## read results
  tmp              <- read.table(paste0("../output_cond//", x), sep="\t", header=T)
  ## add some identifier
  tmp$region_start <- as.numeric(strsplit(x, "\\.")[[1]][5]) 
  tmp$region_end   <- as.numeric(strsplit(x, "\\.")[[1]][6])
  tmp$chr          <- as.numeric(strsplit(x, "\\.")[[1]][4])
  ## indicate no coloc
  if(is.na(tmp$id.ieu[1])){
    tmp$id.ieu <- "no evidence"
  }
  return(tmp)
})
require(plyr)
res.coloc.cond <- do.call(rbind.fill, res.coloc.cond)

## combine to see results
res.coloc.cond <- merge(coloc.cond.regions, res.coloc.cond, by=c("pheno", "chr", "region_start", "region_end"), all = T, suffixes = c(".soma", ".covid"))

##############################################
####            filter results            ####
##############################################

## drop missing runs
res.coloc.cond           <- subset(res.coloc.cond, !is.na(id.ieu) & id.ieu != "no evidence")
## drop bad stats (missing top protein hit)
res.coloc.cond           <- subset(res.coloc.cond, ld.check.sens > .9 & ld.lead.soma > .8)

## flag interesting results
res.coloc.cond$PP.cand   <- ifelse(res.coloc.cond$PP.H4.abf > .8, 1, 0)
res.coloc.cond$LD.cand   <- ifelse(res.coloc.cond$ld.check.top > .8, 1, 0)

## write to file
write.table(res.coloc.cond, "Results.cond.top.coloc.cis.regions.SomaLogic.txt", sep="\t", row.names = F)
