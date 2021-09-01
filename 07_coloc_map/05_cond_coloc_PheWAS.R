#!/usr/bin/env Rscript

## script to run extensive PheWAS analysis for a given regional association
## conditioning on the top hit
## Maik Pietzner 31/10/2020
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/14_PheWAS_2/")

## --> import parameters <-- ##

## aptamer
soma  <- args[1]
## chromosome
chr.s <- as.numeric(args[2])
## start position of the region
pos.s <- as.numeric(args[3])
## end position of the region
pos.e <- as.numeric(args[4])

#-----------------------------------------#
##-- 	       import protein data       --##
#-----------------------------------------#

## load conditional stats
res.soma        <- paste0("input_cond/", soma, ".", chr.s, ".", pos.s, ".", pos.e, ".cma.cojo")
res.soma        <- data.table::fread(res.soma, sep = "\t", header = T, data.table = F)
## do some renaming to align with the results (rename bC to actual beta!)
names(res.soma) <- c("chr", "MarkerName", "pos", "Allele1", "Freq1", "beta_ma", "se_ma", "pval_ma", "TotalSampleSize", "freq_geno", "Effect", "StdErr", "Pvalue")
## add rsID
map             <- data.table::fread("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/Fenland.OMICs.MarkerName.rsID.chr.pos.txt", sep=" ", header=T, data.table=F)
## add to the data
res.soma        <- merge(res.soma, map)
## omit missing values
res.soma        <- subset(res.soma, !is.na(Effect))
rm(map); gc()

## apply MAF filter
res.soma$MAF    <- ifelse(res.soma$Freq1 > .5, 1 - res.soma$Freq1, res.soma$Freq1)

## define break criteria here (min p-value!!)
p.min           <- min(as.numeric(res.soma$Pvalue), na.rm = T)

#-----------------------------------------#
##--   only do if sufficient evidence  --##
#-----------------------------------------#

if(p.min < 1e-6){
  
  #-----------------------------------------#
  ##--       import genotype data        --##
  #-----------------------------------------#
  
  ## read the dosage file
  require(data.table)
  tmp           <- fread(paste0("tmp_input//tmp.", soma, ".", chr.s, ".", pos.s, ".", pos.e, ".dosage"), sep=" ", header=T, data.table=F)
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
  # ## compute full correlation matrix (do only for those needed)
  # ld_mat        <- cor(tmp[,-1])^2
  # rm(tmp); gc()
  
  ## rename
  pheno <- tmp
  rm(tmp); gc()
  
  #-----------------------------------------#
  ##--  define top SNP and get proxies   --##
  #-----------------------------------------#
  
  ## identify top SNP by Z-score 
  top.snp    <- res.soma$MarkerName[which.max(abs(res.soma$Effect/res.soma$StdErr))]
  
  ## get proxies (use rsIDs)
  proxy.snps <- cor(pheno[, tmp.info$id[which(tmp.info$MarkerName == top.snp)]], pheno[, tmp.info$id])^2
  proxy.snps <- data.frame(id=colnames(proxy.snps), R2=t(proxy.snps))
  proxy.snps <- subset(proxy.snps, R2 >= .8)
  ## ease merging
  proxy.snps <- merge(proxy.snps, tmp.info)
  ## ensure overlap with summary stats
  proxy.snps <- subset(proxy.snps, MarkerName %in% res.soma$MarkerName)
  
  #-----------------------------------------#
  ##--  perform PheWAS for all proxies   --##
  #-----------------------------------------#
  
  ## load package to for API to https://gwas.mrcieu.ac.uk/
  require(ieugwasr)
  
  ## query everything (omit FinnGen and BBJ as sample sizes are missing)
  tmp.phewas <- phewas(variants=proxy.snps$rsid, pval=1e-6, batch=c("ebi-a", "eqtl-a", "ieu-a", "ieu-b", "met-a", "met-b", "met-c", "ubm-a", "ukb-b"))
  
  ## define break criteria here as well
  if(nrow(tmp.phewas)>0){
    
    ## take the strongest association for each phenotype
    tmp     <- tmp.phewas[order(tmp.phewas$trait, abs(tmp.phewas$beta/tmp.phewas$se), decreasing = T),]
    tmp$ind <- unlist(lapply(unique(tmp$trait), function(c) seq(1, nrow(subset(tmp, trait == c)))))
    tmp     <- subset(tmp, ind == 1)
    tmp$ind <- NULL
    
    ## test whether any remaining phenotypes may exist
    if(nrow(tmp)>0){
      tmp2 <- subset(tmp.phewas, rsid %in% proxy.snps$rsid & !(trait %in% tmp$trait))
    }else{
      tmp2 <- subset(tmp.phewas, rsid %in% proxy.snps$rsid)
    }
    
    ## do iterative subsetting to keep unqie entries only
    if(nrow(tmp2)>1){
      ## take the strongest association for each phenotype
      tmp2     <- tmp2[order(tmp2$trait, abs(tmp2$beta/tmp2$se), decreasing = T),]
      tmp2$ind <- unlist(lapply(unique(tmp2$trait), function(c) seq(1, nrow(subset(tmp2, trait == c)))))
      tmp2     <- subset(tmp2, ind == 1)
      tmp2$ind <- NULL
    }
    
    ## combine both
    tmp                 <- rbind(tmp, tmp2)
    ## add SNP information
    tmp                 <- merge(proxy.snps[, c("rsid", "id", "MarkerName", "R2")], tmp, by="rsid", suffixes = c(".snp", ".ieu"))
    ## rename and add lead SNP
    tmp$MarkerName.lead <- top.snp
    tmp$rsid.lead       <- tmp.info$rsid[which(tmp.info$MarkerName == top.snp)]
    
    ## add protein stats 
    phewas.results      <- merge(res.soma, tmp, by=c("MarkerName", "rsid"), suffix=c(".pQTL", ".trait"))
    
    #-----------------------------------------#
    ##--    run coloc for all results      --##
    #-----------------------------------------#
    
    require(coloc)
    
    res.coloc <- lapply(1:nrow(phewas.results), function(x){
      
      ## obtain summary stats and trait info
      reg             <- paste0(chr.s, ":", pos.s,"-", pos.e)
      res.trait       <- associations(reg, phewas.results$id.ieu[x])
      ## just to make sure everything will work fine
      res.trait       <- subset(res.trait, !is.na(beta))
      ## make unique, since some data processing error
      res.trait       <- unique(res.trait)
      
      ## get meta information
      # tr.info         <- gwasinfo(phewas.results$id.ieu[x])
      tr.info         <- tryCatch(
        {
          gwasinfo(phewas.results$id.ieu[x])
        }, error=function(e){
          return(NA)
        })
      print(tr.info)
      
      ## merge (this might ignore INDELS)
      res.all              <- merge(res.soma, res.trait, by="rsid")
      ## remove non-biallelelic variants
      ii                   <- table(res.all$rsid)
      res.all              <- subset(res.all, rsid %in% names(ii[ii==1]))
      ## align effect estimates
      res.all$beta.aligned <- ifelse(toupper(res.all$Allele1) == res.all$ea, res.all$beta, -res.all$beta)
      
      #-----------------------------------------#
      ##-- 	         sanity check            --##
      #-----------------------------------------#
      
      ## top signal for SOMAscan in all
      it    <- res.all$MarkerName[which.max(abs(res.all$Effect/res.all$StdErr))]
      ## top signal SOMAscan
      is    <- res.soma$MarkerName[which.max(abs(res.soma$Effect/res.soma$StdErr))]           
      ## get the top SNP for the outcome
      io    <- res.all$MarkerName[which.max(abs(res.all$beta/res.all$se))]
      
      ## keep names
      isnps <- sapply(c(it, is, top.snp, io), function(x) tmp.info$id[which(tmp.info$MarkerName == x)])
      
      ## get the LDs for each of the hits
      it <- cor(pheno[, isnps[1]], pheno[, isnps[3]])^2
      is <- cor(pheno[, isnps[2]], pheno[, isnps[3]])^2
      io <- cor(pheno[, isnps[4]], pheno[, isnps[1]])^2
      
      #-----------------------------------------#
      ##-- 	            run coloc            --##
      #-----------------------------------------#
      
      ## add to ease mapping of LD matrix
      res.all     <- merge(res.all, tmp.info[, c("rsid", "id")], by="rsid", suffix=c(".trait", ".snp"))
      
      ## order by position
      res.all     <- res.all[order(res.all$pos),]
      # print(head(res.all))
      # print(summary(res.all))
      
      ## prepare input
      D1          <- list(beta=res.all$Effect, varbeta=res.all$StdErr^2, 
                          type="quant", 
                          N=max(res.all$TotalSampleSize, na.rm=T), 
                          sdY=1,
                          MAF=res.all$MAF,
                          snp=res.all$id.snp,
                          position=1:nrow(res.all))
      
      ## try out
      if(!("ncase" %in% names(tr.info))){
        print("use quantitative")
        D2          <- list(beta=res.all$beta.aligned, varbeta=res.all$se^2, 
                            type="quant", 
                            # N=tr.info$sample_size,
                            N=max(as.numeric(res.all$n), na.rm=T),
                            sdY=1,
                            MAF=res.all$MAF,
                            snp=res.all$id.snp,
                            position=1:nrow(res.all))
      }else{
        ## binary outcome
        D2          <- list(beta=res.all$beta.aligned, varbeta=res.all$se^2, 
                            type="cc",
                            s=tr.info$ncase/(tr.info$ncontrol+tr.info$ncase), 
                            N=tr.info$sample_size,
                            MAF=res.all$MAF,
                            snp=res.all$id.snp,
                            position=1:nrow(res.all))
      }
      
      ## do naive coloc as well
      naive.coloc <- coloc.signals(D1, D2, method="single", p12=1e-6)
      
      ## add checks to the data
      naive.coloc$summary$ld.lead.soma  <- it
      naive.coloc$summary$ld.check.top  <- io
      naive.coloc$summary$ld.check.sens <- is
      
      ## add the trait id
      naive.coloc$summary$id.ieu        <- phewas.results$id.ieu[x]
      
      #-----------------------------------------#
      ##-- 	        draw selected            --##
      #-----------------------------------------#
      
      if(naive.coloc$summary$PP.H4.abf > .7 | naive.coloc$summary$ld.check.top > .8){
        source("scripts/plot_coloc_results.R")
        png(paste0("graphics_cond//", soma, ".", phewas.results$id.ieu[x], ".", chr.s, ".", pos.s, ".", pos.e, ".png"), width=16, height=8, units="cm", res=100)
        par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=.01, cex.main=.6, font.main=2)
        layout(matrix(c(1,1,2,3),2,2))
        plot.regional.coloc(res.all, naive.coloc$summary, pheno, phewas.results$id.ieu[x], soma, tmp.info)
        dev.off()
      }
      
      
      ## write results to file
      return(naive.coloc$summary)
      
    })
    res.coloc <- do.call(rbind, res.coloc)
    
    #-----------------------------------------#
    ##--  	combine and store the data     --##
    #-----------------------------------------#
    
    ## combine both
    phewas.results       <- merge(phewas.results, res.coloc[, c("id.ieu", "nsnps", "hit1", "hit2", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "ld.lead.soma", "ld.check.top", "ld.check.sens")])
    
    ## add phenotype
    phewas.results$pheno <- soma
    
    ## write to file
    write.table(phewas.results, paste("output_cond//results.phewas", soma, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
    
    
  }else{
    cat("found no evidence for", soma, "on chr", chr.s, "between", pos.s, "and", pos.e, "\n")
    ## write to file
    write.table(data.frame(pheno=soma, id.ieu=NA), paste("output_cond//results.phewas", soma, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
  }
  
}else{
  cat("found no evidence for", soma, "on chr", chr.s, "between", pos.s, "and", pos.e, "\n")
  ## write to file
  write.table(data.frame(pheno=soma, id.ieu=NA), paste("output_cond//results.phewas", soma, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
}

