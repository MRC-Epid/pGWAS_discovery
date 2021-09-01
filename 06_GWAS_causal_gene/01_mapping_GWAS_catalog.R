##############################################
####      GWAS catalog mapping pQTLs      ####
#### Maik Pietzner             25/01/2021 ####
##############################################

rm(list=ls())
setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/19_GWAS_catalog_v2/data")
options(stringsAsFactors = F)
load(".RData")

############################################
####  import results from both studies  ####
############################################

## SOMAscan (already reduced to regional sentinels!)
res.soma  <- read.table("../../05_variant_annotation/output/Regional.sentinels.pGWAS.SomaLogic.annotated.20200813.txt", sep="\t", header=T)
cond.soma <- read.table("../../05_variant_annotation/output/Conditional.results.pGWAS.SomaLogic.annotated.20200813.txt", sep="\t", header=T)

## create file
query     <- data.frame(rsid=unique(c(res.soma$rsid, cond.soma$rsid)))

#######################################################
####           all variants in strong LD           ####
#######################################################

## files for LD
ii <- dir("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/04_ld_files/output/")
ii <- grep("ld$", ii, value=T)

## create LD-matrix for all SNPs
require(data.table)
ld.matrix <- lapply(ii, function(x){
  tmp    <- data.frame(fread(paste0("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/04_ld_files/output/", x), sep=" ", header=T))
  ## do two subsets to get only SNPs of interest
  tmp    <- subset(tmp, SNP_A %in% query$rsid & R2 >= .8)
  print(head(tmp))
  return(tmp)
})
ld.matrix     <- do.call(rbind, ld.matrix)

#################################################
####            GWAS catalogue data          ####
#################################################

## load data from the GWAS catalogue: download 25/01/2021
require(data.table)
gwas.catalogue           <- fread("alternative", sep="\t", header=T, data.table = F)

## rename some columns
names(gwas.catalogue)[8] <- "TRAIT"

## prune GWAS catalogue data
gwas.catalogue           <- subset(gwas.catalogue, !is.na(`OR or BETA`) & is.finite(`OR or BETA`))
## generate risk allele and drop everything w/o this information
gwas.catalogue$riskA     <- sapply(gwas.catalogue$`STRONGEST SNP-RISK ALLELE`, function(x) strsplit(x,"-")[[1]][2])
## delete spaces at the end
gwas.catalogue$riskA     <- trimws(gwas.catalogue$riskA, which = "b")
## drop interaction entries
ii                       <- grep("[0-9]", gwas.catalogue$riskA)
gwas.catalogue           <- gwas.catalogue[-ii,]
## only genome-wide significant ones
gwas.catalogue           <- subset(gwas.catalogue, PVALUE_MLOG > 7.3)
## N = 120,912 entries

## omit some data, e.g. protein data obtained on multiplex platforms
gwas.catalogue           <- subset(gwas.catalogue, !(TRAIT %in% c("Blood protein levels", "Protein biomarker")))
## N = 114,412

## prune some specific studies
gwas.catalogue           <- subset(gwas.catalogue, !(PUBMEDID %in% c(33067605, 27989323, 27989323, 27989323, 29234017, 31320639, 25147954)))

#################################################
####     add information if any cis-pQTL     ####
#################################################

## add ld_matrix
gwas.catalogue.cis <- lapply(1:nrow(gwas.catalogue), function(x){
  
  ## identify all possible mapping pQTLs (careful: many-to-many mappings possible)
  ii <- subset(ld.matrix, SNP_B == gwas.catalogue$SNPS[x])

  ## report only if any
  if(nrow(ii)>0){
    
    ## get mapping pQTLs
    il <- subset(res.soma, rsid %in% ii$SNP_A & cis_trans == "cis")
    tp <- "lead" 
    
    ## in case no lead signal
    if(nrow(il) == 0){
      il <- subset(cond.soma, rsid %in% ii$SNP_A & cis_trans == "cis")
      tp <- "secondary"
    }
    
    ## ensure a cis hit
    if(nrow(il)>0){
      ## take the strongest one as a reference
      if(tp == "lead"){
        pqtl <- il$rsid[which.max(abs(il$Effect/il$StdErr))]
      }else{
        pqtl <- il$rsid[which.max(abs(il$beta.qc.joint/il$se.qc.joint))]
      }
      
      ## get all possible phenotype IDs
      p.id <- paste(unique(il$pheno), collapse = "|")
      ## get all possible protein targets
      p.ta <- paste(unique(il$Target), collapse = "|")
      ## get all possible gene names
      g.n  <- paste(unique(il$Gene.Name.Name), collapse = "|")
      ## ENSEMBL
      e.n  <- paste(unique(il$ensembl_gene_id), collapse = "|")
      
      ## return data frame
      return(data.frame(gwas.catalogue[x,], pQTL=pqtl, R2=ii$R2[which(ii$SNP_A == pqtl)], pheno = p.id,
                        Target=p.ta, protein.gene.name=g.n, protein.gene.ensembl=e.n, type=tp))
    }
  }
})
gwas.catalogue.cis <- do.call(rbind, gwas.catalogue.cis)

#################################################
####  reduce to protein-encoding regions     ####
#################################################

## create colum to indicate matching
gwas.catalogue.cis$match.reported <- apply(gwas.catalogue.cis[, c("REPORTED.GENE.S.", "protein.gene.name")], 1, function(x){
  if(!is.na(x[1]) & !is.na(x[2])){
    if(x[1] == x[2]){
      return("exact_match")
    }else{
      cat("---------------------\n")
      ## search for partial matching (improve search to avoid false negative results)
      rg   <- strsplit(x[1], ", ")[[1]]
      pg   <- strsplit(x[2], "|", fixed=T)[[1]]
      ## test whether any of the protein genes is found in the reported genes
      ii   <- lapply(pg, function(c) c %in% rg)
      ii   <- sum(unlist(ii))
      # x[2] <- paste(unlist(lapply(strsplit(x[2], "|", fixed = T)[[1]], function(c) paste0("^", c, "$"))), collapse = "|")
      # ii   <- grep(x[2], x[1])
      if(ii > 0){
        return("partial_match")
      }else{
        return("no_match")
      }
    }
  }else{
    return("missing")
  }
})

## create common identifier by the GWAS catalog based on proximity
gwas.catalogue.cis$ENSEMBL.GWAS.CATALOG <- apply(gwas.catalogue.cis[, c("UPSTREAM_GENE_ID", "DOWNSTREAM_GENE_ID", "SNP_GENE_IDS")], 1, function(x){
  ## make one common entry
  return(paste(unique(x[which(x!= "")]), collapse = ", "))
})

## compare again
gwas.catalogue.cis$match.ensembl <- apply(gwas.catalogue.cis[, c("ENSEMBL.GWAS.CATALOG", "protein.gene.ensembl")], 1, function(x){
  if(!is.na(x[1]) & !is.na(x[2])){
    if(x[1] == x[2]){
      return("exact_match")
    }else{
      ## search for partial matching
      ii <- grep(x[2], x[1])
      if(length(ii) > 0){
        return("partial_match")
      }else{
        return("no_match")
      }
    }
  }else{
    return("missing")
  }
})


#################################################
####          make locus specific            ####
#################################################

## identify all unqiue protein-encoding regions
tmp              <- unique(gwas.catalogue.cis[, c("pQTL", "pheno", "type")]) ## 590 regions
cis.regions.gwas <- lapply(1:nrow(tmp), function(x){
  
  ## identify all mapping GWAS results
  ii <- subset(gwas.catalogue.cis, pQTL == tmp$pQTL[x] & pheno == tmp$pheno[x])
  
  ## refine the list of reported genes
  jj <- unique(unlist(lapply(ii$REPORTED.GENE.S., function(x) strsplit(x, ", | - ")[[1]])))
  kk <- unique(unlist(lapply(ii$ENSEMBL.GWAS.CATALOG, function(x) strsplit(x, ", | - ")[[1]])))
  
  ## return compressed data set
  return(data.frame(pQTL=tmp$pQTL[x], pheno=tmp$pheno[x], type=tmp$type[x],
                    ## get rid of commas and hyphens before reporting back
                    reported.genes=paste(sort(jj), collapse = "|"),
                    ensembl.gwas=paste(sort(kk), collapse = "|"),
                    match.reported=paste(sort(unique(ii$match.reported)), collapse = "|"),
                    match.ensembl=paste(sort(unique(ii$match.ensembl)), collapse = "|"),
                    protein.gene.name=ii$protein.gene.name[1],
                    protein.gene.ensembl=ii$protein.gene.ensembl[1],
                    n.traits=nrow(ii),
                    gwas.id=paste(unique(ii$STUDY.ACCESSION), collapse = "|"),
                    gwas.traits=paste(unique(ii$TRAIT), collapse = "|"),
                    gwas.mapped.traits=paste(unique(ii$MAPPED_TRAIT), collapse = "|")))
})
cis.regions.gwas <- do.call(rbind, cis.regions.gwas)

## write to file to verify mappings manually
cis.regions.gwas$matched.gene.locus <- sapply(cis.regions.gwas$match.reported, function(x){
  x <- strsplit(x, "|", fixed = T)[[1]]
  if("exact_match" %in% x){
    return("exact_match")
  }else if("partial_match" %in% x){
    return("partial_match")
  }else{
    return("no_match")
  }
})

## write to file
write.table(gwas.catalogue.cis, "Cis.pQTL.mapping.GWAS.catalog.20210125.txt", sep="\t", row.names=F)
write.table(cis.regions.gwas, "Cis.locus.mapping.GWAS.catalog.20210125.txt", sep="\t", row.names=F)
