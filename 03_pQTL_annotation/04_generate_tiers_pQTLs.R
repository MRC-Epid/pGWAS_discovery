############################################################
#### Classification of pQTLs                            ####
#### Maik Pietzner                           20/08/2020 ####
############################################################

rm(list=ls())
setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/07_pQTL_classification/input/")
options(stringsAsFactors = F)
load(".RData")

############################################
####  import results from both studies  ####
############################################

## SOMAscan (already reduced to regional sentinels!)
res.soma  <- read.table("../../05_variant_annotation/output/Regional.sentinels.pGWAS.SomaLogic.annotated.20200813.txt", sep="\t", header=T)
cond.soma <- read.table("../../05_variant_annotation/output/Conditional.results.pGWAS.SomaLogic.annotated.20200813.txt", sep="\t", header=T)

## load label for proteins
require(readxl)
label       <- data.frame(read_excel("../../../COVID_19_targets/proteomics_druggablegenome_secreted_5Dec2019.xlsx"))  
label       <- subset(label, !is.na(SomaId_v4))
label$pheno <- paste0("res_invn_X", gsub("-", "_", label$proteomicsid)) 

############################################
####   query all pQTLs for association  ####
############################################

## create file
query     <- data.frame(rsid=unique(c(res.soma$rsid, cond.soma$rsid)))

## write to file for query on HPC
write.table(query, "pQTL.query.txt", sep="\t", col.names = F, row.names = F, quote = F)

############################################
####          import look up            ####
############################################

res.cross        <- data.table::fread("cross.SOMAmer.lookup.5e8.txt", data.table = F, sep="\t", header = F)
## assign names
names(res.cross) <- c("id", "rsid", "MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "FreqMin", "FreqMax", "Effect",
                      "StdErr", "P.value", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "TotalSampleSize", "chr", "pos")
## generate ID
res.cross$pheno  <- gsub("_Fenland_MA_auto_chrX_filtered.txt.gz", "", res.cross$id)

############################################
####          import GO-terms           ####
############################################

## load mapping GO-terms (from UniProt download)
require(readxl)
uniprot.anno <- data.frame(read_excel("uniprot-yourlist_M202001106746803381A1F0E0DB47453E0216320D57D4CF8.xlsx", 1, guess_max = 5000))

## replace missing entries
uniprot.anno$Gene.ontology..biological.process.[is.na(uniprot.anno$Gene.ontology..biological.process.)] <- "none" 

## create pheno - GO.term mapping
pheno.go.anno        <- label[, c("pheno", "Uniprot_ID", "Target", "TargetFullName")]
pheno.go.anno$GO.bio <- apply(pheno.go.anno, 1, function(x){
  ## get all listed GO-terms in UniProt
  ii <- subset(uniprot.anno, Entry %in% strsplit(x[2], "\\|")[[1]])$Gene.ontology..biological.process.
  ii <- unlist(lapply(ii, function(c) strsplit(c,"; ")[[1]]))
  return(paste(unique(ii), collapse = "; "))
})

## try to identify common GO-terms
pqtl.pleiotropy <- lapply(unique(c(res.soma$MarkerName, cond.soma$MarkerName)), function(x){
  ## identify all associated SOMAmers
  ii <- subset(res.cross, MarkerName == x)$pheno
  print(x)
  if(length(ii) == 0){
    return(data.frame(MarkerName=x, num.SOMAmers=0, GO.term="none", GO.covered=0, GO.perc=0, category="none", pheno.ids="none"))
  }else if(length(ii) == 1){
    return(data.frame(MarkerName=x, num.SOMAmers=1, GO.term="none", GO.covered=0, GO.perc=0, category="specific_protein", pheno.ids=paste(ii, collapse = ", ")))
  }else{
    ## get GO-terms
    tmp <- subset(pheno.go.anno, pheno %in% ii) 
    ## test whether all may refer to the same target or complex
    jj <- table(unlist(lapply(tmp$Uniprot_ID, function(x) strsplit(x, "\\|")[[1]])))
    ## in case of a unique target, i.e. at least one UniProtID common to all SOMAmers
    if(sum(jj == length(ii))>0){
      return(data.frame(MarkerName=x, num.SOMAmers=length(ii), GO.term="none", GO.covered=0, GO.perc=0, category="specific_protein", pheno.ids=paste(ii, collapse = ", ")))
    }else{
      ## test for most common GO-term (or at least x%-coverage)
      jj  <- sort(table(unlist(lapply(tmp$GO.bio, function(x) strsplit(x, "; ")[[1]]))), decreasing = T)
      ## take the frist one
      jj  <- jj[1]
      return(data.frame(MarkerName=x, num.SOMAmers=length(ii), GO.term=names(jj), GO.covered=jj, GO.perc=jj/length(ii),
                        category=ifelse(jj == length(ii), "specific_pathway", "unspecific_effect"), 
                        pheno.ids=paste(ii, collapse = ", ")))
    }
  }
})
pqtl.pleiotropy <- do.call(rbind, pqtl.pleiotropy)
table(pqtl.pleiotropy$category)

## generate another category based on partial overlap
pqtl.pleiotropy$category <- ifelse(pqtl.pleiotropy$GO.perc > .5 & pqtl.pleiotropy$GO.covered > 2 & pqtl.pleiotropy$category != "specific_pathway", "specific_pathway_suggestive", pqtl.pleiotropy$category)

## some checks
pqtl.pleiotropy$signal   <- ifelse(pqtl.pleiotropy$MarkerName %in% res.soma$MarkerName, "sentinel", "secondary")

############################################
####              import GGM            ####
############################################

## import protein network
protein.ggm    <- read.table("Protein.GGM.edge.list.txt", sep="\t", header=T)

## import network-based annotations
ggm.anno       <- read.table("Stats.protein.GGM.txt", sep="\t", header=T)
## edit ID
ggm.anno$pheno <- gsub("SeqId_", "res_invn_X", ggm.anno$pheno)

## do annotation based on protein communities derived from the network
tmp         <- lapply(unique(c(res.soma$MarkerName, cond.soma$MarkerName)), function(x){
  ## identify all associated SOMAmers
  ii <- subset(res.cross, MarkerName == x)$pheno
  print(x)
  if(length(ii) == 0){
    return(data.frame(MarkerName=x, num.SOMAmers=0, community="none", community.n=0, community.c=0, category="none", pheno.ids="none"))
  }else if(length(ii) == 1){
    ## store information
    jj <- ggm.anno$community[which(ggm.anno$pheno == ii)]
    return(data.frame(MarkerName=x, num.SOMAmers=1, 
                      community=ifelse(!is.na(jj), jj, "none"),
                      community.n=1,
                      ## coverage of associated SOMAmers
                      community.c=ifelse(!is.na(jj), 1, 0),
                      category="specific_protein", pheno.ids=paste(ii, collapse = ", ")))
  }else{
    ## get communities
    tmp <- subset(ggm.anno, pheno %in% ii & !is.na(community)) 
    ## if any
    if(nrow(tmp) > 0){
      ## test whether all belong to the same community
      return(data.frame(MarkerName = x, num.SOMAmers = nrow(tmp), 
                        community = paste(unique(tmp$community), collapse = ", "),
                        community.n = length(unique(tmp$community)),
                        ## coverage of associated SOMAmers
                        community.c=nrow(tmp)/length(ii),
                        category=ifelse(length(unique(tmp$community)) == 1, "specific_pathway", "unspecific_effect"), 
                        pheno.ids=paste(ii, collapse = ", ")))
    }else{
      ## test whether all belong to the same community
      return(data.frame(MarkerName = x, num.SOMAmers = 0, 
                        community = "none",
                        community.n = 0,
                        ## coverage of associated SOMAmers
                        community.c=0,
                        category="not_covered", 
                        pheno.ids=paste(ii, collapse = ", ")))
    }
  }

})
tmp <- do.call(rbind, tmp)
table(tmp$category)

############################################
####            combine both            ####
############################################

res.tiers <- merge(pqtl.pleiotropy, tmp, by="MarkerName", suffixes = c(".GOterm", ".GGMnetwork"))

## write results to file
write.table(res.tiers, "Tiers.pQTLs.GOterms.GGMnetwork.txt", sep = "\t", row.names=F)

