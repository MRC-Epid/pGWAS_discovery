####################################################
#### Generation of aptamer network              ####
#### Maik Pietzner                   30/06/2020 ####
####################################################

rm(list=ls())
setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/02_protein_GGM/data/")
options(stringsAsFactors = F)
load(".RData")

############################################
####       import phenotypic data       ####
############################################

#--------------------------#
##-- Fenland Phenotypes --##
#--------------------------#

tmp         <- read.csv("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/COVID_19_targets/PHFENLANDR80001812017_collated12Feb2019.csv")
## create id
tmp$id      <- toupper(tmp$combined_id)
## all names to lower
names(tmp)  <- tolower(names(tmp))

#--------------------------#
##--       SOMAscan     --##
#--------------------------#

## import SOMAscan data
require(readstata13)
sl          <- read.dta13("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/COVID_19_targets/All_AMN_proteinAsAnalyzed-2018-10-04_R205.dta")
## phase 1 only
sl          <- subset(sl, phase == 1)
## N = 12,126

## add protein label
require(readxl)
label       <- data.frame(read_excel("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/COVID_19_targets/proteomics_druggablegenome_secreted_5Dec2019.xlsx"))
## only IDs needed
label       <- subset(label, !is.na(SomaId_v4))
label$pheno <- paste0("SeqId_", gsub("-", "_", label$proteomicsid))

## store targets of interest
target.soma <- label$pheno

## reduce to subset of interest
sl                <- sl[, c("SampleId", target.soma)]
## IvN transform the data
sl[, target.soma] <- apply(sl[, target.soma], 2, function(x){
  qnorm((rank(x ,na.last="keep")-0.5)/sum(!is.na(x)))
})

## add data on testsite
tt          <- data.table::fread("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/COVID_19_targets/17_var_decomp/data/Fenland.SomaLogic.TestSite.txt", sep="\t", header=T, data.table=F)
## edit ID
tt$SampleId <- ifelse(substr(tt$id, 1, 1) == "m", toupper(tt$id), tt$id)
## add to the data
sl          <- merge(sl, tt[,c("SampleId", "TestSite")], by="SampleId")

#-------------------------#
##--    combine all    --##
#-------------------------#

pheno <- merge(tmp, sl, by.x="id", by.y="SampleId")
rm(sl); gc()

#-------------------------#
##--    generate PCs   --##
#-------------------------#

pca.protein <- prcomp(pheno[, target.soma], scale. = T, center = T)
## add to the data
pheno       <- cbind(pheno, pca.protein$x[, paste0("PC", 1:10)])

#-------------------------#
##-- create residuals  --##
#-------------------------#

## correct for test site, age, sex
for(j in 1:length(target.soma)){
  ## define formula to be used
  ff <- paste(target.soma[j], "TestSite + ageattest_dm_attended + sex + PC1 + PC2 + PC3", sep="~")
  ## linear regression model
  ff <- resid(lm(as.formula(ff), data=pheno))
  ## add to the data (keep track of missing values)
  pheno[, paste("resid", target.soma[j], sep="_")]          <- NA
  pheno[names(ff), paste("resid", target.soma[j], sep="_")] <- ff
}

############################################
####  		  reduce redundancy 		####
############################################

## identify proteins with multiple mapping SOMAmers
tmp.lab <- names(which(table(label$Target) > 1))
## 190 --> add correlation of those, 3 proteins are targeted more than twice (DLK1, HSP 70, and Tenascin)
tmp.lab <- subset(label, Target %in% tmp.lab)$pheno
## compute correlation matrix
tmp.cor <- cor(pheno[, paste0("resid_", tmp.lab)])
## convert to list
require(reshape2)
tmp.cor[lower.tri(tmp.cor, diag = T)] <- NA
tmp.cor                               <- melt(tmp.cor)
## omit NA
tmp.cor                               <- subset(tmp.cor, !is.na(value))
## add more specific names
tmp.cor$Var1 <- gsub("resid_", "", tmp.cor$Var1)
tmp.cor$Var2 <- gsub("resid_", "", tmp.cor$Var2)
## add Target names
tmp.cor      <- merge(tmp.cor, label[, c("pheno", "Target", "ApparentKdM")], by.x="Var1", by.y="pheno")
tmp.cor      <- merge(tmp.cor, label[, c("pheno", "Target", "ApparentKdM")], by.x="Var2", by.y="pheno", suffixes = c(".1", ".2"))
## subset to same entries
tmp.cor      <- subset(tmp.cor, Target.1 == Target.2)
## for highly correlated aptamers choose the one with the lower ApparentKdM, ie. better binding
tmp.cor$sel  <- ifelse(tmp.cor$value < .5, "both", 
                       ifelse(tmp.cor$ApparentKdM.1 < tmp.cor$ApparentKdM.2, tmp.cor$Var1, tmp.cor$Var2))
## separate exclusion for DLK1, HSP 70, and Tenascin
tmp.cor      <- subset(tmp.cor, !(Target.1 %in% c("DLK1", "HSP 70", "Tenascin")))

## generate the selection of aptamers
tmp1  <- names(which(table(label$Target) == 1))       
tmp1  <- subset(label, Target %in% tmp1)$pheno
tmp2  <- unique(unlist(subset(tmp.cor, sel == "both")[, c("Var1", "Var2")]))
tmp3  <- tmp.cor$sel[which(tmp.cor$sel != "both")]
t.ggm <- unique(c(tmp1, tmp2, tmp3))
## add the last three individually
t.ggm <- c(t.ggm, "SeqId_6373_54", "SeqId_14237_1", "SeqId_6259_60")

############################################
####          create the network        ####
############################################

## run as a separate job on HPC
tmp <- pheno[, paste0("resid_", t.ggm)]
save(tmp, file="Residuals.GGM.RData")
rm(tmp)

## import the results
protein.ggm <- read.table("SOMAmer.GGM.txt", sep="\t", header=T)

## add more specific names
protein.ggm$var1 <- gsub("resid_", "", protein.ggm$var1)
protein.ggm$var2 <- gsub("resid_", "", protein.ggm$var2)
## add Target names
protein.ggm      <- merge(protein.ggm, label[, c("pheno", "Target", "TargetFullName")], by.x="var1", by.y="pheno")
protein.ggm      <- merge(protein.ggm, label[, c("pheno", "Target", "TargetFullName")], by.x="var2", by.y="pheno", suffixes = c(".1", ".2"))

## save and to do some plotting on Windows
write.table(protein.ggm, "SOMAmer.GGM.annotated.txt", sep="\t", row.names=F)

