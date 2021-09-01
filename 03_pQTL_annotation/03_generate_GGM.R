#!/usr/bin/env Rscript

## script to generate a GGM for SOMAmer data
## Maik Pietzner 16/07/2020
rm(list=ls())

## little options
options(stringsAsFactors = F)

setwd("/rds/project/rjh234/rds-rjh234-mrc-epid/Studies/People/Maik/pGWAS_SomaLogic/02_protein_GGM/data/")

## load the data needed
load("Residuals.GGM.RData")

## get the names of the aptamers
soma <- names(tmp)

## create the network
require("GeneNet")
protein.ggm      <- ggm.estimate.pcor(tmp)
protein.ggm      <- network.test.edges(protein.ggm, plot=F)

## add names
protein.ggm$var1 <- soma[protein.ggm$node1]
protein.ggm$var2 <- soma[protein.ggm$node2]

## subset to significant edges
protein.ggm      <- subset(protein.ggm, pval < .05/nrow(protein.ggm))

## store the final GGM
write.table(protein.ggm, "SOMAmer.GGM.txt", sep="\t", row.names=F)
