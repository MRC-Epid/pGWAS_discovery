#####################################
## function to plot resutls from
## coloc

plot.regional.coloc <- function(res, sum.coloc, pheno, t_n, sl_n, tmp.info){
  
  ## 'res'       -- gwas stats
  ## 'sum.coloc' -- summary from coloc
  ## 'pheno'     -- data set to compute correlations
  ## 't_n'       -- name of the trait
  ## 'sl_n'      -- name for the SOMAmer
  ## 'tmp.info'  -- mapping file for SNPs
  
  ## subset to non-missing SNPs
  ii  <- which(!is.na(res[, "Effect"]) & !is.na(res[, "beta"]))
  res <- res[ii,]
  
  print(dim(sum.coloc))
  
  ## create log10p
  res$log10p.soma  <- -pchisq((res$Effect/res$StdErr)^2, df=1, lower.tail=F, log.p=T)/log(10)
  res$log10p.trait <- -pchisq((res[, "beta"]/res[, "se"])^2, df=1, lower.tail=F, log.p=T)/log(10)
  
  ## identify top SNP in the region
  ii.s   <- which.max(res$log10p.soma)
  ii.o   <- which.max(res$log10p.trait)
  ## choose LD based on protein SNP
  ii     <- ii.s  
  
  ## compute correlation
  ld.mat <- cor(pheno[, tmp.info$id[which(tmp.info$MarkerName == res$MarkerName[ii])]], pheno[, tmp.info$id])^2
  ld.mat <- data.frame(id=colnames(ld.mat), R2=t(ld.mat))
  ## ease merging
  ld.mat <- merge(ld.mat, tmp.info)

  ## add to the data for plotting
  res                   <- merge(res, ld.mat, by="MarkerName", all.x=T)
  res$R2[is.na(res$R2)] <- 0
  res$LD_cat            <- cut(res$R2, c(0, .2, .4, .6, .8, 1), include.lowest = T)
  ## add colour
  LD.colour             <- data.frame(LD_cat=levels(res$LD_cat), LD.col=rev(ggsci::pal_locuszoom()(5)))
  res                   <- merge(res, LD.colour)
  
  ## identify top SNP in the region
  ii                    <- which.max(res$log10p.soma)
  
  ## recolour high-scoring SNPs (assumed to be likely in LD)
  res$LD.col[ii]        <- "#D43F3AFF"
  
  ## load recombination rates
  rec <- data.frame(fread(paste0("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/compare_SomaLogic_Olink/genetic_map/genetic_map_GRCh37_chr", chr.s, ".txt"), sep="\t", header=T))
  ## subset to region of interest
  rec <- subset(rec, Position.bp. >= min(res$pos, na.rm=T) & Position.bp. <= max(res$pos, na.rm=T))
  

  plot(log10p.trait ~ log10p.soma, 
       xlim=c(0,max(res$log10p.soma)*1.05),
       data=res,
       pch=21, cex=.5, lwd=.3,
       bg=res$LD.col, 
       xlab=bquote(.(sl_n)~"GWAS"~~-log[10]("p-value")),
       ylab=bquote(.(t_n)~"GWAS"~~-log[10]("p-value"))
  )
  
  ## highlight top variant somalogic
  jj <- which.max(res$log10p.soma)
  points(res$log10p.soma[jj], res$log10p.trait[jj], pch=24, bg=res$LD.col[jj], cex=.6, lwd=.5)
  ## add text
  text(res$log10p.soma[jj], res$log10p.trait[jj], cex=.4, labels=res$rsid[jj], pos=4, offset = .2, xpd=NA)
  
  ## higlight top variant for trait
  jj <- which.max(res$log10p.trait)
  points(res$log10p.soma[jj], res$log10p.trait[jj], pch=24, bg=res$LD.col[jj], cex=.6, lwd=.5)
  ## add text
  text(res$log10p.soma[jj], res$log10p.trait[jj], cex=.4, labels=res$rsid[jj], pos=2, offset = .2, xpd=NA)
  
  ## add legend for LD
  legend("bottomright", lty=0, pch=22, pt.cex=1.5, pt.bg=ggsci::pal_locuszoom()(5), 
         pt.lwd=.5, title=expression(r^2), cex=.5, legend=rev(levels(res$LD_cat)), bty="n")
  
  ## add legend for Coloc
  legend(ifelse(sum.coloc[1,8] > .5, "topleft", "topright"), lty=0, pch=NA, cex=.5,
         legend=paste(paste0("H", 1:4), "=", sprintf("%.1f", sum.coloc[1,5:8]*100)),
         title="Posterior prob. [%]")
  
  #------------------------------------#
  ##--      simple locus plot       --##
  #------------------------------------#
  
  ## plot the SOMAscan
  par(mar=c(1.5,1.5,1,1.5))
  ## empty plot with recombination rate
  plot(range(res$pos), c(0,100), yaxt="n", ylab="", xlab="genomic position", type="n")
  
  ## add recombination rate
  points(rec$Position.bp., rec$Rate.cM.Mb, type="l", col="grey40", lwd=.2)
  axis(4)
  
  ## now the p-values
  par(new=T)
  plot(log10p.soma ~ pos, data=res, pch=21, cex=.5, bg=res$LD.col, lwd=.2,
       ylab=bquote(.(sl_n)~"GWAS"~~-log[10]("p-value")),
       xlab="genomic position")
  ## highlight SNP LD was based on
  points(res$pos[ii], res$log10p.soma[ii], pch=22, col="#D43F3AFF", lwd=.5)
  print(ii)
  text(res$pos[ii], res$log10p.soma[ii], cex=.4, labels=res$rsid[ii], pos=3, offset = .2, xpd=NA)
  
  
  ## add now trait
  ## empty plot with recombination rate
  plot(range(res$pos), c(0,100), yaxt="n", ylab="", xlab="genomic position", type="n")

  ## add recombination rate
  points(rec$Position.bp., rec$Rate.cM.Mb, type="l", col="grey70", lwd=.2)
  axis(4)
  ## now the p-values
  par(new=T)
  plot(log10p.trait ~ pos, data=res, pch=21, cex=.5, bg=res$LD.col, lwd=.2,
       ylab=bquote(.(t_n)~"GWAS"~~-log[10]("p-value")),
       xlab="genomic position")
  ## highlight SNP LD was based on
  points(res$pos[ii], res$log10p.trait[ii], pch=22, col="#D43F3AFF", lwd=.5)
  text(res$pos[ii], res$log10p.trait[ii], cex=.4, labels=res$rsid[ii], pos=3, offset = .2, xpd=NA)

  
}