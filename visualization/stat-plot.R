statPlot <- function(missingness_file, psc_file){
  
  library(ggplot2)
  library(readr)
  
  # BCFtools Per Sample Counts
  psc<-read_tsv(psc_file, col_names = FALSE)
  names(psc) <- c('psc','id','sample','nRefHom','nNonRefHom', 'nHets', 'nTransitions', 'nTransversions', 'nIndels', 'average_depth', 'nSingletons', 'nHapRef', 'nHapAlt' )
  
  # Transition / Transversion ratio
  
  tstv <- psc$nTransitions/psc$nTransversions
  psc <- cbind(psc, tstv)
  flag_tstv = as.data.frame(ifelse((psc$tstv < (mean(psc$tstv)- 3*sd(psc$tstv)) |psc$tstv > (mean(psc$tstv)+ 3*sd(psc$tstv)) ), 1,NA))
  if (length(flag_tstv)==0){flag_tstv=NA}
  
  
  # Total singletons per sample
  
  singletons <- psc$nSingletons
  psc <- cbind(psc, log10(singletons))
  flag_singletons = as.data.frame(ifelse(psc$`log10(singletons)` > (mean(psc$`log10(singletons)`)+ 3*sd(psc$`log10(singletons)`)), 1,NA))
  if (length(flag_singletons)==0){flag_singletons=NA}
  
  # Total indels per sample
  
  indels <- psc$nIndels
  flag_indels = as.data.frame(ifelse((psc$nIndels < (mean(psc$nIndels)- 3*sd(psc$nIndels)) |psc$nIndels > (mean(psc$nIndels)+ 3*sd(psc$nIndels)) ), 1,NA))
  if (length(flag_indels)==0){flag_indels=NA}
  
  # Total SNVs per sample
  
  snv <- psc$nHets + psc$nNonRefHom
  psc <- cbind(psc, snv)
  flag_snvs = as.data.frame(ifelse((psc$snv < (mean(psc$snv)- 3*sd(psc$snv)) |psc$snv > (mean(psc$snv)+ 3*sd(psc$snv)) ), 1,NA))
  if (length(flag_snvs)==0){flag_snvs=NA}
  
  # Sample depth of coverage 
  
  doc <- psc$average_depth
  flag_doc = as.data.frame(ifelse((psc$average_depth < (mean(psc$average_depth)- 3*sd(psc$average_depth)) |psc$average_depth > (mean(psc$average_depth)+ 3*sd(psc$average_depth)) ), 1,NA))
  if (length(flag_doc)==0){flag_doc=NA}
 
  # Sample Missingness (PLINK)
  
  missingness <- read_table(missingness_file)
  sample_index <- 1:nrow(missingness)
  c <- as.data.frame(cbind(missingness$IID, missingness$F_MISS))
  b <- as.data.frame(cbind(sample_index, missingness$F_MISS))
  flag_miss = as.data.frame(ifelse((b$V2 < (mean(b$V2)- 3*sd(b$V2)) |b$V2 > (mean(b$V2)+ 3*sd(b$V2)) ), 1,NA))
  if (length(flag_miss)==0){flag_miss=NA}
  
  flag <- (cbind(psc$sample, flag_doc, flag_singletons, flag_snvs, flag_indels, flag_miss, flag_tstv))
  names(flag) <- c('sample', 'tstv', 'doc', 'singletons', 'snvs', 'indels', 'miss')
  tstv <-subset(flag, (!is.na(tstv)))[1]
  doc <-subset(flag, (!is.na(doc)))[1]
  singletons <-subset(flag, (!is.na(singletons)))[1]
  snvs <-subset(flag, (!is.na(snvs)))[1]
  indels <-subset(flag, (!is.na(indels)))[1]
  miss <-subset(flag, (!is.na(miss)))[1]

  filter_out <- unique(rbind(tstv, doc, singletons, snvs, indels, miss))
  filter_out <- as.data.frame(filter_out)
  all_samples <- as.data.frame(psc$sample)
  names(all_samples) <- c("sample")
  samples_to_keep <- as.data.frame(setdiff(all_samples, filter_out))

  
  return(samples_to_keep)
  
}


