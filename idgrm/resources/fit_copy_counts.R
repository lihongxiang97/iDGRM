#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))
a <- commandArgs(trailingOnly=TRUE); stopifnot(length(a)==8)
counts <- as.matrix(read.delim(a[[1]],row.names=1,check.names=FALSE))
pairs <- read.delim(a[[2]],check.names=FALSE)
sm <- read.delim(a[[3]],check.names=FALSE,colClasses="character")
out <- a[[4]]; test <- a[[5]]; alpha <- as.numeric(a[[8]])
writeLines(capture.output(sessionInfo()),file.path(out,"R_sessionInfo.txt"))
base <- DESeqDataSetFromMatrix(counts,DataFrame(row.names=colnames(counts)),design=~1)
base <- estimateSizeFactors(base,type="poscounts"); sf <- sizeFactors(base)
write.table(data.frame(sample_id=names(sf),size_factor=unname(sf)),file.path(out,"size_factors.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
if(a[[6]]!="NONE") {
  one <- as.matrix(read.delim(a[[6]],row.names=1,check.names=FALSE))
  two <- as.matrix(read.delim(a[[7]],row.names=1,check.names=FALSE))
} else {
  one <- counts[pairs$gene1,,drop=FALSE];two <- counts[pairs$gene2,,drop=FALSE]
  rownames(one) <- pairs$pair_id;rownames(two) <- pairs$pair_id
}
for(tissue in unique(sm$tissue)) {
  ss <- sm[sm$tissue==tissue,,drop=FALSE]; n <- nrow(ss)
  mat <- cbind(one[pairs$pair_id,ss$sample_id,drop=FALSE],two[pairs$pair_id,ss$sample_id,drop=FALSE])
  rownames(mat) <- pairs$pair_id
  colnames(mat) <- c(paste0("copy1_r",seq_len(n)),paste0("copy2_r",seq_len(n)))
  cd <- DataFrame(biological_library=factor(rep(seq_len(n),2)),
    copy=factor(rep(c("copy1","copy2"),each=n),levels=c("copy1","copy2")),row.names=colnames(mat))
  dds <- DESeqDataSetFromMatrix(mat,cd,design=~biological_library+copy)
  sizeFactors(dds) <- rep(unname(sf[ss$sample_id]),2)
  keep <- rowSums(mat)>=10
  tab <- data.frame(pair_id=pairs$pair_id,baseMean=NA_real_,log2FoldChange=NA_real_,lfcSE=NA_real_,
    stat=NA_real_,pvalue=NA_real_,padj=NA_real_)
  if(sum(keep)>=2) {
    if(test=="wald") fit <- DESeq(dds[keep,],test="Wald",quiet=TRUE,minReplicatesForReplace=Inf)
    else fit <- DESeq(dds[keep,],test="LRT",reduced=~biological_library,quiet=TRUE,minReplicatesForReplace=Inf)
    rr <- as.data.frame(results(fit,contrast=c("copy","copy2","copy1"),alpha=alpha,independentFiltering=FALSE))
    idx <- match(rownames(rr),tab$pair_id)
    for(cc in intersect(names(rr),names(tab)))tab[idx,cc] <- rr[[cc]]
  }
  write.csv(tab,file.path(out,paste0("idgrm.counts.pairs.",tissue,".DESeq2.csv")),row.names=FALSE,quote=FALSE,na="")
}
writeLines("Count fitting complete; final-category FDR not implied",file.path(out,"COMPLETE"))
