prom <- read.table("session_20180305 - interacting_genes - promoter info.txt",head=T,sep="\t",stringsAsFactors=F,comment.char="",quote="")
se <- read.table("20161115_corrected_ES_super_enhancer_with_methylation.txt",head=T,sep="\t",stringsAsFactors=F,comment.char="",quote="")

library(vioplot)

names <- c("2i","Serum","EpiSC")

png("Figure2_promoterMethylation.png",h=6,w=6,res=300,unit="in")
vioplot(na.omit(prom[,6]),na.omit(prom[,8]),na.omit(prom[,10]),col="lightgrey",names=names)
title(main="Promoters",ylab="mCpG/CpG")
axis(1, at=c(1.5),labels=c("ESC"), padj=1.5,tick=F)
dev.off()

png("Figure2_enhancerMethylation.png",h=6,w=6,res=300,unit="in")
vioplot(na.omit(se[,9]),na.omit(se[,15]),na.omit(se[,27]),col="lightgrey",names=names)
title(main="Super Enhancers",ylab="mCpG/CpG")
axis(1, at=c(1.5),labels=c("ESC"), padj=1.5,tick=F)
dev.off()

