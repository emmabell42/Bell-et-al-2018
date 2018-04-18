###################################################################
#
## Analysing the processed in vitro and in vivo RNA-seq from AJ
#
###################################################################
#
setwd("Belletal2018")

library(vioplot) #Load customised vioplot function for multiple colours

expr <- read.table("corrected SE_interacting_genes_expression_vitro_vivo_20180308.txt",sep="\t",head=T,comment.char="",quote="",stringsAsFactors=F)

expr$SE.type[which(expr$SE.type=="PM only")] <- NA

# Offset by 1 to allow log transformation
expr[,5:ncol(expr)] <- expr[,5:ncol(expr)]+1
expr.log <- expr[,5:ncol(expr)]
expr.log <- log(expr.log)

vitro <- c("2i","Serum","EpiSC")
vivo <- c("E3.5","E4.0","E5.5","E6.5")
models <- c("vitro","vivo")

maintained <- which(expr$SE.type=="maintained")
inactivated <- which(expr$SE.type=="inactivated")

se.type <- c("maintained","inactivated")

cols <- c(rgb(128,255,128,maxColorValue=255),rgb(255,128,255,maxColorValue=255))

for(i in 1:length(se.type)){
	type <- se.type[i]
	
	for(j in 1:length(models)){
		
		model <- models[j]
		
				
		pngName <- paste0("Vioplot_",model,"_",type,".png")
		main <- toSentence(type)
		names <- get(model)
		index <- which(expr$SE.type==type)
		
		p <- ""
		
		if(type=="inactivated"){
			
			if(model=="vitro"){
				t <- t.test(expr.log[index,2],expr.log[index,3])
			}
			else {
				t <- t.test(expr.log[index,5],expr.log[index,7])
			}
			p <- t[[3]]
			p <- format(p,digits=3)
			p <- paste0("P = ",p)
		}
		
		png(pngName,w=6,h=6,unit="in",res=300)
		
		if(model=="vitro"){
			vioplot(expr.log[index,1],expr.log[index,2],expr.log[index,3],names=names,col=cols[i])
			axis(1, at=c(1.5),labels=c("ESC"), padj=1.5,tick=F)
		} else {
			vioplot(na.omit(expr.log[index,4]),na.omit(expr.log[index,5]),na.omit(expr.log[index,6]),na.omit(expr.log[index,7]),names=names,col=cols[i])
			axis(1, at=c(1.5, 3.5),labels=c("ICM","Epiblast"), padj=1.5,tick=F)			
		}
		title(main=main,ylab=expression('Log'[10]*'(RPKM)'),xlab=p)
		dev.off()
}
}