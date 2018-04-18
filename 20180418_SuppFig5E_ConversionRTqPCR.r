setwd("Belletal2018/SuppFig5")

cnrq <- read.table("Suppfig5b.txt",sep="\t",head=T,comment.char="",quote="",stringsAsFactors=F)
sumtab <- read.table("Suppfig5b_means.txt",sep="\t",head=T,comment.char="",quote="",stringsAsFactors=F)

targets <- unique(cnrq$Target)
clones <- unique(cnrq$Clone)

conditions <- c("2i","Serum")

for(i in 1:length(targets)){
	
	target <- targets[i]
	
	toplot <- cnrq[which(cnrq$Target==target),]
	average <- sumtab[which(sumtab$Target==target),]
	
	xlim <-c(0,10)
	
	for(j in 1:length(conditions)){
		
		condition <- conditions[j]
		
		ylim <- c(0,max(toplot$CNRQ[which(toplot$Condition==condition)],na.rm=T)*1.1)
		main <- paste(condition,target,sep=" ")
		
		pngName <- paste0(condition,"_",target,"_RelExp.png")
		
		png(pngName,w=6,h=6,unit="in",res=300)
		lp <- plot(CNRQ~Day,data=toplot[which(toplot$Condition==condition & toplot$Clone==clones[1]),],xlim=xlim,ylim=ylim,type="o",pch=20,col="black",main=main,names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Day",ylab="Relative expression",lty=2)
		grid(nx=NA,ny=NULL)
		
		for(k in 2:length(clones)-1){
			points(CNRQ~Day,data=toplot[which(toplot$Condition==condition & toplot$Clone==clones[k]),],type="o",pch=20,col="black",names.arg=NULL,xaxt = 'n',yaxt = 'n',lty=2)
		}
		
		points(RelExp~Day,data=average[which(average$Condition==condition),],type="o",pch=16,names.arg=NULL,xaxt = 'n',yaxt = 'n',col="red",lwd=2)
		axis(1, at=c(0,1,2,3,7,10), labels=c(0,1,2,3,7,"cEpiSC"),tick=T, las=1, cex.axis=1)
		axis(2, at=axTicks(2), cex.axis=0.9, las=2)
		dev.off()
	}
}
