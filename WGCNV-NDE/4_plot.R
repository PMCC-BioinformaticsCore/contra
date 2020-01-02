#USAGE Rscript 4_plot.R <exSummary.f> <wgSummary.f> <sampLabel> <species> <species.chrlen.f> <out.f>
#setwd("W:/Current_Projects/SUPER/NDE_NIMB/ContraOut_TEST/1003-Tumour/WGCNV")
#ex.dat=read.delim("1003-Tumour_exSummary.txt",as.is=T)
#g.dat=read.delim("1003-Tumour_wgSummary.txt",as.is=T)

x=commandArgs()
args=x[match("--args",x)+1:length(x)]

#ex.f="14K0038_exSummary.txt"
#g.f="14K0038_wgSummary.txt"
#samp.lbl="14K0038"
#species="human"
#out.f="plot.pdf"
ex.f=args[1]
g.f=args[2]
samp.lbl=args[3]
species=args[4]
species.chrlen.f=args[5]
out.f=args[6]

plot.max.LR=4


ex.dat=read.delim(ex.f,as.is=T)
g.dat=read.delim(g.f,as.is=T)


if (species == "human") {


#	chrlens=read.delim("hg19_chrlen.txt",as.is=T,header=F)
	chrlens=read.delim(species.chrlen.f,as.is=T,header=F)
	genome.len=sum(as.numeric(chrlens$V2))


	chrlens$chrnum=as.numeric(sub("Y","24",sub("X","23",chrlens$V1)))



	offsetdf=data.frame(chr=1:24,offset=0,stringsAsFactors=F)
	tmp=0
	for (i in 1:nrow(offsetdf)) {
		offsetdf$offset[i]=tmp
		tmp=tmp+chrlens$V2[i]
	}


	rownames(offsetdf)=offsetdf$chr

	addGenomeCoord=function(dat,chr.col,pos.col,pos.offset=0) {
		dat$chrnum=as.numeric(sub("Y","24",sub("X","23",dat[,chr.col])))
		dat$genome.coord=dat[,pos.col]+pos.offset+offsetdf[dat$chrnum,"offset"]
		return (dat)
	}

} else {

	stop("species not implemented")
}


#ex.dat.coord = sub("^[^_]+_","",ex.dat$Exon,perl=T)
#ex.dat$Chromosome=sub(":.*$","",ex.dat.coord)
ex.dat$Chromosome=sub("^.*_","",sub(":.*$","",ex.dat$Exon))

ex.dat.coord = sub("^.*:","",ex.dat$Exon,perl=T)
ex.dat$Start=as.numeric(sub("^.*:","",sub("_.*$","",ex.dat.coord)))
ex.dat$End=as.numeric(sub("^.*_","",ex.dat.coord))
ex.dat$ExonSize=ex.dat$End-ex.dat$Start
wh.ok=is.element(ex.dat$Chromosome,chrlens[,1])
ex.dat2=ex.dat[wh.ok,]

ex.dat3=addGenomeCoord(ex.dat2,"Chromosome","Start",pos.offset=0)
ex.dat3$genome.end=ex.dat3$genome.coord+ex.dat$ExonSize


#seg.dat2=addGenomeCoord(seg.dat,"V1","V2",pos.offset=0)
#seg.dat2$genome.end=seg.dat2$genome.coord+(seg.dat2$V3-seg.dat2$V2)

y.min=-plot.max.LR
y.max=plot.max.LR
preparePlot=function() {
	lbl=samp.lbl

	plot(0,type="n",xlim=c(0,max(ex.dat3$genome.coord)),ylim=c(y.min,y.max),xlab="Genome position",ylab="Log Ratio",main=lbl,xaxt="n")#,xaxs="i")
	plot.chr.lbls.pos=chrlens$V2/2+offsetdf$offset
	plot.chr.lbls=paste("chr",chrlens$V1,sep="")
	axis(side=1,at=plot.chr.lbls.pos,labels=plot.chr.lbls,tick=F,las=3,cex.axis=0.6,line=-0.5)
	abline(v=offsetdf$offset[-1],col="grey88")
}


pdf(out.f,width=10,height=7)
preparePlot()

xs=(ex.dat3$genome.coord+ex.dat3$genome.end)/2
ys=ex.dat3$Av.LR
ys2=ys
ys2[ys>y.max]=y.max
ys2[ys<y.min]=y.min

points(xs,ys2,pch=20,col=rgb(0.2,0.2,0.2,0.1))

wh.outlie.hi=ys>y.max
points(xs[wh.outlie.hi],rep(y.max,sum(wh.outlie.hi)), pch="X",col=rgb(0.2,0.2,0.2,0.9))
wh.outlie.low=ys<y.min
points(xs[wh.outlie.low],rep(y.min,sum(wh.outlie.low)), pch="X",col=rgb(0.2,0.2,0.2,0.9))


#===== EXON
gain.col=rgb(0.6,0.1,0.1,0.3)
wh.gain=ex.dat3$adj.GainLoss=="gain"
xs.gain=xs[wh.gain]
ys.gain=ys2[wh.gain]
points(xs.gain,ys.gain,pch=20,col=gain.col)

loss.col=rgb(0.1,0.6,0.1,0.3)
wh.loss=ex.dat3$adj.GainLoss=="loss"
xs.loss=xs[wh.loss]
ys.loss=ys2[wh.loss]
points(xs.loss,ys.loss,pch=20,col=loss.col)


#===== gene
x.wd = (ex.dat3[nrow(ex.dat3),"genome.end"] - ex.dat3[1,"genome.coord"]) / 200
y.wd = plot.max.LR / 100


gene.gain=g.dat$Gene[g.dat$adj.GainLoss=="gain"]
exon.is.gain=rep(F,nrow(ex.dat3))
for ( g in gene.gain) {
#g="BRAF"
	p1=paste("\\b",g,"(\\b|_)",sep="")
	g.wh=grep(p1,ex.dat3$Exon)
	if (length(g.wh) > 0) {
		dat.ext=ex.dat3[g.wh,]
		LRs=dat.ext$Av.LR
		avg.LR=mean(LRs)
		x.left=dat.ext[1,"genome.coord"]
		x.right=dat.ext[nrow(dat.ext),"genome.end"]
		rect(x.left-x.wd,avg.LR-y.wd,x.right+x.wd,avg.LR+y.wd,col="red",border="black")
		text(x.left,plot.max.LR,g,cex=0.6,col="red",adj=1,srt=90)
	} else {
		print (paste("WARNING: gene not found in data",g))
	}
}

gene.loss=g.dat$Gene[g.dat$adj.GainLoss=="loss"]
exon.is.loss=rep(F,nrow(ex.dat3))
for ( g in gene.loss) {
#g="PTEN"
	p1=paste("\\b",g,"(\\b|_)",sep="")
	g.wh=grep(p1,ex.dat3$Exon)
	if (length(g.wh) > 0) {
		dat.ext=ex.dat3[g.wh,]
		LRs=dat.ext$Av.LR
		avg.LR=mean(LRs)
		x.left=dat.ext[1,"genome.coord"]
		x.right=dat.ext[nrow(dat.ext),"genome.end"]
		rect(x.left-x.wd,avg.LR-y.wd,x.right+x.wd,avg.LR+y.wd,col="green2",border="black")
		text(x.left,-plot.max.LR,g,cex=0.6,col="green2",adj=0,srt=90)
	} else {
		print (paste("WARNING: gene not found in data",g))
	}
}

dev.off()
