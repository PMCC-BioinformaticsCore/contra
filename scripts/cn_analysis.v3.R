# ----------------------------------------------------------------------#
# Copyright (c) 2011, Richard Lupat & Jason Li.
#
# > Source License <
# This file is part of CONTRA.
#
#    CONTRA is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    CONTRA is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with CONTRA.  If not, see <http://www.gnu.org/licenses/>.
#
#
#-----------------------------------------------------------------------#
# Last Updated : 31 Oct 2011 17:00PM


# Parameters Parsing (from Command Line)
options <- commandArgs(trailingOnly = T)
bins = as.integer(options[1])
rd.cutoff = as.integer(options[2])
min.bases = as.integer(options[3])
outf = options[4]
sample.name = options[5]
plotoption = options[6]
actual.bin = as.numeric(options[7])
min_normal_rd_for_call = as.numeric(options[8])
min_tumour_rd_for_call = as.numeric(options[9])
min_avg_cov_for_call = as.numeric(options[10])

if (sample.name == "No-SampleName")
	sample.name = ""

if (sample.name != "")
	sample.name = paste(sample.name, ".", sep="")

# Setup output name
out.f = paste(outf, "/table/", sample.name, "CNATable.", rd.cutoff,"rd.", min.bases,"bases.", bins,"bins.txt", sep="")
pdf.out.f = paste(outf, "/plot/", sample.name, "densityplot.", bins, "bins.pdf", sep="")

# Open and read input files
# cnAverageFile = paste("bin", bins, ".txt", sep="")
cnAverageFile = paste(outf,"/buf/bin",bins,".txt",sep="")
boundariesFile = paste(outf,"/buf/bin",bins,".boundaries.txt",sep="")
#print (cnAverageFile)
cn.average = read.delim(cnAverageFile, as.is=F, header=F)
cn.boundary= read.delim(boundariesFile,as.is=F, header=F)

# Apply thresholds and data grouping
cn.average.aboveTs = cn.average[cn.average$V3>min.bases,]
cn.average.list = as.matrix(cn.average.aboveTs$V4)

# Get the mean and sd for each bins
cn.average.mean = c()
cn.average.sd = c()
cn.average.log= c()

# Density Plots for each bins
if (plotoption == "True"){
	pdf(pdf.out.f)
}
for (j in 1:actual.bin){
	cn.average.nth 	= as.matrix(cn.average.aboveTs[cn.average.aboveTs$V15==j,]$V4)
	cn.coverage.nth	= as.matrix(cn.average.aboveTs[cn.average.aboveTs$V15==j,]$V11)
	boundary.end   	= cn.boundary[cn.boundary$V1==j,]$V2
	boundary.start 	= cn.boundary[cn.boundary$V1==(j-1),]$V2
	boundary.mid	= (boundary.end+boundary.start)/2
	if (plotoption == "True") {
		plot_title = paste("density: bin", bins, sep="")
		#plot(density(cn.average.nth),xlim=c(-5,5), title=plot_title)
		plot(density(cn.average.nth),xlim=c(-5,5))
	}
	cn.average.mean = c(cn.average.mean, mean(cn.average.nth))
#	cn.average.sd 	= c(cn.average.sd, sd(cn.average.nth))
	cn.average.sd 	= c(cn.average.sd, apply(cn.average.nth,2,sd))
	#cn.average.log 	= c(cn.average.log, boundary.mid)
	cn.average.log	= c(cn.average.log, log(mean(cn.coverage.nth),2))
}
if (plotoption == "True"){
	dev.off()
}

# for point outside of boundaries
if (bins > 1) {
        boundary.first  = cn.boundary[cn.boundary$V1==0,]$V2
        boundary.last   = cn.boundary[cn.boundary$V1==bins,]$V2

        b.mean.y2       = cn.average.mean[2]
        b.mean.y1       = cn.average.mean[1]
        b.sd.y2         = cn.average.sd[2]
        b.sd.y1         = cn.average.sd[1]
        b.x2            = cn.average.log[2]
        b.x1            = cn.average.log[1]

        boundary.f.mean = (((b.mean.y2- b.mean.y1)/(b.x2-b.x1))*(boundary.first-b.x1))+b.mean.y1
        boundary.f.sd   = (((b.sd.y2- b.sd.y1)/(b.x2-b.x1))*(boundary.first-b.x1))+b.sd.y1

	if (boundary.f.sd < 0){
		boundary.f.sd = 0
	}

        b.mean.y2       = cn.average.mean[bins]
        b.mean.y1       = cn.average.mean[bins-1]
        b.sd.y2         = cn.average.sd[bins]
        b.sd.y1         = cn.average.sd[bins-1]
        b.x2            = cn.average.log[bins]
        b.x1            = cn.average.log[bins-1]

        boundary.l.mean = (((b.mean.y2- b.mean.y1)/(b.x2-b.x1))*(boundary.last-b.x1))+b.mean.y1
        boundary.l.sd   = (((b.sd.y2- b.sd.y1)/(b.x2-b.x1))*(boundary.last-b.x1))+b.sd.y1

        #cn.average.log         = c(boundary.first, cn.average.log, boundary.last)
        #cn.linear.mean         = c(boundary.f.mean, cn.average.mean, boundary.l.mean)
        #cn.linear.sd           = c(boundary.f.sd, cn.average.sd, boundary.l.sd)

        cn.average.log  = c(boundary.first, cn.average.log)
        cn.linear.mean  = c(boundary.f.mean, cn.average.mean)
        cn.linear.sd    = c(boundary.f.sd, cn.average.sd)

}

# Linear Interpolation
if (bins > 1 ){
        #print(cn.average.log)
        #print(cn.linear.mean)
        #print(cn.linear.sd)
        mean.fn  <- approxfun(cn.average.log, cn.linear.mean, rule=2)
        sd.fn    <- approxfun(cn.average.log, cn.linear.sd, rule=2)
}


# Put the data's details into matrices
ids 		= as.matrix(cn.average.aboveTs$V1)
exons 		= as.matrix(cn.average.aboveTs$V6)
exons.pos 	= as.matrix(cn.average.aboveTs$V5)
gs 		= as.matrix(cn.average.aboveTs$V2)
number.bases	= as.matrix(cn.average.aboveTs$V3)
mean		= as.matrix(cn.average.aboveTs$V4)
sd 		= as.matrix(cn.average.aboveTs$V7)
tumour.rd 	= as.matrix(cn.average.aboveTs$V8)
tumour.rd.ori	= as.matrix(cn.average.aboveTs$V10)
normal.rd	= as.matrix(cn.average.aboveTs$V9)
normal.rd.ori 	= as.matrix(cn.average.aboveTs$V11)
median		= as.matrix(cn.average.aboveTs$V12)
MinLogRatio 	= as.matrix(cn.average.aboveTs$V13)
MaxLogRatio 	= as.matrix(cn.average.aboveTs$V14)
Bin 		= as.matrix(cn.average.aboveTs$V15)
Chr 		= as.matrix(cn.average.aboveTs$V16)
OriStCoordinate = as.matrix(cn.average.aboveTs$V17)
OriEndCoordinate= as.matrix(cn.average.aboveTs$V18)

# Linear Fit
logratios.mean	= mean
logcov.mean	= log2((normal.rd + tumour.rd)/2)
fit.mean	= lm(logratios.mean ~ logcov.mean)
fit.x		= fit.mean$coefficient[1]
fit.y		= fit.mean$coefficient[2]

adjusted.lr	= rep(NA, length(logratios.mean))
for (j in 1:length(logratios.mean)){
	fitted.mean	= fit.x + fit.y * logcov.mean[j]
	adjusted.lr[j]	= logratios.mean[j] - fitted.mean
}

fit.mean2	= lm(adjusted.lr ~ logcov.mean)
fit.mean.a	= fit.mean2$coefficient[1]
fit.mean.b	= fit.mean2$coefficient[2]

fit.mean.fn <- function(x, fit.a, fit.b){
	result = fit.a + fit.b * x
	return (result)
}

# Adjust SD based on the new adjusted log ratios
logratios.sd	= c()
logcov.bins.mean= c()
for (j in 1:actual.bin){
	lr.bins.mean	= as.matrix(adjusted.lr[cn.average.aboveTs$V15==j])
#	logratios.sd	= c(logratios.sd, sd(lr.bins.mean))
	logratios.sd	= c(logratios.sd, apply(lr.bins.mean,2,sd))

	cn.coverage.tumour.nth = as.matrix(cn.average.aboveTs[cn.average.aboveTs$V15==j,]$V8)
	cn.coverage.normal.nth = as.matrix(cn.average.aboveTs[cn.average.aboveTs$V15==j,]$V9)
	cn.coverage.nth	= (cn.coverage.tumour.nth + cn.coverage.normal.nth) /2
	logcov.bins.mean= c(logcov.bins.mean, log2(mean(cn.coverage.nth)))

}

logratios.sd.ori = logratios.sd
if (length(logratios.sd) > 2) {
	logratios.sd 	= logratios.sd[-length(logratios.sd)]
}

logcov.bins.mean.ori = logcov.bins.mean
if (length(logcov.bins.mean) > 2){
	logcov.bins.mean= logcov.bins.mean[-length(logcov.bins.mean)]
}

fit.sd		= lm(log2(logratios.sd) ~ logcov.bins.mean)
fit.sd.a	= fit.sd$coefficient[1]
fit.sd.b	= fit.sd$coefficient[2]

fit.sd.fn <- function(x, fit.a, fit.b){
	result = 2 ^ (fit.mean.fn(x, fit.a, fit.b))
	return (result)
}

# Get the P Values, called the gain/loss
# with average and sd from each bins
pVal.list = c()
gain.loss = c()

for (i in 1:nrow(cn.average.list)){
	#print (i)
	#logratio = cn.average.list[i]
	#logcov	 = log(normal.rd.ori[i],2)
	logratio = adjusted.lr[i]
	logcov	 = logcov.mean[i]
	exon.bin = Bin[i]

	if (length(logratios.sd) > 1){
		#pVal <- pnorm(logratio, fit.mean.fn(logcov, fit.mean.a, fit.mean.b), fit.sd.fn(logcov, fit.sd.a, fit.sd.b))
		pVal <- pnorm(logratio, fit.mean.fn(logcov, fit.mean.a, fit.mean.b), sd.fn(logcov))
	} else {
		pVal <- pnorm(logratio, 0, logratios.sd[exon.bin])
	}

	if (pVal > 0.5){
		pVal = 1-pVal
		gain.loss = c(gain.loss, "gain")
	} else {
		gain.loss = c(gain.loss, "loss")
	}
	pVal.list = c(pVal.list, pVal*2)
}

# Get the adjusted P Values
adjusted.pVal.list = p.adjust(pVal.list, method="BH")

# Write the output into a tab-delimited text files
outdf=data.frame(Targeted.Region.ID=ids,Exon.Number=exons,Gene.Sym=gs,Chr, OriStCoordinate, OriEndCoordinate, Mean.of.LogRatio=cn.average.list, Adjusted.Mean.of.LogRatio=adjusted.lr, SD.of.LogRatio=sd, Median.of.LogRatio=median, number.bases, P.Value=pVal.list ,Adjusted.P.Value=adjusted.pVal.list , gain.loss, tumour.rd, normal.rd, tumour.rd.ori, normal.rd.ori, MinLogRatio, MaxLogRatio, BinNumber = Bin)

#min_normal_rd_for_call=5
#min_tumour_rd_for_call=0
#min_avg_cov_for_call=20
outdf$tumour.rd.ori = outdf$tumour.rd.ori-0.5
outdf$normal.rd.ori = outdf$normal.rd.ori-0.5
wh.to.excl = outdf$normal.rd.ori < min_normal_rd_for_call
wh.to.excl = wh.to.excl | outdf$tumour.rd.ori < min_tumour_rd_for_call
wh.to.excl = wh.to.excl | (outdf$tumour.rd.ori+outdf$normal.rd.ori)/2 < min_avg_cov_for_call
outdf$P.Value[wh.to.excl]=NA
outdf$Adjusted.P.Value[wh.to.excl]=NA


# rounding numbers to 3 digits
outdf$Mean.of.LogRatio <- round(outdf$Mean.of.LogRatio, digits=3)
outdf$Adjusted.Mean.of.LogRatio <- round(outdf$Adjusted.Mean.of.LogRatio, digits=3)
outdf$SD.of.LogRatio <- round(outdf$SD.of.LogRatio, digits=3)
outdf$Median.of.LogRatio <- round(outdf$Median.of.LogRatio, digits=3)

#I dont know if it is smart to round the pValue... but maybe like 6 digits instead?
outdf$P.Value <- round(outdf$P.Value, digits=6)
outdf$Adjusted.P.Value <- round(outdf$Adjusted.P.Value, digits=6)

write.table(outdf,out.f,sep="\t",quote=F,row.names=F,col.names=T)

#Plotting SD
#a.sd.fn  	= rep(fit.sd.a, length(logratios.sd.ori))
#b.sd.fn    	= rep(fit.sd.b, length(logratios.sd.ori))
#sd.after.fit 	= fit.sd.fn(logcov.bins.mean.ori, fit.sd.a, fit.sd.b)
#sd.out.f 	= paste(outf, "/plot/", sample.name, "sd.data_fit.", bins, "bins.txt", sep="")
#sd.outdf 	= data.frame(SD.Before.Fit = logratios.sd.ori, Log.Coverage = logcov.bins.mean.ori, SD.After.Fit = sd.after.fit, a.for.fitting=a.sd.fn, b.for.fitting=b.sd.fn)
#write.table(sd.outdf, sd.out.f,sep="\t", quote=F, row.names=F, col.names=T)


#End of the script
message("End of cn_analysis.R")
#print (i)
