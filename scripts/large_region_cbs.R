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
# Last Updated : 12 October 2011 16:43PM

library(DNAcopy)
options <- commandArgs(trailingOnly = T)
logcov.file = options[1]
param.small.seg = as.numeric(options[2])
param.large.seg = as.numeric(options[3])
param.pval.cutoff = as.numeric(options[4])
param.pvals.size  = as.numeric(options[5])
param.call.low	= as.numeric(options[6])
param.call.high	= as.numeric(options[7])
bufloc		= options[8]

#param.small.seg = 1
#param.large.seg = 25
#param.pval.cutoff = 0.1
#param.pvals.size  = 0.35

#param.call.low	= -0.3
#param.call.high= 0.3


dat = read.delim(logcov.file, as.is=F, header=T)
t.id = dat$Targeted.Region.ID
t.mean = dat$Adjusted.Mean.of.LogRatio
t.chr = dat$Chr

cna.dat <- CNA(t.mean, t.chr, t.id, data.type="logratio")
smooth.cna.dat = smooth.CNA(cna.dat)

kmax.range = c(param.large.seg, param.small.seg)
for (i in kmax.range){
	out.f = paste(bufloc, "/CBS", i ,".txt", sep="")
	out.f2 = paste(logcov.file, ".CBS_", i ,".txt", sep="")
	cna.segment = segment(smooth.cna.dat, kmax = i, verbose = 1)
	#pdf(paste("segment_", "kmax", i, ".pdf", sep=""))
	#for (chr.list in unique(t.chr)){
	#	plotSample(cna.segment, chromlist=chr.list)
	#}
	#dev.off()

	#Get Data
	cna.segment.out <- cna.segment$output
	cna.segment.start = c()
	cna.segment.end	  = c()
	cna.segment.mean  = c()
	cna.segment.pvals.size = c()
	cna.segment.calls = c()
	for (n in 1:nrow(cna.segment.out)){
		cna.start.id	= cna.segment.out[n,]$loc.start
		cna.end.id	= cna.segment.out[n,]$loc.end
		cna.start.coord = dat[dat$Targeted.Region.ID==(cna.start.id),]$OriStCoordinate
		cna.end.coord	= dat[dat$Targeted.Region.ID==(cna.end.id), ]$OriEndCoordinate

		dat.inrange 		= dat[dat$Targeted.Region.ID<=cna.end.id&dat$Targeted.Region.ID>=cna.start.id,]
		cna.segment.logratios 	= dat.inrange$Adjusted.Mean.of.LogRatio
		cna.segment.pvalues   	= dat.inrange$Adjusted.P.Value
		segment.pvals.above	= dat.inrange[dat.inrange$Adjusted.P.Value<=param.pval.cutoff,]$Adjusted.P.Value

		segment.pvals.size	= length(segment.pvals.above)/length(cna.segment.pvalues)
		cna.segment.mean  = c(cna.segment.mean, mean(cna.segment.logratios))
		cna.segment.start = c(cna.segment.start, cna.start.coord)
		cna.segment.end	  = c(cna.segment.end , cna.end.coord)
		cna.segment.pvals.size = c(cna.segment.pvals.size, segment.pvals.size)

		if (segment.pvals.size < param.pvals.size){
			#No-Call
			cna.segment.calls = c(cna.segment.calls, "No")
		} else {
			m = mean(cna.segment.logratios)
			if ((m > param.call.low) && (m < param.call.high)){
				#No-Call
				cna.segment.calls = c(cna.segment.calls, "No")
			} else {
				#Call
				cna.segment.calls = c(cna.segment.calls, "CNV")
			}
		}
	}


	outdf = data.frame(Chr=cna.segment.out$chrom, Target.Start=cna.segment.out$loc.start, Target.End=cna.segment.out$loc.end, NumberOfTargets=cna.segment.out$num.mark, OriStCoordinate=cna.segment.start, OriEndCoordinate=cna.segment.end, CBS.Mean = cna.segment.out$seg.mean, LogRatios = cna.segment.mean, Above.PValues.Cutoff=cna.segment.pvals.size, Calls = cna.segment.calls)

	write.table(outdf, out.f, sep="\t", quote=F, row.names=F, col.names=T)
	write.table(outdf, out.f2, sep="\t", quote=F, row.names=F, col.names=T)

}

small.dat.f = paste(bufloc, "/CBS", param.small.seg ,".txt", sep="")
large.dat.f = paste(bufloc, "/CBS", param.large.seg ,".txt", sep="")

small.dat = read.delim(small.dat.f, as.is=F, header=T)
large.dat = read.delim(large.dat.f, as.is=F, header=T)

large.nocall    = large.dat[large.dat$Calls=="No",]
small.call      = small.dat[small.dat$Calls!="No",]
included.segment= large.dat[large.dat$Calls!="No",]

for (i in 1:nrow(small.call)){
        small.segment   = small.call[i, ]
        chr.small       = small.segment$Chr
        start.small     = small.segment$Target.Start
        end.small       = small.segment$Target.End
        match.large.nocall = large.nocall[large.nocall$Chr==chr.small,]

        match.large.nocall = match.large.nocall[match.large.nocall$Target.End>=start.small,]
        match.large.nocall = match.large.nocall[match.large.nocall$Target.Start<=start.small,]
        merged.hole        = match.large.nocall[match.large.nocall$Target.End>=start.small,]

        if (nrow(merged.hole) > 0){
                included.segment = merge(included.segment, small.segment, all=T)
        }
}

out.f = paste(logcov.file, ".LargeDeletion.txt", sep="")
write.table(included.segment, out.f, sep="\t", quote=F, row.names=F)
