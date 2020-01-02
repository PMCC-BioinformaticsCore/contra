#USAGE Rscript 3_wgsummary.R <scriptPath> <_wg.txt> <_exon.txt> <_wgSummary.txt> <_exSummary.txt> [debug T/F]

# Update 26/09 - changed output columns, _exSummary.txt also
# 28/10 - fixed pnorm to properly use SD not Var...

# command line input based on Jason Li's wgcnv_withBAFplot.R
x=commandArgs()
args=x[match("--args",x)+1:length(x)] # gets the arguments after --args...?

scriptPath = args[1]
f_out = args[4]
f_ex_out = args[5]
f_wg = args[2]
f_exon = args[3]
thresh_wg = args[6]
thresh_exon = args[7]
debug = args[8]

if (debug == "T") {
  print(f_out)
  print(f_ex_out)
  print(f_wg)
  print(f_exon)
  print(thresh_wg)
  print(thresh_exon)
}

# Set the working directory (for local testing)
# setwd("C:/Users/Emma/Desktop/UROP/MyWGCNV")

# Just setting up input and output
# samp_name = "1_12K0038"
# f_out = paste(samp_name, "_wgsummary.txt", sep="")
# f_wg = paste(samp_name, "_wholegene.txt", sep = "")
# f_exon = paste(samp_name, "_exon.txt", sep = "")
# thresh_wg = "thresh_wholegene.txt"
# thresh_exon = "thresh_exon.txt"


# Read in appropriate files. row.names = 1 -> use gene name as row name
data_wg = read.delim(f_wg, as.is = TRUE, row.names="Gene")
data_exon = read.delim(f_exon, as.is = TRUE)
data_thresh_wg = read.delim(thresh_wg, as.is = TRUE, row.names = 1)
data_thresh_exon = read.delim(thresh_exon, as.is = TRUE, row.names = 1)


# Confidence level to use throughout script
mu = 0.05


# Function that, given an exon name, will return a data frame with the following data:
# Exon, AvLR, All_AvLR,All_VarLR, pval, Gain/Loss, adjpval(added later), adjGainLoss (added later)
# Exon Name = gene_startcoord_endcoord
# Extract exon average previously (hard to do now...)
getExonStat = function(exon, ex_av) {
  ex_all_av = data_thresh_exon[exon,]$Cropped.Mean
  ex_all_var = data_thresh_exon[exon,]$Cropped.Var

  # Calculate the pval. If log ratio below average use lower tail, otherwise use upper
 # print(exon)
 # print(ex_av)
  #print(ex_all_av)
  if (ex_av < ex_all_av) {
    tail = TRUE
  } else {
    tail = FALSE
  }
  ex_pval = pnorm(ex_av, ex_all_av, sqrt(ex_all_var), lower.tail = tail)

  ex_gainloss = "."
  # Is it a gain or a loss? Using significance level of mu
  if (ex_pval < mu) {
    # Alright, now decide if its a gain or a loss
    if (ex_av < ex_all_av) {
      ex_gainloss = "loss"
    } else {
      ex_gainloss = "gain"
    }
  }

  ex_outdf = data.frame(Exon = exon, Av.LR = ex_av, Overall.Av.LR = ex_all_av, Overall.Var.LR = ex_all_var,
                        pval = ex_pval, GainLoss = ex_gainloss)

  return(ex_outdf)

}


# Function that, given a gene name, will return a data frame with the following data:
# Gene, AvLR, All_AvLR, All_VarLR, pval, adjpval (added later in process), n.exons, Gain/Loss
# TODO: n.exons.gain, n.exons.loss, n.exons.gain.adj, n.exons.loss.adj
getGeneStat = function(gene) {
  # Pull the first 3 straight out of the data!
  av = data_wg[gene,]
  all_av = data_thresh_wg[gene,]$Cropped.Mean
  all_var = data_thresh_wg[gene,]$Cropped.Var
 # print (av)
  #print (all_av)
  # Calculate the pval. If log ratio below average use lower tail, otherwise use upper
  if (av < all_av) {
    tail = TRUE
  } else {
    tail = FALSE
  }
  pval = pnorm(av, all_av, sqrt(all_var), lower.tail = tail)

  gainloss = "."
  # Is it a gain or a loss? Using significance level of mu
  if (pval < mu) {
    # Alright, now decide if its a gain or a loss
    if (av < all_av) {
      gainloss = "loss"
    } else {
      gainloss = "gain"
    }
  }


  # Now make it into a data frame and return it
  outdf = data.frame(Gene = gene, Av.LR = av, Overall.Av.LR = all_av, Overall.Var.LR = all_var,
                     pval = pval, GainLoss = gainloss)

  # Technically don't have to do this, but do it anyway
  return (outdf)
}

# Do exon processing (to f_ex_out)
app = FALSE
exons = row.names(data_exon)

#Now, we want to process all the genes
for (i in 1:length(exons)) {
  # For each exon
  ex_row = data_exon[i,]
  # Check exon is not NA
  if (is.na(ex_row$Mean.Log.Ratio) == FALSE){

    exon = paste(ex_row$Gene, "_",ex_row$Chr,":",ex_row$Start.Coord, "_",ex_row$End.Coord, sep = "")
#    print(exon)

    if ( ! is.na( data_thresh_exon[exon,"Cropped.Var"]) ) {
	    # Get the data
	    ex_tmpdf = getExonStat(exon, ex_row$Mean.Log.Ratio)
	    # Want to have column names and not append (ie overwrite) for the first gene
	    # But for the rest, want to have no column names (already have!) and append (so we dont lose data)
	    write.table(ex_tmpdf, f_ex_out, col.names = !app, append = app, sep = "\t",
			row.names = FALSE, quote = FALSE)
	    app = TRUE
    }
  }
}

# Now, we want to adjust the pvalues and then overwrite the table with new data
ex_outdf = read.delim(f_ex_out, as.is = TRUE)
ex_outdf$adj.pval = p.adjust(ex_outdf$pval, method = "BH")
write.table(ex_outdf, f_ex_out, col.names = TRUE, append = FALSE, sep = "\t",
            row.names = FALSE, quote = FALSE)
# Then go calculate the adj.GainLoss
ex_outdf = read.delim(f_ex_out, as.is = TRUE)
ex_outdf["adj.GainLoss"] <- '.'
for (i in 1:length(row.names(ex_outdf))) {
  if (ex_outdf[i,]$adj.pval < mu) {
    if (ex_outdf[i,]$Av.LR < ex_outdf[i,]$Overall.Av.LR) {
      ex_outdf[i,]$adj.GainLoss = "loss"
    } else {
      ex_outdf[i,]$adj.GainLoss = "gain"
    }
  }
}
write.table(ex_outdf, f_ex_out, col.names = TRUE, append = FALSE, sep = "\t",
            row.names = FALSE, quote = FALSE)




app = FALSE
genes = row.names(data_wg)
#Now, we want to process all the genes
for (i in 1:length(genes)) {
  # For each gene
  gene = genes[i]
  #print(gene)
  # Check gene actually exists
  #print(is.na(data_wg[gene,]))
  if(is.na(data_wg[gene,]) == FALSE) {

    if (! is.na(data_thresh_wg[gene,]$Cropped.Var) ) {
	    # Get the data
	    tmpdf = getGeneStat(gene)
	    # Want to have column names and not append (ie overwrite) for the first gene
	    # But for the rest, want to have no column names (already have!) and append (so we dont lose data)
	    write.table(tmpdf, f_out, col.names = !app, append = app, sep = "\t",
			row.names = FALSE, quote = FALSE)
	    app = TRUE
    }
  }
}

# Now, we want to adjust the pvalues and then overwrite the table with new data
outdf = read.delim(f_out, as.is = TRUE)
outdf$adj.pval = p.adjust(outdf$pval, method = "BH")
write.table(outdf, f_out, col.names = TRUE, append = FALSE, sep = "\t",
            row.names = FALSE, quote = FALSE)

# Now calculate the exon level information for each gene
# exons.sig.gain, exons.sig.loss, exons.sig.gain.adj, exons.sig.loss.adj, n.exons
# Also do adj.GainLoss here because WHY NOT
outdf = read.delim(f_out, as.is=TRUE)
ex_out = read.delim(f_ex_out, as.is = TRUE)
outdf[c("adj.GainLoss", "exons.sig.gain", "exons.sig.loss", "exons.sig.gain.adj", "exons.sig.loss.adj", "n.exons")] = ".";

for (i in 1:length(genes)) {
  # For each gene
  gene.exons.sig.gain = 0
  gene.exons.sig.loss = 0
  gene.exons.sig.gain.adj = 0
  gene.exons.sig.loss.adj = 0
  gene.n.exons = 0

print (paste("Processing gene ",i,"out of",length(genes)))

  # Go through each exon in f_ex_out
  for (j in 1:length(exons)) {
   # print(genes[i])

      tmp_ex_name = ex_outdf[j,]$Exon
  #print(tmp_ex_name)
  if (!is.na(tmp_ex_name)) { # why is this even a problem??
  if (tmp_ex_name == genes[i] || grepl(genes[i], tmp_ex_name)) {
      # An exon of the gene
      gene.n.exons = gene.n.exons + 1
      if (ex_outdf[j,]$GainLoss == "gain") {
        gene.exons.sig.gain = gene.exons.sig.gain + 1
      } else if (ex_outdf[j,]$GainLoss == "loss") {
        gene.exons.sig.loss = gene.exons.sig.loss + 1
      }
      if (ex_outdf[j,]$adj.GainLoss == "gain") {
        gene.exons.sig.gain.adj = gene.exons.sig.gain.adj + 1
      } else if (ex_outdf[j,]$adj.GainLoss == "loss") {
        gene.exons.sig.loss.adj = gene.exons.sig.loss.adj + 1
      }
    }
    }
  }


  outdf[i,]$exons.sig.gain = gene.exons.sig.gain
  outdf[i,]$exons.sig.loss = gene.exons.sig.loss
  outdf[i,]$exons.sig.gain.adj = gene.exons.sig.gain.adj
  outdf[i,]$exons.sig.loss.adj = gene.exons.sig.loss.adj
  outdf[i,]$n.exons = gene.n.exons

  # adj.GainLoss stuff here
  if (outdf[i,]$adj.pval < mu) {
    if (outdf[i,]$Av.LR < outdf[i,]$Overall.Av.LR) {
      outdf[i,]$adj.GainLoss = "loss"
    } else {
      outdf[i,]$adj.GainLoss = "gain"
    }
  }

}

write.table(outdf, f_out, col.names = TRUE, append = FALSE, sep = "\t",
            row.names = FALSE, quote = FALSE)




closeAllConnections()
