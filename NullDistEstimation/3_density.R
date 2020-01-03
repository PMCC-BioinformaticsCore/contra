#USAGE Rscript 3_wgsummary.R <summary.txt> <out.txt> [debug T/F] [histograms T/F] [(histPath)]
# Last modified 8/10/13

### LAST MODIFIED 25 11 2004 --- NOTE  value "3/4" should be taken from input parameter (proportion of guaranteed germline samples)

# command line input based on Jason Li's wgcnv_withBAFplot.R
x=commandArgs()
args=x[match("--args",x)+1:length(x)] # gets the arguments after --args...?

summaryin= args[1]
threshout= args[2]
debug = args[3]
hists = args[4]

#DEBUG
#setwd("C:/Jason Li 2013/Current Projects/SUPER/NDE_debug")
#summaryin="summary_ex.txt"
#threshout= "thresh_ex_out.txt"
#debug="T"
#hists="T"


if (hists == "T") {
#DEBUG
#histPath = "./hist2.pdf"

  histPath = args[5]
  filename = histPath
  pdf(file = filename)
}



data = read.delim(summaryin, row.names = 1,as.is=T)
rownames = row.names(data)

i = 1
means = c()
vars = c()
left_threshs = c()
right_threshs = c()
cropped_means = c()
cropped_vars = c()

for (i in 1:nrow(data)) {
	###############
print(i)
	row <- data[i,]
	rowname = rownames[i]
	numrow = as.numeric(row)

	# set centre as mean of data
	# trim mean (remove tails)
	mean = mean(numrow, na.rm = TRUE, trim = 0.3)
	sd1=sd(numrow,na.rm=T)

	best_p=1
	best_thresh=NA
	if (length(numrow[!is.na(numrow)]) >= 3) {

		# adjust symmetrically
		# from .75 to 3.5
		#for (i in 5:175) {
		#	thresh = i/50

		n.steps=100
		#num_range = max(numrow) - min(numrow)
		#num_range_step = num_range/100

		#num_range_step = best_thresh / n.steps

		num_range_step = sd1*2*2 / n.steps
		for (i2 in 1:n.steps) {
			thresh = i2 * num_range_step
			#print(paste((mean - thresh), (mean + thresh)))
		  croprow = numrow[numrow >= (mean - thresh) & numrow <= (mean + thresh)]
		  croprow = croprow[!is.na(croprow)]
		  #if (length(croprow) > 3) {
		  if (length(croprow) >  (length(numrow[!is.na(numrow)]) * 3/4) && length(croprow) >= 3) {
		  #	print(croprow)
			test = shapiro.test(croprow)
		#	if (test$p.value > best_p) {
		#		best_p = test$p.value
		#		best_thresh = thresh
		#	}
		#print (test$p.value)
			if (test$p.value < best_p) {
				best_p = test$p.value
				best_thresh = thresh
			}
		  }
		}
		if (is.na(best_thresh)) {
			best_p = shapiro.test(numrow[!is.na(numrow)])$p.value
			best_thresh = max ( max(numrow) - mean, mean - min(numrow))
		}
	}

	if (hists == 'T' && length(numrow[!is.na(numrow)]) > 1) {
	  #x11()
	  plot(density(numrow, na.rm = TRUE),
		main = rowname,
		sub = paste("Mean: ", round(mean,4), " || Threshholds:", round(mean-best_thresh,4), round(mean+best_thresh,4)),
		xlab = paste("Variance:", round(var(numrow, na.rm = TRUE),4),
			" || Cropped variance: ", round(var(numrow[numrow > mean-best_thresh & numrow < mean+best_thresh], na.rm = TRUE),4)
			)
	  )

	  abline(v = mean+best_thresh, col = "red")
	  abline(v = mean-best_thresh, col = "red")
	  abline(v = mean, col = "blue")
	  x = seq(-4,4,length=1000)
	  cropped_sd = sqrt(var(numrow[numrow > mean - best_thresh & numrow < mean + best_thresh], na.rm = TRUE))
	  y = dnorm(x, mean = mean, sd = cropped_sd)
	 # What if we dont crop the  variance?
	 # y = dnorm(x, mean = mean, sd = sqrt(var(numrow, na.rm = TRUE)))
	  #y2 = dnorm(x, mean = mean(numrow[numrow > mean-best_thresh & numrow < mean+best_thresh]), sd = cropped_sd)
	  lines(x,y,type="l", col ="blue")
	  #lines(x,y2,type="l", col = "green")
	  # dev.copy(png, filename)
	}
	# add stuff to vectors
	means <- c(means, mean)
	cropped_means <- c(cropped_means, mean(numrow[numrow > mean-best_thresh & numrow < mean+best_thresh], na.rm = TRUE))
	vars <- c(vars, var(numrow, na.rm = TRUE))
	cropped_vars <- c(cropped_vars, var(numrow[numrow > mean-best_thresh & numrow < mean+best_thresh], na.rm = TRUE))
	left_threshs <- c(left_threshs, mean-best_thresh)
	right_threshs <-c(right_threshs, mean+best_thresh)

}



dataout = data.frame(rownames, means, vars, left_threshs, right_threshs, cropped_means, cropped_vars)

colnames(dataout) <- c("Gene Name", "Mean(trimmed 0.2)", "Variance", "Left Thresh", "Right Thresh",
					"Cropped Mean", "Cropped Var")

write.table(dataout, threshout, quote = FALSE, sep = '\t',
	row.names = FALSE, col.names = TRUE)

if (hists == 'T') {
  dev.off()
}

if (debug == 'T'){
  print(paste("Finished: ", threshout, sep = " "))
}



closeAllConnections()
