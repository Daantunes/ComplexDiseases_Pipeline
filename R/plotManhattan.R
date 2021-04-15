
#
# File: plotManhattan.R
# Created Date: Tuesday March 3rd 2020
# Author: Debora Antunes
# -----
# Last Modified: Thursday, July 9th 2020, 7:19:34 pm
# -----
#

fun = function(var, name){
	#################################################################
	### FUNCTION THAT CREATES A MANHATTAN PLOT
	### var: a dataset with three columns: chromosome number, position of SNP and value
	##################################################################

	### Necessary packages  ###
	#install.packages("qqman")
	library(qqman)
	
	chr = sort(unique(var$"chrNo"))

	if (length(chr) == 24) { 
	labs = c(1:22, "X", "Y")
	} else if (length(chr) == 23 & tail(chr, n=1) == 24) {
	labs = c(1:22, "Y")
	} else if (length(chr) == 23 & tail(chr, n=1) == 23) {
	labs = c(1:22, "X")
	} else {
	labs = c(1:22)
	}

	if (name == '_chi2') {
	limit=c(0,50)
	} else {
	limit=c(0,80)
	}

	print('>>> Saving figure...')
	path=sprintf("../../data/figures/manhattanPlot%s.png", name)
	png(filename=path, width = 1000, height = 500, units = "px")
	manhattan(var, 
				chr = "chrNo", 
				bp = "posNo", 
				p = "val",
				chrlabs = labs,
				ylim=limit, 
				suggestiveline = -log10(0.05), 
				genomewideline = FALSE)
	dev.off()
	}