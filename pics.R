library(tidyverse)


##Site quality###################################################
var_qual <- read_delim("./t2d.lqual", delim = "\t",
           col_names = c("chr", "pos", "qual"), skip = 1)

png(filename="./t2d_total_lqual.png")
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "#a5a5a5", colour = "black", alpha = 0.3)
a + labs(x = "Phred score", y="Density") + theme_light()
dev.off()
png(filename="./t2d_100_lqual.png")
a + labs(x = "Phred score", y="Density") + theme_light() + xlim(0, 100)
dev.off()

print(summary(var_qual$qual))

# ##Variant mean depth###################################################
# var_depth <- read_delim("./t2d.ldepth.mean", delim = "\t",
#            col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
# png(filename="./t2d_total_ldepth_mean.png")
# a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "#a5a5a5", colour = "black", alpha = 0.3)
# a + theme_light()
# dev.off()
# png(filename="./t2d_100_ldepth_mean.png")
# a + theme_light() + xlim(0, 100)
# dev.off()

# print(summary(var_depth$mean_depth))

##Variant missingness###################################################
var_miss <- read_delim("./t2d.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
png(filename="./t2d_lmiss.png")
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "#a5a5a5", colour = "black", alpha = 0.3)
a + labs(x = "Missingness", y="Density") + theme_light()

dev.off()

print(summary(var_miss$fmiss))

##Minor allele frequency###################################################
var_freq <- read_delim("./t2d.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
# find minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
png(filename="./t2d_frq.png")
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "#a5a5a5", colour = "black", alpha = 0.3)
a + labs(x = "Allele Frequency", y="Density") + theme_light()

dev.off()

print(summary(var_freq$maf))


# ##Mean depth per individual###################################################
# ind_depth <- read_delim("./t2d.idepth", delim = "\t",
#                         col_names = c("ind", "nsites", "depth"), skip = 1)
# png(filename="./t2d_idepth.png")                        
# a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "#a5a5a5", colour = "black", alpha = 0.3)
# a + theme_light()
# dev.off()


# ##Proportion of missing data per individual####################################
# ind_miss  <- read_delim("./t2d.imiss", delim = "\t",
#                         col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
# png(filename="./t2d_imiss.png")  
# a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "#a5a5a5", colour = "black", alpha = 0.3)
# a + theme_light()
# dev.off()


##Heterozygosity and inbreeding coefficient per individual######################
ind_het <- read_delim("./t2d.het", delim = "\t",
           col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
png(filename="./t2d_het.png")
a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "#a5a5a5", colour = "black", alpha = 0.3)
a + labs(x = "Wright’s inbreeding coefficient (F)", y="Number of Samples") + theme_light()
dev.off()

##### HWE ##########################################################################


ind_het <- read_delim("./t2d.hwe", delim = "\t",
           col_names = c("chrNo","posNo", "obs", "exp", "pval", "hwe", "b", "c"), skip = 1)
png(filename="./t2d_hwe2.png")
a <- ggplot(ind_het, aes(hwe)) + geom_density(fill = "#a5a5a5", colour = "black", alpha = 0.3)
a + labs(x = "Hardy-Weinberg Equilibrium (HWE) p-value", y="Density") + theme_light() + xlim(0, 0.25)

dev.off()


library(qqman)
ind_hwe <- read_delim("./t2d.hwe", delim = "\t",
           col_names = c("chrNo","posNo", "obs", "exp", "pval", "hwe", "b", "c"), skip = 1,
           col_types = 
  		cols(chrNo = col_double(),
  			posNo = col_double(),
  			obs = col_character(),
  			exp = col_character(),
  			pval = col_double(),
  			hwe=col_double(),
  			b=col_double(),
  			c=col_double()))
	
	chr = sort(unique("ind_hwe$chrNo"))
	print(length(ind_hwe$hwe))
	
	ind_hwe$hwe[is.na(ind_hwe$hwe)] <- 0.0001
	print(sum(ind_hwe$hwe >0.05))
	print(sum(ind_hwe$hwe <=0.05))
	
# [1] 126760
# [1] 121577
# [1] 5183

	if (length(chr) == 24) { 
	labs = c(1:22, "X", "Y")
	} else if (length(chr) == 23 & tail(chr, n=1) == 24) {
	labs = c(1:22, "Y")
	} else if (length(chr) == 23 & tail(chr, n=1) == 23) {
	labs = c(1:22, "X")
	} else {
	labs = c(1:22)
	}

	ind_hwe2 <- na.omit(ind_hwe)


	print('>>> Saving figure...')
	path=sprintf("t2d_hwe.png")
	png(filename=path, width = 1000, height = 500, units = "px")
	manhattan(ind_hwe2, 
				chr = "chrNo", 
				bp = "posNo", 
				p = "hwe",
				chrlabs = labs,
				ylim=c(0,25), 
				suggestiveline = -log10(0.01), 
				genomewideline = FALSE)
	dev.off()
















