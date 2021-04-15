

#
# File: qualityControl.R
# Created Date: Monday, September 28th 2020
# Author: Debora Antunes
# -----
# Last Modified: Wednesday, September 30th 2020, 10:52:14 am
# -----
#

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("GWASTools")
# BiocManager::install("SNPRelate")
library(GWASTools)
library(SNPRelate)
library(plyr)

####### CONVERTING VCF TO GDS ###########
vcf.fn <- "files/data/vcf/cases/All_PT.vcf.gz"
gdsfile <- "files/data/variants/snps.gds"
snpgdsVCF2GDS(vcf.fn, gdsfile, verbose=TRUE)


##### CONVERT CHR VALUES ######
genofile = snpgdsOpen(gdsfile, readonly=FALSE)
chr = read.gdsn(index.gdsn(genofile, "snp.chromosome"))
unique(chr)
chr = revalue(chr, c('X'='23', 'Y'='25'))
unique(chr)
add.gdsn(genofile,'snp.chromosome',chr,replace=TRUE)

# open the GDS file and create a GenotypeData object
gdsfile <- "files/data/variants/snps.gds"
gds <- GdsGenotypeReader(gdsfile)
nscan(gds)
nsnp(gds)
head(getScanID(gds))
head(getSnpID(gds))
tail(getChromosome(gds))
head(getPosition(gds))

########## Missingness and heterozygosity within samples ################
## Missingness
genoData <- GenotypeData(gds)
miss <- missingGenotypeByScanChrom(genoData)
# Examine the results
miss.rate <- t(apply(miss$missing.counts, 1, function(x) {x / miss$snps.per.chr}))
miss.rate <- as.data.frame(miss.rate)
cols <- names(miss.rate) %in% c(1:22, "X", "Y")
boxplot(miss.rate[,cols], 
        main="Missingness by Chromosome", 
        ylab="Proportion Missing", 
        xlab="Chromosome")

## Heterozygosity
het.results <- hetByScanChrom(genoData)
cols <- colnames(het.results) %in% c(1:22, "X", "Y")
boxplot(het.results[,cols], 
        main="Heterozygoty by Chromosome", 
        ylab="Heterozygoty rate", 
        xlab="Chromosome")
close(genoData)


############# Relatedness and IBD Estimation ############
gdsobj = snpgdsOpen(gdsfile, readonly=FALSE)
ibdobj <- snpgdsIBDKING(gdsobj)
length(ibdobj$kinship[ibdobj$kinship < 0.0625])#fourth degree
length(ibdobj$snp.id)