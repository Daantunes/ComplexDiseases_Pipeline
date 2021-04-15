
#
# File: network.R
# Created Date: Tuesday August 18th 2020
# Author: Debora Antunes
# -----
# Last Modified: Thursday, April 15th 2021, 12:08:11 pm
# -----
#

fun = function(interations, rData){
	#################################################################
	### FUNCTION THAT CREATES THE NETWORK AND EXTRACTS THE CENTRALITY METRICS
	### var: a .csv file with the interactions
	### rData: a previously created network
	##################################################################

    
    ### Necessary packages  ###
    # install.packages('igraph')
    # install.packages('files/src/R/dmGWAS_3.0.tar.gz', repos=NULL, type="source")
    library('dmGWAS')

    if(rData == FALSE){
        print('>>> Creating network...')
        geneweight <- read.delim('../../data/proteins/genes.csv',sep=',')
        print(head(geneweight,4))
        network <- read.delim(interations, sep=',')
        print(head(network, 4))
        
        print('>>> Running dmGWAS; It can take days, please be pacient')
        res.list <- dms(network, geneweight, expr1=NULL, expr2=NULL, d=2, r=0.2)

        ############### SAVE INTERMEDIATE FILES
        print('>>> Saving network')  
        save.image("../R/all_data.RData")

    }
    else{
        print('>>> Loading network...')
        load(rData)
    }

    ################ SIZE OF EACH MODULE
    size = c()
    for (i in 1:length(res.list$genesets.clear)){
        a=lengths(res.list$genesets.clear[i], use.names = FALSE)
        size = c(size, a)
    }
    boxplot(size)
    print('>>> Boxplot of sizes of each module; Write \'1\' to save the figure; Write \'0\' to continue')
    x <- scan(nmax=1)
    if(x == 1){
        png(filename="../../data/figures/r_box.png")
        boxplot(size)
        dev.off()
    }

    ############### STATISTICS FOR RES.LIST
    statResults(res.list, top=100)
    print('>>> Statistics for the top 100 modules; Write \'1\' to save the figure; Write \'0\' to continue')
    x <- scan(nmax=1)
    if(x == 1){
        png(filename="../../data/figures/r_stats.png")
        statResults(res.list, top=100)
        dev.off()
    }

    ############## SELECT TOP 
    simpleChoose(res.list, top=1, plot=T)
    print('>>> Top 1 module representation; Press enter to continue')
    x <- scan(nmax=1)
    selected = simpleChoose(res.list, top=50, plot=F)
    g_prot = selected$subnetwork
    name = V(g_prot)$name
    weight = V(g_prot)$weight
    degree_distribuition = c()
    betweenness = c()
    closeness = c()
    print('Number of genes:')
    print(length(name))
    for (i in 1:length(name)){
        prot = name[i]
        d = igraph::degree(g_prot, prot, loops=TRUE, normalized=TRUE)
        b = betweenness(g_prot, v = prot, normalized=TRUE)
        c = closeness(g_prot, v = prot, normalized=TRUE)
        degree_distribuition = c(degree_distribuition, d)
        betweenness = c(betweenness, b)
        closeness = c(closeness, c)
    }
    df = data.frame(name, weight, degree_distribuition, betweenness, closeness)
    print('>>> Saving file of metrics in : data/proteins/r_genes.csv')
    write.csv(df[order(df$w),],"../../data/proteins/r_genes.csv", row.names = FALSE)
    print(df)

}
