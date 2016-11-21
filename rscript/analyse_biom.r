library("phyloseq")
library("vegan")
library("vegan3d")
library("ade4")
setwd("~/Dev/mucobiome/rscript")

# Get biom file 
biomFile = import_biom("xp2.biom", parseFunction=parse_taxonomy_greengenes, parallel = T)


# # Select how to filter otu
taxatype = "Species"
#
#Extract OTU . divide by 2.. This was an error in my pipeline
otu = as.data.frame(otu_table(biomFile)/2)
taxa = as.vector(tax_table(biomFile)[,c(taxatype)])
otu$taxa = taxa

# Aggregate according taxaType
aggdata <-aggregate(otu[,1:ncol(otu)-1], by=list(taxa), FUN=sum, na.rm=TRUE)
rownames(aggdata) = aggdata[,1]
aggdata[,1] = NULL
final = as.data.frame(t(aggdata))
nt#
#
# #Normalize
#
finalN <- scale(final, center=F, scale=colSums(final))
finalN = floor(finalN * 1000000)

#table.value(finalN, csize=0.4)
#
#
acpi = dudi.pca(finalN, scale=F, scannf = F, nf = 2)
scatter(acpi)


# Compute statistics
#
# top10 = sort(colSums(finalN), decreasing = T)[1:10]
# names(top10)
#
# shanon = diversity(finalN,index="shannon")
# alpha = fisher.alpha(finalN)
#
# spa = specaccum(finalN)
# plot(spa)
#
mad = decorana(finalN)
# 
plot(mad)
# points(mad, display = "spec", cex = 0.8, pch = 21, col = "red",bg = "yellow")
# text(mad, display = "sites", cex = 0.7, col = "blue")
# 

