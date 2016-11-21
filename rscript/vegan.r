library("phyloseq")
library("vegan")
library("vegan3d")
library("ade4")
library("plyr")
setwd("~/Dev/mucobiome/rscript")

# Get biom file 
biomFile = import_biom("xp2.biom", parseFunction=parse_taxonomy_greengenes, parallel = T)

# Divide by 2 
biom2 = transform_sample_counts(biomFile, function(OTU) OTU/2 )

#Agglomeration 
biom3 = tax_glom(biom2, "Species")

# # normalize 
normBiom <- transform_sample_counts( biom3, function(x) x/sum(x) )

#normBiom = transform_sample_counts(biom3, function(x) {x / sum(x)})
# 
# # filter 
# filterBiom = filter_taxa(normBiom, function(x) mean(x) > 1e-5, TRUE)
# 
# 

#plot_bar(normBiom, fill = 'Kingdom') 

plot_bar(normBiom, fill = "Species")
