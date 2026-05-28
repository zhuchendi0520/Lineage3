library(devtools)
devtools::install_github("gtonkinhill/fastbaps")
library(fastbaps)
library(ape)
library(ggtree)
library(phytools)
library(ggplot2)
setwd("/Users/zhuchendi/Desktop/")
sparse.data <- import_fasta_sparse_nt("all_strains.fadel-InvMisF5.bak.fa")
sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric")
baps.hc <- fast_baps(sparse.data)
best.partition <- best_baps_partition(sparse.data, baps.hc)
iqtree <- phytools::read.newick("L3.tre")
plot.df <- data.frame(id = colnames(sparse.data$snp.matrix), fastbaps = best.partition, 
                      stringsAsFactors = FALSE)
write.csv(plot.df,"fastbaps.csv")

gg <- ggtree(iqtree)

f2 <- facet_plot(gg, panel = "fastbaps", data = plot.df, geom = geom_tile, aes(x = fastbaps), 
                 color = "blue")
f2
