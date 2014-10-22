library(rphast)
library(ape)

tree <- read.tree('mcc.nh')

nms <- tree$tip.label
abs <- sapply(strsplit(nms, '_'), '[[', 2)
ind <- match(abs, state.abb)
reg <- state.region[ind]
regDNA <- factor(reg,
                 levels=c("Northeast", "South", "North Central", "West"),
                 labels=c('A', 'T', 'C', 'G'))
regDNA <- as.character(regDNA)
pedvMSA <- msa(regDNA, nms)

treeChar <- write.tree(tree)
mod <- phyloFit(pedvMSA, tree=treeChar, no.opt='branches', ninf.sites=1)

# TODO
# 1. stop fitting branch lengths
# 2. add inner loop to do regression



