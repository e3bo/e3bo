library(rphastRegression)
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
mod <- phyloFit(pedvMSA, tree=treeChar, subst.mod='UNREST', no.opt='branches', ninf.sites=1)

#mod$tree <- treeChar
#mod2 <- phyloFit(pedvMSA, init.mod=mod, subst.mod='UNREST', no.opt='branches', ninf.sites=1)

#write.msa(pedvMSA, file='pedv-regions.fasta', format='FASTA')



