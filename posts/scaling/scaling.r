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
#mod <- phyloFit(pedvMSA, tree=treeChar, subst.mod='UNREST', no.opt=c('backgd', 'branches'), ninf.sites=1)

mod <- structure(list(alphabet = "ACGT", backgd = c(0, 0.822784810126582, 
0.0506329113924051, 0.126582278481013), rate.matrix = structure(c(-881.283989747898, 
832.945800602126, 24.1690945728857, 24.1690945728857, 832.945800602126, 
-881.283989747898, 24.1690945728857, 24.1690945728857, 24.1690945728857, 
24.1690945728857, -181.967932737163, 133.629743591391, 24.1690945728857, 
24.1690945728857, 133.629743591391, -181.967932737163), .Dim = c(4L, 
4L), .Dimnames = list(c("0", "1", "2", "3"), NULL)), subst.mod = "UNREST", 
    likelihood = -109.517254528469, eta.coefficients = c(0.2, 0, 0)), .Names = c("alphabet", "backgd", "rate.matrix", 
"subst.mod", "likelihood", "eta.coefficients"), class = "tm")
mod$tree <- treeChar

mod2 <- phyloFit(pedvMSA, init.mod=mod, subst.mod='UNREST', no.opt=c('backgd', 'branches'), ninf.sites=1)

#write.msa(pedvMSA, file='pedv-regions.fasta', format='FASTA')



