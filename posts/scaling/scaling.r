library(phylosim)
library(rphastRegression)


tree <- read.tree('mcc.nh')
root.seq <- NucleotideSequence(length=50)
p<-GTR(rate.params=list("a"=1, "b"=1, "c"=3,"d"=2, "e"=1, "f"=1), base.freqs=c(1,1,1,1)/4)
attachProcess(root.seq,p)
sampleStates(root.seq)

sim <- PhyloSim()
sim$phylo <- tree
sim$rootSeq <- root.seq
Simulate(sim)
saveAlignment(sim,file="sim.fas", skip.internal=TRUE)

nms <- tree$tip.label
pedvMSA <- read.msa("sim.fas")
treeChar <- write.tree(tree)

mod <- structure(list(alphabet = "ACGT", backgd = c(0.25, 0.25, 0.25,
0.25), rate.matrix = structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0), .Dim = c(4L, 4L), .Dimnames = list(c("0", "1", "2", "3"),
NULL)), subst.mod = "UNREST", likelihood = -109.517254528469,
eta.coefficients = c(0.2, 0, 0)), .Names = c("alphabet", "backgd",
"rate.matrix", "subst.mod", "likelihood", "eta.coefficients"), class =
"tm")
mod$tree <- treeChar

mod2 <- phyloFit(pedvMSA, init.mod=mod, subst.mod='UNREST',
                 no.opt=c('backgd', 'branches'), ninf.sites=1)

#write.msa(pedvMSA, file='pedv-regions.fasta', format='FASTA')



