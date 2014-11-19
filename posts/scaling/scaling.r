library(phylosim)
library(rphastRegression)

tree <- read.tree('mcc.nh')
root.seq <- NucleotideSequence(length=10)
p<-GTR(rate.params=list("a"=1, "b"=1, "c"=3,"d"=2, "e"=1, "f"=1), base.freqs=c(1,1,1,1)/4)
attachProcess(root.seq,p)
sampleStates(root.seq)

sim <- PhyloSim()
sim$phylo <- tree
sim$rootSeq <- root.seq
Simulate(sim)
saveAlignment(sim,file="sim10.fas", skip.internal=TRUE)

nms <- tree$tip.label
pedvMSA <- read.msa("sim.fas")
treeChar <- write.tree(tree)

mod <- structure(list(alphabet = "ACGT", backgd = c(0.25, 0.25, 0.25, 
0.25), rate.matrix = structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0), .Dim = c(4L, 4L), .Dimnames = list(c("0", 
"1", "2", "3"), NULL)), subst.mod = "UNREST", likelihood = -109.517254528469, 
    eta.coefficients = c(0.2, 0, 0), design.matrix = structure(c(1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
    1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(12L, 
    3L)), npenalties = 100), .Names = c("alphabet", "backgd", 
"rate.matrix", "subst.mod", "likelihood", "eta.coefficients", 
"design.matrix", "npenalties"), class = "tm")


mod$tree <- treeChar

#mod$design.matrix <- cbind(mod$design.matrix, matrix(runif(n=12*1), nrow=12))
#mod$design.matrix[,-1]
#mod$design.matrix[,-1] <- mod$design.matrix[,-1]/.3892

mod2 <- phyloFit(pedvMSA, init.mod=mod, subst.mod='UNREST',
                 no.opt=c('backgd', 'branches'), ninf.sites=1, design.matrix=mod$design.matrix, log.file='out')

matplot(exp(mod2$penalties), mod2$beta[,-1], type='l')


mod3 <- phyloFit(pedvMSA, subst.mod='UNREST', tree=mod$tree,
                 no.opt=c('backgd', 'branches'), ninf.sites=1, design.matrix=mod$design.matrix)


#write.msa(pedvMSA, file='pedv-regions.fasta', format='FASTA')



