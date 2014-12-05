library(phylosim)

set.seed(12342)

tree <- read.tree('mcc.nh')
root.seq <- NucleotideSequence(length=50)
p<-GTR(rate.params=list("a"=1, "b"=1, "c"=3,"d"=2, "e"=1, "f"=1), base.freqs=c(1,1,1,1)/4)
attachProcess(root.seq,p)
sampleStates(root.seq)

sim <- PhyloSim()
sim$phylo <- tree
sim$rootSeq <- root.seq
Simulate(sim)
saveAlignment(sim,file="sim.fasta", skip.internal=TRUE)
