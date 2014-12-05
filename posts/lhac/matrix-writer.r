set.seed(23412)
M <- cbind(rep(1, 12),
           c(0,0,0,0,0,0,0,0,1,0,0,1),
           c(1,0,0,1,0,0,0,0,0,0,0,0))
M <- cbind(M, matrix(runif(12*(34 - 3), min=-.5, max=.5), nrow=12))
cat(t(M), file='designMat2')
