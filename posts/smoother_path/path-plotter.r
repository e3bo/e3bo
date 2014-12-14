lmd <- scan('/tmp/lmd')
m <- read.csv('/tmp/dat2', header=F)
xmax = lmd[(which(rowSums(m[,-1]) > 0) - 1)[1]] * 1.1
png('path.png', width=800, height=800, pointsize=18)
matplot(lmd, m[seq_along(lmd),-1], type='l', xlab=expression(l[1] ~ "penalty"),
        ylab="Regression coefficients", xlim=c(0,xmax))
dev.off()
