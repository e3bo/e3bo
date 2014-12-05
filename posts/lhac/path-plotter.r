lmd <- scan('/tmp/lmd')
m <- read.csv('/tmp/dat2', header=F)
png('path.png')
matplot(lmd, m[-100,-1], type='l', xlab=expression(l[1] ~ "penalty"),
        ylab="Regression coefficients", xlim=c(5,150))
dev.off()
