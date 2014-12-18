library(Rcpp)

cppFunction('
  NumericMatrix sisPairCt(double lambda, double E, double m0, double n0, int maxSteps, int k){
  // taken in part from http://phylodynamics.blogspot.com/2013/06/continuous-time-markov-chain-ctmc.html
    NumericMatrix res(maxSteps, 5L);
    double t = 0, r1, r2, r3, r4, tot, r, dt;
    int Nss = m0;
    int Nii = n0;
    int Nsi, Nisi, Nssi, Ns, n=0;
    do{
      Nsi = (E - Nss - Nii) / 2;
      Ns = (Nsi + Nss) / k;
      Nisi = Nsi * Nsi / Ns;
      Nssi = Nss * Nsi / Ns;
      r1 = Nsi;
      r2 = r1 + lambda * (Nsi + Nisi);
      r3 = r2 + lambda * Nssi;
      r4 = r3 + Nii;
      tot = r1 + r2 + r3 + r4;
      dt = rexp(1,tot)[0];
      t += dt;
      r = r4 * runif(1)[0];
      if (r < r1) {
        Nss += 2;
      } else if (r < r2){
        Nii += 2;
      } else if (r < r3){
        Nss -= 2;
      } else {
        Nii -= 2;
      }
      res(n,0) = t;
      res(n,1) = Nss;
      res(n,2) = Nii;
      res(n,3) = Nsi;
      res(n,4) = Ns;
      n++;
    } while (n < maxSteps && Nii > 0);
    return res;
  }')

set.seed(1014)

ts <- sisPairCt(lambda=.3, E=2000, m0=0, n0=500, 100000, 4);


N <- 1000
tmpf <- function(x) sisd(R0=x, g=.1, h=0.001, N=N, I0=N/2, maxsteps=100000)

R0 <- seq(from=.2, to=2, len=100)
system.time(ts <- lapply(R0, tmpf))

tmpf <- function(x, burnin=1000) {
    del <- seq_len(1000)
    y <- x[-del, 2]
    mu <- mean(y)
    n <- length(y)
    ss <- var(y)
    c(ac=mean((y[-length(y)] - mu)* (y[-1] - mu) )/ ss,
      ystar=mu/N)
}
sim <- sapply(ts, tmpf)

tmpf <- function(x) infeq(R0=x, recoveryRate=.1, nImports=1, nNodes=N)[c('ac', 'ystar')]
theo <- sapply(R0, tmpf)

sim <- data.frame(t(sim), R0=R0, type='simulation')
theo <- data.frame(t(theo), R0=R0, type='theory')

cmp <- rbind(sim, theo)

plot(ac~R0, data=cmp[cmp$type=='theory', ], type='l',
     ylab='lag-1 autocorr.', xlab='R0 approx.')
points(ac~R0, data=cmp[cmp$type=='simulation',], col=2)

plot(ystar~R0, data=cmp[cmp$type=='theory', ], type='l',
     ylab='equilibrium prevalence', xlab='R0 approx.')
points(ystar~R0, data=cmp[cmp$type=='simulation',], col=2)

#' The theory (lines) and simulation (red circles) agree fairly well
#' even though the theory is for continuous time and simulations were
#' in discrete time. Also, theoretical equilibrium prevalence is the
#' limit as the population size approaches infinity.
