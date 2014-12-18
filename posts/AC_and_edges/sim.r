library(Rcpp)

cppFunction('
  NumericMatrix sisPairCt(double lambda, double E, double m0, double n0, int maxSteps, int k){
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

k=4
E=2000
N <- E/k
lambda = .3
muEq = 1/((lambda*k)^2 + lambda^2*k - lambda)
nuEq = (1 + (1 - 2*lambda*k) * muEq)
NsiEq = E * (1 - muEq - nuEq) / 2
NssEq = E * muEq
NsEq = (NsiEq + NssEq)/k
NiEq = N - NsEq

ts <- sisPairCt(lambda=lambda, E=E, m0=E*muEq, n0=E*nuEq, maxSteps=100000, k=k);
colnames(ts) <- c('time', 'Nss', 'Nii', 'Nsi', 'Ns')
ts <- cbind(ts, Ni=N - ts[, 'Ns'])
yfun <- stepfun(ts[-1, 'time'], ts[, 'Ni'])
domain <- range(knots(yfun))

getAvg <- function(f, lower, upper, ...){
    dx <- 1/(upper - lower)
    tmpf <- function(x) f(x)*dx
    integrate(tmpf, lower=lower, upper=upper, ...)$value
}

mu <- getAvg(yfun, lower=domain[1], upper=domain[2], subdivisions=10000)

devfun <- function(x) yfun(x) - mu

getAutoCov <- function(f, lag, lower, upper, ...){
    dx <- 1/(upper - lower)
    tmpf <- function(x) f(x + lag)*f(x) *dx
    integrate(tmpf, lower=lower, upper=upper, ...)$value
}

aCov <- getAutoCov(devfun, 1, domain[1], domain[2], subdivisions=100000)
sigSq <- getAutoCov(devfun, 0, domain[1], domain[2], subdivisions=100000)

autoCor <- aCov/sigSq
