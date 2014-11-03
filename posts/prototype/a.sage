p1, p2, p3, x1, x2, x3, lambda1, lambda2, lambda3 = var('p1, p2, p3, x1, x2, x3, lambda1, lambda2, lambda3')

p1 = lambda1/(lambda1 + lambda2 + lambda3)
p2 = lambda2/(lambda1 + lambda2 + lambda3)
p3 = lambda3/(lambda1 + lambda2 + lambda3)

l(lambda1, lambda2, lambda3) = x1 * log(p1) + x2 * log(p2) + x3 * log(p3)
H = l.diff(2)

def HDiag(lam, x, i, j):
  sumk = sum(lam)
  c1 = 1/sumk^2
  c2 = -1/lam[i]^2
  ret = x1*0
  for k in range(len(lam)):
     ret += x[k] * c1
  if i == j:
     ret += x[i] * c2
  return ret










