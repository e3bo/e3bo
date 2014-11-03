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

# one step model

p1 = lambda1/(lambda1 + lambda2)
p2 = lambda2/(lambda1 + lambda2)


xg, xs, ps, pg, lam_s, lam_g = var('xg, xs, ps, pg, lam_s, lam_g')

l(lambda1, lambda2) = xg * log(p1) + xs * log(p2)
Hsingle = l.diff(2)

S = lam_s + lam_g
ps = lam_s^2/S^2 + lam_g^2/S^2
pg = 2*lam_s*lam_g/S^2

l(lam_g, lam_s) = xg * log(pg) + xs * log(ps)

Hdouble = l.diff(2)





