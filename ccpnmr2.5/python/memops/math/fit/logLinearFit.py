import math

from memops.math.fit.linearFit import linearFit

def logLinearFit(xs, ys, weights=None, findDev=False):

  n = len(xs)
  assert n == len(ys), 'len(xs) = %d, len(ys) = %d' % (n, len(ys))
  if weights:
    assert n == len(weights), 'len(xs) = %d, len(weights) = %d' % (n, len(weights))

  # only makes sense to fit y of same sign
  # look for positive and negative and fit one with most terms

  cnt1 = cnt2 = 0
  for i, y in enumerate(ys):
    if weights:
      y *= weights[i]
    if y > 0:
      cnt1 += y
    else:
      cnt2 += y

  if cnt1 >= cnt2:
    s = 1
  else:
    s = -1

  xsToFit = []
  ysToFit = []
  if weights:
    weightsToFit = []
  else:
    weightsToFit = None
  for i, y in enumerate(ys):
    t = s * y
    if t > 0:
      xsToFit.append(xs[i])
      ysToFit.append(math.log(t))
      if weights:
        weightsToFit.append(weights[i])

  n = len(xsToFit)
  if n < 2:
    raise MathException('n = %d, need at least 2 points of same sign' % n)

  result = linearFit(xsToFit, ysToFit, weightsToFit, findDev)
  ((a, b), chisq, ysFit) = result[:4]

  # above fit is ax+b, we need A exp(-Bx)
  A = s * math.exp(b)
  B = -a
  ysFit = [A*math.exp(-B*x) for x in xs]

  chisq = 0
  for i, y in enumerate(ys):
    dy = ysFit[i] - y
    chisq += dy*dy
  n = len(xs)
  if n > 2:
    chisq /= n-2
  else:
    chisq = 0

  result = [(A,B), chisq, ysFit]

  if findDev:
    (std_a, std_b) = result[4]
    # below is approximate
    std_A = A * std_b  # factor of A from nonlinear transformation Jacobian
    std_B = std_a
    result.append((std_A, std_B))

  return result

