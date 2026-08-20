import math

from memops.math.gaussJordan import gaussJordan
from memops.math.MathException import MathException

# x[n], y[n] are the data points to be fitted
# w[n] is the weight for the contribution to chisq of the particular point
# a[m] are the coefficients of the nonlinear model
# func is user-provided and returns value and derivative of fitting function at a

# private function
def _findLinearised(x, y, a, func, w=None):

  m = len(a)
  n = len(x)

  alpha = m * [0]
  for j in range(m):
    alpha[j] = m * [0.0]

  beta = m * [0.0]

  yFit = n * [0]
  chisq = 0.0
  for i in range(n):
    (yFit[i], dy_da) = func(x[i], a)
    dy = y[i] - yFit[i]

    for j in range(m):
      wgt = dy_da[j]
      if w:
        wgt *= w[i]
      for k in range(m):
        alpha[j][k] += wgt * dy_da[k]
      beta[j] += wgt * dy

    csq = dy * dy
    if w:
      csq *= w[i]
    chisq += csq

  for j in range(m-1):
    for k in range(j+1, m):
      alpha[j][k] = alpha[k][j]

  return (chisq, alpha, beta, yFit)

# private function
def _nonlinearCovariance(x, y, a, func, lmbda, alpha, beta):

  m = len(a)
  covar = m * [0]
  for j in range(m):
    covar[j] = m * [0]
    for k in range(m):
      covar[j][k] = alpha[j][k]
    covar[j][j] *= 1 + lmbda

  x = beta[:]
  gaussJordan(covar, x)

  return (covar, x)

# private function
def _nonlinearModel(x, y, a, func, w=None, lmbda=None, chisq=None, alpha=None, beta=None, yFit=None):

  m = len(a)

  if lmbda is None:
    lmbda = 0.001
    (chisq, alpha, beta, yFit) = _findLinearised(x, y, a, func, w=w)

  (covar, da) = _nonlinearCovariance(x, y, a, func, lmbda, alpha, beta)

  ap = m * [0]
  for j in range(m):
    ap[j] = a[j] + da[j]

  (newChisq, covar, da, yNewFit) = _findLinearised(x, y, ap, func, w=w)

  if newChisq < chisq:
    a = ap
    lmbda *= 0.1
    chisq = newChisq
    for j in range(m):
      for k in range(m):
        alpha[j][k] = covar[j][k]
    beta = da
    yFit = yNewFit
  else:
    lmbda *= 10.0

  return (a, lmbda, chisq, alpha, beta, yFit)
 
CHISQ_STOP_CRITERION = 0.1
MAX_ITERATION = 20
MAX_CONDITION = 4

def nonlinearFit(x, y, a, func, w=None, noise=None, findDev=False,
                 maxIteration=MAX_ITERATION, maxCondition=MAX_CONDITION):

  n = len(x)
  assert n == len(y), 'len(x) = %d, len(y) = %d' % (n, len(y))

  if n < 1:
    raise MathException('len(x) must be > 0')

  m = len(a)
  if m < 1:
    raise MathException('len(a) must be > 0')

  if not noise:
    noise = 0.1 * max([abs(t) for t in y])  # arbitrary

  chisq_stop_criterion = CHISQ_STOP_CRITERION * noise * noise

  iter = cond = 0
  (a, lmbda, chisq, alpha, beta, yFit) = _nonlinearModel(x, y, a, func, w=w)

  while iter < maxIteration and cond < maxCondition:
    old_chisq = chisq
    (a, lmbda, chisq, alpha, beta, yFit) = _nonlinearModel(x, y, a, func, w=w,
                     lmbda=lmbda, chisq=chisq, alpha=alpha, beta=beta, yFit=yFit)
    if chisq > old_chisq:
      cond = 0
    elif (old_chisq - chisq) < chisq_stop_criterion:
      cond += 1

    iter += 1

  if n > m:
    chisq /= (n-m)
  else:
    chisq = 0

  result = [a, chisq, yFit]

  if findDev:
    (covar, xx) = _nonlinearCovariance(x, y, a, func, alpha, beta)
    aDev = [math.sqrt(chisq*max(covar[i][i], 0)) for i in range(m)]
    result.append(aDev)

  return result

if __name__ == '__main__':

  def func(x, a):

    y = a[0]*x*x + a[1]*x + a[2]
    dy_da = [x*x, x, 1]
    result = (y, dy_da)

    return result

  def printVec(msg, v):

    s = ', '.join(['%.3f' % w for w in v])
    print '%s: %s' % (msg, s)

  import random

  noise = 0.1
  a = [1.0, 1.5, 2.0]
  x = [0.5, 1.0, 1.5, 2.0, 2.5]
  y = [a[0]*w*w+a[1]*w+a[2]+noise*(2*random.random()-1) for w in x]

  printVec('x', x)
  printVec('y', y)
  (a, chisq, yFit) = nonlinearFit(x, y, a, func, noise=noise)
  printVec('a', a)
  printVec('yFit', yFit)
  e2 = sum([(y[i]-yFit[i])**2 for i in range(len(y))])
  print 'chisq = %.3f' % chisq
  print 'e2 = %.3f' % e2

