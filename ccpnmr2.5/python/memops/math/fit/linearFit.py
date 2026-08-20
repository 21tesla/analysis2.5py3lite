import math

from memops.math.MathException import MathException

STT_EPS = 1.0e-10

# fit y = ax+b
# result = [(a,b), chisq, ysFit, optionally (aDev, bDev)]
# where chisq is some of error squared divided by (numPoints-2)
def linearFit(xs, ys, weights=None, findDev=False):

  n = len(xs)
  if n < 2:
    raise MathException('n = %d, need at least 2 points' % n)

  assert n == len(ys), 'len(xs) = %d, len(ys) = %d' % (n, len(ys))
  if weights:
    assert n == len(weights), 'len(xs) = %d, len(weights) = %d' % (n, len(weights))

  if weights:
    s = sx = sy = 0
    for i in range(n):
      w2 = weights[i] * weights[i]
      s += w2
      sx += w2 * xs[i]
      sy += w2 * ys[i]
  else:
    s = float(n)
    sx = sy = 0
    for i in range(n):
      sx += xs[i]
      sy += ys[i]

  a = 0.0
  stt = 0.0
  t = n * [0]
  for i in range(n):
    t[i] = xs[i] - sx/s
    if weights:
      t[i] *= weights[i]
    d = t[i] * ys[i]
    if weights:
      d *= weights[i]
    a += d
    stt += t[i] * t[i]

  if stt < STT_EPS:
    raise MathException('x values all the same (it seems)')

  a /= stt
  b = (sy - a*sx) / s;

  ysFit = n * [0]
  chisq = 0
  for i in range(n):
    ysFit[i] = a*xs[i] + b
    dy = ys[i] - ysFit[i]

    if weights:
      dy *= weights[i] 

    chisq += dy * dy

  if n > 2:
    chisq /= n-2
  else:
    chisq = 0

  result = [(a, b), chisq, ysFit]

  if findDev:
    std_a = 1.0 / stt
    std_b = (1.0 + sx*sx/(s*stt)) / s
    if not weights:
      std_a *= chisq
      std_b *= chisq
    std_a = math.sqrt(std_a)
    std_b = math.sqrt(std_b)
    result.append((std_a, std_b))

  return result

if __name__ == '__main__':

  import random

  def printVec(msg, v):
    s = ', '.join(['%.3f' % w for w in v])
    print '%s: %s' % (msg, s)

  a = 1.4
  b = 2.3
  n = 7
  noise = 0.3
  xs = [0.5*i for i in range(n)]
  ys = [a*x+b+noise*(2*random.random()-1) for x in xs]

  printVec('xs', xs)
  printVec('ys', ys)
  (params, chisq, ysFit) = linearFit(xs, ys)
  printVec('params', params)
  printVec('ysFit', ysFit)
  print 'chisq = %.3f' % chisq

