def calcFit(clazz, xs, ys, weights=None, noise=None, params=None, devMethod=None):

  result = clazz().calcFit(xs, ys, weights, noise, params, devMethod)

  return result

if __name__ == '__main__':

  import math
  import random
  from memops.math.fit.FitExp import FitExp

  def printVec(msg, v):

    s = ', '.join(['%.3f' % w for w in v])
    print '%s: %s' % (msg, s)

  noise = 0.1
  params = [1.8, 0.5]
  xs = [0.5, 1.0, 1.5, 2.0, 2.5]
  ys = [params[0]*math.exp(-params[1]*x)+noise*(2*random.random()-1) for x in xs]

  printVec('actual_params', params)
  printVec('xs', xs)
  printVec('ys', ys)

  (params, chisq, ysFit, paramsDev) = calcFit(FitExp, xs, ys, noise=noise, devMethod='bootstrap')
  printVec('final_params', params)
  printVec('ysFit', ysFit)
  print 'chisq: %.3f' % chisq
  printVec('paramsDev', paramsDev)

