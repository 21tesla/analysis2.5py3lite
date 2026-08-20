import math
import random

from memops.math.MathException import MathException

devMethods = ['normal', 'bootstrap']

# used to model arbitrary functions y = F(x, params)
# and to fit given (x, y) to find best param
class FitClass:

  # could be overridden in subclass
  floatFormat = '%.4f'

  # should be overridden in subclass
  # a name for the function to use in menus
  name = None

  # should be overridden in subclass
  # name of parameters
  paramNames = None

  # should be overridden in subclass
  # textual description of equation
  equation = None

  # should be overridden in subclass
  # a human-readable description of the function
  description = None

  # returns y at (x, params)
  # input:
  #   x float
  #   params tuple/list[float] of size numParams
  # output:
  #   y float
  def getValue(self, x, params):

    raise Exception('must override in subclass')

  # returns ys at (xs, params)
  # input:
  #   xs tuple/list[float]
  #   params tuple/list[float] of size numParams
  # output:
  #   ys list[float]
  def getValues(self, xs, params):

    return [getValue(x, params) for x in xs]

  # returns [params, chisq, ysFit, optionally paramsDev]
  # input:
  #   xs tuple/list[float] of size N (some N)
  #   ys tuple/list[float] of size N
  #   weights tuple/list[float] of size N or None
  #   noise float, which gives an estimate of ys noise, or None
  #   params tuple/list[float] of size numParams, giving initial estimate of params, or None
  #   devMethod string, if want estimate of standard deviation of parameters
  #     if set, must be one of ('normal', 'bootstrap')
  # output:
  #   params list[float] of size numParams
  #   chisq float
  #   ysFit list[float] of size numParams
  #   paramsDev list[float] of size numParams, if calculated
  def calcFit(self, xs, ys, weights=None, noise=None, params=None, devMethod=None):

    raise Exception('must override in subclass')

  # returns paramsDev
  # input:
  #   xs tuple/list[float] of size N (some N)
  #   ys tuple/list[float] of size N
  #   weights tuple/list[float] of size N or None
  #   noise float, which gives an estimate of ys noise, or None
  #   params tuple/list[float] of size numParams, giving initial estimate of params, or None
  #   niters int, number of iterations to run for method
  # output:
  #   paramsDev list[float] of size numParams, if calculated
  # Reference:
  # Bootstrap Methods for Standard Errors, Confidence Intervals and Other
  # Measures of Statistical Accuracy
  # B. Efron and R. Tibshirani
  # Statistical Science, 1986, Vol. 1, No. 1, 54-77
  def bootstrapDev(self, xs, ys, weights=None, noise=None, params=None, niter=1000):

    def samplePoints():
      xxs = []
      yys = []
      n = len(xs)
      for k in range(n):
        m = int(n*random.random())
        m = min(m, n-1)
        xxs.append(xs[m])
        yys.append(ys[m])
      return (xxs, yys)

    if not params:
      params = self.calcFit(xs, ys, weights, noise)[0]

    nparams = len(params)
    paramsAvg = nparams * [0]
    paramsDev = nparams * [0]

    ngood = 0
    for i in range(niter):
      xs2, ys2 = samplePoints()
      try:
        params2 = self.calcFit(xs2, ys2, weights, noise, params)[0]
      except: # fit did not converge
        continue
      ngood += 1
      for j in range(nparams):
        paramsAvg[j] += params2[j]
        paramsDev[j] += params2[j] * params2[j]
      
    if ngood < 2:
      raise MathException('not enough good fits when sampling')

    for j in range(nparams):
      paramsAvg[j] /= ngood
      r = paramsDev[j] - ngood*paramsAvg[j]*paramsAvg[j]
      r = max(r, 0)
      paramsDev[j] = math.sqrt(r/(ngood-1))

    return paramsDev

