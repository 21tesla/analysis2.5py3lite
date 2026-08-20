from memops.math.fit.FitClass import FitClass
from memops.math.fit.nonlinearFit import nonlinearFit

class FitNonlinear(FitClass):

  # default implementation which could be overridden
  # since not very efficient (although might not matter)
  # returns y at (x, params)
  # input:
  #   x float
  #   params tuple/list[float] of size numParams
  # output:
  #   y float
  def getValue(self, x, params):

    result = self.getValueAndDerivative(x, params)[0]

    return result

  # returns (y, dy_dparams) at (x, params)
  # input:
  #   x float
  #   params tuple/list[float] of size numParams
  # output:
  #   y float
  #   dy_dparams list[float] of size numParams
  def getValueAndDerivative(self, x, params):

    raise Exception('must override in subclass')

  # returns initial estimate for params
  # should be overridden in subclass
  # input:
  #   xs tuple/list of size N (some N)
  #   ys tuple/list of size N
  #   weights tuple/list[float] of size N or None
  # output:
  #   params list[float] of size numParams
  def getInitParams(self, xs, ys, weights=None):

    raise Exception('must override in subclass')

  # returns [params, chisq, ysFit, optionally paramsDev]
  def calcFit(self, xs, ys, weights=None, noise=None, params=None, devMethod=None):

    if not params:
      params = self.getInitParams(xs, ys, weights)

    findDev = (devMethod == 'normal')
    result = nonlinearFit(xs, ys, params, self.getValueAndDerivative, weights, noise, findDev)

    if devMethod == 'bootstrap':
      params = result[0]  # start with given answer
      paramsDev = self.bootstrapDev(xs, ys, weights, noise, params)
      result.append(paramsDev)

    return result

