import math

from memops.math.fit.logLinearFit import logLinearFit
from memops.math.fit.FitNonlinear import FitNonlinear

class FitExp(FitNonlinear):

  name = 'Exponential'

  paramNames = ['A', 'B']

  equation = 'y = A exp(-Bx)'

  description = '''This does a fit to an exponentially decaying
function, %s.''' % equation 

  def getInitParams(self, xs, ys, weights=None):

    result = logLinearFit(xs, ys, weights)
    (A, B) = result[0]

    return (A, B)

  def getValue(self, x, params):

    (A, B) = params
    y = A * math.exp(-B*x)

    return y

  def getValueAndDerivative(self, x, params):

    A = params[0]
    B = params[1]
    t = math.exp(-B*x)
    y = A * t
    dy_dparams = (t, -x*y)

    return (y, dy_dparams)

