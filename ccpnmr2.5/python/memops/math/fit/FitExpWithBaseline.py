import math

from memops.math.fit.logLinearFit import logLinearFit
from memops.math.fit.FitNonlinear import FitNonlinear

class FitExpWithBaseline(FitNonlinear):

  name = 'Exponential with Baseline'

  paramNames = ['A', 'B', 'C']

  equation = 'y = A exp(-Bx) + C'

  description = '''This does a fit to an exponentially decaying
function with a baseline, %s.''' % equation

  def getInitParams(self, xs, ys, weights=None):

    result = logLinearFit(xs, ys, weights)
    (A, B) = result[0]

    return (A, B, 0.0)

  def getValue(self, x, params):

    (A, B, C) = params
    y = A * math.exp(-B*x) + C

    return y

  def getValueAndDerivative(self, x, params):

    (A, B, C) = params
    t = math.exp(-B*x)
    y = A * t + C
    dy_dparams = (t, -x*y, 1.0)

    return (y, dy_dparams)

