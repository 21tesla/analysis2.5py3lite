import math

from memops.math.fit.FitNonlinear import FitNonlinear

class FitInversionRecovery(FitNonlinear):

  name = 'Inversion Recovery'

  paramNames = ['equilibriumIntensity', 'R1']

  equation = 'y = A (1 - 2 exp(-Bx))'

  description = '''In an inversion recovery experiment the initial
magnetisation is inverted and so the intensity starts out at the
negative of the equilibrium magnetisation and over time reverts
to the equilibrium.  The eqution used is %s.
Here x is the time, y is the intensity, A is the equilibrium
intensity, B is R1 = 1/T1.''' % equation

  def getInitParams(self, xs, ys, weights=None):

    y_at_xmax = xmax = None
    for n, x in enumerate(xs):
      y = ys[n]
      if xmax is None or x > xmax:
        xmax = x
        y_at_xmax = y

    A = y_at_xmax
    B = 0.0

    return (A, B)

  def getValue(self, x, params):

    (A, B) = params
    y = A * (1 - 2*math.exp(-B*x))

    return y

  def getValueAndDerivative(self, x, params):

    (A, B) = params
    t = 2 * math.exp(-B*x)
    s = 1 - t
    y = A * s

    dy0 = s
    dy1 = A * x * t
    dy_dparams = (dy0, dy1)

    return (y, dy_dparams)

if __name__ == '__main__':

  import random

  fit = FitInversionRecovery()
  A = 1.0
  B = 1.0
  params = (A, B)
  noise = 0.03
  xs = [0.1*x for x in range(1,10)]
  ys = [fit.getValue(x, params) + noise*(2*random.random()-1) for x in xs]
  print 'Ain = %.2f, Bin = %.2f' % params
  print 'xs:', ','.join(['%.3f' % x for x in xs])
  print 'ys:', ','.join(['%.3f' % y for y in ys])

  (params, chisq, ysFit) = fit.calcFit(xs, ys, noise=noise)

  print 'ysFit:', ','.join(['%.3f' % y for y in ysFit])
  print 'Afit = %.2f, Bfit = %.2f' % tuple(params)
  print 'chisq = %.2f' % chisq

