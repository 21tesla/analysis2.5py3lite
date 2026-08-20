import math

from memops.math.fit.FitNonlinear import FitNonlinear

SMALL_X = 1.0e-7

class FitKdMonomerDimerFastExchange(FitNonlinear):

  name = 'Kd for Monomer-Dimer in Fast Exchange'

  paramNames = ['dimerShift-monomerShift', 'Kd', 'monomerShift']

  equation = 'y = A(B+4x-sqrt((B+4x)^2-16x^2))/4x + C'
  description = '''This is for the situation where you have a monomer
and dimer in fast exchange and are measuring the shift for a given peak
group for varying amounts of total concentration.  It is assumed that there
is *not* a known estimate of either the shift for the monomer or the dimer.
It produces an estimate for Kd, and for the shifts of the monomer and
dimer.  The equation used is %s.  Here x is the total concentration, y
is the observed shift, A is the difference between the dimer and the
monomer shifts, B is the Kd and C is the monomer shift.''' % equation

  def getInitParams(self, xs, ys, weights=None):

    # C is limit of y as x --> 0 so take it to be y value at minimum x
    # A+C is limit of y as x --> infty so take it to be y value at maximum x
    # maximum y occurs at Bx around 4.5

    y_at_xmin = xmin = None
    y_at_xmax = xmax = None
    x_at_ymax = ymax = None
    for n, x in enumerate(xs):
      y = ys[n]
      if xmin is None or x < xmin:
        xmin = x
        y_at_xmin = y
      if xmax is None or x > xmax:
        xmax = x
        y_at_xmax = y
      if ymax is None or y > ymax:
        ymax = y
        x_at_ymax = x

    C = y_at_xmin
    A = y_at_xmax - C
    if x_at_ymax > 0:
      B = 4.5 / x_at_ymax
    else:
      B = 1.0

    return (A, B, C)

  def getValue(self, x, params):

    (A, B, C) = params
    if x < SMALL_X:
      y = - A * C
    else:
      t = 1 + B/(4*x)
      y = A * (t - math.sqrt(t-1) - C)

    return y

  def getValueAndDerivative(self, x, params):

    (A, B, C) = params

    x *= B
    x = abs(x)

    if x < SMALL_X:
      s = x * x / 6
      y = A*s + C
      dy0 = s
      dy1 = 2 * A * s / B
    else:
      s = math.sin(x) / x
      t = math.cos(x)
      y = A*(1-s) + C
      dy0 = 1.0 - s
      dy1 = A * (s - t) / B

    dy_dparams = (dy0, dy1, 1.0)

    return (y, dy_dparams)

if __name__ == '__main__':

  import random

  fit = FitKdMonomerDimerFastExchange()
  A = 1.0
  B = 1.0
  C = 1.0
  params = (A, B, C)
  noise = 0.01
  xs = [0.1*x for x in range(1,10)]
  ys = [fit.getValue(x, params) + noise*(2*random.random()-1) for x in xs]
  print 'Ain = %.2f, Bin = %.2f, Cin = %.2f' % params
  print 'xs:', ','.join(['%.3f' % x for x in xs])
  print 'ys:', ','.join(['%.3f' % y for y in ys])

  (params, chisq, ysFit) = fit.calcFit(xs, ys, noise=noise)

  print 'ysFit:', ','.join(['%.3f' % y for y in ysFit])
  print 'Afit = %.2f, Bfit = %.2f, Cfit = %.2f' % tuple(params)
  print 'chisq = %.2f' % chisq

