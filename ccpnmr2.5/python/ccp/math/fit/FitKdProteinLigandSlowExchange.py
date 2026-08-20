import math

from memops.math.fit.FitNonlinear import FitNonlinear

LARGE_FLOAT = 1.0e10

class FitKdProteinLigandSlowExchange(FitNonlinear):

  name = 'Kd for Protein-Ligand in Slow Exchange'

  paramNames = ['Kd']

  equation = 'y = 1 - (P+L+A - sqrt((P+L+A)^2 - 4PL))/2P'

  description = '''This is for the situation where you have a ligand
binding to a protein in slow exchange and have two peaks, one for the
free protein and one for the bound protein.  The fitting here is for
the case where you are considering the ratio of intensities between
the free protein peak intensity in a given experiment and the protein
peak intensity in the case when there is no ligand.  The fit produces
an estimate for Kd.  The equation used is %s.
Here P is the total protein concentration (which often does not vary
between experiments but might), L is the total ligand concentration
(we take x = (P, L)), y is the observed intensity ratio, and A is the
Kd.''' % equation

  def getInitParams(self, xs, ys, weights=None):

    A = 0.5 * xs[0][0]  # half the protein concentration of the first experiment

    return (A,)

  def getValue(self, x, params):

    (P, L) = x
    (A,) = params
    t = P + L + A
    s = t*t - 4*P*L
    s = math.sqrt(s)
    y = 1 - (t - s)/(2*P)

    return y

  def getValueAndDerivative(self, x, params):

    (P, L) = x
    (A,) = params
    if A < 0:
      y = LARGE_FLOAT
      v = -LARGE_FLOAT
    else:
      t = P + L + A
      s = t*t - 4*P*L
      s = math.sqrt(s)
      w = (t - s) / (2*P)
      y = 1 - w
      v = w / s

    dy_dparams = (v,)

    return (y, dy_dparams)

if __name__ == '__main__':

  import random

  fit = FitKdProteinLigandSlowExchange()
  A = 1.0
  params = (A,)
  noise = 0.01
  xs = [(1.0, 0.1*x) for x in range(1,10)]
  ys = [fit.getValue(x, params) + noise*(2*random.random()-1) for x in xs]
  print 'Ain = %.3f' % params
  print 'xs:', ','.join(['%.3f' % x[1] for x in xs])
  print 'ys:', ','.join(['%.3f' % y for y in ys])

  (params, chisq, ysFit) = fit.calcFit(xs, ys, noise=noise)

  print 'ysFit:', ','.join(['%.3f' % y for y in ysFit])
  print 'Afit = %.3f' % tuple(params)
  print 'chisq = %.3f' % chisq

