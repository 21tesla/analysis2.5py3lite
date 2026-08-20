import math

from memops.math.fit.FitNonlinear import FitNonlinear

LARGE_FLOAT = 1.0e10

class FitKdProteinLigandFastExchange(FitNonlinear):

  name = 'Kd for Protein-Ligand in Fast Exchange'

  paramNames = ['boundProteinShift-freeProteinShift', 'Kd']

  equation = 'y = A (P+L+B - sqrt((P+L+B)^2 - 4PL))/2P'

  description = '''This is for the situation where you have a ligand
binding to a protein in fast exchange and are measuring the shift for
a given peak group for varying amounts of ligand and/or protein
concentration.  It is assumed that there is a known estimate for the
shift when there is no ligand (so free protein) but that there is no
known estimate for the shift when all the protein is bound.  It
produces an estimate for Kd, and for the shift of the bound protein.
The equation used is %s.  Here P is the total protein concentration
(which often does not vary between experiments but might), L is the
total ligand concentration (we take x = (P, L)), y is the observed shift
difference between the given experiment and the free protein shift,
A is the difference between the bound and free protein shifts and B
is the Kd.''' % equation

  def getInitParams(self, xs, ys, weights=None):

    A = max(ys)
    B = 0.5 * xs[0][0]  # half the protein concentration of the first experiment

    return (A, B)

  def getValue(self, x, params):

    (P, L) = x
    (A, B) = params
    t = P + L + B
    s = t*t - 4*P*L
    s = math.sqrt(s)
    y = A * (t - s) / (2*P)

    return y

  def getValueAndDerivative(self, x, params):

    (P, L) = x
    (A, B) = params
    if B < 0:
      y = LARGE_FLOAT
      v = 0
      w = -LARGE_FLOAT
    else:
      t = P + L + B
      s = t*t - 4*P*L
      s = math.sqrt(s)
      v = (t - s) / (2*P)
      y = A * v
      w = -y / s

    dy_dparams = (v, w)

    return (y, dy_dparams)

if __name__ == '__main__':

  import random

  fit = FitKdProteinLigandFastExchange()
  A = 1.0
  B = 1.2
  params = (A, B)
  noise = 0.01
  xs = [(1.0, 0.1*x) for x in range(1,10)]
  ys = [fit.getValue(x, params) + noise*(2*random.random()-1) for x in xs]
  print 'Ain = %.2f, Bin = %.2f' % params
  print 'xs:', ','.join(['%.3f' % x[1] for x in xs])
  print 'ys:', ','.join(['%.3f' % y for y in ys])

  (params, chisq, ysFit) = fit.calcFit(xs, ys, noise=noise)

  print 'ysFit:', ','.join(['%.3f' % y for y in ysFit])
  print 'Afit = %.2f, Bfit = %.2f' % tuple(params)
  print 'chisq = %.2f' % chisq

