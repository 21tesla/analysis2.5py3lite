from memops.math.fit.FitLinear import FitLinear
from memops.math.fit.FitLogLinear import FitLogLinear
from memops.math.fit.FitExp import FitExp
from memops.math.fit.FitExpWithBaseline import FitExpWithBaseline

fitClasses = [
  FitLinear,                   # y = Ax + B
  FitLogLinear,                # log(y) = log(A exp(-Bx))
  FitExp,                      # y = A exp(-Bx)
  FitExpWithBaseline,          # y = A exp(-Bx) + C
]

if __name__ == '__main__':

  import random

  # FitLinear test

  sigma = 0.1
  A = 1.0
  B = 2.5
  xs = [1.0, 2.5, 3.7, 8.9, 11.2]
  ys = [A*x + B + sigma*(2*random.random()-1) for x in xs]

  fit = FitLinear()
  ((a,b), chisq, ysFit, (aDev,bDev)) = fit.calcFit(xs, ys, findDev=True)

  print 'FitLinear'
  print 'Input params: %.2f %.2f' % (A, B)
  print 'Fit params: %.2f %.2f' % (a, b)

