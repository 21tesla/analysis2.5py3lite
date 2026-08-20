from memops.math.fit.FitClass import FitClass
from memops.math.fit.linearFit import linearFit

class FitLinear(FitClass):

  name = 'Linear'

  paramNames = ['A', 'B']

  equation = 'y = A x + B'

  description = '''This does a least squares fit on the y values to
determine the best linear fit, %s.''' % equation

  def getValue(self, x, params):

    (A, B) = params
    y = A*x + B

    return y

  def calcFit(self, xs, ys, weights=None, noise=None, params=None, devMethod=None):

    findDev = (devMethod == 'normal')
    result = linearFit(xs, ys, weights, findDev)

    if devMethod == 'bootstrap':
      params = result[0]  # start with given answer
      paramsDev = self.bootstrapDev(xs, ys, weights, noise, params)
      result.append(paramsDev)
    
    return result

