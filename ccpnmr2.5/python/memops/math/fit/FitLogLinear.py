from memops.math.fit.FitClass import FitClass
from memops.math.fit.logLinearFit import logLinearFit

class FitLogLinear(FitClass):

  name = 'Log Linear'

  paramNames = ['A', 'B']

  equation = 'log y = log(A exp(-Bx))'

  description = '''This fits log y versus log x, using a linear, least
squares fit on the log y values, to determine the best fit, %s.''' % equation

  def getValue(self, x, params):

    (A, B) = params
    y = A * exp(-B*x)

    return y

  def calcFit(self, xs, ys, weights=None, noise=None, params=None, devMethod=None):

    findDev = (devMethod == 'normal')
    result = logLinearFit(xs, ys, weights, findDev)

    if devMethod == 'bootstrap':
      params = result[0]  # start with given answer
      paramsDev = self.bootstrapDev(xs, ys, weights, noise, params)
      result.append(paramsDev)
    
    return result
