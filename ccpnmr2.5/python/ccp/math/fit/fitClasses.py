from memops.math.fit.fitClasses import fitClasses, FitLinear

from ccp.math.fit.FitInversionRecovery import FitInversionRecovery
from ccp.math.fit.FitKdMonomerDimerFastExchange import FitKdMonomerDimerFastExchange
from ccp.math.fit.FitKdProteinLigandFastExchange import FitKdProteinLigandFastExchange
from ccp.math.fit.FitKdProteinLigandSlowExchange import FitKdProteinLigandSlowExchange

fitIntensityClasses = fitClasses[:]
fitIntensityClasses.extend([
  FitKdProteinLigandSlowExchange,  # x = (P, L), y = 1 - (P+L+A - sqrt((P+L+A)^2 - 4PL))/2P'
  FitInversionRecovery,            # y = A (1 - 2 exp(-Bx))
])

fitShiftClasses = [
  FitLinear,                       # y = Ax + B
  FitKdMonomerDimerFastExchange,   # y = A (1 + B/(4x) - sqrt((1+B/(4x))^2 - 1) - C)
  FitKdProteinLigandFastExchange,  # x = (P, L), y = A (P+L+B - sqrt((P+L+B)^2 - 4PL))/2P'
]

