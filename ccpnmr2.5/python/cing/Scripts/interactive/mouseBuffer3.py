# python -u $CINGROOT/python/cing/Scripts/interactive/mouseBuffer3.py

try:
    from matplotlib.pyplot import *
except Exception:
    pass

from numpy import ones  # no longer re-exported by modern pyplot
from numpy.random import rand

import cing
from cing.Libs.NTutils import *  #@UnusedWildImport

cing.verbosity = 9

# fake up some data
spread= rand(50) * 100
spread.sort()
center = ones(25) * 50
flier_high = rand(10) * 100 + 100
flier_low = rand(10) * -100
#data =concatenate((spread, center, flier_high, flier_low), 0)
data = spread

# basic plot
wisk_lo, wisk_hi = boxplot(data)
print("wisk_lo, wisk_hi %s %s" % (wisk_lo, wisk_hi))

savefig('/Users/jd/tmp/cingTmp/boxplot_simple.png')


