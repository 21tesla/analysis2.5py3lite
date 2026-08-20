#!/usr/bin/env python
# encoding: utf-8
"""
PeakSeparatorPyMC.py

Created by Daniel O'Donovan on 2010-11-22.
Copyright (c) 2010 University of Cambridge. All rights reserved.
"""

import sys, math

try:
  import numpy as np
  sys.path.append( '/Users/djo35/local/fink/lib/python2.7/site-packages' )
  import pymc

except ImportError:
  print 'Python cannot find PyMC module'
  raise

  # if ndim == 1:
  #   @pymc.deterministic(plot=False)
  #   def peakModel(h=h, sx=s[0], mx=m[0]):
  #     return peakModelNonPyMC(params.peakShape, h, s, m, size, nsignals)
# def peakModelNonPyMC(shape, h, s, m, size, nsignals):
#   """ Peak model function - pymc free """
# 
#   ndim = len( size )
#   x = [ np.zeros(size[i], dtype=np.float32) for i in range(ndim) ]
# 
#   for i in range(ndim):
#     x_range = np.arange( size[i], dtype=np.float32 )
#     if   shape == 3:
#       """ Shape 3 - Gaussian """
#       for n in range( nsignals ):
# 
#         s_array = np.array( s[n][i] )
#         one_over_s_squared = np.array( (-1. / np.power( s_array, 2 ) ) )
#         print one_over_s_squared
#         m_array = np.array( m[n][i] )
#         x_minus_m_squared  = np.array( np.power( x_range - m_array, 2 ) )
#         print x_minus_m_squared
#         x[i] += [h[n] * math.exp(one_over_s_squared[ii] * x_minus_m_squared[ii]) for ii in range(size[i])]
#     
#     elif shape == 4:
#       """ Shape 4 - Lorentzian """
#       for n in range( nsignals ):
#         s_array = np.asarray( s[n][i] )
#         m_array = np.asarray( m[n][i] )
#         x[i] += h[n] * np.power( s_array, 2 ) / (np.power( x_range - m_array, 2 ) + np.power( s_array, 2 ))
# 
#   r = np.ones( shape=[1] )
#   for v in x:
#     r = np.kron( v, r )
#   r = r.reshape( size[::-1] )
#   return r

# @pymc.deterministic(plot=False)
# def peakModel(shape, size, h=h, sx=s[0], sy=s[1], mx=m[0], my=m[1]):
#   ndim = len( size )
#   x_range = np.arange( size[0], dtype=np.float32 )
#   y_range = np.arange( size[1], dtype=np.float32 )
#   data = [ np.zeros(size[i], dtype=np.float32) for i in range(ndim) ]
# 
#   if   shape == 3:
#     """ Shape 3 - Gaussian """
#     for n in range( nsignals ):
#       # OMG this is an annoying way to avoid crappy PyMC bugs
#       data[0] += h[n] * np.exp( (-1. / np.power( sx[n], 2 ) ) * np.power( x_range - mx[n], 2 ) )
#       data[1] += h[n] * np.exp( (-1. / np.power( sy[n], 2 ) ) * np.power( y_range - my[n], 2 ) )
# 
#   elif shape == 4:
#     """ Shape 4 - Lorentzian """
#     for n in range( nsignals ):
#       data[0] += h[n] * np.power( sx[n], 2 ) / (np.power( x_range - mx[n], 2 ) + np.power( sx[n], 2 ))
#       data[1] += h[n] * np.power( sy[n], 2 ) / (np.power( y_range - my[n], 2 ) + np.power( sy[n], 2 ))
# 
#   # this will work in N dim, but was getting other bugs
#   r = np.ones( shape=[1] )
#   for v in data:
#     r = np.kron( v, r )
#   r = r.reshape( size[::-1] )
#   return r

# for i in range(ndim):
#   for n in range(nsignals):
#     # sigma prior
#     @pymc.stochastic( name = 's_%1d_n%1d' % (i, n) )
#     def s_j( value=midSigma[i], lower=params.minSigma[i], upper=params.maxSigma[i] ):
#       return pymc.uniform_like( value, lower, upper )
#     s['%1d_n%1d' % (i, n)] = s_j
# 
#     # mean position prior
#     @pymc.stochastic( name = 'm_%1d_n%1d' % (i, n) )
#     def m_j( value=midPos[i], lower=0., upper=size[i] ):
#       return pymc.uniform_like( value, lower, upper )
#     m['%1d_n%1d' % (i, n)] = m_j


def PeakSeparatorPyMC( params ):

  if params.Ndim != 2:
    print 'PyMC only works with 2d currently'
    return

  ndim = params.Ndim
  size = params.sampleSize

  peakList = params.peakList
  
  # all this just to get the baseLevel (lowest contour level)
  dataSource        = peakList.dataSource
  block_file        = dataSource.block_file
  data              = block_file.getValues( params.sampleStart, params.sampleEnd )

  if params.minAtoms != params.maxAtoms:
    print '&&& PyMC cannot do RJ MCMC - have to specify exact number of peaks'
    print 'Peak separator assuming %d peaks (min peaks)' % params.minAtoms

  nsignals = params.minAtoms
  shape    = params.peakShape

  midSigma = [(max( params.minSigma[i], params.maxSigma[i]) - min( params.minSigma[i], params.maxSigma[i])) / 2. for i in range(ndim)]
  midPos   = [(max( 0., size[i]) - min( 0, size[i])) / 2. for i in range(ndim)]

  params.ClibKeys.sort()
  for key in params.ClibKeys:
    print key, params.__dict__[key]

  # These are all stochastic methods (height, sigma and mean position)
  # hyper-prior
  h  =  pymc.Uniform('h',        lower=params.minHeight,   upper=params.maxHeight,   size=nsignals )
  # priors
  s, m = {}, {}
  for i in range(ndim):
    s['%1d'%i] = pymc.Uniform('s_%1d'%i, lower=0.1, upper=32.,    size=nsignals )
    m['%1d'%i] = pymc.Uniform('m_%1d'%i, lower=0., upper=size[i], size=nsignals )

  x_range = [np.arange( size[i], dtype=np.float32 ) for i in range(ndim) ]

  # likelihood
  @pymc.deterministic(plot=False)
  def peakModel(h=h, s=s, m=m):

    model = [ np.zeros(size[i], dtype=np.float32) for i in range(ndim) ]

    for n in range( nsignals ):

      for i in range(ndim):

        ss      = s['%1d' % i][n]
        mm      = m['%1d' % i][n]

        if   shape == 3:
          # Shape 3 - Gaussian
          model[i] += np.exp( (-1. / np.power( ss, 2 ) ) * np.power( x_range[i] - mm, 2 ) )

        elif shape == 4:
          # Shape 4 - Lorentzian
          model[i] += np.power( ss, 2 ) / (np.power( x_range[i] - mm, 2 ) + np.power( ss, 2 ))

        else:
          print 'Unknown shape %d' % i
          return

      model[0] *= h[n]

    # r = np.outer( model[1], model[0] )
    r = np.ones( shape=[1] )
    for v in model: r = np.kron( v, r )
    r = r.reshape( size[::-1] )
    return r

  D = pymc.Normal('D', mu=peakModel, tau=2, value=data, observed=True)

  mcmcVars = [h, s, m, peakModel, D]

  M = pymc.MCMC( mcmcVars )
  M.sample( iter=5.0E4, burn=1.0E4, thin=1.0E1 )

  results = []

  for n in range(nsignals):

    peak = [None, nsignals]

    peak.append( h.stats()['mean'][n] if nsignals > 1 else h.stats()['mean'] )

    for i in range( ndim ):
    # for i in range( ndim-1, -1, -1 ):
      peak.append( m['%1d' % i].stats()['mean'][n] if nsignals > 1 else m['%1d' % i].stats()['mean'] )
      peak.append( s['%1d' % i].stats()['mean'][n] if nsignals > 1 else s['%1d' % i].stats()['mean'] )

    results.append( peak )

    # height    =   float(sample[2])
    # sigma     = [ float(sample[2*i+2 + 2]) for i in range(tempNdim)]
    # position  = [ float(sample[2*i+1 + 2]) for i in range(tempNdim)]

  return results

