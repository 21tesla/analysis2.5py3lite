#!/usr/bin/env python
# encoding: utf-8
"""
PeakSeparator.py

Created by Daniel O'Donovan on 2008-09-22.
Copyright (c) 2008 University of Cambridge. All rights reserved.
"""

from ccp.api.nmr.Nmr                        import FreqDataDim, SampledDataDim

from ccpnmr.analysis.core.ExperimentBasic   import getPrimaryDataDimRef
from ccpnmr.analysis.core.PeakBasic         import pickPeak, searchPeaks, setManualPeakIntensity
from ccpnmr.analysis.core.UnitConverter     import ppm2pnt, pnt2ppm, pnt2hz, hz2pnt
from ccpnmr.analysis.core.Util              import getMethod

from ccp.api.nmr.Nmr                        import FreqDataDim

from memops.gui.MessageReporter             import showError, showWarning

try:
  from cambridge.c                          import BayesPeakSeparator
except ImportError:
  print 'Error, cannot import BayesPeakSeparator - peak separation will not work.'

from cambridge.bayes.kmeans import kMeans

import math
import os, sys


## Possibly a very poor way to interpret results...
# def getPeaksFromResults( results, verbose=False ):
#   """ Given BayeSys result list, pull out salient peak info """
#   print 'BAD getPeaksFromResults !!! only use for PyMC !!!'
#   npeaks = int( results[-1][1] )
# 
#   peaks = [ results[-(i+1)] for i in range( npeaks )]
# 
#   print peaks
# 
#   return peaks

def getPeaksFromResults( results, verbose=False ):
  """ Given BayeSys result list, pull out salient peak info """

  import numpy as np

  results = np.asarray( results )

  npeaks = int( results[:,1].mean() )

  if npeaks == 1:
    # nice and simple
    peak = [ results[:,i].mean() for i in range( len(results[0]) ) ]
    return [peak]

  else:

    import pickle
    pickle.dump( results, open('results.pkl', 'wb') )

    good_results_i = np.nonzero( results[:,1] == npeaks )[0]
    good_results = results[ good_results_i ][:,2:]

    ( centres, clusters ) = kMeans( good_results, npeaks )

    peaks = []

    for n in range( npeaks ):
      c = np.asarray( clusters[n] )
      peak =  [None, None]
      peak += [ c[:,i].mean() for i in range( len(c[0]) ) ]
      peaks.append( peak )

      print 'peak %2d height %f stdev %f ' % (n, c[:,2].mean(), c[:,2].std())
      open( 'note.csv', 'a' ).write( '%d, %d, %f, %f\n' % ( n, npeaks, c[:,2].mean(), c[:,2].std()) )

    return peaks


def SeparatePeakRoutine(params, peakList, routine='bayesys'):
  """ Run Peak Separate Code and then add resultant peaks to CCPN PeakList (peakList) """

  # check that all C library parameters have been set up before proceeding
  fail = False
  for key in params.ClibKeys:
    if params.__dict__[key] == None:
      print '&&& C Library parameter %s not set, aborting.' % ( key )
      fail = True

  if fail: return
  del( fail )

  # As the Lorentz/Gaussian model not in use
  (params.maxQ, params.minQ) = (1., 1.)
  
  # params.ClibKeys.sort()
  # for key in params.ClibKeys:
  #   print key, params.__dict__[key]

  if params.maxHeight <= params.minHeight:
    print '&&& Peak Separator height mismatch - please report to CCPN'
    return

  if routine == 'bayesys':

    # Maybe a quick parameter check so that no one is doing anything silly ? 
    # Run the BayesPeakSeparator function in C
    results = BayesPeakSeparator.run_bayes(   \
                                    params.dataFile, params.Ndim, params.isBigEndian, params.nPoints,           #  4
                                    params.blockSize,                                                           #  1
                                    params.sampleStart, params.sampleEnd,                                       #  2
                                    params.maxHeight, params.minHeight,                                         #  2
                                    params.maxSigma, params.minSigma,                                           #  2
                                    params.maxQ, params.minQ,                                                   #  2
                                    params.dimWrapped,                                                          #  1
                                    params.peakShape, params.positivePeaks, params.minAtoms, params.maxAtoms,   #  2
                                    params.rate)                                                                #  1
                                                                                                                # 17

  elif routine == 'pymc':

    try:
      from PeakSeparatorPyMC import PeakSeparatorPyMC
    except ImportError:
      print 'Error, cannot import PeakSeparatorPyMC - PyMC peak separation will not work.'
      return

    if params.maxHeight <= params.minHeight:
      print '&&& Peak Separator height mismatch - please report to CCPN'
      return

    if params.Ndim != 2:
      print '&&& Peak Separator PyMC only in two dims currently'
      return

    results = PeakSeparatorPyMC( params )

  if (results == None) or (len(results) == 0):
    print '&&& SeparatePeakRoutine: failed.', results
    return 1

  shapeDict = { 3:'Gaussian', 4:'Lorentzian', 5:'Other' }

  peak_count = 0

  # tempNdim as nDim not correct for pseudo nD spectra
  sampledData = True if False in set( [bool(t) for t in params.isFreqDim] ) else False
  tempNdim = (params.Ndim - 1) if sampledData else params.Ndim

  # for each peak in the returned sample set
  for sample in getPeaksFromResults( results ):

    height    =   float(sample[2])
    sigma     = [ float(sample[2*i+2 + 2]) for i in range(tempNdim)]
    position  = [ float(sample[2*i+1 + 2]) for i in range(tempNdim)]

    peak_count += 1

    cluster_peak_position_point = [0] * params.Ndim

    if params.peakShape == 3: # Gaussian (nd)
      cluster_peak_volume         = height
    if params.peakShape == 4: # Lorentzian (nd)
      cluster_peak_volume         = height
    else:                 # Other (nd)
      cluster_peak_volume         = height

    for i in range( tempNdim ):

      # No longer a cluester peak! cluster_peak_position_point
      cluster_peak_position_point[i]  = float( params.sampleStart[i] ) + position[i] + 1.0

      # Gauss Volume 
      if params.peakShape == 3:
        cluster_peak_volume            *= sigma[i] * math.sqrt( 2.0 * math.pi )

      # Lorentz Volume
      elif params.peakShape == 4:
        cluster_peak_volume            *= sigma[i] * math.pi

      else:
        print '&&& Write code to integrate this shape: %d' % (params.peakShape)

    # plane dim for pseudo 3d Analysis data (titrations / time lapse etc)
    if sampledData:
      for i in range(params.Ndim):
        if not params.isFreqDim[i]:
          cluster_peak_position_point[i] = params.sampleEnd[i]

    # pickPeak from ccpnmr.analysis.core.PeakBasic
    peakListPeak =  pickPeak(peakList, cluster_peak_position_point, unit='point', doFit=False)

    for i, peakDim in enumerate( peakListPeak.sortedPeakDims() ):

      if (sampledData) and (i == (params.Ndim-1)): break

      cluster_peak_sigma              = sigma[i]

      if peakDim.dataDimRef:
        a = pnt2hz( peakDim.position + cluster_peak_sigma, peakDim.dataDimRef )
        b = pnt2hz( peakDim.position,                      peakDim.dataDimRef )

      peakDim.setLineWidth( 2. * abs( a - b ) )

    # Add note to peak in CCPN peak list (so we know it's a BayeSys peak)
    peakListPeak.setDetails( 'PeakSeparator %s' % shapeDict[params.peakShape])
    if sample[0] < 0.:
      merit = 0.
    else:
      merit = 1. / (1. + math.exp( - sample[0] ))
    peakListPeak.setFigOfMerit( merit )

    # Add peak intensities (which have just been calculated)
    setManualPeakIntensity( peakListPeak,              height, intensityType='height' )
    setManualPeakIntensity( peakListPeak, cluster_peak_volume, intensityType='volume' )


###############################################################################
### These functions are for repicking existing peak lists 
def getNeighbourPoints( point, ndim ):
  """ Get neighbour points for given point, up to ndim """

  point = [ int(p) for p in point ]

  points    = []

  for sign in [-1, +1]:

    identity = [ [0.] * ndim for i in range(ndim) ]
    for i in range(ndim): identity[i][i] = sign

    for x in range( ndim ):

      neighbour = identity[x]

      newPoint  = [neighbour[i] + point[i] for i in range( ndim )]

      if len(point) != ndim:
        for i in range( ndim, len(point) ):
          newPoint.append( point[i] )

      points.append( newPoint )

  return points

def getSearchRegion( peak, baseLevel, block_file, addOne=False, isFreqDim=None, ndim=2 ):
  """ Walk around peak to find region containing all peaks in same contoured region
      - and return that region (for searching)
  """

  start = [ dim.position for dim in peak.sortedPeakDims() ]

  # list of points to be checked
  pointsToSearch = getNeighbourPoints( start, ndim )

  # keep track of the good points (above base level)
  pointsAboveLevel = []
  pointsChecked = set()

  # when list of points to be checked is [] we have all points above level
  while pointsToSearch:

    pointToCheck = pointsToSearch.pop()

    pointsChecked.add( tuple(pointToCheck) )

    bfilePointToCheck = pointToCheck[:]
    for dim in range( len( bfilePointToCheck ) ): bfilePointToCheck[dim] = int(bfilePointToCheck[dim]) - 1

    if block_file.getPointValue( tuple(bfilePointToCheck) ) > baseLevel:
      pointsAboveLevel.append( pointToCheck )

      for newPeak in getNeighbourPoints( pointToCheck, ndim ):
        if tuple(newPeak) not in pointsChecked:
          pointsToSearch.append( newPeak )

  # set ridiculous max and min that will be reset
  (maxPoint, minPoint) = ([0] * len(start), [99999] * len(start))
  maxPoint = [    0 for i in range(ndim)]
  minPoint = [99999 for i in range(ndim)]

  # if we have a SampleDataDim, append that
  for i in range( ndim, len(start) ):
    maxPoint.append( int(start[i]) )
    minPoint.append( int(start[i]) )

  for point in pointsAboveLevel:
    for dim in range( ndim ):
      if point[dim] > maxPoint[dim]:
        maxPoint[dim] = point[dim]
      if point[dim] < minPoint[dim]:
        minPoint[dim] = point[dim]

  if addOne:
    if isFreqDim:
      maxPoint = [maxPoint[i] + 1 if isFreqDim[i] else maxPoint[i] for i in range( len(maxPoint) ) ]
      minPoint = [minPoint[i] - 1 if isFreqDim[i] else minPoint[i] for i in range( len(maxPoint) ) ]
    else:
      maxPoint = [maxPoint[i] + 1 for i in range( len(maxPoint) ) ]
      minPoint = [minPoint[i] - 1 for i in range( len(maxPoint) ) ]

  return ( tuple(maxPoint), tuple(minPoint) )

def SeparatePeaksInPeakList( params, HEIGHT_MULTIPLIER=2.5 ):
  """ Given a peak list (params.peakList) repick all peaks into a new list """

  if not params.peakList:
    print '&&& No Peak List set'
    return

  peakList = params.peakList
  
  # all this just to get the baseLevel (lowest contour level)
  dataSource        = peakList.dataSource
  analysisSpectrum  = dataSource.analysisSpectrum
  analysisProject   = analysisSpectrum.analysisProject
  baseLevel         = min(analysisSpectrum.posLevels) * analysisProject.getGlobalContourScale()
  block_file        = dataSource.block_file

  peaksToSearch = [ peak for peak in peakList.sortedPeaks() ]

  # tempNdim as nDim not correct for pseudo nD spectra
  sampledData = True if False in set( [bool(t) for t in params.isFreqDim] ) else False
  tempNdim = (params.Ndim - 1) if sampledData else params.Ndim

  ndim = params.Ndim

  regionsList = []

  # for peak in list of peaks to be searched
  while peaksToSearch:

    sampleNdim = ndim - 1 if sampledData else ndim

    peak    = peaksToSearch.pop()
    region  = getSearchRegion( peak, baseLevel, block_file, ndim=sampleNdim, isFreqDim=params.isFreqDim, addOne=True )

    regionsList.append( region )

  print '\nRepick Peak List - Found %6d regions' % (len(regionsList))

  # remove all but unique regions
  regionsList = list(set(regionsList))

  print 'of which %6d are unique' % (len(regionsList))

  repickList = []

  for region in regionsList:

    ppmRegion       = [ [0.] * ndim for i in range(2) ] # [min x,y,z], [max x,y,z] etc.
    searchPpmRegion = [ [0.] * 2 for i in range(ndim) ] # [min, max] * x, y, z etc.

    # set ppmRegion
    for i in range(ndim):
      if params.isFreqDim[i]:
        # max ppmRegion (ppm is backwards remember)
        maxPpm                = pnt2ppm( region[0][i], params.dataDimRefs[i] )
        minPpm                = pnt2ppm( region[1][i], params.dataDimRefs[i] )
      else:
        # max search region
        maxPpm = region[0][i] - 0.5
        # min search region
        minPpm = region[1][i] + 0.5

      ppmRegion[1][i]       = minPpm
      ppmRegion[0][i]       = maxPpm

      searchPpmRegion[i][1] = minPpm
      searchPpmRegion[i][0] = maxPpm


    peaksInRegion = searchPeaks( [peakList], searchPpmRegion )

    if len( peaksInRegion ) == 0:
      print '&&& Attempted repicking of region with no peaks! ', peaksInRegion
      continue

    repickList.append( [region, ppmRegion, peaksInRegion] )

  newPeakList = peakList.dataSource.newPeakList()
  newPeakList.details = 'Peak Separator rejiggled peak list'

  params.peakList       = newPeakList

  block_file = peakList.dataSource.block_file

  # re-pick into the new peak list based on regions from before
  for ii, repick in enumerate(repickList):

    print 'Repicking region %3d of %3d' % (ii+1, len(repickList))

    region    = repick[0]
    ppmRegion = repick[1]
    npeaks    = int(len(repick[-1]))

    first   = [ int(region[1][i]) for i in range(ndim) ]
    last    = [ int(region[0][i]) for i in range(ndim) ]

    if not params.isFreqDim[-1]:
      first[-1] -= 1

    if npeaks == 0:
      print '&&& Zero peaks found in region, skipping'
      continue

    # these values (may) differ for each run
    params.minAtoms       = npeaks
    params.maxAtoms       = npeaks

    params.maxHeight      = HEIGHT_MULTIPLIER * max( abs( block_file.maxValue(first, last) ), abs( block_file.minValue(first, last) ) )

    params.sampleStart    = first
    params.sampleEnd      = last
    params.sampleSize     = [max( 1, abs( params.sampleEnd[i] - params.sampleStart[i] ) ) for i in range(ndim)]

    params.samplePpmStart = []
    params.samplePpmEnd   = []
    params.samplePpmSize  = []

    # temp just grab the data
    # sampDataDim = dataSource.sortedDataDims()[-1].getPointValues()[last[-1]-1]
    # dump = { 'npeaks':npeaks, 'samp':sampDataDim, 'size':params.sampleSize, 'data':block_file.getValues(first, last) }
    # pickleName = 'double_data_%5.4f.pkl' % sampDataDim
    # print 'pickling ', pickleName
    # import pickle
    # pickle.dump( dump, open(pickleName, 'wb') )

    SeparatePeakRoutine( params, newPeakList )


if __name__ == '__main__':
  print 'PeakSeparator - This should be run from inside '
