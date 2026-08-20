#!/usr/bin/env python
# encoding: utf-8
"""
PeakSeparatorPeakList.py

Created by Daniel O'Donovan on 2010-11-10.
Copyright (c) 2010 University of Cambridge. All rights reserved.

The intention of this file is to update all the Peak Separator 
params that rely on the chosen peak list:

dataFile     Ndim         isBigEndian  nPoints      blockSize
dimWrapped   isFreqDim    dataDimRefs  minHeight

"""

from ccp.api.nmr.Nmr                        import FreqDataDim
from ccpnmr.analysis.core.ExperimentBasic   import getPrimaryDataDimRef


def getPeakListParams( params ):
  """ Get all the parameters needed from the peak list """

  if not params.peakList:
    return

  # These are just handy
  peakList          = params.peakList
  dataSource        = peakList.getDataSource()
  analysisSpectrum  = dataSource.analysisSpectrum
  analysisProject   = analysisSpectrum.analysisProject

  # Now setting parameters (from Peak List)
  params.dataDimRefs = []
  params.isFreqDim   = []

  params.Ndim       = dataSource.numDim                           # set Ndim
  params.dataFile   = dataSource.dataStore.fullPath               # set dataFile

  for i, dataDim in enumerate( dataSource.sortedDataDims() ):

    if isinstance( dataDim, FreqDataDim ):
      dataDimRef = getPrimaryDataDimRef( dataDim )
      params.dataDimRefs.append( dataDimRef )                     # set dataDimRefs
      params.isFreqDim.append( True )                             # set isFreqDim
    else:
      params.dataDimRefs.append( None )                           # set dataDimRefs
      params.isFreqDim.append( False )                            # set isFreqDim

  params.minHeight        = min(analysisSpectrum.posLevels) * analysisProject.getGlobalContourScale()
                                                                  # set minHeight
  params.nPoints          = list(dataSource.dataStore.numPoints)  # set nPoints
  params.blockSize        = list(dataSource.dataStore.blockSizes) # set blockSize
  params.dimWrapped       = [0] * params.Ndim                     # ask wayne about this one...
  params.isBigEndian      =  1 if dataSource.dataStore.isBigEndian else 0
                                                                  # set isBigEndian


