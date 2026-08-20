#!/usr/bin/env python
# encoding: utf-8
"""
PeakSeparatorRegion.py

Created by Daniel O'Donovan on 2010-11-10.
Copyright (c) 2010 University of Cambridge. All rights reserved.

The intention of this file is to update all the Peak Separator 
params that rely on the user chosen region:

sampleStart,    sampleEnd,   sampleSize
maxHeight,

"""

from ccp.api.general.DataLocation       import NumericMatrix
from ccpnmr.analysis.core.UnitConverter import ppm2pnt
from memops.gui.MessageReporter         import showError, showWarning

from PeakSeparatorPeakList              import getPeakListParams

def getRegionParams(params, argServer=None, HEIGHT_MULTIPLIER=2.5):
  """ given region in spec return params """

  # general variables 
  xyz = ['x', 'y', 'z1', 'z2', 'z3']

  if argServer:

    # information from user dragged region
    analysisProject     = argServer.getAnalysisProject()
    window              = argServer.getCurrentWindow()
    ppmPoints           = argServer.parent.currentRegion

    # information from selected peak list
    peakList            = params.peakList
    dataSource          = peakList.dataSource
    analysisSpectrum    = analysisProject.findFirstAnalysisSpectrum(dataSource=dataSource)
    windowPane          = window.findFirstSpectrumWindowPane()
    spectrumWindowView  = windowPane.findFirstSpectrumWindowView(analysisSpectrum=analysisSpectrum)

    if not spectrumWindowView:
      showError( 'Peak Separator', 'The peak list you have chosen is not valid for the current window. Have you selected the correct peak list?')
      return

    if not (spectrumWindowView.isPosVisible or spectrumWindowView.isNegVisible):
      showWarning( 'Peak Separator', 'The spectrum you are picking is not visible in the current window. Have you selected the correct peak list?')

  else:
    print '&&& getRegionParams: No arg server'
    return

  # make sure that dragged region and peak list match!
  if not spectrumWindowView:
    peakListName = '%s:%s:%s' % (spectrum.experiment.name, spectrum.name, peakList.serial)
    showWarning( 'Incorrect Peak List', 
      "Peak list '%s' doesn't match selected spectral region, please select another peak list." % peakListName )
    return

  axisPanel   = windowPane.findFirstAxisPanel(label='z1')
  if axisPanel:
    axisRegion  = axisPanel.findFirstAxisRegion()
    zRegion     = axisRegion.region
  else:
    zRegion = None

  if not isinstance( analysisSpectrum.dataSource.dataStore, NumericMatrix ):
    showError('Not NumericMatrix', 'Peak Separator cannot handle this type of spectra.')
    return

  # This may or may not have been run already - for safety lets do it again
  # set the parameters that are derrived from the peak list
  getPeakListParams( params )

  axisMapping     = [ spectrumWindowView.findFirstAxisMapping( label=xyz[i] ) for i in range(params.Ndim) ]

  analysisDataDim = [ axisMapping[i].analysisDataDim  for i in range(params.Ndim) ]
  dataDim         = [ analysisDataDim[i].dataDim      for i in range(params.Ndim) ]

  localdataDimRef = []

  dataDimOrder    = [dataDim[i].getDim() for i in range(params.Ndim)]

  for i in range( params.Ndim ):

    # tricky ness for kinetic data
    if params.isFreqDim[i]:
      localdataDimRef.append( dataDim[i].findFirstDataDimRef() )
    else:
      localdataDimRef.append( None )

  localPpmPairs   = []

  # extend ppmPoints from z1, z2 etc planes (0,2,1,3,4,5, etc -> x0,x1,y0,y1,z10,z11, etc)
  for i in range( params.Ndim ):
    if i < 2:   # get the dragged region in x,y from argserver (reorder)
      localPpmPairs.append( [max( [ ppmPoints[i], ppmPoints[i+2] ] ), min( [ ppmPoints[i], ppmPoints[i+2] ] )] )
    else:       # get any extra dims z1, z2 from axispanel
      localPpmPairs.append( [ max( zRegion ), min( zRegion ) ] )

  dataDimRef      = []
  PpmPairs        = []

  # re-order local ppm points into actual dim. ordering
  for i in range( params.Ndim ):

    PpmPairs.append( localPpmPairs[ dataDimOrder[i] - 1 ] )

    if localdataDimRef[ dataDimOrder[i] - 1 ]:
      dataDimRef.append( localdataDimRef[ dataDimOrder[i] - 1 ] )
    else:
      dataDimRef.append( None )

  sampDataDim = dataSource.sortedDataDims()

  params.dataDimRef     = dataDimRef[:]

  params.sampleStart    = []
  params.sampleEnd      = []

  params.samplePpmStart = []
  params.samplePpmEnd   = []

  for i in range( params.Ndim ):

    if dataDimRef[i]:

      params.samplePpmStart.append( PpmPairs[i][0] )                              # set ppmStart
      params.samplePpmEnd.append(   PpmPairs[i][1] )                              # set ppmEnd

      params.sampleStart.append( int(ppm2pnt(PpmPairs[i][0], dataDimRef[i] )) )   # set start
      params.sampleEnd.append(   int(ppm2pnt(PpmPairs[i][1], dataDimRef[i] )) )   # set end

    else:

      params.samplePpmStart.append( sampDataDim[i].getPointValues()[ int(PpmPairs[i][0]) -1 ] )
      params.samplePpmEnd.append(   sampDataDim[i].getPointValues()[ int(PpmPairs[i][1]) -1 ] )

      params.sampleStart.append( int( PpmPairs[i][0] ) -1 )                          # set start
      params.sampleEnd.append(   int( PpmPairs[i][1] )    )                          # set end


    if (params.sampleStart[i] < 0) or (params.sampleEnd[i] < 0):
      showError( 'Peak Separator', 'The region you have chosen to pick appears to be invalid. Have you selected the correct peak list?')
      return

  params.sampleSize     = []
  params.samplePpmSize  = []

  for i in range( params.Ndim ): # for higher dim cannot use spectra window

    params.sampleSize.append(    max( 1, abs( params.sampleEnd[i] - params.sampleStart[i] ) ) )
    params.samplePpmSize.append( max( abs( params.samplePpmEnd[i] ), abs( params.samplePpmEnd[i] - params.samplePpmStart[i] ) ) )

  block_file = dataSource.block_file

  # first and last only used in block_file.maxValue etc. Dont use elsewhere
  first = params.sampleStart[:] # (copy)
  last  = params.sampleEnd[:]

  # for i in range( params.Ndim ):
  #   if not params.isFreqDim[i]:
  #     if first[i] == last[i]:
  #       # print first[i], last[i]
  #       first[i] -= 1

  params.maxHeight        = HEIGHT_MULTIPLIER * max(abs(block_file.minValue(first, last)), abs(block_file.maxValue(first, last)))
                                                                                  # set maxHeight

