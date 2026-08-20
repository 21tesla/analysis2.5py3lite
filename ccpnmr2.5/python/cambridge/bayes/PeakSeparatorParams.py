#!/usr/bin/env python
# encoding: utf-8
"""
PeakSeparatorParams.py

Created by Daniel O'Donovan on 2010-11-10.
Copyright (c) 2010 University of Cambridge. All rights reserved.
"""

class PeakSeparatorParams(object):
  """ Just a class to hold the C Library Peak Separator parameters 
      Argument names / numbers reference - needed by C Library
                    # params.dataFile       0
                    # params.Ndim           1
                    # params.isBigEndian    2
                    # params.nPoints        3
                    # params.blockSize      4
                    # params.sampleStart    5
                    # params.sampleEnd      6
                    # params.maxHeight      7
                    # params.minHeight      8
                    # params.maxSigma       9
                    # params.minSigma       10
                    # params.maxQ           11
                    # params.minQ           12
                    # params.dimWrapped     13
                    # params.shape          14
                    # params.pos_peaks      15
                    # params.minAtoms       16
                    # params.maxAtoms       17
                    # params.rate           18

      Others to be set for convenience
  """

  def __init__(self):

    # These parameters chosen in gui
    self.peakList       = None
    self.peakShape      = 3   # 3 is Gaussian
    self.positivePeaks  = 1   # positive peaks ONLY           1
                              # positive AMD negative peaks   0
    self.minAtoms     = 1
    self.maxAtoms     = 1

    self.minSigma     = None
    self.maxSigma     = None

    self.rate         = 0.1

    # These parameters set from chosen peak list
    self.dataFile     = None
    self.Ndim         = None
    self.isBigEndian  = None
    self.nPoints      = None
    self.blockSize    = None
    self.dimWrapped   = None

    self.isFreqDim    = None
    self.dataDimRefs  = None

    # These parameters set from Region select or automated picking
    self.sampleStart  = None
    self.sampleEnd    = None
    self.sampleSize   = None

    self.maxHeight    = None
    self.minHeight    = None

    self.keys = [ key for key in self.__dict__.iterkeys() ]

    self.ClibKeys     = [ 'peakShape', 'positivePeaks', 'minAtoms', 'maxAtoms', \
                          'minSigma', 'maxSigma', 'rate', 'dataFile', 'Ndim', \
                          'isBigEndian', 'nPoints', 'blockSize', 'dimWrapped', \
                          'sampleStart', 'sampleEnd', 'sampleSize', 'maxHeight', 'minHeight' ]

  def getShape(self):
    """ return shape name: shape is only 3, 4, (5) """
    shapes = { 3:'Gauss 2d', 4:'Loren 2d', 5:'G/L   2d'}
    if self.peakShape in shapes.keys():
      return shapes[ self.peakShape ]
    else:
      return 'default'

