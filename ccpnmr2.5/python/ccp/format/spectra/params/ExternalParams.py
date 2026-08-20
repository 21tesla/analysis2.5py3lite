
"""
======================COPYRIGHT/LICENSE START==========================

ExternalParams.py

Copyright (C) 2008 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2008 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""

import re, os

from memops.universal.Io import normalisePath
from memops.universal.Util import isBigEndian

from memops.general.Implementation import ApiError

from ccp.general.Io import getDataStoringFromFilepath

from memops.api import Implementation

isotope_re = re.compile('\d+(\D+)')

# DImension type constants
freqDimType = 'frequency'
fidDimType = 'fid'
sampledDimType = 'sampled'

# Top level parameter tags with their default values
topParamDefs = {'ndim':0, 'head':0, 'blockHead':0, 'dataFile':'', 
                'integer':False, 'big_endian':None, 'swap':False,
                'pulProgName':None, 'pulProgType':None, 'nAxes':0, 
                'acquisitionDim':None,'scale':1.0, 'nbytes':4,
}
    
# Per-dimension parameter tags with their default values
dimParamDefs = {'block':None, 'npts':None, 'sw':1000.0, 'sf':1000.0,
                'refppm':0.0, 'refpt':1.0, 'nuc':'1H', 'isComplex':False,
                'dimType':freqDimType, 'pointValues':None,
                'baseFreq':None, 'numPointsOrig':None, 'pointOffset':0,
}

# per-dimension parameters handled outside normal mechanisms
extraDimParams = ['displayNames', 'scalingFactors']

class ExternalParams:

  format = ''  # should be overridden in subclasses

  
  # how many directories to keep in file path for default DatUrl split
  keepDirectories = 1

  def __init__(self, externalParams=None, extraData=None):
    
    if extraData:
      self.extraData = extraData.copy()
    else:
      self.extraData = {}

    if externalParams:
      self.cloneParams(externalParams)
    else:
      self.makeNewParams()

  def __eq__(self, other):

    if not isinstance(other, ExternalParams):
      return False

    # do not want below because using this for Pipe header comparisons
    #if self.dataFile != other.dataFile:
    #  return False

    if self.ndim != other.ndim:
      return False

    if self.head != other.head:
      return False

    if self.blockHead != other.blockHead:
      return False

    if self.swap != other.swap:
      return False

    if self.integer != other.integer:
      return False

    if self.big_endian != other.big_endian:
      return False

    if self.npts != other.npts:
      return False

    if self.block != other.block:
      return False

    if self.sw != other.sw:
      return False

    if self.sf != other.sf:
      return False

    if self.refpt != other.refpt:
      return False

    if self.refppm != other.refppm:
      return False

    if self.nuc != other.nuc:
      return False

    return True

  def cloneParams(self, externalParams):

    for param in topParamDefs:
      setattr(self, param, getattr(externalParams, param))
    self.topParams = externalParams.topParams.copy()
    
    for param in dimParamDefs:
      setattr(self, param, getattr(externalParams, param)[:])
    self.dimParams = [x.copy() for x in externalParams.dimParams]
    
    
  def makeNewParams(self):

    try:
      self.init()
      self.parseFile()
      self.checkValid()
    except ApiError, e:
      raise ApiError('%s: make sure you selected correct file for %s format' % (e.error_msg, self.format))
    """
    except Exception, e:
      raise ApiError('%s: make sure you selected correct file for %s format' % (e, self.format))
    """

  def __ne__(self, other):

    return not self.__eq__(other)

  def init(self):
  
    for tag, val in topParamDefs.items():
      if not hasattr(self, tag):
        setattr(self, tag, val)
    
    # introduced for easy addition of extra parameters
    self.topParams = {}

  def initDims(self):
    
    # Modified to use nAxes if set (for projection experiments)
    ndim = self.nAxes or self.ndim
    
    # initialise parameters
    for tag, value in dimParamDefs.items():
      setattr(self, tag, ndim*[value])
    
    # special case = dictionaries must be set one by one
    # introduced for easy addition of extra parameters
    self.dimParams = [{} for x in range(ndim)]

  # useful code for setting attributes after generic format parsing
  def setAttrs(self, aPars):
    
    # modified to use nAxes if present (for projection data)
    ndim = aPars['ndim']
    self.ndim = ndim
    self.nAxes = aPars.get('nAxes', ndim)
    self.initDims()

    # TBD: aPars code should do file --> dataFile
    self.dataFile = aPars.get('file', '')
    for key in topParamDefs.keys():
      if aPars.has_key(key):
        setattr(self, key, aPars[key])
    for key in dimParamDefs.keys():
      if aPars.has_key(key):
        ll = aPars[key]
        getattr(self, key)[:len(ll)] = ll

  # introduced for NmrPipe data, where pseudo data seems to be indicated
  # by nuc type ID, which gets changed to None
  def pseudoDataDim(self):

    dims = []
    for dim in range(self.ndim):
      
      if self.dimType[dim] != freqDimType or self.nuc[dim] is None:
        dims.append(dim)
    
    if len(dims) == 1:
      return dims[0]
    else:
      return None

  # introduced for NmrPipe data, where sampled data is determined quasi-manually
  def setSampledDim(self, dim, values):

    if dim < 0 or dim >= self.ndim:
      raise ApiError('dim = %d, must be between 0 and %d' % (dim, self.ndim-1))

    if len(values) != self.npts[dim]:
      raise ApiError('number of values = %d, must be %d' % (len(values), self.npts[dim]))

    self.dimType[dim] = sampledDimType
    self.pointValues[dim] = values
    self.nuc[dim] = None

  # introduced for NmrView data, where extraneous dimensions might be allowed
  # where npts = 1
  def removeTrivialDims(self):

    ndim = self.ndim
    for i in range(ndim-1, -1, -1):
      if self.npts[i] <= 1:
        del self.npts[i]
        del self.block[i]
        del self.refppm[i]
        del self.refpt[i]
        del self.sw[i]
        del self.sf[i]
        del self.nuc[i]
        del self.dimType[i]
        del self.pointValues[i]
        self.ndim = self.ndim - 1

  # introduced for NmrPipe data, where extraneous dimensions might be allowed
  # (although alternatively these are pseudoND experiments)
  def removeNullDims(self):

    ndim = self.ndim
    for i in range(ndim-1, -1, -1):
      if self.nuc[i] is None:
        del self.npts[i]
        del self.block[i]
        del self.refppm[i]
        del self.refpt[i]
        del self.sw[i]
        del self.sf[i]
        del self.nuc[i]
        del self.dimType[i]
        del self.pointValues[i]
        self.ndim = self.ndim - 1

  # introduced as less drastic alternative to removeNullDims
  # only remove dimensions where null nucleus and npts <= 1
  # otherwise set nucleus to standard one
  def fixNullDims(self, ignoreDim=None):

    ndim = self.ndim
    for i in range(ndim-1, -1, -1):
      if i is ignoreDim:
        continue
      if self.nuc[i] is None or self.npts[i] <= 1:
        if self.npts[i] <= 1:
          del self.npts[i]
          del self.block[i]
          del self.refppm[i]
          del self.refpt[i]
          del self.sw[i]
          del self.sf[i]
          del self.nuc[i]
          del self.dimType[i]
          del self.pointValues[i]
          self.ndim = self.ndim - 1
        else:
          self.nuc[i] = self.standardNucleusName()
          self.dimType[i] = freqDimType

  def parseFile(self):

    raise ApiError('must implement parseFile in subclass of ExternalParams')

  def checkValid(self):

    # TBD
    pass

  def standardNucleusName(self, name = ''):

    name = name.upper()
    if (not name):
      name = '1H'
    elif ((name[0] == 'H') or (name == '1H')):
      name = '1H'
    elif ((name[0] == 'C') or (name == '13C')):
      name = '13C'
    elif ((name[0] == 'N') or (name == '15N')):
      name = '15N'
    elif ((name[0] == 'O') or (name == '17O')):
      name = '17O'
    elif ((name[0] == 'F') or (name == '19F')):
      name = '19F'
    elif ((name[0] == 'P') or (name == '31P')):
      name = '31P'
    else:
      name = '1H'

    return name

  def createDataSource(self, experiment, name, dataUrl=None):
    
    # set displayNames if available
    ll = [dd.get('name') for dd in self.dimParams]
    if [x for x in ll if x]:
      # we have some displayNames
      self.displayNames = ll
    
    extraData = self.extraData
    if extraData:
      # set/override data from input extra parameters
      for tag in dimParamDefs.keys() + extraDimParams:
        val = extraData.get(tag)
        if val:
          setattr(self, tag, val)
      
      # set referencing from carppm (carrier pm value)
      carppm = extraData.get('carppm')
      if carppm:
        maxIndex = len(self.numPointsOrig) - 1
        for ii,refppm in enumerate(carppm):
          if self.refpt[ii] == 1.0:
            # reset index to allow for the size arrays being shorter
            # basically use last point by overshoot (project dimensions)
            # skip if refpt !- 1.0, to preserve some user-set values
            jj = min(ii, maxIndex)
            refpt = self.numPointsOrig[jj]/2 +1 - self.pointOffset[jj]
            self.refpt[ii] = refpt
            self.refppm[ii] = refppm
    
    if self.ndim < 1:
      raise ApiError('ndim = 0 (did you choose correct format?)')

    # check that information compatible with experiment
    if experiment.numDim != self.ndim:
      raise ApiError('for now experiment numDim (%d) and dataSource numDim (%d) need to be the same' % (experiment.numDim, self.ndim))

    for n in range(self.ndim):
      # TBD: is this correct? (e.g. do fidDimType have nuc?)
      if self.dimType[n] != freqDimType:
        if self.dimType[n] == fidDimType:
          raise ApiError('dimension %s is FID - not yet supported.' % (n+1))
        else:
          self.nuc[n] = None
    
    # match ExpDimRefs if any, to dims.
    expDims = experiment.sortedExpDims()
    for n in range(self.ndim):
      expDim = expDims[n]

      if self.dimType[n] == freqDimType:
        if expDim.expDimRefs:
          if len(expDim.expDimRefs) == 1:
            expDimRef = expDim.findFirstExpDimRef()
            if (self.nuc[n]):
              if expDimRef.isotopeCodes != (self.nuc[n],):
                print ('Warning: in dim %d experiment has isotopeCodes %s, spectrum has %s, must be the same'
                       % (n+1, expDimRef.isotopeCodes, (self.nuc[n],)))
              if (abs(expDimRef.sf-self.sf[n]) > 1.0e-2):
                print ('Warning: inconsistent sf in dim %d: %3.2f vs %3.2f' 
                       % (n+1, expDimRef.sf, self.sf[n]))
            elif expDimRef.isotopeCodes:
              print ('Warning: in dim %d experiment has isotopeCodes %s, spectrum has no nuc set' 
                     % (n+1, expDimRef.isotopeCodes))
          else:
            # NBNB TBD new code - may not work for all cases
            if (self.nuc[n]):
              expDimRef = None
              expDimRefs = [x for x in expDim.sortedExpDimRefs() 
                            if x.isotopeCodes==(self.nuc[n],)]
              for xdr in expDimRefs:
                if abs(xdr.sf-self.sf[n]) <= 1.0e-2:
                  expDimRef = xdr
                  break
                elif not xdr.sf:
                  expDimRef = xdr
                  
              if expDimRef is None:
                if expDimRefs:
                  expDimRef = expDimRefs[0]
                  print ('Warning: in dim %d no ExpDimRef fits sf %s' 
                         % (n+1, self.nuc[n]), expDimRef.sf, self.sf[n])
                  
                else:
                  expDimRef = expDim.findFirstExpDimRef()
                  print ('Warning: in dim %d no ExpDimRef fits nucleus %s' 
                         % (n+1, self.nuc[n]), expDimRef.sf, self.sf[n])
            else:
              expDimRef = expDim.findFirstExpDimRef()
              print (
               'Warning: in dim %d spectrum has no nuc set' % (n+1)
              )

      else: # non-freq dim
        if expDim.expDimRefs:
          raise ApiError('non-freq dim has experiment dim ref')
    
    if self.dataFile:
      # create DataStore
      
      memopsRoot = experiment.root
      if dataUrl:
        preferDataUrls = [dataUrl]
      else:
        oldDataStores = (x.dataStore for x in experiment.dataSources)
        preferDataUrls = set(x.dataUrl for x in oldDataStores if x)
      
      newDataUrl, filePath = getDataStoringFromFilepath(memopsRoot, 
                            self.dataFile, preferDataUrls=preferDataUrls,
                            keepDirectories=self.keepDirectories)
      
      if dataUrl is not None and newDataUrl is not dataUrl:
        raise ApiError("DataUrl with path %s not valid for file %s" 
                       % (dataUrl.url.path, self.dataFile) )

      dataLocationStore = newDataUrl.dataLocationStore
    
      if self.big_endian is not None:
        bigEndian = self.big_endian
      elif self.swap:
        bigEndian = not isBigEndian()
      else:
        bigEndian = isBigEndian()

      if self.integer:
        numberType = 'int'
      else:
        numberType = 'float'
      
      isComplex = [self.dimParams[ii].get('isComplex') 
                   for ii in range(self.ndim)]
      if None in isComplex:
        isComplex = self.ndim * [False]
      
      # numPoints - necessary in case of projection spectra
      numPoints = self.npts
      if None in numPoints:
        numPoints = numPoints[:numPoints.index(None)]
      
      
      params = {
       'dataUrl':newDataUrl, 'path':filePath, 'numPoints':numPoints,
       'isBigEndian':bigEndian, 'numberType':numberType, 'isComplex':isComplex,
       'headerSize': self.head, 'fileType':self.format, 'nByte':self.nbytes,
      }
 
      if self.format == 'Factorised':
        for tag in ('numShapes', 'isReconstructable', 
                    'isResolved', 'numRecords'):
          if tag in self.topParams:
            params[tag] = self.topParams[tag]
        dataStore = dataLocationStore.newShapeMatrix( **params)
        
        # NBNB TBD assumes there are no 2D shapes. To be changed!
        if len(self.dimParams) == self.topParams['numShapes']:
          # all shapes are 1D
          for dim in range(self.topParams['numShapes']):
            dataStore.newComponentShape(dims = (dim,))
        else:
          raise ApiError("Decompositions must have one shape per dim (for now)")
        
        # add Component elements
        for dd in self.topParams['componentList']:
          dataStore.newComponent(**dd)
        
        
      else:
        # default: blocked binary matrix

        if self.blockHead:
          params['blockHeaderSize'] = self.blockHead
        
      
        # blockSizes - necessary in case of projection spectra
        blockSizes = self.block
        if None in blockSizes:
          blockSizes = blockSizes[:blockSizes.index(None)]

        # Bruker in Topspin 3 can have blockSize = 0 it seems
        for nn, blockSize in enumerate(blockSizes):
          if blockSize == 0:
            blockSizes[nn] = numPoints[nn]
        
        dataStore = dataLocationStore.newBlockedBinaryMatrix(
                        blockSizes=blockSizes, **params)

    else:
      dataStore = None
    
    #set ndims and set up for extra (projection) dims
    ndims = self.ndim
    lastEdrs = expDims[-1].expDimRefs
    if (self.nAxes > ndims and self.dimType[-1] == freqDimType and
        (not lastEdrs or self.nAxes <= ndims + len(lastEdrs))):
      # We have extra (projection) ExpDimRefs, 
      # They fit onto existing, or the last dimension still needs ExpDimRefs
      # Set up to add extra ExpDimRefs to last ExpDim
      expDims.extend((self.nAxes-ndims)*[expDims[-1]])
      ndims = self.nAxes
    
    #set ExpDimRefs and isAcquisition
    hasNoAcqDim = (expDims[0].experiment.findFirstExpDim(isAcquisition=True) 
                   is None)
    
    for n in range(ndims):
      
      if self.dimType[n] == freqDimType:
        expDim = expDims[n]
        lenExpDimRefs = len(expDim.expDimRefs)
 
        if (not lenExpDimRefs or (n >= self.ndim + lenExpDimRefs -1)):
           # The first term ensures an ExpDimRef for every freqDim
           # The second adds extra expDimRefs to the last Dim to use up all
           # read Dims (for projection experiments
           
          if hasNoAcqDim and n == self.acquisitionDim:
            expDim.isAcquisition = True
            hasNoAcqDim = False
 
          if self.nuc[n]:
            xdr = expDim.newExpDimRef(sf=self.sf[n], unit='ppm',
                                      baseFrequency=self.baseFreq[n],
                                      isotopeCodes=[self.nuc[n]])
          else:
            xdr = expDim.newExpDimRef(sf=self.sf[n], unit='None',
                                      baseFrequency=self.baseFreq[n],)
        
    # set displayName
    if hasattr(self, 'displayNames'):
      displayNames = self.displayNames
      if displayNames:
        ii = 0
        for expDim in experiment.sortedExpDims():
          for expDimRef in expDim.sortedExpDimRefs():
            expDimRef.displayName = displayNames[ii]
            ii += 1
    
    dataSource = experiment.newDataSource(name=name, numDim=self.ndim,
                                          dataStore=dataStore, 
                                          dataType='processed',
                                          scale=self.scale)
    scalingFactors = None
    if hasattr(self, 'scalingFactors'):
      ll = self.scalingFactors
      if ll:
        scalingFactors = ll
    
    # Set DataSource.DataDims
    for n in range(self.ndim):
      sw = self.sw[n]
      npts = self.npts[n]
      numPointsOrig = self.numPointsOrig[n] or npts
      expDim = expDims[n]
      if self.dimType[n] == freqDimType:
        freqDataDim = dataSource.newFreqDataDim(dim=n+1, numPoints=npts, 
                                                isComplex=False, 
                                                numPointsOrig=numPointsOrig,
                                                pointOffset=self.pointOffset[n], 
                                                valuePerPoint=sw/float(npts), 
                                                expDim=expDim)
        expDimRef=expDim.findFirstExpDimRef()
        ddr = freqDataDim.newDataDimRef(refPoint=self.refpt[n],
                                        refValue=self.refppm[n],
                                        expDimRef=expDimRef)
        
        # reset to reference on O1
        centrePoint = freqDataDim.numPointsOrig/2 - freqDataDim.pointOffset + 1
        ddr.refValue = ddr.pointToValue(centrePoint)
        ddr.refPoint = centrePoint
        
        if scalingFactors is not None:
          # NB DimensionScaling is set automatically in some cases
          previousDimScaling = freqDataDim.findFirstDimensionScaling(
                                expDimRef=expDimRef)
          if previousDimScaling is None:
            freqDataDim.newDimensionScaling(expDimRef=expDimRef, 
                                    scalingFactors=(scalingFactors[n],))
          else:
            previousDimScaling.scalingFactors = (scalingFactors[n],)
        
      elif self.dimType[n] == fidDimType:
        raise ApiError('fidDimType not supported yet')
      elif self.dimType[n] == sampledDimType:
        # TBD: set conditionVaried somehow
        dataSource.newSampledDataDim(dim=n+1, numPoints=npts, isComplex=False, 
                              pointValues=self.pointValues[n], expDim=expDim)
      else:
        raise ApiError('%s dimType not supported yet' % self.dimType[n])
    
    # Add extra (projection) DataDimRefs to last dataDim
    # NB if ndims > self.ndim the last dimension is a freqDataDim
    # freqDataDim, expDim, and npts come correctly out of the preceding loop.
    expDimRefs = expDim.sortedExpDimRefs()
    ii = 0
    for n in range(self.ndim, ndims):
      ii += 1   # NB we do want to start with expDimRefs[1]
      expDimRef=expDimRefs[ii]
      ddr = freqDataDim.newDataDimRef(refPoint=self.refpt[n],
                                      refValue=self.refppm[n],
                                      expDimRef=expDimRefs[ii],
                                      localValuePerPoint=self.sw[n]/npts)
      # reset to reference on O1
      centrePoint = freqDataDim.numPointsOrig/2 - freqDataDim.pointOffset + 1
      ddr.refValue = ddr.pointToValue(centrePoint)
      ddr.refPoint = centrePoint
      
      if scalingFactors is not None:
        freqDataDim.newDimensionScaling(expDimRef=expDimRef, 
                                scalingFactors=(scalingFactors[n],))
    #  
    return dataSource
