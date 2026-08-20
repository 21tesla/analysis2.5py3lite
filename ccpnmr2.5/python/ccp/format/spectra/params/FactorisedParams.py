
"""
======================COPYRIGHT/LICENSE START==========================

AzaraParams.py

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
import os
from memops.universal.ElementTree import ElementTree as ET
from memops.general.Implementation import ApiError
from ccp.format.spectra.params.ExternalParams import ExternalParams
from ccp.format.spectra.params.ExternalParams import freqDimType, fidDimType
from ccp.util.Software import getMethod
from memops.universal import BlockData

class FactorisedParams(ExternalParams):

  format = 'Factorised'

  def __init__(self, parFile, **kw):

    self.parFile = parFile
    self.dim = -1

    ExternalParams.__init__(self, **kw)

  # ExternalParams requires this to be defined
  def parseFile(self):
    
    # read file and get Decomposition node :
    
    try:
      elementTree = ET.parse(self.parFile)
    except IOError, e:
      raise ApiError(str(e))
      
    topNode = elementTree.getroot()
    if topNode.tag == 'Decomposition':
      # Decomposition is document element
      decomposition = topNode
    else:
      # look for Decompositions directly within document element
      decompositions = list(topNode.findall('Decomposition'))
      if decompositions:
        if len(decompositions) > 1:
          print ('WARNING, found %s Decompositions in file %s\n Using first one'
                 % (len(decompositions),  self.parFile))
        decomposition = decompositions[0]
      
      else:
        print 'WARNING, no Decompositions found in file',  self.parFile
        return
      
    # get top-level parameters from Decomposition node :
    
    self.dataFile = self.parFile
    
    # irrelevant anyway - set to model default to avoid attempts to determine
    self.big_endian = True
    
    self.topParams['nodeTree'] = elementTree
    self.topParams['numShapes'] = int(decomposition.get('nshapes'))
    self.topParams['isReconstructable'] = (decomposition.get('reconstructable') 
                                           in ('True', 'true'))
    self.topParams['isResolved'] = (decomposition.get('resolved') 
                                           in ('True', 'true'))
    pulProgName = decomposition.get('refexperiment')
    if pulProgName is not None:
      self.pulProgName = pulProgName
      self.pulProgType = 'ccpn'
    
    # the following are not meaningful and/or use the defaults:
    # 'head',  'integer', 'numRecords'
    
    
    # get axes in shape number order
    ll = []
    for elem in decomposition.findall('Axis'):
      ll.append((int(elem.get('a')), elem))
    axes = [x[1] for x in sorted(ll)]
    
    # get ndims and initialise
    
    self.ndim = len(axes)
    self.initDims()
 
    # per-dimension parameters (from Axis) :
    
    # 'npts', 'nuc', 'sw', 'sf', 'refppm', 'refpt' 'dimType', 
    # The following is not relevant: 'block'
    # extra parameters 'name', 'isComplex', 'shapeSerial'
    
    for dim,axis in enumerate(axes):
      
      self.npts[dim] = npts = int(axis.get('size'))
      self.sf[dim] = float(axis.get('sfo'))
      swppm = float(axis.get('swppm'))
      self.sw[dim] = swppm * self.sf[dim]
      self.refppm[dim] = refppm = float(axis.get('startppm'))
      self.refpt[dim] = refpt = 1.0
      self.nuc[dim] = self.standardNucleusName(axis.get('nucleus'))
      ss = axis.get('origsize')
      if ss:
        self.numPointsOrig[dim] = numPointsOrig = int(axis.get('origsize'))
      else:
        numPointsOrig = npts
      if axis.get('domain').lower() == 'time':
        self.dimType[dim] = fidDimType
      else:
        self.dimType[dim] = freqDimType
      self.isComplex[dim] = (axis.get('type').lower() == 'complex')
        
      dimParams = self.dimParams[dim]
      dimParams['name'] = axis.get('name')
      dimParams['shapeSerial'] = int(axis.get('a'))
      
      # carrier ppm. Reset referencing, and handle point offset
      carppm = axis.get('carppm')
      if carppm:
        carppm = float(carppm)
        refpt = self.refpt[dim] = refpt + (refppm - carppm) * npts / swppm
        refppm = self.refppm[dim] = carppm
        # NB add 1.5 before int call to get nearest int 
        # (pointOffset is alwayspositive)
        pointOffset = int(numPointsOrig / 2 + 1.5 - refpt)
        
        if numPointsOrig == npts:
          if pointOffset > npts or pointOffset < 0:
            # LikelynumPointsOrig was not set correctly
            pointOffset = 0
        self.pointOffset[dim] = pointOffset
        
    # get Component parameters :
    ll = []
    for elem in decomposition.findall('Component'):
      
      dd = {}
      
      index = elem.get('c')
      if index is None:
        raise ApiError(
         "Incorrect component, lacks 'c' attribute: <Component %s >"
         % elem.attrib
        )
      
      ampl = elem.get('ampl')
      if ampl is None:
        raise ApiError(
         "Incorrect component, lacks 'ampl' attribute: <Component %s >"
         % elem.attrib
        )
      else:
        dd['amplitude'] = float(ampl)
      
      regionId = elem.get('regionid')
      if regionId:
        dd['regionId'] = int(regionId)
      
      dd['status'] = elem.get('status')
      
      dd['annotation'] = elem.get('annotation')
      
      ll.append((index, dd))
    #
    self.topParams['componentList'] = [x[1] for x in sorted(ll)]
        
      
    
  def parsePeaks(self, dataSource, software = None):
    """Parse peaks from file and add them to dataSource new peaklist
    """
    
    decomposition = self.topParams['nodeTree'].getroot()
    nshapes = int(decomposition.get('nshapes'))
    ss = decomposition.get('ncomponents')
    if ss:
      ncomp = int(ss)
    else:
      ncomp = None
    name = decomposition.get('name')
    peakList = dataSource.newPeakList(name=name)
    
    if nshapes != dataSource.numDim:
      raise Exception(
       "Reading peaks for %s:Not implemented yet: nshapes %s != numDims %s"
       % (dataSource, nshapes, dataSource.numDim)
      )
    
    else:
      # one shape per dimension. Add peaks 
      if not software:
        software = getSoftware(dataSource.root, name='unknown', version='unknown')
      heightMethod = getMethod(software, task='find peak height',
                               procedure='Factorised dataSource import')
 
      compFound = 0
      for elem in decomposition:
        if elem.tag == 'Component':
          compFound += 1
          ampl = float(elem.get('ampl'))
          compNo = int(elem.get('c'))
          dd = {}
          lenArray = []
          for ii,ee2 in enumerate(elem):
            if ee2.tag == 'Shape':
              a = int(ee2.get('a'))
              data = None
              for ee3 in ee2:
                if ee3.tag == 'Peaks':
                  data = unpackListElement(ee3)
                  break
              if data:
                dd[a] = data
                lenArray.append(len(data))
              else:
                break
          else:
            # There were peaks picked for all shapes
            # now make peak
            peaks1D = [x[1] for x in sorted(dd.items())]
            nPeaks, cumul = BlockData.cumulativeArray(lenArray)
            for index in range(nPeaks):
              # 26 May 2011: WARNING: arrayOfIndex implementation was broken
              # for Analysis usage so it was fixed but that might break below
              arr = BlockData.arrayOfIndex(index, cumul)
              height = ampl
              peak = peakList.newPeak(componentNumbers=(compNo,))
              peakDims = peak.sortedPeakDims()
              for ii in range(nshapes):
                dd = peaks1D[ii][arr[ii]]
                height *= dd['intensity']
                peakDim = peakDims[ii]
                pos = dd.get('pos')
                if pos is None:
                  posppm = dd['posppm']
                  peakDim.value = posppm
                else:
                  peakDim.position = pos
                pos = dd.get('pos')
                if pos is None:
                  posppm = dd['posppm']
                  peakDim.value = posppm
                else:
                  peakDim.position = pos
 
                # set up dataDimRef NBNB TBD could be neater
                ll = [(x.expDimRef.serial, x) 
                      for x in peakDim.dataDim.dataDimRefs]
                if ll:
                  ll.sort()
                  peakDim.dataDimRef = ll[0][1]
                else:
                  peakDim.dataDimRef = None
 
                # read width
                # NBNB TBD
 
              # NBNB we are not using PeakBasic.pickPeak or any of the setup
              # routines called by it, as we can assume that all peaks found
              # are in the principal region.
              peak.newPeakIntensity(intensityType='height', value=height,
                                    method=heightMethod)
      if compFound != ncomp:
        print ("WARNING, % components expected, % found" % (ncomp, compFound))
        

def unpackListElement(elem):
  """Unpacks list type element (with 'list' column name attribute)
  into a list of dictionari
  """
  result = []
  
  keywds = [x.strip() for x in elem.get('list').split(',')]
  ncols = len(keywds)
  if elem.text:
    data = elem.text.split()
  else:
    data = []
  if len(set(keywds)) == len(keywds):
    # no duplicate keywords
    nlines = len(data) // ncols
    ind = 0
    for ii in range(nlines):
      dd = {}
      result.append(dd)
      for jj in range(ncols):
        xx = data[ind]
        try:
          val = int(xx)
        except:
          try:
            val = float(xx)
          except:
            val = xx
        dd[keywds[jj]] = val
        ind += 1
    
  else:
    raise Exception("duplicate keywords parsing not implemented yet")
  
  #
  return result
