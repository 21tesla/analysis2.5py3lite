"""
======================COPYRIGHT/LICENSE START==========================

Minor.py: Data compatibility handling

Copyright (C) 2007 Rasmus Fogh (CCPN project)
 
=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license.
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

======================COPYRIGHT/LICENSE END============================

To obtain more information about this code:

- CCPN website (http://www.ccpn.ac.uk)

- contact Rasmus Fogh (ccpn@bioc.cam.ac.uk)

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following reference:

===========================REFERENCE START=============================
Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and 
automated software development. Bioinformatics 21, 1678-1684.
===========================REFERENCE END===============================
 
"""

# NB this file will only be used as part of Major upgrades

from memops.general.Implementation import ApiError

def correctData(topObj, delayDataDict, toNewObjDict, mapping=None):
  """ update topObj object tree using information in delayDataDict
  May be used either to postprocess a file load (minor upgrade)
  or as part of an in-memory data transfer (major upgrade)
  
  topObj is the MemopsRoot in the new tree
  toNewObjDict is _ID:newObj for minor 
    and oldObj/oldObjId:newObj for major upgrades
    
    NB this function does nothing without further additions
  """
  
  emptyDict = {}
  emptyList = []
  doGet = delayDataDict.get
  
  
  pName = topObj.packageName
  
  #
  if pName == 'ccpnmr.Analysis':
    # Fix molType and ccpCode (change in code values)
    fixAnalysis(topObj, delayDataDict, toNewObjDict, mapping)
  #
  elif pName == 'ccp.nmr.Nmr':
    # Fix Nmr
    fixNmr(topObj, delayDataDict)
  #
  elif pName == 'memops.Implementation':
    # Add extra repository 'generalData' 
    fixImplementation(topObj, delayDataDict)
  
  # 
  from memops.format.compatibility.upgrade.v_2_0_a2 import Minor as Minor20a2
  Minor20a2.correctData(topObj, delayDataDict, toNewObjDict, mapping)
    
def fixImplementation(topObj, delayDataDict):
  """ Add generalData repository, packageLocator for AnalysisProfile with
  repositories link, and set
  
  NB when this is called PAckageLocators and Repositories
  have already had their child links and crosslinks set.
  """
  from memops.universal import Io as uniIo
  from memops.api.Implementation import Url
  import os.path
  
  emptyDict = {}
  emptyList = []
  doGet = delayDataDict.get
  urlPath = uniIo.normalisePath(os.path.expanduser('~/.ccpn/data'))
  repositories = doGet(topObj).get('repositories')
  for refData in repositories:
    if refData.name == 'refData':
      break
  else:
    raise ApiError("refData repository not found")
  genData = topObj.newRepository(name='generalData', url=Url(path=urlPath))
  topObj.__dict__['repositories']['generalData'] = genData
  
  profileLocator = topObj.newPackageLocator(targetName='ccpnmr.AnalysisProfile',
                                             repositories=(genData, refData))
  topObj.__dict__['packageLocators']['ccpnmr.AnalysisProfile'] = profileLocator
    
def fixNmr(topObj, delayDataDict):
  """ Fix DimensionScaling
  """
  emptyDict = {}
  emptyList = []
  doGet = delayDataDict.get
  
  memopsRoot = topObj.parent
  topObjByGuid = memopsRoot.__dict__.get('topObjects')
  
  for xpr in doGet(topObj, emptyDict).get('experiments', emptyList):
    expDict = doGet(xpr, emptyDict)
    
    # fix NmrExpPrototype mapping for defunct types
    setNmrExpPrototypeLink(xpr, 'refExperiment', topObjByGuid, delayDataDict, 
                           remapPrototypeLink)
    for xd in expDict.get('expDims', emptyList):
      setNmrExpPrototypeLink(xd, 'refExpDim', topObjByGuid, delayDataDict,
                             remapPrototypeLink)
      for xdr in doGet(xd, emptyDict).get('expDimRefs', emptyList):
        setNmrExpPrototypeLink(xdr, 'refExpDimRef', topObjByGuid, delayDataDict,
                               remapPrototypeLink)
    
    # Fix DimensionScalings
    for ds in expDict.get('dataSources', emptyList):
      for dd in doGet(ds, emptyDict).get('dataDims', emptyList):
        if dd.__class__.__name__ != 'SampledDataDim':
          allExpDimRefs = set(
           doGet(dd.expDim, emptyDict).get('expDimRefs', emptyList)
          )
 
          # set up proper structure if not there already
          ddDict = doGet(dd)
          if not ddDict:
            delayDataDict[dd] = ddDict = {}
          dimensionScalings = ddDict.get('dimensionScalings')
          if not dimensionScalings:
            ddDict['dimensionScalings'] = dimensionScalings = []
 
          usedExpDimRefs = set(x.expDimRef for x in dimensionScalings)
          for xdr in allExpDimRefs.difference(usedExpDimRefs):
            # ExpDimRefs that had no DimensionScaling pointing to them
            # Used to mean 'default' now means 'not used'. Change
            scaling = dd.newDimensionScaling(expDimRef=xdr,
                                             scalingFactors=(1.0,))
            dimensionScalings.append(scaling)

def remapPrototypeLink(keyList):
  
  from memops.format.compatibility.upgrade.v_2_0_b2 import Minor as Minor20b2
  
  guid = keyList[0]
    
  # map to different keys
  if guid == 'cam_wb104_2008-01-15-16-06-40_00004':
    # 'H[N[CA[CB]]] defunct
    keyList[0] = 'cam_wb104_2008-01-15-16-06-40_00024'
    keyList[1] = keyList[1] - 2 # refExperiment serial
    
  elif guid == 'cam_wb104_2008-01-15-16-06-40_00033':
    # 'H[N[CO[C[C]]]] defunct
    keyList[0] = 'cam_wb104_2008-01-15-16-06-40_00035'
    refExpMap = {4:1, 6:9, 7:10, 8:11}
    keyList[1] = refExpMap.get(keyList[1],0) # refExperiment serial
  
  Minor20b2.remapPrototypeLink(keyList)


def setNmrExpPrototypeLink(obj, tag, topObjByGuid, delayDataDict,
                           linkmapper):
  """ redirect certain NmrExpPrototype links to other experiments
  """
  doGet = delayDataDict.get
  objDataDict = doGet(obj)
  inDataList = objDataDict.get(tag)
  if inDataList:
    keyList = inDataList[0]
    
    linkmapper(keyList)
    guid = keyList[0]
    
    # set link
    oo = topObjByGuid.get(guid)
    clazz = keyList[-1]['class']
    if (oo is None):
      # NB naughty - _packageName is a private attribute. 
      # But getPackageName is not static
      obj.root.refreshTopObjects(clazz._packageName)
      try:
        oo = topObjByGuid[guid]
      except:
        raise ApiError("""%s.%s: NmrExpPrototype with guid %s not found or loaded"""
                       % (clazz.__name__, tag, guid))
    obj.__dict__[tag] = clazz.getByKey(oo, keyList[1:-1])
    del objDataDict[tag]
    

def fixAnalysis(topObj, delayDataDict, toNewObjDict, mapping):
  """ Fix Analysis package with key changes
  """
  emptyDict = {}
  emptyList = []
  doGet = delayDataDict.get
  
  memopsRoot = topObj.parent
  
  nextSerial = 0
  newData = {}
  topObjDict = doGet(topObj, emptyDict)
  asList = topObjDict.get('analysisSpectra')
  if asList is None:
    asList = topObjDict['analysisSpectra'] = []
  
  for swin in topObjDict.get('spectrumWindows', emptyList):
    for winv in doGet(swin, emptyDict).get('spectrumWindowViews', emptyList):
      dd = doGet(winv)
      dataSource = topObj.nmrProject.findFirstExperiment(
       serial=dd['experimentSerial'][0]).findFirstDataSource(
       serial=dd['dataSourceSerial'][0]
      )
 
      analysisSpectrum = newData.get(dataSource)
      if analysisSpectrum is None:
        analysisSpectrum = topObj.newAnalysisSpectrum(dataSource=dataSource)
        newData[dataSource] = analysisSpectrum
        asList.append(analysisSpectrum)
        delayDataDict[analysisSpectrum] = {
         'analysisPeakLists':[],
         'analysisDataDims':[]
        }
      winv.analysisSpectrum = analysisSpectrum
      aplList = doGet(analysisSpectrum)['analysisPeakLists']
      addList = doGet(analysisSpectrum)['analysisDataDims']
      
      for wpl in doGet(winv, emptyDict).get('windowPeakLists', emptyList):
        dd = doGet(wpl)
        peakList = dataSource.findFirstPeakList(serial=dd['peakListSerial'][0])
 
        analysisPeakList = newData.get(peakList)
        if analysisPeakList is None:
          analysisPeakList = analysisSpectrum.newAnalysisPeakList(peakList=peakList)
          newData[peakList] = analysisPeakList
          aplList.append(analysisPeakList)
        wpl.analysisPeakList = analysisPeakList
      
      for axisMapping in doGet(winv, emptyDict).get('axisMappings', emptyList):
        dd = doGet(axisMapping)
        dataDim = dataSource.findFirstDataDim(dim=dd['dim'][0])
 
        analysisDataDim = newData.get(dataDim)
        if analysisDataDim is None:
          analysisDataDim = analysisSpectrum.newAnalysisDataDim(dataDim=dataDim)
          newData[dataDim] = analysisDataDim
          addList.append(analysisDataDim)
        axisMapping.analysisDataDim = analysisDataDim

  
