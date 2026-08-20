"""
======================COPYRIGHT/LICENSE START==========================

OpenProjectionSpectra.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

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

""" Read Projection Spectra definition file
All files are passed into a single experiment 
"""

from ccpnmr.analysis.popups import OpenSpectrum 

from memops.universal import Io as uniIo
from ccp.general import Util as genUtil

from ccp.general.Io import getDataStoringFromFilepath
from memops.general.Util import copySubTree


supportedFormats = ('Bruker',)

def loadProjDefinitionFile(argServer, inputFileName=None, exptName=None):
  
  # parse input file
  if inputFileName is None:
    inputFileName = argServer.getFile()
  
  if inputFileName:
    fileData = genUtil.parseBrukerSummary(inputFileName)
  else:
    return
  
  # get format specifier
  format = fileData.get('FORMAT')
  if format not in supportedFormats:
    print ("ERROR Unknown format; %s. Supported formats: %s"
           % (format, supportedFormats))
  
  # get OpenSpectrum popup and set up
  top = argServer.parent
  open_spectra = top.popups.get('open_spectrum')
  if open_spectra is None:
    top.openPopup('open_spectrum', OpenSpectrum.OpenSpectrumPopup)
    open_spectra = top.popups.get('open_spectrum')
  
  open_spectra.verifySelect.setSelected(False)
  open_spectra.sharedExpSelect.setSelected(True)
 
  # set up data structures
  extraData = {}
  specNames = []
  fileNames = []
  scalingFactors = []
  useScalingFactors = False
  
  if format == 'Bruker':
    
    keepDirectories = 4

    dimParNameMap = {
     'NUCLEI':'displayNames',
     #'SW_ppm':'swppm',
     #'SW_hz':'sw',
     #'O1_ppm':'carppm',
    }
    
    # set up special data
    numDims = None
    for tagin,tagout in dimParNameMap.items():
      val = fileData.get(tagin)
      if val is not None:
        if numDims is None:
          numDims = len(val)
        elif numDims != len(val):
          print ("ERROR for %s: numDims inconsistent, should be %d" 
                 % (tagin, numDims))
          return
        extraData[tagout] = val
    
    # get file and base directory
    ll = []
    for tt in fileData.items():
      try:
        ind = int(tt[0])
        # this is a string form of an integer. Use it
        ll.append((ind, tt))
      except ValueError:
        # not a string form of an integer - skip
        continue
    
    for xx in sorted(ll):
      tt = xx[1]
      specName, dd = tt
      fileName = dd.get('PATH')
      if fileName:
        specNames.append(str(specName))
        fileNames.append(fileName)
        scalingFactor = dd.get('DEFINITION')
        if scalingFactor:
          if numDims is None:
            # always the case first time for now, but the code might change
            numDims = len(scalingFactor)
          if len(scalingFactor) == numDims:
            useScalingFactors = True
          else:
            scalingFactor = None
        else:
          scalingFactor = None
        scalingFactors.append(scalingFactor)
    
    #indx = list(OpenSpectrum.file_formats).index('Bruker')
    #open_spectra.formatPulldown.setSelectedIndex(indx)
    #open_spectra.chooseFormat(indx,'Bruker')
    open_spectra.formatPulldown.set('Bruker')
    open_spectra.chooseFormat('Bruker')
  else:
    return
 
  # find and set set selection directory
  #baseDir = uniIo.commonSuperDirectory(*fileNames)
  
  startDir=os.path.dirname(inputFileName)
  baseDir, paths = uniIo.suggestFileLocations(fileNames, startDir=startDir)
  
  if baseDir is None:
    print 'WARNING, described files not found. Aborting'
  
  fileNames = [os.path.join(baseDir, path) for path in paths]
 
  defaultExpName = os.path.basename(baseDir)
  open_spectra.fileSelect.changeDir(baseDir)
 
  # set 'Spectra to open' matrix:
  objectList  = []
  textMatrix  = []
 
   # get experiment name
  if exptName is None:
    exptName = argServer.askString(message='New Experiment Name',
                                   default=defaultExpName)
    #exptName = defaultExpName
 
  shiftListName = open_spectra.getShiftLists()[0]
  windowOpt = 1
  for ii, specName in enumerate(specNames):
    path = paths[ii]
    optString = OpenSpectrum.WINDOW_OPTS[windowOpt]
    textMatrix.append([exptName, specName, path, optString, shiftListName])
    objectList.append([exptName, specName, path, windowOpt, shiftListName])

  if len(fileNames) > 1:
    open_spectra.openButton.config(text='Open Spectra')
 
  else:
    open_spectra.openButton.config(text='Open Spectrum')

  open_spectra.scrolledMatrix.update(objectList=objectList,
                                     textMatrix=textMatrix)
    
  if useScalingFactors:
    extraData['scalingFactors'] = scalingFactors[0]
    
  refSpec = open_spectra.openSpectrum(exptName, specNames[0], fileNames[0],
                                      extraData=extraData)
  if not refSpec:
    return
  
  refExp = refSpec.experiment
  
  # set up internal Analysis data
  open_spectra.parent.finishInitSpectrum(refSpec)
  print 'finished opening spectrum', refExp.name, refSpec.name

  refDataStore = refSpec.dataStore
  preferDataUrls = [refDataStore.dataUrl]
  memopsRoot = refSpec.root
  expDimRefs = [xdr for xd in refExp.sortedExpDims() 
                    for xdr in xd.sortedExpDimRefs()]
 
  # 
  for ii, fullPath in enumerate(fileNames):
    if ii == 0:
      # skip first
      continue
 
    if useScalingFactors:
      scalingFactor = scalingFactors[ii]
    else:
      scalingFactor = None
 
    # get new DataLocationStore
    useDataUrl, filePath = getDataStoringFromFilepath(memopsRoot, fullPath,
                          preferDataUrls=preferDataUrls,
                          keepDirectories=keepDirectories)
    
    
    # get data file name from entry path
    if format == 'Bruker':
      from ccp.format.bruker.generalIO import getMatrixFilePath
      urlPath = os.path.join(useDataUrl.url.dataLocation, '')
      dataPath = getMatrixFilePath(os.path.join(urlPath, filePath), 
                                   refDataStore.numDims)
      path = dataPath[len(urlPath):]
    elif format == 'Azara':
      raise NotImplementedError(
             "Azara type projection series not yet implemented"
            )
    else:
      path = filePath
    copyPars = {'dataUrl':useDataUrl, 'path':path, 'nmrDataSources':[]}
    newDataStore = copySubTree(refDataStore, refDataStore.dataLocationStore,
                               topObjectParameters=copyPars)
 
    # get new spectrum
    copyPars = {'name':specNames[ii],'dataStore':newDataStore}
    newSpectrum = copySubTree(refSpec, refExp, topObjectParameters=copyPars)
    
    if useScalingFactors:
      dataDims = newSpectrum.sortedDataDims()
      lastDataDim = dataDims[-1]
      startAt = len(dataDims) - 1
      for jj in range(len(dataDims) - 1, len(scalingFactor)):
        expDimRef = expDimRefs[jj]
        dimScaling = lastDataDim.findFirstDimensionScaling(expDimRef=expDimRef)
        if dimScaling is None:
          lastDataDim.newDimensionScaling(scalingFactors=(scalingFactor[jj],),
                                          expDimRef = expDimRef)
        else:
          dimScaling.scalingFactors = (scalingFactor[jj],)
    # set up internal Analysis data
    open_spectra.parent.visibleSpectra[newSpectrum] = False
    open_spectra.parent.finishInitSpectrum(newSpectrum)
    print 'finished opening spectrum', refExp.name, newSpectrum.name
  
  print
  print 'Projection spectra loaded into experiment', refExp.name
