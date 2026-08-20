
"""
======================COPYRIGHT/LICENSE START==========================

Util.py: code for CCPN data model and code generation framework

Copyright (C) 2012  (CCPN Project)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

"""
from ccpnmr.integrator.core import Util as intUtil
from ccp.lib.nmr.Nmr.DataSource import getIsotopeCodesList



defaultConfiguration = {
  "protocol": {
    "protocolParameters": [ 
      {"name":"fileNameSequence","paramType":"String","value":"protein.seq","code":"FileNameSequence"},
      {"name":"fileFormatSequence","paramType":"String","value":"sparky","code":"FileFormatSequence"},
      {"name":"peakFormat","paramType":"String","value":"nmrview"},
      {"name":"shiftFormat","paramType":"String","value":"xeasy"},
      {"name":"atomNamingSystem","paramType":"String","value":"DIANA"},
      {"name":"skipMultiAssignments","paramType":"Boolean","value":True}

    ]
  }
}


preferredShiftWeights = {
 'hHC':8,
 'cHC':4,
 'NOESY_aliphatic':3,
 'NOESY_all':2,
 'NOESY_aromatic':1,
}

spectrumTypes = {
 'all': 'NOESY_all',
 'aliphatic': 'NOESY_aliphatic',
 'aromatic': 'NOESY_aromatic'
}

def getSpectrumKind(dataSource):
  """ Get spectrum classification
  """
  
  typeMap = {
   ('1H','1H'):'hh',
   ('13C','13C'):'cc',
   ('15N','15N'):'nn',
   ('13C','1H','1H'):'hHC',
   ('15N','1H','1H'):'hHN',
   ('13C','13C','1H'):'cHC',
  }
  
  isotopeCodes = getIsotopeCodesList(dataSource)
  if None in isotopeCodes or [x for x in isotopeCodes if ',' in x]:
    return None
  
  experiment = dataSource.experiment
  # check that spectrum type is through-space
  refExperiment = experiment.refExperiment
  if refExperiment: 
    if refExperiment.nmrExpPrototype.category != 'through-space':
      return None
  elif experiment.findFirstExpTransfer(transferType='through-space')is None:
    return None
  
  #
  return typeMap.get(tuple(sorted(isotopeCodes)))

def adaptNmrCalcRun(nmrCalcRun, protocolName):
  """ Modify nmrCalcRun from generic MultiStructure form to protocol-specific
    Input:
      nmrCalcRun: NmrCalc.Run 
      String protocolName
  """
  
  if protocolName.startswith('UNIO'):
    
    # NBNB TBD Only works for UNIO_CANDID now.
    
    candidNmrCalcRun(nmrCalcRun)
  
  else:
    raise Exception("Protocol %s not recognized" % protocolName)

def candidNmrCalcRun(nmrCalcRun):
  """ Modify nmrCalcRun from generic MultiStructure form to CANDID-specific
    Input:
      nmrCalcRun: NmrCalc.Run 
  """
  
  xx = nmrCalcRun.findFirstRunParameter(name='shiftFormat')
  if xx is None:
    raise Exception("shiftFormat must be set for UNIO_CANDID protocols")
  else:
    shiftFormat = xx.textValue
  xx = nmrCalcRun.findFirstRunParameter(name='peakFormat')
  if xx is None:
    raise Exception("peakFormat must be set for UNIO_CANDID protocols")
  else:
    peakFormat = xx.textValue
  shiftFileSuffix = intUtil.getFileFormatData(shiftFormat, 'shift')['shiftExt']
  peakFileSuffix = intUtil.getFileFormatData(peakFormat, 'peak')['peakExt']
 
  # name : code dictionary for resetting code values
  codeRemapping = {
   'solvent':'SolventType',
   'spectrumType':'PeakType',
   'spectrumKind':'PeakKind',
   'numDataSets':'NoPeakList',
  }
  
  # set proper codes for existing data parameters
  for runParameter in nmrCalcRun.sortedRunParameters():
    newCode = codeRemapping.get(runParameter.name)
    if newCode is not None:
      runParameter.code = newCode
  
  # Reset spectrumType strings
  for runParameter in nmrCalcRun.findAllRunParameters(name='spectrumType'):
    ss = runParameter.textValue
    if ss in spectrumTypes:
      # Spectrum type should be remapped
      runParameter.textValue = spectrumTypes[ss]
  
  # get data sets
  peakListObjs = [x for x in nmrCalcRun.sortedData()
                  if x.className == 'PeakListData']
  
  shiftListObjs = {}
  for datum in nmrCalcRun.findAllData(className='MeasurementListData'):
    measurementList = datum.measurementList
    if measurementList.className == 'ShiftList':
      shiftListObjs[measurementList] = datum
  
  # create data set files and parameters
  numDataSets = 0
  usePeakListObjs = []
  for peakListData in peakListObjs:
    numDataSets += 1
    
    # check  and set spectrum kind
    dataSource = peakListData.dataSource
    spectrumKind = getSpectrumKind(dataSource)
    if spectrumKind is None:
      print 'WARNING, incorrect spectrum type. Removing %s' % peakListData
      for runParameter in peakListData.runParameters:
        runParameter.delete()
      peakListData.delete()
      numDataSets -= 1
      continue
    
    # note in order to select shiftlist later
    usePeakListObjs.append(peakListData)
    
    # Reset connected parameters
    #for runParameter in peakListData.runParameters:
    #  code = runParameter.code
    #  if code:
    #    runParameter.code = code + '[%s]' % numDataSets
    
    nmrCalcRun.newRunParameter(name='spectrumKind', data=peakListData,
                               code=codeRemapping['spectrumKind'],# + '[%s]' % numDataSets),
                               textValue=spectrumKind)
    
    # check and set shiftlist and related parameters
    experiment = dataSource.experiment
    shiftList = experiment.shiftList
    shiftListData = shiftListObjs.get(shiftList)
    if shiftListData is None:
      shiftListData = nmrCalcRun.newMeasurementListData(name='shiftList',
                       measurementList=shiftList)
      shiftListObjs[shiftList] = shiftListData
      shiftFile = intUtil.objectFileName(shiftList, suffix=shiftFileSuffix,
                                         compact=True)
      nmrCalcRun.newRunParameter(name='fileName', textValue=shiftFile, 
                                 data=shiftListData)
      nmrCalcRun.newRunParameter(name='fileFormat', textValue=shiftFormat, 
                                 data=shiftListData)
      useFormat = shiftFormat
    else:
      shiftFile = shiftListData.findFirstRunParameter(name='fileName').textValue
      useFormat = shiftListData.findFirstRunParameter(name='fileFormat').textValue
    
    nmrCalcRun.newRunParameter(name='fileNameShift', textValue=shiftFile, 
                               code='fileNameShift', 
                               data=peakListData)
    nmrCalcRun.newRunParameter(name='fileFormatShift', textValue=shiftFormat, 
                               code='fileFormatShift', 
                               data=peakListData)
    
    fileNamePeakList = intUtil.objectFileName(peakListData.peakList, 
                                              suffix=peakFileSuffix,
                                              compact=True)
    nmrCalcRun.newRunParameter(name='fileName', textValue=fileNamePeakList, 
                               code='fileNamePeak', 
                               data=peakListData)
    nmrCalcRun.newRunParameter(name='fileFormat', textValue=peakFormat, 
                               code='fileFormatPeak', 
                               data=peakListData)
    
  # number of data sets
  nmrCalcRun.newRunParameter(name='numDataSets', 
                             code=codeRemapping['numDataSets'],
                             intValue = numDataSets)
  
  # choose shift list for backbone
  if shiftListObjs:
    if len(shiftListObjs) == 1:
      shiftData = list(shiftListObjs.items())[0][1]
      fileName = shiftData.findFirstRunParameter(name='fileName').textValue
      fileFormat = shiftData.findFirstRunParameter(name='fileFormat').textValue
    else:
      shiftListScore = 0
      for peakListData in usePeakListObjs:
        # check preferred shiftlist
        spectrumKind = peakListData.findFirstRunParameter(
                         name='spectrumKind').textValue
        spectrumType = peakListData.findFirstRunParameter(
                         name='spectrumType').textValue
        shiftFileWeight = (preferredShiftWeights.get(spectrumKind,0) +
                           preferredShiftWeights.get(spectrumType,0))
        if shiftFileWeight > shiftListScore:
          shiftListScore = shiftFileWeight
          fileName = peakListData.findFirstRunParameter(
                       name='fileNameShift').textValue
          fileFormat = peakListData.findFirstRunParameter(
                         name='fileFormatShift').textValue
 
    # main shift list
    nmrCalcRun.newRunParameter(name='fileNameShift', code='fileNameCaCbShift',
                               textValue=fileName)
    nmrCalcRun.newRunParameter(name='fileFormatShift', code='fileFormatCaCbShift',
                               textValue=fileFormat)
 
  
