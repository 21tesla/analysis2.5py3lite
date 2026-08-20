
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
from ccpnmr.analysis.core import AssignmentBasic

from ccp.lib.nmr.Nmr.DataSource import getIsotopeCodesList, getDimCodes


spectrumTypes = {
 'aliphatic': 'C13-ali',
 'aromatic': 'C13-aro'
}



def adaptNmrCalcRun(nmrCalcRun, protocolName):
  """ Modify nmrCalcRun from generic MultiStructure form to protocol-specific
    Input:
      nmrCalcRun: NmrCalc.Run 
      String protocolName
  """
  
  if protocolName.startswith('ASDP'):
    asdpNmrCalcRun(nmrCalcRun)
  
  else:
    raise Exception("Protocol %s not recognized" % protocolName)

def asdpNmrCalcRun(nmrCalcRun):
  """ Modify nmrCalcRun from generic MultiStructure form to ASDP-specific
    Input:
      nmrCalcRun: NmrCalc.Run 
  """
  
  # CCPN params
  xx = nmrCalcRun.findFirstRunParameter(name='shiftFormat')
  if xx is None:
    raise Exception("shiftFormat must be set for ASDP protocols")
  else:
    shiftFormat = xx.textValue
  xx = nmrCalcRun.findFirstRunParameter(name='peakFormat')
  if xx is None:
    raise Exception("peakFormat must be set for ASDP protocols")
  else:
    peakFormat = xx.textValue
  shiftFileSuffix = intUtil.getFileFormatData(shiftFormat, 'shift')['shiftExt']
  peakFileSuffix = intUtil.getFileFormatData(peakFormat, 'peak')['peakExt']
 
  # name : code dictionary for resetting code values
  codeRemapping = {
   #'spectrumType':'x1.type',
   'jobName':'proteinName',
   #'numDataSets':None,
  }
  
  # set proper codes for existing data parameters
  for runParameter in nmrCalcRun.sortedRunParameters():
    name = runParameter.name
    if name in codeRemapping:
      runParameter.code = codeRemapping[name]
  
  # get data sets
  peakListObjs = [x for x in nmrCalcRun.sortedData()
                  if x.className == 'PeakListData']
  
  shiftListObjs = {}
  for datum in nmrCalcRun.findAllData(className='MeasurementListData'):
    measurementList = datum.measurementList
    if measurementList.className == 'ShiftList':
      shiftListObjs[measurementList] = datum
  
  defaultTol = 1.0
  runParam = nmrCalcRun.findFirstRunParameter(name='defaultTolPoints')
  if runParam:
    val = runParam.floatValue
    if val:
      defaultTol = val
  
  minTol = 0.02
  runParam = nmrCalcRun.findFirstRunParameter(name='minTolPpm')
  if runParam:
    val = runParam.floatValue
    if val:
      minTol = val
  
  # create data set files and parameters
  for peakListData in peakListObjs:
    dataSource = peakListData.dataSource
    experiment = dataSource.experiment
    
    dimCodes = getDimCodes(dataSource)
    
    if (None in dimCodes or 
        not experiment.findFirstExpTransfer(transferType='through-space')):
      print ("Removed spectrum %s: not recognised as 2D-4D NOESY" % dataSource)
      for runParameter in peakListData.runParameters:
        runParameter.delete()
      peakListData.delete()
      continue
    
    # peak file name
    fileNamePeakList = intUtil.objectFileName(peakListData.peakList, 
                                              suffix=peakFileSuffix,
                                              compact=True)
    nmrCalcRun.newRunParameter(name='fileName', textValue=fileNamePeakList, 
                               code='fileNamePeak', 
                               data=peakListData)
    nmrCalcRun.newRunParameter(name='fileFormat', textValue=peakFormat, 
                               data=peakListData)
    
    nmrCalcRun.newRunParameter(name='numDim', code='dimension', 
                               intValue=dataSource.numDim, 
                               data=peakListData)
    
    # interchain NOEs not supported
    nmrCalcRun.newRunParameter(name='haveInterChain', code='haveIC', intValue=0, 
                               data=peakListData)
    nmrCalcRun.newRunParameter(name='interChain', code='IC',intValue=0,  
                               data=peakListData)
    
    # Solvent
    runParameter = peakListData.findFirstRunParameter(name='solvent')
    if runParameter.textValue == 'h2o':
      val = 1
    else:
      val= 0
    nmrCalcRun.newRunParameter(name='solventH2O', code='waterFlag',intValue=val,  
                               data=peakListData)
    runParameter.delete()
    
    # column data
    nmrCalcRun.newRunParameter(name='intensityColumn', code='col.intensity',
                               intValue=dataSource.numDim + 4,  
                               data=peakListData)

    isotopeCodes = getIsotopeCodesList(dataSource)
    tolerances = AssignmentBasic.estimateAssignmentTolerances(dataSource,
                                           defPoints=defaultTol, minTol=minTol)

    # Dimension parameters:
    for ii, dataDim in enumerate(dataSource.sortedDataDims()):
      dimCode = dimCodes[ii]
      
      
      # Peak position column column
      tag = 'col.%s' % dimCode
      nmrCalcRun.newRunParameter(name='PositionColumn', code=tag, intValue=ii+2,  
                                 data=peakListData)
      
      # Shift correction (set to default)
      tag = '%s.shift' % dimCode
      nmrCalcRun.newRunParameter(name='ShiftCorrection', code=tag, 
                                 floatValue=0.0, data=peakListData)
      
      # sign (meaning ???) (set to default)
      tag = '%s.sign' % dimCode
      nmrCalcRun.newRunParameter(name=tag, code=tag, intValue=0,  
                                 data=peakListData)
      
      # Assignment tolerance
      tag = '%s.tol' % dimCode
      nmrCalcRun.newRunParameter(name='AssignTolerance', code=tag, 
                                 floatValue=tolerances[ii], data=peakListData)
      
      # sw in ppm
      tag = '%s.sw' % dimCode
      swOrig = 1000.0
      # NBNB program probably assumes that all shifts are within SW.
      # As long as peak positions are correct, swOrig is not useful. 
      # Leave at default
      #for dataDimRef in dataDim.sortedDataDimRefs():
      #  if dataDimRef.expDimRef.measurementType in ('shift','Shift'):
      #    swOrig = dataDimRef.spectralWidthOrig
      
      nmrCalcRun.newRunParameter(name='swPpm', code=tag, 
                                 floatValue=swOrig, data=peakListData)
      
      # dimension type
      isotopeCode = isotopeCodes[ii]
      if isotopeCode == '1H':
        dimType = 'H'
        
      elif isotopeCode == '13C':
        dimType = 'C13'
        specTypeObj = peakListData.findFirstRunParameter(name='spectrumType')
        if specTypeObj:
          dimType = spectrumTypes.get(specTypeObj.textValue) or dimType
        
      else:
        ii = 0
        while isotopeCode[ii] in '0123456789':
          ii += 1
        dimType = isotopeCode[ii:] + isotopeCode[:ii]
      
      tag = '%s.type' % dimCode
      nmrCalcRun.newRunParameter(name='dimensionType', code=tag, 
                                 textValue=dimType, data=peakListData)
    
    # check and set shiftlist and related parameters
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
    
  # main shift list
  if shiftListObjs:
    if len(shiftListObjs) == 1:
      shiftData = list(shiftListObjs.items())[0][1]
      fileName = shiftData.findFirstRunParameter(name='fileName').textValue
      fileFormat = shiftData.findFirstRunParameter(name='fileFormat').textValue
      nmrCalcRun.newRunParameter(name='fileNameShift', code='chemicalShiftFile',
                                 textValue=fileName, data=shiftData)
      nmrCalcRun.newRunParameter(name='fileFormatShift', code='fileFormatShift',
                                 textValue=fileFormat, data=shiftData)
    else:
      raise Exception("All spectra must use same shift list for ASDP")
 
 
  
