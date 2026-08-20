"""
======================COPYRIGHT/LICENSE START==========================

CcpnToAriaXml.py: Extract information from a CCPN project to make an ARIA2
                  XML project file. Requires CCPN & ARIA2 installations.

Copyright (C) 2009 Tim Stevens and Wolfgang Rieping (University of Cambridge)

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

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
Rieping W., Habeck M., Bardiaux B., Bernard A., Malliavin T.E.,
Nilges M.(2007) ARIA2: automated NOE assignment and data
integration in NMR structure calculation. Bioinformatics 23:381-382

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

import os, sys

ARIA_ROOT_DIR = os.environ.get('ARIA2')
FAIL = 'CcpnToAriaXml Failure: '

# # # # # # # CHECK ARIA2 INSTALLATION # # # # # # #

def checkAriaInstallation():

  if not ARIA_ROOT_DIR:
    raise Exception(FAIL + '"ARIA2" environment variable not set')

  elif not os.path.exists(ARIA_ROOT_DIR):
    raise Exception(FAIL + 'ARIA2 path does not exist')

  if ARIA_ROOT_DIR not in sys.path:
    sys.path.append(ARIA_ROOT_DIR)

  try:
    import aria2

  except ImportError:
    raise Exception(FAIL + 'Cannot import ARIA 2 modules:\n' + \
                    'Check installation and ARIA2 environment variable')

  for module in aria2.PATH_MODULES:
    path = os.path.join(ARIA_ROOT_DIR, module)

    if path not in sys.path:
      sys.path.append(path)


  failedModules = []
  for module in aria2.MODULES:
    try:
      __import__(module)

    except ImportError, errMsg:
      failedModules.append(module)
      print errMsg

  if failedModules:
    text = ', '.join(failedModules)
    raise Exception(FAIL + 'Cannot import modules required by ARIA 2:\n%s' % text)

# # # # # # # CHECK CCPN INSTALLATION # # # # # # #

try:
  from memops.api.Implementation import MemopsRoot

except ImportError:
  raise Exception(FAIL + 'Cannot import CCPN modules:\n' + \
                  'Check CCPN installation is on your PYTHONPATH')

from ccp.util.NmrCalc import getObjBooleanParameter, getRunTextParameter, setRunParameter

# # # # # # #

PARAM_ATTR_DICT = {type(1.0):'floatValue',
                   type(1):'intValue',
                   type(True):'booleanValue'}

FILTER_VIOL = 'FilterViol'
KEEP_ASSIGN = 'KeepAssign'
AMBIG_PROTOCOL = 'AmbigProtocol'
USE_IN_CALC = 'UseInCalc'

CNS_EXE = 'CnsExe'
FILE_PREFIX = 'FilePrefix'
WORKING_DIR = 'WorkingDir'
TEMP_DIR = 'TempDir'

ARIA = 'ARIA'

# # # # # # #

def makeAriaProject(ccpnProject, ariaProjectPath, workingDir=None, template=None,
                    tempDir=None, cnsExe=None, calcStore=None, run=None):
  """ Create an ARIA projectXML file form a CCPN MemopsRoot (project) object.
      Options to pass-in specific CCPN objects to determine which ARIA settings
      are to be extracted.

      Input: memops.Implementation.MemopsRoot, string, string, string,
             string, string, molsim.NmrSim.NmrSimStore, molsim.NmrSim.Run
  """
  if not calcStore:
    if run:
      calcStore = run.nmrCalcStore
    else:
      calcStore = ccpnProject.findFirstNmrCalcStore(name=ARIA)

  if not calcStore:
    raise Exception(FAIL + 'No ARIA settings in CCPN project')
    return

  runs = calcStore.sortedRuns()
  if runs and not run:
    run = runs[-1]

  if not run:
    raise Exception(FAIL + 'No ARIA run setup')
    return

  #peakListsData = []
  #constraintData = []
  peakListsData = None
  constraintData = None
  molSystemData = None

  for dataObj in run.sortedInputs():
    className = dataObj.className

    if className == 'PeakListData':
      peakList = dataObj.peakList
      if peakList:
        #peakListsData.append((dataObj, peakList))
        peakListsData = True

    elif className in ('MolSystemData', 'MolResidueData'):
      molSystem = dataObj.molSystem
      if molSystem:
        molSystemData = (dataObj, molSystem)

    elif className == 'ConstraintStoreData':
      nmrConstraintStore = dataObj.nmrConstraintStore
      #if nmrConstraintStore:
      #  constraintLists = dataObj.constraintLists or \
      #                    nmrConstraintStore.sortedConstraintLists()

      #  for constraintList in constraintLists:
      #    constraintData.append((dataObj, constraintList))
      if nmrConstraintStore and (dataObj.constraintLists or 
                                 nmrConstraintStore.constraintLists):
          constraintData=True

  if not molSystemData:
    msg = 'No molecular system selected for ARIA run'
    raise Exception(FAIL + msg)
    return

  if not (constraintData or peakListsData):
    msg = 'No constraint or peak list data selected for ARIA run'
    raise Exception(FAIL + msg)
    return

  if workingDir:
    setRunParameter(run, WORKING_DIR, workingDir)

  if cnsExe:
    setRunParameter(run, CNS_EXE, cnsExe)

  if tempDir:
    setRunParameter(run, TEMP_DIR, tempDir)

  writeAriaProject(run, ariaProjectPath, ARIA_ROOT_DIR, template=template)


def getObjectKeyString(object, delimiter='|'):
  """Descrn: Make an object identifier string for a CCPN data model object
     Inputs: CCPN data model object
     Output: String
  """

  keys = object.getExpandedKey()

  for i in range(len(keys)):
    key = keys[i]

    keyType = type(key)
    if keyType is type([]):
      keys[i] = delimiter.join([str(k) for k in key])
    elif keyType is not type(''):
      keys[i] = str(key)

  return delimiter.join(keys)


# # # # # # # Setup ARIA data structures


def addAriaGeneralData(project, run, externalTemplate):

  cnsExePath = getRunTextParameter(run, CNS_EXE)

  if cnsExePath:
    settings = project.getSettings()

    ccpnProject = run.root
    url = ccpnProject.findFirstRepository(name='userData').url
    ccpnProjPath = url.path

    generalData = {'name':  ARIA+'_'+ccpnProject.name,
                   'file_root': getRunTextParameter(run, FILE_PREFIX) or ccpnProject.name,
                   'working_directory': getRunTextParameter(run, WORKING_DIR) or ccpnProjPath,
                   'temp_root': getRunTextParameter(run, TEMP_DIR) or ccpnProjPath,
                   'run': '%d' % run.serial,}

    settings.update(generalData)

    if os.path.exists(cnsExePath) and not externalTemplate:
      cns = project.getStructureEngine()
      cns.getSettings()['local_executable'] = cnsExePath


def addAriaCcpnData(project, run):

  import aria.DataContainer as DC
  from aria.ariabase import YES, NO

  ccpnProject = run.root
  url = ccpnProject.findFirstRepository(name='userData').url

  ccpn_export = {'export_assignments'        : YES,
                 'export_noe_restraint_list' : 'last',
                 'export_structures'         : YES}

  reporter = project.getReporter()
  reporter['ccpn'].update(ccpn_export)

  d = DC.CCPNData()
  d.reset()

  d['filename'] = url.path

  project.ccpn_model = d

def addAriaMolSystem(project, run):

  import aria.DataContainer as DC

  seq_data = project.getData(DC.DATA_SEQUENCE)[0]
  seq_data.reset()

  dataObj = run.findFirstData(className='MolSystemData', ioRole='input')
  if dataObj:
    chains = dataObj.chains
  else:
    dataObj = run.findFirstData(className='MolResidueData', ioRole='input')
    if dataObj:
      chains = [dataObj.chain]
  
  if not dataObj:
    raise Exception(FAIL + 'No MolSystem in CCPN run')
  else:
    molSystem = dataObj.molSystem

  #chains = dataObj.chains
  
  if not chains:
    chains = molSystem.sortedChains()
  
  chainKeys = [getObjectKeyString(molSystem),] + [ch.code for ch in chains]

  seq_data['format'] = DC.DATA_CCPN
  seq_data['ccpn_id'] = '|'.join(chainKeys)

def addAriaPeakList(project, run):

  import aria.DataContainer as DC
  from aria.ariabase import YES, NO

  for dataObj in run.sortedInputs():
    className = dataObj.className

    if className == 'PeakListData':
      peakList = dataObj.peakList
      if peakList:
        shiftList = peakList.dataSource.experiment.shiftList
        keepAssign = getObjBooleanParameter(dataObj, KEEP_ASSIGN)
        filterViol = getObjBooleanParameter(dataObj, FILTER_VIOL)

        shift_data = DC.ShiftData()
        shift_data.reset()

        shift_data['filename'] = ''
        shift_data['format'] = DC.DATA_CCPN
        shift_data['ccpn_id'] = getObjectKeyString(shiftList)

        peak_data = DC.PeakData() ## using default window size
        peak_data.reset()

        peak_data['filename'] = ''
        peak_data['format'] = DC.DATA_CCPN
        peak_data['ccpn_id'] = getObjectKeyString(peakList)

        spectrum_data = DC.SpectrumData()
        spectrum_data['shifts'] = shift_data
        spectrum_data['peaks'] = peak_data
        spectrum_data['use_assignments'] = keepAssign and YES or NO
        spectrum_data['trust_assigned_peaks'] = filterViol and YES or NO

        project.addData(spectrum_data)

def addAriaConstraints(project, run):

  import aria.DataContainer as DC
  from aria.ariabase import YES, NO

  for dataObj in run.sortedInputs():
    className = dataObj.className

    if className == 'ConstraintStoreData':
      nmrConstraintStore = dataObj.nmrConstraintStore

      if nmrConstraintStore and getObjBooleanParameter(dataObj, USE_IN_CALC):
        constraintLists = dataObj.constraintLists or \
                          nmrConstraintStore.sortedConstraintLists()

        for constraintList in constraintLists:
          className = constraintList.className
          key = getObjectKeyString(constraintList)
          dataDict = None

          if className == 'DistanceConstraintList':
            ambigProtocol = getObjBooleanParameter(dataObj, AMBIG_PROTOCOL)
            filterViol = getObjBooleanParameter(dataObj, FILTER_VIOL)

            if ambigProtocol:
              dataDict = DC.AmbiguousDistanceData()
              dataDict.reset()
              dataDict['filename'] = ''
              dataDict['format'] = DC.DATA_CCPN
              dataDict['ccpn_id'] = key
              dataDict['filter_contributions'] = filterViol and YES or NO
              dataDict['calibrate'] = 'all_iterations_except_first'

            else:
              dataDict = DC.UnambiguousDistanceData()
              dataDict.reset()
              dataDict['filename'] = ''
              dataDict['format'] = DC.DATA_CCPN
              dataDict['ccpn_id'] = key
              dataDict['filter_contributions'] = filterViol and YES or NO
              dataDict['calibrate'] = 'all_iterations_except_first'

          elif className == 'JCouplingConstraintList':
            dataDict = DC.KarplusData()
            dataDict.reset()
            dataDict['filename'] = ''
            dataDict['format'] = DC.DATA_CCPN
            dataDict['ccpn_id'] = key

          elif className == 'DihedralConstraintList':
            dataDict = DC.DihedralData()
            dataDict.reset()
            dataDict['filename'] = ''
            dataDict['format'] = DC.DATA_CCPN
            dataDict['ccpn_id'] = key

          elif className == 'HBondConstraintList':
            dataDict = DC.HBondData()
            dataDict.reset()
            dataDict['filename'] = ''
            dataDict['format'] = DC.DATA_CCPN
            dataDict['ccpn_id'] = key
            dataDict['type'] = 'standard' # Not 'csi'

          if dataDict:
            project.addData(dataDict)

def writeAriaProject(run, filename, rootPath, template=None):

  checkAriaInstallation()

  from aria.AriaXML import AriaXMLPickler

  if not template:
    from aria.ariabase import PROJECT_TEMPLATE
    template = os.path.join(rootPath, 'src', 'py', 'data', PROJECT_TEMPLATE)
    externalTemplate = False
  else:
    externalTemplate = True

  pickler = AriaXMLPickler()

  project = pickler.load_relaxed(template)

  addAriaGeneralData(project, run, externalTemplate)
  addAriaCcpnData(project, run)
  addAriaMolSystem(project, run)
  addAriaPeakList(project, run)
  addAriaConstraints(project, run)

  pickler.dump(project, filename)

# # # # # #  Command line

if __name__ == '__main__':
  checkAriaInstallation()

  cmdArgs = sys.argv[1:]
  nArgs = len(cmdArgs)

  if nArgs < 2:
    print """
      CcpnToAriaXml requires at least two command line arguments:
        CCPN project directory
        Output ARIA XML file name

      Additional arguments may include:
        ARIA working directory (default to current location)
        ARIA temp directory (default to current location)

      Example:

      python CcpnToAriaXml.py /data/ccpnProjDir /home/me/AriaProj.xml /home/me/ariaRuns /temp
      """

    sys.exit(0)

  ccpnProjectDir, ariaProjectPath = cmdArgs[:2]

  workingDir=None
  tempDir=None

  if nArgs > 2:
    workingDir = cmdArgs[2]
  if nArgs > 3:
    tempDir = cmdArgs[3]

  from memops.general.Io import loadProject

  if not os.path.exists(ccpnProjectDir):
    raise Exception(FAIL + 'CCPN project directory does not exist')

  if not os.path.isdir(ccpnProjectDir):
    raise Exception(FAIL + 'CCPN project location is not a directory')

  try:
    ccpnProject = loadProject(ccpnProjectDir)
  except Exception, err:
    raise Exception(FAIL + 'CCPN project failed to load. original error:' + err)


  makeAriaProject(ccpnProject, ariaProjectPath,
                  workingDir=workingDir, tempDir=tempDir)

  print 'Done. Saved ARIA project to: %s' % ariaProjectPath
