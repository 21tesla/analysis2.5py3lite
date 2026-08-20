
"""
======================COPYRIGHT/LICENSE START==========================

Util.py: code for CCPN data model and code generation framework

Copyright (C) 2011  (CCPN Project)

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
import operator, copy

# from memops.general import Io as genIo
from memops.general import Util as genUtil

from ccp.general import Constants as genConstants

from ccpnmr.analysis.core import ConstraintBasic
from ccpnmr.analysis.core import AssignmentBasic

from ccp.lib.nmr.Nmr.DataSource import getDimCodes

from ccpnmr.format.converters import DataFormat
from ccpnmr.format.converters.XEasyFormat import XEasyFormat
from ccpnmr.format.converters.NmrViewFormat import NmrViewFormat
from ccpnmr.format.converters.SparkyFormat import SparkyFormat
from ccpnmr.format.converters.CyanaFormat import CyanaFormat
from ccpnmr.format.converters.CnsFormat import CnsFormat
# from ccpnmr.format.converters.AnsigFormat import AnsigFormat
# from ccpnmr.format.converters.NmrDrawFormat import NmrDrawFormat

from ccpnmr.integrator.plugins.Talos.Io import TalosFormat
from ccpnmr.integrator.plugins.NmrStar.Io import NmrStarWrapper


######################################################################
#
#               plugins and default file and format data
#
######################################################################

# integrator plugin names
pluginDirs = {
 'ROSETTA':'Rosetta',
 'ASDP':'Asdp',
 'UNIO':'Unio',
 'CYANA':'Cyana',
 'ARIA':'Aria',
}

# file data depending on format and file type
# WARNING Do *not* access _fileFormatData directly, use getFileFormatData

# First FormatConverter using formats
_fileFormatData = {
 'CYANA':{
  'handlerClass':CyanaFormat,
  'seqExt':'seq',
  'shiftExt':'prot',
  'peakExt':'peaks',
  'IOkeywords':{}
 },
 'XEASY':{
  'handlerClass':XEasyFormat,
  'seqExt':'seq',
  'shiftExt':'prot',
  'peakExt':'peaks',
  'IOkeywords':{}
 },
 'NMRVIEW':{
  'handlerClass':NmrViewFormat,
  'seqExt':'seq',
  'shiftExt':'ppm',
  'peakExt':'xpk',
  'IOkeywords':{}
 },
 'NMRSTAR':{
  'handlerClass':NmrStarWrapper,
  'seqExt':'bmrb',
  'shiftExt':'bmrb',
  'peakExt':'bmrb',
  'IOkeywords':{}
 },
 'SPARKY':{
  'handlerClass':SparkyFormat,
  'seqExt':'seq',
  'shiftExt':'ppm',
  'peakExt':'peaks',
  'IOkeywords':{}
 },
 'CNS':{
  'handlerClass':CnsFormat,
  'constraintExt':'tbl',
  'IOkeywords':{}
 },
}

mainTags =  [x for x in _fileFormatData.keys() if '_' not in x]

# Set default values.
# NBNB presetResonanceMapping==True is necessary for integrator
# so is compressResonances==False
# NBNB only tested for NmrView peak format
defaultDict = {

 # Mandatory for integrator to work
 'useCcpnChainInfo':True,
 'compressResonances':False,
 'presetResonanceMapping':True,

 # optional - could be changed
 'individualAtomsIfNoSet':True,   # NB consequences may be unpredictable
 'minimalPrompts':True,
 'individualAtoms':False,
 'ignoreChemCompsWithNoSysNames':False,
 'verbose':True,
}
for tag in mainTags:
  dd = _fileFormatData[tag]['IOkeywords']
  for key,val in defaultDict.items():
    if key not in dd:
      dd[key] = val

  # NBNB TBD double check settings ands system for general use

  # NBNB TBD

  # Obsolete. Left here to show how to do peak-specific defaults
  #if 'peakExt' in dd:
  #  # Format has peaks
  #  newTag = tag + '_peak'
  #  _fileFormatData[newTag] = dd2 = copy.deepcopy(dd)
  #  keywds = dd2['IOkeywords']
  #  keywds['compressResonances'] = False


# Special CYANA variants
dd = _fileFormatData['CYANA_peak'] =  copy.deepcopy(_fileFormatData['CYANA'])
dd['IOkeywords']['cyanaFormat'] = True
dd['handlerClass'] = XEasyFormat

dd = _fileFormatData['CYANA_shift'] =  copy.deepcopy(_fileFormatData['CYANA'])
dd['IOkeywords']['cyana21Naming'] = True
dd['handlerClass'] = XEasyFormat

dd = _fileFormatData['CYANA_distance'] =  copy.deepcopy(_fileFormatData['CYANA'])
dd['IOkeywords']['constraintType'] = 'distance'

dd = _fileFormatData['CYANA_dihedral'] =  copy.deepcopy(_fileFormatData['CYANA'])
dd['IOkeywords']['constraintType'] = 'dihedral'

dd = _fileFormatData['CYANA_hBond'] =  copy.deepcopy(_fileFormatData['CYANA'])
dd['IOkeywords']['constraintType'] = 'hBond'

# Non-FormatConverter-using formats
_fileFormatData['TALOS'] = {
 'handlerClass':TalosFormat,
 'shiftExt':'tab',
 'IOkeywords':{'minShiftQuality':0.1}
}



######################################################################
#
#               JSON file tags and maps used for WMS I/O.
#
######################################################################


classGenerators = {
 'MolSystem': 'newMolSystemData',
 'ConstraintStore': 'newConstraintStoreData',
 'PeakList': 'newPeakListData',
 'MeasurementList': 'newMeasurementListData',
 'MolResidue': 'newMolResidueData',
 'ViolationList': 'newViolationListData',
 'Spectrum': 'newSpectrumData',
 'DerivedList': 'newDerivedListData',
 'StructureEnsemble': 'newStructureEnsembleData',
 'External': 'newExternalData',
 'FloatMatrix': 'newFloatMatrixData',
 'SpinSystem':'newData',
 'Tensor':'newTensorData',
 'EnergyTerm':'newEnergyTermData',
}

runDataTags = {
 'MolSystem': frozenset(('serial', 'name', 'code', 'ioRole', 'details',
                         'molSystemCode', 'symmetrySetId','chainCodes')),
 'ConstraintStore': frozenset(('serial', 'name', 'code', 'ioRole', 'details',
                               'constraintStoreSerial',
                               'constraintListSerials')),
 'PeakList': frozenset(('serial', 'name', 'code', 'ioRole', 'details',
                        'experimentSerial', 'dataSourceSerial',
                        'peakListSerial')),
 'MeasurementList': frozenset(('serial', 'name', 'code', 'ioRole', 'details',
                               'measurementListSerial')),
 'MolResidue': frozenset(('serial', 'name', 'code', 'ioRole', 'details',
                          'molSystemCode','chainCode', 'residueSeqIds',
                          'chainCodes')),
 'ViolationList': frozenset(('serial', 'name', 'code', 'ioRole', 'details',
                             'constraintStoreSerial', 'violationListSerial')),
 'Spectrum': frozenset(('serial', 'name', 'code', 'ioRole', 'details',
                        'experimentSerial', 'dataSourceSerial')),
 'DerivedList': frozenset(('serial', 'name', 'code', 'ioRole', 'details',
                           'derivedDataListSerial')),
 'StructureEnsemble': frozenset(('serial', 'name', 'code', 'ioRole', 'details',
                                 'molSystemCode', 'ensembleId','modelSerials')),
 'External': frozenset(('serial', 'name', 'code', 'ioRole', 'details',
                        'dataLocationStoreName', 'dataStoreSerial')),
 'FloatMatrix': frozenset(('serial', 'name', 'code', 'ioRole', 'details',
                           'shape', 'data', 'defaultValue')),
 # TBD:
 #'SpinSystem':(),
 #'Tensor':(),
 #'EnergyTerm':(),
}

nonDataTags = {
 'RunParameter':('name', 'code', 'ioRole',
                 'textValue', 'intValue','floatValue', 'booleanValue'),
 'ParameterGroup':()
}



# Matches WMS  DB as of 11 Oct 2012 (see src/uk/ac/ccpn/wms/admin/InitialiseDb.java)
paramTypeMap = {
'String':'textValue',
'Integer':'intValue',
'Float':'floatValue',
'Boolean':'booleanValue',
'ChainCodes':'molSysChains',
'MeasurementLists':'measurementList',
'Residues':'molResidues',
'Ensembles':'structureModels',
'PeakLists':'peakList',
'ConstraintLists':'constraintLists',
}

paramTypeConverters = {
 'textValue':str,
 'intValue':int,
 'floatValue':float,
 'booleanValue':lambda x: x in ('True','true')

}



######################################################################
#
#               Program data I/O and plugins
#
######################################################################


def getIntegratorPlugin(procName):
  """ Get plugin.
  """
  procName = procName.upper()
  for tag,dirName in pluginDirs.items():
    if tag in procName:
      pluginName = 'ccpnmr.integrator.plugins.' + dirName
      break
  else:
    raise Exception("No plugin found for protocol %s" % procName)
  #
  return __import__(pluginName, globals(), locals(), ['read', 'write', 'Util'])


def getFormatTag(dataFormat, dataType=None):
  """ Get data format tag for internal dictinary access
  """

  formatTag = dataFormat.upper()
  if dataType is not None:
    tag = '_'.join((formatTag,dataType))
    if tag in _fileFormatData:
      formatTag = tag
  #
  return formatTag

def getFileFormatData(format, dataType=None):
  """ Get File format data dictionary
  """

  return _fileFormatData.get(getFormatTag(format, dataType))

class DataIoHandler:
  """ Provides im/exporters for different formats and tasks
  NB nmrCalcRun is required only for some formats
  """

  def __init__(self, project, nmrCalcRun=None):
    self.ioHandlers = {}
    self.project = project
    self.nmrCalcRun = nmrCalcRun

  def resetAll(self):
    self.ioHandlers = {}

  def getIoHandler(self, dataFormat, dataType=None, reset=False):
    """ Get IOhandler of appropriate type for format and variant.
    NB we want the *same* handler object each time for a given handlerClass
    so that we can reuse the same resonanceToAtoms mapping
    """

    formatTag = getFormatTag(dataFormat, dataType)
    handlerClass = _fileFormatData[formatTag]['handlerClass']
    handler = self.ioHandlers.get(handlerClass)
    if handler is not None and not reset:
      return handler

    handler = handlerClass(self.project)
    self.ioHandlers[handlerClass] = handler

    # Special cases:
    if formatTag == 'TALOS':
      # Get residues to write shifts for
      ll = list(self.nmrCalcRun.findAllData(className='MolResidueData'))
      if len(ll) == 1:
        obj =  ll.pop()
        residues = list(obj.residues)
        residues.sort(key=operator.attrgetter('seqId'))
        if len(set(obj.chainCodes)) != 1:
          raise Exception("Run MolResidueData did not have a (single) chain" )
      else:
        raise Exception("Run must have excactly one MolResidueData, %s found" % len(ll))
      #
      if residues:
        handler.residues = residues

    #
    return handler



######################################################################
#
#               RunParameter value handling
#
######################################################################

def getParameterValue(param):
  """ Get value of NmrCalc.RunParameter
  """
  for tag in ('textValue', 'intValue', 'floatValue', 'booleanValue'):
    result = getattr(param,tag)
    if result is not None:
      break
  #
  return result

def setParameterValue(runParameter, value):
  """ Set runParameter value depending on value type.
  Resets "" to None.
  NB does *not* reset previous values
  """

  # NB order is important. True is both a bool and an int; 1 is only an int
  typeTags = [
   (bool, 'booleanValue'),
   (int,'intValue'),
   (float,'floatVaue'),
  ]

  for typ,tag in typeTags:
    if isinstance(value, typ):
      setattr(runParameter, tag, value)
      return

  if value:
    setattr(runParameter, 'textValue', str(value))



def addParamToDict(key, val, dictionary):
  """ Adds Runparameter.value to dictionary
  """

  prev = dictionary.get(key)
  if prev is None:
    dictionary[key] = val
  elif isinstance(prev, list):
    prev.append(val)
  else:
   dictionary[key] = [prev, val]


####################################################################
#
#                   Other code
#
####################################################################


def getNmrCalcIdentifier(nmrCalcRun):
  """ get standard NmrCalc.Run identifier string
  """
  return '%s+%s' % (nmrCalcRun.nmrCalcStore.guid, nmrCalcRun.serial)

def runIsWritable(nmrCalcRun):
  """ check that run can be written to

  NBNB TBD should use run.status at some point
  """
  result = True
  if (nmrCalcRun.findAllData(ioRole='output') or
      nmrCalcRun.findAllRunParameters(ioRole='output')):
    result = False
  #
  return result

def getMeasurementList(nmrCalcRun, measurementType='ShiftList', ioRole='input'):
  """Get MeasurementListData of nmrCalcRun.
  Returns None if the measurementList is missing or ambiguous.
  """
  ll = []
  for obj in nmrCalcRun.findAllData(className='MeasurementListData',
                                    ioRole='input'):
    target = obj.measurementList
    if target and target.className == measurementType:
      ll.append(target)
  if len(ll) == 1:
    result = ll.pop()
  else:
    result = None
  #
  return result

def getConstraintList(nmrCalcRun, measurementType='ChemShiftConstraintList',
                      ioRole='input'):
  """Get constraint list of type measurementType
  from ConstraintListData of nmrCalcRun.
  Returns None if the measurementList is missing or ambiguous.
  """
  result = None

  ll = nmrCalcRun.findAllData(className='ConstraintStoreData',
                              ioRole=ioRole)
  if len(ll) == 1:
    obj = ll.pop()
    targets = obj.findAllConstraintLists(className=measurementType)
    if len(targets) == 1:
      result = targets.pop()
  #
  return result

def parseResidueSelection(selection, molSystem):
  """ convert selection string of form 'A7,A9,A11-88,B' to residue selection.
  NB chain code is optional where there is only one chain.
  """

  residues = []

  selection = selection.strip()
  if selection in ('*',''):
    # use all residues - same as empty
    return residues

  # set up
  chains = {}
  ll = molSystem.sortedChains()
  for chain in ll:
    chains[chain.code] = chain

  if len(ll) == 1:
    singleChain = ll[0]
  else:
    singleChain = None

  # get selected residues
  for ss in selection.split(','):
    ss = ss.strip()

    # get correct chain
    for cCode in chains:
      if ss.startswith(cCode):
        chain = chains[cCode]
        ss = ss[len(cCode):]
        break
    else:
      if singleChain is not None:
        chain = singleChain
      else:
        raise ValueError("Unknown chain code %s in input %s"
                         % (cCode,selection))

    # Add residues
    if not ss:
      # we have an entire chain
      residues.extend(chain.sortedResidues())

    else:
      ll = ss.split('-')
      if len(ll) == 1:
        residues.append(chain.findFirstResidue(seqId=int(ll[0])))
      else:
        for ii in range(int(ll[0]),int(ll[1])+1):
          residues.append(chain.findFirstResidue(seqId=ii))
  #
  return residues

def filterMeasurements(measurements, atomNames=None, residues=None,
                       minFigOfMerit=0.0):
  """ Select measurements where all atom names are in atomNames (if given)
  all residues are in residues (if given)
  and figOfMerit is higher than the minimum

  NB Only accepts single, stereospecifically assigned atoms
  """

  result = []

  if not measurements:
    return result

  xx = measurements[0]
  if hasattr(xx, 'resonances'):
    nAssign = len(xx.resonances)
  else:
    nAssign = 1

  for meas in measurements:

    if meas.figOfMerit >= minFigOfMerit:

      if nAssign == 1:
        resonances = (meas.resonance,)
      else:
        resonances = meas.resonances

      for resonance in resonances:
        resonanceSet = resonance.resonanceSet

        if resonanceSet is not None and len(resonanceSet.resonances) == 1:

          atoms = []
          for atomSet in resonanceSet.atomSets:
            atoms.extend(atomSet.atoms)

          for atom in atoms:
            if atomNames and atom.name not in atomNames:
              break
            if residues and atom.residue not in residues:
              break
          else:
            # Everything OK so far.
            # NB prochirals and methyl groups need to be dealt with here
            if len(atoms) != 1:
              # single assignments only, for now
              break

        else:
          break

      else:
        # it all worked out. Add measurement to result
        result.append(meas)
  #
  return result

def makeSingleConstraintList(constraintStore, measurements,
                             constraintType='ChemShift', name=None,
                             details=None):
  """ Make new one-resonance ConstraintList containing measurements in
  constraintStore
  """

  createFunction = getattr(constraintStore,'new%sConstraintList'
                                           % constraintType)
  newList = createFunction(name=name, details=details)

  measurementLists = set()

  for meas in measurements:
    measurementLists.add(meas.parentList)
    fixedResonance = ConstraintBasic.getFixedResonance(constraintStore,
                                                       meas.resonance)
    createFunction = getattr(newList,'new%sConstraint' % constraintType)
    createFunction(targetValue=meas.value, error=meas.error,
                   resonance=fixedResonance)

  newList.measurementLists = measurementLists
  #
  return newList

def makePairConstraintList(constraintStore, measurements,
                           constraintType='JCoupling', name=None,
                           details=None):
  """ Make new two-reasonance ConstraintList containing measurements in
      constraintStore
  """

  createFunction = getattr(constraintStore,'new%sConstraintList'
                                           % constraintType)
  newList = createFunction(name=name, details=details)

  measurementLists = set()

  for meas in measurements:
    measurementLists.add(meas.parentList)
    fixedResonances = []
    for resonance in meas.resonances:
      fixedResonances.append(ConstraintBasic.getFixedResonance(constraintStore,
                                                              resonance))
    createFunction = getattr(newList, 'new%sConstraint' % constraintType)
    constr = createFunction(targetValue=meas.value, error=meas.error)
    createFunction = getattr(constr, 'new%sConstraintItem' % constraintType)
    createFunction(resonances=fixedResonances)


  newList.measurementLists = measurementLists
  #
  return newList


def constraintListFromData(nmrCalcRun, constraintStore, data, name, code,
                           numResonances, parameterGroup):
  """ NBNB TBD
  """
  return None


def measurementListFromConstraints(nmrCalcRun, name, constraintLists):
  """ NBNB TBD
  """
  return None

def getNmrCalcStore(memopsRoot, nmrProjectName=None, nmrCalcStoreName=None):
  """
  """

  print '### Util getNmrCalcStore', nmrProjectName, nmrCalcStoreName

  result = None

  if nmrProjectName is None:
    nmrProject = memopsRoot.currentNmrProject
  else:
    nmrProject = memopsRoot.findFirstNmrProject(name=nmrProjectName)
  if nmrProject is None:
    raise Exception("No NmrProject with name %s found" % nmrProjectName)

  nmrProjectName = nmrProjectName or nmrProject.name

  if nmrCalcStoreName is None:
    # Take last in sorted NmrCalcStores. NB, sort key is nmrProjectName,name
    for nmrCalcStore in reversed(memopsRoot.sortedNmrCalcStores()):
      if nmrCalcStore.nmrProjectName == nmrProjectName:
        result = nmrCalcStore
        break

  else:
    # Find NmrCalcStore that fits both names.
    result = memopsRoot.findFirstNmrCalcStore(name=nmrCalcStoreName,
                                              nmrProjectName=nmrProjectName)
  if result is None:
    # Nothing found - make a new one
    print '### ALARM making newNmrCalcStore', nmrProjectName, nmrCalcStoreName or 'CCPNauto'
    result = memopsRoot.newNmrCalcStore(nmrProjectName=nmrProjectName,
                                        name=nmrCalcStoreName or 'CCPNauto')
  #
  return result


def makeNmrCalc(memopsRoot, jsonObject):
  """Make NmrCalc instance within memopsRoot, based on info in jsonObj
  Used to make task starting data for WMS.

  Returns nmrCalcId
  """
  #
  # Get or make NmrCalcStore to use
  #
  # First get names
  nmrCalcStoreName = jsonObject['CCPN.NmrCalcStore.name']
  nmrProjectName = jsonObject.get('CCPN.NmrProject.name')

  nmrCalcStore = getNmrCalcStore(memopsRoot, nmrProjectName, nmrCalcStoreName)

  #
  # Now start entering data
  #

  # First make a Run
  wmsProtocolName = jsonObject.get('CCPN.Run.wmsProtocolName')
  run = nmrCalcStore.newRun(wmsProtocolName=wmsProtocolName,
                            details=jsonObject.get('details'))

  nmrCalcId = nmrCalcStore.guid + '+' + str(run.serial)

  # Make Top-level parameters
  for tag in jsonObject:
    if (tag not in classGenerators and tag not in nonDataTags
        and not tag.startswith('CCPN.')):
      param = run.newRunParameter(name=tag)
      setParameterValue(param, jsonObject[tag])

  # Make Data and associated parameters
  try:
    run.root.override = True
    print '### Override On', run.root

    for classTag in sorted(runDataTags):
      tags = runDataTags[classTag]
      if tags:
        generator = getattr(run,classGenerators[classTag])
        for dd in jsonObject.get(classTag, ()):
          print '### Data', run.root.override, run.root, [(x,dd[x]) for x in tags if x in dd]
          parent = generator(**dict((x,dd[x]) for x in tags if x in dd))
          parent.checkValid()
          for tag in dd:
            if tag not in tags:
              # Associated parameter.
              param = run.newRunParameter(code=tag, data=parent)
              setParameterValue(param, dd[tag])
              param.checkValid()

  except:
    raise

  finally:
    print '### Override Off'
    run.root.override = False

  # Make RunParameters
  tags = nonDataTags['RunParameter']
  for dd in jsonObject.get('RunParameter', ()):
    print '### RunParameter', [(x,dd[x]) for x in tags if x in dd]
    par = run.newRunParameter(**dict((x,dd[x]) for x in tags if x in dd))
    dataSerial = dd.get('dataSerial')
    if dataSerial:
      par.data = run.findFirstData(serial=dataSerial)
  '''
  # MolSystemData
  tags = ('name', 'code', 'ioRole', 'details',
          'molSystemCode', 'symmetrySetId','chainCodes')
  for dd in jsonObject.get('MolSystem', ()):
    parent = run.newMolSystemData(**dict((x,dd[x]) for x in tags if x in dd))

  # MolResidueData
  tags = ('name', 'code', 'ioRole', 'details',
          'molSystemCode','chainCode', 'residueSeqIds', 'chainCodes')
  # NBNB You must set either chainCode or chainCodes. See API docs for rules.
  # March 2012 WMS can handle only single-chain sets, so only chainCode
  for dd in jsonObject.get('MolResidue', ()):
    run.newMolResidueData(**dict((x,dd[x]) for x in tags if x in dd))

  # ConstraintStoreData
  tags = ('name', 'code', 'ioRole', 'details',
          'constraintStoreSerial', 'constraintListSerials')
  for dd in jsonObject.get('ConstraintStore', ()):
    run.newConstraintStoreData(**dict((x,dd[x]) for x in tags if x in dd))

  # ViolationListData
  tags = ('name', 'code', 'ioRole', 'details',
          'constraintStoreSerial', 'violationListSerial')
  for dd in jsonObject.get('ViolationList', ()):
    run.newViolationListData(**dict((x,dd[x]) for x in tags if x in dd))

  # PeakListData
  tags = ('name', 'code', 'ioRole', 'details',
          'experimentSerial', 'dataSourceSerial', 'peakListSerial')
  for dd in jsonObject.get('PeakList', ()):
    run.newPeakListData(**dict((x,dd[x]) for x in tags if x in dd))

  # SpectrumData
  tags = ('name', 'code', 'ioRole', 'details',
          'experimentSerial', 'dataSourceSerial')
  for dd in jsonObject.get('Spectrum', ()):
    run.newSpectrumData(**dict((x,dd[x]) for x in tags if x in dd))

  # MeasurementListData
  tags = ('name', 'code', 'ioRole', 'details',
          'measurementListSerial')
  for dd in jsonObject.get('MeasurementList', ()):
    run.newMeasurementListData(**dict((x,dd[x]) for x in tags if x in dd))

  # DerivedListListData
  tags = ('name', 'code', 'ioRole', 'details',
          'derivedDataListSerial')
  for dd in jsonObject.get('DerivedList', ()):
    run.newDerivedListData(**dict((x,dd[x]) for x in tags if x in dd))

  # StructureEnsembleData
  tags = ('name', 'code', 'ioRole', 'details',
          'molSystemCode', 'ensembleId','modelSerials')
  for dd in jsonObject.get('StructureEnsemble', ()):
    run.newStructureEnsembleData(**dict((x,dd[x]) for x in tags if x in dd))

  # ExternalData
  tags = ('name', 'code', 'ioRole', 'details',
          'dataLocationStoreName', 'dataStoreSerial')
  for dd in jsonObject.get('External', ()):
    run.newDerivedListData(**dict((x,dd[x]) for x in tags if x in dd))

  # FloatMatrixData
  tags = ('name', 'code', 'ioRole', 'details',
          'shape', 'data', 'defaultValue')
  for dd in jsonObject.get('FloatMatrix', ()):
    run.newFloatMatrixData(**dict((x,dd[x]) for x in tags if x in dd))

  # SpinSystemData
  # NBNB TBD not implemented. Later, at need, when model is clearer

  # EnergyTerm
  # NBNB Skipped. Later, if needed

  # TensorData
  # NBNB Skipped. Later, if needed

  # ParameterGroup
  # NBNB Skipped. Later, if needed
  '''

  #
  return nmrCalcId



def wmsProjectSummary(memopsRoot, nmrProjectName=None):
  """Get JSON-compatible project summary for use in object selection
  for programs using NmrCalc. Specially designed for WMS. Includes attributes
  needed for making NmrCalc objects, with additional attributes for annotation.
  """

  result = {}

  # get NmrProjects to use
  if nmrProjectName is None:
    nmrProjects = memopsRoot.sortedNmrProjects()
  else:
    obj = memopsRoot.findFirstNmrProject(name=nmrProjectName)
    if obj is None:
      raise KeyError("No ccp.nmr.NmrProject named %s" % nmrProjectName)
    else:
      nmrProjects = [obj]

  # get data that depend on NmrProject(s)
  for nmrProject in nmrProjects:

    # MeasurementList
    for measurementList in nmrProject.sortedMeasurementLists():
      objlist = result.setdefault('MeasurementList', [])
      objlist.append({
       'measurementListSerial' : measurementList.serial,
       'measurementListName' : measurementList.name,
       'details' : measurementList.details,
       'className' : measurementList.className,
      })

    # DerivedDataList
    for derivedDataList in nmrProject.sortedDerivedDataLists():
      objlist = result.setdefault('DerivedList', [])
      objlist.append({
       'measurementListSerial' : derivedDataList.serial,
       'details' : derivedDataList.details,
       'className' : derivedDataList.className,
      })

    # Spectrum and PeakList
    for experiment in nmrProject.sortedExperiments():
      experimentSerial = experiment.serial

      # Spectrum entry
      for dataSource in experiment.sortedDataSources():
        dataSourceSerial = dataSource.serial
        speclist = result.setdefault('Spectrum', [])
        dd = {
         'experimentSerial':experimentSerial,
         'dataSourceSerial':dataSourceSerial,
         'experimentName':experiment.name,
         'dataSourceName':dataSource.name,
         'details':dataSource.details
        }
        speclist.append(dd)

        # PeakList entry
        for peakList in dataSource.sortedPeakLists():
          dd2 = dd.copy()
          objlist = result.setdefault('PeakList', [])
          dd2['details'] = peakList.details
          dd2['peakListSerial'] = peakList.serial
          objlist.append(dd2)

    # ConstraintStore and ViolationList
    for nmrConstraintStore in nmrProject.sortedNmrConstraintStores():
      constraintStoreSerial = nmrConstraintStore.serial

      # ConstraintStoreData
      serials = []
      names = []
      types = []
      for constraintList in nmrConstraintStore.sortedConstraintLists():
        serials.append(constraintList.serial)
        names.append(constraintList.name)
        types.append(constraintList.className)
      if serials:
        objlist = result.setdefault('ConstraintStore', [])
        objlist.append({
         'constraintStoreSerial':constraintStoreSerial,
         'constraintListSerials':serials,
         'constraintListNames':names,
         'constraintTypeNames':types,
        })

      # ViolationList
      for violationList in nmrConstraintStore.sortedViolationLists():
        objlist = result.setdefault('ViolationList', [])
        objlist.append({
         'constraintStoreSerial':constraintStoreSerial,
         'violationListSerial':violationList.serial,
         'details':violationList.details
        })

  # get data that depend on MolSystem
  for molSystem in memopsRoot.sortedMolSystems():
    molSystemCode = molSystem.code

    # MolSystem
    objlist = result.setdefault('MolSystem', [])
    for symSet in [None] + molSystem.sortedMolSystemSymmetrySets():
      dd = {'molSystemCode':molSystemCode,
            'chainCodes':[x.code for x in molSystem.sortedChains()]
      }
      objlist.append(dd)
      if symSet:
        dd['symmetrySetId'] = symSet.symmetrySetId


    # StructureEnsemble
    for structureEnsemble in molSystem.sortedStructureEnsembles():
      objlist = result.setdefault('StructureEnsemble', [])
      objlist.append({
       'molSystemCode':molSystemCode,
       'ensembleId':structureEnsemble.ensembleId,
       'modelSerials':[x.serial for x in structureEnsemble.sortedModels()],
      })

    # MolResidue
    # NB only single-chain sets supported for now
    # The model allows multiple-chain sets.
    for chain in molSystem.sortedChains():
      objlist = result.setdefault('MolResidue', [])
      residues = chain.sortedResidues()
      objlist.append({
       'molSystemCode':molSystemCode,
       'chainCode':chain.code,
       'residueSeqIds':[x.seqId for x in residues],
       'residueCompIds':[x.ccpCode for x in residues],
      })
  #
  return result

def integerListExpression(integers, fieldSep='-', termSep = ','):
  """ Compress list of integers into string of format
  "2-11,15-27,29,31,35-99"
  """

  if not integers:
    return ''

  format = '%s' + fieldSep + '%s'

  integers = list(sorted(integers))
  last = integers[0]
  edges = [last]
  for next in integers[1:]:
    if next - last != 1:
      edges.append(last)
      edges.append(next)
    last = next

  edges.append(integers[-1])

  ll = []
  for ii in range(0,len(edges),2):
    ll.append(format % (edges[ii], edges[ii+1]))
  #
  return termSep.join(ll)


def adaptNmrCalcRun(nmrCalcRun):
  """ Add mandatory, program-independent adaptions fo NmrCalcRun
      nmrCalcRun: NmrCalc.Run
  """

  print 'Starting Util.adaptNmrCalcRun'

  project = nmrCalcRun.root
  nmrProject = nmrCalcRun.nmrCalcStore.nmrProject

  xx = nmrCalcRun.findFirstRunParameter(name='shiftFormat')
  if xx is None:
    shiftFormat = None
  else:
    shiftFormat = xx.textValue
    print '### shiftFormat', shiftFormat
    shiftFileSuffix = getFileFormatData(shiftFormat, 'shift')['shiftExt']
  xx = nmrCalcRun.findFirstRunParameter(name='peakFormat')
  print xx
  if xx is None:
    peakFormat = None
  else:
    peakFormat = xx.textValue
    print '### peakFormat', peakFormat
    peakFileSuffix = getFileFormatData(peakFormat, 'peak')['peakExt']
  xx = nmrCalcRun.findFirstRunParameter(name='restraintFormat')
  if xx is None:
    restraintFormat = None
  else:
    restraintFormat = xx.textValue
    print '### restraintFormat', restraintFormat

  # Get ConstraintStores
  constraintStoreSerials = set(x.constraintStoreSerial
                                 for x in nmrCalcRun.data
                                 if x.className in ('ConstraintStoreData',
                                                    'ViolationListData'))
  if constraintStoreSerials:
    # useSerial = max(constraintStoreSerials)
    mainConstraintStore = project.findFirstNmrConstraintStore(serial=max(constraintStoreSerials))
    # make sure FixedResonances are mapped to resonances(set FixedResonance.resonanceSerial)
    linkFixedResonancesToResonances(mainConstraintStore)
    print '### validating for constraintStore: ', mainConstraintStore.serial
    validateFixedAssignments(mainConstraintStore, nmrProject)

  else:
    mainConstraintStore = project.newNmrConstraintStore(
                                   nmrProject=nmrProject)
  assignStoreData = nmrCalcRun.newConstraintStoreData(
    name='assignment',nmrConstraintStore=mainConstraintStore,
    details='Input assignment state and shifts (as constraints) ' +
            'set by CcpNmr Integrator')


  # copy assignment status to main constraint store
  for resonance in nmrProject.sortedResonances():
    ConstraintBasic.getFixedResonance(mainConstraintStore, resonance)

  # get shift list objects
  shiftListObjs = [x for x in nmrCalcRun.sortedData()
                  if x.className == 'MeasurementListData'
                  and x.measurementList.className == 'ShiftList']
  shiftLists = set(x.measurementList for x in shiftListObjs)

  # set up peak files and add corresponding shift lists
  peakListObjs = [x for x in nmrCalcRun.sortedData()
                  if x.className == 'PeakListData']


  # make file names for constraintLists
  # get constraintList data object
  for xx in nmrCalcRun.sortedData():
    if (xx.className == 'ConstraintStoreData' and xx.ioRole == 'input'
        and xx.name == 'constraintLists'):
      constraintListData = xx
      break
  else:
    constraintListData = None

  if constraintListData and restraintFormat:
    # NBNB file extensions are added later, since they depend on both list and program type.
    for constraintList in constraintListData.constraintLists:
      constraintTag = constraintList.className[:-14]
      constraintTag = constraintTag[0].lower() + constraintTag[1:]
      fileName = '%s_%s.%s' % (constraintTag, constraintList.serial,
                               constraintFileExtension(constraintList, restraintFormat))
      nmrCalcRun.newRunParameter(name='fileName', textValue=fileName,
                                 code='fileName', data=constraintListData)

  # peakFiles = []
  for peakListData in peakListObjs:
    if peakFormat is not None:
      nmrCalcRun.newRunParameter(name='fileFormat', textValue=peakFormat,
                                 data=peakListData)
      fileNamePeakList = objectFileName(peakListData.peakList,
                                        suffix=peakFileSuffix, compact=True)
      nmrCalcRun.newRunParameter(name='fileName', textValue=fileNamePeakList,
                                 code='fileName', data=peakListData)


    if shiftFormat is not None:
      shiftList = peakListData.dataSource.experiment.shiftList
      if shiftList not in shiftLists:
        shiftLists.add(shiftList)
        shiftListData = nmrCalcRun.newMeasurementListData(name='shiftList',
                          measurementList=shiftList)
        shiftListObjs.append(shiftListData)

  # If there are still no shiftLists, take the first non-empty one - we want at least one
  if not shiftListObjs:
    for shiftList in nmrCalcRun.nmrCalcStore.nmrProject.sortedMeasurementLists():
      if shiftList.measurements and shiftList.className == 'ShiftList':
        shiftLists.add(shiftList)
        shiftListObjs.append(nmrCalcRun.newMeasurementListData(name='shiftList',
                             measurementList=shiftList))
        break
    else:
      print('WARNING, no shift lists found')


  if shiftFormat is not None:
    for shiftListData in shiftListObjs:
      nmrCalcRun.newRunParameter(name='fileFormat', textValue=shiftFormat,
                                 data=shiftListData)
      shiftList = shiftListData.measurementList
      shiftFile = objectFileName(shiftList, suffix=shiftFileSuffix,
                                 compact=True)
      nmrCalcRun.newRunParameter(name='fileName', textValue=shiftFile,
                                 code='fileName', data=shiftListData)


  # copy shiftlists to main constraintstore
  resMapDict = mainConstraintStore.quickResonances
  for shiftList in shiftLists:
    #ll = [x for x in mainConstraintStore.constraintLists
    ll = [x for x in mainConstraintStore.findAllConstraintLists(
                                         className='ChemShiftConstraintList')
          if shiftList in x.measurementLists]
    if not ll:
      # make ChemShiftConstraintList
      shiftConstraintList = mainConstraintStore.newChemShiftConstraintList(
        measurementLists=[shiftList], name=shiftList.name, unit=shiftList.unit)
      assignStoreData.addConstraintListSerial(shiftConstraintList.serial)
      for shift in shiftList.sortedMeasurements():
        shiftConstraintList.newChemShiftConstraint(targetValue=shift.value,
          error=shift.error, resonance=resMapDict[shift.resonance.serial])


def validateFixedAssignments(constraintStore, nmrProject):
  """ Check that ConstraintStore assignments are compatible with
  nmrProject asssignments
  """
  if not nmrProject is constraintStore.nmrProject:
    raise Exception("%s is not connected to %s)" % (constraintStore,nmrProject))

  # set up mapping dict
  resonances = {}
  for res in nmrProject.resonances:
     resonances[res.serial] = res

  # check atom compatibility
  unAssigned = {}
  for fres in constraintStore.fixedResonances:
    fresSet = fres.resonanceSet
    fatoms = None
    if fresSet:
      fatoms = set(y for x in fresSet.atomSets for y in x.atoms)
      atoms = None
      if fatoms:
        #unAssigned[fres.serial] = fres
        #res = resonances.get(fres.serial)
        # Weird error. Surely we wannot assume that serials are the same?
        #     resonanceSerial is the correct attribute
        res = resonances.get(fres.resonanceSerial)
        print(res)
        if res:
          resSet = res.resonanceSet
          if resSet:
            atoms = set(y for x in resSet.atomSets for y in x.atoms)
        if atoms:
          if atoms != fatoms:
            raise Exception ('WARNING, incompatibility, assignment %s for %s differs from %s'
                 % (fres.name, fres, res.name))
        else:
          print ('WARNING, incompatibility, extra assignment %s for %s'
                 % (fres.name, fres))
    if not fatoms:
      unAssigned[fres.resonanceSerial] = fres

  # transfer additional assignments
  for res in nmrProject.resonances:
    resSet = res.resonanceSet
    if resSet:
      fres = unAssigned.get(res.serial)
      if fres is not None:
        # transfer assignment
        ConstraintBasic.updateFixedResonance(fres,res)


def objectFileName(ccpnObj, suffix=None, compact=False):
  """Make file name that matches (and identifies CCPN object)
  """

  ll = [str(x) for x in ccpnObj.getExpandedKey()]
  if compact:
    ll[0] = ccpnObj.className
  else:
    ll = [ccpnObj.className] + ll
  result = '_'.join(ll)
  if suffix:
    result = '.'.join((result,suffix))
  #
  return result

def getAmalgamatedTolerances(dataSources, defPoints, minTol,
                             tagOrder=('hx2', 'hx1', 'x1', 'x2')):
  """ Get maximum tolerance over a series of datasources, by dimCode
  Assumes that (h)x1 is the bound dimension if there is only one,
  the acquisition dimension otherwise
  """
  tolerances = []
  dimCodes = []

  for dataSource in dataSources:
    tolerances.append(AssignmentBasic.estimateAssignmentTolerances(dataSource,
                                           defPoints=defPoints, minTol=minTol))
    dimCodes.append(getDimCodes(dataSource))

  # Get tolerance parameter
  tags = ['hx1', 'x1', 'x2', 'hx2']
  dd = {}
  for tag in tags:
    dd[tag] = []

  # match tolerances to bound/free dimensions
  for ii , codes in enumerate(dimCodes):
    tols = tolerances[ii]
    ll = []
    for tag in tags:
      if tag in codes:
        ll.append(tols[codes.index(tag)])
      else:
        ll.append(None)
    if 'x2' in codes and not 'x1' in codes:
      # reverse if hx2 is bound dimension
      ll.reverse()
    for ii,tag in enumerate(tags):
      val = ll[ii]
      if val:
        dd[tag].append(val)
  # get max values and put in order
  result = []
  for tag in ['hx2', 'hx1', 'x1', 'x2']:
    ll = dd[tag]
    if ll:
      result.append(max(ll))
    else:
      break
  #
  return result


def getResonanceAtomMap (namingSystemName, nmrConstraintStore, shiftListObj=None, **keywds):
  """Get map from resonance AND fixedresonances to ResonanceAtom objects,
  compressing atom names depending on parameters"""

  fixedResonances = nmrConstraintStore.sortedFixedResonances()
  resonanceToAtoms = DataFormat.getResonanceAtomMap(namingSystemName,
                                                    fixedResonances, **keywds)

  if shiftListObj and not keywds.get('individualAtoms'):
    compressAssignments(resonanceToAtoms, shiftListObj)

  # Convert from fixed to real resonances

  # make fixedres - res map
  fresMap = {}
  nmrProject = nmrConstraintStore.nmrProject
  for fres,val in resonanceToAtoms.items():
    res = nmrProject.findFirstResonance(serial=fres.resonanceSerial)
    if res is None:
      vv = val[0]
      ss =   '%s.%s.%s - %s' % (vv.chain.code, vv.seqId, vv.atomName or '',
                                vv.atomSetName or '')
      print ("WARNING, no resonance match for FixedResonance %s, (%s) - Skipping"
             % (fres.serial, ss))
    else:
      fresMap[fres] = res

  resonanceToAtoms = {}
  items = list(resonanceToAtoms.items())
  # NB we must make seaprate list as we are modifying teh dictionary during the loop
  for fres, val in items:
    res = fresMap.get(fres)
    if res is not None:
      resonanceToAtoms[res] = ll = []
      for r2a in val:
        newR2a = r2a.clone()
        ll.append(newR2a)
        newR2a.resonance = res
        for tag in ('otherLinkedResonances', 'otherGroupResonances'):
          mm = [fresMap[x] for x in getattr(r2a,tag) if x in fresMap]
          setattr(newR2a, tag, mm)


def compressAssignments(resToAtoms, shiftValueList=None, tolerance=None):
  """ Modify resToAtoms in place
  so that resonanceSets with degenerate or unknown shifts keep their atomSets
  and others have them removed.
  shiftValueList can be either a ChemShiftConstraintList or ShiftList
  """

  if not tolerance:
    tolerance = genConstants.shiftIdentityTolerance

  # set up
  reverseMapping = {}

  # input check
  if shiftValueList.className == 'ChemShiftConstraintList':
    valTag = 'targetValue'
    findFunc = shiftValueList.findFirstConstraint
    for res in resToAtoms:
      if res.className !='FixedResonance':
        raise Exception("ChemShiftConstraintList requires FixedResonance mapping")

  elif shiftValueList.className == 'ShiftList':
    valTag = 'value'
    findFunc = shiftValueList.findFirstMeasurement
    for res in resToAtoms:
      if res.className !='Resonance':
        raise Exception("ShiftList requires FixedResonance mapping")

  else:
    shiftValueList = None
    print 'Warning, compressing shift resonance mapping without shift input:'


  # make reverse mapping
  #dbg=[]
  for res in sorted(resToAtoms, key=operator.attrgetter('serial')):
    r2as = resToAtoms[res]
    for r2a in r2as:
      atomSetName = r2a.atomSetName
      if atomSetName is not None:
        bigChemAtomSet = r2a.getChemAtomSet()
        if bigChemAtomSet:
          key = (r2a.chain.code, r2a.seqId, bigChemAtomSet.name)
          #dbg.append (('\n### PRELIM ', key, r2a.atomName, r2a.resonance.serial,
          #       atomSetName, bigChemAtomSet != r2a.chemAtomSet,
          #       bigChemAtomSet and bigChemAtomSet.name,
          #       r2a.chemAtomSet and r2a.chemAtomSet.name))
          if bigChemAtomSet is not r2a.chemAtomSet:
            # This atom/atomSet uis part of a bigger atomSet
            ll = reverseMapping.get(key)
            if ll is None:
              reverseMapping[key] = [r2a]
            else:
              ll.append(r2a)


  #for tt in sorted(dbg):
  #  print tt

  # loop over reverse maping and modify
  for tt, r2as in sorted(reverseMapping.items()):

    if len(r2as) > 1:

      compress = True

      if shiftValueList is not None:
        objs = [findFunc(resonance=r2a.resonance) for r2a in r2as]
        objs = [x for x in objs if x is not None]
        shifts = [getattr(obj,valTag) for obj in objs]
        shifts = [x for x in shifts if x is not None]
        if len(shifts) > 1:
          diff = max(abs(shift - shifts[0]) for shift in shifts[1:])
          error = max(obj.error or 0.0 for obj in objs)
          if diff > error and diff > tolerance:
            # There are at least two significantly different shifts
            # Do not compress
            compress = False
      #
      if compress:
        #print '\n### COMPR ', tt, r2a.resonance.serial
        for r2a in r2as:
          r2a.useAtomSetName = True
          #print ('---C', r2a.atomName, r2a.atomSetName,
          #              r2a.chemAtom and  r2a.chemAtom.name,
          #              r2a.chemAtomSet and  r2a.chemAtomSet.name,
          #              r2a.useAtomSetName,
          #              [x.serial for x in r2a.otherLinkedResonances],
          #              [x.serial for x in r2a.otherGroupResonances])


      #else:
      #  print '\n### NOCOMP', tt, r2a.resonance.serial
      #  for r2a in r2as:
      #    print ('---N', r2a.atomName, r2a.atomSetName,
      #                  r2a.chemAtom and  r2a.chemAtom.name,
      #                  r2a.chemAtomSet and  r2a.chemAtomSet.name,
      #                  r2a.useAtomSetName,
      #                  [x.serial for x in r2a.otherLinkedResonances],
      #                  [x.serial for x in r2a.otherGroupResonances])
      #    #r2a.atomSetName = None


def setRunParametersFromConfig(nmrCalcRun, confDict, ioRole='input'):
  """ Set or reset RunParameters in nmrCalcRun from confDict.
  confDict can be a jsonObject or a nested dictionary structure
  that matches the protocol description format.
  """

  for dd in confDict['protocol']['protocolParameters']:
    paramType = dd.get('paramType')
    if paramType and 'value' in dd and not 'relatedParameter' in dd:
      name = dd.get('name')
      code = dd.get('code')
      runPar = None
      if name:
        runPar = nmrCalcRun.findFirstRunParameter(name=name, data=None,
                                                  ioRole=ioRole)
      elif code:
        runPar = nmrCalcRun.findFirstRunParameter(code=code, data=None,
                                                  ioRole=ioRole)
      else:
        # Neither name nor code set. Skip.
        continue

      if runPar is None:
        runPar = nmrCalcRun.newRunParameter(name=name, ioRole=ioRole)

      runPar.code = code

      print '### runParam ', name, code, dd['value']

      valueAttr = paramTypeMap[paramType]
      setattr(runPar, valueAttr, dd['value'])

def makeDerivedRun(masterRun):
  """ Make derivedRun as copy from masterRun
  """

  #NB 15 Aug 2013 Rasmus Fogh - Added provisional status and deletion of output

  ll = masterRun.nmrCalcStore.sortedRuns()
  nextSerial = ll[-1].serial + 1
  nmrCalcRun = genUtil.copySubTree(masterRun, masterRun.nmrCalcStore,
                                   maySkipCrosslinks=True,
                                   topObjectParameters={'serial':nextSerial,
                                                        'status':'provisional'})
  # Remove output data
  for ll in (nmrCalcRun.sortedData(), nmrCalcRun.sortedRunParameters()):
    for obj in ll:
      if obj.ioRole == 'output':
        obj.delete()

  # NBNB TODO consider changing this later
  # For now, avoid nested master run (not allowed by model
  while masterRun.masterRun is not None:
    masterRun = masterRun.masterRun

  nmrCalcRun.masterRun = masterRun

  #
  return nmrCalcRun


def constraintFileExtension(constraintList, fileFormat=None, lowerLimit=False):
  """Get extension for constraint file.
  code in seperate function as the problem is ocmplex for definition dictionaries"""

  if 'CYANA' in fileFormat.upper():
    if lowerLimit:
      return 'lol'
    elif constraintList.className == "DihedralConstraintList":
      return 'aco'
    else:
      return 'upl'

  else:
    # Default extension - used whenever nothing specific is known
    return 'tbl'

def linkFixedResonancesToResonances(constraintStore):
  """Match assigned FixedResonances to assigned Resonances
  based on shared atoms assignment and set resonanceSerial to link them"""
  nmrProject = constraintStore.nmrProject
  # atoms2Resonance = []
  atoms2ResonanceSets = []
  for resonanceSets in (constraintStore.sortedFixedResonanceSets(),
                        nmrProject.sortedResonanceSets()):
    # a2Res = {}
    # atoms2Resonance.append(a2Res)
    a2ResSet = {}
    atoms2ResonanceSets.append(a2ResSet)
    for resonanceSet in resonanceSets:
      atoms = frozenset(x for y in resonanceSet.atomSets for x in y.atoms)
      ll = a2ResSet.get(atoms)
      if ll is None:
        a2ResSet[atoms] = [resonanceSet]
      else:
        ll.append(resonanceSet)


  # print('### TOTALS resSets %s; fresSets %s'
  #       % (len(nmrProject.resonanceSets), len(constraintStore.fixedResonanceSets)))

  for atoms, fresonanceSets in atoms2ResonanceSets[0].items():
    resonanceSets = atoms2ResonanceSets[1].get(atoms, ())
    # ll = [sorted(x.serial for y in resonanceSets for x in y.resonances)]
    # fll = [sorted(x.serial for y in fresonanceSets for x in y.resonances)]
    # print("### res %s; fres %s, EQSERIALS %s"
    #       % (ll, fll, ll==fll))

    if len(resonanceSets) == len(fresonanceSets):
      if len(resonanceSets) == 1:
        resonances = list(resonanceSets)[0].sortedResonances()
        fresonances = list(fresonanceSets)[0].sortedResonances()
        fresonance = fresonances[0]
        if len(resonances) == len(fresonances):
          if len(resonances) == 1:
            if fresonance.resonanceSerial is None:
              # print("### FIXED - atoms %s, sets:1-1" %
              #       [x.name for x in atoms])
              fresonance.resonanceSerial = resonances[0].serial
            # elif fresonance.resonanceSerial == resonances[0].serial:
            #   print("### WASOK - atoms %s, sets:1-1" %
            #         [x.name for x in atoms])
            # else:
            #   print("### ASSMISMATCH - atoms %s, sets:1-1" %
            #        [x.name for x in atoms])
    #       else:
    #         print('### MULTIRES OK:%s fres %s fatomset %s res %s atomset %s' %
    #               ([x.name for x in fresonances], not [x for x in fresonances if x.resonanceSerial is None],
    #               [x.name for y in fresonanceSets[0].sortedAtomSets() for x in y.sortedAtoms()],
    #               [x.name for x in resonances],
    #               [x.name for y in resonanceSets[0].sortedAtomSets() for x in y.sortedAtoms()],
    #               ))
    #
    #     else:
    #       print("### RESMISMATCH - atoms %s, res:%s-%s; frnames: %s, names:%s" %
    #             ([x.name for x in atoms], len(resonances), len(fresonances),
    #              [x.name for x in fresonances],
    #              [x.name for x in resonances]
    #             ))
    #   else:
    #
    #     print("### MULTISET - atoms %s, set:%s-%s; " %
    #           ([x.name for x in atoms], len(resonanceSets), len(fresonanceSets)))
    #
    # else:
    #   print("### SETMISMATCH - residue %s atoms %s, set:%s-%s; OKsets %s " %
    #         (list(atoms)[0].residue.seqId, [x.name for x in atoms], len(resonanceSets), len(fresonanceSets),
    #         [not[x for x in y.resonances if x.resonanceSerial is None] for y in fresonanceSets]))


