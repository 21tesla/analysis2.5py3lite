"""Special code for backwards compatibility when reading old release XML

NBNB
NBNB Entire file OBSOLETE after start of version 2.0!!!
NBNB

A converter is defined for each old release supported.
This converter will be passed a project corresponding to an old release,
and will return a project updated to correspond to a new release.

NB when making new compatibility code, the version should be set to 
the first one that is incompatible with older versions,
i.e. the one that has been modified. 
NBNB this was not always observed for the older versions. There are older
upgrades where the version was wset to the last functioning one.
In principle this should be fixed, but the old compatibility code will soon
be history anyway.
The cutoffVersion hack takes care fo the problem meanwhile.
"""

# set up

import types
import os

from memops.general import Constants as genConstants
basicDataTypes = genConstants.baseDataTypeModule


from memops.general.Implementation import ApiError

from memops.general import Version

from memops.universal.Util import semideepcopy

mappingMakers = {}
dataConverters = {}

nonConvertible = 'nonConvertible'

cutoffVersion = Version.getVersion('1.0.207')

#############################################################################
#
# utility functions


def bool2str(value):
  if value:
    return 'true'
  else:
    return 'false'

def str2bool(value):
  if value in ('True', 'true', '1'):
    return True
  elif value in ('False', 'false', '0'):
    return False
  else:
    raise ApiError("String '%s' is not legal for a Boolean" % value)

class DummyClass(object):
  """ dummy class to store temporarily attributes of classes that
  are removed (not just changed).
  
  NB Assumes that class.storage == parent.storage
  """
  def __init__(self, parent, **attrlinks):
    
    # set up
    dataDict = self.__dict__
    dataDict['parent'] = parent
    
    storage = parent.__dict__['storage']
    dataDict['storage'] = storage
    
    # prepare list in parent, and add self to it
    tag = self.__class__.siblingValues
    ll = parent.__dict__.get(tag)
    if ll is None:
      ll = parent.__dict__[tag] = []
    
    ll.append(self)
    
    
    # set attributes
    for key,value in attrlinks.items():
      setattr(self, key, value)
    

#############################################################################
#
# general functions

def getMapping(packageName, fromVersion, currentMapping):
  """ get XML I/O mapping for 'fromVersionString' files into current API
  
  packageName is the package name from 'version' (packages may change names)
  fromVersionString describes the version of the file being read
  currentMapping is the current mapping for the package
  """
  
  return process (
   packageName, fromVersion, currentMapping, mappingMakers, reverse=True
  )
  
  

def postProcess(packageName, fromVersion, headObject):
  """  postprocess data imported from XML  to current API.
  
  packageName is the package name from 'version' (packages may change names)
  fromVersionString describes the version of the file being read
  headObject is the head object for the file being read
  """
  
  return process (packageName, fromVersion, headObject, dataConverters)


def process(packageName, fromVersion, target, infoDict, reverse=False):
  """ applies converter functions belonging to 'packageName' and to 
  version 'fromVersion' or higher versions found to target,
  and returns result. If there is no appropriate func, returns target 
  unchanged.
  
  Checks that packageName either has a func entry or is the name of a
  current package with XML I/O, and throws error otherwise.
  
  infoDict is a {packageName:[(version,func),]} dictionary of lists of tuples
  
  If reverse is False, the functions are applied in the order they appear
  (correct behaviour for conversion functions).
  If reverse is True, the functions are applied in reverse order of appearance
  (correct behaviour for map-modifying functions).
  
  """
  
  currentVersion = genConstants.currentModelVersion
  
  if fromVersion > currentVersion:
    # is processed version compatible with current API version?
    raise ApiError(" Attempt to process %s file with %s  code"
                   % (fromVersion, currentVersion))
                   
  elif fromVersion == currentVersion:
    # no corrections to be made
    return target
  
  funcList = infoDict.get(packageName)
  
  if funcList is None:
    #no funcList for packageName
    
    ll = packageName.split('.')
    ll.insert(1,'xml')
    try:
      # check that a package named PackageName exists and has XML I/O
      dummy = __import__('.'.join(ll))
      

      # it does. Return unchanged
      return target
    
    except ImportError:
      raise ApiError(
       " No XML conversion or I/O found for package %s" %  packageName
      )
  
  else:
    # func found.
    
    funcList = funcList
    funcList2 = []
    
    for vv,func in funcList:
    
      if vv > fromVersion or (vv == fromVersion and vv < cutoffVersion):
        
        if func == nonConvertible:
          return ApiError(
           " %s version %s file not compatible with current XML I/O" 
           % (packageName,target)
          )
        
        else:
          funcList2.append(func)
      
    if reverse:
      funcList2.reverse()
    for func in funcList2:
      # convert
      target =  func(target)
    
    # 
    return target

#############################################################################
#
# specific conversion functions


# NB versions must be  in *increasing* order (newest last) in the lists.
# The lists are sorted automatically but for readability they should be 
# ordered in the file also. 

#########################
#
# ccp.BmrbEntry

# NB package changing name. Data are shared under both names

mapMakers = mappingMakers['ccp.BmrbEntry'] = []
mappingMakers['ccp.Entry'] = mappingMakers['ccp.BmrbEntry']

converters = dataConverters['ccp.BmrbEntry'] = []
dataConverters['ccp.Entry'] = dataConverters['ccp.BmrbEntry']

#
version = Version.getVersion('1.0.b2')

def map_ccp_BmrbEntry_1_0_b2(mapping):
  """ rename 'Entry' to BmrbEntry'
  """
  
  dd = semideepcopy(mapping)
  dd['package']['name'] = 'ccp.Entry'
  
  return dd

mapMakers.append((version, map_ccp_BmrbEntry_1_0_b2))

#
version = Version.getVersion('1.0.210')

def map_ccp_BmrbEntry_1_0_210(mapping):
  """ convert laboratories link (to Organisation) to new one to Group
  """
  
  import ccp.api.People
  
  dd = semideepcopy(mapping)
  dd['Entry']['laboratories'] = {'type':'exo', 'end':'/laboratories>', 
   'tag':'laboratories', 'class':ccp.api.People.Organisation, 'sub':'exoN',
   'key1':{'type':'attrkey', 'end':'/key1', 'tag':'key1', 'xtag':'serial',
    'hicard':1, 'toStr':basicDataTypes.Int.toString,
    'fromStr':basicDataTypes.Int.fromString
   },
   'keyDict':{0: 'key1'}
  }
  # attribute renaming
  dd['Entry']['types']['renames']['laboratories'] = '_tmp_laboratories'
  
  return dd

mapMakers.append((version, map_ccp_BmrbEntry_1_0_210))

def conv_ccp_BmrbEntry_1_0_210(headObject):
  """ convert laboratories link (to Organisation) to new one to Group
  """
  
  for entry in headObject.bmrbEntries:
    if hasattr(entry,'_tmp_laboratories'):
      organisations = entry._tmp_laboratories
      del entry._tmp_laboratories
      groups = []
      for organisation in organisations:
        groups.extend(organisation.groups)
      entry.laboratories = groups
  
  return headObject

converters.append((version, conv_ccp_BmrbEntry_1_0_210))


#########################
#
# ccp.ChemCompCoord

mapMakers = mappingMakers['ccp.ChemCompCoord'] = []
converters = dataConverters['ccp.ChemCompCoord'] = []

#
version = Version.getVersion('1.0.113')

def conv_ccp_ChemCompCoord_1_0_113(headObject):
  """ convert old type descriptor to new
  """
  
  # current conversion
  for xx in headObject.chemCompCoords:
    chemComp = xx.chemComp
    for yy in xx.chemCompVarCoord:
      yy.__dict__['chemCompDescriptor'] = get_new_descriptor_113(
       chemComp, yy.chemCompLinking, yy.chemCompDescriptor
      )
  
  return headObject
  
converters.append((version, conv_ccp_ChemCompCoord_1_0_113)) 

#
version = Version.getVersion('1.0.211')

def map_ccp_ChemCompCoord_1_0_211(mapping):
  """ 
  """
  #
  dd = semideepcopy(mapping)
  
  
  dd['ChemCompVarCoord']['chemCompLinking'] = {
   'type':'text1', 'end':'/chemCompLinking>',
  }
  dd['ChemCompVarCoord']['types']['renames']['chemCompLinking'] = 'linking'
  
  dd['ChemCompVarCoord']['chemCompDescriptor'] = {
   'type':'text1', 'end':'/chemCompDescriptor>',
  }
  dd['ChemCompVarCoord']['types']['renames']['chemCompDescriptor'] = 'descriptor'
  
  #
  return dd

mapMakers.append((version, map_ccp_ChemCompCoord_1_0_211))

def conv_ccp_ChemCompCoord_1_0_211(headObject):
  """ 
  """
  
  
  return headObject
  
converters.append((version, conv_ccp_ChemCompCoord_1_0_211)) 


#########################
#
# ccp.ChemComp

mapMakers = mappingMakers['ccp.ChemComp'] = []
converters = dataConverters['ccp.ChemComp'] = []

#
version = Version.getVersion('1.0.113')

def map_ccp_ChemComp_1_0_113(mapping):
  """ Minor changes
  """
  dd = semideepcopy(mapping)
  dd['StdChemComp']['extendedCode'] = dd['NonStdChemComp']['extendedCode'] = {
   'type':'skip', 'end':'/extendedCode>'
  }
  dd['ChemAtom']['isIonizable'] = {'type':'skip', 'end':'/isIonizable>'}
  #
  return dd

mapMakers.append((version, map_ccp_ChemComp_1_0_113))

#
version = Version.getVersion('1.0.200')

def map_ccp_ChemComp_1_0_200(mapping):
  """ Splitting of ChemAtom-LinkEnd link into two.
  NB 200 may be slightly wrong as version number, but this is a temporary expedient
  and ChemComps, being centrallygenerated, do not exist for most versions.
  """
  dd = semideepcopy(mapping)
  dd['ChemAtom']['linkEnds'] = {'type':'skip', 'end':'/linkEnds>'}
  dd['LinkAtom']['linkEnd'] = {'type':'skip', 'end':'/linkEnd>'}
  dd['LinkEnd']['chemAtoms'] = {'type':'linkN', 'end':'/chemAtoms>'}
  dd['LinkEnd']['linkAtoms'] = {'type':'linkN', 'end':'/linkAtoms>'}
  dd['LinkEnd']['types']['multiLinks'].append('chemAtoms')
  dd['LinkEnd']['types']['multiLinks'].append('linkAtoms')
  #
  return dd

mapMakers.append((version, map_ccp_ChemComp_1_0_200))

def conv_ccp_ChemComp_1_0_200(headObject):
  """ clean up link-splitting for LinkEnd.chemAtoms 
  """
  
  # current conversion
  for linkEnd in headObject.chemComp.linkEnds:
    if hasattr(linkEnd,'chemAtoms'):
      aas = linkEnd.chemAtoms
      if aas:
        setattr(linkEnd,'boundChemAtom',aas[0])
        if len(aas) > 1:
          setattr(linkEnd,'remoteChemAtom',aas[1])
      #
      del linkEnd.chemAtoms
      
    if hasattr(linkEnd,'linkAtoms'):
      aas = linkEnd.linkAtoms
      if aas:
        setattr(linkEnd,'boundLinkAtom',aas[0])
        if len(aas) > 1:
          setattr(linkEnd,'remoteLinkAtom',aas[1])
      #
      del linkEnd.linkAtoms
        
  
  return headObject
  
converters.append((version, conv_ccp_ChemComp_1_0_200))

#
version = Version.getVersion('1.0.201')

def map_ccp_ChemComp_1_0_201(mapping):
  """ Changes to conform to new NmrExpPrototype package
  """
  #
  dd = semideepcopy(mapping)
  
  dd['ChemAtom']['nuclGroupType'] = {'type':'skip', 'end':'/nuclGroupType>'}
  
  #
  return dd

mapMakers.append((version, map_ccp_ChemComp_1_0_201))


version = Version.getVersion('1.0.208')

def map_ccp_ChemComp_1_0_208(mapping):
  """ msdCode reorganised (now derived), renamed varMsdCode
  """
  #
  dd = semideepcopy(mapping)
  
  
  dd['ChemCompVar']['msdCode'] = {'type':'text1', 'end':'/msdCode>',}
  dd['ChemCompVar']['types']['renames']['msdCode'] = 'varMsdCode'
  
  #
  return dd

mapMakers.append((version, map_ccp_ChemComp_1_0_208))


version = Version.getVersion('1.0.211')

def map_ccp_ChemComp_1_0_211(mapping):
  """ renamed waterNmrObs to waterExchangeable, inverting logic
  """
  #
  dd = semideepcopy(mapping)
  
  
  dd['ChemAtom']['waterNmrObs'] = {'type':'attr1',
   'toStr':bool2str,
   'fromStr':str2bool,
    'default':True
  }
  dd['ChemAtom']['types']['renames']['waterNmrObs'] = 'waterExchangeable'
  
  #
  return dd

mapMakers.append((version, map_ccp_ChemComp_1_0_211))

def conv_ccp_ChemComp_1_0_211(headObject):
  """ renamed waterNmrObs to waterExchangeable, inverting logic
  """
  
  # current conversion
  for chemAtom in headObject.chemComp.chemAtoms:
    if hasattr(chemAtom,'waterExchangeable'):
      chemAtom.waterExchangeable = not chemAtom.waterExchangeable
      
  return headObject
  
converters.append((version, conv_ccp_ChemComp_1_0_211))

#########################
#
# ccp.Citation

mapMakers = mappingMakers['ccp.Citation'] = []
converters = dataConverters['ccp.Citation'] = []


version = Version.getVersion('1.0.211')

def map_ccp_Citation_1_0_211(mapping):
  """ removed isPrimary attribute
  """
  #
  dd = semideepcopy(mapping)
  
  dd['BookCitation']['isPrimary'] = {'type':'skip','end':'/isPrimary>'}
  dd['ConferenceCitation']['isPrimary'] = dd['BookCitation']['isPrimary']
  dd['JournalCitation']['isPrimary'] = dd['BookCitation']['isPrimary']
  dd['ThesisCitation']['isPrimary'] = dd['BookCitation']['isPrimary']
  #
  return dd

mapMakers.append((version, map_ccp_Citation_1_0_211))

#########################
#
# ccp.ChemElement

# NB package changing name. Data are shared under both names

mapMakers = mappingMakers['ccp.ChemElement'] = []
mappingMakers['ccp.Nucleus'] = mappingMakers['ccp.ChemElement']

converters = dataConverters['ccp.ChemElement'] = []
dataConverters['ccp.Nucleus'] = dataConverters['ccp.ChemElement']

#
version = Version.getVersion('1.0.104')

def map_ccp_ChemElement_1_0_104(mapping):
  """ rename 'Nucleus' to ChemElement'
  """
  dd = semideepcopy(mapping)
  
  dd['package']['name'] = 'ccp.Nucleus'
  
  dd['Nucleus'] = dd['package']['root']['Nucleus'] = dd['ChemElement']
  dd['Nucleus']['end'] = '/Nucleus>'
  dd['Nucleus']['elementSymbol'] = {'type':'text1', 'end':'/elementSymbol>'}
  dd['Nucleus']['types']['renames']['elementSymbol'] = 'symbol'
  
  return dd

mapMakers.append((version, map_ccp_ChemElement_1_0_104))

def conv_ccp_ChemElement_1_0_104(headObject):
  """ create change case for element symbols
  """
  
  # current conversion
  for chemElement in headObject.chemElements:
    ss = chemElement.symbol
    ss2 = ss[0].upper() + ss[1:].lower()
    chemElement.__dict__['symbol'] = ss2
  
  return headObject
  
converters.append((version, conv_ccp_ChemElement_1_0_104))

#########################
#
# ccp.Method

mapMakers = mappingMakers['ccp.Method'] = []
converters = dataConverters['ccp.Method'] = []

#
version = Version.getVersion('1.0.b3')

def map_ccp_Method_1_0_b3(mapping):
  """ rename 'parameters' to '_tmp_parameters'
  """
  dd = semideepcopy(mapping)
  
  dd['Method']['parameters'] = {'type':'text1', 'end':'/parameters>',}
  dd['Method']['types']['renames']['parameters'] = '_tmp_parameters'
  
  return dd

mapMakers.append((version, map_ccp_Method_1_0_b3))

def conv_ccp_Method_1_0_b3(headObject):
  """ create Parameter child from _tmp_parameters
  """
  from ccp.api.Method import Parameter
  
  # current conversion
  for method in headObject.methods:
    params = method.__dict__.get('_tmp_parameters')
    if params is not None:
      del method._tmp_parameters
      Parameter(method, keyword='parameters', serial=1, value=params)
  
  return headObject
  
converters.append((version, conv_ccp_Method_1_0_b3))


#
version = Version.getVersion('1.0.106')

def map_ccp_Method_1_0_106(mapping):
  """ rename 'keyword' to 'name'
  """
  dd = semideepcopy(mapping)
  
  dd['Parameter']['keyword'] = {'type':'text1', 'end':'/keyword>'}
  dd['Parameter']['types']['renames']['keyword'] = 'name'
  
  return dd

mapMakers.append((version, map_ccp_Method_1_0_106))

def conv_ccp_Method_1_0_106(headObject):
  for method in headObject.methods:
    for parameter in method.parameters:
      if hasattr(parameter,'name') and parameter.name:
        pass
      elif hasattr(parameter,'keyword') and parameter.keyword:
        dd = parameter.__dict__
        dd['name'] = dd['keyword']
        del dd['keyword']
      else:
        parameter.__dict__['name'] = 'parameter'
        
  return headObject
  
converters.append((version, conv_ccp_Method_1_0_106))


#########################
#
# ccp.Molecule

mapMakers = mappingMakers['ccp.Molecule'] = []
converters = dataConverters['ccp.Molecule'] = []

#
version = Version.getVersion('1.0.b3')

def map_ccp_Molecule_1_0_b3(mapping):
  """ remove attributes no longer used
  """
  dd = semideepcopy(mapping)
  
  dd['Polymer']['hasNonStdChemComp'] = {'type':'skip', 'end':'/hasNonStdChemComp>'}
  dd['Polymer']['hasNonStdChirality'] = {'type':'skip', 'end':'/hasNonStdChirality>'}
  dd['Polymer']['hasNonStdLinkage'] = {'type':'skip', 'end':'/hasNonStdLinkage>'}
  dd['Polymer']['hasAmbigResidues'] = {'type':'skip', 'end':'/hasAmbigResidues>'}
  dd['Polymer']['hasAmbigConfs'] = {'type':'skip', 'end':'/hasAmbigConfs>'}
  dd['Polymer']['isAromatic'] = {'type':'skip', 'end':'/isAromatic>'}
 
  dd['NonPolymer']['hasAmbigResidues'] = {'type':'skip', 'end':'/hasAmbigResidues>'}
  dd['NonPolymer']['hasAmbigConfs'] = {'type':'skip', 'end':'/hasAmbigConfs>'}
  dd['NonPolymer']['isAromatic'] = {'type':'skip', 'end':'/isAromatic>'}
  
  dd['MolResidue']['chemCompLoc'] = {'type':'skip','end':'/chemCompLoc>'}
  
  return dd

mapMakers.append((version, map_ccp_Molecule_1_0_b3))

def conv_ccp_Molecule_1_0_b3(headObject):
  """ add default linking
  """
  linkings = ('start','middle','end')
  dd = {}
  for cch in headObject.chemCompHeads:
    mt = cch.molType
    if not mt and hasattr(cch,'_molType'):
      mt = cch._molType
    dd[mt] = None
  if len(dd) == 1:
    molType = dd.keys()[0]
  else:
    molType = None
  for molecule in headObject.molecules:
    if (molType or molecule.molType) in ('protein','DNA','RNA','DNA/RNA'):
      ll = molecule.molResidues
      if ll[0].linking not in linkings:
        ll[0].__dict__['linking'] = 'start'
      if ll[-1].linking not in linkings:
        ll[-1].__dict__['linking'] = 'end'
      for mr in ll[1:-1]:
        if mr.linking not in linkings:
          mr.__dict__['linking'] = 'middle'
  #
  return headObject
  
converters.append((version, conv_ccp_Molecule_1_0_b3))

# Major ChemComp reworking
version = Version.getVersion('1.0.106')

def map_ccp_Molecule_1_0_106(mapping):
  """ 
  """
  dd = semideepcopy(mapping)
  
  # Merge Polymer and NonPolymer into Molecule (which used to be
  # an abstract subclass), and remove a number of attributes. 
  dd['Molecule']['numFreeThiols'] = {'type':'skip', 'end':'/numFreeThiols>'}
  dd['Molecule']['empiricalFormula'] = {'type':'skip', 'end':'/empiricalFormula>'}
  dd['Molecule']['molecularMass'] = {'type':'skip', 'end':'/molecularMass>'}
  dd['Molecule']['formalCharge'] = {'type':'skip', 'end':'/formalCharge>'}
  dd['Molecule']['isParamagnetic'] = {'type':'skip', 'end':'/isParamagnetic>'}
  #
  dd['Polymer'] = dd['package']['root']['Polymer'] = dd['Molecule'].copy()
  dd['Polymer']['end'] = '/Polymer>'
  dd['Polymer']['extinctionCoeff280'] = {'type':'skip', 'end':'/extinctionCoeff280>'}
  dd['Polymer']['polymerDetails'] = {'type':'skip', 'end':'/polymerDetails>'}
  #
  dd['NonPolymer'] = dd['package']['root']['NonPolymer'] = dd['Molecule']
  dd['NonPolymer']['end'] = '/NonPolymer>'
  
  # Rename PolymerSeqBlock to MolSeqFragment
  dd['PolymerSeqBlock'] = dd['package']['root']['PolymerSeqBlock'] = dd['MolSeqFragment']
  dd['PolymerSeqBlock']['end'] = '/PolymerSeqBlock>'
  
  # rename attribute in MolResidue
  dd['MolResidue']['chainLocation'] = {'type':'text1', 'end':'/chainLocation>'}
  dd['MolResidue']['types']['renames']['chainLocation'] = 'linking'
  
  # handle new limits on molType for ChemCompHead
  dd['ChemCompHead']['types']['renames']['molType'] = '_molType'
  
  return dd

mapMakers.append((version, map_ccp_Molecule_1_0_106))

def conv_ccp_Molecule_1_0_106(headObject):
  """ Update Molecule, adding new default values and setting links
  """
  from ccp.api import Molecule
  
  # handle chemCompHeads
  for chemCompHead in headObject.chemCompHeads:
    if hasattr(chemCompHead,'_molType'):
      ss = chemCompHead._molType
      if ss:
        if ss == 'nonpolymer':
          ss = 'other'
        chemCompHead.__dict__['molType'] = ss
      del chemCompHead._molType
        
  # assumes that 1) all molecules either are linear polymers 
  # or have only one residue
  for molecule in headObject.molecules:
    ll = molecule.molResidues
    
    lastnext = None
    for molResidue in ll:
      # assumes 'neutral' is always defined and always the default descrtiptor
      molResidue.__dict__['descriptor'] = 'neutral'
      li = molResidue.linking
      if li in ('start','end','middle'):
        if li != 'start':
          e2 = Molecule.MolResLinkEnd(molResidue,linkCode='prev')
          if lastnext is not None:
            Molecule.MolResLink(molecule,molResLinkEnds=(lastnext,e2))
        if li != 'end': 
          lastnext = Molecule.MolResLinkEnd(molResidue,linkCode='next')
        else:
          lastnext = None
      else:
        lastnext = None
        molResidue.__dict__['linking']='none'
  return headObject
  
converters.append((version, conv_ccp_Molecule_1_0_106))

#
version = Version.getVersion('1.0.113')

def conv_ccp_Molecule_1_0_113(headObject):
  """ convert old type descriptor to new
  """
  
  # current conversion
  for xx in headObject.molecules:
    for yy in xx.molResidues:
      chemComp = yy.chemCompHead.chemComp
      yy.__dict__['descriptor'] = get_new_descriptor_113(
       chemComp, yy.linking, yy.descriptor
      )
  
  return headObject
  
converters.append((version, conv_ccp_Molecule_1_0_113))


#
version = Version.getVersion('1.0.204')

def map_ccp_Molecule_1_0_204(mapping):
  """ rename seqId to serial
  """
  dd = semideepcopy(mapping)
  dd['MolResidue']['seqId'] = dd['MolResidue']['serial'].copy()
  dd['MolResidue']['types']['renames']['seqId'] = 'serial'
  
  return dd

mapMakers.append((version, map_ccp_Molecule_1_0_204))

def conv_ccp_Molecule_1_0_204(headObject):
  """ copy serial to seqCode (taken from seqId)
  """
  
  for molecule in headObject.molecules:
    for molResidue in molecule.molResidues:
      molResidue.seqCode = molResidue.serial
  
  return headObject
  
converters.append((version, conv_ccp_Molecule_1_0_204))


#
version = Version.getVersion('1.0.208')

def conv_ccp_Molecule_1_0_208(headObject):
  """ added isFinalised - should be True for old molecules with chains
  """
  """ Changed to allow for later model change
  for molecule in headObject.molecules:
    if (molecule.chains
      or molecule.testBlueprintComps
      or molecule.blueprintComps
      or molecule.nucTarget
      or molecule.protTarget
      or molecule.refMolComponents
      or molecule.molComponents
      or molecule.refMolCompFeatures
      or molecule.molCompFeatures
    ):
      molecule.__dict__['isFinalised'] = True"""
  for molecule in headObject.molecules:
    if (molecule.chains
      or molecule.nucTargets
      or molecule.protTargets
      or molecule.molComponents
    ):
      molecule.__dict__['isFinalised'] = True
  
  return headObject
  
converters.append((version, conv_ccp_Molecule_1_0_208))


#########################
#
# ccp.MolSystem

mapMakers = mappingMakers['ccp.MolSystem'] = []
converters = dataConverters['ccp.MolSystem'] = []

# major ChemComp update
version = Version.getVersion('1.0.106')
  
def map_ccp_MolSystem_1_0_106(mapping):
  """ 
  """
  dd = semideepcopy(mapping)
  
  # MolSystem attribute changes
  dd['MolSystem']['molSystemTypes'] = {'type':'skip', 'end':'/molSystemTypes>'}
  dd['MolSystem']['numChains'] = {'type':'skip', 'end':'/numChains>'}
  dd['MolSystem']['hasAmbigConfs'] = {'type':'skip', 'end':'/hasAmbigConfs>'}
  dd['MolSystem']['hasAmbigResidues'] = {'type':'skip', 'end':'/hasAmbigResidues>'}
  dd['MolSystem']['paramagnetic'] = {'type':'skip', 'end':'/paramagnetic>'}
  dd['MolSystem']['numFreeThiols'] = {'type':'skip', 'end':'/numFreeThiols>'}
  dd['MolSystem']['molecularMass'] = {'type':'skip', 'end':'/molecularMass>'}
  
  # Chain attribute changes
  dd['Chain']['hasNmrData'] = {'type':'skip', 'end':'/hasNmrData>'}
  
  # Residue attribute changes
  dd['Residue']['chemCompVar'] = {'type':'skip', 'end':'/chemCompVar>'}
  
  # Residue attribute changes
  dd['Atom']['molSystemBonds'] = {'type':'skip', 'end':'/molSystemBonds>'}
  
  # Convert MolSystemBond to MolSystemLink (for later processing)
  dd['MolSystemBond'] = dd['MolSystemLink']
  dd['MolSystemBond']['end'] = '/MolSystemBond>'
  dd['MolSystemBond']['bondType'] = {'type':'skip', 'end':'/bondType>'}
  dd['MolSystemBond']['bondOrder'] = {'type':'skip', 'end':'/bondOrder>'}
  dd['MolSystemBond']['atoms'] = {'type':'linkN', 'end':'/atoms>'}

  return dd
  
  
mapMakers.append((version, map_ccp_MolSystem_1_0_106))


def conv_ccp_MolSystem_1_0_106(headObject):
  """ Update MolSystem
  """
  from ccp.api import MolSystem
  
  for molSystem in headObject.molSystems:
  
    # set up for finding descriptor
    resLinks = {}
    deletedMls = {}
    for ml in molSystem.molSystemLinks:
      for atom in ml.atoms:
        res = atom.residue
        ll = resLinks.get(res)
        if ll is None:
          resLinks[res] = ml
        else:
          ll.append(ml)
  
    # find linking and descriptor for residue
    for chain in molSystem.chains:
      for residue in chain.residues:
      
        # set linking.
        residue.linking = residue.molResidue.linking
        
        # find appropriate descriptor
        chemComp = residue.molResidue.chemComp
        chemCompVars = {}
        for ccv in chemComp.findAllChemCompVars(linking=residue.linking):
          chemCompVars[ccv] = None
        # first select chemCompVars that have all the atoms
        for atom in residue.atoms:
          chemAtom = chemComp.findFirstChemAtom(name=atom.name)
          if chemAtom is None:
            raise ApiError("residue %s atom %s not found in chemComp %s"
             % (residue, atom.name, chemComp)
            )
          else:
            ccvs = chemAtom.chemCompVars
            for ccv in chemCompVars.keys():
              if ccv not in ccvs:
                del chemCompVars[ccv]
        if not chemCompVars:
          raise ApiError("residue %s has no chemCompVar that fits its atoms"
           % (residue)
          )
        # now check which chemCompVars fit the MolSystemLinks
        mls = resLinks.get(residue)
        if mls:
          for ml in resLinks[residue]:
            if not deletedMls.has_key(ml):
              atom = ml.atoms[0]
              if atom.residue is not residue:
                atom = ml.atoms[1]
              le = chemComp.findFirstChemAtom(name=atom.name).linkEnds
              for ccv in chemCompVars.keys():
                if ccv not in le.chemCompVars:
                  if len(chemCompVars) < 2:
                    deletedMls[ml] = None
                    print (
                     "WARNING, MolSystemLink to %s removed: No ChemCompVar found"
                     % residue
                    )
                    ml.delete()
                  else:
                    del chemCompVars[ccv]
        
        # choose one of the ChemCompVars (there will be at least one)
        ccvs = chemCompVars.keys()
        if len(ccvs) == 1:
          descriptor = ccvs[0].descriptor
        else:
          descriptors = [x.descriptor for x in ccvs]
          for ss in ('deprotonated','neutral','protonated'):
            if ss in descriptors:
              descriptor = ss
              break
          else:
            descriptor = descriptors[0]
        
        # finally, set descriptor
        residue.descriptor = descriptor
    
    # set MolSystemLinks and make MolSystemLinkEnds  
    for ml in molSystem.molSystemLinks:
    
      # set dihedralAngle from stereochemistry
      if hasattr(ml,'stereochemistry'):
        ss = ml.stereochemistry
        del ml.stereochemistry
        if ss == 'E':
          ml.dihedralAngle = 180
        elif ss == 'Z':
          ml.dihedralAngle = 0
      
      # find and set up MolSystemLinkEnds from atoms
      residues = []
      linkCodes = []
      for atom in ml.atoms:
        res = atom.residue
        residues.append(res)
        ca = atom.chemAtom
        ll = res.chemCompVar.linkEnds
        linkEnds = [x for x in res.chemCompVar.linkEnds
         if x.chemAtoms[0] is ca
        ]
        if len(linkEnds) == 1:
          linkCodes.append(x.linkCode)
        else:
          break
      else:
        # no break - we found all necessary data
        for ii in (1,2):
          MolSystem.MolSystemLinkEnd(
           residues[ii], linkCode=linkCodes[ii], molSystemLink=ml
          )
        continue
      
      # break in last loop - no link can be made, delete MolSystemLink
      print (
       "WARNING, MolsystemLink could not be made between atoms \n%s and %s"
       % ml.atoms
      )
      ml.delete()
    
    # reading management necessary as createChainFragments creates
    # new objects and only works correstly if isReading == False
    storageDict =  molSystem.storage.__dict__
    reading = storageDict['isReading']
    storageDict['isReading'] = False
    for chain in molSystem.chains:
      chain.createChainFragments()
    storageDict['isReading'] = reading

  return headObject
  
converters.append((version, conv_ccp_MolSystem_1_0_106))


# Minor adjustments
version = Version.getVersion('1.0.113')

def map_ccp_MolSystem_1_0_113(mapping):
  """ 
  """
  dd = semideepcopy(mapping)
  
  dd['MolSystem']['chemExchange'] = dd['MolSystem']['hasChemExchange'].copy()
  dd['MolSystem']['types']['renames']['chemExchange'] = 'hasChemExchange'
  dd['MolSystem']['metalIons'] = {'type':'skip'}
  dd['MolSystem']['organicLigands'] = {'type':'skip'}
  
  return dd
  
mapMakers.append((version, map_ccp_MolSystem_1_0_113))

def conv_ccp_MolSystem_1_0_113(headObject):
  """ Update MolSystem, setting new descriptor
  """
  from ccp.api import MolSystem
  
  for molSystem in headObject.molSystems:
  
    # find descriptor for residue
    for chain in molSystem.chains:
      for residue in chain.residues:
        
        # find appropriate descriptor
        chemComp = residue.molResidue.chemComp
        chemCompVars = {}
        for ccv in chemComp.findAllChemCompVars(linking=residue.linking):
          chemCompVars[ccv] = len(ccv.__dict__['chemAtoms'])
        # first select chemCompVars that have all the atoms
        for atom in residue.atoms:
          chemAtom = chemComp.findFirstChemAtom(name=atom.name)
          if chemAtom is None:
            raise ApiError("residue %s atom %s not found in chemComp %s"
             % (residue, atom.name, chemComp)
            )
          else:
            ccvs = chemAtom.chemCompVars
            for ccv in chemCompVars.keys():
              if ccv not in ccvs:
                del chemCompVars[ccv]
        
        # choose one of the ChemCompVars
        if not chemCompVars:
          raise ApiError("residue %s has no chemCompVar that fits its atoms"
           % (residue)
          )
          
        if len(chemCompVars) == 1:
          descriptor = chemCompVars.keys()[0].descriptor
          
        else:
          # choose Var with fewest atoms.
          ll = [(y,x) for x,y in chemCompVars.items()]
          ll.sort()
          descriptor = ll[0][1].descriptor
        
        # finally, set descriptor
        residue.__dict__['descriptor'] = descriptor

  return headObject
  
converters.append((version, conv_ccp_MolSystem_1_0_113))

#########################
#
# ccp.Nmr

mapMakers = mappingMakers['ccp.Nmr'] = []
converters = dataConverters['ccp.Nmr'] = []

#
version = Version.getVersion('1.0.104')

def map_ccp_Nmr_1_0_104(mapping):
  """ Splitting off NmrConstraints package
  """
  
  dd = semideepcopy(mapping)
  dd2 = semideepcopy(
   __import__('ccp.xml.NmrConstraints',{},{},['mapping']).mapping
  )
  
  # add NmrConstraints in appropriate location
  temp = dd2['package']['root']
  del dd2['package']
  dd.update(dd2)
  dd['NmrConstraintHead']['structures'] = dd['StructureGeneration']['structures']
  dd['StructureGeneration'] = dd['NmrConstraintHead']
  dd['package']['root']['StructureGeneration'] = dd['StructureGeneration']
  global ccp_Nmr_1_0_104_make_NmrConstraintHead
  ccp_Nmr_1_0_104_make_NmrConstraintHead = dd['NmrConstraintHead']['class']
  dd['StructureGeneration']['class'] = map_ccp_Nmr_1_0_104_make_NmrConstraintHead
  dd['StructureGeneration'].update(temp)
  dd['StructureGeneration']['types']['multiExoLinks'].append('structures')
  dd['StructureGeneration']['end'] = '/StructureGeneration>'
  
  # specific changes
  
  # ignore ResStructure and all its contents
  dd['ResStructure']['type'] = 'skip'  
  
  # ignore obsolete links (including obsolete half of many-to-many links)
  # 1) links that are preserved as they are set from the other side
  dd['Experiment']['constraintLists'] = {'type':'skip', 'end':'/constraintLists>'}
  
  dd['Peak']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  
  dd['SubPeak']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  
  # 2) links that are rearranged
  dd['Resonance']['chemShiftConstraints'] = {'type':'skip', 'end':'/chemShiftConstraints>'}
  dd['Resonance']['hBondConstrItems'] = {'type':'skip', 'end':'/hBondConstrItems>'}
  dd['Resonance']['distanceConstrItems'] = {'type':'skip', 'end':'/distanceConstrItems>'}
  dd['Resonance']['dihedralConstraints'] = {'type':'skip', 'end':'/dihedralConstraints>'}
  dd['Resonance']['jCouplingConstrItems'] = {'type':'skip', 'end':'/jCouplingConstrItems>'}
  dd['Resonance']['rdcConstrItems'] = {'type':'skip', 'end':'/rdcConstrItems>'}
  dd['Resonance']['fixedResonanceSets'] = {'type':'skip', 'end':'/fixedResonanceSets>'}
  dd['HBondConstraintItem']['types']['renames']['resonances'] = '_resonances'
  #dd['HBondConstraintItem']['types']['multiLinks'].remove('resonances')
  #dd['HBondConstraintItem']['types']['multiLinks'].append('_resonances')
  dd['DistanceConstraintItem']['types']['renames']['resonances'] = '_resonances'
  #dd['DistanceConstraintItem']['types']['multiLinks'].remove('resonances')
  #dd['DistanceConstraintItem']['types']['multiLinks'].append('_resonances')
  dd['DihedralConstraint']['types']['renames']['resonances'] = '_resonances'
  #dd['DihedralConstraint']['types']['multiLinks'].remove('resonances')
  #dd['DihedralConstraint']['types']['multiLinks'].append('_resonances')
  dd['JCouplingConstraintItem']['types']['renames']['resonances'] = '_resonances'
  #dd['JCouplingConstraintItem']['types']['multiLinks'].remove('resonances')
  #dd['JCouplingConstraintItem']['types']['multiLinks'].append('_resonances')
  dd['RdcConstraintItem']['types']['renames']['resonances'] = '_resonances'
  #dd['RdcConstraintItem']['types']['multiLinks'].remove('resonances')
  #dd['RdcConstraintItem']['types']['multiLinks'].append('_resonances')
  dd['FixedResonanceSet']['types']['renames']['resonances'] = '_resonances'
  #dd['FixedResonanceSet']['types']['multiLinks'].remove('resonances')
  #dd['FixedResonanceSet']['types']['multiLinks'].append('_resonances')
  dd['ChemShiftConstraint']['resonance'] = {'type':'link1',}
  dd['ChemShiftConstraint']['types']['renames']['resonance'] = '_resonance'

  # 3) links that are ignored
  dd['Resonance']['resCoords'] = {'type':'skip', 'end':'/resCoords>'}
  dd['PeakDim']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['SubPeakDim']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['Shift']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['ShiftDifference']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['ShiftAnisotropy']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['HExchRate']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['HExchProtection']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['T1']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['T2']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['T1Rho']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['HomoNoe']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['HeteroNoe']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['JCoupling']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['Rdc']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  dd['DipolarRelaxation']['constraints'] = {'type':'skip', 'end':'/constraints>'}
  
  # 4) link renaming:
  dd['ExpDimRef']['transferFroms'] = dd['ExpDimRef']['expTransfers'].copy()
  dd['ExpDimRef']['types']['renames']['transferFroms'] = 'expTransfers'
  #
  return dd
  
mapMakers.append((version, map_ccp_Nmr_1_0_104))

def map_ccp_Nmr_1_0_104_make_NmrConstraintHead(project,**parameters):
  """ make ContentStorage object for NmrConstraintHead, then create
  and return NmrConstraintHead itself
  """
  
  # make contenteStorage - using sensible defaults
  from memops.api.Implementation import ContentStorage
  # choose a url from teh available ones
  url = project.findFirstUrl(name='data')
  if not url:
    url = project.url
  # make a path, setting name to reflect serial, 
  # or to look different if this is not possible
  ii = parameters['serial']
  path = 'constraintSet%02d/NmrConstraints.xml' % ii
  ii = ii + 1000
  while os.path.exists(os.path.join(url.path,path)):
    path = 'constraintSet%02d/NmrConstraints.xml' % ii
    ii = ii + 1
  cs = ContentStorage(project, url=url, package='ccp.NmrConstraints', path=path)
  cs.__dict__['isReading'] = True
  
  # make and return NmrConstraintHead
  parameters['contentStorage'] = cs
  return ccp_Nmr_1_0_104_make_NmrConstraintHead(project, **parameters)


def conv_ccp_Nmr_1_0_104(headObject):
  from ccp.api.Nmr import StructureGeneration
  from ccp.api.NmrConstraints import ChemShiftConstraintList
  from ccp.api.NmrConstraints import DihedralConstraintList
  from ccp.api.NmrConstraints import FixedResonanceSet
  from ccp.api.NmrConstraints import FixedAtomSet
  
  for conhead in headObject.nmrConstraintHeads:
    
    project = conhead.project
    
    # resolve link to Structures
    strucgen = StructureGeneration(headObject,
     nmrConstraintHead=conhead, serial=conhead.serial
    )
    structures = conhead.__dict__.get('structures')
    if structures is not None:
       strucgen.structures = structures
       del conhead.structures
  
    # resolve links to resonance
    for conlist in conhead.constraintLists:
 
      if isinstance(conlist,ChemShiftConstraintList):
        for con in conlist.constraints:
          con.resonance = conv_ccp_Nmr_1_0_104_FixedResonance(
           conhead, con._resonance
          )
          del con._resonance
      
      elif isinstance(conlist,DihedralConstraintList):
        for con in conlist.constraints:
          con.resonances = [conv_ccp_Nmr_1_0_104_FixedResonance(conhead,x)
           for x in con._resonances
          ]
          del con._resonances
      
      else:
        for con in conlist.constraints:
          for item in con.items:
            item.resonances = [conv_ccp_Nmr_1_0_104_FixedResonance(conhead,x)
             for x in item._resonances
            ]
            del item._resonances
      
    ll = conhead.fixedResonanceSets
    if ll:
      # we have fixedResonanceSets - set link to FixedResonance
      for frset in ll:
        frset.resonances = [conv_ccp_Nmr_1_0_104_FixedResonance(conhead,x)
         for x in frset._resonances
        ]
        del frset._resonances
    else:
      # No FixedResonanceSets - copy over ResonanceSets and AtomSets as needed
      for fres in conhead.fixedResonances:
        resserial = fres.resonanceSerial
        res = project.findFirstResonance(serial=resserial)
        if res is not None:
          resset = res.resonanceSet
          if resset is not None:
            fresset = conhead.findFirstFixedResonanceSet(serial=resset.serial)
            if not fresset:
              fresset = FixedResonanceSet(conhead, serial=resset.serial)
            fres.resonanceSet = fresset
            ll = []
            for atset in resset.atomSets:
              atserial = atset.serial
              fatset = conhead.findFirstFixedAtomSet(serial=atserial)
              if fatset is None:
                fatset = FixedAtomSet(conhead, serial=atserial,
                 name=atset.name, atoms = atset.atoms
                )
              ll.append(fatset)
 
            fresset.atomSets = ll
    
    # unset isReading
    conhead.contentStorage.__dict__['isReading'] = False
  #
  return headObject
          
def conv_ccp_Nmr_1_0_104_FixedResonance(conhead, resonance):
  """ get FixedResonance corresponding to resonance, creating it if necessary
  """
  from ccp.api.NmrConstraints import FixedResonance
  
  fixres = conhead.findFirstFixedResonance(resonanceSerial=resonance.serial)
  
  if fixres is None:
    serial = len(conhead.__dict__['fixedResonances']) + 1
    fixres = FixedResonance(conhead, serial=serial,
     resonanceSerial=resonance.serial, name=resonance.name, 
     isotopeCode=resonance.isotopeCode, details=resonance.details
    )
  #
  return fixres

converters.append((version, conv_ccp_Nmr_1_0_104))

#
version = Version.getVersion('1.0.113')

def map_ccp_Nmr_1_0_113(mapping):
  """ Merged HomoNoe and HeteroNoe into Noe; same with relevant lists
  """
  #
  dd = semideepcopy(mapping)
  
  dd['HomoNoeList'] = dd['package']['root']['HomoNoeList'] = dd['NoeList'].copy()
  dd['HomoNoeList']['end'] = '/HomoNoeList>'
  #
  dd['HeteroNoeList'] = dd['package']['root']['HeteroNoeList'] = dd['NoeList']
  dd['HeteroNoeList']['end'] = '/HeteroNoeList>'
  
  dd['HomoNoe'] = dd['HomoNoeList']['HomoNoe'] = dd['Noe'].copy()
  dd['HomoNoe']['end'] = '/HomoNoe>'
  #
  dd['HeteroNoe'] = dd['HeteroNoeList']['HeteroNoe'] = dd['Noe']
  dd['HeteroNoe']['end'] = '/HeteroNoe>'
  
  # removed attributes
  dd['Resonance']['atomSiteType'] = {'type':'skip', 'end':'/atomSiteType>'}
  dd['ExpDimRef']['relativePhase'] = {'type':'skip', 'end':'/relativePhase>'}
  dd['ExpDimRef']['residueOffsets'] = {'type':'skip', 'end':'/residueOffsets>'}
  
  #
  return dd

mapMakers.append((version, map_ccp_Nmr_1_0_113))

def conv_ccp_Nmr_1_0_113(headObject):
  """ convert old type descriptor to new
  """
  
  # current conversion
  for xx in headObject.resonanceGroups:
    xx.__dict__['descriptor'] = get_new_descriptor_113(
     xx.chemComp, xx.linking, xx.descriptor
    )
  
  return headObject
  
converters.append((version, conv_ccp_Nmr_1_0_113))

#
version = Version.getVersion('1.0.201')

def map_ccp_Nmr_1_0_201(mapping):
  """ Changes to conform to new NmrExpPrototype package
  """
  #
  dd = semideepcopy(mapping)
  
  dd['ExpDimRef']['isLabelled'] = {'type':'skip', 'end':'/isLabelled>'}
  dd['ExpDimRef']['nuclGroupType'] = {'type':'skip', 'end':'/nuclGroupType>'}
  dd['Resonance']['nuclGroupType'] = {'type':'skip', 'end':'/nuclGroupType>'}
  
  #
  return dd
  
mapMakers.append((version, map_ccp_Nmr_1_0_201))


#
version = Version.getVersion('1.0.206')

def conv_ccp_Nmr_1_0_206(headObject):
  """ set links to RefExpDim and RefExpDimRef
  """
  
  # current conversion
  for exp in headObject.nmrExperiments:
    rexp = exp.refExperiment
    if rexp is not None:
      for expDim in exp.expDims:
        rexpDim = rexp.findFirstRefExpDim(dim=expDim.dim)
        expDim.refExpDim = rexpDim
        for expDimRef in expDim.expDimRefs:
          rexpDimRef = rexpDim.findFirstRefExpDimRef(serial=expDimRef.serial)
          expDimRef.refExpDimRef = rexpDimRef
  #
  return headObject
  
converters.append((version, conv_ccp_Nmr_1_0_206))

#
version = Version.getVersion('1.0.207')

def conv_ccp_Nmr_1_0_207(headObject):
  """ clean up PeakDimContribs with duplicate resonances 
  """
  
  for exp in headObject.nmrExperiments:
    for ds in exp.dataSources:
      for pl in ds.peakLists:
        for peak in pl.peaks:
          for peakDim in peak.peakDims:
            dd = {}
            for contrib in peakDim.peakDimContribs:
              if hasattr(contrib,'resonance'):
                key = contrib.resonance
              elif hasattr(contrib,'resonances'):
                key = contrib.resonances
              else:
                continue
              #
              if dd.has_key(key):
                contrib.delete()
              else:
                dd[key] = None
  #
  return headObject
  
converters.append((version, conv_ccp_Nmr_1_0_207))

#
version = Version.getVersion('1.0.212')

def map_ccp_Nmr_1_0_212(mapping):
  """ Changes to conform to new NmrExpPrototype package
  """
  #
  dd = semideepcopy(mapping)
  
  from memops.xml.Implementation import mapping as implMapping
  import memops.api.Implementation
  
  class MyDummyClass(DummyClass):
    siblingValues = 'expSeriesConditions'
  
  # rename SampledDataDim.parameterVaried to conditionVaried
  dd['SampledDataDim']['parameterVaried'] = {'type':'text1', 
                                                  'end':'/parameterVaried>'}
  dd['SampledDataDim']['types']['renames']['parameterVaried'] = 'conditionVaried'
  
  # keep ExpSeries attributes temporarily
  dd['NmrExpSeries']['variableName'] = {'type':'text1', 'end':'/variableName>'}
  dd['NmrExpSeries']['variableUnit'] = {'type':'text1', 'end':'/variableUnit>'}

  dd['NmrExpSeries']['variableDescriptor'] = {'type':'text1', 'end':'/variableDescriptor>'}
  dd['NmrExpSeries']['types']['renames']['variableDescriptor'] = 'details'

  # temporary store old ExpSeriesCondition class
  # slightly modified version of original XML I/O map. APrst may be unnecessary
  dd['ExpSeriesCondition'] = {'type':'class', 'end':'/ExpSeriesCondition>',
   'class':MyDummyClass, 'parent':'nmrExpSeries',
   'types':{'hasElements':True, 'renames':{},
   'singleAttrs': [ 'serial',
    'value',
    'error'
   ],
   'multiAttrs': [],
   'singleTextAttrs': [],
   'headerTextAttrs': [],
   'multiTextAttrs': [],
   'childLinks': [ 'applicationData'
   ],
   'singleLinks': [ 'experiment'
   ],
   'singleLinksLoad': [],
   'multiLinks': [],
   'singleExoLinks': [ 'access'
   ],
   'multiExoLinks': [],
   }
  }
  dd['ExpSeriesCondition']['types']['saving'] = dd['ExpSeriesCondition']['types']

  dd['ExpSeriesCondition']['serial'] = {'type':'attr1',
   'toStr':basicDataTypes.Int.toString,
   'fromStr':basicDataTypes.Int.fromString
  }
  dd['ExpSeriesCondition']['value'] = {'type':'attr1',
   'toStr':basicDataTypes.Float.toString,
   'fromStr':basicDataTypes.Float.fromString
  }
  dd['ExpSeriesCondition']['error'] = {'type':'attr1',
   'toStr':basicDataTypes.Float.toString,
   'fromStr':basicDataTypes.Float.fromString
  }
  dd['ExpSeriesCondition']['experiment'] = {'type':'skip'}
  dd['ExpSeriesCondition']['access'] = {'type':'exo', 'end':'/access>',
   'tag':'access', 'class':memops.api.Implementation.AccessObject, 'sub':'exo1', 'prior':True,
   'key1':{'type':'textkey', 'end':'/key1>', 'tag':'key1',
    'xtag':'name', 'hicard':1},
   'keyDict':{0: 'key1'}
  }
  dd['ExpSeriesCondition']['AppDataBoolean'] = implMapping['AppDataBoolean']
  dd['ExpSeriesCondition']['AppDataDouble'] = implMapping['AppDataDouble']
  dd['ExpSeriesCondition']['AppDataFloat'] = implMapping['AppDataFloat']
  dd['ExpSeriesCondition']['AppDataInt'] = implMapping['AppDataInt']
  dd['ExpSeriesCondition']['AppDataLong'] = implMapping['AppDataLong']
  dd['ExpSeriesCondition']['AppDataString'] = implMapping['AppDataString']
  
  # add handling of Experiment-expSeriesCondition link
  dd['Experiment']['types']['multiLinks'].append('expSeriesConditions')
  dd['Experiment']['expSeriesConditions'] = {'type':'linkN'}
  
  # add 'ExpSeriesCondition' to contents of NmrExpSeries
  dd['NmrExpSeries']['ExpSeriesCondition'] = dd['ExpSeriesCondition']
  
  #
  return dd
  
mapMakers.append((version, map_ccp_Nmr_1_0_212))

def conv_ccp_Nmr_1_0_212(headObject):
  """ clean up PeakDimContribs with duplicate resonances 
  """
  
  for nxs in headObject.nmrExpSeries:
    
    # set conditionNames
    nxs.addConditionName(nxs.variableName)
    
    # convert old ExpSeriesConditions
    if hasattr(nxs,'expSeriesConditions'):
      for xsc in nxs.expSeriesConditions:
        for exp in headObject.nmrExperiments:
          ll = exp.__dict__.get('expSeriesConditions')
          if ll and xsc in ll:
            nxs.addExperiment(exp)
            break
        else:
          raise ("No experiment found for ExpSeriesCondition %s in %s"
           % (xsc.serial, nxs)
          )
        condset = exp.sampleConditionSet
        if condset is None:
          # first calculate correct serial
          serialDict = headObject.__dict__.setdefault('serialDict',{})
          oldserial = serialDict.get('sampleConditionSets',0)
          # then make new condition set
          condset = headObject.newSampleConditionSet(experiments=(exp,),
           details="Generated by XML backwards copmpatibility",
           serial=oldserial+1,
          )
        
        else:
          # check if we already have a condition of this name
          oldCond = condset.findFirstSampleCondition(
                          condition=nxs.variableName
                         )
          if oldCond is None:
            condset.newSampleCondition(condition=nxs.variableName,
                                       unit=nxs.variableUnit, value=xsc.value,
                                       error=xsc.error)
          else:
            #we already have an ExpCondition with this condition 
            if (oldCond.unit == nxs.variableUnit and 
                oldCond.value == xsc.value and oldCond.error == xsc.error):
              
              # the two are compatible, just use the old one
              pass
            
            else:
              # The new value is incompatible with the existing one 
              # warn but use the old one anyway
              print (
               """WARNING experiment %s in %s - condition incompatible with %s
   Using pre-existing SampleCondition instead """
               % (exp,nxs,oldCond)
              )
                                   
      del nxs.expSeriesConditions
    # clean-up
    del nxs.variableName
    del nxs.variableUnit
  
  # further clean-up
  for exp in headObject.nmrExperiments:
    if hasattr(exp, 'expSeriesConditions'):
      del exp.expSeriesConditions
    
  #
  return headObject
  
converters.append((version, conv_ccp_Nmr_1_0_212))
  

#########################
#
# ccp.NmrConstraints

mapMakers = mappingMakers['ccp.NmrConstraints'] = []
converters = dataConverters['ccp.NmrConstraints'] = []

#
version = Version.getVersion('1.0.113')

def map_ccp_NmrConstraints_1_0_113(mapping):

  dd = semideepcopy(mapping)
  dd['Violation']['value'] = {'type':'attr1',
   'toStr':basicDataTypes.Float.toString,
   'fromStr':basicDataTypes.Float.fromString,
  }
  dd['Violation']['types']['renames']['value'] = 'violation'
  dd['Violation']['error'] = {'type':'attr1',
   'toStr':basicDataTypes.Float.toString,
   'fromStr':basicDataTypes.Float.fromString,
  }
  dd['Violation']['types']['renames']['error'] = 'calcValueError'
  
  #
  return dd
  
mapMakers.append((version, map_ccp_NmrConstraints_1_0_113))

def conv_ccp_NmrConstraints_1_0_113(headObject):
  
  for vl in headObject.violationLists:
    for vv in vl.violations:
      if vv.__dict__['violation'] is None:
        vv.violation = 0.0
      val = vl.constraint.__dict__.get('targetValue')
      if val:
        vv.calcValue = vv.violation + val
      else:
        vv.calcValue = 0.0
  
  #
  return headObject
  
converters.append((version, conv_ccp_NmrConstraints_1_0_113))

#
version = Version.getVersion('1.0.212')

def map_ccp_NmrConstraints_1_0_212(mapping):

  dd = semideepcopy(mapping)
  
  dd['Violation']['numViolated'] = {'type':'attr1',
   'toStr':basicDataTypes.Int.toString,
   'fromStr':basicDataTypes.Int.fromString,
  }
  dd['Violation']['types']['renames']['numViolated'] = 'fractionViolated'
  
  #
  return dd
  
mapMakers.append((version, map_ccp_NmrConstraints_1_0_212))

def conv_ccp_NmrConstraints_1_0_212(headObject):
  
  for vl in headObject.violationLists:
    nStructures = 1.0 * len(vl.structures)
    for vv in vl.violations:
      if hasattr(vl,'fractionViolated'):
        vl.fractionViolated = vl.fractionViolated/nStructures
  #
  return headObject
  
converters.append((version, conv_ccp_NmrConstraints_1_0_212))


#########################
#
# ccp.NmrReference

mapMakers = mappingMakers['ccp.NmrReference'] = []
converters = dataConverters['ccp.NmrReference'] = []

#
version = Version.getVersion('1.0.113')

def conv_ccp_NmrReference_1_0_113(headObject):
  """ convert old type descriptor to new
  """
  
  # current conversion
  for xx in headObject.chemCompNmrRefs:
    chemComp = xx.chemComp
    for yy in xx.chemAtomShifts:
      yy.__dict__['chemCompDescriptor'] = get_new_descriptor_113(
       chemComp, yy.chemCompLinking, yy.chemCompDescriptor
      )
  
  return headObject
  
converters.append((version, conv_ccp_NmrReference_1_0_113))


#########################
#
# memops.Implementation

mapMakers = mappingMakers['memops.Implementation'] = []
converters = dataConverters['memops.Implementation'] = []

#
version = Version.getVersion('1.0.b')

def conv_memops_Implementation_1_0_b(headObject):
  """ 
  Removes nmrRefStorage, as this part of the model has changed. 
  Any Nmr Reference data are lost forever, but this part of the model
  was never used in anger in version 1.0b
  """
  
  nmrRefStorage = headObject.findFirstStorage(package = 'ccp.NmrReference')
  
  if nmrRefStorage:
  
    headObject.__dict__['storages'].remove(nmrRefStorage)
  
    url = nmrRefStorage.url
    url.__dict__['storages'].remove(nmrRefStorage)
    
  #
  return headObject
  
converters.append((version, conv_memops_Implementation_1_0_b))

#
version = Version.getVersion('1.0.b2')

def conv_memops_Implementation_1_0_b2(headObject):
  """ rename 'ccp.Entry' storage to 'ccp.BmrbEntry'
  """
  entryStorage = headObject.findFirstStorage(package='ccp.Entry')
  if entryStorage:
    entryStorage.__dict__['package'] = 'ccp.BmrbEntry'
  
  #
  return headObject
  
converters.append((version, conv_memops_Implementation_1_0_b2))

#
version = Version.getVersion('1.0.104')

def conv_memops_Implementation_1_0_104(headObject):
  """ rename 'ccp.Nucleus' storage to 'ccp.ChemElement'
  """
  xx = headObject.findFirstStorage(package='ccp.Nucleus')
  xx.__dict__['package'] = 'ccp.ChemElement'
  #
  return headObject
  
converters.append((version, conv_memops_Implementation_1_0_104))


#########################
#
# ccpnmr.NameMapping

# attribute changing name changing name. 

mapMakers = mappingMakers['ccpnmr.NameMapping'] = []

converters = dataConverters['ccpnmr.NameMapping'] = []

#
version = Version.getVersion('1.0.105')

def map_ccp_NameMapping_1_0_105(mapping):
  """ rename AtomMapping.moleculeCode to moleculeName
  """
  dd = semideepcopy(mapping)
  dd['AtomMapping']['moleculeCode'] = {'type':'text1','end':'/moleculeCode>'}
  dd['AtomMapping']['types']['renames']['moleculeCode'] = 'moleculeName'
  
  return dd

mapMakers.append((version, map_ccp_NameMapping_1_0_105))


#########################
#
# ccpnmr.Analysis

# splitting AtomSetMapping into 
# ChainMapping, ResidueMapping and AtomSetMapping. 

mapMakers = mappingMakers['ccpnmr.Analysis'] = []
converters = dataConverters['ccpnmr.Analysis'] = []

#
version = Version.getVersion('1.0.203')

def map_ccpnmr_Analysis_1_0_203(mapping):
  """ remove AtomSetMapping to force recreation under new model.
  remove StripGroup
  """
  
  dd = semideepcopy(mapping)
  
  # ignore AtomSetMapping and all its contents
  dd['AtomSetMapping']['type'] = 'skip'  
  dd['package']['root']['AtomSetMapping'] = dd['AtomSetMapping']
  
  # ignore StripGroup and all its contents and links to it.
  dd['StripGroup'] = {'type':'skip', 'end':'/StripGroup>'}
  dd['AxisRegionGroup']['stripGroup'] = {'type':'skip', 'end':'/stripGroup>'}
  dd['SpectrumWindow']['stripGroup'] = {'type':'skip', 'end':'/stripGroup>'}
  
  return dd

mapMakers.append((version, map_ccpnmr_Analysis_1_0_203))

#
version = Version.getVersion('1.0.209')

def ccpnmr_Analysis_209_chemAtomSet_linkModifier(objStack, fullKey):
  """ modify exoLink, appending integer 1 
  This allows for change in key of ChemAtomSet, adding subType. 
  """
  result = fullKey[:]
  result.append(1)
  #
  return result

def map_ccpnmr_Analysis_1_0_209(mapping):

  dd = semideepcopy(mapping)
  
  # modify exolink to ChemAtomSet to allow for subType
  dd['AtomSetMapping']['chemAtomSet']['modifierFunc'] = ccpnmr_Analysis_209_chemAtomSet_linkModifier
  #
  return dd
  
mapMakers.append((version, map_ccpnmr_Analysis_1_0_209))

###########################################################################
#
# final processing:

for ll in mappingMakers.values():
  ll.sort()
  
for ll in dataConverters.values():
  ll.sort()


def get_new_descriptor_113(chemComp,linking,descriptor):
  """ get new-format descriptor for upgrade to version 113
  """
  oxidised_CYS = {
   'start':"deprot:H3;link:SG",
   'middle':"link:SG",
   'end':"prot:H'';link:SG",
   'none':"prot:H'';deprot:H3;link:SG",
  }
  protonated_NE2_HIS = {
   'start':"prot:HE2;deprot:H3,HD1",
   'middle':"prot:HE2;deprot:HD1",
   'end':"prot:HE2,H'';deprot:HD1",
   'none':"prot:HE2,H'';deprot:H3,HD1",
  }
  
  if descriptor == 'neutral':
    # use default descriptor, or 'neutral' if there is none
    var = chemComp.findFirstChemCompVar(linking=linking, isDefaultVar = True)
    if var is None:
      # should not happen, but this just might be better than nothing
      return 'neutral'
    else:
      return var.descriptor
  
  elif descriptor == 'oxidized':
    if not (chemComp.molType == 'protein' and chemComp.ccpCode == 'CYS'):
      raise ApiError("(1) compatibility code cannot handle %s"
       % ((chemComp.molType, chemComp.ccpCode,linking,descriptor))
      ) 
    ss = oxidised_CYS.get(linking)
    if ss:
      return ss
    else:
      # should not happen
     return "link:SG"
  
  elif descriptor == 'protonated_NE2':
    if not (chemComp.molType == 'protein' and chemComp.ccpCode == 'HIS'):
      raise ApiError("(2) compatibility code cannot handle %s" 
       % ((chemComp.molType, chemComp.ccpCode,linking,descriptor))
      ) 
    ss = protonated_NE2_HIS.get(linking)
    if ss:
      return ss
    else:
      # should not happen
     return "prot:HE2"
   
  elif descriptor in ('protonated','deprotonated'):
    vvs = chemComp.findAllChemCompVars(linking=linking)
    descrs = []
    hasNeutral = False
    defDescr = None
    for var in vvs:
      ss = var.descriptor
      if ss ==  'neutral':
        hasNeutral = True
      else:
        xx = parse_descriptor_113(ss)
        descrs.append((ss,xx))
        if var.isDefaultVar:
          defDescr = (ss, xx)
 
    if not descrs:
      if hasNeutral:
        # only 'neutral available - use that
        print ( "WARNING, %s.%s.%s descriptor %s forced to 'neutral'"
         % (chemComp.molType, chemComp.ccpCode, linking, descriptor)
        )
        return 'neutral'
 
    elif descriptor == 'deprotonated':
      
      if defDescr and not defDescr[1].get('prot'):
        # default var for linking is deprotonated - use it
        return defDescr[0]
      
      
      elif chemComp.molType == 'protein' and linking == 'end':
        ii = 0
        result = None
        for tag,dd in descrs:
          ll = dd.get('deprot')
          if len(dd) == 2 and ll and dd.get('prot') == ["H''"]:
            if len(ll) > ii:
              ii = len(ll)
              result = tag
        if result:
          return result
 
      else:
        ii = 0
        result = None
        for tag,dd in descrs:
          ll = dd.get('deprot')
          if len(dd) == 1 and ll:
            if len(ll) > ii:
              ii = len(ll)
              result = tag
        if result:
          return result
 
    elif descriptor == 'protonated':
      
      if defDescr and not defDescr[1].get('deprot'):
        # default var for linking is protonated - use it
        return defDescr[0]
        
      elif chemComp.molType == 'protein' and linking == 'start':
        ii = 0
        result = None
        for tag,dd in descrs:
          ll = dd.get('prot')
          if len(dd) == 2 and ll and dd.get('deprot') == ["H3"]:
            if len(ll) > ii:
              ii = len(ll)
              result = tag
        if result:
          return result
 
      else:
        ii = 0
        result = None
        for tag,dd in descrs:
          ll = dd.get('prot')
          if len(dd) == 1 and ll:
            if len(ll) > ii:
              ii = len(ll)
              result = tag
        if result:
          return result
 
    raise ApiError(
     "No appropriate descriptor for %s,%s,%s old descriptor %s"
     % (chemComp.molType, chemComp.ccpCode, linking, descriptor)
    )
  else:
    # descriptor did not fit - may well new and correct, set by previous
    # compatibility software
    
    return descriptor

def parse_descriptor_113(descriptor):
  """ Parse type-113 new descriptor, returning custom dictionary of results
  """
  
  result = {}
  
  if descriptor == 'neutral':
    raise ApiError(
     "parse_descriptor_113 called with improper descriptor 'neutral'"
    )
 
  for ss in descriptor.split(';'):
    tag,ss2 = ss.split(':')
    ats = ss2.split(',')
    result[tag] = ats
  
  return result
  
##########################################################################3
#
#
#  Special conversion functions
#  carried out after reading Implementation package
#
# Designed for operations that need to access several fully loaded packages
# and cannot be called as part of the loading machinery.
#

def getVersion(stream):
  # get version of input stream
  
  # read file and split in elements
  elementStrings = stream.read().split('<')
  
  # find and process _StorageUnit element, 
  for ee in elementStrings:
    if ee[:13] == "_StorageUnit ":
      strings = ee[14:].split('"')
      attrs = {}
      for ii in range(0,len(strings)-1,2):
        attrs[strings[ii][:-1].lstrip()] = strings[ii+1] or None
      
      break
  else:
    raise ApiError(
     loadError + "not (correct) XML. Does not contain <_StorageUnit element"
    )
  
  # check and set up release compatibility
  try:
    fileVersion = attrs['release']
    
  except:
    raise ApiError(loadError + 
     "not (correct) XML. <_StorageUnit element lacks 'release' "
    )
  return Version.getVersion(fileVersion)


#
def specialConversions(project):
  """ Wrapper routine to run special conversions.
  
  NBNB apparently not currrently called from anywhere. Should be.
  """
  
  # converting for Molecule.seqId change starting versin 1.0.204
  ms = project.findFirstStorage(package='ccp.Molecule')
  stream = ms.url.getStream(ms.package, path=ms.path, headObject=project)
  moleculeVersion = getVersion(stream)
 
  if moleculeVersion < Version.getVersion('1.0.204'):
    special_conv_1_0_204(project)


def special_conv_1_0_204(project):
  """ special function to fix Molecule.MolResidue.serial
  
  This functionality requires special handling as multiple packages must
  be loaded before the fix can be carried out.
  
  It is assumed that all Molecules have consecutive serial numbers (after the
  Molecule compatibility function has been executed) and that sequential type
  links (prev-next) are always between consecutive MolResidues. This is
  highly likely to be the case in practice and simplifies the coding.
  The function then resets the serials to start at 1.
  
  The function changes the values of keys (serials) for MolResidue.
  It first loads all packages with linke to MolResidues (to make sure links are
  not broken), and change all attributes in that refer to MolResidue srials 
  (e.g. MolSystem.Residue.seqId). As some of these are keys, the same procedure
  is repeated for the new packages 
  """
  
  molDict = {}
  # find molecules to correct for, and set up
  for molecule in project.molecules:
    ll = [x.serial for x in molecule.molResidues]
    ll.sort()
    if ll and ll[0] != 1:
      corr = ll[0]-1
      molDict[molecule] = corr
      oldSerial = molecule.__dict__['serialDict']['molResidues']
      molecule.__dict__['serialDict']['molResidues'] = oldSerial - corr
  
  if molDict:
    # we need to make corrections
 
    from ccp.api import RefSampleComponent
 
    # load all required packages up front.
    loadPackages = ('ccp.RefSampleComponent','ccp.Annotation',
    'ccp.MolSystem', 'ccp.Nmr', 'ccpnmr.Analysis'
    )
    for storage in project.storages:
      if (storage.package in loadPackages and storage.isStored
       and not storage.isLoaded
      ):
        storage.load()

    # fix ExpBlueprint.BlueprintComponent
    for expBlueprint in project.expBlueprints:
      for blueprintComponent in expBlueprint.blueprintComponents:
        corr = molDict.get(blueprintComponent.molecule)
        if corr is not None:
          nn = blueprintComponent.approxBeginSeqId
          if nn is not None:
            blueprintComponent.approxBeginSeqId = nn - corr
          nn = blueprintComponent.approxEndSeqId
          if nn is not None:
            blueprintComponent.approxEndSeqId = nn - corr

    # fix RefSampleComponent.MolFeature
    for component in project.abstractComponents:
      if isinstance(component,RefSampleComponent.MolComponent):
        for molCompFeature in component.molCompFeatures:
          corr = molDict.get(molCompFeature.molecule)
          if corr is not None:
            nn = molCompFeature.startSeqId
            if nn is not None:
              molCompFeature.startSeqId = nn - corr
            nn = molCompFeature.endSeqId
            if nn is not None:
              molCompFeature.endSeqId = nn - corr

    # fix Molecules
    for molecule,corr in molDict.items():
      for molResidue in molecule.molResidues:
        nn = molResidue.serial
        molResidue.__dict__['serial'] = nn - corr

    # ready to fix ccpnmr.Analysis.ResidueMapping
    # NB the package might not be present
    try:
      chainMappings = project.chainMappings
    except:
      chainMappings = []

    # ready to fix ccpnmr.NameMapping
    # NB the package might not be present
    try:
      nameMappings = project.nameMappings
    except:
      nameMappings = []

    # fix MolSystem.Residue
    for molSystem in project.molSystems:

      molSystemCode = molSystem.code

      for chain in molSystem.chains:
        corr = molDict.get(chain.molecule)
        if corr:
          for residue in chain.residues:
            nn = residue.seqId
            residue.__dict__['seqId'] = nn - corr

          chainCode = chain.code

          # fix Coordinate.Residue
          for structure in molSystem.structures:
            for coordChain in structure.coordChains:
              if coordChain.code == chainCode:
                for residue in coordChain.residues:
                  nn = residue.seqId
                  residue.__dict__['seqId'] = nn - corr

          # fix ccpnmr.Analysis.ResidueMapping
          for chainMapping in chainMappings:
            if chainMapping.chainCode == chainCode:
              for residueMapping in chainMapping.residueMappings:
                nn = residueMapping.seqId
                residueMapping.__dict__['seqId'] = nn - corr

          for nameMapping in nameMapping:

            # fix ccpnmr.NameMapping.ChainMapping
            for chainMapping in nameMapping.chainMappings:
              if chainMapping.chainCode == chainCode:
                nn = chainMapping.firstSeqId
                chainMapping.__dict__['firstSeqId'] = nn - corr

            # fix NameMapping.ChainMapping.AtomMapping
            for atomMapping in nameMapping.atomMapping:
              if (atomMapping.molSystemCode == molSystemCode
               and atomMapping.chainCode == chainCode
              ):
                nn = atomMapping.seqId
                if nn is not None:
                  atomMapping.__dict__['seqId'] = nn - corr
  #
  return project
