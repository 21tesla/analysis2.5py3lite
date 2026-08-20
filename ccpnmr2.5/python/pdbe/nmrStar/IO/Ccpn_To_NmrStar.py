import re

#from ccp.api.molecule import ChemComp
#from ccp.api.molecule import Molecule

from memops.api import Implementation

from ccp.general.Util import findChemCompSysName, findChemCompVarSysName, findAllSysNamesByChemAtomOrSet
from pdbe.nmrStar.IO.Util import getCcpn2NmrStarConstants
from pdbe.nmrStar.IO.Constants import unknownMapping

#
# HERE: ONLY depend on things that are linked to (a) NmrEntry.Entry in the data model!!
#

"""
derivedDataLists
entryMolecules
experiments
laboratories
measurementLists
molSystem
relatedEntries
structureGenerations
structureGroups
study
"""

class WaterMolResidue:
  def __init__(self,serial,seqCode,chemCompVar,ccpCode):
    self.serial      = serial
    self.seqCode     = seqCode
    self.chemCompVar = chemCompVar
    self.ccpCode     = ccpCode

class IonicStrengthCondition:

  def __init__(self,value,parent):

    self.condition = 'ionic_strength'
    self.value     = value
    self.parent    = parent
    self.unit      = 'mM'

class DataCount:

  def __init__(self,name,count = 1):

    self.name  = name
    self.count = count

  def incrementCount(self):

    self.count += 1

class ChemShiftRefList:

  def __init__(self,Id):

    self.Id = Id
    self.chemShiftRefs = []

  def add(self,chemShiftRef):

    if chemShiftRef not in self.chemShiftRefs:
      self.chemShiftRefs.append(chemShiftRef)

#
# Special class to handle individual atom name writing for measurements (shifts)
# Is basically a complex data type that combines several CCPN classes
# Could be extended for other things - follow through how this is used for shift lists to see...
#

class MeasurementByIndividualAtom:

  def __init__(self,measurement,resonance,resonanceToAtom):
  
    #self.serial = serial
    self.measurement = measurement
    self.resonance = resonance
    self.resonanceToAtom = resonanceToAtom

#
# TODO: should probably move over all code from NmrStarExport.py that can be handled here (with self.exportClass have access),
# so that all can be version specific (never know what's going to change). See getShiftAmbiguityCode as example.
#
# Basically means that all code in the exportClass should be totally generic.
#
# TODO: could also remove most passing in of objects ala ('object',self.getSomeInfo), because is set in ccpnVar anyway. But maybe not as clear...
#

class Ccpn_To_NmrStar:

  mappingHelp = \
                \
  """
  
  This help contains information on how to use the dictionary.
  
  General:
  -------
  
    'ccpnLoop': Gets a list of CCPN objects. Has 3 possible forms:
         
         'ccpnObject(.links)'                   Gets list directly
         ('ccpnObject(.links)',self.function)   Gets list via function that uses left-hand list of other ccpn objects as input
         self.function                          Gets list from specific function
    
    'ccpnMap': Defines the main component for mapping.
               Will be defined as a variable (e.g. 'ccpnMap': 'ccpnObjectVariableName', so can be used in the mapping dictionary once defined.
               In conjunction with 'ccpnLoop', will loop over list of CCPN objects from 'ccpnLoop', and set the variable to each in turn.

     A special combined 'ccpnLoop', 'ccpnMap' exists for nested loops that go 'deep' into a set of CCPN objects,
     but still form one loop within NMR-STAR. They are of the form:
     
       'ccpnLoop': [ccpnObject1List,'ccpnObject1.ccpnObjects2','ccpnObject2.ccpnObjects3'],
       'ccpnMap': [('ccpnObject1','corresponding_NmrStarID'),('ccpnObject2','corresponding_NmrStarID'),('ccpnObject3','corresponding_NmrStarID')],
       
       For example:
       
          'ccpnLoop': ['constraintList.constraints','constraint.items'],
          'ccpnMap': [('constraint','Constraint_ID'),('constraintItem','Node_ID')],
          
       Note that the 'corresponding_NmrStarID' can be left empty.
          
       
       The default behaviour in this loop is that within each nested loop the corresponding_NmrStarID is reset when
       the parent loop goes to the next object (so local IDs to each nested loop). This behaviour can be changed by using:
       
       'nmrStarKeyType': 'ccpnObject' or 'global'
       
           Using 'ccpnObject' defines the level at which corresponding_NmrStarIDs will be reset.
           Using 'global', all corresponding_NmrStarIDs are global and nested loops do not matter for the ID count.       
    
    
  
  SaveFrame level:
  ---------------
  
    self.sfDict['saveFrameName'] = {
  
      'ccpnLoop':  When used on this level, multiple saveframes will be created. Each needs a unique name (see 'title')
      'ccpnMap':   Defines the main CCPN variable used for defining the tags in this saveFrame (e.g. 'study' for the study_list saveframe)

      'title':     Defines the (unique) title of the saveframe. 
                   Usually of the form ('ccpnObject.attribute',self.getTitle), as no spaces can occur in the title itself.

      'tags': {    Defines the tag mappings for this saveFrame. The 'tagName's make up the dictionary keys.

      'tables': {  Defines the table (NMR-STAR loops) mappings for this saveFrame. The 'loopName's make up the dictionary keys.
      
      }
  
  
  Table level:
  -----------
  
    'loopName': {
    
        'ccpnLoop': Is required here, as this is a loop and multiple values are possible.
        'ccpnMap':  Same as on saveFrame level.
        
        'tags':     Defines the tag mappings for this NMR-STAR loop
        
        }       
        
    Note that for specific reoccurring 'standard' loops, functions are available to set these up automatically.
    See the 'Common loop/table definitions' header in this class for the ones that are available, e.g.:
        
      'saveFrameName_common_name': self.setCommonNameLoop('listOfCommonNames','commonName'),
      'saveFrameName_keyword': self.setKeywordLoop('listOfKeywords','keyword'),
      'saveFrameName_systematic_name': self.setSysNamesLoop('listOfCcpnSysNameObjects','ccpnSysName'),
   
   Basically the first argument corresponds to the 'ccpnLoop', the second to the 'ccpnMap'
     
  
  Tag level:
  ---------
    
    'tagName': mapping to CCPN
    
    The mapping to CCPN can be of the form:
    
      'ccpnObject.attributeName'
      ('ccpnObject',self.function)
      self.function

  
    A special function that can be used to set variables local to this class is:
    
      'tagName=LOCAL': localVariable or self.variable or self.functionReturningVariable,
      
    Another special function that can be used to set a range of related tags is:
    
      'CUSTOM_actionName': mapping to CCPN
      
      This will call special code inside NmrStarExport (TODO should be here??) to set a range of tags. Heavily used for resonance-atom links.
  
  
  Notes on setting keys:
  ---------------------
  
  The Export script uses the NmrStar dictionary to decide where keys should be set. These can correspond to CCPN objects,
  or to strings. You can specify what an NMR-STAR ID corresponds to with the above system (either in the special combined
  form of 'ccpnLoop', 'ccpnMap', or on the 'tagName' level).
  
  There is also a way to ignore foreign keys inside the NMR-STAR dictionary (but *only* when outside of the saveframe
  that defines the foreign key). See the 'ignoreForeignKeys' variable in this script for a current list.

  """
  
  #
  # List of links from the NMR entry that are used - is in principle CCPN version specific
  #

  nmrEntryLinkList = ['molSystem','structureGenerations','study','authors',
                      'contactPersons','entryMolecules','experiments']
  
  value = None # Dummy to use in lambda functions
    
  ignoreForeignLinks = ['Chem_comp.ID', 'Citation.ID', 'Assembly.ID', 'Entity_assembly.Entity_assembly_name',
                        'Entity_chimera_segment.ID']

  ignorePrimaryKeys = ['Citation_author.Ordinal', 'Entity_poly_seq.Hetero']
  
  # This allows to search for application data with a different saveframe name! Will be mapped...
  # Note that table names also have to be specifically mapped!
  mapAppDataSaveFrames = {'general_distance_constraints':
  
                              ('distance_constraints',
                               {'Distance_constraint_expt':     'Gen_dist_constraint_expt',
                                'Distance_constraint_software': 'Gen_dist_constraint_software',
                                'Dist_constr_software_setting': 'Gen_dist_constraint_software_param',
                                'Dist_constraint_comment_org':  'Gen_dist_constraint_comment_org',
                                'Dist_constraint_parse_err':    'Gen_dist_constraint_parse_err',
                                'Dist_constraint_parse_file':   'Gen_dist_constraint_parse_file',
                                'Dist_constraint_conv_err':     'Gen_dist_constraint_conv_err'
                               } )
                         }
  
  # Dictionary to track chemAtomSysNames for export (to make sure they are PDB_REMED)
  chemAtomSysNames = {}
  
  customSerials = {}

  def __init__(self,exportClass, ccpnVersion = '1.1.a2', nmrStarVersion = '3.1'):
    
    # TODO here have to check if the CCPN and NMR-STAR versions match for the particular Ccpn2Star class!
    self.ccpnVersion = ccpnVersion
    self.nmrStarVersion = nmrStarVersion
    self.exportClass = exportClass
    
    self.constants_ByVersion = getCcpn2NmrStarConstants(self.ccpnVersion,self.nmrStarVersion)

    self.patt = {}
    self.patt['concValue'] = re.compile("^[^0-9\.]*([0-9\.]+)\s*[A-Za-z]+")
    self.patt['concUnit']  = re.compile("^[^0-9\.]*[0-9\.]+\s*([A-Za-z]+)")

    self.patt['confSub']  = re.compile("nmrConformerSubmitted:\'([^\']+)\'")
    self.patt['confCalc'] = re.compile("nmrConformerCalculated:\'([^\']+)\'")
    self.patt['confBest'] = re.compile("nmrStructureRepresent:\'([^\']+)\'")
    self.patt['confSel']  = re.compile("nmrModelSelection:\'([^\']+)\'")

    self.patt['pdbId'] = re.compile("^[A-Za-z\d]{4}$")

    self.patt['chemFormulaSpace'] = re.compile("(\d)([A-Z])")
    self.patt['chemFormulaOne']   = re.compile("([A-Za-z])1( |$)")

    # This is necessary to track custom created objects, to make sure they are not re-created in the
    # second loop of the NMR-STAR export code. If they are a key, then their value will go missing otherwise
    # because the custom object is re-created.
    self.trackCustomObjects = {}
    self.dataListCounts = {}
    self.dataCounts = {}

    self.setSfDict()
    
    # Necessary for subclasses
    self.modifySfDict()

  def printTag(self,tag):

    print 'TAG: [%s]' % tag

    return tag

  def setSerial(self,ccpnObject):
  
    name = self.exportClass.writeStarTableDict['ccpnMap']
  
    if not self.customSerials.has_key(name):
      self.customSerials[name] = 0
      
    self.customSerials[name] += 1
    
    return self.customSerials[name]

  def getStudyKeywords(self,entry):

    study = entry.study

    #print 'STUDY [%s]' % study

    keywords = study.keywords

    #print 'KEYS: [%s]' % keywords

    return keywords

  def getCitations(self,nmrEntry):

    self.authors = []

    self.nmrEntry = nmrEntry

    self.citations = []

    if nmrEntry.primaryCitation not in self.citations:
      self.citations.append(nmrEntry.primaryCitation)

    for citation in nmrEntry.sortedOtherCitations():
      if citation not in self.citations:
        self.citations.append(citation)

    return self.citations

  def getCitationStarClass(self,citation):

    starClass = 'entry citation'

    if citation == self.nmrEntry.primaryCitation:
      starClass = 'entry citation'

    else:
      starClass = 'reference citation'

    return starClass

  def getCitationStarClassLabel(self,citation):

    return self.getCitationStarClass(citation).replace(' ','_')

  def getCitationStarType(self,className):

    starType = 'journal'

    if className == 'JournalCitation':
      starType = 'journal'

    elif className == 'BookCitation':    
      starType = 'book'

    elif className == 'ThesisCitation':
      starType = 'thesis'

    elif className == 'ConferenceCitation':
      starType = 'abstract'

    return starType

  def getCitationAuthors(self,citation):

    #self.citCount = 0

    #keywds = {'application': 'nmrStar',
    #          'keyword':     'CitCount'}

    #appData = citation.findFirstApplicationData(**keywds)

    #if appData:
    #  self.citCount = appData.value

    #self.count = 0
    #self.startFlag = True

    self.authors = None

    if citation:
      self.authors = list(citation.authors)

    #print 'AUTHORS: [%s]' % [auth.familyName for auth in self.authors]

    return self.authors

  def findAuthorIndex(self,author):

    #if self.startFlag:
    #  self.startFlag = False
    #  self.count = 1
    #else:
    #  self.count += 1

    #return self.count

    #print 'NAME: [%s]' % author.familyName

    self.idx = self.authors.index(author)

    #print 'IDX: [%s] ' % idx

    return (self.idx + 1)

  def getMolSystemName(self,molSystem):

    # TODO: better to hook this up to a class so can use self.molSystem?
    #       or nmrEntry.molSystem??!?!

    molSysName = molSystem.name  
    if not molSysName:
      molSysName = molSystem.code

    return molSysName

  def getOrganicLigands(self,chains):

    numOrganicLigands = 0

    for chain in chains:
      molecule = chain.molecule
      if molecule.molType == 'other' and molecule.seqLength == 1:
        molResidue = molecule.findFirstMolResidue()
        if len(molResidue.chemCompVar.chemAtoms) > 1:
          numOrganicLigands += 1

    return numOrganicLigands

  def getMetalIons(self,chains):

    numMetalIons = 0

    for chain in chains:
      molecule = chain.molecule
      if molecule.molType == 'other' and molecule.seqLength == 1:
        molResidue = molecule.findFirstMolResidue()
        if len(molResidue.chemCompVar.chemAtoms) == 1:
          numMetalIons += 1

    return numMetalIons

  def getThiolState(self,chains):
  
    #
    # Enumerated value!
    #

    numFreeThiols = 0
    numBoundThiols = 0

    for chain in chains:
      for residue in chain.residues:
        if residue.ccpCode.title() == 'Cys':
          # TODO THIS IS NOT WORKING CORRECTLY YET!! Need to set link:SG if any disulfides there - now not happening (does deprot:HG)
          if not residue.chemCompVar.descriptor.count('link:SG'):
            numFreeThiols += 1
          else:
            numBoundThiols += 1            
    
    thiolState = 'not present'
    
    if numFreeThiols:
      if not numBoundThiols:
        thiolState = 'all free'
      else:
        thiolState = 'free and bound'
    elif numBoundThiols:
      thiolState = 'all disulfide bound'

    # TODO: not handling this correctly yet - is for bound metal ions, I assume
    #elif ???:
    #  thiolState = 'all other bound'

    return thiolState

  def getMoleculeThiolState(self,molecule):

    numFreeThiols = 0
    numBoundThiols = 0

    for molRes in molecule.molResidues:
      if molRes.ccpCode.title() == 'Cys':
        if not molRes.chemCompVar.descriptor.count('link:SG'):
          numFreeThiols += 1
        else:
          numBoundThiols += 1

    thiolState = 'not present'

    if numFreeThiols:
      if not numBoundThiols:
        thiolState = 'all free'
      else:
        thiolState = 'free and bound'
    elif numBoundThiols:
      thiolState = 'all disulfide bound'

    # TODO: not handling this correctly yet - is for bound metal ions, I assume
    #elif ???:
    #  thiolState = 'all other bound'

    return thiolState

  def getAmbConfStates(self,molecule):

    ambConfState = False

    for chain in molecule.sortedChains():
      #print 'CHAIN: [%s]' % chain
      count = 0
      for chainStateSet in chain.sortedChainStateSets():
        #print 'CHAINSTATESET: [%s]' % chainStateSet
        if chainStateSet.stateSetType in ('protonation', 'isotopomer'):
        #  print 'CONTINUING'
          continue
        #print 'LENGTH: [%s]' % len(chainStateSet.sortedChainStates() )
        count += len(chainStateSet.sortedChainStates() )

      if count > 1:
        #print 'COUNT: [%s]' % count
        ambConfState = True

    #print 'STATE: [%s]' % ambConfState

    return ambConfState

  def getAmbChemCompSites(self,molecule):

    # Can we handle this in data model?

    return False

  def checkType(self,value):

    if type(value) == type(bool):
      return value
    else:
      if value == 'yes':
        return True
      elif value == 'no':
        return False
      else:
        return None

  def negateBoolean(self,boolean):

    boolean = not boolean

    return boolean

  def convertBoolean(self,boolean):

    return self.constants_ByVersion['boolean'][boolean]
    
  def getCoordinateAttr(self,coord,attrName):
    
    return getattr(coord,attrName)

  # TODO: should go elsewhere
  def getTitle(self,value):

    if value:
      value = value.strip()
      value = value.replace(' ', '_')
      value = value.replace('\t', '_') # And others?
      value = value.replace('\n', '_')
    else:
      value = None

    #print 'VAL: [%s]' % value

    return value  

  def toString(self,value):
    return str(value)

  def getName(self,value):

    if value:
      value = value.replace('_', ' ')

    return value

  # TODO: should go elsewhere
  def getLabel(self,value):

    if value:
      value = self.getTitle(value)
      value = "$" + value

    return value  

  def getValue(self,value):

    # Real chemical shifts can be 0.0 for example, so don't set as None if value is 0 in this case.
    # Otherwise - use getNonZeroValue shown below.

    if type(value) != type(float() ):
      return None

    value = '%.3f' % value

    return value

  def getNonZeroValue(self,value):

    if value == 0:
      return None
    elif type(value) != type(float() ):
      return None

    value = '%.3f' % value

    return value

  def getErrorValue(self,value):

    if value == 0:
      return None
    elif type(value) != type(float() ):
      return None

    value = '%.3f' % value

    return value

  def getAssemblyTitle(self,code):

    return 'assembly'

  def getNumChains(self,molSystem):

    foundWater = False

    numChains = 0

    for chain in molSystem.sortedChains():
      if (len(chain.molecule.molResidues) == 1 and
          chain.molecule.findFirstMolResidue().chemComp.ccpCode == 'Hoh'):

        if not foundWater:
          foundWater = True
          numChains += 1

      else:
        numChains += 1

    return numChains

  def getMolSystemChains(self,molSystem):

    self.molSystem = molSystem

    foundWater = False

    self.chains = []

    for chain in molSystem.sortedChains():
      if (len(chain.molecule.molResidues) == 1 and
          chain.molecule.findFirstMolResidue().chemComp.ccpCode == 'Hoh'):

        if not foundWater:
          foundWater = True
          if chain not in self.chains:
            self.chains.append(chain)

      else:
        if chain not in self.chains:
          self.chains.append(chain)

    return self.chains

  def getExptDataReported(self,chain):

    experimentDataObjList = ('datums', 'dipolarRelaxations', 'hExchProtections',
                             'hExchRates', 'isotropicS2s', 'jCouplings', 'noes', 'pkas',
                             'rdcs', 'shiftAnisotropies', 'shiftDifferences', 'shifts',
                             'spectralDensities', 't1Rhos', 't1s', 't2s')
    
    for res in chain.sortedResidues():
      for atom in res.sortedAtoms():
        if not atom.atomSet:
          continue
        for resonanceSet in atom.atomSet.sortedResonanceSets():
          for resonance in resonanceSet.sortedResonances():
            for attrName in experimentDataObjList:
              for item in getattr(resonance, attrName):
                if hasattr(item,'peaks') and item.peaks:
                  return True
                if hasattr(item.parent,'experiments') and item.parent.sortedExperiments():
                  return True

    return False

  def getNonStdBondBool(self,molSystem):

    nonStdBondFlag = 'no'

    for chain in molSystem.sortedChains():
      for mrl in chain.molecule.sortedMolResLinks():
        if not mrl.isStdLinear:
          nonStdBondFlag = 'yes'
          break

    return nonStdBondFlag

  def getDeletedAtoms(self,molSystem):

    self.deletedAtoms = []

    for chain in molSystem.sortedChains():
      for mrl in chain.molecule.sortedMolResLinks():
        if not mrl.isStdLinear:
          for mrle in mrl.sortedMolResLinkEnds():
            if mrle not in self.deletedAtoms:
              self.deletedAtoms.append(mrle)

    return self.deletedAtoms

  def getDeletedAtomName(self,molResLinkEnd):

    res = molResLinkEnd.molResidue.ccpCode.upper()
    link = molResLinkEnd.linkCode

    atomName = None

    # TODO: Need a dictionary of deleted atom names here depending on residue type and linkCode.

    if res == 'CYS' and link == 'SG':
      atomName = 'HG'

    return atomName

  def getNonStdBonds(self,molSystem):

    self.nonStdBonds = []

    for chain in molSystem.sortedChains():
      for mrl in chain.molecule.sortedMolResLinks():
        if not mrl.isStdLinear and mrl not in self.nonStdBonds:
          self.nonStdBonds.append(mrl)

    return self.nonStdBonds

  #def getAllChains(self,molecule):

  #  chains = self.molSystem.findAllChains(molecule = molecule)

  #  return tuple(chains)

  def getBondType(self,molResLink):

    linkEnd1 = molResLink.sortedMolResLinkEnds()[0]
    linkEnd2 = molResLink.sortedMolResLinkEnds()[1]

    res1 = linkEnd1.molResidue.ccpCode.upper()
    res2 = linkEnd2.molResidue.ccpCode.upper()

    link1 = linkEnd1.linkCode
    link2 = linkEnd2.linkCode

    bondType = None

    # TODO: Need a dictionary of bond types here.

    if res1 == 'CYS' and res2 == 'CYS' and link1 == 'SG' and link2 == 'SG':
      bondType = 'disulfide'

    return bondType

  def getBondOrder(self,molResLink):

    linkEnd = molResLink.findFirstMolResLinkEnd().linkEnd

    atomPair = frozenset( (linkEnd.boundChemAtom, linkEnd.boundLinkAtom) )

    chemBond = linkEnd.chemComp.findFirstChemBond(chemAtoms=atomPair)

    bondOrder = None

    # TODO:
    # bondOrder has these choices: 'single', 'double', 'triple', 'aromatic', 'dative', 'singleplanar'
    # what are the corresponding NMR-STAR names?
    
    if chemBond:
      bondOrder = chemBond.bondType

    nmrStarBondOrder = bondOrder

    if bondOrder == 'single':
      nmrStarBondOrder = 'SING'
    elif bondOrder == 'double':
      nmrStarBondOrder = 'DOUB'

    return bondOrder

  def getFirstBondMolResidue(self,molResLink):

    return molResLink.sortedMolResLinkEnds()[0].molResidue

  def getFirstBondMolResID(self,molResLink):

    return molResLink.sortedMolResLinkEnds()[0].molResidue.serial

  def getFirstBondMolResCode(self,molResLink):

    return molResLink.sortedMolResLinkEnds()[0].molResidue.ccpCode.upper()

  def getFirstBondLinkCode(self,molResLink):

    return molResLink.sortedMolResLinkEnds()[0].linkCode

  def getSecondBondMolResidue(self,molResLink):

    return molResLink.sortedMolResLinkEnds()[1].molResidue

  def getSecondBondMolResID(self,molResLink):

    return molResLink.sortedMolResLinkEnds()[1].molResidue.serial

  def getSecondBondMolResCode(self,molResLink):

    return molResLink.sortedMolResLinkEnds()[1].molResidue.ccpCode.upper()

  def getSecondBondLinkCode(self,molResLink):

    return molResLink.sortedMolResLinkEnds()[1].linkCode

  def getDbRefExptMeth(self,dbRef):

    exptMeth = None

    keywds = {'application': 'nmrStar',
              'keyword':     'exptMethod'}

    appData = dbRef.findFirstApplicationData(**keywds)

    if appData:
      exptMeth = appData.value

    return exptMeth

  def getDbRefStructRes(self,dbRef):

    structRes = None

    keywds = {'application': 'nmrStar',
              'keyword':     'structResolution'}

    appData = dbRef.findFirstApplicationData(**keywds)

    if appData:
      structRes = appData.value

    return structRes

  def getMolSysECNumber(self,molSystem):

    ecNumber = None

    molSystemSysName = molSystem.findFirstMolSystemSysName(namingSystem='EC')

    if molSystemSysName:
      ecNumber = molSystemSysName.name

    return ecNumber

  def getMolecules(self,molSystem):

    self.molecules = []

    for chain in molSystem.sortedChains():
      if chain.molecule not in self.molecules:
        self.molecules.append(chain.molecule)

    return self.molecules

  def getMoleculeName(self,molecule):

    molName = molecule.longName

    if not molName:
      molName = molecule.name

    return molName

  def getApproxMolType(self, molecule):

    molType = molecule.molType

    if molType not in ('protein','DNA','RNA','DNA/RNA','carbohydrate'):
      residueTypes = {}
      residueTypeHighest = [None,0]
      for molRes in molecule.sortedMolResidues():
        resType = molRes.molType

        if not residueTypes.has_key(resType):
          residueTypes[resType] = 0

        residueTypes[resType] += 1

        if not residueTypeHighest[0] or residueTypeHighest[1] < residueTypes[resType]:
          residueTypeHighest = [resType,residueTypes[resType] ]

      if residueTypeHighest[0] != 'other':
        if residueTypes.has_key('other') and (
          (residueTypes['other'] < 3 and len(molecule.molResidues) > 5) or
          (residueTypes['other'] * 1.0 / len(molecule.molResidues) < 0.1) ):
          molType = residueTypeHighest[0]

    return molType

  def getMoleculeStarType(self,molecule):

    polymerType = 'non-polymer'

    if len(molecule.molResidues) > 1 and molecule.isStdLinear:
      polymerType = 'polymer'

    elif len(molecule.molResidues) > 1 and not molecule.isStdLinear:

      molType = self.getApproxMolType(molecule)

      if molType in ('protein','DNA','RNA','DNA/RNA','carbohydrate'):
        polymerType = 'polymer'

    return polymerType

  def getMoleculeStarPolymerType(self,molecule):

    molType = self.getApproxMolType(molecule)

    if self.constants_ByVersion['molTypes'].has_key(molType):
      starPolymerType = self.constants_ByVersion['molTypes'][molType]
    else:
      starPolymerType = unknownMapping

    return starPolymerType

  def getChainCodes(self,molecule):

    codes = []

    for chain in molecule.sortedChains():
      if chain.code not in codes:
        codes.append(chain.code)

    codesString = ','.join(codes)

    return codesString

  def getSeqString(self,molecule):

    seqString = molecule.seqString

    seqString = seqString.replace('*', 'X')

    if not seqString or seqString == '':
      seqString = ''
      for molRes in molecule.sortedMolResidues():
        code1Let = molRes.chemComp.code1Letter
        if not code1Let or code1Let == '':
          code1Let = 'X'
        seqString += code1Let

    if len(seqString) > 20:

      i = 0
      newSeqString = '\n'

      for c in seqString:
        i += 1
        newSeqString += c

        if i == 20:
          newSeqString += '\n'
          i = 0

      seqString = newSeqString

    return seqString

  def getMoleculeECNumber(self,molecule):

    ecNumber = None

    molSysName = molecule.findFirstMoleculeSysName(namingSystem='EC')

    if molSysName:
      ecNumber = molSysName.name

    return ecNumber

  def getSeqMutation(self,molecule):

    seqMutation = None

    keywds = {'application': 'nmrStar',
              'keyword':     'Mutation'}

    appData = molecule.findFirstApplicationData(**keywds)

    if appData:
      seqMutation = appData.value

    return seqMutation

  def getSeqFragment(self,molecule):

    seqFragment = None

    keywds = {'application': 'nmrStar',
              'keyword':     'Fragment'}

    appData = molecule.findFirstApplicationData(**keywds)

    if appData:
      seqFragment = appData.value

    return seqFragment

  def getSeqSrcMethod(self,molecule):

    seqSrcMethod = None

    keywds = {'application': 'nmrStar',
              'keyword':     'Src_method'}

    appData = molecule.findFirstApplicationData(**keywds)

    if appData:
      seqSrcMethod = appData.value

    return seqSrcMethod

  # Various old sub-routines to do with cross references
  #   - getDbAlignments is the one in use.

  """
  def getDbEntries(self,entry):

    self.dbEntries = []

    for db in entry.root.sortedDatabases():
      for entry in db.entries:
        self.dbEntries.append(entry)

    return self.dbEntries

  def getDbRefs(self,molecule):

    self.dbRefs = []

    for ali in molecule.sortedAlignments():
      if hasattr(ali, 'dbRef') and ali.dbRef not in self.dbRefs:
        self.dbRefs.append(ali.dbRef)

    return self.dbRefs

  def getComponentDbRefs(self,molecule):

    self.componentDbRefs = []

    for molComp in molecule.sortedMolComponents():
      for compDbRef in molComp.sortedComponentDbRefs():
        if compDbRef not in self.componentDbRefs:
          self.componentDbRefs.append(compDbRef)

    return self.componentDbRefs
  """


  def getDbAlignments(self,molecule):

    self.alignments = []
    codes = []

    # TODO: maybe look for the alignment with the most info
    #   - rather than just the first (probably the one from a PDB file).

    for ali in molecule.sortedAlignments():
      if ali.dbRef:
        code = ali.dbRef.code
      if ali.dbRef.code in codes:
        continue
      if ali not in self.alignments:
        if code:
          codes.append(code)
        self.alignments.append(ali)

    return self.alignments

  def getStartLimitResidue(self,molSeqFragment):

    startLimitResidue = -1

    if hasattr(molSeqFragment, 'limitResidues'):
      startLimitResidue = molSeqFragment.sortedLimitResidues()[0]

    return startLimitResidue

  def getEndLimitResidue(self,molSeqFragment):

    endLimitResidue = -1

    if hasattr(molSeqFragment, 'limitResidues'):
      endLimitResidue = molSeqFragment.sortedLimitResidues()[1]

    return endLimitResidue

  """
  def getPdbSeqStart(self,details):

    start = ''

    if details:
      start = details.split(':')[0]

    return start

  def getPdbSeqEnd(self,details):

    end = ''

    if details:
      end = details.split(':')[1]

    return end
  """

  def getNaturalSources(self,nmrEntry):

    self.naturalSources = []

    for entMol in nmrEntry.sortedEntryMolecules():
      #print 'ENTMOL: [%s]' % entMol
      if entMol.molecule and hasattr(entMol.molecule, 'naturalSource'):
        naturalSource = entMol.molecule.naturalSource
        if naturalSource not in self.naturalSources:
          self.naturalSources.append(naturalSource)

    return self.naturalSources

  def getExperimentalSources(self,nmrEntry):

    self.exptSources = []

    for entMol in nmrEntry.sortedEntryMolecules():
      if entMol.experimentalSource:
        if entMol.experimentalSource not in self.exptSources:
          self.exptSources.append(entMol.experimentalSource)

    return self.exptSources

  def getMolResidues(self,molecule):

    # Tracks waters to 'fake' molResidues, not sure if still necessary with new handling of all waters in same (or a couple of) chain(s).
    # Modified for multiple waters in one chain - assuming that if first molecule in chain is water, so is the rest. (Wim 25/02/11)
    
    molResidues = []

    if molecule.findFirstMolResidue().chemComp.ccpCode == 'Hoh':
      
      if not hasattr(self,'residueToWaterMolRes'):
        self.residueToWaterMolRes = {}
        
      refRes = molecule.findFirstChain().findFirstResidue()
      
      if refRes in self.residueToWaterMolRes:
        for chain in molecule.chains:
          for residue in chain.sortedResidues():
            molResidues.append(self.residueToWaterMolRes[residue])

      else:
        waterMolResCount = 1
        origToNewMolRes = {}
        tempMolResCCV = molecule.findFirstMolResidue().chemCompVar
        
        for chain in molecule.sortedChains():
          for residue in chain.sortedResidues():
            if not residue.molResidue in origToNewMolRes:
              waterMolRes = WaterMolResidue(waterMolResCount,residue.seqCode,tempMolResCCV,tempMolResCCV.chemComp.ccpCode)
              origToNewMolRes[residue.molResidue] = waterMolRes         
              waterMolResCount += 1
            else:
              waterMolRes = origToNewMolRes[residue.molResidue]
  
            self.residueToWaterMolRes[residue] = waterMolRes
            molResidues.append(waterMolRes)
 
    else:
      for molResidue in molecule.sortedMolResidues():
        molResidues.append(molResidue)

    return molResidues

  def getAuthSeqId(self,molResidue):

    authSeqId = None

    if hasattr(molResidue, 'applicationData') and molResidue.applicationData:
      appData = molResidue.findFirstApplicationData(application='nmrStar',keyword='authorSeqCode')

      if appData:
        authSeqId = appData.value

    return authSeqId

  def getActualMolResidue(self,residue):

    if not residue:
      return

    tempMolRes = residue.molResidue

    if tempMolRes.chemComp.ccpCode == 'Hoh':

      chain = residue.chain
      waterMolRes = self.residueToWaterMolRes[residue]

      tempMolRes = waterMolRes

    return tempMolRes
  
  def getModelCoordinate(self,coordAtom):

    coordinate = coordAtom.newCoord(model = self.exportClass.ccpnVar['model'])
    
    return [coordinate]

  def getDepartmentAndInstitution(self,person):

    currentGroup = None

    if hasattr(person,'currentGroup'):
      currentGroup = getattr(person,'currentGroup')

    if not currentGroup:
      currentPinGGroup = person.currentPersonInGroup
      if currentPinGGroup:
        currentGroup = currentPinGGroup.group

    # Should it be this code instead here???

    #currentPersonInGroup = None # or findFirstPersonInGroup() ???

    #if hasattr(person,'currentPersonInGroup'):
    #  currentPersonInGroup = getattr(person,'currentPersonInGroup')
    #  if currentPersonInGroup:
    #    currentGroup = currentPersonInGroup.group

    fullText = None
    
    if currentGroup:
      fullText = currentGroup.name

      #if currentGroup.name != currentGroup.organisation.name:
      #  fullText += "\n" + currentGroup.organisation.name
    
    return fullText

  def getStateProv(self,org):

    stateProv = None

    if not org:
      return None

    keywds = {'application': 'nmrStar',
              'keyword':     'State'}

    appData = org.findFirstApplicationData(**keywds)

    if appData:
      stateProv = appData.value

    return stateProv

  def getStarTitle(self, title):

    starTitle = ''

    if title == 'Dr.' or title == 'Dr':
      starTitle = 'Dr.'
    elif title == 'Mr.' or title == 'Mr':
      starTitle = 'Mr.'
    elif title == 'Prof.' or title == 'Prof' or title == 'Professor':
      starTitle = 'Prof.'
    elif title == 'Mrs.' or title == 'Mrs' or title == 'Madam':
      starTitle = 'Mrs.'
    elif title == 'Miss' or title == 'Ms.' or title == 'Ms':
      starTitle = 'Ms.'
    elif title == 'Sir':
      starTitle = 'Sir' # Not in NmrStar, but better to keep this title!

    return starTitle

  def getMiddleInitials(self,middleInitials):

    tmpmidInit = [ init.strip('.') + '.' for init in middleInitials ]

    return ''.join(tmpmidInit)

  def getLaboratories(self,nmrEntry):

    laboratories = []

    for lab in nmrEntry.sortedLaboratories():
      appData = lab.findFirstApplicationData(application='nmrStar', keyword='PinG Lab')
      if appData and appData.value:
        continue 
      if lab.name == lab.organisation.name:
        if lab not in laboratories:
          laboratories.append(lab)

    return laboratories

  def getContactOrganisations(self,nmrEntry):

    # Could go through nmrEntry.laboratories as well.

    organisations = []

    for cp in nmrEntry.sortedContactPersons():
      for ping in cp.sortedPersonInGroups():
        organisations.append(ping.group.organisation)

    return organisations

  def getDataListCounts(self,nmrEntry):

    measurementListDict = {
      'DipolarRelaxList':    'dipole_dipole_relaxation',
      'HExchProtectionList': 'H_exch_protection_factors',
      'HExchRateList':       'H_exch_rates',
      'JCouplingList':       'coupling_constants',
      'NoeList':             'heteronucl_NOEs',
      'RdcList':             'RDCs',
      'ShiftAnisotropyList': 'chem_shift_anisotropy',
      'ShiftDifferenceList': None,
      'ShiftList':           'assigned_chemical_shifts',
      'T1List':              'heteronucl_T1_relaxation',
      'T1RhoList':           'heteronucl_T1rho_relaxation',
      'T2List':              'heteronucl_T2_relaxation'}

    derivedDataListDict = {
      'DataList':            'other_data_types',
      'IsotropicS2List':     'order_parameters',
      'PkaList':             'pH_titration',
      'SpectralDensityList': 'spectral_density_values'}

    for ml in nmrEntry.sortedMeasurementLists():
      if not self.trackCustomObjects.has_key(ml.className):
        self.trackCustomObjects[ml.className] = {}

      if self.trackCustomObjects[ml.className].has_key(ml):
        continue

      if not ml.className in self.dataListCounts:
        self.dataListCounts[ml.className] = DataCount(measurementListDict[ml.className])
      else:
        self.dataListCounts[ml.className].incrementCount()

    for ml in nmrEntry.sortedMeasurementLists():
      if not self.trackCustomObjects[ml.className].has_key(ml):
        self.trackCustomObjects[ml.className][ml] = self.dataListCounts[ml.className]

    for ddl in nmrEntry.sortedDerivedDataLists():
      if not self.trackCustomObjects.has_key(ddl.className):
        self.trackCustomObjects[ddl.className] = {}

      if self.trackCustomObjects[ddl.className].has_key(ddl):
        continue

      if not ddl.className in self.dataListCounts:
        self.dataListCounts[ddl.className] = DataCount(measurementListDict[ddl.className])
      else:
        self.dataListCounts[ddl.className].incrementCount()

    for ddl in nmrEntry.sortedDerivedDataLists():
      if not self.trackCustomObjects[ddl.className].has_key(ddl):
        self.trackCustomObjects[ddl.className][ddl] = self.dataListCounts[ddl.className]

    for exp in nmrEntry.sortedExperiments():
      for ds in exp.sortedDataSources():
        for pl in ds.sortedPeakLists():
          if not self.trackCustomObjects.has_key(pl.className):
            self.trackCustomObjects[pl.className] = {}

          if self.trackCustomObjects[pl.className].has_key(pl):
            continue

          if not pl.className in self.dataListCounts:
            self.dataListCounts[pl.className] = DataCount('spectral_peak_list')
          else:
            self.dataListCounts[pl.className].incrementCount()

    for exp in nmrEntry.sortedExperiments():
      for ds in exp.sortedDataSources():
        for pl in ds.sortedPeakLists():
          if not self.trackCustomObjects[pl.className].has_key(pl):
            self.trackCustomObjects[pl.className][pl] = self.dataListCounts[pl.className]

    return sorted(self.dataListCounts.values() )

  def getDataCounts(self,nmrEntry):

    dataCountClasses = ('ShiftList', 'NoeList', 'T1List', 'T1RhoList', 'T2List',
                        'HExchProtectionList', 'HExchRateList','JCouplingList', 'RdcList',
                        #'DipolarRelaxList', 'ShiftAnisotropyList', 'ShiftDifferenceList'
                        )

    measurementListNameDict = {
      'HExchProtectionList': 'H exchange protection factors',
      'HExchRateList':       'H exchange rates',
      'JCouplingList':       'coupling constants', # ???
      'NoeList':             'heteronuclear NOE values',
      'RdcList':             'residual dipolar couplings',
      'T1List':              'T1 relaxation values',
      'T1RhoList':           'T1rho relaxation values', # ???
      'T2List':              'T2 relaxation values',
#      'DipolarRelaxList':    
#      'ShiftAnisotropyList': 
#      'ShiftDifferenceList': 
      }

    for ml in nmrEntry.sortedMeasurementLists():

      if ml.className not in dataCountClasses:
        continue

      if ml.className == 'ShiftList':

        for measurement in ml.sortedMeasurements():
          isotopeCode = measurement.resonance.isotopeCode

          if isotopeCode == '1H':
            name = '1H chemical shifts'
          elif isotopeCode == '13C':
            name = '13C chemical shifts'
          elif isotopeCode == '15N':
            name = '15N chemical shifts'
          elif isotopeCode == '31P':
            name = '31P chemical shifts'
          else:
            continue

          if not self.trackCustomObjects.has_key(name):
            self.trackCustomObjects[name] = {}

          if self.trackCustomObjects[name].has_key(measurement):
            continue

          if name not in self.dataCounts:
            self.dataCounts[name] = DataCount(name)
          else:
            self.dataCounts[name].incrementCount()

      else:
        name = measurementListNameDict[ml.className]

        for measurement in ml.sortedMeasurements():

          if not self.trackCustomObjects.has_key(name):
            self.trackCustomObjects[name] = {}

          if self.trackCustomObjects[name].has_key(measurement):
            continue

          if name not in self.dataCounts:
            self.dataCounts[name] = DataCount(name)
          else:
            self.dataCounts[name].incrementCount()

    for ml in nmrEntry.sortedMeasurementLists():
      for measurement in ml.sortedMeasurements():
        if ml.className == 'ShiftList':
          isotopeCode = measurement.resonance.isotopeCode

          if isotopeCode == '1H':
            name = '1H chemical shifts'
          elif isotopeCode == '13C':
            name = '13C chemical shifts'
          elif isotopeCode == '15N':
            name = '15N chemical shifts'
          elif isotopeCode == '31P':
            name = '31P chemical shifts'
          else:
            continue

        else:
          name = measurementListNameDict[ml.className]

        if not self.trackCustomObjects[name].has_key(measurement):
          self.trackCustomObjects[name][measurement] = self.dataCounts[name]

    for ddl in nmrEntry.sortedDerivedDataLists():

      if ddl.className != 'IsotropicS2List':
        continue

      name = 'order parameters'

      for dd in ddl.sortedDerivations():
        for data in dd.sortedDerivedData():

          if not self.trackCustomObjects.has_key(name):
            self.trackCustomObjects[name] = {}

          if self.trackCustomObjects[name].has_key(data):
            continue

          if name not in self.dataCounts:
            self.dataCounts[name] = DataCount(name)
          else:
            self.dataCounts[name].incrementCount()

    for ddl in nmrEntry.sortedDerivedDataLists():

      if ddl.className != 'IsotropicS2List':
        continue

      name = 'order parameters'

      for dd in ddl.sortedDerivations():
        for data in dd.sortedDerivedData():
          if not self.trackCustomObjects[name].has_key(data):
            self.trackCustomObjects[name][data] = self.dataCounts[name]

    return sorted(self.dataCounts.values() )

  def getRelatedDbName(self,dbId):

    dbIdStr = str(dbId)

    relatedDbName = 'BMRB'

    searchObj = self.patt['pdbId'].search(dbIdStr)

    if searchObj and not dbIdStr.isdigit():
      relatedDbName = 'PDB'

    return relatedDbName

  def getNonStdChemCompVarList(self):
  
    # TODO: Could in principle use previous loop to get this info as self.nonStdChemComps or something.
  
    nonStdChemCompVars = []
    
    for molecule in self.molecules:
      for molRes in molecule.sortedMolResidues():
        chemCompVar = molRes.chemCompVar
        if chemCompVar.chemComp.className == 'NonStdChemComp':
          if not chemCompVar in nonStdChemCompVars:
            nonStdChemCompVars.append(chemCompVar)
            
    return nonStdChemCompVars
  
  def getChemCompVarStarType(self,nonStdChemCompVar):

    starType = 'non-polymer'

    if nonStdChemCompVar.chemComp.ccpCode == 'Hoh':

      starType = 'water'

    elif nonStdChemCompVar.chemComp.molType == 'carbohydrate':
            # 'D-saccharide 1,4 and 1,4 linking'
            # 'L-saccharide 1,4 and 1,4 linking'
            # 'D-saccharide 1,4 and 1,6 linking'
            # 'L-saccharide 1,4 and 1,6 linking'
            # 'L-saccharide'
            # 'D-saccharide'

      starType = 'saccharide'

    elif nonStdChemCompVar.chemComp.molType == 'protein':
      # TODO set to D- when available!!
      starType = 'L-peptide '
      linking = nonStdChemCompVar.linking
      if linking != 'none':
        if linking == 'middle':
          starType += 'linking'
        elif linking == 'start':
          starType += 'NH3 amino terminus'
        elif linking == 'end':
          starType += 'COOH carboxy terminus'
     
    elif nonStdChemCompVar.chemComp.molType in ('DNA','RNA'):
      starType = nonStdChemCompVar.chemComp.molType + ' '
      linking = nonStdChemCompVar.linking
      if linking != 'none':
        if linking == 'middle':
          starType += 'linking'
        elif linking == 'start':
          starType += 'OH 5 prime terminus'
        elif linking == 'end':
          starType += 'OH 3 prime terminus'

        # 'other'

    return starType

  def getChemCompFormula(self,formula):

    starFormula = None

    if formula:
      starFormula = self.patt['chemFormulaSpace'].sub('\\1 \\2', formula)

      starFormula = self.patt['chemFormulaOne'].sub('\\1\\2', starFormula)

    return starFormula

  def getCCVTitle(self,chemCompVar):
  
    chemCompSysName = findChemCompSysName('CIF',chemCompVar.chemComp)
    if not chemCompSysName:
      chemCompSysName = findChemCompSysName('CUSTOM',chemCompVar.chemComp) # This for automatically created ChemComps
      if not chemCompSysName:
        chemCompSysName = chemCompVar.chemComp.ccpCode.upper()
  
    return self.getTitle('chem_comp_' +  chemCompSysName)
    #return findChemCompSysName('CIF',chemCompVar.chemComp)

  def getChemCompCif(self,chemComp):
    
    return findChemCompSysName('CIF',chemComp)

  def getChemCompVarCif(self,chemCompVar):
    
    return findChemCompSysName('CIF',chemCompVar.chemComp)
    
  def getPdbCoordAtomName(self,coordAtom):
  
    return self.getPdbAtomName(coordAtom.atom)

  def getPdbAtomName(self,atom):
  
    residue = atom.residue
    atomName = atom.name
    
    # Track in dictionary for speed purposes
    if not self.chemAtomSysNames.has_key(residue):
      chemCompVar = residue.chemCompVar
      chemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemCompVar.chemComp,chemCompVar.chemAtoms,'PDB_REMED')

      self.chemAtomSysNames[residue] = (chemAtomSysNames, {})    
    
    if not self.chemAtomSysNames[residue][1].has_key(atomName):
    
      sysAtomName = None
      
      for casn in self.chemAtomSysNames[residue][0]:
        
        if casn.atomName == atomName:
          # Atom name found, get it and break out of loop
          sysAtomName = casn.sysName
          self.chemAtomSysNames[residue][0].pop(self.chemAtomSysNames[residue][0].index(casn))
          break
    
      # Return CCPN/IUPAC atom name if not found  
      if not sysAtomName:
        sysAtomName = atomName
        
      self.chemAtomSysNames[residue][1][atomName] = sysAtomName
      
    return self.chemAtomSysNames[residue][1][atomName]
    

  def getStarIsotopeLabelling(self,component):

    bmrbLabelDict = {'NatAbun':      'natural abundance',
                     'uni_15N':      '[U-95% 15N]',
                     'uni_15N13C':   '[U-95% 13C; U-95% 15N]',
                     'uni_15N13C2H': '[U-95% 13C; U-95% 15N; U-95% 2H]'}

    molecule = None

    if hasattr(component, 'molecule'):
      molecule = component.molecule

    else:
      return

    labelMix = component.labeledMixture

    if not labelMix:
      return

    labelMol = labelMix.parent

    molLabel = labelMol.findFirstMolLabel()

    firstBmrbLabel = ''
    ccpCodes = []
    resLabFracFlag = False
    atomLabFlag = False

    for resLabel in molLabel.sortedResLabels():
      bmrbLabel = ''

      resLabFrac = resLabel.findFirstResLabelFraction()

      if resLabFrac:
        if atomLabFlag:
          print '  Warning: mixture of labelling types in molecule %s' % component.molecule.name
          break

        resLabFracFlag = True

        bmrbLabel = bmrbLabelDict[resLabFrac.schemeName]
        #print 'BMRB: [%s]' % bmrbLabel

      else:
        if resLabFracFlag:
          print '  Warning: mixture of labelling types in molecule %s' % component.molecule.name
          break

        atomLabFlag = True


        atomLabel13 = resLabel.findFirstAtomLabel(massNumber=13)

        if atomLabel13:
          if bmrbLabel:
            bmrbLabel += '; '

          if atomLabel13.className == 'UniformAtomLabel' and atomLabel13.elementName == 'C':

            if atomLabel13.weight == 0.00:
              bmrbLabel += 'U-%s%s' % (atomLabel13.massNumber, atomLabel13.elementName)
            else:
              percent = atomLabel13.weight*100
              bmrbLabel += 'U-%d%% %s%s' % (percent, atomLabel13.massNumber, atomLabel13.elementName)

          else:

            if atomLabel13.weight == 0.00:
              bmrbLabel += '%s%s' % (atomLabel13.massNumber, atomLabel13.atomName)
            else:
              percent = atomLabel13.weight*100
              bmrbLabel += '%d%% %s%s' % (percent, atomLabel13.massNumber, atomLabel13.atomName)


        atomLabel15 = resLabel.findFirstAtomLabel(massNumber=15)

        if atomLabel15:
          if bmrbLabel:
            bmrbLabel += '; '

          if atomLabel15.className == 'UniformAtomLabel' and atomLabel15.elementName == 'N':

            if atomLabel15.weight == 0.00:
              bmrbLabel += 'U-%s%s' % (atomLabel15.massNumber, atomLabel15.elementName)
            else:
              percent = atomLabel15.weight*100
              bmrbLabel += 'U-%d%% %s%s' % (percent, atomLabel15.massNumber, atomLabel15.elementName)

          else:

            if atomLabel15.weight == 0.00:
              bmrbLabel += '%s%s' % (atomLabel15.massNumber, atomLabel15.atomName)
            else:
              percent = atomLabel15.weight*100
              bmrbLabel += '%d%% %s%s' % (percent, atomLabel15.massNumber, atomLabel15.atomName)


        atomLabel2 = resLabel.findFirstAtomLabel(massNumber=2)

        if atomLabel2:
          if bmrbLabel:
            bmrbLabel += '; '

          if atomLabel2.className == 'UniformAtomLabel' and atomLabel2.elementName == 'H':

            if atomLabel2.weight == 0.00:
              bmrbLabel += 'U-%s%s' % (atomLabel2.massNumber, atomLabel2.elementName)
            else:
              percent = atomLabel2.weight*100
              bmrbLabel += 'U-%d%% %s%s' % (percent, atomLabel2.massNumber, atomLabel2.elementName)

          else:

            if atomLabel2.weight == 0.00:
              bmrbLabel += '%s%s' % (atomLabel2.massNumber, atomLabel2.atomName)
            else:
              percent = atomLabel2.weight*100
              bmrbLabel += '%d%% %s%s' % (percent, atomLabel2.massNumber, atomLabel2.atomName)


      if bmrbLabel:

        if firstBmrbLabel:

          if firstBmrbLabel != bmrbLabel:
            print '  Warning: multiple labels for this molecule %s' % component.molecule.name
            break

        else:
          firstBmrbLabel = bmrbLabel

        if molecule and len(molecule.molResidues) > 1:
          molResidue = component.molecule.findFirstMolResidue(serial=resLabel.resId)

          if molResidue:
            ccpCode = component.molecule.findFirstMolResidue(serial=resLabel.resId).ccpCode

            molType = self.getApproxMolType(molecule)

            if molType in ('DNA', 'RNA', 'DNA/RNA'):

              nucTlcDict = {'A': 'Ade',
                            'C': 'Cyt',
                            'G': 'Gua',
                            'T': 'Thy',
                            'U': 'Ura'}

              ccpCode = nucTlcDict[ccpCode]

            ccpCodes.append(ccpCode)

    if firstBmrbLabel and firstBmrbLabel != 'natural abundance':
      firstBmrbLabel = '[' + firstBmrbLabel + ']'

      if ccpCodes:
        firstCcpCode = ccpCodes[0]
        resFlag = True

        for ccpCode in ccpCodes:
          if ccpCode != firstCcpCode:
            resFlag = False
            break

        if resFlag:
          firstBmrbLabel += '-' + firstCcpCode

    return firstBmrbLabel

  def getCompMolName(self,component):

    compMolName = ''

    if component:
      compMolName = component.details

    if not compMolName and hasattr(component, 'molecule'):
      molecule = component.molecule

      if molecule:
        compMolName = molecule.name

    return compMolName

  def getComponentMolecule(self,component):

    molecule = None

    if hasattr(component, 'molecule'):
      molecule = component.molecule

    return molecule

  def getComponentMoleculeName(self,component):

    molName = None

    if hasattr(component, 'molecule'):
      molecule = component.molecule

      if molecule:
        molName = self.getLabel(molecule.name)

    return molName

  def getComponentConc(self,component):

    conc = component.concentration

    #print 'CONC: [%s]' % conc

    if component.concentrationUnit == 'M': # and component.concDisplayUnit == 'mM':
      conc = conc*1000

    #if component.concentrationUnit == 'M': # and component.concDisplayUnit == 'uM':
    #  conc = conc*1000000

    elif component.concentrationUnit == 'm3/m3': # and component.concDisplayUnit == '%':
      conc = conc*100

    if conc == 0.0:
      appData = component.findFirstApplicationData(application = 'nmrStar', keyword = 'concentration')
      if appData:
        conc = appData.value

        if conc is None or conc == 'None':
          conc = None
      else:
        conc = None

    return conc

  def getComponentUnit(self,component):

    unit = component.concentrationUnit

    if component.concentrationUnit == 'M': # and component.concDisplayUnit == 'mM':
      unit = 'mM' #component.concDisplayUnit

    #if component.concentrationUnit == 'M': # and component.concDisplayUnit == 'uM':
    #  unit = component.concDisplayUnit

    elif component.concentrationUnit == 'm3/m3': # and component.concDisplayUnit == '%':
      unit = '%'

    return unit
    
  def getSamples(self,experiments):

    self.samples = []

    for exp in experiments:
      if exp.sample not in self.samples:
        self.samples.append(exp.sample)

    return self.samples

  def getSampleConcValue(self,details):

    searchObj = self.patt['concValue'].search(details)

    concValue = 0

    if searchObj:
      concValue = searchObj.group(1)

    return concValue

  def getSampleConcUnit(self,details):

    searchObj = self.patt['concUnit'].search(details)

    concUnit = ''

    if searchObj:
      concUnit = searchObj.group(1)

    return concUnit

  def getSampleConditionSetsOld(self,experiments):

    self.sampleCondSets = []

    if not self.trackCustomObjects.has_key('ionicStrength'):
      self.trackCustomObjects['ionicStrength'] = {}

    for exp in experiments:
      if exp.sampleConditionSet not in self.sampleCondSets:
        self.sampleCondSets.append(exp.sampleConditionSet)

      if exp.sample and exp.sample.ionicStrength:
        
        # Make sure doesn't exist already. If not, then make sure that you track it.
        if not self.trackCustomObjects['ionicStrength'].has_key(exp.sampleConditionSet):
          ionicCond = IonicStrengthCondition(exp.sample.ionicStrength, exp.sampleConditionSet)
          self.trackCustomObjects['ionicStrength'][exp.sampleConditionSet] = ionicCond

    return self.sampleCondSets

  def getSampleConditionsOld(self,scs):

    self.sampleConds = []

    issc = None

    if self.trackCustomObjects['ionicStrength'].has_key(scs):
      issc = self.trackCustomObjects['ionicStrength'][scs]
      if issc not in self.sampleConds:
        self.sampleConds.append(issc)
        
    if hasattr(scs, 'sampleConditions'):
      for sc in scs.sortedSsampleConditions():
        if sc not in self.sampleConds:
          self.sampleConds.append(sc)

    return self.sampleConds

  def getSampleConditionSets(self,experiments):

    self.sampleCondSets = []

    for exp in experiments:
      if exp.sampleConditionSet and exp.sampleConditionSet not in self.sampleCondSets:
        self.sampleCondSets.append(exp.sampleConditionSet)

    return self.sampleCondSets

  def getSampleConditions(self,scs):

    self.sampleConds = []

    if hasattr(scs, 'sampleConditions'):
      for sc in scs.sortedSampleConditions():
        if sc not in self.sampleConds:
          self.sampleConds.append(sc)

    return self.sampleConds

  """
  def getSampleConditionSets(self,experiments):

    self.sampleCondSets = []

    for exp in experiments:
      if exp.sampleConditionSet not in self.sampleCondSets:
        self.sampleCondSets.append(exp.sampleConditionSet)

    return self.sampleCondSets
  """

  def getSoftwareTitle(self,software):

    if not software:
      return None

    return self.getTitle(software.name + '_' + software.version)

  def getVendorNames(self,software):

    self.vendorNames = []

    if software.vendorName and software.vendorName not in self.vendorNames:
      self.vendorNames.append(software.vendorName)

    return self.vendorNames

  def getSoftware(self,entry):

    software = []
    softwareNames = []

    for soft in entry.sortedSoftware():
      if soft.name not in softwareNames:
        softwareNames.append(soft.name)
        software.append(soft)

    return software

  def getMethods(self,entry):

    strucGens = entry.sortedStructureGenerations()

    self.methods = []

    for sg in strucGens:
      if sg.method not in self.methods and sg.method is not None and sg.method.name is not None:
        self.methods.append(sg.method)

    software = entry.sortedSoftware()

    for soft in software:
      for meth in soft.sortedMethods():
        if meth not in self.methods and meth is not None and meth.name is not None:
          self.methods.append(meth)

    return self.methods

  def getEntryMolecules(self,nmrEntry):

    self.entryMolecules = []

    for entryMol in nmrEntry.sortedEntryMolecules():
      molType = self.getApproxMolType(entryMol.molecule)
      if molType in ('protein', 'DNA', 'RNA', 'DNA/RNA', 'carbohydrate'):
        if entryMol not in self.entryMolecules:
          self.entryMolecules.append(entryMol)

    return self.entryMolecules

  def getStarConditionType(self,ccpnType):

    self.starType = ccpnType

    if ccpnType == 'ph':
      self.starType = 'pH'

    return self.starType

  def convertUnknown(self,value):

    # TODO - do we need to do this?

    #if value == 'unknown':
    #  value = None

    return value

  def getSpectrometers(self,nmrEntry):

    self.spectrometers = []

    #for spec in self.spectrometers:
    #  print 'SPEC A: [%s]' % spec

    for exp in nmrEntry.sortedExperiments():
      if exp.spectrometer and exp.spectrometer not in self.spectrometers:
        self.spectrometers.append(exp.spectrometer)

    #for spec in self.spectrometers:
    #  print 'SPEC B: [%s]' % spec

    keywds = {'application': 'nmrStar',
              'keyword':     'specFlag',
              'value':       False}

    appData = nmrEntry.findFirstApplicationData(**keywds)

    if appData:
      for instrument in nmrEntry.root.currentInstrumentStore.findAllInstruments(className = 'NmrSpectrometer'):

        if instrument.findFirstApplicationData(**keywds):
          if instrument not in self.spectrometers:
            self.spectrometers.append(instrument)

    #for spec in self.spectrometers:
    #  print 'SPEC C: [%s]' % spec

    return self.spectrometers

  def getNmrProbes(self,experiments):

    self.nmrProbes = []

    for exp in experiments:
      if exp.probe not in self.nmrProbes:
        self.nmrProbes.append(exp.probe)

    return self.nmrProbes

  def getStarExperimentName(self,experiment):

    bmrbName = None

    if experiment.refExperiment:
      refExperiment = experiment.refExperiment

      bmrbSysName = refExperiment.findFirstSystematicName(namingSystem = 'BMRB')

      if bmrbSysName:
        bmrbName = bmrbSysName.name

      if not bmrbName:
        bmrbName = refExperiment.name

        if not bmrbName and refExperiment.nmrExpPrototype:
          bmrbName = refExperiment.nmrExpPrototype.name

          if not bmrbName:
            bmrbName = refExperiment.nmrExpPrototype.synonym

      if bmrbName and not (bmrbSysName and bmrbSysName.name):
        bmrbName = experiment.name + ' (' + bmrbName + ')'

    if not bmrbName:
      bmrbName = experiment.name

    return bmrbName

  def getRawDataFlag(self,experiment):

    rawDataFlag = False

    if experiment.rawData:
      rawDataFlag = True

    #for ds in experiment.sortedDataSources():
    #  if ds.dataType in ('FID', 'part-processed'):
    #    rawDataFlag = True
    #    break

    return rawDataFlag

  def getStarExpState(self,experiment):

    starSampleState = 'isotropic'

    if experiment.sampleState: # Can be 'solid', 'liquid', 'ordered', 'powder', 'crystal'
      if experiment.sampleState == 'liquid':
        starSampleState = 'isotropic'
      else:
        starSampleState = experiment.sampleState # NMR_STAR can be - 'isotropic', 'anisotropic' or 'solid'

    else:
      keywds = {'application': 'nmrStar',
                'keyword':     'sampleState'}

      appData = experiment.findFirstApplicationData(**keywds)

      if appData:
        starSampleState = appData.value

    return starSampleState

  def getStarSampleState(self,sample):

    starSampleType = None

    keywds = {'application': 'nmrStar',
              'keyword':     'sampleState'}

    appData = sample.findFirstApplicationData(**keywds)

    if appData:
      starSampleType = appData.value

    else:
      starSampleType = 'solution'

      sampleState = None

      expt = sample.findFirstNmrExperiment()

      if expt:
        sampleState = expt.sampleState

      if sampleState == 'liquid':
        starSampleType = 'solution'

      elif sampleState: # Can be 'solid', 'liquid', 'ordered', 'powder', 'crystal'
        starSampleType = sampleState # NMR-STAR has many options starting with 'solution' and 'solid'

    return starSampleType

  def getStrucGenSubmittedNew(self,strEns):

    submitted = None

    if strEns and strEns.models:
      submitted = len(strEns.models)

    return submitted

  def getStrucGenCalculated(self,details):

    confCalc = None

    if not details:
      return confCalc

    searchObj = self.patt['confCalc'].search(details)

    if searchObj:
      confCalc = searchObj.group(1)

    if confCalc == 'NULL':
      confCalc = None

    return confCalc

  def getStrucGenSubmitted(self,details):

    confSub = None

    if not details:
      return confSub

    searchObj = self.patt['confSub'].search(details)

    if searchObj:
      confSub = searchObj.group(1)

    if confSub == 'NULL':
      confSub = None

    return confSub

  def getStrucGenSelected(self,details):

    confSel = None

    if not details:
      return confSel

    searchObj = self.patt['confSel'].search(details)

    if searchObj:
      confSel = searchObj.group(1)

    if confSel == 'NULL':
      confSel = None

    return confSel

  def getStrucGenRepresentative(self,details):

    confBest = None

    if not details:
      return confBest

    searchObj = self.patt['confBest'].search(details)

    if searchObj:
      confBest = searchObj.group(1)

    if confBest == 'NULL':
      confBest = None

    return confBest

  def getAllResonanceParents(self,nmrEntry):
  
    nmrProjectsAndConstraintStores = []
    
    for nmrLinks in ('derivedDataLists','experiments','measurementLists'):
      
      nmrLinkedObjects = getattr(nmrEntry,nmrLinks)
      
      for nmrLinkedObject in nmrLinkedObjects:
        if nmrLinkedObject.parent not in nmrProjectsAndConstraintStores:
          nmrProjectsAndConstraintStores.append(nmrLinkedObject.parent)
   
    # TODO CAN I REALLY HANDLE THESE? Also currently ResonanceID is a direct link to the loop, without the _list ID, so can't have multiple ones! 
    #for strucGen in nmrEntry.structureGenerations:
    # 
    #   nmrProjectsAndConstraintStores.append(strucGen.nmrConstraintStore)
        
    return nmrProjectsAndConstraintStores
  
  def getResonances(self,resonanceParent):
  
    if resonanceParent.className == 'NmrProject':
      resonances = resonanceParent.sortedResonances()
    else:
      resonances = resonanceParent.sortedFixedResonances()
      
    self.currentResonances = resonances

    return resonances
  
  def getOrderedResonances(self,constraintItem):
    
    resonances = list(constraintItem.orderedResonances)
    if not resonances:
      resonances = constraintItem.sortedResonances()
    return resonances
    
  def getResonanceSets(self):

    resonanceSets = []
    for resonance in self.currentResonances:
      resSet = resonance.resonanceSet
      if not resSet in resonanceSets:
        resonanceSets.append(resSet)
     
    return resonanceSets
    
  def getResonanceGroups(self):
    
    resonanceGroups = []
    for resonance in self.currentResonances:
      resGroup = resonance.resonanceGroup
      if not resGroup in resonanceGroups:
        resonanceGroups.append(resGroup)
     
    return resonanceGroups
    
  def getResonanceGroupCompID(self,resonanceGroup):
    
    comp_ID = None

    if not resonanceGroup:
      return
    
    if resonanceGroup.residue:
      comp_ID = findChemCompSysName('CIF',resonanceGroup.residue.chemCompVar.chemComp)
    elif resonanceGroup.chemComp:
      comp_ID = findChemCompSysName('CIF',resonanceGroup.chemComp)
    
    return comp_ID

  def getPeakListTitle(self,peakList):

    # TODO - something else here?

    return 'peak_list'

  def getPeakLists(self,nmrEntry):

    self.peakLists = []

    for nmrExp in nmrEntry.sortedExperiments():
      for dataSource in nmrExp.sortedDataSources():
        for peakList in dataSource.sortedPeakLists():
          if peakList not in self.peakLists:
            self.peakLists.append(peakList)

    return self.peakLists

  def getPeakMethods(self,peakList):

    methods = []

    if peakList.fitMethod:
      methods.append(peakList.fitMethod)

    if peakList.intensMethod and peakList.intensMethod not in methods:
      methods.append(peakList.intensMethod)

    return methods

  def getPeakContribs(self,peakList):

    self.peakContribs = []

    for peak in peakList.sortedPeaks():
      for peakContrib in peak.sortedPeakContribs():
        if peakContrib.weight > 0.0 and peakContrib not in self.peakContribs:
          self.peakContribs.append(peakContrib)

    return self.peakContribs

  def getMeasurementTitle(self,measurementList):

    className = measurementList.className

    title = ''

    if className == 'ShiftList' and measurementList.isSimulated == False:
      title = 'assigned_chem_shift_list'
    elif className == 'JCouplingList':
      title = 'coupling_constant_list'
    elif className == 'ShiftAnisotropyList':
      title = 'chem_shift_anisotropy_list' # ???
    elif className == 'ShiftList' and measurementList.isSimulated == True:
      title = 'theoretical_chem_shift_list' # ???
    elif className == 'RdcList':
      title = 'RDC_list'
    elif className == 'DipolarRelaxList':
      title = 'dipole_dipole_relax_list' # ???
    elif className == 'HExchRateList':
      title = 'H_exch_rate_list'
    elif className == 'HExchProtectionList':
      title = 'H_exch_protection_factor_list'
    elif className == 'NoeList':
      title = 'heteronuclear_noe_list'
    elif className == 'T1List':
      title = 'heteronuclear_T1_list'
    elif className == 'T1RhoList':
      title = 'heteronuclear_T1rho_list'
    elif className == 'T2List':
      title = 'heteronuclear_T2_list'

    return title

  def getProjectFileName(self,projName):

    return "From CCPN project: '%s'" % projName

  def getProtonRefFlag(self,chemShiftRefList):

    if 'chemShiftRef' in self.trackCustomObjects and chemShiftRefList in self.trackCustomObjects['chemShiftRef']:
      return self.trackCustomObjects['chemShiftRef'][chemShiftRefList][0]

  def getCarbonRefFlag(self,chemShiftRefList):

    if 'chemShiftRef' in self.trackCustomObjects and chemShiftRefList in self.trackCustomObjects['chemShiftRef']:
      return self.trackCustomObjects['chemShiftRef'][chemShiftRefList][1]

  def getNitrogenRefFlag(self,chemShiftRefList):

    if 'chemShiftRef' in self.trackCustomObjects and chemShiftRefList in self.trackCustomObjects['chemShiftRef']:
      return self.trackCustomObjects['chemShiftRef'][chemShiftRefList][2]

  def getPhosphorusRefFlag(self,chemShiftRefList):

    if 'chemShiftRef' in self.trackCustomObjects and chemShiftRefList in self.trackCustomObjects['chemShiftRef']:
      return self.trackCustomObjects['chemShiftRef'][chemShiftRefList][3]

  def getOtherRefFlag(self,chemShiftRefList):

    if 'chemShiftRef' in self.trackCustomObjects and chemShiftRefList in self.trackCustomObjects['chemShiftRef']:
      return self.trackCustomObjects['chemShiftRef'][chemShiftRefList][4]

  def getChemShiftRefTitle(self,chemShiftRefList):

    index = self.shiftReferenceList.index(chemShiftRefList) + 1
    chemShiftRefList.index = index
    chemShiftRefList.title = 'chemical_shift_reference_' + str(index)

    return chemShiftRefList.title

  def getChemShiftRefLists(self,nmrEntry):

    # TODO: Use trackCustomObjects dictionary here???

    if not hasattr(self, 'shiftReferenceList'):
      self.shiftReferenceList = []

      shiftRefsSetsSet = set()

      for nmrExpt in nmrEntry.sortedExperiments():
        shiftRefsSet = set()

        for shiftRef in nmrExpt.sortedShiftReferences():
          shiftRefsSet.add(shiftRef.serial)

        shiftRefsSetsSet.add(frozenset(shiftRefsSet) )

        shiftRefsSetsList = list(shiftRefsSetsSet)

        for shiftRefsSet1 in shiftRefsSetsList:
          for shiftRefsSet2 in shiftRefsSetsList:
            if len(shiftRefsSet1) > len(shiftRefsSet2):
              if shiftRefsSet2.issubset(shiftRefsSet1) and shiftRefsSet2 in shiftRefsSetsSet:
                shiftRefsSetsSet.remove(shiftRefsSet2)
            elif len(shiftRefsSet2) > len(shiftRefsSet1):
              if shiftRefsSet1.issubset(shiftRefsSet2) and shiftRefsSet1 in shiftRefsSetsSet:
                shiftRefsSetsSet.remove(shiftRefsSet1)

      #print 'SET SET: [%s]' % shiftRefsSetsSet

      shiftRefsIds = []

      for nmrExpt in nmrEntry.sortedExperiments():
        nmrProject = nmrExpt.parent

        shiftRefsSet = set()

        for shiftRef in nmrExpt.sortedShiftReferences():
          shiftRefsSet.add(shiftRef.serial)

        shiftRefsSetUse = None

        for shiftRefsSetFull in shiftRefsSetsSet:
          if shiftRefsSet.issubset(shiftRefsSetFull):
            shiftRefsSetUse = shiftRefsSetFull

        keywds = {'application': 'nmrStar',
                  'keyword':     'ShiftReferenceListId'}
        appData = nmrExpt.findFirstApplicationData(**keywds)

        shiftRefsIdStr = None

        if shiftRefsSetUse:
          if appData:
            shiftRefsIdStr = appData.value
            if shiftRefsSetUse and shiftRefsIdStr != str(shiftRefsSetUse):
              print '  Error: NMR Experiment %s has potentially the wrong chemical shift referencing.' % nmrExpt.name
              shiftRefsIdStr = str(shiftRefsSetUse)

          else:
            shiftRefsIdStr = str(shiftRefsSetUse)
            keywds['value'] = shiftRefsIdStr

            appData = Implementation.AppDataString(**keywds)
            nmrExpt.addApplicationData(appData)

          if shiftRefsIdStr and shiftRefsIdStr not in shiftRefsIds:
            shiftRefsIds.append(shiftRefsIdStr)

            if not self.trackCustomObjects.has_key('chemShiftRef'):
              self.trackCustomObjects['chemShiftRef'] = {}

            if self.trackCustomObjects['chemShiftRef'].has_key(shiftRefsSetUse):
              continue

            self.shiftReferenceList.append(ChemShiftRefList(Id=shiftRefsIdStr) )

            self.trackCustomObjects['chemShiftRef'][shiftRefsSetUse] = self.shiftReferenceList[-1]

            #print 'LIST: [%s]' % shiftRefsIdStr

            for Id in shiftRefsSetUse:
              shiftRef = nmrProject.findFirstShiftReference(serial=Id)
              if not shiftRef:
                print '  Error: NMR Experiment %s does not have chemical shift reference with serial %d' % (nmrExpt, Id)
              else:
                self.shiftReferenceList[-1].add(shiftRef)

    #for shiftRef in self.shiftReferenceList:
    #  print 'REF: [%s]' % shiftRef

    return self.shiftReferenceList

  def getChemShiftRefs(self,chemShiftRefList):

    chemShiftRefData = chemShiftRefList.chemShiftRefs

    protonFlag     = 'no'
    carbonFlag     = 'no'
    nitrogenFlag   = 'no'
    phosphorusFlag = 'no'
    otherFlag      = 'no'

    for chemShiftRef in chemShiftRefData:
      symbol = chemShiftRef.isotope.chemElement.symbol

      if symbol == 'H':
        if (chemShiftRef.isotope.massNumber == 1 and
            chemShiftRef.molName == 'DSS' and
            chemShiftRef.atomGroup == 'methyl protons' and
            chemShiftRef.unit == 'ppm' and
            self.getChemShiftRefValue(chemShiftRef.value) == 0.0 and
            self.getStarReferenceType(chemShiftRef) == 'internal' and
            chemShiftRef.referenceType == 'direct' and
            chemShiftRef.indirectShiftRatio == 1.0):
          protonFlag = 'yes with IUPAC referencing'
        else:
          protonFlag = 'yes'

      elif symbol == 'C':
        if (chemShiftRef.isotope.massNumber == 13 and
            chemShiftRef.molName == 'DSS' and
            chemShiftRef.atomGroup == 'methyl protons' and
            chemShiftRef.unit == 'ppm' and
            self.getChemShiftRefValue(chemShiftRef.value) == 0.0 and
            self.getStarReferenceType(chemShiftRef) == 'n/a' and
            chemShiftRef.referenceType == 'indirect' and
            chemShiftRef.indirectShiftRatio == 0.251449530):
          carbonFlag = 'yes with IUPAC referencing'
        else:
          carbonFlag = 'yes'

      elif symbol == 'N':
        if (chemShiftRef.isotope.massNumber == 15 and
            chemShiftRef.molName == 'DSS' and
            chemShiftRef.atomGroup == 'methyl protons' and
            chemShiftRef.unit == 'ppm' and
            self.getChemShiftRefValue(chemShiftRef.value) == 0.0 and
            self.getStarReferenceType(chemShiftRef) == 'n/a' and
            chemShiftRef.referenceType == 'indirect' and
            chemShiftRef.indirectShiftRatio == 0.101329118):
          nitrogenFlag = 'yes with IUPAC referencing'
        else:
          nitrogenFlag = 'yes'

      elif symbol == 'P':
        if (chemShiftRef.isotope.massNumber == 31 and
            chemShiftRef.molName == 'DSS' and
            chemShiftRef.atomGroup == 'methyl protons' and
            chemShiftRef.unit == 'ppm' and
            self.getChemShiftRefValue(chemShiftRef.value) == 0.0 and
            self.getStarReferenceType(chemShiftRef) == 'n/a' and
            chemShiftRef.referenceType == 'indirect' and
            chemShiftRef.indirectShiftRatio == 0.404808636):
          phosphorusFlag = 'yes with IUPAC referencing'
        else:
          phosphorusFlag = 'yes'

      elif symbol in ('2H', '29Si', '19F', '17O'):
        otherFlag = 'yes'

    if not self.trackCustomObjects.has_key('chemShiftRefFlags'):
      self.trackCustomObjects['chemShiftRefFlags'] = {}

    if not self.trackCustomObjects['chemShiftRefFlags'].has_key(chemShiftRefList):
      self.trackCustomObjects['chemShiftRef'][chemShiftRefList] = (protonFlag, carbonFlag, nitrogenFlag, phosphorusFlag, otherFlag)

    return chemShiftRefData

  def getChemShiftRefValue(self,value):

    if value == float(-999999999.0):
      value = None

    return value

  def getChemShiftRefId(self,chemShiftRef):

    return chemShiftRef.Id

  def getStarReferenceType(self,chemShiftRef):

    refType = chemShiftRef.className

    starRefType = 'n/a'

    if refType == 'InternalShiftReference':
      starRefType = 'internal'
    elif refType == 'ExternalShiftReference':

      keywds = {'application': 'nmrStar',
                'keyword':     'refType'}

      appData = chemShiftRef.findFirstApplicationData(**keywds)

      if appData:
        starRefType = appData.value

      else:
        starRefType = 'external'

    return starRefType

  def getChemShiftLocation(self,chemShiftRef):

    appData = chemShiftRef.findFirstApplicationData(application = 'nmrStar', keyword = 'location')

    location = None

    if hasattr(chemShiftRef, 'location') and chemShiftRef.location:
      location = chemShiftRef.location

    elif appData and appData.value:
      location = appData.value

    return location

  def getChemShiftGeometry(self,chemShiftRef):

    appData = chemShiftRef.findFirstApplicationData(application = 'nmrStar', keyword = 'geometry')

    geometry = None

    if hasattr(chemShiftRef, 'sampleGeometry') and chemShiftRef.sampleGeometry:
      geometry = chemShiftRef.sampleGeometry

    elif appData and appData.value:
      geometry = appData.value

    return geometry

  def getChemShiftAxis(self,chemShiftRef):

    appData = chemShiftRef.findFirstApplicationData(application = 'nmrStar', keyword = 'axis')

    axis = None

    if hasattr(chemShiftRef, 'axis') and chemShiftRef.axis:
      axis = chemShiftRef.axis

    elif appData and appData.value:
      axis = appData.value

    return axis

  def getShiftReferenceListId(self,shiftList):

    valueId = -1
    shiftIdx = None

    if shiftList:
      expt = shiftList.findFirstExperiment()
      if expt:
        appData = expt.findFirstApplicationData(application='nmrStar',keyword='ShiftReferenceListId')
        if appData:
          valueId = appData.value

    for shiftReference in self.shiftReferenceList:
      if valueId == shiftReference.Id:
        shiftIdx = shiftReference
        break

    #print 'IDX: [%s]' % shiftIdx

    return shiftIdx

  def getShiftReferenceListTitle(self,shiftList):

    valueId = -1

    if shiftList:
      expt = shiftList.findFirstExperiment()
      if expt:
        appData = expt.findFirstApplicationData(application='nmrStar',keyword='ShiftReferenceListId')
        if appData:
          valueId = appData.value

    title = None

    for shiftReference in self.shiftReferenceList:
      if valueId == shiftReference.Id:
        title = shiftReference.title
        break

    return title

  def getDataDimsFromDataSource(self,dataSource):
    
    dataDims = []
    for dataDim in dataSource.sortedDataDims():
      if dataDim.className == 'FreqDataDim':
	dataDims.append(dataDim)

    return dataDims

  def getFirstIsotopeName(self,isotopes):

    firstIsotopeName = ''

    if isotopes:
      firstIsotopeName = isotopes[0].chemElement.symbol

    return firstIsotopeName

  def getFirstIsotopeNumber(self,isotopes):

    firstIsotopeNumber = 0

    if isotopes:
      firstIsotopeNumber = isotopes[0].massNumber

    return firstIsotopeNumber

  def getFirstIsotopeCode(self,isotopeCodes):

    firstIsotopeCode = ''

    if isotopeCodes:
      firstIsotopeCode = isotopeCodes[0]

    return firstIsotopeCode

  def getAllIsotropicS2Derivations(self,nmrEntry):

    return self.getAllDerivedDataLists(nmrEntry,'IsotropicS2List')

  def getAllSpectralDensityDerivations(self,nmrEntry):

    return self.getAllDerivedDataLists(nmrEntry,'SpectralDensityList')

  def getAllPkaDerivations(self,nmrEntry):

    return self.getAllDerivedDataLists(nmrEntry,'PkaList')

  def getAllDataDerivations(self,nmrEntry):

    return self.getAllDerivedDataLists(nmrEntry,'DataList')

  def getAllDerivedDataLists(self,nmrEntry,classNames):

    derivedDataLists = []

    for className in classNames:
      dataLists = nmrEntry.findAllDerivedDataLists(className = className)
      for dataList in dataLists:
        derivedDataLists.extend(dataList.sortedDerivations() )

    return derivedDataLists

  def getShiftAmbiguityCode(self,shiftByIndividualAtom):

    shiftList = self.exportClass.ccpnVar['shiftList']
    nmrStarFormat = self.exportClass.ccpnVar['nmrStarFormatIndividual']
    
    chemShiftValue = shiftByIndividualAtom.measurement.value

    ambCode = nmrStarFormat.getShiftAmbiguityCode(chemShiftValue,shiftByIndividualAtom.resonanceToAtom,shiftList)
            
    return ambCode

  def getMeasurementsByIndividualAtom(self,measurements):
  
    if not hasattr(self,'measurementsByIndividualAtoms'):
      self.measurementsByIndividualAtoms = []
    
    measurementsByIndividualAtoms = []
    for measurement in measurements:
      
      # TODO this could be more than one?
      resonance = measurement.resonance

      if not self.exportClass.ccpnVar['nmrStarFormatIndividual'].resonanceToAtoms or resonance not in self.exportClass.ccpnVar['nmrStarFormatIndividual'].resonanceToAtoms:
        continue

      # Need to distinguish because in first loop nmrStarFormat is not yet available
      for resonanceToAtom in self.exportClass.ccpnVar['nmrStarFormatIndividual'].resonanceToAtoms[resonance]:      
        # See if exists already...
        foundMbia = False
        for mbia in self.measurementsByIndividualAtoms:
          if mbia.measurement == measurement and mbia.resonance == resonance and mbia.resonanceToAtom == resonanceToAtom:
            measurementsByIndividualAtoms.append(mbia)
            foundMbia = True
            break
       
        if not foundMbia:
          measurementsByIndividualAtoms.append(MeasurementByIndividualAtom(measurement,resonance,resonanceToAtom) )

    #measurementsByIndividualAtoms.sort(lambda a, b: cmp(a.resonance.serial, b.resonance.serial) )
    
    # sort by assignment
    ll = []
    for xx in measurementsByIndividualAtoms:
      r2a = xx.resonanceToAtom
      ll.append((r2a.chain.code, r2a.seqId, r2a.atomName, 
                 r2a.resonance.serial,xx))
    measurementsByIndividualAtoms = [tt[4] for tt in ll]
    
    # Add to general dict - should probably be divided by measurement list to speed things up if becomes important!      

    for mbia in measurementsByIndividualAtoms:
      if not mbia in self.measurementsByIndividualAtoms:
        self.measurementsByIndividualAtoms.append(mbia)
        
    return measurementsByIndividualAtoms

  def getStarNoeValueType(self,noeValueType):

    starNoeValueType = 'peak height'

    if noeValueType == 'height':
      starNoeValueType = 'peak height'

    elif noeValueType == 'volume':
      starNoeValueType = 'peak integral'

    elif noeValueType == 'contour count':
      starNoeValueType = 'contour count'

    # Also star tag of 'relative intensities'

    return starNoeValueType

  def getStarT1CoherenceType(self,t1List):

    t1CoherenceType = t1List.coherenceType

    t1 = t1List.findFirstMeasurement()

    isotope = t1.resonance.isotope
    if not isotope:
      return

    chemElement = isotope.chemElement
    if not chemElement:
      return

    symbol = chemElement.symbol

    starT1CoherenceType = None

    if t1CoherenceType == 'z':
      if symbol == 'C':
        starT1CoherenceType = 'Cz'
      elif symbol == 'N':
        starT1CoherenceType = 'Nz'
    elif t1CoherenceType == 'zz':
      if symbol == 'C':
        starT1CoherenceType = 'CzHz'  # Not in ADIT-NMR list
      elif symbol == 'N':
        starT1CoherenceType = 'NzHz'

    return starT1CoherenceType

  """
Star tags for main table (and table/tags for representative dihedral angle rmsd's)

_Constraint_stat_list
NOE_dist_averaging_method
NOE_tot_num
RDC_tot_num
Dihedral_angle_tot_num
Protein_dihedral_angle_tot_num (NA?)
NOE_intraresidue_tot_num (ROE?)
NOE_sequential_tot_num
NOE_medium_range_tot_num
NOE_long_range_tot_num
NOE_other_tot_num
RDC_NH_tot_num
RDC_other_tot_num
Protein_phi_angle_tot_num (NA?)
Protein_psi_angle_tot_num
Protein_chi_one_angle_tot_num
Protein_other_angle_tot_num
Protein_ambig_dihedral_tot_num
H_bonds_constrained_tot_num
Constr_def_H_bonds_tot_num
SS_bonds_constrained_tot_num
Constr_def_SS_bonds_tot_num

_Constraint_stat_list_rep
Dihedral_angle_rmsd
Dihedral_angle_rmsd_err
"""
  
  def getStarConstraintType(self,constraintList):

    # Other NMR-STAR types not used here:
    #'dipolar coupling'
    #'protein dihedral angle'
    #'nucleic acid dihedral angle'
    #'other angle'
    #'hydrogen exchange'
    #'line broadening'
    #'pseudocontact shift'
    #'intervector projection angle'
    #'protein peptide planarity'
    #'protein other kinds of constraints'
    #'nucleic acid base planarity'
    #'nucleic acid other kinds of constraints'

    starConstraintType = 'distance'

    if constraintList.className in ('DistanceConstraintList', 'HBondConstraintList'):
      starConstraintType = 'distance'

    elif constraintList.className == 'DihedralConstraintList':
      starConstraintType = 'dihedral angle' # TODO: ADIT-NMR distinguishes protein and nucleic acids - should we?

    elif constraintList.className == 'JCouplingConstraintList':
      starConstraintType = 'coupling constant'

    elif constraintList.className == 'RdcConstraintList':
      starConstraintType = 'rdc' # Could we use 'dipolar coupling' instead?

    elif constraintList.className == 'CsaConstraintList':
      starConstraintType = 'chemical shift anisotropy'

    elif constraintList.className == 'ChemShiftConstraintList':
      starConstraintType = 'chemical shift'

    return starConstraintType

  def getStarSubConstraintType(self,constraintList):

    # Other NMR-STAR subtypes not used here:
    #'Not applicable'
    #'NOE'
    #'NOE buildup'
    #'NOE not seen'
    #'alignment tensor'
    #'chirality'
    #'prochirality'
    #'disulfide bond'
    #'symmetry'
    #'ROE'
    #'peptide'
    #'ring'

    starSubConstraintType = 'Not applicable'

    if constraintList.className == 'DistanceConstraintList':
      starSubConstraintType = 'general distance'

    elif constraintList.className == 'HBondConstraintList':
      starSubConstraintType = 'hydrogen bond'

    return starSubConstraintType
  
  def getPdbCoordinateFileVersion(self,nmrEntry):
    
    appData = nmrEntry.findFirstApplicationData(application='nmrStar',keyword='pdbCoordVers')
    if appData:
      version = appData.value
    else:
      version = None
      
    return version

  def getStarSubSubConstraintType(self,constraintList):

    # Other NMR-STAR subsubtypes not used here:
    #'ambi'
    #'simple'

    starSubSubConstraintType = None

    return starSubSubConstraintType

  def getConstraintListCount(self,constraintList):

    return len(constraintList.constraints)

  def getAllConstraintLists(self,structureGenerations):

    constraintLists = self.getConstraintLists(structureGenerations,
                                              ['DistanceConstraintList','HBondConstraintList',
                                               'DihedralConstraintList','JCouplingConstraintList',
                                               'RdcConstraintList','CsaConstraintList'])
    chemShiftConst13C = self.getChemShiftConstraintLists(structureGenerations,isotopeCode='13C')

    if chemShiftConst13C:
      constraintLists.extend(chemShiftConst13C)

    chemShiftConst1H = self.getChemShiftConstraintLists(structureGenerations,isotopeCode='1H')

    if chemShiftConst1H:
      constraintLists.extend(chemShiftConst1H)

    return constraintLists

  def getDihedralConstraintLists(self,structureGenerations):

    return self.getConstraintLists(structureGenerations,['DihedralConstraintList'])

  def getDistanceConstraintLists(self,structureGenerations):

    return self.getConstraintLists(structureGenerations,['DistanceConstraintList','HBondConstraintList'])

  def getRdcConstraintLists(self,structureGenerations):

    return self.getConstraintLists(structureGenerations,['RdcConstraintList'])

  def getJCouplingConstraintLists(self,structureGenerations):

    return self.getConstraintLists(structureGenerations,['JCouplingConstraintList'])

  def getCACBChemShiftConstraintLists(self,structureGenerations):

    return self.getChemShiftConstraintLists(structureGenerations,isotopeCode='13C')

  def getHChemShiftConstraintLists(self,structureGenerations):

    return self.getChemShiftConstraintLists(structureGenerations,isotopeCode='1H')

  def getOtherConstraintLists(self,structureGenerations):

    # This CCPN constraint list is the only one not covered by NMR-Star.
    return self.getConstraintLists(structureGenerations,['CsaConstraintList'])

  def getChemShiftConstraintLists(self,structureGenerations, isotopeCode):

    constraintLists = []

    if structureGenerations:
      for strucGen in structureGenerations:
        if hasattr(strucGen, 'nmrConstraintStore') and strucGen.nmrConstraintStore:
          constraintLists.extend(list(strucGen.nmrConstraintStore.findAllConstraintLists(className = 'ChemShiftConstraintList', isotopeCode = isotopeCode) ) )

    return constraintLists

  def getConstraintLists(self,structureGenerations,classNames):

    constraintLists = []

    if structureGenerations:
      for strucGen in structureGenerations:
        for className in classNames:
          if hasattr(strucGen, 'nmrConstraintStore') and strucGen.nmrConstraintStore:
            #l = list(strucGen.nmrConstraintStore.findAllConstraintLists(className = className) )
            #if len(l):
            #  print 'NUMBER OF CONSTRAINTS: [' + str(len(l[0].constraints) ) + ']'
            constraintLists.extend(list(strucGen.nmrConstraintStore.findAllConstraintLists(className = className) ) )

    return constraintLists

  def getDistanceConstraintType(self,constraintList):
  
    constraintType = 'NOE'
    
    if constraintList.className == 'HBondConstraintList':
      constraintType = 'hydrogen bond'
    
    else:
      for experiment in constraintList.sortedExperiments():
        if experiment.refExperiment:
          if experiment.refExperiment.name.count('ROESY'):
            constraintType = 'ROE'

    return constraintType

  def getAllConstraintItems(self,constraintList):

    constraintItems = []

    for constraint in constraintList.sortedConstraints():
      for constraintItem in constraint.sortedItems():
        constraintItems.append(constraintItem)

    return constraintItems
    
  def setConstraintTreeNodeMemberId(self,constraintItem):
    
    orderedResonances = self.getOrderedResonances(constraintItem)
    index = orderedResonances.index(self.exportClass.ccpnVar['resonance'])
    
    return (index + 1) 

  def getGeneralConstraintLogic(self,constraint):
    
    constraintLogic = None
    
    if len(constraint.items) > 1:
      constraintLogic = 'OR'
      
    return constraintLogic

  def getFirstAddressLine(self,address):

    index = 0
    
    return self.getArrayEntryIndex(address,index)

  def getSecondAddressLine(self,address):

    index = 1
    
    return self.getArrayEntryIndex(address,index)

  def getThirdAddressLine(self,address):

    index = 2
    
    return self.getArrayEntryIndexToEnd(address,index)

  def getFirstPhoneNumber(self,phoneNums):

    index = 0
    
    return self.getArrayEntryIndex(phoneNums,index)

  def getArrayEntryIndex(self,array,index):

    if type(array) in (type([]),type( () ) ):
      if len(array) > index:
        return array[index]

  def getArrayEntryIndexToEnd(self,array,index):

    if type(array) in (type([]),type( () ) ):
      if len(array) > index:
        return ', '.join(array[index:])

  #################################
  # Common loop/table definitions #
  #################################

  def setKeywordLoop(self,ccpnLoop,ccpnMap):

    keywordDict = {

      'ccpnLoop': ccpnLoop,
      'ccpnMap':  ccpnMap,

      'tags': {
        'Keyword': ccpnMap
      }
    }

    return keywordDict

  def setFunctionLoop(self,ccpnLoop,ccpnMap):

    functionDict = {

      'ccpnLoop': ccpnLoop,
      'ccpnMap':  ccpnMap,

      'tags': {
        'Biological_function': ccpnMap
      }
    }

    return functionDict

  def setCommonNameLoop(self,ccpnLoop,ccpnMap):

    commonNameDict = {

      'ccpnLoop': ccpnLoop,
      'ccpnMap':  ccpnMap,

      'tags': {
        'Name': ccpnMap,
        'Type=LOCAL': 'name'
      }
    }
 
    return commonNameDict

  def setSysNamesLoop(self,ccpnLoop,ccpnMap):

    sysNameDict = {

      'ccpnLoop': ccpnLoop,
      'ccpnMap':  ccpnMap,

      'tags': {
        'Name':          ccpnMap + '.name',
        'Naming_system': ccpnMap + '.namingSystem'
      }
    }

    return sysNameDict

  def setSoftwareLoop(self,ccpnMap):

    software = ccpnMap + '.method.software'
    #softwareName = software + '.name' + '_' + software + '.version'

    method = ccpnMap + '.method'
    methodName = method + '.name'

    softwareDict = {

      'ccpnMap': ccpnMap,

      'tags': {

        'Software_ID':     software,
        #'Software_label':  softwareName,
        'Software_label': (software, self.getSoftwareTitle),
        'Method_ID':       method,
        'Method_label':    methodName,

      }
    }

    return softwareDict

  def setExperimentsLoop(self,ccpnLoop,ccpnMap):

    sample = ccpnMap + '.sample'
    sampleName = sample + '.name'

    experimentDict = {

      'ccpnLoop': ccpnLoop,
      'ccpnMap':  ccpnMap,

      'tags': {

#        'Experiment_ID':   ccpnMap + '.serial',
        'Experiment_name': (ccpnMap, self.getStarExperimentName), # ccpnMap + '.name',
        'Sample_ID':       sample,
        'Sample_label':    sampleName,
        'Sample_state':    (ccpnMap, self.getStarExpState), #ccpnMap + '.sampleState',

      }
    }

    return experimentDict

  def setConExperimentsLoop(self,ccpnLoop,ccpnMap):

    sample = ccpnMap + '.sample'
    sampleName = sample + '.name'

    method = ccpnMap + '.derivationMethod'
    methodName = method + '.name'

    experimentDict = {

      'ccpnLoop': ccpnLoop,
      'ccpnMap':  ccpnMap,

      'tags': {

#        'Experiment_ID':   ccpnMap + '.serial',
        'Experiment_name': (ccpnMap, self.getStarExperimentName), # ccpnMap + '.name',
        'Sample_ID':       sample,
        'Sample_label':    sampleName,
        'Sample_state':    ccpnMap + '.sampleState',
        'Method_ID':       method,
        'Method_label':    methodName,

      }
    }

    return experimentDict

  def setPdbDepositionSite(self,nmrEntry):
    
    depSite = None
    
    appData = nmrEntry.findFirstApplicationData(application='nmrStar',keyword='pdbDepSite')
    
    if appData and appData.value:
      depSite = appData.value
      
    return depSite

  def modifySfDict(self):
  
    # Necessary if want to change specific mapping issue in a subclass.
    pass

  def setSfDict(self):

    self.sfDict = {}

    #######################
    # General information #
    #######################

    """
    self.sfDict['study_list'] = {                                                                             # *** Priority 1 - Not in ADIT-NMR?

      'ccpnMap':                         'study',

      'title':                          ('study.name',self.getTitle),

      'tags': {

        #'Sf_category':                  ['study_list',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],

        },

      'tables': {

        'Study': {                                                                                            # *** Table
        
          'ccpnMap':                     'study',

          'tags': {

            'Name':                      'study.name',
            'Type':                      'study.studyType',
            # TODO studyType needs to be set to one of these:
            #   'Structure analysis', 'Mutant comparison', 'Interactions with different ligands'
            'Details':                   'study.details'

            },
          },

        'Study_keyword':                 self.setKeywordLoop('study.keywords','keyword'),                     # *** Table

        'Study_entry_list': {                                                                                 # *** Table

          'tags': {

            #'Study_ID':                 [None,returnStarInt,'Study.ID',True],
            #'BMRB_accession_code':      [None,lambda x = value: returnStarCode(x,length = 12),None,True],
            #'BMRB_entry_description':   [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Details':                  [None,returnStarString,None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Study_list_ID':            [None,returnStarInt,'Study_list.ID',True],

            },
          },

        }
      }
    """

    self.sfDict['entry_information'] = {                                                                      # *** Priority 1

      'ccpnMap':                         '',

      'title':                          ('nmrEntry.name',self.getTitle),

      'tags': {

        'Title':                         'nmrEntry.title', # TODO: check mapping here!!                       # *** Mandatory
        'NMR_STAR_version=LOCAL':        self.nmrStarVersion,
        'Experimental_method=LOCAL':     'NMR',
        #'Experimental_method':           "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='exptMethod').value",
        'Experimental_method_subtype':   'nmrEntry.entryType',                                                # *** Mandatory - pulldown - auto
        'PDB_coordinate_file_version':   ('nmrEntry',self.getPdbCoordinateFileVersion),
        'Details':                       'nmrEntry.details',                                                  # *** Optional

        'Special_processing_instructions': 'nmrEntry.bmrbProcessing',                                         # *** Optional

        'Version_type':                  "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='versionType').value",
        'Submission_date':               "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='subDate').value",
        'Accession_date':                "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='accDate').value",
        'Origination':                   "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='origin').value",
        #'NMR_STAR_version':              "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='nmrStarVers').value",
        'Original_NMR_STAR_version':     "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='origNmrStarVers').value",
        'Dep_release_code_coordinates':  "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='depRelCoord').value", # *** Mandatory (PDB) - pulldown - auto
        'Dep_release_code_nmr_constraints': "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='depRelConstr').value", # *** Mandatory (PDB) - pulldown - auto
        'Dep_release_code_nmr_exptl':    "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='depRelNmr').value", # *** Mandatory - pulldown - auto
        'Dep_release_code_sequence':     "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='depRelSeq').value", # *** Mandatory - pulldown - auto
        'CASP_target':                   "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='caspTarget').value", # *** Optional (PDB))
        'Update_BMRB_accession_code':    "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='updBmrbCode').value", # *** Optional
        'Replace_BMRB_accession_code':   "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='replBmrbCode').value", # *** Optional
        'Update_PDB_accession_code':     "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='updPdbCode').value",
        'Replace_PDB_accession_code':    "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='replPdbCode').value", # *** Optional (PDB)
        'BMRB_update_details':           "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='bmrbUpdDet').value", # *** Optional
        'PDB_update_details':            "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='pdbUpdDet').value", # *** Optional (PDB)
        'Release_request':               "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='relReq').value",
        'PDB_deposit_site':             ('nmrEntry',self.setPdbDepositionSite),
        'PDB_process_site':              "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='pdbProcSite').value",
        'BMRB_deposit_site':             "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='bmrbDepSite').value",
        'BMRB_process_site':             "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='bmrbProcSite').value",
        'RCSB_annotator':                "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='rcsbAnn').value",
        'Assigned_BMRB_ID':              "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='assBmrbId').value",
        'Assigned_BMRB_deposition_code': "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='assBmrbCode').value",
        'Assigned_PDB_ID':               "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='assPdbId').value",
        'Assigned_PDB_deposition_code':  "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='assPdbCode').value",
        'Date_nmr_constraints':          "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='dateNmrConstr').value",
        'Recvd_author_approval':         "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='recAuthAppr').value",
        'PDB_date_submitted':            "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='pdbDateRec').value",
        'Author_release_status_code':    "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='authRelStatCode').value",
        'Author_approval_type':          "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='authApprType').value",

        },

      'tables': {

        'Entry_proc_cycle': {                                                                                 # *** Table

          'tags': {

            #'Cycle_ID':                 [None,returnStarInt,None,True],
            #'Date_begin_cycle':         [None,returnStarDateTime,None,False],
            #'Date_end_cycle':           [None,returnStarDateTime,None,False],
            #'Details':                  [None,returnStarString,None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

            },
          },

        'Entry_prerelease_seq': {                                                                             # *** Table

          'tags': {

            #'Entity_ID':                [None,returnStarInt,'Entity.ID',True],
            #'Entity_label':             [None,lambda x = value: returnStarLine(x,length = 127),'Entity.Sf_framecode',False],
            #'Seq_one_letter_code':      [None,returnStarString,None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],

            },
          },

        # TODO: This not set for RECOORD - need entries that have this info! Autodep?

        'Contact_person': {                                                                                   # *** Includes mandatory fields

          'ccpnLoop':                    'nmrEntry.sortedContactPersons()',
          'ccpnMap':                     'contactPerson',

          'tags': {

            'Email_address':             'contactPerson.currentPersonInGroup.emailAddress',
            #'Email_address=LOCAL':       'penkett@ebi.ac.uk', # TODO - change!!!                              # *** Mandatory

            'Name_salutation':          ('contactPerson.title', self.getStarTitle),                           # *** Optional - option
            'Given_name':                'contactPerson.givenName',                                           # *** Mandatory
            'Family_name':               'contactPerson.familyName',                                          # *** Mandatory

            #'First_initial':            'contactPerson.firstInitial', # This doesn't exist for contact person?

            'Middle_initials':          ('contactPerson.middleInitials',self.getMiddleInitials),
            'Family_title':              'contactPerson.familyTitle',                                         # *** Optional - option

            'Department_and_institution': ('contactPerson',self.getDepartmentAndInstitution),
            'Mailing_address':           'contactPerson.currentPersonInGroup.mailingAddress', # TODO - set as organisation address info if no mailingAddress
            'Address_1':                ('contactPerson.currentPersonInGroup.group.organisation.addresses',self.getFirstAddressLine), # *** Mandatory
            'Address_2':                ('contactPerson.currentPersonInGroup.group.organisation.addresses',self.getSecondAddressLine), # *** Optional
            'Address_3':                ('contactPerson.currentPersonInGroup.group.organisation.addresses',self.getThirdAddressLine), # *** Optional
            'City':                      'contactPerson.currentPersonInGroup.group.organisation.city',    # *** Mandatory
            'State_province':           ('contactPerson.currentPersonInGroup.group.organisation', self.getStateProv),        # *** Mandatory
            'Country':                   'contactPerson.currentPersonInGroup.group.organisation.country', # *** Mandatory
            'Postal_code':               'contactPerson.currentPersonInGroup.group.organisation.postalCode', # *** Mandatory
            'Phone_number':             ('contactPerson.currentPersonInGroup.phoneNumbers',self.getFirstPhoneNumber), # *** Mandatory
            'FAX_number':                'contactPerson.currentPersonInGroup.faxNumber',                  # *** Optional
            'Role':                      'contactPerson.currentPersonInGroup.position',                   # *** Mandatory - pulldown - auto (Must have at least one principal investigator)
            'Organization_type':         'contactPerson.currentPersonInGroup.group.organisation.organisationType', # *** Mandatory - pulldown - auto

            # Organisation info - could be used as backup:
            #   'phoneNumber', 'faxNumber', 'emailAddress'

            },
          },
              
        'Entry_author': {                                                                                     # *** Includes mandatory fields

          'ccpnLoop':                    'nmrEntry.authors',
          'ccpnMap':                     'author',

          'tags': {

            #'Ordinal':                  [None,returnStarInt,None,True],

            'Given_name':                'author.givenName',                                                  # *** Mandatory
            'Family_name':               'author.familyName',                                                 # *** Mandatory
            #'First_initial':            [None,lambda x = value: returnStarLine(x,length = 15),None,False], # 'author.firstInitial',
            'Middle_initials':          ('author.middleInitials',self.getMiddleInitials),                     # *** Optional
            'Family_title':              'author.familyTitle',                                                # *** Optional - option

            },
          },

        #'SG_project': {                                                                                       # *** Includes mandatory (PDB) fields - also mandatory if selected in entry_interview

        #  'ccpnLoop':                    None,
        #  'ccpnMap':                     None,

        #  'tags': {

        #    'Project_name':             "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='projName').value",
        #    'Full_name_of_center':      "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='centerFullName').value",

        #    },
        #  },

        'Entry_src': {                                                                                        # *** Includes optional fields

          #'ccpnLoop':                   ('nmrEntry', self.getContactOrganisations),
          #'ccpnMap':                     'organisation',

          'ccpnLoop':                   ('nmrEntry', self.getLaboratories),
          'ccpnMap':                     'laboratory',

          'tags': {
            'Project_name=LOCAL':        '.', # TODO - can we force this another way?
            'Organization_initials=LOCAL': '.', # TODO - can we force this another way?
            'Organization_full_name':    'laboratory.name',                                                   # *** Optional - option
            #'Organization_full_name':    'organisation.name',                                                # *** Optional - option

            },
          },

        #'Struct_keywords':               self.setKeywordLoop('study.keywords','keyword'),                    # *** Table

        'Struct_keywords': {                                                                                  # *** Includes mandatory (PDB) fields

          'ccpnLoop':                    'nmrEntry.keywords',
          'ccpnMap':                     'keyword',

          'tags': {

            'Keywords':                  'keyword',                                                           # *** Mandatory (PDB)

            },
          },

        'Data_set': {                                                                                         # *** Table

          'ccpnLoop':                   ('nmrEntry', self.getDataListCounts),
          'ccpnMap':                     'dataListCount',

          'tags': {

            'Type':                      'dataListCount.name',
            'Count':                     'dataListCount.count',
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],

            },
          },

        'Datum': {                                                                                            # *** Table

          'ccpnLoop':                   ('nmrEntry', self.getDataCounts),
          'ccpnMap':                     'dataCount',

          'tags': {

            'Type':                      'dataCount.name',
            'Count':                     'dataCount.count',
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],

            },
          },

        #'Release': {                                                                                          # *** Table

        #  'tags': {

        #  'ccpnLoop':                    None,
        #  'ccpnMap':                     None,

        #    'Release_number':            "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='relNo').value",
        #    'Date':                      "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='date').value",
        #    'Submission_date':           "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='subDateRel').value",
        #    'Type':                      "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='relType').value",
        #    'Author':                    "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='author').value",
        #    'Detail':                    "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='detail').value",
        #    #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],

        #    },
        #  },

        'Related_entries': {                                                                                  # *** Includes optional fields

          'ccpnLoop':                    'nmrEntry.sortedRelatedEntries()',
          'ccpnMap':                     'relatedEntry',

          'tags': {

            #'Database_name':             ('relatedEntry.bmrbId',self.getRelatedDbName),                      # *** Optional - pulldown
            'Database_name':             'relatedEntry.dbName',                                               # *** Optional - pulldown
            'Database_accession_code':   'relatedEntry.dbCode',                                               # *** Optional

            'Relationship':              'relatedEntry.relationship',                                         # *** Optional

            },
          },

        }
      }


    #############
    # Citations #
    #############

    self.sfDict['citations'] = {                                                                              # *** Priority 1

      'ccpnLoop':                       ('nmrEntry',self.getCitations),
      'ccpnMap':                         'citation',

      'title':                          ('citation',self.getCitationStarClassLabel),

      'tags': {

        #'Sf_category':                  ['citations',returnStarString,None,True],
        #'Sf_framecode':                 [None,returnStarString,None,False],                                  # *** Mandatory - set in NmrStarExport - auto?
        #'Entry_ID':                     [None,returnStarString,'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Class':                        ('citation',self.getCitationStarClass),                               # *** Mandatory - pulldown - auto
        'CAS_abstract_code':             'citation.casAbstractCode',
        'MEDLINE_UI_code':               'citation.medlineUiCode',
        'PubMed_ID':                     'citation.pubMedId',                                                 # *** Optional
        'DOI':                           'citation.doi',                                                      # *** Optional
        #'Full_citation':                [None,returnStarString,None,False],
        'Title':                         'citation.title',                                                    # *** Mandatory
        'Status':                        'citation.status',                                                   # *** Mandatory - pulldown - auto
        'Type':                         ('citation.className',self.getCitationStarType),                      # *** Mandatory - pulldown - auto
        'Page_first':                    'citation.firstPage',                                                # *** Optional
        'Page_last':                     'citation.lastPage',                                                 # *** Optional
        'Year':                          'citation.year',                                                     # *** Optional
        'Details':                       'citation.details',
        #'WWW_URL':                      [None,returnStarString,None,False],

        'CONDITIONAL': { # NOTE HERE: can also use 'DEFAULT' value, which is the 'standard' option if all else fails

          ('citation.className',self.getCitationStarType): { # TODO - book chapter, personal communication (no extra fields), internet (WWW_URL - Optional), BMRB only (no extra fields)

            'journal': {

              'Journal_abbrev':          'citation.journalAbbreviation',                                      # *** Mandatory - option
              'Journal_name_full':       'citation.journalFullName',
              'Journal_volume':          'citation.volume',                                                   # *** Optional
              'Journal_issue':           'citation.issue',                                                    # *** Optional
              'Journal_ASTM':            'citation.astm',
              'Journal_ISSN':            'citation.issn',
              'Journal_CSD':             'citation.csd'

              },

            'book': { # 'book chapter' has same fields

              'Book_title':              'citation.bookTitle',                                                # *** Mandatory
              #'Book_chapter_title':     [None,returnStarString,None,False],                                  # *** Mandatory
              'Book_volume':             'citation.volume',                                                   # *** Optional
              'Book_series':             'citation.bookSeries',                                               # *** Optional
              'Book_publisher':          'citation.publisher',                                                # *** Mandatory
              'Book_publisher_city':     'citation.publisherCity',                                            # *** Mandatory
              'Book_ISBN':               'citation.isbn'

              },

            'thesis': {

              'Thesis_institution':      'citation.institution',                                              # *** Mandatory
              'Thesis_institution_city': 'citation.city',                                                     # *** Mandatory
              'Thesis_institution_country': 'citation.country'                                                # *** Mandatory

              },

            'abstract': {

              'Conference_title':        'citation.conferenceTitle',                                          # *** Mandatory
              'Conference_site':         'citation.city',                                                     # *** Mandatory
              'Conference_state_province': 'citation.stateProvince',                                          # *** Mandatory
              'Conference_country':      'citation.country',                                                  # *** Mandatory
              'Conference_start_date':   'citation.startDate',                                                # *** Mandatory
              'Conference_end_date':     'citation.endDate',                                                  # *** Mandatory
              'Conference_abstract_number': 'citation.abstractNumber'                                         # *** Mandatory

              },
            },

          },
        },

      'tables': {

        'Citation_author': {                                                                                  # *** Includes mandatory fields

          'ccpnLoop':                   ('citation',self.getCitationAuthors),
#          'ccpnLoop':                    'citation.authors',
          'ccpnMap':                     'author',

          'tags': {

            'Ordinal':                  ('author',self.findAuthorIndex),
            #'Ordinal':                   'author.serial',
            'Given_name':                'author.givenName',                                                  # *** Mandatory
            'Family_name':               'author.familyName',                                                 # *** Mandatory
            #'First_initial':            [None,lambda x = value: returnStarLine(x,length = 15),None,False], # 'author.firstInitial',
            'Middle_initials':          ('author.middleInitials',self.getMiddleInitials),                     # *** Optional
            'Family_title':              'author.familyTitle',                                                # *** Optional - option
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Citation_ID':              [None,returnStarInt,'Citation.ID',True],

            },
          },

        'Citation_keyword':              self.setKeywordLoop('citation.keywords','keyword'),                  # *** Table

        #'Citation_keyword': {                                                                                # *** Includes optional fields

        #  'tags': {

            #'Keyword':                  [None,returnStarString,None,True],                                   # *** Optional
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Citation_ID':              [None,returnStarInt,'Citation.ID',True],

        #    },
        #  },

        'Citation_editor': {                                                                                  # *** Table

          'ccpnLoop':                    'citation.editors',
          'ccpnMap':                     'editor',

          'tags': {

            #'Ordinal':                  [None,returnStarInt,None,True],
            'Given_name':                'editor.givenName',
            'Family_name':               'editor.familyName',
            #'First_initial':            [None,lambda x = value: returnStarLine(x,length = 15),None,False], # 'editor.firstInitial',
            'Middle_initials':          ('editor.middleInitials',self.getMiddleInitials),
            'Family_title':              'editor.familyTitle',
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Citation_ID':              [None,returnStarInt,'Citation.ID',True],

            },
          },

        }
      }


    ########################
    # Molecule information #
    ########################

    self.sfDict['assembly'] = {                                                                               # *** Priority 1

      'ccpnMap':                         'molSystem',

#      'title':                          ('molSystem.code',self.getTitle),
      'title':                          ('molSystem.code',self.getAssemblyTitle),

      'tags': {

        'Name':                         ('molSystem',self.getMolSystemName),                                  # *** Mandatory
        #'BMRB_ID':                      [None,returnStarString,None,None],
        'Details':                       'molSystem.details',                                                 # *** Optional
        'Number_of_components':         ('molSystem',self.getNumChains),                                      # *** Optional
        'Organic_ligands':              ('molSystem.chains',self.getOrganicLigands),
        'Metal_ions':                   ('molSystem.chains',self.getMetalIons),
        #'Ambiguous_conformational_states': [None,returnStarString,None,None],
        #'Ambiguous_chem_comp_sites':    [None,returnStarString,None,None],
        #'Molecules_in_chemical_exchange': [None,returnStarString,None,None],
        'Paramagnetic':                  'molSystem.isParamagnetic',
        'Thiol_state':                  ('molSystem.chains',self.getThiolState),  # ???
        'Molecular_mass':                'molSystem.molecularMass',                                           # *** Optional
        'Enzyme_commission_number':     ('molSystem',self.getMolSysECNumber),                                 # *** Optional
        'Non_standard_bonds':           ('molSystem',self.getNonStdBondBool),                                 # *** Optional - bool

        },

      'tables': {

        'Assembly_common_name': {                                                                             # *** Table

          'ccpnLoop':                    'molSystem.commonNames',
          'ccpnMap':                     'commonName', # NOTE: this works but is not superelegant!

          'tags': {

            'Name':                      'commonName'

            },
          },

        'Assembly_systematic_name':      self.setSysNamesLoop('molSystem.molSystemSysNames','molSystemSysName'), # *** Table

        'Assembly_type': {                                                                                    # *** Table

          'tags': {

            #'Type':                     [None,returnStarString,None,None],

            },
          },

        'Entity_assembly': {                                                                                  # *** Includes mandatory fields

          #'category':                   'sequence',

          # TODO: Might have to use 'export' chain codes (and seqCodes) here? Or maybe just forget it...
          'ccpnLoop':                   ('molSystem',self.getMolSystemChains),
          'ccpnMap':                     'chain',

          'tags': { # TODO TODO: can multiple CCPN molecules correspond to ONE molecule for BMRB!?!? Does it matter?

            'Entity_ID':                 'chain.molecule',   # TODO; need to force setting this!
            'Entity_assembly_name':     ('chain.molecule.name',self.getName),                                 # *** Mandatory
              #'Entity_label':          ('chain.molecule.name',self.getLabel), # Or should be different structure?
            'Entity_label':              'chain.molecule.name', # Or should be different structure?           # *** Mandatory
            'Asym_ID':                   'chain.code',

            #'Entity_label':             'chain.code', # Or this?
            'Magnetic_equivalence_group_code': 'chain.magnEquivalenceCode',                                   # *** Optional
            'Role':                      'chain.role',                                                        # *** Optional

            'Details':                   'chain.molecule.details',                                            # *** Optional

            'Physical_state':            'chain.physicalState',                                               # *** Mandatory - option - auto
            'Conformational_isomer':    ('chain.conformationalIsomer',self.checkType),                        # *** Mandatory - bool - auto
            'Chemical_exchange_state':  ('chain.chemExchangeState',self.checkType),                           # *** Mandatory - bool - auto

            'Experimental_data_reported': ('chain',self.getExptDataReported),                                 # *** Mandatory - bool

            },
          },

        'Assembly_interaction': {                                                                             # *** Includes optional fields

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],                                      # *** Optional
            #'Entity_assembly_ID_1':     [None,returnStarInt,'Entity_assembly.ID',True],                      # *** Optional
            #'Entity_assembly_ID_2':     [None,returnStarInt,'Entity_assembly.ID',True],                      # *** Optional
            #'Mol_interaction_type':     [None,lambda x = value: returnStarLine(x,length = 127),None,True],   # *** Optional - option
          
            },
          },

        'Chem_comp_assembly': {                                                                               # *** Table

          'tags': {

            #'Entity_assembly_ID':       [None,returnStarString,'Entity_assembly.ID',None],
            #'Entity_ID':                [None,returnStarString,'Entity.ID',None],
            #'Entity_label':             [None,returnStarString,None,None],
            #'Comp_ID':                  [None,returnStarString,'Chem_comp.ID',None],
            #'Chem_comp_label':          [None,returnStarString,None,None],
            #'Comp_index_ID':            [None,returnStarString,'Entity_comp_index.ID',None],
            #'Application_seq_num':      [None,returnStarString,None,None],
            #'Seq_ID':                   [None,returnStarInt,'Entity_poly_seq.Num',None],

            },
          },

        'Atom': {                                                                                             # *** Table

          'tags': {

            #'Entry_atom_ID':            [None,returnStarString,None,None],
            #'Entity_assembly_ID':       [None,returnStarString,'Entity_assembly.ID',None],
            #'Entity_ID':                [None,returnStarString,'Entity.ID',None],
            #'Entity_label':             [None,returnStarString,None,None],
            #'Comp_index_ID':            [None,returnStarString,'Entity_comp_index.ID',None],
            #'Seq_ID':                   [None,returnStarString,'Entity_poly_seq.Num',None],
            #'Comp_ID':                  [None,returnStarString,'Chem_comp.ID',None],
            #'Chem_comp_label':          [None,returnStarString,None,None],
            #'Atom_ID':                  [None,returnStarString,None,None],
            #'Type_symbol':              [None,returnStarString,None,None],
            #'PDB_one_letter_code':      [None,returnStarString,None,None],
            #'PDB_chain_code':           [None,returnStarString,None,None],
            #'PDB_insertion_code':       [None,returnStarString,None,None],
            #'PDB_seq_ID':               [None,returnStarInt,'Entity_poly_seq.Num',None],
            #'PDB_comp_ID':              [None,returnStarString,'Chem_comp.ID',None],
            #'PDB_group':                [None,returnStarString,None,None],
            #'PDB_atom_name':            [None,returnStarString,None,None],
            #'PDB_atom_type':            [None,returnStarString,None,None],

            },
          },

        'Bond': {                                                                                             # *** Mandatory if cross linking/intermolecular bonds are present

          'ccpnLoop':                   ('molSystem',self.getNonStdBonds),
          'ccpnMap':                     'nonStdBond',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],                                      # *** Mandatory
            'Type':                     ('nonStdBond',self.getBondType),                                      # *** Mandatory
            'Order':                    ('nonStdBond',self.getBondOrder),                                     # *** Mandatory
            'Assembly_atom_ID_1':        None, #[None,returnStarInt,'Atom.Assembly_atom_ID',False],
            'Entity_assembly_ID_1':      None, #'nonStdBond.molecule.findFirstChain()', # How do we deal with multiple chains?
            'Entity_assembly_name_1':    None, #'nonStdBond.molecule.name',
            'Entity_ID_1':               None, #'nonStdBond.molecule',
            'Comp_ID_1':                 None, #('nonStdBond',self.getFirstBondMolResCode),
            'Comp_index_ID_1':           None, #('nonStdBond',self.getFirstBondMolResidue),
            'Seq_ID_1':                  None, #('nonStdBond',self.getFirstBondMolResID),
            'Atom_ID_1':                 None, #('nonStdBond',self.getFirstBondLinkCode),
            'Assembly_atom_ID_2':        None, #[None,returnStarInt,'Atom.Assembly_atom_ID',False],
            'Entity_assembly_ID_2':      None, #'nonStdBond.molecule.findFirstChain()', # How do we deal with multiple chains?
            'Entity_assembly_name_2':    None, #'nonStdBond.molecule.name',
            'Entity_ID_2':               None, #'nonStdBond.molecule',
            'Comp_ID_2':                 None, #('nonStdBond',self.getSecondBondMolResCode),
            'Comp_index_ID_2':           None, #('nonStdBond',self.getSecondBondMolResidue),
            'Seq_ID_2':                  None, #('nonStdBond',self.getSecondBondMolResID),
            'Atom_ID_2':                 None, #('nonStdBond',self.getSecondBondLinkCode),
            #'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            'Auth_entity_assembly_name_1': 'nonStdBond.molecule.name',
            'Auth_seq_ID_1':            ('nonStdBond',self.getFirstBondMolResID),                             # *** Mandatory
            'Auth_comp_ID_1':           ('nonStdBond',self.getFirstBondMolResCode),                           # *** Mandatory
            'Auth_atom_ID_1':           ('nonStdBond',self.getFirstBondLinkCode),                             # *** Mandatory
            #'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            'Auth_entity_assembly_name_2': 'nonStdBond.molecule.name',                                        # *** Mandatory
            'Auth_seq_ID_2':            ('nonStdBond',self.getSecondBondMolResID),                            # *** Mandatory
            'Auth_comp_ID_2':           ('nonStdBond',self.getSecondBondMolResCode),                          # *** Mandatory
            'Auth_atom_ID_2':           ('nonStdBond',self.getSecondBondLinkCode),                            # *** Mandatory
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Assembly_ID':              [None,returnStarInt,'Assembly.ID',True],

            },
          },

        'Deleted_atom': {                                                                                     # *** Mandatory if cross linking/intermolecular bonds are present

          'ccpnLoop':                   ('molSystem',self.getDeletedAtoms),
          'ccpnMap':                     'deletedAtom',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            'Assembly_atom_ID':          None, #[None,returnStarInt,'Atom.Assembly_atom_ID',False],
            'Entity_assembly_ID':        None, #[None,returnStarInt,'Entity_assembly.ID',True],
            'Entity_assembly_name':      'deletedAtom.molResidue.parent.name',                                # *** Mandatory
            'Entity_ID':                 None, #[None,returnStarInt,'Entity.ID',True],
            'Entity_label':              None, #[None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',True],
            'Comp_index_ID':             None, #[None,returnStarInt,'Entity_comp_index.ID',True],
            'Seq_ID':                    None, #[None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            'Comp_ID':                   None, #[None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            'Comp_label':                None, #[None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            'Atom_ID':                   None, #[None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Auth_entity_assembly_ID':  [None,lambda x = value: returnStarCode(x,length = 12),None,False],
            'Auth_seq_ID':               'deletedAtom.molResidue.serial',                                     # *** Mandatory
            'Auth_comp_ID':              'deletedAtom.molResidue.ccpCode.upper()',                            # *** Mandatory
            'Auth_atom_ID':             ('deletedAtom', self.getDeletedAtomName),                             # *** Mandatory
            #'Atom_type':                [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Assembly_ID':              [None,returnStarInt,'Assembly.ID',True],

            },
          },

        'Angle': {                                                                                            # *** Table

          'tags': {

            #'Angle_ID':                 [None,returnStarString,None,None],
            #'Angle_name':               [None,returnStarString,None,None],
            #'Entry_atom_ID_1':          [None,returnStarString,'Atom.Entry_atom_ID',None],
            #'Entity_assembly_ID_1':     [None,returnStarString,'Entity_assembly.ID',None],
            #'Entity_ID_1':              [None,returnStarString,'Entity.ID',None],
            #'Entity_label_1':           [None,returnStarString,None,None],
            #'Comp_ID_1':                [None,returnStarString,'Chem_comp.ID',None],
            #'Chem_comp_label_1':        [None,returnStarString,None,None],
            #'Comp_index_ID_1':          [None,returnStarString,'Entity_comp_index.ID',None],
            #'Seq_ID_1':                 [None,returnStarInt,'Entity_poly_seq.Num',None],
            #'Atom_ID_1':                [None,returnStarString,'Chem_comp_atom.Atom_ID',None],
            #'Atom_type_1':              [None,returnStarString,None,None],
            #'Entry_atom_ID_2':          [None,returnStarString,'Atom.Entry_atom_ID',None],
            #'Entity_assembly_ID_2':     [None,returnStarString,'Entity_assembly.ID',None],
            #'Entity_ID_2':              [None,returnStarString,'Entity.ID',None],
            #'Entity_label_2':           [None,returnStarString,None,None],
            #'Comp_ID_2':                [None,returnStarString,'Chem_comp.ID',None],
            #'Chem_comp_label_2':        [None,returnStarString,None,None],
            #'Comp_index_ID_2':          [None,returnStarString,'Entity_comp_index.ID',None],
            #'Seq_ID_2':                 [None,returnStarInt,'Entity_poly_seq.Num',None],
            #'Atom_ID_2':                [None,returnStarString,'Chem_comp_atom.Atom_ID',None],
            #'Atom_type_2':              [None,returnStarString,None,None],
            #'Entry_atom_ID_3':          [None,returnStarString,'Atom.Entry_atom_ID',None],
            #'Entity_assembly_ID_3':     [None,returnStarString,'Entity_assembly.ID',None],
            #'Entity_ID_3':              [None,returnStarString,'Entity.ID',None],
            #'Entity_label_3':           [None,returnStarString,None,None],
            #'Comp_ID_3':                [None,returnStarString,'Chem_comp.ID',None],
            #'Chem_comp_label_3':        [None,returnStarString,None,None],
            #'Comp_index_ID_3':          [None,returnStarString,'Entity_comp_index.ID',None],
            #'Seq_ID_3':                 [None,returnStarInt,'Entity_poly_seq.Num',None],
            #'Atom_ID_3':                [None,returnStarString,'Chem_comp_atom.Atom_ID',None],
            #'Atom_type_3':              [None,returnStarString,None,None],

            },
          },

        'Torsion_angle': {                                                                                    # *** Table

          'tags': {

            #'ID':                       [None,returnStarString,None,None],
            #'Torsion_angle_name':       [None,returnStarString,None,None],
            #'Entry_atom_ID_1':          [None,returnStarString,'Atom.Entry_atom_ID',None],
            #'Entity_assembly_ID_1':     [None,returnStarString,'Entity_assembly.ID',None],
            #'Entity_ID_1':              [None,returnStarString,'Entity.ID',None],
            #'Entity_label_1':           [None,returnStarString,None,None],
            #'Comp_ID_1':                [None,returnStarString,'Chem_comp.ID',None],
            #'Chem_comp_label_1':        [None,returnStarString,None,None],
            #'Comp_index_ID_1':          [None,returnStarString,'Entity_comp_index.ID',None],
            #'Seq_ID_1':                 [None,returnStarInt,'Entity_poly_seq.Num',None],
            #'Atom_ID_1':                [None,returnStarString,'Chem_comp_atom.Atom_ID',None],
            #'Atom_type_1':              [None,returnStarString,None,None],
            #'Entry_atom_ID_2':          [None,returnStarString,'Atom.Entry_atom_ID',None],
            #'Entity_assembly_ID_2':     [None,returnStarString,'Entity_assembly.ID',None],
            #'Entity_ID_2':              [None,returnStarString,'Entity.ID',None],
            #'Entity_label_2':           [None,returnStarString,None,None],
            #'Comp_ID_2':                [None,returnStarString,'Chem_comp.ID',None],
            #'Chem_comp_label_2':        [None,returnStarString,None,None],
            #'Comp_index_ID_2':          [None,returnStarString,'Entity_comp_index.ID',None],
            #'Seq_ID_2':                 [None,returnStarInt,'Entity_poly_seq.Num',None],
            #'Atom_ID_2':                [None,returnStarString,'Chem_comp_atom.Atom_ID',None],
            #'Atom_type_2':              [None,returnStarString,None,None],
            #'Entry_atom_ID_3':          [None,returnStarString,'Atom.Entry_atom_ID',None],
            #'Entity_assembly_ID_3':     [None,returnStarString,'Entity_assembly.ID',None],
            #'Entity_ID_3':              [None,returnStarString,'Entity.ID',None],
            #'Entity_label_3':           [None,returnStarString,None,None],
            #'Comp_ID_3':                [None,returnStarString,'Chem_comp.ID',None],
            #'Chem_comp_label_3':        [None,returnStarString,None,None],
            #'Comp_index_ID_3':          [None,returnStarString,'Entity_comp_index.ID',None],
            #'Seq_ID_3':                 [None,returnStarInt,'Entity_poly_seq.Num',None],
            #'Atom_ID_3':                [None,returnStarString,'Chem_comp_atom.Atom_ID',None],
            #'Atom_type_3':              [None,returnStarString,None,None],
            #'Entry_atom_ID_4':          [None,returnStarString,'Atom.Entry_atom_ID',None],
            #'Entity_assembly_ID_4':     [None,returnStarString,'Entity_assembly.ID',None],
            #'Entity_ID_4':              [None,returnStarString,'Entity.ID',None],
            #'Entity_label_4':           [None,returnStarString,None,None],
            #'Comp_ID_4':                [None,returnStarString,'Chem_comp.ID',None],
            #'Chem_comp_label_4':        [None,returnStarString,None,None],
            #'Comp_index_ID_4':          [None,returnStarString,'Entity_comp_index.ID',None],
            #'Seq_ID_4':                 [None,returnStarInt,'Entity_poly_seq.Num',None],
            #'Atom_ID_4':                [None,returnStarString,'Chem_comp_atom.Atom_ID',None],
            #'Atom_type_4':              [None,returnStarString,None,None],

            },
          },

        'Assembly_segment': {                                                                                 # *** Table

          'tags': {

            #'ID':                       [None,returnStarString,None,None],
            #'Entity_assembly_ID':       [None,returnStarString,'Entity_assembly.ID',None],
            #'Entity_ID':                [None,returnStarString,'Entity.ID',None],
            #'Entity_label':             [None,returnStarString,None,None],
            #'Comp_index_ID':            [None,returnStarString,'Entity_comp_index.ID',None],
            #'Comp_ID':                  [None,returnStarString,'Chem_comp.ID',None],
            #'Chem_comp_label':          [None,returnStarString,None,None],
            #'Seq_ID':                   [None,returnStarInt,'Entity_poly_seq.Num',None],
            #'Atom_ID':                  [None,returnStarString,'Chem_comp_atom.Atom_ID',None],
            #'Entry_atom_ID':            [None,returnStarString,'Atom.Entry_atom_ID',None],

            },
          },

        'Assembly_segment_description': {                                                                     # *** Table

          'tags': {

            #'Assembly_segment_ID':      [None,returnStarString,None,None],
            #'Code':                     [None,returnStarString,None,None],
            #'Details':                  [None,returnStarString,None,None],

            },
          },

        'Assembly_db_link': {                                                                                 # *** Includes optional fields

          'ccpnLoop':                    'molSystem.sortedDbReferences()',
          'ccpnMap':                     'dbReference',

          'tags': {

            'Author_supplied':           "dbReference.findFirstApplicationData(application='nmrStar',keyword='authDbAcc').value",
            'Database_code':             'dbReference.database.name',                                         # *** Optional - option - auto
            'Accession_code':            'dbReference.code',                                                  # *** Optional
            'Entry_mol_name':            'dbReference.name',
            #'Entry_mol_code':           [None,returnStarString,None,None],
            'Entry_experimental_method': ('dbReference',self.getDbRefExptMeth),                               # *** Optional - option
            'Entry_structure_resolution': ('dbReference',self.getDbRefStructRes),                             # *** Optional
            #'Entry_relation_type':      [None,returnStarString,None,None],                                   # *** Optional
            #'Entry_details':            [None,returnStarString,None,None],                                   # *** Optional

            },
          },

        'Assembly_keyword':              self.setKeywordLoop('molSystem.keywords','keyword'),                 # *** Table

        'Assembly_bio_function':         self.setFunctionLoop('molSystem.functions','function'),              # *** Table

        #'Assembly_bio_function': {                                                                           # *** Includes optional fields

        #  'tags': {

            #'Biological_function':      [None,returnStarString,None,None],                                   # *** Optional

        #    },
        #  },

        'Assembly_citation': {                                                                                # *** Table

          'tags': {

            #'Citation_ID':              [None,returnStarString,'Citation.ID',None],
            #'Citation_label':           [None,returnStarString,None,None],

            },
          },

        }
      }


    self.sfDict['entity'] = {                                                                                 # *** Priority 1

      #'category':                       'sequence',

      'ccpnLoop':                       ('molSystem',self.getMolecules),
      'ccpnMap':                         'molecule',  # Defines the main component for mapping...

      'title':                          ('molecule.name',self.getTitle),

      # TODO: Chem_comp.ID can apparently also be set at this stage: if single 'chemcomp'?

      'tags': {

        'Name':                         ('molecule',self.getMoleculeName),
        'Type':                         ('molecule',self.getMoleculeStarType),                                # *** Mandatory - pulldown - auto (Must have at least one polymer)
        'Number_of_monomers':            'molecule.seqLength',
        'Src_method':                   ('molecule',self.getSeqSrcMethod),
        #'Number_of_molecules':          [None,returnStarInt,None,None],
        'Paramagnetic':                  'molecule.isParamagnetic',                                           # *** Mandatory - bool - auto
        'Calc_isoelectric_point':        'molecule.calcIsoelectricPoint',
        #'Formula_weight_exptl':         [None,returnStarFloat,None,None],
        #'Formula_weight_exptl_meth':    [None,returnStarString,None,None],
        'Details':                       'molecule.details',                                                  # *** Optional
        'Thiol_state':                  ('molecule',self.getMoleculeThiolState),                              # *** Mandatory - option
        'Ambiguous_chem_comp_sites':    ('molecule',self.getAmbChemCompSites),                                # *** Mandatory - bool - auto
        'Ambiguous_conformational_states': ('molecule',self.getAmbConfStates),                                # *** Mandatory - bool - auto
        'Nstd_chirality':                'molecule.hasNonStdChirality',                                       # *** Mandatory - bool - auto
        'Nstd_linkage':                 ('molecule.isStdLinear',self.negateBoolean),                          # *** Mandatory - bool - auto

        #'Sf_framecode':                 [None,returnStarString,None,False],                                  # *** Mandatory - set in NmrStarExport - option

        'CONDITIONAL': { # NOTE HERE: can also use 'DEFAULT' value, which is the 'standard' option if all else fails

          ('molecule',self.getMoleculeStarType): {

            'polymer': { # water, aggregate and solvent all have similar extra polymer options for some reason - but not Polymer_type

              'Parent_entity_ID=LOCAL':  None, # Otherwise sets this to current entity_ID. Only seems to be set for non-polymers.
              'Polymer_type':           ('molecule',self.getMoleculeStarPolymerType),                         # *** Mandatory - pulldown
              #'Polymer_type_details':   [None,returnStarString,None,None],
              'Polymer_seq_one_letter_code': ('molecule', self.getSeqString),                       # *** Mandatory
              'Polymer_strand_ID':      ('molecule',self.getChainCodes),                                      # *** Mandatory (PDB)
              'Nstd_monomer':            'molecule.hasNonStdChemComp',                                        # *** Mandatory - bool - auto
              'Polymer_author_defined_seq': "molecule.findFirstApplicationData(application='nmrStar',keyword='authDefSeq').value", # *** Optional
              'Polymer_author_seq_details': 'molecule.seqDetails',                                            # *** Optional

              'Mutation':               ('molecule',self.getSeqMutation),                                     # *** Optional
              'Fragment':               ('molecule',self.getSeqFragment),                                     # *** Optional
              'Formula_weight':          'molecule.molecularMass',                                            # *** Optional
              'EC_number':              ('molecule',self.getMoleculeECNumber),                                # *** Optional

              },

            'non-polymer': {

              'Nonpolymer_comp_ID':    ('molecule.findFirstMolResidue().chemCompVar.chemComp',self.getChemCompCif),
              'Nonpolymer_comp_label': ('molecule.findFirstMolResidue().chemCompVar',self.getCCVTitle),
              'Parent_entity_ID':       'molecule',

              },
            },

          },
        },

      'tables': {

        'Entity_common_name': {                                                                               # *** Table

          'tags': {

            #'Name':                     [None,returnStarString,None,None],
            #'Type':                     [None,returnStarString,None,None],

            },
          },

        #'Entry_ID':                     [None,returnStarString,'Entry.ID',True],
        
        'Entity_systematic_name':        self.setSysNamesLoop('molecule.moleculeSysNames','moleculeSysName'), # *** Table

        #'Entry_ID':                     [None,returnStarString,'Entry.ID',True],

        'Entity_keyword':                self.setKeywordLoop('molecule.keywords','keyword'),                  # *** Table

        'Entity_biological_function':    self.setFunctionLoop('molecule.functions','function'),               # *** Table

        #'Entity_biological_function': {                                                                      # *** Includes optional fields

        #  'tags': {

            #'Biological_function':      [None,returnStarString,None,None],                                   # *** Optional

        #    },
        #  },

        'Entity_db_link': {                                                                                   # *** Includes optional fields - not for non-polymer

          'CONDITIONAL':                ( ('molecule',self.getMoleculeStarType),('polymer',) ),

          'ccpnLoop':                   ('molecule',self.getDbAlignments),
          'ccpnMap':                     'alignment',

          'addCcpnMap':                  {'alignment': [('dbRef','alignment.dbRef')]},

          'tags': {

            'Database_code':             'dbRef.database.name',                                               # *** Optional - option
            'Accession_code':            'dbRef.code',                                                        # *** Optional
            #'Chimera_segment':           'alignment.serial',  # TODO: Is this right???
            'Entry_mol_name':            'dbRef.name',                                                        # *** Optional
            #'Seq_align_begin':           'alignment.dbRefAlignBegin',
            #'Seq_align_end':             'alignment.dbRefAlignEnd',
            'Entry_details':             'dbRef.details',

            'Author_supplied':           "dbRef.findFirstApplicationData(application='nmrStar',keyword='authDbAcc').value",
            #'Entry_mol_code':           [None,returnStarString,None,False],
            #'Entry_experimental_method': [None,returnStarString,None,False],
            #'Entry_structure_resolution': [None,returnStarFloat,None,False],
            #'Entry_relation_type':      [None,returnStarString,None,False],
            #'Chimera_segment':          [None,returnStarInt,None,False],
            #'Seq_query_to_submitted_percent': [None,returnStarFloat,None,False],
            #'Seq_subject_length':       [None,returnStarInt,None,False],
            #'Seq_identity':             [None,returnStarFloat,None,False],
            #'Seq_positive':             [None,returnStarFloat,None,False],
            #'Seq_homology_expectation_val': [None,returnStarFloat,None,False],
            #'Seq_difference_details':   [None,returnStarString,None,False],
            #'Seq_alignment_details':    [None,returnStarString,None,False],
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Entity_ID':                [None,returnStarInt,'Entity.ID',True],

            },
          },

        'Entity_comp_index': {                                                                                # *** Table

          #'category':                   'sequence',

          'ccpnLoop':                   ('molecule',self.getMolResidues),
#          'ccpnLoop':                   'molecule.sortedMolResidues()',
          'ccpnMap':                     'molResidue',  # Defines the main component for mapping...

          'tags': {

            'ID':                        'molResidue.serial',
#            'Auth_seq_ID':               'molResidue.seqCode',
#            'Auth_seq_ID':               "molResidue.findFirstApplicationData(application='nmrStar',keyword='authorSeqCode').value",
            'Auth_seq_ID':              ('molResidue',self.getAuthSeqId),
            'Comp_ID':                  ('molResidue.chemCompVar.chemComp',self.getChemCompCif),

            },
          },

        'Entity_poly_seq': {                                                                                  # *** Table
    
          # TODO Does this only have to be set for hetero groups?
          #'category':                   'sequence',

          'CONDITIONAL':                ( ('molecule',self.getMoleculeStarType),('polymer',) ),

          'ccpnLoop':                    'molecule.sortedMolResidues()',
          'ccpnMap':                     'molResidue',  # Defines the main component for mapping...

          'tags': {

            #'Hetero':                   'molResidue.molType',  # How essential is this? Is this for 'other'?
            'Mon_ID':                   ('molResidue.chemCompVar.chemComp',self.getChemCompCif), # WARNING: THIS an ID but not an integer...
            'Num':                       'molResidue.serial',
            'Comp_index_ID':             'molResidue', # This is link to Entity_comp_index.ID
    
            },
          },

        'Entity_comp_index_alt': {                                                                            # *** Includes optional fields - not for non-polymer
    
          #'CONDITIONAL':               ( ('molecule',self.getMoleculeStarType),('polymer',) ),

          'tags': {

            #'Auth_seq_ID':              [None,lambda x = value: returnStarCode(x,length = 12),None,False],   # *** Optional
            #'Comp_label':               [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
    
            },
          },


        'Entity_fragment': {                                                                                  # *** Table

          'tags': {

            #'ID':                       [None,returnStarString,None,None],
            #'Comp_index_ID_begin':      [None,returnStarString,'Entity_comp_index.ID',None],
            #'Comp_index_ID_end':        [None,returnStarString,'Entity_comp_index.ID',None],
            #'Seq_ID_begin':             [None,returnStarInt,'Entity_poly_seq.Num',None],
            #'Seq_ID_end':               [None,returnStarInt,'Entity_poly_seq.Num',None],
            #'Details':                  [None,returnStarInt,None,None],

            },
          },

        'Entity_citation': {                                                                                  # *** Table

          'tags': {

            #'Citation_ID':              [None,returnStarString,'Citation.ID',None],
            #'Citation_label':           [None,returnStarString,None,None],

            },
          },

        'Entity_chimera_segment': {                                                                           # *** Table

          #'CONDITIONAL':                ( ('molecule',self.getMoleculeStarType),('polymer',) ),

          #'ccpnLoop':                   ('molecule',self.getDbAlignments),
          #'ccpnMap':                     'alignment',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Details':                  [None,returnStarString,None,False],
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Entity_ID':                [None,returnStarInt,'Entity.ID',True],

            'Comp_index_ID_begin':      ('alignment.molSeqFragment',self.getStartLimitResidue),
            'Comp_index_ID_end':        ('alignment.molSeqFragment',self.getEndLimitResidue),
            'Seq_ID_begin':              'alignment.dbRefAlignBegin',
            'Seq_ID_end':                'alignment.dbRefAlignEnd',

            },
          },
        
        }
      }


    #######################
    # Species information #
    #######################

    self.sfDict['natural_source'] = {                                                                         # *** Priority 1

      'ccpnMap':                         'nmrEntry',

      #'title':                          ('naturalSource.organismName',self.getTitle),

      'tags': {

        #'Sf_category':                  ['natural_source',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],

        },

      'tables': {

        'Entity_natural_src': {                                                                               # *** Includes mandatory fields

          'ccpnLoop':                   ('nmrEntry',self.getEntryMolecules),
          'ccpnMap':                     'entryMolecule',

          'addCcpnMap':                  {'entryMolecule': [('naturalSource','entryMolecule.molecule.naturalSource')]},

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            'Entity_ID':                 'entryMolecule.molecule',
            'Entity_label':              'entryMolecule.molecule.name',                                       # *** Mandatory - option
            #'Entity_chimera_segment_ID': [None,returnStarInt,'Entity_chimera_segment.ID',False],
            'NCBI_taxonomy_ID':          'naturalSource.ncbiTaxonomyId',
            'Type':                      'entryMolecule.sourceType',                                          # *** Mandatory - option - auto
            #'Common':                   [None,returnStarString,None,False],
            'Organism_name_scientific':  'naturalSource.scientificName',                                      # *** Mandatory - option
            'Organism_name_common':     ('naturalSource.organismName',self.convertUnknown),
#            'Organism_name_common':      'naturalSource.organismName',
            'Organism_acronym':          'naturalSource.organismAcronym',
            'ICTVdb_decimal_code':       'naturalSource.ictvCode',
            'Superkingdom':              "naturalSource.findFirstApplicationData(application='nmrStar',keyword='superKingdom').value",
            'Kingdom':                   "naturalSource.findFirstApplicationData(application='nmrStar',keyword='kingdom').value",
            'Genus':                     'naturalSource.genus',
            'Species':                   'naturalSource.species',
            'Strain':                    'naturalSource.strain',                                              # *** Optional
            'Variant':                   'naturalSource.variant',                                             # *** Optional
            'Subvariant':                'naturalSource.subVariant',
            'Organ':                     'naturalSource.organ',
            'Tissue':                    'naturalSource.tissue',
            #'Tissue_fraction':          [None,returnStarString,None,False],
            'Cell_line':                 'naturalSource.cellLine',
            'Cell_type':                 'naturalSource.cellType',
            'ATCC_number':               'naturalSource.atccNumber',
            'Organelle':                 'naturalSource.organelle',
            'Cellular_location':         'naturalSource.cellLocation',
            #'Fragment':                 [None,returnStarString,None,False],
            'Fraction':                  'naturalSource.fraction',
            'Secretion':                 'naturalSource.secretion',
            'Plasmid':                   'naturalSource.plasmid',
            'Plasmid_details':           'naturalSource.plasmidDetails',
            'Gene_mnemonic':             'naturalSource.geneMnemonic',                                        # *** Optional
            #'Dev_stage':                [None,returnStarString,None,False],
            'Details':                   'naturalSource.details',                                             # *** Optional
#            'Citation_ID=LOCAL':         None, #[None,returnStarInt,'Citation.ID',False], # Set this to None, otherwise it sets this to the wrong citation
            #'Citation_label':           [None,returnStarString,None,False],
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Entity_natural_src_list_ID': [None,returnStarInt,'Entity_natural_src_list.ID',True],

            },
          },

        'Natural_source_db': {                                                                                # *** Table

          'tags': {

            #'Entity_natural_src_ID':    [None,returnStarInt,'Entity_natural_src.ID',True],
            #'Entity_ID':                [None,returnStarInt,'Entity.ID',True],
            #'Entity_label':             [None,returnStarString,'Entity.Sf_framecode',True],
            #'Entity_chimera_segment_ID': [None,returnStarInt,'Entity_chimera_segment.ID',False],
            #'Database_code':            [None,returnStarString,None,True],
            #'Database_type':            [None,returnStarString,None,True],
            #'Entry_code':               [None,returnStarString,None,True],
            #'Entry_type':               [None,returnStarString,None,True],
            #'ORF_code':                 [None,returnStarString,None,False],
            #'Gene_locus_code':          [None,returnStarString,None,False],
            #'Gene_cDNA_code':           [None,returnStarString,None,False],
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Entity_natural_src_list_ID': [None,returnStarInt,'Entity_natural_src_list.ID',True],

            },
          },

        }
      }


    self.sfDict['experimental_source'] = {                                                                    # *** Priority 1

      'ccpnMap':                         'nmrEntry',

      #'title':                          ('experimentalSource.organismName',self.getTitle),

      'tags': {

        #'Sf_category':                  ['experimental_source',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],

        },

      'tables': {

        'Entity_experimental_src': {                                                                          # *** Includes mandatory fields

          'ccpnLoop':                    'nmrEntry.sortedEntryMolecules()',
          'ccpnMap':                     'entryMolecule',

          'addCcpnMap':                  {'entryMolecule': [('experimentalSource','entryMolecule.experimentalSource')]},

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            'Entity_ID':                 'entryMolecule.molecule',
            'Entity_label':              'entryMolecule.molecule.name',                                       # *** Mandatory - option
            #'Entity_chimera_segment_ID': [None,returnStarInt,'Entity_chimera_segment.ID',False],
            'Production_method':         'entryMolecule.productionMethod',                                    # *** Mandatory - option - auto
            'Host_org_scientific_name':  'experimentalSource.scientificName',                                 # *** Mandatory - option
            'Host_org_name_common':     ('experimentalSource.organismName',self.convertUnknown),
#            'Host_org_name_common':      'experimentalSource.organismName',
            #'Host_org_details':         [None,returnStarString,None,False],
            'Host_org_NCBI_taxonomy_ID': 'experimentalSource.ncbiTaxonomyId',
            'Host_org_genus':            'experimentalSource.genus',
            'Host_org_species':          'experimentalSource.species',
            'Host_org_strain':           'experimentalSource.strain',                                         # *** Optional
            'Host_org_variant':          'experimentalSource.variant',                                        # *** Optional
            'Host_org_subvariant':       'experimentalSource.subVariant',
            'Host_org_organ':            'experimentalSource.organ',
            'Host_org_tissue':           'experimentalSource.tissue',
            #'Host_org_tissue_fraction': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            'Host_org_cell_line':        'experimentalSource.cellLine',
            'Host_org_cell_type':        'experimentalSource.cellType',
            'Host_org_cellular_location': 'experimentalSource.cellLocation',
            'Host_org_organelle':        'experimentalSource.organelle',
            #'Host_org_gene':            [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Host_org_culture_collection': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            'Host_org_ATCC_number':      'experimentalSource.atccNumber',
            'PDBview_host_org_vector_name': "experimentalSource.findFirstApplicationData(application='nmrStar',keyword='pdbVectorName').value",
            #'PDBview_host_org_vector_name': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'PDBview_plasmid_name':     [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            'Vector_type':               'entryMolecule.vectorType',                                          # *** Mandatory - option - auto
            'Vector_name':               'experimentalSource.plasmid',                                        # *** Mandatory
            'Vector_details':            'experimentalSource.plasmidDetails',
            #'Vendor_name':              [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Host_org_dev_stage':       [None,returnStarString,None,False],
            'Details':                   'experimentalSource.details',                                        # *** Optional
#            'Citation_ID=LOCAL':         None, #[None,returnStarInt,'Citation.ID',False], # Set this to None, otherwise it sets this to the wrong citation
            #'Citation_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Entity_experimental_src_list_ID': [None,returnStarInt,'Entity_experimental_src_list.ID',True],

            },
          },

        }
      }


    self.sfDict['chem_comp'] = {                                                                              # *** Priority 1
    
      # TODO check this: shouldn't the entry or molsystem be passed in here??!
      'ccpnLoop':                        self.getNonStdChemCompVarList,
      'ccpnMap':                         'nonStdChemCompVar',

      'addCcpnMap':                     {'nonStdChemCompVar': [('chemComp','nonStdChemCompVar.chemComp')]},
      
      'title':                          ('nonStdChemCompVar',self.getCCVTitle),

      'tags': {

        'ID':                           ('nonStdChemCompVar',self.getChemCompVarCif),
        'Name':                          'nonStdChemCompVar.name',                                            # *** Mandatory - if there are ligands/non-standard residues
        #'BMRB_code':                    '', # TODO could make this part of ref data
        'PDB_code':                      'chemComp.code3Letter',                                              # *** Optional
        #'InCHi_code':                   '',
        'Type':                         ('nonStdChemCompVar',self.getChemCompVarStarType),
        #'PDB_NSTD_flag':                '',
        'Std_deriv_one_letter_code':     'chemComp.code1Letter',
        'Std_deriv_three_letter_code':   'chemComp.stdChemCompCode',
        #'Std_deriv_BMRB_code':          '',
        'Std_deriv_PDB_code':            'chemComp.stdChemComp.code3Letter',
        'Std_deriv_chem_comp_name':      'chemComp.stdChemComp.name',
        'Formal_charge':                 'nonStdChemCompVar.formalCharge',                                    # *** Optional - option
        'Paramagnetic':                  'nonStdChemCompVar.isParamagnetic',                                  # *** Mandatory - if there are ligands/non-standard residues - bool
        'Aromatic':                      'nonStdChemCompVar.isAromatic',                                      # *** Mandatory - if there are ligands/non-standard residues - bool
        'Formula':                      ('nonStdChemCompVar.formula',self.getChemCompFormula),                # *** Optional
        'Formula_weight':                'nonStdChemCompVar.molecularMass',                                   # *** Optional
        #'Image_file_name':              '',                                                                  # *** Optional - option
        #'Image_file_format':            '',                                                                  # *** Optional - option
        #'Topo_file_name':               '',
        #'Topo_file_format':             '',
        #'Struct_file_name':             '',
        #'Struct_file_format':           '',
        #'Stereochem_param_file_name':   '',
        #'Stereochem_param_file_format': '',
        #'Vendor':                       '',
        #'Vendor_product_code':          '',
        'Details':                       'chemComp.details',                                                  # *** Optional
        #'DB_query_date':                '',
        #'DB_last_query_revised_last_date': '',

        },

      'tables': {

        'Chem_comp_common_name':         self.setCommonNameLoop('nonStdChemCompVar.chemComp.commonNames','commonName'), # *** Table


        'Chem_comp_descriptor': {

          'tags': {

            #'Ordinal':                  [None,returnStarInt,None,True],
            #'Descriptor':               [None,lambda x = value: returnStarString(x,length = 1024),None,True],
            #'Type':                     [None,lambda x = value: returnStarLine(x,length = 127),None,True],
            #'Program':                  [None,lambda x = value: returnStarLine(x,length = 127),None,True],
            #'Program_version':          [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

            },
          },

        'Chem_comp_identifier': {

          'tags': {

            #'Ordinal':                  [None,returnStarInt,None,True],
            #'Identifier':               [None,lambda x = value: returnStarString(x,length = 1024),None,True],
            #'Type':                     [None,lambda x = value: returnStarLine(x,length = 127),None,True],
            #'Program':                  [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Program_version':          [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Comp_id':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

            },
          },

        # TODO ChemAtomSysNames!

        'Chem_comp_systematic_name': {                                                                        # *** Table

          'tags': {

            #'Name':                     [None,lambda x = value: returnStarLine(x,length = 1024),None,True],
            #'Naming_system':            [None,lambda x = value: returnStarLine(x,length = 127),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

            },
          },

        'Chem_comp_SMILES': {                                                                                 # *** Table

          'tags': {

            #'Type':                     [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'String':                   [None,lambda x = value: returnStarString(x,length = 1024),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

            },
          },

        'Chem_comp_keyword':             self.setKeywordLoop('nonStdChemCompVar.chemComp.keywords','keyword'), # *** Table

        'Characteristic': {                                                                                   # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Name':                     [None,lambda x = value: returnStarLine(x,length = 127),None,True],
            #'Atom_ID':                  [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Chemical_group':           [None,lambda x = value: returnStarLine(x,length = 127),None,False],

            #'Val':                      [None,returnStarFloat,None,True],
            #'Val_err':                  [None,returnStarFloat,None,True],

            #'Source':                   [None,lambda x = value: returnStarLine(x,length = 127),None,True],
            #'Citation_ID':              [None,returnStarInt,'Citation.ID',True],
            #'Citation_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

            },
          },

        # CJP - missing tags marked with + here:

        'Chem_comp_atom': {                                                                                   # *** Table

          'tags': {

        #  + 'Atom_ID':                  [None,lambda x = value: returnStarAtCode(x,length = 12),None,True],
        #    'PDB_atom_ID':              [None,lambda x = value: returnStarAtCode(x,length = 12),None,False],
        #    'Alt_atom_ID':              [None,lambda x = value: returnStarAtCode(x,length = 12),None,False],
        #    'Auth_atom_ID':             [None,lambda x = value: returnStarCode(x,length = 15),None,False],
        #  + 'Type_symbol':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
        #    'Isotope_number':           [None,returnStarInt,None,False],
        #  + 'Chirality':                [None,lambda x = value: returnStarCode(x,length = 15),None,False],
        #  + 'Charge':                   [None,lambda x = value: returnStarCode(x,length = 15),None,False],
        #    'Partial_charge':           [None,returnStarFloat,None,False],
        #  + 'Oxidation_number':         [None,lambda x = value: returnStarCode(x,length = 15),None,False],
        #  + 'Unpaired_electron_number': [None,lambda x = value: returnStarString(x,length = 3),None,False],
        #    'PDBx_aromatic_flag':       [None,lambda x = value: returnStarCode(x,length = 3),None,False],
        #    'PDBx_leaving_atom_flag':   [None,lambda x = value: returnStarCode(x,length = 3),None,False],
        #    'Substruct_code':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
        #    'Ionizable':                [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
        #    '2D_drawing_coord_x':       [None,returnStarFloat,None,False],
        #    '2D_drawing_coord_y':       [None,returnStarFloat,None,False],
        #    'Model_Cartn_x':            [None,returnStarFloat,None,False],
        #    'Model_Cartn_y':            [None,returnStarFloat,None,False],
        #    'Model_Cartn_z':            [None,returnStarFloat,None,False],
        #    'Details':                  [None,returnStarString,None,False],
        #  + 'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #  + 'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

            },
          },

        'Atom_nomenclature': {                                                                                # *** Table

          'tags': {

            #'Atom_ID':                  [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_name':                [None,lambda x = value: returnStarLine(x,length = 15),None,True],
            #'Naming_system':            [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

            },
          },

        # CJP - missing tags marked with + here:

        'Chem_comp_bond': {                                                                                   # *** Table

          'tags': {

        #    'ID':                       [None,returnStarInt,None,True],
        #    'Type':                     [None,lambda x = value: returnStarLine(x,length = 31),None,True],
        #    'Value_order':              [None,lambda x = value: returnStarLine(x,length = 31),None,True],
        #    'Atom_ID_1':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
        #    'Atom_ID_2':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
        #    'PDB_atom_ID_1':            [None,lambda x = value: returnStarAtCode(x,length = 12),None,False],
        #    'PDB_atom_ID_2':            [None,lambda x = value: returnStarAtCode(x,length = 12),None,False],
        #    'Details':                  [None,returnStarString,None,False],
        #  + 'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #  + 'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

            },
          },

        'Chem_comp_tor': {                                                                                    # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Atom_ID_1':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_ID_2':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_ID_3':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_ID_4':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Details':                  [None,returnStarString,None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

            },
          },

        'Chem_comp_angle': {                                                                                  # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Atom_ID_1':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_ID_2':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_ID_3':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Details':                  [None,returnStarString,None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

            },
          },

        'Chem_comp_db_link': {                                                                                # *** Table

          'tags': {

            #'Author_supplied':          [None,lambda x = value: returnStarYesNo(x,length = 12),None,False],
            #'Database_code':            [None,lambda x = value: returnStarLine(x,length = 127),None,True],
            #'Accession_code':           [None,lambda x = value: returnStarLine(x,length = 15),None,True],
            #'Accession_code_type':      [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Entry_mol_code':           [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Entry_mol_name':           [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Entry_experimental_method': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Entry_relation_type':      [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Entry_details':            [None,returnStarString,None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

            },
          },

        'Chem_comp_citation': {                                                                               # *** Table

          'tags': {

            #'Citation_ID':              [None,returnStarInt,'Citation.ID',True],
            #'Citation_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],

            },
          },

        }
      }   


    self.sfDict['sample'] = {                                                                                 # *** Priority 1

      'ccpnLoop':                       ('nmrEntry.sortedExperiments()', self.getSamples),
      'ccpnMap':                         'sample',

      'title':                          ('sample.name',self.getTitle),

      'tags': {

        #'Sf_category':                  ['sample',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - option - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Type':                         ('sample',self.getStarSampleState),                                   # *** Mandatory - option - auto
        #'Sub_type':                     [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        'Details':                       'sample.details',                                                    # *** Optional
        'Aggregate_sample_number':       "sample.findFirstApplicationData(application='nmrStar',keyword='aggrSamNo').value", # *** Optional - auto
        'Solvent_system':                "sample.findFirstApplicationData(application='nmrStar',keyword='solventSys').value", # *** Mandatory - option
        #'Preparation_date':             [None,returnStarDateTime,None,False],
        #'Preparation_expiration_date':  [None,returnStarDateTime,None,False],
        #'Oriented_sample_prep_protocol': [None,returnStarString,None,False],
        #'Lyophilization_cryo_protectant': [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Storage_protocol':             [None,returnStarString,None,False],
#        'Crystal_grow_method_cit_ID=LOCAL': None, #[None,returnStarInt,'Citation.ID',False], # Set this to None, otherwise it sets this to the wrong citation
#        'Crystal_grow_seeding_cit_ID=LOCAL': None, #[None,returnStarInt,'Citation.ID',False], # Set this to None, otherwise it sets this to the wrong citation
        
        },

      'tables': {

        # Best way to handle SampleComponents...

        'Sample_component': {                                                                                 # *** Includes mandatroy fields

          'ccpnLoop':                    'sample.sortedSampleComponents()',
          'ccpnMap':                     'sampleComponent',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],

            'Mol_common_name':          ('sampleComponent.refComponent',self.getCompMolName),                 # *** Mandatory - option
            'Isotopic_labeling':        ('sampleComponent.refComponent',self.getStarIsotopeLabelling),
            #'Isotopic_labeling':        [None,lambda x = value: returnStarLine(x,length = 127),None,False],
#            'Assembly_ID=LOCAL':         None, #[None,returnStarInt,'Assembly.ID',False], # Set this to None, otherwise it sets this to the wrong assembly
            #'Assembly_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            'Entity_ID':                ('sampleComponent.refComponent',self.getComponentMolecule),
            'Entity_label':             ('sampleComponent.refComponent',self.getComponentMoleculeName), 
            #'Entity_label':             [None,lambda x = value: returnStarLabel(x,length = 127),'Entity.Sf_framecode',False],
            #'Product_ID':               [None,returnStarInt,None,False],
            #'Type':                     [None,lambda x = value: returnStarLine(x,length = 127),None,False],

            'Concentration_val':        ('sampleComponent',self.getComponentConc),                            # *** Mandatory
            'Concentration_val_units':  ('sampleComponent',self.getComponentUnit),                            # *** Mandatory - option - auto
            'Concentration_val_err':     'sampleComponent.concentrationError',                                # *** Optional

            #'Concentration_val_min':    [None,returnStarFloat,None,False],                                   # *** Can be used as a range instead of Concentration_val
            #'Concentration_val_max':    [None,returnStarFloat,None,False],                                   # *** Can be used as a range instead of Concentration_val

            #'Vendor':                   [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Vendor_product_name':      [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Vendor_product_code':      [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],

            },
          },

        # Hacky way to deal with SampleComponents...

        #'Sample_component': {

        #  'ccpnMap':                    'sample',

        #  'tags': {

        #    'Mol_common_name':          None,
        #    'Entity_ID':                None,
        #    'Concentration_val':       ('sample.details',self.getSampleConcValue),
        #    'Concentration_val_units': ('sample.details',self.getSampleConcUnit),

        #    },
        #  },

        }
      }


    self.sfDict['sample_conditions'] = {                                                                      # *** Priority 1

      'ccpnLoop':                       ('nmrEntry.sortedExperiments()', self.getSampleConditionSets),
      'ccpnMap':                         'sampleConditionSet',

      'title':                          ('sampleConditionSet.name',self.getTitle),

      'tags': {

        #'Sf_category':                  ['sample_conditions',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - option - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Details':                       'sampleConditionSet.details',                                        # *** Optional

        },

      'tables': {

        'Sample_condition_variable': {                                                                        # *** Includes mandatory fields

          'ccpnLoop':                   ('sampleConditionSet',self.getSampleConditions),  # Could be sampleConditionSet.getSortedSampleConditions()
          #'ccpnLoop':                  ('nmrEntry.sortedExperiments()', self.getSampleConditions),
          'ccpnMap':                     'sampleCondition',

          'tags': {

            'Type':                     ('sampleCondition.condition',self.getStarConditionType),              # *** Mandatory - pulldown - auto (x4 - temperature, pH, pressure, ionic strength)
            'Val':                      ('sampleCondition.value',self.getNonZeroValue),                       # *** Mandatory - pressure is auto set
            'Val_units':                 'sampleCondition.unit',                                              # *** Mandatory - option - auto (x4 - K, pH, atm, M)
            'Val_err':                   'sampleCondition.error',                                             # *** Optional
            #'Sample_condition_list_ID': 'sampleConditionSet',

            },
          },

        }
      }


    self.sfDict['software'] = {                                                                               # *** Priority 1

      'ccpnLoop':                       ('nmrEntry', self.getSoftware),
      #'ccpnLoop':                        'nmrEntry.sortedSoftware()',
      'ccpnMap':                         'software',

      'title':                          ('software', self.getSoftwareTitle),

      'tags': {

        #'Sf_category':                  ['software',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - option
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],

        'Name':                          'software.name',
        'Version':                      ('software.version', self.convertUnknown),                            # *** Optional
        'Details':                       'software.details',                                                  # *** Optional

        },

      'tables': {

        'Vendor': {                                                                                           # *** Includes mandatory fields

          # TODO: add loop here
          'ccpnLoop':                   ('software', self.getVendorNames),
          'ccpnMap':                     'vendorName',

          'tags': {
            'Name':                      'vendorName',                                                        # *** Mandatory - option
            'Address':                   'software.vendorAddress',                                            # *** Optional
            'Electronic_address':        'software.vendorWebAddress',                                         # *** Optional

            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Software_ID':              [None,returnStarInt,'Software.ID',True],

            },
          },

        'Task': {                                                                                             # *** Includes mandatory fields

          'ccpnLoop':                    'software.tasks',
          'ccpnMap':                     'task',

          'tags': {

            'Task':                      'task',                                                              # *** Mandatory - option

            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Software_ID':              [None,returnStarInt,'Software.ID',True],

            },
          },

        }
      }


    self.sfDict['method'] = {                                                                                 # *** Priority 4 - Not in ADIT-NMR?

      'ccpnLoop':                       ('nmrEntry', self.getMethods),
      'ccpnMap':                         'method',

      'title':                          ('method.name',self.getTitle),

      'tags': {

        #'Sf_category':                  ['method',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Derivation_type':              'method.name',
        'Details':                       'method.name',
        #'Computer_ID':                  [None,returnStarInt,'Computer.ID',False],
        #'Computer_label':               [None,lambda x = value: returnStarLabel(x,length = 127),None,False],

        },
      }


    self.sfDict['NMR_spectrometer'] = {                                                                       # *** Priority 1 - Not in ADIT-NMR?

      'ccpnLoop':                       ('nmrEntry', self.getSpectrometers),
      'ccpnMap':                         'spectrometer',

      'title':                          ('spectrometer.name',self.getTitle),

      'tags': {

        'Name':                          'spectrometer.name',
        'Model':                         'spectrometer.model',
        'Manufacturer':                  'spectrometer.manufacturer.name',
        'Serial_number':                 'spectrometer.serialNumber',
        'Field_strength':                'spectrometer.nominalFreq',
        'Details':                       'spectrometer.details',

        },
      }


    self.sfDict['NMR_spectrometer_list'] = {                                                                  # *** Priority 1

      'ccpnMap':                         'nmrEntry',

      'title':                           None, # TODO - need to do something here and for Sf_framecode

      'tags': {
        #'Sf_category':                  ['NMR_spectrometer_list',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        'Details':                       'nmrEntry.spectrometerListDetails',
        #'ID':                           [None,returnStarInt,None,True],

        },

      'tables': {

        'NMR_spectrometer_view': {                                                                            # *** Includes mandatory fields

          'ccpnLoop':                   ('nmrEntry', self.getSpectrometers),
          'ccpnMap':                     'spectrometer',

          'tags': {

            # Gives four by default all with the same name!!!

            'Name':                      'spectrometer.name',                                                 # *** Mandatory - option - auto
            'Manufacturer':              'spectrometer.manufacturer.name',                                    # *** Mandatory - option
            'Model':                     'spectrometer.model',                                                # *** Mandatory - option
            'Details':                   'spectrometer.details',                                              # *** Optional
            'Field_strength':            'spectrometer.nominalFreq',                                          # *** Mandatory - option
            'Serial_number':             'spectrometer.serialNumber',
            
#            'Citation_ID=LOCAL':         None, #[None,returnStarInt,'Citation.ID',False], # Set this to None, otherwise it sets this to the wrong citation

            },
          },

        }
      }


    self.sfDict['NMR_spectrometer_probe'] = {                                                                 # *** Priority 3 - Not in ADIT-NMR?

      'ccpnLoop':                       ('nmrEntry.sortedExperiments()', self.getNmrProbes),
      'ccpnMap':                         'nmrProbe',

      'title':                          ('nmrProbe.name',self.getTitle),

      'tags': {

        #'Sf_category':                  ['NMR_spectrometer_probe',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],

        'Details':                       'nmrProbe.details',
        'Manufacturer':                  'nmrProbe.manufacturer.name',
        'Model':                         'nmrProbe.model',
        'Serial_number':                 'nmrProbe.serialNumber',
        'Diameter':                      'nmrProbe.diameter',

        #'Rotor_length':                 [None,lambda x = value: returnStarString(x,length = 127),None,False],
        #'Rotor_composition':            [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Internal_volume':              [None,returnStarFloat,None,False],
        #'Spacer_present':               [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],

        },

      'tables': {

        'NMR_probe': {                                                                                        # *** Table

          'ccpnMap':                     'nmrProbe',

          'tags': {
            'Type':                      'nmrProbe.probeType',
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'NMR_spectrometer_probe_ID': [None,returnStarInt,'NMR_spectrometer_probe.ID',True],

            },
          },

        'NMR_spectrometer_probe_citation': {                                                                  # *** Table

          'tags': {

            #'Citation_ID':              [None,returnStarInt,'Citation.ID',True],
            #'Citation_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'NMR_spectrometer_probe_ID': [None,returnStarInt,'NMR_spectrometer_probe.ID',True],

            },
          },

        }
      }


    self.sfDict['experiment_list'] = {                                                                        # *** Priority 1

      'ccpnMap':                         'nmrEntry',

      'title':                           None, # TODO - do something here and for Sf_framecode

      'tags': {

        #'Sf_category':                  ['experiment_list',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Details':                       "nmrEntry.findFirstApplicationData(application='nmrStar',keyword='experimentListDetails').value",

        },

      'tables': {

        'Experiment': {                                                                                       # *** Includes mandatory fields

          'ccpnLoop':                    'nmrEntry.sortedExperiments()',
          'ccpnMap':                     'experiment',

          'tags': {

            'Name':                     ('experiment',self.getStarExperimentName),                            # *** Mandatory - option
            'Raw_data_flag':            ('experiment',self.getRawDataFlag),                                   # *** Mandatory - bool - auto
            'Sample_ID':                 'experiment.sample',
            'Sample_label':             ('experiment.sample.name',self.getLabel),                             # *** Mandatory - option
            'Sample_state':             ('experiment',self.getStarExpState),                                  # *** Mandatory - option - auto
            'Sample_condition_list_ID':  'experiment.sampleConditionSet',
            'Sample_condition_list_label': ('experiment.sampleConditionSet.name',self.getLabel),              # *** Mandatory - option
            'NMR_spectrometer_ID':       'experiment.spectrometer',
            'NMR_spectrometer_label':   ('experiment.spectrometer.name',self.getLabel),                       # *** Optional
            'NMR_spectrometer_probe_ID': 'experiment.probe',
            'NMR_spectrometer_probe_label': ('experiment.probe.name',self.getLabel),
            },
          },

        }
      }


    self.sfDict['NMR_spectrometer_expt'] = {                                                                  # *** Priority 1

      #'ccpnLoop':                       None,
      #'ccpnMap':                        None,

      #'title':                          None,

      'tags': {

        #'Sf_category':                  ['NMR_spectrometer_expt',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Optional - set in NmrStarExport
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Name':                         [None,lambda x = value: returnStarLine(x,length = 127),None,True],
        #'Type':                         [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Sample_volume':                [None,returnStarFloat,None,False],
        #'Sample_volume_units':          [None,lambda x = value: returnStarCode(x,length = 31),None,False],
        #'NMR_tube_type':                [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Sample_spinning_rate':         [None,returnStarFloat,None,False],
        #'Sample_angle':                 [None,returnStarFloat,None,False],
        #'NMR_spectrometer_ID':          [None,returnStarInt,'Sample_condition_list.ID',True],
        #'NMR_spectrometer_label':       [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
        #'NMR_spectrometer_probe_ID':    [None,returnStarInt,'NMR_spectrometer_probe.ID',True],
        #'NMR_spectrometer_probe_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
        #'Carrier_freq_switch_time':     [None,lambda x = value: returnStarString(x,length = 127),None,False],
        #'Software_ID':                  [None,returnStarInt,'Software.ID',True],
        #'Software_label':               [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
        #'Method_ID':                    [None,returnStarInt,'Method.ID',False],
        #'Method_label':                 [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
        #'Pulse_seq_accession_BMRB_code': [None,returnStarInt,None,False],
        #'Details':                      [None,returnStarString,None,False],                                  # *** Optional

        },

      'tables': {

        'NMR_experiment_file': {                                                                              # *** Includes optional fields

          'tags': {

            #'Name':                     [None,lambda x = value: returnStarLine(x,length = 127),None,True],   # *** Optional
            #'Type':                     [None,lambda x = value: returnStarLine(x,length = 127),None,True],   # *** Optional - option
            #'Directory_path':           [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Byte_order':               [None,lambda x = value: returnStarLine(x,length = 31),None,False],
            #'Bytes_per_data_point':     [None,returnStarInt,None,False],
            #'File_header_size':         [None,returnStarInt,None,False],
            #'Record_header_size':       [None,returnStarInt,None,False],
            #'Record_trailer_size':      [None,returnStarInt,None,False],
            #'Compression_algorithm':    [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Details':                  [None,returnStarString,None,False],                                  # *** Optional
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'NMR_spec_expt_ID':         [None,returnStarInt,'NMR_spec_expt.ID',True],

            },
          },

        'NMR_experiment_citation': {                                                                          # *** Includes optional fields

          'tags': {

            #'Citation_ID':              [None,returnStarInt,'Citation.ID',True],
            #'Citation_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'NMR_spec_expt_ID':         [None,returnStarInt,'NMR_spec_expt.ID',True],

            },
          },

        }
      }


    ##############
    # Resonances #
    ##############
    
    self.sfDict['resonance_linker'] = {                                                                       # *** Priority 1 - data - Not in ADIT-NMR?
    
      'ccpnLoop':                       ('nmrEntry',self.getAllResonanceParents),
      'ccpnMap':                         'resonanceParent',
      
      'title':                          ('resonanceParent.name',self.getTitle),

      'tags': {

        #'Details':                      'resonanceParent.details',       # TODO: this does not exist yet!

        },

      'tables': {

        'Resonance': {                                                                                        # *** Table
        
          'ccpnLoop':                   ('resonanceParent',self.getResonances),
          'ccpnMap':                     'resonance',

          'tags': {
                
            'Name':                      'resonance.name',
            'Resonance_set_ID':          'resonance.resonanceSet',   # *TODO*
            'Spin_system_ID':            'resonance.resonanceGroup', # Should these two lines be CCPN objects or IDs?
                               
            },
          }, 
       
        'Resonance_covalent_link': {                                                                          # *** Table
        
          #'ccpnLoop':                   [('resonanceParent',self.getResonances),
          #'ccpnMap':                    'resonance',

          'tags': {

            #.Resonance_ID_1
            #.Resonance_ID_2
                               
            },
          },

        'Resonance_assignment': {                                                                             # *** Table
          
          # TODO: this loop might have problems re-using an existing Atom_set_ID if it was already set earlier? Or OK if specifically set in tags section?
          'ccpnLoop':                    [self.getResonanceSets,'resonanceSet.sortedAtomSets()','atomSet.sortedAtoms()'], # ccpnLoop list
          'ccpnMap':                     [('resonanceSet','Resonance_set_ID'),('atomSet','Atom_set_ID'),('atom','Atom_ID')],

          # This is mechanism to ADD ccpnVars to the equation, so don't have to go through data model every time
          # Necessary to speed things up!!! Also ONLY works with ccpnLoop list as above!!!
          'addCcpnMap':                 {'atom': [('residue','atom.residue'),('chain','residue.chain'),('molResidue',('residue',self.getActualMolResidue) )]},

          'tags': {

            #'Resonance_set_ID':         'resonanceSet',
            #'Assembly_atom_ID':         None,  # TODO: can set this instead of all Label stuff, but then need to define all atoms in the molecular system in loop higher up
            'Entity_assembly_ID':        'chain',
            'Entity_ID':                 'chain.molecule',
            'Comp_index_ID':             'molResidue', # ADD SOMETHING HERE
            'Comp_ID':                  ('molResidue.chemCompVar.chemComp',self.getChemCompCif),
            'Atom_ID':                  ('atom',self.getPdbAtomName),
            'Atom_isotope_number':       'resonanceSet.findFirstResonance().isotope.massNumber',
            #'Atom_set_ID':              'atomSet',
                               
            },
          },
        
        'Spin_system': {                                                                                      # *** Table
        
          # TODO: some info from CCPN not transferred!!!
        
          'ccpnLoop':                    self.getResonanceGroups,
          'ccpnMap':                     'resonanceGroup',

          'tags': {
              
            'Entity_assembly_ID':        'resonanceGroup.residue.chain', # TODO: this could be multiple chains, in principle - how is that handled? Should this be part of PRIMARY?
            #'Entity_ID':                'resonanceGroup.residue.chain.molecule',
            'Comp_index_ID':            ('resonanceGroup.residue',self.getActualMolResidue), # Only relevant if residue set, so should be OK
            'Comp_ID':                  ('resonanceGroup',self.getResonanceGroupCompID),

            },
          },
          
        'Spin_system_link': {                                                                                 # *** Table

          'tags': {

            #'From_spin_system_ID':      [None,returnStarInt,None,False],
            #'To_spin_system_ID':        [None,returnStarInt,None,False],
            #'Offset':                   [None,returnStarInt,None,False],
            #'Type':                     [None,lambda x = value: returnStarString(x,length = 31),None,False],
            #'Selected':                 [None,lambda x = value: returnStarString(x,length = 3),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarString(x,length = 12),None,False],
            #'Resonance_linker_list_ID': [None,returnStarInt,None,False],

            },
          },

        }
      }


    ##############
    # Peak lists #
    ##############

    self.sfDict['spectral_peak_list'] = {                                                                     # *** Priority 1 - data

      'ccpnLoop':                        'nmrEntry.sortedPeakLists()',
      #'ccpnLoop':                       ('nmrEntry' self.getPeakLists),
      'ccpnMap':                         'peakList',

      'title':                          ('peakList',self.getPeakListTitle),

      'addCcpnMap':                     {'peakList': [('dataSource','peakList.dataSource'),('experiment','dataSource.experiment')]},

      'tags': {

        #'Sf_category':                  ['spectral_peak_list',lambda x = value: returnStarCode(x,length = 31),None,False],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],

        'Data_file_name':               ('peakList.root.name',self.getProjectFileName),                       # *** Mandatory
        'Sample_ID':                     'experiment.sample',
        'Sample_label':                 ('experiment.sample.name',self.getLabel),

        'Sample_condition_list_ID':      'experiment.sampleConditionSet',
        'Sample_condition_list_label':  ('experiment.sampleConditionSet.name',self.getLabel),

        'Experiment_ID':                 'experiment',
        'Experiment_name':              ('experiment',self.getStarExperimentName),                            # *** Mandatory - pulldown
        'Details':                       'peakList.details',                                                  # *** Optional
        'Number_of_spectral_dimensions': 'dataSource.numDim',

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'Spectral_dim': {                                                                                     # *** Includes mandatory fields

          # TODO - should multiple isotopes be possible, we need to use dataDimRefs to loop over.

          'ccpnLoop':                    ('dataSource',self.getDataDimsFromDataSource),
          'ccpnMap':                     'freqDataDim',

          'addCcpnMap':                 {'freqDataDim': [('dataDimRef', 'freqDataDim.findFirstDataDimRef()')]},

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],                                      # *** Mandatory

            'Atom_type':                ('dataDimRef.expDimRef.isotopes', self.getFirstIsotopeName),          # *** Mandatory
            'Atom_isotope_number':      ('dataDimRef.expDimRef.isotopes', self.getFirstIsotopeNumber),        # *** Mandatory
            'Spectral_region':          ('dataDimRef.expDimRef.isotopeCodes', self.getFirstIsotopeCode),      # *** Mandatory
            #'Magnetization_linkage_ID': 'experiment.expTransfer',                                            # *** Optional
            'Sweep_width':               'freqDataDim.spectralWidth',                                         # *** Mandatory

            #'Encoding_code':            [None,lambda x = value: returnStarLine(x,length = 31),None,False],   # *** Optional
            #'Encoded_source_dimension_ID': [None,returnStarInt,None,False],                                  # *** Optional
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Spectral_peak_list_ID':    [None,returnStarInt,'Spectral_peak_list.ID',True],

            },
          },

        'Spectral_peak_software': {                                                                           # *** Includes optional fields

          'ccpnLoop':                   ('peakList', self.getPeakMethods),
          'ccpnMap':                     'method',

          'tags': {

            'Software_ID':               'method.software',
            'Software_label':           ('method.software.name',self.getTitle),                               # *** Optional
            'Method_ID':                 'method',
            'Method_label':             ('method.name',self.getTitle),
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Spectral_peak_list_ID':    [None,returnStarInt,'Spectral_peak_list.ID',True],

            },
          },

        'Peak': {                                                                                             # *** Table

          'ccpnLoop':                    'peakList.sortedPeaks()',
          'ccpnMap':                     'peak',

          'tags': {

            'Figure_of_merit':           'peak.figOfMerit',
            'Details':                   'peak.details',

            },
          },

        'Peak_general_char': {                                                                                # *** Table

          'ccpnLoop':                    ['peakList.sortedPeaks()','peak.sortedPeakIntensities()'],
          'ccpnMap':                     [('peak','Peak_ID'),('peakIntensity','Intensity_val')],

          'tags': {

            'Intensity_val':             'peakIntensity.value',
            'Intensity_val_err':         'peakIntensity.error',
            'Measurement_method':        'peakIntensity.intensityType',

            },
          },

        'Peak_char': {                                                                                        # *** Table

          'ccpnLoop':                    ['peakList.sortedPeaks()', 'peak.sortedPeakDims()', 'peak.sortedPeakContribs()'],
          'ccpnMap':                     [('peak','Peak_ID'),('peakDim',''),('peakContrib','Peak_contribution_ID')],

          'tags': {

            'Spectral_dim_ID':           'peakDim.dataDimRef.dataDim',
            #'Peak_contribution_ID':      'peakContrib',
            'Chem_shift_val':           ('peakDim.value',self.getValue),
            'Chem_shift_val_err':       ('peakDim.valueError',self.getErrorValue),
            'Decay_rate_val':            'peakDim.decayRate',
            'Decay_rate_val_err':        'peakDim.decayRateError',
            'Phase_val':                 'peakDim.phase',
            'Phase_val_err':             'peakDim.phaseError',
            'Line_width_val':            'peakDim.lineWidth',
            #'Line_width_val_err':       # not in CCPN

            },
          },

        'Peak_contribution': {                                                                                # *** Table

          'ccpnLoop':                   ('peakList', self.getPeakContribs),
          'ccpnMap':                     'peakContrib',

          'tags': {

            #'ID':                       'peakContrib.serial',
            'Contribution_fractional_val': 'peakContrib.weight',

            },
          },

        'Assigned_peak_chem_shift': {                                                                         # *** Table

          #'ccpnLoop':                    ['peakList.sortedPeaks()', 'peak.sortedPeakDims()', 'peak.sortedPeakContribs()', 'peakContrib.sortedPeakDimContribs()'],
          #'ccpnMap':                     [('peak', 'Peak_ID'),('peakDim', 'Spectral_dim_ID'),('peakContrib','Peak_contribution_ID'),('peakDimContrib','Set_ID'),],

          'ccpnLoop':                    ['peakList.sortedPeaks()', 'peak.sortedPeakDims()', 'peakDim.sortedPeakDimContribs()', 'peakDimContrib.sortedPeakContribs()'],
          'ccpnMap':                     [('peak', 'Peak_ID'),('peakDim', 'Spectral_dim_ID'),('peakDimContrib','Set_ID'),('peakContrib','Peak_contribution_ID')],

          #'addCcpnMap':                {'peakDim': [('resonance', 'peakDim.findFirstPeakDimContrib().resonance')]},
          'addCcpnMap':                 {'peakDimContrib': [('resonance', 'peakDimContrib.resonance')]},

          'tags': {

            'Spectral_dim_ID':           'peakDim.dataDimRef.dataDim',
            #'Peak_contribution_ID':      'peakContrib',
            'Set_ID':                    'peakDimContrib.serial',
            'Val':                      ('peakDim.value',self.getValue),
            #'Figure_of_merit':          'peakContrib.weight',
            'Resonance_ID':              'resonance', #'peakDimContrib.resonance',
            'Atom_chem_shift_ID':        None,

            'CUSTOM_setAtoms':           'resonance', #'peakDimContrib.resonance',

            },
          },

        }
      }


    ##############################
    # Chemical shift referencing #
    ##############################

    self.sfDict['chem_shift_reference'] = {                                                                   # *** Priority 1 - data

      'ccpnLoop':                       ('nmrEntry',self.getChemShiftRefLists),
      'ccpnMap':                         'chemShiftRefList',

      'title':                          ('chemShiftRefList',self.getChemShiftRefTitle),

      'tags': {

        #'Sf_category':                  ['chem_shift_reference',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - option - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                            'chemShiftRefList', #[None,returnStarInt,None,True],
        'Proton_shifts_flag':           ('chemShiftRefList',self.getProtonRefFlag),                           # *** Mandatory - radio
        'Carbon_shifts_flag':           ('chemShiftRefList',self.getCarbonRefFlag),                           # *** Mandatory - radio
        'Nitrogen_shifts_flag':         ('chemShiftRefList',self.getNitrogenRefFlag),                         # *** Mandatory - radio
        'Phosphorus_shifts_flag':       ('chemShiftRefList',self.getPhosphorusRefFlag),                       # *** Mandatory - radio - auto
        'Other_shifts_flag':            ('chemShiftRefList',self.getOtherRefFlag),                            # *** Mandatory - bool - auto
        #'Details':                      [None,returnStarString,None,False],                                  # *** Optional

        },

      'tables': {

        'Chem_shift_ref': {                                                                                   # *** Table

          'ccpnLoop':                   ('chemShiftRefList',self.getChemShiftRefs),
          'ccpnMap':                     'chemShiftRef',

          'tags': {

            'Atom_type':                 'chemShiftRef.isotope.chemElement.symbol',
            'Atom_isotope_number':       'chemShiftRef.isotope.massNumber',
            'Mol_common_name':           'chemShiftRef.molName',
            'Atom_group':                'chemShiftRef.atomGroup',
            #'Concentration_val':        [None,returnStarFloat,None,False],
            #'Concentration_units':      [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Solvent':                  [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Rank':                     [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            'Chem_shift_units':          'chemShiftRef.unit',
            'Chem_shift_val':           ('chemShiftRef.value',self.getChemShiftRefValue),
            'Ref_method':               ('chemShiftRef',self.getStarReferenceType),
            'Ref_type':                  'chemShiftRef.referenceType',
            'Indirect_shift_ratio':      'chemShiftRef.indirectShiftRatio',
            'External_ref_loc':         ('chemShiftRef',self.getChemShiftLocation),
            'External_ref_sample_geometry': ('chemShiftRef',self.getChemShiftGeometry),
            'External_ref_axis':        ('chemShiftRef',self.getChemShiftAxis),
#            'Indirect_shift_ratio_cit_ID': None, #[None,returnStarInt,'Citation.ID',False],
            #'Indirect_shift_ratio_cit_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Ref_correction_type':      [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Correction_val':           [None,returnStarFloat,None,False],
#            'Correction_val_cit_ID':     None, #[None,returnStarInt,'Citation.ID',False],
            #'Correction_val_cit_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Chem_shift_reference_ID':  ('chemShiftRefList',self.getChemShiftRefId),

            },
          },
        }
      }


    ################
    # Measurements #
    ################

    self.sfDict['assigned_chemical_shifts'] = {                                                               # *** Priority 1 - data

      'ccpnLoop':                        "nmrEntry.findAllMeasurementLists(className='ShiftList', isSimulated=False)",
      'ccpnMap':                         'shiftList',
      
      'title':                          ('shiftList',self.getMeasurementTitle),

      'tags': {
      
        # Warning: did not put Labels in! IDs should be fine.

        #'Sf_category':                  ['assigned_chemical_shifts',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Data_file_name':               ('shiftList.root.name',self.getProjectFileName),                      # *** Mandatory
        'Sample_condition_list_ID':      'shiftList.findFirstExperiment().sampleConditionSet',
        'Sample_condition_list_label':   'shiftList.findFirstExperiment().sampleConditionSet.name',           # *** Mandatory - option

        'Chem_shift_reference_ID':      ('shiftList',self.getShiftReferenceListId),
        'Chem_shift_reference_label':   ('shiftList',self.getShiftReferenceListTitle),                        # *** Mandatory - option

        'Chem_shift_1H_err':             "shiftList.findFirstApplicationData(application='nmrStar',keyword='chemShiftErr1H').value", # *** Optional
        'Chem_shift_13C_err':            "shiftList.findFirstApplicationData(application='nmrStar',keyword='chemShiftErr13C').value", # *** Optional
        'Chem_shift_15N_err':            "shiftList.findFirstApplicationData(application='nmrStar',keyword='chemShiftErr15N').value", # *** Optional
        'Chem_shift_31P_err':            "shiftList.findFirstApplicationData(application='nmrStar',keyword='chemShiftErr31P').value", # *** Optional
        'Chem_shift_2H_err':             "shiftList.findFirstApplicationData(application='nmrStar',keyword='chemShiftErr2H').value", # *** Optional
        'Chem_shift_19F_err':            "shiftList.findFirstApplicationData(application='nmrStar',keyword='chemShiftErr19F').value", # *** Optional
        #'Error_derivation_method':      [None,returnStarString,None,False],                                  # *** Optional
        'Details':                       'shiftList.details',                                                 # *** Optional

        # TODO: these would need to be filled in if list can't be handled.
        #   - Could put this in applicationData text for empty shiftList?
        # Also - other measurements with these fields.

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'Chem_shift_experiment':         self.setExperimentsLoop('shiftList.sortedExperiments()','experiment'), # *** Table

        #'Chem_shift_experiment': {                                                                           # *** Includes mandatory fields
        
        #  'ccpnLoop':                   'shiftList.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          ('experiment',self.getStarExperimentName),                            # *** Mandatory - pulldown
            #'Sample_ID':                 'experiment.sample',
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Sample_state':              'experiment.sampleState',
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Assigned_chem_shift_list_ID': [None,returnStarInt,'Assigned_chem_shift_list.ID',True],

        #    },
        #  }, 

        'Systematic_chem_shift_offset': {                                                                     # *** Includes optional fields

          'tags': {

            #'Type':                     [None,lambda x = value: returnStarLine(x,length = 127),None,True],   # *** Optional - option
            #'Atom_type':                [None,lambda x = value: returnStarCode(x,length = 15),None,True],    # *** Optional - option
            #'Atom_isotope_number':      [None,returnStarInt,None,False],                                     # *** Optional

            #'Val':                      [None,lambda x = value: returnStarString(x,length = 15),None,True],  # *** Optional
            #'Val_err':                  [None,lambda x = value: returnStarString(x,length = 15),None,False], # *** Optional

            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Assigned_chem_shift_list_ID': [None,returnStarInt,'Assigned_chem_shift_list.ID',True],

            },
          },

        'Chem_shift_software':           self.setSoftwareLoop('shiftList'),                                   # *** Table

        #'Chem_shift_software': {                                                                             # *** Includes mandatory fields

        #  'ccpnMap':                    'shiftList',

        #  'tags': {

            #'Software_ID':              'shiftList.method.software', #[None,returnStarInt,'Software.ID',True],
            #'Software_label':           'shiftList.method.software.name',                                    # *** Mandatory
            #'Method_ID':                'shiftList.method', #[None,returnStarInt,'Method.ID',False],
            #'Method_label':             'shiftList.method.name',
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Assigned_chem_shift_list_ID': [None,returnStarInt,'Assigned_chem_shift_list.ID',True],

        #    },
        #  },

        'Atom_chem_shift': {                                                                                  # *** Table
        
          'ccpnLoop':                   ('shiftList.sortedMeasurements()',self.getMeasurementsByIndividualAtom),
          'ccpnMap':                     'measurementByIndividualAtom',

          'tags': {
                
            #'ID':                       'shiftByIndividualAtom.serial', # TODO this is a hack - can this be done better?
            'Atom_isotope_number':       'measurementByIndividualAtom.resonance.isotope.massNumber',
            'Val':                      ('measurementByIndividualAtom.measurement.value', self.getValue),
            'Val_err':                  ('measurementByIndividualAtom.measurement.error', self.getErrorValue),
            #'Assign_fig_of_merit':      'measurementByIndividualAtom.measurement.figOfMerit',
            'Ambiguity_code':           ('measurementByIndividualAtom',self.getShiftAmbiguityCode),
            #'Occupancy':
            'Resonance_ID':              'measurementByIndividualAtom.resonance',
            'Details':                   'measurementByIndividualAtom.measurement.details',

            'CUSTOM_setAtoms_individual': 'measurementByIndividualAtom'

            },
          }, 

        # TODO: is this necessary? Should be taken care of by Resonance, no?

        'Ambiguous_atom_chem_shift': {                                                                        # *** Table

          'tags': {

            #'Ambiguous_shift_set_ID':   [None,returnStarInt,None,True],
            #'Atom_chem_shift_ID':       [None,returnStarInt,'Atom_chem_shift.ID',True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Assigned_chem_shift_list_ID': [None,returnStarInt,'Assigned_chem_shift_list.ID',True],

            },
          },

      
        }
      }


    self.sfDict['coupling_constants'] = {                                                                     # *** Priority 1 - data

      'ccpnLoop':                        "nmrEntry.findAllMeasurementLists(className='JCouplingList')",
      'ccpnMap':                         'jCouplingList',

      'title':                          ('jCouplingList',self.getMeasurementTitle),

      'tags': {
      
        #'Sf_category':                  ['coupling_constants',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Data_file_name':               ('jCouplingList.root.name',self.getProjectFileName),                  # *** Mandatory
        'Sample_condition_list_ID':      'jCouplingList.findFirstExperiment().sampleConditionSet',
        'Sample_condition_list_label':   'jCouplingList.findFirstExperiment().sampleConditionSet.name',       # *** Mandatory - option
        'Spectrometer_frequency_1H':     'jCouplingList.sf',                                                  # *** Mandatory - option
        'Details':                       'jCouplingList.details',                                             # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'Coupling_constant_experiment':  self.setExperimentsLoop('jCouplingList.sortedExperiments()','experiment'), # *** Table

        #'Coupling_constant_experiment': {                                                                    # *** Includes mandatory fields
        
        #  'ccpnLoop':                   'jCouplingList.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',False],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Coupling_constant_list_ID': [None,returnStarInt,'Coupling_constant_list.ID',True],

        #    },
        #  }, 

        'Coupling_constant_software':    self.setSoftwareLoop('jCouplingList'),                               # *** Table

        #'Coupling_constant_software': {                                                                      # *** Includes optional fields

        #  'ccpnMap':                    'jCouplingList',

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Coupling_constant_list_ID': [None,returnStarInt,'Coupling_constant_list.ID',True],

        #    },
        #  },

        'Coupling_constant': {                                                                                # *** Table

          'ccpnLoop':                    'jCouplingList.sortedMeasurements()',
          'ccpnMap':                     'measurement',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Code':                     [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Assembly_atom_ID_1':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_1':     [None,returnStarInt,'Entity_assembly.ID',False],
            #'Entity_ID_1':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_1':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_1':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
            #'Comp_ID_1':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_1':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_1':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number_1':    [None,returnStarInt,None,False],
            #'Ambiguity_code_1':         [None,lambda x = value: returnStarCode(x,length = 127),None,False],
            #'Assembly_atom_ID_2':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_2':     [None,returnStarInt,'Entity_assembly.ID',False],
            #'Entity_ID_2':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_2':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_2':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_2':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_2':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_2':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number_2':    [None,returnStarInt,None,False],
            #'Ambiguity_code_2':         [None,lambda x = value: returnStarCode(x,length = 127),None,False],

            'Val':                       'measurement.value',
            'Val_err':                   'measurement.error',

            #'Val_min':                  [None,returnStarFloat,None,False],
            #'Val_max':                  [None,returnStarFloat,None,False],
            #'Resonance_ID_1':           [None,returnStarInt,'Resonance.ID',False],
            #'Resonance_ID_2':           [None,returnStarInt,'Resonance.ID',False],
            #'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID_1':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID_2':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],

            'Details':                   'measurement.details',

            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Coupling_constant_list_ID': [None,returnStarInt,'Coupling_constant_list.ID',True],

            'CUSTOM_setAtoms':           'measurement.sortedResonances()',

            },
          },

        }
      }


    self.sfDict['chem_shift_anisotropy'] = {                                                                  # *** Priority 1 - data
    
      'ccpnLoop':                        "nmrEntry.findAllMeasurementLists(className='ShiftAnisotropyList')",
      'ccpnMap':                         'shiftAnisotropyList',
      
      'title':                          ('shiftAnisotropyList',self.getMeasurementTitle),

      'tags': {
      
        #'Sf_category':                  ['chem_shift_anisotropy',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Data_file_name':               ('shiftAnisotropyList.root.name',self.getProjectFileName),            # *** Mandatory
        'Sample_condition_list_ID':      'shiftAnisotropyList.findFirstExperiment().sampleConditionSet',
        'Sample_condition_list_label':   'shiftAnisotropyList.findFirstExperiment().sampleConditionSet.name', # *** Mandatory - option
        'Spectrometer_frequency_1H':     'shiftAnisotropyList.sf',
        'Val_units':                     'shiftAnisotropyList.unit',                                          # *** Mandatory - option
        'Details':                       'shiftAnisotropyList.details',                                       # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'CS_anistropy_experiment':       self.setExperimentsLoop('shiftAnisotropyList.sortedExperiments()','experiment'), # *** Table

        #'CS_anistropy_experiment': {                                                                         # *** Includes mandatory fields
        
        #  'ccpnLoop':                   'shiftAnisotropyList.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Chem_shift_anisotropy_ID': [None,returnStarInt,'Chem_shift_anisotropy.ID',True],

        #    },
        #  }, 

        'CS_anisotropy_software':        self.setSoftwareLoop('shiftAnisotropyList'),                         # *** Table

        #'CS_anisotropy_software': {                                                                          # *** Includes optional fields

        #  'ccpnMap':                    'shiftAnisotropyList',

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',True],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Chem_shift_anisotropy_ID': [None,returnStarInt,'Chem_shift_anisotropy.ID',True],

        #    },
        #  },

        'CS_anisotropy': {                                                                                    # *** Table
        
          'ccpnLoop':                   ('shiftAnisotropyList.sortedMeasurements()',self.getMeasurementsByIndividualAtom),
          'ccpnMap':                     'measurementByIndividualAtom',

          'tags': {
                
            #'ID':                       'shiftByIndividualAtom.serial', # TODO this is a hack - can this be done better?
            'Atom_isotope_number':       'measurementByIndividualAtom.resonance.isotope.massNumber',
            'Val':                       'measurementByIndividualAtom.measurement.value',
            'Val_err':                   'measurementByIndividualAtom.measurement.error',
            #'Principal_value_sigma_11_val':
            #'Principal_value_sigma_22_val':
            #'Principal_value_sigma_33_val':
            #'Principal_Euler_angle_alpha_val':
            #'Principal_Euler_angle_beta_val':
            #'Principal_Euler_angle_gamma_val':
            #'Occupancy':
            'Resonance_ID':              'measurementByIndividualAtom.resonance',
            'Details':                   'measurementByIndividualAtom.measurement.details',

            'CUSTOM_setAtoms_individual': 'measurementByIndividualAtom'

            },
          }, 

        }
      }


    self.sfDict['theoretical_chem_shifts'] = {                                                                # *** Priority 3 - data

      'ccpnLoop':                        "nmrEntry.findAllMeasurementLists(className='ShiftList', isSimulated=True)",
      'ccpnMap':                         'shiftList',
      
      'title':                          ('shiftList',self.getMeasurementTitle),

      'tags': {
      
        # Warning: did not put Labels in! IDs should be fine.

        #'Sf_category':                  ['theoretical_chem_shifts',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory
        #'Chem_shifts_calc_type_ID':     [None,returnStarInt,'Chem_shifts_calc_type.ID',True],
        #'Chem_shifts_calc_type_label':  [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory
        #'Model_atomic_coordinates_ID':  [None,returnStarInt,'Representative_conformer.ID',True],
        #'Model_atomic_coordinates_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Optional
        #'Chem_shielding_tensor_list_ID': [None,returnStarInt,'Shielding_tensor_list.ID',True],
        #'Chem_shielding_tensor_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Optional
        #'Fermi_contact_spin_density_units': [None,lambda x = value: returnStarCode(x,length = 127),None,False], # *** Optional
        #'Chem_shift_1H_err':            [None,returnStarFloat,None,False],
        #'Chem_shift_2H_err':            [None,returnStarFloat,None,False],
        #'Chem_shift_13C_err':           [None,returnStarFloat,None,False],
        #'Chem_shift_15N_err':           [None,returnStarFloat,None,False],
        #'Chem_shift_19F_err':           [None,returnStarFloat,None,False],
        #'Chem_shift_31P_err':           [None,returnStarFloat,None,False],
      
        'Details':                       'shiftList.details',                                                 # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'Theoretical_chem_shift': {                                                                           # *** Table
        
          'ccpnLoop':                   ('shiftList.sortedMeasurements()',self.getMeasurementsByIndividualAtom),
          'ccpnMap':                     'measurementByIndividualAtom',

          'tags': {
                
            #'ID':                       'shiftByIndividualAtom.serial', # TODO this is a hack - can this be done better?
            'Atom_isotope_number':       'measurementByIndividualAtom.resonance.isotope.massNumber',
            'Val':                       'measurementByIndividualAtom.measurement.value',
            'Val_err':                   'measurementByIndividualAtom.measurement.error',
            #'Fermi_contact_spin_density': [None,lambda x = value: returnStarString(x,length = 127),None,False],
            #'Occupancy':
            #'Details':                  'measurementByIndividualAtom.measurement.details',

            'CUSTOM_setAtoms_individual': 'measurementByIndividualAtom'

            },
          }, 
      
        }
      }


    self.sfDict['Chem_shifts_calc_type'] = {                                                                  # *** Priority 4 - data

      #'ccpnLoop':                       None,
      #'ccpnMap':                        None,

      #'title':                          None,

      'tags': {

        #'Sf_category':                  ['chem_shifts_calc_type',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Calculation_level':            [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory - option
        #'Quantum_mechanical_method':    [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory - option
        #'Quantum_mechanical_theory_level': [None,lambda x = value: returnStarLine(x,length = 127),None,False], # *** Mandatory - option
        #'Quantum_mechanical_basis_set': [None,lambda x = value: returnStarLine(x,length = 31),None,False],   # *** Mandatory
        #'Chem_shift_nucleus':           [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory - option
        #'Modeled_sample_cond_list_ID':  [None,returnStarInt,'Sample_condition_list.ID',True],
        #'Modeled_sample_cond_list_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - option
        #'Chem_shift_reference_ID':      [None,returnStarInt,'Chem_shift_reference.ID',True],
        #'Chem_shift_reference_label':   [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - option
        #'Details':                      [None,returnStarString,None,False],                                  # *** Optional

        },

      'tables': {

        'Chem_shifts_calc_software': {                                                                        # *** Includes optional fields

          'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Chem_shifts_calc_type_ID': [None,returnStarInt,'Chem_shifts_calc_type.ID',True],

            },
          },

        }
      }


    self.sfDict['RDCs'] = {                                                                                   # *** Priority 2 - data

      'ccpnLoop':                        "nmrEntry.findAllMeasurementLists(className='RdcList')",
      'ccpnMap':                         'rdcList',
      
      'title':                          ('rdcList',self.getMeasurementTitle),

      'tags': {

        #'Sf_category':                  ['RDCs',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Data_file_name':               ('rdcList.root.name',self.getProjectFileName),                        # *** Mandatory
        'Sample_condition_list_ID':      'rdcList.findFirstExperiment().sampleConditionSet',
        'Sample_condition_list_label':   'rdcList.findFirstExperiment().sampleConditionSet.name',             # *** Mandatory - option
        'Spectrometer_frequency_1H':     'rdcList.sf',                                                        # *** Mandatory
        'Details':                       'rdcList.details',                                                   # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'RDC_experiment':                self.setExperimentsLoop('rdcList.sortedExperiments()','experiment'), # *** Table

        #'RDC_experiment': {                                                                                  # *** Includes mandatory fields
        
        #  'ccpnLoop':                   'rdcList.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'RDC_list_ID':              [None,returnStarInt,'RDC_list.ID',True],

        #    },
        #  },

        'RDC_software':                  self.setSoftwareLoop('rdcList'),                                     # *** Table

        #'RDC_software': {                                                                                    # *** Includes optional fields

        #  'ccpnMap':                    'rdcList',

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'RDC_list_ID':              [None,returnStarInt,'RDC_list.ID',True],

        #    },
        #  },

        'RDC': {                                                                                              # *** Table

          'ccpnLoop':                    'rdcList.sortedMeasurements()',
          'ccpnMap':                     'measurement',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'RDC_code':                 [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Assembly_atom_ID_1':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_1':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_1':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_1':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_1':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_1':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_1':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_1':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number_1':    [None,returnStarInt,None,False],
            #'Ambiguity_code_1':         [None,lambda x = value: returnStarCode(x,length = 127),None,False],
            #'Assembly_atom_ID_2':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_2':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_2':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_2':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_2':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_2':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_2':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_2':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number_2':    [None,returnStarInt,None,False],
            #'Ambiguity_code_2':         [None,lambda x = value: returnStarCode(x,length = 127),None,False],

            'Val':                       'measurement.value',
            'Val_err':                   'measurement.error',

            #'Val_min':                  [None,returnStarFloat,None,False],
            #'Val_max':                  [None,returnStarFloat,None,False],
            #'Resonance_ID_1':           [None,returnStarInt,'Resonance.ID',False],
            #'Resonance_ID_2':           [None,returnStarInt,'Resonance.ID',False],
            #'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID_1':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID_2':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'RDC_list_ID':              [None,returnStarInt,'RDC_list.ID',True],

            'CUSTOM_setAtoms':           'measurement.sortedResonances()',

            },
          },

        }
      }


    self.sfDict['spectral_density_values'] = {                                                                # *** Priority 3 - data
    
      #'ccpnLoop':                       "nmrEntry.findAllDerivedDataLists(className='SpectralDensityList')",
      'ccpnLoop':                       ('nmrEntry', self.getAllSpectralDensityDerivations),
      'ccpnMap':                         'spectralDensityDerivation',

      'title':                          ('spectralDensityDerivation.details',self.getTitle),

      'tags': {

        #'Sf_category':                  ['spectral_density_values',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory
        #'Sample_condition_list_ID':     [None,returnStarInt,'Sample_condition_list.ID',True],
        #'Sample_condition_list_label':  [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Mandatory - option

        'Details':                       'spectralDensityDerivation.details',                                 # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'Spectral_density_experiment':   self.setExperimentsLoop('spectralDensityDerivation.parent.sortedExperiments()','experiment'), # *** Table

        #'Spectral_density_experiment': {                                                                     # *** Includes mandatory fields
        
        #  'ccpnLoop':                   'spectralDensityDerivation.parent.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional - option
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],    # *** Optional - option
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Spectral_density_list_ID': [None,returnStarInt,'Spectral_density_list.ID',True],

        #    },
        #  },

        'Spectral_density_software':     self.setSoftwareLoop('spectralDensityDerivation.parent'),            # *** Table

        #'Spectral_density_software': {                                                                        # *** Includes optional fields

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',True],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Spectral_density_list_ID': [None,returnStarInt,'Spectral_density_list.ID',True],

        #    },
        #  },

        'Spectral_density': {                                                                                 # *** Table

          'ccpnLoop':                    'spectralDensityDerivation.sortedDerivedData()',
          'ccpnMap':                     'datum',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Assembly_atom_ID':         [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID':       [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID':                [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID':            [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID':                   [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID':                  [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type':                [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number':      [None,returnStarInt,None,False],

            #'W_zero_val':               [None,returnStarFloat,None,False],
            #'W_zero_val_err':           [None,returnStarFloat,None,False],
            #'W_1H_val':                 [None,returnStarFloat,None,False],
            #'W_1H_val_err':             [None,returnStarFloat,None,False],
            #'W_13C_val':                [None,returnStarFloat,None,False],
            #'W_13C_val_err':            [None,returnStarFloat,None,False],
            #'W_15N_val':                [None,returnStarFloat,None,False],
            #'W_15N_val_err':            [None,returnStarFloat,None,False],

            #'Resonance_ID':             [None,returnStarInt,'Resonance.ID',False],
            #'Auth_entity_assembly_ID':  [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID':              [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID':             [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID':             [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Spectral_density_list_ID': [None,returnStarInt,'Spectral_density_list.ID',True],

            'CUSTOM_setAtoms':           'datum.sortedResonances()',

            },
          }, 

        }
      }


    self.sfDict['other_data_types'] = {                                                                       # *** Priority 3 - data

      #'ccpnLoop':                       "nmrEntry.findAllDerivedDataLists(className='DataList')",
      'ccpnLoop':                       ('nmrEntry', self.getAllDataDerivations),
      'ccpnMap':                         'dataDerivation',

      'title':                          ('dataDerivation.details',self.getTitle),

      'tags': {

        #'Sf_category':                  ['other_data_types',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],

        'Name':                          'dataDerivation.name',
        'Definition':                    'dataDerivation.definition',

        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory
        #'Sample_condition_list_ID':     [None,returnStarInt,'Sample_condition_list.ID',True],
        #'Sample_condition_list_label':  [None,lambda x = value: returnStarLabel(x,length = 127),None,True],

        'Details':                       'dataDerivation.details',                                            # *** Mandatory

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'Other_data_experiment':         self.setExperimentsLoop('dataDerivation.parent.sortedExperiments()','experiment'), # *** Table

        #'Other_data_experiment': {                                                                           # *** Includes mandatory fields
        
        #  'ccpnLoop':                   'dataDerivation.parent.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Other_data_type_list_ID':  [None,returnStarInt,'Other_data_type_list.ID',True],

        #    },
        #  }, 

        'Other_data_software':           self.setSoftwareLoop('dataDerivation.parent'),                       # *** Table

        #'Other_data_software': {                                                                              # *** Includes optional fields

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Other_data_type_list_ID':  [None,returnStarInt,'Other_data_type_list.ID',True],

        #    },
        #  },

        'Other_data': {                                                                                       # *** Table

          'ccpnLoop':                    'dataDerivation.sortedDerivedData()',
          'ccpnMap':                     'datum',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Assembly_atom_ID':         [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID':       [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID':                [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID':            [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID':                   [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID':                  [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type':                [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number':      [None,returnStarInt,None,False],

            'Val':                       'datum.value',
            'Val_err':                   'datum.error',

            #'Resonance_ID':             [None,returnStarInt,'Resonance.ID',False],
            #'Auth_entity_assembly_ID':  [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID':              [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID':             [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID':             [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Other_data_type_list_ID':  [None,returnStarInt,'Other_data_type_list.ID',True],

            'CUSTOM_setAtoms':           'datum.sortedResonances()',

            },
          }, 

        }
      }


    self.sfDict['H_exch_rates'] = {                                                                           # *** Priority 2 - data
    
      'ccpnLoop':                        "nmrEntry.findAllMeasurementLists(className='HExchRateList')",
      'ccpnMap':                         'hExchRateList',
      
      'title':                          ('hExchRateList',self.getMeasurementTitle),

      'tags': {

        #'Sf_category':                  ['H_exch_rates',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Data_file_name':               ('hExchRateList.root.name',self.getProjectFileName),                  # *** Mandatory
        'Sample_condition_list_ID':      'hExchRateList.findFirstExperiment().sampleConditionSet',
        'Sample_condition_list_label':   'hExchRateList.findFirstExperiment().sampleConditionSet.name',       # *** Mandatory - option
        'Val_units':                     'hExchRateList.unit',                                                # *** Mandatory - option
        'Details':                       'hExchRateList.details',                                             # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'H_exch_rate_experiment':        self.setExperimentsLoop('hExchRateList.sortedExperiments()','experiment'), # *** Table

        #'H_exch_rate_experiment': {                                                                          # *** Includes mandatory fields
        
        #  'ccpnLoop':                   'HExchRateList.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'H_exch_rate_list_ID':      [None,returnStarInt,'H_exch_rate_list.ID',True],

        #    },
        #  }, 

        'H_exch_rate_software':          self.setSoftwareLoop('hExchRateList'),                               # *** Table

        #'H_exch_rate_software': {                                                                            # *** Includes optional fields

        #  'ccpnMap':                    'hExchRateList',

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'H_exch_rate_list_ID':      [None,returnStarInt,'H_exch_rate_list.ID',True],

        #    },
        #  },

        'H_exch_rate': {                                                                                      # *** Table
        
          'ccpnLoop':                   ('hExchRateList.sortedMeasurements()',self.getMeasurementsByIndividualAtom),
          'ccpnMap':                     'measurementByIndividualAtom',

          'tags': {
                
            #'ID':                       'shiftByIndividualAtom.serial', # TODO this is a hack - can this be done better?
            'Atom_isotope_number':       'measurementByIndividualAtom.resonance.isotope.massNumber',
            'Val':                       'measurementByIndividualAtom.measurement.value',
            'Val_err':                   'measurementByIndividualAtom.measurement.error',
            'Resonance_ID':              'measurementByIndividualAtom.resonance',
            'Details':                   'measurementByIndividualAtom.measurement.details',

            'CUSTOM_setAtoms_individual': 'measurementByIndividualAtom'

            },
          }, 

        }
      }


    self.sfDict['H_exch_protection_factors'] = {                                                              # *** Priority 4 - data
    
      'ccpnLoop':                        "nmrEntry.findAllMeasurementLists(className='HExchProtectionList')",
      'ccpnMap':                         'hExchProtectionList',
      
      'title':                          ('hExchProtectionList',self.getMeasurementTitle),

      'tags': {

        #'Sf_category':                  ['H_exch_protection_factors',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Data_file_name':               ('hExchProtectionList.root.name',self.getProjectFileName),            # *** Mandatory
        'Sample_condition_list_ID':      'hExchProtectionList.findFirstExperiment().sampleConditionSet',
        'Sample_condition_list_label':   'hExchProtectionList.findFirstExperiment().sampleConditionSet.name', # *** Mandatory - option
        #'Std_values_source_cit_ID':     [None,returnStarInt,'Citation.ID',False],
        #'Std_values_source_cit_label':  [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Optional
        'Details':                       'hExchProtectionList.details',                                       # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'H_exch_protection_fact_experiment': self.setExperimentsLoop('hExchProtectionList.sortedExperiments()','experiment'), # *** Table

        #'H_exch_protection_fact_experiment': {                                                               # *** Includes mandatory fields
        
        #  'ccpnLoop':                   'hExchProtectionList.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'H_exch_protection_factor_list_ID': [None,returnStarInt,'H_exch_protection_factor_list.ID',True],

        #    },
        #  }, 

        'H_exch_protect_fact_software':  self.setSoftwareLoop('hExchProtectionList'),                         # *** Table

        #'H_exch_protect_fact_software': {                                                                    # *** Includes optional fields

        #  'ccpnMap':                    'hExchProtectionList',

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'H_exch_protection_factor_list_ID': [None,returnStarInt,'H_exch_protection_factor_list.ID',True],

        #    },
        #  },

        'H_exch_protection_factor': {                                                                         # *** Table
        
          'ccpnLoop':                   ('hExchProtectionList.sortedMeasurements()',self.getMeasurementsByIndividualAtom),
          'ccpnMap':                     'measurementByIndividualAtom',

          'tags': {
                
            #'ID':                       'shiftByIndividualAtom.serial', # TODO this is a hack - can this be done better?
            'Atom_isotope_number':       'measurementByIndividualAtom.resonance.isotope.massNumber',
            'Val':                       'measurementByIndividualAtom.measurement.value',
            'Val_err':                   'measurementByIndividualAtom.measurement.error',
            'Resonance_ID':              'measurementByIndividualAtom.resonance',
            #'Details':                  'measurementByIndividualAtom.measurement.details',

            'CUSTOM_setAtoms_individual': 'measurementByIndividualAtom'

            },
          }, 

        }
      }


    # TODO - implement this table - commented out ccpnMap as don't know how to get from CCPN.

    self.sfDict['homonucl_NOEs'] = {                                                                          # *** Priority 4 - data

      #'ccpnLoop':                       "nmrEntry.findAllMeasurementLists(className='NoeList')",
      #'ccpnMap':                        'noeList',

      #'title':                         ('noeList.name',self.getTitle),

      'tags': {

        #'Sf_category':                  ['homonucl_NOEs',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory
        #'Sample_condition_list_ID':     [None,returnStarInt,'Sample_condition_list.ID',True],
        #'Sample_condition_list_label':  [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Mandatory - option
        #'Homonuclear_NOE_val_type':     [None,lambda x = value: returnStarCode(x,length = 31),None,True],    # *** Mandatory - option
        #'NOE_ref_val':                  [None,returnStarFloat,None,False],                                   # *** Mandatory
        #'NOE_ref_description':          [None,returnStarString,None,False],                                  # *** Optional
        #'Details':                      [None,returnStarString,None,False],                                  # *** Optional
        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        #'Homonucl_NOE_experiment':      self.setExperimentsLoop('noeList.sortedExperiments()','experiment'), # *** Table

        'Homonucl_NOE_experiment': {                                                                          # *** Includes mandatory fields

          'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Homonucl_NOE_list_ID':     [None,returnStarInt,'Homonucl_NOE_list.ID',True],

            },
          },


        #'Homonucl_NOE_software':        self.setSoftwareLoop('noeList'),                                     # *** Table

        'Homonucl_NOE_software': {                                                                            # *** Includes optional fields

          'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Homonucl_NOE_list_ID':     [None,returnStarInt,'Homonucl_NOE_list.ID',True],

            },
          },

        'Homonucl_NOE': {                                                                                     # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Assembly_atom_ID_1':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_1':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_1':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_1':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_1':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_1':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_1':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_1':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number_1':    [None,returnStarInt,None,False],
            #'Assembly_atom_ID_2':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_2':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_2':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_2':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_2':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_2':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_2':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_2':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number_2':    [None,returnStarInt,None,False],

            #'Val':                      [None,returnStarFloat,None,True],
            #'Val_err':                  [None,returnStarFloat,None,False],

            #'Resonance_ID_1':           [None,returnStarInt,'Resonance.ID',False],
            #'Resonance_ID_2':           [None,returnStarInt,'Resonance.ID',False],
            #'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID_1':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID_2':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Homonucl_NOE_list_ID':     [None,returnStarInt,'Homonucl_NOE_list.ID',True],

            #'CUSTOM_setAtoms':          'measurement.sortedResonances()',

            },
          },

        }
      }


    self.sfDict['heteronucl_NOEs'] = {                                                                        # *** Priority 1 - data

      'ccpnLoop':                        "nmrEntry.findAllMeasurementLists(className='NoeList')",
      'ccpnMap':                         'noeList',

      'title':                          ('noeList',self.getMeasurementTitle),

      'tags': {

        #'Sf_category':                  ['heteronucl_NOEs',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Data_file_name':               ('noeList.root.name',self.getProjectFileName),                        # *** Mandatory
        'Sample_condition_list_ID':      'noeList.findFirstExperiment().sampleConditionSet',
        'Sample_condition_list_label':   'noeList.findFirstExperiment().sampleConditionSet.name',             # *** Mandatory - option
        'Spectrometer_frequency_1H':     'noeList.sf',                                                        # *** Mandatory
        'Heteronuclear_NOE_val_type':   ('noeList.noeValueType',self.getStarNoeValueType),                    # *** Mandatory - option
        'NOE_ref_val':                   'noeList.refValue',                                                  # *** Mandatory
        'NOE_ref_description':           'noeList.refDescription',                                            # *** Optional

        'Details':                       'noeList.details',                                                   # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'Heteronucl_NOE_experiment':     self.setExperimentsLoop('noeList.sortedExperiments()','experiment'), # *** Table

        #'Heteronucl_NOE_experiment': {                                                                       # *** Includes mandatory fields
        
        #  'ccpnLoop':                   'noeList.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Heteronucl_NOE_list_ID':   [None,returnStarInt,'Heteronucl_NOE_list.ID',True],
            
        #    },
        #  },

        'Heteronucl_NOE_software':       self.setSoftwareLoop('noeList'),                                     # *** Table

        #'Heteronucl_NOE_software': {                                                                         # *** Includes optional fields

        #  'ccpnMap':                    'noeList',

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Heteronucl_NOE_list_ID':   [None,returnStarInt,'Heteronucl_NOE_list.ID',True],

        #    },
        #  },

        'Heteronucl_NOE': {                                                                                   # *** Table

          'ccpnLoop':                    'noeList.sortedMeasurements()',
          'ccpnMap':                     'measurement',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Assembly_atom_ID_1':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_1':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_1':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_1':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_1':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_1':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_1':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_1':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number_1':    [None,returnStarInt,None,False],
            #'Assembly_atom_ID_2':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_2':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_2':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_2':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_2':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_2':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_2':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_2':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number_2':    [None,returnStarInt,None,False],

            'Val':                       'measurement.value',
            'Val_err':                   'measurement.error',

            #'Resonance_ID_1':           [None,returnStarInt,'Resonance.ID',False],
            #'Resonance_ID_2':           [None,returnStarInt,'Resonance.ID',False],
            #'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID_1':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID_2':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Heteronucl_NOE_list_ID':   [None,returnStarInt,'Heteronucl_NOE_list.ID',True],

            'CUSTOM_setAtoms':           'measurement.sortedResonances()',

            },
          },

        }
      }


    self.sfDict['heteronucl_T1_relaxation'] = {                                                               # *** Priority 1 - data
    
      'ccpnLoop':                        "nmrEntry.findAllMeasurementLists(className='T1List')",
      'ccpnMap':                         't1List',
      
      'title':                          ('t1List',self.getMeasurementTitle),

      'tags': {

        #'Sf_category':                  ['heteronucl_T1_relaxation',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Data_file_name':               ('t1List.root.name',self.getProjectFileName),                         # *** Mandatory
        'Sample_condition_list_ID':      't1List.findFirstExperiment().sampleConditionSet',
        'Sample_condition_list_label':   't1List.findFirstExperiment().sampleConditionSet.name',              # *** Mandatory - option
        'Spectrometer_frequency_1H':     't1List.sf',                                                         # *** Mandatory
        'T1_coherence_type':            ('t1List',self.getStarT1CoherenceType),                               # *** Mandatory - option
        'T1_val_units':                  't1List.unit',                                                       # *** Mandatory - option
        'Details':                       't1List.details',                                                    # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'Heteronucl_T1_experiment':      self.setExperimentsLoop('t1List.sortedExperiments()','experiment'),  # *** Table

        #'Heteronucl_T1_experiment': {                                                                        # *** Includes mandatory fields
        
        #  'ccpnLoop':                   't1List.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Heteronucl_T1_list_ID':    [None,returnStarInt,'Heteronucl_T1_list.ID',True],

        #    },
        #  },

        'Heteronucl_T1_software':        self.setSoftwareLoop('t1List'),                                      # *** Table

        #'Heteronucl_T1_software': {                                                                          # *** Includes optional fields

        #  'ccpnMap':                    't1List',

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Heteronucl_T1_list_ID':    [None,returnStarInt,'Heteronucl_T1_list.ID',True],

        #    },
        #  },

        'T1': {                                                                                               # *** Table
        
          'ccpnLoop':                   ('t1List.sortedMeasurements()',self.getMeasurementsByIndividualAtom),
          'ccpnMap':                     'measurementByIndividualAtom',

          'tags': {
                
            #'ID':                       'shiftByIndividualAtom.serial', # TODO this is a hack - can this be done better?
            'Atom_isotope_number':       'measurementByIndividualAtom.resonance.isotope.massNumber',
            'Val':                       'measurementByIndividualAtom.measurement.value',
            'Val_err':                   'measurementByIndividualAtom.measurement.error',
            'Resonance_ID':              'measurementByIndividualAtom.resonance',
            #'Details':                  'measurementByIndividualAtom.measurement.details',

            'CUSTOM_setAtoms_individual': 'measurementByIndividualAtom'

            },
          }, 

        }
      }


    self.sfDict['heteronucl_T1rho_relaxation'] = {                                                            # *** Priority 1 - data
    
      'ccpnLoop':                        "nmrEntry.findAllMeasurementLists(className='T1RhoList')",
      'ccpnMap':                         't1RhoList',
      
      'title':                          ('t1RhoList',self.getMeasurementTitle),

      'tags': {

        #'Sf_category':                  ['heteronucl_T1rho_relaxation',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Data_file_name':               ('t1List.root.name',self.getProjectFileName),                         # *** Mandatory
        'Sample_condition_list_ID':      't1RhoList.findFirstExperiment().sampleConditionSet',
        'Sample_condition_list_label':   't1RhoList.findFirstExperiment().sampleConditionSet.name',           # *** Mandatory - option
        'Spectrometer_frequency_1H':     't1RhoList.sf',                                                      # *** Mandatory
        'T1rho_coherence_type':          't1RhoList.coherenceType',                                           # *** Mandatory - option
        'T1rho_val_units':               't1RhoList.unit',                                                    # *** Mandatory - option
        'Temp_calibration_method':       't1RhoList.tempCalibMethod',                                         # *** Mandatory - option
        'Temp_control_method':           't1RhoList.tempControlMethod',                                       # *** Mandatory - option

        #'Rex_units':                    [None,lambda x = value: returnStarCode(x,length = 31),None,False],   # *** Mandatory - option

        'Details':                       't1RhoList.details',                                                 # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'Heteronucl_T1rho_experiment':   self.setExperimentsLoop('t1RhoList.sortedExperiments()','experiment'), # *** Table

        #'Heteronucl_T1rho_experiment': {                                                                     # *** Includes mandatory fields
        
        #  'ccpnLoop':                   'T1RhoList.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Heteronucl_T1rho_list_ID': [None,returnStarInt,'Heteronucl_T1rho_list.ID',True],

        #    },
        #  },

        'Heteronucl_T1rho_software':     self.setSoftwareLoop('t1RhoList'),                                   # *** Table

        #'Heteronucl_T1rho_software': {                                                                       # *** Includes optional fields

        #  'tags': {

        #  'ccpnMap':                    't1RhoList',

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Heteronucl_T1rho_list_ID': [None,returnStarInt,'Heteronucl_T1rho_list.ID',True],

        #    },
        #  },


        'T1rho': {                                                                                            # *** Table
        
          'ccpnLoop':                   ('t1RhoList.sortedMeasurements()',self.getMeasurementsByIndividualAtom),
          'ccpnMap':                     'measurementByIndividualAtom',

          'tags': {
                
            #'ID':                       'shiftByIndividualAtom.serial', # TODO this is a hack - can this be done better?
            'Atom_isotope_number':       'measurementByIndividualAtom.resonance.isotope.massNumber',
            'T1rho_val':                 'measurementByIndividualAtom.measurement.value',
            'T1rho_val_err':             'measurementByIndividualAtom.measurement.error',
            #'Rex_val':                  [None,returnStarFloat,None,False],
            #'Rex_val_err':              [None,returnStarFloat,None,False],
            'Resonance_ID':              'measurementByIndividualAtom.resonance',
            #'Details':                  'measurementByIndividualAtom.measurement.details',

            'CUSTOM_setAtoms_individual': 'measurementByIndividualAtom'

            },
          }, 

        }
      }


    self.sfDict['heteronucl_T2_relaxation'] = {                                                               # *** Priority 1 - data
    
      'ccpnLoop':                        "nmrEntry.findAllMeasurementLists(className='T2List')",
      'ccpnMap':                         't2List',
      
      'title':                          ('t2List',self.getMeasurementTitle),

      'tags': {
      
        #'Sf_category':                  ['heteronucl_T2_relaxation',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Data_file_name':               ('t2List.root.name',self.getProjectFileName),                         # *** Mandatory
        'Sample_condition_list_ID':      't2List.findFirstExperiment().sampleConditionSet',
        'Sample_condition_list_label':   't2List.findFirstExperiment().sampleConditionSet.name',              # *** Mandatory - option
        'Spectrometer_frequency_1H':     't2List.sf',                                                         # *** Mandatory
        'T2_coherence_type':             't2List.coherenceType',                                              # *** Mandatory - option
        'T2_val_units':                  't2List.unit',                                                       # *** Mandatory - option
        'Temp_calibration_method':       't2List.tempCalibMethod',                                            # *** Mandatory - option
        'Temp_control_method':           't2List.tempControlMethod',                                          # *** Mandatory - option

        #'Rex_units':                    [None,lambda x = value: returnStarCode(x,length = 31),None,False],   # *** Optional - option

        'Details':                       't2List.details',                                                    # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'Heteronucl_T2_experiment':      self.setExperimentsLoop('t2List.sortedExperiments()','experiment'),  # *** Table

        #'Heteronucl_T2_experiment': {                                                                        # *** Includes mandatory fields
        
        #  'ccpnLoop':                   't2List.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Heteronucl_T2_list_ID':    [None,returnStarInt,'Heteronucl_T2_list.ID',True],

        #    },
        #  }, 

        'Heteronucl_T2_software':        self.setSoftwareLoop('t2List'),                                      # *** Table

        #'Heteronucl_T2_software': {                                                                          # *** Includes optional fields

        #  'ccpnMap':                    't2List',

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Heteronucl_T2_list_ID':    [None,returnStarInt,'Heteronucl_T2_list.ID',True],

        #    },
        #  },

        'T2': {                                                                                               # *** Table
        
          'ccpnLoop':                   ('t2List.sortedMeasurements()',self.getMeasurementsByIndividualAtom),
          'ccpnMap':                     'measurementByIndividualAtom',

          'tags': {
                
            #'ID':                       'shiftByIndividualAtom.serial', # TODO this is a hack - can this be done better?
            'Atom_isotope_number':       'measurementByIndividualAtom.resonance.isotope.massNumber',
            'T2_val':                    'measurementByIndividualAtom.measurement.value',
            'T2_val_err':                'measurementByIndividualAtom.measurement.error',
            #'Rex_val':                  [None,returnStarFloat,None,False],
            #'Rex_err':                  [None,returnStarFloat,None,False],
            'Resonance_ID':              'measurementByIndividualAtom.resonance',
            #'Details':                  'measurementByIndividualAtom.measurement.details',

            'CUSTOM_setAtoms_individual': 'measurementByIndividualAtom'

            },
          }, 

        }
      }


    self.sfDict['dipole_dipole_relaxation'] = {                                                               # *** Priority 2 - data

      'ccpnLoop':                        "nmrEntry.findAllMeasurementLists(className='DipolarRelaxList')",
      'ccpnMap':                         'dipolarRelaxList',
      
      'title':                          ('dipolarRelaxList.name',self.getTitle),

      'tags': {

        #'Sf_category':                  ['dipole_dipole_relaxation',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Data_file_name':               ('dipolarRelaxList.root.name',self.getProjectFileName),               # *** Mandatory
        'Sample_condition_list_ID':      'dipolarRelaxList.findFirstExperiment().sampleConditionSet',
        'Sample_condition_list_label':   'dipolarRelaxList.findFirstExperiment().sampleConditionSet.name',    # *** Mandatory - option
        'Spectrometer_frequency_1H':     'dipolarRelaxList.sf',                                               # *** Mandatory
        'Val_units':                     'dipolarRelaxList.unit',                                             # *** Mandatory - option
        'Details':                       'dipolarRelaxList.details',                                          # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'Dipole_dipole_experiment':      self.setExperimentsLoop('dipolarRelaxList.sortedExperiments()','experiment'), # *** Table

        #'Dipole_dipole_relax_experiment': {                                                                  # *** Includes mandatory fields
        
        #  'ccpnLoop':                   'dipolarRelaxList.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Dipole_dipole_relax_list_ID': [None,returnStarInt,'Dipole_dipole_relax_list.ID',True],

        #    },
        #  }, 

        'Dipole_dipole_relax_software':  self.setSoftwareLoop('dipolarRelaxList'),                            # *** Table

        #'Dipole_dipole_relax_software': {                                                                    # *** Includes optional fields

        #  'tags': {

        #  'ccpnMap':                    'dipolarRelaxList',

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Dipole_dipole_relax_list_ID': [None,returnStarInt,'Dipole_dipole_relax_list.ID',True],

        #    },
        #  },

        'Dipole_dipole_relax': {                                                                              # *** Table

          'ccpnLoop':                    'dipolarRelaxList.sortedMeasurements()',
          'ccpnMap':                     'measurement',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Assembly_atom_ID_1':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_1':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_1':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_1':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_1':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_1':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_1':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_1':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number_1':    [None,returnStarInt,None,False],
            #'Assembly_atom_ID_2':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_2':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_2':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_2':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_2':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_2':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_2':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_2':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number_2':    [None,returnStarInt,None,False],

            'Val':                       'measurement.value',
            'Val_err':                   'measurement.error',

            #'Resonance_ID_1':           [None,returnStarInt,'Resonance.ID',False],
            #'Resonance_ID_2':           [None,returnStarInt,'Resonance.ID',False],
            #'Auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID_1':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID_2':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Dipole_dipole_relax_list_ID': [None,returnStarInt,'Dipole_dipole_relax_list.ID',True],

            'CUSTOM_setAtoms':           'measurement.sortedResonances()',

            },
          },

        }
      }


    self.sfDict['dipole_dipole_cross_correlations'] = {                                                       # *** Priority 3 - data

      #'ccpnLoop':                       None,
      #'ccpnMap':                        None,

      #'title':                          None,

      'tags': {

        #'Sf_category':                  ['dipole_dipole_cross_correlations',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory
        #'Sample_condition_list_ID':     [None,returnStarInt,'Sample_condition_list.ID',True],
        #'Sample_condition_list_label':  [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Mandatory - option
        #'Spectrometer_frequency_1H':    [None,returnStarFloat,None,True],                                    # *** Mandatory
        #'Val_units':                    [None,lambda x = value: returnStarCode(x,length = 31),None,True],    # *** Mandatory - option
        #'Details':                      [None,returnStarString,None,False],                                  # *** Optional
        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

      },

      'tables': {

        #'Cross_correlation_DD_experiment': self.setExperimentsLoop('crossCorDDList.sortedExperiments()','experiment'), # *** Table

        'Cross_correlation_DD_experiment': {                                                                  # *** Includes mandatory fields

          'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Cross_correlation_list_ID': [None,returnStarInt,'Cross_correlation_DD_list.ID',True],

            },
          },

        #'Cross_correlation_DD_software': self.setSoftwareLoop('crossCorDDList'),                             # *** Table

        'Cross_correlation_DD_software': {                                                                    # *** Includes optional fields

          'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Cross_correlation_list_ID': [None,returnStarInt,'Cross_correlation_DD_list.ID',True],

            },
          },

        'Cross_correlation_DD': {                                                                             # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Dipole_1_entry_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Dipole_1_entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',False],
            #'Dipole_1_entity_ID_1':     [None,returnStarInt,'Entity.ID',False],
            #'Dipole_1_comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',False],
            #'Dipole_1_seq_ID_1':        [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
            #'Dipole_1_comp_ID_1':       [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
            #'Dipole_1_atom_ID_1':       [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
            #'Dipole_1_atom_type_1':     [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_1_atom_isotope_number_1': [None,returnStarInt,None,False],
            #'Dipole_1_entry_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Dipole_1_entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',False],
            #'Dipole_1_entity_ID_2':     [None,returnStarInt,'Entity.ID',False],
            #'Dipole_1_comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',False],
            #'Dipole_1_seq_ID_2':        [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
            #'Dipole_1_comp_ID_2':       [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
            #'Dipole_1_atom_ID_2':       [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
            #'Dipole_1_atom_type_2':     [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_1_atom_isotope_number_2': [None,returnStarInt,None,False],
            #'Dipole_2_entry_atom_ID_1': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Dipole_2_entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',False],
            #'Dipole_2_entity_ID_1':     [None,returnStarInt,'Entity.ID',False],
            #'Dipole_2_comp_index_ID_1': [None,returnStarInt,'Entity_comp_index.ID',False],
            #'Dipole_2_seq_ID_1':        [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
            #'Dipole_2_comp_ID_1':       [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
            #'Dipole_2_atom_ID_1':       [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
            #'Dipole_2_atom_type_1':     [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_2_atom_isotope_number_1': [None,returnStarInt,None,False],
            #'Dipole_2_entry_atom_ID_2': [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Dipole_2_entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',False],
            #'Dipole_2_entity_ID_2':     [None,returnStarInt,'Entity.ID',False],
            #'Dipole_2_chem_comp_index_ID_2': [None,returnStarInt,'Entity_comp_index.ID',False],
            #'Dipole_2_seq_ID_2':        [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
            #'Dipole_2_comp_ID_2':       [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
            #'Dipole_2_atom_ID_2':       [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
            #'Dipole_2_atom_type_2':     [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_2_atom_isotope_number_2': [None,returnStarInt,None,False],

            #'Val':                      [None,returnStarFloat,None,True],
            #'Val_err':                  [None,returnStarFloat,None,False],

            #'Resonance_ID':             [None,returnStarInt,'Resonance.ID',False],
            #'Dipole_1_auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_1_auth_seq_ID_1':   [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_1_auth_comp_ID_1':  [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_1_auth_atom_ID_1':  [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_1_auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_1_auth_seq_ID_2':   [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_1_auth_comp_ID_2':  [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_1_auth_atom_ID_2':  [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_2_auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_2_auth_seq_ID_1':   [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_2_auth_comp_ID_1':  [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_2_auth_atom_ID_1':  [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_2_auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_2_auth_seq_ID_2':   [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_2_auth_comp_ID_2':  [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_2_auth_atom_ID_2':  [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Cross_correlation_list_ID': [None,returnStarInt,'Cross_correlation_DD_list.ID',True],

            },
          },

        }
      }


    self.sfDict['dipole_CSA_cross_correlations'] = {                                                          # *** Priority 3 - data

      #'ccpnLoop':                       None,
      #'ccpnMap':                        None,

      #'title':                          None,

      'tags': {

        #'Sf_category':                  ['dipole_CSA_cross_correlations',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory
        #'Sample_condition_list_ID':     [None,returnStarInt,'Sample_condition_list.ID',True],
        #'Sample_condition_list_label':  [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Mandatory - option
        #'Spectrometer_frequency_1H':    [None,returnStarFloat,None,True],                                    # *** Mandatory
        #'Val_units':                    [None,lambda x = value: returnStarCode(x,length = 31),None,True],    # *** Mandatory - option
        #'Details':                      [None,returnStarString,None,False],                                  # *** Optional
        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

      },

      'tables': {

        #'Cross_correlation_D_CSA_experiment': self.setExperimentsLoop('crossCorDCSAList.sortedExperiments()','experiment'), # *** Table

        'Cross_correlation_D_CSA_experiment': {                                                               # *** Includes mandatory fields

          'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Cross_correlation_list_ID': [None,returnStarInt,'Cross_correlation_D_CSA_list.ID',True],

            },
          },

        #'Cross_correlation_D_CSA_software': self.setSoftwareLoop('crossCorDCSAList'),                        # *** Table

        'Cross_correlation_D_CSA_software': {                                                                 # *** Includes optional fields

          'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Cross_correlation_list_ID': [None,returnStarInt,'Cross_correlation_D_CSA_list.ID',True],

            },
          },

        'Cross_correlation_D_CSA': {                                                                          # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Dipole_entry_atom_ID_1':   [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Dipole_entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',False],
            #'Dipole_entity_ID_1':       [None,returnStarInt,'Entity.ID',False],
            #'Dipole_comp_index_ID_1':   [None,returnStarInt,'Entity_comp_index.ID',False],
            #'Dipole_seq_ID_1':          [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
            #'Dipole_comp_ID_1':         [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
            #'Dipole_atom_ID_1':         [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
            #'Dipole_atom_type_1':       [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_atom_isotope_number_1': [None,returnStarInt,None,False],
            #'Dipole_entry_atom_ID_2':   [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Dipole_entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',False],
            #'Dipole_entity_ID_2':       [None,returnStarInt,'Entity.ID',False],
            #'Dipole_comp_index_ID_2':   [None,returnStarInt,'Entity_comp_index.ID',False],
            #'Dipole_seq_ID_2':          [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
            #'Dipole_comp_ID_2':         [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
            #'Dipole_atom_ID_2':         [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
            #'Dipole_atom_type_2':       [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_atom_isotope_number_2': [None,returnStarInt,None,False],
            #'CSA_entry_atom_ID_1':      [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'CSA_entity_assembly_ID_1': [None,returnStarInt,'Entity_assembly.ID',False],
            #'CSA_entity_ID_1':          [None,returnStarInt,'Entity.ID',False],
            #'CSA_comp_index_ID_1':      [None,returnStarInt,'Entity_comp_index.ID',False],
            #'CSA_seq_ID_1':             [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
            #'CSA_comp_ID_1':            [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
            #'CSA_atom_ID_1':            [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
            #'CSA_atom_type_1':          [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'CSA_atom_isotope_number_1': [None,returnStarInt,None,False],
            #'CSA_entry_atom_ID_2':      [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'CSA_entity_assembly_ID_2': [None,returnStarInt,'Entity_assembly.ID',False],
            #'CSA_entity_ID_2':          [None,returnStarInt,'Entity.ID',False],
            #'CSA_comp_index_ID_2':      [None,returnStarInt,'Entity_comp_index.ID',False],
            #'CSA_seq_ID_2':             [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',False],
            #'CSA_comp_ID_2':            [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',False],
            #'CSA_atom_ID_2':            [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',False],
            #'CSA_atom_type_2':          [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'CSA_atom_isotope_number_2': [None,returnStarInt,None,False],

            #'Val':                      [None,returnStarFloat,None,True],
            #'Val_err':                  [None,returnStarFloat,None,False],

            #'Resonance_ID_1':           [None,returnStarInt,'Resonance.ID',False],
            #'Resonance_ID_2':           [None,returnStarInt,'Resonance.ID',False],
            #'Dipole_auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_auth_seq_ID_1':     [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_auth_comp_ID_1':    [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_auth_atom_ID_1':    [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_auth_seq_ID_2':     [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_auth_comp_ID_2':    [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Dipole_auth_atom_ID_2':    [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'CSA_auth_entity_assembly_ID_1': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'CSA_auth_seq_ID_1':        [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'CSA_auth_comp_ID_1':       [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'CSA_auth_atom_ID_1':       [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'CSA_auth_entity_assembly_ID_2': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'CSA_auth_seq_ID_2':        [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'CSA_auth_comp_ID_2':       [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'CSA_auth_atom_ID_2':       [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Cross_correlation_list_ID': [None,returnStarInt,'Cross_correlation_D_CSA_list.ID',True],

            },
          },

        }
      }


    self.sfDict['order_parameters'] = {                                                                       # *** Priority 1 - data
    
      #'ccpnLoop':                       "nmrEntry.findAllDerivedDataLists(className='IsotropicS2List')",
      'ccpnLoop':                       ('nmrEntry', self.getAllIsotropicS2Derivations),
      'ccpnMap':                         'isotropicS2Derivation',

      'title':                          ('isotropicS2Derivation.details',self.getTitle),

      'tags': {

        #'Sf_category':                  ['order_parameters',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory
        #'Sample_condition_list_ID':     [None,returnStarInt,'Sample_condition_list.ID',True],
        #'Sample_condition_list_label':  [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Mandatory - option

        'Tau_e_val_units':               'isotropicS2Derivation.parent.tauEUnit',                             # *** Optional - option
        'Tau_s_val_units':               'isotropicS2Derivation.parent.tauSUnit',                             # *** Optional - option
        'Details':                       'isotropicS2Derivation.details',                                     # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'Order_parameter_experiment':    self.setExperimentsLoop('isotropicS2Derivation.parent.sortedExperiments()','experiment'), # *** Table

        #'Order_parameter_experiment': {                                                                      # *** Includes mandatory fields
        
        #  'ccpnLoop':                   'isotropicS2Derivation.parent.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Order_parameter_list_ID':  [None,returnStarInt,'Order_parameter_list.ID',True],

        #    },
        #  }, 

        'Order_parameter_software':      self.setSoftwareLoop('isotropicS2Derivation.parent'),                # *** Table

        #'Order_parameter_software': {                                                                         # *** Includes optional fields

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Order_parameter_list_ID':  [None,returnStarInt,'Order_parameter_list.ID',True],

        #    },
        #  },

        'Order_param': {                                                                                      # *** Table

          'ccpnLoop':                    'isotropicS2Derivation.sortedDerivedData()',
          'ccpnMap':                     'datum',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Assembly_atom_ID':         [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID':       [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID':                [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID':            [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID':                   [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID':                  [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type':                [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number':      [None,returnStarInt,None,False],

            'Order_param_val':           'datum.value',
            'Order_param_val_fit_err':   'datum.error',

            'Tau_e_val':                 'datum.tauEValue', # TODO - what are the values likely to be?
            'Tau_e_val_fit_err':         'datum.tauEError',

            'Rex_val':                   'datum.rexValue',
            'Rex_val_fit_err':           'datum.rexError',
            'Model_free_sum_squared_errs': 'datum.sumSquaredErrors',
            'Model_fit':                 'datum.modelFit',

            #'Sf2_val':                  [None,returnStarFloat,None,False],
            #'Sf2_val_fit_err':          [None,returnStarFloat,None,False],
            #'Ss2_val':                  [None,returnStarFloat,None,False],
            #'Ss2_val_fit_err':          [None,returnStarFloat,None,False],

            'Tau_s_val':                 'datum.tauSValue',
            'Tau_s_val_fit_err':         'datum.tauSError',

            #'SH2_val':                  [None,returnStarFloat,None,False],
            #'SH2_val_fit_err':          [None,returnStarFloat,None,False],
            #'SN2_val':                  [None,returnStarFloat,None,False],
            #'SN2_val_fit_err':          [None,returnStarFloat,None,False],

            #'Resonance_ID':             [None,returnStarInt,'Resonance.ID',False],
            #'Auth_entity_assembly_ID':  [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID':              [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID':             [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID':             [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Order_parameter_list_ID':  [None,returnStarInt,'Order_parameter_list.ID',True],

            'CUSTOM_setAtoms':           'datum.sortedResonances()',

            },
          }, 

        }
      }


    self.sfDict['pH_titration'] = {                                                                           # *** Priority 4 - data

      #'ccpnLoop':                       "nmrEntry.findAllDerivedDataLists(className='PkaList')",
      'ccpnLoop':                       ('nmrEntry', self.getAllPkaDerivations),
      'ccpnMap':                         'pkaDerivation',

      'title':                          ('pkaDerivation.details',self.getTitle),

      'tags': {

        #'Sf_category':                  ['pH_titration',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory
        #'Sample_condition_list_ID':     [None,returnStarInt,'Sample_condition_list.ID',True],
        #'Sample_condition_list_label':  [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Mandatory - option

        'Expt_observed_param':           'pkaDerivation.findFirstDerivedData().parameterType', # TODO - is this right? # *** Mandatory - option
        'Details':                       'pkaDerivation.details',                                             # *** Optional

        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        'PH_titration_experiment':       self.setExperimentsLoop('pkaDerivation.parent.sortedExperiments()','experiment'), # *** Table

        #'PH_titration_experiment': {                                                                         # *** Includes mandatory fields

        #  'ccpnLoop':                   'pkaDerivation.parent.sortedExperiments()',
        #  'ccpnMap':                    'experiment',

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'PH_titration_list_ID':     [None,returnStarInt,'PH_titration_list.ID',True],

        #    },
        #  },

        'PH_titration_software':         self.setSoftwareLoop('pkaDerivation.parent'),                        # *** Table

        #'PH_titration_software': {                                                                            # *** Includes optional fields

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'PH_titration_list_ID':     [None,returnStarInt,'PH_titration_list.ID',True],

        #    },
        #  },

        'PH_titr_result': {                                                                                   # *** Table

          'ccpnLoop':                    'pkaDerivation.sortedDerivedData()',
          'ccpnMap':                     'datum',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Atm_obs_entry_atom_ID':    [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Atm_obs_entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
            #'Atm_obs_entity_ID':        [None,returnStarInt,'Entity.ID',True],
            #'Atm_obs_comp_index_ID':    [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Atm_obs_seq_ID':           [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Atm_obs_comp_ID':          [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atm_obs_atom_ID':          [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atm_obs_atom_type':        [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atm_obs_atom_isotope_number': [None,returnStarInt,None,False],
            #'Atm_obs_auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Atm_obs_auth_seq_ID':      [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Atm_obs_auth_comp_ID':     [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Atm_obs_auth_atom_ID':     [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Atm_titr_entry_atom_ID':   [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Atm_titr_entity_assembly_ID': [None,returnStarInt,'Entity_assembly.ID',True],
            #'Atm_titr_entity_ID':       [None,returnStarInt,'Entity.ID',True],
            #'Atm_titr_comp_index_ID':   [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Atm_titr_seq_ID':          [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Atm_titr_comp_ID':         [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atm_titr_atom_ID':         [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atm_titr_atom_type':       [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atm_titr_atom_isotope_number': [None,returnStarInt,None,False],
            #'Atm_titr_auth_entity_assembly_ID': [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Atm_titr_auth_seq_ID':     [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Atm_titr_auth_comp_ID':    [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Atm_titr_auth_atom_ID':    [None,lambda x = value: returnStarCode(x,length = 15),None,False],

            'Hill_coeff_val':            'datum.hillCoeff',
            'Hill_coeff_val_fit_err':    'datum.hillCoeffError',
            'High_PH_param_fit_val':     'datum.highPHParam',
            'High_PH_param_fit_val_err': 'datum.highPHParamError',
            'Low_PH_param_fit_val':      'datum.lowPHParam',
            'Low_PH_param_fit_val_err':  'datum.lowPHParamError',

            'PKa_val':                   'datum.value',
            'PKa_val_fit_err':           'datum.error',

            #'PHmid_val':                [None,returnStarFloat,None,False],
            #'PHmid_val_fit_err':        [None,returnStarFloat,None,False],

            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'PH_titration_list_ID':     [None,returnStarInt,'PH_titration_list.ID',True],

            'CUSTOM_setAtoms':           'datum.sortedResonances()',

            },
          },

        }
      }

    self.sfDict['pH_param_list'] = {                                                                          # *** Priority 4 - data

      #'ccpnLoop':                       None,
      #'ccpnMap':                        None,

      #'title':                          None,

      'tags': {

        #'Sf_category':                  ['pH_param_list',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'PH_titration_list_ID':         [None,returnStarInt,'PH_titration_list.ID',True],
        #'PH_titration_list_label':      [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Mandatory - option
        #'Observed_NMR_param':           [None,lambda x = value: returnStarCode(x,length = 31),None,True],    # *** Optional - option
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory
        #'Details':                      [None,returnStarString,None,False],                                  # *** Optional
        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

      },

      'tables': {

        'PH_param': {                                                                                         # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'PH_titr_result_ID':        [None,returnStarInt,'PH_titr_result.ID',True],
            #'PH_val':                   [None,returnStarFloat,None,True],
            #'PH_val_err':               [None,returnStarFloat,None,False],
            #'Observed_NMR_param_val':   [None,returnStarFloat,None,True],
            #'Observed_NMR_param_val_err': [None,returnStarFloat,None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'PH_param_list_ID':         [None,returnStarInt,'PH_param_list.ID',True],

            },
          },

        }
      }


    self.sfDict['D_H_fractionation_factors'] = {                                                              # *** Priority 4 - data

      #'ccpnLoop':                       None,
      #'ccpnMap':                        None,

      #'title':                          None,

      'tags': {

        #'Sf_category':                  ['D_H_fractionation_factors',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],  # *** Mandatory - set in NmrStarExport - auto
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],  # *** Mandatory
        #'Sample_condition_list_ID':     [None,returnStarInt,'Sample_condition_list.ID',True],
        #'Sample_condition_list_label':  [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Mandatory - option
        #'Details':                      [None,returnStarString,None,False],                                  # *** Optional
        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],

        },

      'tables': {

        #'D_H_fract_factor_experiment':  self.setExperimentsLoop('dHFractList.sortedExperiments()','experiment'), # *** Table

        'D_H_fract_factor_experiment': {                                                                      # *** Includes mandatory fields

          'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False], # *** Mandatory - pulldown
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'D_H_fractionation_factor_list_ID': [None,returnStarInt,'D_H_fractionation_factor_list.ID',True],

            },
          },

        #'D_H_fract_factor_software':    self.setSoftwareLoop('dHFractList'),                                 # *** Table

        'D_H_fract_factor_software': {                                                                        # *** Includes optional fields

          'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],  # *** Optional
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'D_H_fractionation_factor_list_ID': [None,returnStarInt,'D_H_fractionation_factor_list.ID',True],

            },
          },

        'D_H_fractionation_factor': {                                                                         # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Assembly_atom_ID':         [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID':       [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID':                [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID':            [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID':                   [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID':                  [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type':                [None,lambda x = value: returnStarCode(x,length = 15),None,True],

            #'Val':                      [None,returnStarFloat,None,True],
            #'Val_err':                  [None,returnStarFloat,None,False],

            #'Resonance_ID':             [None,returnStarInt,'Resonance.ID',False],
            #'Auth_entity_assembly_ID':  [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_seq_ID':              [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID':             [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID':             [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'D_H_fractionation_factor_list_ID': [None,returnStarInt,'D_H_fractionation_factor_list.ID',True],

            },
          },

        }
      }


    ###############
    # Coordinates #
    ###############

    self.sfDict['conformer_family_coord_set'] = {                                                             # *** Priority 1 - data

      'ccpnLoop':                        'nmrEntry.sortedStructureGenerations()',
      'ccpnMap':                         'structureGeneration',
      
      'addCcpnMap':                     {'structureGeneration': [('structureEnsemble','structureGeneration.structureEnsemble')]},

      'title':                          ('structureGeneration.name',self.getTitle),

      'tags': {

        'File_name':                     None,
        'Constraints_PDB_file_ID':       None,
        'PDB_accession_code':            None,
        'Sample_condition_list_ID':      None,
        'Sample_condition_list_label':   None,
        #'Details':                      'structureGeneration.details',

        # CJP - missing tags?

        #'Sf_framecode':                 [None,returnStarString,None,False],
        #'Sample_condition_list_ID':     [None,returnStarInt,'Sample_condition_list.ID',True],
        #'Sample_condition_list_label':  [None,returnStarString,None,True],

        },

      'tables': {

        'Conformer_family_refinement': {                                                                      # *** Includes mandatory (PDB) fields

          'ccpnMap':                     'structureGeneration',

          'tags': {

            'Refine_method':             'structureGeneration.method.name',                                   # *** Mandatory (PDB) - option
            'Refine_details':            'structureGeneration.method.software.details',                       # *** Optional
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

            },
          }, 

        'Conformer_family_software':     self.setSoftwareLoop('structureGeneration'),                         # *** Table

        #'Conformer_family_software': {                                                                        # *** Table

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

        #    },
        #  },

        'Energetic_penalty_function': {                                                                       # *** Table

          'tags': {

            #'Function':                 [None,lambda x = value: returnStarLine(x,length = 127),None,True],
            #'Description':              [None,lambda x = value: returnStarString(x,length = 255),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

            },
          },

        'Conformer_family_coord_set_expt': {                                                                  # *** Table

          'tags': {

        #    CJP - missing tags marked with a +

        #    'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
        #    'Experiment_name':          [None,returnStarString,None,False],
        #  + 'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
        #  + 'Sample_label':             [None,returnStarString,None,True],
        #    'Sample_state':             [None,returnStarString,None,True],
        #  + 'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
        #  + 'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

            },
          },

        'Conf_family_coord_set_constr_list': {                                                                # *** Table

          'tags': {

            #'Constraint_list_category': [None,lambda x = value: returnStarLine(x,length = 127),None,True],
            #'Constraint_list_ID':       [None,returnStarInt,None,True],
            #'Constraint_list_label':    [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

            },
          },

        'Struct_image': {                                                                                     # *** Table

          'tags': {

            #'File_name':                [None,lambda x = value: returnStarLine(x,length = 127),None,True],
            #'File_format':              [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Details':                  [None,returnStarString,None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Conformer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

            },
          },

        'Atom_site': {                                                                                        # *** Table
        
          'ccpnLoop':                    ['structureEnsemble.sortedModels()','structureEnsemble.sortedCoordChains()','coordChain.sortedResidues()','coordResidue.sortedAtoms()',('coordAtom',self.getModelCoordinate)],
          'ccpnMap':                     [('model',''),('coordChain',''),('coordResidue','Label_comp_ID'),('coordAtom',''),('coordinate','ID')],
          'nmrStarKeyType':              'global', # To ensure that atom IDs increase over different 'inside' CCPN loops
                                      # Note that this can be set to, e.g. 'molStructure' to reset ID count (useful for joint primary keys)

          # This is mechanism to ADD ccpnVars to the equation, so don't have to go through data model every time
          # Necessary to speed things up!!! Also ONLY works with ccpnLoop list as above!!!
          'addCcpnMap':                 {'coordResidue': [('residue','coordResidue.residue'),('molResidue',('residue',self.getActualMolResidue)),('chain','residue.chain'),('chemComp','molResidue.chemCompVar.chemComp')],
                        },

          # NOTE: the atom_site ID increases over different models - it is the SOLE primary key for this loop according to Eldon's dictionary!

          'tags': {

            'Model_ID':                  'model.serial',  # This used to be in primary key list, not any more (Wim 11/2009)
            #'Model_site_ID':            None,
            #'Assembly_atom_ID':         None,  # TODO: can set this instead of all Label stuff, but then need to define all atoms in the molecular system in loop higher up
            'Label_entity_assembly_ID':  'chain',
            #'Label_asym_ID':            None, # What's this?
            'Label_entity_ID':           'chain.molecule',
            'Label_comp_index_ID':       'molResidue', # ADD SOMETHING HERE
            'Label_comp_ID':            ('chemComp',self.getChemCompCif),
            'Label_seq_ID':              'molResidue', # .seqCode?
            'Label_atom_ID':            ('coordAtom',self.getPdbCoordAtomName),
            'Type_symbol':               'coordAtom.elementSymbol',
            'Cartn_x':                   'coordinate.x',
            'Cartn_y':                   'coordinate.y',
            'Cartn_z':                   'coordinate.z',
            #'Cartn_x_esd':              None,
            #'Cartn_y_esd':              None,
            #'Cartn_z_esd':              None,
            'Occupancy':                 'coordinate.occupancy',
            #'Occupancy_esd':            None,
            #'Figure_of_merit':          None,
            #'Footnote_ID':              None,
            # This sets the PDB and/or author codes
            'CUSTOM_setApplDataNames':   'coordAtom',
            'Details':                   None,

            },
          },

        'Atom_sites_footnote': {                                                                              # *** Table

          'tags': {

            #'Footnote_ID':              [None,returnStarInt,None,True],
            #'Text':                     [None,returnStarString,None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Confo`rmer_family_coord_set_ID': [None,returnStarInt,'Conformer_family_coord_set.ID',True],

            },
          },

        }
      }


    self.sfDict['conformer_statistics'] = {                                                                   # *** Priority 1

      'ccpnLoop':                        'nmrEntry.sortedStructureGenerations()',
      'ccpnMap':                         'structureGeneration',

      'addCcpnMap':                     {'structureGeneration': [('structureEnsemble','structureGeneration.structureEnsemble')]},

      'title':                          ('structureGeneration.name',self.getTitle),

      'tags': {

        #'Sf_category':                  ['conformer_statistics',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Conformer_ensemble_only':       None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,True],
        'Both_ensemble_and_rep_conformer': None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,True],
        'Representative_conformer_only': None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,True],
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Text_data_format':             [None,lambda x = value: returnStarLine(x,length = 31),None,False],
        #'Text_data':                    [None,returnStarString,None,False],
        #'Original_conformer_stats_file_ID': [None,returnStarInt,None,False],
        #'Conf_family_coord_set_ID':     [None,returnStarInt,'Conformer_family_coord_set.ID',True],
        'Conf_family_coord_set_label':  ('structureGeneration.name',self.getTitle), # BMRB say this should be the same as the name of the conformer_family_coord_set save frame for 'correct' NMR-STAR
        #'Representative_conformer_ID':  [None,returnStarInt,'Representative_conformer.ID',False],
        #'Representative_conformer_label': [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
        'Conformer_calculated_total_num': "structureEnsemble.findFirstApplicationData(application='nmrStar',keyword='calculated').value", # *** Mandatory (PDB)
        'Conformer_submitted_total_num': ('structureEnsemble',self.getStrucGenSubmittedNew),                  # *** Mandatory (PDB)
        'Conformer_selection_criteria':  "structureEnsemble.findFirstApplicationData(application='nmrStar',keyword='criteria').value", # *** Mandatory (PDB)
        'Representative_conformer':      "structureEnsemble.findFirstApplicationData(application='nmrStar',keyword='representative').value", # *** Mandatory (PDB)
        'Rep_conformer_selection_criteria': "structureEnsemble.findFirstApplicationData(application='nmrStar',keyword='repr_criteria').value", # *** Mandatory (PDB) - option
        #'Statistical_struct_param_details': [None,returnStarString,None,False],
        'Details':                       'structureGeneration.details',

        },
      }


    ####################
    # Constraint lists #
    ####################

    self.sfDict['torsion_angle_constraints'] = {                                                              # *** Priority 1 - data - Not in ADIT-NMR?

      'ccpnLoop':                       ('structureGenerations',self.getDihedralConstraintLists),
      'ccpnMap':                         'constraintList',

      'title':                          ('constraintList.name',self.getTitle),

      'tags': {

        'Details':                       'constraintList.details',

        #'Sf_framecode':                 [None,returnStarString,None,False],

        },

      'tables': {

        'Torsion_angle_constraint_software': self.setSoftwareLoop('constraintList'),                          # *** Table

        #'Torsion_angle_constraint_software': {                                                                # *** Table

        #  'tags': {

            #'Software_ID':              [None,returnStarString,'Software.ID',None],
            #'Software_label':           [None,returnStarString,None,None],
            #'Method_ID':                [None,returnStarString,'Method.ID',None],
            #'Method_label':             [None,returnStarString,None,None],

        #    },
        #  },

        'Torsion_angle_constraints_expt': self.setConExperimentsLoop('constraintList.sortedExperiments()','experiment'), # *** Table

        #'Torsion_angle_constraints_expt': {                                                                   # *** Table

        #  'tags': {

        #    CJP - missing tags marked with a +

        #    'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
        #    'Experiment_name':          [None,returnStarString,None,False],
        #    'Method_ID':                [None,returnStarInt,'Method.ID',False],
        #    'Method_label':             [None,returnStarString,None,False],
        #  + 'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
        #  + 'Sample_label':             [None,returnStarString,None,True],
        #    'Sample_state':             [None,returnStarString,None,True],
        #  + 'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
        #  + 'Torsion_angle_constraint_list_ID': [None,returnStarInt,'Torsion_angle_constraint_list.ID',True],

        #    },
        #  },

        'Torsion_angle_constraint': {                                                                         # *** Table

          'ccpnLoop':                   ('constraintList',self.getAllConstraintItems),
          'ccpnMap':                     'constraintItem',

          'tags': {

            # This sets original author resonance info on constraint level
            # Sets either torsion angle name, or original atom names defining the angle
            'CUSTOM_setApplDataNames':   'constraintItem.parent',

            # CCPN attributes - not in NMR-STAR?:
            #targetValue  	Float  	0..1  	Desired value of constrained parameter  
            #error 	Float 	0..1 	Uncertainty (estimated standard deviation) of targetValue  

            'Angle_upper_bound_val':     'constraintItem.upperLimit',
            'Angle_lower_bound_val':     'constraintItem.lowerLimit',

            #'Derivation_code':          [None,returnStarString,None,None],
            #'Source_experiment_ID':     [None,returnStarString,None,None],
            #'Source_experiment_label':  [None,returnStarString,None,None],

            'CUSTOM_setAtoms':           'constraintItem.parent.resonances'

            }
          },

        'TA_constraint_comment_org': {                                                                        # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Comment_text':             [None,returnStarString,None,True],
            #'Comment_begin_line':       [None,returnStarInt,None,False],
            #'Comment_begin_column':     [None,returnStarInt,None,False],
            #'Comment_end_line':         [None,returnStarInt,None,False],
            #'Comment_end_column':       [None,returnStarInt,None,False],
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Torsion_angle_constraint_list_ID': [None,returnStarInt,'Torsion_angle_constraint_list.ID',True],

            }
          },

        'TA_constraint_parse_err': {                                                                          # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Content':                  [None,returnStarString,None,True],
            #'Begin_line':               [None,returnStarInt,None,False],
            #'Begin_column':             [None,returnStarInt,None,False],
            #'End_line':                 [None,returnStarInt,None,False],
            #'End_column':               [None,returnStarInt,None,False],
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Torsion_angle_constraint_list_ID': [None,returnStarInt,'Torsion_angle_constraint_list.ID',True],

            },
          },

        }
      }


    self.sfDict['distance_constraints'] = {                                                                   # *** Priority 1 - data - Not in ADIT-NMR?

      'ccpnLoop':                       ('structureGenerations',self.getDistanceConstraintLists),
      'ccpnMap':                         'constraintList',

      'title':                          ('constraintList.name',self.getTitle),

      'tags': {

        'Constraint_type':              ('constraintList',self.getDistanceConstraintType),
        'Details':                       'constraintList.details',

        #'Averaging_method':             [None,returnStarString,None,None],
        #'Number_monomers':              [None,returnStarInt,None,None],
        #'Pseudoatom_corrections':       [None,returnStarString,None,None],
        #'Potential_function':           [None,returnStarString,None,None],
        #'Temperature':                  [None,returnStarFloat,None,None],
        #'Exponent':                     [None,returnStarInt,None,None],
        #'Switch_distance':              [None,returnStarFloat,None,None],
        #'Asymptote_slope':              [None,returnStarFloat,None,None],
        #'B_high':                       [None,returnStarFloat,None,None],
        #'Ceiling':                      [None,returnStarFloat,None,None],
        #'Negative_offset':              [None,returnStarFloat,None,None],
        #'Scaling_constant':             [None,returnStarFloat,None,None],
        #'Function_detail':              [None,returnStarString,None,None],

        #'Sf_framecode':                 [None,returnStarString,None,False],

        },

      'tables': {

        'Distance_constraint_comment': {                                                                      # *** Table

          'tags': {

            #'Entry_ID':                 [None,returnStarString,'Entry.ID',None],
            #'Distance_constraint_list_ID': [None,returnStarString,'Distance_constraint_list.ID',None],
            #'Comment':                  [None,returnStarString,None,None],

            },
          },

        'Distance_constraint_software':  self.setSoftwareLoop('constraintList'),                              # *** Table

        #'Distance_constraint_software': {                                                                     # *** Table

        #  'tags': {

            #'Entry_ID':                 [None,returnStarString,'Entry.ID',None],
            #'Distance_constraint_list_ID': [None,returnStarString,'Distance_constraint_list.ID',None],
            #'Software_ID':              [None,returnStarString,'Software.ID',None],
            #'Software_label':           [None,returnStarString,None,None],
            #'Method_ID':                [None,returnStarString,'Method.ID',None],
            #'Method_label':             [None,returnStarString,None,None],

        #    },
        #  },

        'Distance_constraint_expt':      self.setConExperimentsLoop('constraintList.sortedExperiments()','experiment'), # *** Table

        #'Distance_constraint_expt': {                                                                         # *** Table

        #  'tags': {

        #    CJP - missing tags marked with a +

        #    'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
        #    'Experiment_name':          [None,returnStarString,None,False],
        #    'Method_ID':                [None,returnStarInt,'Method.ID',False],
        #    'Method_label':             [None,returnStarString,None,False],
        #  + 'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
        #  + 'Sample_label':             [None,returnStarString,None,True],
        #    'Sample_state':             [None,returnStarString,None,True],
        #  + 'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
        #  + 'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

        #    },
        #  },

        'Dist_constraint_tree': {                                                                             # *** Table

          #
          # When multiple ID values are set within the same loop, 
          # then need to define the tags specifically...
          #

          'ccpnLoop':                    ['constraintList.sortedConstraints()','constraint.sortedItems()'], # ccpnLoop list
          'ccpnMap':                     [('constraint','Constraint_ID'),('constraintItem','Node_ID')],
          #'ccpnLoop':                   'constraintList.sortedConstraints()',
          #'ccpnMap':                    'constraint',

          'tags': {

            # Need a tuple of values here - is defined as double key!!

            #'Constraint_ID':            'constraint.serial',
            'CUSTOM_setNodes':           'constraintItem',

            },
          },

        'Dist_constraint': {                                                                                  # *** Table

          #'ccpnLoop':                   ['constraintList.sortedConstraints()','constraint.sortedItems()'], # ccpnLoop list
          #'ccpnMap':                    [('constraint','Tree_node_member_constraint_ID'),('constraintItem','Tree_node_member_node_ID')],
          'ccpnLoop':                    ['constraintList.sortedConstraints()','constraint.sortedItems()',('constraintItem',self.getOrderedResonances)], # ccpnLoop list
          'ccpnMap':                     [('constraint','Tree_node_member_constraint_ID'),('constraintItem','Tree_node_member_node_ID'),('resonance','Constraint_tree_node_member_ID')],

          'tags': {

            #'Contribution_fractional_val': [None,returnStarFloat,None,None],

            #'Constraint_ID':            'constraint.serial',
            'CUSTOM_setAtoms':           'resonance',
            'Constraint_tree_node_member_ID': ('constraintItem',self.setConstraintTreeNodeMemberId),

            },
          },

        'Dist_constraint_value': {                                                                            # *** Table

          'ccpnLoop':                    ['constraintList.sortedConstraints()','constraint.sortedItems()'], # ccpnLoop list
          'ccpnMap':                     [('constraint','Constraint_ID'),('constraintItem','Tree_node_ID')],
          #'ccpnLoop':                   'constraintList.sortedConstraints()',
          #'ccpnMap':                    'constraint',

          'tags': {

            #'Constraint_ID':            [None,returnStarInt,'Dist_constraint_tree.Constraint_ID',True], #'constraint',
            #'Tree_node_ID':             [None,returnStarInt,'Dist_constraint_tree.Node_ID',True], #'constraintItem',
            #'Source_experiment_ID':     [None,returnStarInt,'unknown.ID',False],
            #'Spectral_peak_ID':         [None,returnStarInt,'Peak.ID',False],
            #'Intensity_val':            [None,returnStarFloat,None,False],
            #'Intensity_lower_val_err':  [None,returnStarFloat,None,False],
            #'Intensity_upper_val_err':  [None,returnStarFloat,None,False],

            'Distance_val':              'constraint.targetValue',
            'Distance_lower_bound_val':  'constraint.lowerLimit',
            'Distance_upper_bound_val':  'constraint.upperLimit',

            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

            }
          },

        'Dist_constraint_comment_org': {                                                                      # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Comment_text':             [None,returnStarString,None,True],
            #'Comment_begin_line':       [None,returnStarInt,None,False],
            #'Comment_begin_column':     [None,returnStarInt,None,False],
            #'Comment_end_line':         [None,returnStarInt,None,False],
            #'Comment_end_column':       [None,returnStarInt,None,False],
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

            }
          },

        'Dist_constraint_parse_err': {                                                                        # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Content':                  [None,returnStarString,None,True],
            #'Begin_line':               [None,returnStarInt,None,False],
            #'Begin_column':             [None,returnStarInt,None,False],
            #'End_line':                 [None,returnStarInt,None,False],
            #'End_column':               [None,returnStarInt,None,False],
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

            },
          },

        }
      }


    
    #
    # Added by Wim, 2009/09/01. Test for new distance constraint saveframe.
    #
    
    self.sfDict['general_distance_constraints'] = {                                                           # *** Priority 1 - data - Not in ADIT-NMR?

      'ccpnLoop':                       ('structureGenerations',self.getDistanceConstraintLists),
      'ccpnMap':                         'constraintList',

      'title':                          ('constraintList.name',self.getTitle),

      'tags': {

        'Constraint_type':              ('constraintList',self.getDistanceConstraintType),
        'Details':                       'constraintList.details',

        },

      'tables': {
      
        # TODO Uncomment when inserted in dictionary!

        #'Gen_dist_constraint_comment': {                                                                      # *** Table

        #  'tags': {

        #    },
        #  },

        # TODO: 

        'Gen_dist_constraint_software':  self.setSoftwareLoop('constraintList'),                              # *** Table

        'Gen_dist_constraint_expt':      self.setConExperimentsLoop('constraintList.sortedExperiments()','experiment'), # *** Table

        # Gen_dist_constraint_software_param table NOT handled!

        'Gen_dist_constraint': {                                                                              # *** Table

          #
          # When multiple ID values are set within the same loop, 
          # then need to define the tags specifically...
          #

          'ccpnLoop':                    ['constraintList.sortedConstraints()','generalConstraint.sortedItems()'], # ccpnLoop list
          'ccpnMap':                     [('generalConstraint','ID'),('generalConstraintItem','Member_ID')],

          'tags': {

            'Member_logic_code':        ('generalConstraint',self.getGeneralConstraintLogic),
            'CUSTOM_setAtoms':           'generalConstraintItem.sortedResonances()',

            'Distance_val':              'generalConstraint.targetValue',
            'Distance_lower_bound_val':  'generalConstraint.lowerLimit',
            'Distance_upper_bound_val':  'generalConstraint.upperLimit',

             #_Gen_dist_constraint.Intensity_val
             #_Gen_dist_constraint.Intensity_lower_val_err
             #_Gen_dist_constraint.Intensity_upper_val_err
             #_Gen_dist_constraint.Contribution_fractional_val
    
            },
          },

 
        'Gen_dist_constraint_comment_org': {                                                                  # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Comment_text':             [None,returnStarString,None,True],
            #'Comment_begin_line':       [None,returnStarInt,None,False],
            #'Comment_begin_column':     [None,returnStarInt,None,False],
            #'Comment_end_line':         [None,returnStarInt,None,False],
            #'Comment_end_column':       [None,returnStarInt,None,False],
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

            }
          },

        'Gen_dist_constraint_parse_err': {                                                                    # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Content':                  [None,returnStarString,None,True],
            #'Begin_line':               [None,returnStarInt,None,False],
            #'Begin_column':             [None,returnStarInt,None,False],
            #'End_line':                 [None,returnStarInt,None,False],
            #'End_column':               [None,returnStarInt,None,False],
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Distance_constraint_list_ID': [None,returnStarInt,'Distance_constraint_list.ID',True],

            },
          },
          
        # Gen_dist_constraint_parse_file NOT included!

        }
      }

    self.sfDict['RDC_constraints'] = {                                                                        # *** Priority 1 - data - Not in ADIT-NMR?

      'ccpnLoop':                       ('structureGenerations',self.getRdcConstraintLists),
      'ccpnMap':                         'constraintList',

      'title':                          ('constraintList.name',self.getTitle),

      'tags': {

        #'Dipolar_constraint_calib_method': [None,returnStarString,None,None],
        #'Mol_align_tensor_axial_sym_mol': [None,returnStarFloat,None,None],
        #'Mol_align_tensor_rhombic_mol': [None,returnStarFloat,None,None],
        #'General_order_param_int_motions': [None,returnStarFloat,None,None],
        #'Assumed_H-N_bond_length':      [None,returnStarFloat,None,None],
        #'Assumed_H-C_bond_length':      [None,returnStarFloat,None,None],
        #'Assumed_C-N_bond_length':      [None,returnStarFloat,None,None],

        'Details':                       'constraintList.details',

        #'Sf_framecode':                 [None,returnStarString,None,False],

        },

      'tables': {

        'RDC_constraint_software':       self.setSoftwareLoop('constraintList'),                              # *** Table

        #'RDC_constraint_software': {                                                                          # *** Table

        #  'tags': {

            #'Entry_ID':                 [None,returnStarString,'Entry.ID',None],
            #'RDC_constraint_list_ID':   [None,returnStarString,'RDC_constraint_list.ID',None],
            #'Software_ID':              [None,returnStarString,'Software.ID',None],
            #'Software_label':           [None,returnStarString,None,None],
            #'Method_ID':                [None,returnStarString,'Method.ID',None],
            #'Method_label':             [None,returnStarString,None,None],

        #    },
        #  },

        'RDC_constraint_expt':           self.setConExperimentsLoop('constraintList.sortedExperiments()','experiment'), # *** Table

        #'RDC_constraint_expt': {                                                                              # *** Table

        #  'tags': {

        #    CJP - missing tags marked with a +

        #    'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
        #    'Experiment_name':          [None,returnStarString,None,False],
        #    'Method_ID':                [None,returnStarInt,'Method.ID',False],
        #    'Method_label':             [None,returnStarString,None,False],
        #  + 'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
        #  + 'Sample_label':             [None,returnStarString,None,True],
        #    'Sample_state':             [None,returnStarString,None,True],
        #  + 'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
        #  + 'RDC_constraint_list_ID':   [None,returnStarInt,'RDC_constraint_list.ID',True],

        #    },
        #  },

        'RDC_constraint': {                                                                                   # *** Table

          'ccpnLoop':                    ['constraintList.sortedConstraints()','constraint.sortedItems()'], # ccpnLoop list
          'ccpnMap':                     [('constraint','ID'),('constraintItem','')],

          #
          # TODO: mapping NOT correct - nmrStar has ONLY ONE item per constraint for this!!!
          #

          'tags': {

            # Note: the resonance ordering is in this case handled on the NmrStarExport level!
            'CUSTOM_setAtoms':           'constraintItem.sortedResonances()',

            # CCPN: weight, details

            'RDC_val':                   'constraint.targetValue',
            'RDC_lower_bound':           'constraint.lowerLimit',
            'RDC_upper_bound':           'constraint.upperLimit',
            'RDC_val_err':               'constraint.error',

            #'RDC_val_scale_factor':     [None,returnStarFloat,None,None],
            #'Derivation_code':          [None,returnStarString,None,None],

            }
          },

        'RDC_constraint_comment_org': {                                                                       # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Comment_text':             [None,returnStarString,None,True],
            #'Comment_begin_line':       [None,returnStarInt,None,False],
            #'Comment_begin_column':     [None,returnStarInt,None,False],
            #'Comment_end_line':         [None,returnStarInt,None,False],
            #'Comment_end_column':       [None,returnStarInt,None,False],
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'RDC_constraint_list_ID':   [None,returnStarInt,'RDC_constraint_list.ID',True],

            }
          },

        'RDC_constraint_parse_err': {                                                                         # *** Table

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Content':                  [None,returnStarString,None,True],
            #'Begin_line':               [None,returnStarInt,None,False],
            #'Begin_column':             [None,returnStarInt,None,False],
            #'End_line':                 [None,returnStarInt,None,False],
            #'End_column':               [None,returnStarInt,None,False],
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'RDC_constraint_list_ID':   [None,returnStarInt,'RDC_constraint_list.ID',True],

            },
          },

        }
      }


    self.sfDict['J_three_bond_constraints'] = {                                                               # *** Priority 3 - data - Not in ADIT-NMR?

      'ccpnLoop':                       ('structureGenerations',self.getJCouplingConstraintLists),
      'ccpnMap':                         'constraintList',

      'title':                          ('constraintList.name',self.getTitle),

      'tags': {

        #'Sf_category':                  ['J_three_bond_constraints',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Data_file_format':             [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Constraint_file_ID':           [None,returnStarInt,'Constraint_file.ID',False],
        #'Block_ID':                     [None,returnStarInt,'Constraint_file.Block_ID',False],

        'Details':                       'constraintList.details',

        },

      'tables': {

        'J_three_bond_constraint_software':  self.setSoftwareLoop('constraintList'),                          # *** Table

        #'J_three_bond_constraint_software': {                                                                 # *** Table

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'J_three_bond_constraint_list_ID': [None,returnStarInt,'J_three_bond_constraint_list.ID',True],

        #    },
        #  },

        'J_three_bond_constraint_expt':  self.setConExperimentsLoop('constraintList.sortedExperiments()','experiment'), # *** Table

        #'J_three_bond_constraint_expt': {                                                                     # *** Table

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'J_three_bond_constraint_list_ID': [None,returnStarInt,'J_three_bond_constraint_list.ID',True],

        #    },
        #  },

        'J_three_bond_constraint': {                                                                          # *** Table

          'ccpnLoop':                    ['constraintList.sortedConstraints()','constraint.sortedItems()'],
          'ccpnMap':                     [('constraint','ID'),('constraintItem','')],

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Assembly_atom_ID_1':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_1':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_1':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_1':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_1':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_1':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_1':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_1':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Resonance_ID_1':           [None,returnStarInt,'Resonance.ID',False],
            #'Assembly_atom_ID_2':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_2':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_2':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_2':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_2':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_2':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_2':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_2':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Resonance_ID_2':           [None,returnStarInt,'Resonance.ID',False],
            #'Assembly_atom_ID_3':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_3':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_3':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_3':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_3':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_3':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_3':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_3':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Resonance_ID_3':           [None,returnStarInt,'Resonance.ID',False],
            #'Assembly_atom_ID_4':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_4':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_4':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_4':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_4':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_4':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_4':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_4':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Resonance_ID_4':           [None,returnStarInt,'Resonance.ID',False],

            'Coupling_constant_val':     'constraint.targetValue',
            'Coupling_constant_lower_bound': 'constraint.lowerLimit',
            'Coupling_constant_upper_bound': 'constraint.upperLimit',
            'Coupling_constant_err':     'constraint.error',

            #'Source_experiment_ID':     [None,returnStarInt,'unknown.ID',False],
            #'Auth_asym_ID_1':           [None,lambda x = value: returnStarCode(x,length = 12),None,False],
            #'Auth_seq_ID_1':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_asym_ID_2':           [None,lambda x = value: returnStarCode(x,length = 12),None,False],
            #'Auth_seq_ID_2':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_asym_ID_3':           [None,lambda x = value: returnStarCode(x,length = 12),None,False],
            #'Auth_seq_ID_3':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_3':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_3':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_asym_ID_4':           [None,lambda x = value: returnStarCode(x,length = 12),None,False],
            #'Auth_seq_ID_4':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_4':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_4':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'J_three_bond_constraint_list_ID': [None,returnStarInt,'J_three_bond_constraint_list.ID',True],

            'CUSTOM_setAtoms':           'constraintItem.sortedResonances()',

            },
          },

        }
      }


    self.sfDict['CA_CB_chem_shift_constraints'] = {                                                           # *** Priority 4 - data - Not in ADIT-NMR?

      #'ccpnLoop':                      ('structureGenerations',self.getCACBChemShiftConstraintLists),
      #'ccpnMap':                        'constraintList',

      #'title':                         ('constraintList.name',self.getTitle),

      'tags': {

        #'Sf_category':                  ['CA_CB_chem_shift_constraints',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Data_file_format':             [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Constraint_file_ID':           [None,returnStarInt,'Constraint_file.ID',False],
        #'Block_ID':                     [None,returnStarInt,'Constraint_file.Block_ID',False],

        'Units':                         'constraintList.unit',
        'Details':                       'constraintList.details',

        },

      'tables': {

        'CA_CB_constraint_software': {                                                                        # *** Table

          'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'CA_CB_constraint_list_ID': [None,returnStarInt,'CA_CB_constraint_list.ID',True],

            },
          },

        'CA_CB_constraint_expt': {                                                                            # *** Table

          'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'CA_CB_constraint_list_ID': [None,returnStarInt,'CA_CB_constraint_list.ID',True],

            },
          },

        'CA_CB_constraint': {                                                                                 # *** Table

          #'ccpnLoop':                   ['constraintList.sortedConstraints()','constraint.sortedItems()'],
          #'ccpnMap':                    [('constraint','ID'),('constraintItem','')],

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Assembly_atom_ID_1':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_1':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_1':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_1':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_1':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_1':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_1':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_1':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Resonance_ID_1':           [None,returnStarInt,'Resonance.ID',False],
            #'Assembly_atom_ID_2':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_2':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_2':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_2':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_2':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_2':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_2':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_2':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Resonance_ID_2':           [None,returnStarInt,'Resonance.ID',False],
            #'Assembly_atom_ID_3':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_3':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_3':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_3':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_3':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_3':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_3':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_3':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Resonance_ID_3':           [None,returnStarInt,'Resonance.ID',False],
            #'Assembly_atom_ID_4':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_4':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_4':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_4':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_4':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_4':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_4':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_4':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Resonance_ID_4':           [None,returnStarInt,'Resonance.ID',False],
            #'Assembly_atom_ID_5':       [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID_5':     [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID_5':              [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID_5':          [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID_5':                 [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID_5':                [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID_5':                [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type_5':              [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Resonance_ID_5':           [None,returnStarInt,'Resonance.ID',False],

            #'CA_chem_shift_val':        [None,returnStarFloat,None,True],
            #'CA_chem_shift_val_err':    [None,returnStarFloat,None,True],
            #'CB_chem_shift_val':        [None,returnStarFloat,None,True],
            #'CB_chem_shift_val_err':    [None,returnStarFloat,None,True],

            #'Source_experiment_ID':     [None,returnStarInt,'unknown.ID',False],
            #'Auth_asym_ID_1':           [None,lambda x = value: returnStarCode(x,length = 12),None,False],
            #'Auth_seq_ID_1':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_1':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_asym_ID_2':           [None,lambda x = value: returnStarCode(x,length = 12),None,False],
            #'Auth_seq_ID_2':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_2':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_asym_ID_3':           [None,lambda x = value: returnStarCode(x,length = 12),None,False],
            #'Auth_seq_ID_3':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_3':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_3':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_asym_ID_4':           [None,lambda x = value: returnStarCode(x,length = 12),None,False],
            #'Auth_seq_ID_4':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_4':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_4':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_asym_ID_5':           [None,lambda x = value: returnStarCode(x,length = 12),None,False],
            #'Auth_seq_ID_5':            [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID_5':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID_5':           [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'CA_CB_constraint_list_ID': [None,returnStarInt,'CA_CB_constraint_list.ID',True],

            },
          },

        }
      }


    self.sfDict['H_chem_shift_constraints'] = {                                                               # *** Priority 3 - data - Not in ADIT-NMR?

      'ccpnLoop':                       ('structureGenerations',self.getHChemShiftConstraintLists),
      'ccpnMap':                         'constraintList',

      'title':                          ('constraintList.name',self.getTitle),

      'tags': {

        #'Sf_category':                  ['H_chem_shift_constraints',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Data_file_name':               [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Data_file_format':             [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Constraint_file_ID':           [None,returnStarInt,'Constraint_file.ID',False],
        #'Block_ID':                     [None,returnStarInt,'Constraint_file.Block_ID',False],

        'Units':                         'constraintList.unit',
        'Details':                       'constraintList.details',

        },

      'tables': {

        'H_chem_shift_constraint_software': self.setSoftwareLoop('constraintList'),                           # *** Table

        #'H_chem_shift_constraint_software': {                                                                 # *** Table

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'H_chem_shift_constraint_list_ID': [None,returnStarInt,'H_chem_shift_constraint_list.ID',True],

        #    },
        #  },

        'H_chem_shift_constraint_expt':  self.setConExperimentsLoop('constraintList.sortedExperiments()','experiment'), # *** Table

        #'H_chem_shift_constraint_expt': {                                                                     # *** Table

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'H_chem_shift_constraint_list_ID': [None,returnStarInt,'H_chem_shift_constraint_list.ID',True],

        #    },
        #  },

        'H_chem_shift_constraint': {                                                                          # *** Table

          'ccpnLoop':                    'constraintList.sortedConstraints()',
          'ccpnMap':                     'constraint',

          'tags': {

            #'ID':                       [None,returnStarInt,None,True],
            #'Assembly_atom_ID':         [None,returnStarInt,'Atom.Assembly_atom_ID',False],
            #'Entity_assembly_ID':       [None,returnStarInt,'Entity_assembly.ID',True],
            #'Entity_ID':                [None,returnStarInt,'Entity.ID',True],
            #'Comp_index_ID':            [None,returnStarInt,'Entity_comp_index.ID',True],
            #'Seq_ID':                   [None,returnStarInt,'PDBX_poly_seq_scheme.Seq_ID',True],
            #'Comp_ID':                  [None,lambda x = value: returnStarCode(x,length = 12),'Chem_comp.ID',True],
            #'Atom_ID':                  [None,lambda x = value: returnStarAtCode(x,length = 12),'Chem_comp_atom.Atom_ID',True],
            #'Atom_type':                [None,lambda x = value: returnStarCode(x,length = 15),None,True],
            #'Atom_isotope_number':      [None,returnStarInt,None,False],
            #'Resonance_ID':             [None,returnStarInt,'Resonance.ID',False],

            'Chem_shift_val':            'constraint.targetValue',
            'Chem_shift_val_err':        'constraint.error',

            #'Source_experiment_ID':     [None,returnStarInt,'unknown.ID',False],
            #'Auth_asym_ID':             [None,lambda x = value: returnStarCode(x,length = 12),None,False],
            #'Auth_seq_ID':              [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_comp_ID':             [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Auth_atom_ID':             [None,lambda x = value: returnStarCode(x,length = 15),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'H_chem_shift_constraint_list_ID': [None,returnStarInt,'H_chem_shift_constraint_list.ID',True],

            'CUSTOM_setAtoms':           'constraint.resonance',

            },
          },

        }
      }


    self.sfDict['other_constraints'] = {                                        #                             # *** Priority 3 - data - Not in ADIT-NMR?

      'ccpnLoop':                       ('structureGenerations',self.getOtherConstraintLists),
      'ccpnMap':                         'constraintList',

      'title':                          ('constraintList.name',self.getTitle),

      'tags': {

        #'Sf_category':                  ['other_constraints',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        #'Entity_assembly_ID':           [None,returnStarInt,'Entity_assembly.ID',True],
        #'Entity_ID':                    [None,returnStarInt,'Entity.ID',True],
        #'Type':                         [None,lambda x = value: returnStarLine(x,length = 127),None,True],
        #'Subtype':                      [None,lambda x = value: returnStarLine(x,length = 127),None,True],
        #'File_name':                    [None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'File_format':                  [None,lambda x = value: returnStarLine(x,length = 127),None,True],
        #'Constraint_file_ID':           [None,returnStarInt,'Constraint_file.ID',False],
        #'Block_ID':                     [None,returnStarInt,'Constraint_file.Block_ID',False],
        #'Text':                         [None,returnStarString,None,True],

        'Details':                       'constraintList.details',

        },

      'tables': {

        'Other_constraint_software':     self.setSoftwareLoop('constraintList'),                              # *** Table

        #'Other_constraint_software': {                                                                        # *** Table

        #  'tags': {

            #'Software_ID':              [None,returnStarInt,'Software.ID',True],
            #'Software_label':           [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Other_constraint_list_ID': [None,returnStarInt,'Other_constraint_list.ID',True],

        #    },
        #  },

        'Other_constraint_expt':         self.setConExperimentsLoop('constraintList.sortedExperiments()','experiment'), # *** Table

        #'Other_constraint_expt': {                                                                            # *** Table

        #  'tags': {

            #'Experiment_ID':            [None,returnStarInt,'unknown.ID',True],
            #'Experiment_name':          [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Method_ID':                [None,returnStarInt,'Method.ID',False],
            #'Method_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,False],
            #'Sample_ID':                [None,returnStarInt,'Sample.ID',True],
            #'Sample_label':             [None,lambda x = value: returnStarLabel(x,length = 127),None,True],
            #'Sample_state':             [None,lambda x = value: returnStarLine(x,length = 31),None,True],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Other_constraint_list_ID': [None,returnStarInt,'Other_constraint_list.ID',True],

        #    },
        #  },

        }
      }


    # The next four entries are not read because the ccpnMap key has been commented out.
    # They are here for reference.

    self.sfDict['constraint_statistics'] =  {                                                                 # *** Priority 1

      #'ccpnLoop':                       ('structureGenerations',self.getAllConstraintLists),
      #'ccpnMap':                         'constraintList',

      'ccpnMap':                         'nmrEntry',

      'title':                          ('nmrEntry.name',self.getTitle),

      'tags': {
        #'Entry_ID':                     [None,returnStarString,'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        },

      'tables': {

        'Constraint_file': {                                                                                  # *** Includes optional (PDB) fields

          #'ccpnLoop':                    'constraintList.sortedConstraints()',
          #'ccpnMap':                     'constraint',

          'ccpnLoop':                   ('structureGenerations',self.getAllConstraintLists),
          'ccpnMap':                     'constraintList',

          'tags': {

            #'CUSTOM_setConsStat':        'constraintList', # TODO - does nothing - not in NmrStarExport!

            'ID':                       ('constraintList', self.setSerial),
            'Constraint_filename':      ('constraintList.root.name',self.getProjectFileName),                 # *** Optional (PDB)
            'Software_ID':               'constraintList.method.software',
            'Software_name':             'constraintList.method.software.name',
            'Software_label':            'constraintList.method.software.name',                               # *** Optional (PDB)
            #'Software_label=LOCAL':      'CCPN',
            #'Block_ID':                 [None,returnStarInt,None,False],
            'Constraint_type':          ('constraintList', self.getStarConstraintType),                       # *** Optional (PDB) - option
            'Constraint_subtype':       ('constraintList', self.getStarSubConstraintType),                    # *** Optional (PDB) - option
            'Constraint_subsubtype':    ('constraintList', self.getStarSubSubConstraintType),                 # *** Optional (PDB) - pulldown
            'Constraint_number':        ('constraintList', self.getConstraintListCount),                      # *** Optional (PDB)
            #'Entry_ID':                 [None,returnStarString,'Entry.ID',True],
            #'Constraint_stat_list_ID':  [None,returnStarInt,'Constraint_stat_list.ID',True],

            },
          },

        }
      }


    self.sfDict['org_constr_file_comment'] =  {                                                               # *** Priority 1 - Not in ADIT-NMR?

      'ccpnLoop':                       ('structureGenerations',self.getDistanceConstraintLists),
      #'ccpnMap':                        'constraintList',

      'title':                          ('constraintList.name',self.getTitle),

      'tags': {

        #'Entry_ID':                     [None,returnStarString,'Entry.ID',True],
        #'ID':                           [None,returnStarInt,None,True],
        'Constraint_file_ID':            None, #[None,returnStarInt,'Constraint_file.ID',False],
        'Block_ID':                      None, #[None,returnStarInt,'Constraint_file.Block_ID',False],
        'Details':                       None, #[None,returnStarString,None,False],
        'Comment':                       None, #[None,returnStarString,None,False],

        },
      }


    self.sfDict['entry_interview'] = {                                                                        # *** Priority 1

      #'ccpnMap':                        None,

      'title':                           None,

      'tags': {

        #'Sf_category':                  ['entry_interview',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        'ID':                            None, #[None,returnStarInt,None,True],
        'PDB_deposition':                None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
        'BMRB_deposition':               None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
        'View_mode':                     None, #[None,lambda x = value: returnStarYesNo(x,length = 15),None,False], # *** Mandatory - bool
        'Structural_genomics':           None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], # *** Optional - bool - auto
        'Ligands':                       None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], # *** Mandatory - bool - auto
        'Non_standard_residues':         None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Constraints':                   None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
        'Assigned_chem_shifts':          None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool (deals with assigned_chemical_shifts and chem_shift_reference save frames)
        'Coupling_constants':            None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Chem_shift_anisotropy':         None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Heteronucl_NOEs':               None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Heteronucl_T1_relaxation':      None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Heteronucl_T2_relaxation':      None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Heteronucl_T1rho_relaxation':   None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Order_parameters':              None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Residual_dipolar_couplings':    None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'H_exchange_rate':               None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'H_exchange_protection_factors': None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Spectral_peak_lists':           None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Dipole_dipole_couplings':       None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
        'Quadrupolar_couplings':         None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
        'Homonucl_NOEs':                 None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Dipole_dipole_relaxation':      None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'DD_cross_correlation':          None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Dipole_CSA_cross_correlation':  None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'PKa_value_data_set':            None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto (deals with pH_titration and pH_param_list)
        'D_H_fractionation_factors':     None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Theoretical_chem_shifts':       None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto (deals with theoretical_chem_shifts and chem_shifts_calc_type)
        'Spectral_density_values':       None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto
        'Timedomain_data':               None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
        'Molecular_interactions':        None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
        'Secondary_structure_orientations': None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
        'Metabolite_coordinates':        None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
        'Mass_spec_data':                None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
        'Other_kind_of_data':            None, #[None,lambda x = value: returnStarYesNo(x,length = 3),None,False], *** Mandatory - bool - auto

        },
      }


    self.sfDict['deposited_data_files'] =  {                                                                  # *** Priority 1 - Not as a form in ADIT-NMR

      #'ccpnMap':                        None,

      'title':                           None,

      'tags': {

        #'Sf_category':                  ['deposited_data_files',lambda x = value: returnStarCode(x,length = 31),None,True],
        #'Sf_framecode':                 [None,lambda x = value: returnStarCode(x,length = 127),None,False],
        #'Entry_ID':                     [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
        'ID':                            None, #[None,returnStarInt,None,True],
        'Atomic_coordinate_file_name':   None, #[None,lambda x = value: returnStarLine(x,length = 127),None,False],
        'Atomic_coordinate_file_syntax': None, #[None,lambda x = value: returnStarLine(x,length = 127),None,False],
        'Constraint_file_name':          None, #[None,lambda x = value: returnStarLine(x,length = 127),None,False],
        'Constraint_file_syntax':        None, #[None,lambda x = value: returnStarLine(x,length = 127),None,False],
        #'Precheck_flag':                [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],
        #'Validate_flag':                [None,lambda x = value: returnStarYesNo(x,length = 3),None,False],

        },

      'tables': {

        'Upload_data': {                                                                                      # *** Table

          'tags': {

            #'Data_file_ID':             [None,returnStarInt,None,True],
            #'Data_file_name':           [None,lambda x = value: returnStarLine(x,length = 127),None,True],
            #'Data_file_Sf_category':    [None,lambda x = value: returnStarCode(x,length = 31),None,False],
            #'Data_file_syntax':         [None,lambda x = value: returnStarLine(x,length = 127),None,False],
            #'Entry_ID':                 [None,lambda x = value: returnStarCode(x,length = 12),'Entry.ID',True],
            #'Deposited_data_files_ID':  [None,returnStarInt,'Deposited_data_files.ID',True],

            },
          },

        }
      }

#
# Here subclass based on versions
#

class Ccpn_To_NmrStar_test(Ccpn_To_NmrStar):

  def someNewFunction(self):
  
    # Can here define some new method that is necessary
    pass
  
  def modifyExistingFunction(self):
  
    # Can here modify (reset, in effect) an existing method
    # For example, def getSpectrometers(), in case some changes on CCPN side...
    
    pass
    
  def modifySfDict(self):
    
    # Here can make small changes to the mapping dictionary, if required
    # Major differences would probably be handled better by defining a new setSfDict() method.
    
    # For an example of minor changes, say that the CCPN class name for ShiftList changes to
    # ChemicalShiftList. The original value in the 'assigned_chemical_shifts' saveframe for
    # the 'ccpnLoop' key is:
    
    #"nmrEntry.findAllMeasurementLists(className='ShiftList', isSimulated=False)",
    
    # This would have to be fixed in the following way:

    #self.sfDict['assigned_chemical_shifts']['ccpnLoop'] = "nmrEntry.findAllMeasurementLists(className='ChemicalShiftList', isSimulated=False)",

    print "  Modifying mapping dictionary... "

    pass
 

#
# Version definitions...   
#
# See Constants.py importVersionSep for replacement string for '.' in version numbers (has to be '_', though)
#

Ccpn__1_1_2__To_NmrStar__3_0__  = Ccpn_To_NmrStar
Ccpn__1_1_a2__To_NmrStar__3_0__ = Ccpn_To_NmrStar
Ccpn__1_1_2__To_NmrStar__3_1__  = Ccpn_To_NmrStar
Ccpn__1_1_a2__To_NmrStar__3_1__ = Ccpn_To_NmrStar
Ccpn__1_1_a3__To_NmrStar__3_0__ = Ccpn_To_NmrStar
Ccpn__1_1_a3__To_NmrStar__3_1__ = Ccpn_To_NmrStar

Ccpn__2_0_b3__To_NmrStar__3_0__ = Ccpn_To_NmrStar
Ccpn__2_0_b3__To_NmrStar__3_1__ = Ccpn_To_NmrStar
Ccpn__2_0_b3__To_NmrStar__3_1_test__ = Ccpn_To_NmrStar_test
