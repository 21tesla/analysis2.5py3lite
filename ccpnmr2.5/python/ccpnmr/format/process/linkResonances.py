"""
======================COPYRIGHT/LICENSE START==========================

linkResonances.py: Links resonances to atoms in the data model

Copyright (C) 2008 Wim Vranken (European Bioinformatics Institute)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license
 
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
- PDBe website (http://www.ebi.ac.uk/pdbe/)

- contact Wim Vranken (wim@ebi.ac.uk)
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

import os, sys, string

from memops.universal.Util import drawBox

from ccp.general.Util import setUniqueAppData, findAllSysNamesByChemAtomOrSet
from ccp.general.Util import findChemAtomOrSetSysName
from ccp.general.Util import getDescriptorDict, mergeDescriptorDicts, getDescriptorFromDict

# Import general functions
from ccpnmr.format.general.userInteraction import setupMessageReporter
from ccpnmr.format.general.userInteraction import setupMultiDialog

from ccpnmr.format.general.Util import getAtomSet, updateResonanceNamesDict
from ccpnmr.format.general.Util import getNameInfo, createNonStereoName
from ccpnmr.format.general.Util import getStereoCode, setCcpNmrResonanceName
from ccpnmr.format.general.Util import getNmrIsotope, getApplResNames
from ccpnmr.format.general.Util import returnProchiralAtom, returnProchiralAtomGroup
from ccpnmr.format.general.Util import getResName, getResNameText
from ccpnmr.format.general.Util import SimpleError

from ccpnmr.format.general.TopShared import TopShared

from ccpnmr.format.general.Constants import chemShiftAllowedDiff, namingSystemNames

from ccp.general.Copy import copyResonanceInfo
from ccp.format.general.Constants import defaultSeqInsertCode
from ccp.format.general.Util import getRegularExpressions

import memops.api.Implementation as Implementation
from ccp.api.nmr import Nmr, NmrConstraint
from ccp.api.molecule import ChemComp
from ccpnmr.api import NameMapping

from ccpnmr.format.general.Constants import ccpNmr_kw, assign_kw, ccpNmrMapping_kw
from ccpnmr.format.general.Constants import allResidueAtoms_kw

#####################################
# Classes for easier access to data #
#####################################    

class ResidueAtomToResName:

  """
  
  Fake class similar to AtomToResName, but tracks residue level resonances (for HADDOCK)
  
  """
  def __init__(self,nmrAtoms,resName, verbose = True):
  
    self.verbose = verbose
  
    self.nmrAtoms = list(nmrAtoms)
    self.nmrAtoms.sort()
        
    self.resNameList = []
    if resName:
      self.resNameList.append(resName)

    self.setAllAtomGroups()

  def setAllAtomGroups(self):
  
    #
    # Group the atoms together
    #
  
    self.allAtomGroups = {}
    chemCompVar = self.nmrAtoms[0].residue.chemCompVar
  
    for nmrAtom in self.nmrAtoms:
      
      chemAtomSet = nmrAtom.chemAtom.chemAtomSet
      
      if chemAtomSet:
        chemAtomOrSetName = chemAtomSet.name
      else:
        chemAtomOrSetName = nmrAtom.name
      
      if not self.allAtomGroups.has_key(chemAtomOrSetName):
        self.allAtomGroups[chemAtomOrSetName] = []
      
      self.allAtomGroups[chemAtomOrSetName].append(nmrAtom)

class AtomToResName:

  """
  
  This class keeps track of which atoms will be linked to which resName
  (internal linkResonances use only)

  Since there might be multiple resNames that could be linked to the same
  atoms, this class is an intermediate to hold this information. During the
  actual linking of resonances to atoms the information in here will be used
  to delete/create/copy resonances (or the information linked to them).
  
  """

  def __init__(self,nmrAtoms,chemAtomDict,resName,chemAtomSets,isotopeCode,singleNmrAtom = None, verbose = True):
  
    self.verbose = verbose
  
    self.nmrAtoms = nmrAtoms
    self.nmrAtoms.sort()
    
    self.chemAtomDict = chemAtomDict
    
    self.isotopeCode = isotopeCode
    
    self.resNameList = []
    if resName:
      self.resNameList.append(resName)

    self.eqChemAtomSets = {}
    self.prochChemAtomSets = {}
    self.groupChemAtomSets = {}
    
    self.hasProchChemAtomSets = []
    self.hasEqChemAtomSets = []
    self.hasEmptyChemAtomSets = []
    self.hasGroupChemAtomSets = []

    self.singleAtoms = {}
    self.stereoAtoms = {}
    
    self.setNmrAtomGroups()
    
    chemAtomSets = list(chemAtomSets)
    
    self.setProchiralStatus(chemAtomSets,nmrAtoms[0].residue.chemCompVar)
    
    self.setEquivalentStatus(singleNmrAtom)

    self.setChemAtomSets(chemAtomSets,nmrAtoms)
        
  def addResName(self,resName,chemAtomSets,nmrAtoms,singleNmrAtom = None):

    if resName and resName not in self.resNameList:
      self.resNameList.append(resName)
      self.setChemAtomSets(chemAtomSets,nmrAtoms)
      self.setEquivalentStatus(singleNmrAtom)
    
  def setNmrAtomGroups(self):
  
    #
    # Group the nmr atoms together
    #
  
    tempDict = {}
    self.nmrAtomGroups = []
  
    for nmrAtom in self.nmrAtoms:
      
      chemAtomSet = self.chemAtomDict[nmrAtom].chemAtomSet
      
      if chemAtomSet:
        if not tempDict.has_key(chemAtomSet):
          tempDict[chemAtomSet] = []
        
        tempDict[chemAtomSet].append(nmrAtom)
    
    if tempDict:
      for chemAtomSet in tempDict.keys():
        self.nmrAtomGroups.append(tempDict[chemAtomSet])

  def setChemAtomSets(self,chemAtomSets,nmrAtoms):
  
    resName = self.resNameList[-1]

    self.prochChemAtomSets[resName] = []
    self.eqChemAtomSets[resName] = []
    self.groupChemAtomSets[resName] = []
    
    if self.isProchiral:
    
      #
      # Keep track of stereo assigned ones
      #
        
      if nmrAtoms and len(nmrAtoms) < len(self.nmrAtoms):
        if self.verbose:
          print "  Stereo assigned prochiral %s" % resName
        self.stereoAtoms[resName] = nmrAtoms

    if chemAtomSets:
    
      if chemAtomSets[0].isProchiral:

        self.prochChemAtomSets[resName] = chemAtomSets
        self.hasProchChemAtomSets.append(resName)
        
      else:

        if chemAtomSets[0].isEquivalent == True:
          self.eqChemAtomSets[resName] = chemAtomSets
          self.hasEqChemAtomSets.append(resName)
        
        else:

          #
          # Aromatic (F,Y) or amide (N,Q)
          #
          
          if len(nmrAtoms) > 1:
            
            #
            # Grouped aromatic (HE*,HD*,HD*|HE*)
            # If 'deep' group (HD*|HE*), then split up at this stage
            #
            
            if len(chemAtomSets) == 1 and chemAtomSets[0].chemAtomSets:
              chemAtomSets = chemAtomSets[0].chemAtomSets
            
            self.groupChemAtomSets[resName] = list(chemAtomSets)
            self.hasGroupChemAtomSets.append(resName)
        
    else:
    
      self.hasEmptyChemAtomSets.append(resName)

  def setProchiralStatus(self,chemAtomSets,chemCompVar):

    self.isProchiral = 0
    
    if chemAtomSets:

      if chemAtomSets[0].isProchiral == True:
        self.isProchiral = 1

      elif chemAtomSets[0].chemAtomSet and chemAtomSets[0].chemAtomSet.isProchiral == True:
        self.isProchiral = 1
    
    else:
     
      chemAtomSet = self.chemAtomDict[self.nmrAtoms[0]].chemAtomSet
      
      if chemAtomSet and chemAtomSet.isProchiral == True:
        self.isProchiral = 1

  def setEquivalentStatus(self,singleNmrAtom):
  
    self.isEquivalent = -1

    if singleNmrAtom:
    
      self.singleAtoms[self.resNameList[-1]] = singleNmrAtom
      self.isEquivalent = singleNmrAtom.chemAtom.chemAtomSet.isEquivalent
      
    else:
      
      chemAtomSet = self.nmrAtoms[0].chemAtom.chemAtomSet
      
      if chemAtomSet:
        self.isEquivalent = chemAtomSet.isEquivalent

class ExistingResonances:

  """
  
  This class keeps track of resonances that are already linked to atoms
  (internal linkResonances use only)
  
  """

  def __init__(self):
    
    self.resonances = []
    self.atoms = []
   
  def addResonance(self,resonance):
  
    self.resonances.append(resonance)
    self.atoms.append([])
    
  def getIndex(self,resonance):
  
    return self.resonances.index(resonance)
    
  def addAtom(self,resonance,atom):
  
    self.atoms[self.getIndex(resonance)].append(atom)
    
  def determineMatches(self,nmrAtoms):
     
    self.numResonances = len(self.resonances)

    #
    # Reorganize: is there a full nmrAtom match or a submatch?
    #
    
    self.fullAtomMatch = []
    self.subAtomMatch = []

    for resonance in self.resonances:
    
      if len(self.atoms[self.getIndex(resonance)]) == len(nmrAtoms):
        self.fullAtomMatch.append(resonance)
      else:
        self.subAtomMatch.append(resonance)

  def determineResonanceSets(self):
  
    self.resonanceSets = []
    
    for resonance in self.resonances:
      resonanceSet = resonance.resonanceSet
      
      if self.resonanceSets.count(resonanceSet) == 0:
        self.resonanceSets.append(resonanceSet)
        
  def setResNames(self,assignFormat):
  
    self.resNames = []
    self.ccpNmrResNames = []
    
    for resonance in self.resonances:
      resNameDict = getApplResNames(assignFormat,[resonance], keyword = assign_kw)
      ccpNmrResNameDict = getApplResNames(ccpNmr_kw,[resonance], keyword = assign_kw)
      
      self.resNames.append(resNameDict.keys())      
      self.ccpNmrResNames.append(ccpNmrResNameDict.keys())      

########
# Main #
########

class linkResonances(TopShared):

  """
  
  Linkresonances class, use as follows:
  
  linkResonances = LinkResonances()
  linkResonances.run(project, keywords)
  
    project    Obligatory, a valid CCPN project
    keywords   Keyword arguments, see keywordInfo below
    
  """
  
  keywordInfo = {
                 
    'guiParent':                 (None,  False,"""
        Parent for GUI interface (if set to None, will use text-based interaction)."""),
        
    'messageReporter':           (None,  False,"""
        MessageReporter instance (can be passed in if already created)."""),
        
    'multiDialog':               (None,  False,"""
        MultiDialog instance (can be passed in if already created)."""),
        
    'nameMapping':               (None,  False,"""
        NameMapping instance (can be passed in if already created)."""),   
        
    'assignFormat':              (None,  False,"""
        Name of the format from where the resonances were imported (e.g 'xeasy'). Note that this is not necessarily the same as the atom naming system."""),
        
    'namingSystemName':          (None,  False,"""
        Set the name of the naming system that is used to interprete the atom names (e.g. 'CYANA', 'XPLOR', ...)."""),
        
    'globalStereoAssign':        (False, False,"""
        False: Do not assume anything about stereospecific assignment status of prochiral atoms (depend on user interaction).
        True:  Assume that all prochiral atoms have been stereospecifically assigned."""),
        
    'minimalPrompts':            (False, False,"""
        True:  Use automatic settings as much as possible (this will greatly minimize popups but is less user-oriented).
        False: Interact with the user wherever necessary."""),
        
    'molSystem':                 (None,  False,"""
        A valid CCPN MolSystem object (will be used to get chain information)."""),
        
    'makeMethylsEquivalent':     (True,  False,"""
        True:   Methyl protons are equivalent, e.g. make HG1* if only HG11 given for valine.
        False:  Methyl portons are not equivalent, e.g. in ILV labelling experiments."""),
           
    'setSingleProchiral':        (None,  False,"""
        None:  Do not assume anything on the status of single prochiral atoms (e.g. ASP HB2, but no info on HB3).
        False: Assume that for single prochiral atoms the other prochiral atom carries exactly the same information.
        True:  Assume that for single prochiral atoms the other prochiral atom is not known - information is not copied."""),
        
    'setSinglePossEquiv':        (None,  False,"""
        None:  Do not assume anything on the status of single possibly equivalent atoms (e.g. PHE HD1, but no info on HD2).
        False: Assume that for single possibly equivalent atoms the other atom is equivalent.
        True:  Assume that for single possibly equivalent atoms the other atom is not equivalent."""),
        
    'useCommonNames':            (0,     False,"""
         0:  Do not use 'common names' information.
         1:  Translate 'common names' information; this will translate commonly mistyped atom names using their expected meaning.
        -1:  Ignore 'common names'; this will not link any atoms listed as common names."""),
        
    'useAmbiguity':              (0,     False,"""
         0:  Do not use ambiguity name information.
         1:  Translate ambiguous names; this will translate ambiguous atom names (e.g. HG* for THR) using their expected meaning.
        -1:  Ignore ambiguous names; this will not link any atoms listed as common names."""),
        
    'useIupacMatching':          (0,     False,"""
        0: Do not use IUPAC names for mapping in addition to selected naming system
        1: Use IUPAC names for mapping: if no match is found using the selected naming system the IUPAC names will be tried, but ONLY if there is no chemAtomSet for the match found (this is to avoid problems with HB2/HB3 type mapping).
        2: Use IUPAC names for mapping: if no match is found using the selected naming system the IUPAC names will be tried, and will always be applied if a match found."""),

    'useCcpnAtomMatching':       (False,False,"""
        False: Do not try to map directly to CCPN atom names in addition to selected naming system
        True:  If no match is found, try to map the atom name directly to the CCPN atom (set) names."""),

    'useEmptyNamingSystems':     (True, False,"""
        True:  Always try to use the default naming system, even if no atomSysNames are available (as long as it is referenced to one that is, it will be used).
        False: Do not use the default naming system if it does not have atomSysNames available."""),

    'useLinkResonancePopup':     (True,  False,"""
        True:  If no atom match was found for a particular resonance, display a popup to determine which atom(s) the resonance corresponds to.
        False: Do not ask user to match resonances to atoms if no match was found for a particular resonance."""),

    'forceShiftMerge':           (False, False,"""
        False: Ask user to determine which chemical shift value is correct when there are multiple values for one resonance.
        True:  Will automatically merge ambiguous shift values (takes the average shift value)."""),
                            
    'strucGen':                  (None,  False,"""
        A valid CCPN Nmr.StructureGeneration object. If passed in, will use the related NmrConstraintStore to set the list of fixed resonances."""),
        
    'nmrConstraintStore':        (None,  False,"""
        A valid CCPN NmrConstraint.NmrConstraintStore object. If passed in, will run linkResonances on the list of connected FixedResonances.    """),

    'nmrProject':               (None,  False,"""
        A valid CCPN Nmr.NmrProject object. If passed in, will run linkResonances on the list of connected Resonances.    """),
            
    'fixSingleAtoms':            (True,  False,"""
        False: Will not try to find out if an unrecognized atom set (e.g. HG# for LEU) is available as single atom (HG in this case).
        True:  Will reset unrecognized atom sets to single atom if relevant (e.g. HG# to HG for LEU)."""),
        
    'verbose':                   (True,  False,"""                            
        True:  Print out all detailed information.
        False: Do not print out detailed information.
        """),
        
    'verboseCopy':                (False,  False,"""                            
        True:  Print out detailed information for resonance copying.
        False: Do not print out detailed information for resonance copying.
        """),

    #
    # Advanced keywords for full automation only, use with care! 
    #

    'useResLinkMapping':         (True,  False,"""
        True:  Use stored resonance to atom mappings (from NameMapping).
        False: Do not use stored resonance to atom mappings."""),
        
    'useSingleProchiralMapping': (False, False,"""
        False: Do not use stored single prochiral atom settings (from NameMapping).
        True:  Use stored single prochiral atom settings."""),

    'useSinglePossEquivMapping': (False, False,"""
        False: Do not use stored possible equivalent atom settings (from NameMapping).
        True:  Use stored possible equivalent atom settings."""),

    'forceChainMappings':        (None,  False,"""
        Can be set to a Python list so that chains are mapped automatically.
        Takes the form [[chainCode, formatChainCode, firstSeqId, offset],...].
        Example: [['A','B',1,4]] maps the CCPN chain 'A' to the external format chain 'B'. Residue 1 of the CCPN chain maps to sequence code 5 (1+4) for the relevant chain in the external format."""),
        
    'forceDefaultChainMapping':  (False, False,"""
        Can be set to True to use the default chain mapping (works only for 1 CCPN chain, 1 format chain)"""),

    'ignoreResidues':            (None,  False,"""
        Can be set to a list of CCPN residues for the project - these will be ignored when linking resonances.
        This is useful when linking information from a highly homologous sequence where you want to ignore the residues that do not match.
        Takes the form [ccpnResidue1,ccpnResidue2,...]."""),
        
    'addNameMappings':           (None,  False,"""
        Can be set to a Python dictionary so that specific atom names will be mapped automatically to the internal CCPN atom (set) name.
        Example: {"ALL": [['1H', 'H1'],['2H','H2']], "Thr": [['QG', 'HG2*']]} will map the atom names 1H to H1 and 2H to H2 for all residues, and for Thr residues only will map QG to HG2*."""),
        
    'formatChainCodes':          (None,  False,"""
        Can be set to a dictionary to copy/change the original format chain code to new/existing chain codes (can be necessary for multimers, ...)
        Example: {' ': ['A', 'B']} will change the original chain code ' ' to 'A', and make a copy of all relevant resonances to original chain 'B'."""),
        
    'autoConnectStereo':         (True,  False,"""
        Set to False if unlinked ambiguous resonances should not be automatically connected based on atom names to existing stereospecific resonances."""),
        
    'specificResNameMappings':   (None,  False,"""
        Can be set to a dictionary so that resNames are automatically changed. This is done before any atom name interpretations, so the new identifying string has to conform to the external format names.
        Example: {' .1.H': ' .1.HN'} will convert ONLY the exact resonance identifier string ' .1.H' (from the external format) to ' .1.HN'."""),
        
    'allowPopups':               (True,  False,"""
        True:  Allow user popups to appear.
        False: Do not show any user popups, run everything fully automatically even in case of ambiguities in the data."""),
  }
    
  tempNamingSystemNames = ('PDB_REMED','IUPAC','MSD','PDB')

  checkProtonationStates = {

  #'Asn': ('HD21','HD22'),
  'Asp': ('HD2',),
  #'Gln': ('HE21','HE22'),
  'Glu': ('HE2',),
  'His': ('HD1','HE2'),
  'Lys': ('HZ3',),
  'Arg': ('HH12',),
  'Ser': ('HG',),
  'Thr': ('HG1',),
  'Tyr': ('HH',),
  'Trp': ('HE1',),

  }

  
  LinkResonancesError = SimpleError
  
  def __init__(self,proj, **keywds):

    self.patt = getRegularExpressions()

    self.run(proj,keywds)
    
  def run(self,project,keywds):
    
    #
    # Check keyword information
    #
    
    if not self.setKeywords(keywds):
      return

    #
    # Set up GUI or non-GUI user interaction if not passed in
    #
      
    if not self.messageReporter:
      self.messageReporter = setupMessageReporter(self.guiParent)
      
    if not self.multiDialog:
      self.multiDialog = setupMultiDialog(self.guiParent)
    
    self.status = True

    #
    # Set globalStereoAssign to 1 if everything is stereospecifically assigned
    # Otherwise the script will try to use any information available (eg the nmrStar ambiguitycode)
    #
  
    ##############
    # Initialize #
    ##############

    if not project.molSystems:
      self.messageReporter.showError("Error","Unable to link resonances: need project with molecule data!",parent = self.guiParent)
      self.status = False
      return

    self.project = project
    self.format = self.assignFormat

    if not self.nmrConstraintStore and self.strucGen:
      self.nmrConstraintStore = self.strucGen.nmrConstraintStore
      
    if not self.nmrConstraintStore and not self.nmrProject:
      self.nmrProject = self.project.currentNmrProject
      
    if not self.nmrConstraintStore and not self.nmrProject:
      raise self.LinkResonancesError("No nmrConstraintStore or nmrProject passed in, and no currentNmrProject available. Aborting!")

    
    if self.forceChainMappings:
      self.forceDefaultChainMapping = 0
      
    if not self.ignoreResidues:
      self.ignoreResidues = []
    
    self.resAssignedNames = {}

    self.resUnassignedNames = {}
    self.resUnassignedNameList = []
    self.missingAtomNames = 0

    self.ccpChainDict = {}
    self.formatChains = {}
    self.formatChainCodeDict = {}
    self.chainMappingChain = {}

    self.atomResonances = {}

    self.atomToResNameList = []

    
    self.resonance = None
    self.numResonancesLinked = 0

    self.atomsFound = []
    self.chemAtoms = []
    self.chemAtomSets = None
    
    self.chainCode = None
    self.seqCode = None
    self.chain = None
    self.nmrRes = None
    
    self.ccvSysNames = {}
    
    self.ccNamingSystems = {} # Also track if naming system exists for particular chemComp

    #######################################################
    # Main calls to functions (should go somewhere else?) #
    #######################################################

    #
    # Print message...
    #
    
    print
    print drawBox("Starting linkResonances")
    print

    #
    # Print message...
    #
    
    print
    print drawBox("Stage 1: Determining input format and setting up reference data",indent = "  ")
    print
    
    #
    # Set molSystem: handling only one at a time!
    #
    
    if not self.molSystem:
      self.selectMolSystem()
    
    #
    # Check for nameMapping
    #

    self.checkNameMapping()

    #
    # Do a check on the assignFormat
    #
    
    self.checkAssignFormat()
    
    #
    # Here format chain codes can be changed and/or copied...
    #
    
    self.changeFormatChainCodes(self.formatChainCodes)
    
    #
    # Set up reference data for the resonances.
    # 

    self.setupRefData()
    
    if not self.resUnassignedNameList:
      if not self.minimalPrompts:
        self.messageReporter.showWarning("Warning","No unassigned resonances... aborting linkResonances",parent = self.guiParent)
      return False
      
    #
    # Print message...
    #
    
    print
    print drawBox("Stage 2: Determining chain/residue mapping between format data and CCPN data",indent = "  ")
    print

    #
    # Determine the links between the chains reported for the resonances and the
    # actual chains present in the molSystem
    #

    if not self.setChainLinks(forceChainMappings = self.forceChainMappings):
      self.messageReporter.showWarning("Warning","No chain selected... aborting linkResonances",parent = self.guiParent)
      return False
    
    #
    # Print message...
    #
   
    print
    print drawBox("Stage 3: Determining atom naming system and creating mapping",indent = "  ")
    print

    #
    # Set up resonance groups for resonances without any atom name
    #
    
    if self.createResonanceGroups():
      return True

    #
    # Determine all atom naming systems in library
    # Assuming they are all present on chemComp level
    #
  
    if not self.matchNamingSystems():
      self.messageReporter.showWarning("Warning","No naming system selected... aborting linkResonances",parent = self.guiParent)
      return False
  
    #
    # Check if there is any atom name mapping available, set in application data if so
    #
    
    self.useAtomMapping()
    
    #
    # Use 'common' names if so desired - use with care as this slows down code
    # and can produce wrong results.
    #

    self.applyNameMappings()
       
    #
    # Link resonances to atoms/atom groups 
    #
    
    self.setupAtomToResName()
    
    self.checkAtomToResName()
    
    #
    # Print message...
    #
   
    print
    print drawBox("Stage 4: Assigning resonances (created from format data) to atoms (CCPN data)",indent = "  ")
    print

    #
    # Link resonances to atoms
    #
 
    self.linkResonancesToAtoms()
        
    # TODO: what should be set here?
    self.status = True
    
    return None
  
  def setKeywords(self, keywds, verbose = True):
      
    returnStatus = True
    
    for keyword in self.keywordInfo.keys():

      (default,mandatory,infoText) = self.keywordInfo[keyword]

      if mandatory and not keywds.has_key(keyword):
        
        print "  Error: mandatory linkResonances keyword '%s' missing! Cannot continue." % (keyword)
        returnStatus = False
        break
          
      if keywds.has_key(keyword):
        setattr(self,keyword,keywds[keyword])
        #print "%s fromkeywds:" % functionName,IOkeyword,keywds[IOkeyword]
        del keywds[keyword]

      else:
        setattr(self,keyword,default)
        #print "%s fromdefault:" % functionName,IOkeyword,default

    if keywds and verbose:
      for keywd in keywds.keys():
        print "  Warning: keyword '%s' with value '%s' not recognized." % (keywd,keywds[keywd])
      
    return returnStatus

  def getStatusOtherProchiral(self,atomSets,resName):

    """
    
    User interaction to determine status of other (missing) prochiral
    
    """

    refAtomSet = atomSets[0]
    refAtom = refAtomSet.findFirstAtom()
    otherAtomName = atomSets[1].name
    
    nmrRes = refAtom.residue
    seqCode = nmrRes.seqCode
    ccpCode = nmrRes.molResidue.ccpCode
 
    if self.setSingleProchiral == None:
    
      #
      # First check if mapping available
      #
      
      singleProchiralStatus = None
      resStatusMapping = self.nameMapping.findFirstResStatusMapping(resonanceName = resName)
      
      if self.useSinglePossEquivMapping:
      
      
        if resStatusMapping and resStatusMapping.isSingleProchiral != None:
      
          singleProchiralStatus = resStatusMapping.isSingleProchiral
          if self.verbose:
            print "  Set single prochiral status for '%s' to %d from mapping." % (resName,singleProchiralStatus)


      if singleProchiralStatus == None:

        #
        # User interaction: do the resonances carry the same info? Or was only one found?
        #
    
        interaction = self.multiDialog.SingleProchiralStatus(self.guiParent,refAtomSet,otherAtomName)
        singleProchiralStatus = interaction.singleResonanceStatus
        
        #
        # Set mapping
        #
        
        if not resStatusMapping:
        
          resStatusMapping = NameMapping.ResonanceStatusMapping(self.nameMapping,resonanceName = resName)
          
        resStatusMapping.isSingleProchiral = singleProchiralStatus
      
    else:
      
      singleProchiralStatus = self.setSingleProchiral
 
 
    if self.verbose:

      if singleProchiralStatus == 0:
 
        print "  Adding resonance for prochiral atom/group %s (from residue %s %d, resonance %s)..." % (otherAtomName,ccpCode,seqCode,self.resonance.name)
 
      elif singleProchiralStatus == 1:

        print "  Setting resonance %s (%s %d) to non-stereospecifically assigned..." % (self.resonance.name,ccpCode,seqCode)

    return singleProchiralStatus

  def getStatusResonance(self,refAtom,otherAtomName,resName):

    """
    
    User interaction to determine status of other (missing?) atom
    for (possible equivalent) aromatic/amide atom groups
    
    """
    
    nmrRes = refAtom.residue
    seqCode = nmrRes.seqCode
    ccpCode = nmrRes.molResidue.ccpCode

    if self.setSinglePossEquiv == None:
 
      #
      # First check if mapping available
      #
      
      singleResonanceStatus = None
      resStatusMapping = self.nameMapping.findFirstResStatusMapping(resonanceName = resName)
      
      if self.useSinglePossEquivMapping:
      
      
        if resStatusMapping and resStatusMapping.isSinglePossEquiv != None:
      
          singleResonanceStatus = resStatusMapping.isSinglePossEquiv
          
          if self.verbose:
            print "  Set single possible equivalent status for '%s' to %d from mapping." % (resName,singleResonanceStatus)

      if singleResonanceStatus == None:

        #
        # User interaction: do the resonances carry the same info? Or was only one found?
        #
    
        interaction = self.multiDialog.SingleResonanceStatus(self.guiParent,refAtom,otherAtomName)
        singleResonanceStatus = interaction.singleResonanceStatus
        
        #
        # Set mapping
        #
        
        if not resStatusMapping:
        
          resStatusMapping = NameMapping.ResonanceStatusMapping(self.nameMapping,resonanceName = resName)
          
        resStatusMapping.isSinglePossEquiv = singleResonanceStatus
      
    else:
      
      singleResonanceStatus = self.setSinglePossEquiv
 
    if self.verbose and singleResonanceStatus == 0:
 
      print "  Changing status to equivalent for aromatic/amide atom %s (from residue %s %d)..." % (otherAtomName,ccpCode,seqCode)
 
    return singleResonanceStatus


  def setNamingSystemName(self,namingSystemName):
  
    """
    Sets the selected naming system
    """
  
    self.namingSystemNameSelected = namingSystemName
    
  def checkNameMapping(self):
  
    """
    Checks if a nameMapping exists for this project
    """
  
    if self.nameMapping:
    
      self.format = self.nameMapping.dataFormat

      namingSystemName = self.nameMapping.namingSystem
      
      if namingSystemName == 'Unknown':
        pass

      else:
        self.namingSystemName = namingSystemName
  
  def checkAssignFormat(self):
  
    """
    Checks which original format the resonance names come from (e.g. nmrStar,
    XEasy, ...)
    """
      
    formatsFound = {}
    
    if not self.format:
      
      if self.nmrConstraintStore:
        resonances = self.nmrConstraintStore.fixedResonances
      elif self.nmrProject:
        resonances = self.nmrProject.resonances
        
      for resonance in resonances:
      
         #
         # Only consider non-linked ones
         #
         
         if not resonance.resonanceSet:
         
           assignData = resonance.findAllApplicationData(keyword = assign_kw)
           
           for assignDatum in assignData:
           
             format = assignDatum.application
           
             if formatsFound.has_key(format):
               formatsFound[format] += 1
              
             else:
               formatsFound[format] = 1
      
      #
      # Find the best matching format
      #               
      
      formatNamesList = formatsFound.keys()
      formatNamesList.sort()
      
      for format in formatNamesList:
  
        if formatsFound[format] == len(resonances):
          if not self.format:
            if self.verbose:
              print "  Selecting format %s for atom name translation." % format
            self.format = format
          elif self.verbose:
            # TODO: need Multiple selection really but ignored for now
            print "  WARNING: format %s also fully matches!" % format
      
      if not self.format:
        matchValues = formatsFound.values()
        matchValues.sort()
        
        for format in formatNamesList:
          if formatsFound[format] == matchValues[-1]:
            if not self.format:
              if self.verbose:
                print "  Selecting format %s for atom name translation - has highest match." % format
              self.format = format
            elif self.verbose:
              print "  Warning: format %s also has highest match!" % format
      
  def setupRefData(self):
 
    """
    Set up reference data for the (fixed) resonances.
    
      - if they are linked to atoms already:
           - keep track of them (by nmrAtom reference!!!)
    
      - if they're not linked to atoms:
           - keep track of the chain(s) that occur for them
           - keep track of them (dictionary with format name as key and resonance pointer as value)
    """
 
    if not self.nmrConstraintStore:
      resonanceList = self.nmrProject.sortedResonances()
      self.resonanceClass = Nmr.Resonance
      self.resonanceSetClass = Nmr.ResonanceSet
      self.resonanceParent = self.nmrProject
      
      self.atomSetClass = Nmr.AtomSet
      
      atomSetsList = self.nmrProject.sortedAtomSets()

    else:
      resonanceList = self.nmrConstraintStore.sortedFixedResonances()
      self.resonanceClass = NmrConstraint.FixedResonance
      self.resonanceSetClass = NmrConstraint.FixedResonanceSet
      self.resonanceParent = self.nmrConstraintStore

      self.atomSetClass = NmrConstraint.FixedAtomSet
      
      atomSetsList = self.nmrConstraintStore.sortedFixedAtomSets()
     
    self.atomSetsAtomDict = {}
    for atomSet in atomSetsList:
      atoms = list(atomSet.atoms)
      atoms.sort()
      self.atomSetsAtomDict[atomSet] = atoms

    #
    # Now set up... for either normal or fixed resonances
    #
    
    for resonance in resonanceList:
 
      #
      # Check if the resonance is linked to an atom already
      #
 
      resonanceLinked = 0
      resonanceSet = resonance.resonanceSet
 
      if resonanceSet:
 
        for atomSet in resonanceSet.atomSets:
 
          resonanceLinked = 1
          
          for atom in atomSet.atoms:
            
 
            if not self.atomResonances.has_key(atom):
              self.atomResonances[atom] = []
 
            self.atomResonances[atom].append(resonance)
 
      if not resonanceLinked:

        #
        # If there is no resonance-atom link, set the unassigned resonance names dictionary
        #

        resNames = updateResonanceNamesDict(resonance,self.resUnassignedNames,self.format)

        for resName in resNames[:]:
        
          #
          # Convert name if necessary
          #
        
          if self.specificResNameMappings:
            if resName in self.specificResNameMappings:
              
              newResName = self.specificResNameMappings[resName]
              
              if type(newResName) == type(''):
                newResNames = [newResName]
              else:
                newResNames = newResName
              
              self.resUnassignedNames[resName].pop()
              resNameIndex = resNames.index(resName)
              resNames.pop(resNameIndex)
              
              for newResName in newResNames:

                if not self.resUnassignedNames.has_key(newResName):
                  self.resUnassignedNames[newResName] = []
    
                self.resUnassignedNames[newResName].append(resonance)
                
                if self.verbose:
                  print "  Mapped specific name '%s' to '%s'..." % (resName,newResName)
                
                resNames.insert(resNameIndex,newResName)
        
        #
        # Now use remapped set of resNames to set dictionaries
        #
        
        for resName in resNames:
          
          # 
          # Check if the name can be decomposed or not (no use if not!)
          #
          
          (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resName,verbose = 0)
          
          if seqCode != None:
            # Add to list of unassigned resonance names if not already there
            # Note that this list is just for faster looping: links to the
            # actual resonances are in the resUnassignedNames dict
            if self.resUnassignedNameList.count(resName) == 0:
              self.resUnassignedNameList.append(resName)
 
            # Keep track of reported chainCodes for linking later on
            if not self.formatChains.has_key(chainCode):
              self.formatChains[chainCode] = []
              
            seqCodeKey = (seqCode,seqInsertCode)
 
            # Also keep track of reported seqCodes and seqInsertCodes
            if seqCodeKey not in self.formatChains[chainCode]:
              self.formatChains[chainCode].append(seqCodeKey)
              
            if not atomName:
              self.missingAtomNames = 1

          else:

            #
            # Don't list if has other application data name(s)
            #
            
            if resName or not resonance.findAllApplicationData(keyword = assign_kw):

              # TODO: could have interaction here to select relevant seqCode and
              #       atomName...
              
              if self.verbose:
                print "  Error: Can't handle resonance name '%s'... ignored" % resName
              
            # Delete name again... bit wasteful but unlikely to happen anyway
            self.resUnassignedNames[resName].remove(resonance)
 
            if not self.resUnassignedNames[resName]:
              del self.resUnassignedNames[resName]
 
      else:

        #
        # If there is a link, set the assigned resonance names dictionary
        #

        updateResonanceNamesDict(resonance,self.resAssignedNames,self.format)
 
    self.resUnassignedNameList.sort()
    
    #
    # Also sort seqCodes for resChains
    #
    
    for formatChainCode in self.formatChains.keys():
      self.formatChains[formatChainCode].sort()

  def createResonanceGroups(self):

    """
    If there are valid chainCode/seqCode combinations without atom names, create a resonance group
    """

    if self.missingAtomNames and not self.nmrConstraintStore:
      
      if self.verbose:
        print 'Creating resonance groups for resonances with no atom names...'
        
      skippedResonances = 0
    
      for resName in self.resUnassignedNameList[:]:

        if not self.updateResonance(resName, verbose = 0):
          skippedResonances += 1
          continue

        if not self.atomName and self.nmrRes:
          
          resonanceGroup = self.getResonanceGroup(self.nmrRes)
          resonanceGroupResn = resonanceGroup.resonances
          
          for resonance in self.resUnassignedNames[self.resName]:
            if resonance not in resonanceGroupResn:
              resonanceGroup.addResonance(resonance)
          
          del self.resUnassignedNames[self.resName]
          self.resUnassignedNameList.pop(self.resUnassignedNameList.index(resName))
      
      if len(self.resUnassignedNameList) == skippedResonances:
        self.messageReporter.showWarning("Warning","LinkResonance will not finish because no resonances were found with associated atom names. All resonances were linked to the relevant ResonanceGroup, but not to any atoms. This problem probably occurs because the original imported file does not contain any atom names.",parent = self.guiParent) 
        return True
        
      elif self.verbose:
        self.messageReporter.showWarning("Warning","Some resonances do not have associated atom names. These resonances were linked to the relevant ResonanceGroup.",parent = self.guiParent)
      
    return False

  def getResonanceGroup(self,nmrRes):

    resonanceGroup = self.nmrProject.findFirstResonanceGroup(residue = nmrRes)

    if not resonanceGroup:

      resonanceGroup = Nmr.ResonanceGroup(self.nmrProject,chains = [nmrRes.chain], residue = nmrRes,details = 'Created by formatConverter')
    
    return resonanceGroup
    
  def matchNamingSystems(self):
  
    """
    Finds out which naming system(s) best match the current list of
    resonance names
    """

    if self.namingSystemName:
    
      #
      # Naming system passed in...
      #

      self.setNamingSystemName(self.namingSystemName)
        
    else:
    
      #
      # Set up list of naming systems (use protein chemComps as default for now)
      # Uses all chemComps defined in project... need to access them later
      # on anyway so might as well do it here.
      #
      # Note that added a CUSTOM naming system for automatically created ChemComps
      #

      resNameFound = {'CUSTOM': 0}

      refChemComps = (self.project.findFirstChemComp(molType = 'protein', ccpCode = 'Ala'),
                      self.project.findFirstChemComp(molType = 'DNA', ccpCode = 'T'))
      
      for refChemComp in refChemComps:
        for namingSystem in refChemComp.namingSystems:
          namingSystemName = namingSystem.name
          if not resNameFound.has_key(namingSystemName):
            resNameFound[namingSystemName] = 0

      #
      # Suggest atom naming (based on most matches using CCPN reference lib)
      #

      assigned = 0
      noResidue = 0
      chemCompAtomSysNames = {}

      for resName in self.resUnassignedNameList:

        assigned += self.updateResonance(resName, verbose = self.verbose)

        if self.nmrRes == None:
          noResidue += 1
          continue

        atomName = self.atomName

        #
        # Sysnames now accessible from namingSystem
        #
        
        if not chemCompAtomSysNames.has_key(self.chemComp):
          chemCompAtomSysNames[self.chemComp] = {}
          for namingSystem in self.chemComp.namingSystems:
            
            chemCompAtomSysNames[self.chemComp][namingSystem] = {}
            
            for chemAtom in self.chemComp.chemAtoms:
              atomSysName = findChemAtomOrSetSysName(namingSystem,chemAtom)
              if atomSysName:
                chemCompAtomSysNames[self.chemComp][namingSystem][atomSysName] = (atomSysName.sysName,) + atomSysName.altSysNames

            for chemAtomSet in self.chemComp.chemAtomSets:
              atomSysName = findChemAtomOrSetSysName(namingSystem,chemAtomSet)
              if atomSysName:
                chemCompAtomSysNames[self.chemComp][namingSystem][atomSysName] = (atomSysName.sysName,) + atomSysName.altSysNames

        self.getAtomSysNamesNum(chemCompAtomSysNames[self.chemComp],resNameFound,atomName)

      #
      # GUI interaction: pick naming system. Will have automatically selected highest % one.
      # Don't bother if no connections to any residues
      #

      if noResidue != len(self.resUnassignedNameList):

        if not self.minimalPrompts and self.allowPopups:
          interaction = self.multiDialog.NamingSystem(self.guiParent,resNameFound,assigned)
          namingSystemName = interaction.namingSysName  # namingSystem should be one word
          interaction.destroy()

        else:

          #
          # Pick highest scoring ones
          #
          
          #print resNameFound

          bestNamingSystemNames = [resNameFound.keys()[0]]

          for namingSystemName in resNameFound:

            if resNameFound[namingSystemName] > resNameFound[bestNamingSystemNames[0]]:

              bestNamingSystemNames = [namingSystemName]
          
            elif resNameFound[namingSystemName] == resNameFound[bestNamingSystemNames[0]]:
              if namingSystemName not in bestNamingSystemNames:
                bestNamingSystemNames.append(namingSystemName)

          bestNamingSystemNames.sort()
          
          #
          # Follow hierarchy of naming systems to pick if multiple matches...
          #
          
          namingSystemName = None
          
          if len(bestNamingSystemNames) > 1:
            for namingSystemOrder in (('IUPAC','PDB','GROMOS'),('AQUA','CYANA2.1','PDB_REMED','XPLOR','MSD','DISGEO','SYBYL')):
              for tempNamingSystemName in namingSystemOrder:
                if tempNamingSystemName in bestNamingSystemNames:
                  namingSystemName = tempNamingSystemName
                  break

              if namingSystemName:
                break
          
          if not namingSystemName:
            namingSystemName = bestNamingSystemNames[0]

      else:
        if self.verbose:
          print "  No valid resonances connected to residue... ignoring naming system selection"
        namingSystemName = 'IUPAC'


      if namingSystemName:
        self.setNamingSystemName(namingSystemName)
        self.namingSystemName = namingSystemName # TODO: is this correct?!

        if self.verbose:
          print "  Naming system %s selected" % namingSystemName
          print

      else:
        self.status = False

    #
    # Do some final checking and resetting
    #
  
    if self.useIupacMatching and self.namingSystemName == 'IUPAC':
      self.useIupacMatching = False

    return self.status      


  def getAtomSet(self,nmrAtom):
  
    #
    # Gets the relevant AtomSet or FixedAtomSet
    #
  
    if self.resonanceClass == Nmr.Resonance:
      atomSet = nmrAtom.atomSet
    else:
      atomSet = nmrAtom.findFirstFixedAtomSet(nmrConstraintStore = self.resonanceParent)
    
    return atomSet 
      
  def getAtomSysNamesNum(self,chemCompAtomSysNameDict,resNameFound,atomName):
    
    """
    Look for all atomSysNames from a chemComp
    """
  
    found = 0
    namingSystemNamesMatched = []
   
    for namingSystem in chemCompAtomSysNameDict.keys():
      namingSystemName = namingSystem.name
      for atomSysName in chemCompAtomSysNameDict[namingSystem].keys():
        if atomName in chemCompAtomSysNameDict[namingSystem][atomSysName]:
          if namingSystemName not in namingSystemNamesMatched:
            resNameFound[namingSystemName] += 1
            found += 1
            namingSystemNamesMatched.append(namingSystemName)

    return found

  def setProchiralGroup(self,atomSets, otherResonance = None):
  
    """
    Code to set the resonance-atom link for a prochiral group
    """
    
    #
    # Initialize
    #
    
    resonance = self.resonance
    proj = self.project
    
    
    if not otherResonance:

      #
      # Create a new resonance for prochiral to current resonance
      #

      otherResonance = self.resonanceClass(self.resonanceParent,isotopeCode = resonance.isotopeCode)
      copyResonanceInfo(resonance,otherResonance,nmrShiftSelect = self.nmrShiftSelect,copyTrack = True,verbose = self.verboseCopy)
 
    curResonances = [resonance,otherResonance]
    curResonanceNames = []

    # Stereoassign = 1 assumes all stereospecific assigned...
    # link each resonance straight to an atomSet
    if self.stereoAssign == 1:
 
      for i in range(0,2):
        atomSet = atomSets[i]
        curResonance = curResonances[i]
 
        self.resonanceSetClass(self.resonanceParent,resonances = [curResonance], atomSets = [atomSet])
        curResonanceNames.append(atomSet.name)
 
    # If not, then link them with one resonanceSet to the atomSets
    else:
      self.resonanceSetClass(self.resonanceParent,resonances = curResonances, atomSets = atomSets)
      curResonanceNames.append(createNonStereoName(atomSets[0].name,'a'))
      curResonanceNames.append(createNonStereoName(atomSets[0].name,'b'))
 
    # Set resonance names
    for i in range(0,2):
      atomSet = atomSets[i]
      curResonance = curResonances[i]
      curResonanceName = curResonanceNames[i]
 
      setCcpNmrResonanceName(curResonances[i],atomSet.findFirstAtom(),atomName = curResonanceName)
 
    return otherResonance

  def updateResonance(self,resName, verbose = 1):
 
    """
    This code makes sure all self. variables required by the code are set
    for a resName (and that the resName itself is checked). Is reused in
    several places where code loops over all resNames.
    
    SeqCode is int, atomName is uppercased by this function
    """
    
    nmrResMatches = 1
    
    (chainCode,seqCode,spinSystemId,seqInsertCode,self.atomName) = getNameInfo(resName)

    self.origAtomName = self.atomName
 
    if chainCode != self.chainCode:
    
      if not self.formatChainCodeDict.has_key(chainCode):
        if verbose:
          print "  Warning: ignoring file resonance '%s'... chainCode not mapped" % resName
        return 0

      else:
      
        self.chainMappings = self.formatChainCodeDict[chainCode]
        self.chainCode = chainCode
    
        # Reset seqCode so new nmrRes is assigned in case of multiple chains!
        self.seqCode = ""
        self.seqInsertCode = None

    if seqCode != self.seqCode or seqInsertCode != self.seqInsertCode:

      self.seqCode = seqCode
      self.seqInsertCode = seqInsertCode
      nmrResMatches = 0

      #
      # Loop over all data model chains linked to this file chain:
      # take one that has matching seqCode and/or seqId range
      #
      
      (seqCodeToIdMapping,seqCodeToChainMapping,seqCodes) = self.getChainMappingsSeqCodeInfo(self.chainMappings)
      
      #
      # The seqCode has to fall within the expected lower limit. 
      #
      
      offset = None

      for (firstSeqCodeKey,lastSeqCodeKey) in seqCodes:
        (firstSeqCode,seqInsertCode) = firstSeqCodeKey
        (lastSeqCode,seqInsertCode) = lastSeqCodeKey
        
        if self.seqCode >= firstSeqCode and self.seqInsertCode == seqInsertCode and (lastSeqCode == None or self.seqCode <= lastSeqCode):
          offset = seqCodeToIdMapping[(firstSeqCode,seqInsertCode)] - firstSeqCode
          break
      
      if offset != None:
        
        #
        # If seqId bigger than length of chain it will be automatically ignored.
        # This is already mostly taken care of by new formatLastSeqCode attribute for ChainMapping (wim 06/04/06)
        #
        
        seqId = seqCode + offset

        #
        # The above two statements should ensure that only one nmrRes is selected
        # when one formatChainCode maps onto two (or more) ccp chains. But what in
        # the case of a homodimer with a single list of values that corresponds to
        # one chain but with two ccp chains available? Have to copy all resonances... 
        # is a mess! Use one as a 'reference' chain?!?! TODO TODO but ignore for now
        #
        
        chainMapping = seqCodeToChainMapping[(firstSeqCode,seqInsertCode)]
        chain = self.chainMappingChain[chainMapping]
        nmrRes = chain.findFirstResidue(seqId = seqId)

        if nmrRes:
          if nmrRes not in self.ignoreResidues:
            self.chain = chain
            self.nmrRes = nmrRes
            nmrResMatches += 1
          else:
            self.nmrRes = None
            if verbose:
              resNameText = resName
              if not resNameText:
                resNameText = ','.join([getResNameText(res) for res in self.resUnassignedNames[self.resName]])
            
              print "  Error: ignoring file resonance '%s'... it is linked to a residue that is listed as not valid." % resNameText
            
            return 0
            

    if not nmrResMatches:
      # No match for seqCode within chain, but could be linked to resGroup
      # In any case, can't link this one (or have to do user interaction
      # to link seqCode to nmrResidue...)
      if verbose:
        resNameText = resName
        if not resNameText:
          resNameText = ','.join([getResNameText(res) for res in self.resUnassignedNames[self.resName]])
      
        print "  Error: ignoring file resonance '%s'... no matching residue in data model" % resNameText
      
      self.seqCode = ''
      self.seqInsertCode = None
      self.nmrRes = None
      
      return 0
      
    elif nmrResMatches > 1:
      
      #
      # Also illegal - should only get one match for each.
      #
      
      if verbose:
        print "  Warning: resonance '%s' gives multiple residue matches - ignoring." % resName

      self.seqCode = ''
      self.seqInsertCode = None
      self.nmrRes = None

      return 0

    elif self.nmrRes:
      
      self.resName = resName
      self.stereoAssign = self.globalStereoAssign

      self.chemCompVar = self.nmrRes.chemCompVar
      self.chemComp = self.chemCompVar.chemComp
      self.ccpCode = self.chemComp.ccpCode
      
      # Set local use IUPAC if the names were reset!
      self.localUseIupacMatching = self.useIupacMatching
      for resonance in self.resUnassignedNames[self.resName]:
        applData = resonance.findFirstApplicationData(application = self.format, keyword = 'resetToIupac')
        if applData:
          self.localUseIupacMatching = 2

      return 1
  
    else:
      
      return 0

  def checkStereoStatus(self):
  
    """
    Set the stereospecific assignment status for a resName
    """
  
    if not self.globalStereoAssign:
    
      stereoAssign = None

      #
      # Reset self.stereoAssign based on mapping
      #
    
      resStatusMapping = self.nameMapping.findFirstResStatusMapping(resonanceName = self.resName)
    
      if resStatusMapping and resStatusMapping.isStereo != None:
    
        stereoAssign = resStatusMapping.isStereo
    
      #
      # If no mapping, ONLY reset in case there is ONE resonance
      # and the stereoCode is 1!!
      #
    
      if stereoAssign == None and len(self.resUnassignedNames[self.resName]) == 1:
        resonance = self.resUnassignedNames[self.resName][0]
        stereoCode = getStereoCode(resonance,self.format)
      
        if stereoCode == 1:
          stereoAssign = 1
        
          if not resStatusMapping:

            resStatusMapping = NameMapping.ResonanceStatusMapping(self.nameMapping,resonanceName = self.resName,isStereo = True)

      if stereoAssign != None:
      
        self.stereoAssign = stereoAssign


  def applyNameMappings(self):
     
    """
    Applies direct name mappings (naming system independent)
    """

    if self.useCommonNames or self.useAmbiguity or self.addNameMappings:
    
      tempResUnassignedNameList = self.resUnassignedNameList[:]
       
      for resName in tempResUnassignedNameList:
  
        if not self.updateResonance(resName,verbose = 0):
        
          continue
      
        else:
                   
          if self.useCommonNames:
          
            commonNameMappings = self.chemComp.findAllApplicationData(application = ccpNmr_kw, keyword = 'commonNameMapping')
            self.setApplDataNameMappings(commonNameMappings,self.useCommonNames)
          
          if self.useAmbiguity:
          
            ambiguityMappings = self.chemComp.findAllApplicationData(application = ccpNmr_kw, keyword = 'ambiguityMapping')
            self.setApplDataNameMappings(ambiguityMappings,self.useAmbiguity)
            
          if self.addNameMappings:
          
            ccpCode = self.chemComp.ccpCode
          
            if self.addNameMappings.has_key(ccpCode) or 'ALL' in self.addNameMappings.keys():
            
               nameMappingList = []
               
               if self.addNameMappings.has_key(ccpCode):
                 nameMappingList += self.addNameMappings[ccpCode]
              
               if 'ALL' in self.addNameMappings.keys():
                 nameMappingList += self.addNameMappings['ALL']
            
               for nameMapping in nameMappingList:

                 self.setNameMappings(nameMapping,1)
        
  def setApplDataNameMappings(self,mappings,setStatus):
  
    """
    Sets direct name mappings from application data
    """

    if mappings:
    
      for nameMapping in mappings:

        nameMappings = eval(nameMapping.value)
        
        self.setNameMappings(nameMappings,setStatus)
    
  def setNameMappings(self,nameMappings,setStatus):
  
    """
    Sets the direct name mappings
    """

    if nameMappings[0] == self.atomName:

      #
      # Match. The original resName WILL BE RESET for linking purposes!!
      #
      # If there are multiple mapped names they will also be split up at this stage.
      # It is assumed that the names 'match' each other for whatever namingSystem
      # is in use... the mapping is dealt with at the next stage, and can be reused
      # as long as useCommonNames or useAmbiguity is on.
      #
      
      if self.resUnassignedNameList.count(self.resName) > 0:

        resUnassignedNameIndex = self.resUnassignedNameList.index(self.resName)
        refResonances = self.resUnassignedNames[self.resName]

        #
        # Delete original resName...
        #

        del(self.resUnassignedNames[self.resName])
        self.resUnassignedNameList.pop(resUnassignedNameIndex)

      else:
      
        #print "ERROR: resonance %s not present in unassigned list when setting name mapping..." % self.resName
        return

      #
      # Set new resName info only if setStatus == 1...
      #

      if setStatus == 1:

        for i in range(1,len(nameMappings)):

          newAtomName = nameMappings[i]
          newResName = getResName(self.chainCode,self.seqCode,newAtomName,seqInsertCode = self.seqInsertCode)

          #
          # Set which resonances to use
          #

          if i == 1:
            newResonances = refResonances
            actionText = 'resetting'

          else:
            newResonances = []
            actionText = 'copying'

            for refResonance in refResonances:

              newResonances.append(self.resonanceClass(self.resonanceParent,isotopeCode = 'unknown'))
              copyResonanceInfo(refResonance,newResonances[-1],nmrShiftSelect = self.nmrShiftSelect,verbose = self.verboseCopy)
          
          #
          # Flag the resonance - its name was reset to IUPAC!
          #
          
          for newResonance in newResonances:
            newResonance.addApplicationData(Implementation.AppDataBoolean(application = self.format, keyword = 'resetToIupac', value = 1))

          #
          # Put info in resUnassignedNames dict and list
          #

          if self.resUnassignedNames.has_key(newResName):
            self.resUnassignedNames[newResName].extend(newResonances)
            resNameText = 'existing'

          else:
            self.resUnassignedNames[newResName] = newResonances
            self.resUnassignedNameList.insert(resUnassignedNameIndex,newResName)
            resNameText = 'new'

          if self.verbose:
            print "  Warning: using name %s for %s - %s to %s resName %s." % (newAtomName,self.resName,actionText,resNameText,newResName)

      else:

        if self.verbose:
          print "  Warning: not linking resonance '%s'." % (self.resName)

        #
        # TODO: add error message for nmrStar? Or just leave... will show up as not linked anyway.
        #

  def setupAtomToResName(self):
  
    """
    Prepares information to set up the AtomToResName objects
    """
        
    for resName in self.resUnassignedNameList:
  
      if not self.updateResonance(resName,verbose = 0):
        continue
                    
      #
      # Make sure all resonances with atom names are connected to a resonanceGroup...
      #
      
      if not self.nmrConstraintStore:

        resonanceGroup = self.getResonanceGroup(self.nmrRes)
        resonanceGroupResn = resonanceGroup.resonances
      
        for resonance in self.resUnassignedNames[self.resName]:
          if resonance not in resonanceGroupResn:
            resonanceGroup.addResonance(resonance)
            
      #
      # If it's a residue type constraint, set it here and continue
      #
      
      if self.atomName == allResidueAtoms_kw:
        self.setResidueAtomToResName()
        continue
      
      #
      # Here reset protonation state of a chemCompVar, but ONLY if there is a clear atom name match!
      #
      
      if self.checkProtonationStates.has_key(self.ccpCode):
        atomNames = self.checkProtonationStates[self.ccpCode]
        if self.atomName in atomNames:
          atom = self.nmrRes.findFirstAtom(name = self.atomName)
          if not atom:
            origDesc = self.nmrRes.descriptor
            descDict = getDescriptorDict(origDesc)
            mergeDescriptorDicts(descDict,{'prot': [self.atomName]})
            newDesc = getDescriptorFromDict(descDict)
            self.nmrRes.descriptor = newDesc
            self.nmrRes.newAtom(name = self.atomName)
            self.chemCompVar = self.nmrRes.chemCompVar
            
            print "  Warning: reset protonation state from %s to %s for %s (%s,%d) based on format atom name %s!" % (origDesc,newDesc,self.ccpCode,self.chain.code,self.nmrRes.seqId,self.atomName)
      
      #
      # Check whether the current naming system has any names for this particular chemCompVar
      #
      
      tempNamingSystem = None
      
      if not self.ccNamingSystems.has_key(self.chemComp):
        
        namingSystem = self.chemComp.findFirstNamingSystem(name = self.namingSystemName) 
               
        if namingSystem and self.hasAtomSysNames(namingSystem):
          self.ccNamingSystems[self.chemComp] = namingSystem
        
        else:
          # Need something... look for IUPAC, MSD, PDB in that order...
          for tempNamingSystemName in self.tempNamingSystemNames:
            
            tempNamingSystem = self.chemComp.findFirstNamingSystem(name = tempNamingSystemName)
            if tempNamingSystem and self.hasAtomSysNames(tempNamingSystem):
              self.ccNamingSystems[self.chemComp] = tempNamingSystem
              print "  Warning: using local naming system %s for %s because no %s atom names available!" % (tempNamingSystemName,self.chemComp.ccpCode,self.namingSystemName)
              break
       
      elif self.ccNamingSystems[self.chemComp].name != self.namingSystemName:
        tempNamingSystem = self.ccNamingSystems[self.chemComp]          

      #
      # Check if there is any stereo information available...
      #
      
      self.checkStereoStatus()

      #
      # Set unknownEquivalent - this is to keep track of e.g. HE1
      # atoms for Phe,Tyr
      #

      self.unknownEquivalent = 0

      #
      # Check if there's any atom name mapping available for this resonance (name)
      # Also set the relevant namingSystem
      #
      
      if not self.getAtomNameMapping():
        continue
      
      #
      # Check chemAtomSysNames (use tempNamingSystem if available)
      #
      
      if tempNamingSystem:
        useNamingSystemName = tempNamingSystem.name
      else:
        useNamingSystemName = self.namingSystemName
        
      chemAtomSysNames = self.getAtomSysNames(self.atomName,useNamingSystemName,self.chemCompVar)
      
      #
      # If useIupacMatching is on (and no match on previous),
      # then try to find a match using the IUPAC naming...
      # This is ONLY done for atoms that are NOT part of an atomSet
      # otherwise problems can occur with HB2/HB3!
      #
      
      if self.localUseIupacMatching and not chemAtomSysNames:
        chemAtomSysNames = self.getIupacChemAtomSysNames()

      #
      # Fix problems with single atom names given as atom set name (e.g. HG# for LEU HG)
      #
      
      if self.fixSingleAtoms and not chemAtomSysNames:

        if self.atomName and self.atomName[-1] in '*#%+':

          tryAtomName = self.atomName[:-1]

          chemAtomSysNames = self.getAtomSysNames(tryAtomName,self.namingSystemName,self.chemCompVar)

          if chemAtomSysNames:            
            print "  Warning: reset atom name %s to %s for residue %s - is single atom." % (self.atomName,tryAtomName,self.resName)
      
      """
      #
      # Finally, if nothing found, see if the namingSystem has an atomSetReference. If not, try
      # to use one of the tempNamingSystemNames ones.
      #
      
      if not chemAtomSysNames:
        # TODO: silly that always re-searching for this
        namingSystem = self.chemComp.findFirstNamingSystem(name = useNamingSystemName)
        if not namingSystem.atomSetReference:
          for tempNamingSystemName in self.tempNamingSystemNames:
            chemAtomSysNames = self.getAtomSysNames(self.atomName, tempNamingSystemName, self.chemCompVar)
            
            if chemAtomSysNames:
              print "  Warning: using naming system %s to set atom name %s" % (tempNamingSystemName,self.atomName)
              break
      """
      
      #
      # If match found, look for nmrAtoms automatically
      #
      
      foundChemAtom = False
      
      if chemAtomSysNames:
      
        #
        # Get the matching chemAtoms
        # Also add missing atoms and chemAtomSets (e.g. HD12 reported for Leu)
        #
        
        foundChemAtom = self.getMatchingChemAtoms(chemAtomSysNames)
      
      elif self.useCcpnAtomMatching:
           
        chemAtom = self.chemCompVar.findFirstChemAtom(name = self.atomName)
        
        # TODO this might not be completely correct! Deep chemAtomSets?
        if chemAtom:
          foundChemAtom = True
          self.chemAtoms = [chemAtom]
          self.chemAtomSets = []
        else:
          chemAtomSet = self.chemCompVar.findFirstChemAtomSet(name = self.atomName)
          if chemAtomSet:
            foundChemAtom = True
            self.chemAtomSets = [chemAtomSet]
            self.chemAtoms = chemAtomSet.sortedChemAtoms()
      
      if not chemAtomSysNames and not foundChemAtom:
        print "  Error: No matches at all for resonance name %s (%s)" % (resName,self.resName),
        if self.nmrRes:
          print " - ccpCode %s" % self.nmrRes.ccpCode
        else:
          print
        
      #
      # If no chemAtomSysName or following chemAtom match found,
      # then let user link to relevant atom(s)
      #
      
      if not foundChemAtom:
        
        #
        # Present list of atom names.
        #
        # Also ask user to which level the mapping has to be propagated.
        #
        
        if not self.useLinkResonancePopup or not self.allowPopups or not self.getAtomMatch():
          continue
        
        #
        # Propagate the mapping. If 'Do not link' was selected this
        # function will return None and the loop is continued here.
        #
        
        if not self.setAtomNameMapping():
          continue
                
    
      #
      # Get the actual nmrAtoms corresponding to the chemAtoms that were found
      # 
      
      if not self.getMatchingAtoms():
        continue
        
      #
      # Now have all nmrAtoms, and chemAtomSets (if necessary). Now put info
      # into AtomToResName class.
      #
      
      self.setAtomToResNameList()

  def setAtomToResNameList(self):
  
    """
    Sets up AtomToResName object (or adds information to existing one)
    """
      
    nmrAtoms = []
    chemAtomDict = {}
    keywds = {}
    
    for (nmrAtom,chemAtom) in self.atomsFound:
      nmrAtoms.append(nmrAtom)
      chemAtomDict[nmrAtom] = chemAtom
    
    #
    # Keep track of single but possibly equivalent atoms (e.g. HE1 for F)
    #
    
    if self.unknownEquivalent:
    
      keywds['singleNmrAtom'] = nmrAtoms[0]
    
    nmrAtoms.sort()
    newAtomToResName = 1

    for atomToResName in self.atomToResNameList:
      
      if nmrAtoms == atomToResName.nmrAtoms:
      
        atomToResName.addResName(self.resName,self.chemAtomSets,nmrAtoms,**keywds)
        newAtomToResName = 0
        break
        
    if newAtomToResName:
      isotopeCode = getNmrIsotope(self.project,chemAtom.elementSymbol)
      self.atomToResNameList.append(AtomToResName(nmrAtoms,chemAtomDict,self.resName,self.chemAtomSets,isotopeCode,verbose = self.verbose,**keywds))

  def setResidueAtomToResName(self):
  
    # This will only occur once - things have already been sorted!
    self.atomToResNameList.append(ResidueAtomToResName(self.nmrRes.atoms,self.resName))

  def checkAtomToResName(self):
  
    """
    Checks whether subsets of nmrAtoms exist in AtomToResName instances (slow code)
    """
    
    for i in range(0,len(self.atomToResNameList)-1):
      
      firstSet = self.atomToResNameList[i]
            
      if not firstSet or isinstance(firstSet,ResidueAtomToResName):
        continue

      firstSetNmrAtoms = firstSet.nmrAtoms
      firstSetNumNmrAtoms = len(firstSetNmrAtoms)
      
      for j in range(i+1,len(self.atomToResNameList)):
      
        secondSet = self.atomToResNameList[j]
        
        if not secondSet or isinstance(secondSet,ResidueAtomToResName):
          continue
      
        secondSetNmrAtoms = secondSet.nmrAtoms
        
        #
        # Don't bother if same length: has to be same AtomToResName instance if same
        #
        
        comparisonStatus = -1

        if len(secondSetNmrAtoms) > firstSetNumNmrAtoms:
      
          #
          # First set subset of second?
          #
          
          subSet = firstSet
          fullSet = secondSet
          
          comparisonStatus = self.compareNmrAtoms(firstSetNmrAtoms,secondSetNmrAtoms)
        
        elif len(secondSetNmrAtoms) < firstSetNumNmrAtoms:

          #
          # Second set subset of first?
          #
          
          subSet = secondSet
          fullSet = firstSet
          
          comparisonStatus = self.compareNmrAtoms(secondSetNmrAtoms,firstSetNmrAtoms)
          
          
        if comparisonStatus == 0:
        
          if self.verbose:
            print "  Error: Partial subset detected. Is illegal!"
          
        elif comparisonStatus == 1:
          
          if self.verbose:
            print "  Warning: Full subset detected for %s. Added to full set %s." % (str(subSet.resNameList),
                                                                                     str(fullSet.resNameList))

          #
          # This situation is principle possible for prochirals (some stereo assigned, some not)
          # and for chemAtomSets with isEquivalent set to None (e.g. HE1 reported for Phe,
          # but also HE*).
          #
          # Best way seems to just add the subset to the full one: safest handling. 
          #
                    
          for resName in subSet.resNameList:
                        
            #
            # Have to set all chemAtomSets for this one (if stereo becomes non
            # stereo, if single becomes multi!)
            #
            
            chemAtomSets = [subSet.chemAtomDict[subSet.nmrAtoms[0]].chemAtomSet]
                              
            if chemAtomSets[0].isProchiral:
              if chemAtomSets[0].chemAtomSet:
                chemAtomSetList = list(chemAtomSets[0].chemAtomSet.chemAtomSets)
                otherChemAtomSet = chemAtomSetList[not chemAtomSetList.index(chemAtomSets[0])]
                
                chemAtomSets = list(chemAtomSets[0].chemAtomSet.chemAtomSets)
                  
              else:
                chemAtomSets = []
            
            #
            # Keep track of unknown equivalents
            #
            
            keywds = {}
            
            if subSet.singleAtoms.has_key(resName):
              keywds['singleNmrAtom'] = subSet.singleAtoms[resName]
              chemAtomSets = []
            elif len(subSet.nmrAtoms) == 1:
              keywds['singleNmrAtom'] = subSet.nmrAtoms[0]
              chemAtomSets = []

            fullSet.addResName(resName,chemAtomSets,subSet.nmrAtoms,**keywds)

          self.atomToResNameList[self.atomToResNameList.index(subSet)] = None
          
          if subSet == firstSet:
            
            #
            # Have to end this loop and continue with the main one
            #
            
            break
    
    #
    # Clean up
    #
    
    while(self.atomToResNameList.count(None) > 0):
      
      self.atomToResNameList.pop(self.atomToResNameList.index(None))

    #
    # If there are multiple resNames for the atomToResName, then
    # check whether some are sysName subsets of others... merge
    # if that is the case.
    #

    for atomToResName in self.atomToResNameList:
      
      namingSystem = atomToResName.nmrAtoms[0].chemAtom.chemComp.findFirstNamingSystem(name = self.namingSystemName)
      
      if len(atomToResName.resNameList) > 1:
        
        if atomToResName.hasEqChemAtomSets:
        
          for nmrAtoms in atomToResName.nmrAtomGroups:
          
            chemAtomSet = atomToResName.chemAtomDict[nmrAtoms[0]].chemAtomSet            
            chemAtomSetSysName = findChemAtomOrSetSysName(namingSystem,chemAtomSet)
            
            if not chemAtomSetSysName:
              chemAtomSetSysName = chemAtomSet.name
            else:
              chemAtomSetSysName = chemAtomSetSysName.sysName
            
            atomSetResName = ""
            
            for resName in atomToResName.hasEqChemAtomSets:
              (chainCode,seqCode,spinSystemId,seqInsertCode,atomSetName) = getNameInfo(resName)

              if chemAtomSetSysName == atomSetName:
              
                atomSetResName = resName
                break
                
            if atomSetResName:
            
              for nmrAtom in nmrAtoms:
                
                chemAtom = atomToResName.chemAtomDict[nmrAtom]
                chemAtomSysName = findChemAtomOrSetSysName(namingSystem,chemAtom)
                
                if not chemAtomSysName:
                  chemAtomSysName = chemAtom.name
                else:
                  chemAtomSysName = chemAtomSysName.sysName
                  
                for resName in atomToResName.hasEqChemAtomSets:

                  (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resName)

                  if chemAtomSysName == atomName:
                  
                    if self.verbose:
                      print "  Found equivalent atom group name match for %s - collapsing info on %s" % (resName,atomSetResName)
                    
                    #
                    # Copy associated info
                    #
                    
                    resonance = self.joinResonanceFromResNames([resName],atomToResName.isotopeCode)

                    for refResonance in self.resUnassignedNames[atomSetResName]:
                      copyResonanceInfo(resonance,refResonance,nmrShiftSelect = self.nmrShiftSelect, toResName = atomSetResName,verbose = self.verboseCopy)
                    
                    #
                    # Delete resName and resonance
                    #

                    self.deleteResonance(resonance,[resName])

                    index = atomToResName.hasEqChemAtomSets.index(resName)
                    atomToResName.hasEqChemAtomSets.pop(index)

                    index = atomToResName.resNameList.index(resName)
                    atomToResName.resNameList.pop(index)
                
    """
    set = 0
    for atomToResName in self.atomToResNameList:
      
      print set,len(atomToResName.nmrAtoms)
      
      for resName in atomToResName.resNameList:
        
        print "   ",resName,
        
        chemAtomSets = atomToResName.chemAtomSets[resName]
      
        if chemAtomSets:
          for chemAtomSet in chemAtomSets:
        
            print chemAtomSet.name,
          
        print
      
            
      set += 1
    """

  def compareNmrAtoms(self,subSet,fullSet):
  
    """
    Check whether subSet list is part of fullSet list
    """
  
    matches = 0
    
    for nmrAtom in subSet:
      if nmrAtom in fullSet:
        matches += 1
        
    if matches == len(subSet):
      return 1
    elif matches != 0:
      return 0
    else:
      return -1

  def linkResonancesToAtoms(self):
  
    """
    
    Loop over the all AtomToResName objects and handle the resonances themselves.
    At this stage checks will be made whether resonances with the same atom link
    already exist. If so, the information from the current resonance will be joined
    or linked with the existing resonance(s) (this depend on the exact case).
    
    """
    
    #
    # First sort atomToResName list by actual chain code, residue ID, and atom name
    #
    
    atomToResNameDict = {}
    for atomToResName in self.atomToResNameList:
      nmrAtomKey = [nmrAtom.name for nmrAtom in atomToResName.nmrAtoms]
      nmrAtomKey.sort()
      
      nmrRes = atomToResName.nmrAtoms[0].residue
      
      atomToResNameKey = (nmrRes.chain.code,nmrRes.seqId,tuple(nmrAtomKey))
      atomToResNameDict[atomToResNameKey] = atomToResName
    
    atomToResNameKeys = atomToResNameDict.keys()
    atomToResNameKeys.sort()
    
    self.atomToResNameList = []
    for atomToResNameKey in atomToResNameKeys:
      self.atomToResNameList.append(atomToResNameDict[atomToResNameKey])
      
    allLinkedResonances = []
      
    #
    # Loop over the atomToResNameList... now have to handle resonances themselves.
    #
   
    # This is for output
    oldChainCode = oldSeqId = None
  
    for atomToResName in self.atomToResNameList:
      
      nmrAtoms = atomToResName.nmrAtoms
      self.nmrRes = nmrAtoms[0].residue
      resonanceList = []
      resonanceResidueLink = None
      
      #
      # Print out residue type information
      #
      
      if self.verbose:
        printResidueInfo = False
        
        if oldChainCode != self.nmrRes.chain.code:
          oldChainCode = self.nmrRes.chain.code
          printResidueInfo = True
        if oldSeqId != self.nmrRes.seqId:
          oldSeqId = self.nmrRes.seqId
          printResidueInfo = True
          
        if printResidueInfo:
          print
          print drawBox("Assignments for residue %d, chain '%s' (%s, code %s)" % (self.nmrRes.seqId,self.nmrRes.chain.code,self.nmrRes.ccpCode,self.nmrRes.seqCode), indent = "  ")
          print
              
      #
      # Check if any resonance matches the current self.atomsFound
      #
     
      self.existingResonances = self.findResonancesMatchingAtoms(atomToResName.nmrAtoms)
      
      #
      # Residue level resonances handled here...
      #
      
      if isinstance(atomToResName,ResidueAtomToResName):

        refResName = atomToResName.resNameList[0]
        refResonance = self.resUnassignedNames[refResName][0]
        refResonance.isotopeCode = 'unknown'
        
        if self.existingResonances.resonances:
          # Use existing resonance if available...
          for resonance in self.existingResonances.resonances:
            copyResonanceInfo(resonance,refResonance,nmrShiftSelect = self.nmrShiftSelect,deleteOldInfo = 1,toResName = refResName,verbose = self.verboseCopy)
            for resName in atomToResName.resNameList[1:]:
              self.deleteResonance(resonance,[resName])

        else:
          # Create new one
          atomSets = []
          for chemAtomOrSetName in atomToResName.allAtomGroups:
            nmrAtoms = atomToResName.allAtomGroups[chemAtomOrSetName]
            atomSet = self.getOrCreateAtomSet(chemAtomOrSetName,nmrAtoms)
            atomSets.append(atomSet)

          resonanceSet = self.resonanceSetClass(self.resonanceParent,resonances = [refResonance], atomSets = atomSets)
        
        resonanceList.append(refResonance)
        resonanceResidueLink = self.nmrRes
     
      #
      # Normal atom level resonances handled here...
      #
      
      else:
  
        isotopeCode = atomToResName.isotopeCode
        self.extendResonanceList = []

        #
        # Set up the atom to resonance link(s)
        #

        resNameText = str(atomToResName.resNameList)

        #
        # Prochiral (both Leu HD* ('deep') and Asp HB* ('normal') type)
        #

        if atomToResName.isProchiral:

          if len(nmrAtoms) > 2:

            singleNumberAtoms = 3
            prochiralType = 'deep'
            getAtomSets = self.getDeepAtomSets
            otherChemAtomSets = atomToResName.hasEqChemAtomSets

          else:

            singleNumberAtoms = 1
            prochiralType = 'normal'
            getAtomSets = self.getNormalAtomSets
            otherChemAtomSets = atomToResName.hasEmptyChemAtomSets

          #
          # Set stereo assignment status. Note that if in fact only a
          # HB* (or QB or whatever) given it will NOT stereospecifically
          # assign.
          #

          if atomToResName.stereoAtoms:

            self.stereoAssign = 1
            stereoStatusText = 'stereo'

          else:

            self.stereoAssign = 0
            stereoStatusText = 'non-stereo'


          #
          # Deep prochiral (e.g. Val HG1*, HG2*) or
          # normal prochiral (e.g. Asp HB2, HB3)
          #

          if atomToResName.hasProchChemAtomSets:

            #
            # Prochiral chemAtomSet present (e.g. HG*)
            #

            #
            # Combine all resNames/resonances in this situation and
            # set the remaining resonance as reference
            #

            self.resonance = self.joinResonanceFromResNames(atomToResName.hasProchChemAtomSets,isotopeCode)

            #
            # Create or get two AtomSets and a new resonance.
            #

            atomSets = getAtomSets(atomToResName)

            if not self.linkToExistingResonances(nmrAtoms,atomSets,atomToResName.resNameList,resRequired = 2):

              if self.verbose:
                print "  Set %s as %s %s prochiral (with multi group)" % (resNameText,prochiralType,stereoStatusText)
              otherResonance = self.setProchiralGroup(atomSets)

            if otherChemAtomSets:

              #
              # Is basically situation of resonances already there, have
              # to link these up to them... reset self.existingResonances
              # and use existing code below
              #

              if self.verbose:
                print "    Resetting existing resonances..."

              for nmrAtom in nmrAtoms:

                if not self.atomResonances.has_key(nmrAtom):
                  self.atomResonances[nmrAtom] = []

                #
                # Has to have a link to a resonance by now...
                #

                tempAtomSet = self.getAtomSet(nmrAtom)

                for tempResonanceSet in tempAtomSet.sortedResonanceSets():
                  for tempResonance in tempResonanceSet.sortedResonances():

                    if tempResonance not in self.atomResonances[nmrAtom]:
                      self.atomResonances[nmrAtom].append(tempResonance)

              self.existingResonances = self.findResonancesMatchingAtoms(atomToResName.nmrAtoms)

          #
          # Handle other cases. This HAS to be an if to handle cases
          # where there's a prochiral (e.g. HD*) and equivalent (e.g. HD1*)
          # at the same time
          #
          # This can include a stereo assigned (e.g. HB2) with also a group (e.g. HB*)
          # as resonance
          #

          if len(nmrAtoms) > singleNumberAtoms:

            #
            # Create or get two AtomSets
            #

            atomSets = getAtomSets(atomToResName)

            #
            # Prochiral but no prochiral atom set
            #            

            if len(otherChemAtomSets) > 2:

              #
              # Selection window to group 'a' and 'b'
              # in order to end up with 2 or 1 hasEqChemAtomSets
              #

              interaction = self.multiDialog.ResonanceGroup(self.guiParent,otherChemAtomSets, self.nmrRes, self.allowPopups)
              resNameGroups = interaction.resNameGroups
              interaction.destroy()
              
              if not resNameGroups:
                print "   Error: Not grouping resonances with name %s because automatic grouping failed or no popup allowed." % (otherChemAtomSets)

              #
              # Have to group resonances and reset otherChemAtomSets, use first name as reference...
              #

              otherChemAtomSets = []

              for resNameGroup in resNameGroups:

                if resNameGroup:
                  refResName = resNameGroup[0]
                  refResonance = self.joinResonanceFromResNames([refResName],isotopeCode)

                  otherChemAtomSets.append(refResName)

                  for resName in resNameGroup[1:]:

                    resonance = self.joinResonanceFromResNames([resName],isotopeCode)

                    if self.verbose:
                      print "  Warning: Collapsing group resonances for %s to %s" % (resName,refResName)

                    copyResonanceInfo(resonance,refResonance,nmrShiftSelect = self.nmrShiftSelect,deleteOldInfo = 1,toResName = refResName,verbose = self.verboseCopy)

                    #
                    # Track stereo information
                    #

                    if atomToResName.stereoAtoms.has_key(resName):

                      atomToResName.stereoAtoms[refResName] = atomToResName.stereoAtoms[resName][:]
                      del(atomToResName.stereoAtoms[resName])

                    self.deleteResonance(resonance,[resName])

            if len(otherChemAtomSets) == 2:

              #
              # Assuming that these equivalent sets are different
              # Sorting them by resonance name: this means that the lower
              # number will always get a, the higher the b as non-stereo code!
              #

              resNames = otherChemAtomSets
              resNames.sort()

              if not self.stereoAssign:

                self.resonance = self.joinResonanceFromResNames([resNames[0]],isotopeCode)
                otherResonance = self.joinResonanceFromResNames([resNames[1]],isotopeCode)

                if not self.linkToExistingResonances(nmrAtoms,atomSets,[resNames[0]],resRequired = 2,otherResonance = otherResonance,otherResNames = [resNames[1]]):
                  if self.verbose:
                    print "  Set %s as %s prochiral, non-stereo, two resonances." % (resNameText,prochiralType)
                  self.setProchiralGroup(atomSets, otherResonance = otherResonance)

              else:

                if self.verbose:
                  print "  Linking up single stereo assigned resonances..."

                #
                # Resonances already exist: just link them up
                # Only problem that one might have a stereoAtom but the other not -
                # handling this situation first.
                #

                if len(atomToResName.stereoAtoms) == 1:

                  refResName = atomToResName.stereoAtoms.keys()[0]
                  refResNameIndex = atomToResName.resNameList.index(refResName)
                  nmrAtoms = atomToResName.stereoAtoms[refResName]

                  otherResName = atomToResName.resNameList[not(refResNameIndex)]

                  if len(nmrAtoms) == 1:
                    atomToResName.stereoAtoms[otherResName] = [returnProchiralAtom(nmrAtoms[0],chemCompVar = self.nmrRes.chemCompVar)]

                  else:
                    atomToResName.stereoAtoms[otherResName] = returnProchiralAtomGroup(nmrAtoms,chemCompVar = self.nmrRes.chemCompVar)[0]

                #
                # Now make the actual links
                #

                for resName in resNames:

                  self.resonance = self.joinResonanceFromResNames([resName],isotopeCode)
                  nmrAtoms = atomToResName.stereoAtoms[resName]

                  if not self.linkToExistingResonances(nmrAtoms,atomSets,[resName],resRequired = 1):
                    if self.verbose:
                      print "  Set %s as %s prochiral, stereo, two resonances." % (resName,prochiralType)

                    atomSet = nmrAtoms[0].atomSet
                    resonanceSet = self.resonanceSetClass(self.resonanceParent,resonances = [self.resonance], atomSets = [atomSet])
                    setCcpNmrResonanceName(self.resonance,nmrAtoms[0],atomName = atomSet.name)

                  #
                  # Set resonance list here: need both prochirals
                  #

                  if not self.extendResonanceList:            
                    resonanceList.append(self.resonance)
                  else:
                    resonanceList.extend(self.extendResonanceList)

            elif len(otherChemAtomSets) == 1:

              #
              # Always calling this one a. Also ask if necessary to create b.
              #

              resName = otherChemAtomSets[0]
              self.resonance = self.joinResonanceFromResNames([resName],isotopeCode)

              if not self.stereoAssign:

                singleProchiralStatus = self.getStatusOtherProchiral(atomSets,resName)

                if not singleProchiralStatus:
                  resRequired = 2
                else:
                  resRequired = 1

                if not self.linkToExistingResonances(nmrAtoms,atomSets,atomToResName.resNameList,resRequired = resRequired):

                  if not singleProchiralStatus:

                    if self.verbose:
                      print "  Set %s as %s prochiral, non-stereo, one resonance (copied)." % (resNameText,prochiralType)
                    self.setProchiralGroup(atomSets)

                  else:

                    if self.verbose:
                      print "  Set %s as %s prochiral, non-stereo, one resonance (kept single)." % (resNameText,prochiralType)

                    resonanceSet = self.resonanceSetClass(self.resonanceParent,resonances = [self.resonance], atomSets = atomSets)
                    resonanceName = createNonStereoName(atomSets[0].name,'a')
                    setCcpNmrResonanceName(self.resonance,atomSets[0].findFirstAtom(),atomName = resonanceName)

              else:

                #
                # Resonance already exist: just link it up
                #

                if self.verbose:
                  print "  Linking up single stereo assigned resonance..."

                nmrAtoms = atomToResName.stereoAtoms[resName]

                self.linkToExistingResonances(nmrAtoms,atomSets,[resName],resRequired = 1)

          else:

            #
            # Stereo assigned
            #

            self.resonance = self.joinResonanceFromResNames(atomToResName.resNameList,isotopeCode)

            #
            # Make direct link with correct atomSet OR atom (depending on case)
            #

            atomSets = getAtomSets(atomToResName)

            atomName = atomSets[0].name

            if not self.linkToExistingResonances(nmrAtoms,None,atomToResName.resNameList,resRequired = 1):

              if self.verbose:
                print "  Set %s as %s prochiral (stereo)" % (resNameText,prochiralType)

              self.createResonanceAtomLink(atomName,nmrAtoms)
              setCcpNmrResonanceName(self.resonance,nmrAtoms[0],atomName = atomName)


          if not self.extendResonanceList:            
            resonanceLoopList = [self.resonance]
          else:
            resonanceLoopList = self.extendResonanceList

          for loopResonance in resonanceLoopList:
            if loopResonance not in resonanceList:
              resonanceList.append(loopResonance)


        #
        # Not prochiral
        #

        else:      

          if atomToResName.hasEqChemAtomSets:

            #
            # Equivalent atom sets (methyls)
            #

            self.resonance = self.joinResonanceFromResNames(atomToResName.resNameList,isotopeCode)
            resonanceList.append(self.resonance)

            #
            # Make direct link with correct atomSet
            #

            nmrAtomGroup = atomToResName.nmrAtomGroups[0]
            atomName = atomToResName.chemAtomDict[nmrAtomGroup[0]].chemAtomSet.name

            if not self.linkToExistingResonances(nmrAtoms,None,atomToResName.resNameList,resRequired = 1):

              if self.verbose:
                print "  Set %s as multi group" % resNameText

              self.createResonanceAtomLink(atomName,nmrAtomGroup)
              setCcpNmrResonanceName(self.resonance,nmrAtomGroup[0],atomName = atomName)

          else:

            #
            # First determine if possibly equivalent and set status if necessary
            #

            if atomToResName.singleAtoms and atomToResName.isEquivalent == None:

              #
              # Check for each nmrAtomGroup whether only one relevant
              # nmrAtom available
              #

              singleResonanceStatus = 0

              for nmrAtomGroup in atomToResName.nmrAtomGroups:

                resNamesFound = {}

                for resName in atomToResName.singleAtoms.keys():

                  nmrAtom = atomToResName.singleAtoms[resName]

                  try:
                    index = nmrAtomGroup.index(nmrAtom)
                  except:
                    index = -1

                  if index > -1:
                    if not resNamesFound.has_key(index):

                      resNamesFound[index] = []

                    resNamesFound[index].append(resName)


                if resNamesFound:

                  #
                  # Cleanup: collapse all resonances onto each other for
                  # each nmrAtom
                  #

                  for index in resNamesFound.keys():

                    resonance = self.joinResonanceFromResNames(resNamesFound[index],isotopeCode)

                    for resName in resNamesFound[index][1:]:
                      del(atomToResName.singleAtoms[resName])
                      tIndex = atomToResName.resNameList.index(resName)
                      atomToResName.resNameList.pop(tIndex)

                  #
                  # Continue with checking status
                  #

                  nmrAtom = nmrAtomGroup[index]
                  otherAtomName = nmrAtomGroup[not(index)].name

                  singleResonanceStatus = self.getStatusResonance(nmrAtom,otherAtomName,resNamesFound[index][0])

                  if not singleResonanceStatus:

                    #
                    # Reset refResonance to group in atomToResName
                    #

                    chemAtomSet = atomToResName.chemAtomDict[nmrAtom].chemAtomSet

                    for index in resNamesFound.keys():

                      resName = resNamesFound[index][0]
                      atomToResName.groupChemAtomSets[resName] = [chemAtomSet]
                      atomToResName.hasGroupChemAtomSets.append(resName)
                      del(atomToResName.singleAtoms[resName])

                      if self.verbose:
                        print "  Reset status of resonance %s" % resName


            #
            # For amide-type with no singleAtoms, make sure that are treated
            # correctly (never a single atomSet for e.g. Asn HD11,HD12!!)
            #

            if not atomToResName.singleAtoms and len(atomToResName.nmrAtoms) > 1 and atomToResName.isEquivalent == False:

              #
              # Check for ARG HH*, TYR/PHE HD*|HE* situation
              #

              isGroupedAmide = 1

              if atomToResName.nmrAtoms[0].chemAtom.chemAtomSet:

                deepChemAtomSet = atomToResName.nmrAtoms[0].chemAtom.chemAtomSet.chemAtomSet

                if deepChemAtomSet and deepChemAtomSet.isEquivalent == False and deepChemAtomSet.isProchiral == False:

                  isGroupedAmide = 0


              if isGroupedAmide:

                #
                # They all have to share the same information: copy onto reference, then 
                # use existing or create new resonance
                #

                self.resonance = self.joinResonanceFromResNames(atomToResName.resNameList,isotopeCode)

                refResName = atomToResName.resNameList[0]
                atomSets = self.getNormalAtomSets(atomToResName)

                if not self.linkToExistingResonances(nmrAtoms,atomSets,[refResName],resRequired = 2):
                  if self.verbose:
                    print "  Set %s as grouped amide." % resNameText
                  otherResonance = self.resonanceClass(self.resonanceParent,isotopeCode = isotopeCode)
                  copyResonanceInfo(self.resonance,otherResonance,nmrShiftSelect = self.nmrShiftSelect,copyTrack = True,verbose = self.verboseCopy)
                  resonanceSet = self.resonanceSetClass(self.resonanceParent,resonances = [self.resonance,otherResonance], atomSets = atomSets)

                  #
                  # Set names
                  #

                  resonances = [self.resonance,otherResonance]

                  for i in range(0,2):

                    setCcpNmrResonanceName(resonances[i],atomSets[i].findFirstAtom(),atomName = atomSets[i].findFirstAtom().name)


                if not self.extendResonanceList:            
                  resonanceList.append(self.resonance)
                else:
                  resonanceList.extend(self.extendResonanceList)

                #
                # Clean up
                #

                for resName in atomToResName.resNameList:
                  del(atomToResName.groupChemAtomSets[resName])
                  index = atomToResName.hasGroupChemAtomSets.index(resName)
                  atomToResName.hasGroupChemAtomSets.pop(index)

            #
            # After changing status for possible equivalent atoms continue
            # as normal. This is for both aromatic and amide type (works
            # on isEquivalent = 0 or None)
            #

            if not atomToResName.isEquivalent and len(atomToResName.nmrAtoms) > 1:

              #
              # Aromatic/amide groups are a special case... need to check exactly what is going on
              #

              if atomToResName.singleAtoms:

                #
                # In this case there are single atom/resonances in the list:
                # make a list of nmrAtoms: some might have to be added
                # in case there is a chemAtomSet for the single atoms as
                # well
                #

                if not atomToResName.hasGroupChemAtomSets:

                  nmrAtomList = atomToResName.singleAtoms.values()

                else:

                  nmrAtomList = []

                  groupChemAtomSets = []

                  for groupChemAtomSetList in atomToResName.groupChemAtomSets.values():
                    for groupChemAtom in groupChemAtomSetList:
                      if groupChemAtom not in groupChemAtomSets:
                        groupChemAtomSets.append(groupChemAtom)

                  for nmrAtomGroup in atomToResName.nmrAtomGroups:

                    #
                    # Check if there is a groupChemAtomSet first
                    #

                    nmrAtom = nmrAtomGroup[0]
                    chemAtomSet = atomToResName.chemAtomDict[nmrAtom].chemAtomSet

                    if chemAtomSet not in groupChemAtomSets:
                      continue

                    #
                    # For valid groups, check if there are singleAtoms
                    #

                    for nmrAtom in nmrAtomGroup:

                      if nmrAtom in atomToResName.singleAtoms.values():

                        nmrAtomList.extend(nmrAtomGroup)
                        break

                #
                # Do the processing based on this list
                #

                chemAtomResonances = {}
                chemAtomResNames = {}
                deleteChemAtomSets = []

                for nmrAtom in nmrAtomList:

                  # 
                  # Also add single atoms that are not present
                  #

                  atomName = atomToResName.chemAtomDict[nmrAtom].name
                  setCode = 'Added'
                  self.resonance = None

                  #
                  # Check whether resonance already exists
                  #

                  for singleAtomResName in atomToResName.singleAtoms.keys():

                    if nmrAtom == atomToResName.singleAtoms[singleAtomResName]:

                      self.resonance = self.joinResonanceFromResNames([singleAtomResName],isotopeCode)
                      setCode = 'Set'
                      break

                  #
                  # If not, create new resonance
                  #

                  if not self.resonance:

                    self.resonance = self.resonanceClass(self.resonanceParent,isotopeCode = isotopeCode)
                    singleAtomResName = ''

                  #
                  # Make the link to the atom
                  #

                  if not self.linkToExistingResonances([nmrAtom],None,[singleAtomResName],resRequired = 1):

                    self.createResonanceAtomLink(atomName,[nmrAtom])
                    ccpName = setCcpNmrResonanceName(self.resonance,nmrAtom,atomName = atomName)

                    if not singleAtomResName:
                      singleAtomResName = ccpName

                    if self.verbose:
                      print "  %s '%s' as single aromatic/amide atom" % (setCode,singleAtomResName)

                  #
                  # Keep track
                  #

                  if atomToResName.chemAtomDict.has_key(nmrAtom):
                    chemAtom = atomToResName.chemAtomDict[nmrAtom]
                  else:
                    chemAtom = nmrAtom.chemAtom.chemAtom

                  chemAtomSet = chemAtom.chemAtomSet

                  chemAtomResonances[chemAtom] = self.resonance
                  chemAtomResNames[chemAtom] = singleAtomResName

                  if not self.extendResonanceList:            
                    resonanceList.append(self.resonance)
                  else:
                    resonanceList.extend(self.extendResonanceList)

                  if chemAtomSet not in deleteChemAtomSets:
                    deleteChemAtomSets.append(chemAtomSet)

                #
                # Link grouped atoms to relevant resonances
                #

                tempGroupChemAtomSets = atomToResName.hasGroupChemAtomSets[:]

                for resName in tempGroupChemAtomSets:

                  chemAtomSets = atomToResName.groupChemAtomSets[resName][:]
                  resonance = self.joinResonanceFromResNames([resName],isotopeCode)
                  deleteResonance = 0

                  for chemAtomSet in chemAtomSets:

                    if chemAtomSet in deleteChemAtomSets:
                      index = atomToResName.groupChemAtomSets[resName].index(chemAtomSet)
                      atomToResName.groupChemAtomSets[resName].pop(index)

                      for chemAtom in chemAtomSet.chemAtoms:

                        if chemAtomResonances.has_key(chemAtom):

                          # Possible that only one nmrAtom for this set

                          deleteResonance = 1

                          refResonance = chemAtomResonances[chemAtom]
                          refResName = chemAtomResNames[chemAtom]

                          copyResonanceInfo(resonance,refResonance,nmrShiftSelect = self.nmrShiftSelect,toResName = refResName,copyTrack = True,verbose = self.verboseCopy)

                  #
                  # Clean up
                  #

                  if not atomToResName.groupChemAtomSets[resName]:
                    del(atomToResName.groupChemAtomSets[resName])
                    index = atomToResName.hasGroupChemAtomSets.index(resName)
                    atomToResName.hasGroupChemAtomSets.pop(index)

                    if deleteResonance:
                      self.deleteResonance(resonance,[resName])

              #
              # Check if anything left...
              #

              if atomToResName.hasGroupChemAtomSets:

                #
                # Should be treated as groups. In case F,Y: HD*|HE* always split up into HD*,HE*
                # Since no single atoms listed, base on chemAtomSet(s)
                #

                chemAtomSetNamesResNames = {}

                for resName in atomToResName.hasGroupChemAtomSets:

                  chemAtomSets = atomToResName.groupChemAtomSets[resName]

                  for chemAtomSet in chemAtomSets:

                    chemAtomSetName = chemAtomSet.name

                    if not chemAtomSetNamesResNames.has_key(chemAtomSetName):
                      chemAtomSetNamesResNames[chemAtomSetName] = []

                    if len(chemAtomSets) == 1:

                      #
                      # Case of HD*, ...
                      #

                      chemAtomSetNamesResNames[chemAtomSetName].insert(0,resName)

                    else:

                      #
                      # Case of HD*|HE*
                      #

                      chemAtomSetNamesResNames[chemAtomSetName].append(resName)

                #
                # Set reference info
                #

                chemAtomSetNames = chemAtomSetNamesResNames.keys()

                createNewResonance = 0

                if len(chemAtomSetNames) == 2:

                  if chemAtomSetNamesResNames[chemAtomSetNames[0]] == chemAtomSetNamesResNames[chemAtomSetNames[1]]:

                    #
                    # Have to create a new one if only HD*|HE* type available
                    #

                    createNewResonance = 1

                #
                # Now set links, use first resName resonance in each group as reference
                #

                otherResonances = []

                for i in range(0,len(chemAtomSetNames)):

                  chemAtomSetName = chemAtomSetNames[i]

                  for nmrAtomGroup in atomToResName.nmrAtomGroups:
                    if nmrAtomGroup[0].chemAtom.chemAtomSet.name == chemAtomSetName:
                      break

                  if i == 1 and createNewResonance:

                    #
                    # Special case - only one resonance and HD*|HE* type given (e.g. QR in DIANA)
                    #

                    refResName = ''
                    self.resonance = self.resonanceClass(self.resonanceParent,isotopeCode = resonance.isotopeCode)
                    loopStart = 0
                    setCode = 'Added'

                  else:

                    refResName = chemAtomSetNamesResNames[chemAtomSetName][0]
                    self.resonance = self.joinResonanceFromResNames([refResName],isotopeCode)
                    loopStart = 1
                    setCode = 'Set'

                  #
                  # Here have to check if atomSet already exists... handle differently!!
                  #

                  hasAtomSets = []

                  for nmrAtom in nmrAtomGroup:

                    tempAtomSet = self.getAtomSet(nmrAtom)

                    if tempAtomSet:

                      #
                      # Ignore if same atomSet
                      #

                      if tempAtomSet in hasAtomSets:

                        hasAtomSets.pop(hasAtomSets.index(tempAtomSet))

                      else:

                        hasAtomSets.append(tempAtomSet)

                  #
                  # This means both atoms already have an atomset: use these and split
                  #

                  if hasAtomSets and len(hasAtomSets) == 2:

                    resonance = self.joinResonanceFromResNames(atomToResName.resNameList,isotopeCode)

                    newResonance = self.resonanceClass(self.resonanceParent,isotopeCode = isotopeCode)
                    copyResonanceInfo(resonance,newResonance,copyTrack = True,verbose = self.verboseCopy)

                    allResonances = [resonance,newResonance]

                    for j in range(len(allResonances)):

                      self.resonance = allResonances[j]
                      resonanceList.append(self.resonance)

                      atomSet = hasAtomSets[j]

                      #
                      # Make direct link with correct atomSet
                      #

                      nmrAtom = atomSet.findFirstAtom()
                      atomName = atomToResName.chemAtomDict[nmrAtom].name

                      if not self.linkToExistingResonances([nmrAtom],[atomSet],[resonance.name],resRequired = 1):

                        if self.verbose:
                          print "  Forced %s as single atom (atomSet already present)" % resNameText

                        self.createResonanceAtomLink(atomName,[nmrAtom])
                        setCcpNmrResonanceName(self.resonance,nmrAtom,atomName = atomName)

                  #
                  # No atomSets... just group and get on with it
                  #

                  else:

                    self.createResonanceAtomLink(chemAtomSetName,nmrAtomGroup)
                    ccpName = setCcpNmrResonanceName(self.resonance,nmrAtomGroup[0],chemAtomSetName)

                    if not refResName:
                      refResName = ccpName

                    if self.verbose:
                      print "  %s %s as grouped aromatic/amide atoms" % (setCode,refResName)

                    resonanceList.append(self.resonance)

                    for resName in chemAtomSetNamesResNames[chemAtomSetName][loopStart:]:

                      resonance = self.joinResonanceFromResNames([resName],isotopeCode)

                      if self.verbose:
                        print "  Warning: Collapsing grouped aromatic/amide resonances for %s to %s" % (resName,refResName)

                      copyResonanceInfo(resonance,self.resonance,nmrShiftSelect = self.nmrShiftSelect,toResName = refResName,copyTrack = True,verbose = self.verboseCopy)
                      otherResonances.append([resonance,resName])

                #
                # Delete resonances that are now useless
                #

                for (otherResonance,otherResName) in otherResonances:

                  if otherResonance not in resonanceList:

                    self.deleteResonance(otherResonance,[otherResName])

            else:

              self.resonance = self.joinResonanceFromResNames(atomToResName.resNameList,isotopeCode)
              resonanceList.append(self.resonance)

              #
              # Make direct link with correct atomSet
              #

              nmrAtom = atomToResName.nmrAtoms[0]
              atomName = atomToResName.chemAtomDict[nmrAtom].name

              if not self.linkToExistingResonances(nmrAtoms,None,atomToResName.resNameList,resRequired = 1):

                if self.verbose:
                  print "  Set %s as single atom" % resNameText

                self.createResonanceAtomLink(atomName,[nmrAtom])
                setCcpNmrResonanceName(self.resonance,nmrAtom,atomName = atomName)
 
      #
      # Print out resonance->atom link info for clarity
      #
      
      
      for resonance in resonanceList:
      
        if not resonance:
          if self.verbose:
            print("ERROR: resonance information missing.")
          continue

        resSet = resonance.resonanceSet

        if not resSet:
          if len(resonanceList) > 1 and self.verbose:
            print "ERROR: %s has unassigned resonance in list for printing - ignored" % resNameText
          continue
      
        resonances = list(resSet.sortedResonances())
      
        for tResonance in resonances:
          if tResonance not in allLinkedResonances:
            allLinkedResonances.append(tResonance)
      
        if self.verbose:
          
          atomSets = []

          if not resonanceResidueLink:
            outTextLen = max(len(resonanceList),len(resSet.atomSets),len(resonances))
          else:
            outTextLen = max(len(resonanceList),len(resonances))
            
          outText = outTextLen * ['']
          
          for asToSort in resSet.atomSets:
            asFound = False
            for asIndex in range(len(atomSets)):
              tas = atomSets[asIndex]
              if asToSort.name < tas.name:
                atomSets.insert(asIndex,asToSort)
                asFound = True
                break
            if not asFound:
              atomSets.append(asToSort)
          
          # Resonances

          for i in range(0,len(resonances)):
            resNameDict = getApplResNames(self.format,[resonances[i]], keyword = assign_kw)

            outText[i] += "    Resonance '%s' (%s) " % (resonances[i].name,str(resNameDict.keys()))
            
          # ResonanceSet

          # Setup for output formatting

          maxTextLen = len(outText[0])

          for i in range(1,len(outText)):

            if len(outText[i]) > maxTextLen:
              maxTextLen = len(outText[i])


          for i in range(0,len(outText)):

            textLen = len(outText[i])
            spaces = maxTextLen - textLen

            if i == 0: 
              outText[i] += spaces * ' ' + "-- RS --"
            else:
              if len(resonances) > i:
                startChar = '/'
              else:
                startChar = ' '

              # Special residue case          
              if not resonanceResidueLink and len(resSet.atomSets) > i:
                endChar = "\\"
              else:
                endChar = " "

              outText[i] += spaces * ' ' + startChar + "      " + endChar

          # Atomsets

          # Special residue case          
          if resonanceResidueLink:
            outText[i] += " Residue '%s.%d' -- all atoms" % (resonanceResidueLink.chain.code,resonanceResidueLink.seqCode)
          else:
            
            # Atom set names are already sorted higher up!
                          
            for i in range(0,len(atomSets)):
              outText[i] += " '%s' -- " % (atomSets[i].name)

            # Atoms (full name)

            for i in range(0,len(atomSets)):
              outText[i] += "'"
              atomNames = []
              for atom in atomSets[i].atoms:
                atomNames.append("%s" % atom.name)
              atomNames.sort()

              for atomName in atomNames:
                outText[i] += atomName + ','

              outText[i] = outText[i][:-1] + "' atom(s)"

          print

          for i in range(0,len(outText)):
            if len(resonances) > i or len(atomSets) > i:
              print outText[i]

          print

    self.numResonancesLinked += len(allLinkedResonances)
 
  def getOrCreateAtomSet(self,atomName,nmrAtomGroup):
  
    atomSet = getAtomSet(self.atomSetsAtomDict,atomName,nmrAtomGroup)
    
    if not atomSet:
      atomSet = self.atomSetClass(self.resonanceParent,name = atomName, atoms = nmrAtomGroup)      
      self.atomSetsAtomDict[atomSet] = nmrAtomGroup
       
    return atomSet

  def getDeepAtomSets(self,atomToResName):
  
    """
    Get 'deep' atom sets (e.g. Val HG*)
    """

    atomSets = []        

    for nmrAtomGroup in atomToResName.nmrAtomGroups:
      atomName = atomToResName.chemAtomDict[nmrAtomGroup[0]].chemAtomSet.name
      atomSets.append(self.getOrCreateAtomSet(atomName,nmrAtomGroup))

    return atomSets
            
  def getNormalAtomSets(self,atomToResName):

    """
    Get 'normal' atom sets (e.g. Asp HB*)
    """

    atomSets = []

    for nmrAtom in atomToResName.nmrAtoms:
      atomName = atomToResName.chemAtomDict[nmrAtom].name
      atomSets.append(self.getOrCreateAtomSet(atomName,[nmrAtom]))

    return atomSets

  def joinResonanceFromResNames(self,resNameList,isotopeCode):
    
    """
    Combines ALL resonances that come up in self.resUnassignedNames
    for the list of resNames given!!
    """

    refResName = resNameList[0]
    refResonance = self.resUnassignedNames[refResName][0]
    refResonance.isotopeCode = isotopeCode
    
    for resName in resNameList:
      
      resonances = self.resUnassignedNames[resName][:]

      for resonance in resonances:
      
        if resonance != refResonance:
        
          oldResName = getResNameText(resonance)

          if self.verbose:
            print "  Warning: Collapsing resonances for %s to %s" % (oldResName,refResName)
          copyResonanceInfo(resonance,refResonance,nmrShiftSelect = self.nmrShiftSelect,deleteOldInfo = 1,toResName = refResName,verbose = self.verboseCopy)
          
          self.deleteResonance(resonance,[resName])

    return refResonance

  def linkToExistingResonances(self,nmrAtoms,atomSets,resNames,resRequired = 1,otherResonance = None,otherResNames = None):

    """
    Try to link to existing resonances (very case dependent code)
      
    Returns 1 if linked and all OK
    Returns 0 if not
    """
    
    resNameText = str(resNames)
    
    self.extendResonanceList = []
        
    if self.existingResonances.numResonances > 0:

      resonance = self.resonance
    
      if not otherResonance:
      
        #
        # Have one resonance only as input...
        #
              
        if self.existingResonances.numResonances == 1:
        
          #
          # One resonance already exists
          #
          
          existingResonance = self.existingResonances.resonances[0]
            
          #
          # Copy info from new resonance to existing one
          #
          
          if self.verbose:
            print "  Copying information from %s to existing single resonance..." % (resNameText)
          copyResonanceInfo(resonance,existingResonance,nmrShiftSelect = self.nmrShiftSelect,toResName = str(self.existingResonances.resNames[0]),verbose = self.verboseCopy)
            
          #
          # If there's no full match, link existing resonanceSet to all atomSet(s)
          # (this is a safer approach)
          #
            
          if not self.existingResonances.fullAtomMatch and atomSets:
              
            existingResonance.resonanceSet.atomSets = atomSets
          
          #
          # If need one resonance, delete new resonance
          #

          if resRequired == 1:
            
            if self.verbose:
              print "    Deleting resonance %s..." % (resNameText)
            self.deleteResonance(resonance,resNames)
            self.resonance = existingResonance
            self.addToExtendResonanceList(existingResonance)          # For output at end of script...

          #
          # If need two resonances, link new resonance to same resonanceSet
          #
            
          elif resRequired == 2:
            
            if self.verbose:
              print "    Linking to same resonanceSet..."
            resonance.resonanceSet = existingResonance.resonanceSet
            

          return True


        elif self.existingResonances.numResonances == 2:
            
          #
          # Are the resonances linked to the same resonanceSet?
          #
          # Or, alternatively, is the resonance 'ambiguous' and are the existing resonances stereospecific?
          #
          
          sharedChemAtomSet = None
          
          # This check might not even be necessary, but safer anyway
          if len(self.existingResonances.resonanceSets) == 2:
            for i in range(0,2):
              if sharedChemAtomSet != False and self.existingResonances.atoms[i]:
                chemAtom = self.existingResonances.atoms[i][0].chemAtom
                if chemAtom.chemAtomSet:
                  if sharedChemAtomSet:
                    if sharedChemAtomSet != chemAtom.chemAtomSet:
                      sharedChemAtomSet = False
                  else:
                    sharedChemAtomSet = chemAtom.chemAtomSet
                else:
                  sharedChemAtomSet = False
              else:
                sharedChemAtomSet = False
          
          if (sharedChemAtomSet or len(self.existingResonances.resonanceSets) == 1) and resRequired == 2:

            #
            # If two resonances required and there's only one resonanceSet
            # (with two resonances), then copy the info over and delete the new resonance.
            #
          
            if self.verbose:
              print "  Copying information from %s to two existing resonances..." % (resNameText)

            for i in range(0,len(self.existingResonances.resonances)):
              existingResonance = self.existingResonances.resonances[i]
              copyResonanceInfo(resonance,existingResonance,nmrShiftSelect = self.nmrShiftSelect,toResName = str(self.existingResonances.resNames[i]),copyTrack = True,verbose = self.verboseCopy)
              self.addToExtendResonanceList(existingResonance)
              
            self.deleteResonance(resonance,resNames)
            self.resonance = existingResonance          # For output at end of script...
          
          else:

            #
            # Ambiguous situation: do a selection
            #

            if resRequired == 2:
              multi = 1
            else:
              multi = 0
              
            #
            # Check nmrAtoms - are they clearly linked to one resonance?
            #
            
            selectedResonances = []
            
            for i in range(0,len(self.existingResonances.resonances)):
              nmrAtoms.sort()
              resonanceAtoms = self.existingResonances.atoms[i]
              resonanceAtoms.sort()
              
              if nmrAtoms == resonanceAtoms:
                selectedResonances = [self.existingResonances.resonances[i]]
            
            if not selectedResonances or len(selectedResonances) > 1:
              selectedResonances = self.selectResonances(resNames,multi = multi)
              
            if not selectedResonances:
              print "  Error: no resonance selected... aborting this step."
              return False
            
            for selectedResonance in selectedResonances:
              copyResonanceInfo(resonance,selectedResonance,nmrShiftSelect = self.nmrShiftSelect,copyTrack = True,verbose = self.verboseCopy)
              self.addToExtendResonanceList(selectedResonance)
              
            self.deleteResonance(resonance,resNames)
            self.resonance = selectedResonance          # For output at end of script...
          
          return True
          
        else:
          
          if self.verbose:
            print "  Error: can only handle relinking to one or two resonances."
          return False
          
          # TODO: or do a selection anyway? But could get messy...


      else:
      
        #
        # Two resonances as input... select multi.
        # TODO: is it possible that two as input, one only as linked?
        #       in principle NOT - otherwise impossible to get full atom match... 
        #
        
        resonances = [resonance,otherResonance]
        resNamesList = [resNames,otherResNames]
        
        #
        # Case where you have two resonances to two atoms, and a similar setup already exists...
        #
        
        if len(resonances) == 2 and len(self.existingResonances.resonanceSets) == 1:
        
          if len(self.existingResonances.resonances) == 2:
          
            #
            # Special case: HB2,3 available after HB*. NOTE: will ONLY do this if the resName the same or does not exist for existingResonances!
            #
            # Will now also handle more complex cases where, e.g. HBa, HBb available when setting HB2 HB3 type stuff.
            #

            existingResNames = []

            for tempResNames in self.existingResonances.resNames:
              for tempResName in tempResNames:
                if tempResName not in existingResNames:
                  existingResNames.append(tempResName)

            setConnection = False
            existingResNames.sort()

            if len(existingResNames) <= 1:

              setConnection = True

              if self.verbose:
                if len(existingResNames) == 0:
                  addText = " (no original names available)"
                else:
                  addText = " to %s" % self.existingResonances.resNames[0]
                print "    Warning: special case - autoconnecting resonances %s%s." % (str(resNamesList),addText)

            elif len(existingResNames) == 2:

              matchResNameList = []
              sameResNames = 0
              newResNameList = [resNames,otherResNames]

              for rnIndex in range(0,2):
                tempResNames = self.existingResonances.resNames[rnIndex]
                tempResNames.sort()
                if len(tempResNames) == 1:
                  # Situation where each existing resonance has a resName that can be matched to the new ones (by atom name)
                  (v1,v2,v3,v4,atomName) = getNameInfo(tempResNames[0])
                  if atomName[-1] in '#*%':
                    atomName = atomName[:-1]

                  for nrnIndex in range(0,2):
                    newResNames = newResNameList[nrnIndex]
                    if len(newResNames) == 1:
                      (v1,v2,v3,v4,newAtomName) = getNameInfo(newResNames[0])
                      if newAtomName[-1] in '#*%':
                        newAtomName = newAtomName[:-1]

                      for i in range(2,0,-1):
                        if self.patt['onlydigit'].search(atomName[-i:]) and self.patt['onlydigit'].search(newAtomName[-i:]):
                          if atomName[-i:] == newAtomName[-i:]:
                            matchResNameList.append(nrnIndex)
                          break

                elif tempResNames == existingResNames:
                  # Situation where two resonances with same resNames (e.g. HB2 and HB3)
                  sameResNames += 1

              #
              # See if anything matched...
              #

              if sameResNames == 2:
                setConnection = True
                if self.verbose:
                  print "    Warning: special case - autoconnecting resonances %s with same double name." % (str(resNamesList))

              elif len(matchResNameList) == 2 and matchResNameList[0] != matchResNameList[1]:
                if matchResNameList[1] < matchResNameList[0]:
                  resNamesList.reverse()

                setConnection = True
                if self.verbose:
                  print "    Warning: special case - autoconnecting resonances %s to resonances with matching atomName index." % (str(resNamesList))

            if setConnection:           
              for i in range(0,len(self.existingResonances.resonances)):
                existingResonance = self.existingResonances.resonances[i]
                resonance = resonances[i]
                resNames = resNamesList[i]

                copyResonanceInfo(resonance,existingResonance,nmrShiftSelect = self.nmrShiftSelect,toResName = str(self.existingResonances.resNames[i]),copyTrack = True,verbose = self.verboseCopy)

                self.deleteResonance(resonance,resNames)
                self.addToExtendResonanceList(existingResonance)

              self.resonance = existingResonance          # For output at end of script...           

              return True
            
        
          #
          # Case where you have two resonances to two atoms, but only one existing resonance in similar setup
          #

          elif len(self.existingResonances.resonances) == 1:

            existingResNames = []

            for tempResNames in self.existingResonances.resNames:
              for tempResName in tempResNames:
                if tempResName not in existingResNames:
                  existingResNames.append(tempResName)

            setConnection = False
            existingResNames.sort()

            if len(existingResNames) == 0:

              setConnection = True

              if self.verbose:
                print "    Warning: special case - autoconnecting resonances %s (no existing resNames, creating one new resonance to existing resonanceSet)." % (str(resNamesList))

            elif len(existingResNames) == 1:

              setConnection = True
              matchResNameList = []
              newResNameList = [resNames,otherResNames]

              tempResNames = self.existingResonances.resNames[0]
              tempResNames.sort()
              if len(tempResNames) == 1:
                # Situation where each existing resonance has a resName that can be matched to the new ones (by atom name)
                (v1,v2,v3,v4,atomName) = getNameInfo(tempResNames[0])

                for nrnIndex in range(0,2):
                  newResNames = newResNameList[nrnIndex]
                  if len(newResNames) == 1:
                    (v1,v2,v3,v4,newAtomName) = getNameInfo(newResNames[0])

                    for i in range(2,0,-1):
                      if self.patt['onlydigit'].search(atomName[-i:]) and self.patt['onlydigit'].search(newAtomName[-i:]):
                        if atomName[-i:] == newAtomName[-i:]:
                          matchResNameList.append(nrnIndex)
                        break

              #
              # See if resNames have to be reversed...
              #

              if len(matchResNameList) == 1:
                if matchResNameList[0] == 1:
                  resNamesList.reverse()

              if self.verbose:
                print "    Warning: special case - autoconnecting resonances %s (trying to match to existing resNames, creating one new resonance to existing resonanceSet)." % (str(resNamesList))

            if setConnection:   

              existingResonance = self.existingResonances.resonances[0]
              resonance = resonances[0]
              resNames = resNamesList[0]

              copyResonanceInfo(resonance,existingResonance,nmrShiftSelect = self.nmrShiftSelect,toResName = str(self.existingResonances.resNames[0]),copyTrack = True,verbose = self.verboseCopy)

              self.deleteResonance(resonance,resNames)
              self.addToExtendResonanceList(existingResonance)
              self.resonance = existingResonance          # For output at end of script...

              otherResonance = resonances[1]
              resNames = resNamesList[1]

              otherResonance.resonanceSet = existingResonance.resonanceSet
              self.deleteUnassignedResonance(otherResonance,resNames)     

              return True
            
        
        #
        # Case where there are two resonances (ambiguous), and two existing resonances connected
        # stereospecifically to the right atoms... using chemAtomSysNames to 'backmatch'.
        #
        
        elif self.autoConnectStereo and len(resonances) == 2 and len(self.existingResonances.resonanceSets) == 2 and len(self.existingResonances.resonances) == 2:
          
          matchResNameList = []
          newResNameList = [resNames,otherResNames]
          
          for rnIndex in range(0,2):
            existingResSet = self.existingResonances.resonanceSets[rnIndex]
            if len(existingResSet.atomSets) == 1:
              atomSetName = existingResSet.findFirstAtomSet().name

              for nrnIndex in range(0,2):
                newResNames = newResNameList[nrnIndex]
                if len(newResNames) == 1:
                  (v1,v2,v3,v4,newAtomName) = getNameInfo(newResNames[0])

                  chemAtomSysNames = self.getAtomSysNames(newAtomName,self.namingSystemName,self.nmrRes.chemCompVar)
                  
                  if len(chemAtomSysNames) == 1 and chemAtomSysNames[0].atomName == atomSetName:
                    matchResNameList.append(nrnIndex)
       
          if len(matchResNameList) == 2 and matchResNameList[0] != matchResNameList[1]:
            if matchResNameList[1] < matchResNameList[0]:
              resonances.reverse()
              resNamesList.reverse()
          
            if self.verbose:
              print "    Warning: special case - autoconnecting resonances %s (using atom names of stereospecific existing resonances to connect)." % (str(resNamesList))

            for i in range(0,len(self.existingResonances.resonances)):
              existingResonance = self.existingResonances.resonances[i]
              resonance = resonances[i]
              resNames = resNamesList[i]

              copyResonanceInfo(resonance,existingResonance,nmrShiftSelect = self.nmrShiftSelect,toResName = str(self.existingResonances.resNames[i]),copyTrack = True,verbose = self.verboseCopy)

              self.deleteResonance(resonance,resNames)
              self.addToExtendResonanceList(existingResonance)

            self.resonance = existingResonance          # For output at end of script...           
            
            return True
            
        #
        # All other cases
        #
        
        for i in range(0,len(resonances)):
        
          resonance = resonances[i]
          resNames = resNamesList[i]
            
          selectedResonances = self.selectResonances(resNames,multi = 1)
          
          if selectedResonances:
            for selectedResonance in selectedResonances:
              copyResonanceInfo(resonance,selectedResonance,nmrShiftSelect = self.nmrShiftSelect,copyTrack = True,verbose = self.verboseCopy)
              self.addToExtendResonanceList(selectedResonance)
              
            self.existingResonances.setResNames(self.format)
            self.existingResonances.determineResonanceSets()
            
            self.deleteResonance(resonance,resNames)
            self.resonance = selectedResonance          # For output at end of script...
          
          else:
            if self.verbose:
              print "  Error: no resonances selected."
            return False
             
        return True
        
    return None
  
  def addToExtendResonanceList(self,resonance):
  
    addResonance = True
    
    for tRes in self.extendResonanceList:
      if tRes.resonanceSet == resonance.resonanceSet:
        addResonance = False
        break
        
    if addResonance:
      self.extendResonanceList.append(resonance)
  
  def selectResonances(self,resNames, multi = 0):
    
    while (None in resNames):
      resNames.remove(None)
    
    resName = string.join(resNames,',')
    selectedResonances = []
    
    #
    # Check if there's any mapping for this... using full resName
    #
    
    if self.useResLinkMapping:

      resLinkMapping = self.nameMapping.findFirstResLinkMapping(resonanceName = resName)
    
      if resLinkMapping:
      
        selectedResonanceNames = list(resLinkMapping.selectResonanceNames)
        
        if self.verbose:
          print "  Autoselected resonances '%s' for '%s'" % (str(selectedResonanceNames),resName)
      
        for i in range(0,len(self.existingResonances.resonances)):
        
          existingResonance = self.existingResonances.resonances[i]
          ccpNmrResNames = self.existingResonances.ccpNmrResNames[i]
          
          ccpNmrResNames.sort()
          
          if ccpNmrResNames == selectedResonanceNames:
            selectedResonances.append(existingResonance)
    
    #
    # If nothing found select manually
    #
    
    if not selectedResonances:

      if self.allowPopups:
        interaction = self.multiDialog.ResonanceSelect(self.guiParent,self.existingResonances, self.nmrRes, resName, multi)
        selectedResonances = interaction.selectedResonances
      
        #
        # Set mapping
        #
        
        selectedResonanceNames = []
        
        for selectedResonance in selectedResonances:
          selectedResonanceNames.extend(getApplResNames(ccpNmr_kw,[selectedResonance], keyword = assign_kw).keys())
        
        selectedResonanceNames.sort()
        resLinkMapping = NameMapping.ResonanceLinkMapping(self.nameMapping,resonanceName = resName,selectResonanceNames = selectedResonanceNames)
        
        interaction.destroy()
      
      else:
        tResNames = [ern for ern in self.existingResonances.resNames]
        print "  ERROR: no resonances selected from list %s because no user interaction allowed." % str(tResNames)
      
    
    return selectedResonances

  def deleteResonance(self,resonance,resNameList):           

    for resName in resNameList:
      self.deleteUnassignedResonance(resonance,resName)

    if not resonance.isDeleted:
      resonance.delete()
    elif self.verbose:
      print "  Warning: attempt to delete already deleted %s %s..." % (resonance.className,getResNameText(resonance))

  def deleteUnassignedResonance(self,resonance,resName):
    
    try:
      resonanceIndex = self.resUnassignedNames[resName].index(resonance)

      del(self.resUnassignedNames[resName][resonanceIndex])
    
    except:
      
      pass
    
    try:
    
      if not self.resUnassignedNames[resName]:
        del(self.resUnassignedNames[resName])
    
    except:
    
      pass

  def findResonancesMatchingAtoms(self,nmrAtoms):
  
    #
    # Determine if matching resonance already exists (via the atoms that were found)
    #
    
    existingResonances = ExistingResonances()
    
    for atom in nmrAtoms:

      if self.atomResonances.has_key(atom):
        
        for existingResonance in self.atomResonances[atom]:

          if existingResonances.resonances.count(existingResonance) == 0:
            existingResonances.addResonance(existingResonance)

          existingResonances.addAtom(existingResonance,atom)

    existingResonances.determineMatches(nmrAtoms)
    existingResonances.setResNames(self.format)
    existingResonances.determineResonanceSets()
    
    return existingResonances
  

  def createResonanceAtomLink(self,atomSetName,atoms):
    
    #
    # Create a new atomSet (or get existing one) linked to the atom(s)
    #
    
    atomSet = self.getOrCreateAtomSet(atomSetName,atoms)

    #
    # Create new resonanceSet
    #
    
    self.resonanceSetClass(self.resonanceParent,resonances = [self.resonance], atomSets = [atomSet])

  def getMatchingAtoms(self):

    #
    # Link to the relevant NmrAtoms (has to be done for both the automatic match and the handpicked one)
    #
    
    if not self.chemAtoms:
      return 0
      
    self.atomsFound = self.chemAtoms[:]
    
    for atom in self.nmrRes.atoms:
      
      chemAtom = atom.chemAtom
   
      # See if the chemAtom is in matching chemAtoms,
      # keep track if so... also maintain original order

      if chemAtom in self.chemAtoms:
        self.atomsFound[self.atomsFound.index(chemAtom)] = [atom,chemAtom]
        self.chemAtoms.remove(chemAtom)

    #
    # After automatic mapping or user atom choice, check the results...
    #
    
    chemAtomNames = []

    for chemAtom in self.chemAtoms:
      chemAtomNames.append(chemAtom.name)
        
    if self.atomsFound == self.chemAtoms:
      if self.verbose:
        print "No match with real atoms for chemAtoms %s (%s,%s)" % (str(chemAtomNames),self.seqCode,self.nmrRes.molResidue.ccpCode)
      return 0
      
    elif self.chemAtoms:
      if self.verbose:
        print "Error: ChemAtoms  %s do not match real atoms (%s,%s)." % (str(chemAtomNames),self.seqCode,self.nmrRes.molResidue.ccpCode)
      return 0
      
    else:
      return 1

  def hasAtomSysNames(self,namingSystem,coverage=0.8):
    
    atomSysNamesFound = False
    
    chemComp = namingSystem.chemComp
    chemAtoms = chemComp.sortedChemAtoms()
    chemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemComp,chemAtoms,namingSystem.name)
    
    if len(chemAtomSysNames) * 1.0 / len(chemAtoms) > coverage:
      atomSysNamesFound = True
    
    if not atomSysNamesFound and self.useEmptyNamingSystems:
      atomSysNamesFound = True
    
    return atomSysNamesFound

  def getAtomSysNames(self,atomName,namingSystemName,chemCompVar):
  
    """
    Gets the chemAtomSysName for atomName, namingSystemName and chemCompVar
    
    First look for complete match (but ONLY if that atomName exists for this
    particular chemCompVar! Otherwise problems with, e.g., H*), then loop through
    all names... is now a preferential match searcher to avoid some problems
    with sysNames that are the same (e.g. O for both O' and O for MSI naming...)
    """
    
    testPrint = False
    #if atomName == 'ZN':
    #  testPrint = True
    
    namingSystem = chemCompVar.chemComp.findFirstNamingSystem(name = namingSystemName)
    
    #if testPrint:
    #  print namingSystem
    
    if not self.ccvSysNames.has_key(chemCompVar):
      self.ccvSysNames[chemCompVar] = [[ca.name for ca in chemCompVar.chemAtoms],{}]
    
    if not self.ccvSysNames[chemCompVar][1].has_key(namingSystem):
      self.ccvSysNames[chemCompVar][1][namingSystem] = {}
      
      nsChemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemCompVar.chemComp,chemCompVar.chemAtoms,namingSystemName)
      nsChemAtomSysNames.extend(findAllSysNamesByChemAtomOrSet(chemCompVar.chemComp,chemCompVar.chemAtomSets,namingSystemName))

      #if testPrint:
      #  print nsChemAtomSysNames
                                
      for nsChemAtomSysName in nsChemAtomSysNames:
        self.ccvSysNames[chemCompVar][1][namingSystem][nsChemAtomSysName] = [nsChemAtomSysName.atomName,nsChemAtomSysName.sysName,nsChemAtomSysName.altSysNames]
    
    chemAtomSysNames = []
    
    if atomName in self.ccvSysNames[chemCompVar][0]:
      for casn in self.ccvSysNames[chemCompVar][1][namingSystem].keys():
        (cAtomName,cSysName,cAltSysNames) = self.ccvSysNames[chemCompVar][1][namingSystem][casn]
        
        #if testPrint:
        #  print cAtomName, cSysName, cAltSysNames
          
        if atomName == cAtomName and atomName == cSysName:
          chemAtomSysNames.append(casn)
          #if testPrint:
          #  print "FULLMATCH", atomName, casn

    if not chemAtomSysNames:
    
      for casn in self.ccvSysNames[chemCompVar][1][namingSystem].keys():
        (cAtomName,cSysName,cAltSysNames) = self.ccvSysNames[chemCompVar][1][namingSystem][casn]
        
        #if testPrint:
        #  print cAtomName, cSysName, cAltSysNames

        if atomName == cSysName:
          chemAtomSysNames.append(casn)
          #if testPrint:
          #  print "sysNameMatch", atomName, casn

      if not chemAtomSysNames:
      
        for casn in self.ccvSysNames[chemCompVar][1][namingSystem].keys():
          (cAtomName,cSysName,cAltSysNames) = self.ccvSysNames[chemCompVar][1][namingSystem][casn]

          if atomName in cAltSysNames:
            chemAtomSysNames.append(casn)
            #if testPrint:
            #  print "altSysNameMatch", atomName, casn
    
    return chemAtomSysNames

  def getIupacChemAtomSysNames(self):
  
    """
    Gets the IUPAC name chemAtomSysName
    
    Unless useIupacMatching == 2 this will ONLY return a chemAtomSysName
    that is NOT part of a chemAtomSet - this is to avoid problems with
    HB2/HB3 mixing between naming systems.
    """
  
    chemAtomSysNames = self.getAtomSysNames(self.atomName,'IUPAC',self.chemCompVar)
    matchingChemAtomSysNames = []
        
    if chemAtomSysNames:
    
      for chemAtomSysName in chemAtomSysNames:
    
        chemAtom = self.chemCompVar.findFirstChemAtom(name = chemAtomSysName.atomName, subType = chemAtomSysName.atomSubType)

        if self.localUseIupacMatching == 1 and chemAtom and not chemAtom.chemAtomSet:
          matchingChemAtomSysNames.append(chemAtomSysName)
          
        else:
          matchingChemAtomSysNames.append(chemAtomSysName)
    
    return matchingChemAtomSysNames

  def getAtomNameMapping(self):

    #
    # Check if there's any atom name mapping available for this resonance (name)
    #
  
    firstRes = self.resUnassignedNames[self.resName][0]
    mappingData = firstRes.findFirstApplicationData(application = self.format, keyword = ccpNmrMapping_kw)

    if mappingData:
      
      #
      # Reset atomName if so
      # chemAtomSysNames for IUPAC used to be set in reference data, now not any more!
      #
      
      self.atomName = mappingData.value
      
      if self.atomName:
      
        if self.verbose:
          print "  Mapping resonance '%s' to atomName %s" % (self.resName,self.atomName)
        self.namingSystemName = 'IUPAC'
        return True
        
      else:
      
        #
        # Mapping can also be set to 'Do not link'
        #
        
        if self.verbose:
          print "  Ignoring resonance '%s' - will not be linked" % self.resName
        return False

    else:
      
      #
      # Set original namingSystem if no mapping match
      #
      
      self.namingSystemName = self.namingSystemNameSelected
      return True


  def getAtomMatch(self):
    
    #
    # No match - have to link to the correct nmrAtom by interactive method
    #

    interaction = self.multiDialog.AtomSelect(self.guiParent, self.nmrRes, self.resName, self.resUnassignedNames[self.resName], modal = True)
    chemAtomOrSets = interaction.chemAtomOrSets
    self.propagateMapping = interaction.propagate
        
    interaction.destroy()

    if self.propagateMapping == None:
      if self.verbose:
        print "Can't link resonance name: %s (chemAtom selection dismissed)" % self.resName
      return 0
    
    elif chemAtomOrSets:
      
      if isinstance(chemAtomOrSets[0],ChemComp.ChemAtomSet):
      
        chemAtoms = []
        chemAtomSets = chemAtomOrSets
        chemCompVar = self.nmrRes.chemCompVar

        if len(chemAtomSets) == 1 and chemAtomSets[0].chemAtomSets:
          
          for linkedChemAtomSet in chemAtomSets[0].chemAtomSets:
            chemAtoms.extend(linkedChemAtomSet.chemAtoms)
                  
        else:
        
          for chemAtomSet in chemAtomSets:
            chemAtoms.extend(chemAtomSet.chemAtoms)

      else:
        chemAtoms = chemAtomOrSets
        chemAtomSets = []
      
      self.chemAtoms = chemAtoms
      self.chemAtomSets = chemAtomSets
      self.atomMappingName = chemAtomOrSets[0].name  # TODO: keep track of all names?
    
    else:
    
      #
      # 'Do not link' case - can also be propagated!
      #
      
      self.chemAtoms = None
      self.chemAtomSets = None
      self.atomMappingName = ""
    
    return True

  def useAtomMapping(self):
    
    #
    # Set namingSystem correctly if not set already (no mapping will
    # be available in this case
    #
    
    if self.nameMapping.namingSystem == 'Unknown':

      self.nameMapping.namingSystem = self.namingSystemName
      if self.verbose:
        print "  Warning: set namingSystem for nameMapping from 'Unknown' to %s." % self.namingSystemName
      
    else:
      
      if self.nameMapping.namingSystem != self.namingSystemName:
      
        self.messageReporter.showWarning("Warning","Naming system mismatch!! Using %s for mapping (not %s). This might lead to errors in the mapping..." % (self.nameMapping.namingSystem,self.namingSystemName),parent = self.guiParent)

    #
    # Then get the mapping
    #
    
    for atomMapping in self.nameMapping.atomMappings:
    
      atomMappingName = atomMapping.formatName
      atomName = atomMapping.name
      
      chainCode = atomMapping.chainCode
      moleculeName = atomMapping.moleculeName
      molSystemCode = atomMapping.molSystemCode
      
      #
      # Set chains to be mapped according to above...
      #
      
      if molSystemCode:

        if molSystemCode != self.molSystem.code:
          if self.verbose:
            print "  Warning: molSystem code %s in mapping does not match current molSystem." % molSystemCode

        tempMapChains = self.molSystem.sortedChains()
        
      elif moleculeName:
        
        tempMapChains = []
          
        for chain in self.molSystem.sortedChains():
          if chain.molecule.name == moleculeName:
            tempMapChains.append(chain)
            
      elif chainCode:
      
        tempMapChains = [self.molSystem.findFirstChain(code = chainCode)]
      
      else:
      
        if self.verbose:
          print "  Error: incomplete mapping - ignored."
        continue
      

      seqId = atomMapping.seqId
      ccpCode = atomMapping.ccpCode
          
      #
      # Set atom mapping...
      #
    
      mapChains = list(tempMapChains[:])
    
      for mapChain in tempMapChains:
      
        if not self.ccpChainDict.has_key(mapChain):
          mapChains.pop(mapChains.index(mapChain))
      
      #
      # Check if anything left
      #
      
      if not mapChains:
        
        if self.verbose:
          print "  Error: no chains left to map to - ignored"
        continue


      tempCodes = []
    
      #
      # Generic loop to set the mapping
      #
      
      for mapChain in mapChains:
      
        #
        # First set the residues to be mapped to...
        #
              
        if seqId != None:
          chainResidues = [mapChain.findFirstResidue(seqId = seqId)]

        elif ccpCode != None:
          chainResidues = []
          
          for residue in mapChain.sortedResidues():
            
            if residue.molResidue.ccpCode == ccpCode:
            
              chainResidues.append(residue)
        
        else:
        
          chainResidues = chain.sortedResidues()
        
        #
        # Then use the chain mappings to get the 'format resonance code'
        #
        
        chainMappings = self.ccpChainDict[mapChain]
        mapChainCode = chainMappings[0].formatChainCode
        
        (offset,seqInsertCode,seqIdToCodeMapping,firstSeqIds) = self.getChainMappingsSeqIdInfo(chainMappings)

        for tempNmrRes in chainResidues:
          
          if seqIdToCodeMapping:
            for firstSeqId in firstSeqIds:
              if tempNmrRes.seqId >= firstSeqId:
                offset = seqIdToCodeMapping[firstSeqId][0] - firstSeqId
                seqInsertCode = seqIdToCodeMapping[firstSeqId][1]
                break
          
          formatSeqCode = tempNmrRes.seqId + offset
          resName = getResName(mapChainCode,str(formatSeqCode),atomName,seqInsertCode = seqInsertCode)
          tempCodes.append(resName)

      #
      # Propagate mapping to unassignedNames (also check if it doesn't exist already)
      #
   
      for tempResName in self.resUnassignedNames:
      
        if tempResName in tempCodes:
    
          for tempResonance in self.resUnassignedNames[tempResName]:
          
            mappingData = tempResonance.findFirstApplicationData(application = self.format, keyword = ccpNmrMapping_kw)
          
            # 
            # Mapping exists and is not the same
            #
          
            if mappingData and mappingData.value != atomMappingName:
              if self.verbose:
                print "  Warning: propagating mapping to %s (previous mapping to atomName %s overwritten)." % (tempResName,mappingData.value)
              mappingData.value = atomMappingName
              continue
          
            elif not mappingData:
              if self.verbose:
                print "  Autoset mapping: resonance '%s' to atomName %s" % (resName,atomMappingName)
              mappingData = Implementation.AppDataString(application = self.format, keyword = ccpNmrMapping_kw, value = atomMappingName)
              tempResonance.addApplicationData(mappingData)

  def setAtomMapping(self,chain = None, molecule = None, molSystem = None, ccpCode = None, seqId = None):
  
    keywds = {}
    
    if molSystem:
      keywds['molSystemCode'] = molSystem.code
      
    elif molecule:
      keywds['moleculeName'] = molecule.name

    elif chain:
      keywds['chainCode'] = chain.code
      
    if ccpCode:
      keywds['ccpCode'] = ccpCode
    
    if seqId:
      keywds['seqId'] = seqId
    
    self.atomMapping = self.nameMapping.findFirstAtomMapping(name = self.atomName,
                                                             formatName = self.atomMappingName,
                                                             **keywds)
   
    if not self.atomMapping and self.atomMappingName:
    
      self.atomMapping = NameMapping.AtomMapping(self.nameMapping,name = self.atomName, formatName = self.atomMappingName, **keywds)
      return 1
    
    else:
      return 0

  def setAtomNameMapping(self):
    
    #
    # Set mapping
    #
    
    atomMappingName = self.atomMappingName
    propagateMapping = self.propagateMapping
    resName = self.resName
    chain = self.chain
    atomName = self.atomName
    seqId = self.nmrRes.seqId
    
    #
    # Reset atomName to original name - if get to this stage have to do this...
    #
    
    if self.origAtomName != self.atomName:
      self.atomName = atomName = self.origAtomName

    if atomMappingName:
      if self.verbose:
        print "  New mapping: resonance '%s' to atomName %s" % (resName,atomMappingName)
    else:
      if self.verbose:
        print "  New mapping: resonance '%s' will not be linked." % (resName)
    
    #
    # Propagate only in this chain, in all similar molecules, or in the molSystem?
    # Only use chains in self.chainDict: otherwise will have problems
    # First set up some stuff for setting the atom mapping
    #
    
    mappingKeyWds = {}
    
    if propagateMapping < 5:
      mappingKeyWds['ccpCode'] = self.ccpCode
    
    
    if propagateMapping in [4,7]:
      tempMapChains = chain.molSystem.sortedChains()
      mappingKeyWds['molSystem'] = chain.molSystem
      
    elif propagateMapping in [1,3,6]:
      tempMapChains = list(chain.molSystem.findAllChains(molecule = chain.molecule))
      mappingKeyWds['molecule'] = chain.molecule

      if propagateMapping == 1:
        mappingKeyWds['seqId'] = seqId
      
    else:
      tempMapChains = [chain]
      mappingKeyWds['chain'] = chain

    if not self.setAtomMapping(**mappingKeyWds):
      
      #
      # This means that this particular mapping already set: should be ignored
      #

      return 1
    
    #
    # New atom mapping...
    #
    
    mapChains = list(tempMapChains[:])
    
    for mapChain in tempMapChains:
      
      if not self.ccpChainDict.has_key(mapChain):
        mapChains.pop(mapChains.index(mapChain))
    
    tempCodes = []
    
    #
    # Only this residue is mapped
    #
    
    if propagateMapping == 0:
      tempCodes.append(resName)

    #
    # Set up list to propagate to residues with the same seqId in same molecules
    #
    
    elif propagateMapping == 1:
    
      for mapChain in mapChains:
        
        chainMappings = self.ccpChainDict[mapChain]
        mapChainCode = chainMappings[0].formatChainCode
        
        (offset,seqInsertCode,seqIdToCodeMapping,firstSeqIds) = self.getChainMappingsSeqIdInfo(chainMappings)
        
        for tempNmrRes in mapChain.sortedResidues():
          if tempNmrRes.molResidue.seqCode == seqId:
            if seqIdToCodeMapping:
              for firstSeqId in firstSeqIds:
                if tempNmrRes.seqId >= firstSeqId:
                  offset = seqIdToCodeMapping[firstSeqId][0] - firstSeqId
                  seqInsertCode = seqIdToCodeMapping[firstSeqId][1]
                  break
                  
            formatSeqCode = seqId + offset
            tempCodes.append(getResName(mapChainCode,str(formatSeqCode),atomName,seqInsertCode = seqInsertCode))

    #
    # Set up list to propagate to all residues of this type (in chain or in molSystem)
    #
    
    elif propagateMapping < 5:
    
      for mapChain in mapChains:
      
        chainMappings  = self.ccpChainDict[mapChain]
        mapChainCode = chainMappings[0].formatChainCode
        
        (offset,seqInsertCode,seqIdToCodeMapping,firstSeqIds) = self.getChainMappingsSeqIdInfo(chainMappings)
        
        for tempNmrRes in mapChain.sortedResidues():
          if tempNmrRes.molResidue.ccpCode == self.ccpCode:
            if seqIdToCodeMapping:
              for firstSeqId in firstSeqIds:
                if tempNmrRes.seqId >= firstSeqId:
                  offset = seqIdToCodeMapping[firstSeqId][0] - firstSeqId
                  seqInsertCode = seqIdToCodeMapping[firstSeqId][1]
                  break
                  
            formatSeqCode = tempNmrRes.seqId + offset
            tempCodes.append(getResName(mapChainCode,str(formatSeqCode),atomName,seqInsertCode = seqInsertCode))
      
    #
    # Set up list to propagate to all residues irrespective of type (in chain or in molSystem)
    #
    
    elif propagateMapping >= 5:
      for mapChain in mapChains:

        chainMappings = self.ccpChainDict[mapChain]
        mapChainCode = chainMappings[0].formatChainCode
        
        (offset,seqInsertCode,seqIdToCodeMapping,firstSeqIds) = self.getChainMappingsSeqIdInfo(chainMappings)

        for tempNmrRes in mapChain.sortedResidues():
          if seqIdToCodeMapping:
            for firstSeqId in firstSeqIds:
              if tempNmrRes.seqId >= firstSeqId:
                offset = seqIdToCodeMapping[firstSeqId][0] - firstSeqId
                seqInsertCode = seqIdToCodeMapping[firstSeqId][1]
                break
          formatSeqCode = tempNmrRes.seqId + offset
          tempCodes.append(getResName(mapChainCode,str(formatSeqCode),atomName,seqInsertCode = seqInsertCode))
    
    #
    # Propagate mapping to unassignedNames (also check if it doesn't exist already)
    #
    
    for tempResName in self.resUnassignedNames:
      
      if tempResName in tempCodes:
    
        for tempResonance in self.resUnassignedNames[tempResName]:
          
          mappingData = tempResonance.findFirstApplicationData(application = self.format, keyword = ccpNmrMapping_kw)
          
          # 
          # Mapping exists and is not the same
          #
          
          if mappingData and mappingData.value != atomMappingName:
            if self.verbose:
              print "  Warning: propagating mapping to %s (previous mapping to atomName %s overwritten)." % (tempResName,mappingData.value)
            tempResonance.removeApplicationData(mappingData)
            mappingData = None
          
          if not mappingData:
            mappingData = Implementation.AppDataString(application = self.format, keyword = ccpNmrMapping_kw, value = atomMappingName)
            tempResonance.addApplicationData(mappingData)

    if atomMappingName:
      return 1
    else:
      return 0


  def getChainMappingsSeqIdInfo(self,chainMappings):
  
    #
    # Gets offset (if only one chain mapping for the chain) or
    # detailed mapping info (if multiple mappings for the chain -
    # basically if 'jumps' in the sequence numbering)
    #
    
    seqIdToCodeMapping = {}
    firstSeqIds = []
    seqInsertCode = defaultSeqInsertCode
    offset = 0

    if len(chainMappings) == 1:
      offset = chainMappings[0].formatFirstSeqCode - chainMappings[0].firstSeqId
      seqInsertCode = chainMappings[0].formatSeqInsertCode

    else:
      for chainMapping in chainMappings:
        seqIdToCodeMapping[chainMapping.firstSeqId] = (chainMapping.formatFirstSeqCode,chainMapping.formatSeqInsertCode)
        firstSeqIds.append(chainMapping.firstSeqId)
      firstSeqIds.sort()
      firstSeqIds.reverse()
    
    return (offset,seqInsertCode,seqIdToCodeMapping,firstSeqIds)
    
  def getMatchingChemAtoms(self,chemAtomSysNames):
  
    ccpCode = chemAtomSysNames[0].namingSystem.chemComp.ccpCode
    
    #
    # Do a quick check to see if there are multiple matches...
    #
    
    validChemAtomSysNames = []
    
    for chemAtomSysName in chemAtomSysNames:
      
      chemAtom = self.chemCompVar.findFirstChemAtom(name = chemAtomSysName.atomName, subType = chemAtomSysName.atomSubType)
      
      if chemAtom:
        validChemAtomSysNames.append((chemAtomSysName,chemAtom,'chemAtom'))
      else:
        chemAtomSet = self.chemCompVar.findFirstChemAtomSet(name = chemAtomSysName.atomName, subType = chemAtomSysName.atomSubType)
        
        if chemAtomSet:
                
          #
          # Check if valid chemAtomSet for this particular one...
          # note that matches of type HZ1, HZ2 for HZ* (without HZ3) are also included
          #
          
          validChemAtoms = []
          chemAtomList = chemAtomSet.chemAtoms
          
          for chemAtom in chemAtomList:
            if chemAtom not in self.chemCompVar.chemAtoms:
              validChemAtoms.append(None)
            else:
              validChemAtoms.append(chemAtom)
              
          if None not in validChemAtoms or (len(validChemAtoms) - validChemAtoms.count(None)) >= 2 :
            validChemAtomSysNames.append((chemAtomSysName,chemAtomSet,'chemAtomSet'))
    
    #
    # Prepare everything...
    #
    
    if not validChemAtomSysNames:
      
      atomNames = []
      
      for chemAtomSysName in chemAtomSysNames:
        atomNames.append(chemAtomSysName.atomName)
      
      if self.verbose:
        print "  Error: no match for chemAtomSysName %s, atom(s) %s! (molType %s, ccpCode %s)" % (
      
                                                  chemAtomSysName.sysName,
                                                  str(atomNames),
                                                  chemAtomSysName.namingSystem.chemComp.molType,
                                                  ccpCode)
                                                  
      return False
      
    if len(validChemAtomSysNames) > 1:
    
      if self.verbose:
        print "  Warning: Multiple matches for search on atom name %s (molType %s, ccpCode %s). Picked first valid atom" % (
      
                                                  self.atomName,
                                                  chemAtomSysName.namingSystem.chemComp.molType,
                                                  ccpCode),
    
      for (chemAtomSysName,chemAtomOrSet,chemAtomType) in validChemAtomSysNames:
      
        if isinstance(chemAtomOrSet,ChemComp.ChemAtomSet):
          
          if self.verbose:
            print  "set %s." % chemAtomSysName.atomName
          validChemAtomSysName = (chemAtomSysName,chemAtomOrSet,'chemAtomSet')
          break

      if len(validChemAtomSysNames) > 1:
        if self.verbose:
          print "%s." % validChemAtomSysNames[0][0].atomName
        validChemAtomSysName = validChemAtomSysNames[0]

    else:
     
      validChemAtomSysName = validChemAtomSysNames[0]
      
    #
    # Set all the details for the match...
    #
    
    chemAtoms = []
    
    (chemAtomSysName,chemAtomOrSet,chemAtomType) = validChemAtomSysName
    
    if chemAtomType == 'chemAtom':
    
      #
      # Direct match to an atom
      #
      
      chemAtom = chemAtomOrSet
      chemAtoms.append(chemAtom)
      chemAtomSets = []
        
      #
      # Check if this chemAtom is part of a chemAtomSet...
      #
      
      chemAtomSet = chemAtoms[-1].chemAtomSet
      
      if chemAtomSet:
      
        if not self.stereoAssign and chemAtomSet.isProchiral:

          #
          # For prochiral, not stereo assigned have to use both chemAtoms
          #

          chemAtoms = chemAtomSet.sortedChemAtoms()

          #
          # But keep track of the original match for export later
          #

          for resonance in self.resUnassignedNames[self.resName]:

            # TODO is chemAtomSysName.atomName the correct name to set?
            setUniqueAppData('AppDataString',resonance,self.format,'mappingSingleMatch',chemAtomSysName.atomName)
        
        elif chemAtomSet.isEquivalent == None:
        
          #
          # For chemAtomSets with unknown equivalence status (e.g. aromatic
          # HE1), set special value to keep track. The FIRST chematom
          # in the list is the original one.
          #
          
          chemAtomList = chemAtomSet.sortedChemAtoms()
          
          for otherChemAtom in chemAtomList:
            if chemAtoms[0] != otherChemAtom:
              chemAtoms.append(otherChemAtom)
              
          self.unknownEquivalent = 1
          
        if chemAtomSet.isEquivalent:
        
          #
          # If protons are considered equivalent, need the chemAtomSet unless makeMethylsEquivalent is False
          # This is necessary for checking if only an HG11 is reported for a Val HG1* for example
          # Also use ALL atoms for this set (basically adding them here)
          #
      
          chemAtomSets = [chemAtomSet]
          chemAtomList = chemAtomSet.sortedChemAtoms()
          
          if len(chemAtoms) != len(chemAtomList) and self.makeMethylsEquivalent:
            
            if self.verbose:
              print "  Adding missing equivalent chemAtoms for %s: %s." % (ccpCode,self.resName)

            chemAtoms = chemAtomList
            
          #
          # Also check whether deep prochiral - need to reset values for
          # prochiral, non-stereoassigned: if HG12 only reported for
          # a valine then have to add missing chemAtoms and chemAtomSets
          # Note: code same as bit lower
          #
          
          deepChemAtomSet = chemAtomSet.chemAtomSet
        
          if deepChemAtomSet:
        
            if not self.stereoAssign and deepChemAtomSet.isProchiral and self.makeMethylsEquivalent:
            
              allChemAtomSets = deepChemAtomSet.sortedChemAtomSets()
              otherChemAtomSet = allChemAtomSets[not allChemAtomSets.index(chemAtomSet)]

              #
              # If prochiral and not stereo assigned, need both chemAtomSets
              # and all chemAtoms
              #

              if self.verbose:
                print "  Setting prochiral chemAtomSets for %s: %s." % (ccpCode,self.resName)

              (chemAtoms,chemAtomSets) = self.setDeepProchiral(deepChemAtomSet)

              for resonance in self.resUnassignedNames[self.resName]:
                setUniqueAppData('AppDataString',resonance,self.format,'mappingSingleMatch',chemAtomSet.name)
              
    else:
    
      #
      # For chemAtomSets get all atoms that match
      #
      
      chemAtomSet = chemAtomOrSet
      chemAtomSets = [chemAtomSet]
      
      if chemAtomSet.chemAtomSets:

        #
        # This is a 'deep' match (eg HG* for Valine)
        #

        for tempChemAtomSet in chemAtomSet.sortedChemAtomSets():
          chemAtoms.extend(tempChemAtomSet.chemAtoms)

      else:

        chemAtoms.extend(chemAtomSet.sortedChemAtoms())

        #
        # Check if prochiral (e.g. a Val HG1*)
        #

        deepChemAtomSet = chemAtomSet.chemAtomSet

        if deepChemAtomSet:

          if not self.stereoAssign and deepChemAtomSet.isProchiral:

            allChemAtomSets = deepChemAtomSet.sortedChemAtomSets()
            otherChemAtomSet = allChemAtomSets[not allChemAtomSets.index(chemAtomSet)]
            
            #
            # If prochiral and not stereo assigned, need both chemAtomSets
            # and all chemAtoms
            #

            (chemAtoms,chemAtomSets) = self.setDeepProchiral(deepChemAtomSet)

            for resonance in self.resUnassignedNames[self.resName]:
              setUniqueAppData('AppDataString',resonance,self.format,'mappingSingleMatch',chemAtomSet.name)
  
    self.chemAtomSets = chemAtomSets
    self.chemAtoms = chemAtoms

    return True
 
  def setDeepProchiral(self,deepChemAtomSet):
    
    #
    # TODO: in ccpnmr.format.general.Util?
    #
    
    chemAtomSets = deepChemAtomSet.sortedChemAtomSets()
    chemAtoms = []
    
    for chemAtomSet in chemAtomSets:
      chemAtoms.extend(chemAtomSet.sortedChemAtoms())
      
    return (chemAtoms,chemAtomSets)                    

  def nmrShiftSelect(self,shift,otherShift,toResName):

    shiftValue = shift.value
    shiftError = shift.error

    otherShiftValue = otherShift.value
    otherShiftError = otherShift.error
    
    toRes = otherShift.resonance
    
    if not toResName or toResName == '[]':
      toResName = toRes.name
                
    if not self.forceShiftMerge and abs(shiftValue - otherShiftValue) > chemShiftAllowedDiff and self.allowPopups:

      interaction = self.multiDialog.MeasurementSelect(self.guiParent,toResName,toRes,shiftValue,shift.parentList,'Shift')
      
      newShiftValue = interaction.selection

    else:

      if self.verbose:
        print "      Merging shift values %.3f and %.3f for resonance %s." % (shiftValue,otherShiftValue,toResName)
        print

      newShiftValue = (otherShiftValue + shiftValue) / 2
    
    #
    # Set the shift error
    #
    
    if newShiftValue == shiftValue:
      otherShift.error = shiftError
      
    elif newShiftValue != otherShiftValue:
      otherShift.error = max(shiftError,otherShiftError)
      
  # TODO: pseudoatoms for peakLists are same as set of atoms, but for restraint lists NOT SO!
  #       important to keep that in mind for scripts!!
  
  #
  # TODO: what should be returned?! A list of the linked atoms? Or the resonances?
  #
  

  def changeFormatChainCodes(self,formatChainCodes):
  
    if not formatChainCodes:
      return

    #
    # Get a list of all resonances and the TYPE of resonance
    #

    if not self.nmrConstraintStore:
      resonanceList = self.nmrProject.sortedResonances()
      self.resonanceClass = Nmr.Resonance
      self.resonanceSetClass = Nmr.ResonanceSet
      self.resonanceParent = self.project
      
      self.atomSetClass = Nmr.AtomSet
      
      atomSetsList = self.nmrProject.sortedAtomSets()
      
    else:
      resonanceList = self.nmrConstraintStore.sortedFixedResonances()
      self.resonanceClass = NmrConstraint.FixedResonance
      self.resonanceSetClass = NmrConstraint.FixedResonanceSet
      self.resonanceParent = self.nmrConstraintStore
      
      self.atomSetClass = NmrConstraint.FixedAtomSet
      
      atomSetsList = self.nmrConstraintStore.sortedFixedAtomSets()
      
    self.atomSetsAtomDict = {}
    for atomSet in atomSetsList:
      atoms = list(atomSet.atoms)
      atoms.sort()
      self.atomSetsAtomDict[atomSet] = atoms

    allResNamesDict = {}
 

    for resonance in resonanceList:

      #
      # If there is no resonance-atom link, set the unassigned resonance names dictionary
      #

      updateResonanceNamesDict(resonance,allResNamesDict,self.format)          
    
    allResNamesList = allResNamesDict.keys()
    allResNamesList.sort()
    
    #
    # formatChainCodes = {' ': ['A','B'],...}
    #
    
    formatChainCodeList = formatChainCodes.keys()
    formatChainCodeList.sort()
    
    for formatChainCode in formatChainCodeList:
    
      #
      # Get a list of all resonances for this formatChainCode
      #
      
      origResNamesDict = self.getResNamesByFormatChain(formatChainCode,allResNamesList,allResNamesDict)
      origResNamesList = origResNamesDict.keys()
      origResNamesList.sort()
      
      #
      # Get rid of the application data so it doesn't get copied...
      #
      
      for origResName in origResNamesList:
        for resonance in origResNamesDict[origResName]:

          assignData = resonance.findAllApplicationData(application = self.format, keyword = assign_kw)
          
          for assignDatum in assignData:
            resonance.removeApplicationData(assignDatum)
      
      #
      # Now get all chains that this code has to be mapped to...
      # keep original resonances until all copies made, then remove.
      #
      
      mapToChainList = formatChainCodes[formatChainCode]
      mapToChainList.sort()
      
      for mapToChain in mapToChainList:
      
        mapResNamesDict = self.getResNamesByFormatChain(mapToChain,allResNamesList,allResNamesDict)

        for origResName in origResNamesList:
        
          (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(origResName,verbose = 0)
          
          origResonances = origResNamesDict[origResName]
          matchName = getResName(mapToChain,seqCode,atomName,seqInsertCode = seqInsertCode)
          
          if mapResNamesDict.has_key(matchName):
            
            if self.verbose:
              print " Copying '%s' to '%s'..." % (origResName,matchName)
            matchResonances = mapResNamesDict[matchName]
            
            for origResonance in origResonances:
              for matchResonance in matchResonances:
                copyResonanceInfo(origResonance,matchResonance,toResName = matchName,verbose = self.verboseCopy)

          else:
            
            if self.verbose:
              print " Creating '%s' from '%s'..." % (matchName,origResName)

            for origResonance in origResonances:
              matchResonance = self.resonanceClass(self.resonanceParent,isotopeCode = origResonance.isotopeCode)
              matchResonance.addApplicationData(Implementation.AppDataString(application = self.format, keyword = assign_kw, value = matchName))
              copyResonanceInfo(origResonance,matchResonance,toResName = matchName,verbose = self.verboseCopy)

      #
      # Clean up...
      #
      
      for origResName in origResNamesList:      
        for resonance in origResNamesDict[origResName]:
          resonance.delete()

  def getResNamesByFormatChain(self,formatChainCode,resNamesList,resNamesDict):
  
    tempResNamesDict = {}
    
    for resName in resNamesList:

      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resName,verbose = 0)

      if chainCode == formatChainCode and seqCode != None and atomName:

        tempResNamesDict[resName] = resNamesDict[resName]
    
    return tempResNamesDict
    
