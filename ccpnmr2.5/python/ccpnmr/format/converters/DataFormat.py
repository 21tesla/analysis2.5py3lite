"""
======================COPYRIGHT/LICENSE START==========================

DataFormat.py: Main format conversion script. Contains all generic routines.

Copyright (C) 2008-2012 Wim Vranken (Vrije Universiteit Brussel)

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

import os, string, sys, traceback, operator

#
# Import general functions (normal and Api-dependent)
#

from memops.universal.Util import returnInt, returnInts
from memops.universal.Util import returnFloat
from memops.universal.Util import isBigEndian
from memops.universal.Util import drawBox
from memops.general.Util import returnMemopsWord, returnMemopsLine
from memops.general.Implementation import ApiError
from memops.universal.Io import joinPath, splitPath
from memops.universal.Io import normalisePath
from ccp.general.Constants import refRepositoryName, shiftIdentityTolerance

from ccp.general.Geometry import calcPlanarCoord, convertDegreeRadian
from ccp.general.Io import setDataSourceDataStore, getDataSourceFileName
from ccp.general.Io import getChemComp, getChemCompCoord
from ccp.general.ChemCompOverview import chemCompStdDict

from ccp.general.Util import setUniqueAppData, getDescriptorElementMatches, getDescriptorFromDict
from ccp.general.Util import MakeAtomSets, findAllSysNamesByChemAtomOrSet, findChemAtomSysName
from ccp.general.Util import findChemAtomOrSetSysName
from ccp.general.Util import findChemCompVarSysName, setCurrentStore, getResonancesFromPairwiseConstraintItem
from ccp.general.Util import getAmbigNamingPriority, getChemAtomSetFromAtomSets

from ccp.general.Constants import linkEndDict, linkingEnds, linkCodes
from ccp.general.Constants import linkTorsions, linkBonds
from ccp.general.Constants import standardBackboneAtoms
from ccp.general.Constants import secStrucInfo_kw

from ccp.format.general.Constants import bioPolymerCodes, chainCodeString
from ccp.format.general.Constants import nucleotideList, stableBondTypes
from ccp.format.general.Constants import defaultSeqInsertCode, defaultAltLoc

from ccp.format.general.formatIO import GenericChemComp, GenericAtom

from ccpnmr.format.general.Util import getResName, getIsotopes, getApplResNames
from ccpnmr.format.general.Util import getResNameText, getNameInfo, getShiftRange
from ccpnmr.format.general.Util import getSoftwareDefinitions, createSelection

# TODO: following one will go at one point...
from ccpnmr.format.general.Util import createPpmFreqDataSource

from ccpnmr.format.general.Constants import defaultMolCode, coordinateImportSourceName
from ccpnmr.format.general.Constants import tagSep, IoSetupList
from ccpnmr.format.general.Constants import chemShiftAllowedDiff
from ccpnmr.format.general.Constants import formatNamingDict
from ccpnmr.format.general.Constants import allFormatsDict
from ccpnmr.format.general.Constants import exportChainCodeDict
from ccpnmr.format.general.Constants import code1LetterToCcpCodeDict
from ccpnmr.format.general.Constants import ccpNmr_kw, originalFormat_kw
from ccpnmr.format.general.Constants import assign_kw, peakNum_kw
from ccpnmr.format.general.Constants import importGeneral_kw, atomSerial_kw
from ccpnmr.format.general.Constants import dateTimeStamp_kw, ioCycleTag_kw, seqAssignTag_kw

from ccpnmr.format.general.TopShared import TopShared

# Moved functions:
from ccp.lib.StructureLib import makeEnsemble, checkChemAtomsForDeepSet

#
# Get linkResonances
#

from ccpnmr.format.process.linkResonances import linkResonances

#
# Import user interaction code. Note that a lot of this code will not
# work unless a GUI environment is available (error messages will be
# printed however)
#

from ccpnmr.format.general.userInteraction import setupDataEntry
from ccpnmr.format.general.userInteraction import setupMessageReporter
from ccpnmr.format.general.userInteraction import setupMultiDialog

#
# Import Ccp Api and related
#

import memops.api.Implementation as Implementation
from ccp.api.molecule import ChemComp, MolSystem, Molecule, MolStructure
from ccp.api.general import Instrument, Citation, Affiliation, Method
from ccp.api.nmr import Nmr, NmrConstraint, NmrEntry
from ccpnmr.api import NameMapping

#
# Class ResonanceToAtoms
#
    
class ResonanceToAtoms:

  """

  This class handles export of resonances to atom names. It keeps track of the resonance,
  the corresponding chain and seqId, and the selected atomName (and atomSetName). Other
  resonances linked to this atom are stored in otherLinkedResonances, and prochiral and
  possible grouped resonances (e.g. HD21/22 for Asn) are stored in otherGroupResonances.

  """  

  def __init__(self,resonance,chain,seqId,atomName,atomSetName = None,chemAtomOrSet = None,originalResName = None, atomType = None):
  
    self.resonance = resonance
    self.chain = chain
    self.seqId = returnInt(seqId)
    self.atomName = atomName   # This could be single or for a set!!
    self.atomSetName = atomSetName
    self.atomType = atomType
    
    self.chemAtomSet = None   # Should these be lists?
    self.chemAtom = None
    
    self.useAtomSetName = 0
    
    if chemAtomOrSet:
    
      if chemAtomOrSet.className[-3:] == 'Set':
        self.chemAtomSet = chemAtomOrSet
        
      else:
        self.chemAtom = chemAtomOrSet
    
    self.otherLinkedResonances = []
    self.otherGroupResonances = []
    
    self.originalResName = originalResName
  
  def clone(self):
  
    tags = ('resonance', 'chain', 'seqId', 'atomName',)
    result = ResonanceToAtoms(*(getattr(self,tag) for tag in tags))
    
    for tag in ('chemAtom', 'chemAtomSet', 'originalResName','atomSetName', 
                'useAtomSetName', 'atomType'):
      setattr(result, tag, getattr(self,tag))
    
    for tag in ('otherLinkedResonances', 'otherGroupResonances'):
      setattr(result, tag, list(getattr(self,tag)))
    #
    return result
  
  # WARNING: cannot use this class as pointer in dict if it has __cmp__ method...
  #def __cmp__(self,object):
  
  #  returnCode = 0
  
  #  if isinstance(object,ResonanceToAtoms):
  #    if self.chain.code > object.chain.code:
  #      returnCode = 1
  #    elif self.chain.code < object.chain.code:
  #      returnCode = -1
  #    else:
  #      if self.seqId > object.seqId:
  #        returnCode = 1
  #      elif self.seqId < object.seqId:
  #        returnCode = -1
  #      else:
  #        if self.atomName > object.atomName:
  #          returnCode = 1
  #        elif self.atomName < object.atomName:
  #          returnCode = -1

  #  return returnCode
  
  def getAtomId(self):
     
    """
    Returns chain, seqId, atomName
    """
  
    return [self.chain,self.seqId,self.atomName,self]
  
  def getAtomSetId(self):
  
    """
    Returns chain, seqId, atomSetName
    """

    if self.atomSetName:
      
      return [self.chain,self.seqId,self.atomSetName,self]
    
    else:
    
      return None
    
  def addOtherLinkedResonance(self,resonance):
  
    """
    Adds resonance that is also linked to the chosen atom
    """

    self.otherLinkedResonances.append(resonance)

  def addOtherGroupResonance(self,resonance):
  
    """
    This is for both prochiral and possible groups (e.g. HD21/22, HE1/2, ...)
    """
  
    self.otherGroupResonances.append(resonance)

  def getChemAtomSet(self):
  
    chemAtomSet = self.chemAtomSet
    
    if chemAtomSet and chemAtomSet.chemAtomSet:
      chemAtomSet = chemAtomSet.chemAtomSet

    if not chemAtomSet and self.chemAtom:
      chemAtomSet = self.chemAtom.chemAtomSet
    
    return chemAtomSet

  def getResidue(self):
  
    return self.chain.findFirstResidue(seqId = self.seqId)

#
# Class ResonanceToSpinSystems
#
    
class ResonanceToSpinSystems:

  """

  This class handles export of resonances to based on spin system (ResonanceGroup) organisation

  """  

  def __init__(self,resonance,spinSystem,atomNames,linkedSpinSystems = None):
  
    self.resonance = resonance
    
    self.spinSystem = spinSystem
    self.linkedSpinSystems = linkedSpinSystems
    
    #self.chain = chain

    self.atomNames = atomNames
        
    #self.otherLinkedResonances = []
    #self.otherGroupResonances = []
    
    #self.originalResName = originalResName


#
# Keyword definitions for generic read/write functions defined in DataFormat
#
# Special format-specific keywords are defined in the format-specific class.
#
# First keyword:  function name
# Second keyword: keyword name
# Tuple:          (default,mandatory,explanatory text)
#

IOkeywords = {

  'readSequence': {
    
    'version':                (None,False,'The format version number as string (if relevant, it will be used).'),
    'molSystem':              (None,False,'The new chain will automatically be linked to this MolSystem object.'),
    'molName':                (None,False,'Name for the new molecule (only works if monomer, homodimer, ...).'),
    'numIdenticalChains':     (0,   False,'Number of identical chains for a homodimer, homotrimer, ...\n when reading in a single (monomer) sequence.'),
    'seqCodeStart':           (None,False,'Start sequence codes at this number (monomers only).'),
    'chainCodes':             (None,False,"List of chainCodes to be applied (default 'A','B','C',...)."),
    'useOriginalChainCode':   (None,False,"Use original (author defined) chain code if available. Relevant for archiving formats."),
    'addNameMappings':        (None,False,'Dictionary that can be passed in to map atom names to CCPN atoms.'),
    'autoCreateChemComps':    (False,False,'If set to True (on), unknown chemComps will be automatically created (if atom names are available!).'),
    'ignoreUnknownChemComps': (False,False,'If set to True (on), unknown chemComps will be automatically ignored.'),
    'ignoreSingleOther':      (False,False,"If set to True (on), single residue molecules of type 'other' will be ignored. Useful if, for example, a PDB file with has lots of waters that you want to ignore."), 
    'ignoreResNames':         (None,False,"List of residue names that should be ignored when reading sequence (e.g. ['HOH'])."), 
    'autoLink':               (False,False,'If set to True (on), the sequence elements will be fully automatically determined.'),
    'chemCompMappings':       (None,False,"Dictionary that can be passed in to map chemComp names from the format to CCPN chemComp names. Takes the form: {'formatCode': (molType,ccpCode)}"),
    'chemCompPath':           (None,False,'Path to directory that contains all chemComp information (if locally available).'),
    'saveChemCompLocal':      (True,False,'Save any downloaded chemComps locally. If set to False, will save in refData directory.'),
    'minimalPrompts':         (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.')
  
  },
  
  'writeSequence': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'chains':                   (None,True,'The list of chains to be exported to file.'),
    'resetMapping':             (True,False,'If set to False (off), the existing chain mapping will be used for export\n (otherwise it can be reset).'),
    'useCcpnChainInfo':         (False,False,'If set to True (on), the CCPN chain and sequence codes will be directly used for export.'),
    'useOriginalData':          (False,False,'If set to True (on), original format-specific data (from the import)\n will be used for export.'),
    'forceDefaultChainMapping': (False,False,'If set to True (on), the default chain mapping\n (between the format and the data model)\n will be automatically used.'),
    'noWrite':                  (False,False,'If set to True (on), no file will be written.'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.')
  
  },
  
  'readChemComps': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'molTypes':                 (None,False,'Molecule types for the new chemComps'),
    'ccpCodes':                 (None,False,'Unique CCP codes for the new chemComps'),
    'chemCompPath':             (None,False,'Path to directory that contains all chemComp information (if locally available).'),
    'saveChemComp':             (False,False,'If set to True (on), the new chemComp will be instantly saved after creation.'),
    'saveChemCompLocal':        (True,False,'Save any downloaded chemComps locally. If set to False, will save in refData directory.'),
    'createAtomSets':           (True,False,'If set to False (off), no atomSets will be created for the new chemComp.'),
    'makeNamingSystem':         (None,False,'If set, will create a naming system for all the atom names based on the original name.'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.')
  
  },

  'writeChemComp': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'noWrite':                  (False,False,'If set to True (on), no file will be written.'),
    'copyCoordinates':          (True,False,'If set to False (off), the chemCompCoord files will not be copied if found in the archive or online.'),
    'chemCompVar':              (None,False,'ChemCompVar to write to the external format.'),
    'forceNamingSystemName':    (None,False,'If set, will use this naming system instead of the default one for the atom names. Leave as None for default!'),
    'coordSystem':              (None,False,'Name of the coordinate system to use for getting the atom coordinates.'),
    'writeLinkAtoms':           (True,False,'If set to False (off), will not try to create atoms and coordinates for link atoms'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.')
  
  },

  'readCoordinates': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'molSystem':                (None,False,'The new structures will automatically be linked\n to this MolSystem object.'),
    'strucGen':                 (None,False,'The new structures will automatically be linked\n to this StructureGeneration object.'),
    'structureEnsemble':        (None,False,'The new structures will automatically be linked\n to this StructureEnsemble object.'),
    'maxNum':                   (999,False,'If multiple models are available, the maximum number\n of models that should be read.'),
    'linkAtoms':                (True,False,'If set to False (off), unrecognized coordinate atoms will not be linked.'),
    'matchToAtomNames':         (False,False,'If set to True (on), will try to match the names of the coordinate atoms to the actual atom names if none found for the naming system.'),
    'forceChainMappings':       (None,False,'Dictionary that can be passed in to automate the coordinate\n to data model chain/residue mapping.'),
    'forceDefaultChainMapping': (False,False,'If set to True (on), the default chain mapping\n (between the format and the data model)\n will be automatically used.'),
    'forceReadSequence':        (False,False,'If set to True (on), will read sequence from file as well.'),
    'forceNamingSystemName':    (None,False,'If set, will use this naming system instead of the default one. Leave as None for default!'),
    'addNameMappings':          (None,False,'Dictionary that can be passed in to map coordinate atom names to CCPN atoms.'),
    'autoCreateChemComps':      (False,False,'If set to True (on), unknown chemComps will be automatically created (if atom names are available!).'),
    'ignoreUnknownChemComps':   (False,False,'If set to True (on), unknown chemComps will be automatically ignored (only relevant if sequence read).'),
    'ignoreSingleOther':        (False,False,"If set to True (on), single residue molecules of type 'other' will be ignored. Useful if, for example, a PDB file with has lots of waters that you want to ignore."), 
    'ignoreResNames':           (None,False,"List of residue names that should be ignored when reading coordinates (e.g. ['HOH'])."), 
    'chemCompMappings':         (None,False,"Dictionary that can be passed in to map chemComp names from the format to CCPN chemComp names. Takes the form: {'formatCode': (molType,ccpCode)}"),
    'chemCompPath':             (None,False,'Path to directory that contains all chemComp information (if locally available).'),
    'saveChemCompLocal':        (True,False,'Save any downloaded chemComps locally. If set to False, will save in refData directory.'),
    'trackAtomSerials':         (False,False,'If set to True (on), it will track the atom serials of the first model of the file that is read in.'),
    'resetMapping':             (False,False,'If set to True (on), the existing chain mapping will be reset.'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.')
  
  },

  'writeCoordinates': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'structures':               (None,True,'The list of models from a structure ensemble to be exported to file.'),
    'exportChains':             (None,False,'If set, only the chains in this list will be exported.'), 
    'useOriginalData':          (False,False,'If set to True (on), original format-specific data (from the import)\n will be used for export.'),
    'noWrite':                  (False,False,'If set to True (on), no file will be written.'),
    'resetMapping':             (False,False,'If set to False (off), the existing chain mapping will be used for export\n (otherwise it can be reset).'),
    'useCcpnChainInfo':         (False,False,'If set to True (on), the CCPN chain and sequence codes will be directly used for export.'),
    'forceDefaultChainMapping': (False,False,'If set to True (on), the default chain mapping\n (between the format and the data model)\n will be automatically used.'),
    'forceNamingSystemName':    (None,False,'If set, will use this naming system instead of the default one. Leave as None for default!'),
    #'useFixedAtomSerials':      (None,False,'If set to a format name, will use tracked atom serials for this format if available.'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.')
  
  },

  'readMeasurements': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'measurementType':          (None,True,'The type of measurement list to be created.'),
    'measurementList':          (None,False,'A valid subclass of a MeasurementList (e.g. ShiftList, ...).\n New information will be added to this list\n (and/or old information replaced).'),
    'experiments':              (None,False,'A list of Nmr.Experiment objects.\n The MeasurementList will be linked to these Experiments.'),
    'method':                   (None,False,'A list of Method.Software objects.\n The MeasurementList will be linked to these Softwares.'),
    'measurementListName':      (None,False,'Name for the measurement list.'),
    'unit':                     (None,False,'The unit for measurement values.'),
    'measurementAllowedDiff':   (0.0,False, 'If the difference between two measurements for the same resonance is higher\n than this value, the user will have to select the correct value.'),
    'resonanceMatching':        (True,False,'If set to False (off), the script will not try to match\n the newly created resonances to existing ones (in the Data Model).'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.'),
    'presetResonanceMapping':   (False,False,'If set to True (on), will keep using first (or previously passed-in) Resonance:Atom mapping'),
  
  },
 
  'readHExchProtections':    {   # Specific for proton exchange protection values
  },

  'readHExchRates':    {         # Specific for proton exchange rates    
  },

  'readJCouplings':    {   # Specific for J coupling constants
    
  },

 'readNoes':    {   # Specific for NOEs
    
    #'coherenceType':            (None,False,"Can pass in the coherence type. Is either 'SQ', 'DQ', 'ZQ', 'antiphase', 'crosscorrelation'."),
    'specFreq':                 (None,False,"Spectrometer frequency at which values were recorded, in MHz"),
    #'refValueType':             (None,False,"Type of Intensity measurement"),

  },
  
  'readRdcs':    {   # Specific for RDC measurements
    
  },

  'readShiftAnistropy':    {   # Specific for chemical shift anisotropy
    
    'specFreq':                 (None,False,"Spectrometer frequency at which values were recorded, in MHz"),
    
  },

  'readShifts': {    # Specific to chemical shift reading only
                  
    'autoMergeShifts':          (False,False,'If set to True (on), shifts belonging to the same resonance and shift list will be automatically merged and get the average value.'),
    'measurementAllowedDiff':   (chemShiftAllowedDiff,False, 'If the difference between two measurements for the same resonance is higher\n than this value, the user will have to select the correct value.'),
                  
  }, 
  
  'readT1s':    {   # Specific for T1 relaxation
    
    'coherenceType':            (None,False,"Can pass in the coherence type. Is either 'z' or 'zz'."),
    'specFreq':                 (None,False,"Spectrometer frequency at which values were recorded, in MHz"),
    
  },
  
  'readT1Rhos': {  # Specific for T1 rho relaxation
      
    #'coherenceType':            (None,False,"Can pass in the coherence type. Is either 'SQ', 'DQ', 'ZQ', 'antiphase', 'crosscorrelation'."),
    'specFreq':                 (None,False,"Spectrometer frequency at which values were recorded, in MHz"),
    'tempCalibMethod':          (None,False,"Method used to calibrate the temperature unit on the NMR spectrometer used to collect the experimental data"),
    'tempControlMethod':        (None,False,"Method used to control the temperature of the sample to eliminate heating caused by the long spin-lock pulse"),
    
  },

  'readT2s':    {   # Specific for T2 relaxation
    
    #'coherenceType':            (None,False,"Can pass in the coherence type. Is either 'SQ', 'DQ', 'ZQ', 'antiphase', 'crosscorrelation'."),
    'specFreq':                 (None,False,"Spectrometer frequency at which values were recorded, in MHz"),
    'tempCalibMethod':          (None,False,"Method used to calibrate the temperature unit on the NMR spectrometer used to collect the experimental data"),
    'tempControlMethod':        (None,False,"Method used to control the temperature of the sample to eliminate heating caused by the long spin-lock pulse"),
    
  },
  
  'writeMeasurements': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'measurementList':          (None,True,'The measurement list to be exported to file.'),
    'measurementType':          (None,True,'The type of measurement list to be exported.'),
    'chains':                   (None,False,'A list of MolSystem.Chain objects.\n Only measurements belonging to these chains\n will be exported to file.'),
    'useOriginalData':          (False,False,'If set to True (on), original format-specific data (from the import)\n will be used for export.'),
    'useOriginalResNames':      (False,False,'If set to True (on), original atom names (from the import)\n will be used for export.'),
    'compressResonances':       (True,False,'If set to False (off), the output atom names will not be compressed\n (e.g. HB2 and HB3 will be written out separately,\n even if they share the same value for a measurement).'),
    'individualAtoms':          (False,False,'If set to True (on), atom names will be fully written out (e.g. an Ala HB* is written in individual HB1, HB2, HB3 contributions)'),
    'individualAtomsIfNoSet':   (False,False,'If set to True (on), atom names for an atom set will be fully written out if no correct atomSysName found (e.g. an Ala HB* that has no matching name in the selected naming system is written in individual HB1, HB2, HB3 contributions)'),
    'ignoreChemCompsWithNoSysNames': (False,False,'If set to True (on), residues that do not have the required atom naming system information will not be exported. Default (False) is to use atom names from CCPN instead.'), 
    'noWrite':                  (False,False,'If set to True (on), no file will be written.'),
    'forceDefaultChainMapping': (False,False,'If set to True (on), the default chain mapping\n (between the format and the data model)\n will be automatically used.'),
    'forceNamingSystemName':    (None,False,'If set, will use this naming system instead of the default one. Leave as None for default!'),
    'resetMapping':             (False,False,'If set to False (off), the existing chain mapping will be used for export\n (otherwise it can be reset).'),
    'useCcpnChainInfo':         (False,False,'If set to True (on), the CCPN chain and sequence codes will be directly used for export.'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.'),
    'presetResonanceMapping':   (False,False,'If set to True (on), will keep using first (or previously passed-in) Resonance:Atom mapping'),
    
  },
  
  'writeHExchProtections': {},  # Specific for proton exchange protection values
  'writeHExchRates':    {},    # Specific for proton exchange rates    
  'writeJCouplings':      {},  # Specific for J Couplings
  'writeNoes':      {},  # Specific for NOEs
  'writeRdcs':      {},  # Specific for RDC measurements
  'writeShiftAnistropy': {},   # Specific for chemical shift anisotropy
  'writeShifts':   {},   # Specific to chemical shift writing only
  'writeT1s':      {},   # Specific for T1 relaxation
  'writeT1Rhos':    {},  # Specific for T1 rho relaxation
  'writeT2s':      {},   # Specific for T2 relaxation
 
  'readPeaks': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'dataSource':               (None,False,'A valid Nmr.DataSource object. The new peak list\n will be created using this object.'),
    'dataDimRefs':              (None,False,'A list of Nmr.DataDimRef objects. The peak list\n will be linked to these DataDimRefs.'),
    'peakList':                 (None,False,'A valid Nmr.PeakList object. New information will be \nadded to this list (and/or old information replaced).\n Use with caution!'),
    'specName':                 (None,False,'Name for the peak list.'),
    'resonanceMatching':        (True,False,'If set to False (off), the script will not try to match\n the newly created resonances to existing ones (in the Data Model).'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.'),
    'presetResonanceMapping':   (False,False,'If set to True (on), will keep using first (or previously passed-in) Resonance:Atom mapping'),
  
  },

  'writePeaks': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'peakLists':                (None,True,'A list of valid Nmr.PeakList objects to be exported to file.'),
    'dataDimRefs':              (None,False,'A list of Nmr.DataDimRef objects. The peak list\n will be exported using the order of these dataDimRefs.'),
    'useOriginalNumbers':       (True,False,'If set to False (off), original format-specific peak numbers\n (from the import) will NOT be used for export.'),
    'useOriginalData':          (False,False,'If set to True (on), original format-specific data\n (from the import) will be used for export.'),
    'useOriginalResNames':      (False,False,'If set to True (on), original atom names\n (from the import) will be used for export.'),
    'compressResonances':       (True,False,'If set to False (off), the output atom names will not be compressed\n (e.g. HB2 and HB3 will be written out separately,\n even if their peakDimContrib can be grouped (to HB*)).'),
    'individualAtoms':          (False,False,'If set to True (on), atom names will be fully written out (e.g. an Ala HB* is written in individual HB1, HB2, HB3 contributions)'),
    'individualAtomsIfNoSet':   (False,False,'If set to True (on), atom names for an atom set will be fully written out if no correct atomSysName found (e.g. an Ala HB* that has no matching name in the selected naming system is written in individual HB1, HB2, HB3 contributions)'),
    'ignoreChemCompsWithNoSysNames': (False,False,'If set to True (on), residues that do not have the required atom naming system information will not be exported. Default (False) is to use atom names from CCPN instead.'), 
    'writeAssignments':         (True,False,'If set to False (off), peak assignments are not written to the output file.'),
    'skipMultiAssignments':     (False,False,'If set to True (on), multiple peak assignments are not written to the output file. even for the format(s) taht support this.'),
    'useXeasyDimCodes':         (False,False,'If set to True (on), peakList dimCodes will be set in Xeasy format (NmrView only).'),
    'noWrite':                  (False,False,'If set to True (on), no file will be written.'),
    'forceDefaultChainMapping': (False,False,'If set to True (on), the default chain mapping\n (between the format and the data model)\n will be automatically used.'),
    'forceNamingSystemName':    (None,False,'If set, will use this naming system instead of the default one. Leave as None for default!'),
    'resetMapping':             (False,False,'If set to False (off), the existing chain mapping will be used for export\n (otherwise it can be reset).'),
    'useCcpnChainInfo':         (False,False,'If set to True (on), the CCPN chain and sequence codes will be directly used for export.'),
    'tagIoCycle':               (None,False,'If set, all exported peaks will be marked with their export number.'),
    'forcePositiveVolumes':     (False,False,'If set to True (on), all peaks will be exported with positive volumes/intensities.'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.'),
    'presetResonanceMapping':   (False,False,'If set to True (on), will keep using first (or previously passed-in) Resonance:Atom mapping'),
  
  },

  'readPeakAssignments': {
    
    'version':                (None,False,'The format version number as string (if relevant, it will be used).'),
    'peakList':               (None,False,'A valid Nmr.PeakList object. The assignment information\n will be added to this list (and/or old information replaced).\n Use with caution!'),
    'addMode':                (False,False,'If set to True (on), will add assignment possibilities to the peaks in the peak list.'),
    'overwriteMode':          (False,False,'If set to True (on), will overwrite assignment possibilities\n for the peaks in the peak list.'),
    'useOriginalData':        (False,False,'If set to True (on), original format-specific data\n (from the peak list import) will be used.'),
    'resonanceMatching':      (True,False,'If set to False (off), the script will not try to match the newly\n created resonances to existing ones (in the Data Model).'),
    'minimalPrompts':         (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.'),
    'presetResonanceMapping':   (False,False,'If set to True (on), will keep using first (or previously passed-in) Resonance:Atom mapping'),
  
  },

  'writePeakAssignments': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'peakList':                 (None,True,'The Nmr.PeakList object to be exported to file.'),
    'dataDimRefs':              (None,False,'A list of Nmr.DataDimRef objects. The peak assignments\n will be exported using the order of these dataDimRefs.'),
    'dimMapping':               (None,False,'A list of the dataDimRef indexes that should\n be used for output (obsolete?).'),
    'useOriginalData':          (False,False,'If set to True (on), original format-specific data\n (from the import) will be used for export.'),
    'useOriginalResNames':      (False,False,'If set to True (on), original atom names\n (from the import) will be used for export.'),
    'noWrite':                  (False,False,'If set to True (on), no file will be written.'),
    'resetMapping':             (False,False,'If set to False (off), the existing chain mapping will be used for export\n (otherwise it can be reset).'),
    'useCcpnChainInfo':         (False,False,'If set to True (on), the CCPN chain and sequence codes will be directly used for export.'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.'),
    'presetResonanceMapping':   (False,False,'If set to True (on), will keep using first (or previously passed-in) Resonance:Atom mapping'),
  
  },

  'readConstraints': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'constraintType':           (None,True,'The type of constraint list to be created.'),
    'constraintList':           (None,False,'A valid subclass of a ConstraintList (e.g. DistanceConstraintList, ...). New information will be added to this list (and/or old information replaced).'),
    'constraintListName':       (None,False,'A name for the new constraint list, if one is created.'),
    'nmrConstraintStore':       (None,False,'The new constraint list will automatically be part of this NmrConstraintStore object.'),
    'strucGen':                 (None,False,'The new constraint list will automatically be linked to this StructureGeneration object.'),
    'resonanceMatching':        (True,False,'If set to False (off), the script will not try to match the newly created resonances to existing ones (in the Data Model).'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.'),
    'presetResonanceMapping':   (False,False,'If set to True (on), will keep using first (or previously passed-in) Resonance:Atom mapping'),
  
  },

  'readDistanceConstraints': {},  # Specific for distance constraints
  
  'readDihedralConstraints': {
  
    'newConstraintForDuplicate':  (False,False,'Set to True (on) if you want to create a new constraint when, for a set of resonances that already has a dihedral restraint, a new angle range is found (instead of adding as additional item for the existing one).'),
                              
  },  # Specific for dihedral constraints
  
  'readRdcConstraints': {},       # Specific for RDC constraints
  'readCsaConstraints': {},       # Specific for shift anisotropy constraints
  'readJCouplingConstraints': {}, # Specific for J Coupling constraints
  'readHBondConstraints': {},     # Specific for H-bond constraints

  'writeConstraints': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'constraintList':           (None,True,'The constraint list to be exported to file.'),
    'constraintType':           (None,True,'The type of constraint list to be exported.'),
    'exportAll':                (False,False,'If set to True (on), conversion warnings are exported\n with the constraints (currently NMR-STAR only).'),
    'sortForExport':            (False,False,'If set to True (on), the constraints will be exported sorted by the chain code and sequence ID of the resonances in the first constraint item (or the constraint itself for dihedrals).'),
    'useOriginalData':          (False,False,'If set to True (on), original format-specific data (from the import)\n will be used for export.'),
    'useOriginalResNames':      (False,False,'If set to True (on), original atom names (from the import) will be used for export.'),
    'compressResonances':       (True,False,'If set to False (off), the output atom names will not be compressed\n (e.g. HB2 and HB3 will be written out separately,\n even if their constraint items can be grouped (to HB*)).'),
    'retainOriginalCompression': (True,False,'If set to False (off), the original compression of atom names (based on whether resonances were copied or not) will not be retained.'),
    'individualAtoms':          (False,False,'If set to True (on), atom names will be fully written out (e.g. an Ala HB* is written in individual HB1, HB2, HB3 contributions)'),
    'individualAtomsIfNoSet':   (False,False,'If set to True (on), atom names for an atom set will be fully written out if no correct atomSysName found (e.g. an Ala HB* that has no matching name in the selected naming system is written in individual HB1, HB2, HB3 contributions)'),
    'ignoreChemCompsWithNoSysNames': (False,False,'If set to True (on), residues that do not have the required atom naming system information will not be exported. Default (False) is to use atom names from CCPN instead.'), 
    'noWrite':                  (False,False,'If set to True (on), no file will be written.'),
    'forceDefaultChainMapping': (False,False,'If set to True (on), the default chain mapping\n (between the format and the data model)\n will be automatically used.'),
    'forceNamingSystemName':    (None,False,'If set, will use this naming system instead of the default one. Leave as None for default!'),
    'resetMapping':             (False,False,'If set to False (off), the existing chain mapping will be used for export\n (otherwise it can be reset).'),
    'useCcpnChainInfo':         (False,False,'If set to True (on), the CCPN chain and sequence codes will be directly used for export.'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.'),
    'presetResonanceMapping':   (False,False,'If set to True (on), will keep using first (or previously passed-in) Resonance:Atom mapping'),
  
  },

  'writeDistanceConstraints': {},  # Specific for distance constraints
  'writeDihedralConstraints': {},  # Specific for dihedral constraints
  
  'writeRdcConstraints': {    # Specific for RDC constraints
    'negateValues':           (False,False,'If set to True (on), the RDC values will be negated in the export. This is to avoid problems with programs that use different sign conventions.'),          
  }, 
    
  'writeCsaConstraints': {},       # Specific for shift anisotropy constraints
  'writeJCouplingConstraints': {}, # Specific for J Coupling constraints
  'writeHBondConstraints': {},     # Specific for H-bond constraints

  'readAcqPars': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'experiment':               (None,False,'An Nmr.Experiment object to link the acquisition parameters to.'),
    'fidDataSource':            (None,False,'An Nmr.DataSource object of type FID to link the acquisition parameters to.'),
    'freqDataSource':           (None,False,'An Nmr.DataSource object of type processed to link the acquisition parameters to.'),
    'overwritePars':            (False,False,'Determines whether existing parameters are overwritten. Set to 1 (on) for overwriting.'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.')
  
  },
  
  'readProcPars': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'experiment':               (None,False,'An Nmr.Experiment object to link the processing parameters to.'),
    'expName':                  ('',False,'Name for the experiment to be created'),
    'aliasing':                 (None,False,'Aliasing range for a dimension. Is a dictionary with as key the dimension number,\n as value a tuple (minValue,maxValue).'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.')
  
  },

  'writeProcPars': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'dataDimRefs':              (None,True,'A list of Nmr.DataDimRef objects. The processing parameter information\n will be referenced using these dataDimRefs.'),
    'freqDataSource':           (None,False,'An Nmr.DataSource object of type processed to write the processing parameters from.'),
    'outputFile':               (None,False,'The file the processing script will write its output to.'),
    'noWrite':                  (False,False,'If set to True (on), no file will be written.'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.')
  
  },
  
  'writeProcScript': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'dataDimRefs':              (None,True,'A list of Nmr.DataDimRef objects. The processing parameter information\n will be referenced using these dataDimRefs.'),
    'freqDataSource':           (None,False,'An Nmr.DataSource object of type processed to write the processing parameters from.'),
    'inputFile':                (None,False,'The file the processing script will use as input.'),
    'outputFile':               (None,False,'The file the processing script will write its output to.'),
    'noWrite':                  (False,False,'If set to True (on), no file will be written.'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.')
  
  },
  

  'readPeopleAndCitations': {
    
    'useExistingPersons':       (True,False,'If set to False (off), will not try to match new persons to existing ones.'),
    'useExistingCitations':     (True,False,'If set to False (off), will not try to match new citations to existing ones.'),
    'doNotMakeCcpnObjects':     (False,False,'If set to True (on), will only make the FormatConverter objects and not the CCPN ones'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.')
  
  },

  'readAssignment': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'chain':                    (None,True,'The chain the assignment refers to.'), # TODO this might be specific to MARS, but keeping it for now
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.')
  
  },

  'readProject': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.'),
    'presetResonanceMapping':   (False,False,'If set to True (on), will keep using first (or previously passed-in) Resonance:Atom mapping'),
  
  },
  
  'writeProject': {
    
    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'useCcpnChainInfo':         (False,False,'If set to True (on), the CCPN chain and sequence codes will be directly used for export.'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.'),
    'presetResonanceMapping':   (False,False,'If set to True (on), will keep using first (or previously passed-in) Resonance:Atom mapping'),
  
  },

  #
  # Following are 'combined' read/writes
  #
  
  'readShiftsAndPeaks': {
  
    #
    # Limited options compared to full read for each!
    #

    'version':                  (None,False,'The format version number as string (if relevant, it will be used).'),
    'shiftFile':                (None,True,'Chemical shift file'),
    'measurementListName':      (None,False,'Name for the measurement list.'),
    'peakFile':                 (None,True,'Peak list file'),
    'dataSource':               (None,False,'A valid Nmr.DataSource object.\n The new peak list will be created using this object.'),
    'dataDimRefs':              (None,False,'A list of Nmr.DataDimRef objects.\n The peak list will be linked to these DataDimRefs.'),
    'specName':                 (None,False,'Name for the peak list.'),
    'resonanceMatching':        (True,False,'If set to False (off), the script will not try to match\n the newly created resonances to existing ones (in the Data Model).'),
    'minimalPrompts':           (False,False,'If set to True (on), will use automatic settings as much as possible.\n This will reduce the number of popups that appear but is less user-oriented.'),
    'presetResonanceMapping':   (False,False,'If set to True (on), will keep using first (or previously passed-in) Resonance:Atom mapping'),

  },
  
  #
  # Following are subclass definitions for handling component specific info
  #

  'subClasses': {'Measurements': ['HExchProtections',
                                  'HExchRates',
                                  'JCouplings',
                                  'Noes',
                                  'Rdcs',
                                  'ShiftAnistropy',
                                  'Shifts',
                                  'T1s',
                                  'T1Rhos',
                                  'T2s'],
                                  
                 'Constraints': ['DistanceConstraints',
                                 'DihedralConstraints',
                                 'CsaConstraints',
                                 'RdcConstraints',
                                 'JCouplingConstraints',
                                 'HBondConstraints']}

}
  
#
# Class DataFormat - generic functions
#


class DataFormat(TopShared):

  """
  

  DataFormat class: contains generic code for format import/export

  Some of the functions in here are redefined in the specific format class
  (e.g. XEasyFormat.py)

  
  """

  #
  # FormatConverter Error Class
  #
  
  class FormatConverterError(StandardError):
    
    def __init__(self, value):
      self.value = value
    def __str__(self):
      return repr(self.value)

  def __init__(self,project, guiParent = None, fileName = None, dataEntry = None, messageReporter = None, multiDialog = None, verbose = 0, **keywds):
  
    """
  
    Initialization function. Can also be used to read in project file on
    creation.
    
    Obligatory input:
    
    project:         a valid memops project
    
    
    Keyword input:
    
    guiParent:       usually a Tkinter.Tk() object, but can be None (will then use text-based
                     interaction)
    fileName:        file name for reading a project file
    dataEntry:       dataEntry class object (will be set up if None)
    messageReporter: messageReporter class object (will be set up if None)
    multiDialog:     multiDialog class object (will be set up if None)
    allowPopups:     if set to False, no popups will be displayed (no matter what guiParent is). Not fully functional (05/06/06)
    verbose:         if 1 prints out more messages
    **keywds:        other keywords - these will be passed on to self.readProject()
  
    """

    if project == None or type(project) != Implementation.MemopsRoot:
      raise self.FormatConverterError("Need MemopsRoot object.")

    #
    # Set up GUI or non-GUI user interaction if not passed in
    # TODO: bit nasty as reinitialized every time new format is created (but can create top one
    # and pass in if necessary)
    #
    
    if not dataEntry:
      dataEntry = setupDataEntry(guiParent)
      
    if not messageReporter:
      messageReporter = setupMessageReporter(guiParent)
      
    if not multiDialog:
      multiDialog = setupMultiDialog(guiParent)
        
    self.dataEntry = dataEntry
    self.messageReporter = messageReporter
    self.multiDialog = multiDialog

    #
    # Initialize other stuff
    #

    self.globals = globals()

    self.project = project
    self.guiParent = guiParent
    self.verbose = verbose
    self.status = True
    self.minimalPrompts = True
    
    # Added on 20/05/2009, keeps track of data imported. Can be used to assess whether import successful or not.
    # resetConvertCount resets the counter for whenever a read... function is called. Is set to False when reading
    # multiple data, should then be set back to True at end.
    self.convertCount = {}
    self.resetConvertCount = True
    self.molSystemCreated = False  # This is a hack to count externally created molSystems
    
    # Track newly created resonances. Note that this is a custom variable that needs to be reset
    # by external code (or not) depending on the purpose.
    self.newResonances = []

    self.warnings = {}
    
    self.IOkeywords = IOkeywords    
    self.setFormat()
    self.formatLabel = allFormatsDict[self.format]

    self.setFormatIOs()

    self.setGenericImports()
    
    #
    # Initialise dictionaries to speed things up
    #
    
    self.atomSysNamesByChemAtom = {}
    
    #
    # Hardset popup removal
    #
    
    self.allowPopups = True  # TODO: This is hacky - currently has to be set separately...
    if keywds.has_key('allowPopups'):
      self.allowPopups = keywds['allowPopups']
      del(keywds['allowPopups'])
    
    #
    # Set some keywords...
    #
    
    self.setAssignKeyword()

    #
    # Keeping track of all mapping...
    # Function sets self.nameMapping
    #

    self.getNameMapping()
    
    self.file = None

    if fileName:
      
      #
      # If a fileName is passed in at start (or later on with getFullProject)
      # then read project file and process elements
      #
      
      self.readProject(fileName, **keywds)
    
  def setFormatIOs(self):

    for componentInfo in IoSetupList + [('general',None),('chemShifts',None),('conversionPars',None),('hetNoe',None)]:
      componentType = componentInfo[0]    
      IOType = "%sIO" % componentType
      moduleName = 'ccp.format.%s.%s' % (self.format,IOType)
      #print moduleName,
      try:
        #print "  OK"
        formatIoModule = __import__(moduleName,{},{},[IOType])
        setattr(self,IOType,formatIoModule)
      except:
        #print "  NOGO"
        pass
          
  """

  Subclass independent functions

  """
  
  #
  # PROJECT (dummy defs to get format specific stuff)
  #

  def readProject(self,fileName=None,checkFileExists=True,**keywds):
  
    """
    Dummy def to get to format specific code
    """

    #
    # Set verbose (or not)
    #

    if 'verbose' in keywds.keys():
      self.verbose = keywds.pop('verbose')

    #
    # Get info for full project
    #
    
    if not self.setIOkeywords('readProject', keywds, verbose = False):
      return
  
    if checkFileExists and not self.fileExists(fileName):
      return None

    # Set to False so can collect all info
    self.resetConvertCount = False
    
    try:
      self.getFullProject(fileName,**keywds)
    except IOError, e:
      self.messageReporter.showError("Error","Error reading %s project file: %s" % (self.format,str(e)),parent = self.guiParent)
      return None
      
    self.nameMapping.isOriginalImport = True
    
    # Set back to True so will be reset next time
    self.resetConvertCount = True
    
    self.printConvertCountString(self.verbose)

    #
    # TODO this is rather nasty - should always return an entry?
    #
    
    if hasattr(self,'entry') and self.entry:
      returnValue = self.entry
    else:
      returnValue = self.project
      
    return returnValue
  
  def writeProject(self,fileName = None,**keywds):
  
    """
    
    Goes on to format specific code
    
    """
    
    if not self.setIOkeywords('writeProject', keywds, verbose = False):
      return

    returnValue = self.createFullProject(fileName,**keywds)
    
    return returnValue

  #
  # SEQUENCE
  #
  
  def readSequence(self,fileName = None, verbose = True, **keywds):
  
    """
    
    Generic code for reading and importing sequence information    
    
    Mandatory argument:
    
    fileName:               file name to read the sequence from
    
    
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    # Note that most of these only work on file with ONE sequence (and not on
    # a typical PDB file with ligands/dimers/... - then the info comes from
    # the header)
    #
    
    if not self.setIOkeywords('readSequence', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
    
    self.chains = []
    
    self.nameMapping.isOriginalImport = True
    
    self.setSequenceFileClass()
    
    self.initConvertCount('sequence')

    #
    # Read format file. Note that fileName not necessary if self.file exists (from getProject())
    #

    if fileName:
      
      #
      # If fileName passed in, read the file...
      #
      
      if not self.fileExists(fileName):
        return None

      self.fileName = fileName
      try:
        self.getSequence()
      except IOError, e:
        self.messageReporter.showError("Error","Error reading %s sequence file: %s" % (self.format,str(e)),parent = self.guiParent)
        return None

    elif not self.file:
      self.messageReporter.showError("Error","Need a valid filename to read sequence.",parent = self.guiParent)
      return None
  
    #
    # If there's no reference chemComp data, don't bother
    #
    
    self.setChemCompRepositories()
        
    #
    # Get relevant naming system (necessary for ChemCompSysName and possible ChemAtomSysName searches)
    # 
    
    self.namingSystemName = self.getFormatNamingSystemName()

    #
    # Start processing - first create molecules
    # and keep track of which bits of the self.sequence they belong to...
    #
    
    (createMoleculeDict,molSystemBonds) = self.checkMolecules()

    #
    # Set molSystem
    #

    if not self.molSystem:
      self.checkMolSystem()

    #
    # Then create the chains
    #
   
    self.checkChains(createMoleculeDict,molSystemBonds)

    #
    # And do the spin systems (molecule independent...)
    #
    
    for self.sequence in self.sequenceFile.sequences:
    
      self.checkSpinSystems()
    
    #
    # Print out import info, if necessary (can also be done later on)
    #
    
    self.printConvertCountString(verbose)
    
    return self.chains
    
        
  def writeSequence(self, fileName = None,verbose = True, **keywds):
    
    """
    
    Generic code for exporting and writing sequence information
    
    Mandatory argument:
    
    fileName:               file name to write the sequence to

    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('writeSequence', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
      
    self.fileName = fileName
    self.setSequenceFileClass()
    
    #
    # Assuming that correct chains are always passed in
    #
    
    if not isinstance(self.chains, (list, tuple)):
      
      self.messageReporter.showError("Error","Need valid list of chain(s) to write sequence.",parent = self.guiParent)
      return
    
    else:
    
      for chain in self.chains:
 
        if not isinstance(chain,MolSystem.Chain):
          
          self.messageReporter.showError("Error","Need valid chain object(s) to write sequence. Removed invalid %s" % chain.__class__,parent = self.guiParent)
          self.chains.pop(self.chains.index(chain))
                    
    #
    # Can write multiple chains for formats that support it
    #

    if len(self.chains) == 0:
        
      self.messageReporter.showError("Error","No valid chains passed in.",parent = self.guiParent)
      return
    
    #
    # Get relevant naming system
    # 
    
    self.namingSystemName = self.getFormatNamingSystemName()

    #
    # Set/get the mapping for the chains
    #
    
    self.getChainMapping(resetMapping = self.resetMapping)
    
    #
    # Sort according to exportChainCode, seqCode (not always necessary but doesn't hurt)
    #
    
    self.getSortedChainList()
    
    #
    # Create the format sequence file object
    #
    
    self.createSequenceFile() 
   
    for self.chain in self.chainList:

      self.exportChainCode = self.chainDict[self.chain][0]
      self.seqCode = self.getExportSeqCode(self.chainDict[self.chain][1],self.chain.sortedResidues()[0])
    
      self.createSequence()

      #
      # Get residues and put info in SequenceFileElements
      # Use mapping...
      #
    
      for self.residue in self.chain.sortedResidues():
    
        self.seqCode = self.getExportSeqCode(self.chainDict[self.chain][1],self.residue)
        self.setSequenceFileElements()        
          
    #
    # Write file
    #

    if self.noWrite == False:

      self.sequenceFile.write()
      del self.sequenceFile

    return True
  
  #
  # Coordinates. Allows multiple file input (should all be for the same chains though!)
  #
  
  def readCoordinates(self, fileNames = None,verbose = True, **keywds):

    """
    
    Generic code for reading coordinate files
         
    Mandatory argument:

    fileNames:          file names to read the coordinates from
                        (Warning: this assumes that all files have exactly the same
                         atom naming)
    
    Keyword definitions in IOkeywords. 
        
    """
    
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('readCoordinates', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #

    self.fileNames = fileNames
    self.structures = []    
    self.structureNumber = 0
    self.atomNameDict = {}    # This is used to speed up atom matching
    
    self.setCoordinateFileClass()
    
    self.setChemCompRepositories()

    self.initConvertCount('coordinates')

    #
    # Set single fileName to list...
    # 
    
    if isinstance(self.fileNames, str):
      self.fileNames = [self.fileNames]
      
    #
    # Loop over all fileNames
    #
    
    self.modelSet = None
    self.coordResidues = {}
    self.specialSeqCodeMapping = False  # This is set to True for NMR-STAR imports where the sequence codes are reset - for NRG project

    for self.fileName in self.fileNames:
      
      if self.fileName:
        if not self.fileExists(self.fileName):
          continue

        #
        # Get the coordinates file
        #

        try:
          self.getCoordinates()
        except IOError, e:
          self.messageReporter.showError("Error","Error reading %s coordinates file: %s" % (self.format,str(e)),parent = self.guiParent)
          return None
          
      elif not self.file:
        self.messageReporter.showError("Error","Need a valid filename to read coordinates.",parent = self.guiParent)
        return None   
        
      #
      # Only do below when first file is read (assuming all to be one group!)
      #
      
      if not self.structureNumber:

        #
        # Select molSystem!
        #

        if not self.molSystem:

          self.selectMolSystem()


        #
        # Set up the ccp and external format chain mapping...
        #

        if not self.getChainLinks():
          print "  Error: no complete chain linking found between coordinates and chains. Aborting coordinate reading."
          return False

        #
        # Get relevant naming system
        #
        
        self.namingSystemName = self.getFormatNamingSystemName()

        #
        # Use automatic naming system search. Don't run if forceNamingSystemName given.
        # If  pre-defined naming system for the format exists, still do the search, but only select
        # a different naming system if there is no 100% match.
        #

        if not self.forceNamingSystemName or self.format in ('pseudoPdb',):
          self.autoSelectCoordNamingSystem()


        #
        # (Re)set the name mapping with this naming system
        #

        self.getNameMapping(namingSystemName = self.namingSystemName,  isOriginalImport = True)

        #
        # Select structureEnsemble!
        #

        if not self.structureEnsemble:

          self.selectStructureEnsemble()
      
        #
        # Do coordinate atom mapping only once for an ensemble...
        #
    
        self.mapCoordinateAtoms()

      #
      # Molsystem and chain/residue mapping are now set... do the coordinates
      #     
      
      self.setCoordinates()

    #
    # Print out import info, if necessary (can also be done later on)
    #
    
    self.printConvertCountString(verbose)
  
    self.specialSeqCodeMapping = False  # Reset, just to be sure.
        
    return self.structures

  def writeCoordinates(self, fileName = None, verbose = True, **keywds):
    
    """
    
    Generic code for exporting and writing coordinate information
      
    Mandatory argument:
    
    fileName:               file name to write the structure coordinates to
    
    
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('writeCoordinates', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
        
    self.fileName = fileName
    self.writeKeywds = {}

    self.setCoordinateFileClass()
    
    #
    # Do some checks...
    #
    
    if self.structures:

      if not isinstance(self.structures, (list, tuple)):

        self.messageReporter.showError("Error","Need valid list of structures to write coordinates.",parent = self.guiParent)
        return

      else:

        for structure in self.structures:

          if not isinstance(structure,MolStructure.Model):

            self.messageReporter.showError("Error","Need valid structure object(s) to write sequence. Removed invalid %s" % structure.__class__,parent = self.guiParent)
            self.structures.pop(self.structures.index(structure))    
    else:
    
      self.messageReporter.showError("Error","Need valid list of structures to write coordinates.",parent = self.guiParent)
      return
             
    #
    # Check if anything left...
    #

    if len(self.structures) == 0:
        
      self.messageReporter.showError("Error","No valid structures passed in.",parent = self.guiParent)
      return
       
    #
    # Get all the chains listed for these structures...
    #
    
    self.getCoordinateChains()
    
    #
    # Get/set the mapping for the chains
    #
    
    self.getChainMapping(resetMapping = self.resetMapping)
    
    #
    # Sort according to exportChainCode, seqCode (not always necessary but doesn't hurt)
    #
    
    self.getSortedChainList()
    
    #
    # Get relevant naming system
    #

    self.namingSystemName = self.getFormatNamingSystemName()

    #
    # Create the format coordinate file object
    #
    
    self.createCoordinateFile() 
        
    #
    # Set the coordinates...
    #
    
    self.createCoordinates()
              
    #
    # Write file
    #

    if self.noWrite == False:

      self.coordinateFile.write(**self.writeKeywds)
      del self.coordinateFile
      
    return True
  
  #
  # ChemComp creation
  #
  
  def readChemComps(self, fileName = None,verbose = True, **keywds):

    """
    
    Generic code for reading chemComp information files
         
    Mandatory argument:

    fileName:          file name to read the chemComp information from
    
    Keyword definitions in IOkeywords. 
        
    """
    
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('readChemComps', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
 
    self.chemComps = []       
    
    self.setChemCompFileClass()

    self.initConvertCount('chemComp')
    
    #
    # Get the chemComp info
    #

    if fileName:
    
      if not self.fileExists(fileName):
        return None
    
      self.fileName = fileName

      try:
        self.getChemComps()
      except IOError, e:
        self.messageReporter.showError("Error","Error reading %s chemical compound file: %s" % (self.format,str(e)),parent = self.guiParent)
        return None
    
    #
    # Set info for chemComp(Coord)s
    #
    
    self.setChemCompRepositories()    
    self.setChemCompCoordRepositories()
        
    #
    # If multiple codes available, select which ones to use
    #
    
    self.rawChemCompList = self.selectRawChemComps()
    
    #
    # Loop over information per original chemComp...
    #
    
    for self.rawChemComp in self.rawChemCompList:
    
      self.setupChemCompInfo()      
      self.makeChemComp()

    #
    # Print out import info, if necessary (can also be done later on)
    #
    
    self.printConvertCountString(verbose)
          
    return self.chemComps    
        
  def writeChemComp(self, fileName = None, verbose = True, **keywds):
    
    """
    
    Generic code for exporting and writing chemComp information
      
    Mandatory argument:
    
    fileName:               file name to write the chemCompVar to
    
    
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('writeChemComp', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
        
    self.fileName = fileName
    self.setChemCompFileClass()
    
    self.setChemCompCoordRepositories()
    
    #
    # Do some checks...
    #
    
    #
    #
    # TODO TODO NEED SELECTOR FOR CHEMCOMPVAR, NAMING SYSTEM AND COORDINATE SYSTEM HERE!!!
    #
    #
    
    if self.chemCompVar:

      if not isinstance(self.chemCompVar,ChemComp.ChemCompVar):

            self.messageReporter.showError("Error","Need valid chemCompVar object to write chemComp file (instead of a %s)." % structure.__class__,parent = self.guiParent)

    else:
    
      self.messageReporter.showError("Error","Need valid chemCompVar to write chemComp file.",parent = self.guiParent)
      return
        
    #
    # Set relevant naming system
    #

    self.namingSystemName = self.getFormatNamingSystemName()

    #
    # Create the format chemComp file object
    #
    
    self.createChemCompFile() 
        
    #
    # Get the coordinates and write create the raw chemComp.
    #
    
    self.getChemCompCoord()
    
    if not self.chemCompVarCoord:
      return False
    
    self.createRawChemComp()
              
    #
    # Write file
    #

    if self.noWrite == False:

      self.chemCompFile.write()
      del self.chemCompFile
      
    return True

  #
  # Measurements
  #
  
  def setMeasurementCodes(self,readWrite, keywds, verbose = True):
    
    """
    Defines the name for the constraint list in the Api
    """
  
    if self.measurementType == 'shift':
      self.measurementApiCode = 'Shift'
      self.measurementResNum = 1
      self.formatFileBaseName = 'chemShifts'  # This is the name before the IO bit, e.g. chemShiftsIO in this case
      self.formatFileValueKey = 'chemShifts' # This is the name of the value list in the raw format measurement file, e.g. self.chemShifts
      self.formatFileClassName = 'ChemShift' # This is the name of the File class in the raw format measurement file

    elif self.measurementType[-5:] == 'Relax':
      self.measurementResNum = 1
      self.formatFileBaseName = self.measurementType
      self.formatFileValueKey = "%sValues" % self.measurementType
      self.formatFileClassName = self.measurementType[0].upper() + self.measurementType[1:]
      self.measurementApiCode = self.formatFileClassName[:-5]

    elif self.measurementType[:5] == 'hExch':
      # Todo could combine with above...
      self.measurementResNum = 1
      self.formatFileBaseName = self.measurementType
      self.formatFileValueKey = "%sValues" % self.measurementType
      self.formatFileClassName = self.measurementType[0].upper() + self.measurementType[1:]
      self.measurementApiCode = self.formatFileClassName

    elif self.measurementType in ('shiftAnisotropy',):
      self.measurementResNum = 1
      self.formatFileBaseName = 'csa'
      self.formatFileValueKey = 'csaValues'
      self.formatFileClassName = 'Csa'
      self.measurementApiCode = 'ShiftAnisotropy'
    
    elif self.measurementType in ('noe',):
      self.measurementResNum = 2
      self.formatFileBaseName = 'hetNoe'
      self.formatFileValueKey = 'hetNoeValues'
      self.formatFileClassName = 'HetNoe'
      self.measurementApiCode = 'Noe'
      
    elif self.measurementType in ('jCoupling','rdc'):
      self.measurementResNum = 2
      self.formatFileBaseName = self.measurementType
      self.formatFileValueKey = "%sValues" % self.measurementType
      self.formatFileClassName = self.measurementType[0].upper() + self.measurementType[1:]
      self.measurementApiCode = self.formatFileClassName
             
    else:
      self.messageReporter.showError("Unrecognized type","Unrecognized measurement type %s." % str(self.measurementType))
      return False

    if not self.setIOkeywords(readWrite + self.measurementApiCode + 's', keywds, verbose = verbose):
      return False
  
    return True

  def readChemShifts(self,*args,**keywds):
  
    # Backward compatibility fix

    return self.readShifts(*args, **keywds)

  def readHExchProtections(self, *args, **keywds):
    
    """
    Shell function for calling readMeasurements
    """
    
    keywds['measurementType'] = 'hExchProtection'
    
    if not keywds.has_key('unit'):
      keywds['unit'] = 'ratio'

    self.setIOkeywords('readHExchProtections', keywds, verbose = False)
    
    return self.readMeasurements(*args, **keywds)

  def readHExchRates(self, *args, **keywds):
    
    """
    Shell function for calling readMeasurements
    """
    
    keywds['measurementType'] = 'hExchRate'
    
    if not keywds.has_key('unit'):
      # Other recommended unit is min-1
      keywds['unit'] = 's-1'

    self.setIOkeywords('readHExchRates', keywds, verbose = False)

    return self.readMeasurements(*args, **keywds)

  def readShiftAnistropy(self, *args, **keywds):
    
    """
    Shell function for calling readMeasurements
    """
    
    keywds['measurementType'] = 'shiftAnisotropy'
    
    if not keywds.has_key('unit'):
      keywds['unit'] = 'ppm'

    self.setIOkeywords('readShiftAnistropy', keywds, verbose = False)

    return self.readMeasurements(*args, **keywds)

  def readShifts(self, *args, **keywds):
    
    """
    Shell function for calling readMeasurements
    """
    
    keywds['measurementType'] = 'shift'
    
    if not keywds.has_key('unit'):
      keywds['unit'] = 'ppm'

    self.setIOkeywords('readShifts', keywds, verbose = False)

    return self.readMeasurements(*args, **keywds)

  def readT1s(self, *args, **keywds):
    
    """
    Shell function for calling readMeasurements
    """
    
    keywds['measurementType'] = 't1Relax'

    if not keywds.has_key('unit'):
      keywds['unit'] = 's'

    self.setIOkeywords('readT1s', keywds, verbose = False)

    return self.readMeasurements(*args, **keywds)

  def readT1Rhos(self, *args, **keywds):
    
    """
    Shell function for calling readMeasurements
    """
    
    keywds['measurementType'] = 't1RhoRelax'

    if not keywds.has_key('unit'):
      keywds['unit'] = 's'

    self.setIOkeywords('readT1Rhos', keywds, verbose = False)

    return self.readMeasurements(*args, **keywds)

  def readT2s(self, *args, **keywds):
    
    """
    Shell function for calling readMeasurements
    """
    
    keywds['measurementType'] = 't2Relax'

    if not keywds.has_key('unit'):
      keywds['unit'] = 's'

    self.setIOkeywords('readT2s', keywds, verbose = False)

    return self.readMeasurements(*args, **keywds)

  def readNoes(self, *args, **keywds):
    
    """
    Shell function for calling readMeasurements
    """
    
    keywds['measurementType'] = 'noe'

    self.setIOkeywords('readNoes', keywds, verbose = False)

    return self.readMeasurements(*args, **keywds)

  def readJCouplings(self, *args, **keywds):
    
    """
    Shell function for calling readMeasurements
    """
    
    keywds['measurementType'] = 'jCoupling'

    self.setIOkeywords('readJCouplings', keywds, verbose = False)

    return self.readMeasurements(*args, **keywds)

  def readRdcs(self, *args, **keywds):
    
    """
    Shell function for calling readMeasurements
    """
    
    keywds['measurementType'] = 'rdc'

    self.setIOkeywords('readRdcs', keywds, verbose = False)

    return self.readMeasurements(*args, **keywds)

  #
  # Generic measurement reader...
  #

  def readMeasurements(self,fileName = None, verbose = True, **keywds):
  
    """
    
    Generic code for reading and importing measurement list information    
    
    Mandatory argument:
    
    fileName:               file name to read the measurements from
    
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords (silently because measurement list specific keywds might be available)
    #
    
    if not self.setIOkeywords('readMeasurements', keywds, lockUpdate = True, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #

    setCurrentStore(self.project,'NmrProject')

    if not self.setMeasurementCodes('read', keywds, verbose = verbose):
      return False
    
    self.nameMapping.isOriginalImport = True
    
    self.setMeasurementFileClass()        

    self.initConvertCount('measurements')

    #
    # Read format file. Note that fileName not necessary if self.file exists (from getProject())
    #
    
    if fileName:
    
      #
      # Read file...
      #

      if not self.fileExists(fileName):
        return None
      
      self.fileName = fileName

      self.getMeasurements()
      
      if not self.measurementFile:
        return None
      
    elif not self.file:
      self.messageReporter.showError("Error","Need a valid filename to read %ss." % self.measurementType,parent = self.guiParent)
      return
    
    self.measurementFileValues = getattr(self.measurementFile,self.formatFileValueKey)
       
    #
    # Some further initialization
    #        

    self.ApiMeasurementList = getattr(Nmr,"%sList" % self.measurementApiCode)

    self.ApiMeasurement = getattr(Nmr,"%s" % self.measurementApiCode)
    
    self.setMeasurementListSpecific = getattr(self,'set%sListSpecific' % self.measurementApiCode)
    self.setMeasurements = getattr(self,"set%ss" % self.measurementApiCode)
    self.setMeasurementSpecific = getattr(self,"set%sSpecific" % self.measurementApiCode)
    self.createMeasurement = getattr(self,"create%s" % self.measurementApiCode)

    #
    # Check some values and types
    #
    
    if not self.measurementFile or len(self.measurementFileValues) == 0:
      self.messageReporter.showError("No %ss" % self.measurementType,"No %s values found in file %s: not used" % (self.measurementType,fileName),self.guiParent)
      return None
      
    if self.measurementList:
      if not isinstance(self.measurementList,self.ApiMeasurementList):
        self.messageReporter.showWarning("Invalid %ss" % self.measurementType,"Invalid %s list passed in (%s). Creating new one." % (self.measurementType,self.measurementList.__class__),self.guiParent)
        self.measurementList = None
    
    if self.experiments:

      if not isinstance(self.experiments, (list, tuple)):
        self.messageReporter.showWarning("Invalid experiment list","No list of experiments passed in: are ignored.",self.guiParent)
        self.experiments = None
      
      else:
        for experiment in self.experiments:
          if not isinstance(experiment,Nmr.Experiment):
            self.messageReporter.showWarning("Invalid experiment","Invalid experiment object passed in (%s). Ignored." % experiment.__class__,self.guiParent)
            self.experiments.pop(self.experiments.index(experiment))

    if self.method:

      if not isinstance(self.method, (list, tuple)):
        self.messageReporter.showWarning("Invalid method list","No list of method passed in: are ignored.",self.guiParent)
        self.method = None
      
      else:
        for method in self.method:
          if not isinstance(method,Method.Method):
            self.messageReporter.showWarning("Invalid method","Invalid method object passed in (%s). Ignored." % method.__class__,self.guiParent)
            self.method.pop(self.method.index(method))

    #
    # Check whether measurements already exist in a list (return if so)
    #
    # With respect to resonances:
    #
    # 'resonanceMatching' is on (has a non-zero value):
    #   - Create new resonances & measurements IF they are not found as a direct
    #     name match in the applicationData (otherwise use the exising one)
    #
    # 'resonanceMatching' is off (is None or zero):
    #   - Always create a new resonance/measurement (no checking on exising ones)
    #
    # linkResonance should always be run at one point to clean up the
    # resonances and assign them to atoms
    #
    # Warning: the format specific names are stored as application data!
    # The names of the newly created resonances will be empty
    #

    self.checkMeasurementsAndResonances()

    if self.existingMeasurementList:
      
      #
      # Exactly same list: no user intervention necessary. Just use this one.
      #
    
      return self.existingMeasurementList

    #
    # Create a new measurement list if none given
    #
        
    self.setMeasurementList()
    
    #
    # Now ready to start using measurement...
    #
    
    self.setMeasurements()

    if not self.measurementList.details:
      self.measurementList.details = "Origin %s file %s" % (self.format,self.measurementFile.name)

    else:
      self.measurementList.details += "\nOrigin %s file %s" % (self.format,self.measurementFile.name)

    #
    # Reset in case changed
    #

    self.setAssignKeyword()

    #
    # Print out import info, if necessary (can also be done later on)
    #
    
    self.printConvertCountString(verbose)
    
    #
    # Measurement list is returned - the new resonances should be linked
    # to atoms or to existing resonances with the 'linkResonances' script
    #
     
    return self.measurementList

    
  def writeChemShifts(self,*args,**keywds):
    
    # Backward compatibility fix
  
    return self.writeShifts(*args, **keywds)

  def writeShifts(self, *args, **keywds):
    
    """
    Shell function for calling writeMeasurements
    """
  
    keywds['measurementType'] = 'shift'
   
    self.setIOkeywords('writeShifts', keywds, verbose = self.verbose)
    
    return self.writeMeasurements(*args, **keywds)

  def writeMeasurements(self, fileName = None, verbose = True, **keywds):

    """
    
    Generic code for exporting and writing measurement information to files in a specific format
    
    Mandatory argument:
    
    fileName:               file name to write the measurements to
    
    
    Keyword definitions in IOkeywords. 

    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('writeMeasurements', keywds, lockUpdate = True, verbose = verbose):
      return False
            
    #
    # Initialize other variables...
    #

    if not self.setMeasurementCodes('write', keywds, verbose = verbose):
      return False

    self.fileName = fileName
    self.writeKeywds = {}
    
    self.setMeasurementsExportMode()  # Set the export mode, either 'atoms' or 'spinSystems', depending on format.

    #
    # Initial check
    #

    if not self.measurementList or not isinstance(self.measurementList,getattr(Nmr,"%sList" % self.measurementApiCode)):
      self.messageReporter.showError("Error","No or invalid %s list provided for writing %s: file not written." % (self.measurementType,fileName),parent = self.guiParent)
      return
    
    #
    # Further initialization
    #
    
    self.formatCode = "%s%s" % (self.formatLabel,self.formatFileClassName)
    self.initMeasurementFileResidue = getattr(self,"init%sFileResidue" % self.formatFileClassName)
    self.setMeasurementFileValue = getattr(self,"set%sFileValue" % self.formatFileClassName)

    #
    # If chains were passed in, check if they are valid
    #
    
    if hasattr(self,'chains') and self.chains:
    
      if not isinstance(self.chains, (list, tuple)):
        self.messageReporter.showWarning("Invalid chain list","No list of chains passed in: are ignored.",self.guiParent)
        self.chains = None
      
      else:
        for chain in self.chains:
          if not isinstance(chain,MolSystem.Chain):
            self.messageReporter.showWarning("Invalid chain","Invalid chain object passed in (%s). Ignored." % chain.__class__,self.guiParent)
            self.chains.pop(self.chains.index(chain))
        
        self.molSystem = self.chains[0].molSystem
    
    #
    # Check whether there are unlinked resonances, ask to run linkResonances
    #
    
    self.checkUnlinkedResonances(exportMode = self.exportMode)
    
    #
    # Routine to link measurements to resNames
    #
    
    self.createAtomMeasurements()            

    #
    # Allow chain selection
    #
    
    self.selectChains()

    #
    # Set the mapping for the chains
    #
    
    self.getChainMapping(resetMapping = self.resetMapping)

    #
    # Create external format measurement file
    #
    
    self.setMeasurementFileClass()
    
    self.createMeasurementFile()    

    self.rawMeasurementClass = getattr(self.measurementsIO,self.formatCode)
    self.measurementFileValues = getattr(self.measurementFile,self.formatFileValueKey)

    #
    # Loop over chains and set measurements and other values in measurementFile
    #
    
    self.setMeasurementFileValues()
    
    #
    # Check if atomMeasurements dict empty
    #

    self.checkAtomMeasurements()
    
    #
    # Write file
    #
    
    if self.noWrite == False:
      self.measurementFile.write(**self.writeKeywds) 
      del self.measurementFile

    return True
  
  #
  # PEAKS
  #

  def readPeaks(self,fileName = None, verbose = True, **keywds):
  
    """
    
    Generic code for reading and importing peak list information
    
    Will create new peaklist *except* when peaklist passed in: will then overwrite
    peaks with same peak number, add peaks that are not present, and will delete
    peaks that are not present in the new list. This is *not* the same as creating
    a new list: all links from the original peaklist will be maintained!
    
    Warning: passing in peakList only works for SINGLE peakList: for Ansig
    have to read in each one in the export file separately! Also all peaks 
    have to be linked to same dataDimRefs (or they will be deleted).
      
    Mandatory argument:
    
    fileName:               file name to read the peaks from
    
    
    Keyword definitions in IOkeywords. 

    """    
    
    # 
    # Note: not trivial to make all single inputs above into lists... and is it
    # necessary? Really only valid for Ansig format...
    #
         
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('readPeaks', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
    
    self.setPeakDimOrder()
        
    if not self.dataDimRefs:
      self.dataDimRefs = []
      
    self.peakLists = []
    self.origPeakNumberDict = {}
    
    self.nameMapping.isOriginalImport = True
    
    self.setPeakFileClass()

    self.initConvertCount('peaks')

    setCurrentStore(self.project,'NmrProject')

    if self.peakList:
    
      #
      # If a peaklist was passed in get the values from there
      # if possible.
      #
    
      self.overwrite = 1
      self.dataSource = self.peakList.dataSource
      self.existingPeaks = list(self.peakList.sortedPeaks())
      
      if not self.specName:
        self.specName = self.peakList.name
      else:
        self.peakList.name = self.specName
            
      for peak in self.peakList.sortedPeaks():
        applData = peak.findFirstApplicationData(application = self.format, keyword = peakNum_kw)
        
        if applData:
          self.origPeakNumberDict[applData.value] = peak

    else:
    
      self.overwrite = 0
      self.existingPeaks = []
     
    #
    # Note that fileName not necessary if self.file exists (from getProject())
    #
    
    if fileName:
    
      #
      # Read file...
      #      

      if not self.fileExists(fileName):
        return None

      self.fileName = fileName
      
      try:
        self.getPeaks()
      except IOError, e:
        self.messageReporter.showError("Error","Error reading %s peak list file: %s" % (self.format,str(e)),parent = self.guiParent)
        return None
      
    elif not self.file:
      self.messageReporter.showError("Error","Need a valid filename to read peaks.",parent = self.guiParent)
      return

    #
    # Reference software definitions data
    #

    hasSoftwareDefinitions = False
    
    if self.project.currentMethodStore:
      if self.project.currentMethodStore.findFirstSoftware(name=self.format):
        hasSoftwareDefinitions = True
                  
    if not hasSoftwareDefinitions:
    
      getSoftwareDefinitions(self.project)
      print "Read software definitions."

    #
    # Do some checks on the objects that were passed in
    #
    
    if self.dataSource and not isinstance(self.dataSource,Nmr.DataSource):
    
      self.messageReporter.showError("Error","Need valid dataSource object to read peaks.",parent = self.guiParent)
      return      
    
    if not isinstance(self.dataDimRefs, (list, tuple)):
      
      self.messageReporter.showError("Error","Need valid list of dataDimRefs to read peaks.",parent = self.guiParent)
      return
    
    else:
    
      for dataDimRef in self.dataDimRefs:
      
        if not isinstance(dataDimRef,Nmr.DataDimRef):
          
          self.messageReporter.showError("Error","Need valid dataDimRef object(s) to read peaks.",parent = self.guiParent)
          return
          
    #
    # Create software dicts
    #
    
    self.createSoftwareDict()
      
    #
    # Make dict of resonance names (based on application format if available)
    # Warning: for XEasy this will actually use atomSerials!
    #
    
    self.createResonanceNamesDict(self.project.currentNmrProject.resonances)
    
    #
    # Setup for handling number of spectra & which one(s) to use
    #    

    if not self.setNumberPeakLists():
      return
    
    #
    # Loop over all peaklists that were read in.
    #
    # Note that one format peak list object could have peak lists linked to
    # multiple spectra (e.g. ansig)
    #
    
    for specIndex in range(0,self.numPeakLists):
      
      if not self.specName:
        self.specName = self.peakFile.specNames[specIndex]
        
      self.numDim = self.peakFile.numDims[specIndex]
      
      #
      # Check if peak list has any peaks...
      #
      
      if self.checkPeakListValid():

        #    
        # TODO: use swh, o1, ... for spectrum setup and validation IF available!
        #

        #
        # Get datasource and experiment info (create new one if necessary)
        # 

        self.setDataSource()

        # 
        # Set info for new peaklist
        #

        if not self.peakList:

          peakListName = None

          if not self.minimalPrompts:
            peakListName = self.dataEntry.askString("Enter peak list name","New peak list name",self.specName,self.guiParent)

          if not peakListName:
            peakListName = self.specName

          self.createPeakList(peakListName)
          
        #
        # Format specific things (e.g. Felix format column info selection)
        #
        
        self.setRawPeakFileFormatSpecific()

        #
        #  Select dataDimRefs if not given
        #

        self.selectDataDimRefs()

        #
        #  Edit experiment information (if required)
        #

        self.editDataSourceInformation()

        #
        # Set mapping information for peaklist
        #

        self.setPeakListDimMapping()

        #
        # Set the peaks in the data model and make links to resonances
        #

        self.setPeaks()

        #
        # Keep track of peakLists created for returning
        #

        self.peakLists.append(self.peakList)

        #
        # Delete original peaks left over if overwriting
        #

        self.deleteRemainingPeaks()
      
      elif self.verbose:
      
        print "  Warning: peak list %s doesn not contain any valid peaks - skipped." % self.specName


      #
      # Reset values for next item in loop
      #

      self.peakList = None
      self.specName = None
      self.dataSource = None
      self.dataDimRefs = []


    #
    # Reset in case changed
    #

    self.setAssignKeyword()

    #
    # Print out import info, if necessary (can also be done later on)
    #
    
    self.printConvertCountString(verbose)

    return self.peakLists
    
  def writePeaks(self, fileName = None, verbose = True, **keywds):
        
    """
    
    Generic code for writing and exporting peak lists    
    
    Mandatory argument:
    
    fileName:               file name to write the peak list to
    
    
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('writePeaks', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
   
    self.fileName = fileName
    self.writeKeywds = {}
    
    self.setPeakDimOrder()
    
    self.setPeakFileClass()

    #
    # Hack for passing in single peak lists
    #
    
    if isinstance(self.peakLists,Nmr.PeakList):
      self.peakLists = [self.peakLists]
    
    #
    # Initial check
    #

    if not isinstance(self.peakLists, (list, tuple)):
      self.messageReporter.showError("Error","Invalid list of peaks provided for writing %s: file not written." % fileName,parent = self.guiParent)
      return
      
    else:
      for peakList in self.peakLists:
        if not peakList or not isinstance(peakList,Nmr.PeakList):
          self.messageReporter.showError("Error","Invalid peak list provided for writing %s: deleted." % fileName,parent = self.guiParent)
          self.peakLists.pop(self.peakLists.index(peakList))
          
      if len(self.peakLists) == 0:
        self.messageReporter.showError("Error","No valid peak lists given - aborting." % fileName,parent = self.guiParent)
        return
            
    #
    # Additional check - for some formats (e.g. XEasy) a shift list has to be written
    # out before the assignments can be transferred...
    #

    if not self.peakListAssignmentCheck(self.peakLists):
    
      return

    #
    # If dataDimRefs were passed in, check if they are valid
    # This only works if one peaklist passed in
    #
    
    if self.dataDimRefs and len(self.peakLists) == 1:
    
      if not isinstance(self.dataDimRefs, (list, tuple)):
        self.messageReporter.showWarning("Invalid dataDimRef list","No list of dataDimRefs passed in: are ignored.",self.guiParent)
        self.dataDimRefs = None
      
      else:

        self.dataSource = self.peakLists[0].dataSource
        
        for dataDimRef in self.dataDimRefs:

          if not isinstance(dataDimRef,Nmr.DataDimRef) or dataDimRef.dataDim.dataSource != self.dataSource:
            self.messageReporter.showWarning("Invalid dataDimRef","Invalid dataDimRef object passed in (%s). All ignored." % dataDimRef.__class__,self.guiParent)
            self.dataDimRefs = None
            break

    #
    # Check whether unlinked resonances, ask to run linkResonances
    #
    
    self.checkUnlinkedResonances()
      
    #
    # TODO: for programs that have one peaklist per file, should I write a file
    # per peaklist if multiple provided or just ignore and use first one only?
    # Currently doing latter (is default!)
    #
    
    self.selectPeakLists()    
    
    #
    # Handle resonances linked via peakContribs...
    #
    
    self.createPeakContribAtomLinks()
    
    #
    # Get a list of the chains used...
    #
    
    self.getChainListFromResonanceToAtoms()
 
    #
    # Set/get the mapping for the chains
    #
    
    self.getChainMapping(resetMapping = self.resetMapping)

    # TODO: this ALSO has to be in in loop if multiple peaklists, but only
    # one format peaklist per file allowed. Has to be kept here for ANSIG type files...
    
    self.createPeakFile()

    #
    # Loop over selected peakLists
    #

    for self.peakList in self.peakLists:

      self.dataSource = self.peakList.dataSource
      
      if not self.peakList.peaks:
        print "  Warning: skipping peak list - no peaks available!"        
        continue
      
      # TODO: From here on still have to do for sparky (project), pipp!

      #
      # Get initial info
      #

      self.peakListGetDimMapping()
 
      if not self.dataDimRefs:# or (len(self.peakLists) > 1 and self.minimalPrompts):
      
        #
        # Only works for GUI environment - otherwise have to pass in mapping (or has
        # to be available)
        #
      
        if self.guiParent:
          self.peakDimSelect()
          
        else:
          self.messageReporter.showError("No peakDim mapping","No dataDimRefs passed in: is necessary for writePeaks in non-GUI environment. Aborting.",self.guiParent)
          return
          
      self.setPeakFileInfo()
      
      self.setPeakFilePeaks()
      
      self.dataDimRefs = None
      self.dimMapping = None
    
    #
    # Write file
    #
    
    if self.noWrite == False:
      self.peakFile.write(**self.writeKeywds)
      del self.peakFile
    
    return True
        
  #
  # PEAK ASSIGNMENTS
  #

  def readPeakAssignments(self,fileName = None, verbose = True, **keywds):
    
    """
    
    Generic code for reading and importing (ambiguous) peak assignment information 
    Is right now VERY XEasy specific (still under development)
    
    Mandatory argument:
    
    fileName:               file name to read the peak assignments from
    
    
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('readPeakAssignments', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
    
    self.nameMapping.isOriginalImport = True
    
    self.setPeakAssignmentFileClass()

    self.initConvertCount('peakAssignments')

    #
    # Note that fileName not necessary if self.file exists (from getProject())
    #
    
    if fileName:
    
      #
      # Read file...
      #     

      if not self.fileExists(fileName):
        return None

      self.fileName = fileName
      
      try:
        self.getPeakAssignments() # TODO TODO: check whether file read, exit if not!!
                                  # Also do this for other read/write routines!!
      except IOError, e:
        self.messageReporter.showError("Error","Error reading %s peak assignments file: %s" % (self.format,str(e)),parent = self.guiParent)
        return None
      
    elif not self.file:
      self.messageReporter.showError("Error","Need a valid filename to read peak assignments.",parent = self.guiParent)
      return
      
    #
    # Make dict of resonance names (based on application format if available)
    # Warning: for XEasy this will actually use atomSerials!
    #
    
    self.createResonanceNamesDict(self.project.currentNmrProject.resonances)
        
    # 
    # Set peaklist if not given
    #
     
    self.peakAssignmentCheckPeakList()
     
    #
    # Select which dimensions the information from the peak assignments belongs to...
    # get mapping info from peaklist. If no mapping, then exit (should code
    # in manual linker eventually)
    #
    
    self.peakAssignmentSetDimMapping()

    #
    # Set all the info in the data model
    #
    
    self.peakAssignmentSetData()
    
    #
    # Reset in case changed
    #

    self.setAssignKeyword()

    #
    # Print out import info, if necessary (can also be done later on)
    #
    
    self.printConvertCountString(verbose)

    return self.peakList


    
  def writePeakAssignments(self, fileName = None, verbose = True, **keywds):  

    """
    
    Generic code for writing and exporting peak assignment information.
    (very XEasy specific!)
    
    Mandatory argument:
    
    fileName:               file name to write the peak assignments to
    
    
    Keyword definitions in IOkeywords. 
    
    TODO TODO: add useOriginalNumbers??!
       
    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('writePeakAssignments', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
  
    self.fileName = fileName
    
    self.setPeakDimOrder()
    
    self.setPeakAssignmentFileClass()
    
    #
    # Initial check
    #
    
    if isinstance(self.peakList, (list, tuple)):
      self.peakList = self.peakList[0]
      
    if not self.peakList or not isinstance(self.peakList,Nmr.PeakList):
      self.messageReporter.showError("Error","No or invalid peak list provided for writing %s: file not written." % fileName,parent = self.guiParent)
      return

        
    #
    # Additional check - for some formats (e.g. XEasy) a shift list has to be written
    # out before the assignments can be transferred...
    #
    
    if not self.peakListAssignmentCheck([self.peakList]):    
      return

    #
    # Further initialization
    #

    self.dataSource = self.peakList.dataSource
    self.numDim = self.dataSource.numDim

    #
    # If dataDimRefs were passed in, check if they are valid
    #
    
    if self.dataDimRefs:
    
      if not isinstance(self.dataDimRefs, (list, tuple)):
        self.messageReporter.showWarning("Invalid dataDimRef list","No list of dataDimRefs passed in: are ignored.",self.guiParent)
        self.dataDimRefs = None
      
      else:
        for dataDimRef in self.dataDimRefs:
          if not isinstance(dataDimRef,Nmr.DataDimRef) or dataDimRef.dataDim.dataSource != self.dataSource:
            self.messageReporter.showWarning("Invalid dataDimRef","Invalid dataDimRef object passed in (%s). All ignored." % dataDimRef.__class__,self.guiParent)
            self.dataDimRefs = None
            break

    #
    # Check whether unlinked resonances, ask to run linkResonances
    #
    
    self.checkUnlinkedResonances()

    #
    # Get initial info
    #

    self.peakListGetDimMapping()
    
    if not self.dataDimRefs or not self.minimalPrompts:
    
      if self.guiParent:
        self.peakDimSelect()
        
      else:
        self.messageReporter.showError("No peakDim mapping","No dataDimRefs passed in: is necessary for writePeakAssignments in non-GUI environment. Aborting.",self.guiParent)
        return

    
    if not self.dimMapping:
      self.dimMapping = range(0,len(self.dataDimRefs))
    
    #
    # Create the format peak assignment file
    #
    
    self.createPeakAssignmentFile()

    #
    # Loop over peaks and set assignment information
    #

    for self.peak in self.peakList.sortedPeaks():
    
      self.peakAssignmentsSetAssignmentInfo()
      
      self.peakAssignmentsCreate()

    
    if self.noWrite == False:
      self.peakAssignmentsFile.write()
      del self.peakAssignmentsFile
    
    return True

  #
  # DISTANCE, H-BOND, DIHEDRAL, RDC, SHIFT AND JCOUPLING-type CONSTRAINTS
  #
  
  def setConstraintCodes(self,readWrite, keywds, verbose = True):
    
    """
    Defines the name for the constraint list in the Api
    """
    
    self.constraintFileCode = self.constraintType[0].upper() + self.constraintType[1:]
  
    if self.constraintType == 'distance':
      self.constraintApiCode = 'Distance'
      self.constraintResNum = 2
      
    elif self.constraintType == 'hBond':
      self.constraintApiCode = 'HBond'
      self.constraintResNum = 2
    
    elif self.constraintType == 'dihedral':
      self.constraintApiCode = 'Dihedral'
      self.constraintResNum = 4

    elif self.constraintType == 'jCoupling':
      self.constraintApiCode = 'JCoupling'
      self.constraintResNum = 2

    elif self.constraintType == 'shift':
      self.constraintApiCode = 'ChemShift'
      self.constraintResNum = 1

    elif self.constraintType == 'csa':
      self.constraintApiCode = 'Csa'
      self.constraintResNum = 1

    elif self.constraintType == 'rdc':
      self.constraintApiCode = 'Rdc'
      self.constraintResNum = 2
      
    else:
      self.messageReporter.showError("Unrecognized type","Unrecognized constraint type %s." % str(self.constraintType))
      return False
      
    if not self.setIOkeywords(readWrite + self.constraintFileCode + 'Constraints', keywds, verbose = verbose):
      return False

    return True
  
  def readDistanceConstraints(self, *args, **keywds):
    
    """
    Shell function for calling readConstraints
    """
    
    keywds['constraintType'] = 'distance'

    self.setIOkeywords('readDistanceConstraints', keywds, verbose = False)
    
    return self.readConstraints(*args, **keywds)
    
  def readHBondConstraints(self, *args, **keywds):
    
    """
    Shell function for calling readConstraints
    """
  
    keywds['constraintType'] = 'hBond'

    self.setIOkeywords('readHBondConstraints', keywds, verbose = False)
    
    return self.readConstraints(*args, **keywds)
    
  def readDihedralConstraints(self, *args, **keywds):
    
    """
    Shell function for calling readConstraints
    """
  
    keywds['constraintType'] = 'dihedral'

    self.setIOkeywords('readDihedralConstraints', keywds, verbose = False)
    
    return self.readConstraints(*args, **keywds)
    
  def readJCouplingConstraints(self, *args, **keywds):
    
    """
    Shell function for calling readConstraints
    """
  
    keywds['constraintType'] = 'jCoupling'

    self.setIOkeywords('readJCouplingConstraints', keywds, verbose = False)
    
    return self.readConstraints(*args, **keywds)
    
  def readCsaConstraints(self, *args, **keywds):
    
    """
    Shell function for calling readConstraints
    """
    
    keywds['constraintType'] = 'csa'

    self.setIOkeywords('readCsaConstraints', keywds, verbose = False)
    
    return self.readConstraints(*args, **keywds)
  
  def readRdcConstraints(self, *args, **keywds):
    
    """
    Shell function for calling readConstraints
    """
    
    keywds['constraintType'] = 'rdc'

    self.setIOkeywords('readRdcConstraints', keywds, verbose = False)
    
    return self.readConstraints(*args, **keywds)

  def readConstraints(self, fileName = None, verbose = True, **keywds):
  
    """
    
    Generic code for reading and importing constraint lists    
    Certain functions depend on constraintType and/or original format!
    
    Mandatory argument:
    
    fileName:               file name to read the constraints from
    
    
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('readConstraints', keywds, lockUpdate = True, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
    
    self.nameMapping.isOriginalImport = True
    
    #
    # Set the variable stuff
    #

    if not self.setConstraintCodes('read', keywds, verbose = verbose):
      return False

    self.ApiConstraintList = getattr(NmrConstraint,"%sConstraintList" % self.constraintApiCode)
    self.ApiConstraint = getattr(NmrConstraint,"%sConstraint" % self.constraintApiCode)    
    self.setTypeConstraint = getattr(self,"set%sTypeConstraint" % self.constraintApiCode)

    # This not always the case - ChemShiftConstraint for example
    if hasattr(NmrConstraint,"%sConstraintItem" % self.constraintApiCode):
      self.ApiConstraintItem = getattr(NmrConstraint,"%sConstraintItem" % self.constraintApiCode)
    
    self.setConstraintFileClass()    
    
    self.initConvertCount('constraints')

    #
    # Need a constraint list or an NmrConstraintStore (can also do this via a structure generation) to hang constraints from
    #

    if not self.constraintList: 
      
      if not self.nmrConstraintStore:
        
        #
        # See if a structure generation is available, if not select. TODO should I not just use an nmrConstraintStore here?
        #

        if not self.strucGen:
          
          self.selectStructureGeneration()
            
        #
        # Also create the NmrConstraintStore if not present!
        #
        # Note that this will create a new NmrConstraintStore for each StructureGeneration
        #
        
        if not self.strucGen.nmrConstraintStore:
          nmrConstraintStore = self.project.newNmrConstraintStore(nmrProject = self.project.currentNmrProject)
          self.strucGen.nmrConstraintStore = nmrConstraintStore

        self.nmrConstraintStore = self.strucGen.nmrConstraintStore 

      elif not self.strucGen:
          
       self.selectStructureGeneration()

      #
      # Set constraint list info
      #
      
      if not self.constraintListName:
        self.constraintListName = self.constraintType + ' constraint list'
        
        if not self.minimalPrompts:
  
          newConstraintListName = self.dataEntry.askString("Enter name for constraint list","New %s constraint list name" % self.constraintType,self.constraintListName,self.guiParent)
  
          if newConstraintListName:
            self.constraintListName = newConstraintListName
      
      constraintListPars = {'details': 'Autogenerated by readConstraints.'}
      
      self.constraintList = self.ApiConstraintList(self.nmrConstraintStore,name = self.constraintListName, **constraintListPars)
      self.convertCount[self.mainCode][0] += 1

    else:
      
      if isinstance(self.constraintList,self.ApiConstraintList):
      
        self.nmrConstraintStore = self.constraintList.nmrConstraintStore

      else:
      
        self.messageReporter.showError("Not valid","Constraint list of type %s (should be %s)" % (self.constraintList.className,self.constraintApiCode))
        return False
    
    #
    # Get a list of existing resonance names for this format IF resonanceMatching set to 1
    # In case reading new list and resonances already exist, this should
    # allow them to be linked up straight away
    #
  
    #
    # TODO: problem here with chainCodes... could be different!!
    #
    
    self.createResonanceNamesDict(self.nmrConstraintStore.fixedResonances)
        
    #
    # Note that fileName not necessary if self.file exists (from getProject())
    #
    
    if fileName:

      if not self.fileExists(fileName):
        return None

      self.fileName = fileName
    
      #
      # Get format specific constraint file
      #

      self.getConstraints()
      
    elif not self.file:
      self.messageReporter.showError("Error","Need a valid filename to read constraint lists.",parent = self.guiParent)
      return False
   
    #
    # Now start to actually import constraints into the data model
    # 

    if hasattr(self,'constraintFile') and self.constraintFile:

      #
      # Start analyzing...
      #

      self.setConstraints()

      #
      # Print out import info, if necessary (can also be done later on)
      #
      
      self.printConvertCountString(verbose)

      #
      # Constraint list is returned - the new resonances should be linked
      # to atoms or to existing resonances with the 'linkResonance' script
      # This should also handle the stereochemistry correctly
      #

      return self.constraintList
      
    else:
    
      return None
    
  #
  # Constraint writing....
  #

  def writeDistanceConstraints(self, *args, **keywds):
    
    """
    Shell function for calling writeConstraints
    """
    
    keywds['constraintType'] = 'distance'

    self.setIOkeywords('writeDistanceConstraints', keywds, verbose = self.verbose)
    
    return self.writeConstraints(*args, **keywds)
    
  def writeHBondConstraints(self, *args, **keywds):
    
    """
    Shell function for calling writeConstraints
    """
  
    keywds['constraintType'] = 'hBond'

    self.setIOkeywords('writeHBondConstraints', keywds, verbose = self.verbose)
    
    return self.writeConstraints(*args, **keywds)
    
  def writeDihedralConstraints(self, *args, **keywds):
    
    """
    Shell function for calling writeConstraints
    """
  
    keywds['constraintType'] = 'dihedral'

    self.setIOkeywords('writeDihedralConstraints', keywds, verbose = self.verbose)
    
    return self.writeConstraints(*args, **keywds)
    
  def writeJCouplingConstraints(self, *args, **keywds):
    
    """
    Shell function for calling writeConstraints
    """
  
    keywds['constraintType'] = 'jCoupling'

    self.setIOkeywords('writeJCouplingConstraints', keywds, verbose = self.verbose)
    
    return self.writeConstraints(*args, **keywds)

  def writeCsaConstraints(self, *args, **keywds):
    
    """
    Shell function for calling writeConstraints
    """
  
    keywds['constraintType'] = 'csa'

    self.setIOkeywords('writeCsaConstraints', keywds, verbose = self.verbose)
    
    return self.writeConstraints(*args, **keywds)
   
  def writeRdcConstraints(self, *args, **keywds):
    
    """
    Shell function for calling writeConstraints
    """
  
    keywds['constraintType'] = 'rdc'

    self.setIOkeywords('writeRdcConstraints', keywds, verbose = self.verbose)
    
    return self.writeConstraints(*args, **keywds)

  def writeConstraints(self, fileName = None, verbose = True, **keywds):

    """
    
    Generic code for writing and exporting constraint lists.    
    Certain functions depend on constraintType and/or original format.
    
    Mandatory argument:
    
    fileName:               file name to write the constraint list to
    
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('writeConstraints', keywds, lockUpdate = True, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #

    self.fileName = fileName
    self.writeKeywds = {}
    
    #
    # Set the variable stuff
    #

    if not self.setConstraintCodes('write', keywds, verbose = verbose):
      return False
  
    #
    # Quick check
    #

    if not self.constraintList or self.constraintList.className != self.constraintApiCode + 'ConstraintList':
      self.messageReporter.showError("Error","No or invalid %s constraint list provided for writing %s: file not written." % (self.constraintType,fileName),parent = self.guiParent)
      return
  
    #
    # Check whether unlinked resonances, ask to run linkResonances
    #
  
    self.checkUnlinkedResonances(nmrConstraintStore = self.constraintList.nmrConstraintStore)
            
    #
    # Create external format constraint file object
    #
     
    self.formatCode = "%s%s" % (self.formatLabel,self.constraintApiCode)

    self.setConstraintFileClass()

    self.createConstraintFile()

    self.rawConstraintClass = getattr(self.constraintsIO,"%sConstraint" % self.formatCode)
    self.rawConstraintItemClass = getattr(self.generalIO,"%sConstraintItem" % self.formatLabel)
    self.rawConstraintItemMemberClass = getattr(self.generalIO,"%sConstraintMember" % self.formatLabel)

    self.setConstraintFileConstraints = getattr(self,"set%sConstraintFileConstraints" % self.constraintApiCode)
        
    #
    # Handle resonances
    #
    
    self.createConstraintAtomLinks()
    
    #
    # Get a list of the chains used...
    #
    
    self.getChainListFromResonanceToAtoms()
 
    #
    # Set/get the mapping for the chains
    #
    
    self.getChainMapping(resetMapping = self.resetMapping)
    
    #
    # Setup the constraints and items in the constraintFile
    #
    
    self.setConstraintFileConstraints()
 
    #
    # Write file
    #
    
    if self.noWrite == False:
    
      self.constraintFile.write(**self.writeKeywds)
      del self.constraintFile
        
    return True
    
    
  def readAcqPars(self,fileName = None, verbose = True, **keywds):
     
    """
    
    Generic code for reading and importing acquisition parameters    
    
    Mandatory argument:
    
    fileName:               file name to read the acquisition pars from
    
    
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('readAcqPars', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #


    if not self.fileExists(fileName):
      return None

    self.fileName = fileName
    
    self.setAcqParsDataClass()
        
    #
    # Make sure an experiment is available before reading file...
    # necessary because have to know number of dimensions
    #

    setCurrentStore(self.project,'NmrProject')
    
    if not self.experiment:
    
      if len(self.project.currentNmrProject.experiments) == 1 and self.minimalPrompts:

        self.experiment = self.project.currentNmrProject.sortedExperiments()[0]
        
      elif self.fidDataSource:
        
        self.experiment = self.fidDataSource.experiment
      
      elif self.freqDataSource:
        
        self.experiment = self.freqDataSource.experiment

      elif self.project.currentNmrProject.experiments:

        #
        # Try to select an existing experiment
        #

        (selectionList,selectionDict) = createSelection(self.project.currentNmrProject.sortedExperiments())

        interaction = self.multiDialog.SelectionList(

                           self.guiParent,
                           selectionList,
                           selectionDict = selectionDict,
                           title = "Project '%s': " % self.project.name + 'Select experiment for acquisition parameters',
                           topText = "Format %s, acquisition file %s" % (self.format,fileName),
                           text = "Existing experiments:",
                           dismissText = 'Create new',
                           urlFile = 'SelectExperimentForAcq'

                           )

        if interaction.isSelected:
          self.experiment = interaction.selection
      
      #
      # If nothing was selected, then create a new experiment
      #

      if not self.experiment:

        while (not self.experiment):

          popup = self.multiDialog.ExperimentCreate(self.guiParent, project = self.project)

          if popup:
            self.experiment = popup.experiment
            
            # Kill window if GUI
            if hasattr(popup,'destroy'):
              popup.destroy()
                    
    self.numDim = self.experiment.numDim
  
    #
    # Read format file. Note that fileName not necessary if self.file exists (from getProject())
    #

    if self.fileName:
      
      #
      # If fileName passed in, read the file...
      #
      
      self.fileName = fileName
      
      try:
        self.getAcqPars()
      except IOError, e:
        self.messageReporter.showError("Error","Error reading %s acquisition file: %s" % (self.format,str(e)),parent = self.guiParent)
        return None

    elif not self.file:
      self.messageReporter.showError("Error","Need a valid filename to read acquisition parameters.",parent = self.guiParent)
      return

    #
    # TODO:
    #
    #  - popup for multiple SF selection
    #  - bruker IN selection for 3rd, 4th dim???
    #  - mixed Freq/Fid datasources!
    #
    
    self.setAcqParsGeneral()
            
    self.setSpectrometer(manufacturer = self.formatLabel,
                         model = self.acqParsData.specRefs['1H'].specType,
                         protonFreq = self.acqParsData.specRefs['1H'].baseFreq)
    
    self.setFidDataSource()
    
    self.setExpDimRefs()
    
    # TODO this is removed - is wrong to (re)set frequency data source when reading acquisition pars
    #self.setFreqDataSource()
        
    return (self.fidDataSource,self.freqDataSource)

  def readProcPars(self,fileName = None, verbose = True, **keywds):
     
    """
    
    Generic code for reading and importing parameters for processed spectra
    
    Mandatory argument:
    
    fileName:               file name to read the acquisition pars from
    
    
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('readProcPars', keywds, verbose = verbose):
      return
                    
    #
    # Read format file. Note that fileName not necessary if self.file exists (from getProject())
    #

    if fileName:
      
      #
      # If fileName passed in, read the file...
      #
      
      if not self.fileExists(fileName):
        return None

      self.fileName = fileName

      try:
        self.getProcPars()
      except IOError, e:
        self.messageReporter.showError("Error","Error reading %s processing parameters file: %s" % (self.format,str(e)),parent = self.guiParent)
        return None

    elif not self.file:
      self.messageReporter.showError("Error","Need a valid filename to read processing parameters.",parent = self.guiParent)
      return

    self.numDim = self.procParsFile.fPars['numDim']
    
    #
    # Make sure an experiment is available at this stage.
    #

    if not self.experiment:
    
      if not self.expName and (len(self.project.currentNmrProject.experiments) == 1 and 
          self.numDim  == self.project.currentNmrProject.sortedExperiments()[0].numDim and
          self.minimalPrompts):

        self.experiment = self.project.currentNmrProject.sortedExperiments()[0]

      elif self.project.currentNmrProject.experiments and not self.expName:

        #
        # Try to select an existing experiment
        #
        
        ll = [x for x in self.project.currentNmrProject.sortedExperiments() 
              if x.numDim == self.numDim]
        (selectionList,selectionDict) = createSelection(ll)

        interaction = self.multiDialog.SelectionList(

                           self.guiParent,
                           selectionList,
                           selectionDict = selectionDict,
                           title = "Project '%s': " % self.project.name + 'Select experiment to link processing parameters to',
                           topText = "Format %s, parameter file %s" % (self.format,fileName),
                           text = "Existing experiments:",
                           dismissText = 'Create new',
                           urlFile = 'SelectExperimentForAcq'

                           )

        if interaction.isSelected:
          self.experiment = interaction.selection
      
      #
      # If nothing was selected, then create a new experiment
      #

      if not self.experiment:

        while (not self.experiment):
        
          try:
            (directory,file) = os.path.split(fileName)
          except:
            file = fileName
          
          if self.procParsFile.fPars['nucleus'] and self.procParsFile.fPars['nucleus'][0]:
            infoText = ", isotopes '%s'" % string.join(self.procParsFile.fPars['nucleus'],',')
          elif self.procParsFile.fPars['spectrometerFreq'] and self.procParsFile.fPars['spectrometerFreq'][0]:
            sfreqs = ["%.1f" % sf for sf in self.procParsFile.fPars['spectrometerFreq']]
            infoText = ", spectrometer freqs '%s'" % string.join(sfreqs,',')
          else:
            infoText = ''

          popup = self.multiDialog.ExperimentCreate(self.guiParent,
                                                    numDim = self.numDim,
                                                    project = self.project,
                                                    expName = self.expName,
                                                    topText = "Experiment creation for file '%s'\n(numDim %d%s)" % (file,self.numDim,infoText))

          if popup and not popup.skip:
            self.experiment = popup.experiment
            popup.destroy()
          else:
            return None
    
    self.setProcParsGeneral()
    
    if not self.setProcParsDimInfo():
      # Void experiment could have been created at this stage...
      return None
    
    self.setSpectrometer(protonFreq = self.protonFreq)
         
    return self.dataSource

  def writeProcPars(self,fileName = None, verbose = True, **keywds):
    
    """
    
    Generic code for writing and exporting processing parameters    
    
    Mandatory argument:
    
    fileName:               file name to write the processing parameters to.
        
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    # Warning: for nmrPipe these are the conversion parameters!
    # nmrPipe outputFile should be of style 'auto%03d.fid' for 3D...
    #
    
    if not self.setIOkeywords('writeProcPars', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
    
    self.fileName = fileName
    
    if not self.validDataDimRefs():
      return False
      
    self.experiment = self.dataDimRefs[0].expDimRef.expDim.experiment
    self.rawData = self.experiment.rawData
    self.numDim = self.experiment.numDim
    self.expOrigin = self.rawData.name
    
    # TODO: selector if not present!
    self.inputFile = getDataSourceFileName(self.rawData)
    
    self.createProcParsFile()
    self.setProcPars()
 
    #
    # Write file
    #
    
    if self.noWrite == False:
    
      self.procParsFile.write()
      del self.procParsFile
        
    return True

  def writeProcScript(self,fileName = None, verbose = True, **keywds):
       
    """
    
    Generic code for writing and exporting processing scripts    
    
    Mandatory argument:
    
    fileName:               file name to write the processing script to.
        
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('writeProcScript', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #

    self.fileName = fileName
    
    if not self.validDataDimRefs():
      return False
          
    # azara example
    #
    # writeFile = autoproc.scr (fileName)
    # inputFile = autoproc.par
    # outputFile = spectrum.ft
    # parFile = spectrum.par
    #
    # nmrPipe inputFile should be of style 'auto%03d.fid' for 3D...
    # azara inputFile should be of style 'autoproc.par'
    #
    # nmrPipe outputFile should be of style 'auto%03d.ft' for 3D...
    # azara outputFile should be of style 'spectrum.ft'
    #
    # azara parFile should be of style 'spectrum.par'
    #

    self.experiment = self.dataDimRefs[0].expDimRef.expDim.experiment
    self.rawData = self.experiment.rawData
    self.numDim = self.experiment.numDim
    self.expOrigin = self.rawData.name
        
    self.createProcScriptFile()
    self.setProcScript()
 
    #
    # Write file
    #
    
    if self.noWrite == False:
    
      self.procScriptFile.write()
      del self.procScriptFile
        
    return True

  #
  # Citation/Person creation
  #
  
  def readPeopleAndCitations(self, fileName = None,verbose = True, **keywds):

    """
    
    Generic code for reading Citation and People info from files
         
    Mandatory argument:

    fileName:          file name to read the citation information from
    
    Keyword definitions in IOkeywords. 
        
    """

    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('readPeopleAndCitations', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
 
    self.citations = []       
    self.people = []
      
    #
    # Get the Citation info
    #

    if fileName:
    
      if not self.fileExists(fileName):
        return None
    
      self.fileName = fileName

      try:
        self.getPeopleAndCitations()
      except IOError, e:
        self.messageReporter.showError("Error","Error reading %s citation from file: %s" % (self.format,str(e)),parent = self.guiParent)
        return None
       
    elif not self.file:
      self.messageReporter.showError("Error","Need a valid file name to get citations.",parent = self.guiParent)
      return None

    #
    # Create the citations and people (authors)
    #
    #
    # TODO TODO use Autodep code to look for people!!!
    #

    if not self.doNotMakeCcpnObjects: # For AutoDep to make the objects if it needs them instead
      self.setPeople()
      self.setCitations()

    # TODO which order for people/citation creation

    return (self.people,self.citations)


  #
  # ASSIGNMENT
  #
  
  def readAssignment(self,fileName = None, verbose = True, **keywds):
  
    """
    
    Generic code for reading and importing spin system level assignment information    
    
    Mandatory argument:
    
    fileName:               file name to read the sequence from
    
    
    Keyword definitions in IOkeywords. 
       
    """
      
    #
    # Initialize the keywords
    #
    
    if not self.setIOkeywords('readAssignment', keywds, verbose = verbose):
      return
            
    #
    # Initialize other variables...
    #
    
    self.setAssignmentFileClass()

    #
    # Read format file. Note that fileName not necessary if self.file exists (from getProject())
    #

    if fileName:
      
      #
      # If fileName passed in, read the file...
      #
      
      if not self.fileExists(fileName):
        return None

      self.fileName = fileName
      try:
        self.getAssignment()
      except IOError, e:
        self.messageReporter.showError("Error","Error reading %s assignment file: %s" % (self.format,str(e)),parent = self.guiParent)
        return None

    elif not self.file:
      self.messageReporter.showError("Error","Need a valid filename to read assignment.",parent = self.guiParent)
      return None
  
    self.setAssignments()
    
    return True

  #############################################################
  # Combined functions: use two or more of the above together #
  #############################################################
  
  def readShiftsAndPeaks(self, verbose = True, **keywds):
  
    #
    # Initialize the keywords (silently because measurement list specific keywds might be available)
    #

    # Set to False so can collect all info
    self.resetConvertCount = False
    
    if not self.setIOkeywords('readShiftsAndPeaks', keywds, verbose = verbose):
      return
    
    measurementList = self.readShifts(self.shiftFile, verbose = False, minimalPrompts = self.minimalPrompts, resetIOkeywords = False)
    
    if not measurementList:
      return False
      
    peakLists = self.readPeaks(self.peakFile, verbose = False, minimalPrompts = self.minimalPrompts, resetIOkeywords = False)
    
    if not peakLists:
      return False
    
    for peakList in peakLists:
      if not peakList.dataSource.experiment in measurementList.experiments:
        measurementList.addExperiment(peakList.dataSource.experiment)
    
    if hasattr(self,'peakAssignmentFile') and self.peakAssignmentFile:
      peakList = self.readPeakAssignments(self.peakAssignmentFile, verbose = False, peakList = peakLists[0], minimalPrompts = self.minimalPrompts, resetIOkeywords = False)

    # Set back to True so resets next time
    self.resetConvertCount = True
        
    return (measurementList,peakLists)

  #
  # GENERAL
  #  

  def setIOkeywords(self,functionName, keywds, lockUpdate = False, verbose = True):
      
    returnStatus = True
    
    if keywds.has_key('resetIOkeywords'):
      doReset = keywds['resetIOkeywords']
      del(keywds['resetIOkeywords'])
      #print 'SET RESET VALUE %s' % doReset
    else:
      doReset = True
      
    #print doReset, keywds
    if self.IOkeywords.has_key(functionName):
      
      #self.mapOldKeywords(keywds)
      
      for IOkeyword in self.IOkeywords[functionName].keys():

        (default,mandatory,infoText) = self.IOkeywords[functionName][IOkeyword]
        
        if mandatory and not keywds.has_key(IOkeyword):
          
          if doReset:
            print "  Error: mandatory keyword '%s' for '%s' missing! Cannot continue." % (IOkeyword,functionName)
            returnStatus = False
            break
            
          elif not hasattr(self,IOkeyword):
            print "  Error: mandatory keyword '%s' for '%s' was not set! Cannot continue." % (IOkeyword,functionName)
            returnStatus = False
            break
          
        #print IOkeyword, doReset, default, mandatory

        if keywds.has_key(IOkeyword):
          setattr(self,IOkeyword,keywds[IOkeyword])
          #print "%s fromkeywds:" % functionName,IOkeyword,keywds[IOkeyword]
          del keywds[IOkeyword]

        elif doReset or not hasattr(self,IOkeyword):
          setattr(self,IOkeyword,default)
          #print "%s fromdefault:" % functionName,IOkeyword,default

      if keywds and verbose:
        for keywd in keywds.keys():
          print "  Warning: keyword for '%s' not recognized: '%s'" % (functionName,keywd)

      #
      # Do not reset keywords if recalling method second time - add to keywords (for Measurement and Constraint stuff!)
      # Only set this for the generic calls!
      #
      
      if lockUpdate == True:
        keywds['resetIOkeywords'] = False
    
    else:
      print "  Error: function '%s' not recognized in IOkeywords dictionary!" % (functionName)
      returnStatus = False       
      
    #print functionName
    #print keywds
    #if hasattr(self,'constraintList'):
    #  print self.constraintList
    #print
    
    # Rasmus Fogh 31/7/2012 :
    if (hasattr(self,'presetResonanceMapping') and 
        hasattr(self,'compressResonances')):
      if self.presetResonanceMapping and self.compressResonances:
        raise self.FormatConverterError(
          "presetResonanceMapping and compressResonances are not compatible")
    
    #
    # Print out an identifier for this block of reading/writing
    #
    
    if returnStatus and verbose:
      print
      print drawBox("FormatConverter: executing %s" % functionName)
      print
      
    return returnStatus

  def initConvertCount(self,dataType):
    
    if self.resetConvertCount == True:
      self.convertCount = {}
      
    if dataType == 'sequence':
      
      self.mainCode = 'molecular system(s)'

      if not self.convertCount.has_key(self.mainCode):
        if self.molSystemCreated:
          startCount = 1
        else:
          startCount = 0
        
        self.convertCount[self.mainCode] = [startCount, {'molecule(s)': 0, 'molecule residue(s)': 0, 'chain(s)': 0, 'chain residue(s)': 0}]
      
    elif dataType == 'coordinates':
      
      self.mainCode = 'structure ensemble(s)'
      
      if not self.convertCount.has_key(self.mainCode):
        self.convertCount[self.mainCode] = [0, {'model(s)': 0, 'coordinates': 0}]
      
    elif dataType == 'chemComp':
      
      self.mainCode = 'chemical component(s)'

      if not self.convertCount.has_key(self.mainCode):
        self.convertCount[self.mainCode] = [0, {'atoms': 0, 'bonds': 0}]
      
    elif dataType == 'measurements':

      self.mainCode = '%s list(s)' % self.measurementApiCode
      self.valueCode = '%s values' % self.measurementApiCode
      
      if not self.convertCount.has_key(self.mainCode):
        self.convertCount[self.mainCode] = [0, {self.valueCode: 0}]
      
    elif dataType == 'peaks':
      
      self.mainCode = 'peak list(s)'

      if not self.convertCount.has_key(self.mainCode):
        self.convertCount[self.mainCode] = [0, {'peaks': 0}]
      
    elif dataType == 'peakAssignments':
      
      self.mainCode = 'peak assignments'

      if not self.convertCount.has_key(self.mainCode):
        self.convertCount[self.mainCode] = [0, {}]
      
    elif dataType == 'constraints':

      self.mainCode = '%s constraint list(s)' % self.constraintApiCode
      self.valueCode = '%s constraint values' % self.constraintApiCode
      
      if not self.convertCount.has_key(self.mainCode):
        self.convertCount[self.mainCode] = [0, {self.valueCode: 0}]
  
  def getConvertCountString(self,startIndent = "  "):
    
    infoString = ""
    
    mainCodes = self.convertCount.keys()
    mainCodes.sort()
    
    for mainCode in mainCodes:
      
      mainCount = self.convertCount[mainCode][0]
      
      subCodeDict = self.convertCount[mainCode][1]
      subCodes = subCodeDict.keys()
      subCodes.sort()
      
      subCodeStrings = ["%d %s" % (subCodeDict[subCode],subCode) for subCode in subCodes]
      
      infoString += startIndent + "Created during import: %d %s" % (mainCount,mainCode)
      
      if subCodeStrings:
        infoString += " (with %s)" % ', '.join(subCodeStrings)
      
      infoString += "\n\n"
    
    # Reset to True if temporary switched off (for reading sequence along coordinates)...
    if self.resetConvertCount == 'temporaryFalse':
      self.resetConvertCount = True  
    
    return infoString

  def printConvertCountString(self,verbose):
    
    if verbose and self.resetConvertCount == True:
      print
      print self.getConvertCountString()
      print

  def quit(self):
  
    if self.guiParent:
      self.guiParent.quit()

    else:
      sys.exit()

  def setFormat(self):
  
    self.format = 'generic'
    
  def setGenericImports(self):
    
    pass
  
  def addWarning(self,warning):
    
    if not self.warnings.has_key(warning):
      self.warnings[warning] = 1

    else:
      self.warnings[warning] += 1

  def printWarnings(self, fout = sys.stdout):

    if self.warnings:
    
      warnings = self.warnings.keys()
      warnings.sort()
      
      for warning in warnings:
        fout.write("   Warning: %s (%d times)\n" % (warning,self.warnings[warning]))
      
      # Also reset at this stage.
      self.warnings = {}
  
  def linkResonances(self,**kw):
    
    if not kw.has_key('allowPopups'):
      kw['allowPopups'] = self.allowPopups
      
    if not kw.has_key('assignFormat'):
      kw['assignFormat'] = self.format
  
    self.linkResonancesClass = linkResonances(self.project,guiParent = self.guiParent, **kw)
  
  def checkUnlinkedResonances(self, nmrConstraintStore = None, exportMode = 'atoms', doPopup = True):
  
    """
    
    Check whether more than 5% of total number of resonances are NOT 
    assigned to atom(s). If that's the case, ask user if they
    want to run linkResonances.
    
    """
    
    runLinkResonances = False
  
    #
    # No need to check if original resonance names should be used,
    # or if minimalPrompts in on
    #
    
    if self.useOriginalResNames == 1 or self.minimalPrompts or exportMode == 'spinSystems':
      return runLinkResonances
      
    #
    # Check
    #
  
    numUnlinkedResonances = 0
    
    if not nmrConstraintStore:
      setCurrentStore(self.project,'NmrProject')
      resonances = self.project.currentNmrProject.sortedResonances()
    else:
      resonances = nmrConstraintStore.sortedFixedResonances()

    for resonance in resonances:

      resonanceSet = resonance.resonanceSet

      if not resonanceSet:
        
        numUnlinkedResonances += 1
    
    #
    # This is not sophisticated enough really... it will only pop up
    # if there's > 5% unlinked resonances but this can still be a problem
    #
    
    if numUnlinkedResonances != 0:
      
      #
      # Arbitrary threshold of 5% not linked
      #
      
      if numUnlinkedResonances / len(resonances) > 0.05:
        
        if doPopup:
      
          runLinkResonances = self.messageReporter.showYesNo("Start linkResonances","There are %d resonances with no connection to atoms - do you want to link them first?" % numUnlinkedResonances)
      
          if runLinkResonances:
                          
            popup = self.multiDialog.LinkResonancesSetup(self.guiParent,self.project,nmrConstraintStore)

        else:
          
          runLinkResonance = True

    return runLinkResonances

  """
  
  Subclass dependent functions called by above
  
  """
  
  #
  # Functions that HAVE to be defined in the format classes
  # Do not touch these or it will undoubtedly crash... saveguards
  # against the program running on format areas that haven't been coded yet
  #
  
  #
  # Sequence file reading
  #

  def setSequenceFileClass(self):
    try:
      self.SequenceFileClass = getattr(self.sequenceIO,"%sSequenceFile" % self.formatLabel)    
    except:
      if not self.minimalPrompts:
        self.messageReporter.showWarning("Warning"," No sequence IO available in format %s... aborting import/export." % (self.formatLabel),self.guiParent)
      self.quit()

  def getSequence(self):
  
    """
    
    Dummy code for reading format sequence file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No readSequence available for %s..." % self.formatLabel,self.guiParent)
    self.sequenceFile = None
    self.quit()
    
  def getSequenceGeneric(self):
    
    """
    Generic getSequence function
    """
    
    # TODO HERE: have to figure out what to do if project file read...
    self.sequenceFile = self.SequenceFileClass(self.fileName)
    
    readKeywds = self.getSequenceSetFormatSpecificReadKeywds()
    self.sequenceFile.read(**readKeywds)

    if self.verbose == 1:
      print "Reading sequence from %s file %s" % (self.formatLabel,self.fileName)
    
  def getSequenceSetFormatSpecificReadKeywds(self):
  
    return {}

  #
  # Sequence file writing
  #

  def createSequenceFile(self):
  
    """
    
    Dummy code for creating format sequence file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No writeSequence available for %s..." % self.formatLabel,self.guiParent)
    self.sequenceFile = None
    self.quit()

  def createSequenceFileGeneric(self):
 
    """
    Generic createSequenceFile function
    """

    if self.verbose == 1:
      print "Writing sequence to %s file %s" % (self.formatLabel,self.fileName)
  
    self.sequenceFile = self.SequenceFileClass(self.fileName)

  #
  # Project reading/writing
  #

  def getFullProject(self,fileName = None):
  
    """
    
    Dummy code for reading format project file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No readProject available for %s..." % self.formatLabel,self.guiParent)
    self.quit()
    
  def createFullProject(self,fileName = None):
  
    """
    
    Dummy code for creating format project file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No writeProject available for %s..." % self.formatLabel,self.guiParent)
    self.quit()

  #
  # Coordinate file reading
  #

  def setCoordinateFileClass(self):
    try:
      self.CoordinateFileClass = getattr(self.coordinatesIO,"%sCoordinateFile" % self.formatLabel)
    except:
      if not self.minimalPrompts:
        self.messageReporter.showWarning("Warning"," No coordinate IO available in format %s... aborting import/export." % (self.formatLabel),self.guiParent)
      self.quit()

  def getCoordinates(self, noSequenceRead=False):
  
    """
    
    Dummy code for reading format coordinate file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No readCoordinates available for %s..." % self.formatLabel,self.guiParent)
    self.coordinateFile = None
    self.quit()

  def getCoordinatesGeneric(self, noSequenceRead=False):
    
    """
    Generic getCoordinates function
    """
    
    # TODO HERE: have to figure out what to do if project file read...
    
    self.coordinateFile = self.CoordinateFileClass(self.fileName)
    
    # Set this value in case it is used by the coordinateFile class
    self.coordinateFile.maxModelNum = self.maxNum

    readKeywds = self.getCoordinatesSetFormatSpecificReadKeywds()

    self.coordinateFile.read(maxNum = self.maxNum, **readKeywds)
    
    # Set version if any available
    if not self.version and hasattr(self.coordinateFile,'version'):
      self.version = self.coordinateFile.version
    
    self.getCoordinatesFormatSpecific()
    
    chains = []
    for molSystem in self.project.molSystems:
      for chain in molSystem.chains:
        chains.append(chain)

    if not noSequenceRead and (not chains or self.forceReadSequence):
    
      #
      # Automatically read in sequence if no molSystem...
      # pretend there's a project (bit of a hack)
      #
      
      self.getCoordinatesSetSequenceFile()
      self.file = 'temporary'
      
      #
      # Set special flags that come in from readCoordinates
      #
      
      readSeqKeywds = {'verbose': False}
      
      for readSeqKeywd in ('molSystem','autoCreateChemComps','ignoreUnknownChemComps','ignoreSingleOther','version','chemCompPath','addNameMappings','ignoreResNames','minimalPrompts'):
        
        value = getattr(self,readSeqKeywd)
        readSeqKeywds[readSeqKeywd] = value

      #
      # Create the molecular system in the data model
      #
      #print "SEQUENCE READ", readSeqKeywds
      self.readSequence(**readSeqKeywds)
      
      if self.resetConvertCount == True:
        self.resetConvertCount = 'temporaryFalse'
      
      # Have to reset this here, otherwise won't have coordinate info
      self.initConvertCount('coordinates')
 
      self.file = None
      
      sequenceText = 'and sequence '
    
    else:
      sequenceText = ''

    if self.verbose == 1:
      print "Reading coordinates %sfrom %s file %s" % (sequenceText,self.formatLabel,self.fileName)

  def getCoordinatesSetFormatSpecificReadKeywds(self):

    return {}

  def getCoordinatesFormatSpecific(self):

    pass

  def getCoordinatesSetSequenceFile(self):
    
    self.setSequenceFileClass()
    self.sequenceFile = self.SequenceFileClass(self.fileName)
    
    readKeywds = self.getSequenceSetFormatSpecificReadKeywds()
    self.sequenceFile.read(coordinateFile = self.coordinateFile, **readKeywds)
    
  #
  # Coordinate file writing
  #
  
  def createCoordinateFile(self):
  
    """
    
    Dummy code for writing format coordinate file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No writeCoordinates available for %s..." % self.formatLabel,self.guiParent)
    self.coordinateFile = None
    self.quit()

  def createCoordinateFileGeneric(self):
   
    """
    Generic createCoordinateFile function
    """
  
    if self.verbose == 1:
      print "Writing coordinates to %s file %s" % (self.formatLabel,self.fileName)
  
    self.coordinateFile = self.CoordinateFileClass(self.fileName)
      
    self.createCoordinateFileFormatSpecific()

  def createCoordinateFileFormatSpecific(self):
    
    pass

  #
  # ChemComp file reading
  #

  def setChemCompFileClass(self):
    try:
      self.ChemCompFileClass = getattr(self.chemCompsIO,"%sChemCompFile" % self.formatLabel)
    except:
      if not self.minimalPrompts:
        self.messageReporter.showWarning("Warning"," No chemComp IO available in format %s... aborting import/export." % (self.formatLabel),self.guiParent)
      self.quit()

  def getChemComps(self):
    
    """
    
    Dummy code for reading format chemComp file: is always specifically defined in
    format specific class
    
    """

    self.messageReporter.showWarning("Warning"," No readChemComps available for %s..." % self.formatLabel,self.guiParent) 
    self.chemCompFile = None
    self.quit()
     
  def getChemCompsGeneric(self):
    
    """
    Generic getChemComps function
    """
  
    self.chemCompFile = self.ChemCompFileClass(self.fileName)
    self.chemCompFile.read()

    if self.verbose == 1:
      print "Reading chemComps from %s file %s" % (self.formatLabel,self.fileName)
  
  #
  # ChemComp file writing
  #
  
  def createChemCompFile(self):
  
    """
    
    Dummy code for writing format chemComp file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No writeChemComps available for %s..." % self.formatLabel,self.guiParent)
    self.chemCompFile = None
    self.quit()

  def createChemCompFileGeneric(self):
   
    """
    Generic createChemCompFileGeneric function
    """
  
    if self.verbose == 1:
      print "Writing chemComp to %s file %s" % (self.formatLabel,self.fileName)
  
    self.chemCompFile = self.ChemCompFileClass(self.fileName)

  #
  # Measurement file reading
  #

  def setMeasurementFileClass(self):
    try:
      self.measurementsIO = getattr(self,"%sIO" % self.formatFileBaseName)
      self.MeasurementFileClass = getattr(self.measurementsIO,"%s%sFile" % (self.formatLabel,self.formatFileClassName))
    except:
      if not self.minimalPrompts:
        self.messageReporter.showWarning("Warning"," No IO available for %s in format %s... aborting import/export." % (self.formatFileValueKey,self.formatLabel),self.guiParent)
      self.quit()

  def getMeasurements(self):
  
    """
    
    Dummy code for reading format measurement file: is always specifically defined in
    format specific class
    
    """

    self.messageReporter.showWarning("Warning"," No read%ss available for %s..." % (self.measurementApiCode,self.formatLabel),self.guiParent)     
    self.measurementFile = None
    self.quit()

  def getMeasurementsGeneric(self):
    
    """
    Generic getMeasurements function
    """
    
    #
    # Get a measurement list
    #
    
    try:

      if self.verbose == 1:
        print "Reading %ss from %s file %s" % (self.measurementType,self.formatLabel,self.fileName)

      # TODO HERE: have to figure out what to do if project file read...
      self.measurementFile = self.MeasurementFileClass(self.fileName)
           
      readKeywds = self.getMeasurementsSetFormatSpecificReadKeywds()
      self.measurementFile.read(**readKeywds)
      
      sequenceText = self.getMeasurementsSequenceRead()

      if self.verbose == 1:
        print "Reading %ss%s from %s file %s" % (self.measurementType,sequenceText,self.formatLabel,self.fileName)

    except:

      errorMessage = traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
    
      if self.verbose:
        self.messageReporter.showWarning("Warning"," Cannot read %ss for %s...:\n%s." % (self.measurementType,self.formatLabel,errorMessage),self.guiParent)
        
      self.measurementFile = None
      
      return traceback.format_exception(sys.exc_type,sys.exc_value,sys.exc_info()[2]) 

  def getMeasurementsSetFormatSpecificReadKeywds(self):
    
    return {}

  def getMeasurementsSequenceRead(self):
    
    sequenceText = ''
    
    return sequenceText
            
  #
  # Measurement file writing
  #

  def createMeasurementFile(self):
  
    """
    
    Dummy code for creating format measurement file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No write%ss available for %s..." % (self.measurementApiCode,self.formatLabel),self.guiParent)
    self.measurementFile = None
    self.quit()

  def createMeasurementFileGeneric(self):
    
    """
    Generic createMeasurementsFile function
    """
    
    try:
    
      self.measurementFile = self.MeasurementFileClass(self.fileName)
      
      self.createMeasurementFileFormatSpecific()
     
      if self.verbose == 1:
        print "Writing %ss to %s file %s" % (self.measurementType,self.formatLabel,self.fileName)
  
    except:

      self.messageReporter.showWarning("Warning"," No write%s available for %s..." % (self.measurementApiCode,self.formatLabel),self.guiParent)

  def createMeasurementFileFormatSpecific(self):
  
    pass
  
  #
  # Peak file reading
  #

  def setPeakFileClass(self):
    
    try:
      self.PeakFileClass = getattr(self.peaksIO,"%sPeakFile" % self.formatLabel)
    except:
      if not self.minimalPrompts:
        self.messageReporter.showWarning("Warning"," No peaks IO available in format %s... aborting import/export." % (self.formatLabel),self.guiParent)
      self.quit()
      
  def getPeaks(self):
    
    """
    
    Dummy code for reading format peak file: is always specifically defined in
    format specific class
    
    """

    self.messageReporter.showWarning("Warning"," No readPeaks available for %s..." % self.formatLabel,self.guiParent)     
    self.peakFile = None
    self.quit()
    
  def getPeaksGeneric(self):
  
    if self.verbose == 1:
      print "Reading peak list from %s file %s" % (self.formatLabel,self.fileName)

    # TODO HERE: have to figure out what to do if project file read...
   
    # Read in whole file (even if only experiment wanted later on...)
    
    self.peakFile = self.PeakFileClass(self.fileName)
    self.peakFile.read()
   
  #
  # Peak file writing
  #

  def createPeakFile(self):
  
    """
    
    Dummy code for creating format peak file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No writePeaks available for %s..." % self.formatLabel,self.guiParent)
    self.peakFile = None
    self.quit()

  def createPeakFileGeneric(self):
  
    if self.verbose == 1:
      print "Writing peak list(s) to %s file %s" % (self.formatLabel,self.fileName)
  
    self.peakFile = self.PeakFileClass(self.fileName)
    
    self.createPeakFileFormatSpecific()
    
  def createPeakFileFormatSpecific(self):
    
    pass
  
  #
  # Peak assignments file reading
  #

  def setPeakAssignmentFileClass(self):
    try:
      self.PeakAssignmentFileClass = getattr(self.peakAssignmentsIO,"%sPeakAssignmentFile" % self.formatLabel)
    except:
      if not self.minimalPrompts:
        self.messageReporter.showWarning("Warning"," No chemComp IO available in format %s... aborting import/export." % (self.formatLabel),self.guiParent)
      self.quit()

  def getPeakAssignments(self):
    
    """
    
    Dummy code for reading format peak assignments file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No readPeakAssignments available for %s..." % self.formatLabel,self.guiParent)     
    self.peakAssignmentsFile = None
    self.quit()

  #
  # Peak assignments file writing
  #

  def createPeakAssignmentFile(self):
  
    """
    
    Dummy code for creating format peak assignments file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No writePeakAssignments available for %s..." % self.formatLabel,self.guiParent)
    self.peakAssignmentsFile = None
    self.quit()

  #
  # Constraints file reading
  #
  
  def setConstraintFileClass(self):

    try:
      self.constraintsIO = getattr(self,"%sConstraintsIO" % self.constraintType)
      self.ConstraintFileClass = getattr(self.constraintsIO,"%s%sConstraintFile" % (self.formatLabel,self.constraintFileCode))
    except:
      if not self.minimalPrompts:
        self.messageReporter.showWarning("Warning"," No IO available for %s constraints in format %s... aborting import/export." % (self.constraintType,self.formatLabel),self.guiParent)
      self.quit()
        
  def getConstraints(self):
    
    """
    
    Dummy code for reading format constraints file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No readConstraints available for %s..." % self.formatLabel,self.guiParent)
    self.constraintFile = None
    self.quit()

  def getConstraintsGeneric(self):
    
    """
    Generic getConstraints function
    """   
    try:
    
      self.constraintFile = self.ConstraintFileClass(self.fileName)
      self.constraintFile.read()

      if self.verbose == 1:
        print "Reading %s constraint list from %s file %s" % (self.constraintType,self.formatLabel,self.fileName)

    except:

      errorMessage = traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
      
      if self.verbose:
        self.messageReporter.showWarning("Warning"," Cannot read %s constraints for %s...:\n%s" % (self.constraintApiCode,self.formatLabel,errorMessage),self.guiParent)

      self.constraintFile = None
      
      return traceback.format_exception(sys.exc_type,sys.exc_value,sys.exc_info()[2]) 

  #
  # Constraints file writing
  #

  def createConstraintFile(self):
  
    """
    
    Dummy code for creating format constraints file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No writeConstraints available for %s..." % self.formatLabel,self.guiParent)
    self.constraintFile = None
    self.quit()

  def createConstraintFileGeneric(self):
  
    """
    Generic createConstraintFile function
    """   

    try:

      self.constraintFile = self.ConstraintFileClass(self.fileName)
      
      self.createConstraintFileFormatSpecific()
      
      if self.verbose == 1:
        print "Writing %s constraints to %s file %s" % (self.constraintType,self.formatLabel,self.fileName)
  
    except:
      
      errorMessage = traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
    
      self.messageReporter.showWarning("Warning"," No write%sConstraints available for %s...:\n%s" % (self.constraintApiCode,self.formatLabel,errorMessage),self.guiParent)
      self.quit()

  def createConstraintFileFormatSpecific(self):
  
    pass

  #
  # Acquisition parameters file reading
  #

  def setAcqParsDataClass(self):
    try:
      self.AcqParsDataClass = getattr(self.acqParsIO,"%sAcqParData" % self.formatLabel)
    except:
      if not self.minimalPrompts:
        self.messageReporter.showWarning("Warning"," No chemComp IO available in format %s... aborting import/export." % (self.formatLabel),self.guiParent)
      self.quit()

  def getAcqPars(self):
  
    """
    
    Dummy code for reading format acquisition parameters file(s): is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No readAcqPars available for %s..." % self.formatLabel,self.guiParent)
    self.acqParsData = None
    self.quit()    

  #
  # Processing parameters file reading
  #

  def getProcPars(self):
  
    """
    
    Dummy code for reading format processing parameters file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No readProcPars available for %s..." % self.formatLabel,self.guiParent)
    self.procParsFile = None
    self.quit()
   
  #
  # Processing parameters file writing
  #

  def createProcParsFile(self):
  
    """
    
    Dummy code for creating processing parameters file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No writeProcPars available for %s..." % self.formatLabel,self.guiParent)
    self.procParsFile = None
    self.quit()

    
  #
  # Processing script file writing
  #

  def createProcScriptFile(self):
  
    """
    
    Dummy code for creating processing script file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No writeProcScript available for %s..." % self.formatLabel,self.guiParent)
    self.procScriptFile = None
    self.quit()


  #
  # People and Citation file reading
  #

  def getPeopleAndCitations(self):
    
    """
    
    Dummy code for reading format citation file: is always specifically defined in
    format specific class
    
    """

    self.messageReporter.showWarning("Warning"," No readPeopleAndCitations available for %s..." % self.formatLabel,self.guiParent) 
    self.citationFile = None
    self.quit()

  #
  # Sequence file reading
  #

  def setAssignmentFileClass(self):
    try:
      self.AssignmentFileClass = getattr(self.assignmentIO,"%sAssignmentFile" % self.formatLabel)    
    except:
      if not self.minimalPrompts:
        self.messageReporter.showWarning("Warning"," No assignment IO available in format %s... aborting import/export." % (self.formatLabel),self.guiParent)
      self.quit()

  def getAssignment(self):
  
    """
    
    Dummy code for reading format assignment file: is always specifically defined in
    format specific class
    
    """
  
    self.messageReporter.showWarning("Warning"," No readAssignment available for %s..." % self.formatLabel,self.guiParent)
    self.assignmentFile = None
    self.quit()
    
  def getAssignmentGeneric(self):
    
    """
    Generic getAssignment function
    """
    
    self.assignmentFile = self.AssignmentFileClass(self.fileName)
    
    readKeywds = self.getAssignmentSetFormatSpecificReadKeywds()
    self.assignmentFile.read(**readKeywds)

    if self.verbose == 1:
      print "Reading assignments from %s file %s" % (self.formatLabel,self.fileName)
    
  def getAssignmentSetFormatSpecificReadKeywds(self):
  
    return {}

  #
  # Normal functions (default setup)
  #
    
  def setChemCompRepositories(self):

    #
    # Make sure current project is part of chemComp search (will be saved here if possible)
    #

    packageLocator = self.project.findFirstPackageLocator(targetName = 'ccp.molecule.ChemComp')
    
    repository = self.project.findFirstRepository(name = 'userData')
    repositories = list(packageLocator.repositories)
    
    if not repository in repositories:
      if hasattr(self,'saveChemCompLocal') and self.saveChemCompLocal:
        repositories.insert(0,repository)
      else:
        repositories.append(repository)
      packageLocator.repositories = repositories

  def setChemCompCoordRepositories(self):

    #
    # Make sure current project is part of chemCompCoord search (will be saved here if possible)
    #

    packageLocator = self.project.findFirstPackageLocator(targetName = 'ccp.molecule.ChemCompCoord')
    
    repository = self.project.findFirstRepository(name = 'userData')
    repositories = list(packageLocator.repositories)

    if not repository in repositories:
      if hasattr(self,'saveChemCompLocal') and self.saveChemCompLocal:
        repositories.insert(0,repository)
      else:
        repositories.append(repository)
      packageLocator.repositories = repositories
              
  def setAssignKeyword(self):
    
    """
    
    Reset the assignKeyword. This keyword is used in applicationData to keep
    track of assignments and can be changed in places to 'atomSerial' (for XEasy
    for example)
    
    """
  
    self.assignKeyword = assign_kw

  def checkMeasurementsAndResonances(self):
  
    """
    
    Make a dict of resonance names, and check whether a measurement list that
    is exactly the same exists already
    
    """
    
    setCurrentStore(self.project,'NmrProject')
  
    #
    # Make dict of resonance names (based on application format)
    #
  
    self.resonanceNames = getApplResNames(self.format,self.project.currentNmrProject.resonances)

    #
    # Check if exact same chemical shift list exists already.
    # 
    
    self.checkMeasurementListExists()

    #
    # Remove dict of resonances if resonanceMatching not on
    #
  
    if not self.resonanceMatching:
      self.resonanceNames = {}

  def setMeasurementsExportMode(self):
    
    # Set the export mode, either 'atoms' or 'spinSystems', depending on format.
    self.exportMode = 'atoms'

  def setSingleResonanceMeasurementsGeneric(self):
    
    """
    Set measurements connected to single resonances.
    """
    
    spinSystems = {}
    addSpinSystemInfo = ['resLabel','prevSpinSystemId']
    
    #
    # Assuming that resonances listed only once...
    #
    
    for self.rawMeasurement in self.measurementFileValues:    
    
      # 
      # Assuming only one chain listed in measurement list
      # TODO: if there are cs lists with chainCodes, need to handle that here...
      #

      resName = self.getMeasurementResName(self.rawMeasurement)
      
      if not resName:
      
        continue
      
      makeNewResonance = 0
      
      #
      # Check if this resonance already exists...
      # Note that the resonanceNames dict is empty at this point if resonanceMatching was off
      #
      
      if self.resonanceNames.has_key(resName):
        
        #
        # Link measurement to ALL resonances with this name!
        # If there's any ambiguity a new resonance is created later (solve with linkResonance)
        #

        for resonance in self.resonanceNames[resName]:

          findFunction = getattr(resonance,"findFirst%s" % self.measurementApiCode)
          self.measurement = findFunction(parentList = self.measurementList)

          if self.measurement:
            self.setMeasurementSpecific(resName,resonance)
            
          else:          
            self.measurement = self.createMeasurement(resonance,self.rawMeasurement)

      else:

        #
        # Make a new resonance if not found in resonance dict
        #

        makeNewResonance = 1
        
      if makeNewResonance == 1:
      
        resonance = self.createResonance(resName,origSource = self.rawMeasurement)
        self.measurement = self.createMeasurement(resonance,self.rawMeasurement)
        
        self.resonanceNames[resName] = [resonance]
        
        if self.rawMeasurement.seqCode == None and hasattr(self.rawMeasurement,'spinSystemId'):
          spinSystemId = self.rawMeasurement.spinSystemId
          if not spinSystems.has_key(spinSystemId):
            spinSystems[spinSystemId] = [None,None,[]]
          
            for i in range(0,len(addSpinSystemInfo)):
              attrName = addSpinSystemInfo[i]
              if hasattr(self.rawMeasurement,attrName):
                spinSystems[spinSystemId][i] = getattr(self.rawMeasurement,attrName)
          
          spinSystems[spinSystemId][-1].append(resonance)
          
          # Also set the name and assignNames for the resonance in this case
          resonance.name = resName
          resonance.assignNames = [self.rawMeasurement.atomName]
            
          #
          # Finally also tag the measurement with the original resonance name.
          #
    
          if self.measurement:

            applData = self.measurement.findFirstApplicationData(application = self.format, keyword = assign_kw, value = resName)

            if not applData:
              self.measurement.addApplicationData(Implementation.AppDataString(application = self.format, keyword = assign_kw, value = resName))
            
    #
    # Set spinSystems...
    #
    
    if spinSystems:
      self.createSpinSystems(spinSystems,addSpinSystemInfo)

  def setMultipleResonanceMeasurementsGeneric(self):
    
    """
    Set measurements connected to multiple resonances. Should support spin systems, not tested.
    """
    
    spinSystems = {}
    addSpinSystemInfo = ['resLabel','prevSpinSystemId']
        
    #
    # Loop over the raw measurement data from the format file
    #
    
    for self.rawMeasurement in self.measurementFileValues:
      
      resonances = []
      resNames = []
      
      #
      # Have to loop over info for different resonances
      #

      for resIndex in range(1,self.measurementResNum + 1):
        
        resName = self.getMeasurementResName(self.rawMeasurement,resIndex = resIndex)
      
        if not resName:
          break

        resNames.append(resName)
        
        #
        # Check if this resonance already exists...
        # Note that the resonanceNames dict is empty at this point if resonanceMatching was off
        #
        
        if self.resonanceNames.has_key(resName):
          
          #
          # Link measurement to ALL resonances with this name!
          # If there's any ambiguity a new resonance is created later (solve with linkResonance)
          #
  
          resonances.append(self.resonanceNames[resName])
  
        else:
  
          #
          # Make a new resonance if not found in resonance dict
          #
        
          resonance = self.createResonance(resName,origSource = self.rawMeasurement)
          resonances.append([resonance])

          self.resonanceNames[resName] = [resonance]

          if getattr(self.rawMeasurement,"seqCode%d" % resIndex) == None and hasattr(self.rawMeasurement,'spinSystemId%d' % resIndex):
            spinSystemId = getattr(self.rawMeasurement,"spinSystemId%d" % resIndex)
            if not spinSystems.has_key(spinSystemId):
              spinSystems[spinSystemId] = [None,None,[]]
            
              for i in range(0,len(addSpinSystemInfo)):
                attrName = "%s%d" % (addSpinSystemInfo[i],resIndex)
                if hasattr(self.rawMeasurement,attrName):
                  spinSystems[spinSystemId][i] = getattr(self.rawMeasurement,attrName)
            
            spinSystems[spinSystemId][-1].append(resonance)
            
            # Also set the name and assignNames for the resonance in this case
            resonance.name = resName
            resonance.assignNames = [getattr(self.rawMeasurement,"atomName%d" % resIndex)]

      #
      # Now have the information, create the measurement
      #
      
      if len(resNames) == self.measurementResNum:
        
        resonanceCombs = [[resonance] for resonance in resonances[0]]
        
        self.getResonanceCombinations(resonances[1:],resonanceCombs)
        
        for resonanceComb in resonanceCombs:
          
          self.measurement = self.measurementList.findFirstMeasurement(resonances = set(resonanceComb))

          if self.measurement:
            self.setMeasurementSpecific(resNames,resonanceComb)
          
          # Avoid having same resonance in here - meaningless. Occurs due to problems in input files
          elif resonanceComb.count(resonanceComb[0]) != self.measurementResNum:

            self.measurement = self.createMeasurement(resonanceComb,self.rawMeasurement)
                      
            #
            # Also tag the measurement with the original resonance names.
            #
      
            applData = self.measurement.findFirstApplicationData(application = self.format, keyword = assign_kw, value = str(resNames))

            if not applData:
              self.measurement.addApplicationData(Implementation.AppDataString(application = self.format, keyword = assign_kw, value = str(resNames)))
          
          else:
            print("  Warning: measurement linked to same resonance with name %s - ignored." % resNames[0])
    #
    # Set spinSystems, if any...
    #
    
    if spinSystems:
      self.createSpinSystems(spinSystems,addSpinSystemInfo)
  
  def setHExchProtections(self):
    
    self.setSingleResonanceMeasurementsGeneric()

  def setHExchRates(self):
    
    self.setSingleResonanceMeasurementsGeneric()

  def setJCouplings(self):

    self.setMultipleResonanceMeasurementsGeneric() 

  def setNoes(self):

    self.setMultipleResonanceMeasurementsGeneric() 

  def setRdcs(self):

    self.setMultipleResonanceMeasurementsGeneric() 

  def setShiftAnisotropys(self):
    
    self.setSingleResonanceMeasurementsGeneric()

  def setShifts(self):
    
    self.setSingleResonanceMeasurementsGeneric()
    
  def setT1s(self):

    self.setSingleResonanceMeasurementsGeneric() 
  
  def setT1Rhos(self):

    self.setSingleResonanceMeasurementsGeneric() 

  def setT2s(self):

    self.setSingleResonanceMeasurementsGeneric() 

  def getResonanceCombinations(self,resonanceLists,resonanceCombs):
    
    resonanceList = resonanceLists.pop(0)
    
    for resonanceComb in resonanceCombs:
      resonanceCombBase = resonanceComb[:]
      
      resonanceComb.append(resonanceList[0])
      
      for resonance in resonanceList[1:]:
        resonanceCombs.append(resonanceCombBase[:] + [resonance])
        
    if resonanceLists:
      self.getResonanceCombinations(resonanceLists,resonanceCombs)


  def createSpinSystems(self,spinSystems,addSpinSystemInfo):
    
    #
    # This code not perfect - assumes that spinSystem numbering is 'sequential'. Works for now.
    #
    
    spinSystemIds = spinSystems.keys()
    spinSystemIds.sort()
    
    resGroupsById = {}

    for spinSystemId in spinSystemIds:
      
      if 'resLabel' in addSpinSystemInfo:
        ccpCode = spinSystems[spinSystemId][addSpinSystemInfo.index('resLabel')]
      else:
        ccpCode = None
      
      resGroup = Nmr.ResonanceGroup(self.project.currentNmrProject,name = str(spinSystemId),ccpCode = ccpCode,resonances = spinSystems[spinSystemId][-1],details = 'Set by DataFormat from %s format file' % self.format)
      resGroupsById[spinSystemId] = resGroup

      if 'prevSpinSystemId' in addSpinSystemInfo:
        prevSpinSystemId = spinSystems[spinSystemId][addSpinSystemInfo.index('prevSpinSystemId')]
        if prevSpinSystemId != None:
          prevResGroup = resGroupsById[prevSpinSystemId]
          resGroupProb = Nmr.ResonanceGroupProb(resGroup,linkType = 'sequential', sequenceOffset = -1, possibility = prevResGroup)
      
  def checkMolecules(self):
  
    """
    
    Get the sequence info from the format file object, check whether
    that molecule already exists, and create new one if necessary.
    
    """
    
    #
    # Set up some variables
    #
    
    origSequenceDict = {}
    createMoleculeDict = {}    
    chemCompDict = {}
    
    self.molResToSeqEl = {}

    molType = None
    
    packageLocator = self.project.findFirstPackageLocator(targetName = 'ccp.molecule.ChemComp')
    repositoryList = packageLocator.repositories
        
    #
    # Loop over the sequences
    #
    
    numOrigSequences = len(self.sequenceFile.sequences)
    
    for self.sequence in self.sequenceFile.sequences:
      
      #
      # Check if there's any info available - skip if not
      #

      if not self.sequence.elements and not self.minimalPrompts:
        self.messageReporter.showError("Error","No sequence elements for %s - ignoring." % (self.sequence.molName), parent = self.guiParent)
        self.addWarning("Molecule %s not created because it has no sequence elements." % self.sequence.molName)
        continue
 
      #
      # Check if water, treat differently...
      #
      
      isWater = False
      if len(self.sequence.elements) == 1:
        seqEl = self.sequence.elements[0]
        if hasattr(seqEl,'code3Letter') and seqEl.code3Letter in ('Hoh','HOH'):
          isWater = True

      #
      # Rearrange so exactly same sequences are joined at this stage...
      #
      
      sameMolecule = False
      
      for origMolName in origSequenceDict.keys():
        otherSequence = origSequenceDict[origMolName][0][0]
        
        if isWater:
          otherSeqEl = otherSequence.elements[0]
          if hasattr(otherSeqEl,'code3Letter') and otherSeqEl.code3Letter in ('Hoh','HOH'):
            if otherSequence.originalChainCode == self.sequence.originalChainCode:
              sameMolecule = True

        elif len(self.sequence.elements) == len(otherSequence.elements):
          sameMolecule = True
          for seqIndex in range(len(self.sequence.elements)):
            seqEl = self.sequence.elements[seqIndex]
            otherSeqEl = otherSequence.elements[seqIndex]
            if hasattr(seqEl,'code1Letter') and seqEl.code1Letter:
              if seqEl.code1Letter != otherSeqEl.code1Letter:
                sameMolecule = False
                break
            elif hasattr(seqEl,'code3Letter'):
              if seqEl.code3Letter != otherSeqEl.code3Letter:
                sameMolecule = False
                break

        if sameMolecule:
          break

      #
      # Keep track of sequences and continue...
      #
      
      if sameMolecule:

        if not isWater:
      
          origSequenceDict[origMolName][0].append(self.sequence)
        
          if createMoleculeDict.has_key(origMolName):
            createMoleculeDict[origMolName][1].append((self.sequence,[]))
            for seqIndex in range(len(self.sequence.elements)):
              (otherSeqElement,chemComp,chemCompVar) = createMoleculeDict[origMolName][1][0][1][seqIndex]
              createMoleculeDict[origMolName][1][-1][1].append((self.sequence.elements[seqIndex],chemComp,chemCompVar))

          continue
              
        else:
          
          if self.sequence.originalChainCode == origSequenceDict[origMolName][0][0].originalChainCode:
            
            """
            (otherSeqElement,chemComp,chemCompVar) = createMoleculeDict[origMolName][1][0][1][0]
            createMoleculeDict[origMolName][1][-1][0].elements.append(self.sequence.elements[0])
            createMoleculeDict[origMolName][1][-1][1].append(self.sequence.elements[0],chemComp,chemCompVar)
            """
            
            for seqEl in self.sequence.elements:
              origSequenceDict[origMolName][1].append((seqEl,'other'))
              origSequenceDict[origMolName][0][0].elements.append(seqEl)
              print "   Added water to molecule %s, original code %s" % (origMolName,seqEl.seqCode)
            
            continue
          

      #
      # Check whether a molName is available from the input or the sequence file
      #
      
      molNameIndex = 1
      backupMolName = 'Molecule'
      
      if numOrigSequences == 1 and self.molName:
        baseMolName = self.molName[:76]

      elif self.sequence.molName:
        baseMolName = self.sequence.molName[:76]
        
      else:
        baseMolName = backupMolName + "_%d" % molNameIndex

      if baseMolName:
        baseMolName = baseMolName.strip()

      origMolName = baseMolName
            
      while (origSequenceDict.has_key(origMolName)):
        if origMolName[:len(backupMolName)] == backupMolName:
          origMolName = backupMolName + "_%d" % molNameIndex
        else:
          origMolName = baseMolName + "_%d" % molNameIndex
        molNameIndex += 1
        
      # Modify for water
      if isWater:
        origMolName = "%s_%s" % (origMolName,self.sequence.originalChainCode)
      
      origSequenceDict[origMolName] = ([self.sequence],[])
      
      #
      # Determine if there is mapping available, set molecule if so...
      #
    
      molMappings = self.nameMapping.findAllMoleculeMappings(formatName = origMolName)
      mappingFound = False
      
      if molMappings:
      
        unusedSequenceElements = self.sequence.elements[:]
        molecules = []
      
        for molMapping in molMappings:

          molName = molMapping.name
          molecule = self.project.findFirstMolecule(name = molName)

          if molecule:
            if not molecule in molecules:
              molecules.append(molecule)
          
        if len(molecules) == 1 and molecules[0].seqLength == len(self.sequence.elements):

          print "  Autoset molecule %s from single mapping." % molName
          createMoleculeDict[origMolName] = [molecule,[(self.sequence,[])]]
          
          for seqIndex in range(0,len(self.sequence.elements)):
            seqElement = self.sequence.elements[seqIndex]
            molResidue = molecule.findFirstMolResidue(serial = (seqIndex + 1))
            chemComp = molResidue.chemComp
            chemCompVar = molResidue.chemCompVar
            createMoleculeDict[origMolName][1][0][1].append((seqElement,chemComp,chemCompVar))
            
            self.molResToSeqEl[molResidue] = seqElement
          
          mappingFound = True
            
        else:
            
          print "  Warning: mappings to multiple molecules were found... no automatic mapping done."
      
      #
      # If mapping found, continue
      #
      
      if mappingFound:      
        continue
      

      #
      # If 1 residue only, check if hetero atom/group ('other')
      #

      if len(self.sequence.elements) == 1:
        sequenceElement = self.sequence.elements[0]
        
        if hasattr(sequenceElement,'residueType') and sequenceElement.residueType == 'other':
          origSequenceDict[origMolName][1].append((sequenceElement,'other'))
          print "  Molecule %s set to type 'other' from sequence file information." % origMolName
          
      #
      # If nothing set yet, then do polymer checking...
      #
      
      if not origSequenceDict[origMolName][1]:
            
        #
        # Polymer
        #

        if hasattr(self.sequence,'polymerType') and self.sequence.polymerType != 'unknown':

          print "  Molecule %s set to type %s from sequence file information." % (origMolName,self.sequence.polymerType)
          
          for sequenceElement in self.sequence.elements:
            origSequenceDict[origMolName][1].append([sequenceElement,self.sequence.polymerType])
    
        else:

          #
          # Try determine molecule type (protein,DNA,RNA) from sequence
          #
         
          unknownSequenceElements = []
         
          for sequenceElementIndex in range(0,len(self.sequence.elements)):
            
            sequenceElement = self.sequence.elements[sequenceElementIndex]
            
            if hasattr(sequenceElement,'residueType') and sequenceElement.residueType:
            
              #
              # If info available from sequence, use that...
              #
              
              origSequenceDict[origMolName][1].append([sequenceElement,sequenceElement.residueType])
            
            else:
            
              #
              # Else try to determine what it is...
              #
              
              origSequenceDict[origMolName][1].append([sequenceElement,None])
              unknownSequenceElements.append([sequenceElementIndex,[]])

              if hasattr(sequenceElement,'code1Letter'):

                for bioPolymerType in bioPolymerCodes.keys():

                  if sequenceElement.code1Letter in bioPolymerCodes[bioPolymerType][0]:
                    unknownSequenceElements[-1][-1].append(bioPolymerType)

              elif hasattr(sequenceElement,'code3Letter'):

                #
                # If '3' letter codes, could be both 1-letter one or proprietary (e.g. GUA, ...)
                # Have to translate if proprietary... only one-letter codes used in the end
                #
                
                code3Letter = self.getCode3Letter(sequenceElement)
                #print code3Letter

                for bioPolymerType in bioPolymerCodes.keys():
                  
                  # First check on 3 letter codes
                  if code3Letter in bioPolymerCodes[bioPolymerType][1]:
                    unknownSequenceElements[-1][-1].append(bioPolymerType)
                    #print "1 %s" % bioPolymerType
                  
                  # Then on 1 letter codes, EXCEPT for proteins - this should not happen
                  elif code3Letter in bioPolymerCodes[bioPolymerType][0] and bioPolymerType != 'protein':
                    unknownSequenceElements[-1][-1].append(bioPolymerType)
                    #print "0 %s" % bioPolymerType

                  # Finally a dictionary check...
                  else:
                    for tempCcpCode in chemCompStdDict[bioPolymerType].keys():
                      if code3Letter == tempCcpCode or code3Letter == chemCompStdDict[bioPolymerType][tempCcpCode][1]:
                        unknownSequenceElements[-1][-1].append(bioPolymerType)                     
                        #print "keys %s" % bioPolymerType

          #
          # Try to determine polymer type based on matches for unknown ones...
          #

          if unknownSequenceElements:
           
            bioPolymerMatches = {}
            for bioPolymerType in bioPolymerCodes.keys():
              bioPolymerMatches[bioPolymerType] = 0

            preferredBioPolymer = None

            for (sequenceElementIndex,bioPolymerTypes) in unknownSequenceElements:            
              for bioPolymerType in bioPolymerTypes:
                bioPolymerMatches[bioPolymerType] += 1

            for bioPolymerType in nucleotideList:

              fraction = bioPolymerMatches[bioPolymerType] * 1.0 / len(unknownSequenceElements)
              
              if fraction == 1:
                molType = bioPolymerType
                print "  Molecule type for unknown residues set to %s (unambiguously)" % molType
                break

              elif fraction >= 0.9:
                molType = bioPolymerType
                print "  Partial match for unknown residues with molecule type %s (%.1f %%)" % (molType,fraction * 100)
                break
                
              elif fraction:
                if not preferredBioPolymer:
                  preferredBioPolymer = (fraction,bioPolymerType)
                else:
                  if fraction > preferredBioPolymer[0]:
                    preferredBioPolymer = (fraction,bioPolymerType)

            if not molType:

              bioPolymerType = 'protein'
              fraction = bioPolymerMatches[bioPolymerType] * 1.0 / len(unknownSequenceElements)

              if fraction == 1:
                molType = bioPolymerType
                print "  Molecule type for unknown residues set to %s (unambiguously)" % molType

              elif fraction >= 0.9:
                molType = bioPolymerType
                print "  Partial match for unknown residues with molecule type %s (%.1f %%)" % (molType,fraction * 100)

              elif fraction:
                if not preferredBioPolymer:
                  preferredBioPolymer = (fraction,bioPolymerType)

            #
            # Try to determine polymer type based on matches for unknown ones...
            #

            for (sequenceElementIndex,bioPolymerTypes) in unknownSequenceElements:
              if molType and molType in bioPolymerTypes:
                origSequenceDict[origMolName][1][sequenceElementIndex][-1] = molType

              elif len(bioPolymerTypes) == 1:
                origSequenceDict[origMolName][1][sequenceElementIndex][-1] = bioPolymerTypes[0]
              
              elif preferredBioPolymer and preferredBioPolymer[1] in bioPolymerTypes:
                origSequenceDict[origMolName][1][sequenceElementIndex][-1] = preferredBioPolymer[1]

              else:
                sequenceElement = origSequenceDict[origMolName][1][sequenceElementIndex][0]
                #print "Searching..."
                
                if hasattr(sequenceElement,'code3Letter') and sequenceElement.code3Letter:
                  
                  chemCompFound = False
                  
                  #
                  # First look in local directories, concentrate on protein,DNA,RNA,carbohydrate...
                  #
                  
                  for tmpMolType in bioPolymerCodes.keys() + ['carbohydrate']:
                    file = joinPath('molecule', 'ChemComp' , "%s.%s.xml" % (tmpMolType,sequenceElement.code3Letter.capitalize()))
                    for repository in repositoryList:
                      if os.path.exists(joinPath(repository.url.path,file)):
                        origSequenceDict[origMolName][1][sequenceElementIndex][-1] = tmpMolType
                        chemCompFound = True
                        molType = tmpMolType
                        print "Found in %s, file %s" % (repository.url.path,file)
                        break
                    if chemCompFound:
                      break
                  
                  #
                  # If nothing, *then* look in dictionary... and do special lookup for carbohydrates
                  # because codes are different for CCPN (using Iupac ones)
                  #
                  
                  if not chemCompFound:
 
                    for tmpMolType in bioPolymerCodes.keys() + ['carbohydrate']:
                      for tempCcpCode in chemCompStdDict[tmpMolType].keys():
                        if sequenceElement.code3Letter == tempCcpCode or sequenceElement.code3Letter == chemCompStdDict[tmpMolType][tempCcpCode][1]:
                          origSequenceDict[origMolName][1][sequenceElementIndex][-1] = tmpMolType
                          chemCompFound = True
                          molType = tmpMolType
                          #print "Found in %s in dictionary, molType %s" % (sequenceElement.code3Letter,tmpMolType)
                          break
                        
                    #
                    # Special carbohydrate lookup - TODO TODO MAKE SURE TO HANDLE MULTIPLE CIF CODES (FOR A/B!)
                    #
                    
                    tmpMolType = 'carbohydrate'
                    
                    for carbCcpCode in chemCompStdDict[tmpMolType].keys():
                      (code1letter,cifCode,name) = chemCompStdDict[tmpMolType][carbCcpCode]
                      if sequenceElement.code3Letter == cifCode:
                        origSequenceDict[origMolName][1][sequenceElementIndex][-1] = tmpMolType
                        chemCompFound = True
                        molType = tmpMolType
                        #print "Found in %s in lookuptable, molType %s" % (sequenceElement.code3Letter,molType)
                        break
                 
                  #
                  # Now do 'other'... much slower (many codes)
                  #
                 
                  if not chemCompFound:
                 
                    tmpMolType = 'other'
                    #print "LOOKING FOR OTHER"
                    file = joinPath('molecule', 'ChemComp' , "other.%s.xml" % sequenceElement.code3Letter.capitalize())
                    for repository in repositoryList:
                      if os.path.exists(joinPath(repository.url.path,file)):
                        origSequenceDict[origMolName][1][sequenceElementIndex][-1] = tmpMolType
                        chemCompFound = True
                        molType = tmpMolType
                        break
                  
                    #
                    # If nothing, *then* look in dictionary... and do special lookup for carbohydrates
                    # because codes are different for CCPN (using Iupac ones)
                    #

                    if not chemCompFound:

                      from ccp.general.ChemCompOverview import chemCompOverview
                      
                      for tempCcpCode in chemCompOverview[tmpMolType].keys():
                        if sequenceElement.code3Letter == tempCcpCode or sequenceElement.code3Letter == chemCompOverview[tmpMolType][tempCcpCode][1]:
                          origSequenceDict[origMolName][1][sequenceElementIndex][-1] = tmpMolType
                          molType = tmpMolType
                          chemCompFound = True

    #
    # Now have to find the correct chemComp for each (if possible), and 
    # user intervention if nothing found.
    #
    
    moleculeBonds = {}
    molSystemBonds = {} 
    
    origMolNames = origSequenceDict.keys()
    origMolNames.sort()
    
    refUrl = None
    
    for origMolName in origMolNames:
   
      #
      # If molecule already set by mapping, don't do it...
      #
      
      if createMoleculeDict.has_key(origMolName):
        continue
      
      chemCompDict[origMolName] = []
      oldResidueType = ''
      seqLength = len(origSequenceDict[origMolName][1])
      
      #
      # Option to ignore 'single' residue molecules of type 'other'
      #
      
      if self.ignoreSingleOther and seqLength == 1 and origSequenceDict[origMolName][1][0][1] == 'other':
        continue
      
      for seqElIndex in range(seqLength):
      
        (molSequenceElement,residueType) = origSequenceDict[origMolName][1][seqElIndex]

        seqCode = molSequenceElement.seqCode
        code3Letter = None
        code1Letter = None
        
        for sequence in origSequenceDict[origMolName][0]:
          
          sequenceElement = sequence.elements[seqElIndex]
          
          #
          # Track if any bonds are present...
          #
          
          if sequenceElement.bonds:
            for bondType in stableBondTypes: # Don't use hydrogen or salt bonds...
              
              if sequenceElement.bonds.has_key(bondType):
              
                atomNames = sequenceElement.bonds[bondType].keys()
                atomNames.sort()
                
                for atomName in atomNames:
                  for (bondedSeqEl,bondedAtomName) in sequenceElement.bonds[bondType][atomName]:
                    for tSeqEl in sequence.elements:
                      if bondedSeqEl == tSeqEl:
                        break
                    
                    if bondedSeqEl == tSeqEl:
                      bondDict = moleculeBonds
                    else:
                      bondDict = molSystemBonds
                    
                    if not bondDict.has_key(sequenceElement):
                      bondDict[sequenceElement] = {}
                    
                    if not bondDict[sequenceElement].has_key(atomName):
                      bondDict[sequenceElement][atomName] = []
                      
                    bondDict[sequenceElement][atomName].append((bondedSeqEl,bondedAtomName))
        
        #
        # Check for 1 and 3 letter code. For DNA/RNA the one-letter code is always used (hacked in above)
        #
      
        chemComp = None

        #
        # If mappings available, try to do that first (takes precedence)
        #
        
        #
        # WARNING TODO: currently mixup between uppercased and lowercased three letter codes!!!
        #
        
        if self.chemCompMappings:
          
          if hasattr(sequenceElement,'code3Letter'):
            code = sequenceElement.code3Letter
          elif hasattr(sequenceElement,'code1Letter'):
            code = sequenceElement.code1Letter
          
          if self.chemCompMappings.has_key(code):
            (mappedMolType,mappedCcpCode) = self.chemCompMappings[code]
            chemComp = self.getChemComp(mappedMolType,mappedCcpCode)
            
            if chemComp:
              residueType = mappedMolType
        
        #
        # Do the normal mapping if nothing from the mappings
        #
        
        if not chemComp:

          if residueType in ['protein'] + nucleotideList:

            if hasattr(sequenceElement,'code3Letter') and sequenceElement.code3Letter:

              #
              # DNA and RNA can be listed as three letter code (which is often one letter)
              #

              code3Letter = self.getCode3Letter(sequenceElement)

              chemComp = self.getChemComp(residueType,code3Letter, download = False)

              if chemComp:
                code1Letter = chemComp.code1Letter

              else:

                #
                # HACK HACK This is to reconvert DA, DG type codes...
                #

                if residueType in nucleotideList and code3Letter in bioPolymerCodes[residueType][0]:
                  codeIndex = bioPolymerCodes[residueType][0].index(code3Letter)
                  code3Letter = bioPolymerCodes[residueType][1][codeIndex]

                #
                # A full 3 letter code for DNA/RNA is not the same as the CCP code - reset some information in that case
                #

                if residueType in nucleotideList and code3Letter in bioPolymerCodes[residueType][1]:

                  codeIndex = bioPolymerCodes[residueType][1].index(code3Letter)
                  code1Letter = bioPolymerCodes[residueType][0][codeIndex]
                  # THIS IS A HACK!
                  sequenceElement.code1Letter = code1Letter

                else:
                  code = code3Letter
                  codeText = '3 letter'

            if not chemComp and hasattr(sequenceElement,'code1Letter') and sequenceElement.code1Letter:

              code1Letter = sequenceElement.code1Letter

              if code1LetterToCcpCodeDict.has_key(residueType) and code1LetterToCcpCodeDict[residueType].has_key(code1Letter):

                ccpCode = code1LetterToCcpCodeDict[residueType][code1Letter]

              else:

                ccpCode = code1Letter

              chemComp = self.getChemComp(residueType,ccpCode, download = False)

              if chemComp:
                code3Letter = chemComp.ccpCode
              else:
                code = code1Letter
                codeText = '1 letter'

          elif residueType == 'carbohydrate':

            #
            # Special carbohydrate lookup - TODO MAKE SURE CAN HANDLE MULTIPLE CIF CODES FOR a/b!
            #

            code3Letter = sequenceElement.code3Letter

            for carbCcpCode in chemCompStdDict['carbohydrate']: 
              (code1letter,cifCode,name) = chemCompStdDict['carbohydrate'][carbCcpCode]
              if sequenceElement.code3Letter == cifCode:
                code3Letter = carbCcpCode

            chemComp = self.getChemComp(residueType,code3Letter)

            if chemComp:
              code1Letter = None
            else:
              code = code3Letter
              codeText = '3 letter'

          # This is now changed in case there's an 'other' type residue in the sequence...
          elif residueType == 'other':

            code3Letter = sequenceElement.code3Letter

            chemComp = self.getChemComp(residueType,code3Letter)

            if chemComp:
              code1Letter = None
            else:
              code = code3Letter
              codeText = '3 letter'

          if not chemComp:
            
            if hasattr(sequenceElement,'code3Letter'):
              codeText = '3 letter'
              code = sequenceElement.code3Letter
              
              tryResidueType = None
              # One last try - could be 'other' type in protein, for example
              if hasattr(sequenceElement,'residueType') and sequenceElement.residueType:
                tryResidueType = sequenceElement.residueType
              elif residueType != 'other':
                tryResidueType = 'other'
              
              if tryResidueType:
                chemComp = self.getChemComp(tryResidueType,code)
             
            elif hasattr(sequenceElement,'code1Letter'):
              codeText = '1 letter'
              code = sequenceElement.code1Letter

        #
        # If required, auto create new chemComps. These will have minimal information, only atom names, no bonds, but
        # should be sufficient to link resonance and/or coordinate information.
        #
        
        if not chemComp and self.autoCreateChemComps and hasattr(sequenceElement,'atomNames') and sequenceElement.atomNames:
         
          if not residueType:
            if len(sequence.elements) == 1:
              ccMolType = 'other'
            elif seqElIndex != 0 and sequence.elements[seqElIndex - 1].seqCode != sequence.elements[seqElIndex].seqCode - 1:
              ccMolType = 'other'
            else:
              ccMolType = molType
          else:
            ccMolType = residueType
          
          if hasattr(sequenceElement,'code3Letter'):
            code = sequenceElement.code3Letter.capitalize()
           
          elif hasattr(sequenceElement,'code1Letter'):
            code = sequenceElement.code1Letter
          
          self.rawChemComp = GenericChemComp(code,ccMolType,sequenceElement.atomNames)
          
          for i in range(len(sequenceElement.atomNames)):
            atomName = sequenceElement.atomNames[i]
            atomSerial = i + 1
            
            self.rawChemComp.atoms.append(GenericAtom(atomSerial,atomName))

          self.rawChemCompList = [self.rawChemComp]
          
          #
          # Now use normal chemComp creation approach, fake initialisation
          #
          
          print drawBox("Warning: autocreating chemComp %s" % code, indent = "    ")
          
          # Track original variables, reset afterward
          allowPopups = self.allowPopups
          convertCount = self.convertCount
          mainCode = self.mainCode
          
          # Initialise fake stuff
          self.initConvertCount('chemComp')          
          self.allowPopups = self.createAtomSets = False
          self.saveChemComp = True
          self.molTypes = self.ccpCodes = self.chemComps = []
          self.makeNamingSystem = "CUSTOM"
          
          self.setupChemCompInfo()
          self.makeChemComp()

          # Reset original info
          self.allowPopups = allowPopups
          self.convertCount = convertCount
          self.mainCode = mainCode
          
          #
          # Set new chemComp
          #
          
          chemComp = self.chemComps[0]

        #
        # If nothing found or created, then ask user which chemComp
        #
        
        if not chemComp and not self.ignoreUnknownChemComps and self.allowPopups:
            
          selectNew = True
          if not self.minimalPrompts:
            selectNew = self.messageReporter.showYesNo("Select a chemComp","Unrecognized %s chemComp %s code %s for position %d in %s. Do you want to select a chemComp for this residue?" % (residueType,codeText,code,seqCode,self.fileName), parent = self.guiParent)

          if selectNew:
          
            interaction = self.multiDialog.ChemCompSelection(self.guiParent,self.project, chemCompEntries = ['ChemComp'], origCode = code)

            if interaction and interaction.frame.getSelectedChemComp():

              chemComp = interaction.frame.getSelectedChemComp()
              residueType = chemComp.molType
                              
        #
        # Send warning if no chemComp selected - will be ignored
        #
        
        if not chemComp:
          if not residueType:
            residueTypeText = molType
          else:
            residueTypeText = residueType
          
          if molType in ("protein","DNA","RNA"):
            # Insert an unknown chemComp in these cases to avoid messing up sequence numbering!
            chemComp = self.getChemComp(molType,"Xxx")
            
          if not chemComp:
            print "  Warning: no %s chemComp %s code %s selected. Will be ignored." % (residueTypeText,codeText,code)      
            self.addWarning("No CCPN reference compound could be found for molecular type %s, compound %s with code %s." % (residueType,codeText,code))
            continue

          else:
            print "  Warning: No %s chemComp %s code %s selected, replaced by Xxx chemComp to maintain chain sequence numbering." % (residueTypeText,codeText,code)  
    
        #
        # Keep track of everything...
        #
                
        if oldResidueType != residueType:
            
          if seqElIndex == 0 or not chemCompDict[origMolName]:
            createNew = True
            newResidueType = residueType
          
          else:

            #
            # If there's also a break in seqCode, then start a new segment
            #

            createNew = False

            if seqLength > 1:
              prevSequenceElement = sequence.elements[seqElIndex-1]
              #print oldResidueType, residueType, seqElIndex, prevSequenceElement.seqCode, sequenceElement.seqCode
              if sequenceElement.seqCode - 1 != prevSequenceElement.seqCode:
                createNew = True

            if not createNew:
              newResidueType = 'other'
          
          if createNew:
            chemCompDict[origMolName].append((newResidueType,[]))
          else:
            chemCompDict[origMolName][-1] = (newResidueType,chemCompDict[origMolName][-1][1])        
          
          #
          # If not, then keep on linking...
          #
          
          oldResidueType = residueType
          
        chemCompDict[origMolName][-1][-1].append((sequenceElement,chemComp))

    """
    #
    # Checking code
    #

    for origMolName in origMolNames:
    
      print origMolName
      
      for consistentStretch in chemCompDict[origMolName]:
        
        print "  %s: " % consistentStretch[0],
        
        seqTexts = []
        
        if len(consistentStretch[1]) > 9:
          seqIndexes = [0,1,2,'...',-3,-2,-1]
        
        else:
          seqIndexes = range(0,len(consistentStretch[1]))
          
        for seqIndex in seqIndexes:
          if type(seqIndex) == type(''):
            seqTexts.append(seqIndex)
          else:
            ccHead = consistentStretch[1][seqIndex][1]
            seqTexts.append(ccHead.ccpCode)
            
        print string.join(seqTexts,'-')
    
    print
    """

    #
    # Now reorganize automatically or let the user do it (also selection of chemCompVars!!)
    #
    
    ambiguousSequence = False

    for origMolName in origMolNames:

      #
      # If molecule already set by mapping, don't do it...
      #
      
      if createMoleculeDict.has_key(origMolName):
        continue

      for consistentStretch in chemCompDict[origMolName]:

        if len(chemCompDict[origMolName]) > 1:
          if not self.autoLink:
            ambiguousSequence = True
        
        for consistentStretch in chemCompDict[origMolName]:
          if not consistentStretch[0]:
            # Unknown molType...
            if not self.autoLink:
              ambiguousSequence = True

    #
    # Try automatically...
    #
    
    if not ambiguousSequence and self.minimalPrompts:
      forceUserInteraction = False
    else:
      forceUserInteraction = True    

    #
    # Create default chemCompVar for each consistent stretch...
    #

    for origMolName in origMolNames:
      
      #
      # Possible that already set by mapping information...
      #
      
      if createMoleculeDict.has_key(origMolName):
        continue
      
      origSeqIndex = 0
      
      for sequenceIndex in range(0,len(chemCompDict[origMolName])):

        (molType,sequenceList) = chemCompDict[origMolName][sequenceIndex]
        
        #
        # Set the molecule name... necessary in case sequence was split...
        #
        
        if sequenceIndex == 0:
          molName = origMolName
        else:
          addIndex = sequenceIndex + 1
          molName = "%s_%d" % (origMolName,addIndex)
          if len(sequenceList) == 0:
            seqEl = sequenceList[0][0]
            if hasattr(seqEl,'code3Letter') and seqEl.code3Letter:
              molName = seqEl.code3Letter
          
          while createMoleculeDict.has_key(molName):
            addIndex += 1
            molName = "%s_%d" % (origMolName,addIndex)

        seqLength = len(sequenceList)

        createMoleculeDict[molName] = [None,[]]
        for sequence in origSequenceDict[origMolName][0]:
          createMoleculeDict[molName][1].append((sequence,[]))

        for seqIndex in range(0,seqLength):
        
          if len(sequenceList) == 1:
            linking = 'none'
          elif seqIndex == 0 and not sequence.isCircular:
            linking = 'start'
          elif seqIndex == seqLength - 1 and not sequence.isCircular:
            linking = 'end'
          else:
            linking = 'middle'

          (seqElement,chemComp) = sequenceList[seqIndex]

          if chemComp:
          
            descriptorElements = {}
            origLinking = linking
            
            namingSystem = chemComp.findFirstNamingSystem(name = self.namingSystemName)
            
            if not namingSystem:
              for backupNamingSystemName in ('IUPAC','MSD','CUSTOM'):
                namingSystem = chemComp.findFirstNamingSystem(name = backupNamingSystemName)
                if namingSystem:
                  print "  Warning: selected %s naming system for chemComp %s,%s because %s naming system does not exist." % (backupNamingSystemName,chemComp.molType,chemComp.ccpCode,self.namingSystemName)
                  break
        
            if not namingSystem:
              print "  Warning: no %s naming system or backup naming systems found for chemComp %s,%s... aborting linking/descriptor selection." % (self.namingSystemName,chemComp.molType,chemComp.ccpCode)
              continue
            
            #
            # Also use bond information to set linking/descriptor!
            #
            
            linkedAtomNames = []
            
            for bondDict in (moleculeBonds,molSystemBonds):
            
              if bondDict.has_key(seqElement):
                
                atomNames = bondDict[seqElement].keys()
                
                #
                # Some doublechecking at this stage to get rid of atoms that are contained in prev,next,middle
                #
                
                for (molTypes,removeCcpnAtomNames) in ((('protein',),('C','N')),(('DNA','RNA'),('P',"O3'"))):
                  if chemComp.molType in molTypes:
                    for removeCcpnAtomName in removeCcpnAtomNames:
                      casn = findChemAtomSysName(namingSystem,{'atomName': removeCcpnAtomName})

                      if casn:
                        removeAtomNames = (casn.sysName,) + casn.altSysNames
                      else:
                        removeAtomNames = (removeCcpnAtomName,)
                        
                      for removeAtomName in removeAtomNames:
                        if removeAtomName in atomNames:
                          atomNames.pop(atomNames.index(removeAtomName))
                          del(bondDict[seqElement][removeAtomName])
                                    
                linkedAtomNames.extend(atomNames)
            
            if linkedAtomNames:
            
              linkedAtomNames.sort()
               
              if chemComp.molType in ('protein','DNA','RNA'):
                descriptorElements['link'] = linkedAtomNames
                
              else: # TODO: THIS SHOULD BE CARBOHYDRATE SPECIFIC?!?!
                linking = "link:%s" % ','.join(linkedAtomNames)
                                  
            chemCompVar = None
            
            #
            # Look for a sys name
            #
            
            if hasattr(seqElement,'formatCode') and seqElement.formatCode:
              chemCompSysNames = namingSystem.findAllChemCompSysNames(sysName = seqElement.formatCode)
              
              for chemCompSysName in chemCompSysNames:
                chemCompVar = chemCompSysName.findFirstChemCompVar(linking=linking, isDefaultVar=True) or \
                              chemCompSysName.findFirstChemCompVar(linking=linking)
                
                if chemCompVar:
                  break
                  
            #
            # If atom names are given, use those to try determine Var... bit slower!
            #
            
            if not chemCompVar and hasattr(seqElement,'atomNames') and seqElement.atomNames and namingSystem:
            
              chemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemComp,chemComp.chemAtoms,namingSystem.name)
              chemCompVars = list(chemComp.findAllChemCompVars(linking = linking))
                         
              if not chemCompVars and origLinking != linking:
                print "  Warning: could not find a chemCompVar for (%s,%s), linking '%s'. Switching to linking '%s'." % (chemComp.molType,chemComp.ccpCode,linking,origLinking)
                linking = origLinking
                chemCompVars = list(chemComp.findAllChemCompVars(linking = linking))                

              #
              # Check if can finetune the selection by looking at bonds as well
              #
              
              chemCompVars = getDescriptorElementMatches(chemCompVars,descriptorElements,verbose = self.verbose)
              
              #
              # Set up dictionaries for hits/misses
              #
              
              chemCompVarHits = {}
              chemCompVarMisses = {}
              
              for chemCompVar in chemCompVars: 
                chemCompVarHits[chemCompVar] = 0
                chemCompVarMisses[chemCompVar] = 0
              
              #
              # Also use nameMappings if available
              #
                    
              if self.addNameMappings and self.addNameMappings.has_key(chemComp.ccpCode):
                atomNameMappings = self.addNameMappings[chemComp.ccpCode]
              else:
                atomNameMappings = []
              
              #
              # Now loop over all the names and compare to the vars...
              #
              
              for sysAtomName in seqElement.atomNames:
                atomNames = []
                
                mappedSysAtomName = None
                for (origAtomName,newAtomName) in atomNameMappings:
                  if sysAtomName == origAtomName:
                    mappedSysAtomName = newAtomName
                    break
                
                for chemAtomSysName in chemAtomSysNames:
                  if sysAtomName == chemAtomSysName.sysName or sysAtomName in chemAtomSysName.altSysNames:
                    atomNames.append(chemAtomSysName.atomName)
                  elif mappedSysAtomName and mappedSysAtomName == chemAtomSysName.sysName or mappedSysAtomName in chemAtomSysName.altSysNames:
                    atomNames.append(chemAtomSysName.atomName)
                  
                for atomName in atomNames:
                  for chemCompVar in chemCompVars:
                    if chemCompVar.findFirstChemAtom(name = atomName):
                      chemCompVarHits[chemCompVar] += 1
                    else:
                      chemCompVarMisses[chemCompVar] += 1
                      
              topHit = 0
              selectedChemCompVars = None
              
              #
              # Prioritize selection: 1. highest number of atom hits
              #                       2. lowest number of atom misses
              #                       3. lowest number of atoms (???)
              #
              
              for chemCompVar in chemCompVars:
                if chemCompVarHits[chemCompVar]:
                  if chemCompVarHits[chemCompVar] > topHit:
                    selectedChemCompVars = [chemCompVar]
                    topHit = chemCompVarHits[chemCompVar]
                    lowMiss = chemCompVarMisses[chemCompVar]
                    numAtoms = len(chemCompVar.chemAtoms)
                  elif chemCompVarHits[chemCompVar] == topHit:
                    if chemCompVarMisses[chemCompVar] < lowMiss:
                      selectedChemCompVars = [chemCompVar]
                      lowMiss = chemCompVarMisses[chemCompVar]
                      numAtoms = len(chemCompVar.chemAtoms)
                    elif chemCompVarMisses[chemCompVar] == lowMiss:
                      if len(chemCompVar.chemAtoms) == numAtoms:
                        selectedChemCompVars.append(chemCompVar)
                      elif len(chemCompVar.chemAtoms) < numAtoms:
                        selectedChemCompVars = [chemCompVar]
                        numAtoms = len(chemCompVar.chemAtoms)
                  
              if selectedChemCompVars:
                if len(selectedChemCompVars) > 1:
                  defaultChemCompVars = []
                  for selectedChemCompVar in selectedChemCompVars:
                    if selectedChemCompVar.isDefaultVar == True:
                      defaultChemCompVars.append(selectedChemCompVar)

                  if defaultChemCompVars:
                    selectedChemCompVars = defaultChemCompVars
                    
                if len(selectedChemCompVars) > 1:
                  print "  Warning: %d chemCompVar '%s','%s' matches based on atomnames - picking first match" % (len(selectedChemCompVars),chemComp.ccpCode,linking)

                chemCompVar = selectedChemCompVars[0]
                  
            #
            # Otherwise pick based on presets, or neutral one, or only one available, depending...
            #
            
            if not chemCompVar:
            
              linkingChemCompVars = chemComp.findAllChemCompVars(linking = linking)
              
              if not linkingChemCompVars:
                # If all else fails, try this...
                linkingChemCompVars = chemComp.findAllChemCompVars(linking = 'none')

              if len(linkingChemCompVars) == 1:
                chemCompVar = linkingChemCompVars.pop()
            
              #
              # Check if can finetune the selection by looking at bonds as well
              #
              
              if not chemCompVar:
            
                chemCompVars = getDescriptorElementMatches(linkingChemCompVars,descriptorElements,verbose = False)
              
                if chemCompVars:
                  for ccv in chemCompVars:
                    if ccv.isDefaultVar:
                      chemCompVar = ccv
                      break
                  
                  if not chemCompVar:
                    chemCompVars = list(chemCompVars)
                    chemCompVars.sort()
                    chemCompVar = chemCompVars.pop()

              #
              # Pick default?
              #
              
              if not chemCompVar:
                chemCompVar = chemComp.findFirstChemCompVar(linking = linking, isDefaultVar = True)
              
              #
              # Extra check in case it's a chain but 'other' type - so could be 'none' linking as well
              #
              
              if not chemCompVar:
                
                if molType == 'other':
                  chemCompVar = chemComp.findFirstChemCompVar(linking = 'none', descriptor = 'neutral')
                  if not chemCompVar:
                    linkingChemCompVars = chemComp.findAllChemCompVars(linking = linking)

                    if linkingChemCompVars:

                      if len(linkingChemCompVars) == 1:
                        chemCompVar = linkingChemCompVars[0]
                      else:
                        chemCompVar = chemComp.findFirstChemCompVar(linking = linking, isDefaultVar = True)
                
                #
                # Try neutral one
                #

                if not chemCompVar:
                  chemCompVar = chemComp.findFirstChemCompVar(linking = linking, descriptor = 'neutral')                     

                #
                # If still nothing found... then force user interaction
                #

                if not chemCompVar:
                  forceUserInteraction = True

            if not chemCompVar and not forceUserInteraction:
              print "  ERROR: no chemCompVar %s,%s,%s found while automatically creating molecule - need user interaction." % (chemComp.molType,chemComp.ccpCode,linking)
              forceUserInteraction = True

            for origSequenceNum in range(0,len(origSequenceDict[origMolName][0])):
              sequence = origSequenceDict[origMolName][0][origSequenceNum]
              seqElement = sequence.elements[origSeqIndex]              
              createMoleculeDict[molName][1][origSequenceNum][1].append((seqElement,chemComp,chemCompVar))
          
          else:
            
            resNameText = ""
            
            if hasattr(seqElement,'code3Letter'):
              resNameText += ", resName %s" % seqElement.code3Letter

            print "  Warning: no chemComp selected for %s, seqCode %s%s..." % (origMolName,seqElement.seqCode,resNameText)
          
            for origSequenceNum in range(0,len(origSequenceDict[origMolName][0])):
              sequence = origSequenceDict[origMolName][0][origSequenceNum]
              seqElement = sequence.elements[origSeqIndex]
              createMoleculeDict[molName][1][origSequenceNum][1].append((seqElement,None,None))
              
          origSeqIndex += 1

    #
    # Apply numIdenticalChains, if possible
    # NOTE: will not remove sequences from original list...
    #
    
    if len(origMolNames) == 1 and self.numIdenticalChains:
      sequences = createMoleculeDict[molName][1]
      numSequences = len(sequences)
      (sequence,sequenceList) = sequences[0]
      for i in range(numSequences,self.numIdenticalChains):
        createMoleculeDict[molName][1].append((sequence,sequenceList[:]))
    
    #
    # From here on, only createMoleculeDict and molName is used!
    #
    
    #
    # User interaction to determine exact situation...
    #

    if forceUserInteraction and self.allowPopups:
      
      popup = self.multiDialog.SequenceLink(self.guiParent,self.project,createMoleculeDict)
      
      createMoleculeDict = popup.createMoleculeDict
              
    #
    # TODO: additional check to see whether there's no None? If so,
    # then ask whether to get rid of that chemComp or to abort...
    #

    #
    # Check whether sequence exists, ask whether want
    # to use this molecule setup...
    #

    self.checkMoleculesExists(createMoleculeDict)
    
    #
    # Now create the molecule(s)
    #
    
    molNames = createMoleculeDict.keys()
    molNames.sort()
    
    for molName in molNames:
      
      #
      # If already found, don't bother
      #
      
      if createMoleculeDict[molName][0]:
        self.setMoleculeInfo(createMoleculeDict[molName][0],createMoleculeDict[molName],newMolecule = False)
        continue

      #
      # This is fine... go ahead and make the molecule
      #
      # Try to rename the molecule based on the chemComp if there's only one molResidue...
      #

      if molName == '':
        newMolName = 'Molecule1'
      else:
        newMolName = molName
      
      if newMolName[:8] == 'Molecule' and len(createMoleculeDict[molName][1][0][1]) == 1:
        
        chemComp = createMoleculeDict[molName][1][0][1][0][1]
        if chemComp and chemComp.name:
          tmpMolName = chemComp.name
        
          if not self.project.findFirstMolecule(name = tmpMolName):
            newMolName = tmpMolName        
        
      while self.project.findFirstMolecule(name = newMolName):
        newMolName = self.askMoleculeName(newMolName)
        
      #
      # Set the mapping... can only do this if an original name available
      #
      
      if not molMappings and molName:
      
        molMapping = NameMapping.MoleculeMapping(self.nameMapping, name = molName, formatName = newMolName)
      
      #
      # Create new molecule
      # 

      molecule = self.createMolecule(newMolName,createMoleculeDict[molName])
      createMoleculeDict[molName][0] = molecule
      
      #
      # To set format specific data...
      #
      
      self.setMoleculeInfo(molecule,createMoleculeDict[molName])

      #
      # Create molResidues for this molecule
      #
      
      addText = self.createMolResidues(createMoleculeDict[molName],moleculeBonds)
      
      print "Created %smolecule %s (molType %s, %d chemComps)" % (addText,newMolName,molecule.molType,len(molecule.molResidues))
      
      #molecule.checkAllValid(complete = True)      

    return (createMoleculeDict,molSystemBonds)

  def getChemCompRepository(self):
  
    if self.saveChemCompLocal:
      repositoryName = 'userData'
    else:
      repositoryName = refRepositoryName
 
    return self.project.findFirstRepository(name=repositoryName)
    

  def getChemComp(self,residueType,residueCode, download = True):
    
    if residueType in ('other',):
      
      from ccp.general.ChemCompOverview import chemCompOverview
      
      chemCompInfoDict = chemCompOverview[residueType]
    else:
      chemCompInfoDict = chemCompStdDict[residueType]
    
    # get correct ccpCode
    
    # Make sure three letter codes are capitalized!!
    if len(residueCode) == 3:
      ccpCode = residueCode.capitalize()
    else:
      ccpCode = residueCode
    
    # If not available, see if can find a cifCode match...
    if not ccpCode in chemCompInfoDict.keys():
      for tempCcpCode in chemCompInfoDict.keys():
        if residueCode == chemCompInfoDict[tempCcpCode][1]:
          ccpCode = tempCcpCode
          break
    
    # get ChemComp:
    chemComp = self.project.findFirstChemComp(molType=residueType, 
                                              ccpCode=ccpCode)
    
    if not chemComp:
      # Do a doublecheck for specific molType classes if residueType is 'other'!
      if residueType == 'other' and not chemCompInfoDict.has_key(ccpCode):
        for tempMolType in chemCompStdDict.keys():
          # uses getByNavigation to avoid calling refreshTopObjects again
          chemComp = self.project.getByNavigation(('chemComps',
                                                  (tempMolType,ccpCode)))
          if chemComp:
            break
    
    if not chemComp:
      # this is partially redundant with the first 'findFirst'
      # check in chemCompArchiveDataDir and then check download
      showError = self.messageReporter.showError
      chemComp = getChemComp(self.project, residueType, ccpCode, 
                             download=download, showError=showError,
                             chemCompArchiveDir = self.chemCompPath)    
    #
    return chemComp
  
  def setMoleculeInfo(self,molecule,createMoleculeInfo, newMolecule = True):
    
    """
    Set format-specific information for molecule (application data...)
    """
    
    pass  
  
  def createMolecule(self,molName,createMoleculeInfo):
  
    """
    
    Make a molecule
    
    """
    molecule =  Molecule.Molecule(self.project, name = molName)
    
    self.convertCount[self.mainCode][1]['molecule(s)'] += 1
    
    for (sequence,sequenceList) in createMoleculeInfo[1]:

      if hasattr(sequence,'commonNames') and sequence.commonNames:
        for commonName in sequence.commonNames:
          if commonName:
            molecule.addCommonName(returnMemopsLine(commonName))
    
      if hasattr(sequence,'details') and sequence.details:
        if molecule.details == None:
          molecule.details = ''
        molecule.details = molecule.details + sequence.details
      
    return molecule
  
  def checkMolSystem(self):
  
    """
    
    Find out which molecular system the user wants to link the chain info
    to.
    
    """
    
    if self.minimalPrompts and len(self.project.molSystems) == 1:
      self.molSystem = self.project.findFirstMolSystem()
      
    if not self.molSystem and self.project.molSystems:
    
      (selectionList,selectionDict) = createSelection(self.project.sortedMolSystems())

      interaction = self.multiDialog.SelectionList(

                         self.guiParent,
                         selectionList,
                         selectionDict = selectionDict,
                         title = "Project '%s': " % self.project.name + 'Select molecular system',
                         text = 'Existing molecular system codes:',
                         dismissText = 'Create new molecular system',
                         urlFile = 'SelectMolSystem'

                         )

      if interaction.isSelected:
        self.molSystem = interaction.selection
    
    if not self.molSystem:
    
      if hasattr(self.sequenceFile,'molSystemCode') and self.sequenceFile.molSystemCode:
        molSysCode = self.sequenceFile.molSystemCode
      else:
        molSysCode = 'Molecular system'
          
      #
      # Create molSystem
      #
      # Only ask questions if relevant...
      #
    
      if self.minimalPrompts:
        molSysCodeSet = 1
      else:
        molSysCodeSet = 0

      while (1):
        
        if not molSysCodeSet:

          molSysCode = self.dataEntry.askString("Enter molecular system name","Molecular system name",molSysCode,self.guiParent)
        
        if self.project.findFirstMolSystem(code = returnMemopsWord(molSysCode)):
          self.messageReporter.showError("Error","Molecular system code '%s' already exists in this project - choose another" % returnMemopsWord(molSysCode), parent = self.guiParent)
          molSysCodeSet = 0
        
        elif molSysCode:
        
          self.molSystem = MolSystem.MolSystem(self.project,code = returnMemopsWord(molSysCode), name = molSysCode)
          
          self.convertCount[self.mainCode][0] += 1
          
          break
      
    self.setMolSystemInfo()

  def setMolSystemInfo(self):
    
    """
    Set format-specific information for molecular system (application data...)
    """
    
    pass  

  def checkChains(self, createMoleculeDict,molSystemBonds):
  
    """
    
    See if there are any chain matches. Use them if asked by user, otherwise create
    a new chain for this molSystem.
    
    """

    seqElToChainRes = {}
    
    molTypeOrderList = ('protein','DNA','RNA','DNA/RNA','carbohydrate','other')
    
    origMolNames = createMoleculeDict.keys()
    origMolNamesByMolType = {}

    #
    # Try to use the original chain codes as much as possible...
    #
    
    existingChainCodeList = []

    for origMolName in origMolNames:
      
      #
      # First look for the molTypes of the molecules - try to do proteins, ... first
      #

      molecule = createMoleculeDict[origMolName][0]
      molType = molecule.molType
      
      if molType == None:
        molTypeCount = {'total': 0}
        for molRes in molecule.molResidues:
          if not molTypeCount.has_key(molRes.molType):
            molTypeCount[molRes.molType] = 0
          molTypeCount[molRes.molType] += 1
          molTypeCount['total'] += 1
        if len(molTypeCount.keys()) > 2:
          for tMolType in molTypeOrderList:
            if molTypeCount.has_key(tMolType):
              if molTypeCount[tMolType] * 1.0 / molTypeCount['total'] > 0.9:
                molType = tMolType
                break
      
      if not origMolNamesByMolType.has_key(molType):
        origMolNamesByMolType[molType] = []
      origMolNamesByMolType[molType].append(origMolName)
      
      #
      # Now get the original chain codes
      #

      numSequences = len(createMoleculeDict[origMolName][1])

      for sequenceIndex in range(numSequences):
        (self.sequence,sequenceElements) = createMoleculeDict[origMolName][1][sequenceIndex]
        chainCode = None
        
        if self.sequence and hasattr(self.sequence,'chainCode'):
          chainCode = self.sequence.chainCode
          if not chainCode in existingChainCodeList:
            existingChainCodeList.append(chainCode)
    
    #
    # Remake the origMolNames list
    #
    
    origMolNames = []
    for molType in molTypeOrderList:
      if origMolNamesByMolType.has_key(molType):
        origMolNamesByMolType[molType].sort()
        origMolNames.extend(origMolNamesByMolType[molType])
        del(origMolNamesByMolType[molType])
        
    # For safety, also order anything that is left and not handled by molTypeOrderList
    for molType in origMolNamesByMolType.keys():
      origMolNamesByMolType[molType].sort()
      origMolNames.extend(origMolNamesByMolType[molType])
      
    
    #
    # Now start to create the chains, per molecule
    #    

    for origMolName in origMolNames:
    
      self.molecule = createMoleculeDict[origMolName][0]
      
      if not self.molecule.molResidues:
        print "  Error: molecule %s does not have any molResidues! Chain creation aborted..." % self.molecule.name
        self.addWarning("Molecule %s not created because its sequence elements could not be identified." % self.molecule.name)
        continue
        
      chainLoop = True  # to make sure findChainMatch is not called when multiple chain creation requested...
      numSequences = len(createMoleculeDict[origMolName][1])
      
      for sequenceIndex in range(numSequences):
      
        (self.sequence,sequenceElements) = createMoleculeDict[origMolName][1][sequenceIndex]
                
        validForChain = 0
        
        for sequenceElement in sequenceElements:
          if sequenceElement[1]:
            validForChain = 1
            break
        
        if not validForChain:
          print "Error: no information available for molecule %s! Ignoring..." % self.molecule.name
          continue
          
        chains = []
        if chainLoop:
          chains = self.findChainMatch() # NO PROBLEM

        if chains:

          self.chains.extend(chains)

        else:

          chainLoop = False
          newChain = True
          chainNum = 0

          while (newChain):

            #
            # Set the chainCode (use preset if available)
            #

            formatChainCode = None
            formatFirstSeqCode = sequenceElements[0][0].seqCode
            formatLastSeqCode = sequenceElements[-1][0].seqCode
            
            if hasattr(sequenceElements[0][0],'insertionCode'):
              formatSeqInsertCode = sequenceElements[0][0].insertionCode
            elif hasattr(sequenceElements[0][0],'seqInsertCode'):
              formatSeqInsertCode = sequenceElements[0][0].seqInsertCode
            else:
              formatSeqInsertCode = defaultSeqInsertCode

            chainMapping = None
            
            if self.sequence and hasattr(self.sequence,'chainCode'):

              formatChainCode = self.sequence.chainCode

              #
              # Chain mapping has to be exact!!
              #

              chainMapping = self.nameMapping.findFirstChainMapping(

                                      formatChainCode = formatChainCode,
                                      formatFirstSeqCode = formatFirstSeqCode,
                                      formatLastSeqCode = formatLastSeqCode,
                                      formatSeqInsertCode = formatSeqInsertCode

                                                                   )

              #
              # Set the value from mapping if available
              #

              if chainMapping:

                chainCode = chainMapping.chainCode

              else:

                chainCode = formatChainCode
                
            elif self.useOriginalChainCode and hasattr(self.sequence,'originalChainCode') and self.sequence.originalChainCode:
              
              chainCode = self.sequence.originalChainCode
              formatChainCode = self.sequence.originalChainCode

            elif self.chainCodes:

              chainCode = self.chainCodes[chainNum]

            else:

              chainCode = None

            #
            # Reset if already exists...
            #

            if self.molSystem.findFirstChain(code = chainCode):

              chainCode = None

            #
            # Set a default value if available
            #  

            if not chainCode or (chainCode == defaultMolCode and (not self.minimalPrompts or self.forceNonEmptyChainCode())):
            
              chainCode = None
              characters = ['',] + [x for x in chainCodeString]
              
              for char1 in characters:
                for char2 in characters:
                  for char3 in chainCodeString:
                    tempChainCode = char1+char2+char3
                    if tempChainCode in existingChainCodeList:
                      continue
                    if not self.molSystem.findFirstChain(code=tempChainCode):
                      chainCode = tempChainCode
                      break
                      
                  else:
                    continue
                  break
                   
                else:
                  continue
                break  
                    
              if not self.minimalPrompts and sequenceElements[0][1].ccpCode != 'Hoh':
                chainCode = self.askChainCode(chainCode)

            chainNum += 1
            
            #
            # If originally PDB, get the original chain code...
            #
            
            pdbOneLetterCode = self.getPdbOneLetterCode(chainCode)

            #
            # Create chain
            #

            self.chains.append(self.molSystem.newChain(molecule = self.molecule, code = chainCode, pdbOneLetterCode = pdbOneLetterCode))
            
            self.convertCount[self.mainCode][1]['chain(s)'] += 1
            self.convertCount[self.mainCode][1]['chain residue(s)'] += len(self.chains[-1].residues)
            
            #
            # Set format specific chain information
            #
            
            self.setChainInfo(self.chains[-1])

            #
            # Now two situations possible: either use the mapping to set seqCodes, or
            # start from scratch...
            #

            seqIdToCodeMapping = {}
            seqLen = len(self.chains[-1].residues)

            if chainMapping:

              #
              # Set the seqCodes from the mapping...
              #
              
              chainMappings = self.nameMapping.findAllChainMappings(formatChainCode = formatChainCode)

              for chainMapping in chainMappings:

                seqIdToCodeMapping[chainMapping.firstSeqId] = (chainMapping.formatFirstSeqCode,chainMapping.formatSeqInsertCode)
              
              firstSeqCode = chainMapping.formatFirstSeqCode
              lastSeqCode = '?'

            else:

              #
              # Set the 'difference' between the seqCode in the format file and in the Data Model
              #

              if self.seqCodeStart:
                seqCodeDiff = returnInt(self.seqCodeStart) - sequenceElements[0][0].seqCode

              else:
                seqCodeDiff = 0

              #
              # Now set seqIdToCodeMapping, which tracks if there's any 'jumps'...
              #

              oldSeqCode = None
              oldSeqInsertCode = None
              
              firstSeqCode = None
                            
              for resPos in range(seqLen):

                #
                # Get the seqCode from the file
                #

                residue = self.chains[-1].findFirstResidue(seqId = (resPos + 1))
                
                if not residue:
                  print "  Warning: residue information missing for seqId %d in chain '%s'! Ignored." % (resPos+1,self.chains[-1].code)
                  continue
                  
                seqEl = sequenceElements[resPos][0]

                if not seqEl.seqCode:
                  seqEl.seqCode = 0

                seqCode = seqEl.seqCode + seqCodeDiff
                if hasattr(seqEl,'insertionCode'):
                  seqInsertCode = seqEl.insertionCode
                elif hasattr(seqEl,'seqInsertCode'):
                  seqInsertCode = seqEl.seqInsertCode
                else:
                  seqInsertCode = defaultSeqInsertCode

                if oldSeqCode == None or oldSeqCode + 1 != seqCode or seqInsertCode != oldSeqInsertCode:
                  seqIdToCodeMapping[residue.seqId] = (seqCode,seqInsertCode)

                oldSeqCode = seqCode
                oldSeqInsertCode = seqInsertCode
                
                if firstSeqCode == None:
                  firstSeqCode = seqCode
            
              lastSeqCode = seqCode

            #
            # Now set the seqCodes for the residues in this chain, based on seqIdToCodeMapping
            #

            for resPos in range(0,seqLen):

              residue = self.chains[-1].findFirstResidue(seqId = (resPos + 1))
              
              if residue:
                if seqIdToCodeMapping.has_key(residue.seqId):
                  seqCode = seqIdToCodeMapping[residue.seqId][0]
                  seqInsertCode = seqIdToCodeMapping[residue.seqId][1]

                residue.seqCode = seqCode
                residue.seqInsertCode = seqInsertCode
                
              seqCode += 1
              
              seqEl = sequenceElements[resPos][0]
              
              if not seqElToChainRes.has_key(seqEl):
                seqElToChainRes[seqEl] = []
              
              seqElToChainRes[seqEl].append(residue)
                            
            print "Created chain '%s', start seqCode %s, end seqCode %s, molecule '%s'..." % (self.chains[-1].code,str(firstSeqCode),str(lastSeqCode), self.molecule.name)

            #
            # TODO: create a resonanceGroup by default?!
            #

            #
            # Set mapping... only done the first time a sequence with a particular
            # formatChainCode is read!
            #
            
            if formatChainCode and not chainMapping:

              #
              # Set for each 'chunk' of sequential seqCodes
              #

              seqIds = seqIdToCodeMapping.keys()
              seqIds.sort()
              
              for firstSeqId in seqIds:

                (formatFirstSeqCode,formatSeqInsertCode) = seqIdToCodeMapping[firstSeqId]
                
                chainMapping = self.nameMapping.newChainMapping(
                                    chainCode = chainCode,
                                    formatChainCode = formatChainCode,
                                    formatFirstSeqCode = formatFirstSeqCode,
                                    formatLastSeqCode = formatLastSeqCode,
                                    firstSeqId = firstSeqId,
                                    formatSeqInsertCode = formatSeqInsertCode

                                                     )
                
                if hasattr(self.sequence,'originalChainCode'):
                  chainMapping.addApplicationData(Implementation.AppDataString(application = ccpNmr_kw, keyword = 'originalFormatChainCode', value = self.sequence.originalChainCode))

            #
            # Then check if new chain needs to be created based on the current one...
            #

            if numSequences > sequenceIndex:
              if numSequences > (sequenceIndex+1) and self.sequence == createMoleculeDict[origMolName][1][sequenceIndex+1][0]:
                newChain = False
              elif sequenceIndex > 0 and self.sequence == createMoleculeDict[origMolName][1][sequenceIndex-1][0]:
                newChain = False

            if newChain:
              newChain = self.duplicateChain()
  
    #
    # Now set the linkEnds on the MolSystem level...
    #
    
    self.setMolResLinks(MolSystem,molSystemBonds,self.molSystem,seqElToChainRes)
    
    #
    # Set secondary structure information as application data on chain/residue
    # level. This is 'generic' secondary structure info...
    #
    
    self.setResidueSecStrucInfo(seqElToChainRes)
    
  def forceNonEmptyChainCode(self):
      
    return False
    
  def setResidueSecStrucInfo(self,seqElToChainRes):
  
    """
    Set info on secondary structure (still preliminary!)
    """
    
    for seqEl in seqElToChainRes.keys():
      if seqEl.secStrucInfo:
        #(secStrucType,secStrucSerial,specificInfo) = seqEl.secStrucInfo
        for residue in seqElToChainRes[seqEl]:
          residue.addApplicationData(Implementation.AppDataString(application = self.format,keyword = secStrucInfo_kw,value = str(seqEl.secStrucInfo)))
      
  def setChainInfo(self,chain):
    
    """
    Set format-specific information for chain (application data...)
    """
    
    pass  

  def duplicateChain(self):
    
    if not self.minimalPrompts:
      molecule = self.chains[-1].molecule
      numChains = len(self.molSystem.findAllChains(molecule = molecule))
      newChain = self.messageReporter.showYesNo(
      
                       "Create another chain",
                       "Do you want to create another chain for molecule '%s' in molecular system '%s' (already %d chain(s) present)?" % (molecule.name,self.molSystem.code,numChains)
                       
                       )

    else:
      newChain = False
      
    return newChain
    
  def checkSpinSystems(self):
  
    """
    
    Check if there are any resonanceGroups (spin systems) defined in the
    file - set them if so (and link to correct chain(s))
    
    """
    
    if self.sequence and hasattr(self.sequence,'spinSystems'):
      for spinSystem in self.sequence.spinSystems:
        resGroup = Nmr.ResonanceGroup(self.project.currentNmrProject)
        
        #
        # Also create link to chain so resGroups accessible from there
        #
        
        resGroup.nmrChains = self.chains
        
        if hasattr(spinSystem,'code'):
          resGroup.addApplicationData(Implementation.AppDataInt(application=self.format,keyword='code',value=spinSystem.code))

  def createSequence(self):
  
    """
    
    Set the format sequence. Always format specific
    
    """

    pass

  def setSequenceFileElements(self):
  
    """
    
    Set the format sequence file sequence elements. Always format specific
    
    """

    pass
  
  def checkMeasurementListExists(self):

    """
    
    Check whether contents of measurement file already contained in resonances...
    need to do this because, for example, XEasy shift data often comes in
    .prot/.peaks file pairs
    
    """
  
    self.existingMeasurementList = None
  
    resonanceMatches = {}
    
    #
    # Assuming that resonances listed only once in this file...
    #
    
    for measurement in self.measurementFileValues:
    
      #
      # Ignore if no info available (yet)...
      #
      
      if not hasattr(measurement,'value') or measurement.value is None:
        continue

      #
      # Get the resName(s)
      #
      
      resNames = []
      for resIndex in range(1,self.measurementResNum + 1):
        if self.measurementResNum == 1:
          keywds = {}
        else:
          keywds = {'resIndex': resIndex}
        
        resName = self.getMeasurementResName(measurement,**keywds)
    
        if resName:
          resNames.append(resName) 
    
      if not len(resNames) == self.measurementResNum:
        continue
      
      #
      # Check if this resonance is already linked to an existing MeasurementList
      #
      
      resMeasurementsFound = []
      for resName in resNames:
        resMeasurement = None
        
        if self.resonanceNames.has_key(resName):
          for resonance in self.resonanceNames[resName]:
          
            # TODO: this is not good enough...
            measFindFunc = getattr(resonance,"findAll%ss" % self.measurementApiCode)
            resMeasurements = measFindFunc(value = measurement.value)
            
            # Make sure if multiple resonances that find the SAME measurement!!
            for resMeasurement in resMeasurements:
              measurementExists = False
              for resMeasurementFound in resMeasurementsFound:
                if resMeasurement in resMeasurementFound:
                  measurementExists = True
                  resMeasurementFound.append(resMeasurement)
                
              if not measurementExists:
                resMeasurementsFound.append([resMeasurement])

      #
      # Determine matches, per list
      #
      
      for resMeasurementFound in resMeasurementsFound:
              
        if len(resMeasurementFound) == self.measurementResNum:
      
          existingMeasurementList = resMeasurementFound[0].parentList
          if not resonanceMatches.has_key(existingMeasurementList):
            resonanceMatches[existingMeasurementList] = 0
          resonanceMatches[existingMeasurementList] += 1
 
    #
    # After counting number of occurences, see if there is a full match
    #
 
    for existingMeasurementList in resonanceMatches.keys():
      if resonanceMatches[existingMeasurementList] == len(self.measurementFileValues):
      
        if not self.minimalPrompts:
          if self.messageReporter.showYesNo("Use existing list","A %s list exists with exactly the same values as are being read in: use this list instead?" % self.measurementType):
            continue
            
        self.existingMeasurementList = existingMeasurementList
        return

  def setMeasurementList(self):
  
    """
    
    Create a new measurement list
    
    """
  
    if not self.measurementList:
    
      if not self.measurementListName:
        
        (path,self.measurementListName) = splitPath(self.measurementFile.name)
      
        if not self.minimalPrompts:

          newMeasurementListName = self.dataEntry.askString("Enter %s list name" % self.measurementType,"New %s list name" % self.measurementType,self.measurementListName,self.guiParent)
        
          if newMeasurementListName:
            self.measurementListName = newMeasurementListName
                
      #
      # Set keywords that are always relevant for the list
      #
      
      self.setMeasurementListGeneral()

      #
      # Set any keywords specific for a particular measurement type
      #
      
      self.setMeasurementListSpecific()
      
      #
      # Now create the list
      #

      if not self.measurementListKeywds.has_key('unit') or not self.measurementListKeywds['unit']:
        self.measurementListKeywds['unit'] = 'default'

      #print 'KEYWDS: [%s]' % self.measurementListKeywds

      self.measurementList = self.ApiMeasurementList(self.project.currentNmrProject, name = self.measurementListName, **self.measurementListKeywds)

      self.convertCount[self.mainCode][0] += 1

  def setDataSource(self):

    """
    
    Get datasource and experiment info (create new ones if necessary)
    These popups will always appear - unless dataSource passed in.
    
    """
    
    self.newExperiment = False
    self.newDataSource = False

    if not self.dataSource:
      
      #
      # Select experiment
      #
      
      self.experiment = None
      
      (experimentList,interactionText) = self.getSpecificExperimentList()

      if experimentList and self.guiParent:

        (selectionList,selectionDict) = createSelection(experimentList)

        if selectionList:
          
          for key,val in selectionDict.items():
            if val.name == self.specName:
              selected = key
              break
          else:
            selected = None
            
        
          interaction = self.multiDialog.SelectionList(

                             self.guiParent,
                             selectionList,
                             selectionDict = selectionDict,
                             title = "Project '%s': " % self.project.name + 'Select experiment for peak list',
                             topText = "Format %s, experiment name '%s'" % (self.format,self.specName),
                             text = interactionText,
                             dismissText = 'create new experiment',
                             urlFile = 'SelectExperiment',
                             selected=selected

                             )

          if interaction.isSelected:
            self.experiment = interaction.selection
      
      if self.experiment:

        #
        # Select a datasource
        #
        
        self.selectDataSource(numDim = self.numDim)
      
      #
      # Generic - creates a dataSource if none given, experiment if none given.
      #
      
      if not self.dataSource:
      
        if not self.experiment:
          
          #
          # Create new experiment (code from DataManager, CcpNmrFormatConverter)
          # Have to do this - loop until done so
          #
          
          while (not self.experiment):

            popup = self.multiDialog.ExperimentCreate(self.guiParent,expName = self.specName, numDim = self.numDim, project = self.project)

            if popup:
              self.experiment = popup.experiment
              self.newExperiment = True
        
        # 
        # Create default frequency dataSource for experiment
        #
        
        dataSourceName = self.format
        
        if not self.minimalPrompts:
          newDataSourceName = self.dataEntry.askString("Enter dataSource name","Enter name for processed dataSource (experiment %s)" % self.experiment.name,dataSourceName,self.guiParent)
          
          if newDataSourceName:
            dataSourceName = newDataSourceName
        
        #
        # Create the processed data source. Default values are set here, will be reset if values available from the peak list file
        #

        self.dataSource = createPpmFreqDataSource(self.experiment,dataSourceName,'processed',self.experiment.numDim)
        self.newDataSource = True
    
    else:
    
      self.experiment = self.dataSource.experiment
      
      if self.experiment.numDim != self.numDim:
        self.messageReporter.showWarning("Warning","Mismatch in number of dimensions for peak list and given experiment",self.guiParent)
    
    if not self.dataSource:
      raise self.FormatConverterError("No datasource selected!")
      # TODO: loop over datasource selection (different sub!!?)

  def editDataSourceInformation(self):
  
    if self.newDataSource or self.newExperiment and not self.minimalPrompts:
    
      # TODO: put in a question here if want to edit?
    
      popup = self.multiDialog.AcqProcParsEdit(self.guiParent,experiment = self.experiment)
  
  def setPeakDimOrder(self):
  
    self.peakDimOrder = {
       
       1: [0],
       2: [0,1],
       3: [0,2,1],
       4: [0,2,3,1],
       5: [0,2,3,4,1],
       6: [0,2,3,4,5,1],
       7: [0,2,3,4,5,6,1]
       
      }
       

  def getSpecificExperimentList(self):

    if self.numDim:
      #experimentList = self.project.currentNmrProject.findAllExperiments(numDim = self.numDim)
      experimentList = [x for x in self.project.currentNmrProject.sortedExperiments() 
                        if x.numDim==self.numDim]
      interactionText = "Existing %sD experiment:" % (self.numDim)

    else:
      experimentList = self.project.currentNmrProject.sortedExperiments()
      interactionText = "Existing experiments:"
    
    return (experimentList,interactionText)
    
  """
  
  Generic functions (some to go in general.apiFunctions?)
  
  
  """
  
  def fileExists(self,fileName,verbose = True):
  
    value = 1
    
    if not os.path.exists(fileName):
      if verbose:
        self.messageReporter.showError("Error","File %s does not exist - aborting." % fileName,parent = self.guiParent)
      value = 0
      
    return value

  def createResonance(self,resName,isotopeCode = 'unknown',fixed = False, origSource = None,resLabel=None):
    
    """
    
    Create new resonance (default version). Isotope code should be set by linkResonances
    
    """
    
    if not fixed:
      resonance = self.project.currentNmrProject.newResonance(isotopeCode = isotopeCode)
    else:
      resonance = self.nmrConstraintStore.newFixedResonance(isotopeCode = isotopeCode)
    
    #
    # Set the resonance name in application data. Also check if doesn't
    # exist already and reset if so (can only have one for each format,
    # key or things get out of hand!)
    #
    
    if isinstance(resName, str):

      appDataClass = 'AppDataString'

    elif isinstance(resName, int):

      appDataClass = 'AppDataInt'

    else:

      appDataClass = None
      print "Error: unknown type %s for resName" % str(type(resName))
      
    if appDataClass:

      setUniqueAppData(appDataClass,resonance,self.format,self.assignKeyword,resName)
      
    #
    # Also track residue name, if available
    #
    
    if not resLabel and origSource and hasattr(origSource,'resLabel') and origSource.resLabel:
      resLabel = origSource.resLabel
    
    if resLabel:
      setUniqueAppData('AppDataString',resonance,self.format,'origResLabel',resLabel)

    self.newResonances.append(resonance)

    return resonance

  #
  # Below subdefinitions of self.setMeasurementListGeneral
  #

  def setMeasurementListGeneral(self):
    
    """
    General keywords for measurement list, can be made format specific
    """
    
    self.measurementListKeywds = {}
    
    #
    # Any file details?
    #
    
    if hasattr(self.measurementFile,'details'):
      self.measurementListKeywds['details'] = self.measurementFile.details
    else:
      self.measurementListKeywds['details'] = ""
      
    #
    # Units?
    #

    if hasattr(self.measurementFile,'units'):
      self.measurementListKeywds['unit'] = self.measurementFile.units
    elif hasattr(self.measurementFile,'unit'):
      self.measurementListKeywds['unit'] = self.measurementFile.unit
    elif self.unit:
      self.measurementListKeywds['unit'] = self.unit

    #
    # Experiment and software links
    #
          
    self.setMeasurementListExperiments()
    
    self.setMeasurementListSoftware()
    
  def setMeasurementListExperiments(self):

    #
    # Link new list to experiments passed in.
    # If no experiments passed, select them (if present)
    #
    
    if self.experiments:
      
      self.measurementListKeywds['experiments'] = self.experiments
    
    else:
      
      if not self.minimalPrompts:
        
        setCurrentStore(self.project,'NmrProject')
        
        if self.project.currentNmrProject.experiments:
          
          (selectionList,selectionDict) = createSelection(self.project.currentNmrProject.sortedExperiments())

          interaction = self.multiDialog.MultiSelectionList(
          
                                        self.guiParent,
                                        selectionList,
                                        selectionDict = selectionDict,
                                        title = "Project '%s': " % self.project.name + 'Select experiments linked to %s lists' % self.measurementType,
                                        text = 'Existing experiments:',
                                        dismissText = 'none linked',
                                        urlFile = 'SelectExperimentToMeasurementList'
                                        
                                        )
  
          if interaction.isSelectedList:
            self.measurementListKeywds['experiments'] = interaction.isSelectedList[:]
    
  def setMeasurementListSoftware(self):

    if self.method:
      
      self.measurementListKeywds['method'] = self.method
    
    else:
      
      if not self.minimalPrompts:
        
        setCurrentStore(self.project,'MethodStore')
            
        if self.project.currentMethodStore.methods:
          
          (selectionList,selectionDict) = createSelection(self.project.currentMethodStore.sortedMethods())

          interaction = self.multiDialog.SelectionList(
          
                                        self.guiParent,
                                        selectionList,
                                        selectionDict = selectionDict,
                                        title = "Project '%s': " % self.project.name + 'Select software linked to %s lists' % self.measurementType,
                                        text = 'Existing methods:',
                                        dismissText = 'none linked',
                                        urlFile = 'SelectMethodToMeasurementList' # TODO - where is this?
                                        
                                        )
  
          if interaction.isSelected:
            self.measurementListKeywds['method'] = interaction.selection
    


  #
  # Below subdefinitions of self.setMeasurementListSpecific
  #
        
  def setHExchProtectionListSpecific(self):
    
    pass

  def setHExchRateListSpecific(self):
    
    pass

  def setShiftAnisotropyListSpecific(self):
    
    pass

  def setShiftListSpecific(self):
    
    pass
  
  def setT1ListSpecific(self):
    
    # Note that information that was passed in always get preference over what comes in from file!
    
    if self.coherenceType:
      self.measurementListKeywds['coherenceType'] = self.coherenceType
    elif hasattr(self.measurementFile,'coherenceType'):
      self.measurementListKeywds['coherenceType'] = self.measurementFile.coherenceType
      
    if self.specFreq:
      self.measurementListKeywds['sf'] = self.specFreq
    elif hasattr(self.measurementFile,'specFreq'):
      self.measurementListKeywds['sf'] = self.measurementFile.specFreq
  
  def setT2ListSpecific(self):
    
    # Same info
    #if self.coherenceType:
    #  self.measurementListKeywds['coherenceType'] = self.coherenceType
    #elif hasattr(self.measurementFile,'coherenceType'):
    #  self.measurementListKeywds['coherenceType'] = self.measurementFile.coherenceType
      
    if self.specFreq:
      self.measurementListKeywds['sf'] = self.specFreq
    elif hasattr(self.measurementFile,'specFreq'):
      self.measurementListKeywds['sf'] = self.measurementFile.specFreq
  
    if self.tempCalibMethod:
      self.measurementListKeywds['tempCalibMethod'] = self.tempCalibMethod
    elif hasattr(self.measurementFile,'tempCalibMethod'):
      self.measurementListKeywds['tempCalibMethod'] = self.measurementFile.tempCalibMethod

    if self.tempControlMethod:
      self.measurementListKeywds['tempControlMethod'] = self.tempControlMethod
    elif hasattr(self.measurementFile,'tempControlMethod'):
      self.measurementListKeywds['tempControlMethod'] = self.measurementFile.tempControlMethod

  def setT1RhoListSpecific(self):
    
    self.setT2ListSpecific()
    
  def setNoeListSpecific(self):

    if hasattr(self.measurementFile,'refValue') and self.measurementFile.refValue:
     self.measurementListKeywds['refValue'] = self.measurementFile.refValue
      
    if hasattr(self.measurementFile,'refDesc') and self.measurementFile.refDesc:
     self.measurementListKeywds['refDescription'] = self.measurementFile.refDesc

    #if hasattr(self.measurementFile,'refValueType') and self.measurementFile.refValueType:
    # self.measurementListKeywds['noeValueType'] = self.measurementFile.refValueType

    if hasattr(self.measurementFile,'specFreq') and self.measurementFile.specFreq:
      self.measurementListKeywds['sf'] = self.measurementFile.specFreq

  def setJCouplingListSpecific(self):
    
    if hasattr(self.measurementFile,'specFreq') and self.measurementFile.specFreq:
      self.measurementListKeywds['sf'] = self.measurementFile.specFreq

  def setRdcListSpecific(self):
    
    # Basically the same...
    self.setJCouplingListSpecific()
        
  #
  # Below subdefinitions of self.setMeasurementSpecific
  #
  
  def setSingleResonanceMeasurementSpecific(self,resName,resonance):
    
    if abs(self.measurement.value - self.rawMeasurement.value) > self.measurementAllowedDiff:
      
      interaction = self.multiDialog.MeasurementSelect(self.guiParent,resName,resonance,self.rawMeasurement.value,self.measurementList,self.measurementApiCode)
      self.measurement.value = interaction.selection

  def setHExchProtectionSpecific(self,resName,resonance):
    
    self.setSingleResonanceMeasurementSpecific(resName,resonance)

  def setHExchRateSpecific(self,resName,resonance):
    
    self.setSingleResonanceMeasurementSpecific(resName,resonance)

  def setShiftAnisotropySpecific(self,resName,resonance):
    
    self.setSingleResonanceMeasurementSpecific(resName,resonance)
  
  def setShiftSpecific(self,resName,resonance):
  
    if self.autoMergeShifts or self.minimalPrompts:
      
      if not self.autoMergeShifts:
        print "  Warning: merged chemical shift values %.3f and %.3f for resonance %s" % (self.measurement.value,self.rawMeasurement.value,resName)

      self.measurement.value = (self.measurement.value + self.rawMeasurement.value) /2
            
    else:
  
      self.setSingleResonanceMeasurementSpecific(resName,resonance)
      
  def setT1Specific(self,resName,resonance):
    
    self.setSingleResonanceMeasurementSpecific(resName,resonance)

  def setT1RhoSpecific(self,resName,resonance):
    
    self.setSingleResonanceMeasurementSpecific(resName,resonance)
    
  def setT2Specific(self,resName,resonance):
    
    self.setSingleResonanceMeasurementSpecific(resName,resonance)
    
  def setMultipleResonanceMeasurementSpecific(self,resNames,resonances):
    
    if self.rawMeasurement.value != None and abs(self.measurement.value - self.rawMeasurement.value) > self.measurementAllowedDiff:
      
      print "  Warning: different values for same %s measurement. Keeping first value." % self.measurementApiCode
  
  def setNoeSpecific(self,resNames,resonances):
     
    self.setMultipleResonanceMeasurementSpecific(resNames,resonances)
  
  def setJCouplingSpecific(self,resNames,resonances):
     
    self.setMultipleResonanceMeasurementSpecific(resNames,resonances)

  def setRdcSpecific(self,resNames,resonances):
     
    self.setMultipleResonanceMeasurementSpecific(resNames,resonances)
#
  # Below subdefinitions of self.createMeasurement
  #
  
  def createSingleResonanceMeasurement(self,resonance,rawMeasurement):
    
    """
    
    Create measurement linked to resonance (default version)
    
    """
    
    keywds = {}
    
    for (attrName,keywd) in (('error','valueError'),('figOfMerit','figOfMerit')):
      if hasattr(rawMeasurement,keywd):
        keywds[attrName] = getattr(rawMeasurement,keywd)
    
    if self.thisAssignmentValid():
      self.convertCount[self.mainCode][1][self.valueCode] += 1
      return self.ApiMeasurement(self.measurementList, value = rawMeasurement.value, resonance = resonance, **keywds)
    else:      
      print 'For the residue %s%s assignment probability %.2f is below threshold'\
      %(rawMeasurement.resLabel,rawMeasurement.seqCode,rawMeasurement.figOfMerit)
      return
    
    
  def createHExchProtection(self,resonance,rawMeasurement):
    
    self.createSingleResonanceMeasurement(resonance,rawMeasurement)
    
  def createHExchRate(self,resonance,rawMeasurement):
    
    self.createSingleResonanceMeasurement(resonance,rawMeasurement)
    
  def createShiftAnisotropy(self,resonance,rawMeasurement):
    
    self.createSingleResonanceMeasurement(resonance,rawMeasurement)
    
  def createShift(self,resonance,rawMeasurement):
    
    self.createSingleResonanceMeasurement(resonance,rawMeasurement)
  
  def createT1(self,resonance,rawMeasurement):
    
    self.createSingleResonanceMeasurement(resonance,rawMeasurement)
  
  def createT1Rho(self,resonance,rawMeasurement):
    
    self.createSingleResonanceMeasurement(resonance,rawMeasurement)
  
  def createT2(self,resonance,rawMeasurement):
    
    self.createSingleResonanceMeasurement(resonance,rawMeasurement)

  def createMultipleResonanceMeasurement(self,resonances,rawMeasurement):
    
    """
    
    Create measurement linked to multiple resonances (default version)
    
    """
    
    keywds = {}
    
    for (attrName,keywd) in (('error','valueError'),):
      if hasattr(rawMeasurement,keywd):
        keywds[attrName] = getattr(rawMeasurement,keywd)

    self.convertCount[self.mainCode][1][self.valueCode] += 1

    if not hasattr(rawMeasurement,'value') or rawMeasurement.value == None:
      
      """
      TODO: below bit has to become measurement dependent function!!
      """
      if hasattr(rawMeasurement,'valueMin') and rawMeasurement.valueMin != None:
        minValue = rawMeasurement.valueMin
      else:
        # TODO this always valid?
        minValue = 0.0

      if hasattr(rawMeasurement,'valueMax'):
        maxValue = rawMeasurement.valueMax
        
        if maxValue == None and minValue:
          maxValue = minValue + 4.0  # This is not so great! Need some kind of default error?
        
        value = (minValue + maxValue) / 2.0
        keywds['error'] = maxValue - value
      
      else:
        # Meaningless
        value = None
      
    else:
      value = rawMeasurement.value
  
    return self.ApiMeasurement(self.measurementList, value = value, resonances = resonances, **keywds)

  def createNoe(self,resonances,rawMeasurement):
    
    return self.createMultipleResonanceMeasurement(resonances,rawMeasurement)

  def createJCoupling(self,resonances,rawMeasurement):
    
    return self.createMultipleResonanceMeasurement(resonances,rawMeasurement)

  def createRdc(self,resonances,rawMeasurement):
    
    return self.createMultipleResonanceMeasurement(resonances,rawMeasurement)

  #
  # Other functions relevant for measurements
  #
  
  def getMeasurementResName(self,measurement,resIndex = None):

    argList = []
    for attrName in ('molCode','seqCode','atomName'):
      if resIndex:
        attrName = "%s%d" % (attrName,resIndex)
      argList.append(getattr(measurement,attrName))
      
    kwDict = {}
    for attrName in ('seqInsertCode','spinSystemId'):
      
      if resIndex:
        attrNameIndexed = "%s%d" % (attrName,resIndex)
      else:
        attrNameIndexed = attrName
        
      if hasattr(measurement,attrNameIndexed):
        kwDict[attrName] = getattr(measurement,attrNameIndexed)

    return getResName(*argList,**kwDict)

  #
  # Project level getter
  #

  def getProject(self):
  
    return self.project
  
  #
  # Function to check if molecule exists already
  #

  def checkMoleculesExists(self,createMoleculeDict):
  
    """
    
    Find out whether a molecule with the exact same sequence exists - ask user
    if wants to use this one (unless minimalPrompts is 1 - will automatically
    pick it)
    
    """
    
    molNames = createMoleculeDict.keys()
    molNames.sort()
    
    for molName in molNames:
    
      #
      # Continue if molecule already determined
      #
      
      (molecule,sequenceLists) = createMoleculeDict[molName]
      
      if molecule:
        continue
      
      #
      # Now start checking...
      #
      
      sequenceList = sequenceLists[0][1]
      sequenceChemComps = []
      sequenceChemCompVars = []

      for (seqElement,chemComp,chemCompVar) in sequenceList:
        sequenceChemComps.append(chemComp)
        sequenceChemCompVars.append(chemCompVar)

      for molecule in self.project.sortedMolecules():
        if len(sequenceChemComps) == len(molecule.molResidues):
          molSequenceChemComps = [molResidue.chemComp for molResidue in molecule.sortedMolResidues()]

          if sequenceChemComps == molSequenceChemComps:

            fullMatch = True

            molSequenceChemCompVars = [molResidue.chemCompVar for molResidue in molecule.sortedMolResidues()]

            if sequenceChemCompVars != molSequenceChemCompVars:
              fullMatch = False

            if not fullMatch or not self.minimalPrompts:

              #
              # Set additional text...
              #

              if not fullMatch:
                matchText = '(Warning: not all chemCompVars match for this sequence)\n'
              else:
                matchText = ''

              seqTexts = []

              if len(sequenceChemComps) > 9:
                seqIndexes = [0,1,2,'...',-3,-2,-1]
              else:
                seqIndexes = range(0,len(sequenceChemComps))

              for seqIndex in seqIndexes:
                if isinstance(seqIndex, str):
                  seqTexts.append(seqIndex)
                else:
                  seqTexts.append(sequenceChemComps[seqIndex].ccpCode)

              seqText = string.join(seqTexts,'-')

              useMolecule = self.messageReporter.showYesNo("Use molecule","Sequence:\n%s\ncorresponds to molecule template with name '%s'.\n%sUse this molecule template for sequence setup?" % (seqText,molecule.name,matchText), parent = self.guiParent)

              if useMolecule:
                createMoleculeDict[molName][0] = molecule

            else:
              createMoleculeDict[molName][0] = molecule

  def askMoleculeName(self,molName):
    
    if self.minimalPrompts == 0:
    
      molNameExists = 1
    
      while (molNameExists):
        molName = self.dataEntry.askString('Enter molecule name','Molecule template name',molName,self.guiParent)
      
        if self.project.findFirstMolecule(name = molName):
          self.messageReporter.showError("Error","Molecule template name '%s' already exists in this project - choose another" % molName, parent = self.guiParent)
        elif molName:
          molNameExists = 0
          
    else:
    
      tempMolName = molName
      cnt = 1
      
      while (self.project.findFirstMolecule(name = tempMolName)):       
        
        tempMolName = molName + str(cnt)
        cnt += 1
      
      molName = tempMolName
      
    return molName

  def createMolResidues(self,createMoleculeInfo,moleculeBonds):
    
    addText = ""
    
    molecule = createMoleculeInfo[0]
    sequenceList = createMoleculeInfo[1][0][1]
    numRes = len(sequenceList)
    
    seqElToMolRes = {}
  
    for resPos in range(0,numRes):
    
      #
      # Get the info...
      #
      
      (seqElement,chemComp,chemCompVar) = sequenceList[resPos]
    
      #
      # Create a molResidue
      #
      
      if chemComp and chemCompVar:

        molRes = molecule.newMolResidue(seqCode = resPos + 1, chemComp = chemComp, linking = chemCompVar.linking, descriptor = chemCompVar.descriptor)

        self.convertCount[self.mainCode][1]['molecule residue(s)'] += 1
      
        if hasattr(seqElement,'authorSeqCode') and seqElement.authorSeqCode != None:
          keywds = {'application': self.format,
                    'keyword':     'authorSeqCode',
                    'value':       int(seqElement.authorSeqCode)}

          appData = Implementation.AppDataInt(**keywds)
          molRes.addApplicationData(appData)

        seqElToMolRes[seqElement] = [molRes]
        self.molResToSeqEl[molRes] = seqElement
        
        #
        # Set the linking information...
        #

        if chemCompVar.linking in ['middle','end']:

          prevLink = molRes.findFirstMolResLinkEnd(linkCode = 'prev')
          prevMolRes = molecule.findFirstMolResidue(serial = resPos)
          
          if prevLink and prevMolRes:

            nextLink = prevMolRes.findFirstMolResLinkEnd(linkCode = 'next')
            
            if nextLink:
            
              molResLink = Molecule.MolResLink(molecule,molResLinkEnds = [nextLink,prevLink])

              #
              # Set cis bond info...
              #

              if hasattr(seqElement,'hasCisPeptideBond') and seqElement.hasCisPeptideBond:
                molResLink.dihedralAngle = 0
      
      else:
        
        print " WARNING: missing chemComp or chemCompVar in molecule %s - not created." % molecule.name
    
    #
    # Set MolResLink if circular...
    #
    
    if createMoleculeInfo[1][0][0].isCircular:
      
      prevLink = molecule.sortedMolResidues()[0].findFirstMolResLinkEnd(linkCode = 'prev')
      nextLink = molecule.sortedMolResidues()[-1].findFirstMolResLinkEnd(linkCode = 'next')
      
      if prevLink and nextLink:
        molResLink = Molecule.MolResLink(molecule,molResLinkEnds = [nextLink,prevLink])
        addText = "circular "
        
    #
    # Set other molResBonds...
    #
    
    if moleculeBonds:
    
      self.setMolResLinks(Molecule,moleculeBonds,molecule,seqElToMolRes)        
    
    return addText
      
  def setMolResLinks(self,Package,bondDict,parent,seqElToRes):
  
    #
    # Code to set bonds on either Molecule or MolSystem level!
    #

    seqEls = bondDict.keys()
    seqEls.sort()

    for seqEl in seqEls:
      if seqEl in seqElToRes.keys():

        for residue in seqElToRes[seqEl]:
        
          if Package == MolSystem:
            chemComp = residue.chemCompVar.chemComp
          elif Package == Molecule:
            chemComp = residue.chemComp

          namingSystem = chemComp.findFirstNamingSystem(name = self.namingSystemName)
          
          if not namingSystem:
            print "  Warning: no %s naming system available for %s,%s - aborting molResLink creation." % (self.namingSystemName,chemComp.molType,chemComp.ccpCode)
            continue
          
          chemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemComp,chemComp.chemAtoms,self.namingSystemName)

          sysAtomNames = bondDict[seqEl].keys()
          sysAtomNames.sort()

          for sysAtomName in sysAtomNames:

            chemAtomSysName = None
            for tChemAtomSysName in chemAtomSysNames:
              if sysAtomName == tChemAtomSysName.sysName or sysAtomName in tChemAtomSysName.altSysNames:
                chemAtomSysName = tChemAtomSysName
                break

            if chemAtomSysName:

              chemAtomName = chemAtomSysName.atomName

              if Package == MolSystem:
                prevLink = residue.findFirstMolSystemLinkEnd(linkCode = chemAtomName)

                if not prevLink:
                  ccLinkEnd = residue.molResidue.chemComp.findFirstLinkEnd(linkCode = chemAtomName)

                  if ccLinkEnd:
                    molResLinkEnd = residue.molResidue.findFirstMolResLinkEnd(linkCode = chemAtomName)

                    if not molResLinkEnd:
                      prevLink = MolSystem.MolSystemLinkEnd(residue,linkCode = chemAtomName)
                    else:
                      print "  Error: could not set MolSystem LinkEnd '%s' for residue '%s.%d' because already exists on Molecule level!" % (chemAtomName,residue.chain.code,residue.seqCode)

              elif Package == Molecule:
                prevLink = residue.findFirstMolResLinkEnd(linkCode = chemAtomName)

              if prevLink:
                for (boundSeqEl,boundSysAtomName) in bondDict[seqEl][sysAtomName]:
                  if boundSeqEl in seqElToRes.keys():
                    for boundResidue in seqElToRes[boundSeqEl]:
                    
                      if Package == MolSystem:
                        boundChemComp = boundResidue.chemCompVar.chemComp
                      elif Package == Molecule:
                        boundChemComp = boundResidue.chemComp
                    
                      boundNamingSystem = boundChemComp.findFirstNamingSystem(name = self.namingSystemName)

                      if not boundNamingSystem:
                        print "  Warning: no %s naming system available for %s,%s - aborting molResLink creation." % (self.namingSystemName,boundChemComp.molType,boundChemComp.ccpCode)
                        continue
                      
                      boundChemAtomSysNames = findAllSysNamesByChemAtomOrSet(boundChemComp,boundChemComp.chemAtoms,self.namingSystemName)
                      
                      boundChemAtomSysName = None
                      for tChemAtomSysName in boundChemAtomSysNames:
                        if boundSysAtomName == tChemAtomSysName.sysName or boundSysAtomName in tChemAtomSysName.altSysNames:
                          boundChemAtomSysName = tChemAtomSysName
                          break

                      if boundChemAtomSysName:

                        boundChemAtomName = boundChemAtomSysName.atomName

                        if Package == MolSystem:
                          nextLink = boundResidue.findFirstMolSystemLinkEnd(linkCode = boundChemAtomName)
                          if not nextLink:
                            ccLinkEnd = boundResidue.molResidue.chemComp.findFirstLinkEnd(linkCode = boundChemAtomName)
                            
                            if ccLinkEnd:
                              nextLink = MolSystem.MolSystemLinkEnd(boundResidue,linkCode = boundChemAtomName)

                        elif Package == Molecule:
                          nextLink = boundResidue.findFirstMolResLinkEnd(linkCode = boundChemAtomName)

                        if nextLink:
                        
                          if Package == MolSystem:
                            objectName = 'MolSystemLink'
                            searchDictKey = 'molSystemLinkEnds'
                            infoString = 'MolSystem level '
                            residueInfo = "%s." % residue.chain.code
                            boundResidueInfo = "%s." % boundResidue.chain.code
                            
                          elif Package == Molecule:
                            objectName = 'MolResLink'
                            searchDictKey = 'molResLinkEnds'
                            infoString = 'Molecule %s level ' % residue.molecule.name
                            residueInfo = boundResidueInfo = ''
                          
                          searchFunc = "findFirst%s" % objectName
                          searchDict = {searchDictKey: [nextLink,prevLink]}
                          
                          bondLink = getattr(parent,searchFunc)(**searchDict)

                          if not bondLink:
                            searchDict[searchDictKey].reverse()
                            bondLink = getattr(parent,searchFunc)(**searchDict)

                            if not bondLink:
                              searchDict[searchDictKey].reverse()
                              try:
                                bondLink = getattr(Package,objectName)(parent,**searchDict)                               
                                print "  Set %sbond between (%s%s.%d,%s) and (%s%s.%d,%s)" % (infoString,residueInfo,residue.ccpCode,residue.seqCode,chemAtomName,boundResidueInfo,boundResidue.ccpCode,boundResidue.seqCode,boundChemAtomName)
                              except:
                                print "  Error: failed to set %sbond between (%s%s.%d,%s) and (%s%s.%d,%s)" % (infoString,residueInfo,residue.ccpCode,residue.seqCode,chemAtomName,boundResidueInfo,boundResidue.ccpCode,boundResidue.seqCode,boundChemAtomName)

                        else:
                          print "  Warning: could not set molResLink %s because does not exist for (%s,%s). ChemComp info probably needs updating!" % (boundChemAtomName,boundResidue.molType,boundResidue.ccpCode)

                      else:
                        print "  Warning: could not set molecular bond because no %s sysName found for %s (%s,%s)." % (self.namingSystemName,boundSysAtomName,boundResidue.molType,boundResidue.ccpCode)

                  else:
                    print "  Warning: could not set molecular bond because sequence element (%s,%d) not converted into molResidue." % (boundSeqEl.residueType, boundSeqEl.seqCode)

              else:
                print "  Warning: could not set molResLink %s because does not exist for (%s,%s). ChemComp info probably needs updating!" % (chemAtomName,residue.molType,residue.ccpCode)

            else:
              print "  Warning: could not set molecular bond because no %s sysName found for %s (%s,%s)." % (self.namingSystemName,sysAtomName,residue.molType,residue.ccpCode)

      else:
        print "  Warning: could not find reference dictionary information for (%s,%s) to set bonds." % (seqEl.seqCode,str(seqEl.formatCode))

  def getCode3Letter(self,sequenceElement):

    c3l = sequenceElement.code3Letter

    if c3l:
      c3l = c3l.strip().upper()
    
    return c3l

  def getResidueCode1Letter(self,residue, defaultCode = 'X'):
    
    code1Letter = residue.molResidue.chemComp.code1Letter
    
    if not code1Letter:
      code1Letter = defaultCode

    return code1Letter
  
  def getPdbOneLetterCode(self,currentChainCode):
    
    pdbOneLetterCode = self.getSpecificPdbOneLetterCode(currentChainCode)
     
    if len(currentChainCode) > 1:
      for newCode in chainCodeString:
        if not self.molSystem.findFirstChain(pdbOneLetterCode = newCode):
          pdbOneLetterCode = newCode
          break
      else:
        pdbOneLetterCode = ' '
    
    return pdbOneLetterCode

  def getSpecificPdbOneLetterCode(self,currentChainCode):
  
    return currentChainCode
    
  def findChainMatch(self):
    
    """
    
    Check if a chain with matching sequence exists and ask user whether he
    wants to link info to this one.
    (is this still necessary or should a new chain always be created?)
    
    """
    
    chains = []

    matchingChains = self.molSystem.findAllChains(molecule = self.molecule)

    if matchingChains:
    
      if not self.minimalPrompts:
      
        (selectionList,selectionDict) = createSelection(matchingChains)

        interaction = self.multiDialog.SelectionList(

                             self.guiParent,
                             selectionList,
                             selectionDict = selectionDict,
                             title = "Project '%s': " % self.project.name + 'Select chain',
                             text = 'Existing chain code:',
                             dismissText = 'create new chain',
                             urlFile = 'SelectChain'
                             
                             )
    
        if interaction.isSelected:
          chains.append(interaction.selection)
    
    return chains
    
  def askChainCode(self,chainCode):
  
    """
    
    Get a chain code from the user
    
    """
    
    chainCodeExists = 1
    
    while (chainCodeExists):
      chainCode = self.dataEntry.askString("Enter chain code","New chain code\n(molecule %s)" % self.molecule.name,chainCode,self.guiParent)
      
      if chainCode == "":
        chainCode = defaultMolCode
      
      if self.molSystem.findFirstChain(code = chainCode):
        self.messageReporter.showError("Error","Chain code '%s' already exists in this project - choose another" % chainCode, parent = self.guiParent)
      elif chainCode:
        chainCodeExists = 0
    
    return chainCode

  def createSoftwareDict(self):
    
    """
    
    Create a dict for the methods used in software packages
    
    """
  
    self.methods = {}
    
    if self.project.currentMethodStore:
      for software in self.project.currentMethodStore.software:
        self.methods[software.name] = {}
        for method in software.methods:
          self.methods[software.name][method.procedure] = method

    self.peakPickMethod = 'Automatic and manual'
    
  def createPeakList(self,peakListName):
  
    """
    
    Create a peak list in the data model
    
    """

    peakList = self.dataSource.newPeakList(isSimulated = False)
    peakList.name = peakListName[-80:]
    peakList.details = '%s format, file %s.' % (self.format,self.fileName)

    self.convertCount[self.mainCode][0] += 1
    
    self.peakList = peakList

    self.setFitMethodSerial()
    
  def setFitMethodSerial(self):
  
    self.peakList.fitMethodSerial = self.methods[self.format][self.peakPickMethod]
    
  def setRawPeakFileFormatSpecific(self):
    
    """
    
    For handling format specific things related to the peak list
    
    """
    
    pass
    
  def selectDataDimRefs(self):
  
    """
    
    Select dataDimRefs (for reading peaks and matching to peakDims)
    
    """

    if not self.dataDimRefs:
    
      #
      # Can only do this if there is a GUI environment - need to
      # pass in the info otherwise
      #
    
      if self.guiParent:
      
        #
        # Make a peak selection for reference info (chemical shifts)
        #

        rawPeaks = []
      
        for self.rawPeak in self.peakFile.peaks:
         
          if self.thisPeakValid():

            rawPeaks.append(self.rawPeak)
      
        #
        # Popup a selection window
        #
      
        while (not self.dataDimRefs):

          popup = self.multiDialog.DataDimRefSelectByPeakDim(self.guiParent, self.dataSource, rawPeaks, order = self.peakDimOrder)
      
          if popup.dataDimRefs:
            self.dataDimRefs = popup.dataDimRefs
            
            #
            # If new datasource and reference info available for spectrum from peak list file, set it here...
            #
            
            if self.newDataSource:
              self.setSpectrumInfoFromPeakList()
          
          else:
            self.messageReporter.showError("Error","Have to select dataDimRefs - retry",self.guiParent)
          
      else:

        self.messageReporter.showError("Error","Have to pass in dataDimRefs for reading peaks in non-GUI environment. Aborting.",parent = self.guiParent)
        return
        
        
    if self.numDim == None:
      self.numDim = len(self.dataDimRefs)

  def setSpectrumInfoFromPeakList(self):
    
    """
    Format dependent!
    """
    
    pass

  def setPeakListDimMapping(self):
  
    """
    
    Set the mapping between the ORIGINAL peak list dimension and the peak list
    dimension in the data model (which follows the same order as the experiment
    dimensions)
    
    """
  
    for i in range(0,len(self.dataDimRefs)):
      
      if not self.dataDimRefs[i]:
        continue
      
      dataDimRef = self.dataDimRefs[i]
      peakDimIndex = dataDimRef.dataDim.dim - 1
      
      mapping = str(i) + tagSep + str(peakDimIndex)
      
      setUniqueAppData('AppDataString',self.peakList,self.format,'dimMapping',mapping)

  def setPeaks(self):
  
    """
    
    Generic routine to loop over all valid peak objects from the format file,
    and set values in the data model.
    
    """
  
    for self.rawPeak in self.peakFile.peaks:
    
      if not self.thisPeakValid():
        continue
     
      # This to handle separated assignment items
      self.peakContribs = []
      
      self.setPeak()
      
      self.setPeakExtras()

      self.setPeakIntensity()
      
      self.setPeakDims()
   
  def checkPeakListValid(self):
   
    """
    
    Check whether this peak list has any valid peaks
    
    """
    
    validPeaks = 0
    
    for self.rawPeak in self.peakFile.peaks:
    
      validPeaks += self.thisPeakValid()
      
      if validPeaks:
        break
   
    return validPeaks
   
  def thisPeakValid(self):
  
    """
    
    Check whether this peak is valid. Default is true - is redefined for some formats.
    
    """
  
    return 1
  
  def thisAssignmentValid(self):
  
    """
    
    Check whether this assignment is valid. Default is true - is redefined for some formats.
    
    """
  
    return 1

  def setPeak(self):
  
    """
    
    Creates a peak object inside the data model or uses it to copy information
    onto. Also creates/checks application data with original peak number.
    
    """
    
    newPeak = 1
  
    if self.overwrite:
      
      if self.origPeakNumberDict.has_key(self.rawPeak.num):
        
        peak = self.origPeakNumberDict[self.rawPeak.num]

        self.peak = peak
        
        # TODO: why would this ever fail?
        if peak in self.existingPeaks:
          self.existingPeaks.pop(self.existingPeaks.index(peak))
        
        #
        # Get rid of intensities and peakcontribs
        # TODO: note that peakDims are kept... could in principle create problems with dataDimRefs
        #
        
        for peakIntensity in peak.peakIntensities:
          peakIntensity.delete()
        for peakContrib in peak.sortedPeakContribs():
          peakContrib.delete()
       
        newPeak = 0

    if newPeak:
  
      self.peak = self.peakList.newPeak(details = "Original number %d" % self.rawPeak.num)
      
      self.convertCount[self.mainCode][1]['peaks'] += 1
    
      #
      # Also set original number in applicationdata so can be reused for writing
      #

      self.peak.addApplicationData(Implementation.AppDataInt(application = self.format, keyword = peakNum_kw, value = self.rawPeak.num))
      
  def setPeakIntensity(self):
    
    """
    
    Set intensity objects belonging to peak. No default available - always
    format specific.
    
    """

    pass
      
  def setPeakDims(self):
  
    """
    
    Loop over dimensions for the peak from the format file
    
    """

    for self.rawPeakDimIndex in range(0,self.numDim):
    
      if self.rawPeakDimIndex >= len(self.dataDimRefs) or not self.dataDimRefs[self.rawPeakDimIndex]:
        continue
      
      self.setPeakDim()
      
      self.setPeakDimResonance()
     
  def deleteRemainingPeaks(self):
    
    """
    
    If overwrite mode is on, delete remaining existing peaks that were
    not in new list.
    
    """
  
    for existingPeak in self.existingPeaks:
      applData = existingPeak.findFirstApplicationData(application = self.format, keyword = peakNum_kw)
      if applData:
        addString = " (number %d)" % applData.value
      else:
        addString = ""
      print "  Warning: Deleting original peak%s" % addString
        
      existingPeak.delete()


  def setPeakDim(self):
    
    """
    
    Set peak dimension information for peak. No default available - always
    format specific.
    
    """

    pass
      
  def setPeakDimResonance(self):
  
    """
    
    Set the connection to (existing) resonances.
    
    Not assuming that the peak list is for one chain only... but if no chainCode
    is given by the external reader then assuming that it is all 'one chain'.
    Ambiguity then has to be handled by linkResonances, not by this script.
    
    Note that for XEasy this whole thing is different - have to use atomSerial or
    nothing at all
    
    """
    
    self.getPeakResNames()
    
    for resNameIndex in range(len(self.resNames)):
    
      self.resName = self.resNames[resNameIndex]
      
      #
      # Ignore if None
      #
      
      if self.resName == None:

        continue
    
      #
      # Link peakDim to resonance
      #

      if not self.resonanceNames.has_key(self.resName):
      
        #
        # Create resonance if doesn't exist.
        #
        
        if hasattr(self,'resLabels'):
          resLabel = self.resLabels[resNameIndex]
        else:
          resLabel = None
      
        resonance = self.createResonance(self.resName,resLabel=resLabel)
     
        #
        # Keep track of it...
        #

        resonances = [resonance]
        self.resonanceNames[self.resName] = [resonance]
    
      else:
    
        #
        # Resonance(s) exists, get info from self.resonanceNames
        #
      
        resonances = self.resonanceNames[self.resName]

      #
      # Link resonance(s) to peakDim via peakDimContrib
      # Note that will be linked to peakContrib automatically in this case:
      # no distinction is made between the contributions (by default)
      #
    
      keywds = {}
      if self.peakContribs:
        keywds['peakContribs'] = [self.peakContribs[resNameIndex]]

      for resonance in resonances:
        peakDimContrib = self.peakDim.findFirstPeakDimContrib(resonance = resonance)

        if not peakDimContrib:
          self.peakDim.newPeakDimContrib(resonance = resonance, **keywds)
        elif keywds:
          peakContribs = peakDimContrib.peakContribs
          for peakContrib in keywds['peakContribs']:
            if not peakContrib in peakContribs:
              peakDimContrib.addPeakContrib(peakContrib)
    
  def getPeakResNames(self):
  
    """
    
    Get the format name(s) of the 'resonance' linked to the peak in this dimension
    
    """
  
    self.resNames = []

  def setPeakExtras(self):
    
    """
    
    Set extra attributes belonging to peak. No default available - always
    format specific.
    
    """
    
    pass


  def createResonanceNamesDict(self,resonances):
  
    """
    
    Create a dict with all the resonance names (and corresponding resonances)
    
    """
    
    if self.resonanceMatching:
      self.resonanceNames = getApplResNames(self.format,resonances)
    else:
      self.resonanceNames = {}

  def setNumberPeakLists(self):
  
    """
    
    Determines which spectra to create.
    
      multiple spectra:  - if self.specName, select this one, ignore rest
                         - else loop over all spectra
    
      single spectrum:   - set datasource, ... if available
                         - else loop over spectrum

    """
    
    exitStatus = 1

    self.numPeakLists = len(self.peakFile.specNames)
    
    if self.numPeakLists > 1:
      
      #
      # If specName given, select only this one (if exists)
      #
      
      if self.specName and self.peakFile.specNames.count(self.specName) == 1:
        
        specIndex = self.peakFile.specNames.index(self.specName)
        
        #
        # Reset spectrum list for this file...
        #
        
        self.peakFile.specNames = [self.specName]
        self.peakFile.numDims = [self.peakFile.numDims[specIndex]]
        self.numPeakLists = 1
      
      else:
      
        #
        # Select from the list of spectrum names...
        #
        
        specNameList = self.peakFile.specNames

        if specNameList and not self.minimalPrompts and self.guiParent:
          
          textMessage = "The imported peak list has data for multiple experiments.\nSelect which data in this list "
          
          if self.dataSource:
            textMessage += "belongs to\nexperiment '%s', datasource '%s':" % (self.dataSource.experiment.name,self.dataSource.name)
          else:
            textMessage += "you want to convert to the data model."
        
          popup = self.multiDialog.MultiSelectionList(self.guiParent,specNameList,
                                                      title = "Project '%s': " % self.project.name +'Select experiments from peak list',
                                                      text = textMessage,
                                                      dismissText = 'Select none',
                                                      urlFile = 'SelectExperimentFromPeakList'
                                                      )

          if popup.isSelectedList:
          
            numDims = []

            for specName in popup.isSelectedList:
              specIndex = self.peakFile.specNames.index(specName)
              numDims.append(self.peakFile.numDims[specIndex])
        
            self.peakFile.specNames = popup.isSelectedList[:]
            self.peakFile.numDims = numDims
            self.numPeakLists = len(self.peakFile.specNames)
            
          else:
          
            exitStatus = 0
      
        #
        # Reset these for multiple spectra reading - will be set later on anyway
        #
    
        self.dataSource = None
        self.dataDimRefs = []
    
    return exitStatus
        
  def selectPeakLists(self):
  
    peakList = self.peakLists[0]
  
    self.peakLists = [peakList]
    
    # Do not handle SampledDataDims
    self.numPeakDim = 0
    for dataDim in peakList.dataSource.dataDims:
      if dataDim.className != 'SampledDataDim':
        self.numPeakDim += 1

  def peakListGetDimMapping(self):
    
    """
    
    This routine gets the mapping between the original format peak file dimensions
    and the dataDimRefs: am assuming that these are the same for all peak dimensions
    though. Only used if minimalPrompts is on.
    
    Eventually should have automatic split between peaks with different dataDimRefs?
    
    """
    
    if self.dataDimRefs == None:
      
      # 1D's special case!
      if self.peakList.dataSource.numDim == 1:
        
        tempDataDimRefs = [peak.sortedPeakDims()[0].dataDimRef]
        
      else:
        
        peak = self.peakList.sortedPeaks()[0]
        
        numValidDims = 0
        for dataDim in self.peakList.dataSource.dataDims:
          if dataDim.className != 'SampledDataDim':
            numValidDims += 1
        
        tempDataDimRefs = numValidDims * [None]
      
        for mapping in self.peakList.findAllApplicationData(application = self.format, keyword = 'dimMapping'):
  
          (origPeakFileDim,peakDimIndex) = returnInts(mapping.value.split(tagSep))
          
          try:
            tempDataDimRefs[origPeakFileDim] = peak.sortedPeakDims()[peakDimIndex].dataDimRef
          
          except:
            print "  Error trying to map peakDim %d to original dim %s" % (peakDimIndex,origPeakFileDim)
      
      #
      # Only use info if all went well!
      #
      
      if tempDataDimRefs.count(None) == 0:
      
        self.dataDimRefs = tempDataDimRefs
        

  def peakDimSelect(self):
  
    """
    
    User interaction to get mapping between dimensions in data model peaklist and the
    dimensions in the output format peaklist
    
    """
    
    tempDataDimRefs = None
  
    while (not tempDataDimRefs):
    
      popup = self.multiDialog.PeakDimSelect(self.guiParent, self.peakList, self.format, dataDimRefs = self.dataDimRefs, order = self.peakDimOrder)
    
      if popup.dataDimRefs:
        tempDataDimRefs = popup.dataDimRefs
      else:
        self.messageReporter.showError("Error","Have to select dataDimRefs - retry",self.guiParent)
  
    self.dataDimRefs = tempDataDimRefs
  
  def setPeakFileInfo(self):
  
    pass

  def setPeakFilePeaks(self):
    
    self.warnings = {}
    
    for self.peak in self.peakList.sortedPeaks():
      
      self.setPeakFilePeakNumber()
      
      self.setPeakFilePeakExtras()
      
      self.setPeakFilePeakIntensity()
      
      self.resetPeakFilePeakIntensity()
      
      self.setPeakAssignmentList()
      
      #
      # DataDimRefs are ordered in the sequence for peakFile
      #
      
      for self.dataDimRef in self.dataDimRefs:
        
        self.peakDim = self.peak.findFirstPeakDim(dataDimRef = self.dataDimRef)
        
        if self.peakDim:
          self.setPeakFilePeakDimInfo()
          
        else:
          break

      #
      # Only create a peak if relevant (all right dataDimRefs found)
      #
          
      if self.peakDim:

        self.createPeakFilePeak()

    #
    # Any small problems during write are reported here
    #
    
    self.printWarnings()
  
  def setPeakContribList(self,peakContribList,peakDimDict,peakDimIndex,itemList = None):

    peakDim = self.peak.sortedPeakDims()[peakDimIndex]
    peakDimResToAtoms = peakDimDict[peakDim]

    if not itemList:
      itemList = []

    for resToAtom in peakDimResToAtoms:

      if peakDimIndex == len(self.peak.peakDims) - 1:
        peakContribList.append(itemList + [resToAtom])

      else:
        self.setPeakContribList(peakContribList,peakDimDict,peakDimIndex + 1,itemList = itemList + [resToAtom])

  def setPeakAssignmentList(self):
  
    self.peakAssignmentList = []
      
    #
    # Here do reorganisation...
    #

    if self.writeAssignments:
      
      for peakContrib in self.peak.sortedPeakContribs():
        peakDimDict = {}
        peakContribList = []
          
        for peakDim in self.peak.sortedPeakDims():
          
          peakDimDict[peakDim] = []
          peakDimContribs = peakContrib.findAllPeakDimContribs(peakDim = peakDim)
            
          if peakDimContribs:
            for peakDimContrib in peakDimContribs:
              if self.resonanceToAtoms.has_key(peakDimContrib.resonance):
                resToAtoms = self.resonanceToAtoms[peakDimContrib.resonance]
                for resToAtom in resToAtoms:
                  peakDimDict[peakDim].append(resToAtom)
              else:
                peakDimDict[peakDim].append(None)
          else:
             peakDimDict[peakDim].append(None)
          
        # Now need to make all combinations in a UNIQUE way!
        self.setPeakContribList(peakContribList,peakDimDict,0)
          
        self.peakAssignmentList.extend(peakContribList)
      
    if self.peakAssignmentList:
      
      #
      # TODO: does this handle situations where you have, in each of two different dimensions, a set of prochirals?
      #
      
      numDims = len(self.peakAssignmentList[0])
      for i in range(0,numDims):
        j = 0
        while (j < len(self.peakAssignmentList)):
          contrib = self.peakAssignmentList[j]
          resToAtom = contrib[i]
          
          k = j+1

          while (k < len(self.peakAssignmentList)):
            otherContrib = self.peakAssignmentList[k]
            otherResToAtom = otherContrib[i]
            
            if contrib[:i] == otherContrib[:i] and contrib[i+1:] == otherContrib[i+1:]:

              #
              # ResonanceToAtoms objects are included at end - ignore that bit
              #
              
              if resToAtom == otherResToAtom:
                # Exactly the same - can be removed
                self.peakAssignmentList.pop(k)
                
              elif resToAtom and otherResToAtom:
 
                if resToAtom.getAtomId()[:-1] == otherResToAtom.getAtomId()[:-1]:
                  # Exactly the same - can be removed
                  self.peakAssignmentList.pop(k)

                elif ((self.compressResonances)
                      and resToAtom.getAtomSetId() and otherResToAtom.getAtomSetId() 
                      and resToAtom.getAtomSetId()[:-1] == otherResToAtom.getAtomSetId()[:-1]):
                  resToAtom.useAtomSetName = 1
                  self.peakAssignmentList.pop(k)
                  
                else:
                  k += 1

              else:
                k += 1

            else:
              k += 1

          j += 1
  
  def setPeakFilePeakNumber(self):
    
    self.peakNum = None
    
    #
    # Use original peak number if required.
    #
    
    if self.useOriginalNumbers:
      appData = self.peak.findFirstApplicationData(application = self.format, keyword = peakNum_kw)
      
      if not appData is None:
        self.peakNum = appData.value
      
    if self.peakNum is None:
    
      self.peakNum = self.peak.serial
      
    # Keep track of exported number in case part of IO cycle
    if self.tagIoCycle:
      setTag = str((self.tagIoCycle,self.peakNum))
      self.setIoCycleTag(self.peak,setTag,isUnique = True)

  def setPeakFilePeakExtras(self):
  
    pass

  def setPeakFilePeakIntensity(self):
  
    pass
  
  def resetPeakFilePeakIntensity(self):
    
    if self.forcePositiveVolumes:
      for attrName in ('intensity','volume'):
        if hasattr(self,attrName):
          setattr(self,attrName,abs(getattr(self,attrName)))
    
  def setPeakFilePeakDimInfo(self):
  
    pass
  
  def createPeakFilePeak(self):
   
    pass

  def createAtomMeasurements(self):
    
    self.resonancesToLink = []
    self.atomMeasurements = {}
  
    if self.measurementResNum != 1:
    
      self.messageReporter.showError("Not supported","Output of measurement list with %d resonances not supported yet!" % self.measurementResNum)
      return
      
      
    for measurement in self.measurementList.measurements:
    
      self.resonancesToLink.append(measurement.resonance)
    
    self.forceSpecificNamingSystemName()
    
    if self.exportMode == 'atoms':
      self.getResonanceAtomLinks()
      resonanceLinks = self.resonanceToAtoms
      
    elif self.exportMode == 'spinSystems':
      self.getResonanceSpinSystemLinks()
      resonanceLinks = self.resonanceToSpinSystems      

    #
    # Now create self.atomMeasurements
    #
    
    for resonance in resonanceLinks.keys():
     
      measurement = getattr(resonance,"findFirst%s" % self.measurementApiCode)(parentList = self.measurementList)

      if measurement:
        for resonanceToAtom in resonanceLinks[resonance]:
          self.atomMeasurements[resonanceToAtom] = measurement
       
    #
    # And keep track of original values...
    #

    self.origAtomMeasurements = self.atomMeasurements.copy()
  
  def forceSpecificNamingSystemName(self):
  
    pass
  
  def getChainListFromResonanceToAtoms(self):
  
    """
    
    Sets up a list of chains occuring in the ResonanceToAtoms objects
    
    """
    
    self.chains = []
    
    for resonanceToAtoms in self.resonanceToAtoms.values():
    
      for resonanceToAtom in resonanceToAtoms:
      
        if resonanceToAtom.chain not in self.chains:
        
          self.chains.append(resonanceToAtom.chain)          
    
  def getChainMapping(self,resetMapping = 0):
  
    """
    
    Code to either get the chain mapping if already there or to set it
    This is used in all export functions to make sure the residue information
    is always the same in all exports.
    
    self.chains has to be defined for this function
    
    """

    self.chainDict = {}
    
    if not self.chains:
    
      #
      # Means that no resonances linked to atoms used for output - is fine
      # Only supposed to happen for writePeaks and writeSequence...
      #
    
      return True
    
    #
    # Check if all chains within one molSystem
    #

    self.molSystem = self.chains[0].molSystem
    
    for chain in self.chains[1:]:
      if self.molSystem != chain.molSystem:
        print "  Error: invalid chain list (has to be same molSystem)"
        return False
    
    chainList = self.chains[:]
    
    #
    # Special case - use CCPN information for export (TODO: could result in problems for non-chaincode formats, but hey...)
    #
    
    if self.useCcpnChainInfo:
      for chain in chainList:
        exportChainCode = self.getExportChainCode(chain.code,chain)
        self.chainDict[chain] = [exportChainCode,chain.sortedResidues()[0].seqCode]

      return None
    
    #
    # Find out if mapping available - if resetMapping is off or
    # minimalPrompts is on
    #
    
    if not resetMapping: # used to include  or self.minimalPrompts, but removed (Wim 11/04/05)
     
      for chain in self.chains:
    
        chainMappings = self.nameMapping.findAllChainMappings(chainCode = chain.code)
        chainMapping = None

        for tempChainMapping in chainMappings:
          
          if not chainMapping:
            chainMapping = tempChainMapping
          else:
            if chainMapping.firstSeqId > tempChainMapping.firstSeqId:
              chainMapping = tempChainMapping

        if chainMapping:
        
          #
          # Just in case firstSeqId is not one, recalculate...
          # 
       
          firstSeqCode = chainMapping.formatFirstSeqCode - chainMapping.firstSeqId + 1
      
          exportChainCode = self.getExportChainCode(chainMapping.formatChainCode,chain)
          self.chainDict[chain] = [exportChainCode,firstSeqCode]
          
          chainList.pop(chainList.index(chain))
    
    if chainList:

      #
      # Get info for export window
      # TODO: problem here that if chain mapped previously there is no
      # more consistency checking for seqCodes (if single export chain)
      # Have to build this in as well!
      #
    
      #
      # Some formats (e.g. nmrStar) have fixed seqCodes...
      # others might set mapping automatically if appropriate
      #
      # Exception is when resetMapping is True - can set whatever
      #
      
      mappingChainDict = self.getPresetChainMapping(chainList)
      
      if not mappingChainDict or (resetMapping and not self.minimalPrompts):
      
        #
        # Do manual selection...
        #
      
        requireChainCode = None
    
        if exportChainCodeDict.has_key(self.format):
          requireChainCode = exportChainCodeDict[self.format]
    
        interaction = self.multiDialog.chainExport(self.guiParent,chainList,requireChainCode = requireChainCode)
        
        mappingChainDict = interaction.chainDict
    
      #
      # Set the information...
      #
    
      for chain in mappingChainDict.keys():
        
        (exportChainCode,exportFirstSeqCode) = mappingChainDict[chain]
      
        #
        # TODO: is this correct or should formatChainCode be allowed to be
        # None?
        #
        
        if not exportChainCode:
          exportChainCode = defaultMolCode
          
        chainMappings = self.nameMapping.findAllChainMappings(chainCode = chain.code)

        if not chainMappings:
        
          chainMapping = self.nameMapping.newChainMapping(
                                  
                                  chainCode = chain.code,
                                  firstSeqId = 1,
                                  formatChainCode = exportChainCode,
                                  formatFirstSeqCode = exportFirstSeqCode,
                                  formatSeqInsertCode = defaultSeqInsertCode
              
                                  )

        else:
        
          #
          # Reset chainmapping(s), pick consistent one to reset
          #
          
          lowestSerial = 9999
          for tmpChainMapping in chainMappings:
            if tmpChainMapping.serial < lowestSerial:
              lowestSerial = tmpChainMapping.serial
              chainMapping = tmpChainMapping 
          
          print "  Warning: resetting chain mapping..."
        
          chainMapping.formatChainCode = exportChainCode
          chainMapping.firstSeqId = 1
          chainMapping.formatFirstSeqCode = exportFirstSeqCode
          chainMapping.formatSeqInsertCode = defaultSeqInsertCode
          
          for tmpChainMapping in chainMappings:
            if not tmpChainMapping == chainMapping:
              tmpChainMapping.delete()
   
        exportChainCode = self.getExportChainCode(exportChainCode,chain)
        self.chainDict[chain] = [exportChainCode,exportFirstSeqCode]
  
    return None
  
  def getExportChainCode(self,exportChainCode,chain):
    
    # Returns a one letter code as default (some formats use segId, like CNS, Auremol and CHARMM - this function is redefined there)
    
    if len(exportChainCode) > 1:
      exportChainCode = chain.pdbOneLetterCode
    
    return exportChainCode
  
  def getPresetChainMapping(self,chainList):
  
    return None
  
  def getSortedChainList(self):
  
    #
    # Quick sort routine - according to chain code and seqCode within same chain code
    #
    
    # simplified code. Rasmus 22/7/2011
    self.chainList = [tt[1] for tt in sorted((item[1],item[0]) 
                            for item in self.chainDict.items())]
    
    
  def selectChains(self):
  
    """
    
    Select relevant chains from ones that were found in measurement list.
    Is not working yet. Also note that chains can be passed in!
    
    """  
    
    #
    # TODO: Is this really necessary??!?
    #
    
    #
    # Check if anything left or set to all chains in chemical shift list
    #  
    
    if self.atomMeasurements:
    
      if not self.chains and self.exportMode == 'atoms':

        #
        # Use all chains - note that output order depends on keys list... can change
        #

        resonanceToAtoms = self.atomMeasurements.keys()

        self.chains = [resonanceToAtoms[0].chain]

        for resonanceToAtom in resonanceToAtoms[1:]:

          if resonanceToAtom.chain not in self.chains:

            inserted = 0

            for i in range(0,len(self.chains)):

              if self.chains[i].code > resonanceToAtom.chain.code:

                self.chains.insert(i,resonanceToAtom.chain)
                inserted = 1
                break

            if not inserted:

              self.chains.append(resonanceToAtom.chain)

    else:
   
      #
      # If no measurements, set self.chains to [] so rest of function is ignored...
      #
       
      self.chains = []

  def selectStructureGeneration(self):
    
    # TODO should this be generic code somewhere? self.setNmrProject?
    if not self.project.currentNmrProject:
      nmrProject = self.project.findFirstNmrProject()
      if not nmrProject:
        nmrProject = self.project.newNmrProject(name = self.project.name)
      self.project.currentNmrProject = nmrProject
    
    if self.nmrConstraintStore:
      strucGenList = [self.project.currentNmrProject.findFirstStructureGeneration(nmrConstraintStore = self.nmrConstraintStore)]
    else:
      strucGenList = self.project.currentNmrProject.sortedStructureGenerations()
      
    strucGen = None
    
    if self.minimalPrompts and len(strucGenList) == 1:
    
      strucGen = self.project.currentNmrProject.findFirstStructureGeneration()

    elif strucGenList:
    
      (selectionList,selectionDict) = createSelection(strucGenList)

      interaction = self.multiDialog.SelectionList(

                         self.guiParent,
                         selectionList,
                         selectionDict = selectionDict,
                         title = "Project '%s': " % self.project.name + 'Select structure generation',
                         text = "Structure generations:",
                         dismissText = 'Create new',
                         urlFile = 'SelectStrucGen'
                        )

      if interaction.isSelected:
        strucGen = interaction.selection
          
    if not strucGen:
      strucGenName = 'Structure generation #%d' % (len(self.project.currentNmrProject.structureGenerations) + 1)
      if not self.minimalPrompts:
        strucGenName = self.dataEntry.askString("Enter structure generation name","Name of structure generation",strucGenName,self.guiParent)

      strucGen = Nmr.StructureGeneration(self.project.currentNmrProject, name = strucGenName)
      
    self.strucGen = strucGen

  def selectStructureEnsemble(self):
    
    structureEnsembles = self.molSystem.sortedStructureEnsembles()
    
    # Only do selection if interaction desired - otherwise create new one from scratch
    if not self.minimalPrompts:
    
      self.structureEnsemble = None
      
      (selectionList,selectionDict) = createSelection(structureEnsembles)
    
      if selectionList:
        while (not self.structureEnsemble):

          interaction = self.multiDialog.SelectionList(

                             self.guiParent,
                             selectionList,
                             selectionDict = selectionDict,
                             title = "Project '%s': " % self.project.name + 'Select structure ensemble',
                             text = 'Existing structure ensembles:',
                             urlFile = 'SelectStructureEnsemble')
    
          if interaction.isSelected:
            self.structureEnsemble = interaction.selection
    
    if not self.structureEnsemble:
      
      # NBNB 20/2/09 behaviour changed by Rasmus for global consistency
      #newEnsembleId = 0
      #
      #for tempStrucEnsemble in self.project.structureEnsembles:
      #  if tempStrucEnsemble.ensembleId > newEnsembleId:
      #    newEnsembleId = tempStrucEnsemble.ensembleId
      #
      if structureEnsembles:
        newEnsembleId = max(x.ensembleId for x in structureEnsembles)
      else:
        newEnsembleId = 0
      
      newEnsembleId += 1
      self.structureEnsemble = self.project.newStructureEnsemble(molSystem = self.molSystem, ensembleId = newEnsembleId, atomNamingSystem = self.namingSystemName)
      
      self.convertCount[self.mainCode][0] += 1
      
    if self.strucGen and not self.structureEnsemble.structureGeneration:
      self.structureEnsemble.structureGeneration = self.strucGen

  def selectEntry(self):
  
    entryList = None

    setCurrentStore(self.project,'NmrEntryStore')
    entryList = self.project.currentNmrEntryStore.sortedEntries()

    entry = None

    if entryList:

      (selectionList,selectionDict) = createSelection(entryList)

      interaction = self.multiDialog.SelectionList(

                         self.guiParent,
                         selectionList,
                         selectionDict = selectionDict,
                         title = "Project '%s': " % self.project.name + 'Select entry for project export',
                         text = "Entries:",
                         dismissText = 'Create new',
                         urlFile = 'SelectEntry'
                         
                        )

      if interaction.isSelected:
        entry = interaction.selection
          
    if not entry:
      
      setCurrentStore(self.project,'NmrEntryStore')
    
      entry = NmrEntry.Entry(self.project.currentNmrEntryStore,name = 'allInProject')

      # TODO TODO have to link to selected stuff in project!?!?
      
      entry.experiments = self.project.currentNmrProject.sortedExperiments()
      entry.measurementLists = self.project.currentNmrProject.measurementLists
      
      if self.project.molSystems:
        self.molSystem = None
        self.selectMolSystem()
        entry.molSystem = self.molSystem
      
      entry.structureGenerations = self.project.currentNmrProject.structureGenerations      
      
    self.entry = entry

  def setMeasurementFileValues(self):
  
    """
    
    Loop over the selected chains, get all relevant information and
    set the measurement object and its attributes in the format file object
    
    For spinSystems, get the list of relevant ones and get the atoms for each.
    This can involve reorganisation according to sequential order, ... .
    
    """
    
    self.atomSerial = 0
  
    if self.exportMode == 'atoms':

      #
      # Loop over selected chains
      #
      
      for self.chain in self.chains:
        
        self.chainCode = self.chainDict[self.chain][0]
        
        #
        # Assuming residues to be in order
        #
        
        for self.residue in self.chain.sortedResidues():
          
          self.chemCompVar = self.residue.chemCompVar
  
          #
          # Set some values
          #
          
          self.initMeasurementFileResidue()
          
          #
          # Get a list of resonanceToAtom objects (in correct order)
          #
          
          self.getAtomNamesList()
          
          #
          # Loop over all chemAtoms for this residue and set values for shifts
          #
  
          for self.atomName in self.atomNamesList:
  
            self.atomSerial += 1
            self.setMeasurementFileValue()

    elif self.exportMode == 'spinSystems':
      
      self.getMeasurementSpinSystemList()
      
      self.setSpecificSpinSystemInfo()    

      for self.spinSystem in self.spinSystemList:

        self.getSpinSystemResonanceList()
        
        for self.atomNames in self.atomNamesList:
          
          self.atomSerial += 1
          self.setMeasurementFileValue()
      
  def getAtomNamesList(self):
  
    """
    
    Set up a list of resonanceToAtom objects for this residue
    
    """
    
    #
    # First make a list of (possibly) relevant objects
    #
    
    residueResonanceToAtomList = []
    
    for resonanceToAtom in self.atomMeasurements.keys():
    
      if resonanceToAtom.chain == self.chain and resonanceToAtom.seqId == self.seqId:

        insertIndex = 0
        
        for insertIndex in range(0,len(residueResonanceToAtomList)):
          if residueResonanceToAtomList[insertIndex].atomName > resonanceToAtom.atomName:
            break

        residueResonanceToAtomList.insert(insertIndex,resonanceToAtom)

    #
    # Order the resonanceToAtom objects depending on case (eg XEasy)
    #

    self.atomNamesDict = {}
    self.atomNamesList = []
    
    self.orderByAtomNames(residueResonanceToAtomList)

  def initChemShiftFileResidue(self):
    
    """
    
    Use seqCode of residue to write out to chemical shift file
    Use seqId for identification
    
    """
    
    self.seqId = self.residue.seqId
    self.seqCode = self.getExportSeqCode(self.chainDict[self.chain][1],self.residue)

  def orderByAtomNames(self,residueResonanceToAtomList):
  
    """
    
    Order the list of ResonanceToAtom object (based on atomNamesList)
    Default is order in which atoms occur in reference data.
    
    """
    
    #
    # Only use compressResonances if useOriginalResNames NOT on!!
    #
    
    if self.presetResonanceMapping:
      self.compressPresetResidueResonanceToAtomList(residueResonanceToAtomList)  

    elif (self.compressResonances and not self.useOriginalResNames 
        and not self.individualAtoms):
      
      self.compressResidueResonanceToAtomList(residueResonanceToAtomList)    
    
    #
    # When using original names just order according to atomNames...
    #
    
    if self.useOriginalResNames:
    
      #
      # Exporting resonances in order that were created...
      #
      
      resonanceSerials = {}
    
      for resonanceToAtom in residueResonanceToAtomList:
      
        atomName = resonanceToAtom.atomName
        
        self.atomNamesDict[atomName] = resonanceToAtom
        resonanceSerials[resonanceToAtom.resonance.serial] = atomName
        
      resonanceSerialList= resonanceSerials.keys()
      resonanceSerialList.sort()
      
      for resonanceSerial in resonanceSerialList:
        self.atomNamesList.append(resonanceSerials[resonanceSerial])
    
    else:

      #
      # Order the atom names, purely by alphabetical sorting. This used to be done by the order
      # in which the chemAtoms were returned, but this will soon be irrelevant.
      #
      
      chemComp = self.residue.molResidue.chemComp
      chemCompVar = self.residue.molResidue.chemCompVar
      
      namingSystem = chemComp.findFirstNamingSystem(name = self.namingSystemName)
      
      if not namingSystem:
        print "  Warning: no %s naming system for chemComp %s,%s!" % (self.namingSystemName,chemComp.molType,chemComp.ccpCode)
      
      chemAtomSysNames = []
      chemAtomOrSets = []
      
      for chemAtomOrSetList in (chemCompVar.sortedChemAtoms(), chemCompVar.sortedChemAtomSets()):
        
        for chemAtomOrSet in chemAtomOrSetList:
          chemAtomSysName = findChemAtomOrSetSysName(namingSystem,chemAtomOrSet)
          
          if chemAtomSysName:
            chemAtomSysNames.append(chemAtomSysName)
          
          elif not isinstance(chemAtomOrSet,ChemComp.LinkAtom):
            chemAtomOrSets.append(chemAtomOrSet)
                
      #
      # Look for atomName matches... first try to find full matches...
      #
      # Note that is it necessary to go through both the ChemAtomSysNames and the ChemAtomOrSet
      # lists in case no sysnames are available for a particular atom!!
      #

      recheckChemAtomSysNames = []

      for chemAtomSysName in chemAtomSysNames:
      
        atomName = chemAtomSysName.sysName
        self.atomNamesList.append(atomName)
        
        atomFound = self.getAtomNameMatch(atomName,residueResonanceToAtomList)
    
        if not atomFound:
        
          recheckChemAtomSysNames.append(chemAtomSysName)

      for chemAtomOrSet in chemAtomOrSets:
      
        atomName = chemAtomOrSet.name
        self.atomNamesList.append(atomName)
        
        atomFound = self.getAtomNameMatch(atomName,residueResonanceToAtomList)
    
        # If atom not found, doesn't mean anything... just that it has no measurement
          
      #
      # Could be resonanceToAtom present for IUPAC name (e.g. CD* for CD1, CD2)
      #
        
      for chemAtomSysName in recheckChemAtomSysNames:
          
        atomName = chemAtomSysName.sysName
        atomFound = False

        for resonanceToAtom in residueResonanceToAtomList:

          chemAtomOrSet = chemCompVar.findFirstChemAtom(name = chemAtomSysName.atomName)
          
          if isinstance(chemAtomOrSet,ChemComp.ChemAtom):
          
            chemAtomSet = chemAtomOrSet.chemAtomSet

            if chemAtomSet:
              
              for resonanceToAtom in residueResonanceToAtomList:

                if chemAtomSet.name == resonanceToAtom.atomName:

                  self.atomNamesDict[chemAtomSysName.sysName] = resonanceToAtom
                  atomFound = True
                  break
                             
  def getMeasurementSpinSystemList(self):
    
    spinSystemSerials = {}
    
    self.spinSystemList = []
    self.spinSystemToResonanceDict = {}
    
    for resonanceToSpinSystem in self.atomMeasurements.keys():
      spinSystem = resonanceToSpinSystem.spinSystem
      spinSystemSerials[spinSystem.serial] = spinSystem
      
      if not self.spinSystemToResonanceDict.has_key(spinSystem):
        self.spinSystemToResonanceDict[spinSystem] = []
      self.spinSystemToResonanceDict[spinSystem].append(resonanceToSpinSystem)
      
    spinSystemSerialList = spinSystemSerials.keys()
    spinSystemSerialList.sort()
    
    for sss in spinSystemSerialList:
      self.spinSystemList.append(spinSystemSerials[sss])
  
  def setSpecificSpinSystemInfo(self):
    
    pass
  
  def getSpinSystemResonanceList(self):
    
    self.atomNamesList = []
    self.atomNamesDict = {}
    
    for resonanceToSpinSystem in self.spinSystemToResonanceDict[self.spinSystem]:
      atomNames = tuple(resonanceToSpinSystem.atomNames)
      
      self.atomNamesList.append(atomNames)
      self.atomNamesDict[atomNames] = resonanceToSpinSystem

  def getAtomNameMatch(self,atomName,residueResonanceToAtomList):

    #
    # Look for atomName match...
    #

    atomFound = 0

    for resonanceToAtom in residueResonanceToAtomList:

      if atomName == resonanceToAtom.atomName:
        
        if self.atomNamesDict.has_key(atomName):
          print "  ERROR: Atom name confusion for '%s.%d.%s'!" % (self.chain.code,self.residue.seqId,atomName)
        self.atomNamesDict[atomName] = resonanceToAtom
        residueResonanceToAtomList.pop(residueResonanceToAtomList.index(resonanceToAtom))
        atomFound = 1

        break
    
    return atomFound

  def compressResidueResonanceToAtomList(self,residueResonanceToAtomList):

    self.chemComp = self.residue.molResidue.chemComp
    namingSystem = self.chemComp.findFirstNamingSystem(name = self.namingSystemName)
    
    #
    # Try to 'compress' the resonance names if so desired...
    #

    for resonanceToAtom in residueResonanceToAtomList[:]:

      chemAtomSet = resonanceToAtom.getChemAtomSet()
      chemCompVar = resonanceToAtom.getResidue().chemCompVar

      if chemAtomSet:

        for otherResonanceToAtom in residueResonanceToAtomList:

          if otherResonanceToAtom == resonanceToAtom:
            continue

          otherChemAtomSet = otherResonanceToAtom.getChemAtomSet()

          if chemAtomSet == otherChemAtomSet and self.atomMeasurements.has_key(resonanceToAtom) and self.atomMeasurements.has_key(otherResonanceToAtom):

            measurement = self.atomMeasurements[resonanceToAtom]
            otherMeasurement = self.atomMeasurements[otherResonanceToAtom]

            if (abs(measurement.value - otherMeasurement.value) < shiftIdentityTolerance):
              
              chemAtomSetSysName = self.getChemAtomOrSetSysName((chemAtomSet,namingSystem))

              if chemAtomSetSysName:
                newAtomName = chemAtomSetSysName.sysName
                
                self.setResonanceToAtoms(resonanceToAtom.resonance,
                                         resonanceToAtom.chain,
                                         resonanceToAtom.seqId,
                                         newAtomName)

                newResonanceToAtom = self.resonanceToAtoms[resonanceToAtom.resonance][-1]
                useIndex = residueResonanceToAtomList.index(resonanceToAtom)
                residueResonanceToAtomList.insert(useIndex,newResonanceToAtom)
                self.atomMeasurements[newResonanceToAtom] = measurement
                self.origAtomMeasurements[newResonanceToAtom] = measurement
                
                listAtomNames = [resonanceToAtom.atomName,otherResonanceToAtom.atomName]
                listAtomNames.sort()
                
                print ("  Warning: %s:%s compressed %s and %s to %s name..." 
                       % (self.chemComp.ccpCode, self.residue.seqId,
                          listAtomNames[0], listAtomNames[1], newAtomName))
                          
                residueResonanceToAtomList.pop(residueResonanceToAtomList.index(resonanceToAtom))
                residueResonanceToAtomList.pop(residueResonanceToAtomList.index(otherResonanceToAtom))
                del(self.atomMeasurements[resonanceToAtom])
                del(self.atomMeasurements[otherResonanceToAtom])
                del(self.origAtomMeasurements[resonanceToAtom])
                del(self.origAtomMeasurements[otherResonanceToAtom])

                break

               

  def compressPresetResidueResonanceToAtomList(self,residueResonanceToAtomList):
    """ compress ResidueResonanceToAtomList for preset resonances
    Here shift comparison has been done up front, and useAtomSetName has been set
    """
    
    self.chemComp = self.residue.molResidue.chemComp
    namingSystem = self.chemComp.findFirstNamingSystem(name = self.namingSystemName)
    
    #
    # Try to 'compress' the resonance names if so desired...
    #

    for resonanceToAtom in residueResonanceToAtomList[:]:
      
      #r2a = resonanceToAtom
      #print ('### PRESET', r2a.atomName, r2a.atomSetName, 
      #                r2a.chemAtom and  r2a.chemAtom.name,
      #                r2a.chemAtomSet and  r2a.chemAtomSet.name,
      #                r2a.useAtomSetName, r2a.resonance,  
      #                [x for x in r2a.otherLinkedResonances], 
      #                [x for x in r2a.otherGroupResonances])
      
      if resonanceToAtom.useAtomSetName:

        bigChemAtomSet = resonanceToAtom.getChemAtomSet()
        chemCompVar = resonanceToAtom.getResidue().chemCompVar

        if (bigChemAtomSet 
            and bigChemAtomSet is not resonanceToAtom.chemAtomSet):
              
          chemAtomSetSysName = self.getChemAtomOrSetSysName((bigChemAtomSet,
                                                             namingSystem))

          if chemAtomSetSysName:
            newAtomName = chemAtomSetSysName.sysName
            
            self.setResonanceToAtoms(resonanceToAtom.resonance,
                                     resonanceToAtom.chain,
                                     resonanceToAtom.seqId,
                                     newAtomName, 
                                     chemAtomOrSet=bigChemAtomSet)
            newResonanceToAtom = self.resonanceToAtoms[
                                        resonanceToAtom.resonance][-1]
            
            measurement = self.atomMeasurements.pop(resonanceToAtom)
            del self.origAtomMeasurements[resonanceToAtom]
            self.atomMeasurements[newResonanceToAtom] = measurement
            self.origAtomMeasurements[newResonanceToAtom] = measurement
          
            ii = 0
            while ii < len(residueResonanceToAtomList):
              otherResonanceToAtom = residueResonanceToAtomList[ii]
 
              
              if otherResonanceToAtom == resonanceToAtom:
                residueResonanceToAtomList[ii] = newResonanceToAtom
                ii += 1

              elif (bigChemAtomSet == otherResonanceToAtom.getChemAtomSet()
                  and self.atomMeasurements.has_key(resonanceToAtom)
                  and self.atomMeasurements.has_key(otherResonanceToAtom)):
 
                # We have a match - remove second copy
                del residueResonanceToAtomList[ii]
                del self.atomMeasurements[otherResonanceToAtom]
                del self.origAtomMeasurements[otherResonanceToAtom]
                
                self.setResonanceToAtoms(otherResonanceToAtom.resonance,
                                         otherResonanceToAtom.chain,
                                         otherResonanceToAtom.seqId,
                                         newAtomName, 
                                         chemAtomOrSet=bigChemAtomSet)
              
              else:
                ii += 1
            

               
  def setChemShiftFileValue(self):
  
    """
    
    Always format-specific. Creates the chemical shift object and sets its
    attributes in the format file object
    
    """
    
    pass
  
    
  def checkAtomMeasurements(self):
  
    """
    
    Check if the script was unable to handle some of the measurements/atoms
    
    """

    if self.verbose and self.atomMeasurements != {}:

      print "Unable to process following resonances:"
      
      outputLines = []
      
      for resonanceToAtom in self.atomMeasurements.keys():
      
        measurement = self.atomMeasurements[resonanceToAtom].value
      
        if self.useOriginalResNames:
          outputLines.append("  Original resonance name: '%s', %s %.3f" % (resonanceToAtom.originalResName,self.measurementType,measurement))
        
        else:
          outputLines.append("  Ccp identification: chain %s, seqId %d, atomName %s, %s %.3f" % (resonanceToAtom.chain.code,resonanceToAtom.seqId,resonanceToAtom.atomName,self.measurementType,measurement))
      
      outputLines.sort()
      
      for line in outputLines:
        print line
        
  #
  # Constraint specific
  #
  
  def setConstraints(self):
  
    for self.rawConstraint in self.constraintFile.constraints:
    
      if not self.rawConstraintValid():
        continue
      
      self.setTypeConstraint()
      
  def rawConstraintValid(self):
  
    return 1

  def getConstraintItems(self):
  
    self.items = []
    self.rawItems = []
    
    for item in self.rawConstraint.items:

      if len(item.members) != self.constraintFile.constraintElements:
        print "  %s constraint item with %d members (Id %d)... ignored" % (self.constraintApiCode,len(item.members),self.rawConstraint.Id)
        continue
          
      self.items.append([])
      self.rawItems.append(item)
      
      #
      # WARNING: for Ccp currently *TWO* resonances are linked, except for dihedrals (4) and chemshifts (1)!!!
      #          another exception is J couplings, where in CNS for example 4 members are present, but only 2 are stored in CCPN!
      #
      
      if len(item.members) > self.constraintResNum:
        itemMembers = self.getConstraintItemMembers(item)
      else:
        itemMembers = item.members      
      
      for member in itemMembers:
      
        self.items[-1].append(member)
     
        # Creating resonance(s)... these should be cleaned up later
        # Name of the resonance is format <chainCode>.<seqCode>.<atomName>
        #      or in Ccp                  <chain.code>.<residue.seqCode>.<atom.name>
        
        chainCode = member.chainCode
         
        #
        # Set the resonance name
        #

        resName = getResName(chainCode,member.seqCode,member.atomName,seqInsertCode = member.seqInsertCode)
        
        # Get the resonance (if any)
        # Quickly set a link from the member to the resonance(s)... not good style but easy!
        if self.resonanceNames.has_key(resName):
          member.resonances = self.resonanceNames[resName]
          
        else:
          # Create new FIXED resonance
          resonance = self.createResonance(resName,fixed = True,origSource = member)
          self.resonanceNames[resName] = [resonance]      
          member.resonances = [resonance]
  
  def getConstraintItemMembers(self,item):
    
    # TODO: These are pretty crappy defaults...
    if self.constraintResNum == 1:
      itemMembers = [item.members[0]]
    elif self.constraintResNum == 2:
      itemMembers = [item.members[0],item.members[-1]]
      
    return itemMembers
  
  #
  # Constraint specific (EXCEPT dihedrals and CSA!)
  # 
    
  def setDistanceConstraintItems(self):
  
    self.setConstraintItems()

  def setConstraintItems(self):

    #
    # Now have the resonances... see if there are no constraintItems that
    # correspond to the current resonances
    #
    #
    # Note that this is currently always two (except for chem shift ones but
    # not handled yet)
    #
    
    for i in range(len(self.items)):
      
      (rawItemMember1,rawItemMember2) = self.items[i]
      self.rawConstraintItem = self.rawItems[i]
      
      #
      # This is basically assuming a full OR relationship within the node
      #
      # WARNING: THESE ARE FIXED RESONANCES WITHIN THE NMRCONSTRAINTS PACKAGE!!!
      #
      
      resonances1 = rawItemMember1.resonances
      resonances2 = rawItemMember2.resonances

      for fixedResonance1 in resonances1:
        
        for fixedResonance2 in resonances2:
        
          if fixedResonance1 == fixedResonance2:
            print "  Error: duplicate resonance %s found - constraint not imported" % getResNameText(fixedResonance1)
            self.setConversionError(self.constraintList,self.rawConstraint.Id,1,'Not handling item from distance restraint (original ID %d) - duplicate resonances')
            continue

          fixedResonances = frozenset((fixedResonance1,fixedResonance2)) 
          resonanceMatch = 0
           
          for constrItem in fixedResonance1.findAllPairwiseConstraintItems(className = "%sConstraintItem" % self.constraintApiCode, constraint = self.constraint):
            # Comparing frozensets, order does not matter!              
            if fixedResonances == constrItem.resonances:
              resonanceMatch = 1
              break
           
          if resonanceMatch == 1:
            #print "  constraint item %s to %s already exists... ignored." % ((getApplResNames(constraintFormat,[resonances[0]])).keys(),
            #                                                                          (getApplResNames(constraintFormat,[resonances[1]])).keys())
            pass

          else:
             
            #print "  Doing %s to %s." % ((getApplResNames(self.format,[resonances[0]])).keys(),
            #                             (getApplResNames(self.format,[resonances[1]])).keys())
 
            #
            # Original order is lost in Api unless connected to peakDims...
            #
             
            self.constraintItem = self.ApiConstraintItem(self.constraint, resonances = fixedResonances, firstResonance = fixedResonance1)
            self.setConstraintItemInfo()

    if not self.constraint.items:

      print "  Error: constraint without items created - is now deleted."
      self.constraint.delete()
      self.setConversionError(self.constraintList,self.rawConstraint.Id,1,'Distance restraint (original ID %d) deleted - no constraint items left.')
      self.convertCount[self.mainCode][1][self.valueCode] -= 1
  
  def setConstraintItemInfo(self):
  
    #
    # Format specific info (initiated for nmrStar constraint handling)
    #
  
    pass
  
  #
  # Distance/Hbond constraint specific (import)
  #      

  def setDistanceTypeConstraint(self):
    
    self.setDistanceConstraintPars()
    self.setDistanceConstraint()
    self.setDistanceConstraintInfo()

    self.getConstraintItems()
    self.setDistanceConstraintItems()
    
  def setHBondTypeConstraint(self):
    
    #
    # Assuming this is handled exactly the same as distance constraints
    #
    
    self.setDistanceTypeConstraint()

  def setDistanceConstraintPars(self):
    
    # TODO: this is not really standard...
    
    self.constraintPars = {}
    
    self.setDistanceConstraintLimits()
    
    if hasattr(self.rawConstraint,'weight') and self.rawConstraint.weight != None:
      self.constraintPars['weight'] = self.rawConstraint.weight
    else:
      self.constraintPars['weight'] = 1.0
    
    if hasattr(self.rawConstraint,peakNum_kw) and self.rawConstraint.peakNum != None:
      self.constraintPars['origId'] = self.rawConstraint.peakNum
    else:
      self.constraintPars['origId'] = self.rawConstraint.Id
    
    for origDataKeyword in ('volume','intensity'):
      if hasattr(self.rawConstraint,origDataKeyword):
        value = getattr(self.rawConstraint,origDataKeyword)
        if value != None:
          self.constraintPars['origData'] = value
          break

    #
    # TODO: WHAT ABOUT LINK TO PEAKS? PPM VALUES?!?!
    #       Need a link to the correct peak IF available... or should
    #       this link be created?!
    #

  def setDistanceConstraint(self):
    
    #
    # Check whether original Id exists and set as detail... remove from dict because not API attribute
    #
    
    if self.constraintPars.has_key('origId'):
      if self.constraintPars['origId']:
        self.constraintPars['details'] = "Original Id %d" % self.constraintPars['origId']
      del(self.constraintPars['origId'])
    
    self.constraint = self.ApiConstraint(self.constraintList, **self.constraintPars)
    self.convertCount[self.mainCode][1][self.valueCode] += 1
    
  def setDistanceConstraintInfo(self):
  
    pass

  def setDistanceConstraintLimits(self):
    
    self.constraintPars['upperLimit'] = self.rawConstraint.upperDist
    
    if hasattr(self.rawConstraint,'lowerDist') and self.rawConstraint.lowerDist:
      self.constraintPars['lowerLimit'] = self.rawConstraint.lowerDist
    else:
      self.constraintPars['lowerLimit'] = 1.8
      
    self.constraintPars['targetValue'] = (self.constraintPars['upperLimit'] + self.constraintPars['lowerLimit']) / 2  
  
  #
  # Rdc constraint specific (import)
  #      

  def setRdcTypeConstraint(self):
    
    self.setRdcConstraintPars()
    self.setRdcConstraint() 
    self.setRdcConstraintInfo()

    self.getConstraintItems()
    self.setRdcConstraintItems()
    
  def setRdcConstraintInfo(self):
  
    pass
  
  def setRdcConstraintPars(self):
  
    # TODO: this is not really standard...
    
    self.constraintPars = {}
    
    if hasattr(self.rawConstraint,'value'):
      self.constraintPars['targetValue'] = self.rawConstraint.value
    
    if hasattr(self.rawConstraint,'error'):
      self.constraintPars['error'] = self.rawConstraint.error
    
    if hasattr(self.rawConstraint,'lowerValue'):
      self.constraintPars['lowerLimit'] = self.rawConstraint.lowerValue
    
    if hasattr(self.rawConstraint,'upperValue'):
      self.constraintPars['upperLimit'] = self.rawConstraint.upperValue
  
  def setRdcConstraint(self):
  
    self.setDistanceConstraint()
  
  def setRdcConstraintItems(self):
  
    self.setConstraintItems()

  #
  # Csa constraints (basically same as chemical shift bit lower)
  #      

  def setCsaTypeConstraint(self):
        
    self.setCsaConstraintPars()
    
    # Only ever one item, get it here
    self.getConstraintItems()

    self.setConstraintListElementType()

    self.setCsaConstraint() 
    self.setCsaConstraintInfo()
  
  def setCsaConstraintPars(self):
    
    self.setChemShiftConstraintPars()
  
  def setCsaConstraint(self):
    
    self.setChemShiftConstraint()
  
  def setCsaConstraintInfo(self):
    
    pass
  
  #
  # Chemical shift constraints
  #

  def setChemShiftTypeConstraint(self):
        
    self.setChemShiftConstraintPars()
    
    # Only ever one item, get it here
    self.getConstraintItems()

    self.setConstraintListElementType()

    self.setChemShiftConstraint() 
    self.setChemShiftConstraintInfo()
    
  def setConstraintListElementType(self):
    
    if not self.constraintList.isotopeCode:
      
      rawItemMember = self.items[0][0]
      resonances = rawItemMember.resonances
      
      self.constraintList.isotopeCode = resonances[0].isotopeCode
    
  def setChemShiftConstraintInfo(self):
  
    pass
  
  def setChemShiftConstraintPars(self):
  
    # TODO: this is not really standard...
    
    self.constraintPars = {}
    
    if hasattr(self.rawConstraint,'value'):
      self.constraintPars['targetValue'] = self.rawConstraint.value
    
    if hasattr(self.rawConstraint,'error'):
      self.constraintPars['error'] = self.rawConstraint.error
    
    if hasattr(self.rawConstraint,'lowerValue'):
      self.constraintPars['lowerLimit'] = self.rawConstraint.lowerValue
    
    if hasattr(self.rawConstraint,'upperValue'):
      self.constraintPars['upperLimit'] = self.rawConstraint.upperValue
  
  def setChemShiftConstraint(self):
  
    self.constraint = None
    
    #
    # First get the relevant resonance and check whether a constraint with this resonance exists already
    #
  
    if len(self.items) > 1:
      print "  Warning: %s constraint with multiple items - invalid, only using first one." % self.constraintType
        
    rawItemMember = self.items[0][0]
    self.rawConstraintItem = self.rawItems[0]
          
    resonances = rawItemMember.resonances

    for fixedResonance in resonances:
      
      resonanceExists = self.constraintList.findFirstConstraint(resonance = fixedResonance)
       
      if not resonanceExists:

        #
        # Create the constraint - check whether original Id exists and set as detail... remove from dict because not API attribute
        #
        
        if self.constraintPars.has_key('origId'):
          if self.constraintPars['origId']:
            self.constraintPars['details'] = "Original Id %d" % self.constraintPars['origId']
          del(self.constraintPars['origId'])
          
        self.constraintPars['resonance'] = fixedResonance
        
        self.constraint = self.ApiConstraint(self.constraintList, **self.constraintPars)
        self.convertCount[self.mainCode][1][self.valueCode] += 1
    
    if not self.constraint:

      print "  Error: constraint could not be created."
      self.setConversionError(self.constraintList,self.rawConstraint.Id,1,'Restraint (original ID %d) not created - already exists or not enough information.')  

  #
  # J coupling constraint specific (import)
  #      

  def setJCouplingTypeConstraint(self):
    
    self.setJCouplingConstraintPars()
    self.setJCouplingConstraint() 
    self.setJCouplingConstraintInfo()

    self.setupJCouplingConstraintItems()
    self.getConstraintItems()
    self.setJCouplingConstraintItems()
    
  def setJCouplingConstraintInfo(self):
  
    pass
  
  def setJCouplingConstraintPars(self):
    
    self.constraintPars = {}
    
    if hasattr(self.rawConstraint,'value'):
      self.constraintPars['targetValue'] = self.rawConstraint.value
    
    if hasattr(self.rawConstraint,'error'):
      self.constraintPars['error'] = self.rawConstraint.error
    
    if hasattr(self.rawConstraint,'lowerValue'):
      self.constraintPars['lowerLimit'] = self.rawConstraint.lowerValue
    
    if hasattr(self.rawConstraint,'upperValue'):
      self.constraintPars['upperLimit'] = self.rawConstraint.upperValue
  
  def setJCouplingConstraint(self):
  
    self.setDistanceConstraint()
    
  def setupJCouplingConstraintItems(self):
   
    # Assuming only 2 resonances passed in - in CNS for example this is not the case (4 are present)
    pass

  
  def setJCouplingConstraintItems(self):
  
    self.setConstraintItems()

  #
  # Dihedral constraint specific
  #
    
  def setDihedralTypeConstraint(self):
    
    self.getConstraintItems()
    self.setDihedralConstraint()
    self.setDihedralConstraintItem()

  def setDihedralConstraint(self):
  
    self.setDihedralConstraintPars()

    if hasattr(self.rawConstraint,'origId') and self.rawConstraint.origId:
      self.constraintPars['details'] = "Original Id %d" % self.rawConstraint.origId
      
    #
    # ONLY ONE resonance per item member and ONLY ONE item is allowed for this...
    #
    
    if len(self.items) > 1:
      print "  Warning multiple items for dihedral restraint - only first one used"
    
    rawItemMembers = self.items[0]
    curFixedResonances = []
    
    for i in range(0,self.constraintResNum):
    
      resonances = rawItemMembers[i].resonances
      
      if len(resonances) > 1:
        print "  Warning multiple resonances for dihedral restraint item - only first one used"
    
      fixedResonance = resonances[0]
      
      if fixedResonance not in curFixedResonances:
        curFixedResonances.append(fixedResonance)
      else:
        print "  Error: duplicate resonance %s found - constraint not imported" % getResNameText(fixedResonance)
        self.setConversionError(self.constraintList,self.rawConstraint.Id,1,'Not handling dihedral restraint (original ID %d) - duplicate resonances')
        self.constraint = None
    
    #
    # Only continue if no duplicate resonances found!
    #
    
    if len(curFixedResonances) == self.constraintResNum:
          
      fixedResonancesSorted = curFixedResonances[:]
      fixedResonancesSorted.sort()

      resonanceMatch = False

      for dihedConstr in curFixedResonances[0].findAllDihedralConstraints(parentList = self.constraintList):
        fixedResonancesMatch = list(dihedConstr.resonances)
        fixedResonancesMatch.sort()

        if fixedResonancesSorted == fixedResonancesMatch:
          resonanceMatch = True
          break

      if resonanceMatch and not self.newConstraintForDuplicate:

        #
        # In this case means the dihedral constraint already exists - set it
        #
        print "  Warning: Additional value found for dihedral constraint %s - added as extra item" % (self.rawConstraint.Id)
        self.constraint = self.constraintList.findFirstConstraint(resonances = dihedConstr.resonances)

      else:

        self.constraint = self.ApiConstraint(self.constraintList, resonances = curFixedResonances, **self.constraintPars)
        self.convertCount[self.mainCode][1][self.valueCode] += 1

        self.setDihedralConstraintInfo()
      
  def setDihedralConstraintInfo(self):
    
    #
    # To add application data to dihedral constraint
    #
    
    pass
  
  def setDihedralConstraintPars(self):
    
    # TODO: this is not really standard...
    # Also bit tricky with angles... where +, where -?
    
    self.constraintPars = {}
    
    # TODO: THIS IS NOT CORRECT!! Need to add more pars for Cns-type?!?
    if hasattr(self.rawConstraint,'energyCst') and self.rawConstraint.energyCst:
      self.constraintPars['weight'] = self.rawConstraint.energyCst
    else:
      self.constraintPars['weight'] = 1.0

  def setDihedralConstraintItem(self):
  
    if self.constraint:

      self.setDihedralConstraintItemPars()
      
      if None in (self.constraintPars['upperLimit'],self.constraintPars['lowerLimit']):
        print "  Error: trying to create dihedral constraint item without bounds."
        return
        
      self.constraintItem = self.ApiConstraintItem(self.constraint,**self.constraintPars)
    
      self.setDihedralConstraintItemInfo()
  
  def setDihedralConstraintItemPars(self):

    self.constraintPars = {}
    
    if hasattr(self.rawConstraint,'targetAngle'):
      self.constraintPars['targetValue'] = self.rawConstraint.targetAngle
    
    if hasattr(self.rawConstraint,'devAngle'):
      
      angleErrorRange = self.getDihedralAngleErrorRange()
        
      self.constraintPars['error'] = angleErrorRange
      self.constraintPars['upperLimit'] = self.constraintPars['targetValue'] + self.constraintPars['error']
      self.constraintPars['lowerLimit'] = self.constraintPars['targetValue'] - self.constraintPars['error']

    if hasattr(self.rawConstraint,'upperAngle') and hasattr(self.rawConstraint,'lowerAngle'):
      upperLimit = self.rawConstraint.upperAngle
      lowerLimit = self.rawConstraint.lowerAngle
      
      # Make sure that upper is always larger than lower...
      if upperLimit < lowerLimit:
        upperLimit += 360.00

      self.constraintPars['upperLimit'] = upperLimit
      self.constraintPars['lowerLimit'] = lowerLimit
      
  def getDihedralAngleErrorRange(self):
    
    return self.rawConstraint.devAngle

  def setDihedralConstraintItemInfo(self):
  
    pass

  def createConstraintAtomLinks(self):
    
    self.resonancesToLink = []
    
    if self.constraintType != 'dihedral':
    
      for constraint in self.constraintList.constraints:
        for item in constraint.sortedItems():
          for fixedResonance in item.resonances:
            if fixedResonance not in self.resonancesToLink:
              self.resonancesToLink.append(fixedResonance)
                         
    else:
    
      for constraint in self.constraintList.constraints:
        for fixedResonance in constraint.resonances:
          if fixedResonance not in self.resonancesToLink:
            self.resonancesToLink.append(fixedResonance)
                          
    self.getResonanceAtomLinks()

  def createPeakContribAtomLinks(self):
    
    self.resonancesToLink = []
    
    for peakList in self.peakLists:
  
      for peak in peakList.sortedPeaks():
      
        for peakDim in peak.sortedPeakDims():
        
          for peakDimContrib in peakDim.peakDimContribs:
          
            if isinstance(peakDimContrib,Nmr.PeakDimContrib):

              resonance = peakDimContrib.resonance
        
              if resonance not in self.resonancesToLink:

                self.resonancesToLink.append(resonance)

            else:
            
              print "  Not handling PeakDimContribN items"
    
    # Special cases, like XEASY with CYANA, force a specific naming system
    self.forceSpecificNamingSystemName()
              
    # Now set the info for each resonance with regard to atom, residue, ... name
    self.getResonanceAtomLinks()

  def getFormatNamingSystemName(self, format = None, version = None):
  
    #
    # This can be customized format specifically if desired (selection or such)
    #
    
    if not format:
      format = self.format
    if not version and hasattr(self,'version'):
      version = self.version
    
    namingSystemName = format
    
    if hasattr(self,'forceNamingSystemName') and self.forceNamingSystemName and self.forceNamingSystemName != 'None':
      namingSystemName = self.forceNamingSystemName
    
    elif formatNamingDict.has_key((format,version)):
      namingSystemName = formatNamingDict[(format,version)]

    elif formatNamingDict.has_key((format,None)):
      namingSystemName = formatNamingDict[(format,None)]
      
    if self.verbose:
      print "  Selecting naming system %s." % namingSystemName
  
    return namingSystemName
  
  def peakListAssignmentCheck(self,peakLists):
  
    #
    # Only relevant if assignment not written out in chain.seqcode.atom type form
    # (e.g. for XEasy where atomSerials are used)
    #
    
    return True

  def setResonanceToAtoms(self,resonance,*args,**keywds):
  
    if not self.resonanceToAtoms.has_key(resonance):
    
      self.resonanceToAtoms[resonance] = []
      
    self.resonanceToAtoms[resonance].append(ResonanceToAtoms(resonance,*args,**keywds))
            
  def getResonanceAtomLinks(self):
    
    #for tag in ('compressResonances', 'individualAtomsIfNoSet', 'minimalPrompts', 
    #            'useCcpnChainInfo', 'presetResonanceMapping'):
    #   if (hasattr(self, tag)):
    #     print '###', tag, getattr(self,tag)
    
    # Return if new resonanceToAtoms not needed
    if (hasattr(self, 'resonanceToAtoms') and self.resonanceToAtoms
        and self.presetResonanceMapping):
      # We must use pre-existing mapping
      print '### No map generation - using preset (%s resonances)' % len(self.resonanceToAtoms)
      return
    
    #
    # Initialize some vars
    # 
    
    self.resonanceToAtoms = {}
    preferredAtomNameDict = {}
  
    #
    # Get relevant naming system
    # 
    
    self.namingSystemName = self.getFormatNamingSystemName()
    
    
    #
    # This code checks whether original resonance names should be used.
    # Is now based on namingSystem (format does not matter). Only
    # certain naming systems are allowed for each particular export format!
    #
    
    if self.useOriginalResNames:
    
      self.allowedFormats = []
      self.allowedFormatChainMappings = {}
      
      nameMappings = None
      
      if self.project.currentNameMappingStore:
    
        nameMappings = self.project.currentNameMappingStore.findAllNameMappings(namingSystem = self.namingSystemName, isOriginalImport = True)

        if not nameMappings:

          #
          # Try naming system 'Unknown'. Is a gamble but using original names always is!
          #

          nameMappings = self.project.currentNameMappingStore.findAllNameMappings(namingSystem = 'Unknown', isOriginalImport = True)

          if nameMappings:

            for nameMapping in nameMappings:
              print "  Warning: using original format %s with unknown naming system" % nameMapping.formatName
      
      #
      # Set allowed formats...
      #
      
      if nameMappings:
        
        for nameMapping in nameMappings:
        
          if nameMapping.chainMappings:

            self.allowedFormats.append(nameMapping.formatName)
            self.allowedFormatChainMappings[nameMapping.formatName] = nameMapping.chainMappings
      
      if not self.allowedFormats:
      
        print "  Warning: no allowed existing formats found for %s export in %s naming system - original names cannot be used" % (self.format,self.namingSystemName)
    
    #
    # Here the correct molecular system is chosen: all info has to come from
    # within ONE molecular system (otherwise makes no sense)
    #
    # Only do continue if a molecular system is available, otherwise don't bother
    #
    
    if len(self.project.molSystems) == 0:
      return
    
    self.selectMolSystem(resonanceList = self.resonancesToLink)
    
    #
    # Loop over all resonances to be exported
    #
    
    # Added sorting to make calculation (and output) deterministic
    for resonance in sorted(self.resonancesToLink,
                            key=operator.attrgetter('serial')):
    
      #
      # Set text for identifying resonance in screen output
      #
    
      resonanceNameText = getResNameText(resonance)
    
      #
      # First handle case for reusing original resonance names. Note
      # that some chain mapping must have been done for this to work.
      #
    
      if self.useOriginalResNames:

        matchFound = 0
     
        for allowedFormat in self.allowedFormats:
          
          # TODO: not entirely correct... could in principle be different naming system
          # TODO: also set naming system for each resonance when being linked?!?! But
          #       again defeats part of the purpose: often these are not linked...
          resNames = getApplResNames(self.format,[resonance]).keys()
      
          if resNames:
          
            (seqCodeToIdMapping,seqCodeToChainMapping,seqCodes) = self.getChainMappingsSeqCodeInfo(self.allowedFormatChainMappings[allowedFormat])
            
            for resName in resNames:
        
              (formatChainCode,formatSeqCode,formatSpinSystemId,formatSeqInsertCode,atomName) = getNameInfo(resName)

              for (firstSeqCodeKey,lastSeqCodeKey) in seqCodes:
                (firstSeqCode,seqInsertCode) = firstSeqCodeKey
                (lastSeqCode,seqInsertCode) = lastSeqCodeKey
        
                if formatSeqCode >= firstSeqCode and formatSeqInsertCode == seqInsertCode and (lastSeqCode == None or formatSeqCode <= lastSeqCode):
                  chainMapping = seqCodeToChainMapping[(firstSeqCode,seqInsertCode)]
                  
                  if chainMapping.formatChainCode == formatChainCode:
                    chain = self.molSystem.findFirstChain(code = chainMapping.chainCode)
                    seqLength = len(chain.residues)

                    if chain:
                      matchFound = 1
                      seqId = formatSeqCode - chainMapping.formatFirstSeqCode + chainMapping.firstSeqId
                      self.setResonanceToAtoms(resonance,chain,seqId,atomName,originalResName = resName)

        if not matchFound:
          
          print "  Warning: no original name found for %s. Ignored." % resonanceNameText

      #
      # If no original names are being used for this
      # resonance, or if not treated before, then update self.resonanceToAtoms
      # through interpretation of the resonanceSet-AtomSet-Atom links
      #
            
      if not self.useOriginalResNames and not self.resonanceToAtoms.has_key(resonance):
      
        #
        # First need a resonanceSet
        #

        resonanceSet = None

        if resonance is not None:
          resonanceSet = resonance.resonanceSet
        
        #
        # If there is a resonanceSet, do interpretation
        #
        
        if resonanceSet:
          
          #
          # Keep track of atomSets for this resonanceSet
          #
          
          atomSets = resonanceSet.sortedAtomSets()
          
              
          # Rasmus Fogh 27/6/12 replaces previous code
          
          resSetResonances = resonanceSet.sortedResonances()
          
          # get common residue and naming system
          allAtoms = set(x for y in atomSets for x in y.atoms)
          residues = set(x.residue for x in allAtoms)
          rr = residues.pop()
          if residues:
            # atoms from different residues
            residue = None
          else:
            # one residue only
            residue = rr
            
          namingSystem = rr.chemCompVar.chemComp.findFirstNamingSystem(
                          name=self.namingSystemName)
                    
          if not namingSystem:
            print "  Warning: No %s naming system available for resonance %s" % (self.namingSystemName,getResNameText(resonance)),

            if hasattr(self,'ignoreChemCompsWithNoSysNames') and self.ignoreChemCompsWithNoSysNames:
              print " - not written."
              continue
            
            else:
              print " - trying to use default atom names."
          
          chemAtomSet = getChemAtomSetFromAtomSets(atomSets)
          
          if (len(resSetResonances) == 1 and residue and 
              allAtoms == residue.atoms):
            # Residue-level constraint. Set it.
            self.setResonanceToAtoms(resonance, residue.chain, residue.seqId,
                                     None, atomSetName=None)
          
          elif residue and chemAtomSet and len(atomSets) == 2 and len(resSetResonances) <= 2:
            
            #
            # Case where 1/2 resonances, 2 atomSets, one chemAtomSet
            #
            
            # Rasmus Fogh 27/6/12 end replacement block
            
            seqId = residue.seqId
            chain = residue.chain
 
            exportName = getResName(chain.code,residue.seqCode,chemAtomSet.name)
            
            if (hasattr(self,'compressResonances') and self.compressResonances 
                and self.individualAtomsIfNoSet):
              atomSysName = self.getChemAtomOrSetSysName((chemAtomSet, namingSystem))
              if atomSysName:
                chemAtomSetName = atomSysName.sysName
              else:
                chemAtomSetName = None
                print "  Warning: no chemAtomSetSysName found for %s ... ChemAtomSet not used." % exportName
              
            else:
              atomSysName = self.getChemAtomOrSetSysName((chemAtomSet, namingSystem), 
                                                         exportName, findIupac=True)                        

              if atomSysName:                                                                                  
                chemAtomSetName = atomSysName.sysName                                                          

              else:                                                                                            
                chemAtomSetName = chemAtomSet.name                                                             
                print "  Warning: no chemAtomSetSysName found for %s ... using CCPN name." % exportName   
              
            #print "check: Same chemAtomSet %s (%s) set" % (chemAtomSet.name,chemAtomSet.chemComp.code3Letter)

            #
            # Make a list of chemAtoms or chemAtomSets
            #           
            
            if chemAtomSet.chemAtomSets:
              chemAtomOrSetList = chemAtomSet.sortedChemAtomSets()
            else:
              chemAtomOrSetList = chemAtomSet.sortedChemAtoms()
              
            chemAtomOrSetNames = [caos.name for caos in chemAtomOrSetList]
            #chemAtomOrSetNames.sort() # Objects already sorted by name

            #
            # See which names are 'preferred'
            #
          
            preferredAtomNames = []
            
            # test mapping matches first, to guard against incorrect appdata
            # that gives the same name for both resonances.
            mappingMatches = None
            ll = [x.findFirstApplicationData(keyword='mappingSingleMatch')
                  for x in resSetResonances]
            if None not in ll:
              ll = [x.value for x in ll]
              if len(ll) == len(set(ll)):
                # no duplicates - use mapping data
                mappingMatches = ll
            
            for ii,resonance in enumerate(resSetResonances):
            
              #
              # Already mapped as covalently bound prochiral resonance?
              #
              
              if preferredAtomNameDict.has_key(resonance):
                
                preferredAtomNames.append(preferredAtomNameDict[resonance])
                
              else:
          
                #
                # Check if original mapping was traced
                #
                
                if mappingMatches:
                  preferredAtomNames.append(mappingMatches[ii])
                  

                #mappingSingleMatch = resonance.findFirstApplicationData(keyword = 'mappingSingleMatch')

                #if mappingSingleMatch:

                #  preferredAtomNames.append(mappingSingleMatch.value)

                #
                # If not then check if original name matches anything
                #

                else:
                  index = getAmbigNamingPriority(resonance, chemAtomSet)
                  preferredAtomNames.append(chemAtomOrSetNames[index])
                  
            
            if len(set(preferredAtomNames)) != len(preferredAtomNames):
              raise Exception(
                """Inconsistent assignment: 
resonances %s map onto atom names %s
Cannot map atom names for export.
""" % (([x.serial for x in resSetResonances]), preferredAtomNames))
            #
            # Reorganise...
            #
            
            if chemAtomOrSetList[-1].name == preferredAtomNames[0]:
              chemAtomOrSetList.reverse()
            

            #
            # Now loop over chemAtoms or chemAtomSets
            #
            
            for chemAtomOrSet in chemAtomOrSetList:
              
              chemAtomOrSetName = chemAtomOrSet.name
              
              if chemAtomOrSetName in preferredAtomNames:
                resonanceIndex = preferredAtomNames.index(chemAtomOrSetName)
                resonance = resSetResonances[resonanceIndex]
                del resSetResonances[resonanceIndex]
                del preferredAtomNames[resonanceIndex]                  
              
                #print "        set as preferred %s" % (chemAtomOrSetName)
                     
              else:
                #print "        set as default %s" % (chemAtomOrSetName)
                
                validIndices = range(len(preferredAtomNames))
                
                for tempChemAtomOrSet in chemAtomOrSetList:
                  ss = tempChemAtomOrSet.name
                  if ss in preferredAtomNames:
                    removeIndex = preferredAtomNames.index(ss)
                    validIndices.pop(validIndices.index(removeIndex))
                 
                if validIndices:
                  useIndex = validIndices[0]
                  resonance = resSetResonances.pop(useIndex)
                  del preferredAtomNames[useIndex]
              
              
              #
              # Also use covalently bound if available. Note that uses Analysis appdata tag as well (higher up) if available,
              # but will likely be removed in the future (19/07/2006)
              #
              
              if hasattr(resonance,'covalentlyBound') and resonance.covalentlyBound:  # Temporary fix because FixedResonance has no covalentlyBound
                
                #
                # Only relevant for 'deep' prochiral case (e.g. Valine methyls)...
                #
                
                for covalentResonance in resonance.covalentlyBound:
                
                  covalentResSet = covalentResonance.resonanceSet
                
                  if covalentResSet:
                  
                    for covalentAtomSet in covalentResSet.atomSets:
                  
                      refChemAtom = covalentAtomSet.findFirstAtom().chemAtom
                      
                      if not refChemAtom.chemAtomSet:
                        continue
                                              
                      covalentPreferredAtomName = None

                      for refChemBond in refChemAtom.chemBonds:
                        
                        refChemBondAtoms = set(refChemBond.chemAtoms)
                        refChemBondAtoms.remove(refChemAtom)
                        otherRefChemAtom = refChemBondAtoms.pop()
                        
                        if isinstance(otherRefChemAtom,ChemComp.ChemAtom):

                          if chemAtomOrSet == otherRefChemAtom.chemAtomSet:
                            covalentPreferredAtomName = refChemAtom.name 

                          elif chemAtomOrSet == otherRefChemAtom:
                            covalentPreferredAtomName = refChemAtom.chemAtomSet.name
                            
                          if covalentPreferredAtomName:

                            preferredAtomNameDict[covalentResonance] = covalentPreferredAtomName
                            break

                      if covalentPreferredAtomName:
                        break
              
              #
              # Get the chemAtom(Set)SysNames
              #
              
              tempChemAtomOrSetList = [chemAtomOrSet]
              
              # Decide whether to split this list into individual component Atom(Set)s
              if chemAtomOrSet.className == 'ChemAtomSet':
                splitChemAtomSet = False
                if self.individualAtoms:
                  splitChemAtomSet = 'deep'
                elif self.individualAtomsIfNoSet:
                  if not self.getChemAtomOrSetSysName((chemAtomOrSet, namingSystem)):
                    splitChemAtomSet = 'shallow'
              
                if splitChemAtomSet:
                  # split to direct content of ChemAtom(Set)s
                  tempChemAtomOrSetList = (chemAtomOrSet.sortedChemAtomSets() or
                                           chemAtomOrSet.sortedChemAtoms())
                  if (splitChemAtomSet == 'deep'
                      and tempChemAtomOrSetList[0].className == 'ChemAtomSet'):
                    # split all the way down to individual atoms
                    ii = 0
                    while ii < len(tempChemAtomOrSetList):
                      caos = tempChemAtomOrSetList[ii]
                      if caos.className == 'ChemAtomSet':
                        del tempChemAtomOrSetList[ii]
                        tempChemAtomOrSetList.extend(caos.sortedChemAtomSets())
                        tempChemAtomOrSetList.extend(caos.sortedChemAtoms())
                      else:
                        # we have a ChemAtom here. 
                        # We could probably break, but just to be safe
                        ii += 1
                
              for tempChemAtomOrSet in tempChemAtomOrSetList:
               
                atomSysName = self.getChemAtomOrSetSysName((tempChemAtomOrSet, namingSystem), 
                                                           exportName, findIupac=True)
  
                if not atomSysName:
                
                  # If really nothing, use CCPN atom name
                  refSysName = tempChemAtomOrSet.name
                  print "  Warning: no %s chemAtomSysName found... using CCPN name." % (tempChemAtomOrSet.name)
  
                else:
                  refSysName = atomSysName.sysName
                #print refSysName,chemAtomSetName,chemAtomOrSet
                #print resonance.findFirstApplicationData(application ='Analysis',keyword='ambigProchiralLabel').value
                self.setResonanceToAtoms(resonance,chain,seqId,refSysName,atomSetName=chemAtomSetName, 
                                         chemAtomOrSet=tempChemAtomOrSet, atomType=tempChemAtomOrSet.elementSymbol)
                                
                for otherResonance in resonanceSet.resonances:
                  if otherResonance != resonance:
                    self.resonanceToAtoms[resonance][-1].addOtherGroupResonance(otherResonance)
              
              #
              # Get out of loop if no resonances left
              # TODO: This means that if case 1 resonance via resset to two
              # prochiral atomsets, it will be listed in output as one of those
              # atomsets... not as e.g. HB*. Is that correct?!?!?
              
              if not resSetResonances:
                #print "OUTOFLOOP"
                break
            
          
          elif residue and len(atomSets) == 1 and len(resSetResonances) == 1:
            
            #
            # 1 resonance, 1 atomset. Simplest case.
            # Needs special handling to deal with isopropyl groups
            #
            
            atomSet = atomSets[0]
            resonance = resSetResonances[0]
            
            seqId = residue.seqId
            chain = residue.chain
 
            exportName = getResName(chain.code,residue.seqCode,
                                    atomSet.name or atomSet.findFirstAtom().name)
                                    
            
            chemAtom = atomSet.findFirstAtom().chemAtom
            
            # NB all atoms should ALWAYS point to the same ChemAtomSet, 
            # but it costs nothing ot be careful.
            aSet = set(x.chemAtom.chemAtomSet for x in atomSet.atoms)
            if len(aSet) == 1:
              chemAtomSet = aSet.pop()
            else:
              chemAtomSet = None
            
            if len(atomSet.atoms) == 1:
              
              atomSysName = self.getChemAtomOrSetSysName((chemAtom,namingSystem))
              if atomSysName:
                refSysName = atomSysName.sysName
              else:
                refSysName = chemAtom.name
            
              if chemAtomSet:
                if (hasattr(self,'compressResonances') and self.compressResonances 
                    and self.individualAtomsIfNoSet):
                  chemAtomSetSysName = self.getChemAtomOrSetSysName((chemAtomSet, 
                                                                     namingSystem))
                  if chemAtomSetSysName:
                    chemAtomSetName = chemAtomSetSysName.sysName
                  else:
                    chemAtomSet = None
                    chemAtomSetName = None
                    print "  Warning: no chemAtomSetSysName found for %s ... ChemAtomSet not used." % exportName
                  
                else:
                  chemAtomSetSysName = self.getChemAtomOrSetSysName((chemAtomSet, 
                                                                     namingSystem), 
                                                                    exportName, 
                                                                    findIupac=True)    
                  if chemAtomSetSysName:                                                                                  
                    chemAtomSetName = chemAtomSetSysName.sysName 
                  else:                                                                                            
                    chemAtomSetName = chemAtomSet.name                                                             
                    print "  Warning: no chemAtomSetSysName found for %s ... using CCPN name." % exportName   
  
              else:  
                chemAtomSetName = None
              
            elif chemAtomSet and not self.individualAtoms:
              #assert chemAtomSet and len(atomSet.atoms) > 1
              
              atomSysName = self.getChemAtomOrSetSysName((chemAtomSet,namingSystem))
              if atomSysName:
                refSysName = atomSysName.sysName
              else:
                refSysName = chemAtom.name
              
              deepSet = chemAtomSet.chemAtomSet
              if deepSet:
                chemAtomSet = deepSet
                
              chemAtomSetSysName = self.getChemAtomOrSetSysName((chemAtomSet,
                                                                 namingSystem))
              if chemAtomSetSysName:
                chemAtomSetName = chemAtomSetSysName.sysName
              else:
                chemAtomSetName = chemAtomSet.name
              
            else:
              #assert (individualAtoms or not chemAtomSet) and len(atomSet.atoms) > 1
              # Deal with atoms individualy
              # chemAtomSet should always exist here, but it costs nothing to be prepared.
                                
              chemAtomSetName = None
              if chemAtomSet:
                chemAtomSetSysName = self.getChemAtomOrSetSysName((chemAtomSet,
                                                                   namingSystem))
                if chemAtomSetSysName:
                  chemAtomSetName = chemAtomSetSysName.sysName
                else:
                  chemAtomSetName = chemAtomSet.name
              
              
              for chemAtom in chemAtomSet.chemAtoms:

                exportName = getResName(chain.code,residue.seqCode,chemAtom.name)
                atomSysName = self.getChemAtomOrSetSysName((chemAtom, namingSystem),
                                                           exportName, findIupac=True)
             
                if atomSysName:
                  refSysName = atomSysName.sysName
                else:
                  print "  Warning: no %s sysName for %s... using CCPN name." % (self.namingSystemName,exportName)
                  refSysName = chemAtom.name
                                    
                chemAtomOrSet = chemAtomSet or chemAtom
                self.setResonanceToAtoms(resonance, chain, seqId, refSysName,
                                         atomSetName=chemAtomSetName, 
                                         chemAtomOrSet=chemAtomOrSet, 
                                         atomType=chemAtomOrSet.elementSymbol)
              
              # Done with this resonance
              continue
              
            # set resonanceToAtoms
            chemAtomOrSet = chemAtomSet or chemAtom
            self.setResonanceToAtoms(resonance, chain, seqId, refSysName,
                                     atomSetName=chemAtomSetName, 
                                     chemAtomOrSet=chemAtomOrSet, 
                                     atomType=chemAtomOrSet.elementSymbol)
              
              
          #
          # This is the generic case: ambiguous and not prochiral.
          #
          
          else:
            
            #
            # Treat everything the same: use chemAtomSet if multiple atoms for
            # chemAtomSet, otherwise use chemAtom.
            # Ambiguous so link whatever to whatever!
            # Note that self.resonanceToAtoms is a list of ResonanceToAtoms
            # objects, so in case ambiguity will have multiple members.
            #
            
            #if len(atomSets) == 1:
            #  print "check: Single"
            #else:
            #  print "check: Ambiguous"

            for atomSet in atomSets:
            
              residue = atomSet.findFirstAtom().residue
              seqId = residue.seqId
              chain = residue.chain       
              
              sysNamesLink = []
              chemAtomSet = None
              chemAtomSetName = None
              
              #
              # Check whether sys names for individual atom(s) or for chemAtomSet. Most of this not relevant when exporting all individual atoms.
              #

              if len(atomSet.atoms) > 1:
            
                chemAtomSet = atomSet.findFirstAtom().chemAtom.chemAtomSet
                chemAtomSetName = chemAtomSet.name
                
                #
                # If there is a deep chemAtomSet, try to use this. 
                #
                
                deepSet = chemAtomSet.chemAtomSet
                
                #if chemAtomSet.chemAtomSet and not self.individualAtoms:
                if deepSet and not self.individualAtoms:
                  
                  #
                  # First check if valid...
                  #

                  #chemAtomsFound = checkChemAtomsForDeepSet(chemAtomSet.chemAtomSet,atomSets)
                  chemAtomsFound = checkChemAtomsForDeepSet(deepSet,atomSets)
                  
                  if chemAtomsFound and namingSystem:
                    
                    #
                    # This code might cause some problems for HD*|HE* for aromatics...
                    # Works fine if namingSystem DIANA - will use QR
                    # Never use this name as reference - better to split up!
                    #

                    #chemAtomSetSysName = self.getChemAtomOrSetSysName((chemAtomSet,namingSystem))
                    chemAtomSetSysName = self.getChemAtomOrSetSysName((deepSet,namingSystem))

                    if chemAtomSetSysName:
                      #sysNamesLink = [(chemAtomSetSysName.sysName,chemAtomSet)]
                      sysNamesLink = [(chemAtomSetSysName.sysName,deepSet)]
                      chemAtomSet = None
                      chemAtomSetName = chemAtomSetSysName.sysName
                
                #
                # Find the chemAtomSet sysNames (if still relevant after above)
                #
                
                exportName = getResName(chain.code,residue.seqCode,chemAtomSetName)
                
                if not sysNamesLink and not self.individualAtoms:

                  if self.individualAtomsIfNoSet:
                    atomSysName = self.getChemAtomOrSetSysName((chemAtomSet, namingSystem))
                  else:
                    atomSysName = self.getChemAtomOrSetSysName((chemAtomSet, namingSystem), 
                                                                exportName, findIupac=True)
 
                  if atomSysName:
                    sysNamesLink = [(atomSysName.sysName,chemAtomSet)]
                      
                  else:
                    if chemAtomSet.chemAtomSet and not self.individualAtomsIfNoSet:
                      print ("  Warning: using CCPN name for atomSet '%s'..." 
                             % exportName)
                      sysNamesLink = [(chemAtomSet.name,chemAtomSet)]

                    else:
                      print ("  Warning: no %s sysNames for %s... using individual atom name(s)." 
                             % (self.namingSystemName,exportName))
              
              #
              # If no chemAtomSet sysName(s) were found, do it on the individual atom level...
              #
              
              if not sysNamesLink:
                
                if chemAtomSet:
                  chemAtoms = chemAtomSet.chemAtoms
                else:
                  #chemAtoms = [chemAtom] # NB if we have more than one, we surely want them all
                  chemAtoms = [x.chemAtom for x in atomSet.sortedAtoms()]
                
                for chemAtom in chemAtoms:

                  exportName = getResName(chain.code,residue.seqCode,chemAtom.name)
                  atomSysName = self.getChemAtomOrSetSysName((chemAtom, namingSystem),
                                                              exportName, findIupac=True)
             
                  if atomSysName:
                    sysNamesLink.append((atomSysName.sysName,chemAtom))
                  else:
                    sysNamesLink.append((chemAtom.name,chemAtom))
                    print "  Warning: no %s sysName for %s... using CCPN name." % (self.namingSystemName,exportName)
                    
                #
                # Still set the chemAtomSetName if there is one, for compressing resonances
                #
                
                if (not chemAtomSetName and namingSystem and not self.individualAtoms 
                    and chemAtom.chemAtomSet and hasattr(self,'compressResonances') 
                    and self.compressResonances):
                  
                  chemAtomSetSysName = self.getChemAtomOrSetSysName((chemAtom.chemAtomSet,namingSystem))
                  
                  if chemAtomSetSysName:
                    chemAtomSetName = chemAtomSetSysName.sysName

              sysNamesLink.sort()
              for resonance in resSetResonances:
                for (sysName,chemAtomOrSet) in sysNamesLink:
                  self.setResonanceToAtoms(resonance, chain, seqId, sysName,
                                           atomSetName = chemAtomSetName, 
                                           chemAtomOrSet = chemAtomOrSet, 
                                           atomType = chemAtomOrSet.elementSymbol)
                
            #
            # Set some additional information
            #
            
            for resonance in resSetResonances:
              for otherResonance in resonanceSet.resonances:
                if otherResonance != resonance:
                  self.resonanceToAtoms[resonance][-1].addOtherLinkedResonance(otherResonance)
                
        #
        # If there is no resonanceSet for the resonance, then ignore
        # it: impossible to export (need linkResonances first)
        #
                  
        else:
          
          resonanceNameText = getResNameText(resonance)
          print "  Warning: No resonanceSet for %s with name %s - ignored" % (resonance, resonanceNameText)
    
    #
    # Debug code...
    # 
    
    #for res in self.resonanceToAtoms.keys():
    
    #  for rta in self.resonanceToAtoms[res]:

    #    print rta.resonance
    #    print rta.chain, rta.seqId, rta.atomName, rta.atomSetName
    #    print rta.originalResName
    #    print


  def setResonanceToSpinSystem(self,resonance,*args,**keywds):
  
    if not self.resonanceToSpinSystems.has_key(resonance):
    
      self.resonanceToSpinSystems[resonance] = []
      
    self.resonanceToSpinSystems[resonance].append(ResonanceToSpinSystems(resonance,*args,**keywds))

  def getResonanceSpinSystemLinks(self):
  
    #
    # Get relevant naming system
    # 
    
    self.namingSystemName = self.getFormatNamingSystemName()

    #
    # Initialize some vars
    # 
    
    self.resonanceToSpinSystems = {}
    
    # WARNING: removed originalNames stuff - not relevant for spin systems?
       
    #
    # Loop over all resonances to be exported
    #
    
    # Added sorting to make result (and output) deterministic
    for resonance in sorted(self.resonancesToLink,
                            key=operator.attrgetter('serial')):
      
      assignNames = resonance.assignNames
       
      if not assignNames and resonance.name:
          assignNames = [resonance.name]
      if not assignNames:
        continue

      if not self.resonanceToSpinSystems.has_key(resonance):

        resonanceGroup = resonance.resonanceGroup
        
        if not resonanceGroup and resonance.resonanceProbs:
          for resonanceProb in resonance.resonanceProbs:
            if resonanceProb.weight == 1.0:
              resonanceGroup = resonanceProb.possibility
              break
        
        if not resonanceGroup:
          continue
        
        if self.chains:
          chainFound = True
          if resonanceGroup.chains:
            chainFound = False
            for chain in resonanceGroup.chains:
              if chain in self.chains:
                chainFound = True
                break
          
          if not chainFound and resonanceGroup.residue:
            if resonanceGroup.residue.chain in self.chains:
              chainFound = True
          
          if not chainFound and resonanceGroup.residueProbs:
            for residueProb in resonanceGroup.residueProbs:
              residue = residueProb.possibility
              if residue.chain in self.chains:
                chainFound = True
                break
                      
          if not chainFound:
            continue
        
        bestLink = {}        
        for resonanceGroupProb in resonanceGroup.findAllResonanceGroupProbs(linkType='sequential'):
          
          seqOffset = resonanceGroupProb.sequenceOffset
          weight = resonanceGroupProb.weight
          
          if not bestLink.has_key(seqOffset):
            bestLink[seqOffset] = (0,None)
          
          if weight > bestLink[seqOffset][0]:
            bestLink[seqOffset] = (weight,resonanceGroupProb.possibility)  
        
        linkedResonanceGroups = {}
        for seqOffset in bestLink.keys():
          linkedResonanceGroups[seqOffset] = bestLink[seqOffset][1]
        
        self.setResonanceToSpinSystem(resonance,resonanceGroup,assignNames, linkedSpinSystems = linkedResonanceGroups) #atomSetName = chemAtomSetName, chemAtomOrSet = chemAtomOrSet, atomType = refChemAtom.elementSymbol)
        

  def sortConstraintsForExport(self,constraints,getResonancesFunc):
    
    sortDict = {}
   
    for constraint in constraints:

      sortKeys = []
      
      for resonance in getResonancesFunc(constraint):
        sortKey = ('XXXX',99999999)
      
        if self.resonanceToAtoms.has_key(resonance):
    
          atomId = self.resonanceToAtoms[resonance][0].getAtomId()        
          sortKey = (atomId[0].code,atomId[1])
          
        sortKeys.append(sortKey)
      
      sortKeys = tuple(sortKeys)
      
      if not sortDict.has_key(sortKeys):
        sortDict[sortKeys] = []
      
      sortDict[sortKeys].append(constraint)
    
    sortKeysList = sortDict.keys()
    sortKeysList.sort()
    
    constraints = []
    
    for sortKeys in sortKeysList:
      constraints.extend(sortDict[sortKeys])
      
    return constraints
  
  def getConstraintLevelResonances(self,constraint):
    
    return constraint.resonances    
  
  def getConstraintItemLevelResonances(self,constraint):
    
    resonances = [None,None]
    items = constraint.sortedItems()
    
    if items:
      resonances = getResonancesFromPairwiseConstraintItem(items[0])
    
    return resonances
    
  def setHBondConstraintFileConstraints(self):
  
    """
    Call code for distance constraints: is mostly same. Differences are contained within
    variable self.setRaw... defs
    """
    
    self.setGenericConstraintFileConstraints(constraintApiCode = 'Distance')  

  def setDistanceConstraintFileConstraints(self):
  
    self.setGenericConstraintFileConstraints()

  def setGenericConstraintFileConstraints(self, constraintApiCode = None):
  
    if not constraintApiCode:
      
      constraintApiCode = self.constraintApiCode
    
    setRawConstraintFunc = getattr(self,"setRaw%sConstraint" % constraintApiCode)
    
    """
    Code to set the constraint objects for the export file.
    Only works for output for 2 resonances!
    """

    #
    # Possibly sort by chain code and sequence code of resonances (in order) for first item.
    #
    
    constraints = self.constraintList.sortedConstraints()
    
    if self.sortForExport:
      
      constraints = self.sortConstraintsForExport(constraints,self.getConstraintItemLevelResonances)
    
    #
    # Loop over all constraints, items and resonances
    #
    
    self.resSetNames = 2 * ['']
      
    for self.constraint in constraints:
      
      if self.constraint.items:
        setRawConstraintFunc()

      else:
        print "  Warning: no items for %s constraint %d" % (constraintApiCode.lower(),self.constraint.serial)
        continue
        
      #print self.constraint.serial
      
      itemPairs = []
      
      otherResonance = 2 * ['']
      
      it = 0
      
      constraintItemList = list(self.constraint.sortedItems()[:])
      
      while (it < len(constraintItemList)):
        
        self.item = constraintItemList[it]
        resSetNames = 2 * [None]
        
        #
        # Check whether there are other items that can be combined based on chain, seqId,
        # and atomName/atomSetName
        #
        # First do it for the first resonance, then for the second. Only combine
        # items based on first resonance if the second is exactly the same.
        # Same for second. This way ambiguities should be handled correctly...
        #
        
        itemResonances = getResonancesFromPairwiseConstraintItem(self.item)
        
        for x in range(2):

          #
          # Track original status if necessary...
          #
          
          if not self.compressResonances and self.item.findFirstApplicationData(application = 'ccp', keyword = 'copied', value = x):
          
            isCopied = True
            
          else:
        
            isCopied = False
            
          #
          # Check resonance
          #
          
          fixedResonance = itemResonances[x]
            
          y = not x

          if self.resonanceToAtoms.has_key(fixedResonance):
          
            atomSetId = self.resonanceToAtoms[fixedResonance][0].getAtomSetId()
          
            if len(self.resonanceToAtoms[fixedResonance]) == 1 and atomSetId:
                        
              #
              # Check out the other items...
              #
           
              tempIt = it + 1
            
              while (tempIt < len(constraintItemList)):
            
                tempItem = constraintItemList[tempIt]
            
                matchFixedResonance = getResonancesFromPairwiseConstraintItem(tempItem)[x]
                  
                if self.resonanceToAtoms.has_key(matchFixedResonance):

                  matchAtomSetId = self.resonanceToAtoms[matchFixedResonance][0].getAtomSetId()
                                
                  if matchAtomSetId and len(self.resonanceToAtoms[matchFixedResonance]) == 1 and atomSetId[:-1] == matchAtomSetId[:-1]:
                
                    #
                    # Match for the atomSetId... now check if other side matches exactly
                    #
                  
                    otherFixedResonance = itemResonances[y]                    
                    matchOtherFixedResonance = getResonancesFromPairwiseConstraintItem(tempItem)[y]
                  
                    if otherFixedResonance.resonanceSet and otherFixedResonance == matchOtherFixedResonance:

                      #
                      # If self.compressResonances is default then output depends on
                      # whether an item was duplicated... should retain original 'compression' 
                      #
                      
                      if not self.compressResonances:
  
                        #
                        # Check whether was copied
                        #

                        if not isCopied and tempItem.findFirstApplicationData(application = 'ccp', keyword = 'copied', value = x):

                          isCopied = True

                      #
                      # Compress the items if compressResonances is on or
                      # if it's None and isCopied is True
                      #
                      
                      if self.compressResonances or (self.compressResonances == None and isCopied):
                      
                        otherResidue = otherFixedResonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
                        resSetNames[x] = [atomSetId]
                        chainCodes = [atomSetId[0].code,
                                      otherResidue.chain.code,
                                      matchAtomSetId[0].code]
                        seqIds = [atomSetId[1],
                                  otherResidue.seqId,
                                  matchAtomSetId[1]]


                        print "  %s constraint %d: Removed %s.%d.%s,%s.%d.%s based on %s.%d.%s,%s.%d.%s" % (constraintApiCode,self.constraint.serial,chainCodes[0],seqIds[0],fixedResonance.name,chainCodes[1],seqIds[1],otherFixedResonance.name,chainCodes[2],seqIds[2],matchFixedResonance.name,chainCodes[1],seqIds[1],matchOtherFixedResonance.name)

                        constraintItemList.pop(tempIt)
                        tempIt -= 1
                 
                tempIt += 1

              #
              # Also check that, IF x = 1 and resSetNames[1] is set already,
              # whether tempItem.resonances[y] is in self.resonanceToAtoms
              # and whether its atomSetId matches resSetNames[0]! This is for the
              # case where you have HB1,HB1 HB1,HB2 HB2,HB1 and HB2,HB2 items:
              # the last one would otherwise not be recognized
              #
              # Need to run this loop AFTER all items have been checked (see above)
              #
              
              if x == 1 and resSetNames[1]:
              
                if resSetNames[0] and len(resSetNames[0]) == 1:

                  tempIt = it + 1
            
                  while (tempIt < len(constraintItemList)):
            
                    tempItem = constraintItemList[tempIt]
                    tempItemResonances = getResonancesFromPairwiseConstraintItem(tempItem)
            
                    matchFixedResonance = tempItemResonances[x]

                    if not self.resonanceToAtoms.has_key(matchFixedResonance):
                      tempIt += 1
                      continue
              
                    matchAtomSetId = self.resonanceToAtoms[matchFixedResonance][0].getAtomSetId()
                    matchOtherFixedResonance = tempItemResonances[y]

                    if self.resonanceToAtoms.has_key(matchOtherFixedResonance) and len(self.resonanceToAtoms[matchOtherFixedResonance]) == 1:
                    
                      matchOtherAtomSetId = self.resonanceToAtoms[matchOtherFixedResonance][0].getAtomSetId()
                    
                      if matchAtomSetId and matchOtherAtomSetId and atomSetId[0:2] == matchAtomSetId[0:2] and matchOtherAtomSetId[:-1] == resSetNames[0][0][:-1]:

                        #
                        # If self.compressResonances is default then output depends on
                        # whether an item was duplicated... should retain original 'compression' 
                        #

                        if self.compressResonances == None:

                          #
                          # Check whether was copied
                          #

                          if not isCopied and tempItem.findFirstApplicationData(application = 'ccp', keyword = 'copied', value = x):

                            isCopied = True

                        #
                        # Compress the items if compressResonances is on or
                        # if it's None and isCopied is True
                        #

                        if self.compressResonances or (self.compressResonances == None and isCopied):

                     
                          otherResidue = otherFixedResonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
                          chainCodes = [atomSetId[0].code,
                                        otherResidue.chain.code,
                                        matchAtomSetId[0].code,
                                        otherResidue.chain.code]
                          seqIds = [atomSetId[1],
                                    otherResidue.seqId,
                                    matchAtomSetId[1],
                                    otherResidue.seqId,
                                    ]

                          print "  %s constraint %d: Removed %s.%d.%s,%s.%d.%s based on %s.%d.%s,%s.%d.%s (SET MATCH!)" % (constraintApiCode,self.constraint.serial,chainCodes[0],seqIds[0],fixedResonance.name,chainCodes[1],seqIds[1],otherFixedResonance.name,chainCodes[2],seqIds[2],matchFixedResonance.name,chainCodes[3],seqIds[3],matchOtherFixedResonance.name)

                          constraintItemList.pop(tempIt)
                          tempIt -= 1
                  
                    tempIt += 1
                

            if not resSetNames[x]:
              
              resSetNames[x] = []
              
              for resonanceToAtoms in self.resonanceToAtoms[fixedResonance]:
              
                resSetNames[x].append(resonanceToAtoms.getAtomId())
        
        it += 1

        #
        # Check if matches were found
        #
        
        if None in resSetNames:
        
          textItems = []
          
          for i in range(0,len(resSetNames)):
            
            if resSetNames[i] == None:
            
              resonance = itemResonances[i]

              if resonance:
                textItems.append(getResNameText(resonance))
              else:
                textItems.append('empty')
                
              #
              # Reorganize for exportAll == True (see lower)
              #
                              
              resSetNames[i] = [None]
          
          text = ','.join(textItems)
          
          warningText = "Not handling constraint %d, item %d, resonance(s) '%s' not assigned." % (self.constraint.serial,it,str(text))
          
          if hasattr(self.rawConstraint,'errors'):
            self.rawConstraint.errors.append(warningText)
          
          if self.exportAll == False:

            print "  Warning: %s" % warningText
            continue
          
        #
        # Set the info in the constraintFile - this handles multiple items that result
        # from this particular item (because there is overlap with the other items)!
        #
        
        for atomOrSetId in resSetNames[0]:
        
          if atomOrSetId:
          
            #
            # Use mapping!
            #

            chain = atomOrSetId[0]
            seqId = atomOrSetId[1]
            
            chainCode = self.chainDict[chain][0]
            seqCode = self.getExportSeqCode(self.chainDict[chain][1],chain.findFirstResidue(seqId = seqId))
            atomName = atomOrSetId[2]
          
          else:
          
            chainCode = atomName = ''
            seqCode = -999999
        
          for otherAtomOrSetId in resSetNames[1]:
          
            if [atomOrSetId,otherAtomOrSetId] in itemPairs:
              # DEBUG CODE: should not be necessary any more...
              print "  %s to %s already exists..." % (self.resSetNames[0],self.resSetNames[1])
              continue
            
            if otherAtomOrSetId:
            
              #
              # Use mapping!
              #
              
              otherChain = otherAtomOrSetId[0]
              otherSeqId = otherAtomOrSetId[1]
              
              otherChainCode = self.chainDict[otherChain][0]
              otherSeqCode = self.getExportSeqCode(self.chainDict[otherChain][1],otherChain.findFirstResidue(seqId = otherSeqId))
              otherAtomName = otherAtomOrSetId[2]

            else:

              otherChainCode = otherAtomName = ''
              otherSeqCode = -999999
              
            #
            # Only keep track of itemPairs that don't have None in 'em
            #
            
            if atomOrSetId and otherAtomOrSetId:
            
              itemPairs.append([atomOrSetId,otherAtomOrSetId])
              #print   "       ",self.resSetNames[0] , self.resSetNames[1]
            
            #
            # Always set self.resSetNames
            #
            
            self.resSetNames[0] = getResName(chainCode,seqCode,atomName)
            self.resSetNames[1] = getResName(otherChainCode,otherSeqCode,otherAtomName)
            
            # TODO: this is for nmrStar and Concoord only - should have separate function for this whole bit!!
            self.ccpInfo = 2 * ['']
            self.ccpInfo[0] = atomOrSetId
            self.ccpInfo[1] = otherAtomOrSetId
            
            getattr(self,'setRaw%sConstraintItem' % constraintApiCode)()
            
            # Possible that constraint item creation fails... only relevant if rawConstraintItem exists!
            if not hasattr(self,'rawConstraintItem') or self.rawConstraintItem:
              for funcFormat in ('setRaw%sConstraintItemMembers','setRaw%sConstraintItemInfo'):
                getattr(self,funcFormat % constraintApiCode)()
      
      getattr(self,'setRaw%sConstraintInfo' % constraintApiCode)()

      #
      # Delete constraint if it has no valid items and there is no error support
      #
      
      if hasattr(self.rawConstraint,'errors') and self.rawConstraint.errors:
        delete = False
      else:
        delete = True
      
      getattr(self,'checkRaw%sConstraintItems' % constraintApiCode)(delete = delete)

  def setRawDistanceConstraintInfo(self):
  
    pass

  def setRawDistanceConstraintItemInfo(self):
  
    pass

  def checkRawDistanceConstraintItems(self,delete = True):
  
    return self.checkRawConstraintItems(delete = delete)
    
  def checkRawConstraintItems(self,delete = True):
  
    if not self.rawConstraint.items:
      
      if delete:
        del self.constraintFile.constraints[-1]

      return False
      
    return True
      
  def setRawDistanceConstraintItemMember(self):
    
    pass

  def setRawDistanceConstraintItem(self):
    
    pass

  def setRawDistanceConstraint(self):
  
    pass

  def setDihedralConstraintFileConstraints(self):
  
    """
    Writes out dihedrals. 
    """
    
    #
    # Loop over all constraints, items and resonances
    #
    
    self.setDihedralConstraintFileInfo()
    
    #
    # Possibly sort by chain code and sequence code of resonances (in order) for first item.
    #
    
    constraints = self.constraintList.sortedConstraints()
    
    if self.sortForExport:
      
      constraints = self.sortConstraintsForExport(constraints,self.getConstraintLevelResonances)
    
    #
    # Now loop over the constraints
    #
          
    for self.constraint in constraints:

      if not self.constraint.items:
        print "  Warning: No items for constraint %d" % self.constraint.serial
        continue      
    
      self.resSetNames = self.constraintResNum * ['']
      self.ccpInfo = self.constraintResNum * ['']

      badConstraint = 0
      warningText = ""
      
      constraintResonances = self.constraint.resonances

      for i in range(0,self.constraintResNum):

        resonance = constraintResonances[i]
      
        if self.resonanceToAtoms.has_key(resonance):
        
          atomId = self.resonanceToAtoms[resonance][0].getAtomId()
                  
          # TODO: this is for nmrStar only - should have separate function for this whole bit!!
          self.ccpInfo[i] = atomId

          #
          # Use mapping
          #
          
          chain = atomId[0]
          chainCode = self.chainDict[chain][0]
          
          seqCode = self.getExportSeqCode(self.chainDict[chain][1],chain.findFirstResidue(seqId = atomId[1]))

          self.resSetNames[i] = getResName(chainCode,seqCode,atomId[2])
          
          #
          # TODO: Can't handle ambiguous?
          #

          if len(self.resonanceToAtoms[resonance]) != 1:
            badConstraint = True
            
      #
      # TODO: set warningText from ccp errors as well!!
      #
        
      if badConstraint or '' in self.resSetNames:
        
        textItems = []
        
        for i in range(0,self.constraintResNum):

          if not self.resSetNames[i]:
          
            self.resSetNames[i] = getResName('',-999999,'')
            
            resonance = constraintResonances[i]

            if resonance:
              textItems.append(getResNameText(resonance))
            else:
              textItems.append('empty')
          
        text = string.join(textItems,',')

        warningText = "Not handling dihedral restraint %d, resonance(s) '%s' not assigned." % (self.constraint.serial,text)

        if self.exportAll == False:
        
          badConstraint = True
          print "  Warning: %s" % warningText


      #
      # For dihedrals the items are written out as separate lines in the export...
      # best I can do with this right now (a warning will be given though)
      #

      if len(self.constraint.items) > 1:
        print "  Warning: more than one item for dihedral constraint %d. Written as separate ones." % self.constraint.serial


      for self.item in self.constraint.sortedItems():

        #
        # Set the constraint, items and members info
        # These are created, even if not valid, ONLY if error
        # messages can be written out for them (raw constraint object
        # has .errors attribute)
        #

        self.setRawDihedralConstraint()
        
        delete = True

        if warningText:
        
          if hasattr(self.rawConstraint,'errors'):
          
            self.rawConstraint.errors.append(warningText)
            delete = False

        #
        # Delete constraint if it has no valid items and there's
        # no error support in the output format
        #

        if badConstraint:
        
          if delete:

            del self.constraintFile.constraints[-1]
            
        else:

          self.setRawDihedralConstraintItem()
          self.setRawDihedralConstraintItemMembers()
          self.setRawDihedralConstraintItemInfo()

        self.setRawDihedralConstraintInfo()

  def setDihedralConstraintFileInfo(self):
  
    # To set additional format-specific stuff
  
    pass

  def setRawDihedralConstraintInfo(self):
  
    # To set additional format-specific stuff
    
    pass

  def setRawDihedralConstraintItemInfo(self):
  
    # To set additional format-specific stuff
    
    pass
    
  def setRawDihedralConstraintItemMembers(self):
  
    pass

  def setRawDihedralConstraintItem(self):
    
    pass

  def setRawDihedralConstraint(self):
  
    pass

  def setRdcConstraintFileConstraints(self):
  
    self.setGenericConstraintFileConstraints()

  def setRawRdcConstraintItemMember(self):
    
    pass

  def setRawRdcConstraintItem(self):
    
    pass

  def setRawRdcConstraintItemInfo(self):
    
    pass

  def setRawRdcConstraint(self):
  
    pass

  def setRawRdcConstraintInfo(self):
  
    pass

  def checkRawRdcConstraintItems(self, delete = True):
  
    if not self.rawConstraint.items or len(self.rawConstraint.items) > 1:
      
      if delete:
        del self.constraintFile.constraints[-1]
        print "  Error: only exactly 1 item allowed for rdc constraint %d - deleted." % self.constraint.serial

      return False
      
    return True

  def setConversionError(self,parent,constraintId,errorType,errorText):
    value = [constraintId,errorType,errorText]
    parent.addApplicationData(Implementation.AppDataString(application = self.format, keyword = 'ConversionError', value = str(value)))

  #
  # Defs used in peak assignment
  #  

  def setPeakDimContribs(self,peakAssignDim, peakDims = None, resonances = None):
  
    if len(self.assignCombination) == peakAssignDim:
      #print peakAssignDim, peakDims
      #print resonances
              
      #
      # First check if peakContrib doesn't exist already
      # Note that you don't have to have ALL peak dimensions for this!
      #
      
      checkPeakDimContribs = []

      for i in range(0,len(peakDims)):
      
        resonance = resonances[i]
        
        if self.resonancesPeakContrib[i].has_key(resonance):
          peakDimContrib = self.resonancesPeakContrib[i][resonance]
          checkPeakDimContribs.append(peakDimContrib)
          
      if checkPeakDimContribs:
      
        checkPeakDimContribs.sort()
        
        #
        # TODO: Should really create and maintain a dict for this to keep speed up...
        #
        
        for peakContrib in self.peak.sortedPeakContribs():
        
          matchPeakDimContribs = list(peakContrib.peakDimContribs)
          matchPeakDimContribs.sort()
          
          if checkPeakDimContribs == matchPeakDimContribs:
            
            #
            # Already exists: can return without doing anything
            #
            
            return
    
      #
      # Create a peakContrib
      #
      
      peakContrib = self.peak.newPeakContrib(weight = 1)
      
      #
      # Get or create the peakDimContribs
      #
      
      for i in range(0,len(peakDims)):
      
        resonance = resonances[i]
        
        if not self.resonancesPeakContrib[i].has_key(resonance):
          
          peakDim = peakDims[i]
          
          peakDimContrib = peakDim.newPeakDimContrib(resonance = resonance, peakContribs = [peakContrib])
          
          self.resonancesPeakContrib[i][resonance] = peakDimContrib
          
        else:
        
          peakContrib.addPeakDimContrib(self.resonancesPeakContrib[i][resonance])

    else:
    
      #
      # Continue the loop
      #
  
      peakDim = self.peak.sortedPeakDims()[self.assignDim[peakAssignDim]]
      resName = self.assignCombination[peakAssignDim]

      for resonance in self.peakContribResonances[peakAssignDim][resName]:
      
        self.setPeakDimContribs(peakAssignDim + 1,
                                peakDims = peakDims + [peakDim],
                                resonances = resonances + [resonance])

    
  def createPeakAssignmentResonances(self):
  
    """
    
    Create resonances necessary for peakDimContribs and keep track of them and
    the resNames associated with them
    
    """
  
    #
    # Loop over all possible contributers from format assignment file
    #
    
    for resName in self.peakAssignment.assignments[self.peakAssignDim]:
      
      if not self.resonanceNames.has_key(resName):
    
        #
        # Create resonance if doesn't exist
        #
    
        resonance = self.createResonance(resName)
    
        #
        # Keep track of it...
        #

        resonances = [resonance]
        self.resonanceNames[resName] = [resonance]
    
      else:
    
        #
        # Resonance(s) exists, get info from self.resonanceNames
        #
    
        resonances = self.resonanceNames[resName]

      #
      # Keep track of the resonance(s)...
      #
                
      self.peakContribResonances[-1][resName] = resonances
      
      #
      # Check if a peakDimContrib already exists and keep track if so
      #
      
      for resonance in resonances:
      
        peakDimContrib = self.peakDim.findFirstPeakDimContrib(resonance = resonance)
      
        if peakDimContrib:

          self.resonancesPeakContrib[-1][resonance] = peakDimContrib


  def peakAssignmentCheckPeakList(self):
  
    if not self.peakList and self.guiParent:
      
      peakLists = []
      
      for exp in self.project.currentNmrProject.sortedExperiments():
        for ds in exp.sortedDataSources():
          if ds.numDim == self.peakAssignmentsFile.numDim:
            for pl in ds.sortedPeakLists():
              peakLists.append(pl)

      (selectionList,selectionDict) = createSelection(peakLists)

      interaction = self.multiDialog.SelectionList(

                         self.guiParent,
                         selectionList,
                         selectionDict = selectionDict,
                         title = "Project '%s': " % self.project.name + 'Select peak list to connect to assignment file',
                         text = 'Existing peak lists:',
                         urlFile = 'SelectPeakListForAssignFile'
                         
                         )
    
      if interaction.isSelected:
        self.peakList = interaction.selection
      
      else:
        return None

    else:
    
      #
      # Do check on numDim
      #
    
      if self.peakList.dataSource.numDim != self.peakAssignmentsFile.numDim:
        self.messageReporter.showError("Error","Number of dimensions for peak list and peak assignment file do not match.",self.guiParent)
        return None

    
  def peakAssignmentSetDimMapping(self):
    
    self.assignDim = []
    
    for origDim in self.peakAssignmentsFile.assignmentsDim:

      peakListDimFound = None

      for mapping in self.peakList.findAllApplicationData(application = self.format, keyword = 'dimMapping'):

        (origPeakFileDim,peakListDim) = returnInts(mapping.value.split(tagSep))

        if origPeakFileDim == origDim:
          peakListDimFound = peakListDim
          break
        
      self.assignDim.append(peakListDimFound)
    
  def peakAssignmentSetData(self):
  
    """
    
    Routine to set the ambiguous peak assignments from the peak assignment file in
    the data model
    
    """

    #
    # Once the linking is known, check whether peak number matches.
    # First make temporary dict with original peak number (if any - ignore if not available?!?)
    # Then get the right peak.
    #
    
    peakNumDict = {}
    
    for peak in self.peakList.sortedPeaks():
      
      applData = peak.findFirstApplicationData(application = self.format, keyword = peakNum_kw)
      
      if applData:
        peakNumDict[applData.value] = peak
        
    
    for self.peakAssignment in self.peakAssignmentsFile.peakAssignments:
      
      peakNum = self.peakAssignment.peakNum
      
      #
      # If there is a match, then see if there are already assignments attached
      # Only take the self.assignDims into account
      #

      if peakNumDict.has_key(peakNum):

        self.convertCount[self.mainCode][0] += 1

        self.peak = peakNumDict[peakNum]
        self.peakContribResonances = []    # TODO: could set these up from the start...
        self.resonancesPeakContrib = []
        
        #
        # Check if both relevant dims are empty: if not, then addMode has to be ON
        # before anything is added
        #
                        
        origAssignments = 0
        numPeakAssignDim = len(self.assignDim)
        
        if self.overwriteMode == 0:
        
          for self.peakAssignDim in range(numPeakAssignDim):
            
            if self.assignDim[self.peakAssignDim] != None:

              self.peakDim = self.peak.sortedPeakDims()[self.assignDim[self.peakAssignDim]]
            
              if self.peakDim.peakDimContribs:
  
                origAssignments = 1
        
        #
        # Loop over relevant peakDims
        #
        
        for self.peakAssignDim in range(numPeakAssignDim):
          
          self.peakContribResonances.append({})
          self.resonancesPeakContrib.append({})

          if self.assignDim[self.peakAssignDim] != None:

            self.peakDim = self.peak.sortedPeakDims()[self.assignDim[self.peakAssignDim]]
       
            #
            # If there are assignments present, delete them if overwrite mode is on
            #
          
            if self.peakDim.peakDimContribs:
  
              if self.overwriteMode == 1:
  
                #
                # Delete all exisiting ones
                #
                
                for peakDimContrib in self.peakDim.peakDimContribs:
                  peakDimContrib.delete()
                
                for peakContrib in self.peak.sortedPeakContribs():
                  peakContrib.delete()
            
            #
            # If no assignment present (or just deleted), or if assignments need to be added,
            # then add the new ones, create resonances if necessary
            #
            
            if (not self.peakDim.peakDimContribs and origAssignments == 0) or self.addMode == 1:
  
              self.createPeakAssignmentResonances()
                                
        #
        # If addMode is on, or if there were no original peakdimcontribs, set the peakContribs
        #
        
        if origAssignments == 0 or self.addMode == 1:
        
          for self.assignCombination in self.peakAssignment.assignCombinations:

            self.setPeakDimContribs(0,peakDims = [], resonances = [])
            
  def peakAssignmentsSetAssignmentInfo(self):
  
    pass
    
    
  def peakAssignmentsCreate(self):
  
    pass
        
  def getChainLinks(self):

    self.ccpChainDict = {}
    self.formatChains = {}
    self.formatChainCodeDict = {}
    self.chainMappingChain = {}
    
    self.setCoordChainList()

    for coordChain in self.coordChainList:
      
      chainCode = self.getChainOrRefChainId(coordChain,isChain = True)
      
      # Keep track of reported chainCodes for linking later on
      if not self.formatChains.has_key(chainCode):
        self.formatChains[chainCode] = []

    modelNums = self.coordinateFile.modelCoordinates.keys()
    modelNums.sort()
    
    badChainIdList = []
    
    if modelNums:
      refModel = modelNums[0]

      for coord in self.coordinateFile.modelCoordinates[refModel]:

        chainId = self.getChainOrRefChainId(coord)

        if chainId in self.formatChains:
          seqCodeKey = (coord.seqCode,coord.insertionCode)

          if seqCodeKey not in self.formatChains[chainId]:
            self.formatChains[chainId].append(seqCodeKey)

        elif chainId not in badChainIdList:
          badChainIdList.append(chainId)
    
    else:
      self.messageReporter.showError("Error","No models found in file - cannot continue.",parent = self.guiParent)

    #
    # Also sort seqCodes for resChains, and clean up if necessary
    #

    for formatChainCode in self.formatChains.keys():
      if not self.formatChains[formatChainCode]:
        del(self.formatChains[formatChainCode])
        continue
      self.formatChains[formatChainCode].sort()

    #
    # Link chains
    #
    
    #fck = self.formatChains.keys()
    #fck.sort()
    #for fc in fck:
    #  print fc, self.formatChains[fc]
    
    if badChainIdList:
      print "  Warning: no coordinate chain correspondence for chain(s) with codes %s!" % str(["'%s'" % chainCode for chainCode in badChainIdList])
      
    self.setChainLinks(forceChainMappings = self.forceChainMappings, resetMapping = self.resetMapping)

    #
    # Make sure it all went OK...
    #

    return self.status 

  def setCoordChainList(self):
  
    self.coordChainList = self.coordinateFile.chains
    
  def getChainOrRefChainId(self,coordOrChain,isChain = False):
      
    if hasattr(coordOrChain,'refChainId'):
      chainCode = coordOrChain.refChainId
    else:
      chainCode = coordOrChain.chainId
    
    return chainCode
  
  def autoSelectCoordNamingSystem(self):
    
    testMode = False
    
    if self.verbose:
      print "  Auto selecting naming system for coordinates..."
      
    modelNums = self.coordinateFile.modelCoordinates.keys()
    modelNums.sort()
    
    if not modelNums:
      raise self.FormatConverterError("No models available for determining naming system for coordinates.")
    
    modelNum = modelNums[0]
  
    namingSystemNames = {}
    atomsHandled = [] # Do each atom name only once
    atomSysNameDict = {}
    
    curCoordinates = self.coordinateFile.modelCoordinates[modelNum]

    # This is still necessary - have to map to the existing sequences!
    oldSeqCode = ""
    oldChainId = ""
    oldInsertionCode = ""

    for coordinate in curCoordinates:

      coordAtomName = coordinate.atomName   
      coordAtomId = "%s.%s" % (coordinate.resName,coordAtomName)
      
      if coordAtomId in atomsHandled:
        continue
      
      atomsHandled.append(coordAtomId)

      coordSeqCode = self.getCoordinateSeqCode(coordinate)
      refChainId = self.getChainOrRefChainId(coordinate)
      
      if testMode:
        print coordSeqCode, refChainId

      #
      # If change in chainId, reset variables
      #

      if refChainId != oldChainId:

        oldChainId = refChainId
        
        if self.formatChainCodeDict.has_key(oldChainId):
          chainMappings = self.formatChainCodeDict[oldChainId]
          if testMode:
            print chainMappings
            print self.formatChainCodeDict
        else:
          print "  Error: no chain mapping found for chain '%s'." % oldChainId
          oldChainId = ""
          continue

        #
        # If multiple chainMappings have to select on seqCode!!!!
        #

        (seqCodeToIdMapping,seqCodeToChainMapping,seqCodes) = self.getChainMappingsSeqCodeInfo(chainMappings)

        #
        # Try speed things up a bit...
        #

        if len(seqCodeToIdMapping) > 1:
          multipleChains = True

        else:
          multipleChains = False
          chainMapping = chainMappings[0]
          chain = self.chainMappingChain[chainMapping]
        
        if testMode:
          print "Multiple chains", multipleChains," chainMapping", chainMapping, chain
          print chainMapping.firstSeqId, chainMapping.formatFirstSeqCode
          print chainMapping.parent.formatName, chainMapping.parent.namingSystem
        
        # Make sure that code below gets called as well (could be two subsequent chains with seqCode 1!)
        oldSeqCode = ""

      #
      # If change in residue, reset variables
      #
      
      if coordSeqCode != oldSeqCode or coordinate.insertionCode != oldInsertionCode:

        oldSeqCode = coordSeqCode
        oldInsertionCode = coordinate.insertionCode

        if multipleChains:

          #
          # Have to find correct chainMapping and chain...
          #
          
          chain = None

          for (firstSeqCodeKey,lastSeqCodeKey) in seqCodes:
            (firstSeqCode,seqInsertCode) = firstSeqCodeKey
            (lastSeqCode,seqInsertCode) = lastSeqCodeKey
      
            if oldSeqCode >= firstSeqCode and oldInsertionCode == seqInsertCode and (lastSeqCode == None or oldSeqCode <= lastSeqCode):
              chainMapping = seqCodeToChainMapping[(firstSeqCode,seqInsertCode)]
              chain = self.chainMappingChain[chainMapping]
              break
          
          if not chain:
            oldSeqCode = ''      
            if testMode:
              print "Invalid chain identification based on seqCodes..."
              print seqCodes
              print oldSeqCode, oldInsertionCode
            continue
            
        else:
          
          #
          # Doublecheck
          #
          
          validSeqCode = False
          
          for (firstSeqCodeKey,lastSeqCodeKey) in seqCodes:
            (firstSeqCode,seqInsertCode) = firstSeqCodeKey
            (lastSeqCode,seqInsertCode) = lastSeqCodeKey
      
            if oldSeqCode >= firstSeqCode and oldInsertionCode == seqInsertCode and (lastSeqCode == None or oldSeqCode <= lastSeqCode):
              validSeqCode = True
          
          if not validSeqCode:
            oldSeqCode = ''
            if testMode:
              print "Invalid seq code identification based on seqCodes..."
              print seqCodes
              print oldSeqCode, oldInsertionCode
            continue
          
        seqId = oldSeqCode + chainMapping.firstSeqId - chainMapping.formatFirstSeqCode
        seqLabel = coordinate.resName

        residue = chain.findFirstResidue(seqId = seqId)

        if not residue:
          oldSeqCode = ''
          if testMode:
            print "Can't find residue %s" % residue
          continue

      #
      # Find in which naming systems this atom name exists
      #

      molResidue = residue.molResidue
      chemCompVar = residue.chemCompVar
      chemComp = molResidue.chemComp
      chemCompAtoms = chemComp.sortedChemAtoms()
      
      for namingSystem in chemComp.namingSystems:
        
        if not namingSystemNames.has_key(namingSystem.name):
          namingSystemNames[namingSystem.name] = 0
        
        sysNameKey = (chemComp.ccpCode,namingSystem.name)
        
        if atomSysNameDict.has_key(sysNameKey):
          chemAtomSysNames = atomSysNameDict[sysNameKey]
        else:
          chemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemComp,chemCompAtoms,namingSystem.name)
          atomSysNameDict[sysNameKey] = chemAtomSysNames
          
        atomNameFound = False
        
        for chemAtomSysName in chemAtomSysNames:
          
          if coordAtomName in (chemAtomSysName.sysName,) + chemAtomSysName.altSysNames:
            atomNameFound = True
            break
          
        if atomNameFound:
          if testMode:
            print "  FOUND",sysNameKey, namingSystem.name, coordAtomName
          namingSystemNames[namingSystem.name] += 1
  
    numAtoms = len(atomsHandled)
    
    if testMode:
      print numAtoms
      
      nsn = namingSystemNames.keys()
      nsn.sort()
      
      for ns in nsn:
        print "  %s: %d" % (ns,namingSystemNames[ns])
        
      asns = atomSysNameDict.keys()
      asns.sort()
      
      for asn in asns:
        print "  %s: %s" % (asn, atomSysNameDict[asn])
    
    #
    # Now sort out, and print out list (temporary code)
    #
    
    namingSystemNameList = namingSystemNames.keys()
    namingSystemNameList.sort()
    
    namingSystemHits = namingSystemNames.values()
    namingSystemHits.sort()
    namingSystemHits.reverse()
    
    bestHits = []

    if not namingSystemHits:
      highestMatch = 0
    else:
      highestMatch = namingSystemHits[0] * 100.0 / numAtoms
    
    for namingSystemHit in namingSystemHits:
      for namingSystemName in namingSystemNameList[:]:
        if namingSystemNames[namingSystemName] == namingSystemHit:
          bestHits.append(namingSystemName)
        
          if testMode:
            namingSystemNameList.pop(namingSystemNameList.index(namingSystemName))          
            print namingSystemName, namingSystemHit * 100.0  / numAtoms

      if not testMode:
        break
      
    #
    # Now, depending on situation, reset self.namingSystemName
    #
    # TODO: should probably allow selection here?
    
    if self.namingSystemName not in bestHits and bestHits:
      print "  Resetting naming system to %s for coordinates - has higher atom name match (%.f%%)" % (bestHits[0],highestMatch)
      self.namingSystemName = bestHits[0]
    else:
      print "  Selected naming system %s has %.f%% match with atom names." % (self.namingSystemName,highestMatch)
              
  def mapCoordinateAtoms(self):
    """ June 2011 Rasmus Fogh: Changed to fit new model:
    Map only first model - all models must have same atoms
    Added altLocation to mapping key
    """
  
    #
    # TODO: need to store mapping... also what happens to coordinates
    # if they can't be mapped?
    #
    
    if self.verbose:
      print "   Mapping coordinate atoms..."
  
    modelNums = self.coordinateFile.modelCoordinates.keys()
    modelNums.sort()
    modelNum = modelNums[0]
    
    self.coordAtomMapping = {}   
    self.setMapping = {}
    
    # Track mapping problems, and connect to structureEnsemble.
    mappingErrors = {'chain': [], 'residue': [], 'atom': []}
    
    self.ignoredSeqCodes = []
    
    chemAtomSysNamesByVar = {}
  
    curCoordinates = self.coordinateFile.modelCoordinates[modelNum]

    oldSeqCode = ""
    oldChainId = ""
    oldInsertionCode = ""

    for coordinate in curCoordinates:

      coordAtomName = coordinate.atomName   
      coordSeqCode = self.getCoordinateSeqCode(coordinate)

      refChainId = self.getChainOrRefChainId(coordinate)
      
      if hasattr(coordinate,'altLoc'):
        altLoc = coordinate.altLoc
      else:
        altLoc = defaultAltLoc
      
      
      mappingKey = (refChainId,coordSeqCode,coordinate.insertionCode,coordAtomName,altLoc)
      coordId = "%s.%d%s.%s.%s" % mappingKey 

      #
      # If change in chainId, reset variables
      #

      if refChainId != oldChainId:

        oldChainId = refChainId
        
        if self.formatChainCodeDict.has_key(oldChainId):
          chainMappings = self.formatChainCodeDict[oldChainId]
        else:
          mappingErrors['chain'].append(((oldChainId,),'No mapping'))
          print "  Error: no chain mapping found for chain '%s'." % oldChainId
          oldChainId = None       
          continue

        #
        # If multiple chainMappings have to select on seqCode!!!!
        #

        (seqCodeToIdMapping,seqCodeToChainMapping,seqCodes) = self.getChainMappingsSeqCodeInfo(chainMappings)

        #
        # Try speed things up a bit...
        #

        if len(seqCodeToIdMapping) > 1:
          multipleChains = True

        else:
          multipleChains = False
          chainMapping = chainMappings[0]
          chain = self.chainMappingChain[chainMapping]
        
        # Make sure that code below gets called as well (could be two subsequent chains with seqCode 1!)
        oldSeqCode = ""

      #
      # If change in residue, reset variables
      #
      
      if coordSeqCode != oldSeqCode or coordinate.insertionCode != oldInsertionCode:

        oldSeqCode = coordSeqCode
        oldInsertionCode = coordinate.insertionCode

        if multipleChains:

          #
          # Have to find correct chainMapping and chain...
          #
          
          chain = None
          
          #print oldSeqCode, oldInsertionCode

          for (firstSeqCodeKey,lastSeqCodeKey) in seqCodes:
            (firstSeqCode,seqInsertCode) = firstSeqCodeKey
            (lastSeqCode,seqInsertCode) = lastSeqCodeKey
            
            #print "  ",firstSeqCodeKey
            #print "  ",lastSeqCodeKey
            #print
      
            if oldSeqCode >= firstSeqCode and oldInsertionCode == seqInsertCode and (lastSeqCode == None or oldSeqCode <= lastSeqCode):
              chainMapping = seqCodeToChainMapping[(firstSeqCode,seqInsertCode)]
              chain = self.chainMappingChain[chainMapping]
              break
          
          if not chain:
            print "  Error: no chain/residue mapping for seqCode %d, chainCode '%s' - ignoring." % (oldSeqCode,oldChainId)
            self.ignoredSeqCodes.append((oldChainId,oldSeqCode,oldInsertionCode))
            mappingErrors['residue'].append(((oldChainId,oldSeqCode,oldInsertionCode),'No mapping'))
            oldSeqCode = ''
            
            continue
            
        else:
          
          #
          # Doublecheck
          #
          
          validSeqCode = False
          
          for (firstSeqCodeKey,lastSeqCodeKey) in seqCodes:
            (firstSeqCode,seqInsertCode) = firstSeqCodeKey
            (lastSeqCode,seqInsertCode) = lastSeqCodeKey
      
            if oldSeqCode >= firstSeqCode and oldInsertionCode == seqInsertCode and (lastSeqCode == None or oldSeqCode <= lastSeqCode):
              validSeqCode = True
          
          if not validSeqCode:
            print "  Error: seqCode %d not valid for chainCode '%s' - ignoring." % (oldSeqCode,oldChainId)
            self.ignoredSeqCodes.append((oldChainId,oldSeqCode,oldInsertionCode))
            mappingErrors['residue'].append(((oldChainId,oldSeqCode,oldInsertionCode),'Not a valid sequence code'))
            oldSeqCode = ''

            continue
          
        seqId = oldSeqCode + chainMapping.firstSeqId - chainMapping.formatFirstSeqCode
        seqLabel = coordinate.resName

        residue = chain.findFirstResidue(seqId = seqId)

        if not residue:
          if len(chain.residues) == 1:
            residue = chain.findFirstResidue()
            print "  Warning: single residue chain '%s', picking only residue for seqCode %d (seqId %d)" % (chain.code,oldSeqCode,residue.seqId)
          
          if not residue:

            print "  ERROR: can't handle seqCode %d in chain '%s' - no matching residue! Searched chain %s, seqId %d." % (oldSeqCode,oldChainId,chain.code,seqId)
            self.ignoredSeqCodes.append((oldChainId,oldSeqCode,oldInsertionCode))
            mappingErrors['residue'].append(((oldChainId,oldSeqCode,oldInsertionCode),'Cannot find residue match'))
            oldSeqCode = ''
            continue

        #
        # Get the relevant chemAtomSysNames
        #

        molResidue = residue.molResidue
        chemCompVar = residue.chemCompVar
        chemComp = molResidue.chemComp
        
        namingSystem = chemComp.findFirstNamingSystem(name = self.namingSystemName)

        #
        # If the name is not found for this naming system, try IUPAC as backup
        #
        
        if chemAtomSysNamesByVar.has_key(chemCompVar):
          (chemAtomSysNames,backupChemAtomSysNames) = chemAtomSysNamesByVar[chemCompVar]
        else:
          
          chemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemComp,chemCompVar.sortedChemAtoms(),self.namingSystemName)
          backupChemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemComp,chemCompVar.sortedChemAtoms(),'IUPAC')
          
        if not chemAtomSysNames:
          
          chemAtomSysNames = backupChemAtomSysNames
          backupChemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemComp,chemCompVar.sortedChemAtoms(),'MSD')
          
          if not chemAtomSysNames:
            chemAtomSysNames = backupChemAtomSysNames
            backupChemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemComp,chemCompVar.sortedChemAtoms(),'CUSTOM')
        
        if backupChemAtomSysNames:
          chemAtomSysNames = chemAtomSysNames + backupChemAtomSysNames
          
        if not chemAtomSysNamesByVar.has_key(chemCompVar):
          chemAtomSysNamesByVar[chemCompVar] = (chemAtomSysNames,backupChemAtomSysNames)

        # TODO; could in principle check if cifCode matches!!
        #checkValidResidueName(resid,coordinate.resName)

      #
      # Find the correct atom...
      #
      
      atomName = None
      atomSetName = None

      #
      # Check if the atom is already in the dict
      #

      if self.coordAtomMapping.has_key(mappingKey):

        continue

      #
      # If not then search for it...
      #

      else:

        atom = None

        #
        # See if any mapping was set already
        #

        if self.setMapping.has_key((oldChainId,seqLabel,coordAtomName)):

          atomName = self.setMapping[(oldChainId,seqLabel,coordAtomName)]
          atomStatus = 0

          if not atomName:

            #
            # Do not link this one...
            #

            continue

          else:

            atom = residue.findFirstAtom(name = atomName)  

        #
        # If not then look for sysNames...
        #

        else:

          atomStatus = 1

          searchNames = []

          if self.addNameMappings:
            for nameMappingKey in (residue.ccpCode,'ALL'):
              if self.addNameMappings.has_key(nameMappingKey):
                for (tOrigName,tNewName) in self.addNameMappings[nameMappingKey]:
                  if tOrigName == coordAtomName:
                    searchNames.append(tNewName)
                    break
                if searchNames:
                  break
          
          # Make sure that original atom name comes SECOND after mapped one - this is to avoid problems with
          # HB1/HB2 to HB2/HB3 type mappings for methylenes.
          searchNames.append(coordAtomName)

          #print searchNames

          nameMatchFound = 0

          for searchName in searchNames:
            for chemAtomSysName in chemAtomSysNames:

              #print chemAtomSysName.namingSystem, chemAtomSysName.sysName, chemAtomSysName.altSysNames

              if searchName in (chemAtomSysName.sysName,) + chemAtomSysName.altSysNames:

                nameMatchFound = 1

                atomName = chemAtomSysName.atomName

                #
                # Only break if actually exists...
                #

                atom = residue.findFirstAtom(name = atomName)

                if atom:
                  break

                # TODO THIS COMMENTED OUT - NOT USING SETS IN THIS ANY MORE. SHOULD I?
                #else:
                #  atomSetName = chemAtomSysName.atomName

            if nameMatchFound and atom:
              break
            
          if not nameMatchFound and self.matchToAtomNames:
            for searchName in searchNames:
              atom = residue.findFirstAtom(name = searchName)

        #
        # If there is a sysName match, look for the right atom
        #

        if atomName and not atom:

          #
          # Get the chemAtom...
          #

          chemAtom = chemCompVar.findFirstChemAtom(name = atomName)

          #
          # Get the relevant atom for the chemAtom match
          #

          if chemAtom:   

            atom = residue.findFirstAtom(chemAtom = chemAtom)

        if atom:

          atomStatus = 1
          self.coordAtomMapping[mappingKey] = [chain,residue,atom,atomStatus]

        elif not atomSetName:

          atomStatus = 0
          propagateMapping = None

          #
          # From here on proper mapping...
          #

          if self.linkAtoms:

            if not atomStatus:
              headerLines = []
              headerLines.append("Molecule: '%s'" % chain.molecule.name)
              headerLines.append("Chain code: '%s'" % chain.code)
              interaction = self.multiDialog.CoordinateAtomSelect(self.guiParent,residue,coordAtomName,seqLabel, headerLines = headerLines)

              chemAtomOrSets = interaction.chemAtomOrSets
              propagateMapping = interaction.propagate

          if propagateMapping != None:

            if chemAtomOrSets:

              #
              # Only single atom matches allowed
              #

              chemAtom = chemAtomOrSets[0]
              atom = residue.findFirstAtom(chemAtom = chemAtom)

            else:

              #
              # 'Do not link' case - can also be propagated!
              #

              chemAtom = None
              atom = None

            #
            # Set mapping
            #

            if atom:
              atomName = atom.name
              print "  New mapping: atom '%s' to atom %s" % (coordAtomName,atomName)
              self.coordAtomMapping[mappingKey] = [chain,residue,atom,atomStatus]

            else:
              atomName = None
              print "  New mapping: atom '%s' will not be linked." % (coordAtomName)

            #
            # Propagate for this residue only, for all residues of this type
            # in this chain or for all residues of this type in the file
            #


            if propagateMapping == 1:
              self.setMapping[(oldChainId,seqLabel,coordAtomName)] = atomName

            elif propagateMapping == 2:

              for tempChainId in self.formatChainCodeDict.keys():
                self.setMapping[(tempChainId,seqLabel,coordAtomName)] = atomName

          else:

            print "Error: Not linking atom %s (chemAtom selection dismissed)" % coordId
            self.addWarning("Coordinate atom %s not linked." % coordId)
            mappingErrors['atom'].append((mappingKey,'Cannot find atom match'))
            
            continue
            
    #
    # Now tag the structure ensemble with the mapping error info. Put in as dictionary!
    #
    
    hasMappingErrors = False
    for mappingErrorKey in mappingErrors.keys():
      if mappingErrors[mappingErrorKey]:
        hasMappingErrors = True
        break

    if hasMappingErrors:   
      self.structureEnsemble.addApplicationData(Implementation.AppDataString(application=ccpNmr_kw,keyword ='mappingErrors',value = str(mappingErrors)))
  
  def getCoordinateSeqCode(self,coordinate):
    
    return coordinate.seqCode
          
  def setCoordinates(self):

    modelNums = self.coordinateFile.modelCoordinates.keys()
    modelNums.sort()
    
    #
    # Set index mapping and create coord Atoms (etc.) if necessary
    #
    
    oldSeqCode = ""
    oldSeqInsertCode = ""
    oldChainId = ""
    oldChain = None
    
    indexMapping = {}
    nextIndex = 0
    orderedCoordAtoms = self.structureEnsemble.orderedAtoms
    
    curCoordinates = self.coordinateFile.modelCoordinates[modelNums[0]]
    for coordinate in curCoordinates:
      
      #
      # get mappingKey
      #
      chainId = self.getChainOrRefChainId(coordinate)
      seqCode = self.getCoordinateSeqCode(coordinate)
      seqInsertCode = coordinate.insertionCode
      atomName = coordinate.atomName
      
      if hasattr(coordinate,'altLoc'):
        altLoc = coordinate.altLoc
      else:
        altLoc = defaultAltLoc
      
      mappingKey = (chainId, seqCode, seqInsertCode, atomName, altLoc)
      
      #
      # set up mapping and create coord aAtoms as necessary
      #
      tt =  self.coordAtomMapping.get(mappingKey)
      if tt is not None:
      
        # get data for working
        (chain,residue,atom,atomStatus) = tt
        residueKey = (residue.seqCode,residue.seqInsertCode)

        #
        # Create coordinate chain if required
        #

        if chainId != oldChainId or oldChain != chain:
        
          oldChainId = chainId
          oldChain = chain
          oldSeqCode = ""
          
          cChain = self.structureEnsemble.findFirstCoordChain(code=chain.code)

          if not cChain:
            cChain = self.structureEnsemble.newChain(code=chain.code)
            self.setCoordChainOrigData(cChain,coordinate)
            self.coordResidues[cChain] = {}
          elif not self.coordResidues.has_key(cChain):
            self.coordResidues[cChain] = {}
            for cResidue in cChain.sortedResidues():
              residueKey = (cResidue.seqCode,cResidue.seqInsertCode)
              self.coordResidues[cChain][residueKey] = cResidue
            
        #
        # Create coordinate residue if required
        #

        if seqCode != oldSeqCode or seqInsertCode != oldSeqInsertCode:
        
          oldSeqCode = seqCode
          oldSeqInsertCode = seqInsertCode

          if self.coordResidues[cChain].has_key(residueKey):
            cResidue = self.coordResidues[cChain][residueKey]
          
          else:

            cResidue = cChain.newResidue(seqId = residue.seqId,seqCode = residueKey[0],seqInsertCode = residueKey[1])
            self.setCoordResidueOrigData(cResidue,coordinate,seqCode,seqInsertCode)
            
            self.coordResidues[cChain][residueKey] = cResidue
                       
        #
        # Create coordinate atom if required
        #
        cAtom = cResidue.findFirstAtom(name=atom.name, altLocationCode=altLoc)

        if cAtom is None:
          # No cAtom - make one
          
          if orderedCoordAtoms:
            raise self.FormatConverterError("Cannot create new Atoms in pre-existing Ensemble - %s"
                                            % (mappingKey,))
          
          # make new cAtom
          cAtom = cResidue.newAtom(name=atom.name, altLocationCode=altLoc)
          
          # set index mapping
          indexMapping[mappingKey] = nextIndex
          nextIndex += 1
          
          origAtomName = atomName
          tempAtomStatus = atomStatus
          if hasattr(coordinate,'origAtomName'):
            origAtomName = coordinate.origAtomName
            if origAtomName != atomName:
              tempAtomStatus = 0
          
          self.setCoordAtomOrigData(cAtom,coordinate,origAtomName,tempAtomStatus)
          
          if self.trackAtomSerials:
            setUniqueAppData('AppDataInt',atom,self.format,atomSerial_kw,coordinate.serial)
        
        elif not orderedCoordAtoms:
          #raise self.FormatConverterError(
          # "Duplicate atom record found while setting up new ensemble - %s"
          # % (mappingKey,)
          #)
          #print (
          # "WARNING Duplicate atom record found while setting up new ensemble - %s"
          # % (mappingKey,)
          #)
          pass
           
        else:
          # cAtom found. set mapping
          
          try:
            indexMapping[mappingKey] = orderedCoordAtoms.index(cAtom)
          except ValueError:
            raise self.FormatConverterError(
             "No mapping found for atom - %s"
             % (mappingKey,)
            )
            
    #
    # Create models and add their data
    #
    
    noMappingKeys = {}
    pseudoAtomMessage = 'Warning: Pseudoatom of Q type ignored'
    
    hasBFactors = False
    hasOccupancies = False
    nAtoms = len(indexMapping)
    iiii = 0
    skipcount = 0
    
    for modelNum in modelNums:
      
      indmap = indexMapping.copy()
    
      coordData = [0.0] * (nAtoms*3)
      occData = [0.0] * nAtoms
      bfacData = [0.0] * nAtoms
      
      #
      # Create a new model
      #
    
      self.structureNumber += 1
      print "  Reading model %d..." % self.structureNumber
            
      model = self.structureEnsemble.newModel()
      model.addApplicationData(Implementation.AppDataString(application=ccpNmr_kw, 
                                 keyword=originalFormat_kw, value=self.format))
      self.structures.append(model)

      # loop over coordinates to get data
      # (NB done to process in order and delete temporary data as we go)
      curCoordinates = self.coordinateFile.modelCoordinates.pop(modelNum)
      
      curCoordinates.reverse()
      while curCoordinates:
        coordinate = curCoordinates.pop()
        
        #
        # set mappingKey
        #
      
        chainId = self.getChainOrRefChainId(coordinate)
        seqCode = self.getCoordinateSeqCode(coordinate)
        seqInsertCode = coordinate.insertionCode
        atomName = coordinate.atomName
        
        if hasattr(coordinate,'altLoc'):
          altLoc = coordinate.altLoc
        else:
          altLoc = defaultAltLoc
        
        mappingKey = (chainId, seqCode, seqInsertCode, atomName, altLoc)
        
        # get mapping key
        tt =  self.coordAtomMapping.get(mappingKey)
        
        # get record index
        offset = indmap.get(mappingKey)
        if tt is None or offset is None:
          # set warning records and skip record
          
          skipcount += 1
          
          if (chainId,seqCode,seqInsertCode) in self.ignoredSeqCodes:
            continue
          
          # Special handling of pseudoatomstuff
          if atomName[0] == 'Q':
            mk = pseudoAtomMessage
          else:
            mk = mappingKey
          
          if not mk in noMappingKeys:
            noMappingKeys[mk] = []
          
          noMappingKeys[mk].append(modelNum)
          continue
        
        else:
          
          # delete mapping entry to avoid duplicates
          del indmap[mappingKey]
          # add data to data arrays
          self.convertCount[self.mainCode][1]['coordinates'] += 1
          iiii += 1
          
          tag = 'bFactor'
          if hasattr(coordinate, tag):
            hasBFactors = True
            bfacData[offset] = getattr(coordinate, tag)
          
          tag = 'occupancy'
          if hasattr(coordinate, tag):
            hasOccupancies = True
            occData[offset] = getattr(coordinate, tag)
          
          # coordinates
          ii = offset * 3
          for tag in ('x', 'y', 'z'):
            coordData[ii] = getattr(coordinate, tag)
            ii += 1
          
      # set data in model
      if hasBFactors:
        model.setSubmatrixData('bFactors', bfacData)
      if hasOccupancies:
        model.setSubmatrixData('occupancies', occData)
      model.setSubmatrixData('coordinates', coordData)
      
      self.convertCount[self.mainCode][1]['model(s)'] += 1

    #
    # Print out atoms that didn't work (and where)
    #
    
    if pseudoAtomMessage in noMappingKeys:
      ll = noMappingKeys.pop(pseudoAtomMessage)
      print "  %s in %d instances over all models." % (pseudoAtomMessage,len(ll))
    
    for key, modelNums in sorted(noMappingKeys.items()):
      print ("  Error: no mapping found for %s (models %s)" 
             % (key,','.join(str(x) for x in modelNums)))
    

  def setCoordChainOrigData(self,cChain,coordinate):
    
    chainId = coordinate.chainId
    cChain.addApplicationData(Implementation.AppDataString(application = self.format,keyword = 'originalChainCode',value = chainId))

  def setCoordResidueOrigData(self,cResidue,coordinate,seqCode,seqInsertCode):
    
    cResidue.addApplicationData(Implementation.AppDataInt(application = self.format,keyword = 'originalSeqCode',value = seqCode))
    cResidue.addApplicationData(Implementation.AppDataString(application = self.format,keyword = 'originalSeqInsertCode',value = seqInsertCode))

    if hasattr(coordinate,'resName'):
      cResidue.addApplicationData(Implementation.AppDataString(application = self.format,keyword = 'originalLabel',value = coordinate.resName))

  def setCoordAtomOrigData(self,cAtom,coordinate,atomName,atomStatus):

    if not atomStatus:

      #
      # TODO: here set original atomName as well! But ONLY if not recognized
      # by 'standard' naming system - for those just keep track of the naming
      # system and write them out that way (saves space at least)
      #

      cAtom.addApplicationData(Implementation.AppDataString(application = self.format,keyword = 'originalName',value = atomName))

  def getCoordinateChains(self):
    
    if self.exportChains and isinstance(self.exportChains, (list, tuple, set, frozenset)):
      
      self.chains = list(self.exportChains)
    
    else:
      
      self.chains = []
      
      for model in self.structures:
       
        for coordChain in model.structureEnsemble.coordChains:
        
          if coordChain.chain and coordChain.chain not in self.chains:
          
            self.chains.append(coordChain.chain)
    
  def createCoordinates(self):
    
    attrList1 = ['atom', 'seqCode', 'atomName', 'exportChainCode', 'atomType'] 
    attrList2 = ['origChainCode','origSeqCode','origSeqInsertCode',
                 'origResLabel','origAtomName']
                 
    if self.useOriginalData:
      attrListAll = attrList1 + attrList2
    else:
      attrListAll = attrList1
      
    
    # Sort by serial
    structures = sorted(self.structures, key=operator.attrgetter('serial'))
    
    # get all models into same ensemble. 
    # This ensures that atom names and ordering fit
    aset = set(x.structureEnsemble for x in structures)
    if len(aset) == 1:
      newEnsemble = None
    else:
      newEnsemble = makeEnsemble(structures, replaceStructures=False)
    
    try:
      # necessary to clean up new ensemble in case of error
    
      # set up data structures
      model = structures[0]
      structureEnsemble = model.structureEnsemble
      data = {}
      coordAtoms = structureEnsemble.orderedAtoms
      for tag in attrListAll:
        data[tag] = [None] * len(coordAtoms)
     
      # top level data
      origFormat = None
      if self.useOriginalData:
     
        origNamingSysName = structureEnsemble.atomNamingSystem
        appData = model.findFirstApplicationData(application=ccpNmr_kw, 
                                                 keyword=originalFormat_kw)
        if appData:
          origFormat = appData.value
      
      indices = []
      
      # loop over coordAtoms to get atom data
      for ii,coordAtom in enumerate(coordAtoms):
 
        # get atom, residue and chain
        atom = coordAtom.atom
        if not atom:
          print "   ERROR: no atom linked to coordinate atom %s - skipping" % coordAtom
          continue
        residue = atom.residue
        chain = residue.chain
        if chain not in self.chains:
          continue
        data['atom'][ii] = atom
        
        indices.append(ii)
        
        # get chain info
        data['exportChainCode'][ii] = self.chainDict[chain][0]
        
        # get residue info
        data['seqCode'][ii] = self.getExportSeqCode(self.chainDict[chain][1],
                                                    residue)
        chemCompVar = residue.chemCompVar
        chemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemCompVar.chemComp,
                                                          chemCompVar.chemAtoms,
                                                          self.namingSystemName)
        
        # get atom info
        atomName = None
        for chemAtomSysName in chemAtomSysNames:
          if chemAtomSysName.atomName == atom.name:
            atomName = chemAtomSysName.sysName
            break
        if not atomName:
          atomName = atom.name
          print (
           "  Warning: No chemAtomSysName for atom %s in chain %s, residue %s - using Iupac name" 
           % (atomName,chain.code,residue.seqCode)
          )
        data['atomName'][ii] = atomName
        
        data['atomType'][ii] = atom.chemAtom.elementSymbol
        
        # get original data
        if origFormat is not None:
          
          # chain
          origChainCode = coordChain.findFirstApplicationData(
                           application=origFormat, keyword='originalChainCode'
                          )
          if origChainCode:
            data['origChainCode'][ii] = origChainCode.value
  
          #residue
          origSeqCode = coordResidue.findFirstApplicationData(
                         application=origFormat, keyword='originalSeqCode'
                        )
          if origSeqCode:
            data['origSeqCode'][ii] = origSeqCode.value
            
          origSeqInsertCode = coordResidue.findFirstApplicationData(
                               application=origFormat, 
                               keyword='originalSeqInsertCode'
                              )
          if origSeqInsertCode:
            data['origSeqInsertCode'][ii] = origSeqInsertCode.value
            
          origResLabel = coordResidue.findFirstApplicationData(
                          application=origFormat, keyword='originalLabel'
                         )
          if origResLabel:
            data['origResLabel'][ii] = origResLabel.value
            
          origChemAtomSysNames = findAllSysNamesByChemAtomOrSet(
                                  chemCompVar.chemComp,chemCompVar.chemAtoms,
                                  origNamingSysName
                                 )
          
          #atom
          origAtomName = coordAtom.findFirstApplicationData(
                          application=origFormat, keyword='originalName'
                         )
                         
          if origAtomName:
            origAtomName = origAtomName.value
            
          else:
            for origChemAtomSysName in origChemAtomSysNames:
              if origChemAtomSysName.atomName == atom.name:
                origAtomName = origChemAtomSysName.sysName
                break
            if not origAtomName:
              origAtomName = atomName
          data['origAtomName'][ii] = origAtomName

    
      # do actual setting of coordinate data for all models.
      coordinateSerial = 0
      modelId = 0
      for model in structures:
 
        # model ID and set up
        modelId += 1
        coordinateSerial = self.resetCoordinateSerial(coordinateSerial)
        coordinates = model.coordinates
        bFactors = model.bFactors
        occupancies = model.occupancies
 
        # set input data and call coordinate setup
        for jj in indices:
          coordinateSerial += 1
          offset = jj*3
          atom = data['atom'][jj]
          residue = atom.residue
 
          params = {
          'modelId':modelId,
          'coordinateSerial':coordinateSerial,
          'x':coordinates[offset],
          'y':coordinates[offset+1],
          'z':coordinates[offset+2],
          'bFactor':bFactors[jj],
          'occupancy':occupancies[jj],
          'residue':residue,
          'chain':residue.chain
          }
 
          # set attributes on self for transfer to setRawCoordinate
          # Horrible, horrible, but too much work to redo
          for tag, val in params.items():
            setattr(self, tag, val)
          for tag in attrListAll:
            setattr(self, tag, data[tag][jj])
 
          self.setRawCoordinate()
 
    finally:
      # clean up
      if newEnsemble is not None:
        newEnsemble.delete()
      #raise
        

  def resetCoordinateSerial(self,coordinateSerial):
    
    # Can use this to reset the coordinate serial used in the output file on a per-model basis
    return coordinateSerial
         
  def setRawCoordinate(self):
  
    pass

  #
  # Chemcomp creation stuff...
  #

  def getFormula(self):
  
    atomTypes = {}
    
    for atom in self.rawChemComp.atoms:
    
      if hasattr(atom,'atomType') and atom.atomType:
      
        if not atomTypes.has_key(atom.atomType):
          atomTypes[atom.atomType] = 0

        atomTypes[atom.atomType] += 1
    
    atomTypeList = atomTypes.keys()
    atomTypeList.sort()
    
    self.formula = ''
    
    for atomType in atomTypeList:
      self.formula += "%s%d " % (atomType,atomTypes[atomType])
  
  def setChemCompInfo(self):
    
    defaultCcCode = 'enter_unique_code'
    
    if hasattr(self.rawChemComp,'isChiral'):
      isChiral = self.rawChemComp.isChiral
    else:
      isChiral = 0
      
    if self.molTypes:
      indexNum = self.rawChemCompList.index(self.rawChemComp)
      if indexNum < len(self.molTypes):
        molType = self.molTypes[indexNum]
    elif hasattr(self.rawChemComp,'molType'):
      molType = self.convertRawMolType(self.rawChemComp.molType)      
    else:
      molType = 'other'
    
    if self.ccpCodes:
      indexNum = self.rawChemCompList.index(self.rawChemComp)
      if indexNum < len(self.ccpCodes):
        ccCode = ccName = self.ccpCodes[indexNum]
    elif hasattr(self.rawChemComp,'name'):
      ccCode = ccName = self.rawChemComp.name
    else:
      ccCode = defaultCcCode
      ccName = 'chemComp name'
  
    self.chemCompInfo = [
    
      ('molType',molType),
      ('ccpCode',ccCode),
      ('name',ccName),
      ('merckCode',None),
      ('sigmaAldrichCode',None),
      ('beilsteinCode',None),
      ('casRegCode',None),
      ('details',None),
      ('hasStdChirality',isChiral)
      
    ]
    
    if not self.allowPopups and ccCode != defaultCcCode:
      valuesChecked = True
    else:
      valuesChecked = False
  
    while not valuesChecked:

      interaction = self.multiDialog.ChemCompInfo(self.guiParent,self.formula,len(self.rawChemComp.bonds),self.chemCompInfo)

      if not interaction.updated:
        break
        
      for (attrName,value) in self.chemCompInfo:
        
        if attrName == 'ccpCode':
          
          chemCompPathValid = True

          # Checking all chemComp repositories - would be problem if same ccpCode but different info therein!
          packageLocator = self.project.findFirstPackageLocator(targetName = 'ccp.molecule.ChemComp')
          repositoryList = packageLocator.repositories
 
          for repository in repositoryList:
            curChemCompPath = os.path.join(repository.url.path,'molecule','ChemComp',"%s.%s.xml" % (self.chemCompInfo[0][1],self.chemCompInfo[1][1]))
          
            if os.path.exists(curChemCompPath):
              self.messageReporter.showError("Error","Molecule type '%s', ccpCode '%s' already exists!" % (self.chemCompInfo[0][1],self.chemCompInfo[1][1]),parent = self.guiParent)
              chemCompPathValid = False
              break
          
          if chemCompPathValid:
            valuesChecked = 1
    
  def checkAtomAndBondInfo(self):
  
    atomSerials = {}
    
    self.chemCompAtoms = []
    self.chemCompBonds = []
    
    self.chemCompCoords = {}
    
    atomNames = {}
    
    for rawAtom in self.rawChemComp.atoms:
      
      if hasattr(rawAtom,'name') and rawAtom.name:
        atomName = rawAtom.name
      else:
        if not atomNames.has_key(rawAtom.atomType):
          atomNames[rawAtom.atomType] = 0
        atomNames[rawAtom.atomType] += 1
        atomName = "%s%d" % (rawAtom.atomType,atomNames[rawAtom.atomType])
      
      self.chemCompAtoms.append([
      
          ('name',          atomName),
          ('elementSymbol', rawAtom.atomType),
          ('chirality',     'unknown'),
          ('waterExchangeable',    False),
          ('bondedAtoms',    [])
          
          ])
          
      atomSerials[rawAtom.serial] = self.chemCompAtoms[-1]
      
      coords = []
      for coordId in ('x','y','z'):
        if hasattr(rawAtom,coordId):
          coords.append(getattr(rawAtom,coordId))
      
      if len(coords) == 3:
        self.chemCompCoords[atomName] = coords      
          
    for rawBond in self.rawChemComp.bonds:
      
      #
      # Use x,y,z to determine stereo?
      #
      
      atoms = (atomSerials[rawBond.atomSerial1],atomSerials[rawBond.atomSerial2])
      bondType = self.convertRawBondType(rawBond.bondType)
      
      self.chemCompBonds.append([
      
          ('bondType', bondType),     #   ('single', 'double', 'triple', 'aromatic', 'dative', 'singleplanar')
          ('stereochem', 'unknown'),  #     ('E', 'Z', 'unknown')
          ('atom1', atoms[0]),
          ('atom2', atoms[1])
          
          ])
      
      for i in range(0,2):
        atom = atoms[i]
        bondedAtoms = None
        for (attribute,value) in atom:
          if attribute == 'bondedAtoms':
            bondedAtoms = value
                
        bondedAtoms.append(atoms[not i])
            
    #
    # TODO: Redefine bonds based no new atom name information for protein?
    #        ('H1','N'): 'single',
    #    ('H2','N'): 'single',
    #    ('H','N'): 'single',
    #    ('N','CA'): 'single',
    #    ('CA','C'): 'single',
    #    ('C','O'): 'double',
    #    ('C',"O'"): 'double',
    #    ('C',"O''"): 'singleplanar',
    #    ("O''","H''"): 'single'
    #
    #
    
    #
    # Now user interaction...
    #
    #
    # TODO: for protein, DNA/RNA auto determine and name end sites automatically?
    #    
    
    if self.allowPopups:
      interaction = self.multiDialog.ChemCompAtomBondInfo(self.guiParent,self.chemCompInfo,self.chemCompAtoms,self.chemCompBonds)

    #
    # TODO: what about atomSets?!?!
    #

  def convertRawBondType(self,rawBondType):
    
    return rawBondType

  def convertRawMolType(self,rawMolType):
    
    return rawMolType
          
  def createLinkAndVars(self):
    
    molType = self.chemCompInfo[0][1] 
    
    linkEndInfo = {}
    self.molTypeKey = None
    
    #
    # See if any link end info available
    #
    
    for molTypeTuple in linkEndDict:
      if molType in molTypeTuple:
        linkEndInfo = linkEndDict[molTypeTuple]
        self.molTypeKey = molTypeTuple
        break
    
    linkings = ['none']
    if linkEndInfo:
      linkings.extend(linkingEnds.keys())
      self.setLinkEnds = True
    else:
      self.setLinkEnds = False
      
    #
    # Find out which atoms have to be removed
    #
    
    chemCompByLink = []
  
    for linking in linkings:
    
      atomsRemoved = []
      
      if linkingEnds.has_key(linking):
        linkEndList = linkingEnds[linking]
        
        if len(linkEndList) == 1:
          # Hack: have to switch to other linkEnd to get endAtom info
          linkEndList = [linkingEnds['middle'][not linkingEnds['middle'].index(linkEndList[0])]]
        
        for linkEnd in linkEndList:
          if self.molTypeKey:
            linkAtom = linkEndDict[self.molTypeKey][linkEnd][0]
            
            for atomInfo in self.chemCompAtoms:
              if atomInfo[0][1] == linkAtom:
                otherAtoms = atomInfo[-1][1]
                
                for otherAtom in otherAtoms:
                  if otherAtom[0][1] not in standardBackboneAtoms[self.molTypeKey]:
                    atomsRemoved.append(otherAtom[0][1])
    
      chemCompByLink.append((linking,atomsRemoved))
    
    #
    # Get a name for the descriptor
    #
    
    descriptor = 'neutral'
    
    if not self.minimalPrompts and self.allowPopups:
      descriptor = self.dataEntry.askString("Enter chemCompVar descriptor","Descriptor for the state of this chemCompVar",descriptor,self.guiParent)
    
    self.chemCompVars = [
    
      (descriptor,chemCompByLink)
    
    ]
    
    originalChemCompVar = (descriptor,chemCompByLink)
    
    #
    # Now user interaction...
    #
    
    if self.allowPopups:
      addDescriptor = 1
    else:
      addDescriptor = 0
    
    while (addDescriptor):
    
      descriptors = []
      for chemCompVarInfo in self.chemCompVars:
        descriptors.append(chemCompVarInfo[0])
    
      addDescriptor = self.messageReporter.showYesNo("Add/change descriptor","Do you want to add/change a chemCompVar? '%s' exist already." % string.join(descriptors,', '))

      if addDescriptor:
        
        getDescriptor = 1
        changeDescriptor = 0
        
        while (getDescriptor):
          descriptor = self.dataEntry.askString("Enter chemCompVar descriptor","Give a descriptor name:",descriptor,self.guiParent)
          
          if not descriptor:
            return
          elif descriptor not in descriptors:
            getDescriptor = 0
          else:
            if self.messageReporter.showYesNo("Error","Descriptor '%s' already exists! Do you want to overwrite?" % descriptor,parent = self.guiParent):
              changeDescriptor = 1
              getDescriptor = 0
        
        atomList = []
            
        for atomInfo in self.chemCompAtoms:
          atomList.append(atomInfo[0][1])

        interaction = self.multiDialog.MultiSelectionList(

                                      self.guiParent,
                                      atomList,
                                      title = "Project '%s': " % self.project.name + "Atoms for chemCompVar '%s'" % descriptor,
                                      text = 'Select atoms to remove for this chemCompVar:',
                                      urlFile = ''

                                      )

        if interaction.isSelectedList:
          atomsRemoved = interaction.isSelectedList[:]
        
        #
        # Set it all up...
        #
        
        chemCompByLink = []
        
        for (linking,linkAtomsRemoved) in originalChemCompVar[1]:
          newAtomRemovedList = linkAtomsRemoved[:]
          for atomName in atomsRemoved:
            if atomName not in newAtomRemovedList:
              newAtomRemovedList.append(atomName)
          chemCompByLink.append((linking,newAtomRemovedList))        
        
        if changeDescriptor:
          for i in range(0,len(self.chemCompVars)):
            if descriptor == self.chemCompVars[i][0]:
              self.chemCompVars[i] = (descriptor,chemCompByLink)
        else:
          self.chemCompVars.append((descriptor,chemCompByLink))        
    
    #
    # TODO TODO HAVE TO HANDLE LINKENDS HERE AS WELL!!
    #
    # Need a window to handle descriptors, and add linking for one particular
    # descriptor or all descriptors.
    #
    
  def setupChemCompInfo(self):

    self.getFormula()
    
    self.setChemCompInfo()
    
    self.checkAtomAndBondInfo()
    
    self.createLinkAndVars()
     
  def makeChemComp(self):
    
    #
    # Make the chemComp
    #

    createDict = {}
    for i in range(0,len(self.chemCompInfo)):
      createDict[self.chemCompInfo[i][0]] = self.chemCompInfo[i][1]
  
    chemComp = self.project.newNonStdChemComp(**createDict)     
    self.convertCount[self.mainCode][0] += 1
    
    #
    # Also make a chemCompCoord, if necessary
    #
    
    if self.chemCompCoords:
      keywds = {'sourceName': coordinateImportSourceName,
                'molType':    chemComp.molType,
                'ccpCode':    chemComp.ccpCode,
                'details':    'Generated through FormatConverter import.'}

      chemCompCoord = self.project.newChemCompCoord(**keywds)
      
      chemCompVarToCoord = {}
    
    #
    # Make a naming system, if required
    #
    
    if self.makeNamingSystem:
      namingSystem = chemComp.newNamingSystem(name = self.makeNamingSystem)
    
    #
    # Set the atoms
    #
    
    chemAtoms = []  # Same order as chemCompAtoms!
    
    for atomInfo in self.chemCompAtoms:
      
      createDict = {}
      for i in range(0,len(atomInfo)):
        (keyword,value) = atomInfo[i]
        if keyword not in ('bondedAtoms',):
          createDict[keyword] = value
      
      chemAtom = ChemComp.ChemAtom(chemComp,**createDict)
      chemAtoms.append(chemAtom)
      self.convertCount[self.mainCode][1]['atoms'] += 1
      
      if self.makeNamingSystem:
        namingSystem.newAtomSysName(sysName = chemAtom.name, atomName = chemAtom.name, atomSubType = chemAtom.subType)
    
    #
    # Set the normal bonds and dihedrals
    #
    
    for bondInfo in self.chemCompBonds:
      
      createDict = {}
      for i in range(0,2):
        createDict[bondInfo[i][0]] = bondInfo[i][1]
        
      createDict['chemAtoms'] = []
      for i in range(2,4):
        atomIndex = self.chemCompAtoms.index(bondInfo[i][1])
        chemAtom = chemAtoms[atomIndex]
        createDict['chemAtoms'].append(chemAtom)

      chemBond = ChemComp.ChemBond(chemComp,**createDict)
      self.convertCount[self.mainCode][1]['bonds'] += 1
    
    
    #
    # If there are linkCodes, set them...
    #
    
    if self.setLinkEnds:    
    
      #
      # Set linkCodes, linkAtoms, linkBonds, ends, ...
      #

      for (linkCode,posRange) in linkCodes:
        createLinkAtom = True
        linkEndChemAtomNames = linkEndDict[self.molTypeKey][linkCode]
        for chemAtomName in linkEndChemAtomNames:
          chemAtom = chemComp.findFirstChemAtom(name = chemAtomName)
          if not chemAtom:
            createLinkAtom = False
        
        if createLinkAtom:   
          for pos in range(1,posRange+1):
            name = "%s_%d" % (linkCode,pos)
            linkAtom = ChemComp.LinkAtom(chemComp,name = name)

      #
      # Set the torsions and bonds
      #

      for (torsionName,torsionAtomNames) in linkTorsions[self.molTypeKey]:
        chemAtomList = []

        for chemAtomName in torsionAtomNames:
          chemAtomList.append(chemComp.findFirstChemAtom(name = chemAtomName))

        if None not in chemAtomList and not chemComp.findFirstChemTorsion(chemAtoms = tuple(chemAtomList)):
          linkTorsion = ChemComp.ChemTorsion(chemComp, name = torsionName, chemAtoms = chemAtomList)

      for (bondType,stereoChem,bondAtomNames) in linkBonds[self.molTypeKey]:
        chemAtomList = []

        for chemAtomName in bondAtomNames:
          chemAtomList.append(chemComp.findFirstChemAtom(name = chemAtomName))

        if None not in chemAtomList and not chemComp.findFirstChemBond(chemAtoms = tuple(chemAtomList)):
          linkBond = ChemComp.ChemBond(chemComp, bondType = bondType, chemAtoms = chemAtomList, stereochem = stereoChem)

      #
      # Set the linkEnds
      #

      for (linkCode,posRange) in linkCodes:
        linkEndChemAtomNames = linkEndDict[self.molTypeKey][linkCode]
        linkEndChemAtoms = []
        linkEndLinkAtoms = []
        missingAtoms = []

        for chemAtomName in linkEndChemAtomNames:
          chemAtom = chemComp.findFirstChemAtom(name = chemAtomName)
          linkEndChemAtoms.append(chemAtom)
          if not chemAtom:
            missingAtoms.append(chemAtomName)

        for pos in range(1,posRange+1):
          name = "%s_%d" % (linkCode,pos)
          chemAtom = chemComp.findFirstChemAtom(name = name)
          linkEndLinkAtoms.append(chemAtom)
          if not chemAtom:
            missingAtoms.append(chemAtomName)

        if None not in linkEndChemAtoms + linkEndLinkAtoms:
          
          if not chemComp.findFirstLinkEnd(linkCode = linkCode):
            
            if len(linkEndLinkAtoms) == 2:
              remoteLinkAtom = linkEndLinkAtoms[1]
            else:
              remoteLinkAtom = None
  
            linkEnd = ChemComp.LinkEnd(chemComp, linkCode = linkCode,
                                                 boundChemAtom = linkEndChemAtoms[0],
                                                 remoteChemAtom = linkEndChemAtoms[1],
                                                 boundLinkAtom = linkEndLinkAtoms[0],
                                                 remoteLinkAtom = remoteLinkAtom)
  
            print "  Made linkEnd %s" % linkCode
          
        else:
          print "  Could not make linkEnd %s - %s atoms missing." % (linkCode,str(missingAtoms))

    #
    # Create the vars...
    #
    
    for (descriptor,chemCompByLink) in self.chemCompVars:
      for (linking,atomRemovedList) in chemCompByLink:
        
        isValid = True
        chemAtomList = []

        if linking in linkingEnds.keys():
        
          for linkCode in linkingEnds[linking]:
            linkEnd = chemComp.findFirstLinkEnd(linkCode = linkCode)
            if linkEnd:
              for linkAtom in [linkEnd.boundLinkAtom,linkEnd.remoteLinkAtom]:
                if linkAtom:
                  chemAtomList.append(linkAtom)
            else:
              print "  Error: no link code %s, %s" % (linkCode,descriptor)
              isValid = False
      
        for i in range(len(self.chemCompAtoms)):
          atomName = self.chemCompAtoms[i][0][1]
      
          if atomName not in atomRemovedList:
            chemAtom = chemAtoms[i]
            chemAtomList.append(chemAtom)

        if isValid:
          #TODO: have to set charge, ... stuff on HIGHER level!!
          chemCompVar = ChemComp.ChemCompVar(chemComp, linking = linking,formalCharge = 0, isAromatic = 0, isParamagnetic = 0, descriptor = descriptor, chemAtoms = chemAtomList)
          print "Created %s,%s: %s,%s" % (chemComp.ccpCode,chemComp.molType,descriptor,linking)
        
          if self.chemCompCoords:
            chemCompVarCoord = chemCompCoord.newChemCompVarCoord(linking = chemCompVar.linking, descriptor = chemCompVar.descriptor)
            chemCompVarToCoord[chemCompVar] = chemCompVarCoord

        else:
          print "Could not create %s,%s: %s,%s. Link info missing" % (chemComp.ccpCode,chemComp.molType,descriptor,linking)
          
    #
    # Now create the chemAtomCoords, if necessary
    #
    
    if self.chemCompCoords:
      
      atomNames = self.chemCompCoords.keys()
      for atomName in atomNames:
        chemAtom = chemComp.findFirstChemAtom(name=atomName)
        chemCompVarCoords = []
        for chemCompVar in chemAtom.chemCompVars:
          chemCompVarCoords.append(chemCompVarToCoord[chemCompVar])
        
        (x,y,z) = self.chemCompCoords[atomName]
        
        chemCompCoord.newChemAtomCoord(name=chemAtom.name,subType=chemAtom.subType,chemCompVarCoords=chemCompVarCoords,x=x,y=y,z=z)
    
    #
    # Create atom sets...
    #
    
    if self.createAtomSets:
     
      makeAtomSets = MakeAtomSets()
      
      kywds = {'createSysNames': False}
      if self.makeNamingSystem:
        kywds['createSysNames'] = True
        kywds['namingSystem'] = namingSystem

      makeAtomSets.setChemComps(chemComp, **kywds)
    
    #
    # Save chemComp...
    #
    
    if self.saveChemComp:
      
      repository = self.getChemCompRepository()

      chemComp.saveTo(repository)
      
      if self.chemCompCoords:
        chemCompCoord.saveTo(repository)
    
    self.chemComps.append(chemComp)

  #
  # Writing chemComp to external format...
  #
  
  def getChemCompCoord(self):
    
    self.chemCompVarCoord = None
    
    chemComp = self.chemCompVar.chemComp
    
    localKeywds = {}
    if self.noWrite or not self.copyCoordinates:
      localKeywds['copyFile'] = False

    ccc = getChemCompCoord(self.project, self.coordSystem, chemComp.molType, chemComp.ccpCode, **localKeywds)
      
    if not ccc:
      
      # TODO: need something similar to chemComp XML file copying here?!??!
      pass

    if ccc:
      
      self.chemCompVarCoord = ccc.findFirstChemCompVarCoord(linking = self.chemCompVar.linking, descriptor = self.chemCompVar.descriptor)
  
  
    if not self.chemCompVarCoord:
      print "  Warning: could not find chemCompVarCoord for %s, %s, %s, %s, %s!" % (self.coordSystem,chemComp.molType,chemComp.ccpCode,self.chemCompVar.linking,self.chemCompVar.descriptor)

  def createRawChemComp(self):
    
    self.setRawChemComp()
      
    coords = {}
    linkAtomList = []
    rawAtomSerialDict = {}
    ccpCode = self.chemCompVar.chemComp.ccpCode
    molType = self.chemCompVar.chemComp.molType
    
    namingSystem = self.chemCompVar.chemComp.findFirstNamingSystem(name = self.namingSystemName)
    
    atomSerial = 1
              
    for chemAtom in self.chemCompVar.sortedChemAtoms():
    
      if chemAtom.className == 'LinkAtom':
        linkAtomList.append(chemAtom.name)
        continue

      atomSysName = findChemAtomOrSetSysName(namingSystem, chemAtom)
      
      if not atomSysName:
        print "  Warning: no %s chemAtomSysName for atom %s - using atom name instead." % (self.namingSystemName,chemAtom.name)
        rawAtomName = chemAtom.name
      else:
        rawAtomName = atomSysName.sysName
        
      
      chemAtomCoord = self.chemCompVarCoord.findFirstChemAtomCoord(name = chemAtom.name)
      
      if not chemAtomCoord:
        print " Warning: missing coord for atom %s - set to 0.0" % (chemAtom.name)
        (x,y,z) = (0.0,0.0,0.0)
      else:
        (x,y,z) = (chemAtomCoord.x,chemAtomCoord.y,chemAtomCoord.z)
     
      # TODO; could be more sophisticated... need to fill in nuclGroupType for that though? Or get from Vega type?
      rawAtomType = chemAtom.elementSymbol
      
      # TODO: when available!!
      charge = 0.0
      
      self.setRawChemCompAtom(atomSerial,rawAtomName,x,y,z,rawAtomType,charge,ccpCode)
     
      rawAtomSerialDict[chemAtom.name] = atomSerial
      coords[chemAtom.name] = (x,y,z)

      atomSerial += 1   
    
    #
    # Now loop over the bonds to prep info (have to be able to add to it if link-type atoms added...)
    #
    
    chemBondInfo = []
    
    for chemBond in self.chemCompVar.chemBonds:
      
      chemAtoms = chemBond.sortedChemAtoms()
      chemBondInfo.append(((chemAtoms[0].name,chemAtoms[1].name),chemBond.bondType))
          
    if self.writeLinkAtoms and linkAtomList:

      addAtomData = []
    
      if molType == 'protein':
        
        # TODO this could also be part of reference data!!!
        if 'prev_1' in linkAtomList:
          addAtomData.append(('C_P',(('CA','N','H'),1.32,123.0,'singleplanar')))
          addAtomData.append(('CAP',(('N','C_P','H'),1.53,114.0,'single')))
          addAtomData.append(('O_P',(('N','C_P','CAP'),1.24,125.0,'double')))
        if 'next_1' in linkAtomList:
          addAtomData.append(('N_N',(('CA','C','O'),1.32,114.0,'singleplanar')))
          addAtomData.append(('H_N',(('C','N_N','O'),0.97,123.0,'single')))
          addAtomData.append(('CAN',(('C','N_N','H_N'),1.47,123.0,'single')))
      
      for (linkAtomName,atomData) in addAtomData:
         
        # Also make CH3 out of CA... not ideal! Use 1.080 for C-H bond length, eclipse
        # with respect to C=O bond, same plane as NH (for first one).
        # Other two calculate how? Are out-of-plane by certain amount... need improper type stuff for this
        (x,y,z) = calcPlanarCoord(coords[atomData[0][0]],
                                  coords[atomData[0][1]],
                                  coords[atomData[0][2]],
                                  atomData[1],
                                  convertDegreeRadian(atomData[2]))

        self.setRawChemCompAtom(atomSerial,linkAtomName,x,y,z,linkAtomName[0],0.0,ccpCode)
       
        chemBondInfo.append(((linkAtomName,atomData[0][1]),atomData[3]))
        rawAtomSerialDict[linkAtomName] = atomSerial
        coords[linkAtomName] = (x,y,z)

        print "  Created %s for output." % (linkAtomName)
        
        atomSerial += 1

    #
    # Now set the bonds for real
    #
    
    for (atomNames,bondType) in chemBondInfo:
      
      chemAtomSerials = []
      
      for atomName in atomNames:
        if rawAtomSerialDict.has_key(atomName):
          chemAtomSerials.append(rawAtomSerialDict[atomName])
      
      if len(chemAtomSerials) == 2:
        
        self.setRawChemCompBond(chemAtomSerials[0],chemAtomSerials[1],bondType)

  def setRawChemComp(self):
    
    pass
  
  def setRawChemCompAtom(self,atomSerial,rawAtomName,x,y,z,rawAtomType,charge,chemCompCode):
  
    pass

  def setRawChemCompBond(self,atomSerial1,atomSerial2,bondType):
  
    pass
 
  def makeInitialsForCcpn(self, initials):

    # TODO change this to readable code:
    #   pseudo unix code (x->y): set x = 'A, Z,_@R'; set y = `echo $x | tr 'a-z' 'A-Z' | tr -dc 'A-Z'`; echo $y

    # Was simply:
    #initials = initials.replace(' ','')

    if isinstance(initials, str):
      initials = ','.join(initials)
    return initials.translate(string.maketrans(string.punctuation, ' '*len(string.punctuation))).replace(' ','').upper()

  def setPeople(self):
    
    self.rawPersonToCcpnPerson = {}
    
    for rawPerson in self.peopleAndCitationsFile.persons:
    
      firstName = rawPerson.firstName
      initials = rawPerson.initials

      familyName = rawPerson.familyName

      if rawPerson.familyTitle:
        if rawPerson.familyTitle[:2].upper() in ('JR', 'SR'):
          familyTitle = rawPerson.familyTitle.capitalize()
        else: # TO avoid things like 'Iii'
          familyTitle = rawPerson.familyTitle.upper()
      else:
        familyTitle = None
      
      if rawPerson.title:
        title = rawPerson.title.capitalize()
      else:
        title = None

      if familyName and familyName.upper() == familyName:
        familyName = familyName.capitalize()
      if firstName and firstName.upper() == firstName and not firstName.count('.'):
        firstName = firstName.capitalize()
      
      person = None
      if self.useExistingPersons and self.project.currentAffiliationStore:
        person = self.findVaguePerson(familyName,firstName,initials,familyTitle)

        if person and not person.middleInitials and initials and isinstance(initials, str):
          person.middleInitials = self.makeInitialsForCcpn(initials)

          #print 'FOUND: [%s]' % (person, person.familyName)

      # Make sure things don't get duplicated if similar people available!
      if person in self.people:
        person = None

      if not person:

        keywds = {'familyName': familyName}
        keywds['givenName'] = firstName
        keywds['middleInitials'] = []
        keywds['familyTitle'] = None
        keywds['title'] = None

        if initials:
          if isinstance(initials, str):
            initials = self.makeInitialsForCcpn(initials)

          keywds['middleInitials'] = initials

        if familyTitle:
          keywds['familyTitle'] = familyTitle

        if title:
          if title[-1] == '.':
            title = title[:-1]
          keywds['title'] = title


        setCurrentStore(self.project,'AffiliationStore')
        person = Affiliation.Person(self.project.currentAffiliationStore,**keywds)

        #print 'NEW PERSON: [%s] [%s]' % (person, person.familyName)

        # For setting NmrStar contact information we need to keep track of the person's serial number.

        if hasattr(rawPerson, 'serial') and rawPerson.serial is not None:
          keywds2 = {'application': self.format,
                     'keyword':     'Serial',
                     'value':       rawPerson.serial}

          appData = Implementation.AppDataInt(**keywds2)
          person.addApplicationData(appData)

      self.people.append(person)
      
      if not self.rawPersonToCcpnPerson.has_key(rawPerson):
        self.rawPersonToCcpnPerson[rawPerson] = []
      self.rawPersonToCcpnPerson[rawPerson].append(person)
      
  def findVaguePerson(self,familyName,firstName,initials,familyTitle):

    personsFound = []

    for tperson in self.project.currentAffiliationStore.sortedPersons():
      if tperson.familyName and familyName and tperson.familyName.upper() == familyName.upper():
        if (tperson.givenName and firstName and tperson.givenName.upper() == firstName.upper()) or \
           (firstName and ((tperson.firstInitial and tperson.firstInitial.upper() == firstName[0].upper()) or \
           (tperson.givenName and tperson.givenName[0].upper() == firstName[0].upper()))):

          personsFound.append(tperson)
    
    person = None
    if personsFound:
      for tperson in personsFound:
        if initials is None:
          initials = ''

        if tperson.middleInitials == tuple(initials) and ( (not familyTitle) or tperson.familyTitle == familyTitle):
          person = tperson
          break
      
      if not person:
        if len(personsFound) > 1:
          print "  Warning: match (but no full match) found for person %s %s %s - picking first match." % (firstName,initials,familyName)
        person = personsFound[0]
        
    return person
    
  def setCitations(self):

    for rawCitation in self.peopleAndCitationsFile.citations:    

      title     = rawCitation.title
      status    = rawCitation.status
      #firstPage = rawCitation.firstPage
      #lastPage  = rawCitation.lastPage
      year      = rawCitation.year
      #volume    = rawCitation.volume

      if rawCitation.type == 'journal':
        className = 'JournalCitation'
      elif rawCitation.type == 'thesis':
        className = 'ThesisCitation'
      elif rawCitation.type == 'book':
        className = 'BookCitation'
      elif rawCitation.type in ('conference', 'abstract'):
        className = 'ConferenceCitation'
      else:
        className = 'JournalCitation'

      citation = None

      if self.useExistingCitations and self.project.currentCitationStore:
        citation = self.findVagueCitation(className, title, year)

      # Make sure things don't get duplicated if similar citations available!
      if citation in self.citations:
        citation = None

      if not citation:
        keywds = {}

        if rawCitation.type == 'journal':
          className = 'JournalCitation'

          for (ccpnKeywd,attrName) in (('journalAbbreviation','pubShortName'),('journalFullName','pubLongName'),('issue','issue')):
            attrValue = getattr(rawCitation,attrName)
            if attrValue != None:
              keywds[ccpnKeywd] = attrValue

        else:

          if rawCitation.type == 'thesis':
            className = 'ThesisCitation'

            keywds['title'] = rawCitation.title

          elif rawCitation.type == 'book':
            className = 'BookCitation'

            keywds['title'] = rawCitation.title

          elif rawCitation.type in ('conference', 'abstract'):
            className = 'ConferenceCitation'

          if rawCitation.type == 'book':
            keywds['publisher'] = rawCitation.publisher
            keywds['publisherCity'] = rawCitation.publisherPlace
            keywds['bookTitle'] = rawCitation.bookTitle
            keywds['bookSeries'] = rawCitation.bookSeries

          elif rawCitation.type == 'thesis':
            keywds['institution'] = rawCitation.publisher
            keywds['city'] = rawCitation.publisherPlace            
            keywds['country'] = rawCitation.country

          elif rawCitation.type in ('abstract', 'conference'):
            keywds['city'] = rawCitation.publisherPlace
            keywds['conferenceTitle'] = rawCitation.bookTitle
            keywds['country'] = rawCitation.country

          if rawCitation.editors:
            keywds['editors'] = []
            for rawPerson in rawCitation.editors:

              # Have to avoid duplication of people with same name... is possible (see PDB entry 1fO3)
              for person in self.rawPersonToCcpnPerson[rawPerson]:
                if not person in keywds['editors']:
                  keywds['editors'].append(person)
                  if hasattr(rawPerson, 'serial') and rawPerson.serial is not None:

                    keywds2 = {'application': self.format,
                               'keyword':     'serialEdit',
                               'value':       rawPerson.serial}

                    appData = Implementation.AppDataInt(**keywds2)
                    person.addApplicationData(appData)

                  break

        #
        # Set authors
        #

        if rawCitation.authors:
          keywds['authors'] = []
          for rawPerson in rawCitation.authors:
            for person in self.rawPersonToCcpnPerson[rawPerson]:
              if not person in keywds['authors']:
                keywds['authors'].append(person)

                if hasattr(rawPerson, 'serial') and rawPerson.serial is not None:

                  keywds2 = {'application': self.format,
                             'keyword':     'serialCit',
                             'value':       rawPerson.serial}

                  appData = Implementation.AppDataInt(**keywds2)
                  person.addApplicationData(appData)

                break

        if rawCitation.keywords:
          keywds['keywords'] = []

          for keyword in rawCitation.keywords:
            if not keyword in keywds['keywords']:
              keywds['keywords'].append(keyword.keywd)

        #
        # Set title
        #

        # TODO remove returnMemopsLine() when data type changed!!

        if title:
          title = rawCitation.title.replace('\n', ' ')
          title = title.strip()

        keywds['title'] = title

        #
        # Set whether primary ref
        # TODO THIS IS PART OF ENTRY NOW!!! Not set on citation level 

        #keywds['isPrimary'] = rawCitation.isPrimary

        #
        # Published?
        #

        if rawCitation.status == 'published':
          keywds['status'] = 'published'

          if rawCitation.firstPage:
            keywds['firstPage'] = rawCitation.firstPage

          if rawCitation.year:
            keywds['year'] = returnInt(rawCitation.year)

          if rawCitation.volume and className not in('ThesisCitation','ConferenceCitation'):
            keywds['volume'] = rawCitation.volume

          if rawCitation.lastPage:
            keywds['lastPage'] = rawCitation.lastPage

          #if className == 'BookCitation':
          #  keywds['bookSeries'] = curInfo['attrs']['seriesInfo']

        elif rawCitation.status in ('submitted', 'accepted', 'in preparation', 'to be published', 'in press'):
          keywds['status'] = rawCitation.status

        else:
          keywds['status'] = 'to be published'

        if rawCitation.pubMedId:
          keywds['pubMedId'] = rawCitation.pubMedId

        if rawCitation.doi:
          keywds['doi'] = rawCitation.doi

        if rawCitation.medline:
          keywds['medlineUiCode'] = rawCitation.medline

        if rawCitation.details:
          keywds['details'] = rawCitation.details

        #
        # Search for reference, if not found create
        #

        setCurrentStore(self.project,'CitationStore')

        citations = self.project.currentCitationStore.findAllCitations(className = className)

        citation = None

        for tcit in citations:
          allMatch = True
          for keywd in keywds:
            # TODO: How do I handle these correctly? Need to sort...
            if keywd in ('authors','editors'):
              if list(getattr(tcit,keywd)).sort() != keywds[keywd].sort():
                allMatch = False
            if getattr(tcit,keywd) != keywds[keywd]:
              allMatch = False
              break

          if allMatch:
            citation = tcit
            break

        if not citation:
          citation = getattr(Citation,className)(self.project.currentCitationStore,**keywds)

      self.citations.append(citation)
          
      # TODO set this?
      #self.citation.setMolSystems([self.molSystem])

  def setKeywords(self):

    for keyword in self.peopleAndCitationsFile.keywords:
      keyword = rawCitation.keywd


  def findVagueCitation(self, className, title, year):

    citationsFound = []

    citations = self.project.currentCitationStore.findAllCitations()

    for tcitation in citations:
      if tcitation.className == className and (tcitation.title and title and tcitation.title.upper() == title.upper() ):
        #and (tcitation.status and status and tcitation.status.upper() == status.upper() ):
        citationsFound.append(tcitation)

    citation = None

    if citationsFound:
      for tcitation in citationsFound:
        if tcitation.year == year:
          citation = tcitation
          break
      
      if not citation:
        if len(citationsFound) > 1:
          print "  Warning: match (but no full match) found for citation %s %s %s - picking first match." % (title, status, year)
        citation = citationsFound[0]
        
    return citation

  #
  # Assignment file stuff
  #
  
  def setAssignments(self):
    
    self.setSpecificAssignmentInfo()
    
    for assignment in self.assignmentFile.assignments:
      
      if self.assignmentValid(assignment):
        
        resonanceGroup = self.project.currentNmrProject.findFirstResonanceGroup(serial = assignment.spinSystemId)
        residue = self.chain.findFirstResidue(seqId = assignment.seqId)
        
        if resonanceGroup and residue:
          resonanceGroup.residue = residue
        else:
          print "  Error: could not set assignment:",
          if not resonanceGroup:
            print " resonance group with serial %s not found." % (str(assignment.spinSystemId))
          elif not residue:
            print " residue with seqId %d not found in chain '%s'." % (assignment.seqId,self.chain.code)
            
  def setSpecificAssignmentInfo(self):
    
    pass
        
  def assignmentValid(self,assignment):
    
    # Has to be format-specific
    
    return True
    
  #
  # Selection code...
  #

  def selectRawChemComps(self):
    
    rawChemCompList = self.chemCompFile.chemComps
    
    if len(self.chemCompFile.chemComps) > 1 and not self.minimalPrompts:
      
      ccDict = {}
      ccNames = []
      
      for chemComp in self.chemCompFile.chemComps:
        if hasattr(chemComp,'name') and chemComp.name:
          ccName = chemComp.name
        else:
          atomNameList = [atom.name for atom in chemComp.atoms]
          atomNameList.sort()
          ccName = string.join(atomNameList,',')
        
        ccNames.append(ccName)
        ccDict[ccName] = chemComp
          
      popup = self.multiDialog.MultiSelectionList(self.guiParent,ccNames,selectionDict=ccDict,
                                                  title = "Project '%s': " % self.project.name +'Select chemComps to import',
                                                  text = "Select the chemComps you want to import from this file.",
                                                  dismissText = 'Select none'
                                                  )

      if popup.isSelectedList:   
        rawChemCompList = popup.isSelectedList
    
    return rawChemCompList
    
  def selectDataSource(self, dataType = None, numDim = None, title = 'Select data source for experiment'):
  
    if not dataType:
    
      dataSources = self.experiment.sortedDataSources()
      
    else:
    
      kw = {}
      
      if numDim:
        kw['numDim'] = numDim
      if dataType:
        kw['dataType'] = dataType
    
      dataSources = self.experiment.findAllDataSources(**kw)
  
    if dataSources:

      (selectionList,selectionDict) = createSelection(dataSources)

      interaction = self.multiDialog.SelectionList(

                         self.guiParent,
                         selectionList,
                         selectionDict = selectionDict,
                         title = "Project '%s': " % self.project.name + title,
                         text = "Data sources for experiment %s:" % self.experiment.name,
                         dismissText = 'create new datasource',
                         urlFile = 'SelectDataSource'
                         
                         )

      if interaction.isSelected:
        self.dataSource = interaction.selection

  def setFidDataSource(self):
    
    if not self.fidDataSource:
      
      self.dataSource = None
      self.selectDataSource(dataType = 'FID',title = 'Select FID data source for experiment')
      
      if not self.dataSource:
        
        # self.formatLabel was now redundant fileType
        self.dataSource = Nmr.DataSource(self.experiment,
                                         name = self.format,
                                         numDim = self.numDim,
                                         dataType = 'FID')
        
            
      self.fidDataSource = self.dataSource
      
    self.experiment.rawData = self.fidDataSource
    
    self.setFidDataDims()
    
    self.setFidDataSourceSpecific()
    
  def setFidDataDims(self):
    
    #
    # Defaults for folding, complex
    #
    
    isComplex = 1

    for i in range(0,self.numDim):
    
      expDim = self.experiment.findFirstExpDim(dim = (i+1))
      fidDataDim = self.fidDataSource.findFirstDataDim(dim = expDim.dim)
      
      values = {'isComplex':     isComplex,
                'valuePerPoint': self.acqParsData.fidDimRefs[i].valuePerPoint,
                'numPoints':     self.acqParsData.fidDimRefs[i].numPoints,
                'numPointsValid':self.acqParsData.fidDimRefs[i].numPointsValid,
                'unit':          's'
                }

      if fidDataDim:
        
        for valueKey in values:
          if self.overwrite or not getattr(fidDataDim,valueKey):
            setattr(fidDataDim,valueKey,values[valueKey])
          
      else:
      
        values.update({'expDim': expDim, 'dim':expDim.dim})
        fidDataDim = Nmr.FidDataDim(self.fidDataSource, **values)
                
  def setExpDimRefs(self):

    folding = 0

    for i in range(0,self.numDim):
    
      expDim = self.experiment.findFirstExpDim(dim = (i+1))

      for expDimRef in expDim.expDimRefs:
      
        nucleus = expDimRef.isotopeCodes[0]
        
        if self.acqParsData.specRefs.has_key(nucleus):
        
          expDimRef.sf = self.acqParsData.specRefs[nucleus].specFreq
          expDimRef.baseFrequency = self.acqParsData.specRefs[nucleus].baseFreq
          expDimRef.folding = folding
          
        else:
        
          print "  Warning: no nucleus information for expDimRef %d" % (expDimRef.dim)

  def selectSpectrometer(self):
  
    if self.project.currentInstrumentStore:

      spectrometerList = self.project.currentInstrumentStore.findAllInstruments(className = 'NmrSpectrometer')
  
      if spectrometerList:
  
        (selectionList,selectionDict) = createSelection(spectrometerList)
  
        interaction = self.multiDialog.SelectionList(
  
                           self.guiParent,
                           selectionList,
                           selectionDict = selectionDict,
                           title = "Project '%s': " % self.project.name + 'Select spectrometer for experiment',
                           text = "Spectrometers available",
                           dismissText = 'Create new spectrometer',
                           urlFile = 'SelectSpectrometer')
  
        if interaction.isSelected:
          self.spectrometer = interaction.selection

  def setSpectrometer(self, manufacturer = 'Bruker', model = 'unknown', protonFreq = 500.00): 

    #
    # Look for spectrometer (and set and/or create)
    #

    if not self.experiment.spectrometer:
    
      self.spectrometer = None
    
      self.selectSpectrometer()
      
      if not self.spectrometer:
        
        #
        # TODO: create a Bruker People.Organisation object to link manufacturer to!
        #
        
        values = [manufacturer,model,protonFreq]
        valueInfo = ['Manufacturer','Model','Proton freq']
      
        if not self.minimalPrompts or model == 'unknown':
    
          interaction = self.multiDialog.EntryList(

                                 self.guiParent,
                                 values,
                                 valueInfo,
                                 title = 'Select spectrometer for experiment %s' % self.experiment.name,
                                 topText = 'Set spectrometer information')
      
          if interaction.values:
            manufacturer = interaction.values[0]
            model = interaction.values[1]
            protonFreq = returnFloat(interaction.values[2])

        nominalFreq = str((returnInt((protonFreq + 5)/ 10)) * 10)
        
        setCurrentStore(self.project,'InstrumentStore')        
        self.spectrometer = Instrument.NmrSpectrometer(self.project.currentInstrumentStore, name = manufacturer + '_' + model, model = model, nominalFreq = nominalFreq,protonFreq = protonFreq)

      self.experiment.spectrometer = self.spectrometer

  def setFreqDataSource(self):
    
    if not self.freqDataSource:
      
      self.dataSource = None
      self.selectDataSource(dataType = 'processed', title = 'Select frequency data source for experiment')
      
      if not self.dataSource:
      
        self.dataSource = Nmr.DataSource(self.experiment,
                                         name = importGeneral_kw,
                                         numDim = self.numDim,
                                         dataType = 'processed')
 
      self.freqDataSource = self.dataSource
      
    self.setFreqDataDims()
    
  def setFreqDataDims(self):

    for i in range(0,self.numDim):

      expDim = self.experiment.findFirstExpDim(dim = (i+1))
      freqDataDim = self.freqDataSource.findFirstDataDim(dim = expDim.dim)
      
      #

      # NOTE THAT isComplex default is set to FALSE! Assuming imaginary data will be deleted...
      #
      
      values = {'isComplex':       False,
                'valuePerPoint':   self.acqParsData.freqDimRefs[i].valuePerPoint,
                'numPoints':       self.acqParsData.freqDimRefs[i].numPoints,
                'numPointsOrig':   self.acqParsData.freqDimRefs[i].numPointsOrig,
                'phase0':          self.acqParsData.freqDimRefs[i].phase[0],
                'phase1':          self.acqParsData.freqDimRefs[i].phase[1],
                'unit':            'Hz'
                }

      if freqDataDim:
        
        for valueKey in values:
        
          setattr(freqDataDim,valueKey,values[valueKey])
          
      else:
      
        values.update({'expDim': expDim, 'dim':expDim.dim})
        freqDataDim = Nmr.FreqDataDim(self.freqDataSource, **values)

  def setAcqParsGeneral(self):
  
    self.experiment.numScans = self.acqParsData.numScans
    
    self.setAcqParsSpecific()
    
  def setAcqParsSpecific(self):
    
    #
    # Always format specific
    #
    
    pass
    
  def setFidDataSourceSpecific(self):
    
    #
    # Always format specific
    #
    
    pass
    
  def setFidDataSourceLocation(self,expFidName):

    if not self.fidDataSource.dataStore or self.overwrite:
      (dataUrlPath,localPath) = splitPath(self.fileName)

      #TODO this was originally in here - why? fileName = joinPath(fileDir,expFidName)
      dataStore = setDataSourceDataStore(self.fidDataSource, dataUrlPath, localPath)
      
      # TODO SET OTHER INFO ON RAW DATA!!

  def setProcParsGeneral(self):
     
    #
    # This code taken from analysis' ExternalParams...
    #

    if (self.experiment.numDim != self.procParsFile.fPars['numDim']):
      raise self.FormatConverterError('Experiment and dataSource numDim needs to be the same')
          
    # Set name of the data source
    if self.fileName:
      (fileDir,dataSourceName) = os.path.split(self.fileName)
    else:
      dataSourceName = self.format 
      
    self.dataSource = Nmr.DataSource(self.experiment, name=dataSourceName,
                                     numDim=self.numDim, dataType='processed')
    
    fPars = self.procParsFile.fPars

    # TODO probably need to split up inputFile into logical components depending on format... in fpars?
    if fPars.has_key('inputFile') and fPars['inputFile']:
      (dataUrlPath,localPath) = os.path.split(fPars['inputFile'])
            
    dataStore = setDataSourceDataStore(self.dataSource, dataUrlPath, localPath)
    dataStore.nByte = 4
    
    #
    # Set rest info from fPars, if possible
    #
    
    if fPars.has_key('bigEndian') and fPars['bigEndian'] is not None:
      dataStore.isBigEndian = fPars['bigEndian']
    elif fPars.has_key('swap') and fPars['swap']:
      dataStore.isBigEndian = not isBigEndian()
    else:
      dataStore.isBigEndian = isBigEndian()
      
    if fPars.has_key('integer') and fPars['integer']:
      dataStore.numberType = 'int'
    else:
      dataStore.numberType = 'float'
      
    if fPars.has_key('headerSize') and fPars['headerSize']:
      dataStore.headerSize = fPars['headerSize']
  
  def setProcParsDimInfo(self):
  
    #
    # At this stage should have expDimRefs... link them up to the ones
    # that were read in...
    #
    
    expDimRefs = None
    fPars = self.procParsFile.fPars

    self.dataDimRefs = [] # Keep track for possible later reference
    self.protonFreq = 500.0
    
    # Also track variables to set dataSource.dataStore NumericMatrix
    dataStoreInfo = {'blockSizes': [], 'isComplex': [], 'numPoints': []}
    
    while (not expDimRefs):

      interaction = self.multiDialog.ProcParToExpDim(self.guiParent,self.experiment,fPars)
      
      if not hasattr(interaction,'expDimRefs') or not interaction.expDimRefs:
        return None
      else:
        expDimRefs = interaction.expDimRefs
        
    for dim in range(0,len(expDimRefs)):
    
      expDimRef = expDimRefs[dim]
      expDim = expDimRef.expDim
      
      #
      # Reset expDimRef values...
      #
      
      expDimRef.sf = self.procParsFile.fPars['spectrometerFreq'][dim]
      
      if self.procParsFile.fPars['nucleus'][dim]:
        expDimRef.isotopeCodes = [self.procParsFile.fPars['nucleus'][dim]]
      elif expDimRef.refExpDimRef:
        expDimRef.isotopeCodes = [ass.isotopeCode for ass in expDimRef.refExpDimRef.expMeasurement.atomSites]
      
      if self.aliasing and self.aliasing.has_key(dim):
        expDimRef.minAliasedFreq = self.aliasing[dim][0]
        expDimRef.maxAliasedFreq = self.aliasing[dim][1]
      
      #
      # Set freqDataDim and dataDimRef
      #
      
      if self.procParsFile.fPars['nucleus'][dim] == '1H':
        self.protonFreq = self.procParsFile.fPars['spectrometerFreq'][dim]
        
      isComplex = False # TODO WHY!?!
            
      freqDataDim = Nmr.FreqDataDim(self.dataSource,dim=expDim.dim,
                                           numPoints=fPars['numPoints'][dim],
                                           isComplex=isComplex,
                                           numPointsOrig=fPars['numPoints'][dim],
                                           valuePerPoint=fPars['spectralWidth'][dim]/float(fPars['numPoints'][dim]),
                                           expDim=expDim)
                                           
      dataDimRef = Nmr.DataDimRef(freqDataDim, refPoint=fPars['refPoint'][dim],
                                      refValue=fPars['refPpm'][dim],
                                      expDimRef=expDimRef)
      
      # Set info for dataSource.dataStore.
      dataStoreInfo['blockSizes'].append(fPars['blockSize'][dim])
      dataStoreInfo['isComplex'].append(isComplex)
      dataStoreInfo['numPoints'].append(fPars['numPoints'][dim])
                                      
      self.dataDimRefs.append(dataDimRef)
      
    #
    # Set dataSource.dataStore info. Should be in correct order.
    #
    
    for attrName in dataStoreInfo.keys():
      setattr(self.dataSource,attrName,dataStoreInfo[attrName])

    return 1
  
  #
  # Procpar writing
  #
  
  def setProcPars(self): 
    
    #
    # TODO: Set general stuff!
    #
    
    self.setProcParsSpecific()
    
  def setProcParsSpecific(self):
  
    pass

  def setProcScript(self): 
    
    #
    # TODO: Set general stuff!
    #
    
    self.setProcScriptSpecific()

  def setProcScriptSpecific(self):
  
    pass

  def getComplexFactor(self,fidDataDim):
      
    #
    # Handle complex/real number of points
    #
    
    if fidDataDim.isComplex == True:
      return 2 # If complex multiple num points by 2 to get actual number of points!
      
    else:
      return 1

  def validDataDimRefs(self):
    
    #
    # Do some checks on the dataDimRefs that were passed in...
    #
    
    if isinstance(self.dataDimRefs, tuple):
      self.dataDimRefs = list(self.dataDimRefs)

    if isinstance(self.dataDimRefs, list):
      
      for dataDimRef in self.dataDimRefs:
        if not isinstance(dataDimRef,Nmr.DataDimRef):
          self.messageReporter.showError("Error","Need valid dataDimRef object(s) to write processing parameters.",parent = self.guiParent)
          return False

    else:
      self.messageReporter.showError("Error","Need list of dataDimRef object(s) to write processing parameters.",parent = self.guiParent)
      return False
      
    return True


  def getChemAtomOrSetSysName(self, searchKey, exportName=None, findIupac=False):
    
    if self.atomSysNamesByChemAtom.has_key(searchKey):
      atomSysName = self.atomSysNamesByChemAtom[searchKey]
      
    else:
      (chemAtomOrSet, namingSystem) = searchKey
      atomSysName = findChemAtomOrSetSysName(namingSystem, chemAtomOrSet,
                                             exportName, findIupac=findIupac)
                                             
      if atomSysName:
        self.atomSysNamesByChemAtom[searchKey] = atomSysName
    
    return atomSysName
  
  def getChemAtomSysNamesFromSet(self,chemAtomSet,chemCompVar,findIupac = False):
  
    chemAtomSysNameList = []
    
    chemAtoms = list(chemAtomSet.chemAtoms)
    
    if not chemAtoms:
      chemAtoms = []
      for chemAtomSet in chemAtomSet.chemAtomSets:        
        chemAtoms.extend(list(chemAtomSet.chemAtoms))
    
    namingSystem = chemCompVar.chemComp.findFirstNamingSystem(name = self.namingSystemName)
    if findIupac:
      iupacNamingSystem = chemCompVar.chemComp.findFirstNamingSystem(name = 'IUPAC')
    
    for chemAtom in chemAtoms:
    
      chemAtomSysName = self.getChemAtomOrSetSysName((chemAtom,namingSystem))
      
      if findIupac and not chemAtomSysName:
      
        chemAtomSysName = self.getChemAtomOrSetSysName((chemAtom,iupacNamingSystem))
      
      if chemAtomSysName:
        chemAtomSysNameList.append(chemAtomSysName.sysName)
      else:
        chemAtomSysNameList.append(None)
    
    return chemAtomSysNameList

  def findChemCompVarSysName(self,namingSystemName,chemCompVar):

    sysName = findChemCompVarSysName(namingSystemName,chemCompVar)
    
    return sysName

  def getExportSeqCode(self,firstSeqCode,residue):
    
    if self.useCcpnChainInfo:
      seqCode = residue.seqCode
    else:
      seqCode = firstSeqCode + residue.seqId - 1
   
    return seqCode 

  def getResidueLabelOneLetter(self,resonanceToAtom):
  
    residue = resonanceToAtom.chain.findFirstResidue(seqId = resonanceToAtom.seqId)
    residueLabel = residue.molResidue.chemComp.code1Letter
    
    if not residueLabel:
      residueLabel = 'X'
      
    return residueLabel
    
  def getResidueLabelThreeLetter(self,resonanceToAtom):
  
    residue = resonanceToAtom.chain.findFirstResidue(seqId = resonanceToAtom.seqId)
    residueLabel = residue.molResidue.chemComp.code3Letter
    
    if not residueLabel:
      residueLabel = 'XXX'
      
    return residueLabel

  def getMultiChainFormatPresetChainMapping(self,chainList):
  
    mappingChainDict = {}
    
    #
    # Only set this automatically if self.minimalPrompts is on
    # Now using seqCode... should be seqId here?
    #
        
    if self.minimalPrompts or (hasattr(self,'forceDefaultChainMapping') and self.forceDefaultChainMapping):
            
      for chain in chainList:

        mappingChainDict[chain] = (chain.code,chain.sortedResidues()[0].seqCode)

    return mappingChainDict  

  def getSingleChainFormatPresetChainMapping(self,chainList):
  
    mappingChainDict = {}
    
    #
    # This value determines how far apart the seqCodes for different
    # chains have to be spaced! Could become user-set...
    #
    
    linkResidueLength = 50
    
    #
    # Only try to set this automatically if self.minimalPrompts is on
    #
        
    if self.minimalPrompts:
      
      seqCodeRange = []      
      setMapping = 1
      
      for chain in chainList:
      
        validRange = 1
        
        #
        # Assuming seqCodes to be ordered from start to end...
        # no overlap is allowed
        #
        
        lowerSeqCode = chain.sortedResidues()[0].seqCode
        upperSeqCode = chain.sortedResidues()[-1].seqCode
        
        for (tempLowerSeqCode,tempUpperSeqCode) in seqCodeRange:
        
          if (lowerSeqCode >= tempLowerSeqCode and lowerSeqCode <= tempUpperSeqCode) or \
             (upperSeqCode >= tempLowerSeqCode and upperSeqCode <= tempUpperSeqCode):
          
            validRange = 0
            break

        if validRange:
          
          seqCodeRange.append((lowerSeqCode-linkResidueLength,upperSeqCode+linkResidueLength))
        
        else:
        
          setMapping = 0
          break
    
      if setMapping:
      
        #
        # Set the same as the seqCodes...
        #
      
        for chain in chainList:
        
          mappingChainDict[chain] = (None,chain.sortedResidues()[0].seqCode)
          
    #
    # If forceDefaultChainMapping is on, then set them automatically...
    # can only be used when writing a sequence.
    #
    
    if not mappingChainDict and hasattr(self,'forceDefaultChainMapping'):
    
      if self.forceDefaultChainMapping:
      
        seqCode = 0
      
        for chain in chainList:
        
          for tCode in range(1,1001,linkResidueLength):
            if seqCode < tCode:
              seqCode = tCode
              break
        
          mappingChainDict[chain] = (None,seqCode)
          seqCode += len(chain.residues) + linkResidueLength

    return mappingChainDict  

  def setGenericTag(self,label,object,setTag,isUnique = True, AppDataClass = Implementation.AppDataString):
    
    if setTag:
      appData = object.findFirstApplicationData(application = self.format, keyword = label)
      
      if isUnique:
        object.removeApplicationData(appData)
        appData = None
        
      if not appData:
        object.addApplicationData(AppDataClass(application = self.format, keyword = label, value = setTag))

  def setDateTimeStamp(self,*args,**keywds):
    
    self.setGenericTag(dateTimeStamp_kw,*args,**keywds)

  def setIoCycleTag(self,*args,**keywds):
     
    self.setGenericTag(ioCycleTag_kw,*args,**keywds)

  def setSeqAssignTag(self,*args,**keywds):
     
    self.setGenericTag(seqAssignTag_kw,*args,**keywds)


  
  
def filteredNameInfo(resonanceToAtoms, mode=None):
  """ get list of chain,seqId,atomName tuples
  free of duplicates, compressed according to compression mode:
  'compress', 'preset' or None
  Used in FormatConverter
  """
  
  resonanceToAtoms = list(resonanceToAtoms)
  
  result = list()
  
  if mode in ('compress', 'preset'):
    # Merge resonances from shared atomSets
    # if mode is 'preset' only do it if r2a.useAtomSet
    
    chemAtomSets = [r2a.getChemAtomSet() for r2a in resonanceToAtoms]
    
    #for r2a in resonanceToAtoms:
    #  print ('   +++ ', r2a.resonance.serial, r2a.seqId, r2a.atomName, r2a.atomSetName,
    #                  r2a.chemAtom and  r2a.chemAtom.name,
    #                  r2a.chemAtomSet and  r2a.chemAtomSet.name,
    #                  r2a.useAtomSetName,
    #                  [x.serial for x in r2a.otherLinkedResonances],
    #                  [x.serial for x in r2a.otherGroupResonances])
    
    keyMap = {}
    
    ii = 0
    while ii < len(resonanceToAtoms):
      r2a = resonanceToAtoms[ii]
      atomName = r2a.atomName
      if mode == 'compress' or r2a.useAtomSetName:
        chemAtomSet = chemAtomSets[ii]
        if chemAtomSet:
          chain = r2a.chain
          seqId = r2a.seqId
          jj = ii + 1
          while jj < len(resonanceToAtoms):
            otherR2a = resonanceToAtoms[jj]
            if (otherR2a.chain == chain and otherR2a.seqId == seqId
                and chemAtomSets[jj] == chemAtomSet):
              # resonanceToAtoms share a chemAtomSet. Merge
              if r2a.atomSetName:
                # use atomSetName instead of atomName
                atomName = r2a.atomSetName
              # remove duplicates
              del resonanceToAtoms[jj]
              del chemAtomSets[jj]
            else:
              jj += 1
      key = (r2a.chain, r2a.seqId, atomName)
      tt = keyMap.get(key)
      if (tt is None or tt[0].atomSetName == atomName or 
          not tt[0].useAtomSetName or not tt[0].getChemAtomSet()):
        # favour pseudoatoms by selecting later, or incompletely set r2a 
        keyMap[key] = (r2a, atomName)
      ii += 1
    
    result = [keyMap[x] for x in sorted(keyMap)]
    
  else:
    # no compression, just transform, remove duplicates, and sort
    result = [(r2a, r2a.atomName) for r2a in resonanceToAtoms]
  
  #
  # transform to list, sort, and return
  return sorted(result)

      
def cleanAtomNameLists (resNames, resAtomSetNames):
  """Remove duplicate assignments
  Used in FormatConverter
  """
  debug = 'Q' in ''.join(str(x) for x in resAtomSetNames)
  
  if debug:
    print '\n@@@1-', resNames
    print '@@@1 ', resAtomSetNames
  
  resNames = list(resNames)
  resAtomSetNames = list(resAtomSetNames)
  
  # First remove individual atoms if their atomSetNames are also present
  ii = 0
  while ii < len(resNames):
    ss = resNames[ii]
    jj = ii +1
    while jj < len(resNames):
      if resAtomSetNames[jj] == ss:
        del resNames[jj]
        del resAtomSetNames[jj]
      else:
        jj += 1
    ii +=1
  
  if debug:
    print '\n@@@2-', resNames
    print '@@@2 ', resAtomSetNames
  
  # Next remove duplicates that share an AtomSetName
  ii = 0
  while ii < len(resNames):
    resAtomSetName = resAtomSetNames[ii]
    if resAtomSetName is not None:
      jj = ii+1
      while jj < len(resNames):
        if resAtomSetNames[jj] == resAtomSetName:
          del resNames[jj]
          del resAtomSetNames[jj]
          resNames[ii] = resAtomSetName
        else:
          jj += 1
    ii += 1
  
  if debug:
    print '\n@@@3-', resNames
    print '@@@3 ', resAtomSetNames
    
  # Finally remove duplicate names
  ii = 0
  while ii < len(resNames):
    ss = resNames[ii]
    jj = ii +1
    while jj < len(resNames):
      if resNames[jj] == ss:
        del resNames[jj]
        del resAtomSetNames[jj]
      else:
        jj += 1
    ii +=1
  
  if debug:
    print '\n@@@4-', resNames
    print '@@@4 ', resAtomSetNames
  
  #
  return resNames, resAtomSetNames


def getResonanceAtomMap(namingSystemName, resonances,
                        individualAtomsIfNoSet=True, individualAtoms=False,
                        ignoreChemCompsWithNoSysNames=False,
                        compressResonances=False):
  """ get resonance: ResonanceToAtoms map for specific naming system
  resonances is list or tuple of resonances or fixedResonances
  Stand-alone replacement for DataFormat.getResonanceAtomLinks,
  minus original atom names and sysname caching
  """
  
  print '### starting getResonanceAtomMap'
    
  # set up
  resonanceToAtoms = {}
  atomSysNameMap = {}
  
  preferredAtomNameDict = {}
  
  if not resonances:
    return resonanceToAtoms
  
  # set high level variables
  project = resonances[0].root
  
  usedResonanceSets = set()
  
  for resonance in resonances: 
          
    # Make assignment map for resonance
    resonanceSet = resonance.resonanceSet
    if resonanceSet:
      # We have assignment, so dd to map (otherwise skip)
      
      # Only do each ResonanceSet once
      if resonanceSet in usedResonanceSets:
        continue
      else:
        usedResonanceSets.add(resonanceSet)
      
      # set some variables
      resSetResonances = resonanceSet.sortedResonances()
      atomSets = resonanceSet.sortedAtomSets()
      chemAtomSet = getChemAtomSetFromAtomSets(atomSets)
      
      # get common residue and naming system
      allAtoms = set(x for y in atomSets for x in y.atoms)
      residues = set(x.residue for x in allAtoms)
      residue = residues.pop()
      exportName = getResName(residue.chain.code, residue.seqCode,
                     (chemAtomSet or atomSets[0].findFirstAtom()).name)
      
      namingSystem = residue.chemCompVar.chemComp.findFirstNamingSystem(
                      name=namingSystemName)
      if not namingSystem:
        print ("  Warning: No %s naming system available for resonance %s" 
               % (namingSystemName, exportName))

        if ignoreChemCompsWithNoSysNames:
          print " - not written."
          continue
        else:
          print " - trying to use default atom names."
          
      if residues:
        # There were atoms from different residues - set residue back to None
        residue = None
      
      
      
      # Start setting mappings:
      
      if (len(resSetResonances) == 1 and residue and 
          allAtoms == residue.atoms):
        # Residue-level constraint. Set it.
        
        resonanceToAtoms[resonance] = ll =[]
        ll.append(ResonanceToAtoms(resonance, residue.chain, residue.seqId, 
                                   None, atomSetName=None))



      elif (residue and chemAtomSet and len(atomSets) == 2 
            and len(resSetResonances) <= 2):
        # Case where 1/2 resonances, 2 atomSets, one chemAtomSet
        # Typically prochiral pair
        
        # set up parameters
        
        # get chemAtomSet atomSysName and update atomSysNameMap
        if chemAtomSet in atomSysNameMap:
          # NB Do not use .get. Might be in map with value None
          atomSysName = atomSysNameMap[chemAtomSet]
        else:
          atomSysName = findChemAtomOrSetSysName(namingSystem, chemAtomSet, 
                                                 exportName, 
                                                 findIupac=not individualAtomsIfNoSet)
                                                 #findIupac=not compressIfHasSet)
          atomSysNameMap[chemAtomSet] = atomSysName
        
        # get chemAtomSetName
        if atomSysName:
          chemAtomSetName = atomSysName.sysName
        elif individualAtomsIfNoSet:               
          chemAtomSetName = None
          print ("  Warning: no %s sysNames found for %s (%s). ChemAtomSet not used." 
                 % (namingSystemName, chemAtomSet.name, exportName))
        else:                                          
          chemAtomSetName = chemAtomSet.name
          print ("  Warning 1: no %s sysNames found for %s (%s)... using CCPN name." 
                 % (namingSystemName, chemAtomSet.name, exportName))
          
        # Make list of chemAtomSet contents (chemAtoms or chemAtomSets)
        chemAtomOrSetList = (chemAtomSet.sortedChemAtomSets() or 
                             chemAtomSet.sortedChemAtoms())

        # Get their atom names in preferred order
        preferredAtomNames = []
        for reson in resSetResonances:
        
          ss = preferredAtomNameDict.get(reson)
          if ss is not None:
            # Already mapped as covalently bound prochiral resonance
            preferredAtomNames.append(ss)
            
          else:
            # set preference order by standard heuristics
            ii = getAmbigNamingPriority(reson, chemAtomSet)
            preferredAtomNames.append(chemAtomOrSetList[ii].name)
              
        # Reverse chemAtomOrSetList if necessary, to make it fit preferred order
        if chemAtomOrSetList[-1].name == preferredAtomNames[0]:
          chemAtomOrSetList.reverse()
        

        # Now loop over chemAtoms or chemAtomSets and do the work
        assignResonances = list(resSetResonances)  # copy, must be modified
        for chemAtomOrSet in chemAtomOrSetList:
          
          # select resonance for ChemAtomOrSet
          chemAtomOrSetName = chemAtomOrSet.name
          if chemAtomOrSetName in preferredAtomNames:
            # if name in list take out corresponding resonance for processing
            useIndex = preferredAtomNames.index(chemAtomOrSetName)
                 
          else:
            # take first resonance that does *not* fit a prefered name position
            validIndices = range(len(preferredAtomNames))
            for caos in chemAtomOrSetList:
              ss = caos.name
              if ss in preferredAtomNames:
                removeIndex = preferredAtomNames.index(ss)
                validIndices.pop(validIndices.index(removeIndex))
             
            useIndex = validIndices[0]
          #
          reson = assignResonances.pop(useIndex)
          del preferredAtomNames[useIndex]


          # Add to preferredAtomNameDict based on covalently bound resonances
          # Only relevant for 'deep' prochiral case (e.g. Valine methyls)...
          for covalentResonance in reson.covalentlyBound:
            # NB from July 2012, FixedResonances have covalentlyBound too
          
            covalentResSet = covalentResonance.resonanceSet
            if covalentResSet:
            
              covalentPreferredAtomName = None
              for covalentAtomSet in covalentResSet.atomSets:
              
                refChemAtom = covalentAtomSet.findFirstAtom().chemAtom
                if refChemAtom.chemAtomSet:
                
                  for refChemBond in refChemAtom.chemBonds:
 
                    atSet = set(refChemBond.chemAtoms)
                    atSet.remove(refChemAtom)
                    otherRefChemAtom = atSet.pop()
                    
                    if hasattr(otherRefChemAtom,'chemAtomSet'):
                      # the 'hasattr' required to weed out LinkAtoms
 
                      if chemAtomOrSet is otherRefChemAtom.chemAtomSet:
                        covalentPreferredAtomName = refChemAtom.name

                      elif chemAtomOrSet is otherRefChemAtom:
                        covalentPreferredAtomName = refChemAtom.chemAtomSet.name
 
                      if covalentPreferredAtomName:
                        preferredAtomNameDict[covalentResonance] = covalentPreferredAtomName
                        break

                  if covalentPreferredAtomName:
                    break
          
          
          # Set up to get the chemAtom(Set)SysNames
          tempChemAtomOrSetList = [chemAtomOrSet]
          # Decide whether to split this list into individual component Atom(Set)s
          if chemAtomOrSet.className == 'ChemAtomSet':
            splitChemAtomSet = False
            if individualAtoms:
              # switch says always split
              splitChemAtomSet = 'deep'
            elif individualAtomsIfNoSet:
              if not findChemAtomOrSetSysName(namingSystem, chemAtomOrSet):
                # Split, as naming system has no name for AtomSet
                splitChemAtomSet = 'shallow'
            
            if splitChemAtomSet:
              # split to direct content of ChemAtom(Set)s
              tempChemAtomOrSetList = (chemAtomOrSet.sortedChemAtomSets() or
                                       chemAtomOrSet.sortedChemAtoms())
              if (splitChemAtomSet == 'deep' 
                  and tempChemAtomOrSetList[0].className == 'ChemAtomSet'):
                # split all the way down to individual atoms
                ii = 0
                while ii < len(tempChemAtomOrSetList):
                  caos = tempChemAtomOrSetList[ii]
                  if caos.className == 'ChemAtomSet':
                    del tempChemAtomOrSetList[ii]
                    tempChemAtomOrSetList.extend(caos.sortedChemAtomSets())
                    tempChemAtomOrSetList.extend(caos.sortedChemAtoms())
                  else:
                    # we have a ChemAtom here. We could probably break, 
                    # but just to be safewe continue till the end
                    ii += 1
          
          # get sysName and set mapping
          for tempChemAtomOrSet in tempChemAtomOrSetList:
            atomSysName = findChemAtomOrSetSysName(namingSystem, 
                                                   tempChemAtomOrSet, 
                                                   exportName,
                                                   findIupac=True)
            if atomSysName:
              refSysName = atomSysName.sysName
            
            else:
              # If really nothing, use CCPN atom name
              refSysName = tempChemAtomOrSet.name
              print ("  Warning 2: no %s sysNames found for %s (%s)... using CCPN name." 
                     % (namingSystemName, tempChemAtomOrSet.name, exportName))
                         
            # Finally set mapping
            r2a = ResonanceToAtoms(reson, residue.chain, residue.seqId, 
                                    refSysName, atomSetName=chemAtomSetName,
                                    chemAtomOrSet=tempChemAtomOrSet,
                                    atomType=tempChemAtomOrSet.elementSymbol)
            ll = resonanceToAtoms.get(reson)
            if ll:
              ll.append(r2a)
            else:
              resonanceToAtoms[reson] = [r2a]
            # Add otherResonances
            for otherResonance in resSetResonances:
              if otherResonance != reson:
                r2a.addOtherGroupResonance(otherResonance)

          if not assignResonances:
            break
        
          
      elif residue and len(atomSets) == 1 and len(resSetResonances) == 1:
        # NBNB Added and not tested in situ - RHF 12/11/13
        
        #
        # 1 resonance, 1 atomset. Simplest case.
        # Needs special handling to deal with isopropyl groups
        #
        
        atomSet = atomSets[0]
        resonance = resSetResonances[0]
        
        seqId = residue.seqId
        chain = residue.chain

        exportName = getResName(chain.code,residue.seqCode,
                                atomSet.name or atomSet.findFirstAtom().name)
                                
        
        chemAtom = atomSet.findFirstAtom().chemAtom
        
        # NB all atoms should ALWAYS point to the same ChemAtomSet, 
        # but it costs nothing ot be careful.
        aSet = set(x.chemAtom.chemAtomSet for x in atomSet.atoms)
        if len(aSet) == 1:
          chemAtomSet = aSet.pop()
        else:
          chemAtomSet = None
        
        if len(atomSet.atoms) == 1:
          
          atomSysName = findChemAtomOrSetSysName(namingSystem, chemAtom)
          if atomSysName:
            refSysName = atomSysName.sysName
          else:
            refSysName = chemAtom.name
        
          if chemAtomSet:
            if compressResonances and individualAtomsIfNoSet:
              chemAtomSetSysName = findChemAtomOrSetSysName(namingSystem, 
                                                            chemAtomSet)
              if chemAtomSetSysName:
                chemAtomSetName = chemAtomSetSysName.sysName
              else:
                chemAtomSet = None
                chemAtomSetName = None
                print "  Warning: no chemAtomSetSysName found for %s ... ChemAtomSet not used." % exportName
              
            else:
              chemAtomSetSysName = findChemAtomOrSetSysName(namingSystem, 
                                                            chemAtomSet, 
                                                            exportName, 
                                                            findIupac=True)    
              if chemAtomSetSysName:                                                                                  
                chemAtomSetName = chemAtomSetSysName.sysName 
              else:                                                                                            
                chemAtomSetName = chemAtomSet.name                                                             
                print "  Warning: no chemAtomSetSysName found for %s ... using CCPN name." % exportName   

          else:  
            chemAtomSetName = None
          
        elif chemAtomSet and not individualAtoms:
          #assert chemAtomSet and len(atomSet.atoms) > 1
          
          atomSysName = findChemAtomOrSetSysName(namingSystem, chemAtomSet)
          if atomSysName:
            refSysName = atomSysName.sysName
          else:
            refSysName = chemAtom.name
          
          deepSet = chemAtomSet.chemAtomSet
          if deepSet:
            chemAtomSet = deepSet
            
          chemAtomSetSysName = findChemAtomOrSetSysName(namingSystem, 
                                                        chemAtomSet)
          if chemAtomSetSysName:
            chemAtomSetName = chemAtomSetSysName.sysName
          else:
            chemAtomSetName = chemAtomSet.name
          
        else:
          #assert (individualAtoms or not chemAtomSet) and len(atomSet.atoms) > 1
          # Deal with atoms individualy
          # chemAtomSet should always exist here, but it costs nothing to be prepared.
                            
          chemAtomSetName = None
          if chemAtomSet:
            chemAtomSetSysName = findChemAtomOrSetSysName(namingSystem, chemAtomSet)
            if chemAtomSetSysName:
              chemAtomSetName = chemAtomSetSysName.sysName
            else:
              chemAtomSetName = chemAtomSet.name
          
          
          for chemAtom in chemAtomSet.chemAtoms:

            exportName = getResName(chain.code,residue.seqCode,chemAtom.name)
            atomSysName = findChemAtomOrSetSysName(namingSystem, chemAtom, 
                                                   exportName, findIupac=True)
         
            if atomSysName:
              refSysName = atomSysName.sysName
            else:
              print ("  Warning: no %s sysName for %s... using CCPN name." % 
                     (namingSystemName,exportName))
              refSysName = chemAtom.name
                                
            chemAtomOrSet = chemAtomSet or chemAtom
            r2a = ResonanceToAtoms(resonance, chain, seqId, refSysName,
                                   atomSetName=chemAtomSetName, 
                                   chemAtomOrSet=chemAtomOrSet, 
                                   atomType=chemAtomOrSet.elementSymbol)
            ll = resonanceToAtoms.get(resonance)
            if ll:
              ll.append(r2a)
            else:
              resonanceToAtoms[resonance] = [r2a]
          
          # Done with this resonance
          continue
          
        # set resonanceToAtoms
        chemAtomOrSet = chemAtomSet or chemAtom
        r2a = ResonanceToAtoms(resonance, chain, seqId, refSysName,
                               atomSetName=chemAtomSetName, 
                               chemAtomOrSet=chemAtomOrSet, 
                               atomType=chemAtomOrSet.elementSymbol)
        ll = resonanceToAtoms.get(resonance)
        if ll:
          ll.append(r2a)
        else:
          resonanceToAtoms[resonance] = [r2a]
                                 

      else:
        # This is the generic case: will be clear if only one atomSet,
        # otherwise very ambiguous
        
        
        #
        # Treat everything the same: use chemAtomSet if multiple atoms for
        # chemAtomSet, otherwise use chemAtom.
        # Ambiguous so link whatever to whatever!
        #

        # Loop over atomSets and do the work
        for atomSet in atomSets:
          
          # set up
          sysNamesLink = []
          chemAtomSet = None
          chemAtomSetName = None
          residue = atomSet.findFirstAtom().residue
          chain = residue.chain       
          
          
          #
          # Check whether sys names for individual atom(s) or for chemAtomSet. 
          #

          if len(atomSet.atoms) > 1 and not individualAtoms:
            
            chemAtomSet = atomSet.findFirstAtom().chemAtom.chemAtomSet
            chemAtomSetName = chemAtomSet.name
            if chemAtomSet.chemAtomSet:
              # If there is a deep chemAtomSet, try to use this. 

              if (namingSystem and 
                  checkChemAtomsForDeepSet(chemAtomSet.chemAtomSet, atomSets)):
                
                #
                # This code might cause some problems for HD*|HE* for aromatics...
                # Works fine if namingSystem DIANA - will use QR
                # Never use this name as reference - better to split up!
                #
                chemAtomSetSysName = findChemAtomOrSetSysName(namingSystem, 
                                                              chemAtomSet)
                if chemAtomSetSysName:
                  chemAtomSet = None
                  chemAtomSetName = chemAtomSetSysName.sysName
                  sysNamesLink = [(chemAtomSetName, chemAtomSet)]
            
            exportName = getResName(chain.code, residue.seqCode,
                                           chemAtomSetName)
            
            # Find the chemAtomSet sysNames (if still relevant after above)
            if not sysNamesLink:

              atomSysName = findChemAtomOrSetSysName(namingSystem, 
                              chemAtomSet, exportName,
                              findIupac=not individualAtomsIfNoSet)
              if atomSysName:
                sysNamesLink = [(atomSysName.sysName,chemAtomSet)]
                  
              else:
                if chemAtomSet.chemAtomSet and not individualAtomsIfNoSet:
                  print ("  Warning 3: no %s sysNames for %s ... using CCPN name ."
                         % (namingSystemName, exportName), chemAtomSetName)
                  sysNamesLink = [(chemAtomSet.name,chemAtomSet)]

                else:
                  chemAtomSetName = None
                  print ("  Warning: no %s sysNames for %s... using individual atom name(s)."
                         % (namingSystemName, exportName))
          
          #
          # If no chemAtomSet sysName(s) were found, do it on the individual atom level...
          #
          
          if not sysNamesLink:
            
            if chemAtomSet:
              chemAtoms = chemAtomSet.chemAtoms
            else:
              chemAtoms = [x.chemAtom for x in atomSet.sortedAtoms()]
            
            for chemAtom in chemAtoms:

              exportName = getResName(chain.code, residue.seqCode,
                                             chemAtom.name)  
              atomSysName = findChemAtomOrSetSysName(namingSystem, 
                              chemAtom, exportName, findIupac=True)
         
              if atomSysName:
                sysNamesLink.append((atomSysName.sysName, chemAtom))
              else:
                sysNamesLink.append((chemAtom.name,chemAtom))
                print ("  Warning 4: no %s sysName for %s... using CCPN name." 
                       % (namingSystemName, exportName), chemAtomSetName)
                
            #
            # Still set the chemAtomSetName if there is one, for compressing resonances
            # 
            
            if (not chemAtomSetName and namingSystem and not individualAtoms
                and chemAtom.chemAtomSet):
              
              chemAtomSetSysName = findChemAtomOrSetSysName(namingSystem, 
                                                            chemAtom.chemAtomSet)
              if chemAtomSetSysName:
                chemAtomSetName = chemAtomSetSysName.sysName

          sysNamesLink.sort()
          for reson in resSetResonances:
            for (sysName,chemAtomOrSet) in sysNamesLink:
              # Finally set mapping
              
              r2a = ResonanceToAtoms(reson, chain, residue.seqId, sysName, 
                                     atomSetName=chemAtomSetName,
                                     chemAtomOrSet=chemAtomOrSet, 
                                     atomType=chemAtomOrSet.elementSymbol)
              ll = resonanceToAtoms.get(reson)
              if ll:
                ll.append(r2a)
              else:
                resonanceToAtoms[reson] = [r2a]
        #
        # Set some additional information
        #
        
        for reson in resSetResonances:
          for otherResonance in resSetResonances:
            if otherResonance != reson:
              resonanceToAtoms[reson][-1].addOtherLinkedResonance(otherResonance)
            
              
    else:
      # If there is no resonanceSet for the resonance, then ignore it
      print ("  Warning: No assignment for resonance %s (%s) - ignored" 
             % (resonance.serial, resonance.name))
  
  
  #
  # Debug code...
  # 
  
  #printList = []
  #tag = 'resonanceSerial'
  #for res in sorted(resonanceToAtoms, key=operator.attrgetter(tag)):
  #  ll = [(rta.chain and rta.chain.code, rta.seqId, rta.atomName, rta.atomSetName) 
  #        for rta in resonanceToAtoms[res]]
  #  ll2 = [ll, getattr(res,tag), res.name, len(ll)]
  #  printList.append(ll2)
  
  #printList.sort()
  
  #for ll in printList:
  #  print   '+++', ll[1], ll[2], ll[3], ['%s.%s.%s.%s' % tt for tt in ll[0]]
  
  #
  return resonanceToAtoms
  


