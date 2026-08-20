"""
======================COPYRIGHT/LICENSE START==========================

Util.py: General utility functions

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

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
- PDBe website (http://www.ebi.ac.uk/pdbe/)

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

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================
"""

import re, sys, copy, os

from memops.api import Implementation

from memops.universal.MessageReporter import messageReporter

from ccp.api.molecule import ChemComp
from ccp.general.Geometry import calcTorsionAngleDegrees

from ccp.general.Constants import voidAtomSysName
from ccp.general.Constants import ccpCodeToCode1LetterDict


#
# Generally useful functions
#

def setCurrentStore(project, className, linkName = None, createNew = True):

  currentStore = None

  if className:

    currentName = "current%s" % (className)

    if not linkName:
      if className.endswith('y'):
        linkName = "sorted%sies" % className[:-1]

      else:
        linkName = "sorted%ss" % className

    currentStore = getattr(project,currentName)

    if not currentStore:
      storeList = getattr(project,linkName)()
      if storeList:
        setattr(project,currentName,storeList[0])
      elif createNew:
        # this automatically sets project.current..Store
        if className == 'Classification':
          newStore = getattr(project,"new%s" % className)(namingSystem=project.name)
        else:
          newStore = getattr(project,"new%s" % className)(name=project.name)

      currentStore = getattr(project,currentName)

  return currentStore

def findAllSysNamesByChemAtomOrSet(chemComp, chemAtomOrSetList, namingSystemName):

  #
  # Returns all sysNames for a chemComp based on a list of chemAtoms or chemAromSets.
  # Returns None if no direct matches to the namingSystem
  #

  namingSystem = chemComp.findFirstNamingSystem(name=namingSystemName)
  chemAtomSysNames = []

  if namingSystem:
    for chemAtom in chemAtomOrSetList:
      casn = findChemAtomOrSetSysName(namingSystem, chemAtom)

      if casn:
        chemAtomSysNames.append(casn)

  return chemAtomSysNames

def findChemAtomSysName(namingSystem,searchDict):

  return findAtomSysName(namingSystem,searchDict,'atomReference')

def findChemAtomSetSysName(namingSystem,searchDict):

  return findAtomSysName(namingSystem,searchDict,'atomSetReference')

def findAtomSysName(namingSystem,searchDict,refType):
  
  casn = None
  
  if namingSystem:

    casn = namingSystem.findFirstAtomSysName(**searchDict)
  
    if not casn:
      newNamingSystem = getattr(namingSystem,refType)
      if newNamingSystem:
        casn = findAtomSysName(newNamingSystem,searchDict,refType)
      else:
        casn = None
    
    if casn and casn.sysName == voidAtomSysName:
      casn = None

  return casn


def findChemAtomOrSetSysName(namingSystem, chemAtomOrSet, warnAtomId=None,
                             findIupac=False):
  """Input: ChemComp. ChemAtom or ChemAtomSet
            ChemComp.NamingSystem
            Boolean findIUPAC. Check IUPAC if nothing found in namingSystem?
     Output: ChemComp.AtomSysName matching chemAtomOrSet
  """
  
  # initialise
  result = None
  
  # set link name for finding next naming system depending on input type
  if chemAtomOrSet.className == 'ChemAtomSet':
    nextSysTag = 'atomSetReference'
  else:
    nextSysTag = 'atomReference'
  
  # Look through successive namingSystems for a match
  name = chemAtomOrSet.name
  subType = chemAtomOrSet.subType
  ns = namingSystem
  while bool(ns and not result):
  
    result = ns.findFirstAtomSysName(atomName=name, atomSubType=subType)
    ns = getattr(ns, nextSysTag)
  
  if findIupac and (not result or result.sysName == voidAtomSysName):
    # Now try IUPAC instead
    
    if warnAtomId:
      if namingSystem:
        print ("  Warning: no %s atomSysName available for '%s', trying IUPAC."
               % (namingSystem.name, warnAtomId))
      else:
        print ("  Warning: no namingSystem found for '%s', trying IUPAC."
               % warnAtomId)

    ns = chemAtomOrSet.chemComp.findFirstNamingSystem(name='IUPAC')
    if ns:
      result = ns.findFirstAtomSysName(atomName=name, atomSubType=subType)
  
  if result and result.sysName == voidAtomSysName:
    result = None
  #
  return result
  
def findChemCompSysName(namingSystemName,chemComp):

  sysName = None
  namingSystem = chemComp.findFirstNamingSystem(name = namingSystemName)
  if namingSystem:
    for ccsn in namingSystem.chemCompSysNames:
      if not ccsn.specificChemCompVars:
        sysName = ccsn.sysName

  return sysName

def findChemCompVarSysName(namingSystemName,chemCompVar):

  sysName = None
  namingSystem = chemCompVar.chemComp.findFirstNamingSystem(name = namingSystemName)
  if namingSystem:
    for ccsn in namingSystem.chemCompSysNames:
      if chemCompVar in ccsn.specificChemCompVars:
        sysName = ccsn.sysName
        break

    if not sysName:
      sysName = findChemCompSysName(namingSystemName,chemCompVar.chemComp)

  return sysName

def getResonancesFromPairwiseConstraintItem(pairwiseConstraintItem):

  resonanceList = pairwiseConstraintItem.orderedResonances
  if not resonanceList:
    resonanceList = pairwiseConstraintItem.sortedResonances()

  return list(resonanceList)

def setUniqueAppData(appDataClass,object,application,keyword,value):

  """
  Checks whether a unique application/keyword application data object is already available
  for the current object. If the value of this application data is different from the new
  value it will reset it. If no application data was found it will create a new one.
  """

  appData = object.findFirstApplicationData(application = application, keyword = keyword)

  if appData:
    object.removeApplicationData(appData)

  appData = getattr(Implementation,appDataClass)(application = application,keyword = keyword,value = value)
  object.addApplicationData(appData)

def getOtherAtom(refAtom,bond):

  """
  Returns atom at other end of bond. Returns None if not present
  """

  bondedAtoms = bond.sortedChemAtoms()

  otherAtom = None

  if bondedAtoms[0] == refAtom:
    otherAtom = bondedAtoms[1]
  elif bondedAtoms[1] == refAtom:
    otherAtom = bondedAtoms[0]

  return otherAtom

def getAtomPositionType(atomName):

  """
  For protein atoms and CCPN atom names, returns whether atom is sidechain or backbone.
  """

  if atomName in ('N','CA','C','O','HA','HA2','HA3','H','H1','H2','H3',"O'","O''","H''"):
    atomType = 'bb'
  else:
    atomType = 'sc'

  return atomType

def getSortContactDist(residueKey1,residueKey2,maxRangeCheck):

  """
  Based on input of two CCPN MolSystem Residues, determines the type of
  sorting that has to be done and the relevant contactDist (0 intra, 1 sequential, ... , -1 long range, -2 intermolecular).

  This is useful for distance constraint handling.
  """

  sortType = 'ccpCode'
  if residueKey1[0] == residueKey2[0]:
    seqDiff = abs(residueKey1[1] - residueKey2[1])
    if seqDiff <= maxRangeCheck:
      sortType = 'seqId'
      contactDist = seqDiff
    else:
      contactDist = -1
  else:
    contactDist = -2

  return (sortType,contactDist)

def getDescriptorDict(descriptor):

  """
  Creates a dictionary from a chemCompVar.descriptor
  """

  ddict = {}

  if descriptor and descriptor not in ('neutral','none'):

    comps = descriptor.split(';')
    for comp in comps:
      (key,atomNames) = comp.split(':')
      atomList = atomNames.split(',')
      ddict[key] = atomList

  return ddict

def mergeDescriptorDicts(mainDict,addDict):

  """
  Merges dictionaries describing chemCompVar.descriptors
  """

  for addKey in addDict.keys():
    
    # First make sure to remove atoms from existing info!
    for atomName in addDict[addKey]:
      for mainKey in mainDict.keys():
        if atomName in mainDict[mainKey]:
          mainDict[mainKey].pop(mainDict[mainKey].index(atomName))
    
    if mainDict.has_key(addKey):
      for atomName in addDict[addKey]:
        if atomName not in mainDict[addKey]:
          mainDict[addKey].append(atomName)
    else:
      mainDict[addKey] = addDict[addKey]

def getDescriptorFromDict(mainDict):

  """
  Creates a valid chemCompVar.descriptor from a dictionary describing it (see getDescriptorDict())
  """

  descKeys = []

  for tempKey in ['prot','deprot','link']:
    if tempKey in mainDict.keys():
      descKeys.append(tempKey)

  descStringList = []

  for descKey in descKeys:
    if mainDict[descKey]:
      mainDict[descKey].sort()
      descStringList.append("%s:%s" % (descKey,','.join(mainDict[descKey])))

  if descStringList:
    descriptor = ';'.join(descStringList)
  else:
    descriptor = 'neutral'

  return descriptor

def getDescriptorElementMatches(chemCompVars,descriptorElements,verbose = False):

  if descriptorElements and chemCompVars:
    matchingChemCompVars = []
    for tccv in chemCompVars:
      descriptorDict = getDescriptorDict(tccv.descriptor)
      allMatch = True
      for descriptorKey in descriptorElements.keys():
        if descriptorDict.has_key(descriptorKey):
          for atomName in descriptorElements[descriptorKey]:
            if atomName not in descriptorDict[descriptorKey]:
              allMatch = False
              break
        else:
          allMatch = False
          break

      if allMatch:
        matchingChemCompVars.append(tccv)

    if matchingChemCompVars:
      chemCompVars = matchingChemCompVars
    elif verbose:
      chemComp = chemCompVars[0].chemComp
      print "  Error: could not select chemCompVars from (%s,%s) with descriptor elements %s." % (chemComp.molType,chemComp.ccpCode,getDescriptorFromDict(descriptorElements))

  return chemCompVars

def createMoleculeTorsionDict(molSystem):

  """

  Function to get torsion angle definitions per molecule for a molecular system.

  Returns a dictionary organized in the following way (molecule, molRes are CCPN objects):

  moleculeTorsions[molecule][molRes] = [[torsion,[(molRes1,chemAtom1),(molRes2,chemAtom2),(molRes3,chemAtom3),(molRes4,chemAtom4)]],...]

  """

  moleculeTorsions = {}

  for chain in molSystem.chains:

    molecule = chain.molecule

    if not moleculeTorsions.has_key(molecule):

      moleculeTorsions[molecule] = {}

      for molRes in molecule.molResidues:

        moleculeTorsions[molecule][molRes] = []

        chemCompVar = molRes.chemCompVar

        for chemTorsion in chemCompVar.chemTorsions:

          moleculeTorsions[molecule][molRes].append([chemTorsion,[]])
          currentAtoms = moleculeTorsions[molecule][molRes][-1][-1]

          for chemAtom in chemTorsion.chemAtoms:

            if isinstance(chemAtom,ChemComp.ChemAtom):
              currentAtoms.append((molRes,chemAtom))

            else:

              if chemAtom.boundLinkEnd:
                linkCode = chemAtom.boundLinkEnd.linkCode
                distanceType = 'bound'
              else:
                linkCode = chemAtom.remoteLinkEnd.linkCode
                distanceType = 'remote'

              molResLinkEnd = molRes.findFirstMolResLinkEnd(linkCode = linkCode)
              molResLink = molResLinkEnd.molResLink

              if not molResLink:
                print "  Warning: could not create torsion information for %s.%s.%s - no linking to %s residue." % (molecule.name,molRes.seqCode,chemTorsion.name,linkCode)
                continue

              molResLink_LinkEnds = molResLink.sortedMolResLinkEnds()

              if molResLink_LinkEnds[0] == molResLinkEnd:
                otherMolResLinkEnd = molResLink_LinkEnds[1]
              else:
                otherMolResLinkEnd = molResLink_LinkEnds[0]

              otherMolRes = otherMolResLinkEnd.molResidue
              otherChemCompVar = otherMolRes.chemCompVar
              otherLinkEnd = otherChemCompVar.findFirstLinkEnd(linkCode = otherMolResLinkEnd.linkCode)

              otherChemAtom = getattr(otherLinkEnd,distanceType + 'ChemAtom')

              currentAtoms.append((otherMolRes,otherChemAtom))

  return moleculeTorsions

def getTorsions(structureGeneration,moleculeTorsions = None):

  """

  Function to get torsion angles from a structureGeneration.

  Returns a dictionary organized in the following way (structure, chain,
  residue are CCPN objects):

  torsionAngles[structure][chain][residue][torsionName] = torsionAngle

  """

  torsionAngles = {}
  
  structureEnsemble = structureGeneration.structureEnsemble
  if structureEnsemble:

    molSys = structureEnsemble.molSystem

    if not moleculeTorsions:
      moleculeTorsions = createMoleculeTorsionDict(molSys)

    for structure in structureEnsemble.models:

      torsionAngles[structure] = {}

      for coordChain in structure.coordChains:

        chain = coordChain.chain
        molecule = chain.molecule

        torsionAngles[structure][chain] = {}

        if moleculeTorsions.has_key(molecule):

          molResDict = {}

          for coordRes in coordChain.residues:
            molResDict[coordRes.residue.molResidue] = coordRes

          for coordRes in coordChain.residues:

            residue = coordRes.residue
            molRes = residue.molResidue

            torsionAngles[structure][chain][residue] = {}

            if moleculeTorsions[molecule].has_key(molRes):
              for (torsion,torsionChemAtoms) in moleculeTorsions[molecule][molRes]:

                torsionName = torsion.name

                torsCoords = []
                for (torsionMolRes,torsionChemAtom) in torsionChemAtoms:
                  if molResDict.has_key(torsionMolRes):
                    torsionCoordRes = molResDict[torsionMolRes]

                    torsionCoordAtom = torsionCoordRes.findFirstAtom(name = torsionChemAtom.name)
                    coords = torsionCoordAtom.coords[0]
                    torsCoords.append((coords.x,coords.y,coords.z))

                  else:
                    torsCoords.append(None)

                if None in torsCoords:
                  print "  Error: unable to calculate torsion angle %s (molSys '%s', chain '%s', seqCode %s)" % (torsionName,molSys.name,chain.code,residue.seqCode)
                  continue

                torsionAngle = calcTorsionAngleDegrees(torsCoords[0],torsCoords[1],torsCoords[2],torsCoords[3])

                if not torsionAngles[structure][chain][residue].has_key(torsionName):
                  torsionAngles[structure][chain][residue][torsionName] = torsionAngle
                else:
                  print "  Error: duplicate torsion angle %s (molSys '%s', chain '%s', seqCode %s). Set value %.2f, also found value %.2f." % (torsionName,molSys.name,chain.code,residue.seqCode,torsionAngles[structure][chain][residue][torsionName],torsionAngle)

  return torsionAngles                
                

def createOneLetterSequence(molType,sequence):
  
  """
  Creates a one-letter sequence string from a list of sequence labels.
  
  Input:  molType:  Molecular type
          sequence: List/tuple of (seqId,seqLabel) elements

  Output: String with one letter codes
  
  """

  if ccpCodeToCode1LetterDict.has_key(molType):
    shortSeq = ""
    for seqEl in sequence:
      shortCode = 'X'
      if seqEl[1]:
        ccpCode = seqEl[1].capitalize()
        if ccpCodeToCode1LetterDict[molType].has_key(ccpCode):
          shortCode = ccpCodeToCode1LetterDict[molType][ccpCode]

      shortSeq += shortCode

  else:
    shortSeq = 'X' * len(sequence)
    
  return shortSeq


'''
#
# ReferenceShifts class creates a dictionary of random coil shift values based on a molecule
#

class ReferenceShifts:

  """
  Can be used to create a dictionary of random coil shift values based on a molecule

  Use as:

  refShifts = ReferenceShifts(molecule, nmrRefStorePath = refPath) # If nmrRefStorePath not set, will use default
  valueDict = refShifts.calculateAll()

  """

  def __init__(self, molecule, nmrRefStorePath = None):
    
    from ccp.general.Io import getNmrReferenceStoragePath
  
    self.molecule = molecule
    self.sourceDict = {'pescador': ('mean','corrected'),
                       'merutka':  ('mean',)}

    if not nmrRefStorePath:
      self.nmrRefStorePath = getNmrReferenceStoragePath()
    else:
      self.nmrRefStorePath = nmrRefStorePath

    self.initialize()

  def initialize(self):

  
    self.nmrRefStores = {}
    self.moleculeCcpCodes = []
    self.chemCompNmrRefDict = {}

    for molResidue in self.molecule.molResidues:
      self.moleculeCcpCodes.append(molResidue.ccpCode)
      molResKey = self.getMolResKey(molResidue)
      if not self.nmrRefStores.has_key(molResKey):
        nmrRefStore = self.molecule.root.findFirstNmrReferenceStore(molType=molResKey[0],ccpCode=molResKey[1])        
        self.nmrRefStores[molResKey] = nmrRefStore

  def setupChemCompNmrRefDict(self,sourceName):

    if not self.chemCompNmrRefDict.has_key(sourceName):
      self.chemCompNmrRefDict[sourceName] = {}
      for molResKey in self.nmrRefStores.keys():
        self.chemCompNmrRefDict[sourceName][molResKey] = self.nmrRefStores[molResKey].findFirstChemCompNmrRef(sourceName = sourceName)

  def calculateAll(self):

    valueDict = {}

    sourceList = self.sourceDict.keys()
    sourceList.sort()

    for sourceName in sourceList:
      self.setupChemCompNmrRefDict(sourceName)
      valueDict[sourceName] = {}
      for calcType in self.sourceDict[sourceName]:
        valueDict[sourceName][calcType] = []
        for molResidue in self.molecule.molResidues:
          valueDict[sourceName][calcType].append([molResidue.seqId,[]])
          for chemAtom in molResidue.chemCompVar.chemAtoms:
            shiftValue = self.getMeanShiftValue(molResidue,sourceName,calcType,chemAtom.name)
            if shiftValue != None:
              valueDict[sourceName][calcType][-1][-1].append((chemAtom.name,shiftValue))

    return valueDict

  def getMolResKey(self,molResidue):

    return (molResidue.molType,molResidue.ccpCode)

  def getChemAtomShifts(self,molResidue,chemCompNmrRef,atomName):

    chemCompVar = molResidue.chemCompVar
    chemAtomShifts = chemCompNmrRef.findAllChemAtomShifts(chemCompLinking = chemCompVar.linking,
                                                          chemCompDescriptor = chemCompVar.descriptor,
                                                          atomName = atomName)

    if not chemAtomShifts:
      chemAtomShifts = chemCompNmrRef.findAllChemAtomShifts(chemCompLinking = chemCompVar.linking,
                                                            chemCompDescriptor = 'any',
                                                            atomName = atomName)

      if not chemAtomShifts:
        chemAtomShifts = chemCompNmrRef.findAllChemAtomShifts(chemCompLinking = 'any',
                                                              chemCompDescriptor = 'any',
                                                              atomName = atomName)
    return chemAtomShifts

  def getMeanShiftValue(self,molResidue,sourceName,calcType,atomName):

    if sourceName not in self.sourceDict:
      print "  Error: unknown source %s!" % sourceName
      return

    shiftValue = None

    molResKey = self.getMolResKey(molResidue)
    chemCompNmrRef = self.chemCompNmrRefDict[sourceName][molResKey]

    if chemCompNmrRef:

      chemAtomShifts = self.getChemAtomShifts(molResidue,chemCompNmrRef,atomName)

      if chemAtomShifts:

        if len(chemAtomShifts) > 1:
          print "  Error: multiple reference shifts available... using first one."

        chemAtomShift = chemAtomShifts[0]

        #
        # Pescador data
        #

        if sourceName == 'pescador':

          shiftValue = chemAtomShift.meanValue

          if calcType == 'corrected':

            i = molResidue.seqId

            for j in [-2,-1,1,2]:

              if 0 < (j+i) <= len(self.moleculeCcpCodes):

                chemAtomShiftCorr = chemAtomShift.findFirstChemAtomShiftCorr(seqOffset = j, ccpCode = self.moleculeCcpCodes[i+j-1])
                shiftValue += chemAtomShiftCorr.value

        #
        # Merutka data
        #

        elif sourceName == 'merutka':

          shiftValue = chemAtomShift.randomCoilValue

    return shiftValue
'''

#
# ChemComp stuff
#

def getAtomPatterns(chemAtomList):

  """
  Returns regular expression dict for recognizing atom groups
  """

  atomPatterns = {}

  # Do all chemAtoms for chemAtomList
  for chemAtom in chemAtomList:
    chemAtomName = chemAtom.name

    # Assuming IUPAC name: only using names that have a digit
    # at the end.
    if chemAtomName[-1] in '0123456789':

      found = 0

      for atomPattern in atomPatterns:
        if atomPattern.search(chemAtomName):
          # Match
          found = 1
          atomPatterns[atomPattern].append(chemAtom)
          break

      if found == 0:
        try:
          atomPattern = re.compile("^" + chemAtomName[0:-1] + "\d$")
          atomPatterns[atomPattern] = [chemAtom]
        except:
          print "ERROR: bad atom name '%s'" % chemAtomName
          #raise

  return atomPatterns

class MakeAtomSets:

  """

  Class to make automatically make AtomSets for a new chemComp

  """

  def checkAtomPattern(self,atom,foundAtomPattern):

    for atomPattern in self.atomPatterns:
      if atomPattern.search(atom.name):
        if not foundAtomPattern.has_key(atomPattern):
          foundAtomPattern[atomPattern] = [atom]
        else:
          foundAtomPattern[atomPattern].append(atom)

  def setChemComps(self,chemComp,createSysNames = True,namingSystem=None,verbose = False):

    if verbose:
      from memops.universal.Util import drawBox
      print drawBox("Creating atom sets",indent = "  ")

    self.chemComp = chemComp
    self.namingSystem = namingSystem
    self.atomSetsCreated = []
    self.atomsCreated = []
    self.errorMessages = []

    self.createSysNames = createSysNames

    for chemCompVar in self.chemComp.chemCompVars:

      self.setChemCompVar(chemCompVar)
    
    for errorMessage in self.errorMessages:
      print errorMessage
    print


    if self.atomsCreated:
      self.atomsCreated.sort()
      print "    Created chemAtoms: %s" % str(self.atomsCreated)
        
    if self.atomSetsCreated:
      self.atomSetsCreated.sort()
      print "    Created chemAtomSets: %s" % str(self.atomSetsCreated)

  def setChemCompVar(self,chemCompVar):

    # See if IUPAC exists, use that if possible... otherwise default to MSD
    if self.namingSystem:
      namingSystem = self.namingSystem
    else:
      for namingSystemName in ('IUPAC','MSD'):
        namingSystem = self.chemComp.findFirstNamingSystem(name = namingSystemName)
        if namingSystem:
          break

    if not namingSystem:
      namingSystem = self.chemComp.newNamingSystem(name = 'IUPAC')
      
    # Alternate set names - Q/M naming
    alternateNamingSystem = self.chemComp.findFirstNamingSystem(name = 'AQUA')
    
    if not alternateNamingSystem:
      alternateNamingSystem = self.chemComp.newNamingSystem(name = 'AQUA')
      
    #
    # First look for atom name patterns - this code can be useful in naming the chemAtomSet later on
    # Because names are often not consistent, it cannot be used to find the chemAtomSets though.
    #

    self.atomPatterns = getAtomPatterns(chemCompVar.chemAtoms)
    self.atomPatternList = self.atomPatterns.keys()
    self.atomPatternList.sort()

    #
    # Now check which atoms (of the same element) are connected to each other via ONLY ONE other atom.
    #
    
    atomConnects = {}
    
    # TODO should this happen on a chemCompVar level??!?! Probably yes? Not trivial without
    # certain amount of curation...
    atomList = self.chemComp.sortedChemAtoms()
    #atomConnectSet = {}

    for i in range(0,len(atomList)-1):
      chemAtom1 = atomList[i]
    
      #
      # Get atom at other end of bond
      #
        
      for bond1 in chemAtom1.sortedChemBonds():
        chemAtomCentral = getOtherAtom(chemAtom1,bond1)
    
        #if not atomConnectSet.has_key(chemAtomCentral):
        # atomConnectSet[chemAtomCentral] = False
          
        #
        # Loop over other atoms connected to this chemAtomCentral
        #
        
        for centralBond in chemAtomCentral.sortedChemBonds():
          chemAtom2 = getOtherAtom(chemAtomCentral,centralBond)
          
          if chemAtom1 == chemAtom2:
            continue
               
          # The joint sets have to be of the the same element and subType! And in the same list..
          if chemAtom2.elementSymbol != chemAtom1.elementSymbol or chemAtom2.subType != chemAtom1.subType:
            continue
    
          if not atomConnects.has_key(chemAtomCentral):
            atomConnects[chemAtomCentral] = []
          
          setOfChemAtoms = set((chemAtom1,chemAtom2))
          
          setAvailable = False
          for existingSetOfChemAtoms in atomConnects[chemAtomCentral]:
            if setOfChemAtoms in existingSetOfChemAtoms:
              setAvailable = True
              break
          
          if not setAvailable: 
            
            tmpChemAtomList = list(setOfChemAtoms)
            
            chemAtomAdded = False
            for k in range(len(tmpChemAtomList)):
              chemAtom = tmpChemAtomList[k]
              for existingSetOfChemAtoms in atomConnects[chemAtomCentral]:
                if chemAtom in existingSetOfChemAtoms:
                  chemAtomAdded = True
                  existingSetOfChemAtoms.add(tmpChemAtomList[not tmpChemAtomList.index(chemAtom)])

            if not chemAtomAdded:
              atomConnects[chemAtomCentral].append(setOfChemAtoms)
    
    #
    # Make a sorted list so always handles in same order
    #

    chemAtomCentralNameList = []
    chemAtomCentralNames = {}

    for chemAtomCentral in atomConnects.keys():

      chemAtomCentralNameList.append(chemAtomCentral.name)
      chemAtomCentralNames[chemAtomCentral.name] = chemAtomCentral

    chemAtomCentralNameList.sort()


    #
    # Get a list of the bonds relevant for this chemCompVar
    #

    chemCompVarBonds = chemCompVar.chemBonds

    # Loop over connected ones
    #
    # Check if part of a single atomPattern. If not, check if at least two match
    # with one (and use that one then)
    #
    # If connected ones are protons:
    #   - if two connected ones, and four bonds to central atom, are prochiral
    #   - if two connected ones, and three bonds to central atom, are not equivalent or prochiral
    #   - if three connected ones, are equivalent

    for chemAtomCentralName in chemAtomCentralNameList:

      chemAtomCentral = chemAtomCentralNames[chemAtomCentralName]
      
      #print
      #print "CENTRAL ATOM:",chemAtomCentral.name
      
      #
      # Re-sort so that starts with largest chemAtomLists. Better for consistency
      #
      
      setOfChemAtomsLengths = {}
      for setOfChemAtoms in atomConnects[chemAtomCentral]:
        setOfChemAtomsLen = len(setOfChemAtoms)
        if not setOfChemAtomsLen in setOfChemAtomsLengths.keys():
          setOfChemAtomsLengths[setOfChemAtomsLen] = []
          
        if not setOfChemAtoms in setOfChemAtomsLengths[setOfChemAtomsLen]:
          setOfChemAtomsLengths[setOfChemAtomsLen].append(setOfChemAtoms)
        
      sortedSetOfChemAtomsLengths = setOfChemAtomsLengths.keys()
      sortedSetOfChemAtomsLengths.sort()
      sortedSetOfChemAtomsLengths.reverse()
      
      sortedChemAtomTuples = []
      for setOfChemAtomsLen in sortedSetOfChemAtomsLengths:
        for setOfChemAtoms in setOfChemAtomsLengths[setOfChemAtomsLen]:
          sortedChemAtomTuples.append(tuple(setOfChemAtoms))
                
      #
      # Now go over lists
      #

      for chemAtomTuple in sortedChemAtomTuples:
        #print "  START SET:", ', '.join(["%s.%d" % (ca.name,ca.subType) for ca in chemAtomTuple])
        numConnect = 0
        for chemBond in chemAtomCentral.chemBonds:
          if chemBond in chemCompVarBonds:
            numConnect += 1

        foundAtomPattern = {}

        # Loop over connected atoms
        for chemAtom in chemAtomTuple:
          # Check if matches atom Pattern
          self.checkAtomPattern(chemAtom,foundAtomPattern)

        # No matches with any pattern
        if len(foundAtomPattern) == 0:
          message = "      ERROR: no matches for names %s." % [ca.name for ca in chemAtomTuple]
          if not message in self.errorMessages:
            self.errorMessages.append(message)

        # Matches multiple patterns: pick highest frequency one
        elif len(foundAtomPattern) > 1:
          highNum = 0
          for atomPattern in foundAtomPattern.keys():
            if len(foundAtomPattern[atomPattern]) > highNum:
              chemAtomTuple = foundAtomPattern[atomPattern]
              highNum = len(foundAtomPattern[atomPattern])

        # Protons
        if chemAtomTuple[0].elementSymbol == 'H':

          chemAtomSetName = chemAtomTuple[0].elementSymbol + chemAtomCentralName[1:] + '*'
          
          if chemCompVar.chemComp.findFirstChemAtom(name = chemAtomSetName):
            print "      ERROR: Unable to create ChemAtomSet, atom with name '%s' already exists!" % chemAtomSetName
            continue

          isProchiral = False

          if len(chemAtomTuple) == 2:

            isEquivalent = False

            if numConnect == 4:
              # Prochiral set
              isProchiral = True

            elif numConnect == 3 or numConnect == 2:
              # Not equivalent, not prochiral set
              pass

            else:
              print "  Error invalid number of connecting bonds (%d) to central atom %s" % (numConnect,chemAtomCentral.name)

          elif len(chemAtomTuple) == 3:
            # Equivalent set
            isEquivalent = True

          else:
            print "  Error invalid number of matched protons."
            isEquivalent = False

          #
          # See if the chemAtomSet already exists
          #

          chemAtomSet = self.chemComp.findFirstChemAtomSet(name = chemAtomSetName, chemAtoms = set(chemAtomTuple))
          #print "    PROTONS: " + ', '.join(["%s.%d" % (ca.name,ca.subType) for ca in chemAtomTuple]) + " CAS:" + str(chemAtomSet)
          if not chemAtomSet:
              
            existingChemAtomSet = self.chemComp.findFirstChemAtomSet(name=chemAtomSetName)
            if existingChemAtomSet:
              isSubSet = False
              if set(chemAtomTuple) in existingChemAtomSet.chemAtoms or existingChemAtomSet.chemAtoms in set(chemAtomTuple):
                isSubSet = True
              
              if not isSubSet:
                print "    ERROR: invalid new chemAtomSet being created, %s not a subset of %s - cowardly refusing this." % ([ca.name for ca in chemAtomTuple],[ca.name for ca in existingChemAtomSet.sortedChemAtoms()])
                continue

            #print "      NEW CAS %s:" % chemAtomSetName, chemCompVar.linking, chemCompVar.descriptor, ', '.join(["%s.%d" % (ca.name,ca.subType) for ca in chemAtomTuple])
 
            subType = 0
            
            for tmpChemAtomSet in self.chemComp.findAllChemAtomSets(name = chemAtomSetName):
              if set([ca.name for ca in chemAtomTuple]) == set([ca.name for ca in tmpChemAtomSet.chemAtoms]):
                subType = tmpChemAtomSet.subType - 1
              elif tmpChemAtomSet.subType > subType:
                subType = tmpChemAtomSet.subType

            subType += 1

            #
            # Create new chemAtoms for new subType...
            #

            if subType > 1:

              #print "    Creating new subtypes"

              # This list includes the current chemCompVar!
              changeChemCompVars = []
              for otherChemCompVar in self.chemComp.chemCompVars:
                makeChange = True

                # First check if current chemCompVar has all the relevant atoms
                for chemAtom in chemAtomTuple:
                  if not chemAtom in otherChemCompVar.chemAtoms:
                    makeChange = False
                    break

                # Then check if, if these chemAtoms are part of a chemAtomSet, this chemAtomSet is NOT valid for this chemCompVar
                if chemAtom.chemAtomSet:
                  if chemAtom.chemAtomSet in otherChemCompVar.chemAtomSets:
                    makeChange = False

                if makeChange:
                  changeChemCompVars.append(otherChemCompVar)

              newChemAtomList = []

              for chemAtom in chemAtomTuple:
                newChemAtom = self.chemComp.findFirstChemAtom(name = chemAtom.name, subType = subType)
                if not newChemAtom:
                  newChemAtom = self.chemComp.newChemAtom(name = chemAtom.name, subType = subType, elementSymbol = chemAtom.elementSymbol, waterExchangeable = chemAtom.waterExchangeable)
                newChemAtomList.append(newChemAtom)
                
                self.atomsCreated.append((newChemAtom.name,subType))

                for chemBond in chemAtom.chemBonds:
                  otherChemAtom = getOtherAtom(chemAtom,chemBond)
                  
                  existingChemBond = self.chemComp.findFirstChemBond(chemAtoms = set((otherChemAtom,newChemAtom)))
                  if not existingChemBond:
                    self.chemComp.newChemBond(chemAtoms = [otherChemAtom,newChemAtom], bondType = chemBond.bondType)

                for changeChemCompVar in changeChemCompVars:
                  changeChemCompVar.removeChemAtom(chemAtom)
                  changeChemCompVar.addChemAtom(newChemAtom)

            else:
              newChemAtomList = chemAtomTuple

            chemAtomSet = self.chemComp.findFirstChemAtomSet(name = chemAtomSetName, subType = subType)
            if not chemAtomSet:
              #print newChemAtomList
              chemAtomSet = self.chemComp.newChemAtomSet(name = chemAtomSetName, subType = subType, isEquivalent = isEquivalent, isProchiral = isProchiral, chemAtoms = newChemAtomList)

            if chemAtomSetName not in self.atomSetsCreated:
              self.atomSetsCreated.append((chemAtomSetName,subType))

            if self.createSysNames:

              keywds = {'atomName': chemAtomSetName,
                        'atomSubType': chemAtomSet.subType,
                        'sysName': chemAtomSetName}

              if not namingSystem.findFirstAtomSysName(**keywds):
                chemAtomSetSysName = namingSystem.newAtomSysName(**keywds)
              
                self.setAlternateSysName(alternateNamingSystem,newChemAtomList,keywds)
              
          else:
            #print "  ERROR: double atomSet name for %s... ignored" % chemAtomSetName
            pass

    # TODO: CRAP code... can't be bothered to make it better
    # hetero's have to be done after protons 'cos otherwise code
    # for atomset linking won't work

    # If connected ones are NOT protons:
    #   - Forget equivalent ones: problems with Trp for example:
    #       HH2 to HZ3,HZ2, ... is very similar to Tyr, Phe...
    #   - if two connected ones, and four bonds to central atom, are considered prochiral

    for chemAtomCentralName in chemAtomCentralNameList:

      chemAtomCentral = chemAtomCentralNames[chemAtomCentralName]

      numConnect = len(chemAtomCentral.chemBonds)

      for setOfChemAtoms in atomConnects[chemAtomCentral]:
        
        chemAtomList = list(setOfChemAtoms)

        foundAtomPattern = {}

        # Loop over connected atoms
        for chemAtom in chemAtomList:
          # Check if matches atom Pattern
          self.checkAtomPattern(chemAtom,foundAtomPattern)

        # No matches with any pattern
        if len(foundAtomPattern) == 0:
          message = "      ERROR: no matches for set of atom names %s." % [ca.name for ca in chemAtomList]
          if not message in self.errorMessages:
            self.errorMessages.append(message)

        # Matches multiple patterns: pick highest frequency one
        elif len(foundAtomPattern) > 1:
          highNum = 0
          for atomPattern in foundAtomPattern.keys():
            if len(foundAtomPattern[atomPattern]) > highNum:
              chemAtomList = foundAtomPattern[atomPattern]
              highNum = len(foundAtomPattern[atomPattern])

        # Heteroatoms
        if chemAtomList[0].elementSymbol != 'H':

          # Only try to find prochiral ones - forget about equivalent
          # (gets too complicated)

          chemAtomSetName = chemAtomList[0].name[:-1] + '*'

          if len(chemAtomList) == 2:
            if numConnect == 4:
              # Do some additional checking (protons only)
              chemAtomProtons = {}

              noProtons = 0

              # Get list of all connected protons
              for chemAtom in chemAtomList:
                chemAtomProtons[chemAtom] = []

                # Get connected protons
                for bond in chemAtom.chemBonds:

                  otherAtom = getOtherAtom(chemAtom,bond)

                  # Only get protons
                  if otherAtom.elementSymbol == 'H':
                    chemAtomProtons[chemAtom].append(otherAtom)

                # If any has no protons, just ignore this list...
                if not chemAtomProtons[chemAtom]:
                  noProtons = 1

              # If any has no protons, just ignore this chemAtomCentral...
              if noProtons == 1:
                continue

              else:
                # Set a new pattern and check if all protons match
                atomPattern = re.compile("^" + chemAtomProtons[chemAtom][0].name[0:-2] + "\d\d$")
                protonMatch = 0
                protonAll = 0
                chemAtomProtonSetNames = {}

                for chemAtom in chemAtomList:
                  for chemAtomProton in chemAtomProtons[chemAtom]:
                    protonAll += 1

                    if atomPattern.search(chemAtomProton.name):
                      protonMatch += 1
                      # Set pattern match for this one...
                      chemAtomProtonSetNames[chemAtomProton.name[:-1] + '*'] = 1

                # All matches, check if proton set codes exist and go ahead with set
                if protonMatch == protonAll:

                  noProtonSet = 0
                  protonSetNumber = -1

                  for chemAtomProtonSetName in chemAtomProtonSetNames.keys():
                    chemAtomProtonSet = self.chemComp.findFirstChemAtomSet(name = chemAtomProtonSetName)

                    if not chemAtomProtonSet:
                      noProtonSet = 1

                    else:

                      if protonSetNumber == -1:
                        protonSetNumber = len(chemAtomProtonSet.chemAtoms)

                      # Simple check for prochirality:
                      # if number of protons on groups is not the same then not prochiral
                      # At least takes care of Ile...
                      if protonSetNumber != len(chemAtomProtonSet.chemAtoms):
                        noProtonSet = 1

                  if noProtonSet != 1:

                    # Prochiral set
                    isEquivalent = False
                    isProchiral = True

                    # Set heteroAtom prochiral (eg CG*)

                    #print "      NEW CAS (1) %s:" % chemAtomSetName, chemCompVar.linking, chemCompVar.descriptor, ', '.join(["%s.%d" % (ca.name,ca.subType) for ca in chemAtomList])

                    chemAtomSet = self.chemComp.findFirstChemAtomSet(name = chemAtomSetName)
                    if not chemAtomSet:
                      chemAtomSet = self.chemComp.newChemAtomSet(name = chemAtomSetName, isEquivalent = isEquivalent, isProchiral = isProchiral, chemAtoms = chemAtomList)

                      if chemAtomSetName not in self.atomSetsCreated:
                        self.atomSetsCreated.append((chemAtomSetName,chemAtomSet.subType))

                      if self.createSysNames:

                        keywds = {'atomName': chemAtomSetName,
                                  'atomSubType': chemAtomSet.subType,
                                  'sysName': chemAtomSetName}

                        if not namingSystem.findFirstAtomSysName(**keywds):
                          chemAtomSetSysName = namingSystem.newAtomSysName(**keywds)
    
                          self.setAlternateSysName(alternateNamingSystem,chemAtomList,keywds)

                    # Set proton pseudoatom (eg HG*)
                    chemAtomProtonSets = []
                    for chemAtomProtonSetName in chemAtomProtonSetNames.keys():
                      chemAtomProtonSet = self.chemComp.findFirstChemAtomSet(name = chemAtomProtonSetName)
                      chemAtomProtonSets.append(chemAtomProtonSet)

                    chemAtomSetName = chemAtomProtonSet.name[0:-2] + "*"

                    #
                    # Make sure that the atomSet used has the highest number of chemAtoms
                    #

                    chemAtomSet = self.chemComp.findFirstChemAtomSet(name = chemAtomSetName)

                    #print "      NEW CAS (2) %s:" % chemAtomSetName, chemCompVar.linking, chemCompVar.descriptor, ', '.join(["%s.%d" % (cas.name,cas.subType) for cas in chemAtomProtonSets])

                    if chemAtomSet and not chemAtomSet.chemAtoms and len(chemAtomSet.chemAtomSets) < len(chemAtomProtonSets):
                      print "  Warning: deleting chemAtomSet %s with less atomSets (%d) than new chemAtomSet (%d)" % (chemAtomSetName,len(chemAtomSet.chemAtomSets),len(chemAtomProtonSets))
                      chemAtomSet.delete()

                      if len(chemAtomProtonSets) > 2:
                        isProchiral = False

                    #
                    # Now go ahead...
                    #

                    if not self.chemComp.findFirstChemAtomSet(name = chemAtomSetName):

                      if len(chemAtomProtonSets) >= 2:
                        chemAtomSet = self.chemComp.newChemAtomSet(name = chemAtomSetName, isEquivalent = isEquivalent, isProchiral = isProchiral, chemAtomSets = chemAtomProtonSets)

                        if chemAtomSetName not in self.atomSetsCreated:
                          self.atomSetsCreated.append((chemAtomSetName,chemAtomSet.subType))

                        if self.createSysNames:

                          keywds = {'atomName': chemAtomSetName,
                                    'atomSubType': chemAtomSet.subType,
                                    'sysName': chemAtomSetName}

                          if not namingSystem.findFirstAtomSysName(**keywds):
                            chemAtomSetSysName = namingSystem.newAtomSysName(**keywds)
                            
                            chemAtomList = []
                            for chemAtomProtonSet in chemAtomProtonSets:
                              chemAtomList.extend(chemAtomProtonSet.sortedChemAtoms())

                            self.setAlternateSysName(alternateNamingSystem,chemAtomList,keywds)

                      else:
                        print "  Warning: invalid prochiral center %s! Has %d protons... not created." % (chemAtomSetName,len(chemAtomProtonSets))

                    else:
                      pass
                      #print "  Error: double chemAtomSet name %s... ignored" % chemAtomSetName

  def setAlternateSysName(self,alternateNamingSystem,chemAtomList,keywds):
    
    #print chemAtomList,keywds
    #print
    #
    # Set AQUA names, only for H connected to C and C.
    #

    chemAtomSetName = keywds['atomName']
    setAlternateName = False

    if chemAtomSetName[0] == 'H':
      refChemAtom = chemAtomList[0]
      chemBond = refChemAtom.findFirstChemBond()
      centralAtom = getOtherAtom(refChemAtom,chemBond)
      if centralAtom.elementSymbol == 'C':
        setAlternateName = True
    elif chemAtomSetName[0] == 'C':
      setAlternateName = True
                
    if setAlternateName:
      if len(chemAtomList) == 2:
        nameType = 'Q'
      elif len(chemAtomList) == 3:
        nameType = 'M'
      elif len(chemAtomList) == 6:
        nameType = 'QQ'
      else:
        nameType = None
        
      if nameType:
        keywds['sysName'] = nameType + chemAtomSetName[1:-1]
         
        if not alternateNamingSystem.findFirstAtomSysName(**keywds):
          alternateChemAtomSysName = alternateNamingSystem.newAtomSysName(**keywds)

def createNewDescriptors(chemComp,newDict,removeAtomNames,renameAtoms,saveChemComp=True):

  """

  Creates a new set of descriptors for a particular chemComp.

  chemComp is the chemComp to modify

  newDict is the new descriptor (e.g. {'link': ['ND2']})

  removeAtomNames is a list of atom names that have to be removed for this new descriptor

  renameAtoms is a dictionary with atoms that have to be renamed (e.g. {'HD22': 'HD2'})

  """

  chemComp.setIsModifiable(True)

  molType = chemComp.molType
  ccpCode = chemComp.ccpCode

  allCcvs = chemComp.chemCompVars

  ccvs = getDescriptorElementMatches(allCcvs,newDict)

  if ccvs != allCcvs:
    print "  Exiting: chemCompVars already exist with '%s'" % (getDescriptorFromDict(newDict))
    sys.exit()

  #
  # Check whether atoms exist and get a list of atoms that are bound to the ones for the new code...
  #

  for removeAtomName in removeAtomNames:
    ca = chemComp.findFirstChemAtom(name = removeAtomName)
    if not ca:
      removeAtomNames.pop(removeAtomNames.index(removeAtomName))
      print "  Removing chemAtom %s from 'remove' list - does not exist." % removeAtomName

  boundChemAtomNames = {}
  linkEnd = None

  for ndkey in newDict.keys():
    atomNames = newDict[ndkey]
    for atomName in atomNames:
      linkingChemAtom = chemComp.findFirstChemAtom(name = atomName)

      if not linkingChemAtom:
        print "  Removing chemAtom %s - does not exist." % atomName
        newDict[ndkey].pop(newDict[ndkey].index(atomName))

      else:
        if not boundChemAtomNames.has_key(atomName):
          boundChemAtomNames[atomName] = []

        for bond in linkingChemAtom.chemBonds:
          otherChemAtom = getOtherAtom(linkingChemAtom,bond)
          boundChemAtomNames[atomName].append(otherChemAtom.name)

    if not newDict[ndkey]:
      del(newDict[ndkey])

    elif ndkey == 'link':

      #
      # Create linkEnds
      #

      linkedAtomName = linkingChemAtom.name

      linkEnd = chemComp.findFirstLinkEnd(linkCode = linkedAtomName)

      if not linkEnd:
        boundChemAtom = linkingChemAtom
        boundLinkAtom = getLinkAtom(boundChemAtom)
        linkEnd = ChemComp.LinkEnd(chemComp,linkCode = linkedAtomName, boundChemAtom = boundChemAtom, boundLinkAtom = boundLinkAtom)
        print "   Made linkEnd %s" % linkedAtomName



  if not newDict:
    print "  Exiting: no linking elements left...  ."
    sys.exit()

  #
  # Decide how to set link: info
  #

  if chemComp.isLinearPolymer:
    linksInDescriptor = True
  else:
    linksInDescriptor = False

  #
  # Now start the actual merging...
  #

  newCodeList = []
  ccvBaseList = {}
  ignoredCcvList = []

  for ccv in allCcvs:

    linkCodeText = "'%s','%s'" % (ccv.linking,ccv.descriptor)
    infoTexts = []

    newCodes = {'descriptor': copy.deepcopy(newDict), 'linking': {}}

    if linksInDescriptor:
      linkingDict = None
    else:
      # TODO: what is right way to handle this? link: goes in descriptor, and linking only used for protein/DNA/RNA, not ever for
      # anything else?
      linkingDict = getDescriptorDict(ccv.linking)
      if newCodes['descriptor'].has_key('link'):
        newCodes['linking']['link'] = newCodes['descriptor']['link'][:]
        del(newCodes['descriptor']['link'])

    descDict = getDescriptorDict(ccv.descriptor)

    #
    # If this ccv has prot: for one of the 'removeatoms', then don't bother!!!
    #

    if removeAtomNames and descDict.has_key('prot'):
      protonatedByRemoveAtom = False
      for removeAtomName in removeAtomNames:
        if removeAtomName in descDict['prot']:
          protonatedByRemoveAtom = True
          break
      if protonatedByRemoveAtom:
        ignoredCcvList.append(linkCodeText)
        continue

    #
    # For non-linear polymers, remove any new link: codes that don't make sense for this one.
    # For linear polymers, do the same but use the descDict
    #

    changeDict = None
    newCodesDict = None

    if linkingDict and newCodes['linking']:
      cleanupDictByLink(ccv,linkingDict,newCodes['linking'],boundChemAtomNames)
    elif not linkingDict and newCodes['descriptor'].has_key('link'):
      cleanupDictByLink(ccv,descDict,newCodes['descriptor'],boundChemAtomNames)


    #
    # If this ccv has deprot for the removeAtoms, then get rid of that in descDict...
    #

    if removeAtomNames and descDict.has_key('deprot'):
      for removeAtomName in removeAtomNames:
        if removeAtomName in descDict['deprot']:
          descDict['deprot'].pop(descDict['deprot'].index(removeAtomName))
      if not descDict['deprot']:
        del(descDict['deprot'])

    #
    # Add descDict and linkingDict to what's now the new code...
    #

    if linkingDict:
      mergeDescriptorDicts(newCodes['linking'],linkingDict)
    elif not linkingDict and not newCodes['linking']:
      newCodes['linking'][ccv.linking] = None

    if descDict:
      mergeDescriptorDicts(newCodes['descriptor'],descDict)

    #
    # Keep track of relevant chemAtoms for this one...
    #

    calist = list(ccv.chemAtoms)

    heavyAtomsNeedRemoving = False
    heavyAtomRemoved = False
    for removeAtomName in removeAtomNames:
    
      if removeAtomName[0] != 'H':
        heavyAtomsNeedRemoving = True
    
      ca = ccv.findFirstChemAtom(name = removeAtomName)
      if ca:
        calist.pop(calist.index(ca))
        
        if ca.elementSymbol != 'H':
          heavyAtomRemoved = True
        
        infoTexts.append("Removed atom %s" % ca.name)
        
    # Here check if  heavy atoms needed removing; only valid if they were!
    if heavyAtomsNeedRemoving and not heavyAtomRemoved:
      ignoredCcvList.append("%s, aborting creation, heavy atoms not removed" % linkCodeText)
      continue

    for renameAtomName in renameAtoms.keys():
      ca = ccv.findFirstChemAtom(name = renameAtomName)
      if ca:
        calist.pop(calist.index(ca))
        infoTexts.append("Removed atom %s" % ca.name)
        cbonds = ca.chemBonds
        newAtomName = renameAtoms[renameAtomName]
        newCa = ccv.parent.findFirstChemAtom(name = newAtomName)

        if not newCa:
          newCa = ChemComp.ChemAtom(ccv.parent,name = newAtomName, elementSymbol = ca.elementSymbol, chirality = ca.chirality, waterNmrObs = ca.waterNmrObs)
          for cb in cbonds:
            catoms = list(cb.chemAtoms)
            caindex = catoms.index(ca)
            catoms.pop(caindex)
            catoms.insert(caindex,newCa)
            ChemComp.ChemBond(ccv.parent,chemAtoms = catoms, bondType = cb.bondType)
          infoTexts.append("Created atom %s" % newCa.name)

        infoTexts.append("Added atom %s" % newCa.name)
        calist.append(newCa)

    #
    # Keep track of final results
    #

    ccvBaseList[linkCodeText] = [newCodes,infoTexts]
    newCodeList.append((newCodes,calist,ccv))

    #
    # TODO TODO CHEMCOMPCOORDS!!!!!  WILL HAVE TO DO SEPARATELY FOR protein,ASN!!
    #

  #
  # If nothing found, then just make a new ccv based on input
  #

  if not newCodeList:
    newCodes = {'descriptor': newDict.copy(), 'linking': {}}

    if linksInDescriptor:
      linkingDict = None
    else:
      if newCodes['descriptor'].has_key('link'):
        newCodes['linking']['link'] = newCodes['descriptor']['link'][:]
        del(newCodes['descriptor']['link'])

    # TODO TODO ADD VALID ATOM LIST HERE!!

    newCodeList.append(newCodes)


  #
  # Print info for new ones
  #

  ignoredCcvList.sort()

  print "Ignored ccv with linking, descriptor:"
  for lct in ignoredCcvList:
    print "   %s" % lct

  linkCodeTexts = ccvBaseList.keys()
  linkCodeTexts.sort()

  print
  print "Creating ccv with linking, descriptor based on original:"
  for lct in linkCodeTexts:
    linkDict = ccvBaseList[lct][0]['linking']
    descDict = ccvBaseList[lct][0]['descriptor']

    if len(linkDict.keys()) == 1 and not linkDict[linkDict.keys()[0]]:
      linking = linkDict.keys()[0]
    else:
      linking = getDescriptorFromDict(linkDict)

    descriptor = getDescriptorFromDict(descDict)

    print "   '%s','%s' (based on %s)" % (linking,descriptor,lct)

    for infoText in ccvBaseList[lct][1]:
      print "     %s" % infoText

  #
  # Then make them
  #

  for (newCodes,calist,refCcv) in newCodeList:

    if linkEnd:
      calist.append(linkEnd.boundLinkAtom)

    linkDict = newCodes['linking']
    descDict = newCodes['descriptor']

    if len(linkDict.keys()) == 1 and not linkDict[linkDict.keys()[0]]:
      linking = linkDict.keys()[0]
    else:
      linking = getDescriptorFromDict(linkDict)

    descriptor = getDescriptorFromDict(descDict)

    ccv = ChemComp.ChemCompVar(chemComp,chemAtoms = calist,
                                        linking = linking,
                                        descriptor = descriptor,
                                        formalCharge = refCcv.formalCharge,
                                        isParamagnetic = refCcv.isParamagnetic,
                                        isAromatic = refCcv.isAromatic)

  #
  # Do final check and save
  #

  if chemComp.isModified:

    chemComp.checkAllValid(complete = True)
    print "  Modified & valid..."

    if saveChemComp and messageReporter.showYesNo('Confirm','Do you really want to create these new chemComps?'):

      chemComp.save()


def cleanupDictByLink(ccv,linkingDict,newCodeDict,boundChemAtomNames):

  """

  Auxilary to createNewDescriptors

  """

  codeKey = 'link'

  if newCodeDict.has_key('link'):
    atomNames = newCodeDict[codeKey][:]
    for atomName in atomNames:

      codeKeyIndex = newCodeDict[codeKey].index(atomName)

      if linkingDict.has_key(codeKey) and atomName in linkingDict[codeKey]:
          newCodeDict[codeKey].pop(codeKeyIndex)
          continue

      if boundChemAtomNames.has_key(atomName) and not ccv.findFirstChemAtom(name = atomName):
        for boundChemAtomName in boundChemAtomNames[atomName]:
          if boundChemAtomName in linkingDict[codeKey]:
            newCodeDict[codeKey].pop(codeKeyIndex)

    if not newCodeDict[codeKey]:
      del(newCodeDict[codeKey])


def getLinkAtom(chemAtom,bondType = 'single'):

  linkAtomName = "%s_1" % chemAtom.name
  linkAtom = chemAtom.parent.findFirstChemAtom(className = 'LinkAtom', name = linkAtomName)
  if not linkAtom:
    linkAtom = ChemComp.LinkAtom(chemAtom.parent,name = linkAtomName)
    ChemComp.ChemBond(chemAtom.parent,chemAtoms = [chemAtom,linkAtom], bondType = bondType)

  return linkAtom

def getResAtomObjectMapping(resonances,chainCodeFilter=None):

  """

  Function that maps resonances to residues and actual atom names (not sets)

  In case of prochiral atoms (two resonances, two atomSets) will only return
  one of the two atomSets

  """

  resMapping = {}

  #
  # First make quick link for resonance -> atom
  #

  for resonance in resonances:

    #
    # Only a link from the resonance to an atom if there is a resonanceSet...
    #

    resonanceSet = resonance.resonanceSet

    if resonanceSet:

      resonanceSetResonances = list(resonanceSet.sortedResonances())
      atomSets = resonanceSet.sortedAtomSets()

      if len(resonanceSetResonances) == len(atomSets) == 2:
        # Prochirals... assuming they are anyway
        rsrSerialDict = {}
        for rsr in resonanceSetResonances:
          rsrSerialDict[rsr.serial] = rsr
        rsrSerials = rsrSerialDict.keys()
        rsrSerials.sort()

        asNameDict = {}
        for ass in atomSets:
          asNameDict[ass.name] = ass
        asNames = asNameDict.keys()
        asNames.sort()

        resonanceIndex = rsrSerials.index(resonance.serial)
        atomSets = [atomSets[resonanceIndex]]

      residue = atomSets[0].findFirstAtom().residue
      
      if chainCodeFilter and residue.chain.code not in chainCodeFilter:
        continue

      #
      # Go over the atomSets and atoms...
      #

      atomList = []

      for atomSet in atomSets:

        refAtom = atomSet.findFirstAtom()
        curResidue = refAtom.residue

        #
        # Check if all is OK (no mapping to different residues)
        #

        if curResidue != residue:
          print "  ERROR two residues to same resonance!"
          atomList = []
          break

        for atom in atomSet.sortedAtoms():
          atomList.append(atom)


      #
      # Do some sorting...
      #

      resMapping[resonance] = [residue,tuple(atomList)]

  return (resMapping)

#
# Code for getting secondary structure information for a residue
#

def getResidueSsCode(residue, defaultSsCode = None):

  #
  # Warning: this is currently PDB and application data specific!
  #
  # Also should only use 'C' as ssCode default if there is ss info available for a fragment...
  #

  from ccp.general.Constants import pdbToDsspNamesDict, secStrucInfo_kw

  ssCode = defaultSsCode
  ssAppData = residue.findFirstApplicationData(keyword = secStrucInfo_kw)

  if ssAppData:
    try:
      secStrucInfo = eval(ssAppData.value)
      if secStrucInfo[0] == 'sheet':
        secStruc = 'sheet'
      elif secStrucInfo[0] == 'helix':
        secStruc = secStrucInfo[2]['class']
      elif secStrucInfo[0] == 'turn':
        secStruc = 'turn'
    except:
      # This is to handle CSI type data...
      secStruc = ssAppData.value

    ssCode = pdbToDsspNamesDict[secStruc]

  return ssCode

#
# Code for getting pseudoatom corrections
#

def getPseudoCorrectionsWuthrich(resonances):

  #
  # Returns dictionary with corrections for pseudoatoms based on Wuthrich 1986
  #

  from ccp.general.Constants import pseudoAtomCorrectionsWuthrich

  #
  # TODO also need to check whether two items that together are pseudo (HB2/HB3 type stuff)
  #

  pseudoCorrections = {}

  for resonance in resonances:

    #
    # Only a link from the resonance to an atom if there is a resonanceSet...
    #

    if resonance.resonanceSet:

      atomSets = list(resonance.resonanceSet.sortedAtomSets())
      residue = atomSets[0].findFirstAtom().residue

      if pseudoAtomCorrectionsWuthrich.has_key(residue.ccpCode):
        corrDict = pseudoAtomCorrectionsWuthrich[residue.ccpCode]
      else:
        continue
      
      for atomSet in atomSets:
      
        chemAtomSet = atomSet.findFirstAtom().chemAtom.chemAtomSet

        if chemAtomSet:
        
          if corrDict.has_key(chemAtomSet.name):
            if not pseudoCorrections.has_key(resonance):
              pseudoCorrections[resonance] = {}
            
            if not pseudoCorrections[resonance].has_key(chemAtomSet.name):
              pseudoCorrections[resonance][chemAtomSet.name] = corrDict[chemAtomSet.name]
         
          if chemAtomSet.chemAtomSet:
            chemAtomSetName = chemAtomSet.chemAtomSet.name
            if corrDict.has_key(chemAtomSetName):
              pseudoCorrections[resonance][chemAtomSetName] = corrDict[chemAtomSetName]

  return (pseudoCorrections)

def getPseudoCorrections(resonances):

  #
  # Returns dictionary with corrections for pseudoatoms based on reference data CCPN
  # These corrections are based on the values used in the CONCOORD software.
  #

  pseudoCorrections = {}

  for resonance in resonances:

    #
    # Only a link from the resonance to an atom if there is a resonanceSet...
    #

    if resonance.resonanceSet:

      atomSets = resonance.resonanceSet.sortedAtomSets()

      allAtoms = []

      for atomSet in atomSets:
        for atom in atomSet.sortedAtoms():
          allAtoms.append(atom)

      refAtom = allAtoms[0]
      refChemAtom = refAtom.chemAtom

      chemAtomSets = []

      if refChemAtom.chemAtomSet:

        if refChemAtom.chemAtomSet.distCorr:
          chemAtomSets.append((refChemAtom.chemAtomSet,list(refChemAtom.chemAtomSet.chemAtoms)))

        if refChemAtom.chemAtomSet.chemAtomSet:
          if refChemAtom.chemAtomSet.chemAtomSet.distCorr:
            chemAtoms = []
            for cas in refChemAtom.chemAtomSet.chemAtomSet.chemAtomSets:
              chemAtoms.extend(list(cas.chemAtoms))
            chemAtomSets.append((refChemAtom.chemAtomSet.chemAtomSet,chemAtoms))

      pseudoCorrection = 0.0

      for (chemAtomSet,chemAtoms) in chemAtomSets:
        allAtomsCopy = allAtoms[:]
        for atom in allAtoms:
          chemAtom = atom.chemAtom
          if chemAtom.chemAtomSet:
            if chemAtom.chemAtomSet == chemAtomSet:
              allAtomsCopy.pop(allAtomsCopy.index(atom))
              chemAtoms.pop(chemAtoms.index(chemAtom))
            if chemAtom.chemAtomSet.chemAtomSet and chemAtom.chemAtomSet.chemAtomSet == chemAtomSet:
              allAtomsCopy.pop(allAtomsCopy.index(atom))
              chemAtoms.pop(chemAtoms.index(chemAtom))

        if not allAtomsCopy and not chemAtoms:
          pseudoCorrections[resonance] = chemAtomSet.distCorr
          break

  return (pseudoCorrections)


def parseValueString(txt):
  """ parse and evaluate comma-separated list of subvalues
  subvalues are either numeric (all eval'able)
  or string (which are stripped)
  Empty strings are replaced by None
  """
 
  ll = [x.strip() for x in txt.split(',')]
  try:
    # treat as numeric
    vals = [eval(x) for x in ll]
  except:
    # not numeric - treat as string
    vals = [x or None for x in ll]
 
  if len(vals) == 1:
    # single value
    return vals[0]
  else:
    # list of values
    return vals


def parseBrukerSummary(fileName, **kw):
  """ parse Bruker summary experiment description file
      - a general keyword-value format.
  
  - starting and trailing spaces are ignored
  - Lines starting with '#' (possibly preceded by whitespace) are ignored
    NB - this does *not* hold for continuation lines
  - lines ending in backcslash "\", possibly followed by whitespace,
    are joined with the succeding line (including leading whitespace)
  - The following line types are read:
    - <keyword> = <value>
    - <keyword> = <tag>:<value> <tag>:<value> ...
  - values can be either single subvalues, or comma-separated lists of
     subvalues
  - subvalues can be either numeric or string; all subvalues in a value
    must fall in the same category. If 'eval' succeeds for all subvalues
    the 'eval'ed values are used ('numeric'), otherwise the raw strings are used.
  - keywords, tags, values and subvalues are stripped before processing. 
  - keywords and tags must be alphanumeric.
  - tags must be preceded by whitespace and followed immediately by 
    ': ' (colon, whitespace).
  """
  
  try:
    fp = open(fileName)
  except:
    print "Could not read file", fileName
    raise
  
  tagcolon = re.compile("\s+(\w+):\s+")  # Regular Expression for tag expression:
  commentline = re.compile("\s*#")       # Regular Expression for comment:
  dataline = re.compile("^\s*(\w+)\s*=") # Regular Expression for data line:
  
  result = {}
  
  # merge continuation lines
  lines = []
  doAppend = False
  for line in fp:
    
    line = line.rstrip()
  
    nextAppend = (line[-1] == '\\')
    if nextAppend:
      line = line[:-1] + ' '
      
    if doAppend:
      lines[-1] = lines[-1] + line
    else:
      lines.append(line)
      
    doAppend = nextAppend
  
  # extract keywords, tags, and values
  for line in lines:
    
    if commentline.match(line):
      # comment line - skip
      continue
  
    ll = dataline.split(line)
    if len(ll) == 1:
      # no match - skip
      continue
    
    else:
      # correct data line
      keyword = ll[1].strip()
      txt = ll[2]
      
      matches = tagcolon.split(txt)
      if len(matches) <>1:
        # mathes. contains tag:value pairs
        
        result[keyword] = dd = {}
        
        ii = 1
        while ii < len(matches):
          tag = matches[ii]
          val = parseValueString(matches[ii+1])
          dd[tag] = val
          ii += 2
          
      else:
        # no match - just keyowrd=value
        result[keyword] = parseValueString(txt)
  #
  return result



def getAmbigNamingPriority(resonance, chemAtomSet=None):
  """Determine if an ambigous prochiral Resonance or FixedResonance 
  (non-stereospecifically assigned) should be assigned to the first
  (return 0) or the second name (return 1) alphabetically
  
  ChemAtomSet must match the assignment of the Resonance. If not passed in it is
  generated from the atomSets
  
  The function gives reproducible results for the same assignment state,
  regardless of shift values.
  Results for proton and carbon of CH(CH3)2 groups are consistent.
  If Resonance.name are set according to Analysis shift value convention, the
  resultant priorities are reproduced.
  If resonance.name matches a ChemAtomOrSet exactly, this will be matched
  
  .. describe:: Input
  
  Nmr.Resonance or NmrConstraint.FixedResonance
  ChemComp.ChemAtomSet (optional)  
  
  .. describe:: Output
  
  0/1/None Index of preferred name in alphabetical list of alternative names,
  """
  
  # set up
  priority = None
  resonanceSet = resonance.resonanceSet
  resonanceName = resonance.name
  if resonanceSet:
  
    # get otherResonance
    ll = [x for x in resonanceSet.resonances if x is not resonance]
    if ll:
      otherResonance = ll[0]
    else:
      otherResonance = None
    
    if resonance.isotopeCode != '1H':
      # heavy atom resonance. 
      # First try using bound protons(needed e.g. for CH(CH3)2 groups)
      resonance2 = resonance.findFirstCovalentlyBound(isotopeCode='1H')
      if resonance2 is not None:
        # get priority of bound proton
        priority = getAmbigNamingPriority(resonance2)
      
      if priority is None and otherResonance is not None:
        resonance3 = otherResonance.findFirstCovalentlyBound(isotopeCode='1H')
        if resonance3 is not None:
          # get priority from bound proton of partner
          priority = getAmbigNamingPriority(resonance3)
          if priority is not None:
            priority = 1-priority
    
    if priority is None:
      # now use the resonance itself
      
      if otherResonance is None:
        if chemAtomSet is None:
          chemAtomSet = getChemAtomSetFromAtomSets(resonanceSet.atomSets)
        if chemAtomSet is not None:
          ll = chemAtomSet.sortedChemAtomSets() or chemAtomSet.sortedChemAtoms()
          names = [x.name for x in ll]
          if resonanceName in names:
            # Priority from exact match to ChemAtomOrSet.name
            priority = names.index(resonanceName)
        
        if priority is None:
          if resonanceName and (resonanceName.endswith('b') 
                                or resonanceName.endswith('b*')):
            # Get priority from Analysis-type assignment names
            priority = 1
          else:
            # Only resonance available, and nothing to work with. Use first name
            priority = 0
            
      
      else:
        # We have two resonances
        # sort on name, then serial to get priority
        if ( (resonanceName, resonance.serial) < 
            (otherResonance.name, otherResonance.serial)):
          priority = 0
        else:
          priority = 1
  
  #
  return priority


def getChemAtomSetFromAtomSets(atomSets):
  """ Get a ChemAtomSet that matches all atoms in atomSets
  
  .. describe:: Input
  
  List of Nmr.AtomSet or NmrConstraint.FixedAtomSet
  
  .. describe:: Output
  
  ChemComp.ChemAtomSet (or None)
  
  """
  
  atoms = []
  for atomSet in atomSets:
    atoms.extend(atomSet.atoms)
  
  chemAtoms = [x.chemAtom for x in atoms]
  if None in chemAtoms:
    return None
  
  chemAtomSets = set(x.chemAtomSet for x in chemAtoms)
  if None in chemAtomSets:
    return None
  
  nChemAtoms = 0
  for chemAtomSet in chemAtomSets:
    nChemAtoms += len(chemAtomSet.chemAtoms)
  if nChemAtoms != len(atoms):
    return None
    
  
  while len(chemAtomSets) > 1:
    chemAtomSets = set(x.chemAtomSet for x in chemAtomSets)
    if None in chemAtomSets:
      return None
  else:
    return chemAtomSets.pop()


def integerStringSortKey(key):
  """return sort key so that strings starting with an integer sort as if by integer"""

  result = list(key)
  for ii,val in enumerate(result):
    if isinstance(val, str):
      vv = val.strip()
      if vv and vv.isdigit():
        result[ii] = '%30s' % vv

  return result
