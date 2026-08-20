"""
======================COPYRIGHT/LICENSE START==========================

stereoAssignmentSwap.py: Swaps/combines stereo assignments based on restraints and coordinates

Copyright (C) 2005-2011 Wim Vranken (Vrije Universiteit Brussel) Jurgen Doreleijers (CMBI, Radboud Universiteit)

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

from ccp.format.nmrStar.bmrb.File import File
from ccp.general.Geometry import getDistanceFromCoordinates
from ccp.general.Util import getResonancesFromPairwiseConstraintItem
from ccpnmr.format.general.Util import getResNameText
from memops.api import Implementation
import math
import os




#
# Dictionary creation (could be useful for other scripts!)
#

def createResonanceAtomAndAtomSetDict(resonances):

  #
  # Creates a fully ambiguous resonance to atom dict: all atoms
  # 'on the other side' are included.
  #

  resAtomDict = {}
  resAtomSetDict = {}

  for resonance in resonances:

    resAtomDict[resonance] = []
    resAtomSetDict[resonance] = []

    if resonance.resonanceSet:
      for atomSet in resonance.resonanceSet.sortedAtomSets():
        resAtomSetDict[resonance].append(atomSet)
        for atom in atomSet.sortedAtoms():
          resAtomDict[resonance].append(atom)

    resAtomDict[resonance] = tuple(resAtomDict[resonance])
    resAtomSetDict[resonance] = tuple(resAtomSetDict[resonance])

  return resAtomDict,resAtomSetDict

def createResAtomSwapDict(resAtomSetDict,compareWithWattos=False):

  #
  # Create a list of the atoms if the assignments were 'swapped' for prochirals
  #

  resAtomSwapDict = {}
  prochiralResonancesDict = {}

  for resonance in resAtomSetDict:
    atomSetList = resAtomSetDict[resonance]

    #
    # Only stereo if corresponds to one atomSet only...
    #

    if len(atomSetList) == 1:

      refAtom = atomSetList[0].findFirstAtom()
      if not refAtom.chemAtom:
        print "  Error: no chemAtom for atom %s in residue %s, %s" % (refAtom.name,refAtom.residue.molType, refAtom.residue.ccpCode)
        continue
#      if compareWithWattos:
#          resName = refAtom.residue.ccpCode.upper()
#          if resName == 'PHE' or resName == 'TYR':
##              print 'DEBUG: skipping Phe/Tyr for comparison with Wattos.'
#              continue

      refChemAtomSet = refAtom.chemAtom.chemAtomSet

      if refChemAtomSet:

        atomList = getAtomsOtherProchiral(refAtom,prochiralResonancesDict,resonance)

        if atomList:
          resAtomSwapDict[resonance] = tuple(atomList)

  return (resAtomSwapDict,prochiralResonancesDict)

#
# General functions working with above dictionaries...
#

def getAtomsOtherProchiral(atom,prochiralResonancesDict,resonance):

  atomList = []

  chemAtom = atom.chemAtom
  chemAtomSet = chemAtom.chemAtomSet

  ccpCode = chemAtom.chemComp.ccpCode

  if chemAtomSet:

    otherChemAtoms = []
    prochiralChemAtomSet = None

    nmrConstraintStore = resonance.parent
    atomSet = atom.findFirstFixedAtomSet(nmrConstraintStore = nmrConstraintStore)

    #
    # Normal prochiral
    #

    if len(atomSet.atoms) == 1 and chemAtomSet.isProchiral:
      sortedChemAtomsList = list(chemAtomSet.sortedChemAtoms())
      chemAtomIndex = sortedChemAtomsList.index(chemAtom)
      prochiralChemAtomSet = chemAtomSet
      otherChemAtoms.append(sortedChemAtomsList[not chemAtomIndex])

    #
    # 'Deep' prochiral (e.g. HD1*, HD2* for Leu)
    #

    elif len(atomSet.atoms) == 3 and chemAtomSet.chemAtomSet and chemAtomSet.chemAtomSet.isProchiral:
      sortedChemAtomSetsList = list(chemAtomSet.chemAtomSet.sortedChemAtomSets())
      chemAtomSetIndex = sortedChemAtomSetsList.index(chemAtomSet)
      prochiralChemAtomSet = chemAtomSet.chemAtomSet
      otherChemAtomSet = sortedChemAtomSetsList[not chemAtomSetIndex]
      otherChemAtoms.extend(list(otherChemAtomSet.sortedChemAtoms()))

    #
    # Special cases.
    #

    elif (ccpCode == 'Asn' and chemAtomSet.name == 'HD2*') or \
         (ccpCode == 'Gln' and chemAtomSet.name == 'HE2*') or \
         (ccpCode in ('Phe','Tyr') and chemAtomSet.name in ('HD*','HE*')) or \
         (ccpCode in ('G','A','C','G','T','U','DUR') and chemAtomSet.name == "H2'*") or \
         (ccpCode in ('A') and chemAtomSet.name == 'H6*') or \
         (ccpCode in ('C') and chemAtomSet.name == 'H4*'):

      sortedChemAtomsList = list(chemAtomSet.sortedChemAtoms())
      chemAtomIndex = sortedChemAtomsList.index(chemAtom)
      prochiralChemAtomSet = chemAtomSet
      otherChemAtoms.append(sortedChemAtomsList[not chemAtomIndex])
   
    """
    elif (ccpCode in ('Phe','Tyr') and chemAtomSet.chemAtomSet and chemAtomSet.chemAtomSet.name == 'HD*|HE*'):
      sortedChemAtomSetsList = list(chemAtomSet.chemAtomSet.sortedChemAtomSets())
      chemAtomSetIndex = sortedChemAtomSetsList.index(chemAtomSet)
      prochiralChemAtomSet = chemAtomSet.chemAtomSet
      otherChemAtomSet = sortedChemAtomSetsList[not chemAtomSetIndex]
      otherChemAtoms.extend(list(otherChemAtomSet.sortedChemAtoms()))
    """


    if prochiralChemAtomSet:

      residue = atom.residue
      for otherChemAtom in otherChemAtoms:
        atomList.append(residue.findFirstAtom(chemAtom = otherChemAtom))

      seqId = residue.seqId
      chainCode = residue.chain.code
      prochiralKey = (chainCode,seqId,prochiralChemAtomSet)
      prochiralResonancesDict[resonance] = prochiralKey

  return atomList


#
# Is this generic?!??!
#

def getDistConstrSumAvg(distConstr,resAtomDict,atomCoordDict,factor,resAtomSwapDict,prochiralResonancesDict,swap = False):

  #
  # Calculate sum averaged distance (Nilges et al., Proteins 17, 297-309, 1993)
  #

  sideCodes = ['left','right']
  avgSum = {}
  avgDist = {}
  prochiralAvgSum = {}
  prochiralContribs = {}

  #
  # Do a 'left' and a 'right' component: this is to try and handle swaps when
  # both 'sides' of a restraint item have a prochiral group.
  #

  for sideCode in sideCodes:
    avgSum[sideCode] = 0
    avgDist[sideCode] = 0

  for sideCode in sideCodes:
    prochiralAvgSum[sideCode] = {}
    prochiralContribs[sideCode] = {}

  #
  # Now loop over all the items
  #

  for item in distConstr.items:

    #
    # Do first resonance
    #

    setLocalSumContribs(0,sideCodes,item,resAtomDict,swap,resAtomSwapDict,atomCoordDict,prochiralResonancesDict,prochiralAvgSum,avgSum)

    #
    # Do second resonance
    #

    setLocalSumContribs(1,sideCodes,item,resAtomDict,swap,resAtomSwapDict,atomCoordDict,prochiralResonancesDict,prochiralAvgSum,avgSum)

  #
  # Finally set the results...
  #

  for sideCode in sideCodes:

    if avgSum[sideCode]:
      avgDist[sideCode] = math.pow(1.0 / avgSum[sideCode], factor)

      for prochiralKey in prochiralAvgSum[sideCode].keys():

        #
        # Calculate the average contribution to the distance (for violation calc if necessary)
        #
        # Done according to Nilges, O'Donoghue in Meth. Enzym. 32 (1998) 107-139 (see page 113)
        #

        tempProchiralAvgSum = prochiralAvgSum[sideCode][prochiralKey][0]

        prochiralContribs[sideCode][prochiralKey] = (tempProchiralAvgSum / avgSum[sideCode], prochiralAvgSum[sideCode][prochiralKey][1])

  return avgDist, prochiralContribs

def setLocalSumContribs(resIndex,sideCodes,item,resAtomDict,swap,resAtomSwapDict,atomCoordDict,prochiralResonancesDict,prochiralAvgSum,avgSum):

  sideCode = sideCodes[resIndex]

  itemResonances = getResonancesFromPairwiseConstraintItem(item)

  resonance = itemResonances[resIndex]
  atomList = resAtomDict[resonance]

  if swap and resAtomSwapDict.has_key(resonance):
    atomList = resAtomSwapDict[resonance]

  for atom in atomList:

    if atomCoordDict.has_key(atom):

      otherResonance = itemResonances[not resIndex]
      otherAtomLists = [resAtomDict[otherResonance]]

      if resAtomSwapDict.has_key(otherResonance):
        otherAtomLists.append(resAtomSwapDict[otherResonance])

      avgLocalSums = [0.0,0.0]

      for i in range(0,len(otherAtomLists)):

        otherAtomList = otherAtomLists[i]

        for otherAtom in otherAtomList:
          if atom != otherAtom and atomCoordDict.has_key(otherAtom) and atomCoordDict[atom] and atomCoordDict[otherAtom]:

            avgSumContrib = 1.0 / math.pow(getDistanceFromCoordinates(atomCoordDict[atom],atomCoordDict[otherAtom]),6)
            avgLocalSums[i] += avgSumContrib

      #
      # Pick the largest average local sum: is the smallest distance!!
      #

      if avgLocalSums[0] < avgLocalSums[1]:
        avgLocalSum = avgLocalSums[1]
      else:
        avgLocalSum = avgLocalSums[0]

      avgSum[sideCode] += avgLocalSum

      if avgLocalSum:
        if resonance in prochiralResonancesDict.keys():
          prochiralKey = prochiralResonancesDict[resonance]

          if not prochiralAvgSum[sideCode].has_key(prochiralKey):
            prochiralAvgSum[sideCode][prochiralKey] = [0.00,[]]
          prochiralAvgSum[sideCode][prochiralKey][0] += avgLocalSum
          prochiralAvgSum[sideCode][prochiralKey][1].append(item)

#
# Make swapping itself a class with attributes for common variables...
#

class StereoAssignmentSwapCheck:

  def __init__(self,nmrConstraintStore,structureEnsemble,verbose = False):

    self.nmrConstraintStore = nmrConstraintStore
    self.structureEnsemble = structureEnsemble

    self.distanceConstraintLists = []
    for constraintList in self.nmrConstraintStore.sortedConstraintLists():
      if constraintList.className == 'DistanceConstraintList':
        self.distanceConstraintLists.append(constraintList)

    self.verbose = verbose

  def checkSwapsAndClean(self,method = 'SUM_AVERAGING', violationCodes = None, swapFraction = 0.75, deassignAll = False):

    if not self.distanceConstraintLists or not self.structureEnsemble or not self.structureEnsemble.models:
      print "Error: no constraint lists or no structures available! Aborting..."
      return

    print
    print "Checking swap status and cleaning prochiral groups in constraint lists..."
    print

    #
    # Make this self attribute, just in case
    #

    self.deassignAll = deassignAll

    #
    # Initialize... see parameters above for swapFraction
    #
    # Set a dictionary with violationCodes (what is a large violation?)
    #

    if violationCodes:
      self.violationCodes = violationCodes

    else:
      self.violationCodes = {}
      self.violationCodes['xl'] = {'violation': 2.0, 'fraction': 0.00001}
      self.violationCodes['l'] =  {'violation': 1.0, 'fraction': 0.5}

    #
    # Set up a dict of resonances
    # TODO: Should this be done in cleanStereoAssignments?!?
    #       Should I make this a 'violationHandler' class??!?!
    #

    (self.resAtomDict,self.resAtomSetDict) = createResonanceAtomAndAtomSetDict(self.distanceConstraintLists[0].parent.fixedResonances)
    if self.verbose:
      print "Made resAtomDict, resAtomSetDict"

    (self.resAtomSwapDict,self.prochiralResonancesDict) = createResAtomSwapDict(self.resAtomSetDict)
    if self.verbose:
      print "Made resAtomSwapDict,prochiralResonancesDict"

    structureViolations = []

    """
    infoStrings = []
    for resonance in resAtomSwapDict.keys():
      atoms = resAtomDict[resonance]
      infoString = "%3d.%s" % (atoms[0].residue.seqId, atoms[0].atomSet.name)
      atoms = resAtomSwapDict[resonance]
      infoString += "-> swap is %3d.%s" % (atoms[0].residue.seqId, atoms[0].name)
      infoStrings.append(infoString)
    infoStrings.sort()
    for infoString in infoStrings:
      print infoString
    """

    #
    # Set the factor for calculating violations
    #

    if method == 'SUM_AVERAGING':
      factor = 1.0/6.0

    #
    # Loop over the structures
    #

    self.models = self.structureEnsemble.sortedModels()
    #self.models = self.models[:3]
    for model in self.models:

      self.prochiralViolationDict = self.createProchiralViolationDict()

      #violationList = distanceConstraintLists[0].structureGeneration.newViolationList(molStructures = [model])

      totalViols = 0

      #
      # Set up dict for coordinates
      #

      self.atomCoordDict = {}

      # No need for sorting here...
      for cChain in self.structureEnsemble.coordChains:
        for cRes in cChain.residues:
          for cAtom in cRes.atoms:
            if cAtom.atom:

              #
              # TODO: cannot handle multiple coords for one atom!!
              #

              self.atomCoordDict[cAtom.atom] = model.findFirstCoord(atom = cAtom)

      if self.verbose:
        print "Made atomCoordDict for model %d" % model.serial

      #
      # Go over the distance constraints
      #

      for distanceConstraintList in self.distanceConstraintLists:

        for distConstr in distanceConstraintList.sortedConstraints():

          # TODO: extend to other methods?

          if method == 'SUM_AVERAGING':

            #
            # Calculate sum averaged distance (Nilges et al., Proteins 17, 297-309, 1993)
            #

            (avgDist,prochiralContribs) = getDistConstrSumAvg(distConstr,self.resAtomDict,self.atomCoordDict,factor,self.resAtomSwapDict,self.prochiralResonancesDict)
            totalViols += self.setViolationDict(avgDist,prochiralContribs,distConstr,'orig')

            #
            # Now do the same but for the swap...
            #

            if prochiralContribs:
              (avgDist,prochiralContribs) = getDistConstrSumAvg(distConstr,self.resAtomDict,self.atomCoordDict,factor,self.resAtomSwapDict,self.prochiralResonancesDict,swap = True)
              self.setViolationDict(avgDist,prochiralContribs,distConstr,'swap')

      structureViolations.append(self.prochiralViolationDict)

      if self.verbose:
        print "Total violations %d" % totalViols

    #
    # Check whether original or swap state was the best...
    #
    # TODO: might have to do this in a specific order, based on total number of restraints it's involved in, then number of unambiguous restraints...
    # see http://nmr.cmbi.ru.nl/~jd/wattos/doc/Wattos/Soup/Constraint/AssignStereo.html#doAssignStereo%28float,%20float,%20float,%20float,%20float,%20float,%20java.lang.String%29
    #

    prochiralSwaps = {}

    for i in range(len(structureViolations)):
      prochiralViolationDict = structureViolations[i]

      for prochiralKey in prochiralViolationDict['orig']['total']:

        #
        # Calculate total violation over all structures for swapping...
        #

        if prochiralViolationDict['orig']['total'][prochiralKey] > prochiralViolationDict['swap']['total'][prochiralKey]:
          if not prochiralSwaps.has_key(prochiralKey):
            prochiralSwaps[prochiralKey] = 0
          prochiralSwaps[prochiralKey] += 1


    #
    # Finally make the changes in the data model where appropriate
    #

    cutoff = swapFraction * len(structureViolations)
    self.infoText = {'deassign': [], 'swap': []}

    totalStructures = len(self.models)

    for prochiralKey in prochiralViolationDict['orig']['total']:

      #
      # Set the swapstatus
      #

      if prochiralSwaps.has_key(prochiralKey) and prochiralSwaps[prochiralKey] > cutoff:
        swapStatus = 'swap'
        percent = prochiralSwaps[prochiralKey] * 100 / len(structureViolations)
        self.infoText['swap'].append("Swapping stereo assignment status chain %s, residue %3d group %s: %.3f%%" % (prochiralKey[0],prochiralKey[1],prochiralKey[2].name,percent))

      else:
        swapStatus = 'orig'

      #
      # Check for constraints with large violations based on the swap state...
      #

      self.constraintItemsReset = []

      for violationCode in self.violationCodes.keys():

        constraintsViolated = {}
        cutoffFraction = self.violationCodes[violationCode]['fraction']

        #
        # Again look over all structures
        #

        for prochiralViolationDict in structureViolations:
          for distConstr in prochiralViolationDict[swapStatus][violationCode][prochiralKey].keys():
            if not constraintsViolated.has_key(distConstr):
              constraintsViolated[distConstr] = [0.0,[]]
            constraintsViolated[distConstr][0] += 1.0

            for constrItem in prochiralViolationDict[swapStatus][violationCode][prochiralKey][distConstr][1]:
              if constrItem not in constraintsViolated[distConstr][1]:
                constraintsViolated[distConstr][1].append(constrItem)

        #
        # Then recalculate and add a constraint item if necessary
        #

        distConstraints = constraintsViolated.keys()
        distConstraints.sort()

        for distConstr in distConstraints:

          fractionViolated = constraintsViolated[distConstr][0] / totalStructures
          #print prochiralKey, violationCode, fractionViolated
          if fractionViolated >= cutoffFraction:

            prochiralResonances = []
            for resonance in self.prochiralResonancesDict.keys():
              if self.prochiralResonancesDict[resonance] == prochiralKey:
                prochiralResonances.append(resonance)

            #
            # Have to make a new prochiral resonance if there's only one!!
            #

            if len(prochiralResonances) == 1:

              otherProchiralResonance = self.nmrConstraintStore.newFixedResonance(isotopeCode = prochiralResonances[0].isotopeCode)
              print "NEW1 resonance %d" % otherProchiralResonance.serial
              otherAtoms = self.resAtomSwapDict[prochiralResonances[0]]
              otherAtomSet = otherAtoms[0].findFirstFixedAtomSet(nmrConstraintStore=self.nmrConstraintStore)
              if not otherAtomSet:
                otherAtomSet = self.nmrConstraintStore.newFixedAtomSet(atoms = otherAtoms)

              otherProchiralResonanceSet = self.nmrConstraintStore.newFixedResonanceSet(resonances = [otherProchiralResonance], atomSets = [otherAtomSet]) #@UnusedVariable

              prochiralResonances.append(otherProchiralResonance)
              self.prochiralResonancesDict[otherProchiralResonance] = prochiralKey

            violatedItems = constraintsViolated[distConstr][1]

            # Always reset the violated restraint
            self.resetConstraintItems(violatedItems,prochiralResonances,prochiralKey,violationCode,fractionViolated)

            # If using deassignAll option, also deassign all other *distance* constraints
            if self.deassignAll:

              # TODO/Note: should do this for other types of constraints?! H bonds at least...
              for prochiralResonance in prochiralResonances:
                allDistConstrItems = []
                for constrItem in prochiralResonance.sortedPairwiseConstraintItems():
                  if constrItem.className == "DistanceConstraintItem":
                    allDistConstrItems.append(constrItem)
                self.resetConstraintItems(allDistConstrItems,prochiralResonances, prochiralKey,violationCode,fractionViolated,basedOnOtherConstraint=distConstr.serial)

              break

      #
      # Reset based on swapStatus
      #

      if swapStatus == 'swap':

        #
        # Switch the assignments...
        #

        prochiralResonances = []
        for resonance in self.prochiralResonancesDict.keys():
          if self.prochiralResonancesDict[resonance] == prochiralKey:
            prochiralResonances.append(resonance)

        if len(prochiralResonances) == 2:

          resSet1 = prochiralResonances[0].resonanceSet
          atomSet1 = resSet1.sortedAtomSets()[0]
          resSet2 = prochiralResonances[1].resonanceSet
          atomSet2 = resSet2.sortedAtomSets()[0]

          resSet1.addAtomSet(atomSet2)
          resSet1.removeAtomSet(atomSet1)
          resSet2.addAtomSet(atomSet1)
          resSet2.removeAtomSet(atomSet2)

        else:
          resSet = prochiralResonances[0].resonanceSet
          atomSet1 = resSet.sortedAtomSets()[0]

          otherAtoms = self.resAtomSwapDict[prochiralResonances[0]]

          otherAtomSet = otherAtoms[0].findFirstFixedAtomSet(nmrConstraintStore=self.nmrConstraintStore)
          if not otherAtomSet:
            otherAtomSet = self.nmrConstraintStore.newFixedAtomSet(atoms = otherAtoms)

          resSet.addAtomSet(otherAtomSet)
          resSet.removeAtomSet(atomSet1)

    #
    # Print out info - TODO need this as XML?
    #

    infoTypes= self.infoText.keys()
    infoTypes.sort()

    for infoType in infoTypes:
      for line in self.infoText[infoType]:
        print line
      print

  def resetConstraintItems(self,constraintItems,prochiralResonances,prochiralKey,violationCode,fractionViolated,basedOnOtherConstraint=None,verbose=True):
    
    deassign = True
    
    # Make new resonance if required...
    if len(prochiralResonances) == 1:
      prochiralResonance = prochiralResonances[0]
      resSet = prochiralResonance.resonanceSet
      
      # Only do this if it's not an aromatic HD*/HE*!
      if len(resSet.atomSets) != 1 or len(resSet.sortedAtomSets()[0].atoms) != 2:
            
        otherAtoms = self.resAtomSwapDict[prochiralResonance]
        atomSet = self.nmrConstraintStore.newFixedAtomSet(atoms=otherAtoms)
        resSet.addAtomSet(atomSet)
              
        otherProchiralResonance = self.nmrConstraintStore.newFixedResonance(name=otherAtoms[0].name,isotopeCode=prochiralResonance.isotopeCode,resonanceSet=resSet)
        print "NEW2 resonance %d" % otherProchiralResonance.serial
        for ats in resSet.sortedAtomSets():
          if ats == atomSet:
            tt = "NEW"
          else:
            tt = "EXI"
          print "  %s AS" % tt,ats,ats.sortedAtoms()
  
        atomList = []
        for atomSet in resSet.sortedAtomSets():
          for atom in atomSet.sortedAtoms():
            atomList.append(atom)
        
        # Reset both resonances!
        self.resAtomDict[otherProchiralResonance] = tuple(atomList)
        self.resAtomDict[prochiralResonance] = tuple(atomList)
        
        prochiralResonances.append(otherProchiralResonance)
    
      else:
        
        # No point in deassigning if already HE* and HD* for aromatics!
        deassign = False
    
    
    if deassign:

      # Now deassign all constraints...
      for constraintItem in constraintItems:
  
        #
        # Don't redo the item if it was already reset..
        #
  
        if (prochiralKey,constraintItem) in self.constraintItemsReset:
          continue
  
        distConstr = constraintItem.constraint
        resonances = getResonancesFromPairwiseConstraintItem(constraintItem)
        newResonances = None
  
        for resonance in prochiralResonances:
          if resonance in resonances:
            otherResonanceIndex = not resonances.index(resonance)
            otherResonance = resonances[otherResonanceIndex]
  
            otherProchiralResonance = prochiralResonances[not prochiralResonances.index(resonance)]
  
            # Put these the wrong way around on purpose!
            if otherResonanceIndex == 1:
              newResonances = [otherResonance,otherProchiralResonance]
            else:
              newResonances = [otherProchiralResonance,otherResonance]
  
            break
  
        if newResonances and newResonances[0] != newResonances[1]:
  
          constraintExists = distConstr.findFirstItem(resonances = tuple(newResonances))
  
          if not constraintExists:
            newResonances.reverse()
  
            constraintExists = distConstr.findFirstItem(resonances = tuple(newResonances))
  
            if not constraintExists:
            
              className = distConstr.className
              constructor = getattr(distConstr,'new%sItem' % className)
  
              constructor(resonances = newResonances)
  
              prochiralText = "%s.%d.%s" % (prochiralKey[0],prochiralKey[1],prochiralKey[2].name)
  
              if verbose:
                if basedOnOtherConstraint:
                  infoLine = "  + based on constraint %d, deassigned constraint %d to prochiral %s.\n" % (basedOnOtherConstraint,distConstr.serial,prochiralText)
                  infoLine += "     --> added new item for '%s' to '%s'" % (getResNameText(newResonances[0]),getResNameText(newResonances[1]))
                else:
                  infoLine = "Deassigned constraint %d to prochiral %s: violation > %.1f in %.1f%% of structures.\n" % (distConstr.serial,prochiralText,self.violationCodes[violationCode]['violation'],fractionViolated * 100)
                  infoLine += "     --> added new item for '%s' to '%s'" % (getResNameText(newResonances[0]),getResNameText(newResonances[1]))
  
                self.infoText['deassign'].append(infoLine)
  
              self.constraintItemsReset.append((prochiralKey,constraintItem))

  def createProchiralViolationDict(self):

    prochiralViolationDict = {}

    for swapState in ['orig','swap']:
      prochiralViolationDict[swapState] = {}
      for violType in ['total','xl','l']:
        prochiralViolationDict[swapState][violType] = {}
        for prochiralKey in self.prochiralResonancesDict.values():
          if violType == 'total':
            prochiralViolationDict[swapState][violType][prochiralKey] = 0.0
          else:
            prochiralViolationDict[swapState][violType][prochiralKey] = {}

    return prochiralViolationDict

  def setViolationDict(self,avgDist,prochiralContribs,distConstr,swapState):

    totalViols = 0

    for sideCode in avgDist.keys():

      viol = 0

      if avgDist[sideCode]:
        if distConstr.lowerLimit != None and avgDist[sideCode] < distConstr.lowerLimit:
          viol = distConstr.lowerLimit - avgDist[sideCode]
        elif distConstr.upperLimit != None and avgDist[sideCode] > distConstr.upperLimit:
          viol = avgDist[sideCode] - distConstr.upperLimit

      if viol:

        totalViols += 1

        #
        # TODO: Violation creation not really necessary here - but need such a function
        #
        #violation = violationList.newViolation(constraint = distConstr, value = viol)

        for prochiralKey in prochiralContribs[sideCode]:
          violContrib = prochiralContribs[sideCode][prochiralKey][0] * viol
          self.prochiralViolationDict[swapState]['total'][prochiralKey] += violContrib

          for violationCode in self.violationCodes.keys():
            violLimit = self.violationCodes[violationCode]['violation']
            if violContrib > violLimit:
              #print prochiralKey, swapState, violationCode, violContrib, viol
              prochiralViolation = self.prochiralViolationDict[swapState][violationCode][prochiralKey]
              if not prochiralViolation.has_key(distConstr):
                prochiralViolation[distConstr] = [0,[]]
              prochiralViolation[distConstr][0] += 1
              prochiralViolation[distConstr][1].extend(prochiralContribs[sideCode][prochiralKey][1])

    if totalViols:
      totalViols = 1

    return totalViols

class StereoAssignmentCleanup(StereoAssignmentSwapCheck):

  """
  Updated version of above, based on Jurgen's algorithm in Wattos, with input from Geerten Vuister.

  This algorithm:

    1. Finds a list of 'triplets' (stereospecifically assigned prochirals or equivalents)

    2. Sorts this list according to -1- total number of restraints involved
                                    -2- number of restraints with unique assignments
    3. For each triplet T:

        1. Find set of restraints S containing T
        2. Calculates E, energy of S for the different models
        3. Calculates Eflip, energy of S with swapped stereospecific assignment for T
        4. If Eflip < E for >XX% of the models then
             swap the assignment for T
        5. Deassign the triplet T (for ALL constraints) if a single restraint is violated more than 2.0 angstrom in any model, or more than 1.0 angstrom in more than half of models (see 'violationCodes')

    For comparison the Wattos code for swapping (brackets removed for brevity):
        if (totalEnergy[0] <= totalEnergy[1]) {
            General.showDebug("criterium not met: totalEnergy[0] > totalEnergy[1]: " + totalEnergy[0] + " and " + totalEnergy[1]);
            return true;
        if (percentageModelFavoured < model_criterium) {
            General.showDebug("criterium not met: percentageModelFavoured >= model_criterium");
            return true;
        if (energyDifference < energy_abs_criterium) {
            General.showDebug("criterium not met: energyDifference >= energy_abs_criterium");
            return true;
        if (energyDifferencePercentage < energy_rel_criterium) {
            General.showDebug("criterium not met: energyDifferencePercentage >= energy_rel_criterium");
            return true;
        swapped = true;
        
    Difference between Wattos and FC code:
    - FC counts more ambis than Wattos does when there are deassigned atoms on the 'other' side that in other restraints are SSA.
        See comment 58 by Wim at http://code.google.com/p/cing/issues/detail?id=231    
  """

  # created by $CINGROOT/python/cing/Database/Scripts/extractPseudoAtomNameMap.py
  # removed variants by hand. Need not be complete.
  mapCcpn2IupacPseudo = {
    "ALA,HB*": "MB",
    "ARG,HB*": "QB",
    "ARG,HG*": "QG",
    "ARG,HD*": "QD",
    "ARG,HH1*": "QH1",
    "ARG,HH2*": "QH2",
    "ASN,HB*": "QB",
    "ASN,HD2*": "QD",
    "ASP,HB*": "QB",
    "CYS,HB*": "QB",
    "GLN,HB*": "QB",
    "GLN,HG*": "QG",
    "GLN,HE2*": "QE",
    "GLU,HB*": "QB",
    "GLU,HG*": "QG",
    "GLY,HA*": "QA",
    "HIS,HB*": "QB",
    "ILE,HG2*": "MG",
    "ILE,HG1*": "QG",
    "ILE,HD1*": "MD",
    "LEU,HB*": "QB",
    "LEU,HD1*": "MD1",
    "LEU,HD2*": "MD2",
    "LEU,HD*": "QD",
    "LYS,HB*": "QB",
    "LYS,HG*": "QG",
    "LYS,HD*": "QD",
    "LYS,HE*": "QE",
    "LYS,HZ*": "QZ",
    "MET,HB*": "QB",
    "MET,HG*": "QG",
    "MET,HE*": "ME",
    "PHE,HB*": "QB",
    "PHE,HD*": "QD",
    "PHE,HE*": "QE",
    "PHE,HD*|HE*": "QR",
    "PRO,HB*": "QB",
    "PRO,HG*": "QG",
    "PRO,HD*": "QD",
    "SER,HB*": "QB",
    "THR,HG2*": "MG",
    "TRP,HB*": "QB",
    "TYR,HB*": "QB",
    "TYR,HD*": "QD",
    "TYR,HE*": "QE",
    "TYR,HD*|HE*": "QR",
    "VAL,HG1*": "MG1",
    "VAL,HG2*": "MG2",
    "VAL,HG*": "QG",
    "DA,H5'*": "Q5'",
    "DA,H2'*": "Q2'",
    "DA,H6*": "Q6",
    "A,H5'*": "Q5'",
    "A,H6*": "Q6",
    "DC,H5'*": "Q5'",
    "DC,H2'*": "Q2",
    "DC,H4*": "Q4",
    "C,H5'*": "Q5'",
    "C,H4*": "Q4",
    "DG,H5'*": "Q5'",
    "DG,H2'*": "Q2'",
    "DG,H2*": "Q2",
    "G,H5'*": "Q5'",
    "G,H2*": "Q2",
    "DT,H5'*": "Q5'",
    "DT,H2'*": "Q2'",
    "DT,H7*": "M7",
    "U,H5'*": "Q5'",
    "HOH,H*": "QH",
    }

  compareWithWattos = False

  VIOLATION_CODE_REPORTINGS_STR = 'reporting s' # essentially zero.
  VIOLATION_CODE_REPORTINGL_STR = 'reporting l' # just for getting _Stereo_assign.Multi_mdl_crit_count
  VIOLATION_CODE_REPORTINGX_STR = 'reporting x' # just for getting _Stereo_assign.Single_mdl_crit_count
  VIOLATION_CODE_REPORTING_LIST = [ VIOLATION_CODE_REPORTINGS_STR, VIOLATION_CODE_REPORTINGL_STR, VIOLATION_CODE_REPORTINGX_STR ]
  REQUIRES_DEASSIGNMENT_STR = 'requiresDeassignment'

  SWAP_TYPE_ORG = 'original'

  """
  Return True on error
  """

  def storeToAppData(self, star_text ):
    appData1 = Implementation.AppDataString(value=star_text, application='FormatConverter', keyword='stereoAssignmentCorrectionsFile')
    self.nmrConstraintStore.addApplicationData(appData1)
    print "Added STAR file to application data within nmrConstraintStore"

  def checkSwapsAndClean( self,              # For comparison the NRG tags and defaults on March 2nd, 2011 are presented.
                          energy_abs_criterium = 0.1,         # _Stereo_assign_list.Crit_abs_e_diff      0.100
                          energy_rel_criterium = 0.0,         # _Stereo_assign_list.Crit_rel_e_diff      0.000
                          swapPercentage = 75.0,              # _Stereo_assign_list.Crit_mdls_favor_pct  75.0
                          singleModelCutoff = 1.0,            # _Stereo_assign_list.Crit_sing_mdl_viol   1.000 (inclusive)
                          multiModelCutoff = 0.5,             # _Stereo_assign_list.Crit_multi_mdl_viol  0.500 (inclusive)
                          multiModelPercentageCutoff = 50.0,  # _Stereo_assign_list.Crit_multi_mdl_pct   50.0  (inclusive)
                          method = 'SUM_AVERAGING',           # TODO: code others.
                          outputFileName = 'stereo_assign.str', # will be written to current directory if not an absolute path. Ignored if output type is custom
                          debug = False,                      # Print debug info?
                          useLowestAromaticViolation = False, # Check for lowest violation for single HD1/2 HE1/2 distance constraint items
                          outputType = 'NMRSTAR'              # Will write out NMR-STAR file. Can also be 'custom', will then only print info
                          ):
    """
    Return True on error.
    """
    if not self.distanceConstraintLists or not self.structureEnsemble or not self.structureEnsemble.models:
      print "Error: no constraint lists or no structures available! Aborting..."
      return True

    #
    # Initialize... see parameters above for swapPercentage
    #
    # Set a dictionary with violationCodes (what is a large violation?)
    #
    #      smallFloat = 0.000000000001 # same for cutoff distance and fraction

    negativeFraction = -999.9 # fraction set to always happen as it's under cut off.

    self.violationCodes = {}
    self.violationCodes['xl'] =                               {'violation': singleModelCutoff, 'fraction': negativeFraction}
    self.violationCodes['l']  =                               {'violation': multiModelCutoff,  'fraction': multiModelPercentageCutoff/100.}
    self.violationCodes[self.VIOLATION_CODE_REPORTINGX_STR] = {'violation': singleModelCutoff, 'fraction': negativeFraction}
    self.violationCodes[self.VIOLATION_CODE_REPORTINGL_STR] = {'violation': multiModelCutoff,  'fraction': negativeFraction}
    self.violationCodes[self.VIOLATION_CODE_REPORTINGS_STR] = {'violation': 0.0,               'fraction': negativeFraction}


    # JFD changed indentation here so that below statement is always executed.
    # Order in which they are checked, if found will abort so xl violation is prioritized
    self.violationCodeList = ['xl','l',
        self.VIOLATION_CODE_REPORTINGS_STR,
        self.VIOLATION_CODE_REPORTINGL_STR,
        self.VIOLATION_CODE_REPORTINGX_STR ]
    for violationCode in self.violationCodeList:
        if not self.violationCodes.has_key(violationCode):
            print 'ERROR: expected violationCode [%s] in StereoAssignmentCleanup.violationCodes ' % violationCode
            return True
#        print 'DEBUG: self.violationCode[%s] : %s' % ( violationCode, str(self.violationCodes[violationCode]))

    #
    # Initialise some variables
    #

    self.useLowestAromaticViolation = useLowestAromaticViolation

    #
    # Set the factor for calculating violations
    #

    self.method = method
    if self.method == 'SUM_AVERAGING':
      self.factor = 1.0/6.0

    #
    # Initialise resonance and 'triplet' information
    #

    print
    print "Checking swap status and cleaning prochiral groups in constraint lists..."
    print

    (self.resAtomDict,self.resAtomSetDict) = createResonanceAtomAndAtomSetDict(self.distanceConstraintLists[0].parent.fixedResonances)
    if self.verbose:
      print "Made resAtomDict, resAtomSetDict"

    # resAtomSwapDict is list of atoms associated with a resonance, prochiralResonancesDict links to (chainCode,seqId,prochiralChemAtomSet) tuple
    (self.resAtomSwapDict,self.prochiralResonancesDict) = createResAtomSwapDict(self.resAtomSetDict,compareWithWattos=self.compareWithWattos)
    if self.verbose:
      print "Made resAtomSwapDict,prochiralResonancesDict"

    self.triplets = {}

    # Generate a list of triplets, only for ones that have resonances - rest is dealt with later on.
    resList = self.prochiralResonancesDict.keys()
    resList.sort()

    for res in resList:
      atomTuple = self.resAtomDict[res]
      prochiralKey = self.prochiralResonancesDict[res]

      if not self.triplets.has_key(prochiralKey):
        self.triplets[prochiralKey] = {}

      if not self.triplets[prochiralKey].has_key(atomTuple):
        self.triplets[prochiralKey][atomTuple] = []

      self.triplets[prochiralKey][atomTuple].append(res)

    #
    # Now prioritise the triplets...
    #

    prochiralPriority = {}
    self.prochiralConstraints = {}

    prochiralKeys = self.triplets.keys()
    prochiralKeys.sort()
    Triplet_count = len(prochiralKeys)
    if Triplet_count < 1:
        print "WARNING: expected at least one triplet. Are there SSA distance restraints available?"
        return
    invalidTripletCount = 0 # Like 1a24 1    185    LEU    CD* that is invalid and can easily be recognized because it gets no involved restraints.
    for prochiralKey in prochiralKeys:
      #print prochiralKey
      atomTuples = self.triplets[prochiralKey].keys()
      atomTuples.sort()
      connectedConstraints = []
      unambiguousStereoConstraints = []  # These are constraints where there is no additional stereo ambiguity in the constraint items involving the prochiral
      allResonancesSet = set()

      otherItems = {}

      for atomTuple in atomTuples:
        #print "",atomTuple,triplets[prochiralKey][atomTuple]
        for resonance in self.triplets[prochiralKey][atomTuple]:
          allResonancesSet.add(resonance) # Note will not add the same item twice, so this is fine!
          for constraintItem in resonance.pairwiseConstraintItems:
            constraint = constraintItem.constraint
            if not otherItems.has_key(constraint):
              otherItems[constraint] = {}

            # Track other resonance in the item for filtering out fully ambiguous restraints
            orderedResonances = list(constraintItem.orderedResonances)
            otherResonance = orderedResonances[not orderedResonances.index(resonance)]
            if otherResonance not in otherItems[constraint]: # Use this now for future Python3 compatibility
              otherItems[constraint][otherResonance] = set()
            otherItems[constraint][otherResonance].add(resonance)

            if constraint.className in ('DistanceConstraint','HBondConstraint'):
              if constraint not in connectedConstraints:
                connectedConstraints.append(constraint)
                # So only 'unambiguous' if the 'other' resonance in the item has a resonance assignment, is assigned to one atomSet, and is prochiral (so could be deassigned)
                if otherResonance.resonanceSet and len(otherResonance.resonanceSet.atomSets) == 1 and otherResonance in self.prochiralResonancesDict:
                  #if self.resAtomDict[resonance][0].residue.seqId == 48:
                  #  print self.resAtomDict[resonance], self.resAtomDict[otherResonance], otherResonance.resonanceSet.atomSets
                  unambiguousStereoConstraints.append(constraint)
                else:
                  pass
#                  print 'DEBUG: ambi in %s:\n    %s' % (prochiralKey, ccpnDistanceRestraintToString(constraint)) # JFD doesn't know how to easily show atoms here.

      #
      # Clean up restraints so that constraints that are already fully ambiguous for the prochiral resonances (and they point to exactly the same resonances) are not included in the list to check..
      #

      if len(allResonancesSet) > 1:
        for constraint in otherItems:
          allMatch = True
          for otherResonance in otherItems[constraint]:
            if allResonancesSet != otherItems[constraint][otherResonance]:
              allMatch = False
              break

          if allMatch:
            if constraint in connectedConstraints:
              connectedConstraints.pop(connectedConstraints.index(constraint))
            if constraint in unambiguousStereoConstraints:
              unambiguousStereoConstraints.pop(unambiguousStereoConstraints.index(constraint))

      #
      # Set their priority
      #

      chainIdCcpn = prochiralKey[0]
      resIdCcpn = prochiralKey[1]
      chemAtomSetName = prochiralKey[2].name
      priorityKey = (len(connectedConstraints),len(unambiguousStereoConstraints),chainIdCcpn,resIdCcpn,chemAtomSetName)
#      print "DEBUG: priorityKey:", priorityKey
      if not prochiralPriority.has_key(priorityKey):
        prochiralPriority[priorityKey] = []

      prochiralPriority[priorityKey].append(prochiralKey)

      connectedConstraints.sort()
      self.prochiralConstraints[prochiralKey] = connectedConstraints

    
    #
    # Sort by priority and reorganise...
    #
    
    priorityKeys = prochiralPriority.keys()

    ## custom sort needs to return an int.
    def tripletComparator(x, y):
        if x[0] != y[0]:
            return x[0] - y[0] # ascending connectedConstraints
#        if not self.compareWithWattos:
        if x[1] != y[1]:
            return y[1] - x[1] # ascending unambiguousStereoConstraints
        if x[2] != y[2]:
            if x[2] < y[2]: # descending chainIdCcpn character
                return  1
            else:
                return -1
        resIdX = int(x[3])
        resIdY = int(y[3])
        if resIdX != resIdY:
            return resIdY - resIdX # descending resIdCcpn
        if x[4] != y[4]:
            if x[4] < y[4]: # descending chemAtomSetName
                return  1
            else:
                return -1
        return 0
    # end def

    priorityKeys.sort(cmp=tripletComparator)
    priorityKeys.reverse()

    if debug:
        for pk in priorityKeys:
          for pck in prochiralPriority[pk]:
            print "pck: ", pck
            for at in self.triplets[pck].keys():
              print "  at, self.triplets[pck][at]: ",at, self.triplets[pck][at]
            print

    #
    # Now calculate the total 'energy' for each constraint, and track whether there are any serious violations
    #
    # The 'energy' is the sum of the squared violations (over all models and restraints).
    #

    self.createAtomCoordDict() # This is static, fine to keep like this!

    # Corresponds to the indexes of avgLocalSums

    self.swapTypes = [self.SWAP_TYPE_ORG,'swapped']
    self.constraintItemsReset = []

    #
    # First only do swapping...
    #

    swapInfo = {}
    orgMaxViolation = {}
    orgViolationSingleModelCriteriumCount = {}
    orgViolationMultiModelCriteriumCount = {}

    Swap_count = 0 # Using captial to distinguish from original FC and use exact same as Wattos.
    Deassign_count = 0
    Total_e_low_states = 0.0
    Total_e_high_states = 0.0
    tripletIdx = 0
    for priorityKey in priorityKeys:
      for prochiralKey in prochiralPriority[priorityKey]:
        tripletIdx += 1
        if debug:
          print prochiralKey

        (prochiralViolationInfo,allConstraintItems) = self.checkProchiralKeyConstraints(prochiralKey,debug)

        # Find max violation of original assignment
        orgMaxViolation[                        prochiralKey] = 0.0
        orgViolationSingleModelCriteriumCount[  prochiralKey] = 0
        orgViolationMultiModelCriteriumCount[   prochiralKey] = 0
        violResultTupleList = prochiralViolationInfo[self.SWAP_TYPE_ORG][self.REQUIRES_DEASSIGNMENT_STR]
        for violationCode, violationList in violResultTupleList:
            if violationCode == self.VIOLATION_CODE_REPORTINGS_STR: # Includes any possible violation.
                orgMaxViolation[prochiralKey] = max( orgMaxViolation[prochiralKey], max(violationList)) # a list of violations
            elif violationCode == self.VIOLATION_CODE_REPORTINGX_STR:
                orgViolationSingleModelCriteriumCount[prochiralKey] += self.numModels - violationList.count(0.0)
            elif violationCode == self.VIOLATION_CODE_REPORTINGL_STR:
                orgViolationMultiModelCriteriumCount[prochiralKey]  += self.numModels - violationList.count(0.0)
        # end for violation results

        #
        # Now check whether needs to be swapped
        #

        doSwapCount = 0.0
        totalEnergyHighState = 0.0 # actual high state will be determined after next loop. For now assume state 0 (unswapped)
        totalEnergyLowState = 0.0
        for modelIndex in range(self.numModels):
          energyHighState = prochiralViolationInfo[self.swapTypes[0]]['energy'][modelIndex]
          energyLowState = prochiralViolationInfo[self.swapTypes[1]]['energy'][modelIndex]

#          totalEnergyDiff = prochiralViolationInfo[self.swapTypes[0]]['energy'][modelIndex] - prochiralViolationInfo[self.swapTypes[1]]['energy'][modelIndex] # this is a bug? Needs to be cumulative over models.
          totalEnergyHighState += energyHighState
          totalEnergyLowState += energyLowState
          if energyHighState > energyLowState: # swapping needed because for this model the assumption on the unswapped being the highest energy state was correct
            doSwapCount += 1.0
#          print "DEBUG: tripletIdx,modelIndex,energyHighState,energyLowState: %s" % str((tripletIdx,modelIndex,energyHighState,energyLowState))
        # end for model loop
        swappedFavouredFraction = doSwapCount / self.numModels

        # Adapted from Wattos
        totalEnergyHighState /= self.numModels # For criteria it's important to use one that can be compared over entries. Ensemble size should not influence result.
        totalEnergyLowState /= self.numModels
        if totalEnergyHighState < totalEnergyLowState: # Get this right before deciding on swapping.
          tmpEnergy = totalEnergyHighState
          totalEnergyHighState = totalEnergyLowState
          totalEnergyLowState = tmpEnergy
        # end if
        energyDifference = totalEnergyHighState - totalEnergyLowState # guaranteed positive or zero
        totalEnergyDiff = energyDifference # FC name
        percentageModelFavoured = 100.0 * swappedFavouredFraction
        if totalEnergyHighState > 0.0: # Strange in Wattos code there's no safety on totalEnergyHighState being zero. Added here.
            energyDifferencePercentage = 100.0 * energyDifference / totalEnergyHighState
        else:
            energyDifferencePercentage = 0.0
            if energyDifference > 0.0:
                energyDifferencePercentage = 100.0
        # end if/else

        # If any criteria is not met then the assignment will be maintained.
        swapAssignment = False
        if totalEnergyHighState <= totalEnergyLowState:
            msg = "criterium not met: totalEnergyHighState > totalEnergyLowState: %.3f and %.3f" % ( totalEnergyHighState, totalEnergyLowState )
        elif percentageModelFavoured < swapPercentage:
            msg = "criterium not met: percentageModelFavoured >= swapPercentage: %.1f %.1f" % ( percentageModelFavoured, swapPercentage)
        elif energyDifference < energy_abs_criterium: # If diff is close to zero do nothing.
            msg = "criterium not met: energyDifference >= energy_abs_criterium: %.3f and %.3f" % ( energyDifference, energy_abs_criterium )
        elif energyDifferencePercentage < energy_rel_criterium:
            msg = "criterium not met: energyDifferencePercentage >= energy_rel_criterium: %.1f %.1f" % ( energyDifferencePercentage, energy_rel_criterium)
        else:
            swapAssignment = True
        # end if/else
        if not swapAssignment:
            print "DEBUG maintaining tripletIdx %s because %s" % ( tripletIdx, msg)
        else:
            print "DEBUG swapping    tripletIdx %s" % tripletIdx
        # end if
        finalSwapType = self.swapTypes[0]
        favouredPercent = (1 - swappedFavouredFraction) * 100.0
        if swapAssignment:
          finalSwapType = self.swapTypes[1]
          favouredPercent = 100.0 - favouredPercent
          Swap_count += 1

        Total_e_low_states += totalEnergyLowState
        Total_e_high_states += totalEnergyHighState
        swapInfo[prochiralKey] = (swapAssignment,finalSwapType,energyDifferencePercentage,totalEnergyDiff, totalEnergyHighState, totalEnergyLowState,
                                  favouredPercent,swappedFavouredFraction,tripletIdx)


        #
        # Now make changes in CCPN... deassignment gets priority over swapping.
        #

        if swapAssignment:

          prochiralResonances = []
          for resList in self.triplets[prochiralKey].values():
            for resonance in resList:
              if not resonance in prochiralResonances:
                prochiralResonances.append(resonance)

          #
          # Switch the assignments...
          #
          
          if debug:
            print
            print "SWAPPING", prochiralResonances
            print

          if len(prochiralResonances) == 2:

            resSet1 = prochiralResonances[0].resonanceSet
            atomSet1 = resSet1.sortedAtomSets()[0]
            resSet2 = prochiralResonances[1].resonanceSet
            atomSet2 = resSet2.sortedAtomSets()[0]

            resSet1.addAtomSet(atomSet2)
            resSet1.removeAtomSet(atomSet1)
            resSet2.addAtomSet(atomSet1)
            resSet2.removeAtomSet(atomSet2)

            # Reset some dictionaries as well - note that resAtomSwapDict gives atoms of the *other* prochiral, so below is correct!
            atomTuple1 = tuple(atomSet1.sortedAtoms())
            atomTuple2 = tuple(atomSet2.sortedAtoms())

            self.resAtomSwapDict[prochiralResonances[0]] = atomTuple2
            self.resAtomSwapDict[prochiralResonances[1]] = atomTuple1

            # Reset triplets info
            self.triplets[prochiralKey] = {}
            self.triplets[prochiralKey][atomTuple1] = [prochiralResonances[1]]
            self.triplets[prochiralKey][atomTuple2] = [prochiralResonances[0]]

          elif len(prochiralResonances) == 1:
            resSet = prochiralResonances[0].resonanceSet
            atomSet1 = resSet.sortedAtomSets()[0]

            otherAtoms = self.resAtomSwapDict[prochiralResonances[0]]

            otherAtomSet = otherAtoms[0].findFirstFixedAtomSet(nmrConstraintStore=self.nmrConstraintStore)
            if not otherAtomSet:
              otherAtomSet = self.nmrConstraintStore.newFixedAtomSet(atoms = otherAtoms)
              
            if otherAtomSet != atomSet1:
              resSet.addAtomSet(otherAtomSet)
              atomSet1.removeResonanceSet(resSet)

              # Reset some dictionaries as well - note that resAtomSwapDict gives atoms of the *other* prochiral, so below is correct!
              atomTuple1 = tuple(atomSet1.sortedAtoms())
              
            else:
              # Same atomSet, possible for HD1/2 HE1/2 aromatics
              atomTuple1 = otherAtoms

            self.resAtomSwapDict[prochiralResonances[0]] = atomTuple1

            # Reset triplets info
            self.triplets[prochiralKey] = {}
            self.triplets[prochiralKey][atomTuple1] = []
            self.triplets[prochiralKey][otherAtomSet] = [prochiralResonances[0]]

    #
    # Then do deassigning. and track info for final printout...
    #

    finalList = {}

    self.swapTypes = [self.SWAP_TYPE_ORG] # Swapped not necessary any more
    priorityCount = 0

    for priorityKey in priorityKeys:
      priorityCount += 1
      for prochiralKey in prochiralPriority[priorityKey]:

        if debug:
          print prochiralKey

        (prochiralViolationInfo,allConstraintItems) = self.checkProchiralKeyConstraints(prochiralKey,debug=debug)

        #
        # Now check whether needs to be deassigned
        #

        finalSwapType = self.SWAP_TYPE_ORG

        numViol = {}
        deassign = False

        violResultTupleList = prochiralViolationInfo[finalSwapType][self.REQUIRES_DEASSIGNMENT_STR]
        for violationCodeToTest in self.violationCodeList:
          if violationCodeToTest in self.VIOLATION_CODE_REPORTING_LIST:
              continue
          if deassign:
              continue
          fractionByViolationCode = self.violationCodes[violationCodeToTest]['fraction']
#          numViol[violationCodeToTest] = 0
          for violationCode, violationList in violResultTupleList:
              if violationCodeToTest != violationCode:
                  continue
              # Look for every violationCodeToTest (a large single model cutoff and a smaller multi model cutoff) if fraction is met.
              numViol = self.numModels - violationList.count(0.0)
              fractionFound = ( 1.0 * numViol ) / self.numModels
              if fractionFound >= fractionByViolationCode: # inclusive
                  if debug:
                      print "DEBUG:    DEASSIGNING BASED ON %s %s" % (violationCode, str(prochiralViolationInfo[finalSwapType][self.REQUIRES_DEASSIGNMENT_STR]))
                  deassign = True
                  Deassign_count += 1
                  break # no need to look at other potentially qualifying restraints
          # end for
        # end for violationCodeToTest

        # Retrieve the swap info...
        (swapAssignment,finalSwapType,energyDifferencePercentage,totalEnergyDiff, totalEnergyHighState, totalEnergyLowState,
            favouredPercent,swappedFavouredFraction,tripletIdx) = swapInfo[prochiralKey]

        chainCode = prochiralKey[0]
        seqId = prochiralKey[1]
        chemAtomSetName = prochiralKey[2].name
        ccpCode = prochiralKey[2].chemComp.ccpCode
        totalConstraints = priorityKey[0]
        ambiguousConstraints = priorityKey[1]

        maximum_violation                      = orgMaxViolation[                      prochiralKey]
        violation_single_model_criterium_count = orgViolationSingleModelCriteriumCount[prochiralKey]
        violation_multi_model_criterium_count  = orgViolationMultiModelCriteriumCount[ prochiralKey]

        # chainCode, seqId, ccpCode, chemAtomSetName, swapAssignment, favouredPercent, totalEnergyDiff, totalConstraints, unambiguousStereoConstraints, deassign, numVeryLargeViol, numLargeViol
#        dummyIdxForComparisonWithWattos = '1' # TODO: reset to sensible output. chainCode
#        mapChainId2Idx = { 'A': '1', 'B': '2', 'C': '3' }
#        if mapChainId2Idx.has_key(chainCode):
#            dummyIdxForComparisonWithWattos = mapChainId2Idx[chainCode]
        pseudoNameKey = '%s,%s' % (ccpCode.upper(), chemAtomSetName)
        iupacPseudo = chemAtomSetName
        if self.mapCcpn2IupacPseudo.has_key(pseudoNameKey):
            iupacPseudo = self.mapCcpn2IupacPseudo[ pseudoNameKey ]
        lineItem = "%1s %4d %5s %-10s"  % ( chainCode, seqId, ccpCode.upper(), iupacPseudo )
        lineItem += " %3d %-3s %7.1f %7.1f %6.1f"  % ( tripletIdx, booleanPythonToJavaStr(swapAssignment), favouredPercent, energyDifferencePercentage, totalEnergyDiff )
        lineItem += " %6.1f %6.1f %3d"  % ( totalEnergyHighState, totalEnergyLowState, totalConstraints )
        lineItem += " %3d"  % ( ambiguousConstraints )
        lineItem += " %-5s %7.3f"  % ( booleanPythonToJavaStr(deassign), maximum_violation )
        lineItem += " %3d %3d"  % ( violation_single_model_criterium_count, violation_multi_model_criterium_count)
        
        if totalConstraints:
            finalList[(chainCode,seqId,chemAtomSetName)] = lineItem
        else:
            print "warning skipping triplet without restraints: " + lineItem
            invalidTripletCount += 1
        # end if
        
        #
        # Now make changes in CCPN... deassignment gets priority over swapping.
        #


        if deassign:

          violationCode = 'xxx'
          fractionViolated = 0.00

          prochiralResonances = []
          for resList in self.triplets[prochiralKey].values():
            for resonance in resList:
              if not resonance in prochiralResonances:
                prochiralResonances.append(resonance)

          self.resetConstraintItems(allConstraintItems,prochiralResonances, prochiralKey,violationCode,fractionViolated,verbose=False)

    #
    # Print out for checking
    #
    
    if outputType == 'custom':
    
      print """# Columns below (* means new):
#   1 chainCode
#   2 seqId
#   3 ccpCode
#   4 chemAtomSetName
#   5 priority (1 was handled first)
#   6 swapAssignment
#   7 favouredPercent (so for the swapped state if swapped!)
#   8 energyDifferencePercentage (*)
#   9 totalEnergyDiff ensemble averaged
#  10 totalEnergyHighState ensemble averaged (*)
#  11 totalEnergyLowState ensemble averaged (*)
#  12 totalConstraints
#  13 ambiguousConstraints (optional)
#  14 deassign
#  15 maximumViolation (pre processing)
#  16 numVeryLargeViol (post processing TODO: check)
#  17 numLargeViol (post processing TODO: check)
"""

    finalIds = finalList.keys()
    finalIds.sort()

    meat = ''

    for finalId in finalIds:
      if outputType == 'custom':
        print finalList[finalId]
      else:
        meat += str( finalList[finalId] ) + '\n'

    
    #
    # NMR-STAR Wattos type output
    #
    #    meat = """
    #  A    4   Met  HB*         82 False  100.0    0.000   2   0 False   0.000   0   0
    #  A    5   Arg  HD*         81 False  100.0    0.000   4   2 False   0.000   0   0
    #  A    6   Leu  HB*         23 False   90.0   14.328  26   7 True    1.812  11   0
    #
    #       1   6 LEU QB  22 no   90.0  78.6  8.803 11.204  2.402 26 10 yes 2.200  11  11
    #       1   6 LEU QD   8 no    5.0   0.0  0.000  1.649  1.649 34 14 yes 1.651  19  22
    #       1   9 GLU QG  96 no  100.0   0.0  0.000  0.000  0.000 10  0 no  0.000   0   0
    #"""

    if outputType == 'NMRSTAR':

      # Let's do the same with a STAR table.
      if invalidTripletCount:
          print "Warning: found triplets without restraints."
      validTripletCount = Triplet_count - invalidTripletCount
      if validTripletCount < 1:
          print "Error: found no triplets with restraints."
          return True
      validTripletCount2 = len(finalIds) # double check.
      if validTripletCount != validTripletCount2:
          print "Error: found number of triplets with restraints %d but number of report list %d" % ( validTripletCount, validTripletCount2)
#          return True
          
      Swap_percentage = ( 100.0 * Swap_count ) / validTripletCount
      Deassign_percentage = ( 100.0 * Deassign_count ) / validTripletCount
      Model_count = self.numModels
      Crit_abs_e_diff = energy_abs_criterium
      Crit_rel_e_diff = energy_rel_criterium
      Crit_mdls_favor_pct = swapPercentage
      Crit_sing_mdl_viol = self.violationCodes['xl']['violation']
      Crit_multi_mdl_viol = self.violationCodes['l']['violation']
      Crit_multi_mdl_pct = self.violationCodes['l']['fraction'] * 100.0

      header = """data_entry


  save_assign_stereo
      _Stereo_assign_list.Sf_category          stereo_assignments
      _Stereo_assign_list.Triplet_count        %s
      _Stereo_assign_list.Swap_count           %s
      _Stereo_assign_list.Swap_percentage      %.1f
      _Stereo_assign_list.Deassign_count       %s
      _Stereo_assign_list.Deassign_percentage  %.1f
      _Stereo_assign_list.Model_count          %s
      _Stereo_assign_list.Total_e_low_states   %.1f
      _Stereo_assign_list.Total_e_high_states  %.1f
      _Stereo_assign_list.Crit_abs_e_diff      %.1f
      _Stereo_assign_list.Crit_rel_e_diff      %.1f
      _Stereo_assign_list.Crit_mdls_favor_pct  %.1f
      _Stereo_assign_list.Crit_sing_mdl_viol   %.3f
      _Stereo_assign_list.Crit_multi_mdl_viol  %.3f
      _Stereo_assign_list.Crit_multi_mdl_pct   %.1f""" % (
          validTripletCount,
          Swap_count,
          Swap_percentage,
          Deassign_count,
          Deassign_percentage,
          Model_count,
          Total_e_low_states,
          Total_e_high_states,
          Crit_abs_e_diff,
          Crit_rel_e_diff,
          Crit_mdls_favor_pct,
          Crit_sing_mdl_viol,
          Crit_multi_mdl_viol,
          Crit_multi_mdl_pct
      )


      explanations = """
      _Stereo_assign_list.Details
;

Description of the tags in this list:
*  1 * NMR-STAR 3 administrative tag
*  2 * NMR-STAR 3 administrative tag
*  3 * NMR-STAR 3 administrative tag
*  4 * Number of triplets (atom-group pair and pseudo)
*  5 * Number of triplets that were swapped
*  6 * Percentage of triplets that were swapped
*  7 * Number of deassigned triplets
*  8 * Percentage of deassigned triplets
*  9 * Number of models in ensemble
* 10 * Energy of the states with the lower energies summed for all triplets (Ang.**2) ensemble averaged
* 11 * Energy of the states with the higher energies summed for all triplets (Ang.**2) ensemble averaged
* 12 * Item 9-8
* 13 * Criterium for swapping assignment on the absolute energy difference (Ang.**2)
* 14 * Criterium for swapping assignment on the relative energy difference (Ang.**2)
* 15 * Criterium for swapping assignment on the percentage of models favoring a swap
* 16 * Criterium for deassignment on a single model violation (Ang.)
* 17 * Criterium for deassignment on a multiple model violation (Ang.)
* 18 * Criterium for deassignment on a percentage of models
* 19 * this tag

Description of the tags in the table below:
*  1 * Chain identifier (can be absent if none defined)
*  2 * Residue number
*  3 * Residue name
*  4 * Name of pseudoatom representing the triplet
*  5 * Ordinal number of assignment (1 is assigned first)
*  6 * 'yes' if assignment state is swapped with respect to restraint file
*  7 * Percentage of models in which the assignment with the lowest
        overall energy is favored
*  8 * Percentage of difference between lowest and highest overall energy
        with respect to the highest overall energy
*  9 * Difference between lowest and highest overall energy ensemble averaged
* 10 * Energy of the highest overall energy state (Ang.**2) ensemble averaged
* 11 * Energy of the lowest overall energy state (Ang.**2) ensemble averaged
* 12 * Number of restraints involved with the triplet. The highest ranking
        triplet on this number, is assigned first (optional)
* 13 * Number of restraints involved with the triplet that are ambiguous
        besides the ambiguity from this triplet
* 14 * 'yes' if restraints included in this triplet are deassigned
* 15 * Maximum unaveraged violation before deassignment (Ang.)
* 16 * Number of violated restraints above threshold for a single model
        before deassignment (given by Single_mdl_crit_count)
* 17 * Number of violated restraints above threshold for a multiple models
        before deassignment (given by Multi_mdl_crit_count)
;


      loop_
         _Stereo_assign.Chain_ID
         _Stereo_assign.Comp_index_ID
         _Stereo_assign.Comp_ID
         _Stereo_assign.Pseudo_Atom_ID
         _Stereo_assign.Num
         _Stereo_assign.Swapped
         _Stereo_assign.Models_favoring_pct
         _Stereo_assign.Energy_difference_pct
         _Stereo_assign.Energy_difference
         _Stereo_assign.Energy_high_state
         _Stereo_assign.Energy_low_state
         _Stereo_assign.Constraint_count
  """
  #    if not self.compareWithWattos:
      explanations += "       _Stereo_assign.Constraint_ambi_count\n"
      # end if
      explanations += """         _Stereo_assign.Deassigned
         _Stereo_assign.Violation_max
         _Stereo_assign.Single_mdl_crit_count
         _Stereo_assign.Multi_mdl_crit_count

"""

      footer = """    stop_

  save_

  """


      star_text = header + explanations + meat + footer

      starFile = File()
      if starFile.read(text=star_text):
          print "Error: reading STAR text by STAR api."
          return True
      if starFile.check_integrity():
          print "Error: STAR text failed integrity check."
          return True
      starFile.filename = outputFileName
      if starFile.write():
          print "Error: writing file %" % outputFileName
          return True
      if not os.path.exists(outputFileName):
          print "Error: failed to find STAR file %s" % outputFileName
          return True
#      print "Written meta data to STAR file: %s" % outputFileName      # already printed by write()
          
      
      self.storeToAppData( star_text )
    # end def

  """
  Returns: (prochiralViolationInfo, allConstraintItems)
                prochiralViolationInfo[swapType][self.REQUIRES_DEASSIGNMENT_STR].append((violationCode,violationInfo[violationCode]))
    So result[0][swapType][self.REQUIRES_DEASSIGNMENT_STR] is a list of tuples:
        violationCode,violationInfo[violationCode]
    With violationInfo[violationCode] being filled like:
         violationInfo[violationCode][modelIndex] = viol
    So the get the max back out we need to:
  """
  def checkProchiralKeyConstraints(self,prochiralKey,debug=0):

    constraints = self.prochiralConstraints[prochiralKey]
    atomTuples = self.triplets[prochiralKey].keys()

    allConstraintItems = []
    prochiralViolationInfo = {}
    for swapType in self.swapTypes:
      prochiralViolationInfo[swapType] = {'energy': [0.0] * self.numModels, self.REQUIRES_DEASSIGNMENT_STR: []}

    #if prochiralKey[1] == 25 and prochiralKey[2].name == 'HB*':
    #  debug = True
    #else:
    #  debug = False

    for constraint in constraints:
    
      #if constraint.serial == 1617:
      #  debug = True
      #else:
      #  debug = False

      if debug:
        print constraint
      lowerLimit = constraint.lowerLimit
      upperLimit = constraint.upperLimit

      if self.method == 'SUM_AVERAGING':

        # The first index tracks the current situation, the second the swapped one
        avgLocalSums = [[0.0] * self.numModels,[0.0] * self.numModels]

        # Check if this constraint was already treated, ignore if so
        itemAtomsList = []
        
        # This on required in case of HE1/HE2 in separate item and using lowest violation (e.g. for HE1) in one item, then use other violation to other atom (e.g. HE2) if the same
        # item reoccurs.
        checkAromaticLowestViolations = {}

        if self.useLowestAromaticViolation:

          #
          # First check whether there's any Phe/Tyr HD1/2 HE1/2 items - if so, will treat these differently when calculating distances (smallest distance taken)
          # This step here is required in case there are constraint items to both the HD1 and HD2 - in that case the special treatment is not required
          #

          for item in constraint.items:

            itemResonances = getResonancesFromPairwiseConstraintItem(item)

            resonance = itemResonances[0]
            otherResonance = itemResonances[1]

            atomList = self.resAtomDict[resonance]
            otherAtomList = self.resAtomDict[otherResonance]

            tmpList = [(resonance,atomList),(otherResonance,otherAtomList)]
            for (tmpResonance, tmpAtomList) in tmpList:
              if self.prochiralResonancesDict.has_key(tmpResonance) and prochiralKey != self.prochiralResonancesDict[tmpResonance] and len(tmpAtomList) == 1:
                otherProchiralChemAtomSet = self.prochiralResonancesDict[tmpResonance][2]
                if otherProchiralChemAtomSet.name in ('HD*','HE*') and otherProchiralChemAtomSet.chemComp.ccpCode in ('Phe','Tyr'):
                  otherInfo = tmpList[not tmpList.index((tmpResonance,tmpAtomList))]
                  otherLowestViolationAtom = self.resAtomSwapDict[tmpResonance]

                  checkLowestViolationItemExists = (otherInfo[1],otherLowestViolationAtom)

                  setItem = True
                  for tmpItem in checkAromaticLowestViolations.keys():
                    if checkLowestViolationItemExists in checkAromaticLowestViolations[tmpItem][0]:
                      # Remove and ignore this one if it exists
                      del(checkAromaticLowestViolations[tmpItem])
                      setItem = False

                  if setItem:
                    # So set this, and if there is an item to the other aromatic (e.g. if this one is HD1 and there's an item to HD2) then remove it again
                    checkAromaticLowestViolations[item] = ((otherInfo[1],tmpAtomList),otherLowestViolationAtom,itemResonances.index(tmpResonance))
        
        #
        # Now do all items
        #

        for item in constraint.items:

          allConstraintItems.append(item)

          itemResonances = getResonancesFromPairwiseConstraintItem(item)

          resonance = itemResonances[0]
          otherResonance = itemResonances[1]

          #if resonance.serial == 1617 or otherResonance.serial == 1617:
          #  debug = True
          #else:
          #  debug = False
          
          if debug:
            print "   R     ",resonance, resonance.resonanceSet.atomSets
            print "   OtherR",otherResonance, otherResonance.resonanceSet.atomSets

          itemAtoms = (self.resAtomDict[resonance],self.resAtomDict[otherResonance])
          if itemAtoms in itemAtomsList:
            continue
          else:
            itemAtomsList.append(itemAtoms)

          atomLists = [self.resAtomDict[resonance]]
          otherAtomLists = [self.resAtomDict[otherResonance]]

          checkLowestViolationAtomListDict = {} #@UnusedVariable

          # Do the swap, only one of these should ever apply, so either resonance or otherResonance...
          atomListsAdded = 0
          for (tmpResonance, tmpAtomLists) in ((resonance,atomLists),(otherResonance,otherAtomLists)):
            if self.prochiralResonancesDict.has_key(tmpResonance) and prochiralKey == self.prochiralResonancesDict[tmpResonance]:

              # WARNING/TODO: what if HE1 only? See 2knr!
              addAtomList = self.resAtomSwapDict[tmpResonance]

              if addAtomList:
                tmpAtomLists.append(addAtomList)
                atomListsAdded += 1
              else:
                print prochiralKey, atomTuples, tmpAtomLists[0]
                sys.exit() #@UndefinedVariable
                              
          # end for
          if atomListsAdded == 2:
            print "  Warning: restraint item between prochiral resonances of same triplet %s.%d.%s. Ignoring... ." % (prochiralKey[0],prochiralKey[1],prochiralKey[2].name)
            continue

          if debug:
            print prochiralKey, atomListsAdded
            print "  AL", atomLists
            print "  OAL", otherAtomLists

          localSumsIndex = 0

          for atomList in atomLists:
            for otherAtomList in otherAtomLists:
              for atom in atomList:
                if self.atomCoordDict.has_key(atom):
                  for otherAtom in otherAtomList:

                    if atom == otherAtom:
                      print "  Warning: restraint between same atom."
                      continue

                    checkAromaticAtomInfo = None
                    if item in checkAromaticLowestViolations.keys():
                      checkAromaticAtomInfo = [atom,otherAtom]
                      
                      (_checkInfo,replaceAtom,atomIndex) = checkAromaticLowestViolations[item]
                      checkAromaticAtomInfo[atomIndex] = replaceAtom[0]
                      
                    # Get the distances from the coordinates
                    if self.atomCoordDict.has_key(otherAtom):
                      if self.atomCoordDict[atom] and self.atomCoordDict[otherAtom]:
                        for i in range(self.numModels):
                          avgSumContrib = 1.0 / math.pow(getDistanceFromCoordinates(self.atomCoordDict[atom][i],self.atomCoordDict[otherAtom][i]),6)
                          
                          # Check if another sum contrib needs to be calculated..
                          otherAvgSumContrib = 0.0
                          if checkAromaticAtomInfo:
                            otherAvgSumContrib = 1.0 / math.pow(getDistanceFromCoordinates(self.atomCoordDict[checkAromaticAtomInfo[0]][i],self.atomCoordDict[checkAromaticAtomInfo[1]][i]),6)
                          
                          # Take shortest distance (highest avgSumContrib)
                          avgLocalSums[localSumsIndex][i] += max(avgSumContrib,otherAvgSumContrib)
                          
                      else:
                        print "  Warning: empty coordinate for atom %s or %s" % (atom,otherAtom)
                    else:
                      print "  Warning: no coordinate for atom %s" % otherAtom
                else:
                  print "  Warning: no coordinate for atom %s" % atom

              localSumsIndex += 1
          # end for atomList

          # Make sure item info added if not a prochiral item
          if localSumsIndex == 1:
            for i in range(self.numModels):
              avgSumContrib = avgLocalSums[0][i]
              avgLocalSums[localSumsIndex][i] += avgSumContrib
        # at indent 9 end for item in constraint.items:

        # Now check for big violations, calculate energies, track for original and swapped situation. Is fine like this
        # because either original or swapped situation above.

        for i in range(len(self.swapTypes)):
          swapType = self.swapTypes[i]
          violationInfo = {}
          for violationCode in self.violationCodeList:
            violationInfo[violationCode] = [0.0] * self.numModels

          if debug:
            print "  ",swapType

          for modelIndex in range(self.numModels):
            # Possible that some data missing or not linked, ignore in that case.
            if not avgLocalSums[i][modelIndex]:
              continue

            avgDist = math.pow(1.0 / avgLocalSums[i][modelIndex], self.factor)

            if debug:
              print "    ",modelIndex, avgDist

            viol = None
            if lowerLimit and avgDist < lowerLimit:
              viol = lowerLimit - avgDist
            elif upperLimit and avgDist > upperLimit:
              viol = avgDist - upperLimit

            # If violated, add to overall prochiral energy and check whether it's a bad violation...
            if viol:
              if debug:
                print "    -> VIOL",viol,math.pow(viol,2)
              prochiralViolationInfo[swapType]['energy'][modelIndex] += math.pow(viol,2)

              for violationCode in self.violationCodeList:
                violLimit = self.violationCodes[violationCode]['violation']
                if viol > violLimit:
                  violationInfo[violationCode][modelIndex] = viol
#                  break # JFD notes: report for all violationCode instead so commented out here.
              # end for violationCode
            # end if viol
          # end for over models
#          print "DEBUG: Now check whether deassignment required for this particular restraint"
          for violationCode in self.violationCodeList:
            violationCount = self.numModels - violationInfo[violationCode].count(0.0)
            violationFraction = violationCount * 1.0 / self.numModels
            maxViol = max(violationInfo[violationCode]) #@UnusedVariable
#            if violationCode == self.VIOLATION_CODE_REPORTINGS_STR and maxViol > 0.2:
#            print "DEBUG: violationCode, maxViol, violationCount, violationFraction, self.violationCodes[violationCode]['fraction'] :", violationCode, maxViol, violationCount, violationFraction, self.violationCodes[violationCode]['fraction']
            if violationCount and (violationFraction > self.violationCodes[violationCode]['fraction']):
              prochiralViolationInfo[swapType][self.REQUIRES_DEASSIGNMENT_STR].append((violationCode,violationInfo[violationCode]))
          # end for violationCode
          if debug:
            print "  ", prochiralViolationInfo
            print
          # end if debug
        # at indent 9 end for over swapTypes
      # at indent 7 end if self.method
    # at indent 5 end for constraint
    return (prochiralViolationInfo, allConstraintItems)
  # end def


  def createAtomCoordDict(self):

    self.models = self.structureEnsemble.sortedModels()
    self.atomCoordDict = {}

    #self.models = self.models[:2]
    # Warning: might fail if not all atoms in all models.
    self.numModels = len(self.models)

    for cChain in self.structureEnsemble.coordChains:
      for cRes in cChain.residues:
        for cAtom in cRes.atoms:
          atom = cAtom.atom
          if atom:
            # Warning: not handling multiple coords for one atom!!
            self.atomCoordDict[atom] = []
            for model in self.models:
              self.atomCoordDict[atom].append(model.findFirstCoord(atom = cAtom))

def booleanPythonToJavaStr( b ):
    if b:
        return 'yes'
    return 'no'


def ccpnDistanceRestraintToString(ccpnConstraint):
    """
    """
    lowerLimit = None
    upperLimit = None
    if hasattr(ccpnConstraint, 'lowerLimit'):
        lowerLimit = ccpnConstraint.lowerLimit
        upperLimit = ccpnConstraint.upperLimit
    result = 'Constraint [%d]: [%s] - [%s]' % (ccpnConstraint.serial, lowerLimit, upperLimit)

    for constItem in ccpnConstraint.sortedItems():
        atomList = []
        # Sometimes there may also be an ordered<Class> method.
        resonanceList = None
        if hasattr(constItem, 'orderedResonances'): # dihedrals don't have this.
            resonanceList = constItem.orderedResonances
        # Otherwise, use the usual sorted<Class> method.
        if not resonanceList:
            resonanceList = constItem.sortedResonances()
        # Here, resonanceList should always have 2 resonances.

        resonanceListLength = len(resonanceList)
        assert(resonanceListLength == 2) # During a regular run (not with -O option given to python interpreter) this might cause a exception being thrown.
        if resonanceListLength != 2:
            print "ERROR: expected a pair but found number: %d for ccpnConstraint %s" % (resonanceListLength, ccpnConstraint)
            return None
        for resonance in resonanceList:
            resAtomList = []
            resonanceSet = resonance.resonanceSet
            if resonanceSet:
                for atomSet in resonanceSet.sortedAtomSets():
                    # atom set is a group of atoms that are in fast exchange and therefore are not assigned to individually (e.g. methyl group).
                    for atom in atomSet.sortedAtoms():
                        resAtomList.append('%d.%s' % (
                            atom.residue.seqCode, atom.name))
            else:
                print "WARNING: No resonanceSet (means unassigned) for ccpnConstraint %s" % ccpnConstraint
            resAtomList.sort()
            resAtomString = ','.join(resAtomList)
            atomList.append(resAtomString)
        result += '  [%s] - [%s]' % (atomList[0], atomList[1])
    # end for
    return result
# end def

if __name__ == '__main__':

  #inFile = '/ebi/msd/nmrqual/reference/AartProject/db700/monomers/2ezm/2ezm.xml'
  #inFile = '/ebi/msd/nmrqual/reference/AartProject/db700/monomers/1n6u/ccpn/1n6u.xml'
  #inFile = '/ebi/msd/nmrqual/reference/AartProject/db700/monomers/1bvm/1bvm.xml'
  #inFile = '/ebi/msd/nmrqual/reference/AartProject/db700/monomers/1d8v/1d8v.xml'
  #inFile = '/ebi/msd/nmrqual/reference/AartProject/db700/monomers/1eza/1eza.xml'
  #inFile = '/ebi/msd/nmrqual/reference/AartProject/db700/monomers/1tih/1tih.xml'
  #inFile = '/ebi/msd/nmrqual/reference/AartProject/db700/monomers/1i6g/1i6g.xml'
  #inFile = '/ebi/msd/nmrqual/reference/AartProject/db700/monomers/1jvr/1jvr.xml'
  inFile = '/ebi/msd/nmrqual/reference/AartProject/db700/monomers/1a6x/1a6x.xml'
  #inFile = '/ebi/msd/nmrqual/reference/AartProject/db700/monomers/1afi/1afi.xml'
  #inFile = '/ebi/msd/nmrqual/reference/AartProject/db700/monomers/1ah9/1ah9.xml'
  #inFile = '/ebi/msd/nmrqual/reference/AartProject/db700/monomers/1b3c/1b3c.xml'

  # Read xml file
  from memops.format.xml import XmlIO
  project = XmlIO.loadProjectFile(inFile)

  strucGen = project.currentNmrProject.sortedStructureGenerations()[1]

  sssc = StereoAssignmentSwapCheck(strucGen.nmrConstraintStore,strucGen.structureEnsemble,verbose = True)

  sssc.checkSwapsAndClean()
  sssc.checkSwapsAndClean()
