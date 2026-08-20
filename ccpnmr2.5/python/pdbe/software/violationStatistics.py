
"""
======================COPYRIGHT/LICENSE START==========================

ViolationStatistics.py: Script to generate distance constraint violation stats from CCPN project

Copyright (C) 2007 Wim Vranken (European Bioinformatics Institute)

=======================================================================

Can only be distributed and used within the CCPNGRID project, no other
use or distribution allowed.

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

import os, sys, math, string

# TODO this should be replaced by other functions in order to work for all!
# NB meanwhile moved into functions def to allow load, at least.
#from rpy import r

from memops.universal.Util import drawBox, getRms
from memops.general.Io import loadProject

from ccp.general.Util import getAtomPositionType, getSortContactDist
from ccp.general.Util import getResidueSsCode

# Angle and coordinate stuff
from memops.universal.Geometry import angleStats
from memops.universal.Geometry import angleNormaliseAroundZero
from ccp.general.Geometry import getDistanceFromCoordinates
from ccp.general.Geometry import calcTorsionAngleDegrees
from ccp.general.Geometry import calcAngleViolation

from pdbe.software.Util import HtmlPage, ResonanceCoordinateHandler, ContactOccurrenceHandler
from pdbe.software.Constants import styleSheet

contactDistToType = {

  0: 'intra',
  1: 'seq',
 -1: 'long',
 -2: 'inter'

}

#
# Violation analysis specific code....
#

class InfoTypes:

  dictInfo = {
    
      'contact': 'contactTypes',
      'atom':    'atomTypes',
      'secStruc':'ssTypes',
      'residues':'residues'
      
      }

   # Order in which produced is found in ViolationStatistics!!
   
class ViolationInfo(InfoTypes):

  keysInList = ['secStruc','residues']

  def __init__(self,violation,contactTypes,atomTypes,ssTypes,residues):
  
    self.violation = violation
    self.contactTypes = contactTypes
    self.atomTypes = atomTypes
    self.ssTypes = ssTypes
    self.residues = residues
    
  def setInfo(self,violationDict):

    violationDict['all'].append(self.violation)
    
    for dictKey in self.dictInfo.keys():

      otherDataDict = getattr(self,self.dictInfo[dictKey])
      dictSecondKeyList = otherDataDict.keys()

      # Can't really use item contributions for a violation, so using the whole thing...
      for dictSecondKey in dictSecondKeyList:
        if not violationDict[dictKey].has_key(dictSecondKey):
          violationDict[dictKey][dictSecondKey] = [[],0,0]
        violationDict[dictKey][dictSecondKey][0].append(self.violation)
        violationDict[dictKey][dictSecondKey][1] += otherDataDict[dictSecondKey] # Total item contrib
        violationDict[dictKey][dictSecondKey][2] += 1 # Total violations
      
    return self.violation

class ConstraintCount(InfoTypes):

  def __init__(self):
  
    self.infoDict = {}
    
    for dictKey in self.dictInfo.keys():
      self.infoDict[dictKey] = {}

  def addValues(self,contactTypes,atomTypes,ssTypes,residues):
  
    self.setListValue('contact',contactTypes)
    self.setListValue('atom',atomTypes)
    
    for ssType in ssTypes.keys():
      self.setValue('secStruc',ssType,ssTypes[ssType])
    
    residueKeys = {}
    for residue in residues.keys():
      residueKey = (residue.chain.code,residue.seqCode,residue.ccpCode)
      self.setValue('residues',residueKey,residues[residue])
      residueKeys[residueKey] = residues[residue]
    
    return residueKeys
    
  def setListValue(self,dictKey,dictSecondKeyDict):

    dictSecondKeys = dictSecondKeyDict.keys()
    
    for dictSecondKey in dictSecondKeys:
      self.setValue(dictKey,dictSecondKey,dictSecondKeyDict[dictSecondKey])
      
    return dictSecondKeys
  
  def setValue(self,dictKey,dictSecondKey,itemContrib):

    if not self.infoDict[dictKey].has_key(dictSecondKey):
      self.infoDict[dictKey][dictSecondKey] = [0,0]
    self.infoDict[dictKey][dictSecondKey][0] += itemContrib
    self.infoDict[dictKey][dictSecondKey][1] += 1   # Also track directly the number of constraints this element occurs, for use in violation analysis!
    
    return dictSecondKey
      
class ViolationStatistics(ResonanceCoordinateHandler,ContactOccurrenceHandler):

  # TODO should be customisable from GUI or something!

  infoTags =    ("Number","Total","Percent","Rms","Average","Median")  
  infoTypeKeyList = ('contact','atom','secStruc','residues')
  
  largeViolations = {'distance': 0.5, 'dihedral': 10.0}
  smallViolations = {'distance': 0.1, 'dihedral': 5.0}
  
  output = {'ensemble': {}, 'model': {}}
  
  getStatisticsExecuted = False
  
  # TODO: get this from contact occurrence stuff (if available)?
  maxRangeCheck = 6
          
  HtmlPage = HtmlPage
  styleSheet = styleSheet

  def __init__(self,nmrProject,strucGen = None, excludeStructures = None, useContactOccurrence = True, chainCodeFilter = None):
  
    self.projectName = nmrProject.name
    self.dataType = 'orig'
    
    if not excludeStructures:
      self.excludeStructures = []
    
    if not strucGen:
      strucGen = nmrProject.findFirstStructureGeneration()
    self.strucGen = strucGen
    
    self.structureList = strucGen.structureEnsemble.sortedModels()
    
    self.useContactOccurrence = useContactOccurrence
    
    self.chainCodeFilter = chainCodeFilter
    
  def getSecondaryStructureInfo(self):
  
    #
    # Get the original data with sec struc info
    #
    
    origResidueList = self.createOrigResidueList()

    #
    # Get the number of chains for non-original info
    #
    
    chainCodes = []
    if self.dataType != 'orig':
      for residue in self.assignedResidues:
        if residue.chain.code not in chainCodes:
          chainCodes.append(residue.chain.code)

    #
    # Get the info
    #
    
    self.residueSecStruc = {}
    origChainCodes = []

    for residue in origResidueList:
      
      ssCode = getResidueSsCode(residue, defaultSsCode = 'C')
      
      residueKey = self.getResidueKey(residue)
  
      self.residueSecStruc[residueKey] = ssCode  
      
      if not residueKey[0] in origChainCodes:
        origChainCodes.append(residueKey[0])
      
    #
    # Some after-interpretation
    #
    
    if len(origChainCodes) == 1 and len(chainCodes) == 1 and origChainCodes[0] != chainCodes[0]:
      for residueKey in self.residueSecStruc.keys()[:]:
        newResidueKey = (chainCodes[0],residueKey[1])
        self.residueSecStruc[newResidueKey] = self.residueSecStruc[residueKey]
        del(self.residueSecStruc[residueKey])

  def createOrigResidueList(self):

    return self.assignedResidues

  def getResidueKey(self,residue):
  
    return (residue.chain.code,residue.seqId)
  
  def getStatistics(self,testMode):
    
    if self.getStatisticsExecuted:
      return
  
    invSixth = -1.0 / 6.0

    ncs = self.strucGen.nmrConstraintStore
    fixedResonances = ncs.fixedResonances

    self.setAssignedAtomsAndResidues(fixedResonances,chainCodeFilter=self.chainCodeFilter)
    
    if self.useContactOccurrence:
      self.setupContactOccurrenceInfo(ncs)
        
    #
    # Get secondary structure info
    #
    
    self.getSecondaryStructureInfo()

    #
    # Set up reference info for structure coords
    #
        
    if testMode:
      self.structureList = self.structureList[:2]
    else:
      self.structureList = self.structureList[:]
      
    self.numStructures = len(self.structureList)

    if self.numStructures == 1:
      print "  Warning: only one structure used in analysis!"
    
    self.validStructures = self.numStructures
    if self.excludeStructures:
      for strucIndex in range(self.numStructures):
        if strucIndex in self.excludeStructures:
          self.validStructures -= 1
          
    #
    # Create a dictionary with coordinate info for the relevant atoms
    #
    
    self.createCoordAtomInfoDict()

    #
    # Get info on atoms close to atoms involved in constraints with unassigned stuff...
    #

    #
    # Now go over the constraint lists, get a list of atoms close to the 'missing' one...
    #

    #
    # TODO NEED TO COMBINE RESTRAINTS WITH SAME INFO FIRST AS WELL!!?!?! Not crucial at this stage..
    #    
    
    self.constraintCount = {}
    violatedConstraints = {}
    overallViolations = {}
    modelViolations = {}

    for strucIndex in range(self.numStructures):
      modelViolations[strucIndex] = {}
    
    allDistanceConstraintLists = list(ncs.findAllConstraintLists(className = 'DistanceConstraintList')) + list(ncs.findAllConstraintLists(className = 'HBondConstraintList'))

    if allDistanceConstraintLists:

      print "--> Found distance restraints for analysis"

      overallViolations['upper'] = {}
      overallViolations['lower'] = {}

      for strucIndex in range(self.numStructures):
        for violType in overallViolations.keys():
          modelViolations[strucIndex][violType] = {}

      numConstraints = 0

      self.constraintCount['distance'] = ConstraintCount()

      for dcl in allDistanceConstraintLists:

        for violType in overallViolations.keys():
          overallViolations[violType][dcl] = {}
          for strucIndex in modelViolations.keys():
            modelViolations[strucIndex][violType][dcl] = {}

        violatedConstraints[dcl] = []

        dclConstraints = dcl.sortedConstraints()

        for constraint in dclConstraints:

          numConstraints += 1

          constraintItems = constraint.sortedItems()
          
          if not constraintItems:
            continue

          #lowerDist = constraint.lowerLimit
          lowerDist = 1.0 # TODO have to check how low this can get between two methyls!
          upperDist = constraint.upperLimit
          
          if upperDist == None:
            upperDist = constraint.targetValue
          
          if upperDist == None:
            print ' No value for distance limit, ignoring...'
            continue

          distPerStruc = {}
          for strucIndex in range(self.numStructures):
            distPerStruc[strucIndex] = []
          hasDistance = False

          atomCombs = []
          contactTypes = {}
          atomTypes = {}
          ssTypes = {}
          residues = {}
          
          # Make sure to handle ambiguous constraints correctly!!
          itemContrib = 1.0 / len(constraintItems)

          for item in constraintItems:

            resonances = item.orderedResonances
            
            # Can't count on above returning something...
            if not resonances:
              resonances = list(item.resonances)

            allAtomList = []
            itemResidueKeys = []
            atomContactType = []

            for resIndex in range(2):
              resonance = resonances[resIndex]
              if self.resObjectMapping.has_key(resonance):
                (residue,atomList) = self.resObjectMapping[resonance]
                allAtomList.append(atomList)
                itemResidueKeys.append((residue.chain,residue.seqId))
                atomRefName = atomList[0].name
                atomType = getAtomPositionType(atomRefName)
                atomContactType.append(atomType)

                if not residue in residues.keys():
                  residues[residue] = 0
                residues[residue] += itemContrib / 2.0  # Here also divided in half, so each constraint counts for half for each of the residues involved...

                residueKey = self.getResidueKey(residue)
                ssCode = 'C'
                if residueKey in self.residueSecStruc.keys():
                  ssCode = self.residueSecStruc[residueKey]
                if not ssCode in ssTypes.keys():
                  ssTypes[ssCode] = 0
                ssTypes[ssCode] += itemContrib / 2.0  # Note that this is per resonances, so need to divide in half...

            atomContactType.sort()
            atomContactType = tuple(atomContactType)
            if not atomContactType in atomTypes.keys():
              atomTypes[atomContactType] = 0
            atomTypes[atomContactType] += itemContrib

            if len(allAtomList) != 2:
              continue

            (sortType,contactDist) = getSortContactDist(itemResidueKeys[0],itemResidueKeys[1],self.maxRangeCheck)

            if contactDistToType.has_key(contactDist):
              contactType = contactDistToType[contactDist]
            else:
              contactType = 'i+%d' % contactDist
            
            if not contactType in contactTypes.keys():
              contactTypes[contactType] = 0
            contactTypes[contactType] += itemContrib

            # Should in principle do sets, but this will do for now.
            for atom in allAtomList[0]:
              if self.coordAtomInfo.has_key(atom):
                for otherAtom in allAtomList[1]:
                  if self.coordAtomInfo.has_key(otherAtom):

                    atomComb = set((atom,otherAtom))

                    # Only count each contribution once!!
                    if atomComb in atomCombs:
                      continue

                    atomCombs.append(atomComb)

                    for strucIndex in range(self.numStructures):
                      if strucIndex in self.excludeStructures:
                        continue
                      coord = self.coordAtomInfo[atom][strucIndex]
                      otherCoord = self.coordAtomInfo[otherAtom][strucIndex]

                      if coord and otherCoord: 
                        distance = getDistanceFromCoordinates(coord,otherCoord)
                        distPerStruc[strucIndex].append(distance)

          residues = self.constraintCount['distance'].addValues(contactTypes,atomTypes,ssTypes,residues)

          if distPerStruc:
            overallDistSum = 0.0
            violated = False

            for strucIndex in range(self.numStructures):
              distSum = 0.0
              for distance in distPerStruc[strucIndex]:
                if strucIndex in self.excludeStructures:
                  continue
                distSum += math.pow(distance,-6)

              if not distSum:
                continue

              #if strucIndex == 0 and dcl.serial == 1 and constraint.serial in (733,):
              #  print constraint.serial
              #  print distSum
              #  print distPerStruc[strucIndex]
              #  print atomCombs

              distance = math.pow(distSum,invSixth)

              if distance > upperDist:
                modelViolations[strucIndex]['upper'][dcl][constraint] = ViolationInfo(distance - upperDist,contactTypes,atomTypes,ssTypes,residues)
                violated = True

              elif distance < lowerDist:
                modelViolations[strucIndex]['lower'][dcl][constraint] = ViolationInfo(lowerDist - distance,contactTypes,atomTypes,ssTypes,residues)
                violated = True

              # This is correct
              overallDistSum += distSum

            if overallDistSum:
              overallDistSum = overallDistSum / self.validStructures
              avgDist = math.pow(overallDistSum,invSixth)

              if avgDist > upperDist:
                overallViolations['upper'][dcl][constraint] = ViolationInfo(avgDist - upperDist,contactTypes,atomTypes,ssTypes,residues)
                violated = True

              elif avgDist < lowerDist:
                overallViolations['lower'][dcl][constraint] = ViolationInfo(lowerDist - avgDist,contactTypes,atomTypes,ssTypes,residues)              
                violated = True

            if violated and not constraint in violatedConstraints[dcl]:
              violatedConstraints[dcl].append(constraint)       
            
    #
    # Now do dihedrals
    #
    
    allDihedralConstraintLists = list(ncs.findAllConstraintLists(className = 'DihedralConstraintList'))

    if allDihedralConstraintLists:
    
      print "--> Found dihedral restraints for analysis"
    
      overallViolations['dihedral'] = {}
      for strucIndex in range(self.numStructures):
        if not modelViolations.has_key(strucIndex):
          modelViolations[strucIndex] = {}
        modelViolations[strucIndex]['dihedral'] = {}

      numDihedralConstraints = 0

      self.constraintCount['dihedral'] = ConstraintCount()

      for dhcl in allDihedralConstraintLists:

        for violType in overallViolations.keys():
          overallViolations['dihedral'][dhcl] = {}
          for strucIndex in modelViolations.keys():
            modelViolations[strucIndex]['dihedral'][dhcl] = {}

        violatedConstraints[dhcl] = []

        dhclConstraints = dhcl.sortedConstraints()

        for constraint in dhclConstraints:

          numDihedralConstraints += 1

          #
          # Get relevant coords and such...
          #

          resonances = constraint.resonances

          allAtomList = []
          dihedralPerStruc = {}

          contactTypes = {}
          atomTypes = {}
          ssTypeCount = {}
          ssTypes = {}  # This used for info only, not in getting contact occurrence!?
          ssTypeList = []
          residues = {}

          for resIndex in range(4):
            resonance = resonances[resIndex]
            if self.resObjectMapping.has_key(resonance):
              (residue,atomList) = self.resObjectMapping[resonance]

              # Only relevant if one atom!
              if len(atomList) != 1:
                continue

              allAtomList.append(atomList[0])

              if residue not in residues.keys():
                residues[residue] = 0
              residues[residue] += 1

              residueKey = self.getResidueKey(residue)
              
              ssKey = 'C'
              if residueKey in self.residueSecStruc.keys():
                ssKey = self.residueSecStruc[residueKey]
  
              if not ssTypeCount.has_key(ssKey):
                ssTypeCount[ssKey] = 0
              ssTypeCount[ssKey] += 1
              ssTypeList.append(ssKey)

          #
          # Only use if all atoms there (and single resonance->atom connection)
          #

          if len(allAtomList) != 4:
            print "  Atoms missing for dihedral angle, only recognized %s - ignored." % str(allAtomList)
            continue

          #
          # Do some reinterpretation of secondary structure type info...
          #

          for ssKey in ssTypeCount.keys():
            if ssTypeCount[ssKey] in (3,4):
              ssType = ssKey
              break
            else:
              ssType = ssTypeList[-1]
            
          if ssType not in ssTypes.keys():
            ssTypes[ssType] = 0
          ssTypes[ssType] += 1

          #
          # Now get dihedral angle per structure
          #

          for strucIndex in range(self.numStructures):

            if strucIndex in self.excludeStructures:
              continue

            angleCoords = []
            for atom in allAtomList:
              if self.coordAtomInfo.has_key(atom):
                coord = self.coordAtomInfo[atom][strucIndex]                
                if coord:
                  angleCoords.append((coord.x,coord.y,coord.z))
                else:
                  angleCoords.append(None)

            if None not in angleCoords and len(angleCoords) == 4:
              angleDegrees = calcTorsionAngleDegrees(angleCoords[0],angleCoords[1],angleCoords[2],angleCoords[3])
              # RETURNS VALUE BETWEEN -180.0 AND 180.0!!!!
              dihedralPerStruc[strucIndex] = angleDegrees

          #
          # Now take the average angle - TODO should really only do this if the circular variance
          # is below a cutoff value, otherwise this doesn't make sense... 
          #

          if len(dihedralPerStruc) > 1:
            (averageAngle,devAngle) = angleStats(dihedralPerStruc.values(),units = 'degrees')
            averageAngle = angleNormaliseAroundZero(averageAngle)
          elif len(dihedralPerStruc) == 1:
            averageAngle = dihedralPerStruc.values()[0]
            devAngle = 0.0

          residues = self.constraintCount['dihedral'].addValues(contactTypes,atomTypes,ssTypes,residues)

          if dihedralPerStruc:

            violated = False

            constraintItems = constraint.sortedItems()

            constraintLimits = []
            # Normalise upper and lower limit to be between 180 and -180
            for constraintItem in constraintItems:
              upperLimit = angleNormaliseAroundZero(constraintItem.upperLimit)
              lowerLimit = angleNormaliseAroundZero(constraintItem.lowerLimit)
              constraintLimits.append((upperLimit,lowerLimit))

            #
            # Structure level violations
            #            

            violated = False

            for strucIndex in range(self.numStructures):
              if strucIndex in self.excludeStructures:
                continue

              angleDegrees = dihedralPerStruc[strucIndex]

              finalViolationDegrees = 0.0

              for (upperLimit,lowerLimit) in constraintLimits:

                violationDegrees = calcAngleViolation(angleDegrees,upperLimit,lowerLimit)

                if violationDegrees:
                  if not finalViolationDegrees or violationDegrees < finalViolationDegrees:
                    finalViolationDegrees = violationDegrees

                else:
                  # Constraint satisfied
                  finalViolationDegrees = violationDegrees
                  break

              if finalViolationDegrees:
                modelViolations[strucIndex]['dihedral'][dhcl][constraint] = ViolationInfo(finalViolationDegrees,contactTypes,atomTypes,ssTypes,residues)
                violated = True

            #
            # Also check ensemble level (though pretty meaningless if wide spread of dihedral angles)
            #

            finalViolationDegrees = 0.0

            for (upperLimit,lowerLimit) in constraintLimits:

              violationDegrees = calcAngleViolation(averageAngle,upperLimit,lowerLimit)

              if violationDegrees:
                if not finalViolationDegrees or violationDegrees < finalViolationDegrees:
                  finalViolationDegrees = violationDegrees

              else:
                # Constraint satisfied
                finalViolationDegrees = violationDegrees
                break

            if finalViolationDegrees:
              overallViolations['dihedral'][dhcl][constraint] = ViolationInfo(finalViolationDegrees,contactTypes,atomTypes,ssTypes,residues)

            #
            # Keep track if constraint violated
            #

            if violated and not constraint in violatedConstraints[dhcl]:
              violatedConstraints[dhcl].append(constraint)       

    #
    # Set up dictionaries to handle violations
    #
    
    specificKeys = InfoTypes().dictInfo.keys()

    self.violTypes = overallViolations.keys()
    self.violTypes.sort()
    
    # Ignore lower for now... not very interesting
    if 'lower' in self.violTypes:
      self.violTypes.pop(self.violTypes.index('lower'))
    
    allViolations = {}
    allModelViolations = {}
    for violType in self.violTypes:
      allViolations[violType] = {}
      allModelViolations[violType] = {}
      for strucIndex in range(self.numStructures):
        allModelViolations[violType][strucIndex] = {}
      
      allViolations[violType]['all'] = []
      for strucIndex in range(self.numStructures):
        allModelViolations[violType][strucIndex]['all'] = []
      
      for specificKey in specificKeys:
        allViolations[violType][specificKey] = {}
        for strucIndex in range(self.numStructures):
          allModelViolations[violType][strucIndex][specificKey] = {}
          
    #
    # Set up info for output
    #
    
    modelRms = {}
    
    self.output['constraintInfo'] = {}
    self.output['ensemble']['constraintLists'] = {}
    self.output['models'] = {}
    
    for strucIndex in range(self.numStructures):
      self.output['models'][strucIndex] = {}
      self.output['models'][strucIndex]['constraintLists'] = {}

    #
    # Set up info for distance constraint lists...
    #
    
    clKeys = {'distance': [], 'dihedral': []}
    clDict = {}
    self.clTypeDict = {}

    for (clList,clType) in ((allDistanceConstraintLists,'distance'),(allDihedralConstraintLists,'dihedral')):
      for cl in clList:
        clKeys[clType].append(cl.serial)
        clDict[cl.serial] = cl
        self.clTypeDict[cl.serial] = clType

      clKeys[clType].sort()

    for clType in clKeys.keys():
      for clKey in clKeys[clType]:

        cl = clDict[clKey]
        constraintsNumber = len(cl.constraints)
        
        if not constraintsNumber:
          print "  Error: %s constraint list %d has no constraints - ignored." % (clType,clKey)
          continue

        self.output['constraintInfo'][clKey] = {}
        
        self.output['ensemble']['constraintLists'][clKey] = {}
        self.output['ensemble']['constraintLists'][clKey]['total'] = constraintsNumber
        self.output['ensemble']['constraintLists'][clKey]['violations'] = {}
        self.output['ensemble']['constraintLists'][clKey]['constraints'] = {}
        self.output['ensemble']['violations'] = {}

        for strucIndex in range(self.numStructures):
          self.output['models'][strucIndex]['constraintLists'][clKey] = {}
          self.output['models'][strucIndex]['constraintLists'][clKey]['total'] = constraintsNumber
          self.output['models'][strucIndex]['constraintLists'][clKey]['violations'] = {}
          self.output['models'][strucIndex]['constraintLists'][clKey]['constraints'] = {}
          self.output['models'][strucIndex]['violations'] = {}


        for violType in self.violTypes:

          # TODO HACK THIS HAS TO BE BETTER!
          if not self.isValidType(violType,clType):
            continue
            
          violationNumber = len(overallViolations[violType][cl])  
          self.output['ensemble']['constraintLists'][clKey]['violations'][violType] = {

                'number': violationNumber,
                'percent': violationNumber * 100.0 / constraintsNumber

                }

          for strucIndex in range(self.numStructures):
            violationNumber = len(modelViolations[strucIndex][violType][cl])  
            self.output['models'][strucIndex]['constraintLists'][clKey]['violations'][violType] = {

                'number': violationNumber,
                'percent': violationNumber * 100.0 / constraintsNumber

                }

          consKeys = []
          consDict = {}

          for cons in violatedConstraints[cl]:
            consKeys.append(cons.serial)
            consDict[cons.serial] = cons

          consKeys.sort()

          for consKey in consKeys:
            cons = consDict[consKey]

            violation = None
            numViolations = 0
            hasLargeViolation = False

            if overallViolations[violType][cl].has_key(cons):

              overallViolation = overallViolations[violType][cl][cons].setInfo(allViolations[violType])

              if overallViolation:
                if not self.output['ensemble']['constraintLists'][clKey]['constraints'].has_key(consKey):
                  self.output['ensemble']['constraintLists'][clKey]['constraints'][consKey] = {}
                self.output['ensemble']['constraintLists'][clKey]['constraints'][consKey][violType] = overallViolation

            for strucIndex in range(self.numStructures):
              if strucIndex in self.excludeStructures:
                continue
              if modelViolations[strucIndex][violType][cl].has_key(cons):

                violation = modelViolations[strucIndex][violType][cl][cons].setInfo(allModelViolations[violType][strucIndex])

                if violation:
                  numViolations += 1
                  if not self.output['models'][strucIndex]['constraintLists'][clKey]['constraints'].has_key(consKey):
                    self.output['models'][strucIndex]['constraintLists'][clKey]['constraints'][consKey] = {}
                  self.output['models'][strucIndex]['constraintLists'][clKey]['constraints'][consKey][violType] = violation

                  if violation > self.largeViolations[clType]:
                    hasLargeViolation = True

            if violation:
            
              if clType == 'distance':

                if not self.output['constraintInfo'][clKey].has_key(consKey):

                  # Set generic constraint info only once...

                  self.output['constraintInfo'][clKey][consKey] = {

                    'contactTypes': ','.join(contactTypes),
                    'atomTypes':    str(atomTypes),
                    'items':        [],

                  }

                  for item in cons.items:

                    resonances = item.orderedResonances

                    # Can't count on above returning something...
                    if not resonances:
                      resonances = item.sortedResonances()

                    resTexts = self.getResonanceTexts(resonances)

                    #
                    # Set contact occurrence for each item if required...
                    # 
                    
                    contactOccurrence = 0.0
                    
                    if self.useContactOccurrence:

                      resInfo = []

                      for i in range(2):
                        resonance = resonances[i]
                        if self.resMapping.has_key(resonance):
                          (residue,atomNameTuple) = self.resMapping[resonance]
                          # Can happen
                          if atomNameTuple.count(None):
                            continue
                          resInfo.append(self.residueInfo[residue] + (atomNameTuple,))

                      if len(resInfo) == 2:
                        ssCodes = []
                        
                        # Note: was originally self.residueSecStrucDict.
                        if self.residueSecStruc:
                          for residueKey in (resInfo[0][1],resInfo[1][1]):
                            ssCode = None
                            if self.residueSecStruc.has_key(residueKey):
                              ssCode = self.residueSecStruc[residueKey]
                            ssCodes.append(ssCode)
                            
                        ssCodes = tuple(ssCodes)

                        # Or should getContactOccurrence be a separate bit of code after all?
                        (contactOccurrence,atomNameTuple,averageDist) = self.getContactOccurrence(resInfo[:], ssCodes = ssCodes, contactOccurrenceDefault = 0.0)
                      
                    self.output['constraintInfo'][clKey][consKey]['items'].append((resTexts,contactOccurrence))
                      
                  self.output['constraintInfo'][clKey][consKey]['items'].sort()
                  
                #
                # Set violation information info
                #

                self.output['constraintInfo'][clKey][consKey][violType] = {

                  'limit':             getattr(cons,"%sLimit" % violType),
                  'numViolations':     numViolations,
                  'hasLargeViolation': hasLargeViolation

                }
                
              elif clType == 'dihedral':

                if not self.output['constraintInfo'][clKey].has_key(consKey):

                  # Set generic constraint info only once...

                  self.output['constraintInfo'][clKey][consKey] = {'atoms': []}

                  resonances = cons.resonances

                  resTexts = self.getResonanceTexts(resonances)

                  self.output['constraintInfo'][clKey][consKey]['atoms'] = resTexts

                #
                # Set violation information info
                #

                self.output['constraintInfo'][clKey][consKey][violType] = {

                  'limits':            [],
                  'numViolations':     numViolations,
                  'hasLargeViolation': hasLargeViolation

                }
                
                
                for item in cons.items:
                  self.output['constraintInfo'][clKey][consKey][violType]['limits'].append((item.upperLimit,item.lowerLimit))
                self.output['constraintInfo'][clKey][consKey][violType]['limits'].sort()
                
    #
    # Analyse overall violations
    #
    
    for clType in clKeys.keys():
    
      for violType in self.violTypes:
      
        if violType == 'dihedral':
          tempNum = numDihedralConstraints
        else:
          tempNum = numConstraints

        if not self.isValidType(violType,clType):
          continue

        self.output['ensemble']['violations'][violType] = {}      
        numViolations = len(allViolations[violType]['all'])
        self.output['ensemble']['violations'][violType]['all'] = self.getListStats((allViolations[violType]['all'],numViolations,numViolations),(tempNum,tempNum))

        for infoTypeKey in self.constraintCount[clType].infoDict.keys():
          self.output['ensemble']['violations'][violType][infoTypeKey] = self.getSpecificListStats(self.constraintCount[clType].infoDict[infoTypeKey],allViolations[violType][infoTypeKey],infoTypeKey)

        for strucIndex in range(self.numStructures):
          if strucIndex in self.excludeStructures or 'violations' not in self.output['models'][strucIndex].keys():
            continue
          self.output['models'][strucIndex]['violations'][violType] = {}
          numViolations = len(allModelViolations[violType][strucIndex]['all'])
          self.output['models'][strucIndex]['violations'][violType]['all'] = self.getListStats((allModelViolations[violType][strucIndex]['all'],numViolations,numViolations),(tempNum,tempNum))

          # Track RMS for sorting structure order in violation overview
          if violType == 'upper':
            rmsValue = self.output['models'][strucIndex]['violations'][violType]['all'][3]
            if not modelRms.has_key(rmsValue):
              modelRms[rmsValue] = []
            modelRms[rmsValue].append(strucIndex)

          for infoTypeKey in self.constraintCount[clType].infoDict.keys():
            self.output['models'][strucIndex]['violations'][violType][infoTypeKey] = self.getSpecificListStats(self.constraintCount[clType].infoDict[infoTypeKey],allModelViolations[violType][strucIndex][infoTypeKey],infoTypeKey)
     
    #
    # Sort strucs by RMS
    #
    
    self.output['modelOrder'] = []
    
    rmsValues = modelRms.keys()
    rmsValues.sort()
    
    for rmsValue in rmsValues:
      self.output['modelOrder'].extend(modelRms[rmsValue])
    
    self.getStatisticsExecuted = True

  def isValidType(self,violType,clType):            

    if violType == clType == 'dihedral':
      validType = True
    elif clType == 'distance' and violType != 'dihedral':
      validType = True
    else:
      validType = False
    
    return validType 

  def getResonanceTexts(self,resonances):

    resTexts = []

    for resonance in resonances:
      if resonance in self.resObjectMapping.keys():
        (residue,atomList) = self.resObjectMapping[resonance]
        atomSets = {}
        atomNames = []

        for atom in atomList:
          atomNames.append(atom.name)
          if atom.chemAtom.chemAtomSet:
            chemAtomSet = atom.chemAtom.chemAtomSet
            atomSetName = chemAtomSet.name
            if not atomSets.has_key(atomSetName):
              atomSets[atomSetName] = [[],len(chemAtomSet.chemAtoms)]
            atomSets[atomSetName][0].append(atom.name)

        for atomSetName in atomSets.keys():
          # Only use atomSet name if all atoms found!
          if len(atomSets[atomSetName][0]) == atomSets[atomSetName][1]:
            for atomName in atomSets[atomSetName][0]:
              atomNames.pop(atomNames.index(atomName))
            atomNames.append(atomSetName)

        atomNames.sort()
        resTexts.append("%s.%d.%s" % (residue.chain.code,residue.seqCode,string.join(atomNames,',')))

      else:
        resTexts.append("resonance %s" % resonance.name)

    return resTexts

  def getInfoTypeText(self,infoTypeKey,infoTypeValue):

    #if infoTypeKey in ('contact',):
    #  infoTypeText = string.join(infoTypeValue,',')

    if infoTypeKey in ('contact','secStruc',):
      infoTypeText = infoTypeValue

    elif infoTypeKey == 'atom':
      atomTypeTexts = []
      for atomType in infoTypeValue:
        atomTypeTexts.append(string.join(atomType,'-'))
      infoTypeText = string.join(atomTypeTexts,',')

    elif infoTypeKey == 'residues':
      infoTypeText = "%s.%d.%s" % (infoTypeValue[0],infoTypeValue[1],infoTypeValue[2])

    return infoTypeText
      
  def getListStats(self,valueList,constrValues):
  
    from rpy import r
    
    if valueList[0]:
      
      contribCount = constrValues[0]
      totalCount = constrValues[1]
      
      violList = valueList[0]
      violContribCount = valueList[1]
      violTotalCount = valueList[2]
      
      statSummary = r.summary(violList)
      mean = statSummary['Mean']
      median = statSummary['Median']
      
      rms = getRms(valueList[0],total = totalCount)

      percent = violContribCount * 100.0 / contribCount
      infoValues = (violContribCount,contribCount,percent,rms,mean,median)

    else:
      infoValues = (0.0,constrValues[0],0.0,0.0,0.0,0.0)
      
    return infoValues

  def getSpecificListStats(self,constraintCountInfoDict,violationDict,infoTypeKey):
  
    infoValuesList = []

    infoTypeValueList = constraintCountInfoDict.keys()
    infoTypeValueList.sort()

    for infoTypeValue in infoTypeValueList:
    
      if violationDict.has_key(infoTypeValue):
        violationList = violationDict[infoTypeValue]
      else:
        violationList = ([],None)
    
      infoValues = self.getListStats(violationList,constraintCountInfoDict[infoTypeValue])
      infoTypeText = self.getInfoTypeText(infoTypeKey,infoTypeValue)
      
      infoValuesList.append((infoTypeText,infoValues))
      
    return infoValuesList

  def getSummaryText(self):

    self.getStatistics(False)

    dclKeys = self.output['ensemble']['constraintLists'].keys()
    dclKeys.sort()

    finalData = {}

    for dclKey in dclKeys:
      clType = self.clTypeDict[dclKey]

      # Only handle distance at the moment

      if clType != 'distance':
        continue

      consKeys = self.output['constraintInfo'][dclKey].keys()
      consKeys.sort()

      if dclKey not in finalData:
        finalData[dclKey] = {}

      for violType in self.violTypes:
        if clType == 'distance' and violType == 'dihedral':
          continue
        elif clType == 'dihedral' and violType == 'upper':
          continue

        #print 'DCL: [%s] CLTYPE: [%s] VIOLTYPE: [%s]\n' % (dclKey, clType, violType)

        for consKey in consKeys:
          #if consKey not in finalData[dclKey]:
          #  finalData[dclKey][consKey] = None

          #print "%4d " % consKey,
          consInfo = self.output['constraintInfo'][dclKey][consKey]

          sumVal   = 0.0
          sumSqVal = 0.0
          violMax  = 0.0

          count   = 0
          count_0 = 0
          count_1 = 0
          count_3 = 0
          count_5 = 0

          for strucIndex in self.output['modelOrder']:
            viol = self.getViolationData(self.output['models'][strucIndex]['constraintLists'][dclKey]['constraints'],consKey,violType,clType)

            if viol > violMax:
              violMax = viol

            sumVal   += viol
            sumSqVal += viol*viol

            count += 1
            if viol > 0.5:
              count_5 += 1

            if viol > 0.3:
              count_3 += 1

            if viol > 0.1:
              count_1 += 1

            if viol:
              count_0 += 1

          viol = self.getViolationData(self.output['ensemble']['constraintLists'][dclKey]['constraints'],consKey,violType,clType)

          #if viol:
          #  print "%.3f " % viol,
          #else:
          #  print "     -",

          ovMean = sumVal/count
          #violMean = sumVal/count_0

          sd = math.sqrt( (count*sumSqVal - sumVal*sumVal)/(count * (count-1) ) )
          #print "%.3f %.3f %.3f %3d %3d %3d %3d %3d" % (ovMean, sd, violMax, count, count_0, count_1, count_3, count_5)

          finalData[dclKey][consKey] = ["%.3f" % viol, "%.3f" % ovMean, "%.3f" % sd, "%.3f" % violMax, count, count_0, count_1, count_3, count_5]

    return finalData

  def getViolationData(self,consDict,consKey,violType,clType):

    violation = 0.000
  
    if consDict.has_key(consKey) and consDict[consKey].has_key(violType):
      violation = consDict[consKey][violType]
      #print "%.3f " % violation,

    return violation
  
  ######################
  # File/screen output #
  ######################
  
  def writeText(self,fout,verbose = True, testMode = False):
  
    self.getStatistics(testMode)

    #
    # TODO: some of this info could be considered generic, no?
    #

    infoTagFormat = "%6s  "
    infoFormats = ("  %5.1f  ","  %5.1f  "," %5.1f%% "," %6.3f "," %6.3f "," %6.3f ") 

    typeTagFormat = "      %-15s:"
    
    self.fout = fout
  
    #
    # Below is screen/file output - TODO replace print by fout.write()!!
    #
    
    for violType in self.violTypes:

      fout.write("\n")
      fout.write(drawBox("%s bound violations overview" % (violType.capitalize())))
      fout.write("\n" * 2)
      
      dclKeys = self.output['ensemble']['constraintLists'].keys()
      dclKeys.sort()
      
      for dclKey in dclKeys:
      
        clType = self.clTypeDict[dclKey]
        
        if not self.isValidType(violType,clType):
          continue
      
        dclDict = self.output['ensemble']['constraintLists'][dclKey]

        fout.write("  Constraint list %d (%d constraints, ensemble violations %d)" % (dclKey,dclDict['total'],dclDict['violations'][violType]['number']))
        fout.write("\n" * 2)
        
        fout.write("    %-12s" % "Model:")
        for strucIndex in range(self.numStructures):
          fout.write("  %3d  " % (strucIndex + 1))

        fout.write("\n" * 2)

        fout.write("    %-12s" % "Violations:")
        for strucIndex in range(self.numStructures):
          violDict = self.output['models'][strucIndex]['constraintLists'][dclKey]['violations'][violType]
          fout.write("  %3d  " % (violDict['number']))
        
        # Track large/small violations
        numLargeViolations = []
        numSmallViolations = []
        for strucIndex in range(self.numStructures):
          #self.output['ensemble']['constraintLists'][dclKey]['constraints']
          numLarge=numSmall=0
          for consKey in self.output['models'][strucIndex]['constraintLists'][dclKey]['constraints'].keys():
            violation = self.output['models'][strucIndex]['constraintLists'][dclKey]['constraints'][consKey][violType]
            if violation > self.largeViolations[clType]:
              numLarge += 1
            elif violation > self.smallViolations[clType]:
              numSmall += 1
          numLargeViolations.append(numLarge)
          numSmallViolations.append(numSmall)
          

        fout.write("\n")
        fout.write("    %-12s" % "    Large:")
        for strucIndex in range(self.numStructures):          
          fout.write("  %3d  " % numLargeViolations[strucIndex])
          
        fout.write("\n")
        fout.write("    %-12s" % "    Small:")
        for strucIndex in range(self.numStructures):          
          fout.write("  %3d  " % numSmallViolations[strucIndex])


        fout.write("\n" * 2)
        fout.write("    %-12s " % "Percent:")
        for strucIndex in range(self.numStructures):
          violDict = self.output['models'][strucIndex]['constraintLists'][dclKey]['violations'][violType]
          fout.write(" %5.1f%%" % (violDict['percent']))

        fout.write("\n" * 2)

        if verbose:
        
          consKeys = dclDict['constraints'].keys()
          consKeys.sort()
          
          for consKey in consKeys:
          
            if dclDict['constraints'][consKey].has_key(violType):
              
              overallViolation = dclDict['constraints'][consKey][violType]
              consInfo = self.output['constraintInfo'][dclKey][consKey]
 
              fout.write("    %d: %.3f  (%s  %s)\n" % (consKey,overallViolation,consInfo['contactTypes'],consInfo['atomTypes']))

              for (resTexts,contactOccurrence) in consInfo['items']:
                if self.useContactOccurrence:
                  addText = " (%6.3f)" % contactOccurrence
                else:
                  addText = ""
                
                fout.write("      %-20s - %-20s%s\n" % (resTexts[0],resTexts[1],addText))

              fout.write("\n")

          sys.__stdout__.flush()

        fout.write("\n")

      fout.write("\n")

      fout.write("  Total violations overview:")
      fout.write("\n" * 2)

      fout.write("    Ensemble level:")
      fout.write("\n" * 2)

      if self.output['ensemble']['violations'][violType]['all'][0]:

        # Print header
        fout.write(typeTagFormat % "")
        for infoTag in self.infoTags:
          fout.write(infoTagFormat % infoTag)
 
        fout.write("\n" * 2)

        self.outputValues(typeTagFormat,"Total",self.output['ensemble']['violations'][violType]['all'],infoFormats)

        fout.write("\n")

        for infoTypeKey in self.infoTypeKeyList:
          # These are not relevant for dihedral constraints
          if violType == 'dihedral' and infoTypeKey in ('contact','atom'):
            continue

          self.outputSpecificValues(self.output['ensemble']['violations'][violType][infoTypeKey],typeTagFormat,infoFormats)
          
        fout.write("\n")

      else:

        print "    Nothing to report - all good."

      fout.write("\n")
      fout.write("    Per model level:")
      fout.write("\n" * 2)

      # Print header
      fout.write(typeTagFormat % "")
      for infoTag in self.infoTags:
        fout.write(infoTagFormat % infoTag)

      fout.write("\n" * 2)
      
      self.minMaxValues = {'min': [None] * (len(infoFormats) ), 'max': [None] * (len(infoFormats) )}

      for strucIndex in range(self.numStructures):
        self.outputValues(typeTagFormat,"Model %d" % strucIndex,self.output['models'][strucIndex]['violations'][violType]['all'],infoFormats,trackMinMax = True)

      fout.write("\n")

      self.outputValues(typeTagFormat,"Minimum",self.minMaxValues['min'],infoFormats)
      self.outputValues(typeTagFormat,"Maximum",self.minMaxValues['max'],infoFormats)

      fout.write("\n" * 2)

  def outputValues(self,typeTagFormat,tag,infoValues,infoFormats,trackMinMax = False):

    self.fout.write(typeTagFormat % tag)

    for i in range(len(infoValues)):
      self.fout.write(infoFormats[i] % infoValues[i])
      
      if trackMinMax:
        for trackType in ('min','max'):
          if self.minMaxValues[trackType][i] == None:
            self.minMaxValues[trackType][i] = infoValues[i]
          else:
            if trackType == 'min' and self.minMaxValues[trackType][i] > infoValues[i]:
              self.minMaxValues[trackType][i] = infoValues[i]
            elif trackType == 'max' and self.minMaxValues[trackType][i] < infoValues[i]:
              self.minMaxValues[trackType][i] = infoValues[i]
      
    self.fout.write("\n")

  def outputSpecificValues(self,valuesList,typeTagFormat,infoFormats):

    for (infoTypeText,infoValues) in valuesList:
      self.outputValues(typeTagFormat,infoTypeText,infoValues,infoFormats)

    self.fout.write("\n")



  #####################
  # HTML PAGES output #
  #####################
  
  def getFullSideMenuList(self,sideMenuItems):
    
    fullSideMenuList = []
    
    for sideMenuItem in sideMenuItems:
      fullSideMenuList.append((sideMenuItem,'%s.html' % sideMenuItem))
      fullSideMenuList.append([])

      subMenuList = [

        ("Ensemble analysis","%s/ensemble.html"),
        ("Model analysis","%s/models.html"),

      ]

      for (subMenuItem,subMenuLink) in subMenuList:
        fullSideMenuList[-1].append((subMenuItem,subMenuLink % sideMenuItem))
 
    return fullSideMenuList

  def writeHtml(self,saveDir, sideMenuItems = None, testMode = False):
  
    self.getStatistics(testMode)

    self.infoFormats = ("%5.1f","%5.1f","%5.1f%%","%6.3f","%6.3f","%6.4f") 
  
    self.baseName = self.dataType
    
    if not sideMenuItems:
      sideMenuItems = []
 
    if not self.baseName in sideMenuItems:
      sideMenuItems.append(self.baseName)

    self.fullSideMenuList = self.getFullSideMenuList(sideMenuItems)
    
    #
    # Directory creation and handling...
    #
    
    saveDir = os.path.join(saveDir,'html')
    if not os.path.exists(saveDir):
      os.mkdir(saveDir)
    
    subDirName = self.projectName
    
    htmlDir = os.path.join(saveDir,subDirName)
    if not os.path.exists(htmlDir):
      os.mkdir(htmlDir)
    
    self.mainSaveDir = htmlDir  # TODO this to be used further on for subpages
    
    htmlDetailsDir = os.path.join(htmlDir,self.baseName)
    if not os.path.exists(htmlDetailsDir):
      os.mkdir(htmlDetailsDir)
    
    #
    # Set up graph stuff
    #

    self.graphDir = "graphs"
    self.xValues = range(0,self.validStructures + 1)
    
    #
    # Write out main overview page
    #
      
    self.createMainPage(htmlDir)    

    #
    # Details for ensemble and models...
    #
    
    self.createDetailsPage(htmlDetailsDir,'ensemble',self.output['ensemble'])
    
    #
    # Info per model
    #
    
    modelNames = []
    for strucIndex in range(self.numStructures):
      modelName = 'model_%d' % (strucIndex + 1)
      modelNames.append(modelName)
      if strucIndex not in self.excludeStructures:
        self.createDetailsPage(htmlDetailsDir,modelName,self.output['models'][strucIndex])

    self.createModelTopPage(htmlDetailsDir,'models',modelNames)
    
    #
    # Info per constraint list
    #
    
    constraintListNames = []
 
    dclKeys = self.output['ensemble']['constraintLists'].keys()
    dclKeys.sort()

    for dclKey in dclKeys:
     constraintListName = "constraintList_%d" % dclKey
     self.createConstraintListPage(htmlDetailsDir,constraintListName,dclKey)
     constraintListNames.append(constraintListName)


  def createMainPage(self,htmlDir):
    
    #
    # Initialize (generic)
    #
    
    localGraphDir = self.getLocalGraphDir(htmlDir)

    htmlBaseName = self.getHtmlBaseName(htmlDir)

    #
    # Initialize (non generic)
    #

    ignoreInfoTags = ('Total')  # This ignored for totals - is pointless to display since same value over and over again

    #
    # Create pages
    # 

    mainPage = self.HtmlPage(os.path.join(htmlDir,"%s.html" % (self.baseName)), htmlBaseName = "Void", styleSheet = self.styleSheet)

    mainPage.setupHtml("Violation analysis for %s, data type '%s'" % (self.projectName,self.dataType))
    mainPage.htmlBaseLevel = 0 # Hack! TODO get this sorted!

    mainPage.createFullSideMenu(self.fullSideMenuList, activeItem = (self.baseName,))

    mainPage.mainPageHtml("Summary of all violations")
    
    
    for violType in self.violTypes:
    
      graphInfo = self.initialiseGraphInfo(violType,ignoreInfoTags,self.baseName,localGraphDir)

      mainPage.setEmptyRow(colspan = 99)

      # Top level row
      mainPage.writeTableHeader('%s limit violations (out of %d constraints)' % (violType.capitalize(),self.output['ensemble']['violations'][violType]['all'][1]), colspan = 99)
      
      # Write column header info
      mainPage.addMainTableRow()
      mainPage.addMainTableColumn("Level",addText = ' class="subheading"')
      for infoTag in self.infoTags:
        if infoTag not in ignoreInfoTags:
          graphLink = graphInfo['links'][infoTag]
          mainPage.addMainTableColumn("<a href='%s'>%s</a>" % (graphLink,infoTag),addText = ' class="subheading"')
      mainPage.closeMainTableRow()

      mainPage.addMainTableRow()
      mainPage.addMainTableColumn("<a href='%s'>Ensemble</a>" % (os.path.join(self.baseName,'ensemble.html')), addText = ' class="subheading"')

      # Write data values for ensemble
      infoValues = self.output['ensemble']['violations'][violType]['all']
      for j in range(len(infoValues)):
        infoTag = self.infoTags[j]
        if infoTag not in ignoreInfoTags:
          mainPage.addMainTableColumn(self.infoFormats[j] % infoValues[j], addText = ' class="centered"')
          graphInfo['data'][infoTag].append(infoValues[j])
      mainPage.closeMainTableRow()
      
      mainPage.setEmptyRow(colspan = 99)

      # Write data values for models
      for strucIndex in range(self.numStructures):
      
        strucNum = strucIndex + 1

        mainPage.addMainTableRow()

        if strucIndex not in self.excludeStructures:
          
          #
          # Structure included in analysis
          #
          
          mainPage.addMainTableColumn("<a href='%s'>Model %s</a>" % (os.path.join(self.baseName,'model_%d.html' % strucNum),strucNum), addText = ' class="subheading"')

          infoValues = self.output['models'][strucIndex]['violations'][violType]['all']
          for j in range(len(infoValues)):
            infoTag = self.infoTags[j]
            if infoTag not in ignoreInfoTags:
              mainPage.addMainTableColumn(self.infoFormats[j] % infoValues[j], addText = ' class="centered"')
              graphInfo['data'][infoTag].append(infoValues[j])
        
        else:
          
          #
          # Structure not included in analysis
          #
          
          mainPage.addMainTableColumn("Model %s" % strucNum, addText = ' class="subheading"')

          infoValues = ("Not included",) + ("n/a",) * (len(self.infoTags) - 1)
          for j in range(len(self.infoTags)):
            infoTag = self.infoTags[j]
            if infoTag not in ignoreInfoTags:
              mainPage.addMainTableColumn("%s" % infoValues[j], addText = ' class="centered"')
        
        mainPage.closeMainTableRow()
      
      #
      # Create graphs...
      #
      
      for infoTag in graphInfo['data'].keys():
        title = "%s values for %s (%s limit)" % (infoTag,self.baseName,violType)
        self.writePlot(graphInfo['fileNames'][infoTag],self.xValues,graphInfo['data'][infoTag],title,"Models (0 is ensemble)",infoTag)
    
    mainPage.closeMainTable()
    mainPage.finishHtml()    

  def createDetailsPage(self,htmlDetailsDir,htmlPageName,infoDict):

    #
    # Initialize (generic)
    #
    
    localGraphDir = self.getLocalGraphDir(htmlDetailsDir)

    htmlBaseName = self.getHtmlBaseName(htmlDetailsDir)

    #
    # Initialize (non generic)
    #

    ignoreInfoTags = tuple()
    ignoreInfoTypeKeyValues = ('residues',)
    ignoreInfoTypeKeyGraphs = ('contact','atom','secStruc')

    #
    # Create page
    # 

    htmlPage = self.HtmlPage(os.path.join(htmlDetailsDir,"%s.html" % (htmlPageName)), htmlBaseName = htmlBaseName, styleSheet = self.styleSheet)

    htmlPage.setupHtml("Violation analysis for %s, data type '%s': %s level" % (self.projectName,self.dataType,htmlPageName))

    htmlPage.createFullSideMenu(self.fullSideMenuList, activeItem = (self.baseName,))

    htmlPage.mainPageHtml("Details of violations")


    for violType in self.violTypes:

      htmlPage.setEmptyRow(colspan = 99)

      # Top level row
      htmlPage.writeTableHeader('%s limit violations' % (violType.capitalize()), colspan = 99)
      
      #
      # Write constraint list level info
      #

      # Write column info
      htmlPage.addMainTableRow()
      htmlPage.addMainTableColumn("Constraint list",addText = ' class="subheading"')
      for infoTag in ('Number','Percent','Total'):
        htmlPage.addMainTableColumn(infoTag,addText = ' class="subheading"')
      htmlPage.closeMainTableRow()

      dclKeys = infoDict['constraintLists'].keys()
      dclKeys.sort()
      
      for dclKey in dclKeys:
      
        clType = self.clTypeDict[dclKey]
        
        if not self.isValidType(violType,clType):
          continue
      
        constraintListName = "constraintList_%d" % dclKey
        dclDict = infoDict['constraintLists'][dclKey]
        
        htmlPage.addMainTableRow()
        htmlPage.addMainTableColumn("<a href='%s.html'>%s</a>" % (constraintListName,dclKey), addText = ' class="subheading"')
        
        
        htmlPage.addMainTableColumn("%5.1f" % dclDict['violations'][violType]['number'], addText = ' class="centered"')
        htmlPage.addMainTableColumn("%5.1f%%" % dclDict['violations'][violType]['percent'], addText = ' class="centered"')
        htmlPage.addMainTableColumn("%5.1f" % dclDict['total'], addText = ' class="centered"')

        htmlPage.closeMainTableRow()
     
      #
      # Write information type level info
      #

      for mainInfoTypeKey in self.infoTypeKeyList:
      
        # These are not relevant for dihedral constraints
        if violType == 'dihedral' and mainInfoTypeKey in ('contact','atom'):
          continue
        
        mainInfoTypeText = mainInfoTypeKey.capitalize()
        
        if mainInfoTypeKey not in ignoreInfoTypeKeyGraphs:
          graphBaseName = "%s_%s" % (self.baseName,mainInfoTypeKey)
          graphInfo = self.initialiseGraphInfo(violType,ignoreInfoTags,graphBaseName,localGraphDir)

        htmlPage.setEmptyRow(colspan = 99)

        # Write column header info
        htmlPage.addMainTableRow()
        htmlPage.addMainTableColumn(mainInfoTypeText,addText = ' class="subheading"')
        for infoTag in self.infoTags:
          if infoTag not in ignoreInfoTags:
            if mainInfoTypeKey not in ignoreInfoTypeKeyGraphs:
              graphLink = graphInfo['links'][infoTag]
              infoTagHtml = "<a href='%s'>%s</a>" % (graphLink,infoTag)
            else:
              infoTagHtml = infoTag
            htmlPage.addMainTableColumn(infoTagHtml,addText = ' class="subheading"')
        htmlPage.closeMainTableRow()

        # Write data values for ensemble
        valuesList = infoDict['violations'][violType][mainInfoTypeKey]
        
        for (infoTypeText,infoValues) in valuesList:
          if mainInfoTypeKey not in ignoreInfoTypeKeyValues:
            htmlPage.addMainTableRow()
            htmlPage.addMainTableColumn(infoTypeText, addText = ' class="subheading"')
          for j in range(len(infoValues)):
            infoTag = self.infoTags[j]
            if infoTag not in ignoreInfoTags:
              if mainInfoTypeKey not in ignoreInfoTypeKeyValues:
                htmlPage.addMainTableColumn(self.infoFormats[j] % infoValues[j], addText = ' class="centered"')
              if mainInfoTypeKey not in ignoreInfoTypeKeyGraphs:
                graphInfo['data'][infoTag].append(infoValues[j])
          htmlPage.closeMainTableRow()

        #
        # Create graphs...
        #

        if mainInfoTypeKey not in ignoreInfoTypeKeyGraphs:
          for infoTag in graphInfo['data'].keys():
            title = "%s values for %s for %s data (%s limit)" % (infoTag,self.baseName,mainInfoTypeText,violType)
            self.writePlot(graphInfo['fileNames'][infoTag],range(1,len(valuesList)+1),graphInfo['data'][infoTag],title,mainInfoTypeText,infoTag,lab = (int(len(valuesList)/10),10,2))

    htmlPage.closeMainTable()
    htmlPage.finishHtml()    

  def createModelTopPage(self,htmlDetailsDir,htmlPageName,modelNames):

    #
    # Initialize (generic)
    #

    htmlBaseName = self.getHtmlBaseName(htmlDetailsDir)

    #
    # Create page
    # 

    htmlPage = self.HtmlPage(os.path.join(htmlDetailsDir,"%s.html" % (htmlPageName)), htmlBaseName = htmlBaseName, styleSheet = self.styleSheet)

    htmlPage.setupHtml("Violation analysis for %s, data type '%s': links to models" % (self.projectName,self.dataType))

    htmlPage.createFullSideMenu(self.fullSideMenuList, activeItem = (self.baseName,))

    htmlPage.mainPageHtml("List of models")

    for modelName in modelNames:

      htmlPage.addMainTableRow()
      htmlPage.addMainTableColumn("<a href='%s.html'>%s</a>" % (modelName,modelName),addText = ' class="centered"')
      htmlPage.closeMainTableRow()

    htmlPage.closeMainTable()
    htmlPage.finishHtml()    
  

  def createConstraintListPage(self,htmlDetailsDir,htmlPageName,dclKey):

    #
    # Initialize (generic)
    #
    
    #localGraphDir = self.getLocalGraphDir(htmlDetailsDir)
    htmlBaseName = self.getHtmlBaseName(htmlDetailsDir)

    #
    # Initialize (non generic)
    #
    
    clType = self.clTypeDict[dclKey]

    #
    # Create page
    # 

    htmlPage = self.HtmlPage(os.path.join(htmlDetailsDir,"%s.html" % (htmlPageName)), htmlBaseName = htmlBaseName, styleSheet = self.styleSheet)

    htmlPage.setupHtml("Violation analysis for %s, data type '%s': %s" % (self.projectName,self.dataType,htmlPageName))

    htmlPage.createFullSideMenu(self.fullSideMenuList, activeItem = (self.baseName,))

    htmlPage.mainPageHtml("Details of violations")

    for violType in self.violTypes:
    
      if not self.isValidType(violType,clType):
        continue

      htmlPage.setEmptyRow(colspan = 99)

      # Top level row
      htmlPage.writeTableHeader('%s limit violations' % (violType.capitalize()), colspan = 99)
      
      #
      # Write constraint list level info
      #

      # Write column info
      htmlPage.addMainTableRow()
      htmlPage.addMainTableColumn("Constraint",addText = ' class="subheading" rowspan=2')
      
      if clType == 'distance':
        htmlPage.addMainTableColumn("Limit",addText = ' class="subheading" rowspan=2')
        htmlPage.addMainTableColumn("Items",addText = ' class="subheading" rowspan=2 colspan=2')
        if self.useContactOccurrence:
          htmlPage.addMainTableColumn("Occurrence",addText = ' class="subheading" rowspan=2 colspan=1')
        
      elif clType == 'dihedral':  
        htmlPage.addMainTableColumn("Atoms",addText = ' class="subheading" rowspan=2')
        htmlPage.addMainTableColumn("Limits",addText = ' class="subheading" rowspan=2 colspan=2')
      
      htmlPage.addMainTableColumn("<a href='ensemble.html'>Ensemble</a>",addText = ' class="subheading" rowspan=2')
      htmlPage.addMainTableColumn("Models",addText = ' class="subheading" colspan=%d' % self.validStructures)
      htmlPage.closeMainTableRow()

      htmlPage.addMainTableRow()

      for strucIndex in self.output['modelOrder']:
        if strucIndex in self.excludeStructures:
          continue
        modelIndex = strucIndex + 1
        htmlPage.addMainTableColumn("<a href='%s'>%s</a>" % ("model_%d.html" % modelIndex,modelIndex),addText = ' class="subheading"')
      htmlPage.closeMainTableRow()


      # TODO might want to put this info back in here...
      #dclDict = 

      #numConstraints = dclDict['total']
      #numViolations = dclDict['violations'][violType]

      #for strucIndex in range(self.numStructures):
      #  dclDictModel = self.output['models'][strucIndex]['constraintLists'][dclKey]
      #  numViolations = dclDictModel['violations'][violType]

      consKeys = self.output['constraintInfo'][dclKey].keys()
      consKeys.sort()
          
      for consKey in consKeys:
        consInfo = self.output['constraintInfo'][dclKey][consKey]

        fractionViolated = consInfo[violType]['numViolations'] * 1.0 / self.validStructures
        if consInfo[violType]['hasLargeViolation'] or fractionViolated > 0.8:
          (prefix,suffix) = ('<font color="#FF0000">',"</font>")
        else:
          (prefix,suffix) = ("","")

        htmlPage.addMainTableRow()
        htmlPage.addMainTableColumn(consKey, addText = ' class="subheading"')
        
        if clType == 'distance':
          htmlPage.addMainTableColumn("%s%.2f%s" % (prefix,consInfo[violType]['limit'],suffix), addText = ' class="centered"')

          htmlResTexts = [[],[]]
          occurrenceTexts = []
          for (resTextList,contactOccurrence) in consInfo['items']:
            occurrenceTexts.append("%6.3f" % contactOccurrence)
            for i in range(len(resTextList)):
              htmlResTexts[i].append("%s%s%s" % (prefix,resTextList[i],suffix))

          for i in range(2):
            htmlResText = '<BR>'.join(htmlResTexts[i])
            htmlPage.addMainTableColumn(htmlResText, addText = ' class="centered"')
          
          if self.useContactOccurrence:
            htmlOccText = "<BR>".join(occurrenceTexts)
            htmlPage.addMainTableColumn(htmlOccText, addText = ' class="centered"')
        
        elif clType == 'dihedral':
          htmlPage.addMainTableColumn("%s%s%s" % (prefix,string.join(consInfo['atoms'],'-'),suffix), addText = ' class="centered"')

          for limits in consInfo[violType]['limits']:
            for limit in limits:
              htmlLimitText = "%s%.1f%s" % (prefix,limit,suffix)
              htmlPage.addMainTableColumn("%s <BR>" % htmlLimitText, addText = ' class="centered"')
        
        # Coloured by severity of violation!
        self.writeViolationColumn(htmlPage,self.output['ensemble']['constraintLists'][dclKey]['constraints'],consKey,violType,clType)
        
        for strucIndex in self.output['modelOrder']:
          if strucIndex in self.excludeStructures:
            continue
          self.writeViolationColumn(htmlPage,self.output['models'][strucIndex]['constraintLists'][dclKey]['constraints'],consKey,violType,clType)
        
        htmlPage.closeMainTableRow()

    htmlPage.closeMainTable()
    htmlPage.finishHtml()    

  def writePlot(self,fileName,xValues,yValues,main,xlab,ylab, lab = None):

    from rpy import r
    
    if not lab:
      lab = r.c(len(xValues),len(yValues),2)

    r.bitmap(fileName, res=200)
    
    ylim = (0,r.max(yValues))
    xlim = (r.min(xValues),r.max(xValues))

    # TODO: use matplotlib instead of this? Should be higher quality...
    r.plot(xValues,yValues, type = "h", main = main, xlab = xlab, ylab=ylab, xlim = xlim, ylim = ylim, lwd = 4, pch = 6, lab = lab)

    # Turn off this print - otherwise end up with a load of different threads.
    r.dev_off()
  
  def getLocalGraphDir(self,htmlDir):

    localGraphDir = os.path.join(htmlDir,self.graphDir)
    if not os.path.exists(localGraphDir):
      os.mkdir(localGraphDir)

    return localGraphDir

  def getHtmlBaseName(self,htmlDir):

    if self.mainSaveDir == htmlDir:
      htmlDir = self.baseName
    else:
      htmlDir = htmlDir.replace("%s/" % self.mainSaveDir,'')
  
    return htmlDir

  def initialiseGraphInfo(self,violType,ignoreInfoTags,graphBaseName,localGraphDir):
      
    graphInfo = {'data': {}, 'fileNames': {}, 'links': {}}

    for infoTag in self.infoTags:
      if infoTag not in ignoreInfoTags:
        graphInfo['data'][infoTag] = []

        graphName = "%s_%s_%s.png" % (graphBaseName,violType,infoTag)
        graphInfo['links'][infoTag] =  os.path.join(self.graphDir,graphName)
        graphInfo['fileNames'][infoTag] = os.path.join(localGraphDir,graphName)
    
    return graphInfo      
        
  def writeViolationColumn(self,htmlPage,consDict,consKey,violType,clType):
  
    if consDict.has_key(consKey) and consDict[consKey].has_key(violType):
      violation = consDict[consKey][violType]
      htmlFormat = "%.2f"
      if violation > self.largeViolations[clType]:
        htmlFormat = '<font color="#FF0000">%.2f</font>'
      elif violation > self.smallViolations[clType]:
        htmlFormat = '<font color="#880000">%.2f</font>'
      
      htmlPage.addMainTableColumn(htmlFormat % violation, addText = ' class="centered"')

    else:
      htmlPage.addMainTableColumn("-", addText = ' class="centered"')

if __name__ == '__main__':
  
  import sys
  
  if len(sys.argv) > 1:
    projDir = sys.argv[1]
  else:
    projDir ="/Users/wim/workspace/stable/all/data/coco/1qnd/coco/ccpn"
      
  ccpnProject = loadProject(projDir)
  nmrProject = ccpnProject.currentNmrProject
  if not nmrProject:
    nmrProject = ccpnProject.findFirstNmrProject()
  strucGen = nmrProject.findFirstStructureGeneration()
  
  violationStatistics = ViolationStatistics(nmrProject, strucGen = strucGen) # Can also pass in structures to exclude (by index)
  violationStatistics.writeHtml("local/test",testMode = False)
