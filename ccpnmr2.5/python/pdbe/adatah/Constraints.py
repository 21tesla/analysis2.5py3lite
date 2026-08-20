
from ccpnmr.format.process.stereoAssignmentSwap import StereoAssignmentSwapCheck

from ccpnmr.format.general.Util import getResNameText, getAtomNameText

from ccp.general.Copy import createNewNmrConstraintStoreFromCopy

# For recalibration stuff - maybe should go in there
import sys, math, os

from memops.api import Implementation
from ccp.api.nmr import NmrConstraint

from pdbe.software.Util import ResonanceCoordinateHandler

from ccp.general.Util import getPseudoCorrectionsWuthrich, getPseudoCorrections
from ccp.general.Util import getAtomPositionType, getResonancesFromPairwiseConstraintItem
from ccp.general.Geometry import getDistanceFromCoordinates

class ConstraintsHandler(ResonanceCoordinateHandler):
  
  """
  For use with dataHandler classes, can also be used as standalone if self. info set correctly
  """
  
  def readDistanceConstraintFiles(self,formatName,distanceConstraintFiles,nmrConstraintStore):

    keywds = {}

    if self.presets.has_key('readDistanceConstraints'):

      scriptPresets = self.presets['readDistanceConstraints']

      if scriptPresets.has_key('keywds'):

        keywds = scriptPresets['keywds']
        

    for distanceConstraintFile in distanceConstraintFiles:
      (pathName,baseName) = os.path.split(distanceConstraintFile)
      self.formatObjectDict[formatName].readDistanceConstraints(distanceConstraintFile, constraintListName = baseName, nmrConstraintStore = nmrConstraintStore, minimalPrompts = 1, **keywds)
  
      print "  Read %s file in %s format..." % (distanceConstraintFile,formatName)
  
                
  def createSwapCheckedNmrConstraints(self,origNmrConstraintStore,
                                           structureEnsemble,
                                           newStrucAnalName = "Swap checked constraints",
                                           newStrucGenName = None,
                                           numSwapCheckRuns = 2,
                                           saveChanges = False,
                                           deassignAll=False):

    """
    Input:
    
    origNmrConstraintStore is original constraint store - this will be used as base for the new one
    structureEnsemble is the structure ensemble with the coordinates used for swap checking.
    
    self.ccpnProject, self.nmrProject, self.entry have to be set.
    
    WARNING: if newStrucGenName given, will reset self.strucGen, and add self.strucGen to self.entry.structureGenerations!
    
    """

    newNmrConstraintStore = createNewNmrConstraintStoreFromCopy(self.ccpnProject,
                                                                self.nmrProject,
                                                                origNmrConstraintStore,
                                                                structureEnsemble,
                                                                newStrucGenName = newStrucGenName,
                                                                newStrucAnalName = newStrucAnalName)
    
    if newStrucGenName:
      self.strucGen = newNmrConstraintStore.findFirstStructureGeneration()
      returnObject = self.strucGen
    else:
      self.strucAnal = newNmrConstraintStore.findFirstStructureAnalysis()
      returnObject = self.strucAnal

    #
    # Connect to the 'reprotonated' coordinate files from CNS
    #
    # TODO: eventually this will be the cleaned up files from the E-MSD!
    #

    #
    # Reset the connections for aromatic protons - this is hardcoded.
    #

    self.resetAromaticAssignments(newNmrConstraintStore,structureEnsemble)

    #
    # Run the 'swapcheck' code
    #

    self.swapCheck(newNmrConstraintStore,structureEnsemble,numSwapCheckRuns,deassignAll=deassignAll)
    
    #
    # Add to existing entry
    #
    
    if newStrucGenName:
      self.entry.addStructureGeneration(self.strucGen)
    else:
      self.entry.addStructureAnalysis(self.strucAnal)
    
    #
    # Also save to XML file if necessary
    #
    
    if saveChanges:
      newNmrConstraintStore.save()
      
    #
    # Give back the structure generation/analysis
    #
    
    return returnObject

  def resetAromaticAssignments(self,nmrConstraintStore,structureEnsemble):

    """
    Input:
    
    nmrConstraintStore
    structureEnsemble
    
    Will reset badly set aromatic assigments
    
    """

    print "\n### Resetting aromatic assignments ###\n"
    
    if not structureEnsemble:
      print "  Error no structureEnsemble available - aborting"
      return
      
    models = structureEnsemble.sortedModels()

    if not models:
      print "  Error no structureEnsemble models available for structure ensemble - aborting"
      return
    
    refMolStructure = models[0]
    molSystem = structureEnsemble.molSystem

    for chain in molSystem.sortedChains():
      for aromResName in ['Phe','Tyr']:
        aromResidues = chain.findAllResidues(ccpCode = aromResName)

        for aromRes in aromResidues:
          for aromProtonGroup in [('HD1','HD2'),('HE1','HE2')]:
            for i in range(0,2):
              aromProtonName = aromProtonGroup[i]
              aromProton = aromRes.findFirstAtom(name = aromProtonName)

              atomSet = aromProton.findFirstFixedAtomSet(nmrConstraintStore = nmrConstraintStore)

              #
              # Assuming that only one resonanceSet for the atomSet!
              #

              if atomSet and atomSet.resonanceSets:

                otherAromProtonName = aromProtonGroup[not i]
                otherAromProton = aromRes.findFirstAtom(name = otherAromProtonName)

                otherAtomSet = otherAromProton.findFirstFixedAtomSet(nmrConstraintStore = nmrConstraintStore)

                if otherAtomSet:
                  if otherAtomSet != atomSet:

                    if otherAtomSet.resonanceSets and otherAtomSet.resonanceSets != atomSet.resonanceSets:
                      resonanceSet = atomSet.sortedResonanceSets()[0]
                      otherResonances = otherAtomSet.sortedResonanceSets()[0].sortedResonances()
                      otherAtomSet.sortedResonanceSets()[0].delete()
                      for otherResonance in otherResonances:
                        resonanceSet.addResonance(otherResonance)
                        print "  Warning: moved aromatic resonance %s to joined resonanceSet..." % getResNameText(otherResonance)


                    atomSet.addAtom(otherAromProton)
                    otherAtomSet.delete()
                    print "  Warning: made joined atomSet for atoms '%s' and '%s'..." % (getAtomNameText(aromProton),getAtomNameText(otherAromProton))
                    break

                  else:
                    break

                else:
                  atomSet.addAtom(otherAromProton)
                  print "  Warning: added atom '%s' to joined atomSet..." % (getAtomNameText(otherAromProton))


  def swapCheck(self,nmrConstraintStore,structureEnsemble,numSwapCheckRuns,deassignAll=False):

    """
    Input:
    
    nmrConstraintStore
    structureEnsemble
    
    numSwapCheckRuns: number of times this swap check is performed. 2 should be enough
   
    """

    print "\n### Checking stereo swaps and deassignment ###"

    swapCheck = StereoAssignmentSwapCheck(nmrConstraintStore,structureEnsemble,verbose = True)

    violationCodes = {'xl': {'violation': 1.0, 'fraction': 0.00001},
                      'l': {'violation': 0.5, 'fraction': 0.5}}

    for swapCheckRun in range(numSwapCheckRuns):
      swapCheck.checkSwapsAndClean(violationCodes=violationCodes, deassignAll=deassignAll)

    print
    print

  # TODO: is reasonably independent, could be ripped out and put somewhere else... is that useful though?
  def recalibrateConstraints(self,
                             nmrConstraintStore,
                             structures,
                             pseudo = 'Wuthrich',
                             saveResults = True,
                             ensembleAverage = 'median', # Can be 'median', 'average', or 'NOE'. Last one is r-6 average
                             classify = [],
                             fout = sys.stdout,
                             correlationMethod = 'spearman'):
  
    """
    
    Recalibrates a set of distance constraint lists based on a set of coordinates.
    
    Input:
    
     args:
      nmrConstraintStore: A CCPN Nmr.NmrConstraintStore (NOTE THIS WILL BE MODIFIED IF SAVERESULTS IS TRUE!)
      structures:         A list of CCPN MolSystem.Models
     
     keywds:
      pseudo:       'Wuthrich' to use Wuthrich pseudo atom corrections
                    'Generic' to use generic pseudo atom corrections (based on CONCOORD ones)
      saveResults:  Boolean. Write out the results (or not)
      classify:     Set to ['atomType'] for using bb,sc atom classifications in calibration.
      fout:         File handle to write result information to.
  
    Output:
    
      None (except for writing outFile and possibly new data to CCPN project)
        
    """

    from pdbe.analysis.Stats import getCorrelation
      
    #
    # Initialise some stuff
    #
    
    invSixth = - 1.0 / 6.0    
    fixedResonances = nmrConstraintStore.sortedFixedResonances()
    
    #
    # Check if pseudo correction info OK
    #
    
    if pseudo == 'Wuthrich':
      pseudoCorrections = getPseudoCorrectionsWuthrich(fixedResonances)
    elif pseudo == 'Generic':
      pseudoCorrections = getPseudoCorrections(fixedResonances)
    else:
      print "  ERROR: Unrecognized pseudo correction system '%s' - aborting."
      return
    
    #
    # Write header
    #
    
    if classify:
      classifyString = ", classification by '%s'" % str(classify)
    else:
      classifyString = ""
    
    headerLine = " # Recalibrating constraints, based on '%s', ensemble averaging '%s'%s. #\n" % (pseudo,ensembleAverage,classifyString)
    headerFrameLine = " %s\n" % ("#" * (len(headerLine)-1))
    
    fout.write("\n")
    fout.write(headerFrameLine)  
    fout.write(headerLine)
    fout.write(headerFrameLine)
    fout.write("\n")
  
    #
    # Set up dict for resonance->atom links
    #
    
    self.setAssignedAtomsAndResidues(fixedResonances) # From ResonanceCoordinateHandler class!
      
    #
    # Set up reference info for structure coords
    #
    
    self.structureList = list(structures)
    
    numStructures = len(self.structureList)    
    if numStructures == 1:
      fout.write("  Warning: only one structure used in analysis!\n")
  
    self.createCoordAtomInfoDict() # From ResonanceCoordinateHandler class!
        
    #
    # Now go over the constraint lists, handle 1 by 1...
    #
    
    allDistanceConstraintLists = []
    
    for constraintList in nmrConstraintStore.sortedConstraintLists():
      if constraintList.className == 'DistanceConstraintList':
        allDistanceConstraintLists.append(constraintList)
    
    for dcl in allDistanceConstraintLists:
    
      dclConstraints = dcl.sortedConstraints()
      constraintKeys = []
      
      distances = []
      targetDistances = []
      targetCorrected = []
      upperDistances = []
      upperCorrected = []
  
      distanceClasses = []
      pseudoDistances = []
      pseudoTargetDistances = []
      pseudoTargetCorrectedDistances = []
      pseudoUpperDistances = []
      pseudoUpperCorrectedDistances = []
      
      if not classify:
        coordinateVolumeSum = {'all': 0.0}
      else:
        coordinateVolumeSum = {}
  
      for constraint in dclConstraints:
      
        distPerStruc = []
        atomCombs = []
        atomTypes = []
        
        for strucIndex in range(numStructures):
          distPerStruc.append([])
        
        hasDistance = False
        pseudoCorrection = self.getPseudoCorrection(constraint,pseudoCorrections)
        
        """           
        if dcl.serial == 9 and constraint.serial in [6,7,8,9]:
          print "  ", constraint.serial
          #print "    ",constraint.sortedItems()
          #print "    ", distPerStruc
          #print "    ", ["%s.%d.%s-%s.%d.%s" % (atomComb[0].residue.ccpCode,atomComb[0].residue.seqId,atomComb[0].name,atomComb[1].residue.ccpCode,atomComb[1].residue.seqId,atomComb[1].name)  for atomComb in atomCombs]
        """           

        for item in constraint.sortedItems():
          (resonance,otherResonance) = getResonancesFromPairwiseConstraintItem(item)
          
          if self.resObjectMapping.has_key(resonance) and self.resObjectMapping.has_key(otherResonance):
            (residue,atomList) = self.resObjectMapping[resonance]
            (otherResidue,otherAtomList) = self.resObjectMapping[otherResonance]
            
            """           
            if dcl.serial == 9 and constraint.serial in [6,7,8,9]:
              print "    ", residue, atomList
              print "    ", otherResidue,otherAtomList
            """           

            # Set the atom type
            if 'atomType' in classify:
              atomContactType = []
  
              for atomRefName in (atomList[0].name,otherAtomList[0].name):
                atomType = getAtomPositionType(atomRefName)
                atomContactType.append(atomType)
  
              atomContactType.sort()
              atomContactType = tuple(atomContactType)
              if not atomContactType in atomTypes:
                atomTypes.append(atomContactType)
  
            # Get distance info...
            for atom in atomList:
              if self.coordAtomInfo.has_key(atom):
                for otherAtom in otherAtomList:
                  if atom == otherAtom:
                    continue
                    
                  atomComb = (atom,otherAtom)
                  
                  # Only count each contribution once!!
                  if atomComb in atomCombs or (atomComb[1],atomComb[0]) in atomCombs:
                    continue
                  
                  atomCombs.append(atomComb)
                    
                  if self.coordAtomInfo.has_key(otherAtom):
                    hasDistance = True
                    for strucIndex in range(numStructures):
                      coord = self.coordAtomInfo[atom][strucIndex]
                      otherCoord = self.coordAtomInfo[otherAtom][strucIndex]
                      
                      if coord and otherCoord: 
                        distance = getDistanceFromCoordinates(coord,otherCoord)
                        distPerStruc[strucIndex].append(distance)
           
        if not hasDistance:
          continue
        
        #
        # Take average based on median, straight average or NOE average (r-6)
        #
        
        overallDistSum = 0.0
        tempDistTotal = 0.0
        numDistZero = 0
        tempDistanceList = []
  
        for strucIndex in range(numStructures):
          distSum = 0.0
          for distance in distPerStruc[strucIndex]:
            distSum += math.pow(distance,-6)
         
          if distSum:
            if ensembleAverage == 'NOE':
              # This is correct
              overallDistSum += distSum
            else:
              tempDist = math.pow(distSum,invSixth)
              tempDistanceList.append(tempDist)
              tempDistTotal += tempDist
          else:
            numDistZero += 1
              
        # TODO: here aim 'lower' than actual target distance? Aymeric can do this himself though...
        # Not implemented currently
        numDistances = numStructures - numDistZero
        
        if not numDistances:
          continue
        
        if ensembleAverage == 'NOE':
          overallDistSum = overallDistSum / numDistances
          avgDist = math.pow(overallDistSum,invSixth)
  
        else:
          numDistances = len(tempDistanceList)
          # TODO ARE THERE GOOD BUILTIN PYTHON FUNCIONS FOR THIS?
          if ensembleAverage == 'median':
            tempDistanceList.sort()
            middleIndex = int(numDistances/2)
            if numDistances % 2 == 0:
              avgDist = (tempDistanceList[middleIndex-1] + tempDistanceList[middleIndex]) / 2
            else:
              avgDist = tempDistanceList[middleIndex]
          elif ensembleAverage == 'average':
            avgDist = tempDistTotal / numDistances
          else:
            avgDist = 0.0
  
          overallDistSum = math.pow(avgDist,-6)
          
        """           
        if dcl.serial == 9 and constraint.serial in [6,7,8,9]:
          print "  ", avgDist
          print
        """           

        #
        # Now set the information
        #
          
        if not classify:
          classType = 'all'
        elif 'atomType' in classify:
          atomTypes.sort()
          atomTypes = tuple(atomTypes)
          classType = atomTypes
          # TODO Need to add other types if necessary, also combine info...
          if not coordinateVolumeSum.has_key(classType):
            coordinateVolumeSum[classType] = 0.0
  
        coordinateVolumeSum[classType] += overallDistSum
  
        distances.append(avgDist)
        constraintKeys.append(constraint.getFullKey())
  
        if pseudoCorrection:
          pseudoDistances.append(avgDist)
        
        #
        # Get the restraint distance info
        #
                
        if constraint.targetValue:
          targetValue = constraint.targetValue
        else:
          targetValue = 0.0
          
        distanceClasses.append(classType)
        targetDistances.append(targetValue)        
        targetCorrected.append(targetValue - pseudoCorrection)
        
        if constraint.upperLimit:
          upperLimit = constraint.upperLimit
        else:
          upperLimit = 0.0
          
        upperDistances.append(upperLimit)
        upperCorrected.append(upperLimit - pseudoCorrection)
  
        if pseudoCorrection:
          pseudoTargetDistances.append(targetValue)
          pseudoTargetCorrectedDistances.append(targetValue - pseudoCorrection)
          pseudoUpperDistances.append(upperLimit)
          pseudoUpperCorrectedDistances.append(upperLimit - pseudoCorrection)
          
          #print residue.ccpCode, atomList[0].name, otherResidue.ccpCode, otherAtomList[0].name, avgDist, upperLimit, pseudoCorrection
      
      
      #
      # TODO: also check number of values that are available? If only a couple (out of 1000s), then don't use that one!
      #
      
      fout.write("Constraint list %s, with %d distances.\n" % (dcl.serial, len(distances)))
      if not distances:
        fout.write("  No distances - ignored.\n")
        continue
  
      fout.write("  Correlations by %s method:\n" % correlationMethod)
      
      targetCorr = getCorrelation(distances,targetDistances, correlationMethod = correlationMethod)
      fout.write("    All target distances, as is: %.4f.\n" % targetCorr)
      upperCorr = getCorrelation(distances,upperDistances, correlationMethod = correlationMethod)    

      """           
      print distances[0], distances[-1], upperDistances[0], upperDistances[-1]
      if dcl.serial == 9:
        print
        print distances
        print
      """  
               
      fout.write("    All upper distances, as is:  %.4f.\n" % upperCorr)
      
      if targetCorr > upperCorr:
        selectText = "Selecting target distances as best option"
        useDistances = targetDistances
        distType = 'target'
      
      elif 0 < upperCorr < 1:
  
        pseudoUpperCorr = 0.0
        
        if pseudoDistances:    
          pseudoUpperCorr = getCorrelation(distances,upperCorrected, correlationMethod = correlationMethod)    
          fout.write("    All upper distances, pseudo corrected:  %.4f\n" % pseudoUpperCorr)
        
        if not pseudoUpperCorr or pseudoUpperCorr < upperCorr:
          selectText = "Selecting upper bounds as best option"
          useDistances = upperDistances
          distType = 'upper'
        else:
          selectText = "Selecting pseudo corrected upper bounds as best option"
          useDistances = upperCorrected
          distType = 'pseudo'
      
      else:
        
        # TODO THIS IS NOT GOOD - just copy over original list as is? Or what?!?!?
        # Info will now go missing - check where this happens!
      
        fout.write("  Ignoring list - no valid distances.\n\n")
        continue
  
      fout.write("%s\n\n\n" % selectText)
      
      newConstraintList = NmrConstraint.DistanceConstraintList.getByKey(nmrConstraintStore, (dcl.serial,))
      newConstraintList.details = selectText
      
      newConstraintList.addApplicationData(
      
          Implementation.AppDataString(application = 'ConstraintsHandler', keyword = 'distType', value = distType)
          
          )
      
      #
      # Now create a new list, with recalculated target distances
      #
      
      volumes = []
      volumeSum = {}
      for classType in coordinateVolumeSum.keys():
        volumeSum[classType] = 0.0
      
      for distIndex in range(len(useDistances)):
        useDistance = useDistances[distIndex]
        classType = distanceClasses[distIndex]
        if useDistance:
          volume = math.pow(useDistance,-6)
        else:
          volume = 0.0
        volumeSum[classType] += volume
        volumes.append((volume,classType))
  
      correctionFactor = {}
      for classType in volumeSum.keys():
        correctionFactor[classType] = coordinateVolumeSum[classType] / volumeSum[classType]
        
        if type(classType) == type(""):
          classTypeString = classType
        else:
          classTypeString = ','.join(["%s-%s" % classTypeItem for classTypeItem in classType])

        newConstraintList.addApplicationData(
      
          Implementation.AppDataFloat(application = 'ConstraintsHandler', keyword = 'correctionFactor_%s' % classTypeString, value = correctionFactor[classType])
          
          )

      
      for i in  range(0,len(volumes)):
        (volume,classType) = volumes[i]
        
        if volume:
          newDistance = math.pow(volume * correctionFactor[classType],invSixth)
        else:
          # This is the default value
          newDistance = 1.8
          
        # Consistency check
        if newDistance < 1.8:
          newDistance = 1.8
               
        #
        # Change the relevant constraint...
        #
        
        constraintKey = constraintKeys[i][1:]
        newConstraint = NmrConstraint.DistanceConstraint.getByKey(nmrConstraintStore, constraintKey)
  
        newConstraint.targetValue = newDistance  
        
        # Do some consistency checks...
        if not newConstraint.weight:
          newConstraint.weight = 1.0
        
        if not newConstraint.lowerLimit or not newConstraint.upperLimit:
          correction = 0.125 * (newConstraint.targetValue ** 2)
  
          if not newConstraint.lowerLimit:
            lowerLimit = newConstraint.targetValue  - correction
            if lowerLimit < 1.8:
              lowerLimit = 1.8
            newConstraint.lowerLimit = lowerLimit
          
          if not newConstraint.upperLimit:
            newConstraint.upperLimit = newConstraint.targetValue + correction
        
    if saveResults:
      # TODO will this take care of everything?
      nmrConstraintStore.root.saveModified()
        
    return nmrConstraintStore


  def applyPseudoCorrection(self,
                            nmrConstraintStore,
                            pseudo = 'Wuthrich',
                            distanceType = 'upper', # Can also be target
                            saveResults = True,
                            fout = sys.stdout):
  
    """
    
    Applies pseudoatom corrections to upperbounds or target distances.
    
    Input:
    
     args:
      nmrConstraintStore: A CCPN Nmr.NmrConstraintStore (NOTE THIS WILL BE MODIFIED IF SAVERESULTS IS TRUE!)
     
     keywds:
      pseudo:       'Wuthrich' to use Wuthrich pseudo atom corrections
                    'Generic' to use generic pseudo atom corrections (based on CONCOORD ones)
      distanceType: 'upper' to apply corrections to upper limits
                    'target' to apply corrections to target distance values
      saveResults:  Boolean. Write out the results (or not)
      fout:         File handle to write result information to.
  
    Output:
    
      None
        
    """
        
    #
    # Check if pseudo correction info OK
    #
    
    fixedResonances = nmrConstraintStore.sortedFixedResonances()
    
    if pseudo == 'Wuthrich':
      pseudoCorrections = getPseudoCorrectionsWuthrich(fixedResonances)
    elif pseudo == 'Generic':
      pseudoCorrections = getPseudoCorrections(fixedResonances)
    else:
      print "  ERROR: Unrecognized pseudo correction system '%s' - aborting."
      return
        
    self.setAssignedAtomsAndResidues(fixedResonances) # From ResonanceCoordinateHandler class!

    #
    # Now go over the constraint lists, handle 1 by 1...
    #
    
    allDistanceConstraintLists = nmrConstraintStore.findAllConstraintLists(className = 'DistanceConstraintList')
    
    for dcl in allDistanceConstraintLists:
      
      numberCorrections = 0
    
      dclConstraints = dcl.sortedConstraints()
  
      for constraint in dclConstraints:
              
        pseudoCorrection = self.getPseudoCorrection(constraint,pseudoCorrections)

        if pseudoCorrection:

          if distanceType == 'upper' and constraint.upperLimit:
            constraint.upperLimit = constraint.upperLimit - pseudoCorrection
            numberCorrections += 1
          elif distanceType == 'target' and constraint.targetValue:
            constraint.targetValue = constraint.targetValue - pseudoCorrection
            numberCorrections += 1
              
      print "Distance constraint list %d: corrected %d constraints (out of %d)" % (dcl.serial,numberCorrections,len(dclConstraints))
             
    if saveResults:
      nmrConstraintStore.root.saveModified()

  def getPseudoCorrection(self,constraint,pseudoCorrections):

    constraintItems = constraint.sortedItems()

    #
    # Can only really do pseudocorrection on unambiguous item,
    # or on multiple items that designate an atom set...
    #

    resonanceList = []

    for consItem in constraintItems:
      resonanceList.append(getResonancesFromPairwiseConstraintItem(consItem))

    atomsForPseudoCorrection = []

    residueList = []

    #print constraint.serial

    for i in range(2):

      chemAtomSets = []
      deepChemAtomSets = []
      atomNames = []

      useAtomSetName = None

      if len(resonanceList) > 1:
        checkForSets = True
      else:
        checkForSets = False

      for resonances in resonanceList: 
        resonance = resonances[i]
        
        if self.resObjectMapping.has_key(resonance):
          (residue,atomList) = self.resObjectMapping[resonance]
        else:
          useAtomSetName = checkForSets = False
          atomNames = []
          break

        if len(atomList) > 1:
          checkForSets = True

        residueList.append(residue)

        atomNames.append(atomList[0].name)

        chemAtomSet = deepChemAtomSet = None

        chemAtom = atomList[0].chemAtom
        if chemAtom.chemAtomSet:
          chemAtomSet = chemAtom.chemAtomSet.name

          if chemAtom.chemAtomSet.chemAtomSet:
            deepChemAtomSet = chemAtom.chemAtomSet.chemAtomSet.name

        chemAtomSets.append(chemAtomSet)
        deepChemAtomSets.append(deepChemAtomSet)

      # First check if got multiple 'deep' chemAtomSets
      if not deepChemAtomSets.count(None) and checkForSets:
        if deepChemAtomSets.count(deepChemAtomSets[0]) == len(deepChemAtomSets):
          useAtomSetName = deepChemAtomSets[0]
      if not chemAtomSets.count(None) and checkForSets:
        if chemAtomSets.count(chemAtomSets[0]) == len(chemAtomSets):
          # Note that this will reset HG* to HG1*, for example, if this is the only
          # one that occurs!
          useAtomSetName = chemAtomSets[0]
      if (not checkForSets or not useAtomSetName) and atomNames and atomNames.count(atomNames[0]) == len(atomNames):
        useAtomSetName = atomNames[0]

      #print "    %d" % i, chemAtomSets, deepChemAtomSets, atomNames

      atomsForPseudoCorrection.append(useAtomSetName)

    #print "  ",atomsForPseudoCorrection

    if residueList and len(residueList) == residueList.count(residueList[0]):
      intraResidue = True
      #print "  intrares, %s" % residueList[0].ccpCode
    else:
      intraResidue = False
      #print "  interres, %s-%s" % (residueList[0].ccpCode,residueList[-1].ccpCode)

    #
    # Now check if have to apply correction...
    #

    pseudoCorrection = 0.0

    if None not in atomsForPseudoCorrection and len(resonanceList) == 2:

      for i in range(2):
        resonance = resonanceList[0][i]
        if pseudoCorrections.has_key(resonance):

          if type(pseudoCorrections[resonance]) == type(0.0):
            pseudoCorrection += pseudoCorrections[resonance]

          # Now handle Wuthich exceptions...
          elif pseudoCorrections[resonance].has_key(atomsForPseudoCorrection[i]):
            corrInfo = pseudoCorrections[resonance][atomsForPseudoCorrection[i]]

            if intraResidue and corrInfo[1].has_key(atomsForPseudoCorrection[not i]):
              pseudoCorrection += corrInfo[1][atomsForPseudoCorrection[not i]]
            else:
              pseudoCorrection += corrInfo[0]

    return pseudoCorrection
