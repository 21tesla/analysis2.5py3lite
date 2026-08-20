"""cyana2ccpnConverter.py

Uses CyanaParser from Cing to import cyana data

Uses Project from which the calculation originated or creates a new one if desired

Initialises experiment in a new project or clones existing experiments

Creates a new peak list for each new/cloned experiment by picking peaks defined in the CYANA output peak files
and manually sets the intensities

Stores original assignment state as a Chemical Shift Restraint List and creates new shiftlist to reflect CYANA assignment
state.

Assigns analysis peak lists using a cing-to-CCPN resonance dictionary generated using CYANA assignment state

Imports restraints into the CCPN project and assigns them using CYANA assignment state

Imports structure(s) into the CCPN project using a CCPN routine

"""

import CyanaParser.CyanaParser as cp
from ccpnmr.integrator.core import Util as intUtil
from ccpnmr.analysis.core.PeakBasic import pickPeak, setManualPeakIntensity, copyPeakList
from ccpnmr.analysis.core.AssignmentBasic import assignResToDim, assignAtomsToRes
from ccpnmr.analysis.core.ConstraintBasic import makeStructureGeneration, getFixedResonance
from ccp.lib.StructureIo import getStructureFromFile
from ccpnmr.analysis.core.ExperimentBasic import cloneExperiment, setExperimentShiftList
from ccpnmr.format.converters.DataFormat import getResonanceAtomMap
import os
import json



def importFromCyana(nmrCalcRun, targetDir):

  c = cp.CyanaParser()

  rootProject = nmrCalcRun.root
  os.chdir(targetDir)

  peakFiles = []
  uplFiles = []
  acoFiles = []
  rdcFiles = []
  cyaFiles = []
  noaFiles = []
  txtFiles = []
  dataSources = []
  ccpnConfigFile = open("Properties.ccpn.json")
  ccpnConfig = json.load(ccpnConfigFile)
  configFile = open("cyana2ccpn.json")
  config = json.load(configFile)
  runId = ccpnConfig['CCPN.nmrCalcId'].split('+')[1]


  for peakList in ccpnConfig['PeakListData']:
    peakListName = peakList['fileName'].split('.peaks')[0]+'-cycle7.peaks'
    peakFiles.append(peakListName)
  for distanceRestraintList in config["distanceRestraints"]:
    uplFiles.append(distanceRestraintList)
  for dihedralRestraintList in config["dihedralRestraints"]:
    acoFiles.append(dihedralRestraintList)
  for rdcRestraintList in config["rdcRestraints"]:
    rdcFiles.append(rdcRestraintList)
  for cyaFile in config["cyaFiles"]:
    cyaFiles.append(cyaFile)
  for noaFile in config["noaFiles"]:
    noaFiles.append(noaFile)
  seqFile = ccpnConfig['RunParameter']['fileNameSequence']
  originalProtFile = ccpnConfig['MeasurementListData'][0]['fileName']

  finalProtFile = originalProtFile.split(".prot")[0] + "-final.prot"

  c.parse(seqFile=seqFile, originalProtFile=originalProtFile,
          finalProtFile=finalProtFile, peakFiles=peakFiles,
         rdcFiles=rdcFiles,
          acoFiles=acoFiles, cyaFiles=cyaFiles,noaFiles=noaFiles)
  nmrProject = rootProject.currentNmrProject
  nmrConstraintStore = rootProject.newNmrConstraintStore(nmrProject=nmrProject)
  if rootProject.currentAnalysisProject is not None:
    AnalysisProject = rootProject.currentAnalysisProject
  else:
    AnalysisProject=rootProject.newAnalysisProject(name="analysisProject",nmrProject=nmrProject)
  shiftListName = c.finalProtFile
  molSystem = rootProject.findFirstMolSystem()
  newShiftList = nmrProject.newShiftList(name=shiftListName)


  atomToResonanceMap = createAtomtoResonanceMap(nmrProject)
  resonanceDictionaries = assignResonances(c.resonances, atomToResonanceMap, molSystem, newShiftList,
                                           AnalysisProject, nmrConstraintStore)
  loadPdb(config['pdbFile'], molSystem)

  # newPeakLists = []
  for peakList in ccpnConfig['PeakListData']:
    peakListName = peakList['fileName'].split('.peaks')[0]
    experimentSerial = int(peakList['fileName'].split('_')[1])
    ccpnExperiment = nmrProject.findFirstExperiment(serial=experimentSerial)
    cloneName = ccpnExperiment.name + '_run' + str(runId)
    clonedExperiment = cloneExperiment(ccpnExperiment, cloneName)
    spectrum = clonedExperiment.findFirstDataSource()
    dataSources.append(spectrum)
    setExperimentShiftList(clonedExperiment, newShiftList)
    peakList = spectrum.newPeakList()
    spectrum.activePeakList = peakList

    for cingPeakList in c.peakLists:

      if cingPeakList.name.split('-')[0] == peakListName:
        peakList.name = cingPeakList.name
        pickPeaksFromCing(cingPeakList, peakList, resonanceDictionaries['resonanceDict'])
        nmrCalcRun.newPeakListData(name=cingPeakList.name, ioRole='output', peakList=peakList)
        # newPeakLists.append(peakList)


    # for newPeakList in newPeakLists:
    assignedPeakList = spectrum.newPeakList()
    unassignedPeakList = spectrum.newPeakList()
    assignedPeakList.details = 'assigned'
    unassignedPeakList.details = 'unassigned'
    splitAssignedUnassigned(peakList, assignedPeakList, unassignedPeakList)

  loadChemShiftRestraints(c.chemicalShiftRestraints, originalProtFile, nmrConstraintStore, atomToResonanceMap,
                          resonanceDictionaries['fixedResonanceDict'], molSystem)
  if len(c.distanceRestraintLists) != 0:
    loadDistanceRestraints(c.distanceRestraintLists, nmrProject, nmrCalcRun, nmrConstraintStore,
                           resonanceDictionaries['fixedResonanceDict'],
                           resonanceDictionaries['cingFixedResonanceDict'],molSystem, AnalysisProject)

  if len(c.violationLists) != 0:
    violatedPeaks = loadViolatedDistanceRestraints(c.violationLists, nmrProject, nmrCalcRun, nmrConstraintStore,
                           resonanceDictionaries['fixedResonanceDict'],
                           resonanceDictionaries['cingFixedResonanceDict'],molSystem, AnalysisProject)
    for peakList,peaks in violatedPeaks.iteritems():
      spectrum = peakList.getDataSource()
      newPeakList = spectrum.newPeakList()
      newPeakList.details = 'violated'
      copyPeakList(peakList, newPeakList, peaks=peaks)


  if len(c.dihedralRestraintLists) != 0:
    loadDihedralRestraints(c.dihedralRestraintLists,nmrConstraintStore, nmrCalcRun,
                           resonanceDictionaries['cingFixedResonanceDict'], molSystem)

  if len(c.rdcRestraintLists) != 0:
    loadRdcRestraints(c.rdcRestraintLists,nmrConstraintStore,
                      resonanceDictionaries['cingFixedResonanceDict'], molSystem)

  pluginModule = intUtil.getIntegratorPlugin(ccpnConfig["CCPN.Run.wmsProtocolName"])
  pluginModule.read.read(nmrCalcRun, targetDir)
  print "Import Complete"
  return dataSources




def splitAssignedUnassigned(peakList, assignedPeakList, unassignedPeakList):

  assignedPeaks = []
  unassignedPeaks = []
  for peak in peakList.peaks:
    contribs =  peak.sortedPeakContribs()
    if len(contribs) > 0:
      assignedPeaks.append(peak)
    else:
      unassignedPeaks.append(peak)

  copyPeakList(peakList, assignedPeakList, peaks=assignedPeaks)
  copyPeakList(peakList, unassignedPeakList, peaks=unassignedPeaks)


def createAtomtoResonanceMap(nmrProject):
  atomToResonanceMap = {}
  resonanceToAtomMap = getResonanceAtomMap('IUPAC', nmrProject.sortedResonances())
  for k, v in resonanceToAtomMap.iteritems():
    atomTuple = (v[0].chain, v[0].seqId, v[0].atomName)
    atomToResonanceMap[atomTuple] = k
  return atomToResonanceMap


def assignResonances(resonances, atomToResonanceMap,molSystem, ccpnShiftList, AnalysisProject, nmrConstraintStore):

  print "assigning resonances"
  atomTupleDict = {}
  resonanceDict = {}
  fixedResonanceDict = {}
  cingFixedResonanceDict = {}

  for cingResonance in resonances:

    ccpnResidueArray = cingResonance.atom.residue.typeIdentifier.translate('CCPN').split()
    ccpnResidueName = ccpnResidueArray[1]
    if cingResonance.stereo == True:
      ### find the ccpnResonance from which the cingResonance originated and assign it.


      ringAtoms = ['CD1', 'CD2', 'CE1', 'CE2', 'HD1', 'HD2', 'HE1', 'HE2']
      ccpnRealAtomName = cingResonance.atom.typeIdentifier.translate('CCPN')
      if ccpnResidueName == 'Phe' and ccpnRealAtomName in ringAtoms:
        ccpnRealAtomName = ccpnRealAtomName[:2] + "*"
      if ccpnResidueName == 'Tyr' and ccpnRealAtomName in ringAtoms:
        ccpnRealAtomName = ccpnRealAtomName[:2] + "*"
        ccpnMappingAtomName = ccpnRealAtomName[0] + ccpnRealAtomName[1:].lower()
      if len(ccpnRealAtomName) == 1:
        ccpnMappingAtomName = ccpnRealAtomName
      else:
        ccpnMappingAtomName = ccpnRealAtomName[0] + ccpnRealAtomName[1:].lower()

      ccpnChainCode = molSystem.findFirstChain(code=cingResonance.atom.chainId)
      ccpnResNum = cingResonance.atom.sequenceId
      ccpnResidue = molSystem.findFirstChain().findFirstResidue(seqCode=ccpnResNum)
      ccpnResidueMapping = AnalysisProject.findFirstChainMapping(
        molSystemCode=molSystem.code,chainCode=cingResonance.atom.chainId).findFirstResidueMapping(
        seqId=ccpnResidue.seqId)
      ccpnAtomMapping = ccpnResidueMapping.findFirstAtomSetMapping(name=ccpnMappingAtomName)

      atomTuple = (ccpnChainCode, ccpnResidue.seqId, ccpnRealAtomName)
      if atomTuple not in atomTupleDict:
        atomTupleDict[atomTuple] = cingResonance
      resonance = atomToResonanceMap[atomTuple]
      assignAtomsToRes(ccpnAtomMapping.atomSets, resonance)

      # except KeyError:
      #   pass



    else:
      ### find the ccpnResonance from which the cingResonance originated and assign it.
      ringAtoms = ['CD1', 'CD2', 'CE1', 'CE2', 'HD1', 'HD2', 'HE1', 'HE2']
      ccpnRealAtomName = cingResonance.atom.typeIdentifier.translate('CCPN')
      if ccpnResidueName == 'Phe' and ccpnRealAtomName in ringAtoms:
        ccpnRealAtomName = ccpnRealAtomName[:2] + "*"

      if ccpnResidueName == 'Tyr' and ccpnRealAtomName in ringAtoms:
        ccpnRealAtomName = ccpnRealAtomName[:2] + "*"

      ccpnChainCode = molSystem.findFirstChain(code=cingResonance.atom.chainId)
      ccpnResNum = cingResonance.atom.sequenceId
      ccpnResidue = molSystem.findFirstChain().findFirstResidue(seqCode=ccpnResNum)
      atomTuple = (ccpnChainCode, ccpnResidue.seqId, ccpnRealAtomName)
      if atomTuple not in atomTupleDict:
        atomTupleDict[atomTuple] = cingResonance

  for atomTuple in atomTupleDict:
    cingResonance = atomTupleDict[atomTuple]
    try:
      resonance = atomToResonanceMap[atomTuple]
      shift=ccpnShiftList.newShift(value=cingResonance.value, resonance=resonance)
      if resonance is not None:
        resonanceDict[cingResonance] = resonance
        fixedResonance = getFixedResonance(nmrConstraintStore, resonance)
        fixedResonanceDict[resonance] = fixedResonance
        cingFixedResonanceDict[cingResonance.atom] = fixedResonance

    except KeyError:
      pass
  return {"resonanceDict":resonanceDict, "fixedResonanceDict":fixedResonanceDict,
            "cingFixedResonanceDict":cingFixedResonanceDict}

  print "Resonances assigned"


def pickPeaksFromCing(cingPeakList, peakList, resonanceDict):


  for cingPeak in cingPeakList:
    peakList.root.override = True
    peak = pickPeak(peakList, cingPeak.positions, unit='ppm', doFit=False,
                    serial=cingPeak.xeasyIndex)
    peakList.root.override = False
    setManualPeakIntensity(peak, value=cingPeak.height.value, intensityType="height")
    if len(peak.peakDims) == len(cingPeak.resonances):
      for j, peakDim in enumerate(peak.sortedPeakDims()):
        cingResonance = cingPeak.resonances[j]

        resonance = resonanceDict.get(cingResonance)
        if resonance is not None:
          assignResToDim(peakDim, resonance)

    else:
      # original, fixed
      for j, peakDim in enumerate(peak.sortedPeakDims()):
        numAssignments = len(cingPeak.resonances) / len(peak.peakDims)
        ccpnResonances = []
        indx = j
        for k in range(numAssignments):
          cingResonance = cingPeak.resonances[indx]
          resonance = resonanceDict.get(cingResonance)
          ccpnResonances.append(resonance)
          indx += len(peak.peakDims)
        for res in ccpnResonances:
          if res is not None:
            assignResToDim(peakDim, res)

  print "Peaks Picked"

def loadPdb(filename, molSystem):
  structure = getStructureFromFile(molSystem, filename, doWarnings=False)
  makeStructureGeneration(structure)
  print "Structures loaded"

def loadDistanceRestraints(cingRestraintList, nmrProject, nmrCalcRun, nmrConstraintStore,
                           fixedResonanceDict, cingFixedResonanceDict,molSystem, AnalysisProject):
  for cingList in cingRestraintList:
    name = cingList.name
    distanceRestraintList = nmrConstraintStore.newDistanceConstraintList(name=name)
    for cingRestraint in cingList:

      restraint = distanceRestraintList.newDistanceConstraint(upperLimit=float(cingRestraint.upper), lowerLimit=cingRestraint.lower)
      for atoms in cingRestraint.atomPairs:
        ccpnFixedResonances = [[], []]
        for i, atom in enumerate(atoms):
          if "Q" in atom.atomId:
            ccpnPseudoAtom = atom.typeIdentifier.translate('CCPN')
            ccpnPseudoAtomMappingName = ccpnPseudoAtom[0] + ccpnPseudoAtom[1:].lower()
            ccpnResNum = atom.sequenceId
            ccpnResidue = molSystem.findFirstChain().findFirstResidue(seqCode=ccpnResNum)
            ccpnChainCode = atom.chainId
            ccpnResidueMapping = AnalysisProject.findFirstChainMapping(molSystemCode=molSystem.code,
                                                                            chainCode=ccpnChainCode
                                                                      ).findFirstResidueMapping(seqId=ccpnResidue.seqId)
            ccpnPseudoAtomMapping = ccpnResidueMapping.findFirstAtomSetMapping(name=ccpnPseudoAtomMappingName)
            resonances = []
            for atomSet in ccpnPseudoAtomMapping.atomSets:
              for resonanceSet in atomSet.resonanceSets:
                resonances.extend(resonanceSet.resonances)
            resonances = list(set(resonances))
            for res in resonances:
              ccpnFixedResonances[i].append(fixedResonanceDict[res])


          else:
            try:
              ccpnFixedResonances[i].append(cingFixedResonanceDict[atom])
            except KeyError:
              pass


        for fr1 in ccpnFixedResonances[0]:
          for fr2 in ccpnFixedResonances[1]:
            try:
              restraint.newDistanceConstraintItem(resonances=[fr1, fr2])
            except:
              pass
      ccpnPeaks=set()
      for ccpnExperiment in nmrProject.experiments:
        try:
          spectrum = ccpnExperiment.findFirstDataSource()
          peak = spectrum.findFirstPeakList(name=cingRestraint.peakList).findFirstPeak(serial=int(cingRestraint.peak))
          ccpnPeaks.add(peak)
        except (AttributeError, TypeError):
          pass
      if None not in ccpnPeaks:
        restraint.setPeaks(ccpnPeaks)
    nmrCalcRun.newConstraintStoreData(name=name, ioRole='output',constraintStoreSerial=nmrConstraintStore.serial)


  print "Distance restraints assigned"


def loadViolatedDistanceRestraints(cingRestraintList, nmrProject, nmrCalcRun, nmrConstraintStore,
                                   fixedResonanceDict, cingFixedResonanceDict,molSystem, AnalysisProject):
  violatedPeaks = {}
  for cingList in cingRestraintList:
    name = cingList.name+"_violated"

    distanceRestraintList = nmrConstraintStore.newDistanceConstraintList(name=name)
    for cingRestraint in cingList:

      restraint = distanceRestraintList.newDistanceConstraint(upperLimit=float(cingRestraint.upper),lowerLimit=cingRestraint.lower)
      for atoms in cingRestraint.atomPairs:
        ccpnFixedResonances = [[], []]
        for i, atom in enumerate(atoms):
          if "Q" in atom.atomId:
            ccpnPseudoAtom = atom.typeIdentifier.translate('CCPN')
            ccpnPseudoAtomMappingName = ccpnPseudoAtom[0] + ccpnPseudoAtom[1:].lower()
            ccpnResNum = atom.sequenceId
            ccpnResidue = molSystem.findFirstChain().findFirstResidue(seqCode=ccpnResNum)
            ccpnChainCode = atom.chainId
            ccpnResidueMapping = AnalysisProject.findFirstChainMapping(molSystemCode=molSystem.code,
                                                                            chainCode=ccpnChainCode
                                                                      ).findFirstResidueMapping(seqId=ccpnResidue.seqId)
            ccpnPseudoAtomMapping = ccpnResidueMapping.findFirstAtomSetMapping(name=ccpnPseudoAtomMappingName)
            resonances = []
            for atomSet in ccpnPseudoAtomMapping.atomSets:
              for resonanceSet in atomSet.resonanceSets:
                resonances.extend(resonanceSet.resonances)
            resonances = list(set(resonances))
            for res in resonances:
              ccpnFixedResonances[i].append(fixedResonanceDict[res])


          else:
            try:
              ccpnFixedResonances[i].append(cingFixedResonanceDict[atom])
            except KeyError:
              pass

        for fr1 in ccpnFixedResonances[0]:
          for fr2 in ccpnFixedResonances[1]:
            restraint.newDistanceConstraintItem(resonances=[fr1, fr2])
      ccpnPeaks=set()
      for ccpnExperiment in nmrProject.experiments:
        try:
          spectrum = ccpnExperiment.findFirstDataSource()
          peakList = spectrum.findFirstPeakList(name=cingRestraint.peakList)
          peak = peakList.findFirstPeak(serial=int(cingRestraint.peak))
          ccpnPeaks.add(peak)
          if peak is not None:
            ll = violatedPeaks.get(peakList)
            if ll is None:
              violatedPeaks[peakList] = [peak]
            else:
              if peak not in ll:
                ll.append(peak)
        except (AttributeError, TypeError):
          pass
      if None not in ccpnPeaks:
        restraint.setPeaks(ccpnPeaks)

    nmrCalcRun.newConstraintStoreData(name=name, ioRole='output',constraintStoreSerial=nmrConstraintStore.serial)

  print "Violations assigned"

  return violatedPeaks



def loadDihedralRestraints(cingRestraintList, nmrConstraintStore, nmrCalcRun, cingFixedResonanceDict, molSystem):
  for cingList in cingRestraintList:
    name = cingList.name
    dihedralRestraintList = nmrConstraintStore.newDihedralConstraintList(name=name)
    for cingRestraint in cingList:
      fixedResonances = []
      #~ ccpnAtoms = []
      for i, atom in enumerate(cingRestraint.atoms):
        if atom in cingFixedResonanceDict:
          fixedResonances.append(cingFixedResonanceDict[atom])

        else:
          if atom.atomId[0] == 'C':
            isotopeCode = '13C'
          if atom.atomId[0] == 'N':
            isotopeCode = '15N'
          ccpnResNum = atom.sequenceId
          ccpnResidue = molSystem.findFirstChain().findFirstResidue(seqCode=ccpnResNum)
          ccpnAtoms = ccpnResidue.findFirstAtom(name=str(atom.atomId))
          fixedRes = nmrConstraintStore.newFixedResonance(name=atom.atomId, isotopeCode=isotopeCode)
          fixedResonances.append(fixedRes)
          atomSet = nmrConstraintStore.newFixedAtomSet(atoms=(ccpnAtoms,))
          nmrConstraintStore.newFixedResonanceSet(resonances=(fixedRes,), atomSets=(atomSet,))
      dihedralConstraint = dihedralRestraintList.newDihedralConstraint(resonances=fixedResonances)
      dihedralConstraint.newDihedralConstraintItem(upperLimit=cingRestraint.upper, lowerLimit=cingRestraint.lower)
    nmrCalcRun.newConstraintStoreData(name=cingRestraintList.name, ioRole='output',
                                        constraintStoreSerial=nmrConstraintStore.serial)


  print "Dihedral restraints assigned"


def loadRdcRestraints(cingRestraintList, nmrConstraintStore, cingFixedResonanceDict, molSystem):
  for cingList in cingRestraintList:
    name = cingList.name
    rdcRestraintList = nmrConstraintStore.newRdcConstraintList(name=name)
    for cingRestraint in cingList:
      fixedResonances = []
      #~ ccpnAtoms = []
      for i, atom in enumerate(cingRestraint.atomPairs):
        if atom in cingFixedResonanceDict:
          fixedResonances.append(cingFixedResonanceDict[atom])

        else:
          if atom.atomId[0] == 'C':
            isotopeCode = '13C'
          if atom.atomId[0] == 'N':
            isotopeCode = '15N'
          ccpnResNum = atom.sequenceId
          ccpnResidue = molSystem.findFirstChain().findFirstResidue(seqCode=ccpnResNum)
          ccpnAtoms = ccpnResidue.findFirstAtom(name=str(atom.atomId))
          fixedRes = nmrConstraintStore.newFixedResonance(name=atom.atomId, isotopeCode=isotopeCode)
          fixedResonances.append(fixedRes)
          atomSet = nmrConstraintStore.newFixedAtomSet(atoms=(ccpnAtoms,))
          nmrConstraintStore.newFixedResonanceSet(resonances=(fixedRes,), atomSets=(atomSet,))
      rdcConstraint = rdcRestraintList.newRdcConstraint(resonances=fixedResonances)
      rdcConstraint.newRdcConstraintItem(targetValue=cingRestraint.value, error=cingRestraint.error)

  print "RDC restraints assigned"

def loadChemShiftRestraints(cingRestraintList,originalProtFile,nmrConstraintStore,
                            atomToResonanceMap,fixedResonanceDict, molSystem):
  chemShiftRestraintList = nmrConstraintStore.newChemShiftConstraintList(name=originalProtFile)
  for cingRestraint in cingRestraintList:
    ringAtoms = ['CD1', 'CD2', 'CE1', 'CE2', 'HD1', 'HD2', 'HE1', 'HE2']
    ccpnResidueName = cingRestraint.atom.residueId[0] + cingRestraint.atom.residueId[1:].lower()
    ccpnRealAtomName = cingRestraint.atom.atomId

    if ccpnResidueName == 'Phe' and ccpnRealAtomName in ringAtoms:
      ccpnRealAtomName = ccpnRealAtomName[:2] + "*"

    if ccpnResidueName == 'Tyr' and ccpnRealAtomName in ringAtoms:
      ccpnRealAtomName = ccpnRealAtomName[:2] + "*"

    ccpnResNum = cingRestraint.atom.sequenceId
    ccpnChainCode = molSystem.findFirstChain(code=cingRestraint.atom.chainId)

    try:
      atomTuple = (ccpnChainCode, ccpnResNum, ccpnRealAtomName)
      resonance = atomToResonanceMap[atomTuple]
      fixedRes = fixedResonanceDict[resonance]
      chemShiftConstraint = chemShiftRestraintList.newChemShiftConstraint(resonance=fixedRes,
                                                                        targetValue=cingRestraint.value,
                                                                        error=cingRestraint.error)
    except KeyError:
      pass