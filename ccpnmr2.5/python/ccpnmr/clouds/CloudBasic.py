"""
  These functions control the entire CLOUDS protocol

======================COPYRIGHT/LICENSE START==========================

CloudBasic.py: Part of the CcpNmr Clouds program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
import re

from math import sqrt

from os import listdir

from memops.universal.Util import returnInt, returnFloat

code1LetterDict = {'Ala':'A','Cys':'C','Asp':'D','Glu':'E','Phe':'F','Gly':'G',
                   'His':'H','Ile':'I','Lys':'K','Leu':'L','Met':'M','Asn':'N',
                   'Gln':'Q','Arg':'R','Ser':'S','Thr':'T','Val':'V','Trp':'W',
                   'Tyr':'Y','Asx':'B','Glx':'Z','Hyp':'O','Pro':'P'
                  }
                  
def filterPeakListWithCloud(argServer, peakList=None, pattern='t_intra_0\d+.pdb'):

  from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes

  if peakList is None:
    peakList = argServer.getPeakList()

  fileNames = getFileNamesFromPattern(pattern, '.')
  project   = peakList.root
  clouds    = getCloudsFromFile(fileNames, project)
  isotopes  = getSpectrumIsotopes(peakList.dataSource)

  hDims = []
  for i in range(len(isotopes)):
    if isotopes[i] == '1H':
      hDims.append(i)

  for peak in peakList.peaks:
    peakDims = peak.sortedPeakDims()

    hDim1 = peakDims[hDims[0]]
    hDim2 = peakDims[hDims[1]]
    
    contribs1 = hDim1.peakDimContribs
    contribs2 = hDim2.peakDimCOntribs
    
    if contribs1 and contribs2:
      resonances1 = [c.resonance for c in contribs1]   
      resonances2 = [c.resonance for c in contribs2]   

      coords1 = getResonanceCoords(clouds, resonances1)
      coords2 = getResonanceCoords(clouds, resonances2)
      dist    = getEnsembleCoordsDist(coord1,coord2)
      
def getResonanceCoords(clouds, resonances):

  ensemble = []
  for cloud in clouds:
    coords = []
    for resonance in resonances:
      coord = cloud.get(resonance)
      if coord:
        coords.append(coord)

    if coords:
      coord = averageCoord(coords)
      ensemble.append(coord)

  return ensemble

def getEnsembleCoordsDist(ensemble1, ensemble2):

  n    = len(ensemble1) 
  dist = 0.0

  for j in range(n):
    coords1 = ensemble1[j]
    coords2 = ensemble2[j]
    sum = 0.0
    for i in range(3):
      diff = coords2[i]-coords1[i]
      sum += diff * diff
    
    dist += sqrt(sum)
  dist /= float(n)
    
  return sqrt(sum)

def averageCoord(coordList):

  sum = [0.0, 0.0, 0.0]
  for coord in coordList:
    for i in range(3):
      sum[i] += coord[i]

  N = float(len(coordList))
  for i in range(3):
    sum[i] /= N

  return sum


def readCloudFile(fileName, hydrogenOnly=True):

  file = open(fileName, 'r')
  coords = []
  atoms  = []
  shifts = []
  serials = []
  
  for line in file.readlines():
    key = line[0:6].strip()
    if key == 'ATOM':
      a = line[13:16].strip()
      x = returnFloat(line[30:38])
      y = returnFloat(line[38:46])
      z = returnFloat(line[46:54])
      n = int(float(line[54:60].strip() or 0.0) * 100)
      s = float(line[60:66].strip() or 0.0)/10.0
            
      if hydrogenOnly:
        if (a[0]=='H') or (a == 'CH3'):
          coords.append([x,y,z])
          shifts.append(s)
          atoms.append(a)
          serials.append(n)
      else:
        coords.append([x,y,z])
        shifts.append(s)
        atoms.append(a)
        serials.append(n)

  return (coords, serials, shifts, atoms)

def getFileNamesFromPattern(pattern, directory):

  from re import match
  from memops.gui.MessageReporter import showWarning
  
  fileNames = []
  for fileName in listdir(directory):
    match = re.match(pattern, fileName)
    if match:
      fileNames.append(fileName)
  
  if not fileNames:
    showWarning('Failure','No file names were matched with pattern %s' % pattern)
  
  return fileNames
  

def getCloudsFromFile(fileNames, project, hydrogenOnly=False):
  # For multiple filenames the mean of the ensemble is considered

  shifts   = None
  serials  = []
  clouds   = []
  ensemble = []
  
  n = 0
  for fileName in fileNames:
    coords, serials, shifts, atomTypes = readCloudFile(fileName, hydrogenOnly=hydrogenOnly)
    ensemble.append(coords)
    clouds.append( {} )
    n += 1
  
  if shifts is None:
    # no matching files
    return []

  
  resonances = connectResonances(project, serials, shifts, atomTypes)

  for resonance in resonances:
    #resonance.name = getResonanceType(resonance)

    for i in range(n):
      clouds[i][resonance] = ensemble[i].pop(0)
  
  return clouds


def getResonanceType(resonance):

  name = 'r'
  if resonance.resonanceSet:
    index = list(resonance.resonanceSet.resonances).index(resonance)
    atoms = list(resonance.resonanceSet.atomSets)[index].atoms

    if len(atoms) == 3:
      name = 'MTH'
      
    else:  
      name = atoms[0].name
      if name == 'H1':
        name = 'H'
 
  elif resonance.assignNames:
    name = resonance.assignNames[0]
  
  return name

def writeTypedPdbCloud(atomCoordList, pdbFileName, resonances):

  fp = open(pdbFileName, 'w')
  N  = len(atomCoordList)
  if hasattr(atomCoordList[0],'x'):
    for n in range(N):
      a = atomCoordList[n]
      name = getResonanceType(resonances[n])
      serial = resonances[n].serial
      if resonances[n].shifts:
        shift = resonances[n].findFirstShift().value
      else:
        shift = 0.0
      
      fp.write('ATOM  %5d  %-3.3s HYD %5d    %8.3f%8.3f%8.3f%6.2f%6.2f\n' \
               % (n+1, name, n+1, a.x, a.y, a.z, float(serial)/100.0, shift*10.0))
  else:
    for n in range(N):
      [x,y,z] = atomCoordList[n]
      name = getResonanceType(resonances[n])
      serial = resonances[n].serial
      if resonances[n].shifts:
        shift = resonances[n].findFirstShift().value
      else:
        shift = 0.0
      
      fp.write('ATOM  %5d  %-3.3s HYD %5d    %8.3f%8.3f%8.3f%6.2f%6.2f\n' \
               % (n+1, name, n+1, x, y, z, float(serial)/100.0, shift*10.0))
    
  fp.close()
  


def connectResonances(project, serials, shifts, atomTypes, isotopeCodes=('1H',), shiftList=None):

  shiftLists = project.currentNmrProject.findAllMeasurementLists(className='ShiftList')
  
  if not shiftList:
    shiftListIndex = 0
    if not shiftLists:
      shiftList = project.newShiftList(unit='ppm')
    else:
      shiftList = shiftLists[0]
  else:
    shiftListIndex = list(shiftLists).index(shiftList)

  isotopes = { 'H':'1H','N':'15N','C':'13C' }
    
  shift1Dict = {}
  shift2Dict = {}
  shift3Dict = {}
  serialDict = {}
  for resonance in project.currentNmrProject.resonances:
    if resonance.shifts and (shiftListIndex < len(resonance.shifts)):
      shift = resonance.shifts[shiftListIndex].value
      shift1Dict['%6.4f' % shift] = resonance
      shift2Dict['%5.3f' % shift] = resonance
      shift3Dict['%4.2f' % shift] = resonance
      
    serialDict[resonance.serial] = resonance
  
  Cnew = 0
  Cfound = 0
  resonances = []
  for i in range(len(serials)):
    serial = serials[i]
    shift  = shifts[i]
  
    resonance = serialDict.get(serial)
    
    if resonance is None:
      print 'Serial %d missing' % serial
      resonance = shift1Dict.get('%6.4f' % shift)
    if resonance is None:
      resonance = shift2Dict.get('%5.3f' % shift)
    if resonance is None:
      resonance = shift3Dict.get('%4.2f' % shift)
    if resonance is None:
      Cnew += 1
      isotope   = isotopes.get(atomTypes[i][0],'1H')
      resonance = project.newResonance(isotopeCode=isotope,assignNames=[atomTypes[i],])
      if shift is not None:
        newShift  = shiftList.newShift(value=shift,resonance=resonance)
    
    else:
      Cfound += 1
    
    if not resonance.assignNames:
      resonance.assignNames = [atomTypes[i],]
    
    resonances.append(resonance)

  
  print "Total: %d, Found: %d, New: %d" % (len(resonances),Cfound,Cnew)
  return resonances

def makeConstraintsFromStructure(argServer, structure=None, atomRoots=None, threshold=5.5, tolerance=0.1, adcThreshold=7.0, atomTypes=('H',), adcAtoms=('H',)):
  # find all NMR visible atom pairs within a specified distance to 
  # generate distance constraints

  if not structure:
    structure = argServer.getStructure()

  from ccpnmr.analysis.core.AssignmentBasic import newResonance
  from ccpnmr.analysis.core.StructureBasic  import getAtomSetsDistance
  from ccpnmr.analysis.core.ConstraintBasic import getFixedResonance, makeNmrConstraintHead

  notVisible = {'Lys':["H''","HZ1","HZ2","HZ3"],
                'Arg':["H''","HH11","HH12","HH21","HH22","CZ"],
                'Ser':["H''","HG",],
                'Cys':["H''","HG",],
                'His':["H''","HD1","CG"],
                'Thr':["H''","HG1",],
                'Tyr':["H''","HH","CG"],
                'Phe':["H''","CG"],
                'Glu':["H''","HE2",],
                'Asp':["H''","HD2",],
                'Trp':["H''","CG","CD2","CE2"],
               }
  
  project = structure.root
  chain   = structure.findFirstCoordChain().chain
  
  atomSetDict    = {}
  anchorAtomSets = {}
  amideAtomSets  = []
  resonanceDict  = {}
  adcAtomSets    = {}
  
  print "Getting atomSets and resonances"
  for residue in chain.residues:
    atomSetDict[residue] = []
    exclude = notVisible.get(residue.ccpCode) or ["H''",]
    for atom in residue.atoms:
      if atom.name[0] not in atomTypes:
        continue
      if atom.name in exclude:
        continue
        
      atomSet = atom.atomSet
      if atomSet:
        if resonanceDict.get(atomSet) is None:
          if atomSet.resonanceSets:
            resonanceSet = atomSet.findFirstResonanceSet()
            i = list(resonanceSet.atomSets).index(atomSet)
            if i >= len(resonanceSet.resonances):
              i = 0
            resonance = resonanceSet.resonances[i]
          
          else:
            resonance = newResonance(project, isotopeCode='1H')

          resonanceDict[atomSet] = resonance
          atomSetDict[residue].append(atomSet)
      
        if atom.name in ('H','H1'):
          amideAtomSets.append( atomSet )
        elif (residue.ccpCode == 'PRO') and (atom.name == 'HD2'):  
          amideAtomSets.append( atomSet )
 
        if atomRoots:
          if atom.name in atomRoots:
            anchorAtomSets[atomSet] = 1
        else:
          anchorAtomSets[atomSet] = 1
        
        if adcAtoms and (atom.name in adcAtoms):
          adcAtomSets[atomSet] = 1
          
  print "Calculating distances"
  resonanceDists = []
  resonanceDists2 = []
  for amide1 in amideAtomSets:
    residue1 = amide1.finsFirstAtom().residue
    print "  %d %s" % (residue1.seqCode, residue1.ccpCode)
    
    for amide2 in amideAtomSets:
      residue2 = amide2.findFirstAtom().residue
      if amide1 is amide2:
        amideDist = 0.0
      else:
        amideDist = getAtomSetsDistance((amide1,), (amide2,), structure)

      if (amideDist is not None) and (amideDist <= 25.0):
        for atomSet1 in atomSetDict[residue1]:
          if anchorAtomSets.get(atomSet1) is None:
            continue
        
          for atomSet2 in atomSetDict[residue2]:
            if atomSet1 is atomSet2:
              continue
          
            dist = getAtomSetsDistance((atomSet1,), (atomSet2,), structure)
            if dist and dist <= threshold:
              resonanceDists.append( (dist,resonanceDict[atomSet1],resonanceDict[atomSet2]) )
      
      if adcAtoms:
        for atomSet1 in atomSetDict[residue1]:
          if anchorAtomSets.get(atomSet1) is None:
            continue
          if adcAtomSets.get(atomSet1) is None:
            continue
            
          for atomSet2 in atomSetDict[residue2]:
            if adcAtomSets.get(atomSet2) is None:
              continue
            if atomSet1 is atomSet2:
              continue

            dist = getAtomSetsDistance((atomSet1,), (atomSet2,), structure)
            if dist >= adcThreshold:
              resonanceDists2.append( (dist,resonanceDict[atomSet1],resonanceDict[atomSet2]) )

  print "Generating constraints"
  if resonanceDists:
    constraintHead = makeNmrConstraintHead(project)
    constraintList = constraintHead.newDistanceConstraintList()
    
    for dist, resonance1, resonance2, in resonanceDists:
      if resonance1 is resonance2: 
        continue
      fixedResonance1 = getFixedResonance(constraintHead, resonance1)
      fixedResonance2 = getFixedResonance(constraintHead, resonance2)
      delta           = tolerance * dist
      minDist         = max(0.0,dist - delta)
      maxDist         = dist + delta
      
      constraint = constraintList.newDistanceConstraint(weight=1.0, origData=dist, targetValue=dist,upperLimit=maxDist, lowerLimit=minDist, error=delta)
      item       = constraint.newDistanceConstraintItem(resonances=[fixedResonance1,fixedResonance2])

    if adcAtoms:
      constraintList2 = constraintHead.newDistanceConstraintList()
      # make ADCs
      for dist, resonance1, resonance2, in resonanceDists2:
        if resonance1 is resonance2: 
          continue
        dist = 75.0
        fixedResonance1 = getFixedResonance(constraintHead, resonance1)
        fixedResonance2 = getFixedResonance(constraintHead, resonance2)
        minDist         = 5.0
        maxDist         = 150
        delta           = maxDist - minDist
 
        constraint = constraintList2.newDistanceConstraint(weight=1.0, origData=dist, targetValue=dist,upperLimit=maxDist, lowerLimit=minDist, error=delta)
        item       = constraint.newDistanceConstraintItem(resonances=[fixedResonance1,fixedResonance2])

  print "Done"

def getMeanCloud(clouds):

  N = len(clouds)
  
  meanCloud = {}

  for resonance in clouds[0].keys():
    x,y,z = 0.0,0.0,0.0
    for i in range(N):
      x += clouds[i][resonance][0]
      y += clouds[i][resonance][1]
      z += clouds[i][resonance][2]
    x /= M
    y /= M
    z /= M
    meanCloud[resonance] = (x,y,z)
      
  return meanCloud

