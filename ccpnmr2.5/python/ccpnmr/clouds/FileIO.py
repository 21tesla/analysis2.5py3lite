"""
Simple IO routines for PDB style clouds coordinates files

======================COPYRIGHT/LICENSE START==========================

FileIO.py: Part of the CcpNmr Clouds program

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
from memops.universal.Util import returnInt, returnFloat

def readPdbCloud(pdbFileName):

  fileHandle = open(pdbFileName, 'r')

  # assumes order of atoms

  coordList = []

  for line in fileHandle.readlines():
    key = line[0:6].strip()
    if key == 'ATOM':
      #serial  = returnInt(line[6:11])
      #seqCode = returnInt(line[22:26])
      x = returnFloat(line[30:38])
      y = returnFloat(line[38:46])
      z = returnFloat(line[46:54])
      coordList.append([x,y,z])

  return coordList

def readTypedPdbCloud(pdbFileName):

  fileHandle = open(pdbFileName, 'r')

  # assumes order of atoms

  coordList = []
  atomList = []

  for line in fileHandle.readlines():
    key = line[0:6].strip()
    if key == 'ATOM':
      #serial  = returnInt(line[6:11])
      #seqCode = returnInt(line[22:26])
      a = line[13:15].strip()
      x = returnFloat(line[30:38])
      y = returnFloat(line[38:46])
      z = returnFloat(line[46:54])
      coordList.append([x,y,z])
      atomList.append(a)

  return (coordList, atomList)

def writePdbCloud(atomCoordList, pdbFileName):

  fp = open(pdbFileName, 'w')
  N  = len(atomCoordList)
  if hasattr(atomCoordList[0],'x'):
    for n in range(N):
      a = atomCoordList[n]
      fp.write('ATOM  %5d  H   HY1 %5d    %8.3f%8.3f%8.3f%6.2f%6.2f\n' \
               % (n+1, n+1, a.x, a.y, a.z, 1.0, 1.0))
  else:
    for n in range(N):
      [x,y,z] = atomCoordList[n]
      fp.write('ATOM  %5d  H   HY1 %5d    %8.3f%8.3f%8.3f%6.2f%6.2f\n' \
               % (n+1, n+1, x, y, z, 1.0, 1.0))
    
  fp.close()

def writeStructureCloud(structure, pdbFileName, coordIndex=0):

  atomSetDict = {}
  resonances = []
  atomCoordList = []
  for chain in structure.coordChains:
    for residue in chain.residues:
      for atom in residue.atoms:
        atomSet = atom.atom.atomSet
        if atomSet:
          coord = atom.coords[coordIndex]
          if atomSetDict.get(atomSet) is None:
            atomSetDict[atomSet] = []
          
          ll = [coord.x, coord.y, coord.z]
          atomSetDict[atomSet].append( ll ) 
  
  for atomSet in atomSetDict.keys():
    if atomSet.resonanceSets:

      x = 0.0
      y = 0.0
      z = 0.0
      n = 0.0
      for coord in atomSetDict[atomSet]:
        x += coord[0]
        y += coord[1]
        z += coord[2]
        n += 1.0
 
      x /= n
      y /= n
      z /= n
      atomCoordList.append((x,y,z))
    
      resonance = None
      for resonanceSet in atomSet.resonanceSets:
        if len(resonanceSet.resonances) == 1:
          resonance = resonanceSet.findFirstResonance()
          break
      
      if not resonance:
        resonanceSet = atomSet.findFirstResonanceSet()
        index = list(resonanceSet.atomSets).index(atomSet)
        if index >= len(resonanceSet.resonances):
          resonance = resonanceSet.resonances[-1]
        else:
          resonance = resonanceSet.resonances[index]
  
      resonances.append(resonance)

  writeTypedPdbCloud(atomCoordList, pdbFileName, resonances)  

def writeTypedPdbCloud(atomCoordList, pdbFileName, resonances):

  fp = open(pdbFileName, 'w')
  N  = len(atomCoordList)
  if hasattr(atomCoordList[0],'x'):
    for n in range(N):
      a = atomCoordList[n]
      name = resonances[n].name
      if resonances[n].shifts:
        shift = resonances[n].findFirstShift().value
      else:
        shift = 0.0
      
      fp.write('ATOM  %5d  %-3.3s HY1 %5d    %8.3f%8.3f%8.3f%6.2f%6.4f\n' \
               % (n+1, name, n+1, a.x, a.y, a.z, 1.0, shift))
  else:
    for n in range(N):
      [x,y,z] = atomCoordList[n]
      name = resonances[n].name
      if resonances[n].shifts:
        shift = resonances[n].finsFirstShift().value
      else:
        shift = 0.0
      
      fp.write('ATOM  %5d  %-3.3s HY1 %5d    %8.3f%8.3f%8.3f%6.2f%6.4f\n' \
               % (n+1, name, n+1, x, y, z, 1.0, shift))
    
  fp.close()
  
