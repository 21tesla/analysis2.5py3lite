LICENSE = """
======================COPYRIGHT/LICENSE START==========================

StructureBasic.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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


# NBNB entire file is DEPRECATED! and contents have been copied elsewhere.
# Use ccp.lib.StructureIo, StructureLib or DataConvertLib instead

# NB If you continue to use this, please make any bug fixes in ccp.lib as well.

from math import sqrt, cos, sin, atan2
from os import path, mkdir

#from ccp.general.Geometry import calcTorsionAngleRadians, calcTorsionAngleDegrees
#from ccp.general.Io import getStdChemComps
#from ccp.general.Io import getChemComp
#from ccp.util.Molecule import makeChain, addMolResidues, makeMolecule

#from ccpnmr.analysis.core.MoleculeBasic import findMatchingChains, getLinkedResidue

try:
  from memops.gui.MessageReporter import showOkCancel, showYesNo
  from memops.gui.MessageReporter import showError, showWarning
except ImportError:
  from memops.universal.MessageReporter import showOkCancel, showYesNo
  from memops.universal.MessageReporter import showError, showWarning

from memops.universal.Util import returnInt, returnFloat
from memops.universal.Geometry import matrixMultiply
                  
BACKBONE_ATOMS = {'protein':('N','C','CA'),
                  'RNA':("OP1","P","O3'","O5'","C3'","C4'","C5'"),
                  'DNA':("OP1","P","O3'","O5'","C3'","C4'","C5'"),
                  'carbohydrate':("C1","C2","C3","C4","C5","C6","C7","C8"),}


# priority order of naming systems
namingSystemPriorityOrder = ['PDB','IUPAC','PDB_REMED','BMRB','XPLOR',
                             'CYANA2.1','DIANA', 'GROMOS','MSD','SYBYL',
                             'UCSF','AQUA','DISGEO','DISMAN', 'MOLMOL','MSI']

TWOPI = 6.2831853071795864


# Functions moved elsewhere
#from ccp.lib.StructureIo import makeStructureDictFromPdb # no longer needed.
#from ccp.lib.StructureIo import makeStructureDictFromRoughPdb # no longer needed.
from ccp.lib.StructureIo import makePdbFromStructure
from ccp.lib.StructureIo import makeStructureEnsemble, getStructureFromFile

from ccp.lib.DataConvertLib import getBestNamingSystem, getBestNamingSystemCC
from ccp.lib.DataConvertLib import getBestChemComp
from ccp.lib.DataConvertLib import getBestMolType
from ccp.lib.DataConvertLib import findMatchingMolSystemAtom

from ccp.lib.StructureLib import makeEmptyEnsembleCopy, copyModelToEnsemble, makeEnsemble
from ccp.lib.StructureLib import alignStructures, compareEnsembles
from ccp.lib.StructureLib import getResiduePhiPsi, getAtomSetCoords
from ccp.lib.StructureLib import getAtomSetsDihedral, getAtomSetsDistance


# Obsolete, unused, and incorrect
def alignCoordinates(coords1, coords0, allCoords, W):
  """Align on two sets of coordinates (may be sub sets of whole)
             and then update all coordinates given this new position.
             Coordinates are weighted for alignment.
  .. describe:: Input
  
  List of List of Floats (x,y,z), List of List of Floats (reference x,y,z),
             List of List of Floats (x,y,z), List of Floats (weights)

  .. describe:: Output

  List of List of Floats (aligned x,y,z),
             List of List of Floats (all updated x,y,z),
             Float (fitting score)
  NBNB consider refactoring to avoid use of Coords
  """
  from ccp.c.StructUtil import alignCoordinates as cAlignCoordinates
  
  print ("DEPRECATED, function alignCoordinates should not be used")
  
  rMat, error = cAlignCoordinates(coords1, coords0)
  # coords1 are modified in-place

  for n in range(len(allCoords)):
    [allCoords[n][0],allCoords[n][1],allCoords[n][2]] = matrixMultiply(rMat,allCoords[n])
    
  return coords1, allCoords, error


# Obsolete and unused
def getMeanStrucCoords(structureCoords):
  """Find the mean position of input coordinates of multiple structures.
  .. describe:: Input
  
  List of List of List of Floats (x,y,z per atom per structure)

  .. describe:: Output

  List of List of Floats (mean x,y,z per atom)
  """
  print ("DEPRECATED, function getMeanStrucCoords should not be used")

  meanCoords = []
  N = len(structureCoords)
  
  for coord in range(len(structureCoords[0])):
    sx = 0
    sy = 0
    sz = 0
    for s in range(N):
      sx  += structureCoords[s][coord][0] 
      sy  += structureCoords[s][coord][1] 
      sz  += structureCoords[s][coord][2]

    meanCoords.append([sx/N,sy/N,sz/N])

  return meanCoords



# Obsolete and unused
def getRmsd(structureCoords):
  """Find the root mean square devatation for a list of structures
             (their coordinates in lists)
  .. describe:: Input
  
  List of List of List of Floats (x,y,z per atom per structure)

  .. describe:: Output

  Float
  """
  
  print ("DEPRECATED, function getRmsd should not be used")

  N  = 0
  d2 = 0
  for coord in range(len(structureCoords[0])):
    d2i = 0
    Ni = 0
    for s1 in range(len(structureCoords)-1):
      for s2 in range(s1+1,len(structureCoords)):
        dx  = structureCoords[s1][coord][0] - structureCoords[s2][coord][0]
        dy  = structureCoords[s1][coord][1] - structureCoords[s2][coord][1]
        dz  = structureCoords[s1][coord][2] - structureCoords[s2][coord][2]
        d2  += (dx*dx) + (dy*dy) + (dz*dz)
        d2i += (dx*dx) + (dy*dy) + (dz*dz)
        Ni += 1
        N  += 1
    d2i = sqrt(d2i/Ni)
    #print "RMSD>", coord, d2i

  return sqrt(d2/N)


# Obsolete and unused
def centerCoordinates(coords):
  """More the centre of mass of input coordinates to the origin.
  .. describe:: Input
  
  List of List of Floats (x,z,y coords)

  .. describe:: Output

  List of List of Floats (x,z,y coords)
  NBNB consider refactoring to avoid use of Coords
  """
  print ("DEPRECATED, function centerCoordinates should not be used")

  if coords:
    centreOfMass = getCentreOfMass(coords)
    vector = [-c for c in centreOfMass]
    moveCoords(coords, vector)
    return coords


# Obsolete and unused
def centerStructures(structures):
  """Move the input structures' centre of mass  to the origin.
  .. describe:: Input
  
  List of MolStructure.StructureEnsembles

  .. describe:: Output

  None
  NBNB consider refactoring to avoid use of Coords
  """
  print ("DEPRECATED, function centerCoordinates should not be used")

  for structure in structures:
    coords = getStructureCoordinates(structure)
    centerCoordinates(coords)
 
# Obsolete and unused
def moveCoords(coords, vector):
  """Move data model coordinates according to a vactor.
  .. describe:: Input
  
  List of MolStructure.Coords, List of Floats (translation vector)

  .. describe:: Output

  List of MolStructure.Coords
  NBNB consider refactoring to avoid use of Coords
  """
  print ("DEPRECATED, function moveCoords should not be used")

  (dx,dy,dz) = vector
  for coord in coords:
    coord.x += dx
    coord.y += dy
    coord.z += dz

  return coords

# Obsolete and unused
def getMeanCoords(coordList):
  """Find the mean position of data model coordinates.
  .. describe:: Input
  
  List of MolStructure.Coords

  .. describe:: Output

  List of Floats (x,y,z coords of mean)
  """
  print ("DEPRECATED, function getMeanCoords should not be used")

  N = len(coordList)
  x = 0
  y = 0
  z = 0
  for coord in coordList:
    x += coord.x
    y += coord.y
    z += coord.z
    
  if N > 0:
    x /= N
    y /= N
    z /= N
   
  return (x,y,z) 

# Obsolete and unused
def getCentreOfMass(coords):
  """Find the center of mass of data model coordinates. (Uss atomic weights)
  .. describe:: Input
  
  List of MolStructure.Coords

  .. describe:: Output

  List of Floats (x,y,z centre of mass)
  NBNB consider refactoring to avoid use of Coords
  """
  print ("DEPRECATED, function getCentreOfMass should not be used")

  (cx,cy,cz) = [0,0,0]

  M = 0
  for coord in coords:
    atom = coord.atom
    mass = atom.atom.chemAtom.chemElement.mass
    #coord = atom.findFirstCoord()  # surely you should use coord RHF June 2011
    cx += mass * coord.x
    cy += mass * coord.y
    cz += mass * coord.z
    M  += mass 
    
  cx = cz/M
  cy = cz/M
  cz = cz/M

  return (cx,cy,cz)

