"""
======================COPYRIGHT/LICENSE START==========================

coordinatesIO.py: I/O for GROMOS coordinate files

Copyright (C) 2012 Wim Vranken (Vrije Universiteit Brussel)

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

# Import general functions
from memops.universal.Util import returnInt, returnFloat, returnFloats
from ccp.format.gromos.generalIO import GromosFile

from ccp.format.general.Constants import defaultMolCode, defaultSeqInsertCode

#####################
# Class definitions #
#####################

class GromosCoordinateFile(GromosFile):

  def initialize(self):
  
    self.modelCoordinates = {}
    
    self.chains = []
    
    self.chainCodes = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789" # Let's hope this is enough!

  def read(self, verbose=False, maxNum = 999, ignoreResNames = None):
  
    """
MD of 2 waters, t= 0.0
    6
    1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
    1WATER  HW2    2   0.190   1.661   1.747  0.8085  0.3191 -0.7791
    1WATER  HW3    3   0.177   1.568   1.613 -0.9045 -2.6469  1.3180
    2WATER  OW1    4   1.275   0.053   0.622  0.2519  0.3140 -0.1734
    2WATER  HW2    5   1.337   0.002   0.680 -1.0641 -1.1349  0.0257
    2WATER  HW3    6   1.326   0.120   0.568  1.9427 -0.8216 -0.0244
   1.82060   1.82060   1.82060

Lines contain the following information (top to bottom):

    title string (free format string, optional time in ps after 't=')
    number of atoms (free format integer)
    one line for each atom (fixed format, see below)
    box vectors (free format, space separated reals), values: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y), the last 6 values may be omitted (they will be set to zero). Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0. 

This format is fixed, ie. all columns are in a fixed position. Optionally (for now only yet with trjconv) you can write gro files with any number of decimal places, the format will then be n+5 positions with n decimal places (n+1 for velocities) in stead of 8 with 3 (with 4 for velocities). Upon reading, the precision will be inferred from the distance between the decimal points (which will be n+5). Columns contain the following information (from left to right):

    residue number (5 positions, integer)
    residue name (5 characters)
    atom name (5 characters)
    atom number (5 positions, integer)
    position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places) 

Note that separate molecules or ions (e.g. water or Cl-) are regarded as residues. If you want to write such a file in your own program without using the GROMACS libraries you can use the following formats:

C format
    "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" 
Fortran format
    (i5,2a5,i5,3f8.3,3f8.4) 
    """

    if verbose:
      print "Reading %s coordinate file %s" % (self.format,self.name)
      
    if not ignoreResNames:
      ignoreResNames = ('SOL',)

    #
    # Read coordinates and other info
    #
    
    currentChainCode = None
    oldResNum = None

    fin = open(self.name, 'rU')
    lines = fin.readlines()
    fin.close()
    
    self.title = lines[0]
    self.numAtoms = returnInt(lines[1])
    self.boxCoords = returnFloats(lines[-1].split())
    
    modelNum = 1
    self.modelCoordinates[modelNum] = []
    
    for line in lines[2:-1]:

      resNum = returnInt(line[:5])
      resName = line[5:10].strip()
      atomName = line[10:15].strip()
      atomSerial = returnInt(line[15:20])
      
      # Don't read solvent by default, can be changed...
      if resName in ignoreResNames:
        continue
      
      coordAndVelocities = line[20:].split()
      
      (x,y,z) = returnFloats(coordAndVelocities[:3])
      
      if len(coordAndVelocities) > 3:
        velocities = returnFloats(coordAndVelocities[3:])
      else:
        velocities = []
        
      # Initialise this at start of file reading
      if currentChainCode == None:
        currentChainCode = self.chainCodes[0]
        oldResNum = resNum
        self.chains.append(GromosChain(currentChainCode))
        
      # New chain code if residue number jump
      if resNum < oldResNum or oldResNum + 1 < resNum:
        currentChainCode = self.chainCodes[self.chainCodes.index(currentChainCode) + 1] 
        oldResNum = resNum
        self.chains.append(GromosChain(currentChainCode))
                  
      self.modelCoordinates[modelNum].append(GromosCoordinate(atomSerial,atomName,resName,currentChainCode,resNum,x,y,z,velocities = velocities))
      
      oldResNum = resNum

    fin.close()

  """
  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing %s coordinate file %s" % (self.format,self.name)

    fout = open(self.name,'w')

    #
    # TODO: CURRENTLY NO HEADER!
    #
    
    multipleModels = 0
    
    if len(self.modelCoordinates) > 1:
      multipleModels = 1
      
    modelNums = self.modelCoordinates.keys()
    modelNums.sort()
      
    for modelNum in modelNums:
        
      oldChainId = self.modelCoordinates[modelNum][0].chainId

      if multipleModels:
      
        fout.write("MODEL %7s" % modelNum + self.newline)
        
      for coord in self.modelCoordinates[modelNum]:
      
        if oldChainId != coord.chainId:
          fout.write("TER" + self.newline)
      
        atomText = 'ATOM'
        atomFormat = " %-3s"
       
        if self.patt['onlydigit'].search(coord.atomName[0]) or len(coord.atomName) == 4:
          
          atomFormat = "%-4s"
      
        lineFormat = "%-6s%5d " + atomFormat + "%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s"
      
        fout.write(lineFormat %
        
                    (atomText,coord.serial,coord.atomName,'',coord.resName,coord.chainId,
                     coord.seqCode,coord.insertionCode,coord.x,coord.y,coord.z,1.00,0.00,coord.segId,'',''))
                     
        fout.write(self.newline)
            
      fout.write("TER" + self.newline)
      
      if multipleModels:
      
        fout.write("ENDMDL" + self.newline)

    fout.close()
  """

class GromosChain:

  def __init__(self,chainCode):

    self.chainCode = self.chainId = chainCode

class GromosCoordinate:

  def __init__(self,atomSerial,atomName,resName,chainCode,resNum,x,y,z, velocities = None):
  
    self.serial = atomSerial
    self.atomName = atomName
    self.resName = resName
    self.seqCode = resNum
    self.x = x
    self.y = y
    self.z = z
    self.insertionCode = defaultSeqInsertCode
    
    self.velocities = velocities

    self.chainCode = self.chainId = chainCode
    
    # On-the-fly mapping, hacky
    if resName == 'ILE' and atomName == 'CD':
      self.atomName = 'CD1'
    
    elif atomName == 'OC1':
      self.atomName = "O"
    elif atomName == 'OC2':
      self.atomName = 'OXT'
    
    
