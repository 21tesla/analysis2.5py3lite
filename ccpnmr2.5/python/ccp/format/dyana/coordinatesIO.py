"""
======================COPYRIGHT/LICENSE START==========================

coordinatesIO.py: Cartesian coordinate (.cor) I/O for Diana/Dyana files

Copyright (C) 2005-2010 Wim Vranken (European Bioinformatics Institute)

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

from ccp.format.dyana.generalIO import DyanaGenericFile

from memops.universal.Util import returnInt, returnFloat

from ccp.format.general.Constants import defaultSeqInsertCode

"""

Cartesian coordinate files in DG (.cor) format are used by DYANA for the input of start
conformations and for the output of the final conformations after the minimization. The
same format is used in both cases for the data lines: (6X,A5,I6,1X,A5,3F11.4)
corresponding to the atom name, the residue number and name, and the x-, y- and
z-coordinates of the atom in Angstrom units. For compatibility with other programs, the
first three lines are always comment lines even if they do not start with the # character;
further comment lines are not allowed. The three comment lines indicate the structure
index of the conformation and the final value of its target function, the program version
used, and the number of residues and atoms listed in the coordinate file, respectively. 

Optionally, the Cartesian atomic coordinates may be followed by the covalent
connectivities. In this case the format is (6X,A5,I6,1X,A5,3F11.4,4I6).

Note that on input, coordinates are only used to calculate all dihedral angles; the
structure will be rebuilt in DYANA according to the standard geometry obtained from the
residue library file. Therefore, conformations may be significantly changed if the
Cartesian coordinates do not imply exactly the bond lengths, bond angles and chirality of
the standard geometry! The same applies for Cartesian coordinates where the dihedral angle
of the peptide bonds are not exactly in the cis or trans conformation as defined in the
sequence input file.

An example output Cartesian coordinate file (in part) from DYANA follows:


 
Structure 001 from DYANA, f=4.953E-01 
Number of residues:    20 Number of atoms:     373 
    1 N         1 ALA       1.3450     0.0000     0.0000 
    2 H         1 ALA       1.8550     0.0001     0.8833 
    3 CA        1 ALA       2.1151    -0.0001    -1.2323 
    4 HA        1 ALA       1.4767    -0.1864    -2.0959 
    5 CB        1 ALA       2.8635     1.3270    -1.3717 
    6 HB1       1 ALA       3.4429     1.3320    -2.2949 
    7 HB2       1 ALA       3.5395     1.4630    -0.5277 
    8 HB3       1 ALA       2.1406     2.1427    -1.3822 
    9 QB        1 ALA       3.0410     1.6459    -1.4016 
   10 C         1 ALA       3.1856    -1.0908    -1.1595 
   11 O         1 ALA       2.9768    -2.0748    -0.4350 
   12 N         2 ARG       4.2924    -0.8983    -1.8990 
   13 H         2 ARG       4.3848    -0.0579    -2.4697 
   14 CA        2 ARG       5.3828    -1.8587    -1.9176 
   15 HA        2 ARG       6.0643    -1.5746    -2.7194 
   16 CB        2 ARG       6.1324    -1.8862    -0.5842 
   17 HB1       2 ARG       5.9446    -2.8408    -0.0927 
   18 HB2       2 ARG       5.7584    -1.0815     0.0488 
   19 QB        2 ARG       5.8516    -1.9612    -0.0220 
   ...
   
"""

#####################
# Class definitions #
#####################
     
class DyanaCoordinateFile(DyanaGenericFile):
  # Information on file level
  def initialize(self):
  
    self.modelCoordinates = {}
    self.chains = []

  def read(self, maxNum = None, verbose = 0):

    if verbose == 1:
      print "Reading %s coordinate file %s" % (self.format,self.name)
    
    #
    # One structure per file...
    #
    
    modelNum = 1
    self.modelCoordinates[modelNum] = []
    self.chains.append(DyanaChain(self.defaultMolCode))
    currentRefChainId = self.defaultMolCode
    oldSerial = None
    coordinateAddedToNewChain = False
    modelAdded = False
    
    #
    # Start reading file
    #
    
    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      cols = line.split()
      
      coordFound = False

      if len(cols) == 0 or self.patt['hash'].search(line):
        pass
        
      elif self.patt[self.format + "NewModel"].search(line):
      
        newModelNum = returnInt(self.patt[self.format + "NewModel"].search(line).group(1))
        
        # Just safeguard...
        if modelNum != newModelNum:
          modelAdded = True
          modelNum = newModelNum
          self.modelCoordinates[modelNum] = []
     
      elif self.patt[self.format + 'CoordinateAtomLine'].search(line):
        
        serial = line[5:11]
        atomName = line[12:17]
        resName = line[17:20].strip()
        chainCode = line[21:22]
        seqCode = line[22:26]
        insertionCode = line[26:27]
        x = line[30:38]
        y = line[38:46]
        z = line[46:54]

        coordFound = True

      elif self.patt[self.format + 'CoordinateInfoLine'].search(line):
        
        #(6X,A5,I6,1X,A5,3F11.4)
        serial = line[0:6]
        atomName = line[6:11]
        seqCode = line[11:17]
        insertionCode = line[17:18]
        resName = line[18:23].strip()
        x = line[23:34]
        y = line[34:45]
        z = line[45:56]
        chainCode = self.defaultMolCode
        
        coordFound = True
      
      #
      # Set the coordinate info if found... 
      #
      
      if coordFound:
      
        if oldSerial == None:
          oldSerial = serial
        
        if oldSerial > serial and not modelAdded:
          # This is for the 'old' format...
          modelNum += 1
          self.modelCoordinates[modelNum] = []
        
        modelAdded = False
        
        # Remove linker residue coordinates - irrelevant for CCPN, and cannot be stored.
        if resName in self.linkerResidueCodes:
          print "  Warning: Removed %s linker residue" % resName
          
          if coordinateAddedToNewChain:
            currentRefChainId = self.chainCodesString[len(self.chains) - 1]
            self.chains.append(DyanaChain(currentRefChainId))
            coordinateAddedToNewChain = False
            
        else:
        
          # Now add the coordinate info...
          coordinate = DyanaCoordinate(serial,atomName,seqCode,resName,chainCode,x,y,z,insertionCode = insertionCode,refChainId=currentRefChainId,verbose=False)
          if not None in (coordinate.x,coordinate.y,coordinate.z):
            self.modelCoordinates[modelNum].append(coordinate)
            coordinateAddedToNewChain = True

          oldSerial = serial
        
      line = fin.readline()

    fin.close()


  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing %s coordinate file %s" % (self.format,self.name)

    fout = open(self.name,'w')
    
    if len(self.modelCoordinates) > 1:
      print "Warning: can only handle 1 model per file for %s! Writing 1st model in list..." % self.format
    
    modelNums = self.modelCoordinates.keys()
    modelNums.sort()
    modelNum = modelNums[0]
    
    numResidues = 0
    numAtoms = 0
    
    oldSeqCode = None
    
    for coord in self.modelCoordinates[modelNum]:
    
      if coord.seqCode != oldSeqCode:
        numResidues += 1
        oldSeqCode = coord.seqCode
      
      numAtoms += 1
    
    fout.write("Structure 001 from CCPN, f=0.000" + self.newline)
    fout.write("Number of residues: %5d Number of atoms: %7d " % (numResidues,numAtoms) + self.newline)

    lineFormat = "%5d %-5s%6d%1s%-5s%11.4f%11.4f%11.4f "

    #
    # Write coordinates...
    #
    
    for coord in self.modelCoordinates[modelNum]:

      fout.write(lineFormat % (coord.serial,
                               coord.atomName,
                               coord.seqCode,
                               coord.insertionCode,
                               coord.resName,
                               coord.x,
                               coord.y,
                               coord.z))

      fout.write(self.newline)

    fout.close()
    
class DyanaChain:

  def __init__(self,chainId):

    self.chainId = chainId

class DyanaCoordinate:

  def __init__(self,serial,atomName,seqCode,resName,chainId,x,y,z,insertionCode = defaultSeqInsertCode,refChainId=None,verbose=False):
  
    self.serial = returnInt(serial,verbose=verbose)
    self.atomName = atomName.strip()
    self.seqCode = returnInt(seqCode,verbose=verbose)
    self.resName = resName.strip()
    self.chainId = chainId
    self.x = returnFloat(x,default=None,verbose=verbose)
    self.y = returnFloat(y,default=None,verbose=verbose)
    self.z = returnFloat(z,default=None,verbose=verbose)
    self.segId = ''
    self.insertionCode = insertionCode
    
    if refChainId:
      self.refChainId = refChainId
    else:
      self.refChainId = chainId
