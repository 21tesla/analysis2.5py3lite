#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

coordinatesIO.py: I/O for CHARMM coordinate files

Copyright (C) 2005 Wim Vranken (European Bioinformatics Institute)

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

import os, string

# Import general functions
from memops.universal.Util import returnInt, returnFloat
from ccp.format.charmm.generalIO import CharmmGenericFile
from memops.universal.Io import getTopDirectory

from ccp.format.general.Util import getRegularExpressions
from ccp.format.general.Constants import defaultMolCode, defaultSeqInsertCode

#####################
# Class definitions #
#####################

class CharmmCoordinateFile(CharmmGenericFile):

  def initialize(self):
  
    self.modelCoordinates = {}
    
    self.patt.update(getRegularExpressions('pdb'))
    
    self.chains = []

  def read(self,maxNum = 999):

    atomRead = 1
    modelNum = ''
    oldSegId = None

    #
    # Read coordinates and other info
    #

    fin = open(self.name, 'rU')
    line = fin.readline()

    while line:
      
      # Is this relevant for CHARMM???
      if self.patt['pdbNewModel'].search(line): # Not necessarily there!! Is empty if so
        modelNum = returnInt(self.patt['pdbNewModel'].search(line).group(1))

        if modelNum > maxNum:
          atomRead = 0
          break

        self.modelCoordinates[modelNum] = []

      elif atomRead == 1 and self.patt[self.format + 'NumberAtoms'].search(line):
        
        numberAtoms = returnInt(line)
        
      elif atomRead == 1 and self.patt[self.format + 'AtomLine'].search(line):

        # Case one model only
        if modelNum == '':
          modelNum = 1
          self.modelCoordinates[modelNum] = []

        serial = line[:5]
        seqCode = line[5:10]
        resName = line[10:16]
        atomName = line[16:20]
        x = line[20:30]
        y = line[30:40]
        z = line[40:50]
        # charmm USES SegId rather than chainId!!
        segId = line[50:55]
        segSeqCode = line[55:60]
        tempFactor = line[60:]
        
        #
        # Parsing consistency check...
        #
        
        coordinate = CharmmCoordinate(serial,atomName,resName,segId,seqCode,segSeqCode,x,y,z,verbose=False)
        
        if None not in (coordinate.x,coordinate.y,coordinate.z):
          self.modelCoordinates[modelNum].append(coordinate)

        if segId != oldSegId:
          
          oldSegId = segId
          
          foundChain = 0
          
          for chain in self.chains:
            if chain.segId == segId:
              foundChain = 1
              break
              
          if not foundChain:  

            self.chains.append(CharmmChain(segId))
            
        # TODO: add in occupancy, tempFactor, ... fields if needed!

      line = fin.readline()

    fin.close()

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing %s coordinate file %s" % (self.format,self.name)

    fout = open(self.name,'w')

    #
    # TODO: only writing out one model at a time per file...
    #
    
    modelNums = self.modelCoordinates.keys()
    modelNums.sort()
    
    modelCoords = self.modelCoordinates[modelNums[0]]
    
    fout.write("%5d" % len(modelCoords))
    fout.write(self.newline)

    for coord in modelCoords:

      lineFormat = "%5s%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-5d%9.5f"

      fout.write(lineFormat %

                  (coord.serial,coord.seqCode,coord.resName,coord.atomName,coord.x,
                   coord.y,coord.z,coord.segId,coord.segSeqCode,0.00))

      fout.write(self.newline)
            
    fout.close()

class CharmmChain:

  def __init__(self,segId):

    self.segId = segId
    self.chainId = string.strip(segId)
    if not self.chainId:
      self.chainId = defaultMolCode

class CharmmCoordinate:

  def __init__(self,serial,atomName,resName,segId,seqCode,segSeqCode,x,y,z, insertionCode = defaultSeqInsertCode,verbose=False):
  
    self.serial = returnInt(serial,verbose=verbose)
    self.atomName = string.strip(atomName)
    self.resName = string.strip(resName)
    self.seqCode = returnInt(seqCode,verbose=verbose)
    self.segSeqCode = returnInt(segSeqCode,verbose=verbose)
    self.x = returnFloat(x,default=None,verbose=verbose)
    self.y = returnFloat(y,default=None,verbose=verbose)
    self.z = returnFloat(z,default=None,verbose=verbose)
    self.insertionCode = insertionCode

    self.chainId = self.segId = string.strip(segId)
    if not self.chainId:
      self.chainId = defaultMolCode

###################
# Main of program #
###################

if __name__ == "__main__":  

  file = os.path.join(getTopDirectory(),'../../reference/charmm/ref03_min.crd')
  charmmFile = CharmmCoordinateFile(file)
  
  charmmFile.read(maxNum = 1)

  for mc in charmmFile.modelCoordinates:
  
    for coo in charmmFile.modelCoordinates[mc]:
      print coo.serial, coo.resName, coo.chainId
      
  for chain in charmmFile.chains:
    print chain, chain.chainId

  charmmFile.name = 'local/test.coord'  
  charmmFile.write()
  
