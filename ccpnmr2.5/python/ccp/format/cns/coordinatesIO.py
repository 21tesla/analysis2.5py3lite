"""
======================COPYRIGHT/LICENSE START==========================

coordinatesIO.py: I/O for XPLOR/CNS coordinate files

Copyright (C) 2005-2008 Wim Vranken (European Bioinformatics Institute)

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

import os

# Import general functions
from memops.universal.Util import returnInt, returnFloat
from ccp.format.cns.generalIO import CnsGenericFile
from memops.universal.Io import getTopDirectory

from ccp.format.general.Util import getRegularExpressions
from ccp.format.general.Constants import defaultMolCode, defaultSeqInsertCode
from ccp.format.general.Constants import defaultAltLoc

#####################
# Class definitions #
#####################

class CnsCoordinateFile(CnsGenericFile):

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
      
      if self.patt['pdbNewModel'].search(line): # Not necessarily there!! Is empty if so
        modelNum = returnInt(self.patt['pdbNewModel'].search(line).group(1))

        if modelNum > maxNum:
          atomRead = 0
          break

        self.modelCoordinates[modelNum] = []

      elif atomRead == 1 and self.patt['pdbAllAtom'].search(line):

        # Case one model only
        if modelNum == '':
          modelNum = 1
          self.modelCoordinates[modelNum] = []

        serial = line[6:11]
        atomName = line[12:16]
        altLoc = line[16:17]
        resName = line[17:20]
        seqCode = line[22:26]
        insertionCode = line[26:27]
        x = line[30:38]
        y = line[38:46]
        z = line[46:54]
        occupancy = line[54:60]
        tempFactor = line[60:66]
        
        # CNS USES SegId rather than chainId!!
        #chainId = line[21:22]
        segId = line[72:76]

        self.modelCoordinates[modelNum].append(CnsCoordinate(
         serial,atomName,resName,segId,seqCode,x,y,z,
         insertionCode=insertionCode, altLoc=altLoc)
        )

        if segId != oldSegId:
          
          oldSegId = segId
          
          foundChain = 0
          
          for chain in self.chains:
            if chain.segId == segId:
              foundChain = 1
              break
              
          if not foundChain:  

            self.chains.append(CnsChain(segId))
            
        # TODO: add in occupancy, tempFactor, ... fields if needed!

      line = fin.readline()

    fin.close()

  def write(self,endStatement = 'END',verbose = False,writeString=False):

    if verbose == 1:
      print "Writing %s coordinate file %s" % (self.format,self.name)

    # Provision for writing out a string...
    if writeString:
      from memops.universal.Io import RedirectToString
      fout = RedirectToString()
    else:
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
          oldChainId = coord.chainId
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
            
      fout.write(endStatement + self.newline)
      
      if multipleModels:
      
        fout.write("ENDMDL" + self.newline)

    #
    # Make sure string version is accessible
    #
    
    if writeString:
      self.coordFileString = fout.string
   
    fout.close()    

class CnsChain:

  def __init__(self,segId):

    self.segId = segId
    self.chainId = segId.strip()
    if not self.chainId:
      self.chainId = defaultMolCode

class CnsCoordinate:

  def __init__(self,serial,atomName,resName,segId,seqCode,x,y,z, 
               insertionCode=defaultSeqInsertCode, altLoc=defaultAltLoc, chainId = None):
  
    self.serial = returnInt(serial)
    self.atomName = atomName.strip()
    self.resName = resName.strip()
    self.seqCode = returnInt(seqCode)
    self.x = returnFloat(x)
    self.y = returnFloat(y)
    self.z = returnFloat(z)
    self.insertionCode = insertionCode
    self.altLoc = altLoc

    self.segId = segId.strip()
    
    if chainId == None:
      self.chainId = self.segId
      if not self.chainId:
        self.chainId = defaultMolCode
    else:
      self.chainId = chainId
    
###################
# Main of program #
###################

if __name__ == "__main__":  

  #file = os.path.join(getTopDirectory(),'../../reference/xplor/neuhaus.pdb')
  #file = os.path.join(getTopDirectory(),'../../reference/ccpNmr/davidLMB/041029/complex_01.pdb')
  file = os.path.join(getTopDirectory(),'../../reference/ccpNmr/aartUtrecht/041118/1a7f/coord_reprotonated/1a7f_1.pdb')
  cnsFile = CnsCoordinateFile(file)
  
  cnsFile.read(maxNum = 1)

  for mc in cnsFile.modelCoordinates:
  
    for coo in cnsFile.modelCoordinates[mc]:
      print coo.serial, coo.resName, coo.chainId
      
  for chain in cnsFile.chains:
    print chain, chain.chainId

  cnsFile.name = 'local/test.coord'  
  cnsFile.write()
  
