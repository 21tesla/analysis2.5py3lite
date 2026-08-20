#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for Auremol sequence files

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
from memops.universal.Util import returnInt
from ccp.format.auremol.generalIO import AuremolFile
from memops.universal.Io import getTopDirectory
from ccp.format.general.formatIO import Sequence, SequenceElement
from ccp.format.general.Constants import defaultSeqInsertCode

#####################
# Class definitions #
#####################
      
class AuremolSequenceFile(AuremolFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.sequences = []

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s sequence file %s" % (self.format,self.name)

    self.sequences.append(AuremolSequence())

    seqCode = 1
    codeLen = 3
    seqRead = 0

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      cols = line.split()

      if len(cols) == 0 or self.patt['hash'].search(line):
        pass

      elif cols[0] == 'NAME:':
      
        self.sequences[-1].molName = cols[1]

      elif cols[0] == 'CODE:':
      
        codeLen = returnInt(cols[1])
        
        if codeLen not in (1,3):
          print "Error: invalid sequence code letter length %d." % codeLen
          return

      elif cols[0] == 'NUMBER:':
      
        seqLen = returnInt(cols[1])
      
      elif cols[0] == 'OFFSET:':
    
        seqCode += returnInt(cols[1])

      elif cols[0] == 'START:':
      
        seqRead = 1

      elif cols[0] == 'END:':
      
        seqRead = 0
        
      elif seqRead:
      
        for codeLetter in cols:
        
          if codeLetter == '-':
            break
            
          if codeLen == 3:
          
            self.sequences[-1].elements.append(AuremolSequenceElement(seqCode,code3Letter = codeLetter))
           
          elif codeLen == 1:
          
            self.sequences[-1].elements.append(AuremolSequenceElement(seqCode,code1Letter = codeLetter))

          seqCode += 1

      else:

        print "Warning: unknown keyword %s!" % cols[0]

      line = fin.readline()

    fin.close()

  def readFromCoordinates(self,coordinateFile, verbose = 0):
  
    if verbose == 1:
      print "Extracting %s sequence from coordinate file %s" % (self.format,coordinateFile.name)

    seqCode = ""
    seqInsertCode = defaultSeqInsertCode
    segId = -1
    residueName = ""
    
    modelNums = coordinateFile.modelCoordinates.keys()
    modelNums.sort()

    for coordinate in coordinateFile.modelCoordinates[modelNums[0]]:

      if segId != coordinate.segId:

        #
        # New sequence
        #

        self.sequences.append(AuremolSequence(chainCode = coordinate.segId))
        segId = coordinate.segId
        seqCode = ""
        seqInsertCode = defaultSeqInsertCode

      if seqCode != coordinate.seqCode or seqInsertCode != coordinate.insertionCode:

        #
        # New residue/item
        # 
        
        seqCode = coordinate.seqCode
        seqInsertCode = coordinate.insertionCode

        fullSeqCode = str(seqCode) + seqInsertCode

        self.sequences[-1].elements.append(AuremolSequenceElement(fullSeqCode,code3Letter = coordinate.resName))

"""
  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing %s sequence file %s" % self.name

    if len(self.sequences) > 1:
      print "Warning: multiple sequences - writing to same file."        

    fout = open(self.name,'w')

    for sequence in self.sequences:

      #
      # Write three letter codes (one per line)
      #

      for residue in sequence.elements:

        fout.write(" %3s" % (string.upper(residue.code3Letter)))
        fout.write(self.newline)

    fout.close()
"""
    
class AuremolSequence(Sequence):

  def setFormatSpecific(self,*args,**keywds):
  
    pass
 
class AuremolSequenceElement(SequenceElement):

  def setFormatSpecific(self,code3Letter = None, code1Letter = None):
    
    if code3Letter:
      self.code3Letter = string.upper(code3Letter)
      
    elif code1Letter:
      self.code1Letter = string.upper(code1Letter)

###################
# Main of program #
###################
                                                      
if __name__ == "__main__":

  files = [('../../reference/auremol/csp.seq',None),
           ('../../reference/auremol/WERNERCSP14_7_IUPAC.pdb','coords')
          ]

  for (inFile,readType) in files:
    
    sequenceFile = AuremolSequenceFile(os.path.join(getTopDirectory(),inFile))
    
    if not readType:
      sequenceFile.read(verbose = 1)
    else:
      from coordinatesIO import AuremolCoordinateFile
      coordFile = AuremolCoordinateFile(os.path.join(getTopDirectory(),inFile))
      coordFile.read(verbose = 1)
      sequenceFile.readFromCoordinates(coordFile,verbose = 1)

    for seq in sequenceFile.sequences:
      for seqel in seq.elements:
        if hasattr(seqel,'code3Letter'):
          print seqel.seqCode, seqel.code3Letter
        else:
          print seqel.seqCode, seqel.code1Letter

