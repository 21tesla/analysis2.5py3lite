"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for Diana/Dyana sequence files

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

from memops.universal.Util import returnInt

from ccp.format.general.formatIO import Sequence, SequenceElement
from ccp.format.general.Constants import defaultSeqInsertCode

# Diana/Dyana/Cyana data...
from ccp.format.cyana.cyanaLibParser import CyanaLibrary

#
# TODO: is this supposed to be the same format as XEasy?!?
#

#####################
# Class definitions #
#####################
      
class DyanaSequenceFile(DyanaGenericFile):
  # Information on file level
  def initialize(self):
  
    self.sequences = []

    self.cyanaLib = CyanaLibrary(version = self.version)
    
  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s sequence file %s" % (self.format,self.name)
    
    #
    # Parse the file
    #
    
    sequenceCols = []

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      cols = line.split()

      if len(cols) == 0 or self.patt['hash'].search(line):
        pass

      else:
        for col in cols:
          if not self.patt['hash'].search(col):
            sequenceCols.append(col)
          else:
            break

      line = fin.readline()

    fin.close()
    
    #
    # Loop through the sequence info from the file
    #
    # Will now remove the following (and start a new chain):
    #
    #
    # PL: Linker from a protein amino acid residue to a linker residue
    # LL: Linker residue with a virtual bond length of 1 A
    # LL2: Linker residue with a virtual bond length of 2 A
    # LL5: Linker residue with a virtual bond length of 5 A
    # LP: Linker from a linker residue to a protein amino acid residue 
    #

    seqCode = 0
    seqColNum = 0
    
    lineErrors = []

    self.sequences.append(DyanaSequence())
    chainResiduesAdded = False

    while (seqColNum < len(sequenceCols)):

      seqCol = sequenceCols[seqColNum]

      if not self.patt['onlydigit'].search(seqCol):
      
        isLinker = False
      
        # Removing linker residues - not for CCPN...
        if seqCol in self.linkerResidueCodes:
          if chainResiduesAdded:
            # Create a new chain!
            chainCode = self.chainCodesString[len(self.sequences)-1]
            self.sequences.append(DyanaSequence(chainCode=chainCode))
            print "  Warning: started new sequence '%s' based on linker residues." % chainCode
            chainResiduesAdded = False
            
          isLinker = True
        
        #
        # Check whether recognized...
        #
        
        if not self.cyanaLib.findResLabel(seqCol) and not self.cyanaLib.findResLabel(seqCol[1:]):
          lineErrors.append("  Warning: %s sequence element not recognized by standard CYANA library." % seqCol)

        #
        # Check if next element is a number
        #

        if seqColNum + 1 < len(sequenceCols) and self.patt['onlydigit'].search(sequenceCols[seqColNum + 1]):

          seqColNum += 1
          seqCode = returnInt(sequenceCols[seqColNum])

        else:

          seqCode += 1

        #
        # Set the data, if valid
        #
        
        if not isLinker:

          self.sequences[-1].elements.append(DyanaSequenceElement(seqCode,seqCol))
          self.sequences[-1].elements[-1].setFormatCode(self.sequences[-1].elements[-1].code3Letter)
         
          chainResiduesAdded = True

      else:

        lineErrors.append("  Error: %s sequence element %s doesn't fit format." % (self.format,seqCol))

      seqColNum += 1
      
    #
    # Validity check
    #
    
    fileReadOk = True
    
    numLineErrors = len(lineErrors)
    
    # Assume that if have a certain amount of errors, this file is not the right format...
    if numLineErrors > min((seqColNum * 0.1),5):
      self.sequences[-1].elements = []
      fileReadOk = False
    else:
      for lineError in lineErrors:
        print lineError
    
    return fileReadOk

  def readFromCoordinates(self,coordinateFile, verbose = 0):
  
    if verbose == 1:
      print "Extracting %s sequence from coordinate file %s" % (self.format,coordinateFile.name)

    seqCode = ""
    seqInsertCode = defaultSeqInsertCode
    
    chainCode = None
    
    modelNums = coordinateFile.modelCoordinates.keys()
    modelNums.sort()

    for coordinate in coordinateFile.modelCoordinates[modelNums[0]]:

      if chainCode != coordinate.refChainId:

        self.sequences.append(DyanaSequence(chainCode=coordinate.refChainId))
        chainCode = coordinate.refChainId


      if seqCode != coordinate.seqCode or seqInsertCode != coordinate.insertionCode:

        #
        # New residue/item
        # 
        
        seqCode = coordinate.seqCode
        seqInsertCode = coordinate.insertionCode

        fullSeqCode = str(seqCode) + seqInsertCode

        self.sequences[-1].elements.append(DyanaSequenceElement(fullSeqCode,coordinate.resName))
        self.sequences[-1].elements[-1].setFormatCode(coordinate.resName)

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing %s sequence file %s" % (self.format,self.name)

    if len(self.sequences) > 1:
      print "Warning: multiple sequences - writing to same file."        

    #
    # TODO TODO: have to fill in sequence gaps with linker residues!?>!
    #
    fout = open(self.name,'w')

    for sequence in self.sequences:

      #
      # Write three letter codes (one per line)
      #

      for residue in sequence.elements:

        if residue.hasCisPeptideBond:
          resText = 'c' + residue.formatCode
        else:
          resText = residue.formatCode

        fout.write("%-4s %4d" % (resText,residue.seqCode))
        fout.write(self.newline)

    fout.close()


DyanaSequence = Sequence

class DyanaSequenceElement(SequenceElement):
  
  def setResidueCode(self,*args):
    
    code3Letter = args[0]
    
    #
    # Determine if cis
    #
    
    if code3Letter[0] == 'c':
      self.setCisPeptideBond()
      code3Letter = code3Letter[1:]
    
    #
    # Get rid of +,- (currently anyway)
    #  
    
    self.formatCode = code3Letter.upper()
    
    if len(code3Letter) > 3:
      code3Letter = code3Letter[:3]
    
    self.code3Letter = code3Letter.upper()
