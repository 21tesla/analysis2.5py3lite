"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for CSI chemical shift files

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

import string

# Import general functions
from memops.universal.Util import returnFloats, returnFloat
from memops.universal.Util import returnInt
from ccp.format.csi.generalIO import CsiGenericFile

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class CsiChemShiftFile(CsiGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.chemShifts = []
    
    self.chemShiftsBySeqCodeAndAtomName = {}
    self.seqCodeLabels = {}
    
    self.headerAtomNames = []
        
  def read(self,isOutputFile = False, headerCols = None, verbose = 0):
  
    if verbose == 1:
      print "Reading %s chemical shift list %s" % (self.format,self.name)

    fin = open(self.name, 'rU')

    headerCols = []
    
    line = fin.readline()
    
    while line:

      if self.patt['%sComment' % self.format].search(line) or self.patt['emptyline'].search(line):

        line = fin.readline()
        continue
              
      #
      # Get the info... 
      #
      
      cols = line.split()
      
      if not headerCols:
        headerCols = line.split()
        colLen = len(headerCols)
        self.headerAtomNames = headerCols[2:]
        
      elif len(cols) == colLen:
        
        seqCode = returnInt(cols[0],default=None,verbose=False)
        resLabel = cols[1]
        
        if seqCode != None:
          self.seqCodeLabels[seqCode] = resLabel
          
          for i in range(2,colLen):
            atomName = headerCols[i]
            value = returnFloat(cols[i],default=None,verbose=False)
            
            # Ignore if value is None - is unlikely for these values
            if value != None:
              
              # Need to handle GLY differently...
              if resLabel == 'G' and atomName == 'HA':
                atomName = atomName + '*'
              elif atomName == 'CO':
                atomName = 'C'
  
              self.chemShifts.append(CsiChemShift(self,value,atomName,seqCode,resLabel,self.defaultMolCode))

      line = fin.readline()

    fin.close()

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing %s chemical shift list %s" % (self.format,self.name)


    fout = open(self.name,'w')

    #
    # Write out chem shifts
    #

    fout.write(">CSI input file written by FormatConverter" + self.newline)
    fout.write(self.newline)
    
    # Write header    
    fout.write("%-4s%-3s" % ("#","AA"))
    for headerAtomName in self.headerAtomNames:
      fout.write("   %-4s   " % headerAtomName)
    fout.write(self.newline)
    
    # This should be ints
    seqCodes = self.chemShiftsBySeqCodeAndAtomName.keys()
    seqCodes.sort()
    
    if seqCodes:
      
      seqCodeRange = range(seqCodes[0],seqCodes[-1]+1)
      
      for seqCode in seqCodeRange:
        
        if not self.seqCodeLabels.has_key(seqCode):
          print "  Warning: no CSI output for sequence code %d - no information available." % seqCode
          continue
        
        resLabel = self.seqCodeLabels[seqCode]
        fout.write("%-4s%-3s" % (str(seqCode),resLabel))
  
        for atomName in self.headerAtomNames:
        
          if self.chemShiftsBySeqCodeAndAtomName.has_key(seqCode) and \
             self.chemShiftsBySeqCodeAndAtomName[seqCode].has_key(atomName):
            
            value = self.chemShiftsBySeqCodeAndAtomName[seqCode][atomName].value
          
          else:
            
            value = 0.0
          
          fout.write(" %8.2f " % (value))
          
        fout.write(self.newline)
    
    else:
      
      print "  Error: no sequence codes set for %s export. Aborting." % self.format

class CsiChemShift:

  def __init__(self,parent,value,atomName,seqCode,resLabel,molCode):
  
    self.value = value
    self.atomName = atomName
    #self.atomType = atomType
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.molCode = molCode
    self.resLabel = resLabel
    
    if not parent.chemShiftsBySeqCodeAndAtomName.has_key(self.seqCode):
      parent.chemShiftsBySeqCodeAndAtomName[self.seqCode] = {}
    
    parent.chemShiftsBySeqCodeAndAtomName[self.seqCode][self.atomName] = self
