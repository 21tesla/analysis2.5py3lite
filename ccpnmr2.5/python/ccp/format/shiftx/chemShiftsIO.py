"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for Shiftx chemical shift files

Copyright (C) 2009 Wim Vranken (European Bioinformatics Institute)

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
from ccp.format.shiftx.generalIO import ShiftxGenericFile

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################


class ShiftxChemShift:

  def __init__(self,parent,value,atomName,seqCode,resLabel,molCode):
  
    self.resLabel = resLabel
    self.value = value
    self.atomNameOrig = atomName
    self.atomName = self.checkAtomName(atomName)
    #self.atomType = atomType
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.molCode = molCode
    
    if not parent.chemShiftsBySeqCodeAndAtomName.has_key(self.seqCode):
      parent.chemShiftsBySeqCodeAndAtomName[self.seqCode] = {}
    
    if not parent.chemShiftsBySeqCodeAndAtomName[self.seqCode].has_key(self.atomNameOrig):
      parent.chemShiftsBySeqCodeAndAtomName[self.seqCode][self.atomNameOrig] = self

  def checkAtomName(self,atomName):

    # Need to handle GLY differently...
    if self.resLabel == 'G' and atomName == 'HA':
      atomName = "QA"
    
    return atomName

class ShiftxChemShiftFile(ShiftxGenericFile):
  """
  Information on file level
  """
  
  shiftClass = ShiftxChemShift
  
  def initialize(self):
  
    self.chemShifts = []
    
    self.chemShiftsBySeqCodeAndAtomName = {}
    self.seqCodeLabels = {}
    
  def read(self,isOutputFile = False, headerCols = None, verbose = 0):
  
    if verbose == 1:
      print "Reading %s chemical shift list %s" % (self.format,self.name)

    fin = open(self.name, 'rU')

    headerCols = []
    colLen = 0
    
    line = fin.readline()
    
    while line:

      if self.patt['emptyline'].search(line):

        line = fin.readline()
        continue
              
      #
      # Get the info... 
      #
      
      cols = line.split()
      
      if cols[0].upper().count('NUM'):
        headerCols = line.split()
        colLen = len(headerCols)
        
      elif len(cols) == colLen and not cols[0].count("-"):
        
        badShifts = False
        
        # These are fishy shifts, need to remove the *!
        if cols[0][0] == '*':
          cols[0] = cols[0][1:]
          badShifts = True
        
        seqCode = returnInt(cols[0],default=None,verbose=False)
        resLabel = cols[1]
        
        if seqCode != None:
          self.seqCodeLabels[seqCode] = resLabel
          
          for i in range(2,colLen):
            atomName = headerCols[i]
            value = returnFloat(cols[i],default=0.00,verbose=False)
            
            # Ignore if value is 0.00. This is not good but has to be done because format is set up this way.
            if value != 0.00:

              self.chemShifts.append(self.shiftClass(self,value,atomName,seqCode,resLabel,self.defaultMolCode))

      line = fin.readline()

    fin.close()

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing %s chemical shift list %s" % (self.format,self.name)


    fout = open(self.name,'w')

    #
    # Write out chem shifts
    #

    # This should be ints
    seqCodes = self.chemShiftsBySeqCodeAndAtomName.keys()
    seqCodes.sort()

    if seqCodes:

      seqCodeRange = range(seqCodes[0],seqCodes[-1]+1)

      # Write backbone specific info
      fout.write("  NUM RES   HA     H       N        CA      CB       C" + self.newline)
      fout.write("--- --- ------ ------ -------- ------- ------- --------" + self.newline)

      for seqCode in seqCodeRange:

        if not self.seqCodeLabels.has_key(seqCode):
          print "  Warning: no Shiftx output for sequence code %d - no information available." % seqCode
          continue

        resLabel = self.seqCodeLabels[seqCode]
        fout.write(" %-6s%-2s" % (str(seqCode),resLabel))

        for atomName in ("HA","H","N","CA","CB","C"):

          #if resLabel == 'G' and atomName == 'HA':
          #  searchAtomName = atomName + '*'
          #else:
          #  searchAtomName = atomName

          if self.chemShiftsBySeqCodeAndAtomName.has_key(seqCode) and \
             self.chemShiftsBySeqCodeAndAtomName[seqCode].has_key(atomName):

            value = self.chemShiftsBySeqCodeAndAtomName[seqCode][atomName].value

          else:

            value = 0.0

          fout.write(" %6.4f" % (value))

        fout.write(self.newline)
      
      fout.write(self.newline)
      
      # Write out all atom info
      fout.write(" NUM RES  ")
      headerAtomNames = ("H","HA","HB","HB2","HB3","HD1","HD2","HD21","HD22","HD3","HE","HE1","HE2","HE21","HE22","HE3","HG","HG1","HG12","HG13","HG2","HG3","HZ")
      for headerAtomName in headerAtomNames:
        fout.write("%-5s" % headerAtomName)
      fout.write(self.newline)
      
      
      for seqCode in seqCodeRange:
      
        if not self.seqCodeLabels.has_key(seqCode):
          continue
          
        resLabel = self.seqCodeLabels[seqCode]
        fout.write(" %-6s%-2s" % (str(seqCode),resLabel))

        for atomName in headerAtomNames:

          #if resLabel == 'G' and atomName == 'HA':
          #  searchAtomName = atomName + '*'
          #else:
          #  searchAtomName = atomName

          if self.chemShiftsBySeqCodeAndAtomName.has_key(seqCode) and \
             self.chemShiftsBySeqCodeAndAtomName[seqCode].has_key(atomName):

            value = self.chemShiftsBySeqCodeAndAtomName[seqCode][atomName].value

          else:

            value = 0.0

          fout.write(" %4.2f" % (value))

        fout.write(self.newline)
        
    
    else:
      
      print "  Error: no sequence codes set for %s export. Aborting." % self.format
