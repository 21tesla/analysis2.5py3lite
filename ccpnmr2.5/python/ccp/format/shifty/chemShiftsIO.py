"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for Shifty chemical shift files

Copyright (C) 2011-2012 Wim Vranken (Vrije Universiteit Brussel)

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
from ccp.format.shifty.generalIO import ShiftyGenericFile
from ccp.format.shiftx.chemShiftsIO import ShiftxChemShiftFile, ShiftxChemShift

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class ShiftyChemShift(ShiftxChemShift):

  def checkAtomName(self,atomName):

    # Need to handle GLY differently...
    if self.resLabel == 'G' and atomName == 'HA':
      atomName = "QA"
    elif atomName == 'CO':
      atomName = 'C'
    elif atomName == 'HN':
      atomName = 'H'
      
    return atomName
  
class ShiftyChemShiftFile(ShiftyGenericFile,ShiftxChemShiftFile):

  shiftClass = ShiftyChemShift
    
  def read(self, verbose=False):
  
    if verbose:
      print "Reading %s chemical shift list %s" % (self.format,self.name)

    fin = open(self.name, 'rU')
    lines = fin.readlines()
    fin.close()

    headerCols = []
    
    for line in lines:

      if self.patt['emptyline'].search(line):
        continue
              
      #
      # Get the info... 
      #
      
      cols = line.split()
      
      if cols[0].count('NUM'):
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
            value = returnFloat(cols[i],default=None,verbose=False)
            
            # Ignore if value is 0.00. This is stupid but has to be done because format is set up this way.
            if value != 0.00:

              self.chemShifts.append(self.shiftClass(self,value,atomName,seqCode,resLabel,self.defaultMolCode))

    fin.close()

  def write(self,verbose=False):

    if verbose:
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

      # Write only backbone specific info
      fout.write("#NUM AA   HA     CA      CB       CO        N       HN    " + self.newline)
      #fout.write("#NUM AA   HA     H       HN        CA      CB       CO" + self.newline)

      for seqCode in seqCodeRange:

        if not self.seqCodeLabels.has_key(seqCode):
          print "  Warning: no %s output for sequence code %d - no information available." % (self.format,seqCode)
          continue

        resLabel = self.seqCodeLabels[seqCode]
        fout.write(" %-6s%-2s" % (str(seqCode),resLabel))

        #for atomName in ("HA","H","HN","CA","CB","CO"):
        for atomName in ("HA","CA","CB","CO","N","HN"):

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
          
    else:
      
      print "  Error: no sequence codes set for %s export. Aborting." % self.format
