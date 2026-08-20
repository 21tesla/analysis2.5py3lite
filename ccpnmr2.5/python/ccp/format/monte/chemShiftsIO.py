#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for Monte chemical shift files

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
from ccp.format.monte.generalIO import MonteGenericFile

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class MonteChemShiftFile(MonteGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.chemShifts = []
    
  def read(self,isOutputFile = False, headerCols = None, verbose = 0):
    
    #
    # TODO TODO also need a 'readFinal' one to read the FINAL information!!
    #           set a flag in MonteFormat.py to select which cs info to read...
    
    """
    The first column in any data entry is a peak identification number. This is a four
    digit number with one decimal place. The last digit may be incremented, for example,
    to indicated a single parent amide resonance that may give rise to more than one
    carbon or proton cross peaks, i.e. from minor protein conformations or overlapping
    resonance signals.
    
    The second entry is an optional comment of up to eight characters. If there is a
    number in this field, the software will assume that this is a tentative assignment.
    This number is not used in the Monte Carlo protocol; however, is carried over to the
    output files. Using this field to add/change tentative assignment is useful in
    comparing output from multiple Monte runs.
    
    Subsequent entries in each row are chemical shift frequencies and may be separated by
    spaces or comas. The first and second columns must contain the parent amide nitrogen
    and amide proton shifts, in either order (defined by the .par file). The following
    chemical shifts may be in any order, but must be consistent throughout the table. The
    order is user defined in the .par file. (see section on Parameter File)
    
    Any chemical shifts at the very end of the list may be left blank if the shift is
    unknown. However, a unknown chemical shift field positioned within a series of known
    shifts must be entered 0.00 to hold its place (see example peak 101.0 above, no shift
    is identified with the intra residue CA).    
    """

    if verbose == 1:
      print "Reading %s chemical shift list %s" % (self.format,self.name)

    fin = open(self.name, 'rU')

    line = fin.readline()
    
    spinSystemId = 0
    resLabel = oldResLabel = None

    while line:

      if self.patt['%sComment' % self.format].search(line):

        if not isOutputFile and not self.chemShifts and not headerCols:

          #
          # Get atom info from first line...
          #
          
          headerCols = line.split()
          headerCols.pop(0)

        line = fin.readline()
        continue

      if self.patt['emptyline'].search(line):
        line = fin.readline()
        continue
        
      #
      # Make sure header info is available - otherwise no point
      #
      
      if not headerCols:
        raise "Error: no header column information available. Try reading .par file!"
        return
      
      #
      # Get the info... should really come for .par file!!
      #
      
      cols = line.split()
      
      infoCode = None
      
      if not isOutputFile:
      
        stripId = returnFloat(cols.pop(0))

        #
        # NOt necessarily info string available...
        #

        if self.patt['onlyFloat'].search(cols[0]):
          seqCode = None
          resLabel = None

        else:
          assignment = cols.pop(0)

          searchAssignment = self.patt['%sAssignment' % self.format].search(assignment)

          resLabel = searchAssignment.group(1)
          seqCode = searchAssignment.group(2)
      
      else:
      
        seqCode = cols.pop(0)
        if seqCode[-1] in '+':
          seqCode = seqCode[:-1]
          infoCode = seqCode[-1]
          
        oldResLabel = resLabel
        resLabel = cols.pop(0)
        stripId = returnFloat(cols.pop(0))
        voidCol = cols.pop(0)
      
      #
      # Set up info for atoms...
      #
      
      if not seqCode or seqCode == '?':
        seqCode = None
        spinSystemId = spinSystemId + 2
      else:
        seqCode = returnInt(seqCode)

      if len(cols) == 1:
        cols = cols.split(',')

      values = returnFloats(cols)

      for i in range(0,len(values)):
        atomId = headerCols[i]
        value = values[i]
        
        if value == 0.0:
          continue
        
        atomSearch = self.patt['%sAtomInfo' % self.format].search(atomId)
        
        atomName = atomSearch.group(1)
        atomPlace = atomSearch.group(2)
        
        if atomName == 'HA1':
          nextAtomValue = values[i+1]
          if nextAtomValue == 0.00:
            atomName = 'HA'
        
        curSeqCode = seqCode
        curResLabel = None
        
        if seqCode == None:
          curSpinSystemId = spinSystemId
          prevSpinSystemId = spinSystemId - 1
        else:
          curSpinSystemId = None
          prevSpinSystemId = None
        
        if atomPlace == '(i-1)' or atomPlace == '-1':

          if seqCode != None:
            curSeqCode = seqCode - 1
          else:
            curSpinSystemId = spinSystemId - 1
            prevSpinSystemId = None
          
          if not isOutputFile:
            curResLabel = resLabel
          else:
            curResLabel = oldResLabel
        
        elif isOutputFile:
          curResLabel = resLabel

        self.chemShifts.append(MonteChemShift(value,atomName,curSeqCode,curSpinSystemId,stripId,curResLabel,self.defaultMolCode, infoCode = infoCode, prevSpinSystemId = prevSpinSystemId))

      line = fin.readline()

    fin.close()

  def write(self,verbose = 0):
    
    #
    # TODO TODO set up for writing spin system info for input!!
    #
    
    print "NOT WORKING YET"
    return

    if verbose == 1:
      print "Writing %s chemical shift list %s" % (self.format,self.name)


    fout = open(self.name,'w')

    #
    # Write out chem shifts
    #

    for chemShift in self.chemShifts:

      fout.write("%3d%s%-5s %8.3f %d" % (chemShift.seqCode,
                    	 self.seqNameSep,
                    	 chemShift.atomName,
                    	 chemShift.value,
                    	 chemShift.unknownCode))
      fout.write(self.newline)

class MonteChemShift:

  def __init__(self,value,atomName,seqCode,spinSystemId,stripId,resLabel,defaultMolCode, infoCode = None, prevSpinSystemId = None):
  
    self.value = returnFloat(value)
    self.atomName = atomName
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.molCode = defaultMolCode
    self.resLabel = resLabel
    
    self.spinSystemId = spinSystemId
    self.prevSpinSystemId = prevSpinSystemId
    self.stripId = stripId
    self.infoCode = infoCode
