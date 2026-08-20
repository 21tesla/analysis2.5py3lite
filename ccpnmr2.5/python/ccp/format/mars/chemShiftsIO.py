"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for Mars chemical shift files

Copyright (C) 2007 Wim Vranken (European Bioinformatics Institute)

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
from memops.universal.Util import returnFloat
from memops.universal.Util import returnInt
from ccp.format.mars.generalIO import MarsGenericFile

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class MarsChemShiftFile(MarsGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.chemShifts = []
    self.spinSystemLinks = {}
    
    self.backboneAssignNames = ('H', 'N', 'CO', 'CO-1', 'CA', 'CA-1', 'CB', 'CB-1', 'HA', 'HA-1')
  
  def findChemShift(self,spinSystemId,atomName):
  
    csFound = None

    # IUPAC name for CO is C, need to compensate for this if written from MarsFormat.py
    if atomName == 'CO':
      searchAtomNames = ('C', atomName)
    else:
      searchAtomNames = (atomName,)
  
    for chemShift in self.chemShifts:
      if chemShift.spinSystemId == spinSystemId and chemShift.atomName in searchAtomNames:
        csFound = chemShift

    return csFound
        
  def read(self, headerCols = None, verbose = 0):
    
    if verbose == 1:
      print "Reading %s chemical shift list %s" % (self.format,self.name)

    fin = open(self.name, 'rU')

    line = fin.readline()
    
    spinSystemId = 1

    while line:

      if self.patt['%sAtomNameHeader' % self.format].search(line):

        if not headerCols:

          #
          # Get atom header info...
          #
          
          headerCols = line.split()

        line = fin.readline()
        continue

      if self.patt['emptyline'].search(line):
        line = fin.readline()
        continue
      
      #
      # Make sure header info is available - otherwise no point
      #
      
      if not headerCols:
        raise self.FileParseError("No header column information available. Cannot read file!")
        return
      
      #
      # Get the info... this is all Sparky-ish
      #
      
      cols = line.split()
      
      curSsLabel = cols.pop(0)
      
      shiftValueStrings = cols     
      
      for i in range(0,len(shiftValueStrings)):
      
        atomId = headerCols[i]
        value = shiftValueStrings[i]
        
        if value == '-':
          continue
          
        value = returnFloat(value)
        
        atomSearch = self.patt['%sAtomInfo' % self.format].search(atomId)
        
        atomName = atomSearch.group(1)
        atomPlace = atomSearch.group(2)
        
        if atomPlace and atomPlace == '-1':
          curSpinSystemId = spinSystemId - 1
          prevSpinSystemId = None
        
        else:        
          curSpinSystemId = spinSystemId
          prevSpinSystemId = spinSystemId - 1
          
        self.chemShifts.append(MarsChemShift(value,atomName,curSpinSystemId,curSsLabel,self.defaultMolCode, prevSpinSystemId = prevSpinSystemId))

      self.spinSystemLinks[spinSystemId] = {-1: (spinSystemId - 1)}
      
      spinSystemId += 2

      line = fin.readline()

    fin.close()

  def write(self,verbose = 0):
    
    if verbose == 1:
      print "Writing %s chemical shift list %s" % (self.format,self.name)

    fout = open(self.name,'w')
    
    #
    # Determine relevant header items and organise for output
    #
    
    spinSystemAtomList = {}
    spinSystemLabels = {}
    
    atomNamesFound = []
    
    spinSystemIds = self.spinSystemLinks.keys()
    
    for spinSystemId in spinSystemIds:
      
      spinSystemAtomList[spinSystemId] = []
      
      prevSpinSystemId = None
      
      if self.spinSystemLinks[spinSystemId].has_key(-1) and self.spinSystemLinks[spinSystemId][-1] != None:
        prevSpinSystemId = self.spinSystemLinks[spinSystemId][-1]
      for atomName in self.backboneAssignNames:
        chemShift = None
        if atomName[-2:] == '-1':
          if prevSpinSystemId != None:
            chemShift = self.findChemShift(prevSpinSystemId,atomName[:-2])
        else:
          chemShift = self.findChemShift(spinSystemId,atomName)
          
          if chemShift and not spinSystemLabels.has_key(spinSystemId):
            if chemShift.ssLabel:
              spinSystemLabels[spinSystemId] = chemShift.ssLabel
            else:
              spinSystemLabels[spinSystemId] = "%d" % spinSystemId
              
        spinSystemAtomList[spinSystemId].append(chemShift)
        
        if chemShift and atomName not in atomNamesFound:
          atomNamesFound.append(atomName)
          
    #
    # Only print out relevant info
    #

    headerColumns = []
    
    for atomName in self.backboneAssignNames:
      if atomName in atomNamesFound:
        headerColumns.append(atomName)
      else:
        headerColumns.append(None)
        
    #
    # Print out header
    #
    
    fout.write("         ")
    for header in headerColumns:
      if header:
        fout.write("   %-4s   " % header)

    fout.write(self.newline)

    #
    # Write out chem shifts
    #
    
    spinSystemIds.sort()
    
    for spinSystemId in spinSystemIds:
      # Only write out relevant columns
      fout.write("%-8s " % spinSystemLabels[spinSystemId])
      for i in range(len(headerColumns)):
        if headerColumns[i]:
          chemShift = spinSystemAtomList[spinSystemId][i]
          
          if chemShift:
            fout.write("%8.3f  " % chemShift.value)
          else:
            fout.write("    -     ")
      
      fout.write(self.newline)

class MarsChemShift:

  def __init__(self,value,atomName,spinSystemId,ssLabel,defaultMolCode, prevSpinSystemId = None):
  
    self.value = returnFloat(value)
    self.atomName = atomName
    self.molCode = defaultMolCode
    self.ssLabel = ssLabel
    
    self.spinSystemId = spinSystemId
    self.prevSpinSystemId = prevSpinSystemId
    
    self.seqCode = None      # Dummy value
    self.seqInsertCode = ''  # Dummy value
