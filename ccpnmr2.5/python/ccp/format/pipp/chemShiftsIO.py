"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for PIPP chemical shift files

Copyright (C) 2005-2007 Wim Vranken (European Bioinformatics Institute)

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
from memops.universal.Util import returnFloat, returnInt
from ccp.format.pipp.generalIO import PippGenericFile

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class PippChemShiftFile(PippGenericFile):

  """
  Information on file level
  """
  
  def initialize(self):
  
    self.chemShifts = []
    self.experiments = []
    
    self.vars = {}

  def readLine(self,fin):
  
    line = fin.readline()
    
    while (1):

      if len(line.split()) == 0 or self.patt['hash'].search(line):
        
        line = fin.readline()
        continue
      
      break
    
    return line
    
  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s chemical shift list %s" % (self.format,self.name)

    fin = open(self.name, 'rU')
    
    experimentRead = 0
    shiftRead = 0

    #
    # Start reading, look for header line
    #

    line = fin.readline()

    while line:
    
      cols = line.split()

      if len(cols) == 0 or self.patt['hash'].search(line):
        
        pass

      else:
      
        if cols[0] in ['SHIFT_FL_FRMT','FIRST_RES_IN_SEQ','RES_ID','SPIN_SYSTEM_ID','RES_TYPE','PREV_TYPE','HETEROGENEITY']:
          
          self.vars[cols[0]] = cols[1]

          if cols[0] == 'SPIN_SYSTEM_ID':
            
            experimentRead = 0
            shiftRead = 1
            currentAtomNames = []
        
        elif cols[0] == 'EXP_PEAK_PICK_TBL':
          
          line = self.readLine(fin)
          labels = line.split('  ')
          
          labelInd = 0
          
          while(labelInd < len(labels)):
            if labels[labelInd] == '':
              labels.pop(labelInd)
            else:
              labels[labelInd] == labels[labelInd].strip()
              labelInd += 1
          
          experimentRead = 1
        
        elif cols[0] == 'END_RES_DEF':
        
          shiftRead = 0
        
        elif shiftRead:
        
          searchObj = self.patt[self.format + 'Shift'].search(line)
          
          if not searchObj:
            searchObj = self.patt[self.format + 'ShiftNoAss'].search(line)
            peakAssign = None
          else:
            peakAssign = searchObj.group(3)
          
          if searchObj and self.vars['RES_ID'] != '*':
          
            atomNameMark = None
            atomName = searchObj.group(1)
          
            if atomName.count(':'):
              (atomName,atomNameMark) = atomName.strip().split(':')
            
            if atomNameMark != '-1':
              chemShift = PippChemShift(self)
              chemShift.initialiseFromParent(atomName,searchObj.group(2),peakAssign,currentAtomNames)
              self.chemShifts.append(chemShift)
        
        elif experimentRead:
          
          self.experiments.append({})

          for i in range(0,len(labels)):
            # Possible that end ones (?) are not present...
            if i < len(cols):
              self.experiments[-1][labels[i]] = cols[i]
            else:
              self.experiments[-1][labels[i]] = None

      line = fin.readline()

    fin.close()

  def write(self,verbose = 0):
    
    if verbose == 1:
      print "Writing %s chemical shift list %s" % (self.format,self.name)

    fout = open(self.name,'w')
    
    shiftsByResidue = {}
    residueTypes = {}
    for chemShift in self.chemShifts:
      if not shiftsByResidue.has_key(chemShift.seqCode):
        shiftsByResidue[chemShift.seqCode] = {}
        residueTypes[chemShift.seqCode] = chemShift.resLabel
      atomNameKey = (chemShift.atomName,) + tuple(chemShift.otherAtomNames)
      shiftsByResidue[chemShift.seqCode][atomNameKey] = chemShift
      
    
    seqCodes = shiftsByResidue.keys()
    seqCodes.sort()
    
    for ssIndex in range(len(seqCodes)):
      seqCode = seqCodes[ssIndex]
      
      fout.write("RES_ID\t\t%d" % seqCode + self.newline)
      fout.write("RES_TYPE\t%s" % residueTypes[seqCode] + self.newline)
      fout.write("SPIN_SYSTEM_ID\t%d" % (ssIndex + 1) + self.newline)
   
      #
      # Write out chem shifts
      #
      
      atomNameKeys = shiftsByResidue[seqCode].keys()
      atomNameKeys.sort()
      
      for atomNameKey in atomNameKeys:
        
        atomNameText = string.join(atomNameKey,'|')
        chemShift = shiftsByResidue[seqCode][atomNameKey]

        fout.write("    %-7s  %8.4f" % (atomNameText,chemShift.value))
        fout.write(self.newline)
      
      fout.write("END_RES_DEF" + self.newline)
      fout.write(self.newline)

class PippChemShift:

  def __init__(self,parent):

    self.parent = parent
    
  def setValue(self,attrName,attrValue):
    
    setattr(self,attrName,attrValue)
  
  def initialiseFromParent(self,atomName,value,peakAss,currentAtomNames):
    
    self.resLabel = self.parent.vars['RES_TYPE']
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(returnInt(self.parent.vars['RES_ID']))
    
    #
    # Handle atomNames
    #
    
    atomNames = atomName.split('|')
    for i in range(0,len(atomNames)):
      atomNames[i] = atomNames[i].strip()
      
    self.atomName = None
    
    for atomName in atomNames:
      if atomName not in currentAtomNames:
        self.atomName = atomName
        atomNames.pop(atomNames.index(atomName))
        self.otherAtomNames = atomNames
        currentAtomNames.append(atomName)
        break
    
    if not self.atomName:
      self.atomName = atomNames[0].strip()
      self.otherAtomNames = atomNames[1:]
      
    #
    # Set ambiguity if possible
    #
    
    self.ambCode = None
    
    if not self.otherAtomNames:
      self.ambCode = 1
    
    else:
      nonStereo = 0
      for otherAtomName in self.otherAtomNames:
        if otherAtomName[:-1] == self.atomName[:-1]:
          nonStereo += 1
      if nonStereo == len(self.otherAtomNames):
        self.ambCode = 2
    
    #
    # Handle shift values
    #
    
    self.value = returnFloat(value)
    self.valueError = 0.0
    
    #
    # Handle assignments
    #
    
    self.peakAssignments = []
    self.numAssignments = 0
    
    if peakAss:
      self.peakAssignments = peakAss.split()
      self.numAssignments = len(self.peakAssignments)
    
    self.molCode = self.parent.defaultMolCode
