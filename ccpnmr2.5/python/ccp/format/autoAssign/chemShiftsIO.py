#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for AutoAssign chemical shift files

Copyright (C) 2006 Wim Vranken (European Bioinformatics Institute)

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
from ccp.format.autoAssign.generalIO import AutoAssignGenericFile

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class AutoAssignChemShiftFile(AutoAssignGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.chemShifts = []
    self.seqCodes = []
    
  def read(self, verbose = 0):

    if verbose == 1:
      print "Reading %s chemical shift list %s" % (self.format,self.name)

    fin = open(self.name, 'rU')
    
    atomCols = []

    line = fin.readline()
    
    while line:
      
      cols = line.split()
     
      if cols:
     
        seqCode1Or3LetterAndCodeSearch = self.patt[self.format + 'seqCode1Or3LetterAndCode'].search(cols[0])
        
        if seqCode1Or3LetterAndCodeSearch:
          
          resCode = seqCode1Or3LetterAndCodeSearch.group(1)
          #if len(resCode) == 1:
          #  codeType = 'code1Letter'
          #else:
          #  codeType = 'code3Letter'
          
          seqCode = returnInt(seqCode1Or3LetterAndCodeSearch.group(2))
          
          if not self.seqCodes or self.seqCodes[-1] != (seqCode,resCode):
            self.seqCodes.append((seqCode,resCode))
  
          spinSystemInfoSearch = self.patt[self.format + 'spinSystemInfo'].search(line)
          
          if spinSystemInfoSearch:
            spinSystemInfo = spinSystemInfoSearch.group(1)
            line = line.replace(spinSystemInfo,"")
            cols = line.split()
            (spinSystemId,peakId) = spinSystemInfo[1:-1].split()
          else:
            spinSystemId = peakId = None
           
         
          for colNum in range(1,len(cols)):
  
            value = returnFloat(cols[colNum])
            curSeqCode = seqCode
            curResCode = resCode
            iCode = 'i'
            
            lineIndex = line.find(cols[colNum])
            line = line.replace(cols[colNum]," " * len(cols[colNum]),1)
            
            atomName = None
            
            for (tAtomName,atomLineIndex) in atomCols:
              for i in range(0,3):
                if lineIndex == atomLineIndex - i:
                  atomName = tAtomName
                  break                
            
            if atomName.count('-1') > 0:
              # Check if it's already there...
              
              atomName = atomName[:-2]
              for i in range(len(self.chemShifts)-1,-1,-1):
                chemShift = self.chemShifts[i]
                if chemShift.seqCode == (seqCode - 1) and chemShift.atomName == atomName:
                  #print "Resetting", seqCode -1, atomName,value
                  atomName = None
                  chemShift.value =  (chemShift.value + value) / 2
                  chemShift.allValues['i+1'] = value
                  chemShift.valueError = abs(chemShift.value - chemShift.allValues['i+1'])
              
              if not atomName:
                continue
              
              (curSeqCode,curResCode) = self.seqCodes[-2]
              iCode = 'i+1'
                  
            self.chemShifts.append(AutoAssignChemShift(value,atomName,curSeqCode,spinSystemId,peakId,curResCode,self.defaultMolCode,iCode = iCode))
            #print "set",curSeqCode,curResCode,atomName,value
            
        elif cols[0] == 'AA':
         
          for atomName in cols[1:]:
            useAtomName = atomName
            if atomName == 'N15':
              useAtomName = 'N'
            atomCols.append((useAtomName,line.find(atomName)))
            line = line.replace(atomName," " * len(atomName),1)
 
      line = fin.readline()

    fin.close()

  def write(self,verbose = 0):
    
    print "Not relevant"

class AutoAssignChemShift:

  def __init__(self,value,atomName,seqCode,spinSystemId,peakId,resLabel,defaultMolCode,iCode = 'i'):
  
    self.value = returnFloat(value)
    self.valueError = 0.0
    self.atomName = atomName
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.molCode = defaultMolCode
    self.resLabel = resLabel
    
    self.allValues = {iCode: self.value}
    
    self.spinSystemId = spinSystemId
    self.peakId = peakId
