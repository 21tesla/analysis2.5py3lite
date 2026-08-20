"""
======================COPYRIGHT/LICENSE START==========================

generalIO.py: General I/O information for AMBER files

Copyright (C) 2010 Wim Vranken (European Bioinformatics Institute)

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
from memops.universal.Util import returnInt

from ccp.format.general.formatIO import FormatFile

from ccp.format.general.Constants import defaultMolCode
from ccp.format.general.Util import getSeqAndInsertCode

from memops.universal.Util import returnFloat
from memops.universal.Util import returnInt

#####################
# Class definitions #
#####################

class AmberGenericFile(FormatFile):

  def setGeneric(self):
    
    self.format = 'amber'
    self.defaultMolCode = defaultMolCode
    
  def readConstraints(self):
       
    #
    # Start reading...
    #

    fin = open(self.name, 'rU')
    lines = fin.readlines()
    fin.close()
    
    #
    # Parse data...
    #
          
    constraintNum = 0
    constraintLines = []
    
    for line in lines:

      if self.patt['emptyline'].search(line) or self.patt['hash'].search(line):
        continue
      
      newConstraint = False
      
      restraintStartSearch = self.patt[self.format + 'RestraintStart'].search(line)
      if restraintStartSearch:
        line = line.replace(restraintStartSearch.group(1),'')
        constraintNum += 1
        newConstraint = True
        
      restraintEndSearch = self.patt[self.format + 'RestraintEnd'].search(line)
      if restraintEndSearch:
        line = line.replace(restraintEndSearch.group(1),'')
              
      if newConstraint:
        constraintLines.append(line.strip())
      else:
        constraintLines[-1] += line.strip()
    
    #
    # Put info in dictionary
    #
    
    constraintInfoList = []

    for constraintLine in constraintLines:
      
      constraintCols = constraintLine.split(',')
      constraintContent = {}
      curConstraintKey = None

      for constraintCol in constraintCols:
      
        constraintCol = constraintCol.strip()
        
        if not constraintCol:
          continue
        
        if constraintCol.count('='):
          (constraintKey,constraintValue) = constraintCol.split('=')
          curConstraintKey = constraintKey.strip()
          constraintContent[curConstraintKey] = []
        else:
          constraintValue = constraintCol
        
        constraintValue = constraintValue.strip()
        if constraintValue.count('.'):
          value = returnFloat(constraintValue)
        else:
          value = returnInt(constraintValue)
          
        constraintContent[curConstraintKey].append(value)
      
      constraintInfoList.append(constraintContent)
    
    return constraintInfoList

class AmberGenericConstraint:

  def __init__(self,parent,Id):
    
    self.parent = parent
    self.Id = returnInt(Id)
    
    self.items = []

  def addItem(self,coordinateAtoms):
     
    item = AmberConstraintItem()
    self.items.append(item)
    
    for coordinateAtom in coordinateAtoms:
      
      item.members.append(AmberConstraintMember(coordinateAtom.chainId,coordinateAtom.resName,coordinateAtom.seqCode,coordinateAtom.atomName))
    
class AmberConstraintItem:

  def __init__(self):
    
    self.members = []
    
class AmberConstraintMember:

  def __init__(self,chainCode,resLabel,seqCode,atomName):
    
    self.chainCode = chainCode
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.atomName = atomName
    self.resLabel = resLabel
