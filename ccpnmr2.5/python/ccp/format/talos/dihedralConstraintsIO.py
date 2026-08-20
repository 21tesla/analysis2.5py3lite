"""
======================COPYRIGHT/LICENSE START==========================

dihedralConstraintsIO.py: I/O for Talos dihedral constraints (output) file

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

# Import general functions
from memops.universal.Util import returnFloat
from memops.universal.Util import returnInt
from ccp.format.talos.generalIO import TalosGenericFile, TalosGenericDataItem

from ccp.format.general.Util import getSeqAndInsertCode
from ccp.format.general.Constants import defaultMolCode

#####################
# Class definitions #
#####################
    
class TalosDihedralConstraintFile(TalosGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.constraints = []
    self.constraintElements = 4

    self.constraintNum = 1 
    
    #
    # Reference data setup
    #
    
    self.refData = [['RESID', '%4d','seqCode',returnInt,True],
                   	['RESNAME','%s','resLabel',str,True], 
                   	['PHI','%8.3f','value',returnFloat,True],
                   	['PSI','%8.3f','value',returnFloat,True],
                   	['DPHI','%8.3f','error',returnFloat,True],
                   	['DPSI','%8.3f','error',returnFloat,True],
                   	['DIST','%8.3f','dist',returnFloat,True],
                   	['S2','%5.3f','s2',returnFloat,False],
                   	['COUNT','%2d','count',returnInt,True],
                    ['CS_COUNT','%2d','chemShiftCount',returnInt,False],
                   	['CLASS','%s','className',str,True],
                   ]
    

  def printInfo(self,action):
    
    print "%s %s dihedral constraint list %s" % (action,self.format,self.name)

  def setVarsLine(self,varsDict):
    
    for angleName in ['PSI','PHI']:
      
      refData = []
      
      for refDataItem in self.refData:
        if refDataItem[0] in ['RESID','RESNAME','DIST','COUNT','CLASS']:
          refData.append(refDataItem)
        elif refDataItem[0][-3:] == angleName:
          refData.append(refDataItem)
 
      if varsDict[angleName] != '9999.000':
        self.constraints.append(TalosDihedralConstraint(self,angleName))
        self.constraints[-1].setVarsDict(varsDict,refData = refData)
        self.constraintNum += 1

  def getVarItems(self):
  
    self.varItems = self.constraints

class TalosDihedralConstraint(TalosGenericDataItem):
 
  def __init__(self,parent,name):
  
    self.parent = parent
    self.items = []
    self.name = name

  def setItemSpecificVars(self):
  
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(self.seqCode)
    self.Id = self.parent.constraintNum
    
    self.targetAngle = self.value
    self.devAngle = self.error
    
    self.setAtomMembers()
    
  def setAtomMembers(self):
    
    self.items.append(TalosConstraintItem())
    
    if self.name == 'PHI':
      
      refAtoms = (('C',-1),('N',0),('CA',0),('C',0))
      
    elif self.name == 'PSI':
      
      refAtoms = (('N',0),('CA',0),('C',0),('N',1))
    
    for (atomName,location) in refAtoms:
            
      self.items[-1].members.append(TalosConstraintMember(self.seqCode + location,self.seqInsertCode,atomName))
      
class TalosConstraintItem:

  def __init__(self):
    
    self.members = []
    
class TalosConstraintMember:

  def __init__(self,seqCode,seqInsertCode,atomName):
    
    self.chainCode = defaultMolCode
    self.seqCode = seqCode
    self.seqInsertCode = seqInsertCode
    self.atomName = atomName
