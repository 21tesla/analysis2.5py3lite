#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

rdcConstraintsIO.py: I/O for nmrStar rdc constraints saveframes

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

import os

from memops.universal.Io import getTopDirectory

from ccp.format.nmrStar.generalIO import NmrStarConstraintFile
from ccp.format.nmrStar.generalIO import NmrStarFile
from ccp.format.nmrStar.generalIO import GenericConstraint

#
# For DIANA torsion angle name information
#

from ccp.format.cyana.cyanaLibParser import CyanaLibrary

from ccp.format.general.Constants import defaultMolCode
from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################


class NmrStarFile(NmrStarFile):

  def initialize(self, version = '2.1.1'):
  
    self.constraintFiles = []
    self.files = self.constraintFiles

    if not self.version:
      self.version = version
    
    if self.version == '3.1':
      self.saveFrameName = 'RDC_constraints'    
    else:
      self.saveFrameName = 'residual_dipolar_couplings'

    self.DataClassFile = NmrStarRdcConstraintFile
    
    self.components = [self.saveFrameName]
    self.setComponents()

  def read(self,verbose = 0):
  
    self.readComponent(verbose = verbose) 

class NmrStarRdcConstraintFile(NmrStarConstraintFile):

  def initialize(self,parent,saveFrame = None):

    self.constraints = []
    
    self.constraintElements = 2
    
    self.cyanaLibUsed = 0
    self.cyanaLib = CyanaLibrary()
    
    self.saveFrame = saveFrame
         
    self.parent = parent
    self.version = parent.version

    if self.saveFrame:
      self.parseSaveFrame()

  def parseSaveFrame(self):
    
    if not self.checkVersion():
      return

    #
    # Set up reference info
    #
    
    tableTagNames = {}
    
    if self.version == '2.1.1':

      tableName = '_Residual_dipolar_coupling_ID'
      tableTagNames['ID'] = tableName

      tableTagNames['atomPresent'] = "_Atom_two_atom_name"
      tableTagNames['position'] = ['one','two']
  
      tableTagNames['code'] = '_Residual_dipolar_coupling_code'
      tableTagNames['value'] = '_Residual_dipolar_coupling_value'
      
      #
      # Hack - both of these occur in the files for some reason
      #
      
      if self.saveFrame.tables.has_key(tableName):

        constraintTableTags = self.saveFrame.tables[tableName].tags

        # This value not necessarily available...
        if constraintTableTags.has_key('_Residual_dipolar_coupling_value_error'):
          tableTagNames['error'] = '_Residual_dipolar_coupling_value_error'
        elif constraintTableTags.has_key('_Residual_dipolar_coupling_error'):
          tableTagNames['error'] = '_Residual_dipolar_coupling_error'

      tableTagNames['lowerValue'] = '_Residual_dipolar_coupling_lower_bound_value'
      tableTagNames['upperValue'] = '_Residual_dipolar_coupling_upper_bound_value'
      
      # Same as in dihedralConstraintsIO!
      tableTagNames['chainCode'] = "_Atom_%s_mol_system_component_name"
      tableTagNames['seqCode'] = "_Atom_%s_residue_seq_code"
      tableTagNames['atomName'] = "_Atom_%s_atom_name"
      tableTagNames['resLabel'] = "_Atom_%s_residue_label"
      
    elif self.version[0] == '3':
    
      tableName = '_RDC_constraint'
      tableTagNames['ID'] = 'ID'

      tableTagNames['position'] = ['1','2']
  
      #tableTagNames['code'] = '_Residual_dipolar_coupling_code'
      tableTagNames['value'] = 'RDC_val'
      tableTagNames['error'] = 'RDC_val_err'
      tableTagNames['lowerValue'] = 'RDC_lower_bound'
      tableTagNames['upperValue'] = 'RDC_upper_bound'
 
      if self.version == '3.0':
        tableTagNames['atomPresent'] = "Label_atom_ID_2"
        tableTagNames['chainCode'] = "Label_entity_assembly_ID_%s"
        tableTagNames['seqCode'] = "Label_comp_index_ID_%s"
        tableTagNames['atomName'] = "Label_atom_ID_%s"
        tableTagNames['resLabel'] = "Label_comp_ID_%s"

      elif self.version == '3.1':
        tableTagNames['atomPresent'] = "Atom_ID_2"
        tableTagNames['chainCode'] =   "Entity_assembly_ID_%s"
        tableTagNames['seqCode'] =     "Comp_index_ID_%s"
        tableTagNames['atomName'] =    "Atom_ID_%s"
        tableTagNames['resLabel'] =    "Comp_ID_%s"
 
      if self.saveFrame.tables.has_key(tableName) and not self.saveFrame.tables[tableName].tags[tableTagNames['seqCode'] % 1][0]:
        tableTagNames['atomPresent'] = "Auth_atom_ID_2"
        
        if self.version == '3.0':
          tableTagNames['chainCode'] = "Auth_segment_code_%s"
        elif self.version == '3.1':
          tableTagNames['chainCode'] = "Auth_asym_ID_%s"
        
        tableTagNames['seqCode'] = "Auth_seq_ID_%s"
        tableTagNames['atomName'] = "Auth_atom_ID_%s"
        tableTagNames['resLabel'] = "Auth_comp_ID_%s"
    
    #
    # Could be missing...
    #

    if self.saveFrame.tables.has_key(tableName):
      
      constraintTableTags = self.saveFrame.tables[tableName].tags

      for i in range(0,len(constraintTableTags[tableTagNames['ID']])):

        rdcConstraint = NmrStarRdcConstraint()
        rdcConstraint.setData(constraintTableTags,i,tableTagNames)

        self.constraints.append(rdcConstraint)

        rdcConstraint.setErrors(self.saveFrame.tables[tableName].tagErrors,tableTagNames,i)

        if rdcConstraint.Id and not constraintTableTags[tableTagNames['atomPresent']][i]:

          #
          # Only one atom given - find other one... also set code to something for
          # output later on (hack for writing constraint nmrStar files for Jurgen)
          #

          rdcConstraint.setAtomMembers(constraintTableTags,i,self.cyanaLib,tableTagNames)
          rdcConstraint.code = True

          self.cyanaLibUsed = 1

        else:

          rdcConstraint.nodes.append(NmrStarRdcConstraintItem())

          for position in tableTagNames['position']:

            rdcConstraint.nodes[-1].members.append(NmrStarRdcConstraintMember())
            rdcConstraint.nodes[-1].members[-1].addStarInfo(constraintTableTags,i,position,tableTagNames)

    #
    # Comments
    #
        
    self.parseCommentsLoop()

class NmrStarRdcConstraint(GenericConstraint):

  def __init__(self,Id = None):
    
    self.Id = Id
    self.nodes = []   # THIS IS A HACK! Not really a node...
    self.errors = []
    self.code = None
  
  def setData(self,constraintLogicTableTags,i,tableTagNames):

    self.Id = constraintLogicTableTags[tableTagNames['ID']][i]

    self.value = constraintLogicTableTags[tableTagNames['value']][i]
    
    if tableTagNames.has_key('error'):
      self.error = constraintLogicTableTags[tableTagNames['error']][i]
    else:
      self.error = None
    
    if constraintLogicTableTags.has_key(tableTagNames['lowerValue']):
      self.lowerValue = constraintLogicTableTags[tableTagNames['lowerValue']][i]
      
    if constraintLogicTableTags.has_key(tableTagNames['upperValue']):
      self.upperValue = constraintLogicTableTags[tableTagNames['upperValue']][i]
    
    if tableTagNames.has_key('code') and constraintLogicTableTags.has_key(tableTagNames['code']):
      self.code = constraintLogicTableTags[tableTagNames['code']][i]
    
  def setAtomMembers(self,constraintTableTags,i,cyanaLib,tableTagNames):
   
    #
    # Need to set two members... for DIANA-CYANA only one available.
    # Using DIANA library to find out where the 'other' atom is...
    #
  
    chainCode = constraintTableTags[tableTagNames['chainCode'] % tableTagNames['position'][0]][i]

    if not chainCode:
      chainCode = defaultMolCode

    seqCode = constraintTableTags[tableTagNames['seqCode'] % tableTagNames['position'][0]][i]
    resLabel = constraintTableTags[tableTagNames['resLabel'] % tableTagNames['position'][0]][i]
    refAtomName = constraintTableTags[tableTagNames['atomName'] % tableTagNames['position'][0]][i]
    
    refAtom = cyanaLib.findAtom(resLabel,refAtomName)

    if not refAtom:
      print "  Error: atom %s (%s %s, chain '%s') not found in Cyana library!" % (refAtomName,resLabel,seqCode,chainCode)
    
    elif refAtom.bondedAtomSerials.count(0) != 3:
      print "  Error: invalid single atom %s (%s %s, chain '%s'). No or multiple bonded atoms" % (refAtomName,resLabel,seqCode,chainCode)
    
    else:
    
      self.nodes.append(NmrStarRdcConstraintItem())

      self.nodes[-1].members.append(NmrStarRdcConstraintMember(chainCode,seqCode,refAtomName,resLabel))

      atomSerial = refAtom.bondedAtomSerials[0]     
      bondedAtom = cyanaLib.findAtomBySerial(resLabel,atomSerial)
    
      self.nodes[-1].members.append(NmrStarRdcConstraintMember(chainCode,seqCode,bondedAtom.name,resLabel))

class NmrStarRdcConstraintItem:

  def __init__(self):
    
    self.Id = 1
    self.members = []    
    
class NmrStarRdcConstraintMember:

  def __init__(self,chainCode = None,seqCode = None, atomName = None, resLabel = None):
    
    self.chainCode = chainCode
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.atomName = atomName
    self.resLabel = resLabel
    
    self.defaultMolCode = defaultMolCode
    
  def addStarInfo(self,constraintTableTags,i,position,tableTagNames):

    # Tag can be missing in 2.1.1 files...
    if constraintTableTags.has_key(tableTagNames['chainCode'] % position):
      self.chainCode = constraintTableTags[tableTagNames['chainCode'] % position][i]
    else:
      self.chainCode = None

    if not self.chainCode:
      self.chainCode = self.defaultMolCode

    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(constraintTableTags[tableTagNames['seqCode'] % position][i])
    self.atomName = constraintTableTags[tableTagNames['atomName'] % position][i]
    self.resLabel = constraintTableTags[tableTagNames['resLabel'] % position][i]

###################
# Main of program #
###################

if __name__ == "__main__":  

  files = [#['../reference/ccpNmr/jurgenBmrb/1jwe/restraints.star','2.1.1'],
           #['../reference/ccpNmr/jurgenBmrb/1j7p/restraints.star','2.1.1']
           ['../reference/ccpNmr/jurgenBmrb/1b4c/1b4c.str','3.0']
          ]
  
  for (file,version) in files:
    
    file = os.path.join(getTopDirectory(), file)
    
    nmrStarFile = NmrStarFile(file,version = version)

    nmrStarFile.read(verbose = 1)

    for constraintFile in nmrStarFile.constraintFiles:
      for constraint in constraintFile.constraints:

        print constraint.Id, constraint.value, constraint.error, constraint.lowerValue, constraint.upperValue

        for item in constraint.nodes:
          mlist = []
          for member in item.members:
            mlist.append([member.seqCode,member.atomName])
          print "   " + str(mlist)
