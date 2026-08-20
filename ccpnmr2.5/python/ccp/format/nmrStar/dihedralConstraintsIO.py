#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

dihedralConstraintsIO.py: I/O for nmrStar dihedral constraints saveframe

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

from ccp.general.Constants import code1LetterToCcpCodeDict

from ccp.format.general.Constants import defaultMolCode, bioPolymerCodes
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
    
    self.saveFrameName = 'torsion_angle_constraints'
    self.DataClassFile = NmrStarDihedralConstraintFile
    
    self.components = [self.saveFrameName]
    self.setComponents()

  def read(self,verbose = 0):
  
    self.readComponent(verbose = verbose) 

class NmrStarDihedralConstraintFile(NmrStarConstraintFile):

  def initialize(self,parent,saveFrame = None):

    self.constraints = []
    
    self.constraintElements = 4
    
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
    # Parse the constraints
    #
    
    tableTagNames = {}
    
    if self.version == '2.1.1':

      tableName = '_Constraint_ID'
      tableTagNames['constraintID'] = tableName
      tableTagNames['atomPresent'] = "_Atom_one_atom_name"
      tableTagNames['position'] = ['one','two','three','four']
      
      tableTagNames['name'] = '_Angle_name'
      tableTagNames['upperAngle'] = '_Angle_upper_bound_value'
      tableTagNames['lowerAngle'] = '_Angle_lower_bound_value'
      tableTagNames['energyCst'] = '_Force_constant_value'
      tableTagNames['funcExp'] = '_Potential_function_exponent'
 
      tableTagNames['chainCode'] = "_Atom_%s_mol_system_component_ID"
      tableTagNames['seqCode'] = "_Atom_%s_residue_seq_code"
      tableTagNames['atomName'] = "_Atom_%s_atom_name"
      tableTagNames['resLabel'] = "_Atom_%s_residue_label"

    elif self.version[0] == '3':
    
      tableName = '_Torsion_angle_constraint'
      tableTagNames['constraintID'] = 'ID'
      tableTagNames['position'] = ['1','2','3','4']
      
      tableTagNames['name'] =       'Torsion_angle_name'
      tableTagNames['upperAngle'] = 'Angle_upper_bound_val'
      tableTagNames['lowerAngle'] = 'Angle_lower_bound_val'
      
      if self.version == '3.0':
        tableTagNames['energyCst'] =  'Force_constant_value'
        tableTagNames['funcExp'] =    'Potential_function_exponent'
      
        #
        # TODO: NOTE that original input here is different from output!!!!
        # need to take that into account when mapping!!
        #

        tableTagNames['atomPresent'] = "Label_atom_ID_1"
        tableTagNames['chainCode'] = "Label_entity_assembly_ID_%s"
        tableTagNames['seqCode'] = "Label_comp_index_ID_%s"
        tableTagNames['atomName'] = "Label_atom_ID_%s"
        tableTagNames['resLabel'] = "Label_comp_ID_%s"
      
      elif self.version == '3.1':

        tableTagNames['atomPresent'] = "Atom_ID_1"
        tableTagNames['chainCode'] =   "Entity_assembly_ID_%s"
        tableTagNames['seqCode'] =     "Comp_index_ID_%s"
        tableTagNames['atomName'] =    "Atom_ID_%s"
        tableTagNames['resLabel'] =    "Comp_ID_%s"
      
      if self.saveFrame.tables.has_key(tableName) and not self.saveFrame.tables[tableName].tags[tableTagNames['seqCode'] % 1][0]:
        tableTagNames['atomPresent'] = "Auth_atom_ID_1"
        
        if self.version == '3.0':
          tableTagNames['chainCode'] = "Auth_segment_code_%s"
        elif self.version == '3.1':
          tableTagNames['chainCode'] = "Auth_asym_ID_%s"
          
        tableTagNames['seqCode'] = "Auth_seq_ID_%s"
        tableTagNames['atomName'] = "Auth_atom_ID_%s"
        tableTagNames['resLabel'] = "Auth_comp_ID_%s"

    if self.saveFrame.tables.has_key(tableName):

      constraintTableTags = self.saveFrame.tables[tableName].tags

      for i in range(0,len(constraintTableTags[tableTagNames['constraintID']])):

        dihedralConstraint = NmrStarDihedralConstraint()
        dihedralConstraint.setData(constraintTableTags,i,self.cyanaLib,tableTagNames)

        self.constraints.append(dihedralConstraint)

        dihedralConstraint.setErrors(self.saveFrame.tables[tableName].tagErrors,tableTagNames,i)

        if dihedralConstraint.name and not constraintTableTags[tableTagNames['atomPresent']][i]:

          #
          # Note that following will set item and members automatically
          #

          dihedralConstraint.setAtomMembers(constraintTableTags,i,tableTagNames)
          self.cyanaLibUsed = 1

        else:

          dihedralConstraint.nodes.append(NmrStarDihedralConstraintItem())

          for position in tableTagNames['position']:

            dihedralConstraint.nodes[-1].members.append(NmrStarDihedralConstraintMember())
            dihedralConstraint.nodes[-1].members[-1].addStarInfo(constraintTableTags,i,position,tableTagNames)
        
    self.parseCommentsLoop()
      
class NmrStarDihedralConstraint(GenericConstraint):

  def __init__(self, Id = None):
  
    self.Id = Id
    self.nodes = []   # THIS IS A HACK! Not really a node...
    self.errors = []
    
    self.name = None

  def setData(self,constraintLogicTableTags,i,cyanaLib,tableTagNames):
    
    self.cyanaLib = cyanaLib
    self.defaultMolCode = defaultMolCode

    self.Id = constraintLogicTableTags[tableTagNames['constraintID']][i]
    self.name = constraintLogicTableTags[tableTagNames['name']][i]
    self.upperAngle = constraintLogicTableTags[tableTagNames['upperAngle']][i]
    self.lowerAngle = constraintLogicTableTags[tableTagNames['lowerAngle']][i]
    
    for attrName in ('energyCst','funcExp'):
      if tableTagNames.has_key(attrName):
        setattr(self,attrName,constraintLogicTableTags[tableTagNames[attrName]][i])
    
  def setAtomMembers(self,constraintTableTags,i,tableTagNames):
   
    #
    # TODO: this is rather nasty: assuming it's DIANA if dihedral angle name given.
    #       should really be checking type... leave this for now
    #

    chainCode = constraintTableTags[tableTagNames['chainCode'] % tableTagNames['position'][0]][i]
    
    if not chainCode:
      chainCode = defaultMolCode

    # Can only handle real numbers here - no insert codes allowed...
    (seqCode,seqInsertCode) = getSeqAndInsertCode(constraintTableTags[tableTagNames['seqCode'] % tableTagNames['position'][0]][i])
    resLabel = constraintTableTags[tableTagNames['resLabel'] % tableTagNames['position'][0]][i]
    
    refAngle = self.cyanaLib.findAngle(resLabel,self.name)
    
    # Special casing RNA/DNA!
    if not refAngle:
      # Assuming RNA
      if len(resLabel) == 1 and resLabel[0] in code1LetterToCcpCodeDict['RNA'].keys():
        molType = 'RNA'
      # Assuming DNA
      elif len(resLabel) == 2 and resLabel[0] == 'D' and resLabel[1] in code1LetterToCcpCodeDict['DNA'].keys():
        molType = 'DNA'
      else:
        print "    Unrecognized dihedral angle '%s' for residue label '%s' - ignoring..." % (self.name,resLabel)
        return
        
      for i in range(len(bioPolymerCodes[molType][0])):
        if resLabel == bioPolymerCodes[molType][0][i]:
          if molType == 'DNA':
            refCode = bioPolymerCodes[molType][1][i]
          else:
            refCode = "R%s" % bioPolymerCodes[molType][1][i]
            
          #print "Looking for %s for %s..." % (refCode,self.name)
          refAngle = self.cyanaLib.findAngle(refCode,self.name)
          if refAngle:
            break
    
    self.nodes.append(NmrStarDihedralConstraintItem())
    
    if refAngle:
      for atom in refAngle.atoms:

        if atom.location:
          curSeqCode = seqCode + atom.location
          curResLabel = ""

        else:
          curSeqCode = seqCode
          curResLabel = resLabel

        self.setAtomMember(self.nodes[-1],chainCode,curSeqCode,atom.name,curResLabel)

    else:
      print "  Error: no reference angle found for name %s..." % self.name
      for i in range(4):
        curSeqCode = seqCode
        curResLabel = resLabel
        self.setAtomMember(self.nodes[-1],chainCode,curSeqCode,'',curResLabel)
    
  def setAtomMember(self,node,chainCode,seqCode,atomName,resLabel):
  
    node.members.append(NmrStarDihedralConstraintMember(chainCode,seqCode,atomName,resLabel))

class NmrStarDihedralConstraintItem:

  def __init__(self):
    
    self.members = []    

#
# TODO: MERGE with RDC constraints!!
#

class NmrStarDihedralConstraintMember:

  def __init__(self,chainCode = None,seqCode = None, atomName = None, resLabel = None):
    
    self.chainCode = chainCode
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.atomName = atomName
    self.resLabel = resLabel
    
    self.defaultMolCode = defaultMolCode
    
  def addStarInfo(self,constraintTableTags,i,position,tableTagNames):

    self.chainCode = constraintTableTags[tableTagNames['chainCode'] % position][i]
    
    if not self.chainCode:
      self.chainCode = self.defaultMolCode

    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(constraintTableTags[tableTagNames['seqCode'] % position][i])
    self.atomName = constraintTableTags[tableTagNames['atomName'] % position][i]
    self.resLabel = constraintTableTags[tableTagNames['resLabel'] % position][i]

###################
# Main of program #
###################

if __name__ == "__main__":  

  files = [#['../reference/nmrStar/1afi.dihed.restraints','2.1.1'],
           #['../reference/ccpNmr/aartUtrecht/1byx/restraints.star','2.1.1'],
           #['../reference/ccpNmr/aartUtrecht/1bct/restraints.star','2.1.1'],
           #['../reference/ccpNmr/aartUtrecht/1d3z/restraints.star','2.1.1'],
           #['../reference/ccpNmr/aartUtrecht/1bct/restraints.star','2.1.1'],
           #['../reference/ccpNmr/aartUtrecht/1bsh/restraints.star','2.1.1'],
           #['../reference/ccpNmr/aartUtrecht/1fht/restraints.star','2.1.1'],
           #['../reference/ccpNmr/jurgenBmrb/1d8b/restraints.star','2.1.1']
           #['../reference/ccpNmr/jurgenBmrb/1n8x/1n8x.str','3.0']
           ['../reference/ccpNmr/jurgenBmrb/1ao9/1ao9.str','3.0']
          ]
  
  for (file,version) in files:
    
    file = os.path.join(getTopDirectory(), file)
    
    nmrStarFile = NmrStarFile(file,version = version)

    nmrStarFile.read(verbose = 1)

    for constraintFile in nmrStarFile.constraintFiles:
      #print constraintFile.comments
      for constraint in constraintFile.constraints:

        print constraint.Id, constraint.name, constraint.lowerAngle, constraint.upperAngle

        if constraint.errors:
          print constraint.errors

        for item in constraint.nodes:
          mlist = []
          for member in item.members:
            mlist.append([member.seqCode,member.seqInsertCode,member.atomName,member.resLabel])
          print "   " + str(mlist)
