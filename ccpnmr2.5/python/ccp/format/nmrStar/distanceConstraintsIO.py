"""
======================COPYRIGHT/LICENSE START==========================

distanceConstraintsIO.py: I/O for nmrStar distance constraints saveframes

Copyright (C) 2005-2008 Wim Vranken (European Bioinformatics Institute)

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

import os, string

from memops.universal.Io import getTopDirectory

# Import general functions
from ccp.format.nmrStar.util import returnStarInt

from ccp.format.nmrStar.generalIO import NmrStarConstraintFile
from ccp.format.nmrStar.generalIO import NmrStarFile
from ccp.format.nmrStar.generalIO import GenericConstraint

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
    
    # Note: this will create problems for general_distance_constraints. TODO FIX!! 
    self.saveFrameName = 'distance_constraints'
    self.DataClassFile = NmrStarDistanceConstraintFile
    
    self.components = [self.saveFrameName]
    self.setComponents()

  def read(self,verbose = 0):
  
    self.readComponent(verbose = verbose) 

class NmrStarDistanceConstraintFile(NmrStarConstraintFile):

  def initialize(self,parent,saveFrame = None):

    self.constraints = []
    
    self.constraintElements = 2
    
    self.saveFrame = saveFrame
    
    self.parent = parent
    self.version = parent.version
    
    if self.saveFrame:
      self.parseSaveFrame()

  def parseSaveFrame(self):

    findConstraint = {}
       
    if not self.checkVersion():
      return

    #
    # Nodes & their logic relationships
    #
    
    tableTagNames = {}
    
    isGeneralType = False
    
    if self.version == '2.1.1':

      tableName = '_Constraint_ID'
      tableTagNames['constraintID'] = tableName
      tableTagNames['nodeID'] = '_Constraint_tree_node_ID'
      tableTagNames['logicOperation'] = '_Constraint_tree_logic_operation'
      tableTagNames['downNodeID'] = '_Constraint_tree_down_node_ID'
      tableTagNames['rightNodeID'] = '_Constraint_tree_right_node_ID'
      
    elif self.version == '3.0':
    
      tableName = '_Dist_constraint_tree'
      tableTagNames['constraintID'] = 'ID'
      tableTagNames['nodeID'] = 'Node_ID'
      tableTagNames['logicOperation'] = 'Logic_operation'
      tableTagNames['downNodeID'] = 'Down_node_ID'
      tableTagNames['rightNodeID'] = 'Right_node_ID'

    elif self.version == '3.1':
    
      tableName = '_Dist_constraint_tree'
      if self.saveFrame.tables.has_key(tableName):
        tableTagNames['constraintID'] = 'Constraint_ID'
        tableTagNames['nodeID'] = 'Node_ID'
        tableTagNames['logicOperation'] = 'Logic_operation'
        tableTagNames['downNodeID'] = 'Down_node_ID'
        tableTagNames['rightNodeID'] = 'Right_node_ID'
      else:
        isGeneralType = True
        tableName = "_Gen_dist_constraint"
        tableTagNames['constraintID'] = 'ID'

        tableTagNames['memberID'] = 'Member_ID'
        tableTagNames['memberLogic'] = 'Member_logic_code'
        
        tableTagNames['chainCode'] = 'Entity_assembly_ID_%s'
        tableTagNames['seqCode'] = 'Comp_index_ID_%s'
        tableTagNames['atomName'] = 'Atom_ID_%s'
        tableTagNames['resLabel'] = "Comp_ID_%s"
        tableTagNames['target'] = 'Distance_val'
        tableTagNames['lowerBound'] = 'Distance_lower_bound_val'
        tableTagNames['upperBound'] = 'Distance_upper_bound_val'
        tableTagNames['intensity'] = 'Intensity_val'
    
    #
    # Could be missing...
    #
    
    if self.saveFrame.tables.has_key(tableName):
      
      constraintTableTags = self.saveFrame.tables[tableName].tags

      for i in range(len(constraintTableTags[tableTagNames['constraintID']])):

        Id = constraintTableTags[tableTagNames['constraintID']][i]

        if not findConstraint.has_key(Id):
          self.constraints.append(NmrStarDistanceConstraint())
          if not isGeneralType:
            self.constraints[-1].setData(constraintTableTags,i,tableTagNames)
          findConstraint[Id] = self.constraints[-1]
          
        elif not isGeneralType:
          self.constraints[-1].addNode(constraintTableTags,i,tableTagNames)

        # Always have to add a node for general type
        if isGeneralType:
          self.constraints[-1].nodes.append(NmrStarDistanceConstraintItem(Id = constraintTableTags[tableTagNames['memberID']][i]))

        self.constraints[-1].setErrors(self.saveFrame.tables[tableName].tagErrors,tableTagNames,i)
        
        # Deal with things immediately, all info available.
        if isGeneralType:
          
          # TODO also handle logic here?
          
          node = self.constraints[-1].nodes[-1]
          
          if node.Id != '1':
            node.logicOperation = constraintTableTags[tableTagNames['memberLogic']][i]
  
          for index in ('1','2'):
            
            tmpTableTagNames = {'memberID': 'Member_ID'} # Note: this is reset!
            
            for tableTagName in tableTagNames.keys():
              if tableTagNames[tableTagName].count("%s"):
                tmpTableTagNames[tableTagName] = tableTagNames[tableTagName] % index
            
            # Reset member ID - the memberId in gen_dist is really the node ID (or item)
            node.addMember(constraintTableTags,i,tmpTableTagNames)
            node.members[-1].Id = index
          
          node.addDistance(constraintTableTags,i,tableTagNames)

      #
      # Node members
      #
      
      if not isGeneralType:

        tableTagNames = {}
  
        if self.version == '2.1.1':
  
          tableName = '_Constraint_tree_node_member_constraint_ID'
          tableTagNames['constraintID'] = tableName
          tableTagNames['nodeID'] = '_Constraint_tree_node_member_node_ID'
  
          tableTagNames['memberID'] = '_Constraint_tree_node_member_ID'
          tableTagNames['fracVal'] = '_Contribution_fractional_value'
          tableTagNames['chainCode'] = '_Mol_system_component_code'
          tableTagNames['seqCode'] = '_Residue_seq_code'
          tableTagNames['atomName'] = '_Atom_name'
          tableTagNames['resLabel'] = "_Atom_residue_label"
  
        elif self.version == '3.0':
  
          tableName = '_Dist_constraint'
          tableTagNames['constraintID'] = 'Dist_constraint_tree_ID'
          tableTagNames['nodeID'] = 'Tree_node_member_node_ID'
  
          tableTagNames['memberID'] = 'Constraint_tree_node_member_ID'
          tableTagNames['fracVal'] = 'Contribution_fractional_val'
          
          tableTagNames['chainCode'] = 'Label_entity_assembly_ID'
          tableTagNames['seqCode'] = 'Label_comp_index_ID'
          tableTagNames['atomName'] = 'Label_atom_ID'
          tableTagNames['resLabel'] = "Label_comp_ID"
  
          # 3.0 or higher only...
          tableTagNames['chemCompID'] = 'Label_comp_ID'
  
        elif self.version == '3.1':
  
          tableName = '_Dist_constraint'
          tableTagNames['constraintID'] = 'Tree_node_member_constraint_ID'
          tableTagNames['nodeID'] = 'Tree_node_member_node_ID'
          tableTagNames['memberID'] = 'Constraint_tree_node_member_ID'
  
          #tableTagNames['fracVal'] = 'Contribution_fractional_val'
          
          tableTagNames['chainCode'] = 'Entity_assembly_ID'
          tableTagNames['seqCode'] = 'Comp_index_ID'
          tableTagNames['atomName'] = 'Atom_ID'
          tableTagNames['resLabel'] = "Comp_ID"
  
          # 3.0 or higher only...
          tableTagNames['chemCompID'] = 'Comp_ID'
  
        
        if self.version[0] == '3' and not self.saveFrame.tables[tableName].tags[tableTagNames['atomName']][0]:
            
          if self.version == '3.0':
            tableTagNames['chainCode'] = 'Auth_segment_code'
          elif self.version == '3.1':        
            tableTagNames['chainCode'] = 'Auth_asym_ID'
  
          tableTagNames['seqCode'] = 'Auth_seq_ID'
          tableTagNames['atomName'] = 'Auth_atom_ID'
          tableTagNames['resLabel'] = "Auth_comp_ID"
  
          # 3.0 or higher only...
          tableTagNames['chemCompID'] = 'Auth_comp_ID'
  
        nodeMembersTableTags = self.saveFrame.tables[tableName].tags
  
        for i in range(len(nodeMembersTableTags[tableTagNames['constraintID']])):
  
          Id = nodeMembersTableTags[tableTagNames['constraintID']][i]
  
          if not findConstraint.has_key(Id):
            print "  Error in nmrStar file: constraint %s is not defined but is referenced." % Id
  
          else:
            treeNodeId = nodeMembersTableTags[tableTagNames['nodeID']][i]
            constraint = findConstraint[Id]
            node = constraint.findNode(treeNodeId)
            node.addMember(nodeMembersTableTags,i,tableTagNames)
  
            constraint.setErrors(self.saveFrame.tables[tableName].tagErrors,tableTagNames,i)
  
        #
        # Distances (can be extended!)
        #
  
        tableTagNames = {}
  
        if self.version == '2.1.1':
  
          tableName = '_Distance_constraint_ID'
          tableTagNames['constraintID'] = tableName
          tableTagNames['nodeID'] = '_Distance_constraint_tree_node_ID'
  
          tableTagNames['target'] = '_Distance_value'
          tableTagNames['lowerBound'] = '_Distance_lower_bound_value'
          tableTagNames['upperBound'] = '_Distance_upper_bound_value'
  
        elif self.version == '3.0':
  
          tableName = '_Dist_constraint_value'
          tableTagNames['constraintID'] = 'Constraint_ID'
          tableTagNames['nodeID'] = 'Tree_node_ID'
  
          tableTagNames['target'] = 'Distance_val'
          tableTagNames['weight'] = 'Weight'
          tableTagNames['lowerBound'] = 'Distance_lower_bound_val'
          tableTagNames['upperBound'] = 'Distance_upper_bound_val'
          tableTagNames['intensity'] = 'Intensity_val'
  
        elif self.version == '3.1':
  
          tableName = '_Dist_constraint_value'
          tableTagNames['constraintID'] = 'Constraint_ID'
          tableTagNames['nodeID'] = 'Tree_node_ID'
  
          tableTagNames['target'] = 'Distance_val'
          tableTagNames['lowerBound'] = 'Distance_lower_bound_val'
          tableTagNames['upperBound'] = 'Distance_upper_bound_val'
          tableTagNames['intensity'] = 'Intensity_val'
  
        distanceTableTags = self.saveFrame.tables[tableName].tags
  
        for i in range(len(distanceTableTags[tableTagNames['constraintID']])):
  
          Id = distanceTableTags[tableTagNames['constraintID']][i]
  
          if not findConstraint.has_key(Id):
            print "  Error in nmrStar file: constraint %s is not defined but is referenced." % Id
  
          else:
            treeNodeId = distanceTableTags[tableTagNames['nodeID']][i]
            constraint = findConstraint[Id]
            node = constraint.findNode(treeNodeId)
            node.addDistance(distanceTableTags,i,tableTagNames)
  
            constraint.setErrors(self.saveFrame.tables[tableName].tagErrors,tableTagNames,i)

    #
    # Comments
    #
        
    self.parseCommentsLoop()

class NmrStarDistanceConstraint(GenericConstraint):

  def __init__(self, Id = None):
  
    self.Id = Id
    self.nodes = []
    self.errors = []
  
  def setData(self,constraintTableTags,i,tableTagNames):
  
    self.Id = constraintTableTags[tableTagNames['constraintID']][i]
    self.addNode(constraintTableTags,i,tableTagNames)
  
  def addNode(self,constraintTableTags,i,tableTagNames):
  
    self.nodes.append(NmrStarDistanceConstraintItem())
    self.nodes[-1].setData(constraintTableTags,i,tableTagNames)
    
  def addEmptyNode(self, Id = None):
  
    self.nodes.append(NmrStarDistanceConstraintItem(Id = Id))

  def findNode(self,treeNodeId):
  
    treeNodeId = returnStarInt(treeNodeId)
    for node in self.nodes:
      if treeNodeId == node.Id:
        return node
        
    print "  TreenodeID %d not found!" % treeNodeId
    
    return None

class NmrStarDistanceConstraintItem:

  #
  # Formerly TreeNode
  #

  def __init__(self, Id = None, constraintType = ""):
  
    self.Id = Id
    self.members = []
    self.type = constraintType
    self.logicOperation = None
    self.downId = None
    self.rightId = None   
  
  def setData(self,constraintLogicTableTags,i,tableTagNames):
  
    self.Id = constraintLogicTableTags[tableTagNames['nodeID']][i]
    self.logicOperation = constraintLogicTableTags[tableTagNames['logicOperation']][i]
    self.downId = constraintLogicTableTags[tableTagNames['downNodeID']][i]
    self.rightId = constraintLogicTableTags[tableTagNames['rightNodeID']][i]
  
  def addMember(self,nodeMembersTableTags,i,tableTagNames):
    #
    # Relationship between multiple members with same ID in node always OR!!
    #
    self.members.append(NmrStarDistanceConstraintMember())
    self.members[-1].setData(nodeMembersTableTags,i,tableTagNames)
    
  def addEmptyMember(self,Id = None):
  
    self.members.append(NmrStarDistanceConstraintMember(Id = Id))
    
  def addDistance(self,distanceTableTags,i,tableTagNames):
  
    #
    # TODO: COULD link the class attribute name to the tableTagNames dict keys...
    # then run loop. Also have room for function for special cases!!!
    #
          
    self.type = 'Distance'
    self.target = distanceTableTags[tableTagNames['target']][i]
    self.lowerBound = distanceTableTags[tableTagNames['lowerBound']][i]
    self.upperBound = distanceTableTags[tableTagNames['upperBound']][i]
    
    self.intensity = distanceTableTags[tableTagNames['intensity']][i]
    
    if tableTagNames.has_key('weight'):
      self.weight = distanceTableTags[tableTagNames['weight']][i]
      
class NmrStarDistanceConstraintMember:

  #
  # Formerly TreeNodeMember
  #
  
  def __init__(self,Id = None, chainCode = None,seqCode = None, atomName = None, resLabel = None):
  
    self.Id = Id
    self.fracVal = None
    self.chainCode = chainCode
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.atomName = atomName
    self.resLabel = resLabel

  def setData(self,nodeMembersTableTags,i,tableTagNames):

    self.Id = nodeMembersTableTags[tableTagNames['memberID']][i]
    
    if tableTagNames.has_key('fracVal'):
      self.fracVal = nodeMembersTableTags[tableTagNames['fracVal']][i]
    
    self.chainCode = nodeMembersTableTags[tableTagNames['chainCode']][i]
    if not self.chainCode:
      self.chainCode = defaultMolCode
    
    self.resLabel = nodeMembersTableTags[tableTagNames['resLabel']][i]
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(nodeMembersTableTags[tableTagNames['seqCode']][i])
    
    # This is a hack, is possible in NMR-STAR file from CNS/XPLOR data!
    if self.seqCode == None and self.chainCode:
      self.seqCode = 1
    
    self.atomName = nodeMembersTableTags[tableTagNames['atomName']][i]

    if tableTagNames.has_key('chemCompID'):
      self.chemCompID = nodeMembersTableTags[tableTagNames['chemCompID']][i]


###################
# Main of program #
###################

if __name__ == "__main__":  

  files = [#['../reference/nmrStar/1e91.restraints','2.1.1'],
           #['../reference/ccpNmr/aartUtrecht/1byx/restraints.star','2.1.1'],
           #['../reference/ccpNmr/aartUtrecht/1bsh/restraints.star','2.1.1'],
           #['../reference/ccpNmr/jurgenBmrb/1d8b/restraints.star','2.1.1']
           #['../reference/ccpNmr/jurgenBmrb/1d8b/1d8b.str','3.0']
           ['../reference/nmrStar/tmp/tmp/1ak7/restraints_nonredun.str','3.0']
           #['../reference/ccpNmr/jurgenBmrb/1ao9/1ao9.str','3.0']
          ]  
  for (file,version) in files:
    
    file = os.path.join(getTopDirectory(), file)
    
    nmrStarFile = NmrStarFile(file,version = version)

    nmrStarFile.read(verbose = 1)

    for constraintFile in nmrStarFile.constraintFiles:
      """
      for comment in constraintFile.comments:
        print comment
    """

      for constraint in constraintFile.constraints:
        if constraint.errors:
          print constraint.errors
        for node in constraint.nodes:
          """
          if hasattr(node,'weight') and node.weight:
            print node.weight,str(node.weight) # nmrStar 3.0 ONLY!!
          """
          if hasattr(node,'target'):
            print constraint.Id, node.Id, node.target, node.lowerBound, node.upperBound
          mlist = []
          for member in node.members:
            mlist.append([member.Id,member.seqCode,member.atomName])
          print "   " + str(mlist)
          """
    
          """
