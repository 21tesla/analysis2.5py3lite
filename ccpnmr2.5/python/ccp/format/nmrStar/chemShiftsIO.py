"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for nmrStar chemical shift saveframe

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

import os

from ccp.format.nmrStar.generalIO import NmrStarGenericFile
from ccp.format.nmrStar.generalIO import NmrStarFile

from ccp.format.general.Util import getSeqAndInsertCode
from ccp.format.general.Constants import defaultMolCode

#####################
# Class definitions #
#####################

class NmrStarFile(NmrStarFile):

  def initialize(self, version = '2.1.1'):
  
    self.chemShiftFiles = []
    
    self.files = self.chemShiftFiles
    
    if not self.version:
      self.version = version
    
    self.saveFrameName = 'assigned_chemical_shifts'
    self.DataClassFile = NmrStarChemShiftFile
   
    self.components = [self.saveFrameName]
    self.setComponents()

  def read(self,verbose = 0):

    if os.path.exists(self.name):
      text = open(self.name,'r').read()
      if not text.count(self.saveFrameName) and text.count("loop_"):
        
         # Also set version - can't do this from info
        if text.count('_Atom_chem_shift'):
          sfCatText = "   _Assigned_chem_shift_list.Sf_category       assigned_chemical_shifts"
        else:
          sfCatText = "   _Saveframe_category               assigned_chemical_shifts"
        
        # Try to fix the missing saveframe, and set this as the text.
        loopIndex = text.index("loop_")     
        
        # Also try to catch saveframe indicator
        if text[:loopIndex].count("save_"):
          saveFrameStart = ""
        else:
          saveFrameStart = "save_shifts\n"
           
        saveFrameStart += "%s\n  loop_" % sfCatText
        text = text.replace("loop_",saveFrameStart,1)
        
        if not self.patt[self.format + 'EndSaveTag'].search(text):
          text += "save_\n"
        
        self.text = text
        print self.text
        print "Warning: fixing input file to be correct NMR-STAR. This might not work."
        
    self.readComponent(verbose = verbose) 

class NmrStarChemShiftFile(NmrStarGenericFile):

  """
  Information on file level
  """
  
  def initialize(self,parent,saveFrame = None):

    # 2.1 doesn't use this attribute.
 
    self.attrToTagMappings = (
      
      ('details','Details',None),
      ('chemShiftRefId','Chem_shift_reference_ID',None),

    )

    self.chemShifts = []
    
    self.saveFrame = saveFrame
    
    self.parent = parent
    self.version = parent.version

    if self.saveFrame:
      self.parseSaveFrame()

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    if self.version == '2.1.1':
      # In case is missing in loop_ only shift files.
      if self.saveFrame.tags.has_key('_Mol_system_component_name'):
        molCode = self.saveFrame.tags['_Mol_system_component_name']
      else:
        molCode = defaultMolCode
      
      if '_Atom_shift_assign_ID' not in self.saveFrame.tables.keys():
        # Possible apparently
        return
      
      chemShiftTableTags = self.saveFrame.tables['_Atom_shift_assign_ID'].tags
      numChemShifts = len(chemShiftTableTags['_Atom_shift_assign_ID'])
    else:
      if self.attrToTagMappings:

        for (attrName,tagName,default) in self.attrToTagMappings:
          if self.saveFrame.tags.has_key(tagName):
            attrValue = self.saveFrame.tags[tagName]
          else:
            attrValue = default
        
          if not hasattr(self,attrName) or not self.attrName:
            setattr(self,attrName,attrValue)

      molCode = None
      chemShiftTableTags = self.saveFrame.tables['_Atom_chem_shift'].tags
      if chemShiftTableTags.has_key('ID'):
        numChemShifts = len(chemShiftTableTags['ID'])
      elif chemShiftTableTags.has_key('Comp_ID'):
        numChemShifts = len(chemShiftTableTags['Comp_ID'])
      else:
        numChemShifts = 0
        

    for i in range(0,numChemShifts):
      
      tmpMolCode = molCode
      if not tmpMolCode:
        # TODO could in principle find entity_assembly name (or code)
        # What is the best way to handle this?
        if chemShiftTableTags.has_key('Entity_assembly_ID'):
          tmpMolCode = str(chemShiftTableTags['Entity_assembly_ID'][i])
        else:
          tmpMolCode = defaultMolCode

      self.chemShifts.append(NmrStarChemShift(tmpMolCode,self))       
      self.chemShifts[-1].setData(chemShiftTableTags,i)

    tableName = '_Chem_shift_experiment'
    self.setMeasureExperiments(tableName)

    tableName = '_Chem_shift_software'
    self.setMeasureSoftwares(tableName)

class NmrStarChemShift:

  def __init__(self,molCode,parent):

    self.molCode = molCode[:65]  # Make sure this is truncated! Used as resonance name later on...
    self.parent = parent
    self.seqCode = None
  
  def setData(self,chemShiftTableTags,i):
      
    # Values already formatted...

    if self.parent.version == '2.1.1':
      assignList = [['Id',       '_Atom_shift_assign_ID',None],
                    ['seqCode',  '_Residue_seq_code',None],
                    ['resLabel', '_Residue_label',None],
                    ['atomName', '_Atom_name',None],
                    ['atomType', '_Atom_type',None],
                    ['value',    '_Chem_shift_value',None],
                    ['valueError', '_Chem_shift_value_error',0.0],
                    ['ambCode',    '_Chem_shift_ambiguity_code',1],
                    ['ambCode',    '_Chem_shift_ambiguity_type',1], # This apparently also occurs sometimes...
                    ['authorSeqCode', '_Residue_author_seq_code',None]
                   ]
    else:
      # TODO can get a lot more out...
      assignList = [['Id',            'ID',None],
                    ['seqCode',       'Seq_ID',None],
                    ['backupSeqCode', 'Comp_index_ID',None],
                    ['resLabel',      'Comp_ID',None],
                    ['atomName',      'Atom_ID',None],
                    ['atomType',      'Atom_type',None],
                    ['value',         'Val',None],
                    ['valueError',    'Val_err',0.0],
                    ['figOfMerit',    'Assign_fig_of_merit',None],
                    ['ambCode',       'Ambiguity_code',1],
                    ['details',       'Details',None],
                    ['authorSeqCode', 'Author_seq_ID',None]
                   ]
      
    
    for (attrName,tagName,default) in assignList:

      if chemShiftTableTags.has_key(tagName):
        if chemShiftTableTags[tagName][i] != None:
          setattr(self,attrName,chemShiftTableTags[tagName][i])
        else:
          setattr(self,attrName,default)
    
    # Little hack - not sure what is best in the end...
    if not self.seqCode and hasattr(self,'backupSeqCode') and self.backupSeqCode:
      self.seqCode = self.backupSeqCode
          
    # For completeness...
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(self.seqCode)

    # This is a hack, is possible in older NMR-STAR files! Use resLabel as chain ID
    if self.seqCode == None and hasattr(self,'resLabel') and self.resLabel:
      self.seqCode = 1
      self.chainCode = self.resLabel
