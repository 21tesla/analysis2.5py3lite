"""
======================COPYRIGHT/LICENSE START==========================

jCouplingIO.py: I/O for nmrStar J coupling saveframe

Copyright (C) 2008 Chris Penkett (European Bioinformatics Institute)

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

- contact Chris Penkett (penkett@ebi.ac.uk)
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

from ccp.format.nmrStar.generalIO import NmrStarGenericFile
from ccp.format.nmrStar.generalIO import NmrStarFile

from ccp.format.general.Util import getSeqAndInsertCode
from ccp.format.general.Constants import defaultMolCode

#####################
# Class definitions #
#####################

class NmrStarFile(NmrStarFile):

  def initialize(self, version = '3.1'):
  
    self.jCouplingFiles = []
    
    self.files = self.jCouplingFiles
    
    if not self.version:
      self.version = version
    
    self.saveFrameName = 'coupling_constants'
    self.DataClassFile = NmrStarJCouplingFile
   
    self.components = [self.saveFrameName]
    self.setComponents()

  def read(self,verbose = 0):

    self.readComponent(verbose = verbose) 


class NmrStarJCouplingFile(NmrStarGenericFile):

  """
  Information on file level
  """
  
  def initialize(self,parent,saveFrame = None):

    
    # Warning: should in principle be version specific, this is 3.1
    self.attrToTagMappings = (
      
        ('details','Details',None),
        ('specFreq','Spectrometer_frequency_1H',None),
    
    )

    self.jCouplingValues = []
    
    self.saveFrame = saveFrame
    
    self.parent = parent
    self.version = parent.version

    if self.saveFrame:
      self.parseSaveFrame()

  def parseSaveFrame(self):

    if not self.checkVersion():
      return
    
    #
    # Read the table with the measurements
    #
  
    tmpMolCode1 = tmpMolCode2 = None

    if self.version == '2.1.1':
      if self.saveFrame.tags.has_key('_Mol_system_component_name'):
        tmpMolCode1 = tmpMolCode2 = self.saveFrame.tags['_Mol_system_component_name']
        
      if '_Coupling_constant_code' not in self.saveFrame.tables.keys():
        return
      
      jCouplingTableTags = self.saveFrame.tables['_Coupling_constant_code'].tags
      numJCouplingValues = len(jCouplingTableTags['_Coupling_constant_code'])
      
    else:
      if self.attrToTagMappings:

        for (attrName,tagName,default) in self.attrToTagMappings:
          if self.saveFrame.tags.has_key(tagName):
            attrValue = self.saveFrame.tags[tagName]
          else:
            attrValue = default
        
          if not hasattr(self,attrName) or not self.attrName:
            setattr(self,attrName,attrValue)

      tableName = '_Coupling_constant'

      if self.saveFrame.tables.has_key(tableName):

        jCouplingTableTags = self.saveFrame.tables[tableName].tags
        numJCouplingValues = len(jCouplingTableTags['ID'])

    for i in range(numJCouplingValues):
            
      if not tmpMolCode1:
        if jCouplingTableTags.has_key('Entity_assembly_ID_1'):
          tmpMolCode1 = str(jCouplingTableTags['Entity_assembly_ID_1'][i])
        else:
          tmpMolCode1 = defaultMolCode
    
      if not tmpMolCode2:
        if jCouplingTableTags.has_key('Entity_assembly_ID_2'):
          tmpMolCode2 = str(jCouplingTableTags['Entity_assembly_ID_2'][i])
        else:
          tmpMolCode2 = defaultMolCode

      self.jCouplingValues.append(NmrStarJCoupling(tmpMolCode1,tmpMolCode2,self))       
      self.jCouplingValues[-1].setData(jCouplingTableTags,i)  

    tableName = '_Coupling_constant_experiment'
    self.setMeasureExperiments(tableName)

    tableName = '_Coupling_constant_software'
    self.setMeasureSoftwares(tableName)

class NmrStarJCoupling:

  def __init__(self,molCode1,molCode2,parent):

    self.molCode1 = molCode1
    self.molCode2 = molCode2
    self.parent = parent
  
  def setData(self,jCouplingTableTags,i):

    if self.parent.version == '2.1.1':
      assignList = [['type',      '_Coupling_constant_code',None],
                    ['seqCode1',  '_Atom_one_residue_seq_code',None],
                    ['resLabel1', '_Atom_one_residue_label',None],
                    ['atomName1', '_Atom_one_name',None],
                    ['seqCode2',  '_Atom_two_residue_seq_code',None],
                    ['resLabel2', '_Atom_two_residue_label',None],
                    ['atomName2', '_Atom_two_name',None],
                    ['value',     '_Coupling_constant_value',None],
                    ['valueMin',  '_Coupling_constant_min_value',None],
                    ['valueMax',  '_Coupling_constant_max_value',None],
                    ['valueError','_Coupling_constant_value_error',0.0]
                   ]
    else:  
      assignList = [['Id',             'ID',None],
                    ['seqCode1',       'Seq_ID_1',None],
                    ['resLabel1',      'Comp_ID_1',None],
                    ['atomName1',      'Atom_ID_1',None],
                    ['atomType1',      'Atom_type_1',None],
                    ['seqCode2',       'Seq_ID_2',None],
                    ['resLabel2',      'Comp_ID_2',None],
                    ['atomName2',      'Atom_ID_2',None],
                    ['atomType2',      'Atom_type_2',None],
                    ['value',          'Val',None],
                    ['valueMin',       'Val_min',None],
                    ['valueMax',       'Val_max',None],
                    ['valueError',     'Val_err',0.0],
                    #['valueMin',       'Val_min',0.0],
                    #['valueMax',       'Val_max',0.0],
                    ['authorSeqCode1', 'Author_seq_ID_1',None],
                    ['authorSeqCode2', 'Author_seq_ID_2',None]
                    ]

    for (attrName,tagName,default) in assignList:

      if jCouplingTableTags.has_key(tagName):
        if jCouplingTableTags[tagName][i] != None:
          setattr(self,attrName,jCouplingTableTags[tagName][i])
        else:
          setattr(self,attrName,default)
          
    # For completeness...
    (self.seqCode1,self.seqInsertCode1) = getSeqAndInsertCode(self.seqCode1)
    (self.seqCode2,self.seqInsertCode2) = getSeqAndInsertCode(self.seqCode2)

###################
# Main of program #
###################

if __name__ == "__main__":  

  files = ['/homes/penkett/project/nmrStar/files/new2/bmr15834_3.str']
  
  for file in files:
    
    #file = os.path.join(getTopDirectory(), file)
    
    nmrStarFile = NmrStarFile(file, version='3.1')

    nmrStarFile.read(verbose = 1)

    for jCouplingFile in nmrStarFile.jCouplingFiles:
      for jCoupling in jCouplingFile.jCouplingValues:
        print jCoupling.Id, jCoupling.seqCode1, jCoupling.resLabel1, jCoupling.atomName1, jCoupling.seqCode2, jCoupling.resLabel2, jCoupling.atomName2, jCoupling.value, jCoupling.valueError
