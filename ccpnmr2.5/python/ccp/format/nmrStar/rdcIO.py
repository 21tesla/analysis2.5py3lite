"""
======================COPYRIGHT/LICENSE START==========================

rdcIO.py: I/O for nmrStar RDC saveframe

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

  def initialize(self, version = '2.1.1'):
  
    self.rdcFiles = []
    
    self.files = self.rdcFiles
    
    if not self.version:
      self.version = version
    
    if self.version[0] == '3':
      self.saveFrameName = 'RDCs'
    else:
      self.saveFrameName = 'residual_dipolar_couplings'
      
    self.DataClassFile = NmrStarRdcFile
   
    self.components = [self.saveFrameName]
    self.setComponents()

  def read(self,verbose = 0):

    self.readComponent(verbose = verbose) 


class NmrStarRdcFile(NmrStarGenericFile):

  """
  Information on file level
  """
  
  def initialize(self,parent,saveFrame = None):

    
    # Warning: should in principle be version specific, this is 3.1
    self.attrToTagMappings = (
      
        ('details','Details',None),
        ('specFreq','Spectrometer_frequency_1H',None),
    
    )

    self.rdcValues = []
    
    self.saveFrame = saveFrame
    
    self.parent = parent
    self.version = parent.version

    if self.saveFrame:
      self.parseSaveFrame()

  def parseSaveFrame(self):

    if not self.checkVersion():
      return
    
    #
    # Set version-specific info
    #
        
    if self.version == '2.1.1':
      
      rdcTableTags = self.saveFrame.tables['_Residual_dipolar_coupling_ID'].tags
      numRdcs = len(rdcTableTags['_Residual_dipolar_coupling_ID'])

    else:
      
      if self.attrToTagMappings:

        for (attrName,tagName,default) in self.attrToTagMappings:
          if self.saveFrame.tags.has_key(tagName):
            attrValue = self.saveFrame.tags[tagName]
          else:
            attrValue = default
        
          if not hasattr(self,attrName) or not self.attrName:
            setattr(self,attrName,attrValue)

      
      rdcTableTags = self.saveFrame.tables['_RDC'].tags
      numRdcs = len(rdcTableTags['ID'])
          
    #
    # Read the table with the measurements
    #

    for i in range(0,numRdcs):

      if rdcTableTags.has_key('Entity_assembly_ID_1'):
        molCode1 = str(rdcTableTags['Entity_assembly_ID_1'][i])
        molCode2 = str(rdcTableTags['Entity_assembly_ID_2'][i])

      elif rdcTableTags.has_key('_Atom_one_mol_system_component_name') and rdcTableTags['_Atom_one_mol_system_component_name'][i]:
        molCode1 = rdcTableTags['_Atom_one_mol_system_component_name'][i]
        molCode2 = rdcTableTags['_Atom_two_mol_system_component_name'][i]

      else:
        molCode1 = molCode2 = defaultMolCode

      self.rdcValues.append(NmrStarRdc(molCode1,molCode2,self))       
      self.rdcValues[-1].setData(rdcTableTags,i)

    tableName = '_RDC_experiment'
    self.setMeasureExperiments(tableName)

    tableName = '_RDC_software'
    self.setMeasureSoftwares(tableName)

class NmrStarRdc:

  def __init__(self,molCode1,molCode2,parent):

    self.molCode1 = molCode1
    self.molCode2 = molCode2
    self.parent = parent
  
  def setData(self,rdcTableTags,i):

    if self.parent.version == '2.1.1':
      assignList = [['Id',             '_Residual_dipolar_coupling_ID',None],
                    ['seqCode1',       '_Atom_one_residue_seq_code',None],
                    ['resLabel1',      '_Atom_one_residue_label',None],
                    ['atomName1',      '_Atom_one_atom_name',None],
                    ['atomType1',      '_Atom_one_atom_type',None],
                    ['seqCode2',       '_Atom_two_residue_seq_code',None],
                    ['resLabel2',      '_Atom_two_residue_label',None],
                    ['atomName2',      '_Atom_two_atom_name',None],
                    ['atomType2',      '_Atom_two_atom_type',None],
                    ['value',          '_Residual_dipolar_coupling_value',None],
                    ['valueError',     '_Residual_dipolar_coupling_value_error',0.0],
                    ['valueMin',       '_Residual_dipolar_coupling_min_value',0.0],
                    ['valueMax',       '_Residual_dipolar_coupling_max_value',0.0],
                    #['authorSeqCode1', 'Author_seq_ID_1',None],
                    #['authorSeqCode2', 'Author_seq_ID_2',None]
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
                    ['valueError',     'Val_err',0.0],
                    #['valueMin',       'Val_min',0.0],
                    #['valueMax',       'Val_max',0.0],
                    ['authorSeqCode1', 'Author_seq_ID_1',None],
                    ['authorSeqCode2', 'Author_seq_ID_2',None]
                    ]

    for (attrName,tagName,default) in assignList:

      if rdcTableTags.has_key(tagName):
        if rdcTableTags[tagName][i] != None:
          setattr(self,attrName,rdcTableTags[tagName][i])
        else:
          setattr(self,attrName,default)
          
    # For completeness...
    (self.seqCode1,self.seqInsertCode1) = getSeqAndInsertCode(self.seqCode1)
    (self.seqCode2,self.seqInsertCode2) = getSeqAndInsertCode(self.seqCode2)

###################
# Main of program #
###################

if __name__ == "__main__":  

  files = ['/homes/penkett/project/nmrStar/files/new/bmr15655_3.str']
  
  for file in files:
    
    #file = os.path.join(getTopDirectory(), file)
    
    nmrStarFile = NmrStarFile(file, version='3.1')

    nmrStarFile.read(verbose = 1)

    for rdcFile in nmrStarFile.rdcFiles:
      for rdc in rdcFile.rdcValues:
        print rdc.Id, rdc.seqCode1, rdc.resLabel1, rdc.atomName1, rdc.seqCode2, rdc.resLabel2, rdc.atomName2, rdc.value, rdc.valueError
