"""
======================COPYRIGHT/LICENSE START==========================

orderParam.py: I/O for nmrStar order parameters saveframe

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

#####################
# Class definitions #
#####################

class NmrStarFile(NmrStarFile):

  def initialize(self, version = '3.1'):

    self.orderParamFiles = []
    
    self.files = self.orderParamFiles
    
    if not self.version:
      self.version = version
    
    if self.version[0] == '3':
      self.saveFrameName = 'order_parameters'
    else:
      self.saveFrameName = 'S2_parameters'
      
    self.DataClassFile = NmrStarOrderParamFile
   
    self.components = [self.saveFrameName]
    self.setComponents()

  def read(self,verbose = 0):

    self.readComponent(verbose = verbose) 

class NmrStarOrderParamFile(NmrStarGenericFile):

  """
  Information on file level
  """
  
  def initialize(self,parent,saveFrame = None):
    
    # Warning: should in principle be version specific, this is 3.1
    self.attrToTagMappings = (
      
        ('details','Details',None),

    )
  
    self.orderParamValues = []
    
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
      
      orderParamTableTags = self.saveFrame.tables['_Residue_seq_code'].tags
      numOrderParamValues = len(orderParamTableTags['_Residue_seq_code'])

    else:
      
      if self.attrToTagMappings:

        for (attrName,tagName,default) in self.attrToTagMappings:
          if self.saveFrame.tags.has_key(tagName):
            attrValue = self.saveFrame.tags[tagName]
          else:
            attrValue = default

          if not hasattr(self,attrName) or not self.attrName:
            setattr(self,attrName,attrValue)          

      tableName = '_Order_param'

      orderParamTableTags = self.saveFrame.tables[tableName].tags
      numOrderParamValues = len(orderParamTableTags['ID'])
      
    
    #
    # Read the table with the measurements
    #

    for i in range(numOrderParamValues):

      if self.version == '2.1.1':
        tmpMolCode = str(self.saveFrame.tags['_Mol_system_component_name'])      
      else:
        tmpMolCode = str(orderParamTableTags['Entity_assembly_ID'][i])

      self.orderParamValues.append(NmrStarOrderParam(tmpMolCode,self))       
      self.orderParamValues[-1].setData(orderParamTableTags,i)
    
    #
    # Set additional info (V3 only)
    #
    
    tableName = '_Order_parameter_experiment'
    self.setMeasureExperiments(tableName)

    tableName = '_Order_parameter_software'
    self.setMeasureSoftwares(tableName)

class NmrStarOrderParam:

  def __init__(self,molCode,parent):

    self.molCode = molCode
    self.parent = parent
  
  def setData(self,orderParamTableTags,i):

    if self.parent.version == '2.1.1':

      assignList = [['Id',            'SET', i+1],
                    ['seqCode',       '_Residue_seq_code',None],
                    ['resLabel',      '_Residue_label',None],
                    ['atomName',      '_Atom_name',None],
                    ['s2Value',       '_S2_value',None],
                    ['s2Error',       '_S2_value_fit_error',0.0],
                    ['tauEValue',     '_Tau_e_value',None],
                    ['tauEError',     '_Tau_e_value_fit_error',0.0],
                    ]
    
    else:
    
      assignList = [['Id',            'ID',None],
                    ['seqCode',       'Seq_ID',None],
                    ['resLabel',      'Comp_ID',None],
                    ['atomName',      'Atom_ID',None],
                    ['atomType',      'Atom_type',None],
                    ['tauEValue',     'Tau_e_val',None],
                    ['tauEError',      'Tau_e_val_fit_err',0.0],
                    ['authorSeqCode', 'Author_seq_ID',None]
                    ]

    for (attrName,tagName,default) in assignList:

      if orderParamTableTags.has_key(tagName):
        if orderParamTableTags[tagName][i] != None:
          setattr(self,attrName,orderParamTableTags[tagName][i])
        else:
          setattr(self,attrName,default)

      elif tagName == 'SET':
        setattr(self,attrName,default)
          
    # For completeness...
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(self.seqCode)

###################
# Main of program #
###################

if __name__ == "__main__":  

  files = ['/data/work/archives/bmrb/data/bmr15536/bmr15536_21.str']
  
  for file in files:
    
    #file = os.path.join(getTopDirectory(), file)
    
    nmrStarFile = NmrStarFile(file, version='2.1.1')

    nmrStarFile.read(verbose = 1)

    for orderParamFile in nmrStarFile.orderParamFiles:
      for orderParam in orderParamFile.orderParamValues:
        print orderParam.Id, orderParam.seqCode, orderParam.resLabel, orderParam.atomName, orderParam.s2Value, orderParam.tauEValue, orderParam.tauEError
