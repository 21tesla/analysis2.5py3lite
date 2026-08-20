#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

t2RelaxIO.py: I/O for nmrStar T2 relaxation saveframe

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
  
    self.t2RelaxFiles = []
    
    self.files = self.t2RelaxFiles
    
    if not self.version:
      self.version = version
    
    self.saveFrameName = 'heteronucl_T2_relaxation'
    self.DataClassFile = NmrStarT2RelaxFile
   
    self.components = [self.saveFrameName]
    self.setComponents()

  def read(self,verbose = 0):

    self.readComponent(verbose = verbose) 


class NmrStarT2RelaxFile(NmrStarGenericFile):

  """
  Information on file level
  """
  
  def initialize(self,parent,saveFrame = None):

    # Warning: should in principle be version specific, this is 3.1
    self.attrToTagMappings = (
      
        ('details','Details',None),
        ('units','T2_val_units',None),
        ('specFreq','Spectrometer_frequency_1H',None),
        ('tempCalibMethod','Temp_calibration_method',None),
        ('tempControlMethod','Temp_control_method', None),
        ('coherenceType','T2_coherence_type',None),
    
    )

    self.t2RelaxValues = []
    
    self.saveFrame = saveFrame
    
    self.parent = parent
    self.version = parent.version

    if self.saveFrame:
      self.parseSaveFrame()

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    #
    # Set saveframe (measurement list) level information.
    #

    for (attrName,tagName,default) in self.attrToTagMappings:
      if self.saveFrame.tags.has_key(tagName):
        attrValue = self.saveFrame.tags[tagName]
      else:
        attrValue = default
        
      if not hasattr(self,attrName) or not self.attrName:
        setattr(self,attrName,attrValue)
    
    #
    # Read the table with the measurements
    #

    tableName = '_T2'

    if self.saveFrame.tables.has_key(tableName):
      t2RelaxTableTags = self.saveFrame.tables[tableName].tags
      numT2RelaxValues = len(t2RelaxTableTags['ID'])

      for i in range(0,numT2RelaxValues):
      
        tmpMolCode = str(t2RelaxTableTags['Entity_assembly_ID'][i])

        self.t2RelaxValues.append(NmrStarT2Relax(tmpMolCode,self))       
        self.t2RelaxValues[-1].setData(t2RelaxTableTags,i)

    tableName = '_Heteronucl_T2_experiment'
    self.setMeasureExperiments(tableName)

    tableName = '_Heteronucl_T2_software'
    self.setMeasureSoftwares(tableName)

class NmrStarT2Relax:

  def __init__(self,molCode,parent):

    self.molCode = molCode
    self.parent = parent
  
  def setData(self,t2RelaxTableTags,i):
      
    assignList = [['Id',            'ID',None],
                  ['seqCode',       'Seq_ID',None],
                  ['resLabel',      'Comp_ID',None],
                  ['atomName',      'Atom_ID',None],
                  ['atomType',      'Atom_type',None],
                  ['value',         'T2_val',None],
                  ['valueError',    'T2_val_err',0.0],
                  ['rexValue',      'Rex_val',None],
                  ['rexValueError', 'Rex_err',0.0],
                  ['authorSeqCode', 'Author_seq_ID',None]
                  ]

    for (attrName,tagName,default) in assignList:

      if t2RelaxTableTags.has_key(tagName):
        if t2RelaxTableTags[tagName][i] != None:
          setattr(self,attrName,t2RelaxTableTags[tagName][i])
        else:
          setattr(self,attrName,default)
          
    # For completeness...
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(self.seqCode)

###################
# Main of program #
###################

if __name__ == "__main__":  

  files = ['/homes/penkett/project/nmrstar/files/bmr15230_3_pub.str']
  
  for file in files:
    
    #file = os.path.join(getTopDirectory(), file)
    
    nmrStarFile = NmrStarFile(file, version='3.1')

    nmrStarFile.read(verbose = 1)

    for t2RelaxFile in nmrStarFile.t2RelaxFiles:
      for t2Relax in t2RelaxFile.t2RelaxValues:
        print t2Relax.Id, t2Relax.seqCode, t2Relax.resLabel, t2Relax.atomName, t2Relax.value, t2Relax.valueError, t2Relax.rexValue, t2Relax.rexValueError
