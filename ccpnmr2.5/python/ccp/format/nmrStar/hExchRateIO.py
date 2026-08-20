#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

hExchRateIO.py: I/O for nmrStar proton exchange rates saveframe

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
  
    self.hExchRateFiles = []
    
    self.files = self.hExchRateFiles
    
    if not self.version:
      self.version = version
    
    self.saveFrameName = 'H_exch_rates'
    self.DataClassFile = NmrStarHExchRateFile
   
    self.components = [self.saveFrameName]
    self.setComponents()

  def read(self,verbose = 0):

    self.readComponent(verbose = verbose) 


class NmrStarHExchRateFile(NmrStarGenericFile):

  """
  Information on file level
  """
  
  def initialize(self,parent,saveFrame = None):

    self.attrToTagMappings = (
      
        ('details','Details',None),
        ('units','Val_units',None),
    
    )

    self.hExchRateValues = []
    
    self.saveFrame = saveFrame
    
    self.parent = parent
    self.version = parent.version

    if self.saveFrame:
      self.parseSaveFrame()

  def parseSaveFrame(self):

    if not self.checkVersion():
      return

    for (attrName,tagName,default) in self.attrToTagMappings:
      if self.saveFrame.tags.has_key(tagName):
        attrValue = self.saveFrame.tags[tagName]
      else:
        attrValue = default
        
      if not hasattr(self,attrName) or not self.attrName:
        setattr(self,attrName,attrValue)

    tableName = '_H_exch_rate'

    if self.saveFrame.tables.has_key(tableName):
      hExchRateTableTags = self.saveFrame.tables[tableName].tags
      numHExchRateValues = len(hExchRateTableTags['ID'])

      for i in range(0,numHExchRateValues):
      
        tmpMolCode = str(hExchRateTableTags['Entity_assembly_ID'][i])

        self.hExchRateValues.append(NmrStarHExchRate(tmpMolCode,self))       
        self.hExchRateValues[-1].setData(hExchRateTableTags,i)

    tableName = '_H_exch_rate_experiment'
    self.setMeasureExperiments(tableName)

    tableName = '_H_exch_rate_software'
    self.setMeasureSoftwares(tableName)

class NmrStarHExchRate:

  def __init__(self,molCode,parent):

    self.molCode = molCode
    self.parent = parent
  
  def setData(self,hExchRateTableTags,i):
      
    assignList = [['Id',            'ID',None],
                  ['seqCode',       'Seq_ID',None],
                  ['resLabel',      'Comp_ID',None],
                  ['atomName',      'Atom_ID',None],
                  ['atomType',      'Atom_type',None],
                  ['value',         'Val',None],
                  ['valueError',    'Val_err',0.0],
#                  ['valueMin',      'Val_min',None],
#                  ['valueMax',      'Val_max',0.0],
                  ['details',       'Details',None],
                  ['authorSeqCode', 'Author_seq_ID',None]
                  ]

    for (attrName,tagName,default) in assignList:

      if hExchRateTableTags.has_key(tagName):
        if hExchRateTableTags[tagName][i] != None:
          setattr(self,attrName,hExchRateTableTags[tagName][i])
        else:
          setattr(self,attrName,default)
          
    # For completeness...
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(self.seqCode)

###################
# Main of program #
###################

if __name__ == "__main__":

  # No test data

  files = ['/homes/penkett/project/nmrstar/files/bmr15230_3_pub.str']
  
  for file in files:
    
    #file = os.path.join(getTopDirectory(), file)
    
    nmrStarFile = NmrStarFile(file, version='3.1')

    nmrStarFile.read(verbose = 1)

    for hExchRateFile in nmrStarFile.hExchRateFiles:
      for hExchRate in hExchRateFile.hExchRateValues:
        print hExchRate.Id, hExchRate.seqCode, hExchRate.resLabel, hExchRate.atomName, hExchRate.value, hExchRate.valueError
