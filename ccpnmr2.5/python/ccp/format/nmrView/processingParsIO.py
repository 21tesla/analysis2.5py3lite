"""
======================COPYRIGHT/LICENSE START==========================

processingParsIO.py: NmrPipe spectrum parameter reader

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

from memops.universal.Io import getTopDirectory, normalisePath, splitPath, joinPath
from memops.universal.Constants import dirsep
from memops.universal.Util import returnInt, returnFloat

from ccp.format.nmrView.generalIO import NmrViewGenericFile

#####################
# Class definitions #
#####################

class NmrViewProcessingParsFile(NmrViewGenericFile):

  def initialize(self):

    self.scriptName = 'processingParsIO.py'
    
    self.fPars = {}
    self.fPars['bigEndian'] = True
    self.fPars['numberType'] = 'float'
    #(void,self.fPars['inputFile']) = os.path.split(self.name)
    
    self.allPars = {
    
      'mainPars': [
    
        ('header',3,'Header size and header block size'),
        ('felix',1,'Felix data'),
        ('vnmr',1,'VNMR data'),
        ('dim',-1,'Number of dimensions'),
        ('posneg',2,'????'),
        ('lvl',2,'????'),
        ('scale',2,'????'),
        ('rdims',2,'????'),
        ('datatype',2,'Integer type?'),
        ('byteorder',2,'Swap data? Endianness?')
    
       ],
    
      'dimPars': [
      
        ('sw',3,'Spectral width'),
        ('sf',3,'Spectrometer frequence'),
        ('ref',4,'Reference ppm and point'),
        ('label',3,'Label')
    
       ]
    }
    
    self.allParKeys = []
    self.parKeys = {}
    
    for parType in self.allPars.keys():
      parKeys = []
      for parInfo in self.allPars[parType]:
        parKeys.append(parInfo[0])
        
      self.allParKeys.extend(parKeys)
      self.parKeys[parType] = parKeys
    
  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s processing parameters file %s" % (self.format,self.name)

    try:
      fin = open(self.name, 'rU')
    except IOError, e:
      raise str(e)

    lines = fin.readlines()

    fin.close()

    for line in lines:

      n = line.find('!')
      if (n >= 0):
        line = line[:n]

      line = line.strip()

      if (line):

        fields = line.split()
        self.parseFields(fields)

  def parseFields(self, fields):
  
    #
    # This is taken almost directly from Wayne's AzaraParams.py file
    #

    key = fields[0]
    assert key in self.allParKeys, 'key = %s' % key
    
    parType = numFields = None

    for tempParType in self.allPars.keys():
      if key in self.parKeys[tempParType]:
        for parInfo in self.allPars[tempParType]:
          if key == parInfo[0]:
            (origName,numFields,comment) = parInfo
            parType = tempParType
            break
            
      if parType:
        break
      
    if not numFields:
      return

    if numFields != -1:
      assert len(fields) == numFields, 'len(fields) = %d, should be %d' % (len(fields),numFields)

    
    if parType == 'mainPars':
      
      """
        ('felix',1,'Felix data'),
        ('vnmr',1,'VNMR data'),
        
        ('posneg',2,'????'),
        ('lvl',2,'????'),
        ('scale',2,'????'),
        ('rdims',2,'????'),
      """
      
      if key == 'header':
        self.fPars['headerSize'] = returnInt(fields[1])
        self.fPars['headerBlockSize'] = returnInt(fields[2])

      elif key == 'dim':
        self.numDim = returnInt(fields[1])
        self.fPars['numDim'] = self.numDim
        
        fieldNum = 2
        for fParName in ('numPoints','blockSize'):
          self.fPars[fParName] = [''] * self.numDim
          for dimIndex in range(self.numDim):
            self.fPars[fParName][dimIndex] = returnInt(fields[fieldNum])
            fieldNum += 1

      elif key == 'byteorder':
        if returnInt(fields[1]):
          self.fPars['bigEndian'] = False

      # TODO not clear what this parameter means!!
      elif key == 'datatype':
        if returnInt(fields[1]):
          self.fPars['numberType'] = 'int'
      
      elif numFields == 1:
        self.fPars[key] = 1
      
      elif numFields == 2:
        self.fPars[key] = fields[1]
      
    elif parType == 'dimPars':
    
      dimIndex = returnInt(fields[1]) - 1
      
      if key in ('sw','sf'):
           
        if key == 'sw':
          fParName = 'spectralWidth'
        elif key == 'sf':
          fParName = 'spectrometerFreq'

        self.setDimValue(fParName,dimIndex,returnFloat(fields[2]))
        
      elif key == 'ref':
        
        self.setDimValue('refPpm',dimIndex,returnFloat(fields[2]))
        self.setDimValue('refPoint',dimIndex,returnFloat(fields[3]))

      elif key == 'label':
        
        self.setDimValue('label',dimIndex,fields[2])
        
        # TODO: this is a bit nasty...
        if 'H' in fields[2] or 'D' in fields[2]:
          self.setDimValue('nucleus',dimIndex,'1H')
        elif 'C' in fields[2]:
          self.setDimValue('nucleus',dimIndex,'13C')
        elif 'N' in fields[2]:
          self.setDimValue('nucleus',dimIndex,'15N')          
      
  def setDimValue(self,fParName,dimIndex,value):

    if not self.fPars.has_key(fParName):
      self.fPars[fParName] = [''] * self.numDim

    self.fPars[fParName][dimIndex] = value

###################
# Main of program #
###################
    
if __name__ == "__main__":

  files = ['../../reference/ccpNmr/atomCrystal/dl_nhsqc.par',
          ]

  
  for inFile in files:
    
    procParFile = NmrViewProcessingParsFile(os.path.join(getTopDirectory(),inFile))

    procParFile.read(verbose = 1)
    
    print procParFile.fPars
