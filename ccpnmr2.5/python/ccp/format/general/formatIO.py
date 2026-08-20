
"""
======================COPYRIGHT/LICENSE START==========================

formatIO.py: Generic superclass for format classes in subdirectories

Copyright (C) 2005-2009 Wim Vranken (European Bioinformatics Institute)

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

from memops.universal.Util import returnInt

from ccp.format.general.Util import getRegularExpressions, getSeqAndInsertCode
from ccp.format.general.Constants import defaultMolCode

class FileParseError(StandardError):

  def __init__(self, value):
    
    self.value = value
        
  def __str__(self):
    
    if self.value.count("\n"):
      returnStr = self.value
    else:
      returnStr = repr(self.value)
    
    return returnStr

class FormatFileError(FileParseError):
  
  pass

class FormatFile:

  FileParseError = FileParseError
  FormatFileError = FormatFileError

  def __init__(self,name,*args,**keywds):
  
    self.name = name
    self.newline = os.linesep
    
    self.version = None
    if 'version' in keywds.keys():
      self.version = keywds['version']
      del(keywds['version'])

    self.setGeneric()
    self.patt = getRegularExpressions(self.format)
    
    self.initialize(*args,**keywds)

  def setGeneric(self):
  
    self.format = None

class Citation:

  def __init__(self,parent,type,title,status,pubShortName,volume,issue,firstPage,year,publisherPlace,publisher = None,country = None,isPrimary = False,lastPage = None, pubLongName = None, pubMedId = None, details = None, bookTitle = None, bookSeries = None, medline = None, doi = None):

    self.parent = parent

    self.type = type
    self.title = title
    self.status = status
    
    self.pubShortName = pubShortName
    self.pubLongName = pubLongName

    self.pubMedId = pubMedId
    self.medline = medline
    self.doi = doi

    self.publisherPlace = publisherPlace
    self.publisher = publisher
    self.country = country
    
    self.volume = volume
    self.issue = issue
    self.firstPage = firstPage
    self.lastPage = lastPage
    self.year = year

    self.details = details
    self.isPrimary = isPrimary

    self.bookTitle = bookTitle
    self.bookSeries = bookSeries
    
    self.authors = []
    self.editors = []
    self.keywords = []
    
  def setAuthor(self,person):
    self.authors.append(person)
  
  def setEditor(self,person):
    self.editors.append(person)

  def setKeyword(self,person):
    self.keywords.append(person)
 
class Person:
  
  def __init__(self,parent,firstName,initials,familyName,familyTitle = None, title = None, serial = None):
    self.parent = parent
    
    self.firstName = firstName
    self.initials = initials    
    self.familyName = familyName
    self.familyTitle = familyTitle
    self.title = title
    
    self.serial = serial

class Keyword:

  def __init__(self,parent,keywd):

    self.parent = parent
    self.keywd = keywd

class Sequence:

  """
  Information on sequences stored in file
  """

  def __init__(self,molName = None, chainCode = defaultMolCode, *args, **keywds):

    if molName:
      self.molName = molName.strip()
    else:
      self.molName = molName
      
    self.chainCode = chainCode
    
    self.originalChainCode = chainCode
    
    self.elements = []
    
    self.isCircular = False
    
    self.secStrucInfo = {}
    
    self.setFormatSpecific(*args,**keywds)

  def setFormatSpecific(self,*args,**keywds):
  
    pass

class SequenceElement:
  
  """
  Information for an element in the sequence (residue)
  """
  
  def __init__(self,seqCode,*args,**keywds):
  
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    
    self.formatCode = None
    self.atomNames = []
    self.hasCisPeptideBond = False
    
    self.bonds = {}
    
    self.setResidueCode(*args)
    
    self.setFormatSpecific(*args,**keywds)
    
    self.secStrucInfo = None

  def setResidueCode(self,*args):
  
    if len(args) >= 1:
      self.code3Letter = args[0].upper()

  def setFormatSpecific(self,*args,**keywds):
  
    pass

  def setFormatCode(self,formatCode):
     
     self.formatCode = formatCode

  def addAtomName(self,atomName):
  
    self.atomNames.append(atomName)

  def setCisPeptideBond(self):
  
    self.hasCisPeptideBond = True
    
  def setBond(self,bondType,atomName,bondedSeqEl,bondedAtomName):
  
    if not self.bonds.has_key(bondType):
      self.bonds[bondType] = {}
    
    if not self.bonds[bondType].has_key(atomName):
      self.bonds[bondType][atomName] = []
    
    if not (bondedSeqEl,bondedAtomName) in self.bonds[bondType][atomName]:
      self.bonds[bondType][atomName].append((bondedSeqEl,bondedAtomName))
      print "  Set %s bond from %s.%s - %s.%s" % (bondType,self.seqCode,atomName,bondedSeqEl.seqCode,bondedAtomName)
    
  def setSecStrucInfo(self,secStrucType,secStrucSerial,specificInfo):
    
    self.secStrucInfo = (secStrucType,secStrucSerial,specificInfo)

class SpinSystem:

  """
  Information about a spin system
  """

  def __init__(self,spinSysCode,*args,**keywds):

    self.code = returnInt(spinSysCode)
    
    self.setFormatSpecific(*args,**keywds)

  def setFormatSpecific(self,*args,**keywds):
  
    pass

class GenericProcessingParsFile:
  
  """
  
  Generic processing parameters stuff
  
  """

  dimTypeFrequency = 'frequency'
  dimTypeFid = 'fid'
  dimTypeSampled = 'sampled'
  dimTypeNotProcessed= 'not processed'
  dimTypes = (dimTypeFrequency, dimTypeFid, dimTypeSampled, dimTypeNotProcessed)

  #
  # Processing parameters stuff - mapping between Wayne's stuff and the DataFormat stuff...
  #

  parNameMapping = {
  
      'ndim': ('numDim','Number of dimensions'),
      'file': ('inputFile','Name of input FID'),
      'head': ('headerSize','Header size in words'),
      'integer': ('integer','Data is integer (not floating point)'),
      'swap': ('swap','Swap data'),
      'big_endian': ('bigEndian','Data is big endian'),
      'varian': ('varian','Varian data'),      
      'dim': ('dim','Dimension'),
      'npts': ('numPoints','Number of real points'),
      'block': ('blockSize','Block size'),
      'sw': ('spectralWidth','Spectral width'),
      'sf': ('spectrometerFreq','Spectrometer frequence'),
      'refppm': ('refPpm','Reference ppm'),
      'refpt': ('refPoint','Reference point'),
      'nuc': ('nucleus','Nucleus'),
      'dimType': ('dimType','Dimension type'),
      'sampledValues': ('sampledValues','Sampled values'),
  }

  def initGeneric(self):  

    self.aPars = {}
    self.fPars = {}

    self.allParKeys = []
    self.parKeys = {}
    
    for parType in self.allPars.keys():
      parKeys = []
      for parInfo in self.allPars[parType]:
        parKeys.append(parInfo[0])
        if parType in ('mainPars','dimPar') and parInfo[1]:
          self.setValue(parInfo,parInfo[3],initialize = 1)
        
      self.allParKeys.extend(parKeys)
      self.parKeys[parType] = parKeys

  def initDims(self):
  
    if self.aPars.has_key('ndim') and self.aPars['ndim']:
      numDim = self.aPars['ndim']
    else:
      numDim = self.fPars['numDim']
      
    self.aPars['dimType'] = numDim * [self.dimTypeFrequency]

    for parInfo in self.allPars['dimPars']:
      
      self.setValue(parInfo,numDim * [None])

      for i in range(numDim):
        self.setListValue(parInfo,i,parInfo[3],initialize = 1)

  def setValue(self,parInfo,value,initialize = 0):
  
    if parInfo[1]:
      self.aPars[parInfo[1]] = value
    
    
    if not initialize and self.parNameMapping.has_key(parInfo[1]):
      if type(value) == type([]):
        value = value[:]
      self.fPars[self.parNameMapping[parInfo[1]][0]] = value
    
    return value
    
  def setListValue(self,parInfo,index,value, initialize = 0):
    
    
    if parInfo[1]:
      self.aPars[parInfo[1]][index] = value
      
    if not initialize and self.parNameMapping.has_key(parInfo[1]):
      # Reset to default 
      # - necessary for AXNUC/NUC1 settings from processingParsIO.py
      if value is None and parInfo[-1]:
        value = parInfo[-1]
        
      self.fPars[self.parNameMapping[parInfo[1]][0]][index] = value
    
    return value

  def findParInfo(self, pars, key):

    for parInfo in pars:
      if parInfo[1] == key:
        return parInfo

    return None

#
# Generic chemcomps, only for use in DataFormat.py at the moment (10/2009). Wim
#

class GenericChemComp:

  def __init__(self,name,molType,atomNames):

    self.name = name
    self.molType = molType    
    self.atomNames = atomNames

    self.atoms = []
    self.bonds = []

class GenericAtom:

  def __init__(self,serial,atomName):
  
    self.serial = serial
    self.name = atomName
    
    # Dirty determination of element symbol, could probably be done more consistently but not that relevant
    atomNameNoNumbers = ""
    
    for i in range(len(atomName)):
      
      if atomName[i] not in '0123456789':
        atomNameNoNumbers += atomName[i]
        
    #
    # Try to set this now...
    #
    
    elementSymbol = None

    if atomNameNoNumbers: 
      if atomName[0] in 'CHNOPSF':
        elementSymbol = atomName[0]
      else:
        elementSymbol = atomName
    
    self.atomType = elementSymbol
      
    #self.otherCoords = []
    
  #def addInfo(self,atomInfo):
  #  
  #  self.otherCoords.append({})
  #  
  #  for attrName in ('x','y','z'):
  #    self.otherCoords[-1][attrName] = atomInfo[attrName]
"""
class GenericBond:

  def __init__(self,parent,atomSerial1,atomSerial2,bondType,bondStatus):
  
    self.parent = parent
        
    self.atomSerial1 = atomSerial1
    self.atomSerial2 = atomSerial2
    
    self.bondType = bondType

    self.bondStatus = bondStatus
"""
