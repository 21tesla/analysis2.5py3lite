#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

cyanaLibParser.py: Parser for Cyana reference information

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

import string

from memops.universal.Util import returnInt, returnInts
from memops.universal.Io import joinPath

from ccp.format.general.Util import getRegularExpressions
from ccp.general.Io import getDataPath

#
# This parser reads the CYANA reference information so it can be used
# for DIANA/DYANA/CYANA file parsing.
#
# The CYANA reference information file is the standard library file from the
# CYANA distribution.
#

#####################
# Class definitions #
#####################

class CyanaLibrary:

  def __init__(self, version = None):
  
    self.patt = getRegularExpressions()
  
    self.readComponents = ['STARTGRP','HEADRES','RESIDUE','TAILRES']
  
    self.refInfo = {}
    
    if version in ['2.1']:
      fileName = 'cyana2.1.lib'
    else:
      fileName = 'cyana.lib'
      
    self.read(getDataPath('ccp', 'cyana', fileName))
    
    self.setAtomPosition()
    
    self.setAngleAtoms()
  
  def read(self,file):
  
    fin = open(file, 'rU')
    
    name = None
  
    line = fin.readline()
  
    while line:
    
      cols = string.split(line)
  
      if not cols:
      
        pass
  
      elif cols[0] in self.readComponents:
        
        infoType = cols[0]
        name = cols[1]
        
        if not self.refInfo.has_key(infoType):
        
          self.refInfo[infoType] = []
        
        cyanaRefInfo = CyanaRefInfo(name)
        
        self.refInfo[infoType].append(cyanaRefInfo)
      
      elif name and self.patt['onlydigit'].search(cols[0]):
        
        if len(cols) == 10:
        
          #
          # Angle info
          #
          
          cyanaRefInfo.setAngle(cols[0],cols[1],cols[5:9])
          
        elif len(cols) == 13:
        
          #
          # Atom info
          #
          
          cyanaRefInfo.setAtom(cols[0],cols[1],cols[2],cols[8:12])
          
      else:
        name = None
        
        
      line = fin.readline()
  
  def setAngleAtoms(self):
    
    for refInfo in self.refInfo.keys():
    
      for cyanaRefInfo in self.refInfo[refInfo]:
      
        for angle in cyanaRefInfo.angles:
        
          for atomSerial in angle.atomSerials:
          
            for atom in cyanaRefInfo.atoms:
              
              atomFound = 0
            
              if atom.serial == atomSerial:
              
                angle.atomNames.append(atom.name)
                angle.atoms.append(atom)
                
                atomFound = 1
                break
              
            if not atomFound:
              
              print "  No match for atom serial %s, angle %s (%s)" % (atomSerial,angle.name,cyanaRefInfo.name)

  def setAtomPosition(self):
  
    preAtoms = ['C','O',"C3'","O3'"]
    postAtoms = ['N','P']
    
    #
    # Determines the position of (backbone) atoms (-1,0,1)
    #
    
    for refInfo in self.refInfo.keys():
    
      for cyanaRefInfo in self.refInfo[refInfo]:
      
        for i in range(0,len(cyanaRefInfo.atoms)-1):
        
          for j in range(i+1,len(cyanaRefInfo.atoms)):
          
            if cyanaRefInfo.atoms[i].name == cyanaRefInfo.atoms[j].name:
            
              #
              # Set the location for these (i or j depending on case)
              #
              
              if cyanaRefInfo.atoms[i].name in preAtoms:
              
                cyanaRefInfo.atoms[i].location = -1
              
              elif cyanaRefInfo.atoms[j].name in postAtoms:
              
                cyanaRefInfo.atoms[j].location = 1
                
              else:
                pass
                #print "No location can be set for ", cyanaRefInfo.name, cyanaRefInfo.atoms[j].name

  def getAtomOrderDict(self,verbose = 0):
  
    atomOrderDict = {}
    
    for cyanaRefInfo in self.refInfo['RESIDUE']:
    
      cyanaCode = string.strip(cyanaRefInfo.name)
      
      atomOrderDict[cyanaCode] = []
      
      for atom in cyanaRefInfo.atoms[2:-2]:
      
        #
        # Ignore oxygens
        #
        
        #if atom.name[0] in ['O']:
        
        #  continue
        
        atomOrderDict[cyanaCode].append(atom.name)
    
    if verbose:
      resNames = atomOrderDict.keys()
      resNames.sort()
      for resName in resNames:
        print resName
        print "   ",atomOrderDict[resName]
        print
        
    return atomOrderDict

  def show(self):
  
    for refInfo in self.refInfo.keys():
    
      print refInfo
    
      for cyanaRefInfo in self.refInfo[refInfo]:
      
        print "  ",cyanaRefInfo.name
      
        if cyanaRefInfo.angles:
          
          print
          print "    Angles:"

          for angle in cyanaRefInfo.angles:
            
            print "      ",angle.name,angle.atomNames
        
      
        if cyanaRefInfo.atoms:
          
          print
          print "    Atoms:"
          
          text =  "      "

          for atom in cyanaRefInfo.atoms:
            
            text += "%s (%s)," % (atom.name,atom.location)
            
          print text[:-1]

      print
  
  def findResLabel(self,resLabel):
    
   for refInfo in self.refInfo.keys():
      for cyanaRefInfo in self.refInfo[refInfo]:
        if cyanaRefInfo.name == resLabel:
          return True
    
   return False
      
  def findAngle(self,resLabel,angleName):
  
    for refInfo in self.refInfo.keys():
    
      for cyanaRefInfo in self.refInfo[refInfo]:
      
        if cyanaRefInfo.name == resLabel:
      
          for angle in cyanaRefInfo.angles:
          
            if angle.name == angleName:
            
              return angle

    return None

  def findAngleFromAtoms(self,resLabel,atomNames):
  
    reversedAtomNames = atomNames[:]
    reversedAtomNames.reverse()
  
    for refInfo in self.refInfo.keys():
    
      for cyanaRefInfo in self.refInfo[refInfo]:
      
        if cyanaRefInfo.name == resLabel:
      
          for angle in cyanaRefInfo.angles:
          
            if angle.atomNames == atomNames or angle.atomNames == reversedAtomNames:
            
              return angle
    
    return None

      
  def findAtom(self,resLabel,atomName):
  
    for refInfo in self.refInfo.keys():
    
      for cyanaRefInfo in self.refInfo[refInfo]:
      
        if cyanaRefInfo.name == resLabel:
      
          for atom in cyanaRefInfo.atoms:
          
            if atom.name == atomName:
            
              return atom
              
    return None

  def findAtomBySerial(self,resLabel,atomSerial):
  
    for refInfo in self.refInfo.keys():
    
      for cyanaRefInfo in self.refInfo[refInfo]:
      
        if cyanaRefInfo.name == resLabel:
      
          for atom in cyanaRefInfo.atoms:
          
            if atom.serial == atomSerial:
            
              return atom
              
    return None

class CyanaRefInfo:

  def __init__(self,name):
  
    self.name = name

    self.angles = []
    self.atoms = []
    
  def setAngle(self,serial,name,atomSerials):
  
    self.angles.append(CyanaAngle(serial,name,atomSerials))
    
  def setAtom(self,serial,name,atomType,bondedAtomSerials):
  
    self.atoms.append(CyanaAtom(serial,name,atomType,bondedAtomSerials))
  
  def __cmp__(self,other):
  
    if not isinstance(other, CyanaRefInfo):
      return cmp(id(self), other)

    if self.name < other.name:
      return -1
    elif self.name == other.name:
      return 0
    elif self.name > other.name:
      return 1
      
class CyanaAngle:

  def __init__(self,serial,name,atomSerials):
  
    self.serial = returnInt(serial)
    self.name = name
    self.atomSerials = returnInts(atomSerials)
    self.atomNames = []
    self.atoms = []

class CyanaAtom:

  def __init__(self,serial,name,atomType,bondedAtomSerials):
  
    self.serial = returnInt(serial)
    self.name = name
    self.atomType = atomType
    self.location = 0
    
    self.bondedAtomSerials = returnInts(bondedAtomSerials)
    
########
# Main #
########

if __name__ == "__main__":

  cyanaLib = CyanaLibrary()
  
  cyanaLib.show()
  cyanaLib.getAtomOrderDict(verbose=1)
  #print cyanaLib.findAtom('ASN','CA')
  #print cyanaLib.findAtomBySerial('ASN',15)
