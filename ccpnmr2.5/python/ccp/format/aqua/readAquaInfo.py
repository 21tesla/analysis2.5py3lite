#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

readAquaInfo.py: Reads information from Aqua library files

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

import os, string, re

from ccp.general.Io import getDataPath

from memops.universal.Io import joinPath


patt = {}
patt['dollar'] = re.compile("^\s*\$")

xplorCodes = ['%','*','+']

#############
# Functions #
#############

def readAquaFormats(file,format,chemCompCodes):

  fin = open(file, 'rU')

  # Read in reference information (quick and dirty)
  line = fin.readline()
  
  while line:
    cols = line.split()
    
    if len(cols) == 0 or patt['dollar'].search(line):
      pass
    
    #
    # Using - and + codes as well
    #

    elif cols[0] == 'def':
      chemCompCode = cols[1]
        
      formatAtomName = cols[3]
      aquaAtomName = cols[4]
      
      if not chemCompCodes.has_key(chemCompCode):
        chemCompCodes[chemCompCode] = {}
      
      if not chemCompCodes[chemCompCode].has_key(aquaAtomName):
        chemCompCodes[chemCompCode][aquaAtomName] = {}
      
      if not chemCompCodes[chemCompCode][aquaAtomName].has_key(format):
        chemCompCodes[chemCompCode][aquaAtomName][format] = [formatAtomName]
        
        # HACK FOR XPLOR ONLY!
        # For aqua xplor names ending in hash '#',
        # will add same name but with '*','%','+' (see xplorCodes)
        # These should go in altSysNames in the data model
        if format == 'xplor':
          findHash = string.find(formatAtomName,'#')
          if findHash == len(formatAtomName) - 1:
            for xplorCode in xplorCodes:
              newName = string.replace(formatAtomName,'#',xplorCode)
              chemCompCodes[chemCompCode][aquaAtomName][format].append(newName)
              
            # For 'deep' matches, add ++, **, %% and ##
            if chemCompCodes[chemCompCode][aquaAtomName].has_key('diana'):
              if chemCompCodes[chemCompCode][aquaAtomName]['diana'][0][:2] == 'QQ':
                for xplorCode in xplorCodes + ['#']:
                  newName = string.replace(formatAtomName,'#',xplorCode)
                  newName += xplorCode
                  chemCompCodes[chemCompCode][aquaAtomName][format].append(newName)
              
      else:
        if formatAtomName not in chemCompCodes[chemCompCode][aquaAtomName][format]:
          chemCompCodes[chemCompCode][aquaAtomName][format].append(formatAtomName)

    line = fin.readline()

  return chemCompCodes
   
def readAquaDefs(file,info,chemCompCodes):
  
  fin = open(file, 'rU')
  
  # Read in reference information (quick and dirty)
  line = fin.readline()
  
  while line:
    cols = line.split()
 
    if len(cols) == 0 or patt['dollar'].search(line):
      pass
    
    #
    # Using - and + codes as well
    #
    
    elif cols[0] == 'def':
      chemCompCode = cols[1]

      aquaAtomName = cols[3]
      aquaCorrespondingAtoms = cols[5:]
      
      chemCompCodes[chemCompCode][aquaAtomName][info] = aquaCorrespondingAtoms
    
    line = fin.readline()

  return chemCompCodes

def getAquaAtomInfo(chemCompCodes = {}):

  aquaRefDir = getDataPath('ccp', 'aqua')
  files = os.listdir(aquaRefDir)

  for file in files:
    if file[0:8] == 'AtomLIB-':
      format = file[8:]
      
      readAquaFormats(joinPath(aquaRefDir,file),format,chemCompCodes)
  
  for file in files:
    if file[-8:] == 'LIB-aqua':
      info = file[0:-8]
      readAquaDefs(joinPath(aquaRefDir,file),info,chemCompCodes)
  
  return chemCompCodes
  
###################
# Main of program #
###################

if __name__ == "__main__":  

  chemCompCodes = getAquaAtomInfo()

  
  chemCompCodeList = chemCompCodes.keys()
  chemCompCodeList.sort()

  for chemCompCode in chemCompCodeList:
    print chemCompCode

    atomNameList = chemCompCodes[chemCompCode].keys()
    atomNameList.sort()

    for aquaAtomName in atomNameList:
      
      print "  " + aquaAtomName
      
      formatList = chemCompCodes[chemCompCode][aquaAtomName].keys()
      formatList.sort()
      
      for format in formatList:
        print "  " + str(chemCompCodes[chemCompCode][aquaAtomName][format]) + " (%s)" % format
