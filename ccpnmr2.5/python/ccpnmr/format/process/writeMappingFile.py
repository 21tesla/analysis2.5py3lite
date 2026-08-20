#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

writeMappingFile.py: Writes a resonance-atom mapping file

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

import os, sys, time, string

from memops.universal.Util import drawBox
from memops.universal.Util import returnInt


from ccpnmr.format.general.Util import getApplResNames
from ccpnmr.format.general.Util import getNameInfo

from ccpnmr.format.general.Constants import newline, allFormatsDict
#
# Write out ccpn mapping per original file format
#

def writeMappingFile(project,outFile,originalFormats = None):
  
  if originalFormats:
    checkFormatList = originalFormats
    
  else:
    checkFormatList = allFormatsDict.keys()
  
  #
  # First check if there's any mapping info at all...
  #
  
  formatList = []   
  allResonances = project.currentNmrProject.sortedResonances()
 
  for format in checkFormatList:
    resonancesDict = getApplResNames(format,allResonances)
    
    if resonancesDict != {}:
      formatList.append(format)


  if not formatList:
  
    #
    # No resonances linked
    #
    
    return False
    
  formatList.sort()


  fout = open(outFile,'w')
  
  #
  # Write top header with time/user info
  #
  
  fout.write("#" + newline)
  fout.write("# CcpNmr atom name mapping file." + newline)
  fout.write("#" + newline)
  fout.write("# Contains list of original names from the imported file(s) and the" + newline)
  fout.write("# atom names they were mapped to in the data model." + newline)
  fout.write("#" + newline)
  fout.write("# File written on %s." % time.ctime()+ newline)
  fout.write("#" + newline * 2)
    
  for format in formatList:
  
    resonancesDict = getApplResNames(format,allResonances)
    resonanceNames = []
    
    #
    # Do some sorting...
    #
    
    resonancesSortDict = {}
    
    for resonanceName in resonancesDict.keys():
      
      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resonanceName)
      
      if seqCode == None:
        continue
        
      seqCode = returnInt(seqCode)

      if not resonancesSortDict.has_key(chainCode):
        resonancesSortDict[chainCode] = {}
      
      if not resonancesSortDict[chainCode].has_key(seqCode):
        resonancesSortDict[chainCode][seqCode] = {}

      if not resonancesSortDict[chainCode][seqCode].has_key(seqInsertCode):
        resonancesSortDict[chainCode][seqCode][seqInsertCode] = {}
      
      resonancesSortDict[chainCode][seqCode][seqInsertCode][atomName] = resonanceName
    
    chainCodeList = resonancesSortDict.keys()
    chainCodeList.sort()
    
    for chainCode in chainCodeList:
    
      seqCodeList = resonancesSortDict[chainCode].keys()
      seqCodeList.sort()
    
      for seqCode in seqCodeList:
      
        seqInsertCodeList = resonancesSortDict[chainCode][seqCode].keys()
        seqInsertCodeList.sort()
        
        for seqInsertCode in seqInsertCodeList:
    
          atomNameList = resonancesSortDict[chainCode][seqCode][seqInsertCode].keys()
          atomNameList.sort()

          for atomName in atomNameList:

            resonanceNames.append(resonancesSortDict[chainCode][seqCode][seqInsertCode][atomName])
    
    #
    # Write it all out
    #
    
    if resonancesDict != {}:
    
      fout.write(drawBox("Origin file format: %s" % format,liner = "*"))
      fout.write(newline)
      
      fout.write("  %-15s    %-15s  %s" % ('Original name','Mapped to atoms','Molecules:chains:residues'))
      fout.write(newline * 2)
      
      for resonanceName in resonanceNames:
        
        resonances = resonancesDict[resonanceName]
        atomNames = []
        
        resCodes = []
        moleculeNames = []
        chainCodes = []
        
        for resonance in resonances:
          atomName = ""
          
          resonanceSet = resonance.resonanceSet
        
          if resonanceSet:
          
            for atomSet in resonanceSet.atomSets:
            
              #
              # Get some info out of data model for clarity...
              #
            
              residue = atomSet.findFirstAtom().residue
              resCode = residue.molResidue.ccpCode + str(residue.seqCode) + string.strip(residue.seqInsertCode)
              chainCode = residue.chain.code
              moleculeName = residue.chain.molecule.name
                
              if resCode not in resCodes:
                resCodes.append(resCode)
                
              if chainCode not in chainCodes:
                chainCodes.append(chainCode)
                
              if moleculeName not in moleculeNames:
                moleculeNames.append(moleculeName)
              
              #
              # Get name of atom(s) this resonance was mapped to
              #
            
              atomName += atomSet.name + "/"
          
            if atomName:
              atomName = atomName[:-1]
              
              if atomName not in atomNames:
                atomNames.append(atomName)
              
                    
        if atomNames:
        
          moleculeText = string.join(moleculeNames,',')
          chainText = string.join(chainCodes,',')
          resCodesText = string.join(resCodes,',')
          
          addText = moleculeText + ':' + chainText + ':' + resCodesText
        
          fout.write("  %-15s -> %-15s  %s" % ("'" + resonanceName + "'",string.join(atomNames,','),addText))

        else:
          fout.write("  %-15s    %-15s" % ("'" + resonanceName + "'",'NOT LINKED'))
        
        fout.write(newline)
        
        
  return True
  
###################
# Main of program #
###################

if __name__ == "__main__":

  inFile = sys.argv[1]
  outFile = sys.argv[2]
  
  # Read xml file
  from memops.format.xml import XmlIO
  proj = XmlIO.loadProjectFile(inFile)
  
  writeMappingFile(proj,outFile)
