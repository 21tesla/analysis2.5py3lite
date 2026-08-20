#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

makeChemShiftSelections.py: Makes a selection of chemical shifts and writes them out

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


from ccpnmr.format.general.Constants import newline, chemShiftClasses

# Get format classes
from ccpnmr.format.converters.NmrStarFormat import NmrStarFormat

#############
# Functions #
#############

def makeChemShiftSelections(guiParent,chains,chemShiftList,verbose = 0):
  
  #
  # TODO: only NON-EXCHANGABLE hydrogens for sidechain!!! Mark them in ref data?
  #
    
  #
  # Get the info from the chemical shift list and the relevant chains
  #
  
  if not chains or not chemShiftList:
    print "  ERROR need chains or chemShiftList!"
    return False
  
  #
  # Set up for all chains...
  #
   
  results = {}

  for chain in chains:
  
    if chain.molecule.molType not in chemShiftClasses.keys():
      continue
  
    chemShiftClassMolType = chemShiftClasses[chain.molecule.molType]
  
    for residue in chain.sortedResidues():
    
      atomSetsDone = []
    
      for atom in residue.sortedAtoms():
      
        atomType = atom.chemAtom.elementSymbol

        if atomType not in chemShiftClassMolType.keys():
          continue

        atomChemShiftClasses = chemShiftClassMolType[atomType]

        atomName = atom.name
        atomRecognized = 0
        
        for chemShiftClass in atomChemShiftClasses.keys():

          (findAtomNames,ignoreAtomNames) = atomChemShiftClasses[chemShiftClass]
          
          atomFound = 0

          if atomName in findAtomNames:

            atomFound = 1

          elif not findAtomNames:

            if atomName not in ignoreAtomNames:
            
              chemAtom = atom.chemAtom
              chemAtomSet = chemAtom.chemAtomSet
              
              if chemAtomSet and (chemAtomSet.isEquivalent or chemAtomSet.isEquivalent == None):
                
                if chemAtomSet.name not in atomSetsDone:
                  atomSetsDone.append(chemAtomSet.name)
                  atomFound = 1
                
                atomRecognized = 1
                               
              else:
                atomFound = 1
          
          if atomFound:
          
            atomRecognized = 1
          
            if not results.has_key(chemShiftClass):
              results[chemShiftClass] = [0,0]
      
            results[chemShiftClass][0] += 1


        if not atomRecognized:
          print " Chain atom name %s not found" % atomName            
        
  #
  # Do per chain... but keep track of all.
  #
  
  format = NmrStarFormat(chains[0].root,guiParent,verbose = 0)
  format.getFormatNamingSystemName = getFormatNamingSystemName
  
  for chain in chains:
  
    if chain.molecule.molType not in chemShiftClasses.keys():
      continue
  
    chemShiftClassMolType = chemShiftClasses[chain.molecule.molType]
  
    format.writeShifts('voidFileNotWritten',measurementList = chemShiftList,chains = [chain],noWrite = True, verbose = False, individualAtoms = False)
    
    for chemShift in format.measurementFile.chemShifts:
      
      atomType = chemShift.atomType
      
      if atomType not in chemShiftClassMolType.keys():
        continue
        
      atomChemShiftClasses = chemShiftClassMolType[atomType]
      atomRecognized = 0
      
      for chemShiftClass in atomChemShiftClasses.keys():
      
        (findAtomNames,ignoreAtomNames) = atomChemShiftClasses[chemShiftClass]
        
        atomFound = 0
        
        if chemShift.atomName in findAtomNames:
          
          atomFound = 1
        
        elif not findAtomNames:
          
          if chemShift.atomName not in ignoreAtomNames:
            atomFound = 1
        
        if atomFound:
          atomRecognized = 1
          results[chemShiftClass][1] += 1
          
      if not atomRecognized:
        print " Chemical shift atom name %s not found" % atomName            
            
    del format.measurementFile

  #
  # Print results
  #
  
  if verbose:

    chemShiftResultClasses = results.keys()
    chemShiftResultClasses.sort()

    for chemShiftClass in chemShiftResultClasses:
      percentage = results[chemShiftClass][1] * 100.0 / results[chemShiftClass][0]
      print "  %-20s: %6.2f found (%d/%d total)" % (chemShiftClass,percentage,results[chemShiftClass][1],results[chemShiftClass][0])

  return results

def getFormatNamingSystemName():

  return 'IUPAC'


  
