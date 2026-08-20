
"""
======================COPYRIGHT/LICENSE START==========================

generalIO.py: General I/O information for Tripos Mol2 format files

Copyright (C) 2006-2009 Wim Vranken (European Bioinformatics Institute)

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

# Import general functions
from memops.universal.Util import returnInt, returnFloat

from ccp.format.general.formatIO import FormatFile
from ccp.format.general.Constants import defaultMolCode

#####################
# Class definitions #
#####################

class Mol2GenericFile(FormatFile):

  def setGeneric(self):
    
    self.format = 'mol2'
    self.defaultMolCode = defaultMolCode
    
    self.readDict = {
    
      'MOLECULE': {
           1: ('molName',),
           2: ('numAtoms','numBonds','numSubst','numFeat','numSets'),
           3: ('molType',),# One of SMALL, BIOPOLYMER, PROTEIN, NUCLEIC_ACID, SACCHARIDE
           4: ('chargeType',),# NO_CHARGES, DEL_RE, GASTEIGER, GAST_HUCK, HUCKEL, PULLMAN, GAUSS80_CHARGES, AMPAC_CHARGES, MULLIKEN_CHARGES, DICT_ CHARGES, MMFF94_CHARGES, USER_CHARGES
           5: ('molStatus',),
           6: ('molComment',)
                   },
                   
      'ATOM': {
      
           0: ('atomId',returnInt),
           1: ('atomName',string.strip),
           2: ('x',returnFloat),
           3: ('y',returnFloat),
           4: ('z',returnFloat),
           5: ('atomType',self.getAtomType),
           6: ('substId',returnInt),       # Substructure ID
           7: ('substName',string.strip),
           8: ('charge',returnFloat),
           9: ('atomStatus',string.strip)  #DSPMOD, TYPECOL, CAP, BACKBONE, DICT, ESSENTIAL, WATER and DIRECT.

                },

      'BOND': {
      
           0: ('bondId',returnInt),
           1: ('atom1Id',returnInt),
           2: ('atom2Id',returnInt),
           3: ('bondType',string.strip),  #1 = single, 2 = double, 3 = triple, am = amide, ar = aromatic, du = dummy, un = unknown, nc = not connected
           4: ('bondStatus',string.strip) #TYPECOL, GROUP, CAP, BACKBONE, DICT and INTERRES
           
               },

      'SUBSTRUCTURE': {
      
           0: ('substId',returnInt),
           1: ('substName',string.strip),
           2: ('rootAtomId',returnInt),
           3: ('substType',string.strip), # temp, perm, residue, group or domain
           4: ('dictType',returnInt),     # dictionary associated with structure?
           5: ('chain',string.strip),  
           6: ('chainSubType',string.strip),   # Chain subtype
           7: ('interBonds',returnInt),        # The number of inter substructure bonds
           8: ('substStatus',string.strip),    #LEAF, ROOT, TYPECOL, DICT, BACKWARD and BLOCK
           9: ('comment',string.strip),  
           
               },
    }
    
    self.infoDict = {}

  def initialize(self):
    
    pass
    
  def readGeneric(self,verbose = 0):

    if verbose == 1:
      print "Reading %s file %s" % (self.format,self.name)
      
    #
    # Read all information...
    #

    fin = open(self.name, 'rU')
    line = fin.readline()
    
    triposSection = None

    while line:
    
      if self.patt['emptyline'].search(line):
        line = fin.readline()
        continue
    
      triposTagSearch = self.patt[self.format + 'TriposTag'].search(line)
      
      if triposTagSearch:
        triposSection = triposTagSearch.group(1).strip()      
        sectionLine = 0
        if not self.infoDict.has_key(triposSection):
          self.infoDict[triposSection] = {}
      
      elif triposSection:
      
        sectionLine += 1

        if triposSection == 'MOLECULE':
        
          if sectionLine in (1,3,4,5,6):
            infoName = self.readDict[triposSection][sectionLine][0]
            self.infoDict[triposSection][infoName] = line.strip()
          
          elif sectionLine == 2:
            numbers = line.split()
            for i in range(len(numbers)):
              infoName = self.readDict[triposSection][sectionLine][i]
              self.infoDict[triposSection][infoName] = returnInt(numbers[i])
           
        elif triposSection in ('ATOM','BOND','SUBSTRUCTURE'):
        
          if not self.infoDict[triposSection]:
            self.infoDict[triposSection] = []
          
          self.infoDict[triposSection].append({})

          atomCols = line.split()
          for i in range(len(atomCols)):
            (infoName,returnFunc) = self.readDict[triposSection][i]
            self.infoDict[triposSection][-1][infoName] = returnFunc(atomCols[i])

      line = fin.readline()
      
    fin.close()
  
  #
  # Specific functions
  #
  
  def getAtomType(self,atomTypeInfo):
    
    # Can be of type H.ar
    atomTypeData = atomTypeInfo.split('.')
    
    return atomTypeData[0]

  def findSubstAtoms(self,totalAtoms,rootAtomId,substName):

    atomNames = []
    atomIds = {}
    substAtoms = []

    for i in range(totalAtoms):
      
      atomInfo = self.infoDict['ATOM'][i]

      if atomInfo['atomId'] == rootAtomId and atomInfo['substName'] == substName:
      
        #
        # Here determine whether all atom names the same
        # If so, renumber as reorganising the atom info!!
        #
      
        addAtomIndex = {}        
        atomNamesByElement = {}
        
        iTmp = i
        while (atomInfo['substName'] == substName):
          
          atomType = atomInfo['atomType']
          atomName = atomInfo['atomName']
          
          if not atomNamesByElement.has_key(atomType):
            atomNamesByElement[atomType] = {}
          if not atomNamesByElement[atomType].has_key(atomName):
            atomNamesByElement[atomType][atomName] = 0
            
          atomNamesByElement[atomType][atomName] += 1
        
          iTmp += 1
          if iTmp < totalAtoms:
            atomInfo = self.infoDict['ATOM'][iTmp]
          else:
            break
            
        for atomType in atomNamesByElement.keys():
          addAtomIndex[atomType] = 0
          if len(atomNamesByElement[atomType]) == 1:
            atomName = atomNamesByElement[atomType].keys()[0]
            if atomNamesByElement[atomType][atomName] > 1:
              addAtomIndex[atomType] = 1
                      
        #
        # Reset before starting data reorganisation
        #
        
        atomInfo = self.infoDict['ATOM'][i]

        while (atomInfo['substName'] == substName):
          if addAtomIndex[atomInfo['atomType']]:
            newAtomName = "%s%d" % (atomInfo['atomName'],addAtomIndex[atomInfo['atomType']])
            atomInfo['atomName'] = newAtomName
            addAtomIndex[atomInfo['atomType']] += 1
          elif atomInfo['atomName'] in atomNames:
            for i in range(10000):
              newAtomName = "%s_%d" % (atomInfo['atomName'],i)
              if newAtomName not in atomNames:
                break
          else:
            newAtomName = atomInfo['atomName']
           
          substAtoms.append(atomInfo)
          atomNames.append(newAtomName)
          atomIds[atomInfo['atomId']] = atomInfo
          i += 1

          if i < totalAtoms:
            atomInfo = self.infoDict['ATOM'][i]
          else:
            break

        break
    
    return (substAtoms,atomNames,atomIds)
