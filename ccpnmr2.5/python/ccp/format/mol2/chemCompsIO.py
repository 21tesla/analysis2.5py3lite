"""
======================COPYRIGHT/LICENSE START==========================

chemCompsIO.py: I/O for Tripos mol2 format atom/molecule information file

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


# Import generic file reader
from ccp.format.mol2.generalIO import Mol2GenericFile

#####################
# Class definitions #
#####################

class Mol2ChemCompFile(Mol2GenericFile):

  def initialize(self):
  
    self.chemComps = []

  def read(self,mol2File = None,verbose = 0):

    if not mol2File:
      self.readGeneric()
    else:
      self.infoDict = mol2File.infoDict
    
    if verbose == 1:
      print "Getting chemComp information."
    
    totalAtoms = len(self.infoDict['ATOM'])
    totalBonds = len(self.infoDict['BOND'])
    
    #
    # Get info for the chemComps that are present in this file.
    # If not given, assume one chemComp only and set default info.
    #
    
    if self.infoDict.has_key('SUBSTRUCTURE'):
      substInfoList = self.infoDict['SUBSTRUCTURE']
    else:
      substInfo = {
          
          'substId':    1,
          'substName':  self.infoDict['ATOM'][0]['substName'],
          'rootAtomId': 1,
          'chain':      'A',
          
                  }
      
      substInfoList = [substInfo]

    #
    # Loop over the different chemComps in the file
    #
    
    for substInfo in substInfoList:
      
      substId = substInfo['substId']
      substName = substInfo['substName']
      rootAtomId = substInfo['rootAtomId']
      
      #
      # Find the atoms
      #
      
      (substAtoms,atomNames,atomIds) = self.findSubstAtoms(totalAtoms,rootAtomId,substName)
    
      if not substAtoms:
        for attrType in ('chain','chainSubType'):
          attrName = substInfo['chain']
          if verbose:
            print "  Warning: no atoms found for substructure name '%s' - trying %s name '%s'" % (substName,attrType,attrName)
          (substAtoms,atomNames,atomIds) = self.findSubstAtoms(totalAtoms,rootAtomId,attrName)
          
          if substAtoms:
            break
         
      #
      # Find the bonds
      #
      
      substBonds = []
      bondTypes = []

      for i in range(totalBonds):
        bondInfo = self.infoDict['BOND'][i]

        if bondInfo['atom1Id'] in atomIds.keys() and bondInfo['atom2Id'] in atomIds.keys():
          substBonds.append(bondInfo)
          bondTypes.append(bondInfo['bondType'])
          
      #
      # Check if already exists. It is possible that the same name occurs with different
      # atom names!
      #
      # TODO: check if this is valid!
      #
      
      atomNames.sort()
      bondTypes.sort()
      
      chemComp = None
      onlyNameMatch = False
      
      for tChemComp in self.chemComps:
        if tChemComp.name == substName:
          onlyNameMatch = True
          if tChemComp.atomNames == atomNames and tChemComp.bondTypes == bondTypes:
            chemComp = tChemComp
            onlyNameMatch = False
            break

      if onlyNameMatch:
        print "  Warning: %s substructure %s appears with different atom names and/or bond types." % (self.format,substName)
          
      if not chemComp:
        chemComp = Mol2ChemComp(self,substName,atomNames,bondTypes,substInfo)
        self.chemComps.append(chemComp)
        
        for atomInfo in substAtoms:

          self.chemComps[-1].atoms.append(Mol2Atom(self.chemComps[-1],atomInfo['atomId'],atomInfo['atomName'],atomInfo))

        for bondInfo in substBonds:
          
          if bondInfo.has_key('bondStatus'):      
            bondStatus = bondInfo['bondStatus']  
          else:
            bondStatus = None
            
          self.chemComps[-1].bonds.append(Mol2Bond(self.chemComps[-1],bondInfo['atom1Id'],bondInfo['atom2Id'],bondInfo['bondType'],bondInfo))
      
      else:
        for atomInfo in substAtoms:
          for atom in chemComp.atoms:
            if atom.name == atomInfo['atomName']:
              atom.addInfo(atomInfo)

  def write(self,verbose = 0):
    
    if verbose == 1:
      print "Writing chemComp information to %s file %s." % (self.format,self.name)
    
    #
    # Only write one chemComp at a time...
    #
    
    if len(self.chemComps) > 1:
      print "  Warning: can only write out one chemComp at a time (multiple ones given)."
    
    chemComp = self.chemComps[0]
    
    #
    # Open output file
    #
    
    fout = open(self.name,'w')
    
    #
    # Write out header
    #
    
    fout.write('#\n')
    fout.write('# Mol2 file written by ccpNmr FormatConverter.\n')
    fout.write('#\n')
    
    if chemComp.overallCharge != None:
      fout.write('#\n')
      fout.write('# Overall charge: %d\n' % chemComp.overallCharge)
      fout.write('#\n')
      
    
    fout.write('@<TRIPOS>MOLECULE\n')
    fout.write('%s\n' % chemComp.name)
    fout.write('%5d %5d %5d %5d %5d\n' % (len(chemComp.atoms),len(chemComp.bonds),1,0,0))
    fout.write('%s\n' % chemComp.molType)
    fout.write('%s\n' % 'USER_CHARGES') # TODO: make this settable from reference data info?
    fout.write('\n\n') 

    #
    # Write atoms
    #

    fout.write('@<TRIPOS>ATOM\n')
    
    for atom in chemComp.atoms:
      
      fout.write("%7d %-8s %9.4f %9.4f %9.4f %-4s %6d %-7s %9.4f %s\n" % (
      
                   atom.serial,
                   atom.name,
                   atom.x,
                   atom.y,
                   atom.z,
                   atom.atomType,
                   atom.substId,  # This should always be for single chemComp...
                   atom.substName,
                   atom.charge,
                   atom.atomStatus
      
                  ))

    #
    # Write bonds
    #

    fout.write('@<TRIPOS>BOND\n')
    
    for i in range(0,len(chemComp.bonds)):

      bond = chemComp.bonds[i]
      
      fout.write("%5d %5d %5d %s\n" % (
      
                  (i+1),
                  bond.atomSerial1,
                  bond.atomSerial2,
                  bond.bondType
                  
                  ))
          
class Mol2ChemComp:

  def __init__(self,parent,name,atomNames,bondTypes,substInfo,molType = 'SMALL',overallCharge = None):

    self.parent = parent

    self.name = name
    self.molType = molType
    
    self.atomNames = atomNames
    self.bondTypes = bondTypes
    
    self.overallCharge = overallCharge
    
    # TODO: this is not entirely correct - what if writing this file?
    for attrName in substInfo.keys():
      if attrName not in ('substId','substName','rootAtomId'):
        setattr(self,attrName,substInfo[attrName])

    self.atoms = []
    self.bonds = []

class Mol2Atom:

  def __init__(self,parent,atomSerial,atomName,atomInfo):
  
    self.parent = parent
    
    self.serial = atomSerial
    self.name = atomName
    
    for attrName in atomInfo.keys():
      if attrName not in ('atomId','atomName'):
        setattr(self,attrName,atomInfo[attrName])
        
    self.otherCoords = []
    
  def addInfo(self,atomInfo):
    
    self.otherCoords.append({})
    
    for attrName in ('x','y','z'):
      self.otherCoords[-1][attrName] = atomInfo[attrName]

class Mol2Bond:

  def __init__(self,parent,atomSerial1,atomSerial2,bondType,bondStatus):
  
    self.parent = parent
        
    self.atomSerial1 = atomSerial1
    self.atomSerial2 = atomSerial2
    
    self.bondType = bondType

    self.bondStatus = bondStatus

