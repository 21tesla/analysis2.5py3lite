"""
======================COPYRIGHT/LICENSE START==========================

CosmosFormat.py: Contains functions specific to Cosmos conversions.

Copyright (C) 2012 Wim Vranken (Vrije Universiteit Brussel)

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

import os, traceback, sys, copy

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccpnmr.format.general.Util import getNameInfo, getResName

from ccpnmr.format.general.Constants import allResidueAtoms_kw

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['writeProject']['models'] = (None,False,'The structure models with coordinates to be exported.')
IOkeywords['writeProject']['shiftList'] = (None,False,'The chemical shift list for the exported file.')
IOkeywords['writeProject']['constraintList'] = (None,False,'The distance constraint list for the exported file.')
IOkeywords['writeProject']['projectType'] = ('COD',True,'The type of project file. Use COO for coordinates/shifts, COD for constraints/shifts.')

class CosmosFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'cosmos'
    
    # Necessary for project COO write
    self.atomToCoordinateSerial = {}
    
    # Necessary for project COD write
    #self.resonancesDoneForShiftValue = {}

  def setGenericImports(self):
    
    #self.getSequence = self.getSequenceGeneric
    #self.createSequenceFile = self.createSequenceFileGeneric

    self.createCoordinateFile = self.createCoordinateFileGeneric
    
    self.createMeasurementFile = self.createMeasurementFileGeneric
    
    self.createConstraintFile = self.createConstraintFileGeneric

  #
  # Deviations from generic import stuff
  #

  def createFullProject(self,fileName,models=None,shiftList=None,constraintList=None,projectType='COO'):

    self.projectType = projectType

    if self.verbose == 1:
      print "Writing %s %s project from file %s" % (self.formatLabel,self.projectType,fileName)
   
    self.file = self.projectIO.CosmosProjectFile(fileName)
   
    if self.projectType == 'COO':
      self.writeCoordinates(fileName, structures=models, noWrite=True, resetIOkeywords=False, useCcpnChainInfo=True)
      self.file.coordinateFile = self.coordinateFile
      
    if self.projectType == 'COD' and constraintList:
      self.writeDistanceConstraints(fileName,constraintList=constraintList,noWrite=True,resetIOkeywords=False,resetMapping=False,individualAtoms=True,compressResonances=False,useCcpnChainInfo=True,retainOriginalCompression=False)
      self.file.constraintFile = self.constraintFile
    
    if shiftList:
      self.writeShifts(fileName, measurementList=shiftList, noWrite=True, resetMapping=False, resetIOkeywords=False, individualAtoms=True, useCcpnChainInfo=True)
      self.file.chemShiftFile = self.measurementFile
    
    self.file.write(fileType=self.projectType)
    
    return self.file
    
  #
  # Functions different to default functions in DataFormat
  #

  def setRawCoordinate(self):
  
    if not self.coordinateFile.modelCoordinates.has_key(self.modelId):

      self.coordinateFile.modelCoordinates[self.modelId] = []
    
    chemCompVar = self.residue.chemCompVar
    chemComp = chemCompVar.chemComp
    
    resName = chemComp.ccpCode.upper()
    atomName = self.atomName
    
    atomNumber = self.atom.chemAtom.chemElement.atomNumber
 
    modelCoordinate = self.coordinatesIO.CosmosCoordinate(self.coordinateSerial,
                                                          self.chain.code,
                                                          self.residue.seqId,        
                                                          atomName,
                                                          atomNumber,
                                                          resName,
                                                          self.x,
                                                          self.y,
                                                          self.z)
      
    self.coordinateFile.modelCoordinates[self.modelId].append(modelCoordinate)
    
    if not self.atomToCoordinateSerial.has_key(self.chain):
      self.atomToCoordinateSerial[self.chain] = {}
    if not self.atomToCoordinateSerial[self.chain].has_key(self.residue):
      self.atomToCoordinateSerial[self.chain][self.residue] = {}
    if not self.atomToCoordinateSerial[self.chain][self.residue].has_key(self.atom):
      self.atomToCoordinateSerial[self.chain][self.residue][self.atom] = []
      
    self.atomToCoordinateSerial[self.chain][self.residue][self.atom].append(self.coordinateSerial)
    
  def setChemShiftFileValue(self):
    
    if self.atomName in self.atomNamesDict.keys():
    
      resonanceToAtom = self.atomNamesDict[self.atomName]
      measurement = self.origAtomMeasurements[resonanceToAtom]
      
      residueLabel = self.getResidueLabelThreeLetter(resonanceToAtom)
      
      coordAtomSerial = None
      atom = self.residue.findFirstAtom(name = self.atomName)
      
      if hasattr(self,'atomToCoordinateSerial') and self.atomToCoordinateSerial.has_key(self.chain) and \
         self.atomToCoordinateSerial[self.chain].has_key(self.residue) and \
         self.atomToCoordinateSerial[self.chain][self.residue].has_key(atom):
         
        coordAtomSerial = self.atomToCoordinateSerial[self.chain][self.residue][atom][0]
      
      #
      # Set the value
      #
      
      self.measurementFileValues.append(
      
               self.rawMeasurementClass(measurement.value,
                                        atomSerial = coordAtomSerial
                                        ))
      
      self.measurementFileValues[-1].addAssignment(self.atomName,residueLabel.upper(),resonanceToAtom.chain.code,resonanceToAtom.seqId)
      
      #
      # Add other assignments for COD file
      #
      
      if self.projectType == 'COD':  
        for resonance in [resonanceToAtom.resonance] + resonanceToAtom.otherLinkedResonances + resonanceToAtom.otherGroupResonances:
          for otherResonanceToAtom in self.atomMeasurements.keys():
            if otherResonanceToAtom == resonanceToAtom:
              # Don't write out own info again!
              continue
            if resonance == otherResonanceToAtom.resonance:
              otherAtomName = otherResonanceToAtom.atomName
              self.measurementFileValues[-1].addAssignment(otherAtomName, self.getResidueLabelThreeLetter(otherResonanceToAtom).upper(),otherResonanceToAtom.chain.code,otherResonanceToAtom.seqId)
              
              # This might be a problem?
              if resonance == resonanceToAtom.resonance and otherAtomName in self.atomNamesDict.keys():
                del(self.atomNamesDict[otherAtomName])
                
      #
      # Remove self anyway... might be a problem for ambiguous assignments
      #
      
      elif self.projectType == 'COO':
        if resonanceToAtom in self.atomMeasurements.keys():
          del(self.atomMeasurements[resonanceToAtom])
              
  def setRawDistanceConstraint(self):
  
    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraint.serial))

    self.rawConstraint = self.constraintFile.constraints[-1]
 
    self.rawConstraint.setDistanceData(self.constraint.upperLimit,self.constraint.targetValue)
              
  def setRawDistanceConstraintItem(self):
    
    self.rawConstraint.items.append(self.rawConstraintItemClass())
    self.rawConstraintItem = self.rawConstraint.items[-1]
    
  def setRawDistanceConstraintItemMembers(self):
    
    for i in range(2):
            
      resonanceToAtom = self.ccpInfo[i][3]
      residue = resonanceToAtom.getResidue()

      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])
      
      if atomName == allResidueAtoms_kw:
        atomName = None
      
      # Using residue.seqCode - OK since useCcpnChainInfo = True for all exports
      self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass(chainCode,residue.seqId,atomName,residue.ccpCode.upper()))

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)
