#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

CnsFormat.py: Contains functions specific to Cns/XPLOR conversions.

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

import copy

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccp.format.general.Constants import defaultSeqInsertCode

from ccpnmr.format.general.Util import getNameInfo

from ccpnmr.format.general.Constants import distanceConstraintDefaultLowerLimit
from ccpnmr.format.general.Constants import defaultMolCode

from ccpnmr.format.general.Constants import volume_kw, height_kw, allResidueAtoms_kw

#
# Add some information to IOkeywords...
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['writeCoordinates']['forceExportSegId'] = (None,False,'Force this seg ID to be set for ALL chains.')
IOkeywords['writeCoordinates']['forceExportChainId'] = (None,False,'Force this chain ID to be set for ALL chains.')
IOkeywords['writeRdcConstraints']['refFrameSegId'] = (defaultMolCode,False,'The segment Id for the reference frame.')
IOkeywords['writeRdcConstraints']['refFrameResId'] = (999,False,'The residue code for the reference frame.')
IOkeywords['writeRdcConstraints']['writeWeight'] = (False,False,'Include restraint weight in exported file.')
IOkeywords['writeRdcConstraints']['filterAtomNames'] =  (None,False,'If set (to a comma separated list of CCPN atom name strings), only those atoms will be exported. Can be used to filter out, for example, only HN-N RDCs so a different weight can be set for them.')
IOkeywords['writeCsaConstraints']['refFrameSegId'] = (defaultMolCode,False,'The segment Id for the reference frame.')
IOkeywords['writeCsaConstraints']['refFrameResId'] = (999,False,'The residue code for the reference frame.')
IOkeywords['writeCsaConstraints']['writeWeight'] = (False,False,'Include restraint weight in exported file.')
IOkeywords['writeDistanceConstraints']['oneLinePerItem'] = (False,False,'Write out as one line per constraint item.')
IOkeywords['writeDistanceConstraints']['usePeakInfo'] = (True,False,'Use any information about peaks attached to the constraints and write out in ARIA style.')
IOkeywords['writeHBondConstraints']['oneLinePerItem'] = (False,False,'Write out as one line per constraint item.')

class CnsFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'cns'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric 

    self.getCoordinates = self.getCoordinatesGeneric
    self.createCoordinateFile = self.createCoordinateFileGeneric
    
    self.getMeasurements = self.getMeasurementsGeneric
    self.createMeasurementFile = self.createMeasurementFileGeneric

    self.getConstraints = self.getConstraintsGeneric
    self.createConstraintFile = self.createConstraintFileGeneric

  #
  # Deviations from generic import stuff
  #

  #
  # Functions different to default functions in DataFormat
  #

  def duplicateChain(self):
  
    #
    # All chain info should be contained in coordinate file...
    #
    
    return 0   

  def setDistanceConstraintLimits(self):
  
    self.constraintPars['targetValue'] = self.rawConstraint.targetDist
    self.constraintPars['upperLimit'] = self.constraintPars['targetValue'] + self.rawConstraint.plusDist
    self.constraintPars['lowerLimit'] = self.constraintPars['targetValue'] - self.rawConstraint.minusDist

  def createConstraintFileFormatSpecific(self):
 
    if hasattr(self,'oneLinePerItem') and self.oneLinePerItem and self.constraintType in ('distance','hBond'):
      self.writeKeywds['oneLinePerItem'] = self.oneLinePerItem
      
    if hasattr(self,'writeWeight') and self.constraintType == 'rdc':
      self.writeKeywds['includeWeight'] = self.writeWeight

    if hasattr(self,'filterAtomNames') and self.filterAtomNames and self.constraintType == 'rdc':
      filterAtomNameList = []
      for filterAtomName in self.filterAtomNames.strip().split(','):
        filterAtomNameList.append(filterAtomName.strip())

      print "  Warning: filtering RDC restraint list on atom names %s" % str(filterAtomNameList)
      self.filterAtomNames = filterAtomNameList

  def setupJCouplingConstraintItems(self):
   
    self.constraintFile.constraintElements = 2
    
    # 4 resonances from CNS, only 2 in data model... rearranging this here.
    for itemIndex in range(len(self.rawConstraint.items)):
      
      firstMember = self.rawConstraint.items[itemIndex].members[0]
      lastMember = self.rawConstraint.items[itemIndex].members[-1]
      
      # This is a hack, but should work.
      self.rawConstraint.items[itemIndex].members = [firstMember,lastMember]    

  def getConstraintItemMembers(self,item):
    
    # TODO: These are pretty crappy defaults... at least will work for RDC now.
    if self.constraintResNum == 1:
      if self.constraintApiCode == 'Csa':
        itemMembers = [item.members[-1]]
      else:
        itemMembers = [item.members[0]]
    elif self.constraintResNum == 2:
      if self.constraintApiCode == 'Rdc':
        itemMembers = [item.members[-2],item.members[-1]]
      else:
        itemMembers = [item.members[0],item.members[-1]]
      
    return itemMembers
  
  def setRawDistanceConstraint(self):
  
    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraint.serial,self.constraint.serial))

    self.rawConstraint = self.constraintFile.constraints[-1]
    
    #
    # Have to get some sensible values out (even if only list with
    # upper distances read in!)
    #
    
    lowerLimit = self.constraint.lowerLimit    
    if lowerLimit == None:
      if self.constraint.targetValue and self.constraint.targetValue < distanceConstraintDefaultLowerLimit:
        lowerLimit = self.constraint.targetValue
      else:
        lowerLimit = distanceConstraintDefaultLowerLimit
    
    upperLimit = self.constraint.upperLimit
    if upperLimit == None:
      if self.constraint.targetValue:
        upperLimit = self.constraint.targetValue
      else:
        upperLimit = lowerLimit
        
    targetValue = self.constraint.targetValue
    if targetValue == None:
      targetValue = (lowerLimit + upperLimit) / 2
    
    minusDist = targetValue - lowerLimit
    plusDist = upperLimit - targetValue
    
    self.rawConstraint.setDistanceData(targetValue,minusDist,plusDist)
    
    if self.constraint.weight:
      self.rawConstraint.weight = self.constraint.weight
    
    # Also try to set peak information. Assumes only linked to one peak list (this can be done better!)
    peaks = self.constraint.sortedPeaks()
    if peaks and self.usePeakInfo:
      
      # TODO THIS should be done better. Based on peak list or something... .
      if len(peaks) > 1:
        print "  Warning: multiple peaks available for constraint - picking first available."
      
      peak = peaks[0]
      
      self.rawConstraint.peakNum = peak.serial
      
      volume = 0.0
      for peakIntensity in peak.peakIntensities:
        for intKw in (volume_kw,height_kw):
          if peakIntensity.intensityType == intKw:
            volume = peakIntensity.value
            break
                
      self.rawConstraint.peakVol = volume
      
      # TODO AGAIN this should be set beforehand so correct dims selected manually or automatically...
      peakDims = peak.sortedPeakDims()
      ppms = []
      
      for peakDim in peakDims:
        if '1H' in peakDim.dataDimRef.expDimRef.isotopeCodes:
          if len(ppms) < 2:
            ppms.append(peakDim.value)
          else:
            print "  Warning: more than 2 valid dimensions for creating ppm values for constraint list."

      self.rawConstraint.ppms = ppms
          
  def setRawDistanceConstraintItem(self):
    
    self.rawConstraint.items.append(self.rawConstraintItemClass())
    self.rawConstraintItem = self.rawConstraint.items[-1]
    
  def setRawDistanceConstraintItemMembers(self):
    
    for i in range(2):
            
      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])
      
      if atomName == allResidueAtoms_kw:
        atomName = None
      
      self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass(chainCode,seqCode,atomName))

  def setRawRdcConstraint(self):
  
    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraint.serial,self.constraint.serial))

    self.rawConstraint = self.constraintFile.constraints[-1]
    
    #
    # Have to get some sensible values out
    #
    
    lowerLimit = self.constraint.lowerLimit    
    upperLimit = self.constraint.upperLimit        
    targetValue = self.constraint.targetValue
    if targetValue == None:
      targetValue = (lowerLimit + upperLimit) / 2
    
    if self.negateValues:
      targetValue = -targetValue
        
    self.rawConstraint.setRdcData(targetValue,self.constraint.error,self.constraint.weight)
   
  def setRawRdcConstraintItem(self):

    validConstraintItem = True
    if self.filterAtomNames:
      if self.ccpInfo[0][2] not in self.filterAtomNames or self.ccpInfo[1][2] not in self.filterAtomNames:
        validConstraintItem = False
    
    if validConstraintItem:
      self.rawConstraint.items.append(self.rawConstraintItemClass())
      self.rawConstraintItem = self.rawConstraint.items[-1]
    else:
      self.rawConstraintItem = None
    
  def setRawRdcConstraintItemMembers(self):
    
    for coordSystem in ['OO','Z','X','Y']:
      self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass(self.refFrameSegId,self.refFrameResId,coordSystem))
    
    for i in range(0,2):
            
      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])
      
      self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass(chainCode,seqCode,atomName))

  def setRawDihedralConstraint(self):

    #
    # Have to get some sensible values out (even if only list with
    # upper distances read in!)
    #
    
    targetAngle = self.item.targetValue
    lowerLimit = self.item.lowerLimit
    upperLimit = self.item.upperLimit
    
    if not targetAngle:
      
      #
      # Problem here: which direction from upper to lower?
      # Assuming that smallest value comes first (is this then lowerLimit?!)
      #
      
      if upperLimit < lowerLimit and upperLimit < 0:

        upperLimit += 360
          
      if upperLimit < lowerLimit:
      
        print "  Warning: upper limit lower than lower limit... cannot be handled."
        return
        
      else:
        targetAngle = (upperLimit - lowerLimit) / 2 + lowerLimit
      
      devAngle = upperLimit - targetAngle
    
    else:

      if self.item.error:
        devAngle = self.item.error

      else:
        
        #
        # Not entirely correct but will have to do!
        #
          
        devAngle = targetAngle - lowerLimit
        
        
    if self.constraint.weight:
      energyCst = self.constraint.weight
    else:
      energyCst = 2.0
      
    funcExp = 2  # TODO: get this from applicationData?
        
  
    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraint.serial,self.constraint.serial))

    self.rawConstraint = self.constraintFile.constraints[-1]

    self.rawConstraint.setAngleData(targetAngle,devAngle,energyCst,funcExp)

    
  def setRawDihedralConstraintItem(self):
    
    self.rawConstraint.items.append(self.rawConstraintItemClass())
    self.rawConstraintItem = self.rawConstraint.items[-1]
    
  def setRawDihedralConstraintItemMembers(self):
    
    for i in range(self.constraintResNum):
            
      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])
      
      self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass(chainCode,seqCode,atomName))

  def getPresetChainMapping(self,chainList):
  
    return self.getMultiChainFormatPresetChainMapping(chainList)
    

  def setRawCoordinate(self):
  
    if not self.coordinateFile.modelCoordinates.has_key(self.modelId):

      self.coordinateFile.modelCoordinates[self.modelId] = []
    
    chemCompVar = self.residue.chemCompVar
    chemComp = chemCompVar.chemComp
    
    resName = None
    if self.useOriginalData and self.origResLabel != None:
      resName = self.origResLabel

    # Try to get a name from a naming system...
    if not resName:
      resName = self.findChemCompVarSysName('XPLOR',chemCompVar)
    if not resName:
      resName = self.findChemCompVarSysName('CIF',chemCompVar)
    
    # Make sure something gets set
    if not resName:
      if chemComp.code3Letter:
        resName = chemComp.code3Letter
      else:
        resName = chemComp.ccpCode
      
    # Note: chain code will be set as seg ID, unless forceChainId set!
    if self.forceExportSegId != None:
      segId = self.forceExportSegId
    elif self.useOriginalData and self.origChainCode != None:
      segId = self.origChainCode
    else:
      segId = self.exportChainCode
      
    if self.forceExportChainId != None:
      chainId = self.forceExportChainId
    else:
      chainId = None # Will be set to seg ID
  
    if self.useOriginalData and self.origSeqCode != None:
      seqCode = self.origSeqCode
    else:
      seqCode = self.seqCode
 
    if self.useOriginalData and self.origSeqInsertCode != None:
      seqInsertCode = self.origSeqInsertCode
    else:
      seqInsertCode = defaultSeqInsertCode
 
    if self.useOriginalData and self.origAtomName != None:
      atomName = self.origAtomName
    else:
      atomName = self.atomName
 
    #
    # TODO: setting chainId but is really segId... does CNS use the chainID position?
    #

    modelCoordinate = self.coordinatesIO.CnsCoordinate(self.coordinateSerial,
                                                  atomName,
                                                  resName,
                                                  segId,
                                                  seqCode,
                                                  self.x,
                                                  self.y,
                                                  self.z,
                                                  insertionCode = seqInsertCode,
                                                  chainId = chainId)
      
    self.coordinateFile.modelCoordinates[self.modelId].append(modelCoordinate)

    
  def setChemShiftFileValue(self):
    
    if self.atomNamesDict.has_key(self.atomName):
    
      resonanceToAtom = self.atomNamesDict[self.atomName]
      measurement = self.origAtomMeasurements[resonanceToAtom]
       
      self.measurementFileValues.append(
      
               self.rawMeasurementClass(measurement.value,
                                        self.atomName,
                                        self.seqCode,
                                        defaultMolCode))
      
      if self.atomMeasurements.has_key(resonanceToAtom):
        del(self.atomMeasurements[resonanceToAtom])
                           
  def getExportChainCode(self,exportChainCode,chain):
    
    # Returns original chain code - this format uses segId so OK if longer than 1 char
    
    return exportChainCode
