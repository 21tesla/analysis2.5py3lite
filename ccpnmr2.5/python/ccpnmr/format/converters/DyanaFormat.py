#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

DyanaFormat.py: Contains functions specific to Diana/Dyana conversions.

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

import string, copy, traceback, sys

from memops.universal.Io import joinPath, splitPath 

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccp.format.general.Constants import defaultSeqInsertCode

from ccp.general.Util import createMoleculeTorsionDict
from ccp.general.Util import getResonancesFromPairwiseConstraintItem
from ccp.general.Geometry import calcTorsionAngleDegrees

from ccpnmr.format.general.Constants import distanceConstraintDefaultLowerLimit, allResidueAtoms_kw
from ccpnmr.format.general.Util import getNameInfo

#
# Add some information to IOkeywords...
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readDistanceConstraints']['lowerLimitFile'] = (None,False,'The distance constraints lower limit file.')
IOkeywords['readHBondConstraints']['lowerLimitFile'] = (None,False,'The distance constraints lower limit file.')
IOkeywords['writeDistanceConstraints']['lowerLimitFile'] = (None,False,'The distance constraints lower limit file.')
IOkeywords['writeHBondConstraints']['lowerLimitFile'] = (None,False,'The distance constraints lower limit file.')
IOkeywords['writeRdcConstraints']['singleAtom'] = (False,False,'Write out RDC constraints with single atom name (old format).')

class DyanaFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'dyana'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    self.createSequenceFile = self.createSequenceFileGeneric
    
    self.getCoordinates = self.getCoordinatesGeneric
    self.createCoordinateFile = self.createCoordinateFileGeneric 

  #
  # Deviations from generic import stuff
  #

  def getCoordinatesSetSequenceFile(self):
    
    self.setSequenceFileClass()
    self.sequenceFile = self.SequenceFileClass(self.fileName)
    self.sequenceFile.readFromCoordinates(self.coordinateFile)

  def getConstraints(self):
    
    #
    # Format slightly different for constraints....
    #
    
    try:
      
      addKeywords = self.getLowerLimitFile()
    
      self.constraintFile = self.ConstraintFileClass(self.fileName, **addKeywords)
      self.constraintFile.read()

      if self.verbose:
        print "Reading %s constraint list from %s file %s" % (self.constraintType,self.formatLabel,self.fileName)

    except:

      errorMessage = traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
      self.messageReporter.showWarning("Warning"," Cannot read %s constraints for %s...:\n%s" % (self.constraintApiCode,self.formatLabel,errorMessage),self.guiParent)
      self.constraintFile = None
      
      return traceback.format_exception(sys.exc_type,sys.exc_value,sys.exc_info()[2]) 

  def createConstraintFile(self):
  
    try:

      addKeywords = self.getLowerLimitFile(messageType = 'write')

      self.constraintFile = self.ConstraintFileClass(self.fileName, **addKeywords)
    
      if self.verbose:
        print "Writing %s constraints to %s file %s" % (self.constraintType,self.formatLabel,self.fileName)
      
      self.setSpecificConstraintFileWriteInfo()
  
    except:

      self.messageReporter.showWarning("Warning"," No write%sConstraints available for %s..." % (self.constraintApiCode,self.formatLabel),self.guiParent)
      self.quit()

  #
  # Functions different to default functions in DataFormat
  #

  def getLowerLimitFile(self,messageType = 'read'):
    
    addKeywds = {}
    
    if self.constraintType in ['distance','hBond']:
    
      (path,fileName) = splitPath(self.fileName)
      fileName = string.replace(fileName,'upl','lol')
      fileName = joinPath(path,fileName)
      
      #
      # TODO: could here in principle check if file exists and read automatically?
      #
  
      getLowerDist = False

      if not self.lowerLimitFile and not self.minimalPrompts:

        getLowerDist = self.messageReporter.showYesNo("%s lower limit constraint file" % string.capitalize(messageType),
                                                      "Do you want to %s the accompanying lower limit constraint file?" % messageType)

        if getLowerDist:

          component = '%sConstraints' % self.constraintType
          interaction = self.multiDialog.FileName(self.guiParent,file = fileName, component = component, selectionText = 'Select file', format = self.format)
          self.lowerLimitFile = interaction.file

      if hasattr(self,'lowerLimitFile'):
        addKeywds['lowerDistanceFile'] = self.lowerLimitFile
    
    return addKeywds
  
  def setSpecificConstraintFileWriteInfo(self):

    if self.constraintType == 'rdc' and hasattr(self,'singleAtom'):
      self.writeKeywds['singleAtom'] = self.singleAtom
    
  def createSequence(self):
  
    self.sequenceFile.sequences.append(self.sequenceIO.DyanaSequence(molName = self.chain.molecule.name))
    self.sequence = self.sequenceFile.sequences[-1]
    
  def setSequenceFileElements(self):
    
    #
    # Trying to get a sysName for this chemCompVar...
    #
    
    chemCompVar = self.residue.chemCompVar
    namingSystem = chemCompVar.chemComp.findFirstNamingSystem(name = self.namingSystemName)
   
    chemCompSysName = chemCompVar.findFirstSpecificSysName(namingSystem = namingSystem)
   
    if not chemCompSysName:
      chemCompSysName = chemCompVar.findFirstChemCompSysName(namingSystem = namingSystem)
     
    if chemCompSysName:
      residueCode = chemCompSysName.sysName
    else:
      residueCode = self.residue.ccpCode.upper()
    
    self.sequence.elements.append(self.sequenceIO.DyanaSequenceElement(self.seqCode,residueCode))
 
    #
    # Trying to see if it's cis...
    #
    
    molRes = self.residue.molResidue
    linkEnd = molRes.findFirstMolResLinkEnd(linkCode = 'prev')
    
    if linkEnd and linkEnd.molResLink and linkEnd.molResLink.dihedralAngle == 0:
      self.sequence.elements[-1].setCisPeptideBond()

  def setRawDistanceConstraint(self):
  
    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraint.serial))

    self.rawConstraint = self.constraintFile.constraints[-1]
    
    lowerLimit = self.constraint.lowerLimit
    
    if lowerLimit < 0.0:
      lowerLimit = distanceConstraintDefaultLowerLimit
    
    self.rawConstraint.setDistanceData(self.constraint.upperLimit,lowerLimit)
    
  def setRawDistanceConstraintItem(self):
    
    self.rawConstraint.items.append(self.rawConstraintItemClass())
    self.rawConstraintItem = self.rawConstraint.items[-1]
    
  def setRawDistanceConstraintItemMembers(self):
    
    itemResonances = getResonancesFromPairwiseConstraintItem(self.item)
    
    for i in range(0,2):
      
      resonance = itemResonances[i]

      resLabel = self.getResonanceResLabel(resonance)
      
      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])

      if atomName == allResidueAtoms_kw:
        print "  Error: not handling residue level constraints for %s." % self.format
        continue
      
      self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass(chainCode,seqCode,resLabel,atomName))

  def setDihedralConstraintFileInfo(self):
  
    #
    # Get all torsions for reference
    #
  
    self.moleculeTorsions = createMoleculeTorsionDict(self.molSystem)

  def setRawDihedralConstraint(self):
    
    lowerLimit = self.item.lowerLimit
    upperLimit = self.item.upperLimit
    
    keywds = {}
        
    if self.constraint.weight:
      keywds['energyCst'] = self.constraint.weight
    
    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraint.serial,None,self.constraintFile.cyanaLib))

    self.rawConstraint = self.constraintFile.constraints[-1]

    self.rawConstraint.setAngleData(lowerLimit,upperLimit,**keywds)

  def setRawDihedralConstraintItem(self):
    
    self.rawConstraint.items.append(self.rawConstraintItemClass())
    self.rawConstraintItem = self.rawConstraint.items[-1]
    
  def setRawDihedralConstraintItemMembers(self):
    
    #
    # For Dyana, set the angle name (and derive correct angle name from atom names if none found)
    #
    
    (chain,seqId,atomName,resToAtom) = self.ccpInfo[0]
    molecule = chain.molecule
    
    if not self.moleculeTorsions.has_key(molecule):
      if self.verbose:
        print "  Error: unknown molecule '%s'!" % (molecule.name)    
      return
    
    #
    # Find out which residue this torsion angle is for and make a list of chemAtoms...
    #
    
    seqIds= {}
    chemAtoms = []
    
    constraintResonances = self.constraint.resonances
    
    for i in range(0,self.constraintResNum):
      (chain,seqId,atomName,resToAtom) = self.ccpInfo[i]

      if not seqIds.has_key(seqId):
        seqIds[seqId] = 0
      seqIds[seqId] += 1
      
      chemAtoms.append(resToAtom.chemAtom)
    
      resLabel = self.getResonanceResLabel(constraintResonances[i])
      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])
      
      self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass(chainCode,seqCode,resLabel,atomName))

    if len(seqIds) > 1:
      seqId = -99999999999
      maxNum = 0
      for tSeqId in seqIds.keys():
        if seqIds[tSeqId] > maxNum or (seqIds[tSeqId] == maxNum and tSeqId > seqId):
          maxNum = seqIds[tSeqId]
          seqId = tSeqId
          
    else:
      seqId = seqIds.keys()[0]
    
    residue = chain.findFirstResidue(seqId = seqId)
    
    if not residue:
      if self.verbose:
        print "  Error: could not identify residue with seqId %d in chain '%s'!" % (seqId,chain.code)
      return
    
    molRes = residue.molResidue
    
    if not molRes or not self.moleculeTorsions[molecule].has_key(molRes):
      if self.verbose:
        print "  Error: no molresidue for residue with seqId %d in chain '%s'!" % (seqId,chain.code)
      return
    
    #
    # Now get the relevant torsion angles for this molecule, molResidue, and determine which one it is...
    #
    # Note that (in torsion terms) H-N-CA-HA is the same as N-H-HA-CA! (described by _inout)
    #
    
    matchPatterns = {'forward': {}, 'backward': {}, 'forward_inout': {}, 'backward_inout': {}}
    matchInfo = (('forward',(0,1,2,3)),('backward',(3,2,1,0)),('forward_inout',(1,0,3,2)),('backward_inout',(2,3,0,1)))
    torsion = None

    for (torsion,atomList) in self.moleculeTorsions[molecule][molRes]:
       
      #
      # Check matching atoms...
      #
      
      for (matchType,matchIndexes) in matchInfo:
      
        match = ""
        
        # AtomList sometimes not complete!
        if len(atomList) == self.constraintResNum:
          for i in range(0,self.constraintResNum):
            caIndex = matchIndexes[i]
            if chemAtoms[caIndex] == atomList[i][1]:
              match += "1"
            else:
              match += "0"
            
        if match == '1111':
          # All is fine...
          matchPatterns = {}
          break
        else:
          matchPatterns[matchType][match] = torsion
          
      if not matchPatterns:
        break
        
    #
    # Try to see if the angle can be converted...
    #
    
    if matchPatterns or not torsion:
      #print molRes.ccpCode
      #for ca in chemAtoms:
      #  print ca.name,
      #print 
    
      torsion = None

      #for (matchType,matchIndexes) in matchInfo:       
      #  patterns = matchPatterns[matchType].keys()
      #  patterns.sort()
      #  for pattern in patterns:
      #    print "  ",matchType,pattern,
      #    for ca in matchPatterns[matchType][pattern].chemAtoms:
      #      print ca.name,
      #    print

      for pattern in ("0111","1110","0110"):
        for (matchType,matchIndexes) in matchInfo:          
          if matchPatterns[matchType].has_key(pattern):
            torsion = matchPatterns[matchType][pattern]
            
            if matchType[:8] == 'backward':
              chemAtoms.reverse()
            
            if matchType[-5:] == 'inout':
              tca = chemAtoms[:]
              chemAtoms[0] = tca[1]
              chemAtoms[1] = tca[0]
              chemAtoms[2] = tca[3]
              chemAtoms[3] = tca[2]
              
            break
      
      #
      # Have to work 'locally' within chemComp for this... use prev/next coords as well
      #
      
      if torsion:

        # TODO: need something similar to chemComp XML file copying here?!??!
        ccc = None
        for coordSysName in ('ideal','pdb'):
          ccc = self.project.findFirstChemCompCoord(sourceName = coordSysName, molType = molRes.molType, ccpCode = molRes.ccpCode) # , showError = self.messageReporter.showError)
          if ccc:
            break
        
        # TODO: really need to calculate angle between two sets of 3 atoms, with 2 atoms shared (??)
        angles = {}
        #print torsion.name

        for (angleType,caList) in (('ref',torsion.chemAtoms),('actual',chemAtoms)):
          
          coords = []
          angles[angleType] = None
          
          for ca in caList:
            if not ca:
              coords.append(None)
            else:
              #print ca.name,
              coord = ccc.findFirstChemAtomCoord(chemAtom = ca)
              if coord:
                coords.append((coord.x,coord.y,coord.z))
              else:
                coords.append(None)
                
          #print
          #print coords
          
          if None not in coords:
            angles[angleType] = calcTorsionAngleDegrees(coords[0], coords[1], coords[2], coords[3])
          else:
            #
            # Special case... O used instead of CA for OMEGA
            #
            if torsion.name == 'OMEGA' and angleType == 'actual':
              if coords[2] and coords[3]:
                angles['actual'] = angles['ref'] - 180.00
                
                if chemAtoms[3].name == 'H':
                  angles['actual'] += 180.00
        
        if None not in angles.values():
          
          #
          # Have all reference info, now reset the value...
          #
          
          #print angles
          
          # TODO: IS THIS ALWAYS CORRECT??? CHECK WITH WAYNE!!
          angleDiff = angles['ref']-angles['actual']
                    
          lowerLimit = self.rawConstraint.lowerAngle + angleDiff
          upperLimit = self.rawConstraint.upperAngle + angleDiff
          
          if lowerLimit > 360.0 and upperLimit > 360.0:
            lowerLimit -= 360.0
            upperLimit -= 360.0
            
          #print angleDiff, lowerLimit, upperLimit
          
          self.rawConstraint.setAngleData(lowerLimit,upperLimit)
          
        else:
          torsion = None
    
    if torsion:
      namingSystem = torsion.chemComp.findFirstNamingSystem(name = self.namingSystemName)
      
      if namingSystem:
        torsionSysName = torsion.findFirstSysName(namingSystem = namingSystem)
        if torsionSysName:
          self.rawConstraint.name = torsionSysName.sysName
        elif self.verbose:
          print "  Error: could not find a %s sysName for %s, angle name %s" % (self.namingSystemName,molRes.ccpCode,torsion.name)
      elif self.verbose:
        print "  Error: could not find a %s naming system for chemComp %s" % (self.namingSystemName,molRes.ccpCode)
    else:
      if self.verbose:
        if None not in chemAtoms:
          print "  Error: could not define a torsion angle name for %s, atoms %s-%s-%s-%s." % (molRes.ccpCode,chemAtoms[0].name,chemAtoms[1].name,chemAtoms[2].name,chemAtoms[3].name)
        else:
          print "  Error: could not define a torsion angle name for %s because of missing chemAtoms" % (molRes.ccpCode)
    
  def setRawRdcConstraint(self):
  
    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraint.serial,self.constraintFile.cyanaLib))

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
        
    self.rawConstraint.setRdcData(targetValue,error = self.constraint.error,energyCst = self.constraint.weight)
   
  def setRawRdcConstraintItem(self):
    
    self.rawConstraint.items.append(self.rawConstraintItemClass())
    self.rawConstraintItem = self.rawConstraint.items[-1]
    
  def setRawRdcConstraintItemMembers(self):
        
    itemResonances = getResonancesFromPairwiseConstraintItem(self.item)
    
    for i in range(0,2):
            
      (chainCode,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])
      resLabel = self.getResonanceResLabel(itemResonances[i])
      
      self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass(chainCode,seqCode,resLabel,atomName))

  def setRawCoordinate(self):
  
    if not self.coordinateFile.modelCoordinates.has_key(self.modelId):

      self.coordinateFile.modelCoordinates[self.modelId] = []
    
    chemCompVar = self.residue.chemCompVar
    chemComp = chemCompVar.chemComp
    cifCode = self.findChemCompVarSysName('CIF',chemCompVar)
    
    if self.useOriginalData and self.origResLabel != None:
      resName = self.origResLabel
    elif cifCode:
      resName = cifCode
    elif chemComp.code3Letter:
      resName = chemComp.code3Letter
    else:
      resName = chemComp.ccpCode.upper()
      
    
    if self.useOriginalData and self.origChainCode != None:
      chainCode = self.origChainCode
    else:
      chainCode = self.exportChainCode
  
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
    # TODO: set chainId? Is not segId... does CNS use the chainID position?
    #

    modelCoordinate = self.coordinatesIO.DyanaCoordinate(self.coordinateSerial,
                                                    atomName,
                                                    seqCode,
                                                    resName,
                                                    '',
                                                    self.x,
                                                    self.y,
                                                    self.z,
                                                    insertionCode = seqInsertCode)
      
    self.coordinateFile.modelCoordinates[self.modelId].append(modelCoordinate)

  def getResonanceResLabel(self,resonance):
  
   resLabel = None
  
   chemCompVar = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue.chemCompVar
   namingSystem = chemCompVar.chemComp.findFirstNamingSystem(name = self.namingSystemName)
   
   chemCompSysName = chemCompVar.findFirstSpecificSysName(namingSystem = namingSystem)
   
   if not chemCompSysName:
     chemCompSysName = chemCompVar.findFirstChemCompSysName(namingSystem = namingSystem)
     
   if chemCompSysName:
     resLabel = chemCompSysName.sysName
   else:
     resLabel = chemCompVar.chemComp.ccpCode.upper()
   
   return resLabel

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)
