"""
======================COPYRIGHT/LICENSE START==========================

XEasyFormat.py: Contains functions specific to XEasy conversions.

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

import traceback, copy

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccp.general.Util import setUniqueAppData, findAllSysNamesByChemAtomOrSet
from memops.universal.Util import returnInt

from ccpnmr.format.general.Util import getApplResNames
from ccpnmr.format.general.Constants import volume_kw, height_kw, peakNum_kw

from ccp.lib.nmr.Nmr.DataSource import getXEasyDimCodes

import memops.api.Implementation as Implementation
import ccp.api.molecule.ChemComp as ChemComp
import ccp.api.nmr.Nmr as Nmr

#
# Add some information to IOkeywords...
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readPeaks']['cyanaFormat'] = (None,False,'Read as CYANA format (with multiple assignments).')
IOkeywords['readPeaks']['forceHeight'] = (None,False,'Force import of XEASY volumes as heights in CCPN.')
IOkeywords['writePeaks']['cyanaFormat'] = (None,False,'Write as CYANA format (with multiple assignments).')
IOkeywords['writePeaks']['integrationMethod'] = ('volume',False,'Type of integration method to use in output file.')
IOkeywords['readShifts']['readAllSerials'] = (None,False,'Read all atom serials, even if no shifts available.')
IOkeywords['writeShifts']['cyana21Naming'] = (False,False,'Use CYANA 2.1 atom names.')
IOkeywords['readShiftsAndPeaks']['peakAssignmentFile'] = (None,False,'Peak assignment file.')
IOkeywords['readShiftsAndPeaks']['cyanaFormat'] = (None,False,'Read as CYANA format (with multiple assignments).')

class XEasyFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'xeasy'
    self.IOkeywords = IOkeywords
    
    self.shiftListKeyword = 'shiftList'
    self.atomSerialKeyword = 'atomSerial'

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    self.createSequenceFile = self.createSequenceFileGeneric

    #self.getPeaks = self.getPeaksGeneric
    self.createPeakFile = self.createPeakFileGeneric
    
    self.getMeasurements = self.getMeasurementsGeneric
    self.createMeasurementFile = self.createMeasurementFileGeneric

  #
  # Deviations from generic import stuff
  #

  def createMeasurementFileFormatSpecific(self):
  
    #
    # This is specific for XEasy only! Is only necessary if resonances were 'compressed'
    #
    
    self.atomSerialDict = {}

  # TODO: spinsystems!!
  
  def getMeasurementsSetFormatSpecificReadKeywds(self):
    
    readKeywds = {}
    
    if self.measurementApiCode == 'Shift' and self.readAllSerials:
      readKeywds = {'readAll': True}
    
    return readKeywds
  
  def getPeaks(self):
  
    if self.verbose == 1:
      print "Reading peak list from %s file %s" % (self.formatLabel,self.fileName)

    # TODO HERE: have to figure out what to do if project file read...
    # TODO: Xeasy peak colour, userCode do not fit right now.
  
    # Assuming only one peaklist per file  

    self.peakFile = self.PeakFileClass(self.fileName)
    self.peakFile.read(cyanaFormat = self.cyanaFormat)

  def createPeakFileFormatSpecific(self):

    self.peakFile.setSpectrumInfo(self.fileName,self.numPeakDim)
    self.writeKeywds['cyanaFormat'] = self.cyanaFormat
  
    
  
  def getPeakAssignments(self):
  
    if self.verbose == 1:
      print "Reading peak assignments from %s file %s" % (self.formatLabel,self.fileName)

    # TODO HERE: have to figure out what to do if project file read...
    # TODO: Xeasy peak colour, userCode do not fit right now.
  
    # Assuming only one peaklist per file  

    self.peakAssignmentsFile = self.PeakAssignmentFileClass(self.fileName)
    self.peakAssignmentsFile.read()
    
  def createPeakAssignmentFile(self):
  
    if self.verbose == 1:
      print "Writing peak assignments to %s file %s" % (self.formatLabel,self.fileName)
  
    self.peakAssignmentsFile = self.PeakAssignmentFileClass(self.fileName)
    self.peakAssignmentsFile.setDimensions(self.numDim)

    # TODO TODO get this sorted correctly!! Will probably fall over
    # if information from other formats used.
    # Select relevant dimensions for XEasy .assign file
    # ONLY last two selected (this is a hack - only two
    # should be passed in/used for XEasy)
    #
  
    for dataDimRefIndex in self.dimMapping[-2:]:
      self.peakAssignmentsFile.assignmentsDim.append(dataDimRefIndex)

  #
  # Functions different to default functions in DataFormat
  #
    
  def createSequence(self):
  
    self.sequenceFile.sequences.append(self.sequenceIO.XEasySequence(molName = self.chain.molecule.name))
    self.sequence = self.sequenceFile.sequences[-1]
    
  def setSequenceFileElements(self):
      
      self.sequence.elements.append(self.sequenceIO.XEasySequenceElement(self.seqCode,self.residue.molResidue.chemComp.ccpCode))

  def forceSpecificNamingSystemName(self):
  
    if (hasattr(self,'cyana21Naming') and self.cyana21Naming) or (hasattr(self,'cyanaFormat') and self.cyanaFormat):
      self.forceNamingSystemName = self.getFormatNamingSystemName(format = 'cyana', version = '2.1') 
      print "  Warning: setting naming system to %s." % self.forceNamingSystemName

  def createShift(self,resonance,chemShift):
    
    #
    # For xeasy also use chemShift to create atomSerial application data.
    #
    # Warning: if there is already an atomSerial available for this resonance
    # it will be overwritten!!
    #
    
    applData = resonance.findFirstApplicationData(application = self.format, keyword = self.atomSerialKeyword)

    if applData:
      resonance.removeApplicationData(applData)
    
    resonance.addApplicationData(Implementation.AppDataInt(application = self.format, keyword = self.atomSerialKeyword, value = chemShift.atomSerial))
    
    #
    # Create shift linked to resonance - NOTE changed so can read in all atom serials, even if no shift value given!
    #
    
    if chemShift.value:

      self.convertCount[self.mainCode][1][self.valueCode] += 1
      shift = self.measurementList.newShift(value = chemShift.value, error = chemShift.valueError, resonance = resonance)

    else:
     
      shift = None
     
    return shift

  def setPeakIntensity(self):
  
    self.intMethod = self.peakFile.translateIntMethod(self.format,self.rawPeak.intCode)
    
    if self.intMethod != 'Not integrated' or self.rawPeak.volume != None:
        
      if self.forceHeight or self.intMethod == 'Not integrated':
        
        # PeakIntensity attributes
        peakInt = self.peak.newPeakIntensity(value = self.rawPeak.volume, method = self.methods['unknown']['Intensity'])
        intensityType = height_kw
        
      else:
   
        # PeakIntensity attributes
        # Assuming that always volume for xeasy (all integration)
        peakInt = self.peak.newPeakIntensity(value = self.rawPeak.volume, method = self.methods[self.format][self.intMethod])
        intensityType = volume_kw
      
      peakInt.error = self.rawPeak.volumeError
      peakInt.intensityType = intensityType

  def setPeakDim(self):
  
    dataDimRef = self.dataDimRefs[self.rawPeakDimIndex]

    self.peakDim = self.peak.findFirstPeakDim(dim = dataDimRef.dataDim.dim)

    self.peakDim.dataDimRef = dataDimRef
    
    self.peakDim.value = self.rawPeak.ppm[self.rawPeakDimIndex]

  def createResonanceNamesDict(self,resonances):
    
    #
    # For XEasy am using atomSerials instead of resNames for matching...
    #
    
    self.assignKeyword = self.atomSerialKeyword
    
    if self.resonanceMatching:
      self.resonanceNames = getApplResNames(self.format, resonances, keyword = self.assignKeyword)
    
    else:
      self.resonanceNames = {}
      
    if self.resonanceNames == {}:
      self.messageReporter.showWarning("Warning","No atom serial numbers available for xeasy peak assignment linking - assignments cannot be transferred",self.guiParent)
      
  def setPeakExtras(self):
  
    if self.cyanaFormat and self.rawPeak.ambiguousAssign:
      
      #
      # Make peakContribs for this one...
      #
      
      numContribs = len(self.rawPeak.ambiguousAssign) + 1
      
      for i in range(numContribs):
        self.peakContribs.append(Nmr.PeakContrib(self.peak))

  def getPeakResNames(self):
    
    #
    # Assignments put in application data as numbers (need .prot file for atom linking)
    # Warning: Keyword here is 'atomSerial'
    #
    
    if self.rawPeak.assign[self.rawPeakDimIndex] != None:
    
      resName = self.rawPeak.assign[self.rawPeakDimIndex]
      
      #
      # Ignore peaks with atomSerial 0
      #
      
      if resName == 0:
        self.resNames = []
        
      else:
        self.resNames = [resName]
        
        if self.cyanaFormat and self.rawPeak.ambiguousAssign:
        
          for ambAssignment in self.rawPeak.ambiguousAssign:
            
            resName = ambAssignment[self.rawPeakDimIndex]
            
            if resName == 0:
              self.resNames.append(None)
            else:
              self.resNames.append(resName)
    
    else:
    
      self.resNames = []
    
  def peakListAssignmentCheck(self,peakLists):
  
    #
    # Check whether shift file was written... if not then do it!
    #
    
    checkResonances = []
    
    for peakList in peakLists:
  
      for peak in peakList.sortedPeaks():
      
        for peakDim in peak.sortedPeakDims():
        
          for peakDimContrib in peakDim.peakDimContribs:
          
            if isinstance(peakDimContrib,Nmr.PeakDimContrib):

              resonance = peakDimContrib.resonance
        
              if resonance not in checkResonances:

                checkResonances.append(resonance)
    
    atomSerialsFound = 0
         
    for resonance in checkResonances:
      if self.getResonanceAtomSerial(resonance,verbose = 0):
        atomSerialsFound += 1
    
    if checkResonances:
      percentageFound = float(atomSerialsFound) * 100 / len(checkResonances)
    else:
      percentageFound = 100
    
    continueWriting = True
    
    if not self.minimalPrompts and percentageFound < 90:
    
      continueWriting = self.messageReporter.showYesNo("Write shift list","Only %.1f percent of the resonances has an atom serial - assignments cannot be transferred before a shift list is written. Do you want to continue?" % percentageFound)
    
    
    return continueWriting
    
  def setPeakFilePeakExtras(self):
  
    # TODO: Xeasy peak colour, userCode do not fit right now.
    
    self.colour = 1              # TODO: application data?
    self.userCode = "T"          # TODO: application data?
    self.assign = []             # Can only do this if atomSerial assigned...
    self.ambiguousAssign = []    # Can only do this if atomSerial assigned and cyanaFormat option selected.
    self.ppm = []
    
  def setPeakFilePeakIntensity(self):
  
    self.volume = 0.0
    self.volumeError = 0.0
    self.intCode = '-'

    intensity = self.peak.findFirstPeakIntensity(intensityType = self.integrationMethod)
    
    if intensity:
      
      self.volume = intensity.value
      
      if intensity.error != None:
        self.volumeError = intensity.error
      else:
        self.volumeError = 0.0
        
      intMethod = intensity.method.procedure
      
      #
      # XEasy specific codes...
      # Now assuming 'automatic' as default if a volume is available
      #
      
      self.intCode = self.peakFile.translateIntMethod('ccpn',intMethod)
  
  def getAtomSerial(self,resToAtom):

    atomSerial = 0

    if resToAtom:
      if resToAtom.useAtomSetName:
        atomName = resToAtom.atomSetName
      else:
        atomName = resToAtom.atomName    

      atomSerialDictKey = (resToAtom.chain,resToAtom.seqId,atomName)

      if self.atomSerialDict.has_key(atomSerialDictKey):
        atomSerial = self.atomSerialDict[atomSerialDictKey]
      else:
        self.addWarning("Key %s not found in self.atomSerialDict." % str(atomSerialDictKey))
  
    return atomSerial
    
  def setPeakFilePeakDimInfo(self):
    
    self.ppm.append(self.peakDim.getValue())
    
    atomSerial = 0
    
    #
    # Transfer assignment for compressresonances - only works if cs list written earlier!
    # 
    
    assignLen = len(self.peakAssignmentList)
    
    if hasattr(self,'atomSerialDict') and assignLen >= 1:
    
      peakDimIndex = list(self.peak.sortedPeakDims()).index(self.peakDim)

      if self.cyanaFormat or assignLen == 1:
        resToAtom = self.peakAssignmentList[0][peakDimIndex]
        atomSerial = self.getAtomSerial(resToAtom)

      #
      # Do ambiguous ones if possible.
      #

      if not self.cyanaFormat:

        self.addWarning("Multiple assignments for XEasy peak. Have to write .assign file.")

      else:

        for i in range(1,assignLen):

          resToAtom = self.peakAssignmentList[i][peakDimIndex]
          otherAtomSerial = self.getAtomSerial(resToAtom)

          if len(self.ambiguousAssign) < assignLen - 1:
            self.ambiguousAssign.append([])

          self.ambiguousAssign[i-1].append(otherAtomSerial)
        
    self.assign.append(atomSerial)      

  def peakAssignmentsSetAssignmentInfo(self):
  
    # TODO: redundant - could use self.peakAssignmentsList!!!
    # Get all assignments
    #
    
    self.assignCombinations = []
    peakDims = list(self.peak.sortedPeakDims())

    for peakContrib in self.peak.sortedPeakContribs():
    
      peakDimContribs = self.numDim * [None]
      
      #
      # Get all peakDimContribs connected to this peakContrib
      # and organize them according to peakDim
      #
    
      for peakDimContrib in peakContrib.peakDimContribs:
      
        if isinstance(peakDimContrib,Nmr.PeakDimContribN):

          print "   Warning: cannot handle multiple resonances for one peakDim contribution. Ignored."
      
        else:
          
          peakDimIndex = peakDims.index(peakDimContrib.peakDim)
          peakAssignmentIndex = self.dataDimRefs.index(self.peak.sortedPeakDims()[peakDimIndex].dataDimRef)
          
          atomSerial = self.getResonanceAtomSerial(peakDimContrib.resonance)
          
          if not peakDimContribs[peakAssignmentIndex]:
            peakDimContribs[peakAssignmentIndex] = []
            
          peakDimContribs[peakAssignmentIndex].append(atomSerial)
    
      #
      # Now loop over the peakDimContribs and make all possible combinations
      #
        
      self.setAssignCombinations(peakDimContribs,0,[])

 
  def peakAssignmentsCreate(self):
    
    if self.assignCombinations:
    
      appData = self.peak.findFirstApplicationData(application = self.format, keyword = peakNum_kw)
      
      if appData:
        peakNum = appData.value
      else:
        peakNum = self.peak.serial

      peakAssignment = self.peakAssignmentsIO.XEasyPeakAssignment(self.peakAssignmentsFile,peakNum)
      
      # 
      # set the assignCombinations
      #      

      for assignCombination in self.assignCombinations:
        
        combi = [assignCombination[self.peakAssignmentsFile.assignmentsDim[0]],assignCombination[self.peakAssignmentsFile.assignmentsDim[1]]]
        peakAssignment.assignCombinations.append(combi)

      self.peakAssignmentsFile.peakAssignments.append(peakAssignment)
  
  def setAssignCombinations(self,peakDimContribs,peakAssignmentIndex,currentContrib):
  
    """
    
    Recursive def. Loops over all possible combinations for peakDimContribs
    and sets them in self.assignCombinations in the order of the FORMAT
    peak assignment file!!!
    
    """
  
    if len(peakDimContribs) == peakAssignmentIndex:
    
      self.assignCombinations.append(currentContrib)

    else:
    
      if peakDimContribs[peakAssignmentIndex]:
    
        for atomSerial in peakDimContribs[peakAssignmentIndex]:
      
          self.setAssignCombinations(peakDimContribs,peakAssignmentIndex + 1, currentContrib + [atomSerial])
      
      else:
      
          self.setAssignCombinations(peakDimContribs,peakAssignmentIndex + 1, currentContrib + [0])

  def getResonanceAtomSerial(self,resonance,verbose = 1):
    
    atomSerial = 0
    
    #
    # First look if there's an atomSerial from a cs list that was written
    #
    
    applData = resonance.findFirstApplicationData(application = self.format, keyword = self.shiftListKeyword)
    
    if not applData:
     
      #
      # Then look if there's original (read in from XEasy file) numbers
      #
    
      applData = resonance.findFirstApplicationData(application = self.format, keyword = self.atomSerialKeyword)
    
    if not applData:
      if verbose:
        self.addWarning("No atomSerial for XEasy peak %d. Try writing XEasy shift file first." % self.peak.serial)
    
    else:
      atomSerial = applData.value
      
    return atomSerial

  def createPeakFilePeak(self):
    
    self.peakFile.peaks.append(self.peaksIO.XEasyPeak(self.peakNum,self.ppm,self.colour,self.userCode,self.volume,self.volumeError,self.intCode,self.assign))
    
    if self.cyanaFormat:
      self.peakFile.peaks[-1].ambiguousAssign = self.ambiguousAssign[:]

  def setPeakFileInfo(self):    
    
    
    dimCodes = getXEasyDimCodes(self.peakList.dataSource)
    
    self.peakFile.dimCodes = dimCodes

    #
    # Also do a check for which intensity type is available, reset if not good
    #
    
    intensityTypes = []
    
    for peak in self.peakList.peaks: 
      for intensity in peak.findAllPeakIntensities():
        if not intensity.intensityType in intensityTypes:
          intensityTypes.append(intensity.intensityType)
    
    if self.integrationMethod not in intensityTypes and len(intensityTypes) == 1:
      print "  Warning: resetting peak integration method to %s - no %s values available." % (intensityTypes[0],self.integrationMethod)
      self.integrationMethod = intensityTypes[0]
      
  def orderByAtomNames(self,residueResonanceToAtomList):

    #
    # Only use compressResonances if useOriginalResNames NOT on!!
    #
    

    
    if self.presetResonanceMapping:
      self.compressPresetResidueResonanceToAtomList(residueResonanceToAtomList)  

    elif self.compressResonances and not self.useOriginalResNames:
      
      self.compressResidueResonanceToAtomList(residueResonanceToAtomList)    
   
    #
    # Order according to applData 'atomOrder' linked to XEasy chemAtomSysNames
    #
    
    chemAtomOrder = {}

    #
    # Use 'middle', 'neutral' chemComp for reference info
    #
   
    chemCompVar = self.residue.chemCompVar
    
    if not chemCompVar:
      chemCompVar = self.residue.molResidue.chemComp.findFirstChemCompVar(linking = 'middle', isDefaultVar = True)
      
      if not chemCompVar:
        print "  ERROR: no chemCompVar found for XEasy - might be non-standard residue."
        return
        
    chemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemCompVar.chemComp,chemCompVar.chemAtoms,self.namingSystemName)
    chemAtomSetSysNames = findAllSysNamesByChemAtomOrSet(chemCompVar.chemComp,chemCompVar.chemAtomSets,self.namingSystemName)
    namingSystem = chemCompVar.chemComp.findFirstNamingSystem(name = self.namingSystemName)

    for chemAtomSysName in chemAtomSysNames + chemAtomSetSysNames:
    
      chemAtomOrSet = chemCompVar.findFirstChemAtom(name = chemAtomSysName.atomName, subType = chemAtomSysName.atomSubType)
        
      if not chemAtomOrSet:
      
        chemAtomOrSet = chemCompVar.findFirstChemAtomSet(name = chemAtomSysName.atomName, subType = chemAtomSysName.atomSubType)
        
      if not chemAtomOrSet:
        print "  Warning: ignoring atom with name %s for %s export." % (chemAtomSysName.atomName,self.format)
        continue
      
      applData = None
      if namingSystem:     
        serialKeyword = "atomSerial_%s.%d" % (chemAtomOrSet.name,chemAtomOrSet.subType)
        applData = namingSystem.findFirstApplicationData(application= self.format, keyword=serialKeyword)
      
      # TODO THIS IS A SLIGHT HACK TO AVOID UNNECESSARY PRINTOUTS
      if not applData:
        if chemAtomSysName.sysName not in ['HN','O','H',"O'","O''"]:
          print "  Warning: no XEASY serial number for %s ('%s','%s'), sysName %s" % (self.residue.ccpCode,chemCompVar.linking,chemCompVar.descriptor,chemAtomSysName.sysName)
        continue
        
      serialOrder = applData.value
      
      chemAtomOrder[serialOrder] = chemAtomSysName.sysName
      
      #
      # Already set the atomNamesDict...
      #
      
      if self.useOriginalResNames:
      
        for resonanceToAtom in residueResonanceToAtomList:
        
          if chemAtomSysName.sysName == resonanceToAtom.atomName:
      
            self.atomNamesDict[chemAtomSysName.sysName] = resonanceToAtom
            break
      
      else:
      
        atomFound = 0
      
        for resonanceToAtom in residueResonanceToAtomList:
        
          if chemAtomSysName.sysName == resonanceToAtom.atomName:
          
            self.atomNamesDict[chemAtomSysName.sysName] = resonanceToAtom
            atomFound = 1
            
            break        
          
        if not atomFound:

          #
          # Could be resonanceToAtom present for IUPAC name (e.g. CD* for CD1, CD2)
          #
          
          if isinstance(chemAtomOrSet,ChemComp.ChemAtom):
            chemAtomSet = chemAtomOrSet.chemAtomSet

            if chemAtomSet:

              for resonanceToAtom in residueResonanceToAtomList:

                if chemAtomSet.name == resonanceToAtom.atomName:

                  self.atomNamesDict[chemAtomSysName.sysName] = resonanceToAtom
                  break        
 
    #
    # Make sure atomNames are in XEasy order...
    #
    
    atomsIndexList = chemAtomOrder.keys()
    atomsIndexList.sort()
        
    for atomIndex in atomsIndexList:
    
      self.atomNamesList.append(chemAtomOrder[atomIndex])
    
  def setChemShiftFileValue(self):

    #
    # Warning: chainCode is not actually used here... 
    #
    
    value = 999.000
    error = 0.000
    
    if self.atomNamesDict.has_key(self.atomName):
      
      resonanceToAtom = self.atomNamesDict[self.atomName]
      
      chemShift = self.origAtomMeasurements[resonanceToAtom]
      
      value = chemShift.value
      error = chemShift.error
      
      if self.atomMeasurements.has_key(resonanceToAtom):
      
        del(self.atomMeasurements[resonanceToAtom])

      #
      # Set serial number as application data for resonance
      # TODO: this is now only done when writing cs file: otherwise
      # serials don't have meaning anyway!!!
      #
      
      setUniqueAppData('AppDataInt',chemShift.resonance,self.format,self.shiftListKeyword,self.atomSerial)
       
    self.measurementFile.chemShifts.append(
    
              self.rawMeasurementClass(self.measurementFile,
                                       self.atomSerial,
                                       value,
                                       error,
                                       self.atomName,
                                       self.seqCode))
    
    
    self.atomSerialDict[(self.chain,self.residue.seqId,self.atomName)] = self.atomSerial
    
    # TODO: also have to include spin system related info:
    #       see 'xeasySpinSysAtoms' for XEasy!!!

  def setPeakDimOrder(self):
  
    self.peakDimOrder = {
       
       1: [0],
       2: [1,0],
       3: [1,2,0],
       4: [1,2,3,0],
       5: [1,2,3,4,0],
       6: [1,2,3,4,5,0]
       
      }

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)
