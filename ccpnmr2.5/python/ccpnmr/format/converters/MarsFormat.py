"""
======================COPYRIGHT/LICENSE START==========================

MarsFormat.py: Contains functions specific to MARS conversions.

Copyright (C) 2007 Wim Vranken (European Bioinformatics Institute)

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

import os, sys, copy, string

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords
from ccpnmr.format.general.Constants import newline
from ccp.general.ChemCompOverview import chemCompStdDict
from memops.universal.Io import normalisePath

#from ccpnmr.format.general.Util import getResName

import ccp.api.nmr.Nmr as Nmr

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readAssignment']['reliabilityCutoff'] = ('H',True,'The cutoff for including the assignment (H/M/L, with L being the lowest reliability).')
IOkeywords['writeMeasurements']['ignoreFixedConnection'] = (False,False,'Write out a fixed connections file for the spin systems.')
IOkeywords['writeProject']['chain'] = (None,True,'The chain to be analyzed in MARS.')
IOkeywords['writeProject']['shiftList'] = (None,True,'The shift list to use for the MARS input.')
IOkeywords['writeProject']['inputDirectory'] = (None,True,'The directory to write the MARS input files to.')
IOkeywords['writeProject']['useFixedConnection'] = (True,False,'Write out a fixed connections file for the spin systems.')
IOkeywords['writeProject']['useFixedAssignments'] = (True,False,'Write out a fixed assignments file for the spin systems.')

class MarsFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'mars'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    self.createSequenceFile = self.createSequenceFileGeneric
    
    self.getMeasurements = self.getMeasurementsGeneric
    self.createMeasurementFile = self.createMeasurementFileGeneric
    
    self.getAssignment = self.getAssignmentGeneric

  #
  # Deviations from generic import stuff
  #

  """
  def getFullProject(self,fileName, verbose = 0):
  
    #
    # TODO TODO REALLY should only read output _final file - although there should be an option
    # to read the spin system information from the input file only.
    #
  
    if self.verbose == 1:
      print "Reading %s project from file %s" % (self.formatLabel,fileName)
   
    self.file = self.projectIO.MonteProjectFile(fileName)
    self.file.read(verbose = verbose)
    
    (projectDir,projectFileName) = os.path.split(fileName)
    
    csFiles = {}
    
    # TODO: need to be able to pass in these file paths/names!!
    for (fileType,fileFormat,isOutputFile) in (('input',"%s.cs",0),('output',"%s_final",1)):
      for subDir in (self.file.subDirPath,''):
        csFileName = os.path.join(projectDir,subDir,fileFormat % self.file.rootFileName)
        if os.path.exists(csFileName):
          csFile = self.chemShiftsIO.MonteChemShiftFile(csFileName)
          csFile.read(headerCols = self.file.headerCols,isOutputFile = isOutputFile, verbose = verbose)
          csFiles[fileType] = csFile
          break
    
    #
    # Rearrange if necessary
    #
    
    if len(csFiles.keys()) == 1:
      self.measurementFile = csFiles.values()[0]

    elif len(csFiles.keys()) == 2:
      
      assignedStripIds = []
      for chemShift in csFiles['output'].chemShifts:
        if chemShift.stripId not in assignedStripIds:
          assignedStripIds.append(chemShift.stripId)

      for chemShift in csFiles['input'].chemShifts:
        if chemShift.stripId not in assignedStripIds:
          csFiles['output'].chemShifts.append(chemShift)
          
      self.measurementFile = csFiles['output']
      
    self.readShifts(minimalPrompts = self.minimalPrompts)
    
  """
  
  def createFullProject(self,fileName):

    if self.verbose == 1:
      print "Writing %s project from file %s" % (self.formatLabel,fileName)
      
    self.file = self.projectIO.MarsProjectFile(fileName)
    
    (fileDir,fileName) = os.path.split(fileName)
    fileDir = normalisePath(fileDir)
    marsInputDir = normalisePath(self.inputDirectory)
    if marsInputDir != fileDir and marsInputDir.count(fileDir):
      marsInputDir = marsInputDir.replace(fileDir,'')
      if marsInputDir[0] == os.sep:
        marsInputDir = marsInputDir[1:]
        
    if not os.path.exists(marsInputDir):
      os.mkdir(marsInputDir)
    
    #
    # Set tag for export/import cycle
    #
    
    #self.setDateTimeStamp(self.project,setTag,isUnique = False)
    #self.setDateTimeStamp(self.chain,setTag,isUnique = True) # Only one export per chain!
    
    #    
    # FASTA type sequence file
    #
    
    sequenceFileName = 'fastaSequence.tab'
    self.file.files['sequence'] = os.path.join(marsInputDir,sequenceFileName)
    self.writeSequence(os.path.join(self.inputDirectory,sequenceFileName), chains = [self.chain], minimalPrompts = self.minimalPrompts)

    #
    # Chemical shift tab file - also used to set up spin system info
    #
    
    chemShiftFile = 'cs.tab'
    self.file.files['csTab'] = os.path.join(marsInputDir,chemShiftFile)
    
    if self.useFixedConnection:
      ignoreFixedConnection = True
    else:
      ignoreFixedConnection = False
    
    # Need to keep track of info in measurement file    
    self.writeShifts(os.path.join(self.inputDirectory,chemShiftFile), measurementList = self.shiftList, chains = [self.chain], ignoreFixedConnection = ignoreFixedConnection, noWrite = True, minimalPrompts = self.minimalPrompts)
    self.measurementFile.write(**self.writeKeywds)
    
    #
    # Fixed connection info - TODO make a ccp/format level file for this? Or not worth it?
    #

    if self.useFixedConnection:
      fixedConnectionFile = 'fix_con.tab'
      self.file.files['fixConn'] = os.path.join(marsInputDir,fixedConnectionFile)
      
      fixConnOut = open(os.path.join(self.inputDirectory,fixedConnectionFile),'w')
      
      # Ignore info written in chemical shift file - is already filtered to not include chain-linked spin systems
      ignoreSsSerials = []
      for ssSerial in self.measurementFile.spinSystemLinks.keys():
        ignoreSsSerials.append(ssSerial)
        ignoreSsSerials.append(self.measurementFile.spinSystemLinks[ssSerial][-1])

      # Now go over spins systems, and print out links
      serialConnections = []
      for spinSystem in self.spinSystemList:
        ssSerial = spinSystem.serial        
        if ssSerial in ignoreSsSerials:
          continue

        # Only report spin systems that are part of the chemical shift tab file!!
        if not self.spinSystemToPrevId.has_key(spinSystem) or not self.spinSystemToPrevId[spinSystem]:
          continue
        
        refResonanceToSpinSystem = self.spinSystemToResonanceDict[spinSystem][0]
        
        serialKey = None
        if refResonanceToSpinSystem.linkedSpinSystems.has_key(1):
          serialKey = (spinSystem.serial,refResonanceToSpinSystem.linkedSpinSystems[1].serial)
        elif refResonanceToSpinSystem.linkedSpinSystems.has_key(-1):
          serialKey = (refResonanceToSpinSystem.linkedSpinSystems[-1].serial,spinSystem.serial)
        
        if serialKey and serialKey not in serialConnections:
          serialConnections.append(serialKey)
      
      for serialKey in serialConnections:
        fixConnOut.write("%-4d  %-4d" % (serialKey[0],serialKey[1]) + newline)
      fixConnOut.close()
    
    #
    # Fixed assignment info - TODO make a ccp/format level file for this? Or not worth it?
    #
    
    if self.useFixedAssignments:
      fixedAssignmentFile = 'fix_ass.tab'
      self.file.files['fixAss'] = os.path.join(marsInputDir,fixedAssignmentFile)
      
      fixAssOut = open(os.path.join(self.inputDirectory,fixedAssignmentFile),'w')
      
      # Go over spins systems, and print out info if available
      ssInfo = []
      for spinSystem in self.spinSystemList:
        
        # Only report spin systems that are part of the chemical shift tab file!!
        if not self.spinSystemToPrevId.has_key(spinSystem) or not self.spinSystemToPrevId[spinSystem]:
          continue
        
        seqIdList = []
        resCodeList = []
        
        # First use 'clean' information
        if spinSystem.residue:
          if spinSystem.residue.chain == self.chain:
            seqIdList.append(spinSystem.residue.seqId)
          else:
            # Don't use if residue in wrong chain
            continue
                  
        else:
          if spinSystem.residueProbs:
            for residueProb in spinSystem.residueProbs:
              residue = residueProb.possibility
              if residue.chain == self.chain:
                seqIdList.append(residue.seqId)
              
            # Don't use if all residues are in wrong chain
            if not seqIdList:
              continue
        
          elif spinSystem.ccpCode:
            resCodeList.append(self.getCode1Letter(spinSystem.ccpCode))

          else:
            for residueTypeProb in spinSystem.residueTypeProbs:
              code1Letter = self.getCode1Letter(residueTypeProb.possibility.ccpCode)
              if not code1Letter in resCodeList:
                resCodeList.append(code1Letter)
            
        if seqIdList or resCodeList:
          ssInfo.append((spinSystem.serial,seqIdList,resCodeList))
      
      for (ssSerial,seqIdList,resCodeList) in ssInfo:
        fixAssOut.write("%-4d  " % ssSerial)
        if seqIdList:
          for seqId in seqIdList:
            fixAssOut.write("%d " % seqId)
        if resCodeList:
          fixAssOut.write("%s" % string.join(resCodeList,''))
        fixAssOut.write(newline)
 
      fixAssOut.close()
          
    #
    # Finish up
    #
              
    self.file.write()
    
    return self.file
  
    # TODO: secondary structure, PDB/RDC stuff...

    
  #
  # Functions different to default functions in DataFormat
  #
  
  def getCode1Letter(self,ccpCode):
    # Always protein
    return chemCompStdDict['protein'][ccpCode][0]

  def setSpecificAssignmentInfo(self):
    
    reliabilities = ['H','M','L']
    
    self.allowedReliabilities = reliabilities[0:reliabilities.index(self.reliabilityCutoff)]
  
  def assignmentValid(self,assignment):
    
    if assignment.reliability in self.allowedReliabilities:
      return True

  def setMeasurementsExportMode(self):
    
    # Set the export mode, either 'atoms' or 'spinSystems', depending on format.
    self.exportMode = 'spinSystems'

  def createSequence(self):
  
    self.sequenceFile.sequences.append(self.sequenceIO.MarsSequence(molName = self.chain.molecule.name))
    self.sequence = self.sequenceFile.sequences[-1]

  def setSequenceFileElements(self):
        
    self.sequence.elements.append(self.sequenceIO.MarsSequenceElement(self.seqCode,self.getResidueCode1Letter(self.residue)))

  def setSpecificSpinSystemInfo(self):
  
    self.spinSystemToPrevId = {}
    
    for resonanceToSpinSystem in self.atomMeasurements.keys():
      
      spinSystem = resonanceToSpinSystem.spinSystem
      
      if not self.spinSystemToPrevId.has_key(spinSystem):
        self.spinSystemToPrevId[spinSystem] = None
      
      if resonanceToSpinSystem.linkedSpinSystems:
        
        if resonanceToSpinSystem.linkedSpinSystems.has_key(-1):
          prevSpinSystem = resonanceToSpinSystem.linkedSpinSystems[-1]
          prevResonanceToSpinSystem = self.findResonanceToSpinSystem(prevSpinSystem)
        elif resonanceToSpinSystem.linkedSpinSystems.has_key(1):
          prevSpinSystem = spinSystem
          spinSystem = resonanceToSpinSystem.linkedSpinSystems[1]          
          prevResonanceToSpinSystem = resonanceToSpinSystem
          resonanceToSpinSystem = self.findResonanceToSpinSystem(spinSystem)
        else:
          continue

        self.spinSystemToPrevId[spinSystem] = prevSpinSystem.serial
        
        #
        # TODO HERE: if cs.tab is main file for identifing spin systems, then ALWAYS need to write them out here!
        # Also then check if spinSystem.serial in self.measurementFile for fix_ass and fix_conn file before writing anything.
        #
        
        # Only write out spin systems that are not in a 'chain' - this is a pretty rough way of doing things, but should work
        if prevResonanceToSpinSystem and resonanceToSpinSystem and (self.ignoreFixedConnection or (not resonanceToSpinSystem.linkedSpinSystems.has_key(1) and not prevResonanceToSpinSystem.linkedSpinSystems.has_key(-1))):
          self.measurementFile.spinSystemLinks[spinSystem.serial] = {-1: prevSpinSystem.serial}

  def findResonanceToSpinSystem(self,spinSystem):
    
    resonanceToSpinSystemFound = None
    
    for resonanceToSpinSystem in self.atomMeasurements.keys():
      if resonanceToSpinSystem.spinSystem == spinSystem:
        resonanceToSpinSystemFound = resonanceToSpinSystem
        break
    
    return resonanceToSpinSystemFound

  def setChemShiftFileValue(self):
  
    self.setMeasurementFileGenericValue()

  def setMeasurementFileGenericValue(self):

    if self.atomNamesDict.has_key(self.atomNames):
    
      resonanceToSpinSystem = self.atomNamesDict[self.atomNames]
      measurement = self.origAtomMeasurements[resonanceToSpinSystem]
      
      for atomName in self.atomNames:
        self.measurementFileValues.append(
        
          self.rawMeasurementClass(
          
            measurement.value,
            atomName,
            self.spinSystem.serial,
            None,
            None,
            prevSpinSystemId = self.spinSystemToPrevId[self.spinSystem]
            
          )
        )
      
      if self.atomMeasurements.has_key(resonanceToSpinSystem):
        del(self.atomMeasurements[resonanceToSpinSystem])

  def getPresetChainMapping(self,chainList):
  
    return self.getMultiChainFormatPresetChainMapping(chainList)
