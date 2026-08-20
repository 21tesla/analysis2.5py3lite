"""
======================COPYRIGHT/LICENSE START==========================

NmrStarFormat.py: Contains functions specific to NmrStar conversions.

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

import traceback, sys, os, time, string, copy, re

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccpnmr.format.general.Util import getResName
from ccpnmr.format.general.Util import getNameInfo
from ccpnmr.format.general.Util import getApplResNames
from ccpnmr.format.general.Constants import assign_kw

from memops.universal.Util import returnInt, returnFloat
from memops.general.Util import returnMemopsLine, returnMemopsWord

from ccp.format.general.Constants import defaultSeqInsertCode
from ccp.format.general.Util import getSeqAndInsertCode
from ccp.general.Util import getResonancesFromPairwiseConstraintItem

from ccp.general.Util import setCurrentStore

from ccpnmr.format.general.Constants import tagSep, defaultMolCode, allResidueAtoms_kw

from ccpnmr.format.general.Util import getOriginalData, setOriginalData, getOriginalDataValue

import memops.api.Implementation as Implementation
import ccp.api.nmr.Nmr as Nmr
import ccp.api.molecule.Molecule as Molecule
import ccp.api.molecule.ChemComp as ChemComp

# Specific for nmrStar!
import ccp.api.nmr.NmrEntry as NmrEntry

#
# Additional IOkeywords definitions
#

IOkeywords = copy.deepcopy(IOkeywords)

IOkeywords['readCoordinates']['disableAuthorSeqCode'] = (None,False,'If set, will not use author sequence codes when mapping coordinates (basically turns off use of author/PDB sequence codes for mapping)')
IOkeywords['readMeasurements']['readAllSaveFrames'] = (None,False,'If set, all saveframes containing measurement data will be read in into a joined chemical shift list.')
IOkeywords['writeProject']['entry'] = (None,True,'The entry to be exported to the NMR-STAR file.')
IOkeywords['writeProject']['useOriginalData'] = (False,False,'Use original data for export (if available).')
IOkeywords['writeProject']['exportAll'] = (False,False,'Export all information (if available).')
IOkeywords['writeMeasurements']['individualAtoms'] = (True,False,'If set to False (off), atom names will not be fully written out (e.g. an Ala HB* is written out instead of individual HB1, HB2, HB3 contributions)')

#
# TODO: nmrStar is a special case... still handling on same level
# as other formats for now.
# This means that many nmrStar scripts are superfluous though...
#

class NmrStarFormat(DataFormat):

  ccpn2NmrStarMolTypes = {
        
        'protein':      'polypeptide(L)',
        'DNA':          'polydeoxyribonucleotide',
        'RNA':          'polyribonucleotide',
        'DNA/RNA':      'DNA/RNA hybrid',
        'carbohydrate': 'carbohydrates',
        'other':        'other'
                        
                        }
                        
  chainIdToCcpn = {}
  starIdToCcpnObject = {}
  ccpnClassToStarId = {}

  count = 0

  def setFormat(self):
  
    self.format = 'nmrStar'
    self.IOkeywords = IOkeywords

  def getSequence(self):
  
    # TODO HERE: have to figure out what to do if project file read...

    nmrStarFile = self.sequenceIO.NmrStarFile(self.fileName)
    nmrStarFile.read()
    
    #
    # Some special gynmastics to make single access to sequences...
    #
    
    if nmrStarFile.sequenceFiles:
    
      self.sequenceFile = nmrStarFile.sequenceFiles[0]

      for sequenceFile in nmrStarFile.sequenceFiles[1:]:
        self.sequenceFile.dbLinks.extend(sequenceFile.dbLinks)
        for sequence in sequenceFile.sequences:
          self.sequenceFile.sequences.append(sequence)

    self.setCcpnMolTypes()
    
    if self.verbose == 1:
      print "Reading sequence from %s file %s" % (self.formatLabel,self.fileName)

  def setMoleculeInfo(self,molecule,createMoleculeInfo, newMolecule = True):

    #print 'CREATEMOL: [%s]' % createMoleculeInfo[1][0][0].dbLinks
    
    sequence = createMoleculeInfo[1][0][0]

    if hasattr(sequence,'details') and sequence.details:
      molecule.details = sequence.details

    if hasattr(sequence,'seqDetails') and sequence.seqDetails:
      molecule.seqDetails = sequence.seqDetails

    if hasattr(sequence,'EC_number') and sequence.EC_number:
      molecule.newMoleculeSysName(name=sequence.EC_number, namingSystem='EC')

    if hasattr(sequence,'Mutation') and sequence.Mutation:
      keywds = {'application': self.format,
                'keyword':     'Mutation',
                'value':       sequence.Mutation}

      appData = Implementation.AppDataString(**keywds)
      molecule.addApplicationData(appData)

    if hasattr(sequence,'Fragment') and sequence.Fragment:
      keywds = {'application': self.format,
                'keyword':     'Fragment',
                'value':       sequence.Fragment}

      appData = Implementation.AppDataString(**keywds)
      molecule.addApplicationData(appData)

    if hasattr(sequence,'Src_method') and sequence.Src_method:
      keywds = {'application': self.format,
                'keyword':     'Src_method',
                'value':       sequence.Src_method}

      appData = Implementation.AppDataString(**keywds)
      molecule.addApplicationData(appData)

    if hasattr(sequence,'authDefSeq') and sequence.authDefSeq:
      keywds = {'application': self.format,
                'keyword':     'authDefSeq',
                'value':       sequence.authDefSeq}

      appData = Implementation.AppDataString(**keywds)
      molecule.addApplicationData(appData)

    if hasattr(sequence,'ambChemCompSites') and sequence.ambChemCompSites:
      keywds = {'application': self.format,
                'keyword':     'ambChemCompSites',
                'value':       sequence.ambChemCompSites}

      appData = Implementation.AppDataString(**keywds)
      molecule.addApplicationData(appData)

    if hasattr(sequence, 'ambConfStates') and sequence.ambConfStates:
      keywds = {'application': self.format,
                'keyword':     'ambConfStates',
                'value':       sequence.ambConfStates}

      appData = Implementation.AppDataString(**keywds)
      molecule.addApplicationData(appData)

    for dbLink in sequence.dbLinks:
      dbName = dbLink.dbName
      
      # Can sometimes be empty - problem in NMR-STAR file, ignore
      if not dbName:
        print "  Warning: sequence DB link without name in original file."
        dbName = 'unknown'
      
      db = self.project.findFirstDatabase(name = dbName)

      if not db:
        db = self.project.newDatabase(name = dbName)

      dbAcc = dbLink.dbAcc
      dbMolName = dbLink.dbMolName
      dbDetails = dbLink.dbDetails

      if dbMolName:
        dbMolName = returnMemopsWord(dbMolName)
      else:
        dbMolName = dbAcc

      dbEntry = db.findFirstEntry(name = dbAcc)

      aln = None

      if not dbEntry:
        dbEntry = db.newEntry(name = dbAcc, code = dbMolName, details = dbDetails)

      else:
        aln = molecule.findFirstAlignment(dbRef = dbEntry)

      if dbEntry:
        keywds = {'application': self.format,
                  'keyword':     'dbName',
                  'value':       dbName}

        appData = Implementation.AppDataString(**keywds)
        dbEntry.addApplicationData(appData)

      if not aln:
        aln = molecule.newAlignment(dbRef = dbEntry)

    for bioFunction in sequence.bioFunctions:
      if bioFunction.bioFunc not in molecule.functions:
        if len(bioFunction.bioFunc) < 80:
          molecule.addFunction(bioFunction.bioFunc)
        else:
          bioFuncDetailsText = "Function: %s" % bioFunction.bioFunc
          if molecule.details:
            molecule.details += bioFuncDetailsText
          else:
            molecule.details = bioFuncDetailsText
            
    # Possible that not set, in case reading only molecule component...
    if hasattr(self,'starIdToCcpnObject'):
  
      foreignIdTagName = ('Entity','ID')
  
      if not self.starIdToCcpnObject.has_key(foreignIdTagName):
        self.starIdToCcpnObject[foreignIdTagName] = []
  
      (sequence,sequenceList) = createMoleculeInfo[1][0]
      self.starIdToCcpnObject[foreignIdTagName].append((sequence.entityId,molecule))

  def setMolSystemInfo(self):

    #print 'DICT: [%s] [%s]' % (dir(self.molSystem), dir(self.sequenceFile) )

    if self.sequenceFile.details:
      self.molSystem.details = self.sequenceFile.details

    if hasattr(self.sequenceFile,'EC_number') and self.sequenceFile.EC_number:
      self.molSystem.newMolSystemSysName(name=self.sequenceFile.EC_number, namingSystem='EC')

    for dbLink in self.sequenceFile.dbLinks:

      dbName = dbLink.dbName
      db = self.project.findFirstDatabase(name = dbName)

      if not db:
        db = self.project.newDatabase(name = dbName)

      dbAcc = dbLink.dbAcc

      dbEntry = db.findFirstEntry(name = dbAcc)

      if not dbEntry:
        dbEntry = db.newEntry(name = dbAcc, code = dbAcc)

        if dbLink.exptMethod:
          keywds = {'application': self.format,
                    'keyword':     'exptMethod',
                    'value':       dbLink.exptMethod}

          appData = Implementation.AppDataString(**keywds)
          dbEntry.addApplicationData(appData)

        if dbLink.structResolution:
          keywds = {'application': self.format,
                    'keyword':     'structResolution',
                    'value':       dbLink.structResolution}

          appData = Implementation.AppDataFloat(**keywds)
          dbEntry.addApplicationData(appData)

      self.molSystem.addDbReference(dbEntry)

    for bioFunction in self.sequenceFile.bioFunctions:
      if bioFunction.bioFunc not in self.molSystem.functions:
        self.molSystem.addFunction(bioFunction.bioFunc)

  def getMeasurements(self, sfnum = 0):
    
    #
    # Get a measurement list
    #
    
    try:
    
      self.nmrStarFile = self.measurementsIO.NmrStarFile(self.fileName)

      self.nmrStarFile.read()

      # Handle chemical shift file exception here. Not good.
      if self.measurementType == 'shift':
        measurementFilesAttr = "chemShiftFiles"
      elif self.measurementType == 'noe':
        measurementFilesAttr = "hetNoeFiles"
      else:
        measurementFilesAttr = "%sFiles" % self.measurementType
      
      measurementFiles = getattr(self.nmrStarFile,measurementFilesAttr)
            
      if measurementFiles:
 
        if self.verbose == 1:
          print "Reading %s list(s) from %s file %s" % (self.measurementType,self.formatLabel,self.fileName)

        if len(measurementFiles) > 1:
          if sfnum >= len(measurementFiles):
            sfnum = 0
          
          if self.readAllSaveFrames:
            self.measurementFile = measurementFiles[0]
            
            for otherSfNum in range(1,len(measurementFiles)):
              otherMeasurementFile = measurementFiles[otherSfNum]
              self.measurementFile.chemShifts.extend(otherMeasurementFile.chemShifts)
          
          else:
            self.messageReporter.showWarning("Warning","Multiple %s saveframes found - reading saveframe %d. Import whole project to get all saveframes!" % (self.measurementType,sfnum),self.guiParent)
            self.measurementFile = measurementFiles[sfnum]

        else:
          self.measurementFile = measurementFiles[sfnum]
        
      else:
      
        self.messageReporter.showWarning("Warning"," No %s in %s file %s..." % (self.formatFileValueKey,self.formatLabel,self.fileName),self.guiParent)
        self.measurementFile = None

    except:
      raise
      errorMessage = "        %s" % traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
    
      self.messageReporter.showWarning("Warning"," Cannot read %s for %s...:\n%s" % (self.formatFileValueKey,self.formatLabel,errorMessage),self.guiParent)
      self.measurementFile = None

      return traceback.format_exception(sys.exc_type,sys.exc_value,sys.exc_info()[2]) 

    # TODO HERE: have to figure out what to do if project file read...
  

    # WARNING: confusing names:
    #    seqCode in nmrStar is CCPN seqId
    #    authorSeqCode in nmrStar is CCPN seqCode!
    # Could also in principle check authorSeqCode with seqCode set for CCPN residues...

  def createMeasurementFile(self):
  
    try:
    
      self.nmrStarFile = self.measurementsIO.NmrStarFile(self.fileName)

      self.measurementFile = self.MeasurementFileClass(self.fileName,self.nmrStarFile)
     
      if self.verbose == 1:
        print "Writing %ss to %s file %s" % (self.measurementType,self.formatLabel,self.fileName)
  
    except:

      self.messageReporter.showWarning("Warning"," No write%s available for %s..." % (self.measurementApiCode,self.formatLabel),self.guiParent)

  def getConstraints(self, sfnum = 0):
    
    #
    # Format slightly different for constraints....
    # TODO: have to be able to read multiple ones!!
    #
    
    try:
    
      nmrStarFileClass = getattr(self.constraintsIO,"%sFile" % (self.formatLabel))
      self.nmrStarFile = nmrStarFileClass(self.fileName,version=self.version)
      self.nmrStarFile.read()
            
      if self.nmrStarFile.constraintFiles:
 
        if self.verbose == 1:
          print "Reading %s constraint list(s) from %s file %s" % (self.constraintType,self.formatLabel,self.fileName)

        if len(self.nmrStarFile.constraintFiles) > 1:
          if sfnum >= len(self.nmrStarFile.constraintFiles):
            sfnum = 0
          
          self.messageReporter.showWarning("Warning","Multiple %s constraint saveframes found - reading saveframe %d. Import whole project to get all saveframes!" % (self.constraintApiCode,sfnum),self.guiParent)

        self.constraintFile = self.nmrStarFile.constraintFiles[sfnum]
        
      else:
      
        self.messageReporter.showWarning("Warning"," No %s constraints in %s file %s..." % (self.constraintApiCode,self.formatLabel,self.fileName),self.guiParent)
        self.constraintFile = None

    except:

      errorMessage = traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
      self.messageReporter.showWarning("Warning"," Cannot read %s constraints for %s...:\n%s" % (self.constraintApiCode,self.formatLabel,errorMessage),self.guiParent)
      self.constraintFile = None
      
      return traceback.format_exception(sys.exc_type,sys.exc_value,sys.exc_info()[2]) 

  def createConstraintFile(self):
  
    try:
          
      if self.file:
        refFile = self.file
      else:
        refFile = self.refFile

      constraintFileClass = getattr(self.constraintsIO,"%sConstraintFile" % (self.formatCode))
      self.constraintFile = constraintFileClass(self.fileName,refFile)

      self.generalIO.NmrStarConstraintItem = getattr(self.constraintsIO,"%s%sConstraintItem" % (self.formatLabel,self.constraintApiCode))
      self.generalIO.NmrStarConstraintMember = getattr(self.constraintsIO,"%s%sConstraintMember" % (self.formatLabel,self.constraintApiCode))
     
      if self.verbose == 1:
        print "Writing %s constraints to %s file %s" % (self.constraintType,self.formatLabel,self.fileName)
  
    except:
      self.messageReporter.showWarning("Warning"," No write%sConstraints available for %s..." % (self.constraintApiCode,self.formatLabel),self.guiParent)

  def getCoordinates(self, noSequenceRead=False):
    
    self.coordinateFile = self.coordinatesIO.NmrStarFile(self.fileName)
    self.coordinateFile.maxModelNum = self.maxNum

    self.coordinateFile.read()
        
    if self.verbose == 1:
      print "Reading coordinates from %s file %s" % (self.formatLabel,self.fileName)
      
  def getCoordinateSeqCode(self,coordinate):
    
    seqCode = None

    # Bit dangerous, but necessary to handle joint NMR-STAR files - this is because the sequence codes are
    # handled by format, and because we're trying to preserve author info in the NMR-STAR file...
    # Can be turned off by 'disableAuthorSeqCode'.
    if not self.disableAuthorSeqCode:
      if hasattr(coordinate,'authorSeqCode') and coordinate.authorSeqCode != None:
        seqCode = returnInt(coordinate.authorSeqCode)
      elif hasattr(coordinate,'pdbSeqCode') and coordinate.pdbSeqCode != None:
        seqCode = returnInt(coordinate.pdbSeqCode)
        
      if seqCode != None:
        self.specialSeqCodeMapping = True
    
    if seqCode == None:
      seqCode = coordinate.seqCode
      
    return seqCode

  def createCoordinateFile(self):
  
    if self.verbose == 1:
      print "Writing coordinates to %s file %s" % (self.formatLabel,self.fileName)
  
    self.coordinateFile = self.coordinatesIO.NmrStarCoordinateFile(self.fileName,self.file)

  def getFullProject(self,fileName,entry = None, molSystem = None, strucGen = None, components = [], title = 'None', details = 'None', version = '2.1.1', **keywds):
  
    if self.verbose == 1:
      print "Reading %s project from file %s" % (self.formatLabel,fileName)
    
    # Little hack - should come from readProject
    if not self.version:
      self.version = version
   
    # TODO: this could possibly cause problems!!
    if not self.file:
      self.file = self.projectIO.NmrStarProjectFile(fileName,version = self.version)
    
    status = self.file.read(**keywds)

    if not status:
      if self.verbose:
        print "Invalid NMR-STAR project file, cannot read, aborting."
      return False
    
    #
    # Track mapping...
    #
    
    self.starIdToCcpnObject = {}

    self.ccpnClassToStarId = {}
    
    #
    # Extensions specifically for nmrStar reading.
    # Functions only inside this Format!!
    #
    
    #
    # Create an entry...
    #
    
    self.createEntry(entry,title,details,self.version)
   
    #
    # Process elements - call functions 'as if' reading real file, but pass
    # on objects created during project import
    #
    
    self.fileName = fileName
    
    #
    # Read molecule/molSystem information
    #
    
    chains = None
    
    if not components or ('sequence' in components):

      if self.file.sequenceFiles:
        self.sequenceFile = self.file.sequenceFiles[0]
        if len(self.file.sequenceFiles) > 1:
          for sequenceFile in self.file.sequenceFiles[1:]:
            for sequence in sequenceFile.sequences:
              self.sequenceFile.sequences.append(sequence)

        self.setCcpnMolTypes()

        localKeywds = {}

        for localKeywd in ('chemCompPath','chemCompPaths','chemCompMappings','ignoreUnknownChemComps','useOriginalChainCode', 'autoCreateChemComps', 'ignoreResNames'):
          if keywds.has_key(localKeywd):
            localKeywds[localKeywd] = keywds[localKeywd]

        chains = self.readSequence(molSystem = molSystem, minimalPrompts = self.minimalPrompts, version = self.version, **localKeywds)
        # TODO: need to set sequence tag stuff!!!

        self.entry.molSystem = self.molSystem
        
    #
    # Read coordinates - moved this after sequence to avoid mapping issues
    #
    
    if not components or ('coordinates' in components):
      
      localKeywds = {}
      
      for localKeywd in ('chemCompPath','chemCompPaths','linkAtoms','forceChainMappings','forceDefaultChainMapping','addNameMappings','ignoreUnknownChemComps','maxNum', 'autoCreateChemComps', 'ignoreResNames'):
        if keywds.has_key(localKeywd):
          localKeywds[localKeywd] = keywds[localKeywd]
      
      for self.coordinateFile in self.file.coordinateFiles:
        
        # Hardset the linking if sequence info also from project...
        if chains:
          self.setReadCoordinatesChainMapping(chains,localKeywds)

        if self.coordinateFile.modelCoordinates:

          self.readCoordinates(molSystem = molSystem, strucGen = strucGen, minimalPrompts = self.minimalPrompts, version = self.version, **localKeywds)

          if self.strucGen and self.strucGen not in self.entry.structureGenerations:
            self.entry.addStructureGeneration(self.strucGen)

    #
    # Read constraints
    #

    if not components or ('constraints' in components):

      for self.constraintFile in self.file.distanceConstraintFiles:
        self.readDistanceConstraints(strucGen = strucGen, minimalPrompts = self.minimalPrompts, version = self.version)
        self.setConstraintSfApplData()

        if self.strucGen not in self.entry.structureGenerations:
          self.entry.addStructureGeneration(self.strucGen)

      for self.constraintFile in self.file.hBondConstraintFiles:
        self.readHBondConstraints(strucGen = strucGen,minimalPrompts = self.minimalPrompts, version = self.version)
        self.setConstraintSfApplData()

        if self.strucGen not in self.entry.structureGenerations:
          self.entry.addStructureGeneration(self.strucGen)

      for self.constraintFile in self.file.dihedralConstraintFiles:
        self.readDihedralConstraints(strucGen = strucGen,minimalPrompts = self.minimalPrompts, version = self.version)
        self.setConstraintSfApplData()

        if self.strucGen not in self.entry.structureGenerations:
          self.entry.addStructureGeneration(self.strucGen)

      for self.constraintFile in self.file.rdcConstraintFiles:
        self.readRdcConstraints(strucGen = strucGen,minimalPrompts = self.minimalPrompts, version = self.version)
        self.setConstraintSfApplData()

        if self.strucGen not in self.entry.structureGenerations:
          self.entry.addStructureGeneration(self.strucGen)

    #
    # Read people/citations
    #
    
    if not components or ('peopleAndCitations' in components):
      for self.peopleAndCitationsFile in self.file.peopleAndCitationsFiles:
        
        # Ignore if empty - is possible!
        if not self.peopleAndCitationsFile.persons and not self.peopleAndCitationsFile.citations:
          continue
        
        (persons,citations) = self.readPeopleAndCitations(minimalPrompts = self.minimalPrompts)

        if self.peopleAndCitationsFile.personInGroups:
          self.setPersonInGroups(persons)

        if self.peopleAndCitationsFile.infoType == 'contactPersons':
          for person in persons:
            self.entry.addContactPerson(person)
        elif self.peopleAndCitationsFile.infoType == 'authors':
          for person in persons:
            self.entry.addAuthor(person)
        else:
          # Is citation in this case... only one per saveframe
          citation = citations[0]
          if self.peopleAndCitationsFile.isPrimary:
            self.entry.primaryCitation = citation
          else:
            self.entry.addOtherCitation(citation)

    #
    # Read project/entry level information, only when reading project for 3.1 and up files
    #
    # NOTE: this is completely NmrStar specific for now (02/2008 Wim). Some of this might have to go into generic FC class...
    #
    # NOTE: for the priority order in which to handle these, see ccp.format.nmrStar.projectIO!!!!
    #
    # TODO: make sure to link the information to the entry (or other CCPN objects) - either here
    # if already read in, or further down. Might want to rearrange order of below to reflect best
    # way to merge NMR-STAR data with setup of CCPN objects.
    #
    
    if not components:

      # Entry level information can be set here - links to self.entry!
      if self.file.entryInformation:
        self.setEntryInformation()

      # Studies - not really used by the BMRB???
      for self.studyFile in self.file.studyFiles:
        self.readStudies()

      # Instruments
      for self.instrumentFile in self.file.instrumentFiles:
        self.readInstruments()

      # Probes
      for self.probeFile in self.file.probeFiles:
        self.readProbes()

      # Software - TODO might want to pre-create a lot of this as reference data based on NMR-STAR dictionaries!
      for self.softwareFile in self.file.softwareFiles:
        self.readSoftware()

      # Methods
      for self.methodFile in self.file.methodFiles:
        self.readMethods()

      # Ensemble stats
      for self.confStatFile in self.file.confStatFiles:
        self.readConfStats()

      # Sample information
      for self.sampleFile in self.file.sampleFiles:
        self.readSamples()

      # Sample conditions
      for self.sampleConditionFile in self.file.sampleConditionFiles:
        self.readSampleConditions()
      
      # NMR experiments
      for self.nmrExperimentFile in self.file.nmrExperimentFiles:
        self.readNmrExperiments() # TODO should this link into self.readAcq/ProcPars?

        self.setNmrSpectrometerAppData()

      # Chemical shift referencing
      if not self.file.chemShiftRefFiles:
        keywds = {'application': self.format,
                  'keyword':     'chemShiftRefs',
                  'value':       'empty'}

        appData = Implementation.AppDataString(**keywds)
        self.entry.addApplicationData(appData)

      else:
        for self.chemShiftRefFile in self.file.chemShiftRefFiles:
          self.readChemShiftRefs()

      # Natural/experimental source
      for self.sourceFile in self.file.sourceFiles:
        self.readSourceFile()

    if not components or ('coordinates' in components):

      # TODO - this is a bit of a hack as it assumes there is only one structure ensemble in the star file
      #   and we then re-use 'self.coordinateFile' from above (probably the last one in the 'for' loop
      #   that contains self.readCoordinates above).

      # It sets that structure ensemble to the first structure generation it can find
      #   or the structure generation passed into getFullProject
      #   and then looks for method/software info. in the 'self.coordinateFile'
      #   to add.

      self.setStructureGenerationEnsemble(strucGen)

    #
    # Read measurement data (TODO: update for j couplings, ...)
    #

    if not components or ('measurements' in components):
      
      for (nmrStarFileListName,dataFormatReadFuncName) in (('chemShiftFiles','readShifts'),
                                                           ('hExchProtectionFiles','readHExchProtections'),
                                                           ('hExchRateFiles','readHExchRates'),
                                                           ('hetNoeFiles','readNoes'),
                                                           ('jCouplingFiles','readJCouplings'),
                                                           ('rdcFiles','readRdcs'),
                                                           #('shiftAnisotropyFiles','readShiftAnisotropy'),
                                                           #('orderParamFiles','readIsotropicS2s'),
                                                           ('t1RelaxFiles','readT1s'),
                                                           ('t1RhoRelaxFiles','readT1Rhos'),
                                                           ('t2RelaxFiles','readT2s') ):

        for self.measurementFile in getattr(self.file,nmrStarFileListName):

          readMeasurementsKeywds = {}

          if dataFormatReadFuncName == 'readShifts':
            for readShiftKeywd in ('autoMergeShifts',):
              if keywds.has_key(readShiftKeywd):
                readMeasurementsKeywds[readShiftKeywd] = keywds[readShiftKeywd]

          readFunc = getattr(self, dataFormatReadFuncName) 
          readFunc(minimalPrompts = self.minimalPrompts, version = self.version, **readMeasurementsKeywds)

          if self.measurementList:
            
            appData = Implementation.AppDataString(application=self.format,keyword='saveFrameTitle',value=self.measurementFile.saveFrame.title)
            self.measurementList.addApplicationData(appData)
            
            if hasattr(self.measurementFile, 'details') and self.measurementFile.details:
              self.measurementList.setDetails(self.measurementFile.details)
              
            if dataFormatReadFuncName == 'readShifts':
              tagNames = ('chemShiftErr1H','chemShiftErr13C','chemShiftErr15N','chemShiftErr31P','chemShiftErr2H','chemShiftErr19F')

              for tagName in tagNames:
                value = getattr(self.measurementFile,tagName)
                if value is None:
                  continue

                keywds = {'application': self.format,
                          'keyword':     tagName,
                          'value':       str(value)}

                appData = Implementation.AppDataString(**keywds)
                self.measurementList.addApplicationData(appData)

              if 'ChemicalShiftReference' in self.ccpnClassToStarId:
                for (chemShiftRefId,chemShiftRef) in self.starIdToCcpnObject[self.ccpnClassToStarId['ChemicalShiftReference'] ]:

                  if self.measurementList.experiments:
                    for nmrExpt in self.measurementList.experiments:
                     if self.measurementFile.chemShiftRefId == chemShiftRefId:
                        if chemShiftRef and chemShiftRef not in nmrExpt.shiftReferences:
                          nmrExpt.addShiftReference(chemShiftRef)


            if self.measurementList not in self.entry.measurementLists:
              self.entry.addMeasurementList(self.measurementList)

    #
    # TODO Here have to custom set links between CCPN objects based on tracking between NMR-STAR IDs and
    # actual CCPN objects (e.g. _Chem_shift_software.Software_ID, set link between MeasurementList and Software)
    #
    
    #self.setRemainingLinks()
    
    #
    # Avoid reading project again, also reset some values
    #
    
    self.file = None
    self.strucGen = None
    
    return self.entry

  def setStructureGenerationEnsemble(self, strucGen):

    # TODO - this is a bit messy, but means any emsembles in NMR-STAR files will get a link through to the entry.

    strucEns = self.project.currentStructureEnsemble or self.project.findFirstStructureEnsemble()

    if strucEns:
      if not strucGen:
        strucGen = self.project.currentNmrProject.findFirstStructureGeneration()

      if not strucGen:
        strucGen = self.project.currentNmrProject.newStructureGeneration(name='sgEnsemble', generationType='refinement')

      strucEns.setStructureGeneration(strucGen)

      appData = strucEns.findFirstApplicationData(application=self.format, keyword='sgDetails')

      if appData:
        strucGen.details = appData.value

    if strucGen:
      coordMethod = self.setCoordinatesMethod(strucGen)

      if coordMethod:
        strucGen.method = coordMethod

        if coordMethod.software and coordMethod.software not in self.entry.software:
          self.entry.addSoftware(coordMethod.software)

      if strucGen not in self.entry.structureGenerations:
        self.entry.addStructureGeneration(strucGen)

  def setCoordinatesMethod(self, strucGen):

    setCurrentStore(self.project,'MethodStore')

    coordMethod = None

    if not (hasattr(self, "coordinateFile") and self.coordinateFile):
      return coordMethod

    if 'Software' in self.ccpnClassToStarId:
      for (softId,software) in self.starIdToCcpnObject[self.ccpnClassToStarId['Software'] ]:
        if hasattr(self.coordinateFile, "structSofts") and self.coordinateFile.structSofts:
          for soft in self.coordinateFile.structSofts:
            if soft.softId == softId:
              if 'Method' in self.ccpnClassToStarId:
                for (methodId,method) in self.starIdToCcpnObject[self.ccpnClassToStarId['Method'] ]:
                  if soft.methodId == methodId:
                    method.setSoftware(software)

                    coordMethod = method
                    break
              else:
                methName = soft.methodLabel
                coordMethod = self.project.currentMethodStore.findFirstMethod(name=methName)
                if not coordMethod:
                  coordMethod = self.project.currentMethodStore.newMethod(name=meth.methodLabel)

                if not coordMethod.software:
                  coordMethod.setSoftware(software)

    if not coordMethod:
      if 'Method' in self.ccpnClassToStarId:
        for (methodId,method) in self.starIdToCcpnObject[self.ccpnClassToStarId['Method'] ]:
          if hasattr(self.coordinateFile, "structMeths") and self.coordinateFile.structMeths:
            for meth in self.coordinateFile.structMeths:
              if meth.methodLabel == method.name:
                soft = self.project.currentMethodStore.findFirstSoftware(method=method)
                if not soft:
                  soft = self.project.currentMethodStore.newSoftware(name=method.name,
                                                                     version='any',
                                                                     details=meth.methodDetails,
                                                                     methods=[method,])
                coordMethod = method
                break

      else:
        if hasattr(self.coordinateFile, "structMeths") and self.coordinateFile.structMeths:
          meth = self.coordinateFile.structMeths[0]
          coordMethod = self.project.currentMethodStore.findFirstMethod(name=meth.methodLabel)

          if not coordMethod:
            coordMethod = self.project.currentMethodStore.newMethod(name=meth.methodLabel)

          soft = self.project.currentMethodStore.findFirstSoftware(method=coordMethod)
          if not soft:
            soft = self.project.currentMethodStore.newSoftware(name=coordMethod.name,
                                                               version='any',
                                                               details=meth.methodDetails,
                                                               methods=[coordMethod,])

    return coordMethod
              
  def setPersonInGroups(self, persons):

    setCurrentStore(self.project,'AffiliationStore')

    for person in persons:

      keywds = {'application': self.format,
                'keyword':     'Serial'}

      appData = person.findFirstApplicationData(**keywds)

      serial = appData.value

      for persInGrp in self.peopleAndCitationsFile.personInGroups:
        if serial == persInGrp.serial:

          affStore = self.project.currentAffiliationStore

          name = persInGrp.dept

          if not name:
            name = persInGrp.addr2

            if not name:
              name = persInGrp.addr1

              if not name:
               name = 'Default'

          if name:
            name = name.replace('\n', ' ')

            if len(name) >= 80:
              name = name[:80]

          orgn = affStore.findFirstOrganisation(name = name)

          if not orgn:

            addresses = []

            if persInGrp.addr1:
              addresses.append(persInGrp.addr1)
            if persInGrp.addr2:
              addresses.append(persInGrp.addr2)
            if persInGrp.addr3:
              addresses.append(persInGrp.addr3)

            keywds2 = {'name':             returnMemopsLine(name),
                       'addresses':        addresses,
                       'city':             persInGrp.city,
                       'country':          persInGrp.country,
                       'postalCode':       persInGrp.postCode,
                       'organisationType': persInGrp.orgType}

            orgn = affStore.newOrganisation(**keywds2)

            if persInGrp.stateProv:
              keywds5 = {'application': self.format,
                         'keyword':     'State',
                         'value':       persInGrp.stateProv}

              appData = Implementation.AppDataString(**keywds5)
              orgn.addApplicationData(appData)

          group = orgn.findFirstGroup(name = name)

          if not group:
            group = orgn.newGroup(name = name)

          personInGroup = person.findFirstPersonInGroup(emailAddress = persInGrp.emailAddr)

          if not personInGroup:
            mailAddress = None

            if persInGrp.mailAddr:
              mailAddress = persInGrp.mailAddr
            elif addresses:
              mailAddress = ', '.join(addresses)

            keywds3 = {'emailAddress':   persInGrp.emailAddr,
                       'mailingAddress': mailAddress,
                       'position':       persInGrp.role,
                       'faxNumber':      persInGrp.faxNo,
                       'group':          group}

            if persInGrp.phoneNo:
              keywds3['phoneNumbers'] = [persInGrp.phoneNo]

            personInGroup = person.newPersonInGroup(**keywds3)

            if not group in self.entry.laboratories:
              self.entry.addLaboratory(group)

              keywds4 = {'application': self.format,
                         'keyword':     'PinG Lab',
                         'value':       group.name}

              appData = Implementation.AppDataString(**keywds4)
              group.addApplicationData(appData)

          if personInGroup:
            person.currentPersonInGroup = personInGroup

  #
  # setMeasurementList functions redefined from DataFormat.py
  #
  
  def setMeasurementListExperiments(self):

    self.measurementListKeywds['experiments'] = []

    if 'Experiment' in self.ccpnClassToStarId:
      for (nmrExptId,nmrExpt) in self.starIdToCcpnObject[self.ccpnClassToStarId['Experiment'] ]:
        if self.measurementFile.measureExpts:
          for expt in self.measurementFile.measureExpts:
            if expt.exptId == nmrExptId and nmrExpt not in self.measurementListKeywds['experiments']:
              self.measurementListKeywds['experiments'].append(nmrExpt)

  def setMeasurementListSoftware(self):

    setCurrentStore(self.project,'MethodStore')

    self.measurementListKeywds['method'] = None

    if 'Software' not in self.ccpnClassToStarId:
      return

    self.measurementListKeywds['method'] = []

    software = None

    if 'Software' in self.ccpnClassToStarId:
      for (softId,software) in self.starIdToCcpnObject[self.ccpnClassToStarId['Software'] ]:
        if self.measurementFile.measureSofts:
          for soft in self.measurementFile.measureSofts:
            if soft.softId == softId:
              if 'Method' in self.ccpnClassToStarId:
                method = None
                for (methodId,method) in self.starIdToCcpnObject[self.ccpnClassToStarId['Method'] ]:
                  if soft.methodId == methodId:
                    method.setSoftware(software)

                    if method not in self.measurementListKeywds['method']:
                      self.measurementListKeywds['method'].append(method)
              else:
                methName = soft.methodLabel
                if not methName or methName in ('$?', '?', '$.', '.', '$', ''):
                  methName = software.name
                else:
                  if methName[0] == '$':
                    methName = methName[1:]
                  if methName.count('_'):
                    methName = methName.replace('_', ' ')

                method = self.project.currentMethodStore.findFirstMethod(name = methName)
                if not method:
                  method = self.project.currentMethodStore.newMethod(name = methName)
                method.setSoftware(software)

              if method and method not in self.measurementListKeywds['method']:
                self.measurementListKeywds['method'].append(method)
                
    if len(self.measurementListKeywds['method']) > 0:
      self.measurementListKeywds['method'] = self.measurementListKeywds['method'][0]
    else:
      del self.measurementListKeywds['method']

  def setReadCoordinatesChainMapping(self,chains,localKeywds):

    if not localKeywds.has_key('forceChainMappings'):
      
      #
      # First check if any mappings are available...
      #
  
      self.getNameMapping()
  
      #
      # Only set forcing of chain stuff if no chain mappings exist!
      #
          
      forceChainMappings = []
          
      for chain in chains:
        if not self.nameMapping.findAllChainMappings(chainCode = chain.code):
          forceChainMappings.append([chain.code,chain.code,1,0])
          localKeywds['disableAuthorSeqCode'] = True
  
      localKeywds['forceChainMappings'] = forceChainMappings
  
  def createFullProject(self,fileName, headerComment = '', peakKeyWdList = None):
  
    if self.verbose == 1:
      print "Writing %s project to file %s" % (self.formatLabel,fileName)
   
    #
    # Writing ONLY version 3.0!!
    #
    
    self.file = self.projectIO.NmrStarProjectFile(fileName, version = '3.0')
    
    if not self.entry:
      self.selectEntry()
      
    self.starIDs = {}
    self.starLabels = {}
    
    self.constraintListIDs = {}  # NOTE: this is for writing out star when NO original file available!!
    
    self.createFileCharsSf(headerComment)
        
    self.createAssemblyAndEntitySf()

    self.createCommentsSf()

    strucGens = self.entry.structureGenerations
    
    if strucGens:
    
      for strucGen in strucGens:
      
        if strucGen.nmrConstraintStore:
      
          for constraintList in strucGen.nmrConstraintStore.constraintLists:

            self.constraintApiCode = string.replace(constraintList.className,'ConstraintList','')
            self.constraintType = string.lower(self.constraintApiCode[0]) + self.constraintApiCode[1:]

            #if self.constraintType != 'distance':
            #  continue

            writeFunc = getattr(self,"write%sConstraints" % self.constraintApiCode)
            writeFunc(fileName, constraintList = constraintList,noWrite = True, resetIOkeywords = False)

            createSfFunc = getattr(self,"create%sConstraintsSf" % self.constraintApiCode)
            createSfFunc(constraintList)

        #
        # Write out structures...
        #

        self.writeCoordinates(fileName,structures = list(strucGen.molStructures),noWrite = True, resetIOkeywords = False)
        self.createCoordinatesSf()

    self.file.write()
    
    return self.file

  def createCoordinatesSf(self):
  
    if not hasattr(self,'coordinateFile') or not self.coordinateFile:
      print "  Error: cannot create coordinate saveframe (no file object created)"
      return
  
    sfName = 'conformer_family_coord_set'
    
    #
    # Generic saveframe setup
    #
    
    if not self.starIDs.has_key(sfName):
    
      coordFamilyID = 1
    
    else:
  
      coordFamilyID = self.starIDs[sfName][-1][-1] + 1
      
    #
    # Setup saveframe
    #

    title = sfName + '_' + str(coordFamilyID)
    
    self.starIDs[sfName] = [[self.coordinateFile,coordFamilyID]]
  
    self.setupSaveFrame(sfName,title)
    self.setSfTag('ID',self.starIDs[sfName][-1][-1])

    #
    # Loop _Dist_constraint_tree
    #

    tableName = '_Atom_site'

    self.setupTable(tableName)
    
    # THESE KEYS ARE INTEGERS - leave as is!
    modelNumbers = self.coordinateFile.modelCoordinates.keys()
    modelNumbers.sort()
    
    for modelNumber in modelNumbers:
    
      coordinateID = 1
    
      for coordinate in self.coordinateFile.modelCoordinates[modelNumber]:

        self.setTableTag('Conformer_family_coord_set_ID',coordFamilyID)
        self.setTableTag('Model_ID',modelNumber)
        self.setTableTag('ID',coordinateID)

        chain = coordinate.chain

        for (tempChain,chainID) in self.starIDs['_Entity_assembly']:
          if chain == tempChain:
            self.setTableTag('Label_entity_assembly_ID',chainID)
            tempChain = None
            break

        if tempChain:
          print "  Error: label not found for chain %s!!!" % (chain.code)

        chainLabel = coordinate.chainCode
        moleculeNames = self.chainLabelMoleculeNames[chainLabel]

        for (tempMoleculeNames,moleculeID) in self.starIDs['entity']:

          if moleculeNames == tempMoleculeNames:
            self.setTableTag('Label_entity_ID',moleculeID)
            tempMoleculeNames = None
            break

        if tempMoleculeNames:
          print "  Error: label not found for molecules %s!!!" % (str(tempMoleculeNames))

        #
        # Now get the chemComp based on the sequence code
        #

        for (IDs,chemCompCode) in self.starIDs['_Entity_poly_seq']:
          (tempEntityID,tempSeqId,tempSeqCode) = IDs

          if moleculeID == tempEntityID and coordinate.seqCode == tempSeqCode:

            self.setTableTag('Label_comp_index_ID',tempSeqId)
            self.setTableTag('Label_comp_ID',chemCompCode)
            
            chemCompCode = None
            break

        if chemCompCode:
          print "  Error: label not found for molecule %s, seqCode %d!!!" % (moleculeID,coordinate.seqCode)


        self.setTableTag('Label_atom_ID',coordinate.atomName)

        #
        # Get original info (PROBLEMS if copied!!)
        #

        (origChainCode,origSeqCode,origSeqInsertCode,origResLabel,origAtomName) = ('','',defaultSeqInsertCode,None,None)

        if hasattr(coordinate,'origData') and coordinate.origData:

          (origChainCode,origSeqCode,origSeqInsertCode,origResLabel,origAtomName) = coordinate.origData

        #
        # Set other info
        #
        
        # HACK FOR CHAIN CODE TO NULL!
        if origChainCode == defaultMolCode:
          origChainCode = None
        
        # Print out empty origSeqCode as . (do not force string conversion)
        origSeqInsertCode = string.strip(origSeqInsertCode)
        if origSeqCode == None and not origSeqInsertCode:
          origSeqCodeValue = None
        else:
          if origSeqCode == None:
            origSeqCode = ''
          origSeqCodeValue = str(origSeqCode) + origSeqInsertCode      

        self.setTableTag('Auth_segment_code',origChainCode)
        self.setTableTag('Auth_seq_ID', origSeqCodeValue)
        self.setTableTag('Auth_comp_ID',origResLabel)
        self.setTableTag('Auth_atom_ID',origAtomName)

        self.setTableTag('Type_symbol',coordinate.atomType)
        self.setTableTag('Cartn_x',coordinate.x)
        self.setTableTag('Cartn_y',coordinate.y)
        self.setTableTag('Cartn_z',coordinate.z)
        self.setTableTag('PDB_extracted_Bfactor_col',coordinate.bFactor)
        
        coordinateID += 1
    
  def createFileCharsSf(self,headerComment):
    
    saveFrameName = 'file_characteristics'
  
    self.setupSaveFrame(saveFrameName,'global_MR_file_characteristics')
    
    if not headerComment:
    
      headerComment = """
#
# This NMR-STAR file was created by the CcpNmr FormatConverter software
# developed by the PDBe group at the European Bioinformatics Institute (EBI).
# See http://www.ccpn.ac.uk/ or http://www.ebi.ac.uk/pdbe/docs/NMR/main.html
# for more information on this software and the CCPN project.
#
    """
    
    self.saveFrame.comment = headerComment
      
    applData = self.entry.findFirstApplicationData(application = self.format, keyword = saveFrameName)
  
    if applData:
    
      sfKey = applData.value
      keyword = saveFrameName + '_tags_' + sfKey
      sfTagValues = eval(self.entry.findFirstApplicationData(application = self.format, keyword = keyword).value)
      
      for i in range(1,len(self.file.sfDict[saveFrameName].tags)):
        
        tagInfo = self.file.sfDict[saveFrameName].tags[i]
        tagName = tagInfo[0]

        tagValue = sfTagValues[i-1]
    
        self.setSfTag(tagName,tagValue)
      
      self.originalFileName = self.saveFrame.tags[self.prefix + 'MR_file_name']
      
      #
      # Overwrite conversion date
      #

      self.setSfTag('Conversion_date',time.strftime("%Y-%m-%d"))
      
      #
      # Write out loops from application data
      #
      
      self.setLoopsFromApplData(self.entry,saveFrameName,sfKey)

  def createCommentsSf(self):
  
    saveFrameName = 'MR_file_comment'
    
    applData = self.entry.findAllApplicationData(application = self.format, keyword = saveFrameName)
    
    for applDatum in applData:
    
      self.setSfFromApplData(applDatum,saveFrameName,self.entry)
        
  def createAssemblyAndEntitySf(self):
  
    molSystem = self.entry.molSystem
    
    if not molSystem:
      if self.project.molSystems:
        molSystem = self.project.sortedMolSystems()[0]
      else:
        return None
    
    entitySfName = 'entity'
    chemCompSfName = 'chem_comp'
    
    self.chains = list(molSystem.sortedChains())
    self.chainDict = {}
    self.getChainMapping(resetMapping = True)
    
    #
    # Set up some variables for export...
    #
    
    self.chainLabelChains = {}
    self.chainLabels = []
    self.chainLabelMolecules = {}
    self.chainLabelMoleculeNames = {}

    for chain in self.chains:
    
      if not self.chainDict.has_key(chain):
      
        continue
      
      [exportChainCode,exportFirstSeqCode] = self.chainDict[chain]
      
      if exportChainCode not in self.chainLabels:
        self.chainLabels.append(exportChainCode)
        self.chainLabelChains[exportChainCode] = []
        self.chainLabelMolecules[exportChainCode] = []
        self.chainLabelMoleculeNames[exportChainCode] = []
        
      #
      # Sort by exportFirstSeqCode
      #
            
      insertPos = len(self.chainLabelChains[exportChainCode])

      for i in range(0,len(self.chainLabelChains[exportChainCode])):
        
        tempChain = self.chainLabelChains[exportChainCode][i]
        
        if int(exportFirstSeqCode) < int(self.chainDict[tempChain][1]):
          insertPos = i
          break
      
      self.chainLabelChains[exportChainCode].insert(insertPos,chain)
      self.chainLabelMolecules[exportChainCode].insert(insertPos,chain.molecule)
      self.chainLabelMoleculeNames[exportChainCode].insert(insertPos,chain.molecule.name)
      
    #
    # Assembly saveframe
    #
    
    sfName = 'assembly'
    title = sfName + '_' + string.replace(molSystem.code,' ','_')
    
    self.starIDs[sfName] = [[molSystem,1]]
  
    self.setupSaveFrame(sfName,title)
    self.setSfTag('ID',self.starIDs[sfName][-1][-1])
    
    #
    # Hack to get name for molSystem
    #
    
    molSysName = molSystem.name  
    if not molSysName:
      molSysName = molSystem.code
      if not molSysName:
        molSysName = 'None'
        
    self.setSfTag('Name',molSysName)
  
    chainID = 1
    tableName = '_Entity_assembly' # TODO: use lowercase, no _ name?
    
    self.setupTable(tableName)
    
    self.starIDs[tableName] = []
    self.starLabels[tableName] = []
    self.starIDs[entitySfName] = []
    self.starLabels[entitySfName] = []
    
    tempEntityID = 0
    
    for chainLabel in self.chainLabels:
      for chain in self.chainLabelChains[chainLabel]:

        #
        # New chain for export
        #

        self.starIDs[tableName].append([chain,chainID])
        self.starLabels[tableName].append([chain,chainLabel])
      
      #
      # Look if nmrStar molecule/entity already exists
      #
      
      entityID = None
      entityLabel = None
      
      moleculeNames = self.chainLabelMoleculeNames[chainLabel]

      for [tempMoleculeNames,tempEntityID] in self.starIDs[entitySfName]:
        if moleculeNames == tempMoleculeNames:
          entityID = tempEntityID
          
          for [tempMoleculeNames,entityLabel] in self.starLabels[entitySfName]:
            if moleculeNames == tempMoleculeNames:
              break
       
      if not entityID:
        entityID = tempEntityID+1
        entityLabel = '$molecule_' + str(entityID)
        
        self.starIDs[entitySfName].append([moleculeNames,entityID]) # TODO: better to always have ID, then get label from there!
        self.starLabels[entitySfName].append([moleculeNames,entityLabel])
      
      self.setTableTag(string.capitalize(sfName) + '_ID',self.starIDs[sfName][-1][-1])
      self.setTableTag('ID',chainID)
      self.setTableTag('Entity_ID',entityID)
      #self.setTableTag('Entity_assembly_name',chainLabel)
      self.setTableTag('Entity_label',entityLabel)

      chainID += 1

    #
    # Entity saveframe
    #

    entityLabelsDone = []
    
    self.starIDs[chemCompSfName] = []
    self.starIDs['_Entity_poly_seq'] = []

    for chainLabel in self.chainLabels:
     
      moleculeNames = self.chainLabelMoleculeNames[chainLabel]

      for (tempMolNames,entityID) in self.starIDs[entitySfName]:

        if moleculeNames == tempMolNames:

          for (tempMolNames,entityLabel) in self.starLabels[entitySfName]:
            if moleculeNames == tempMolNames:
              break
              
          break
      
      #
      # Check if already printed.
      #

      if entityLabel in entityLabelsDone:
        continue
        
      else:
        entityLabelsDone.append(entityLabel)
        
      #
      # Set up saveframe
      #
      
      self.setupSaveFrame(entitySfName,entityLabel[1:])
      self.setSfTag('ID',entityID)
     
      #
      # Map polymer type
      #

      molecules = self.chainLabelMolecules[chainLabel]
      
      isPolymer = 0
      
      for molecule in molecules:

        if molecule.isStdLinear:
          isPolymer = 1
          break
          
      #
      # Basically using type of first encountered polymer chain... tough luck if
      # people map badly for output
      #
      
      if isPolymer:
      
        self.setSfTag('Type','polymer')
        
        #
        # TODO: is this mapping correct?!
        #
        
        if self.ccpn2NmrStarMolTypes.has_key(molecule.molType):
          molType = self.ccpn2NmrStarMolTypes[molecule.molType]
        else:
          molType = molecule.molType
        
        self.setSfTag('Pol_type',molType)
        self.setSfTag('Seq_length',molecule.seqLength)
        
        #
        # Get 1-letter code, set 'X' if non standard.
        #
        
        seqCode1LetterString = ""
        
        for molecule in molecules:
        
          for molRes in molecule.sortedMolResidues():
        
            chemComp = molRes.chemComp
          
            if isinstance(chemComp,ChemComp.StdChemComp):
          
              seqCode1LetterString += chemComp.code1Letter
          
            else:
          
              seqCode1LetterString += 'X'
        
        self.setSfTag('Seq',seqCode1LetterString)

      else:
      
        self.setSfTag('Type','non-polymer')


      #
      # _Entity_comp_index loop
      #
        
      tableName = '_Entity_comp_index'
      tableLoopID = 1  # TODO: use this for all parent loop IDs!!

      self.setupTable(tableName)
      
      # TODO: also set Num in dict?
      
      for molecule in molecules:

        for molRes in molecule.sortedMolResidues():
        
          chemComp = molRes.chemComp

          chemCompCode = None

          for (tempChemComp,tempChemCompCode) in self.starIDs[chemCompSfName]:
            if chemComp == tempChemComp:
              chemCompCode = tempChemCompCode
              break

          chemCompVar = molRes.chemCompVar
          cifCode = self.findChemCompVarSysName('CIF',chemCompVar)

          if not chemCompCode:
            self.starIDs[chemCompSfName].append([chemComp,cifCode])

          self.setTableTag(string.capitalize(entitySfName) + '_ID',entityID)
          self.setTableTag('Num',tableLoopID) # Essentially seqId
          self.setTableTag('Comp_ID',cifCode)

          tableLoopID += 1

      #
      # _Entity_poly_seq loop
      #
      
      tableName = '_Entity_poly_seq'
      tableLoopID = 1

      self.setupTable(tableName) 
 
      for i in range(0,len(molecules)):
       
        molecule = molecules[i]
        
        chain = self.chainLabelChains[chainLabel][i]
        [exportChainCode,seqCode] = self.chainDict[chain]

        ccIndex = 1
        
        for molRes in molecule.sortedMolResidues():

          chemComp = molRes.chemComp

          chemCompCode = None

          for (tempChemComp,tempChemCompCode) in self.starIDs[chemCompSfName]:
            if chemComp == tempChemComp:
              chemCompCode = tempChemCompCode
              break

          chemCompVar = molRes.chemCompVar
          cifCode = self.findChemCompVarSysName('CIF',chemCompVar)

          self.starIDs[tableName].append([[entityID,ccIndex,seqCode],cifCode])

          self.setTableTag(string.capitalize(entitySfName) + '_ID',entityID)
          self.setTableTag('Num',tableLoopID) # Refers to previous loop
          # Using molecule mapping...
          self.setTableTag('Comp_index_num',ccIndex)
          self.setTableTag('Mon_ID',cifCode)

          tableLoopID += 1
          ccIndex += 1
          seqCode += 1

    #
    # Chemcomp saveframe
    #

    for [chemComp,cifCode] in self.starIDs[chemCompSfName]:
      
      if isinstance(chemComp,ChemComp.NonStdChemComp):

        self.setupSaveFrame(chemCompSfName,cifCode)
        self.setSfTag('ID',cifCode)
        self.setSfTag('Type','non-polymer')

  def createHBondConstraintsSf(self,constraintList):
  
    self.createDistanceConstraintsSf(constraintList, subType = 'hydrogen bond')

  def createDistanceConstraintsSf(self,constraintList, subType = None):
  
    sfName = 'distance_constraints'
    
    if not self.constraintListIDs.has_key(sfName):
      self.constraintListIDs[sfName] = 1
    
    conversionErrors = []
    
    (constraintListID,sfKey) = self.setupConstraintSaveFrame(sfName,constraintList, subType = subType)
    
    #
    # Loop _Dist_constraint_tree
    #

    tableName = '_Dist_constraint_tree'

    #
    # Keep track of constraintIDs separately - need to be sequential
    #

    constraintID = 0
    constraintIDs = {}

    self.setupTable(tableName)
    
    #
    # Add application data conversion errors...
    #

    appDataConversionErrors = constraintList.findAllApplicationData(application = self.format, keyword = 'ConversionError')

    if appDataConversionErrors:
      for convError in appDataConversionErrors:
        (constraintId,errorType,errorText) = eval(convError.value)
        conversionErrors.append(['1','$' + self.saveFrame.title,constraintId,errorType,errorText])

    for constraint in self.constraintFile.constraints:      
      
      if constraint.errors:
      
        #
        # This can be made much more informative based on 
        # generic error messages from linkResonances!
        #
        # TODO: also list when there's items removed (because not
        # linked to atom) but other items available!!
        #
        
        constraint.errors.reverse()
        
        for errorText in constraint.errors:
        
          conversionErrors.append(['1','$' + self.saveFrame.title,constraint.Id,'1',errorText])
        
      #
      # Continue if no nodes...
      #
            
      if not constraint.nodes or not constraint.nodes[0].Id:
        
        if not constraint.errors:
          print "  nmrStar write error: removed constraint %d in list %d" % (constraint.Id,constraintListID)
        
        continue
      
      elif constraint.errors:
      
        continue
      
      else:
      
        constraintID += 1
        constraintIDs[constraint.Id] = constraintID
               
      #
      # Loop over the nodes
      #

      for node in constraint.nodes:
      
        self.setTableTag('Constraints_ID',constraintListID)
      
        self.setTableTag('ID',constraintIDs[constraint.Id])
      
        self.setTableTag('Node_ID',node.Id)
        self.setTableTag('Down_node_ID',node.downId)
        self.setTableTag('Right_node_ID',node.rightId)
        self.setTableTag('Logic_operation',node.logicOperation)
    
    #
    # Loop _Dist_constraint
    #

    tableName = '_Dist_constraint' # TODO: use lowercase, no _ name?

    self.setupTable(tableName)
    
    for constraint in self.constraintFile.constraints:
    
      if not constraintIDs.has_key(constraint.Id):
        continue
      
      for node in constraint.nodes:
      
        for member in node.members:
      
          self.setTableTag('Constraints_ID',constraintListID)
      
          self.setTableTag('Dist_constraint_tree_ID',constraintIDs[constraint.Id])    
          self.setTableTag('Tree_node_member_node_ID',node.Id)
          self.setTableTag('Contribution_fractional_val',member.fracVal)
          self.setTableTag('Constraint_tree_node_member_ID',member.Id)

          chainID = self.getEntityAssemblyID(member.chain)
          self.setTableTag('Label_entity_assembly_ID',chainID)
                  
          moleculeID = self.getEntityID(member.chainCode)
          self.setTableTag('Label_entity_ID',moleculeID)
          
          (compID,compIndex) = self.getCompIDAndCompIndex(member.seqCode,moleculeID)          
          self.setTableTag('Label_comp_index_ID',compIndex)
          self.setTableTag('Label_comp_ID',compID)
          
          self.setTableTag('Label_atom_ID',member.atomName)
          
          #
          # Get original info (PROBLEMS if copied!!)
          #

          (nmrStarCode,nmrStarSeqID,seqInsertCode,nmrStarAtomID,nmrStarLabel) = ('','',defaultSeqInsertCode,None,None)

          nmrStarName = member.nmrStarName
          
          if nmrStarName:
          
            (nmrStarCode,nmrStarSeqID,void,seqInsertCode,nmrStarAtomID) = getNameInfo(nmrStarName)
            
            if not nmrStarCode or nmrStarCode == defaultMolCode:
              nmrStarCode = None
              
            nmrStarLabel = member.nmrStarLabel
            
          else:
            print "  Error: constraint %d has missing original name" % (constraintIDs[constraint.Id])
          
          # TODO: additional setting of original stuff here?
          
          self.setTableTag('Auth_segment_code',nmrStarCode)
          self.setTableTag('Auth_seq_ID', str(nmrStarSeqID) + string.strip(seqInsertCode))
          self.setTableTag('Auth_comp_ID',nmrStarLabel)
          self.setTableTag('Auth_atom_ID',nmrStarAtomID)

    #
    # Loop _Dist_constraint_value
    #

    tableName = '_Dist_constraint_value'

    self.setupTable(tableName)
    
    attrNames = ['target','lowerBound','upperBound']
    attrDict = {}
    
    for constraint in self.constraintFile.constraints:
      
      if not constraintIDs.has_key(constraint.Id):
        continue

      for attrName in attrNames:
        attrDict[attrName] = getattr(constraint.nodes[0],attrName)
      
      for node in constraint.nodes:
      
        if not node.Id or node.logicOperation:
          continue

        for attrName in attrNames:
        
          if hasattr(node,attrName):

            value = getattr(node,attrName)
          
            if value:
              attrDict[attrName] = value
                  
        self.setTableTag('Constraints_ID',constraintListID)
        self.setTableTag('Constraint_ID',constraintIDs[constraint.Id])
        self.setTableTag('Tree_node_ID',node.Id)

        self.setTableTag('Source_experiment_ID',None)
        self.setTableTag('Spectral_peak_ID',None)
        self.setTableTag('Intensity_val',None)
        self.setTableTag('Intensity_lower_val_err',None)
        self.setTableTag('Intensity_upper_val_err',None)

        self.setTableTag('Distance_val',attrDict['target'])
        self.setTableTag('Distance_lower_bound_val',attrDict['lowerBound'])
        self.setTableTag('Distance_upper_bound_val',attrDict['upperBound'])

        if hasattr(node,'weight'):
          self.setTableTag('Weight',node.weight)
        else:
          self.setTableTag('Weight',None)

        self.setTableTag('Spectral_peak_ppm_1',None)
        self.setTableTag('Spectral_peak_ppm_2',None)

    #
    # Comments loop (only one that is saved in application data)
    #
    
    self.setLoopsFromApplData(constraintList,sfName,sfKey)

    #
    # Conversion errors loops (if any found)
    #
    
    self.setConversionErrorsLoop(sfName,constraintListID,conversionErrors)


  def getEntityAssemblyID(self,chain):

    if chain:
      for (tempChain,chainID) in self.starIDs['_Entity_assembly']:
        if chain == tempChain:
          tempChain = None
          break

      if tempChain:
        print "  Error: label not found for chain %s!!!" % (chain.code)

    else:
      chainID = ''
      
    return chainID
                  
  def getEntityID(self,chainLabel):
          
    if chainLabel:
      moleculeNames = self.chainLabelMoleculeNames[chainLabel]

      #molecule = chain.molecule        
      #self.starLabels[entitySfName].append([moleculeNames,entityLabel])
      for (tempMoleculeNames,moleculeID) in self.starIDs['entity']:

        if moleculeNames == tempMoleculeNames:
          tempMoleculeNames = None
          break

      if tempMoleculeNames:
        print "  Error: label not found for molecules %s!!!" % (str(tempMoleculeNames))

    else:

      moleculeID = ''
      
    return moleculeID

  def getCompIDAndCompIndex(self,seqCode,moleculeID):

    compIndex = compID = ''
    
    #
    # Set Label_comp_index_ID to export sequence code
    #

    if seqCode != -999999:

      #
      # Now get the chemComp based on the sequence code
      #

      for (IDs,chemCompCode) in self.starIDs['_Entity_poly_seq']:
        (tempEntityID,tempSeqId,tempSeqCode) = IDs

        if moleculeID == tempEntityID and seqCode == tempSeqCode:

          compID = chemCompCode
          chemCompCode = None
          compIndex = tempSeqId
          break

      if chemCompCode:
        print "  Error: label not found for molecule %s, seqCode %d!!!" % (moleculeID,seqCode)
        compId = compIndex = ''
    
    return (compID, compIndex)
          
  def createDihedralConstraintsSf(self,constraintList):
  
    sfName = 'torsion_angle_constraints'
    
    if not self.constraintListIDs.has_key(sfName):
      self.constraintListIDs[sfName] = 1

    conversionErrors = []
    
    (constraintListID,sfKey) = self.setupConstraintSaveFrame(sfName,constraintList)

    #
    # Loop _Torsion_angle_constraint
    #

    tableName = '_Torsion_angle_constraint' 

    self.setupTable(tableName)
    
    tableInfo = self.findTableInfo(sfName,tableName)
     
    tableLoopID = 1
    
    #
    # Add application data conversion errors...
    #

    appDataConversionErrors = constraintList.findAllApplicationData(application = self.format, keyword = 'ConversionError')

    if appDataConversionErrors:
      for convError in appDataConversionErrors:
        (constraintId,errorType,errorText) = eval(convError.value)
        conversionErrors.append(['1','$' + self.saveFrame.title,constraintId,errorType,errorText])

    for constraint in self.constraintFile.constraints:

      if constraint.errors:
      
        #
        # This can be made much more informative based on 
        # generic error messages from linkResonances!
        #
        # TODO: also list when there's items removed (because not
        # linked to atom) but other items available!!
        #
        
        constraint.errors.reverse()
        
        for errorText in constraint.errors:
        
          conversionErrors.append(['1','$' + self.saveFrame.title,constraint.Id,'1',errorText])
        
      #
      # Continue if no nodes...
      #
            
      if not constraint.nodes:
        
        if not constraint.errors:
          print "  nmrStar write error: removed constraint %d in list %d" % (constraint.Id,constraintListID)
        
        continue

      elif constraint.errors:
      
        continue


      for node in constraint.nodes:
      
        #
        # Keep track of tags for one line, then write out in same order as
        # in reference data...
        #
      
        tableTags = {}
      
        tableTags['Constraints_ID'] = constraintListID
        tableTags['ID'] = tableLoopID

        for i in range(0,len(node.members)):
          
          member = node.members[i]
          
          chainID = self.getEntityAssemblyID(member.chain)
          tableTags['Label_entity_assembly_ID_%d' % (i+1)] = chainID
            
          moleculeID = self.getEntityID(member.chainCode)
          tableTags['Label_entity_ID_%d' % (i+1)] = moleculeID
                    
          (compID,compIndex) = self.getCompIDAndCompIndex(member.seqCode,moleculeID)          
          tableTags['Label_comp_index_ID_%d' % (i+1)] = compIndex
          tableTags['Label_comp_ID_%d' % (i+1)] = compID
          
          tableTags['Label_atom_ID_%d' % (i+1)] = member.atomName
      
          #
          # Get original info
          #
          
          (nmrStarCode,nmrStarSeqID,seqInsertCode,nmrStarLabel,nmrStarAtomID) = ('','',defaultSeqInsertCode,None,None)
          
          if constraint.nmrStarNames:
          
            try:

              nmrStarName = constraint.nmrStarNames[i]
              (nmrStarCode,nmrStarSeqID,void,seqInsertCode,nmrStarAtomID) = getNameInfo(nmrStarName)
              
              if not nmrStarCode or nmrStarCode == defaultMolCode:
                nmrStarCode = None
                
              nmrStarLabel = constraint.nmrStarLabels[i]
              
            except:
            
              print "  Error: constraint %d has missing original name" % (constraint.Id)
          
          else:
            
            if hasattr(member,'origSeqCode'):
              nmrStarSeqID = member.origSeqCode
                      
            if hasattr(member,'origResLabel'):
              nmrStarLabel = member.origResLabel
          
          tableTags['Auth_segment_code_%d' % (i+1)] = nmrStarCode
          tableTags['Auth_seq_ID_%d' % (i+1)] = str(nmrStarSeqID) + string.strip(seqInsertCode)
          tableTags['Auth_comp_ID_%d' % (i+1)] = nmrStarLabel
          tableTags['Auth_atom_ID_%d' % (i+1)] = nmrStarAtomID
          
        #
        # Add angle and other information
        #
        
        tableTags['Torsion_angle_name'] = constraint.name
        tableTags['Angle_upper_bound_val'] = constraint.upperAngle
        tableTags['Angle_lower_bound_val'] = constraint.lowerAngle
        tableTags['Force_constant_value'] = constraint.energyCst
        tableTags['Potential_function_exponent'] = constraint.funcExp

        #
        # Set the information
        #

        for i in range(0,len(tableInfo)):
          
          tagName = tableInfo[i][0]
          
          if tableTags.has_key(tagName):

            value = tableTags[tagName]

          else:
            
            #
            # Ignore missing tag names?
            # Or is this an error?
            #
            
            continue
          
          self.setTableTag(tagName,value)
        
      tableLoopID += 1

    #
    # Comments loop (only one that is saved in application data)
    #
    
    self.setLoopsFromApplData(constraintList,sfName,sfKey)

    #
    # Conversion errors loops (if any found)
    #
    
    self.setConversionErrorsLoop(sfName,constraintListID,conversionErrors)

  def createRdcConstraintsSf(self,constraintList):
  
    sfName = 'residual_dipolar_couplings'
    
    if not self.constraintListIDs.has_key(sfName):
      self.constraintListIDs[sfName] = 1

    conversionErrors = []
    
    (constraintListID,sfKey) = self.setupConstraintSaveFrame(sfName,constraintList)

    #
    # Loop _RDC_constraint
    #

    tableName = '_RDC_constraint'

    self.setupTable(tableName)
    
    tableInfo = self.findTableInfo(sfName,tableName)
     
    tableLoopID = 1
        
    #
    # Add application data conversion errors...
    #

    appDataConversionErrors = constraintList.findAllApplicationData(application = self.format, keyword = 'ConversionError')

    if appDataConversionErrors:
      for convError in appDataConversionErrors:
        (constraintId,errorType,errorText) = eval(convError.value)
        conversionErrors.append(['1','$' + self.saveFrame.title,constraintId,errorType,errorText])

    for constraint in self.constraintFile.constraints:
      
      if constraint.errors:
      
        #
        # This can be made much more informative based on 
        # generic error messages from linkResonances!
        #
        # TODO: also list when there's items removed (because not
        # linked to atom) but other items available!!
        #
        
        constraint.errors.reverse()
        
        for errorText in constraint.errors:
        
          conversionErrors.append(['1','$' + self.saveFrame.title,constraint.Id,'1',errorText])
        
      #
      # Continue if no nodes...
      #
            
      if not constraint.nodes:
        
        if not constraint.errors:
          print "  nmrStar write error: removed constraint %d in list %d" % (constraint.Id,constraintListID)
        
        continue

      elif constraint.errors:
      
        continue


      for node in constraint.nodes:
      
        #
        # This should only be 1 node for RDCs...
        #
        # Keep track of tags for one line, then write out in same order as
        # in reference data...
        #
      
        tableTags = {}
      
        tableTags['Constraints_ID'] = constraintListID
        tableTags['ID'] = tableLoopID
       
        for i in range(0,len(node.members)):
          
          member = node.members[i]
          
          chainID = self.getEntityAssemblyID(member.chain)
          tableTags['Label_entity_assembly_ID_%d' % (i+1)] = chainID
            
          moleculeID = self.getEntityID(member.chainCode)
          tableTags['Label_entity_ID_%d' % (i+1)] = moleculeID
                    
          (compID,compIndex) = self.getCompIDAndCompIndex(member.seqCode,moleculeID)          
          tableTags['Label_comp_index_ID_%d' % (i+1)] = compIndex
          tableTags['Label_comp_ID_%d' % (i+1)] = compID
          
          tableTags['Label_atom_ID_%d' % (i+1)] = member.atomName
      
          #
          # Get original info (PROBLEMS if copied!!)
          # TODO: THIS BIT SIMILAR TO DISTANCE CONSTRAINTS! SHould
          # be separate function from above@!!!
          #
          # WARNING: is different because possible that no label for second member
          # (in case of DIANA files...)
          #
          
          (nmrStarCode,nmrStarSeqID,seqInsertCode,nmrStarAtomID,nmrStarLabel) = ('','',defaultSeqInsertCode,None,None)
          
          if member.nmrStarName:
          
            try:

              (nmrStarCode,nmrStarSeqID,void,seqInsertCode,nmrStarAtomID) = getNameInfo(member.nmrStarName)
              
              if not nmrStarCode or nmrStarCode == defaultMolCode:
                nmrStarCode = None
                
              nmrStarLabel = member.nmrStarLabel
              
            except:
            
              print "  Error: constraint %d has missing original name" % (constraint.Id)
          
          tableTags['Auth_segment_code_%d' % (i+1)] = nmrStarCode
          tableTags['Auth_seq_ID_%d' % (i+1)] = str(nmrStarSeqID) + string.strip(seqInsertCode)
          tableTags['Auth_comp_ID_%d' % (i+1)] = nmrStarLabel
          tableTags['Auth_atom_ID_%d' % (i+1)] = nmrStarAtomID

          
        #
        # Add other information
        #
        
        #tableTags['name'] = constraint.name
        tableTags['RDC_val'] = constraint.value
        tableTags['RDC_lower_bound'] = constraint.lowerValue
        tableTags['RDC_upper_bound'] = constraint.upperValue
        tableTags['RDC_val_err'] = constraint.error

        #
        # Set the information
        #

        for i in range(0,len(tableInfo)):
          
          tagName = tableInfo[i][0]
          
          if tableTags.has_key(tagName):

            value = tableTags[tagName]

          else:
            
            #
            # Ignore missing tag names?
            # Or is this an error?
            #
            
            continue
          
          self.setTableTag(tagName,value)
        
      tableLoopID += 1

    #
    # Comments loop (only one that is saved in application data)
    #
    
    self.setLoopsFromApplData(constraintList,sfName,sfKey)

    #
    # Conversion errors loops (if any found)
    #
    
    self.setConversionErrorsLoop(sfName,constraintListID,conversionErrors)

  def setupConstraintSaveFrame(self,sfName,constraintList,subType = None):
    
    #
    # Generic constraint saveframe setup
    #
    
    if not self.starIDs.has_key(sfName):
    
      constraintListID = self.constraintListIDs[sfName]
    
    else:
  
      constraintListID = self.starIDs[sfName][-1][-1] + 1

    #
    # Set the tags from input file info if available...
    #
    
    applData = constraintList.findFirstApplicationData(application = self.format, keyword = sfName)

    if applData:
      
      sfKey = applData.value
      self.setSfFromApplData(applData,sfName,constraintList)
      
      constraintListID = returnInt(self.saveFrame.tags[self.saveFrame.prefix + tagSep + 'ID'])
      
    else:
      
      print "  Warning: no original data for saveframe %s" % (sfName)
      
      title = sfName + '_' + str(constraintListID)
      self.setupSaveFrame(sfName,title)
      self.setSfTag('ID',constraintListID)
      
      if sfName == 'distance_constraints':
        constrType = 'distance'
      elif sfName == 'torsion_angle_constraints':
        constrType = "dihedral angle"
      elif sfName == 'residual_dipolar_couplings':
        constrType = "dipolar coupling"
      else:
        constrType = None
        
      self.setSfTag('Type',constrType)
      self.setSfTag('Subtype',subType)
      self.constraintListIDs[sfName] += 1
      sfKey = None
   
    self.starIDs[sfName] = [[constraintList,constraintListID]]
    
    return (constraintListID,sfKey)

  def setupSaveFrame(self,sfName,title):
  
    self.saveFrame = self.file.setupSaveFrame(sfName,title)
    
    self.prefix = self.saveFrame.prefix + tagSep
    
    self.setSfTag('Sf_category',sfName)
    
  def setSfTag(self,name,value):

    self.saveFrame.setTag(self.prefix + name, value)
    
  def setupTable(self,tableName):

    self.table = self.saveFrame.setupTable(tableName)
    
    self.prefix = tableName + tagSep
    
  def setTableTag(self,name,value):

    self.table.setTag(self.prefix + name, value)
    
  #
  # Functions different to default functions in DataFormat
  #
  
  def createMolecule(self,molName,createMoleculeInfo):
        
    molecule = self.project.newMolecule(name = molName)
    self.convertCount[self.mainCode][1]['molecule(s)'] += 1

    #print 'CREATE: [%s]' % createMoleculeInfo

    if hasattr(self.sequenceFile,'commonName'):
      molecule.longName = self.sequenceFile.commonName
        
    # TODO; abbreviation doesn't fit in data model...
    
    if hasattr(self.sequenceFile,'otherNames'):
      # TODO: this not correct... how are multiple names stored in star?
      for commonName in self.sequenceFile.otherNames:
        molecule.addCommonName(returnMemopsLine(commonName))
   
    return molecule
    
  def setChainInfo(self,chain):
  
    #
    # Use this to track ID later if mapping to this chain from within the NMR-STAR file...
    #
    
    setOriginalData(self.format,chain,self.sequence,'Id')

    if self.sequence.details:
      chain.details = self.sequence.details

    if self.sequence.role:
      chain.role = self.sequence.role

    if self.sequence.physState:
      chain.physicalState = self.sequence.physState

    if self.sequence.confIsomer:
      chain.conformationalIsomer = self.sequence.confIsomer

    if self.sequence.chemExchState:
      chain.chemExchangeState = self.sequence.chemExchState

    if self.sequence.magnEquivCode:
      chain.magnEquivalenceCode = self.sequence.magnEquivCode

    if self.sequence.waterMoleculeCompIndexId:
      setOriginalData(self.format,chain,self.sequence,'waterMoleculeCompIndexId')
      setOriginalData(self.format,chain,self.sequence,'originalChainCode')

    setCurrentStore(self.project,'NmrProject')
    appData = chain.molecule.findFirstApplicationData(application = self.format, keyword = 'ambConfStates', value = 'yes')

    if chain.molecule.findFirstApplicationData(application = self.format, keyword = 'ambConfStates', value = 'yes'):
      chainStateSet = self.project.currentNmrProject.newChainStateSet(chain = chain, stateSetType = 'conformation')

      chainStateA = chainStateSet.newChainState(name = 'unknown conformation A', details = 'From NMRSTAR file')
      chainStateB = chainStateSet.newChainState(name = 'unknown conformation B', details = 'From NMRSTAR file')

      #chain.molecule.removeApplicationData(appData) # Can't do as all chains for this molecule need to be tagged.

    #
    # Also track star ID to CCPN object for chain - is entity assembly level
    #
    
    if hasattr(self,'starIdToCcpnObject'):
      
      foreignIdTagName = ('Entity_assembly','ID')

      if not self.starIdToCcpnObject.has_key(foreignIdTagName):
        self.starIdToCcpnObject[foreignIdTagName] = []
  
      self.starIdToCcpnObject[foreignIdTagName].append((self.sequence.Id,chain))

  def getChainOrRefChainId(self,coordOrChain,isChain = False):
  
    chainCode = None
    
    if isChain:

      for chain in self.molSystem.sortedChains():
        
        starChainId = getOriginalDataValue(self.format,chain,'Id')
        coordOrChainId = coordOrChain.chainId
        tmpChainCode = chain.code
                 
        if starChainId == coordOrChainId:          
          self.chainIdToCcpn[starChainId] = tmpChainCode
          coordOrChain.refChainId = tmpChainCode 
          chainCode = tmpChainCode
          
          # Water is handled differently... need to use original chain code.
          # TODO WARNING: this might possibly mess up things if there is overlap with the CCPN chain codes, but should
          # not be the case... .
          if coordOrChain.isWaterChain:
          
            waterMoleculeCompIndexId = getOriginalDataValue(self.format,chain,'waterMoleculeCompIndexId')
            originalChainCode = getOriginalDataValue(self.format,chain,'originalChainCode') # This is the asymId in entity_assembly!
    
            starChainId = (starChainId,waterMoleculeCompIndexId)
            coordOrChainId = (coordOrChainId,coordOrChain.firstIndexId)
            chainCode = originalChainCode
          
          break
   
    else:
      entityAssemblyId = coordOrChain.entityAssemblyId
      if coordOrChain.compId == 'HOH':
        coordOrChainId = (entityAssemblyId,coordOrChain.indexId)
      else:
        coordOrChainId = entityAssemblyId
      
      if entityAssemblyId in self.chainIdToCcpn:
        chainCode = self.chainIdToCcpn[entityAssemblyId]
    
    return chainCode

  def createShift(self,resonance,chemShift):
    
    keywds = {}
    
    #
    # For nmrStar also set ambiguity code if listed
    #
    if hasattr(chemShift,'ambCode'):
      if chemShift.ambCode != None:
        resonance.addApplicationData(Implementation.AppDataInt(application = self.format, keyword = 'ambiguityCode', value = returnInt(chemShift.ambCode)))
    
    #
    # Shift error could be missing...
    #
    
    if hasattr(chemShift,'valueError'):

      keywds['error'] = chemShift.valueError
      
    #
    # Create shift linked to resonance
    #
    
    if chemShift.value != None:
      shift = self.measurementList.newShift(value = chemShift.value, resonance = resonance, **keywds)
      self.convertCount[self.mainCode][1][self.valueCode] += 1
    else:
      shift = None

    return shift
    
  def setChemShiftFileValue(self):
        
    if self.atomNamesDict.has_key(self.atomName):
    
      resonanceToAtom = self.atomNamesDict[self.atomName]
      chemShift = self.origAtomMeasurements[resonanceToAtom]
    
      self.measurementFileValues.append(self.rawMeasurementClass(defaultMolCode,self.measurementFile))

      rawChemShift = self.measurementFile.chemShifts[-1]
      rawChemShift.Id = self.atomSerial
      rawChemShift.seqCode = self.seqCode
      rawChemShift.resLabel = self.residue.ccpCode
      rawChemShift.atomName = self.atomName
      rawChemShift.atomType = resonanceToAtom.resonance.resonanceSet.findFirstAtomSet().findFirstAtom().chemAtom.elementSymbol
      rawChemShift.value = chemShift.value
      rawChemShift.valueError = chemShift.error
      
      #
      # Not entirely correct but attempt at least...
      # Need to do more complicated cases as well...
      #
      
      rawChemShift.ambCode = self.getShiftAmbiguityCode(chemShift.value,resonanceToAtom,self.measurementList)
      
      #rawChemShift.authorSeqCode = 
      
      if self.atomMeasurements.has_key(resonanceToAtom):
        del(self.atomMeasurements[resonanceToAtom])

  def setShiftListSpecific(self):
    
    #
    # Track the _Systematic_chem_shift_offset loop, if present
    #
    
    tableName = '_Systematic_chem_shift_offset'
    
    if self.measurementFile.saveFrame.tables.has_key(tableName):

      sfKey = self.measurementFile.saveFrame.tags['ID']
      saveFrameName = self.measurementFile.saveFrame.tags['Sf_category']

      self.setSfTablesApplData(self,self.measurementList,saveFrameName,sfKey,tableList = [tableName])
 
  def setT1ListSpecific(self):

    if self.coherenceType:
      self.measurementListKeywds['coherenceType'] = self.coherenceType
    elif hasattr(self.measurementFile,'coherenceType'):
      
      nmrStarCoherenceType = self.measurementFile.coherenceType
      ccpnCoherenceType = ""
      elementSymbols = ""
      slen = 0
      
      while slen < len(nmrStarCoherenceType):
        elementSymbols    += nmrStarCoherenceType[slen]
        ccpnCoherenceType += nmrStarCoherenceType[slen + 1]
        slen += 2  
      
      self.measurementListKeywds['coherenceType'] = ccpnCoherenceType

      print "  Warning: ignoring element symbols for T1 list! Are %s" % elementSymbols
      
    if self.specFreq:
      self.measurementListKeywds['sf'] = self.specFreq
    elif hasattr(self.measurementFile,'specFreq'):
      self.measurementListKeywds['sf'] = self.measurementFile.specFreq    

  def setT2ListSpecific(self):

    """
    if self.coherenceType:
      self.measurementListKeywds['coherenceType'] = self.coherenceType
    elif hasattr(self.measurementFile,'coherenceType'):
      
      nmrStarCoherenceType = self.measurementFile.coherenceType
      ccpnCoherenceType = ""
      elementSymbols = ""
      slen = 0
      
      #while slen < len(nmrStarCoherenceType):
      #  elementSymbols    += nmrStarCoherenceType[slen]
      #  ccpnCoherenceType += nmrStarCoherenceType[slen + 1]
      #  slen += 2  
      
      #self.measurementListKeywds['coherenceType'] = ccpnCoherenceType

      self.measurementListKeywds['coherenceType'] = nmrStarCoherenceType
    """

    if self.specFreq:
      self.measurementListKeywds['sf'] = self.specFreq
    elif hasattr(self.measurementFile,'specFreq'):
      self.measurementListKeywds['sf'] = self.measurementFile.specFreq

    if self.tempCalibMethod:
      self.measurementListKeywds['tempCalibMethod'] = self.tempCalibMethod
    elif hasattr(self.measurementFile,'tempCalibMethod'):
      self.measurementListKeywds['tempCalibMethod'] = self.measurementFile.tempCalibMethod

    if self.tempControlMethod:
      self.measurementListKeywds['tempControlMethod'] = self.tempControlMethod
    elif hasattr(self.measurementFile,'tempControlMethod'):
      self.measurementListKeywds['tempControlMethod'] = self.measurementFile.tempControlMethod

  def setNoeListSpecific(self):

    if hasattr(self.measurementFile,'refValueType') and self.measurementFile.refValueType:
     refValueType = self.measurementFile.refValueType

     if refValueType == 'peak_height':
       self.measurementListKeywds['noeValueType'] = 'height'
     elif refValueType == 'peak_integral':
       self.measurementListKeywds['noeValueType'] = 'volume'
     elif refValueType == 'contour_count':
       self.measurementListKeywds['noeValueType'] = 'contour count'
     elif refValueType == 'relative_intensities':
       pass # Not an option in CCPN

    if hasattr(self.measurementFile,'refValue') and self.measurementFile.refValue:
     self.measurementListKeywds['refValue'] = self.measurementFile.refValue
      
    if hasattr(self.measurementFile,'refDesc') and self.measurementFile.refDesc:
     self.measurementListKeywds['refDescription'] = self.measurementFile.refDesc

    if hasattr(self.measurementFile,'specFreq') and self.measurementFile.specFreq:
      self.measurementListKeywds['sf'] = self.measurementFile.specFreq

  def getShiftAmbiguityCode(self,chemShiftValue,resonanceToAtom,shiftList):
        
    ambCode = 1
    
    chemAtom = resonanceToAtom.chemAtom
    
    # Current sensitivity, could be changed, but to what?
    chemShiftValueRef = int(chemShiftValue * 1000)
    
    # Check prochiral and aromatic chemAtoms
    if chemAtom:
      chemAtomSet = chemAtom.chemAtomSet
      if chemAtomSet:
        # Get info for 'deep' chemAtomSets
        if chemAtomSet.isEquivalent and chemAtomSet.chemAtomSet:
          chemAtomSet = chemAtomSet.chemAtomSet
      
        if chemAtomSet.isProchiral and resonanceToAtom.otherGroupResonances:
          refResonance = resonanceToAtom.resonance
      
          for resonance in resonanceToAtom.otherGroupResonances:          
          
            shift = resonance.findFirstShift(parentList = shiftList)
        
            if shift and int(shift.value * 1000) != chemShiftValueRef:
              ambCode = 2
              break
          
        elif chemAtomSet.isProchiral == False and chemAtomSet.isEquivalent == None:
          # Should be aromatic!
          ambCode = 3
          
    return ambCode

  def rawConstraintValid(self):
    
    #
    # Note: this has to work on ALL types of constraints!
    #
    
    if not self.rawConstraint.nodes:
      print "  Error: no nodes for constraint!"
      return 0

    if len(self.rawConstraint.nodes) > 1:
      
      #
      # Multiple nodes: need to work with downnode
      #
      
      firstNode = self.rawConstraint.nodes[0]
      logicOperator = firstNode.logicOperation
      #downNode = firstNode.downId
      
      if logicOperator != 'OR':
        print "  Currently handling only 'OR' noded nmrStar distance constraints..."
        return 0
      
      for node in self.rawConstraint.nodes[1:]:
        if node.downId != None:
          print "  Currently handling only simple 'OR' relationships..."
          return 0
          
      self.rawConstraint.items = self.rawConstraint.nodes[1:]
      
    else:
    
      self.rawConstraint.items = self.rawConstraint.nodes[:]
    
    # TODO: This is a hack...
    self.rawConstraintNode = self.rawConstraint.nodes[-1]

    return 1
    
  def getConstraintItems(self):
  
    #
    # Have to do a hack to get specific code in for distances/hbonds
    #
    
    if self.constraintType in ['distance','hBond']:
      self.getDistanceConstraintItems()
    else:
      DataFormat.getConstraintItems(self)

  def getDistanceConstraintItems(self):
  
    self.items = []
    self.rawItems = []
    
    for item in self.rawConstraint.items:
      
      #
      # For nmrStar multiple members with same Id mean OR relationship
      #
      
      memberDict = {}
      
      for member in item.members:
        
        if not memberDict.has_key(member.Id):

          memberDict[member.Id] = []
          
        memberDict[member.Id].append(member)
        
      
      if len(memberDict) != self.constraintFile.constraintElements:
        print "  Warning can only handle exactly 2 member Ids for one node - ignored node"
        continue
      
      
      memberIds = memberDict.keys()
      members = ['','']
      
      for members[0] in memberDict[memberIds[0]]:
        for members[1] in memberDict[memberIds[1]]:
          
          self.items.append([members[0],members[1]])
          self.rawItems.append(item)
             
          # Creating resonance(s)... these should be cleaned up later
          # Name of the resonance is format <chainCode>.<seqCode>.<atomName>
          #      or in CCPN                 <chain.code>.<residue.seqCode>.<atom.name>
          
          for i in range(0,2):
            member = members[i]
            
            chainCode = member.chainCode

            # Set the resonance name
            resName = getResName(chainCode,member.seqCode,member.atomName,seqInsertCode = member.seqInsertCode)
        
            # Get the resonance (if any)
            # Quickly set a link in the rawNodeMember class to the resonance(s)... not good style but easy!
            if self.resonanceNames.has_key(resName):
              member.resonances = self.resonanceNames[resName]
            else:
              # Create new resonance
              resonance = self.createResonance(resName, fixed = True, origSource = member)
              self.resonanceNames[resName] = [resonance]      
              member.resonances = [resonance]
    
  def setConstraintItemInfo(self):
    
    getItemNumFunc = getattr(self,"get%sConstraintItemNumResonances" % self.constraintApiCode)
    numResonances = getItemNumFunc()
    
    if numResonances == 2:
      itemResonances = getResonancesFromPairwiseConstraintItem(self.constraintItem)
    else:
      itemResonances = self.constraintItem.sortedResonances()
    
    for i in range(0,numResonances):
    
      resonance = itemResonances[i]
      
      applData = resonance.findFirstApplicationData(application = self.format, keyword = assign_kw)
      
      if applData:
        self.constraintItem.addApplicationData(Implementation.AppDataString(application = self.format, keyword = 'origAssign', value = applData.value))
        
      setOriginalData(self.format,self.constraintItem,self.rawConstraintItem.members[i],'resLabel',isUnique = False)
    
    setOriginalData(self.format,self.constraintItem,self.rawConstraintItem,'Id',isUnique = False)

  def getDistanceConstraintItemNumResonances(self):
  
    return len(self.constraintItem.resonances)
    
  def getHBondConstraintItemNumResonances(self):
  
    return self.getDistanceConstraintItemNumResonances()

  def getRdcConstraintItemNumResonances(self):
  
    numResonances = len(self.constraintItem.resonances)
  
    if self.rawConstraint.code:
    
      return (numResonances - 1)
  
    else:

      return numResonances

  def setDistanceConstraintInfo(self):
  
    #
    # TODO: make this optional?!?
    #
    
    setOriginalData(self.format,self.constraint,self.rawConstraintNode,'weight')
    setOriginalData(self.format,self.constraint,self.rawConstraint,'errors')  

  def setDihedralConstraintInfo(self):
  
    #
    # Need separate treatment if angle name given... TODO WARNING this
    # only works for DYANA type stuff!! Should 'decompose' the angle
    # name only on the DataFormat level and treat angle names like
    # atom names (get naming system, ...)?
    #
    
    self.constraint.details = "Original Id %d" % self.rawConstraint.Id
    
    if self.rawConstraint.name:
    
      setOriginalData(self.format,self.constraint,self.rawConstraint,'name')
      
      #
      # Have to loop to get the original seqCode, resLabel: possible
      # that not set for a residue!
      #
      
      for refMember in self.rawConstraint.nodes[0].members:
        
        if refMember.resLabel:

          setOriginalData(self.format,self.constraint,refMember,'seqCode')
          setOriginalData(self.format,self.constraint,refMember,'resLabel')
          break
    
    else:
      
      constraintResonances = self.constraint.resonances
    
      for i in range(0,len(constraintResonances)):
  
        resonance = constraintResonances[i]
      
        applData = resonance.findFirstApplicationData(application = self.format, keyword = assign_kw)
        
        if applData:
          self.constraint.addApplicationData(Implementation.AppDataString(application = self.format, keyword = 'origAssign', value = applData.value))
        
        setOriginalData(self.format,self.constraint,self.rawConstraint.nodes[0].members[i],'resLabel',isUnique = False)
    
    #
    # TODO: make this optional?!?
    #

    setOriginalData(self.format,self.constraint,self.rawConstraint,'energyCst')
    setOriginalData(self.format,self.constraint,self.rawConstraint,'funcExp')
    
  def setDihedralConstraintItemInfo(self):
  
    setOriginalData(self.format,self.constraintItem,self.rawConstraint,'errors') 
  
  def setRdcConstraintInfo(self):
  
    setOriginalData(self.format,self.constraint,self.rawConstraint,'errors')  
    setOriginalData(self.format,self.constraint,self.rawConstraint,'Id',isUnique = False)
    
  def setJCouplingConstraintInfo(self):
    
    setOriginalData(self.format,self.constraint,self.rawConstraint,'errors')  
    setOriginalData(self.format,self.constraint,self.rawConstraint,'Id',isUnique = False)
    
  def setDistanceConstraintLimits(self):
    
    self.constraintPars['upperLimit'] = self.rawConstraintNode.upperBound
    self.constraintPars['lowerLimit'] = self.rawConstraintNode.lowerBound
    self.constraintPars['targetValue'] = self.rawConstraintNode.target
    
    if hasattr(self.rawConstraintNode,'intensity'):
      value = getattr(self.rawConstraintNode,'intensity')
      if value != None:
        self.constraintPars['origData'] = value

  def setRawDistanceConstraint(self):
    
    self.constraintFile.constraints.append(self.rawConstraintClass(Id = self.constraint.serial))

    self.rawConstraint = self.constraintFile.constraints[-1]
    
    #
    # Have to create first node and set data there...
    #
    
    self.rawConstraint.addEmptyNode()
    node = self.rawConstraint.nodes[-1]
    
    node.target = self.constraint.targetValue
    node.lowerBound = self.constraint.lowerLimit
    node.upperBound = self.constraint.upperLimit
    node.weight = self.constraint.weight
    
  def setRawDistanceConstraintItemInfo(self):
  
    if self.useOriginalData:
    
      getOriginalData(self.format,self.constraint,self.rawConstraint.nodes[-1],'weight')
    
  def setRawDistanceConstraintInfo(self):
  
    if self.useOriginalData:
    
      getOriginalData(self.format,self.constraint,self.rawConstraint,'errors')

  def setRawDistanceConstraintItem(self):
      
    node = self.rawConstraint.nodes[-1]
  
    if node.Id:
    
      if node.Id == 1:
      
        # Add a node, move members from first node to this one
        # and set logical operation and down_node to this...
        
        self.rawConstraint.addEmptyNode(Id = 2)
        self.rawConstraint.nodes[1].members = self.rawConstraint.nodes[0].members[:]
        self.rawConstraint.nodes[0].members = []
        self.rawConstraint.nodes[0].downId = 2
        self.rawConstraint.nodes[0].logicOperation = 'OR'
        node = self.rawConstraint.nodes[1]
      
      newId = node.Id + 1
      self.rawConstraint.nodes[-1].rightId = newId
        
      self.rawConstraint.addEmptyNode(Id = newId)
      self.rawConstraintItem = self.rawConstraint.nodes[-1]
    
    else:
    
      node.Id = 1
      self.rawConstraintItem = node
    
  def setRawDistanceConstraintItemMembers(self):
  
    nmrStarNames = self.item.findAllApplicationData(application = self.format,keyword = 'origAssign')
    nmrStarLabels = self.item.findAllApplicationData(application = self.format,keyword = 'origResLabel')
    
    if nmrStarNames and len(nmrStarNames) != 2:
      print "  Error: constraint %d has invalid number of names %s" % (self.rawConstraint.Id,str(nmrStarNames))
      return
      
    for i in range(0,2):
            
      (chainCode,seqCode,void,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])

      if atomName == allResidueAtoms_kw:
        atomName = None
      
      self.rawConstraintItem.addEmptyMember(Id = i + 1)
      member = self.rawConstraintItem.members[-1]
      
      member.chainCode = chainCode
      member.seqCode = seqCode
      member.atomName = atomName
      
      # TODO: hack again... adding chain, residue info...
      
      if self.ccpInfo[i]:
        member.chain = self.ccpInfo[i][0]
      else:
        member.chain = None
      
      if nmrStarNames:
        member.nmrStarName = nmrStarNames[i].value
        
        try:
          member.nmrStarLabel = eval(nmrStarLabels[i].value)
        except:
          member.nmrStarLabel = nmrStarLabels[i].value

      else:
        member.nmrStarName = None
        member.nmrStarLabel = None
      
  def setRawDihedralConstraint(self):

    #
    # Have to get some sensible values out (even if only list with
    # upper distances read in!)
    #
    
    lowerLimit = self.item.lowerLimit
    upperLimit = self.item.upperLimit
     
    #
    # TODO: these have to be transferred correctly!!!
    #
    
    if self.constraint.weight:
      energyCst = self.constraint.weight
    else:
      energyCst = None
      
    funcExp = None
        
    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraint.serial))

    self.rawConstraint = self.constraintFile.constraints[-1]

    self.rawConstraint.upperAngle = upperLimit
    self.rawConstraint.lowerAngle = lowerLimit
    self.rawConstraint.energyCst = energyCst
    self.rawConstraint.funcExp = funcExp
    
  def setRawDihedralConstraintItemInfo(self):
     
    #
    # Keep track of original names, use angle name if necessary
    # Plenty of hacking going on in here...
    #
    
    hasAngle = getOriginalData(self.format,self.constraint,self.rawConstraint,'name')

    self.rawConstraint.nmrStarNames = []
    self.rawConstraint.nmrStarLabels = []

    if hasAngle and self.useOriginalData:
    
      refMember = self.rawConstraint.nodes[0].members[0]
    
      getOriginalData(self.format,self.constraint,refMember,'seqCode',setFormatName = 'origSeqCode')
      getOriginalData(self.format,self.constraint,refMember,'resLabel',setFormatName = 'origResLabel')

    else:
    
      #
      # This can't go to the resonance level because those are reorganized
      # when running linkResonances...
      #
        
      for applData in self.constraint.findAllApplicationData(application = self.format,keyword = 'origAssign'):

        self.rawConstraint.nmrStarNames.append(applData.value)
        
      for applData in self.constraint.findAllApplicationData(application = self.format,keyword = 'origResLabel'):

        try:
          value = eval(applData.value)
        except:
          value = applData.value

        self.rawConstraint.nmrStarLabels.append(value)

  def setRawDihedralConstraintInfo(self):    

    #
    # Reset to original data...
    #
    
    if self.useOriginalData:
    
      getOriginalData(self.format,self.constraint,self.rawConstraint,'energyCst')
      getOriginalData(self.format,self.constraint,self.rawConstraint,'funcExp')
      getOriginalData(self.format,self.constraint,self.rawConstraint,'errors')
      getOriginalData(self.format,self.item,self.rawConstraint,'errors')
  
  def setRawDihedralConstraintItem(self):
    
    self.rawConstraint.nodes.append(self.rawConstraintItemClass())
    self.rawConstraintItem = self.rawConstraint.nodes[-1]
    
  def setRawDihedralConstraintItemMembers(self):
  
    #
    # TODO: add original names?!?! Not there yet?!?
    #
    
    for i in range(0,self.constraintResNum):
      
      (chainCode,seqCode,void,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])
      
      self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass(chainCode,seqCode,atomName))

      if self.ccpInfo[i]:
        self.rawConstraintItem.members[-1].chain = self.ccpInfo[i][0]
      else:
        self.rawConstraintItem.members[-1].chain = None
      
  def setRawRdcConstraint(self):
        
    self.constraintFile.constraints.append(self.rawConstraintClass(self.constraint.serial))

    self.rawConstraint = self.constraintFile.constraints[-1]

    self.rawConstraint.value = self.constraint.targetValue
    self.rawConstraint.error = self.constraint.error
    self.rawConstraint.lowerValue = self.constraint.lowerLimit
    self.rawConstraint.upperValue = self.constraint.upperLimit

    if self.negateValues:
      for attrName in ('value','lowerValue','upperValue'):
        setattr(self.rawConstraint,attrName,-getattr(self.rawConstraint,attrName))

  def setRawRdcConstraintItem(self):
    
    self.rawConstraint.nodes.append(self.rawConstraintItemClass())
    self.rawConstraintItem = self.rawConstraint.nodes[-1]

  def setRawRdcConstraintInfo(self):
    
    #
    # Reset to original data...
    #
    
    if self.useOriginalData:
    
      getOriginalData(self.format,self.constraint,self.rawConstraint,'errors')
    
  def setRawRdcConstraintItemMembers(self):
       
    #
    # Keep track of original names
    #

    nmrStarNames = self.item.findAllApplicationData(application = self.format,keyword = 'origAssign')
    nmrStarLabels =  self.item.findAllApplicationData(application = self.format,keyword = 'origResLabel')
           
    if nmrStarNames and len(nmrStarNames) > 2:
      
      #
      # Note that this is different from distances: it is possible that only one
      # original member was listed (e.g. DYANA format)
      #
    
      print "  Error: rdc constraint %d has invalid number of names %s" % (self.rawConstraint.Id,str(nmrStarNames))

    for i in range(0,self.constraintResNum):
            
      (chainCode,seqCode,void,seqInsertCode,atomName) = getNameInfo(self.resSetNames[i])
      
      self.rawConstraintItem.members.append(self.rawConstraintItemMemberClass(chainCode,seqCode,atomName))

      if self.ccpInfo[i]:
        self.rawConstraintItem.members[-1].chain = self.ccpInfo[i][0]
      else:
        self.rawConstraintItem.members[-1].chain = None
     
      member = self.rawConstraintItem.members[-1]
     
      if nmrStarNames and i < len(nmrStarNames):

        member.nmrStarName = nmrStarNames[i].value
        
        try:
          member.nmrStarLabel = eval(nmrStarLabels[i].value)
        except:
          member.nmrStarLabel = nmrStarLabels[i].value

      else:
      
        member.nmrStarName = None
        member.nmrStarLabel = None

  def checkRawConstraintItems(self,delete = True):
    
    if not self.rawConstraint.nodes:
    
      if delete:
        del self.constraintFile.constraints[-1]
        
      return False
      
    return True

  def checkRawRdcConstraintItems(self, delete = True):
  
    if len(self.rawConstraint.nodes) > 1:
      
      if delete:
        del self.constraintFile.constraints[-1]
        print "  Error: no multiple items allowed for rdc constraint %d - deleted." % self.constraint.serial

      return False
      
    return True
  
  def setRawCoordinate(self):
  
    if not self.coordinateFile.modelCoordinates.has_key(self.modelId):

      self.coordinateFile.modelCoordinates[self.modelId] = []
      
    self.coordinateFile.modelCoordinates[self.modelId].append(self.coordinatesIO.NmrStarCoordinate())
    
    modelCoordinate = self.coordinateFile.modelCoordinates[self.modelId][-1]

    modelCoordinate.Id = self.coordinateSerial
    modelCoordinate.modelId = self.modelId
    modelCoordinate.atomName = self.atomName
    modelCoordinate.chainCode = self.exportChainCode
    modelCoordinate.seqCode = self.seqCode

    modelCoordinate.atomType = self.coordAtom.elementSymbol

    modelCoordinate.x = self.x
    modelCoordinate.y = self.y
    modelCoordinate.z = self.z
    modelCoordinate.bFactor = self.bFactor
    
    #
    # Hacks for writing...
    #
    
    modelCoordinate.chain = self.chain
    
    if self.useOriginalData:
      
      modelCoordinate.origData = (self.origChainCode,self.origSeqCode,self.origSeqInsertCode,self.origResLabel,self.origAtomName)

  def getSpecificChainMappingsSeqCodeInfo(self,chainMappings):
  
    """
    This function RESETS TopShared function - always use main ID information here, not whatever is 'format' from NMR-STAR.
    The information in relation to codes is quite confused for historical reasons... 
    """

    if not self.forceChainMappings:
      firstSeqCodeKey = (1,' ')
      lastSeqCodeKey = (99999,' ')
      
      seqCodeToIdMapping = {firstSeqCodeKey: 1}
      seqCodes = [(firstSeqCodeKey,lastSeqCodeKey)]
      seqCodeToChainMapping = {firstSeqCodeKey: None}
      
      seqCodes.sort()
      seqCodes.reverse()
      
      specificInfo = (seqCodeToIdMapping,seqCodeToChainMapping,seqCodes)

    else:
      
      specificInfo = None
    
    return specificInfo

  def setCoordChainOrigData(self,cChain,coordinate):
    
    if coordinate.authorChainId != None:
      cChain.addApplicationData(Implementation.AppDataString(application = self.format,keyword = 'originalChainCode',value = coordinate.authorChainId))

    if coordinate.pdbChainId != None:
      cChain.addApplicationData(Implementation.AppDataString(application = 'pdb',keyword = 'originalChainCode',value = coordinate.pdbChainId))
      
    if hasattr(coordinate,'pdbxChainId') and coordinate.pdbxChainId != None:
      cChain.addApplicationData(Implementation.AppDataString(application = 'pdbx',keyword = 'originalChainCode',value = coordinate.pdbxChainId))

  def setCoordResidueOrigData(self,cResidue,coordinate,seqCode,seqInsertCode):
    
    if coordinate.authorSeqCode != None:
      (seqCode,seqInsertCode) = getSeqAndInsertCode(coordinate.authorSeqCode)
      cResidue.addApplicationData(Implementation.AppDataInt(application = self.format,keyword = 'originalSeqCode',value = seqCode))
      cResidue.addApplicationData(Implementation.AppDataString(application = self.format,keyword = 'originalSeqInsertCode',value = seqInsertCode))

    if coordinate.authorResName != None:
      cResidue.addApplicationData(Implementation.AppDataString(application = self.format,keyword = 'originalResLabel',value = coordinate.authorResName))

    # Also set PDB specific info, if available...
    if coordinate.pdbSeqCode != None:

      (seqCode,void) = getSeqAndInsertCode(coordinate.pdbSeqCode)
      cResidue.addApplicationData(Implementation.AppDataInt(application = 'pdb',keyword = 'originalSeqCode',value = seqCode))
      
      if coordinate.pdbSeqInsertCode != None:
        cResidue.addApplicationData(Implementation.AppDataString(application = 'pdb',keyword = 'originalSeqInsertCode',value = coordinate.pdbSeqInsertCode))

    if coordinate.pdbResName != None:
      cResidue.addApplicationData(Implementation.AppDataString(application = 'pdb',keyword = 'originalResLabel',value = coordinate.pdbResName))

  def setCoordAtomOrigData(self,cAtom,coordinate,atomName,atomStatus):
  
    authorAtomName = coordinate.authorAtomName

    if authorAtomName != None and coordinate.atomId != authorAtomName:

      cAtom.addApplicationData(Implementation.AppDataString(application = self.format,keyword = 'originalName',value = authorAtomName))

    if coordinate.pdbAtomName:

      cAtom.addApplicationData(Implementation.AppDataString(application = 'pdb',keyword = 'originalName',value = coordinate.pdbAtomName))

  #
  #
  # NMRSTAR FORMAT SPECIFIC CODE HERE
  #
  #
      
  def createEntry(self,entry,title,details,version):
    
    #
    # TODO: set name, title correctly (should be BMRB ID if available!)
    #
    #
    # TODO: read in first (entry) frame and check for info from there if
    #       real BMRB file...
    #
    
    if not entry:
      
      setCurrentStore(self.project,'NmrEntryStore')

      fileName = '.'.join(os.path.basename(self.file.name).split('.')[:-1])

      entryName = returnMemopsLine(fileName)
      entryTitle = returnMemopsLine(title)
      entryTitle = entryTitle.strip()

      self.entry = self.project.currentNmrEntryStore.findFirstEntry(name=entryName)

      if not self.entry:
        self.entry = NmrEntry.Entry(self.project.currentNmrEntryStore,name = entryName, title = entryTitle,details = details)
    else:
      self.entry = entry

    if version[0] == '3':

      #
      # Specific for constraint entry: make sure everything is tracked!
      # TODO: should NOT be in here ?@!?!
      #

      #
      # Block definition in constraint file
      #
      
      if version == '3.0':
        saveFrameName = 'file_characteristics'
      elif version == '3.1':
        saveFrameName = 'constraint_statistics'
        
      self.setAllSfApplData(self.entry,saveFrameName,'ID')

      #
      # File level comment blocks
      #

      if version == '3.0':
        saveFrameName = 'MR_file_comment'
      elif version == '3.1':
        saveFrameName = 'org_constr_file_comment'

      self.setAllSfApplData(self.entry,saveFrameName,'ID')

  #
  # Below all code that could go to DataFormat.py, but currently NMR-STAR project reading (3.1 and up) specific!
  #
  
  def setEntryInformation(self):

    myEntry = self.file.entryInformation.entry
    
    #
    # Only reset if title and details if they have the default values (so new entry created for reading this project!)
    #

    if self.entry.title == 'None':
      self.entry.setTitle(myEntry.title)
    if self.entry.details == 'None':
      self.entry.setDetails(myEntry.details)
    if self.entry.entryType == 'None':
      self.entry.setEntryType(myEntry.exptMethodSub)
    if self.entry.bmrbProcessing is None:
      self.entry.setBmrbProcessing(myEntry.specProcInstr)

    tagNames = ('exptMethod', 'versionType','subDate','accDate','origin','nmrStarVers','origNmrStarVers',
                'depRelCoord','depRelConstr','depRelNmr','depRelSeq','caspTarget','updBmrbCode','replBmrbCode','updPdbCode',
                'replPdbCode','bmrbUpdDet','pdbUpdDet','relReq','pdbDepSite','pdbProcSite','bmrbDepSite','bmrbProcSite',
                'rcsbAnn','assBmrbId','assBmrbCode','assPdbId','assPdbCode','dateNmrConstr','recAuthAppr','pdbDateRec',
                'authRelStatCode','authApprType')

    for tagName in tagNames:
      value = getattr(myEntry,tagName)
      if value is None:
        continue

      keywds = {'application': self.format,
                'keyword':     tagName,
                'value':       str(value)}

      appData = Implementation.AppDataString(**keywds)
      self.entry.addApplicationData(appData)

    for entrySrc in self.file.entryInformation.entrySrcs:
      newOrganisation = self.project.currentAffiliationStore.newOrganisation(name = entrySrc.orgFullName)

      newGroup = newOrganisation.newGroup(name = entrySrc.orgFullName)

      self.entry.addLaboratory(newGroup)

      #for contactPerson in self.entry.sortedContactPersons():
      #  for ping in contactPerson.sortedPersonInGroups():
      #    addresses = ping.mailingAddress.split(', ')

      #    if entrySrc.orgFullName in addresses:
      #      ping.setGroup(newGroup)
      #      break

    for relEntry in self.file.entryInformation.relatedEntries:

      if relEntry.dbName not in ('PDB', 'BMRB'):
        continue

      newRelEntry = self.entry.newRelatedEntry(dbCode = relEntry.dbAcc,
                                               dbName = relEntry.dbName,
                                               relationship = relEntry.reln)

    if self.file.entryInformation.structKeywds:
      for structKeywd in self.file.entryInformation.structKeywds:
        if structKeywd.keywd not in self.entry.keywords:
          self.entry.addKeyword(structKeywd.keywd)

    #for release in self.file.entryInformation.releaseInfo:
      
    #  tagNames2 = ('entryId','relNo','date','subDateRel','relType','author','detail')

    #  for tagName in tagNames2:
    #    value = getattr(release,tagName)
    #    if value is None:
    #      continue

    #    keywds2 = {'application': self.format,
    #               'keyword':     tagName,
    #               'value':       str(value)}

    #    appData = Implementation.AppDataString(**keywds2)
    #    self.entry.addApplicationData(appData)

    #i = 0

    #for strucGenom in self.file.entryInformation.strucGenomInfo:

    #  i += 1

    #  keywds4 = {'application': self.format,
    #               'keyword':     tagName,
    #               'value':       'SG' + str(i)}

    #  tagNames3 = ('projName','centerFullName')

    #  for tagName in tagNames3:
    #    value = getattr(strucGenom,tagName)
    #    if value is None:
    #      continue

    #    keywds5 = {'application': self.format,
    #               'keyword':     tagName,
    #               'value':       str(value + str(i) )}

    #    appData = Implementation.AppDataString(**keywds5)
    #    self.entry.addApplicationData(appData)

  def readStudies(self):
    
    #
    # TODO: Might want to copy over some information if entry was given and already has study, if necessary?
    #
    
    if not self.entry.study:

      saveFrameName = self.studyFile.saveFrame.name
  
      foreignIdTagName = (saveFrameName,'ID')
  
      if not self.starIdToCcpnObject.has_key(foreignIdTagName):
        self.starIdToCcpnObject[foreignIdTagName] = []
  
        self.ccpnClassToStarId['Study'] = foreignIdTagName
  
      study = self.studyFile.studies[0]
  
      newStudy = self.project.currentNmrEntryStore.newStudy(name      = study.name,
                                                            studyType = study.type)
  
      self.starIdToCcpnObject[foreignIdTagName].append((study.Id,newStudy))
  
      self.entry.setStudy(newStudy)
  
      for studyKeywd in self.studyFile.studyKeywds:
        newStudy.addKeyword(studyKeywd.keywd)


  def readInstruments(self):

    setCurrentStore(self.project,'InstrumentStore')
    setCurrentStore(self.project,'AffiliationStore')

    tableName = self.instrumentFile.tableName

    foreignIdTagName = (tableName,'ID')

    if self.instrumentFile.specListDetails:
      self.entry.spectrometerListDetails = self.instrumentFile.specListDetails

    if not self.starIdToCcpnObject.has_key(foreignIdTagName):
      self.starIdToCcpnObject[foreignIdTagName] = []

      self.ccpnClassToStarId['NmrSpectrometer'] = foreignIdTagName

    for instrument in self.instrumentFile.instruments:

      #print 'INSTR: [%s]' % (instrument)

      manufacturer = self.project.currentAffiliationStore.findFirstOrganisation(name = instrument.manufacturer.capitalize())

      if not manufacturer:
        manufacturer = self.project.currentAffiliationStore.newOrganisation(name = instrument.manufacturer.capitalize())

      newNmrSpec = self.project.currentInstrumentStore.newNmrSpectrometer(name         = instrument.name,
                                                                          model        = instrument.model,
                                                                          nominalFreq  = instrument.field_strength,
                                                                          serialNumber = instrument.serialNumber,
                                                                          details      = instrument.details,
                                                                          manufacturer = manufacturer)

      newNmrSpec.protonFreq = returnFloat(instrument.field_strength)

      keywds = {'application': self.format,
                'keyword':     'specFlag',
                'value':       True}

      if not self.entry.findFirstApplicationData(**keywds):
        appData = Implementation.AppDataBoolean(**keywds)
        self.entry.addApplicationData(appData)

      keywds = {'application': self.format,
                'keyword':     'specFlag',
                'value':       False}

      if not newNmrSpec.findFirstApplicationData(**keywds):
        appData = Implementation.AppDataBoolean(**keywds)
        newNmrSpec.addApplicationData(appData)

      self.starIdToCcpnObject[foreignIdTagName].append((instrument.Id,newNmrSpec))


  def readProbes(self):

    setCurrentStore(self.project,'InstrumentStore')
    setCurrentStore(self.project,'AffiliationStore')

    tableName = self.probeFile.tableName

    foreignIdTagName = (tableName,'ID')

    if not self.starIdToCcpnObject.has_key(foreignIdTagName):
      self.starIdToCcpnObject[foreignIdTagName] = []

      self.ccpnClassToStarId['NmrProbe'] = foreignIdTagName

    for probe in self.probeFile.probes:

      manufacturer = self.project.currentAffiliationStore.findFirstOrganisation(name = probe.manufacturer.capitalize())

      if not manufacturer:
        manufacturer = self.project.currentAffiliationStore.newOrganisation(name = probe.manufacturer.capitalize())

      name = manufacturer.name + '-' + probe.model

      newNmrProbe = self.project.currentInstrumentStore.newNmrProbe(name         = name,
                                                                    model        = probe.model,
                                                                    serialNumber = probe.serialNumber,
                                                                    details      = probe.details,
                                                                    probeType    = probe.pType,
                                                                    diameter     = probe.diameter,
                                                                    manufacturer = manufacturer)

      self.starIdToCcpnObject[foreignIdTagName].append((probe.Id,newNmrProbe))


  def readSoftware(self):

    setCurrentStore(self.project,'MethodStore')

    saveFrameName = self.softwareFile.saveFrame.name

    foreignIdTagName = (saveFrameName,'ID')

    if not self.starIdToCcpnObject.has_key(foreignIdTagName):
      self.starIdToCcpnObject[foreignIdTagName] = []

      self.ccpnClassToStarId['Software'] = foreignIdTagName

    for software in self.softwareFile.software:

      keywds = {'name':    software.name,
                'details': software.details}

      if software.version:
        keywds['version'] = software.version

      else:
        keywds['version'] = 'unknown' # TODO - what do we set version to as is required by CCPN

      if software.venName:
        keywds['vendorName'] = software.venName

      if software.venAddress:
        keywds['vendorAddress'] = software.venAddress

      if software.venEAddress:
        keywds['vendorWebAddress'] = software.venEAddress

      if software.tasks:

        if software.tasks.count(','):
          keywds['tasks'] = software.tasks.split(',')
        else:
          keywds['tasks'] = [software.tasks]

      newSoft = self.project.currentMethodStore.findFirstSoftware(**keywds)

      if newSoft:
        print "  Warning: software %s already exists! Cannot create again." % software.name
      else:
        newSoft = self.project.currentMethodStore.newSoftware(**keywds)

        self.starIdToCcpnObject[foreignIdTagName].append((self.softwareFile.softwareId,newSoft))

        self.entry.addSoftware(newSoft)


  def readMethods(self):

    setCurrentStore(self.project,'MethodStore')

    saveFrameName = self.methodFile.saveFrame.name

    foreignIdTagName = (saveFrameName,'ID')

    if not self.starIdToCcpnObject.has_key(foreignIdTagName):
      self.starIdToCcpnObject[foreignIdTagName] = []

      self.ccpnClassToStarId['Method'] = foreignIdTagName

    for method in self.methodFile.method:

      keywds = {'name': method.details}

      newMethod = self.project.currentMethodStore.findFirstMethod(**keywds)

      if newMethod:
        print "  Warning: method %s already exists! Cannot create again." % method.details
      else:
        newMethod = self.project.currentMethodStore.newMethod(**keywds)

        self.starIdToCcpnObject[foreignIdTagName].append((self.methodFile.methodId,newMethod))


  def readConfStats(self):

    # TODO - get the correct StructureEnsemble here.

    strucEns = self.project.currentStructureEnsemble or self.project.findFirstStructureEnsemble()

    if not strucEns:
      return

    saveFrameName = self.confStatFile.saveFrame.name

    for confStat in self.confStatFile.confStat:

      # TODO: find the right StructureEnsemble using this ID.
      #strucEnsId = confStat.ensId

      for (appDataName, tagName) in ( ('calculated','calculated'),
                                      ('criteria','criteria'),
                                      ('representative','representative'),
                                      ('repr_criteria','repr_criteria'),
                                      ('sgDetails','details') ):

        if not hasattr(confStat, tagName):
          continue

        tagValue = getattr(confStat, tagName)

        keywds = {'application': self.format,
                  'keyword':     appDataName}

        appData = strucEns.findFirstApplicationData(**keywds)

        if appData:
          for appData in strucEns.findAllApplicationData(**keywds):
            strucEns.removeApplicationData(appData)

        keywds['value'] = str(tagValue)

        appData = Implementation.AppDataString(**keywds)
        strucEns.addApplicationData(appData)

  def readNmrExperiments(self):
    # TODO should this link into self.readAcq/ProcPars?
    setCurrentStore(self.project,'NmrProject')
    setCurrentStore(self.project,'AffiliationStore')

    tableName = self.nmrExperimentFile.tableName

    tagNames = ('experimentListDetails',)

    for tagName in tagNames:
      value = getattr(self.nmrExperimentFile,tagName)
      if value is None:
        continue

      keywds = {'application': self.format,
                'keyword':     tagName,
                'value':       str(value)}

      appData = Implementation.AppDataString(**keywds)
      self.entry.addApplicationData(appData)

    foreignIdTagName = (tableName,'ID')

    if not self.starIdToCcpnObject.has_key(foreignIdTagName):
      self.starIdToCcpnObject[foreignIdTagName] = []

      #print 'ID: [%s]' % str(foreignIdTagName)

      self.ccpnClassToStarId['Experiment'] = foreignIdTagName

    for nmrExperiment in self.nmrExperimentFile.experiments:
      
      # Possible in old files
      if not nmrExperiment.name:
        continue

      newExpt = self.project.currentNmrProject.newExperiment(name   = nmrExperiment.name,
                                                             numDim = 2)

      self.starIdToCcpnObject[foreignIdTagName].append((nmrExperiment.Id,newExpt))

      self.entry.addExperiment(newExpt)

      if nmrExperiment.sampleState:
        newExpt.sampleState = nmrExperiment.sampleState

      #tagNames = ('sampleState',)

      #for tagName in tagNames:
      #  value = getattr(nmrExperiment,tagName)
      #  if value is None:
      #    continue

      #  keywds = {'application': self.format,
      #            'keyword':     tagName,
      #            'value':       str(value)}

      #  appData = Implementation.AppDataString(**keywds)
      #  newExpt.addApplicationData(appData)

      if 'NmrSpectrometer' in self.ccpnClassToStarId:
        for (specId, instrument) in self.starIdToCcpnObject[self.ccpnClassToStarId['NmrSpectrometer'] ]:
          if nmrExperiment.specId == specId:
            newExpt.setSpectrometer(instrument)

            keywds = {'application': self.format,
                      'keyword':     'specFlag',
                      'value':       False}

            appData = instrument.findFirstApplicationData(**keywds)

            if instrument.findFirstApplicationData(**keywds):
              instrument.removeApplicationData(appData)

              keywds['value'] = True

              #print 'Adding App data'

              appData = Implementation.AppDataBoolean(**keywds)
              instrument.addApplicationData(appData)

            break

      if 'NmrProbe' in self.ccpnClassToStarId:
        for (probeId, probe) in self.starIdToCcpnObject[self.ccpnClassToStarId['NmrProbe'] ]:
          if nmrExperiment.probeId == probeId:
            newExpt.setProbe(probe)
            break

      if 'Sample' in self.ccpnClassToStarId:
        for (sampleId,sample) in self.starIdToCcpnObject[self.ccpnClassToStarId['Sample'] ]:
          if nmrExperiment.sampleId == sampleId:
            newExpt.setSample(sample)
            break
 
      if 'SampleConditionSet' in self.ccpnClassToStarId:
        for (sampleCondListId,samCondSet) in self.starIdToCcpnObject[self.ccpnClassToStarId['SampleConditionSet'] ]:
        
          if nmrExperiment.sampleCondListId == sampleCondListId:
            newExpt.setSampleConditionSet(samCondSet)
            break

  def setNmrSpectrometerAppData(self):

    specFlag = True

    setCurrentStore(self.project,'InstrumentStore')

    for instrument in self.entry.root.currentInstrumentStore.findAllInstruments(className = 'NmrSpectrometer'):

      keywds = {'application': self.format,
                'keyword':     'specFlag',
                'value':       False}

      appData = instrument.findFirstApplicationData(**keywds)

      #print 'IN APP DATA: [%s] [%s]' % (instrument, appData)

      if appData:
        specFlag = False
        break

    if not specFlag:
      keywds = {'application': self.format,
                'keyword':     'specFlag',
                'value':       True}

      appData = self.entry.findFirstApplicationData(**keywds)

      if self.entry.findFirstApplicationData(**keywds):
        self.entry.removeApplicationData(appData)

        keywds['value'] = specFlag # Must be False

        appData = Implementation.AppDataBoolean(**keywds)
        self.entry.addApplicationData(appData)

  def parseBmrbLabelName(self, bmrbLabelName):

    specificResidue = None

    resLabPatt = re.compile("\]\-([A-Za-z\d]+)$")

    searchObj = resLabPatt.search(bmrbLabelName)

    if searchObj:
      specificResidue = searchObj.group(1)

    bmrbLabelName = resLabPatt.sub(']', bmrbLabelName)

    #print 'SPEC: [%s] [%s]' % (bmrbLabelName, specificResidue)

    if bmrbLabelName.count(';'):
      labels = bmrbLabelName.split(';')

    else:
      labels = [bmrbLabelName]

    labelInfo = []

    for label in labels:
      label = label.strip('[] \t')

      #print 'LABEL: [%s]' % label

      if label == 'natural abundance':
        return None # pass ???

      else:
        labelPatt = re.compile("^(U\-)?(\d+%)? ?(\d+)([HCN])([A-Z]*)$")

        searchObj = labelPatt.search(label)

        if searchObj:
          uniform  = searchObj.group(1)
          percent  = searchObj.group(2)
          mass     = searchObj.group(3)
          atomType = searchObj.group(4)
          specificAtom = searchObj.group(5)

          labelInfo.append( (uniform, percent, mass, atomType, specificResidue, specificAtom) )

    return labelInfo

  def makeUniformLabels(self, resLabel, atomType, mass, percent):

    isotopeMasses = {'C': '13',
                     'N': '15',
                     'H': '2'}

    isotopeNuclei = {'C': 'carbon',
                     'N': 'nitrogen',
                     'H': 'hydrogen'}

    if mass != isotopeMasses[atomType]:
      print '  Warning: labelled %s isotope does not have mass of %s' % (isotopeNuclei[atomType], isotopeMasses[atomType])
      return

    if not percent:
      weight = 0.00
    else:
      weight = float(percent[:-1])/100.0

    if weight < 0 and weight > 1:
      print '  Warning: degree of isotope labelling not in the correct range (0-100%)'
      return

    uniformAtomLabel2 = None

    uniformAtomLabel1 = resLabel.newUniformAtomLabel(elementName=atomType, massNumber=int(mass), weight=weight)
    if weight != 1.00:
      uniformAtomLabel2 = resLabel.newUniformAtomLabel(elementName=atomType, massNumber=int(mass)-1, weight=1.00-weight)

    return (uniformAtomLabel1, uniformAtomLabel2)

  def makeSpecAtomLabels(self, resLabel, atomType, mass, percent, specificAtom):

    isotopeMasses = {'C': '13',
                     'N': '15',
                     'H': '2'}

    isotopeNuclei = {'C': 'carbon',
                     'N': 'nitrogen',
                     'H': 'hydrogen'}

    if mass != isotopeMasses[atomType]:
      print '  Warning: labelled %s isotope does not have mass of %s' % (isotopeNuclei[atomType], isotopeMasses[atomType])
      return

    if not percent:
      weight = 0.00
    else:
      weight = float(percent[:-1])/100.0

    if weight < 0 and weight > 1:
      print '  Warning: degree of isotope labelling not in the correct range (0-100%)'
      return

    specAtomLabel2 = None

    specAtomLabel1 = resLabel.newSingleAtomLabel(atomName=atomType + specificAtom, massNumber=int(mass), weight=weight)
    if weight != 1.00:
      specAtomLabel2 = resLabel.newSingleAtomLabel(atomName=atomType + specificAtom, massNumber=int(mass)-1, weight=1.00-weight)

    return (specAtomLabel1, specAtomLabel2)

  def doLabels(self, labelInfo, bmrbLabelName, molLabel, molecule=None, molRes=None):

    schemeNameDict = {'natural abundance':                'NatAbun',
                      '[U-95% 15N]':                      'uni_15N',
                      '[U-95% 13C; U-95% 15N]':           'uni_15N13C',
                      '[U-95% 13C; U-95% 15N; U-95% 2H]': 'uni_15N13C2H'}

    if molRes is None:
      molResId = 1
    else:
      molResId = molRes.serial

    resLabel = molLabel.newResLabel(resId=molResId)

    if bmrbLabelName in schemeNameDict:
      resLabelFrac = resLabel.newResLabelFraction(schemeName=schemeNameDict[bmrbLabelName])

      #print 'RES: [%s] [%s] [%s] [%s]' % (resLabel, resLabelFrac, resLabelFrac.schemeName, resLabelFrac.getLabelingScheme() )

    else:
      uniformLabels = None

      if labelInfo:
        for label in labelInfo:

          (uniform, percent, mass, atomType, specificResidue, specificAtom) = label

          ccpCode = ''

          if molecule and molecule.molType in ('DNA', 'RNA', 'DNA/RNA') and molRes:

            nucTlcDict = {'A': 'Ade',
                          'C': 'Cyt',
                          'G': 'Gua',
                          'T': 'Thy',
                          'U': 'Ura'}

            ccpCode = nucTlcDict[molRes.ccpCode]

          #print "DATA: [%s] [%s] [%s] [%s] [%s] [%s] [%s] [%s]" % (uniform, percent, mass, atomType, specificResidue, molRes, ccpCode, specificAtom)

          if specificResidue is not None and specificResidue != ccpCode:
            continue

          if uniform == 'U-':
            uniformLabels = self.makeUniformLabels(resLabel, atomType, mass, percent)

            if uniformLabels is None:
              continue

          elif specificAtom:
            specAtomLabels = self.makeSpecAtomLabels(resLabel, atomType, mass, percent, specificAtom)

          else:
            continue

      #print 'RES2: [%s] [%s] [%s]' % (resLabel, uniformLabels[0], uniformLabels[1])


  def readSamples(self):
    
    setCurrentStore(self.project,'SampleStore')
    setCurrentStore(self.project,'Classification')
    setCurrentStore(self.project,'RefSampleComponentStore')

    if not hasattr(self, 'sampleCategory'):
      self.sampleCategory = self.project.currentClassification.newSampleCategory(name = self.format)

    saveFrameName = self.sampleFile.saveFrame.name
    foreignIdTagName = (saveFrameName,'ID')

    if not self.starIdToCcpnObject.has_key(foreignIdTagName):
      self.starIdToCcpnObject[foreignIdTagName] = []

      self.ccpnClassToStarId['Sample'] = foreignIdTagName

    for sample in self.sampleFile.samples:

      #print 'SAMPLE: [%s]' % sample.name

      newSam = self.project.currentSampleStore.newSample(name = sample.name,
                                                         sampleCategories = [self.sampleCategory],
                                                         details = sample.details)

      tagNames = ('sampleState','solventSys','aggrSamNo')

      for tagName in tagNames:
        value = getattr(sample,tagName)
        if value is None:
          continue

        keywds = {'application': self.format,
                  'keyword':     tagName,
                  'value':       str(value)}

        appData = Implementation.AppDataString(**keywds)
        newSam.addApplicationData(appData)

      self.starIdToCcpnObject[foreignIdTagName].append((self.sampleFile.sampleId,newSam))

      tableName = self.sampleFile.tableName
      foreignIdTagName = (tableName,'ID')

      if not self.starIdToCcpnObject.has_key(foreignIdTagName):
        self.starIdToCcpnObject[foreignIdTagName] = []

        self.ccpnClassToStarId['AbstractComponent'] = foreignIdTagName

      for sampleComponent in self.sampleFile.sampleComponents:

        #print 'SAM COMP: [%s]' % sampleComponent.name

        molecule = None

        moleculeForeignKey = ('Entity','ID') # Could be done in a more sensible way?
        if self.starIdToCcpnObject.has_key(moleculeForeignKey):
          for (starId,ccpnObject) in self.starIdToCcpnObject[moleculeForeignKey]:
            if starId == sampleComponent.molId:
              molecule = ccpnObject
              break

        #if not molecule:
        #  continue

        self.count += 1

        #print 'COUNT: [%s]' % self.count

        newRefComp = self.project.currentRefSampleComponentStore.findFirstComponent(name = sampleComponent.name) # + ':' + sample.name

        if not newRefComp:
          newRefComp = self.project.currentRefSampleComponentStore.newMolComponent(name = sampleComponent.name, details =  sampleComponent.name) # Need to set right molecule, took out this from name: + ':' + sample.name

        if molecule:

          if not molecule.isFinalised:
            molecule.isFinalised = True

          newRefComp.setMolecule(molecule)

        bmrbLabelName = sampleComponent.label

        if bmrbLabelName:
          newLabMol = self.project.findFirstLabeledMolecule(name=newRefComp.name)

          if not newLabMol:
            newLabMol = self.project.newLabeledMolecule(name=newRefComp.name)

          #print 'REF NAME: [%s]' % newRefComp.name

          molLabel = newLabMol.newMolLabel()

          labMix = newLabMol.newLabeledMixture()
          newRefComp.setLabeledMixture(labMix)

          molLabelFrac = labMix.newMolLabelFraction(molLabel=molLabel)

          labelInfo = self.parseBmrbLabelName(bmrbLabelName)

          if molecule:
            for molRes in molecule.sortedMolResidues():
              self.doLabels(labelInfo, bmrbLabelName, molLabel, molecule, molRes)

          else:
            self.doLabels(labelInfo, bmrbLabelName, molLabel)

        concentration = 0.0

        if sampleComponent.conc:
          numberPatt = re.compile("^[\d\.]+$")

          searchObj = numberPatt.search(sampleComponent.conc)

          if searchObj:
            concentration = returnFloat(sampleComponent.conc)

        concentrationUnit = sampleComponent.conc_unit

        displayUnit = concentrationUnit

        if concentrationUnit == 'mM':
          concentration = concentration / 1000.0
          concentrationUnit = 'M'
          displayUnit = 'mM'

        if concentrationUnit == 'uM':
          concentration = concentration / 1000000.0
          concentrationUnit = 'M'
          displayUnit = 'uM'

        elif concentrationUnit == '%':
          concentration = concentration / 100.0
          concentrationUnit = 'm3/m3' # This may not be correct, but generally %'s refer to volumes
          displayUnit = '%'
          
        elif concentrationUnit == 'mg/ml':
          concentration = concentration
          concentrationUnit = 'kg/m3' # This may not be correct, but generally %'s refer to volumes
          displayUnit = 'mg/mls'

        elif concentrationUnit == 'x':
          print "  Warning: cannot handle unit 'x' for sample concentration - ignoring info."
          continue

        sampComp = newSam.newSampleComponent(refComponent       = newRefComp,
                                             concentration      = concentration,
                                             concentrationError = sampleComponent.conc_err,
                                             concentrationUnit  = concentrationUnit,
                                             concDisplayUnit    = displayUnit)

        if concentration == 0.0:
          sampComp.addApplicationData(Implementation.AppDataString(
            application = self.format,
            keyword = 'concentration',
            value = str(sampleComponent.conc)
            ) )

        self.starIdToCcpnObject[foreignIdTagName].append((sampleComponent.Id,newRefComp))


  def readSampleConditions(self):
    
    setCurrentStore(self.project,'NmrProject')

    newSampleConditionSet = self.project.currentNmrProject.newSampleConditionSet(name = self.sampleConditionFile.sampleConditionSetName, details = self.sampleConditionFile.sampleConditionSetDetails)

    # Map between star ID and CCPN object! Only necessary for some objects, not all...
    # TODO should probably use nmrStarDict to set foreign keys automatically...
    saveFrameName = self.sampleConditionFile.saveFrame.name
    foreignIdTagName = (saveFrameName,'ID') # Could do something better?

    if not self.starIdToCcpnObject.has_key(foreignIdTagName):
      self.starIdToCcpnObject[foreignIdTagName] = []

      self.ccpnClassToStarId['SampleConditionSet'] = foreignIdTagName

    self.starIdToCcpnObject[foreignIdTagName].append((self.sampleConditionFile.sampleConditionSetId,newSampleConditionSet))
    
    for sampleCondition in self.sampleConditionFile.sampleConditions:
      if not sampleCondition.type:
        continue

      samCondSet = newSampleConditionSet.findFirstSampleCondition(condition = sampleCondition.type,
                                                                  value     = sampleCondition.value,
                                                                  error     = sampleCondition.error,
                                                                  unit      = sampleCondition.units)

      if not samCondSet:

        samCondSet = newSampleConditionSet.newSampleCondition(condition = sampleCondition.type,
                                                              value     = sampleCondition.value,
                                                              error     = sampleCondition.error,
                                                              unit      = sampleCondition.units)

        #print 'SCD: [%s] [%s] [%s] [%s]' % (samCondSet.condition,samCondSet.value,samCondSet.error,samCondSet.unit)
        #print 'SCD2: [%s] [%s] [%s] [%s]' % (sampleCondition.type,sampleCondition.value,sampleCondition.error,sampleCondition.units)


  def readChemShiftRefs(self):

    setCurrentStore(self.project,'NmrProject')
    setCurrentStore(self.project,'ChemElementStore')

    saveFrameName = self.chemShiftRefFile.saveFrame.name
  
    foreignIdTagName = (saveFrameName,'ID')
  
    if not self.starIdToCcpnObject.has_key(foreignIdTagName):
      self.starIdToCcpnObject[foreignIdTagName] = []
  
      self.ccpnClassToStarId['ChemicalShiftReference'] = foreignIdTagName

    for chemShiftRef in self.chemShiftRefFile.chemShiftReferences:

      newChemShiftRef = None
      refType = chemShiftRef.refType
      if not refType:
        # Possibly missing in older files
        continue

      if chemShiftRef.refMeth == 'internal':

        if chemShiftRef.chemShiftVal is None:
          chemShiftRef.chemShiftVal = -999999999.0

        newChemShiftRef = self.project.currentNmrProject.newInternalShiftReference(atomGroup = chemShiftRef.atomGroup,
                                                                                   isotopeCode = '%s%s' % (chemShiftRef.isotopeNum, chemShiftRef.atomType),
                                                                                   molName = chemShiftRef.molName,
                                                                                   referenceType = refType,
                                                                                   value = chemShiftRef.chemShiftVal,
                                                                                   unit = chemShiftRef.chemShiftUnits,
                                                                                   indirectShiftRatio = chemShiftRef.indirectShiftRatio)

        if chemShiftRef.location:
          newChemShiftRef.addApplicationData(Implementation.AppDataString(
            application = self.format,
            keyword = 'location',
            value = str(chemShiftRef.location)
            ) )

        if chemShiftRef.geometry:
          newChemShiftRef.addApplicationData(Implementation.AppDataString(
            application = self.format,
            keyword = 'geometry',
            value = str(chemShiftRef.geometry)
            ) )

        if chemShiftRef.axis:
          newChemShiftRef.addApplicationData(Implementation.AppDataString(
            application = self.format,
            keyword = 'axis',
            value = str(chemShiftRef.axis)
            ) )

      elif chemShiftRef.refMeth == 'external':

        if chemShiftRef.chemShiftVal is None:
          chemShiftRef.chemShiftVal = -999999999.0

        newChemShiftRef = self.project.currentNmrProject.newExternalShiftReference(atomGroup = chemShiftRef.atomGroup,
                                                                                   isotopeCode = '%s%s' % (chemShiftRef.isotopeNum, chemShiftRef.atomType),
                                                                                   molName = chemShiftRef.molName,
                                                                                   referenceType = refType,
                                                                                   value = chemShiftRef.chemShiftVal,
                                                                                   unit = chemShiftRef.chemShiftUnits,
                                                                                   indirectShiftRatio = chemShiftRef.indirectShiftRatio,
                                                                                   location = chemShiftRef.location,
                                                                                   sampleGeometry = chemShiftRef.geometry,
                                                                                   axis = chemShiftRef.axis)

      else:

        if chemShiftRef.chemShiftVal is None:
          chemShiftRef.chemShiftVal = -999999999.0

        # Can only do Internal or External - do as External and set chemShiftRef.refMeth as application data

        newChemShiftRef = self.project.currentNmrProject.newExternalShiftReference(atomGroup = chemShiftRef.atomGroup,
                                                                                   isotopeCode = '%s%s' % (chemShiftRef.isotopeNum, chemShiftRef.atomType),
                                                                                   molName = chemShiftRef.molName,
                                                                                   referenceType = refType,
                                                                                   value = chemShiftRef.chemShiftVal,
                                                                                   unit = chemShiftRef.chemShiftUnits,
                                                                                   indirectShiftRatio = chemShiftRef.indirectShiftRatio,
                                                                                   location = chemShiftRef.location,
                                                                                   sampleGeometry = chemShiftRef.geometry,
                                                                                   axis = chemShiftRef.axis)

        newChemShiftRef.addApplicationData(Implementation.AppDataString(
          application = self.format,
          keyword = 'refType',
          value = str(chemShiftRef.refMeth)
          ) )

      self.starIdToCcpnObject[foreignIdTagName].append((self.chemShiftRefFile.chemShiftRefId,newChemShiftRef))

      keywds = {'application': self.format,
                'keyword':     'ShiftReferenceListId',
                'value':       self.chemShiftRefFile.chemShiftRefId}
      appData = Implementation.AppDataInt(**keywds)
      newChemShiftRef.addApplicationData(appData)

  def readSourceFile(self):
    
    setCurrentStore(self.project,'Taxonomy')

    tableName = self.sourceFile.tableName
    foreignIdTagName = (tableName,'ID')

    if not self.starIdToCcpnObject.has_key(foreignIdTagName):
      self.starIdToCcpnObject[foreignIdTagName] = []

      self.ccpnClassToStarId['NaturalSource'] = foreignIdTagName

    exptSource = None
    natSource = None

    for source in self.sourceFile.sources:

      if source.sourceType == 'experimental':
        exptSource = source
      elif source.sourceType == 'natural':
        natSource = source

      if not source.orgName:
        if not source.orgSci:
          source.orgSci = source.orgName = 'unknown'
        else:
          source.orgName = 'unknown'

      newSource = self.project.currentTaxonomy.newNaturalSource(organismName   = source.orgName,
                                                                scientificName = source.orgSci,
                                                                strain         = source.strain,
                                                                variant        = source.variant,
                                                                organ          = source.organ,
                                                                tissue         = source.tissue,
                                                                cellLine       = source.cell_line,
                                                                cellType       = source.cell_type,
                                                                atccNumber     = source.atccNo,
                                                                plasmid        = source.plasmid,
                                                                ncbiTaxonomyId = source.ncbiTaxId,
                                                                details        = source.details,
                                                                species        = source.species,
                                                                genus          = source.genus,
                                                                geneMnemonic   = source.geneMnemonic,
                                                                organelle      = source.organelle)

      if source.sourceType == 'natural' and newSource:
        tagNames = ('superKingdom','kingdom')

        for tagName in tagNames:
          value = getattr(natSource,tagName)
          if value is None:
            continue

          keywds = {'application': self.format,
                    'keyword':     tagName,
                    'value':       str(value)}

          appData = Implementation.AppDataString(**keywds)
          newSource.addApplicationData(appData)

      elif source.sourceType == 'experimental' and newSource:
        tagNames = ('pdbVectorName',)

        for tagName in tagNames:
          value = getattr(exptSource,tagName)
          if value is None:
            continue

          keywds = {'application': self.format,
                    'keyword':     tagName,
                    'value':       str(value)}

          appData = Implementation.AppDataString(**keywds)
          newSource.addApplicationData(appData)
      
      self.starIdToCcpnObject[foreignIdTagName].append((source.Id,newSource))

      entryMol = None

      foreignIdTagNameMol = ('Entity','ID')

      if foreignIdTagNameMol in self.starIdToCcpnObject:

        if exptSource:
          for (molId,molecule) in self.starIdToCcpnObject[foreignIdTagNameMol]:

            if exptSource.entId == molId:
              prodMethod = exptSource.prodMethod

              if prodMethod is None:
                prodMethod = 'recombinant organism'  # TODO: we are forcing this - is this right?

              entryMol = self.entry.findFirstEntryMolecule(molecule=molecule)

              if not entryMol:
                entryMol = self.entry.newEntryMolecule(productionMethod=prodMethod,
                                                       molecule=molecule,
                                                       experimentalSource=newSource)

              else:
                entryMol.setProductionMethod(prodMethod)
                entryMol.setExperimentalSource(newSource)

              vecType = exptSource.vecType

              if vecType is not None:
                entryMol.setVectorType(vecType)
              elif prodMethod in ('recombinant organism', 'recombinant technology'):
                entryMol.setVectorType('plasmid') # TODO: we are forcing this value as well

              break

        if natSource:
          for (molId,molecule) in self.starIdToCcpnObject[foreignIdTagNameMol]:
            if natSource.entId == molId:
              if not molecule.naturalSource:
                molecule.setNaturalSource(newSource)

              entryMol = self.entry.findFirstEntryMolecule(molecule=molecule)

              if not entryMol:
                entryMol = self.entry.newEntryMolecule(productionMethod='recombinant organism', molecule=molecule)

              orgType = natSource.orgType

              if orgType is not None:
                entryMol.setSourceType(orgType)

      #if self.sourceFile.source == 'natural':
      #  self.entry.setNaturalSource(newSource)
      #elif self.sourceFile.source == 'experimental':
      #  self.entry.setExperimentalSource(newSource)

  #
  # Application data handling for NMR-STAR...
  #
  
  def setSequenceSfApplData(self):
  
    #
    # TODO: SET and MATCH by name!! Use IDs for all this stuff!!
    #
    
    pass
    """
    MRblockKey = 'MR_file_block_position'
    
    if self.constraintFile.saveFrame.tags.has_key(MRblockKey):

      sfKey = self.constraintFile.saveFrame.tags[MRblockKey]
      self.saveFrame = self.constraintFile.saveFrame
      saveFrameName = self.saveFrame.tags['Sf_category']

      self.setSfApplData(self.constraintList,saveFrameName,sfKey,self.saveFrame.title)
      self.setSfTagsApplData(self.constraintList,saveFrameName,sfKey)
    
      #
      # Also add in comment and original parse error loops
      #

      tableName = self.file.sfDict[saveFrameName].preTagShort + 'constraint_comment_org'
      self.setSfTablesApplData(self.constraintList,saveFrameName,sfKey,tableList = [tableName])

      tableName = self.file.sfDict[saveFrameName].preTagShort + 'constraint_org_file_parse_err'
      self.setSfTablesApplData(self.constraintList,saveFrameName,sfKey,tableList = [tableName])
    """
    
  def setConstraintSfApplData(self):
  
    if self.version[0] == '3':
      if self.version == '3.0':
        blockKey = 'MR_file_block_position'
      elif self.version == '3.1':
        blockKey = 'Block_ID'

    else:
      return
    
    if self.constraintFile.saveFrame.tags.has_key(blockKey):

      sfKey = self.constraintFile.saveFrame.tags[blockKey]
      self.saveFrame = self.constraintFile.saveFrame
      saveFrameName = self.saveFrame.tags['Sf_category']

      self.setSfApplData(self.constraintList,saveFrameName,sfKey,self.saveFrame.title)
      self.setSfTagsApplData(self.constraintList,saveFrameName,sfKey)
    
      #
      # Also add in comment and original parse error loops
      #
      
      tableName = self.findTableByEndText(self.file.sfDict[saveFrameName],'constraint_comment_org')
      self.setSfTablesApplData(self.constraintList,saveFrameName,sfKey,tableList = [tableName])

      if self.version == '3.0':
        loopName = 'constraint_org_file_parse_err'
      elif self.version == '3.1':
        loopName = 'constraint_parse_err'
        
      tableName = self.findTableByEndText(self.file.sfDict[saveFrameName],loopName)
      self.setSfTablesApplData(self.constraintList,saveFrameName,sfKey,tableList = [tableName])
  
  def setAllSfApplData(self,parent,saveFrameName,sfKeyName):
    
    if self.file.sfs.has_key(saveFrameName):
    
      for self.saveFrame in self.file.sfs[saveFrameName]:
        
        if self.saveFrame.tags.has_key(sfKeyName):
        
          sfKey = self.saveFrame.tags[sfKeyName]
          sfTitle = self.saveFrame.title

          self.setSfApplData(parent,saveFrameName,sfKey,sfTitle)

          self.setSfTagsApplData(parent,saveFrameName,sfKey)

          self.setSfTablesApplData(parent,saveFrameName,sfKey)
       
  def setSfApplData(self,parent,saveFrameName,sfKey,sfTitle):
   
    parent.addApplicationData(Implementation.AppDataString(
          
          application = self.format,
          keyword = saveFrameName,
          value = str(sfKey)
                                ))

    parent.addApplicationData(Implementation.AppDataString(
          
          application = self.format,
          keyword = saveFrameName + '_title_' + str(sfKey),
          value = sfTitle
                                ))

  def setSfTablesApplData(self,parent,saveFrameName,sfKey,tableList = None):
  
    #
    # Default is to store all loops if tableList not given
    #
    
    tableNames = []

    for tableInfo in self.file.sfDict[saveFrameName].tables:
    
      tableName = tableInfo[0]

      if tableList:
        if tableName not in tableList:
          continue

      if self.saveFrame.tables.has_key(tableName):
        
        #
        # Also set a list of tables in application data...
        #
        
        tableNames.append(tableName)
        
        #
        # TODO: THIS IS A HACK! Need label to determine 'right' column
        # Can this come from nmrStar dictionary?!? Could in principle
        # also be multiple values...
        #
               
        tableRef = tableInfo[1][0][0]
        #if tableRef == 'Constraints_ID':
        #  tableRef = tableInfo[1][1][0]

        # Also track tag names
        parent.addApplicationData(Implementation.AppDataString(

          application = self.format,
          keyword = saveFrameName + '_table_' + str(sfKey) + tableName + '_tagNames' ,
          value = str(self.saveFrame.tables[tableName].tagNames)
                                     ))
        
        tableTags = self.saveFrame.tables[tableName].tags
        
        for i in range(0,len(tableTags[tableRef])):
        
          #tagKey = str(tableTags[tableRef][i])
          tagKey = str(i+1)
          tagList = []

          for tagInfo in tableInfo[1]:

            tagName = tagInfo[0]
        
            if tableTags.has_key(tagName):
      
              tagList.append(tableTags[tagName][i])
        
            else:
      
              tagList.append(None)

          # Track tag values
          parent.addApplicationData(Implementation.AppDataString(

            application = self.format,
            keyword = saveFrameName + '_table_' + str(sfKey) + tableName + '_' + str(tagKey),
            value = str(tagList)
                                       ))
    
    #
    # Set tablenames for this sf...
    #
                      
    if tableNames:

      parent.addApplicationData(Implementation.AppDataString(

        application = self.format,
        keyword = saveFrameName + '_tables',
        value = str(tableNames)
        
        ))
                
  def setSfTagsApplData(self,parent,saveFrameName,sfKey):
  
    #
    # Order in which tags stored in list is same as in sfDict (ccp.format.nmrStar.generalIO)
    # Leave out saveframe category (is known)
    #
    
    tagList = []
    tagNames = []
    
    for tagInfo in self.file.sfDict[saveFrameName].tags[1:]:
    
      tagName = tagInfo[0]
      tagNames.append(tagName)
      
      if self.saveFrame.tags.has_key(tagName):
      
        tagList.append(self.saveFrame.tags[tagName])
        
      else:
      
        tagList.append(None)
    
    parent.addApplicationData(Implementation.AppDataString(
          
          application = self.format,
          keyword = saveFrameName + '_tagNames_' + str(sfKey),
          value = str(tagNames)
                                ))
    
    parent.addApplicationData(Implementation.AppDataString(
          
          application = self.format,
          keyword = saveFrameName + '_tags_' + str(sfKey),
          value = str(tagList)
                                ))

  def setSfFromApplData(self,applData,saveFrameName,parent):
  
      sfKey = applData.value
      
      keyword = saveFrameName +  '_title_' + sfKey
      sfTitle = parent.findFirstApplicationData(application = self.format, keyword = keyword).value

      self.setupSaveFrame(saveFrameName,sfTitle)
  
      keyword = saveFrameName + '_tags_' + sfKey

      sfTagValues = eval(parent.findFirstApplicationData(application = self.format, keyword = keyword).value)
          
      for i in range(1,len(self.file.sfDict[saveFrameName].tags)):
        
        tagInfo = self.file.sfDict[saveFrameName].tags[i]
        tagName = tagInfo[0]

        tagValue = sfTagValues[i-1]
    
        self.setSfTag(tagName,tagValue)

  def setLoopsFromApplData(self,parent,saveFrameName,sfKey):
  
    applData = parent.findFirstApplicationData(application = self.format, keyword = saveFrameName + '_tables')

    if applData:

      tableNames = eval(applData.value)

      for tableName in tableNames:

        #
        # Get reference info
        #

        tableInfo = self.findTableInfo(saveFrameName,tableName)
        
        if not tableInfo:

          print "  Warning: no reference data for saveframe %s, table %s" % (saveFrameName,tableName)
          continue

        index = 1
        self.setupTable(tableName)

        keyword = saveFrameName + '_table_' + str(sfKey) + tableName + '_' + str(index)
        applData = parent.findFirstApplicationData(application = self.format, keyword = keyword)

        while (applData):

          tableTagValues = eval(applData.value)

          for i in range(0,len(tableTagValues)):

            tagName = tableInfo[i][0]
            tagValue = tableTagValues[i]

            self.setTableTag(tagName,tagValue)

          #
          # Get next line
          #

          index += 1
          keyword = saveFrameName + '_table_' + str(sfKey) + tableName + '_' + str(index)
          applData = parent.findFirstApplicationData(application = self.format, keyword = keyword)

  def findTableInfo(self,saveFrameName,tableName):

    for (tempTableName,tempTableInfo) in self.file.sfDict[saveFrameName].tables:

      if tableName == tempTableName:

        return tempTableInfo
        
    return None

  def setConversionErrorsLoop(self,saveFrameName,constraintListID,conversionErrors):
        
    if conversionErrors:
    
      tableName = self.findTableByEndText(self.file.sfDict[saveFrameName],'constraint_parse_file')

      self.setupTable(tableName)
      self.setTableTag('Constraints_ID',constraintListID)
      self.setTableTag('ID','1')
      
      if hasattr(self,'originalFileName'):
        name = self.originalFileName
      else:
        name = self.entry.name
      
      self.setTableTag('Name',name)

      self.setupTable(tableName)

      tableInfo = self.findTableInfo(saveFrameName,tableName)

      for i in range(len(conversionErrors)):

        conversionError = conversionErrors[i]

        self.setTableTag('Constraints_ID',constraintListID)
        self.setTableTag('ID',i+1)

        for j in range(0,len(conversionError)):

          self.setTableTag(tableInfo[j+2][0],conversionError[j])
  
  def findTableByEndText(self,saveFrame,tableEndText):
  
    tableName = None
    
    for tableInfo in saveFrame.tables:
      if tableInfo[0].endswith(tableEndText):
        tableName = tableInfo[0]
        break
        
    return tableName

  def setCcpnMolTypes(self):
      
    #
    # Reset original polymerType!
    #

    for sequence in self.sequenceFile.sequences:
      if hasattr(sequence,'polymerType'):
        if sequence.polymerType in self.ccpn2NmrStarMolTypes.values():
          for molType in self.ccpn2NmrStarMolTypes.keys():
            if self.ccpn2NmrStarMolTypes[molType] == sequence.polymerType:
              sequence.polymerType = molType
        # Hack, otherwise too much hassle
        elif sequence.polymerType == 'polypeptide(D)':
          sequence.polymerType = 'protein'

  def getPresetChainMapping(self,chainList):
  
    mappingChainDict = {}
    
    for chain in chainList:
    
      mappingChainDict[chain] = (chain.code,1)
    
    return mappingChainDict
