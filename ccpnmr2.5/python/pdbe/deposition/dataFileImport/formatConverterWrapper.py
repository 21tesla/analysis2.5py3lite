"""
======================COPYRIGHT/LICENSE START==========================

formatConverterWrapper.py: Code wrapper to handle import of files for PDB/BMRB deposition

Copyright (C) 2010 Wim Vranken (PDBe, EBI)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.

A copy of this license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wim@ebi.ac.uk
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

import sys, os

from memops.api import Implementation

from ccpnmr.format.general.Conversion import FormatConversion

from ccpnmr.format.process.sequenceCompare import SequenceCompare

# Needs to connect to GUI, this for testing
#
#
# ALWAYS one sequence/coordinate file, one chemical shift file!
#


class DepositionImportError(StandardError):

  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

class FormatConverterWrapper:

  """
  
  
  Restrictions:

  - Will not work with complex carbohydrates or multiple ligands connected to each other.
  
  
  """
  
  coordinateFormatNames = ['auremol','charmm','cns','cyana','dyana','molmol','nmrStar','pdb','pseudoPdb']

  formatNameLists = {
  
    'shifts': 

       # 'cosmos',
       ['auremol','autoAssign','cns','csi','mars','monte','nmrStar','nmrView','pipp','pistachio','pronto','shiftx','sparky','talos','xeasy'],
  
    'coordinates':
    
       # 'cosmos'
       coordinateFormatNames,
       
    'chemComps':
    
       ['mol2','pdb'], # Only these, need atom names! Pdb overlap is OK - different parser.
       
    'sequence':
  
       # 'mol2',
       coordinateFormatNames[:] + ['ansig','aria','autoAssign','csi','fasta','mars','monte','nmrView','pipp','pistachio','pronto','shiftx','sparky','talos','xeasy'],
    
    'distanceConstraints':
    
       ['cns','cyana'],
       
    'rdcConstraints':
    
       ['cns','cyana'],
       
    'dihedralConstraints':
    
       ['cns','cyana'],
       
  }

  def __init__(self, ccpnProjectName=None, ccpnProject=None, guiRoot=None):
  
    # TODO: this could start from existing project, for example when reading in second chain!
    self.formatConversion = FormatConversion(ccpnProjectName=ccpnProjectName, ccpnProject=ccpnProject,
                                             silentRead=True,guiRoot=guiRoot)
    
    self.formatNamesAtImport = {}

    self.importReturns = {}
    self.importConversion = {}
    self.importSuccess = {}

    self.newResonances = {}
    
    self.sequenceComparison = SequenceCompare()
  
  def readChemicalShiftFile(self,formatName,filePath):
  
    dataType = 'shifts'
    
    return self.readFile(dataType,formatName,filePath)

  def readSequenceFile(self,formatName,filePath):
  
    dataType = 'sequence'
    
    return self.readFile(dataType,formatName,filePath)

  def readLigandFile(self,formatName,filePath):
  
    # TODO doesn't work yet in Conversion.py!
    dataType = 'chemComps'
    
    return self.readFile(dataType,formatName,filePath)

  def readCoordinateFile(self,formatName,filePath):
  
    dataType = 'coordinates'

    return self.readFile(dataType,formatName,filePath)

  def readFile(self,dataType,formatName,filePath,addKeywords=None):
    
    # Track name of format name used for import of this file type
    self.formatNamesAtImport[dataType] = formatName
    
    if not addKeywords:
      addKeywords = {}
    
    # Hacks to handle coordinates - handling these files is a mess in the formatConversion class!
    preparseFilePath = filePath
    if dataType == 'coordinates':
      addKeywords['autoCreateChemComps'] = True
      preparseFilePath = [filePath]

    (fileRead,fileInformation) = self.formatConversion.preparseFile(dataType,formatName,preparseFilePath)

    self.importReturns[dataType] = self.formatConversion.importFile(dataType,formatName,filePath,addKeywords=addKeywords)

    self.ccpnObjectOrList = self.importReturns[dataType]
    self.conversionInfo = self.formatConversion.conversionInfo
    self.conversionSuccess = self.formatConversion.conversionSuccess

    # Set info for export to NMR-STAR
    self.setNmrStarExportInfo(formatName)
    
    # Track newly created resonances
    if dataType not in ('sequence','coordinates','chemComps'):
      self.addNewResonances(formatName)
      
    return (fileRead,fileInformation)

  def setNmrStarExportInfo(self,importFormatName):
  
    resonances = self.formatConversion.getFormatClass(importFormatName).newResonances

    #
    # Copy over original assignment info for NMR-STAR export!
    #
    
    for resonance in resonances:
      for appData in resonance.findAllApplicationData(application=importFormatName):
        newAppData = None
        if appData.keyword == 'assign':
          newAppData = Implementation.AppDataString(application='nmrStar',keyword='origAssign',value=appData.value)
        elif appData.keyword == 'origResLabel':
          newAppData = appData
          
        if newAppData:
          resonance.addApplicationData(newAppData)

  def addNewResonances(self,importFormatName):
  
    newResonances = self.formatConversion.getFormatClass(importFormatName).newResonances
    
    if newResonances:
    
      resonanceParent = newResonances[0].parent
      
      if resonanceParent not in self.newResonances.keys():
        self.newResonances[resonanceParent] = {}
      if importFormatName not in self.newResonances[resonanceParent].keys():
        self.newResonances[resonanceParent][importFormatName] = []
      
  
      for newResonance in newResonances:
        if newResonance not in self.newResonances[resonanceParent][importFormatName]:
          self.newResonances[resonanceParent][importFormatName].append(newResonance)
    
    
  def determineFileInfo(self, filePath):
    """ Analyse file and file name for readable file and return info dictionary
    
    NB currently only for restrains and coordinates as in the dataTypes parameter
    Problem is that sequence reading will interpret almost any text as a sequence and try to
    download ChemComps name '<DocType' etc.
    Rasmus Fogh 26/6/2013
    """
    
    constraintTypes = ('distanceConstraints', 'dihedralConstraints', 
                       'rdcConstraints', )
    coordinateFormats = ('pdb', 'pseudoPdb')
    ignoreExtensions = ('.xml',)
    
    
    result = {}
    
    # Set up
    fileName = os.path.basename(filePath)
    head, ext = os.path.splitext(fileName)
  
    result['name'] = head
    
    # lower case to simplify tests below
    head = head.lower()
    ext = ext.lower()
    
    # set dataTypes and formatNames from extension
    dataTypes = None
    formatNames = None
    
    if ext in ('.coord', '.pdb'):
      dataTypes = ('coordinates',)
      formatNames = coordinateFormats
      
    elif ext in ('.upl', '.lol'):
      dataTypes = ('distanceConstraints',)
      formatNames = ('cyana',)
 
    elif ext == '.aco':
      dataTypes = ('dihedralConstraints',)
      formatNames = ('cyana',)
 
    elif ext == '.tbl':
      formatNames = ('cns',)
 
      if 'dihe' in head:
        dataTypes = ('dihedralConstraints',)
      
      elif 'rdc' in head:
         dataTypes = ('rdcConstraints',)
         
      else:
        for tag in ('ambig', 'unambig', 'noe', 'dist'):
          if tag in head:
            dataTypes = ('distanceConstraints',)
            break
        else:
          dataTypes = constraintTypes
    
    else:
      # extension not recognised
      # do nothing for now
      # Later add shifts, sequences, peaks ... here
      pass
      
    # Check if any of the proposed dataType,formtName combos work
    fileRead = None
    if dataTypes and formatNames:
      # extension shows which format(s) to check
      for dataType in dataTypes:
        for formatName in formatNames:
          fileRead,fileInformation = self.formatConversion.preparseFile(
                                            dataType, formatName, filePath)
          if fileRead:
            break
        else:
          continue
        break
        
    if fileRead:
      # success above. Set result
        result['dataType'] = dataType
        result['formatName'] = formatName
          
    #
    return result
      
  
  def determineFormatNamesForFile(self,dataType,filePath):
  
    formatNameSuggestions = self.formatConversion.determineFormatNamesForFile(dataType,filePath, formatNameList = self.formatNameLists[dataType])
  
    return formatNameSuggestions
    
  def linkAllResonancesToAtoms(self):
  
    linkingInfo = {}
  
    # Links everything, will only work for monomers at this stage
    resonanceParents = self.newResonances.keys()
    resonanceParents.sort()
        
    for resonanceParent in resonanceParents:
    
      importFormatNames = self.newResonances[resonanceParent].keys()
      importFormatNames.sort()
      
      for importFormatName in importFormatNames:
        
        origUnlinked = -999
        
        """
        TODO: Deal with problem of files with different chain codes; could run first on distance stuff (residue
        info), then force same mapping for other chaincodes... best to make sure original format is fixed though?
        """
        forceChainMappings = self.linkResonancesToSequence(resonanceParent=resonanceParent,importFormatName=importFormatName,allowMultipleFormatChains=True)
        
        # This tries to apply a chain mapping from another format to this one; bit dangerous
        if not forceChainMappings:
          if existingForceChainMappings:
            forceChainMappings = self.setForceChainMappingsFromPreviousMapping(existingForceChainMappings,resonanceParent=resonanceParent,importFormatName=importFormatName)
        else:
          existingForceChainMappings = forceChainMappings
        
        if forceChainMappings:
        
          chainMappingResetAttempted = False
        
          while True:

            self.formatConversion.linkResonances(forceChainMappings=forceChainMappings,setSingleProchiral=False,setSinglePossEquiv=False)        
            numResonancesLinked = self.formatConversion.numResonancesLinked
            linkingInfo[(resonanceParent,importFormatName)] = numResonancesLinked

            if numResonancesLinked['unlinked'] == 0:
              break
            elif self.formatConversion.numResonancesLinked['origUnlinked'] == origUnlinked:
              # Here trying to apply a working chain mapping to one where there is not enough information, within the same format.
              if existingForceChainMappings:
                forceChainMappings = self.setForceChainMappingsFromPreviousMapping(existingForceChainMappings,resonanceParent=resonanceParent,importFormatName=importFormatName)
                
              if chainMappingResetAttempted or not forceChainMappings:
                break
              else:
                print("  WARNING: trying new chain mapping {}".format(str(forceChainMappings)))
                chainMappingResetAttempted = True
              
            origUnlinked = self.formatConversion.numResonancesLinked['origUnlinked']

    return linkingInfo
        
  def setForceChainMappingsFromPreviousMapping(self, existingForceChainMappings, resonanceParent=None, importFormatName=None):
  
    self.setSequenceComparisonFormatFileInfo(resonanceParent,importFormatName)
    
    # Make sure to exclude format chain codes in the existing mapping; are picked up again in above method
    existingMappedFormatChainCodes = []
    for (ccpnChainCode,existingFormatChainCode,ccpnSeqId,offset) in existingForceChainMappings:
      if existingFormatChainCode not in existingMappedFormatChainCodes:
        existingMappedFormatChainCodes.append(existingFormatChainCode)

    # Now try to find out whether a new chain mapping can be created based on the existing one.
    forceChainMappings= []

    formatChainCodes = self.sequenceComparison.formatFileResidueDict.keys()
    for formatChainCode in formatChainCodes[:]:
      if formatChainCode in existingMappedFormatChainCodes:
        formatChainCodes.pop(formatChainCodes.index(formatChainCode))

    if len(formatChainCodes) == 1:
      for (ccpnChainCode,existingFormatChainCode,ccpnSeqId,offset) in existingForceChainMappings:
        forceChainMappings.append((ccpnChainCode,formatChainCodes[0],ccpnSeqId,offset))
    else:
      print("   Warning: cannot use previous mapping, multiple chain codes in format.")
      
    return forceChainMappings
  
  def setSequenceComparisonFormatFileInfo(self,resonanceParent,importFormatName):
  
    #
    # New resonances for shift list are tracked by FC
    #
    
    if not resonanceParent:
      if self.newResonances.keys():
        resonanceParent = self.newResonances.keys()[0]
    
    assert resonanceParent in self.newResonances.keys(), "No resonance parent class defined!"
    
    if not importFormatName:
      if self.newResonances[resonanceParent].keys():
        importFormatName = self.newResonances[resonanceParent].keys()[0]

    resonances = self.newResonances[resonanceParent][importFormatName]
    self.numNewResonances = len(resonances)
    
    #print importFormatName, resonances

    #
    # Get info from the resonances
    #
    
    self.sequenceComparison.getFormatFileInformation(resonances,importFormatName)

  def linkResonancesToSequence(self, chain=None, resonanceParent=None, importFormatName=None, allowMultipleFormatChains=False):
  
    """
    Note: taken from ccpnmr.format.process.matchResonToMolSys
    """
    
    forceChainMappings = {}
    
    
    self.setSequenceComparisonFormatFileInfo(resonanceParent,importFormatName)
    
    #
    # Now generate the sequence information (one-letter codes for polymer, three letter otherwise) for the chain(s) in CCPN
    #
    
    if not chain:
      
      # 2013JUn26 Rasmus Fogh
      # Changed to allow use of pre-exiating MolSystem, e.g. if reading into existing project
    
      #if 'sequence' not in self.importReturns.keys() and 'coordinates' not in self.importReturns.keys():
      #  raise DepositionImportError('Sequence information missing, cannot continue!')

      if 'sequence' in self.importReturns.keys():
        chains = self.importReturns['sequence']
      else:
        chains = self.formatConversion.ccpnProject.currentMolSystem.sortedChains()

      # Rasmus Fogh addition
      if not chains:
        raise DepositionImportError('Sequence information missing, cannot continue!')
      
      # TODO REMOVE THIS WHEN BECOMES POSSIBLE!
      if len(chains) > 1:
        raise DepositionImportError('Multiple chains created during import, cannot continue!')
    
    else:
      
      chains = [chain]

    self.sequenceComparison.createCcpnChainInformation(chains)
        
    #
    # Now generate the sequence (one-letter codes) for the information from the chemical shift file
    #
    
    if not allowMultipleFormatChains and len(self.sequenceComparison.formatFileChainDict) > 1:
      raise DepositionImportError('Multiple format chain codes in imported chemical shift file, cannot continue!')
    
    self.sequenceComparison.createFormatFileChainInformation()
      
    #
    # Now run the comparison...
    #
    
    forceChainMappings = self.sequenceComparison.compareFormatFileToCcpnInfo()

    print "\n*** Chain mappings set by alignment information ***\n"
    print forceChainMappings
      
    #
    # Reset list of resonances for FormatClass!
    #

    self.formatConversion.getFormatClass(importFormatName).newResonances = []

    
    return forceChainMappings

if __name__ == '__main__':

  """
  
  AT IMPORT:
  
  - either upload a sequence (protein/DNA/RNA or known ligand code) or a .mol2/PDB file (with atom names matching shifts)
    
  - then load shifts, ONLY for this chain
  
    --> for both, give list of valid format names on upload
  
  - repeat if necessary (with tabs?)
  
  - when finalised, give option to load coordinates (?). Ligand codes HAVE TO MATCH!
  
  
  
  
  
  - RENAME TO SOMETHING CLEARER!
  
  - Have to do preparse, use FormatConversion code, to validate file
     --> Popup with nice description if not good!

  - Might need to work on linkResonances popups - are the texts clear enough? Can be better?

  - Have to work off some crap files, INCLUDE NMR-STAR for sure as well!!
  
  
  """
  
  pass
