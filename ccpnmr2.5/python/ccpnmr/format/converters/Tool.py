
"""
======================COPYRIGHT/LICENSE START==========================

Tool.py: Tool class specified here simplifies use of format conversion code

Copyright (C) 2005-2010 Wim Vranken (European Bioinformatics Institute)

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

from memops.universal.Util import returnInt
from memops.universal.Io import splitPath

from memops.general.Util import returnMemopsWord
from memops.general.Io import loadProject

from ccpnmr.format.general.Util import createExperiment, getRefExpFromOldExpType
from ccpnmr.format.general.Util import createPpmFreqDataSource
from ccpnmr.format.general.Constants import defaultMolCode, allFormatsDict
from ccpnmr.format.process.makeChemShiftsFromPeakLists import createChemShifts

from memops.api import Implementation
from ccp.api.nmr import Nmr
from ccp.api.molecule import MolSystem

class Tool:

  """
  
  Is essentially a wrapper around all the ...Format.py files in this
  directory. Allows for easier import/creation of data model objects,
  as well as easier writing.
  
  See end of this file for an example on how to use this class
  
  """

  def __init__(self, guiParent = None):
    
    self.fileNameDict = {}
    self.formatDict = {}

    self.chains = []

    #
    # Trying to load GUI environment: is used by default
    # Bit nasty to put import here but easier.
    #
    
    if not guiParent:

      try:

        import Tkinter
        self.guiParent = Tkinter.Tk()
  
      except:

        self.guiParent = None
    
    else:
    
      self.guiParent = guiParent

  def setupProject(self,dirPath,projectName):
  
    self.projectDirPath = dirPath

    try:
      self.ccpnProject = loadProject(dirPath)
      print "Read existing CCPN project..."
      
    except:
      print "Creating new CCPN project..."
      
      projectName = returnMemopsWord(projectName)
      self.ccpnProject = Implementation.MemopsRoot(name = projectName)
      self.nmrProject = self.ccpnProject.newNmrProject(name = self.ccpnProject.name)
      
      if not os.path.exists(dirPath):
        os.mkdir(dirPath)
      
      repository = self.ccpnProject.findFirstRepository(name = 'userData')
      
      if not repository.url.path == dirPath:
        url = Implementation.Url(path=dirPath) 
        repository.url = url
    
    # TODO: check name of project: if not same CRASH OUT!
          
    return self.ccpnProject
    
  def setupMolSystem(self,molSysName):
  
    self.molSystem = self.ccpnProject.findFirstMolSystem(code = molSysName)
    
    if not self.molSystem:

      self.molSystem = MolSystem.MolSystem(self.ccpnProject,code = molSysName)

  def setupFormat(self, label):
  
    if self.formatDict.has_key(label):
      self.format = self.formatDict[label]
      
    else:
      formatModule = __import__('ccpnmr.format.converters.%sFormat' % label,{},{},['%sFormat' % label])
      formatClass = getattr(formatModule,'%sFormat' % label)
      self.format = formatClass(self.ccpnProject,self.guiParent,verbose = 1)
      self.formatDict[label] = self.format
    
    return self.format
    
  def readSequence(self,fileName,
                        moleculeName,
	    moleculeType,
	    moleculeSegID,
	    moleculeFirstResidueNumber):
	    	    
    #
    # First do a simple check to see if molecule/chain already exists (with the relevant segID!)
    #
    
    molecule = self.ccpnProject.findFirstMolecule(name = moleculeName)
    
    if molecule:
      chain = self.molSystem.findFirstChain(code = moleculeSegID, molecule = molecule)
    else:
      chain = None

    if not chain:
    
      chains = self.format.readSequence(fileName,
                                        molSystem = self.molSystem,
		molName = moleculeName,
                                        seqCodeStart = moleculeFirstResidueNumber,
                                        chainCodes = [moleculeSegID],
                                        minimalPrompts = 1)
      
      self.chains.extend(chains)
    
    else:
      
      if not chain in self.chains:
        self.chains.append(chain)

    return self.chains
    
  def readPeakList(self,fileName,refExperimentOldName,spectrumName,peakListDimOrder):
    
    if self.fileNameDict.has_key(fileName):
      self.peakList = self.fileNameDict[fileName]
      
    else:
    
      refExperiment = getRefExpFromOldExpType(self.ccpnProject,refExperimentOldName)
      
      self.experiment = self.ccpnProject.currentNmrProject.findFirstExperiment(name = spectrumName, refExperiment = refExperiment)
      
      if not self.experiment:
        self.experiment = createExperiment(self.ccpnProject, spectrumName, refExperiment)
      
      self.dataSource = self.experiment.findFirstDataSource(dataType = 'processed')

      if not self.dataSource:
        self.dataSource = createPpmFreqDataSource(self.experiment,spectrumName,'processed',self.experiment.numDim)
      
      self.peakList = self.dataSource.findFirstPeakList(name = spectrumName)
        
      #
      # Determine the peak list relationship
      # In ccp data model setup assuming that:
      #   - first dim is acquisition
      #   - first hetero dim belongs to acquisition dim
      #

      protonDims = []
      heteroDims = []
      
      for dataDim in self.dataSource.sortedDataDims():
      
        if '1H' in dataDim.findFirstDataDimRef().expDimRef.isotopeCodes:
    
          protonDims.append(dataDim.findFirstDataDimRef())
        
        else:
      
          heteroDims.append(dataDim.findFirstDataDimRef())
      
      dataDimRefs = self.experiment.numDim * ['']
      
      for dim in range(0,len(peakListDimOrder)):

        peakListDim = returnInt(peakListDimOrder[dim])
      
        if peakListDim:

          if dim >= 2:
            dataDimRef = heteroDims[dim - 2]
          else:
            dataDimRef = protonDims[dim - 2]
        
          dataDimRefs[peakListDim-1] = dataDimRef
      
      self.peakLists = self.format.readPeaks(fileName,
                                            dataSource = self.dataSource,
	                                    dataDimRefs = dataDimRefs,
		                            specName = spectrumName,
		                            peakList = self.peakList,
	                                    minimalPrompts = 1)
    
    if self.experiment not in self.chemShiftList.experiments:
      self.chemShiftList.addExperiment(self.experiment)

    return self.peakLists

  def readPeakListAssignments(self,fileName):
    
    if self.fileNameDict.has_key(fileName):
    
      self.peakList = self.fileNameDict[fileName]
      
    else:
    
      self.peakList = self.format.readPeakAssignments(fileName,
                                                      self.peakList)
		              
      self.fileNameDict[fileName] = self.peakList


    return self.peakList

  def readShifts(self,fileName,chemShiftListName):
    
    #
    # ALWAYS read chemical shift list BEFORE peak list!
    #
    
    if self.fileNameDict.has_key(fileName):
      self.chemShiftList = self.fileNameDict[fileName]
      
    else:
      
      self.chemShiftList = self.ccpnProject.currentNmrProject.findFirstMeasurementList(className = 'ShiftList', name = chemShiftListName)
      
      self.chemShiftList = self.format.readShifts(fileName,
                                                  measurementListName = chemShiftListName,
		                                                measurementList = self.chemShiftList,
 	                                                minimalPrompts = 1)

    return self.chemShiftList
  
  def createChemShiftsFromPeakList(self,chemShiftListName):
  
    shiftList = createChemShifts([self.peakList], guiParent = self.guiParent, shiftListName = chemShiftListName)
    
    return shiftList
    
  def linkResonances(self,namingSystemName,**keywds):
    
    #
    # TODO: set nameMapping in advance?
    #
    
    if namingSystemName == 'auto' or namingSystemName == '':
      namingSystemName = None

    self.format.linkResonances(namingSystemName = namingSystemName,**keywds)

  def writeXmlProject(self):
  
    self.ccpnProject.saveModified()


if __name__ == '__main__':

  ccpnConv = Tool()
  
  #mainDir = "/ebi/msd/nmrqual/"
  mainDir = "/Users/wim/"

  #
  # CCPN project:
  #
  #   Loads an existing project if fileName corresponds to ccpn project xml file
  #   Creates a new project if fileName unknown or wrong format
  #
  # Note that you should only create a new project file when starting a set
  # of structure calculations for a protein: otherwise information about
  # previous peak lists and assignments will be lost.
  #

  projectName = "My_project"   # This info can come from XML file...
  dirPath = os.path.join(mainDir,"workspace/stable/all/python/ccpnmr/format/converters/local/testTool/myProject")

  ccpnConv.setupProject(dirPath,projectName)

  #
  # Define molSystem! No white spaces allowed in this name.
  #

  molSysName = "MyMolecularSystem"

  ccpnConv.setupMolSystem(molSysName)

  #
  # Define sequence (TODO: also need to read 'non-standard' PDB for residues!)
  #

  sequenceFormatName = "XEasy"

  sequenceFileName = os.path.join(mainDir,"reference/xeasy/werner.seq")
  moleculeName = "My molecule"
  moleculeType = "protein"
  moleculeSegID = "A"
  moleculeFirstResidueNumber = 1
  #number chains? Or interprete that from input??

  ccpnConv.setupFormat(sequenceFormatName)

  ccpnConv.readSequence(sequenceFileName,
                        moleculeName,
	    moleculeType,
	    moleculeSegID,
	    moleculeFirstResidueNumber)

  #
  # Define experiments/spectra/peaklists
  #

  dataFilesFormatName = "XEasy"

  dataFilesFormat = ccpnConv.setupFormat(dataFilesFormatName)

  # Has to be one of: noesy.hh, noesy_hsqc_HNH.hhn,
  # noesy_hsqc_HCH.hhc, noesy_hsqc_HNNH.hhnn,
  # noesy_hsqc_HCNH.hhcn (?) (always use hsqc here - doesn't
  # matter if different experiment type)
  # Could also use a dict for names (noesy's only)
  
  spectrumType = "noesy_hsqc_HNH.hhn" 
  spectrumName = "noehsqc"

  chemicalShiftFileName = os.path.join(mainDir,"reference/xeasy/werner.noehsqc.prot")

  ccpnConv.readShifts(chemicalShiftFileName,spectrumName)

  peakListFileName = os.path.join(mainDir,"reference/xeasy/werner.noehsqc.peaks")
  peakListProton1 = 3
  peakListProton2 = 2
  peakListHetero1 = 1
  peakListHetero2 = None

  peakListDimOrder = [peakListProton1,peakListProton2,peakListHetero1,peakListHetero2]

  ccpnConv.readPeakList(peakListFileName,
                        spectrumType,
  	    spectrumName,
  	    peakListDimOrder)
  
  
  
  spectrumType = "noesy_hsqc_HCH.hhc" 
  spectrumName = "hmqcnoe"

  chemicalShiftFileName = os.path.join(mainDir,"reference/xeasy/werner.hmqcnoe.prot")

  ccpnConv.readShifts(chemicalShiftFileName,spectrumName)

  peakListFileName = os.path.join(mainDir,"reference/xeasy/werner.hmqcnoe.peaks")
  peakListProton1 = 3
  peakListProton2 = 2
  peakListHetero1 = 1
  peakListHetero2 = None

  peakListDimOrder = [peakListProton1,peakListProton2,peakListHetero1,peakListHetero2]

  ccpnConv.readPeakList(peakListFileName,
                        spectrumType,
  	    spectrumName,
	    peakListDimOrder)


  namingSystemName = 'auto'

  #ccpnConv.linkResonances(namingSystemName)

  ccpnConv.writeXmlProject()
  
