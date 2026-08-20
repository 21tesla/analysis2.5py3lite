"""
======================COPYRIGHT/LICENSE START==========================

AnsigFormat.py: Contains functions specific to Ansig conversions.

Copyright (C) 2008-2011 Wim Vranken (Vrije Universiteit Brussel)

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

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords, cleanAtomNameLists
from ccpnmr.format.general.Util import getResName
from ccpnmr.format.general.Util import getResNameText

from ccpnmr.format.general.Constants import volume_kw


import ccp.api.nmr.Nmr as Nmr

from ccp.format.general.Util import getSeqAndInsertCode

import string, copy

# This is necessary for handling the parameter files!
from ccpnmr.format.converters.AzaraFormat import AzaraFormat

#
# Format specific keywords...
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readProject']['peakFile'] = (None,False,'The peak file in export or storage format associated with this project.')
IOkeywords['writePeaks']['forceNumDim'] = (None,False,'Force the number of dimensions in the output file.')

class AnsigFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'ansig'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    self.createSequenceFile = self.createSequenceFileGeneric
    
    self.getPeaks = self.getPeaksGeneric
    self.createPeakFile = self.createPeakFileGeneric
    
  #
  # Deviations from generic import stuff
  #

  def createPeakFileFormatSpecific(self):
    
    if self.forceNumDim:
      numDim = self.forceNumDim
    else:
      numDim = self.numPeakDim

    self.peakFile.setSpectrumInfo('Export',0,numDim)

  def getFullProject(self,fileName,peakFile = None,peakKeyWdList = None, sequenceKeyWds = None):
  
    if self.verbose == 1:
      print "Reading %s project from file %s" % (self.formatLabel,fileName)
   
    self.file = self.projectIO.AnsigProjectFile(fileName)
    self.file.read()
    
    #
    # Ask if peakFile has to be read if not given...
    #
    
    if not peakFile and not self.minimalPrompts:
    
      getPeakFile = self.messageReporter.showYesNo("Read a peak file","Do you want to read in an associated peak file (in export or storage format)?")
      
      if getPeakFile:
      
        interaction = self.multiDialog.FileName(self.guiParent,component = 'peaks', format = self.format, selectionText = "Select peak file", title = 'Select peak file')
        peakFile = interaction.file
    
    #
    # Process elements - call functions 'as if' reading real file, but pass
    # on objects created during project import
    #
        
    if self.file.files['sequence']:
      
      if sequenceKeyWds:
        keywds = sequenceKeyWds
      else:
        keywds = {}
        
      self.readSequence(fileName = self.file.files['sequence'], minimalPrompts = self.minimalPrompts, **keywds)
    
    if self.file.spectraFile:
    
      spectraList = self.file.spectraFile.spectra
      azaraFormat = AzaraFormat(self.project,self.guiParent)

      for spectrumInfo in spectraList:
      
        #
        # TODO: use 'solvent',  'temperature', 'pH' to set up sample info!
        #
        
        specName = spectrumInfo['name']
        if hasattr(spectrumInfo,'type'):
          specType = spectrumInfo['type']
        else:
          specType = None
        numDim = spectrumInfo['numDim']
        parFile = spectrumInfo['parFile']
        
        #
        # Use ppmRange if there is an aliased dimension...
        # assuming that dimension order in ANSIG file same as in Azara one.
        #
        
        aliasing = {}
        
        for dim in range(0,len(spectrumInfo['dims'])):
          dimName = spectrumInfo['dims'][dim]
          if spectrumInfo['dimInfo'][dimName].has_key('aliased'):
            nucleus = spectrumInfo['dimInfo'][dimName]['nucleus']
            if self.file.spectraFile.ppmRange.has_key(nucleus):
              aliasing[dim] = self.file.spectraFile.ppmRange[nucleus]
        
        if not self.fileExists(parFile, verbose = False):
          print "  Warning: parameter file %s does not exist - reading skipped" % parFile
          continue
            
        dataSource = azaraFormat.readProcPars(fileName = parFile, expName = specName, aliasing = aliasing, minimalPrompts = self.minimalPrompts)
        
        if dataSource:
        
          #
          # Reset some values. DO NOT USE the mapping information - 
          #

          #dataDimRefs = azaraFormat.dataDimRefs        
          azaraFormat.experiment = None
          
          # Force this if it was set...
          if self.file.spectraFile.binaryFormat == 'bigendian':
            dataSource.isBigEndian = True

          #
          # Now try to read in the relevant peak information
          #

          if peakFile:
            # dataDimRefs = dataDimRefs,
            peakLists = self.readPeaks(fileName = peakFile, dataSource = dataSource, specName = specName, minimalPrompts = self.minimalPrompts)
    
    #
    # Avoid reading project again
    #
    
    self.file = None
    
    return self.project
    

  #
  # Functions different to default functions in DataFormat
  #

  def createSequence(self):
  
    self.sequenceFile.sequences.append(self.sequenceIO.AnsigSequence(molName = self.chain.molecule.name, chainCode = self.exportChainCode))
    self.sequence = self.sequenceFile.sequences[-1]

  def setSequenceFileElements(self):
  
    #
    # TODO: get 'details' information... from resGroup in principle.
    #
        
    self.sequence.elements.append(self.sequenceIO.AnsigSequenceElement(self.seqCode,self.residue.molResidue.chemComp.ccpCode))
    
  def thisPeakValid(self):
    
    #
    # Check whether peak from correct spectrum
    # 

    if self.specName != self.rawPeak.specName:
      return 0
    
    else:
      return 1

  def setPeakIntensity(self):
        
    # PeakIntensity attributes
    if self.rawPeak.volume != 0:
      peakInt = Nmr.PeakIntensity(self.peak,value = self.rawPeak.volume, method = self.methods[self.format]['Volume'])

      # Assuming always volume for ansig
      peakInt.intensityType = volume_kw

  def setPeakDim(self):
  
    dataDimRef = self.dataDimRefs[self.rawPeakDimIndex]

    self.peakDim = self.peak.findFirstPeakDim(dim = dataDimRef.dataDim.dim)

    self.peakDim.dataDimRef = dataDimRef
    
    self.peakDim.value = self.rawPeak.ppm[self.rawPeakDimIndex]

  def setPeakExtras(self):

    if self.rawPeak.deleted == 1:
      # TODO: is this the correct way to 'delete' a peak?
      self.peak.delete = True
    
  def getPeakResNames(self):
  
    if self.rawPeak.assignResNum[self.rawPeakDimIndex] != '':
      (seqCode,seqInsertCode) =  getSeqAndInsertCode(self.rawPeak.assignResNum[self.rawPeakDimIndex])
      self.resNames =  [getResName(' ',seqCode,self.rawPeak.assignAtom[self.rawPeakDimIndex],seqInsertCode = seqInsertCode)]
    
    else:
      self.resNames = []
      
  def selectPeakLists(self):
    
    #
    # Use all peaklists for Ansig (can handle multiple per file)
    # Have to find maximum number of dims first
    #
    
    self.numPeakDim = 0
    
    for peakList in self.peakLists:
      
      numPeakDim = 0
      for dataDim in peakList.dataSource.dataDims:
        if dataDim.className != 'SampledDataDim':
          numPeakDim += 1

      if self.numPeakDim < numPeakDim:
        self.numPeakDim = numPeakDim

  def setPeakFileInfo(self):

    self.numDim = len(self.peakList.sortedPeaks()[0].peakDims)
    self.peakFile.numDims.append(self.numDim)
    
    #
    # Limit to 12 chars... should really have popup or something
    #

    specName = None
    
    if self.peakList.name and not self.peakFile.specNames.count(self.peakList.name):
      specName = self.peakList.name[:12]

    elif self.minimalPrompts:
      specName = origSpecName = "%s_%d.%d" % (self.peakList.dataSource.experiment.name[:8],self.peakList.dataSource.serial,self.peakList.serial)
      startId = 1
      while (self.peakFile.specNames.count(specName)):
        specName = "%s_%d" % (origSpecName[:10],startId)
        startId += 1
        
    else:
      if self.peakList.name:
        expName = self.peakList.name
      else:
        expName = self.peakList.dataSource.experiment.name
      specName = self.dataEntry.askString("Enter experiment name","Ansig peak list export: Please enter a name for the peaklist for experiment '%s':" % expName,expName,self.guiParent)
      specName = specName[:12]
      
    while (self.peakFile.specNames.count(specName)):
      specName = self.dataEntry.askString("Enter experiment name","Ansig peak list export: name '%s' exists already. Please enter another name:" % specName,specName,self.guiParent)
      specName = specName[:12]
      
    self.peakFile.specNames.append(specName)
    self.specName = specName

  def setPeakFilePeakIntensity(self):
  
    self.volume = 0.0
   
    if self.peak.peakIntensities:
      
      # Only use obvious volume. Ignore otherwise
      
      for peakIntensity in self.peak.peakIntensities:
        if peakIntensity.intensityType == volume_kw:
          self.volume = peakIntensity.value

  def setPeakFilePeakExtras(self):

    self.assignResNum = [] 
    self.assignRes = []
    self.assignAtom = []
    
    self.ppm = []

  def setPeakFilePeakDimInfo(self):
    
    self.ppm.append(self.peakDim.getValue())

    #
    # Transfer assignment
    #

    peakDimContribs = self.peakDim.sortedPeakDimContribs()
    
    code3Letter = seqCode = atomName = ''
    
    if self.writeAssignments and peakDimContribs:
      
      if isinstance(peakDimContribs[0],Nmr.PeakDimContribN):
        print "   Warning: cannot handle multiple resonances for one peakDim contribution. Ignored."
        return
        
      atomSetIds = []
      atomIds = []
        
      for peakDimContrib in peakDimContribs:

        resonance = peakDimContrib.resonance
        
        if resonance.resonanceSet:
          code3Letter = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue.molResidue.chemComp.ccpCode
          code3Letter = string.capitalize(code3Letter)

        else:
          code3Letter = ""
      
        if self.resonanceToAtoms.has_key(resonance):

          if self.useOriginalResNames:
          
            if len(self.resonanceToAtoms[resonance]) > 1:
              print "   Warning: Only first original resName taken for peak %d." % self.peak.serial

            resonanceToAtom = self.resonanceToAtoms[resonance][0]

            atomIds.append([resonanceToAtom.chain.code,resonanceToAtom.seqId,resonanceToAtom.atomName,None,code3Letter])
            
          else:
          
            for resonanceToAtom in self.resonanceToAtoms[resonance]:
                
              #
              # Always keep track of single ones...
              #
              
              # Modified because cleanAtomNameLists needs atomId and atomSetId that compare
              atomId = resonanceToAtom.getAtomId()
              atomId[-1] = code3Letter
              
              if atomId and atomIds.count(atomId) == 0:
                atomIds.append(atomId)
                
                #
                # Check atomSet: use this if the same for all...
                #
 
                atomSetId = resonanceToAtom.getAtomSetId()

                if atomSetId:
                  atomSetId[-1] = code3Letter
                else:
                  atomSetId = None
                
                atomSetIds.append(atomSetId)
                
      if len(atomIds) > 1 and self.compressResonances:
        # Replace atomId with single atomSetId for duplicates
        atomIds, atomSetIds = cleanAtomNameLists(atomIds, atomSetIds)

      chain = None
      
      if atomIds:
        chain, seqId, atomName, code3Letter = atomIds[0]
        
        
        if len(atomIds) > 1:
          print "   Warning: Multiple assignment possibilities for peak %d - only first used." % self.peak.serial

        
      else:
      
        print "    Warning: no atom link for resonance %s" % getResNameText(resonance)

      #
      # Use chain and seqcode mapping
      #

      if chain:
      
        #chainCode = self.chainDict[chain][0]
        seqCode = self.getExportSeqCode(self.chainDict[chain][1],chain.findFirstResidue(seqId = seqId))

    self.assignResNum.append(str(seqCode))  
    self.assignRes.append(code3Letter) 
    self.assignAtom.append(atomName) 


  def createPeakFilePeak(self):
  
    self.peakFile.npeaks += 1
    
    self.peakFile.peaks.append(self.peaksIO.AnsigPeak(self.peakNum,0,self.ppm,self.volume,self.assignResNum,self.assignRes,self.assignAtom,self.numDim,self.specName))


  def setPeakDimOrder(self):
  
    self.peakDimOrder = {
       
       1: [0],
       2: [1,0],
       3: [1,2,0],
       4: [1,3,2,0],
       5: [1,4,3,2,0],
       6: [1,5,4,3,2,0],
       7: [1,6,5,4,3,2,0]
       
      }

  def getSpecificExperimentList(self):
  
    experimentList = self.project.currentNmrProject.sortedExperiments()
    interactionText = "Existing experiments:"
    
    return (experimentList,interactionText)
      
  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)
