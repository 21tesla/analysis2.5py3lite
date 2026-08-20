
"""
======================COPYRIGHT/LICENSE START==========================

AriaXmlFormat.py: Contains functions specific to Aria XML format conversions.

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
import ccp.format.aria.xmlreader as xmlreader
import sys, string

from ccpnmr.format.general.Constants import defaultMolCode
from ccpnmr.format.general.Util import getRefExpFromOldExpType
from ccpnmr.format.converters.Tool import Tool

from memops.api import Implementation

from ccpnmr.format.general.Constants import ccpNmr_kw, isAriaInput_kw

class AriaXmlFormat:

  def __init__(self,fileName):

    contentHandler = xmlreader.XMLContentHandler()
    xmlFileReader = xmlreader.XMLReader(contentHandler)

    self.xmlFile = xmlFileReader.load(fileName)

    self.checkMandatory()
    
    self.ccpNmrConversion()
    
    #self.writeXmlFile()
    
  def writeXmlFile(self):
  
    self.ccpNmrConv.writeXmlProject()
    
  def tagAriaInput(self,object):
  
    appData = object.findFirstApplicationData(application = ccpNmr_kw, keyword = isAriaInput_kw)
    
    if not appData:
      object.addApplicationData(Implementation.AppDataBoolean(application = ccpNmr_kw, keyword = isAriaInput_kw, value = True))
    
  def ccpNmrConversion(self):
  
    ccpNmrConv = Tool()
    
    ccpNmrConv.setupProject(self.xmlFile.conversion.project.filename,
                            self.xmlFile.conversion.project.name)
    
    ccpNmrConv.setupMolSystem(self.xmlFile.conversion.molecular_system.name)
    
    self.tagAriaInput(ccpNmrConv.molSystem)
    
    #
    # Molecules
    #
    
    molecules = self.xmlFile.conversion.molecular_system.molecule
    
    if type(molecules) != type([]):
      molecules = [molecules]
      
    for molecule in molecules:
  
      ccpNmrConv.setupFormat(molecule.sequence_format)
      
      if not hasattr(molecule,'segID'):
        segID = defaultMolCode
      else:
        segID = molecule.segID
      
      ccpNmrConv.readSequence(molecule.sequence_filename,
                              molecule.name,
	          molecule.type,
	          segID,
	          molecule.first_residue_number)
    
    #
    # Spectra
    #
    
    spectra = self.xmlFile.conversion.spectrum
    
    if type(spectra) != type([]):
      spectra = [spectra]
      
    for spectrum in spectra:
  
      ccpNmrConv.setupFormat(spectrum.data_files_format)
      
      if spectrum.chemical_shifts.filename:
        shiftList = ccpNmrConv.readShifts(spectrum.chemical_shifts.filename,spectrum.name)
      else:
        shiftList = ccpNmrConv.createChemShiftsFromPeakList(spectrum.name)

      self.tagAriaInput(shiftList)

      
      peakListDimOrder = [spectrum.cross_peaks.proton1,
                          spectrum.cross_peaks.proton2,
	                         spectrum.cross_peaks.hetero1,
	                         spectrum.cross_peaks.hetero2]
                          
                        
      # TODO TODO Should the names not directly reflect the ones from the data model?!?! Better...  
      # or have alternative system for this at least
      refExperiment = getRefExpFromOldExpType(ccpNmrConv.ccpnProject,spectrum.type)
      
      if not refExperiment:
        print "  Error unknown experiment type %s... not created." % spectrum.type

      peakLists = ccpNmrConv.readPeakList(spectrum.cross_peaks.filename,
                                         refExperiment,
  	                                      spectrum.name,
  	                                      peakListDimOrder)
					 
      self.tagAriaInput(peakLists[0])
	              
      if spectrum.cross_peaks_assignments.filename:
        ccpNmrConv.readPeakListAssignments(spectrum.cross_peaks_assignments.filename)

    #
    # Link resonances to atoms
    #

    ccpNmrConv.linkResonances(self.xmlFile.conversion.project.naming_system,setSinglePossEquiv = 0, setSingleProchiral = 0, useIupacMatching = 1, minimalPrompts = 1)
    
    self.ccpNmrConv = ccpNmrConv
    
  def checkMandatory(self):  

    blocks = ['project','molecular_system','spectrum']
    
    subBlocks = {'project': ['name','filename','naming_system'],
    
                 'molecular_system': ['name',
	                 {'molecule': ['name','type','first_residue_number','sequence_filename','sequence_format']}
                                    ],
	                
                 'spectrum': ['name',
                              'type',
	          'data_files_format',
	          {'chemical_shifts': ['filename']},
	          {'cross_peaks': ['filename','proton1','proton2','hetero1','hetero2']},
	          {'cross_peaks_assignments': ['filename']}
                             ]
                }
    
    for blockName in blocks:
      
      if not hasattr(self.xmlFile.conversion,blockName):
        
        self.exitError(blockName)
      
      else:
        
        xmlItems = getattr(self.xmlFile.conversion,blockName)
        
        if type(xmlItems) == type([]):
          for xmlItem in xmlItems:
            self.checkSubBlocks(xmlItem,subBlocks[blockName],blockName)
        
        else:
          self.checkSubBlocks(xmlItems,subBlocks[blockName],blockName)
        

  def checkSubBlocks(self,xmlItem,subBlockList,parentBlockName):

    for subBlock in subBlockList:
    
      if type(subBlock) == type({}):
      
        subBlockName = subBlock.keys()[0]
      
        if not hasattr(xmlItem,subBlockName):
    
          self.exitError(parentBlockName + '.' + subBlockName)
        
        else:
          
          xmlSubItems = getattr(xmlItem,subBlockName)
          
          if type(xmlSubItems) == type([]):
            for xmlSubItem in xmlSubItems:
              self.checkSubBlocks(xmlSubItem,subBlock[subBlockName],subBlockName)
          
          else:
            self.checkSubBlocks(xmlSubItems,subBlock[subBlockName],subBlockName)
        
      else:
      
        if not hasattr(xmlItem,subBlock):
    
          self.exitError(parentBlockName + '.' + subBlock)
    
  def exitError(self,attrName):
  
    print "Invalid conversion xml file: need %s information" % (attrName)
    sys.exit()


if __name__ == '__main__':

  #AriaXmlFormat('local/testAria/ariaConvTest.xml')
  ariaXmlFormat = AriaXmlFormat('local/testAria/ariaConvTest2.xml')

  ccpnProject = ariaXmlFormat.ccpNmrConv.ccpnProject
  
  print ccpnProject
  
  #To write...
  #ariaXmlFormat.writeXmlFile()
