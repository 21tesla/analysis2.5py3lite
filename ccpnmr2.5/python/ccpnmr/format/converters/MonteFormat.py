#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

MonteFormat.py: Contains functions specific to Monte conversions.

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

import os, traceback, sys

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccpnmr.format.general.Util import getResName

from ccpnmr.format.general.Constants import height_kw

import ccp.api.nmr.Nmr as Nmr

class MonteFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'monte'

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    self.createSequenceFile = self.createSequenceFileGeneric
    
    self.getMeasurements = self.getMeasurementsGeneric

  #
  # Deviations from generic import stuff
  #

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
    
    
  #
  # Functions different to default functions in DataFormat
  #

  def createSequence(self):
  
    self.sequenceFile.sequences.append(self.sequenceIO.MonteSequence(molName = self.chain.molecule.name))
    self.sequence = self.sequenceFile.sequences[-1]
    
  def setSequenceFileElements(self):
  
    self.sequence.elements.append(self.sequenceIO.MonteSequenceElement(self.seqCode,self.residue.molResidue.chemComp.ccpCode))

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)
