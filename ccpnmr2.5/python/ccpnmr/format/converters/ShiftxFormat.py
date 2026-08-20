"""
======================COPYRIGHT/LICENSE START==========================

ShiftxFormat.py: Contains functions specific to ShiftX conversions.

Copyright (C) 2009 Wim Vranken (European Bioinformatics Institute)

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

from ccpnmr.format.general.Constants import defaultMolCode

import ccp.api.nmr.Nmr as Nmr

class ShiftxFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'shiftx'
    
  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric
    
    self.getMeasurements = self.getMeasurementsGeneric
    self.createMeasurementFile = self.createMeasurementFileGeneric

  #
  # Deviations from generic import stuff
  #
    
  #
  # Functions different to default functions in DataFormat
  #
    
  def initChemShiftFileResidue(self):
    
    """
    
    Different from DataFormat because have to set self.seqCodeLabel dictionary
    
    Use seqCode of residue to write out to chemical shift file
    Use seqId for identification
    
    """
    
    self.seqId = self.residue.seqId
    self.seqCode = self.getExportSeqCode(self.chainDict[self.chain][1],self.residue)
    
    residueLabel = self.residue.molResidue.chemComp.code1Letter
    
    if not residueLabel:
      residueLabel = 'X'
    
    self.measurementFile.seqCodeLabels[self.seqCode] = residueLabel

  def setChemShiftFileValue(self):
    
    if self.atomNamesDict.has_key(self.atomName):
    
      resonanceToAtom = self.atomNamesDict[self.atomName]
      residueLabel = self.getResidueLabelOneLetter(resonanceToAtom)

      # Need some special handling for HA values.
      searchMatch = False
      if residueLabel == 'G' and self.atomName[:2] == 'HA':
        shiftxAtomName = 'HA'
        if self.atomName in ('HA2','HA3'):
          searchMatch = True
      elif self.atomName == 'H*':
        shiftxAtomName = 'H'
        searchMatch = True
      else:
        shiftxAtomName = self.atomName
      
      #
      # Export all atoms (not necessarily written!)
      #

      measurement = self.origAtomMeasurements[resonanceToAtom]

      if searchMatch:
        searchMatch = False
        for rawMeasurement in self.measurementFileValues:
          if rawMeasurement.seqCode == self.seqCode and rawMeasurement.atomName == shiftxAtomName:
            rawMeasurement.value = (rawMeasurement.value + measurement.value) / 2
            searchMatch = True
            break

      if not searchMatch:

        self.measurementFileValues.append(

               self.rawMeasurementClass(self.measurementFile,
                                        measurement.value,
                                        shiftxAtomName,
                                        self.seqCode,
                                        residueLabel,
                                        defaultMolCode))
        
      if self.atomMeasurements.has_key(resonanceToAtom):
        del(self.atomMeasurements[resonanceToAtom])
                           
  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)
