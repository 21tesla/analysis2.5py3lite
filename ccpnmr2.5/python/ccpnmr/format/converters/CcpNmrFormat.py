"""
======================COPYRIGHT/LICENSE START==========================

CcpNmrFormat.py: Contains functions specific to CcpNmr Analysis tab file conversions.

Copyright (C) 2010 Wim Vranken (European Bioinformatics Institute)

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

import copy, re

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccpnmr.format.general.Constants import defaultMolCode
from ccpnmr.format.general.Constants import volume_kw, height_kw, ccpNmr_kw
from ccpnmr.format.general.Util import getResName, getResNameText

#
# Add some information to IOkeywords...
#

#IOkeywords = copy.deepcopy(IOkeywords)

class CcpNmrFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'ccpNmr'
    self.IOkeywords = IOkeywords
    
    self.assignPatt = re.compile("([A-Za-z]*)(\d+)([A-Z][a-z]*)([^a-z]+)")

  def setGenericImports(self):
    
    #self.getSequence = self.getSequenceGeneric
        
    self.getPeaks = self.getPeaksGeneric
    
    self.getMeasurements = self.getMeasurementsGeneric

  #
  # Deviations from generic import stuff
  #
    
  #
  # Functions different to default functions in DataFormat
  #

  def setPeakIntensity(self):
  
    # PeakIntensity attributes
    # TODO: Should probably be using some CcpNmr Analysis value here?
    
    if self.rawPeak.volume != 0:
      peakInt = self.peak.newPeakIntensity(value = self.rawPeak.volume, method = self.methods[self.format]['Volume'])
      peakInt.intensityType = volume_kw

    if self.rawPeak.intensity != 0:
      peakInt = self.peak.newPeakIntensity(value = self.rawPeak.intensity, method = self.methods[self.format]['Intensity'])
      peakInt.intensityType = height_kw

  def setPeakExtras(self):
  
    if self.peakFile.hasAssignItems:
      
      #
      # Find maximum number (could in principle be different?)
      #
      
      self.peakContribs.append(self.peak.newPeakContrib())

  def setPeakDim(self):
  
    # TODO: rawPeak.shape[i]
    dataDimRef = self.dataDimRefs[self.rawPeakDimIndex]

    self.peakDim = self.peak.findFirstPeakDim(dim = dataDimRef.dataDim.dim)

    self.peakDim.dataDimRef = dataDimRef

    self.peakDim.value = self.rawPeak.ppm[self.rawPeakDimIndex]
    
    if self.rawPeak.width and self.rawPeak.width[self.rawPeakDimIndex]:
      self.peakDim.decayRate = self.rawPeak.width[self.rawPeakDimIndex]
      
  def getPeakResNames(self):
  
    self.resNames = []
    self.resLabels = []
    
    origAssignName = self.rawPeak.assign[self.rawPeakDimIndex] 

    if origAssignName not in [None,'?','null','']:
    
      # Spin systems, not handled although in principle possible.
      if origAssignName.count("{"):
        pass
      else:
        
        # TODO Problems: No chain code, what if get for example 1CC5 for DNA or RNA? No label separator.
        assignSearch = self.assignPatt.search(origAssignName)
        
        if assignSearch:
          if assignSearch.group(1):
            chainCode = assignSearch.group(1)
          else:
            chainCode = defaultMolCode
        
          self.resNames.append(getResName(chainCode,assignSearch.group(2),assignSearch.group(4)))
          self.resLabels.append(assignSearch.group(3))
          
        else:
          print "   No decomposition possible for %s assignment" 

  def getPresetChainMapping(self,chainList):
    # TODO: should be able to handle multiple!
    return self.getSingleChainFormatPresetChainMapping(chainList)
