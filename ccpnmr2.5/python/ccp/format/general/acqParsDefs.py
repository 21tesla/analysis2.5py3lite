
"""
======================COPYRIGHT/LICENSE START==========================

acqParsDefs.py: General definitions for acquisition parameter reading

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
from memops.universal.Util import returnInt
from memops.universal.Util import returnFloat, returnFloats
from memops.universal.Util import getUpperPowerTwo

from ccp.format.general.formatIO import FormatFile

##########################################################
# Class definitions for general use in acqParsIO scripts #
##########################################################

class SpectrometerReference:

  def __init__(self,specFreq,baseFreq):
  
    self.specFreq = specFreq
    self.baseFreq = baseFreq

  def setSpecType(self,specType):
  
    if specType:
      self.specType = specType
      
    else:
      self.specType = 'Unknown'
        
class FidDimReference:

  def __init__(self,numPoints,numPointsValid,valuePerPoint,nucleus):
  
    self.numPoints = returnInt(numPoints)
    self.numPointsValid = returnInt(numPointsValid)
    self.valuePerPoint = returnFloat(valuePerPoint)
    self.nucleus = nucleus
    
class FreqDimReference:

  def __init__(self,numPoints,numPointsOrig,valuePerPoint,nucleus,phase = None):
  
    self.numPoints = returnInt(numPoints)
    self.numPointsOrig = returnInt(numPointsOrig)
    self.valuePerPoint = returnFloat(valuePerPoint)
    self.nucleus = nucleus
    
    if phase:
      self.phase = returnFloats(phase)
    else:
      self.phase = 2 * [0.0]

  def setDataTypes(self):
  
    self.numPoints = returnInt(self.numPoints)
    self.numPointsOrig = returnInt(self.numPointsOrig)
    self.valuePerPoint = returnFloat(self.valuePerPoint)
    self.phase = returnFloats(self.phase)


#############################
# Functions for general use #
#############################

class GenericAcqParData(FormatFile):

  def getSpecRefs(self,pars,nullValue,specTypeRef = None):

    specRefs = {}

    for refs in self.refList:
      noValue = 0
      for ref in refs:
        if not pars.has_key(ref):
          noValue = 1
        elif pars[ref].values[0] == nullValue:     # This could be a problem, but seems to work allright
          noValue = 1

      if noValue == 0:
        (specFreqCode,specFreqOffsetCode,nucleusCode) = refs

        nucleus = self.convertNucleusCode(pars[nucleusCode].values[0])

        # Check if already a frequency with this nucleus: is an error!
        # TODO: can make SpectrometerReference values lists to accomodate all values...
        if specRefs.has_key(nucleus) > 0:
          print "  ERROR: %s has two transmitters channels with nucleus %s. Second one ignored." % (self.fileDir,nucleus)
          continue

        specFreq = pars[specFreqCode].values[0]
        baseFreq = specFreq - (pars[specFreqOffsetCode].values[0] / 1000000.0) 
        #print basefreq, sfreq, pars[sfreqOffsetCode].values[0]

        # Additional check: specFreq can be ~150 MHz, specFreqOffset ~ -50000000 Hz (so baseFreq ~ 200MHz)
        # Resetting to normal values this way... 
        if int(baseFreq/10) != int(specFreq/10):       
          print "  Error: Exp %s bad referencing: baseFreq %.3f, specFreq %.3f. Reset specFreq to baseFreq" % (self.fileDir,baseFreq,specFreq)
          specFreq = baseFreq

        specRefs[nucleus] = SpectrometerReference(specFreq,baseFreq)

      else:
        #print "  Warning: no value for acquisition parameter data %s..." % str(refs)
        pass

    #
    # Set other spectrometer data
    #

    specType = None

    if specTypeRef:
      if pars.has_key(specTypeRef):
        specType = pars[specTypeRef].values[0]

    for nucleus in specRefs.keys():
      specRefs[nucleus].setSpecType(specType)

    self.specRefs = specRefs

  def getFreqDimInfo(self,valppoint,npoints):

    if npoints > 1:
      npointsFreq = getUpperPowerTwo(npoints) * 2 # Make this default to the next power of two - better for FT

    else:
      print "  WARNING: %s has only %d number of points for FID" % (self.fileDir,npoints)
      npointsFreq = 2

    valppointFreq = (1 / valppoint) / npointsFreq

    npointsorig = npointsFreq

    return npointsFreq,npointsorig,valppointFreq

  
  def convertNucleusCode(self, nucleusCode):
    """ Convert nucleus code to accepted style. e.g. '1H'
    """

    searchObj = self.patt['nucleusLetterDigit'].search(nucleusCode)

    if searchObj:
      # Switch letter to back for 'H1' style nuclei
      result = searchObj.group(2) + searchObj.group(1)
    else:
      result = nucleusCode
    #
    return result
