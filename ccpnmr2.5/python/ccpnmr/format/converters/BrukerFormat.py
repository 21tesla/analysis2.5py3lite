

"""
======================COPYRIGHT/LICENSE START==========================

BrukerFormat.py: Contains functions specific to Bruker conversions.

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
import string, copy

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccpnmr.format.general.Constants import tagSep

from ccp.general.Util import setUniqueAppData

from ccp.format.bruker.acqParsHelpReader import BrukerParHelp

import memops.api.Implementation as Implementation

#
# Bruker specific functions (should these be here?!?!)
#

IOkeywords = copy.deepcopy(IOkeywords)
IOkeywords['readAcqPars']['switchRefDimInfo'] = (True,False,'Switch acqu2s/3s info when reading\n(so it ends up as information for a different dimension).')

format = 'bruker'

def getDMXtypeSpectrometer(spectrometer):

  if spectrometer.model in ['DMX','DRX','Avance']:
      return spectrometer.model
      
  return None

def getBrukerDQDdata(rawData):

  dqdData = {}
  
  if rawData.findFirstApplicationData(application = format, keyword = 'acqMode', value = 'DQD'):
    
    for applData in rawData.findAllApplicationData(application = format):
    
      if applData.keyword[:8] == ('acqMode' + tagSep):
      
        (void,par) = string.split(applData.keyword,tagSep)
        dqdData[par] = applData.value
      
  return dqdData

#
# Real class
#

class BrukerFormat(DataFormat):

  def setFormat(self):
  
    self.format = format
    self.IOkeywords = IOkeywords
    
    #
    # Only done once...
    #
    
    self.acqParsHelp = BrukerParHelp('help')

  def getAcqPars(self):
    
    self.acqParsData = self.AcqParsDataClass(self.fileName,self.numDim,self.acqParsHelp.tags)
    self.acqParsData.get(verbose = self.verbose)

  def getProcPars(self):
  
    self.procParsFile = self.processingParsIO.BrukerProcessingParData(self.fileName,self.acqParsHelp.tags)
    self.procParsFile.get(verbose = self.verbose)

    if self.verbose == 1:
      print "Reading processing pars from %s file %s" % (self.formatLabel,self.fileName)
  
  #
  #
  # Here format specific stuff
  #
  #
  
  def setFidDataSourceSpecific(self):
  
    #
    # TODO: need to check if file there - selection if not
    #
  
    if self.numDim >=2:

      expFidName = 'ser'
      
    else:

      expFidName = 'fid'         # Bruker 1D's are called fid!
    
    
    self.setFidDataSourceLocation(expFidName)

    #
    # Set some extra data if DQD spectra...
    #
    
    if self.acqParsData.acqMode:

      setUniqueAppData('AppDataString',self.fidDataSource,self.format,'acqMode',self.acqParsData.acqMode)
        
      for par in self.acqParsData.acqModeValues:
        
        setUniqueAppData('AppDataString',self.fidDataSource,self.format,'acqMode' + tagSep + par,str(self.acqParsData.acqModeValues[par]))

  def setAcqParsSpecific(self):
  
    #
    # Switch dimensions to match the nucleus setup if necessary...
    # currently only works for 3D spectra.
    #
    
    if self.numDim == 3:

      nucleusDim2 = self.acqParsData.fidDimRefs[1].nucleus
      nucleusDim3 = self.acqParsData.fidDimRefs[2].nucleus

      if nucleusDim3:
      
        expDim2 = self.experiment.findFirstExpDim(dim = 2)
        expDim3 = self.experiment.findFirstExpDim(dim = 3)

        if self.matchNucleus(expDim2,nucleusDim3) and self.matchNucleus(expDim3,nucleusDim2):
          
          print "  Warning: switched Bruker data acqu2s/3s info for experiment '%s'!" % (self.experiment.name)
          self.switchRefDims()
            
        elif self.matchNucleus(expDim2,nucleusDim2) and self.matchNucleus(expDim3,nucleusDim3):
          
          if self.verbose:
            print "Bruker dimension order fits acqu<n>s files"
            
            
        else:
          print "  Warning: %s experiment '%s' might have dimension order problems... isotopes do not fit acqu<n>s files" % (self.formatLabel,self.experiment.name)
          
          if self.switchRefDimInfo:
            self.switchRefDims()
            print "  Warning: default switch for Bruker data acqu2s/3s info for experiment '%s'!" % (self.experiment.name)
  
  def switchRefDims(self):
  
    (tFidDimRef,tFreqDimRef) = (self.acqParsData.fidDimRefs[1],self.acqParsData.freqDimRefs[1])

    self.acqParsData.fidDimRefs[1] = self.acqParsData.fidDimRefs[2]
    self.acqParsData.freqDimRefs[1] = self.acqParsData.freqDimRefs[2]

    self.acqParsData.fidDimRefs[2] = tFidDimRef
    self.acqParsData.freqDimRefs[2] = tFreqDimRef
    
  def matchNucleus(self,expDim,matchIsotope):

    hits = 0

    for expDimRef in expDim.getExpDimRefs():
      
      if len(expDimRef.isotopeCodes) > 1:
      
        print "  Multiple isotopes for same expDimRef. Picking first one."

      isotope = expDimRef.isotopeCodes[0]

      if string.upper(isotope) == string.upper(matchIsotope):
        hits += 1

    return hits
