"""
======================COPYRIGHT/LICENSE START==========================

AmberFormat.py: Contains functions specific to AMBER conversions.

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

import copy

from ccpnmr.format.converters.DataFormat import DataFormat, IOkeywords

from ccp.format.general.Constants import defaultSeqInsertCode

from ccpnmr.format.general.Constants import distanceConstraintDefaultLowerLimit
from ccpnmr.format.general.Constants import defaultMolCode

#
# Add some information to IOkeywords...
#

IOkeywords = copy.deepcopy(IOkeywords)

IOkeywords['readDistanceConstraints']['coordFile'] = (None,True,'The coordinate file related to this restraint list. Required for converting atom serials to atom names.')
IOkeywords['readHBondConstraints']['coordFile'] =    (None,True,'The coordinate file related to this restraint list. Required for converting atom serials to atom names.')
IOkeywords['readDihedralConstraints']['coordFile'] = (None,True,'The coordinate file related to this restraint list. Required for converting atom serials to atom names.')

class AmberFormat(DataFormat):

  def setFormat(self):
  
    self.format = 'amber'
    self.IOkeywords = IOkeywords

  def setGenericImports(self):
    
    self.getSequence = self.getSequenceGeneric

    self.getCoordinates = self.getCoordinatesGeneric

  #
  # Deviations from generic import stuff
  #
  
  def getConstraints(self):
    
    try:
    
      self.constraintFile = self.ConstraintFileClass(self.fileName)
      self.constraintFile.read(coordFile=self.coordFile,verbose=self.verbose)

      if self.verbose:
        print "Reading %s constraint list from %s file %s" % (self.constraintType,self.formatLabel,self.fileName)

    except:

      errorMessage = traceback.format_exception_only(sys.exc_type,sys.exc_value)[-1]
      self.messageReporter.showWarning("Warning"," Cannot read %s constraints for %s...:\n%s" % (self.constraintApiCode,self.formatLabel,errorMessage),self.guiParent)
      self.constraintFile = None
      
      raise
      
      return traceback.format_exception(sys.exc_type,sys.exc_value,sys.exc_info()[2]) 
  
  #
  # Functions different to default functions in DataFormat
  #

  def duplicateChain(self):
  
    #
    # All chain info should be contained in coordinate file...
    #
    
    return False

  def getPresetChainMapping(self,chainList):
  
    return self.getSingleChainFormatPresetChainMapping(chainList)    
