
"""
======================COPYRIGHT/LICENSE START==========================

generalIO.py: General I/O information for mmCIF files

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

from ccp.format.general.formatIO import FormatFile
from ccp.format.general.Constants import defaultMolCode

# Sans parser code
from ccp.format.mmCif.sans.lexer import STARLexer
from ccp.format.mmCif.sans.CifParser import parser

# Custom data dictionary for reading mmCIF
from ccp.format.mmCif.Util import DataDictionary_mmCIF

# TODO USE THIS -> look at PDB 
from ccp.format.pdb.cifCodeRedirect import redirectDict

#####################
# Class definitions #
#####################

class MMCIFGenericFile(FormatFile):

  def setGeneric(self):
    
    self.format = 'mmCif'
    self.defaultMolCode = defaultMolCode
    
    self.version = None

class MMCIFFile(MMCIFGenericFile):

  """
  
  mmCIF file handling. Uses code from Dimitri Maziuk (BMRB) for parsing files, the information is then
  extracted here to be forwarded to MMCIFFormat.py code (in ccpnmr/format/converters/)
  
  """
  
  def readGeneric(self, verbose=False):
  
    if verbose:
      print "  Reading %s file %s..." % (self.format,self.name)

    fin = open(self.name, 'rU')
    lexer = STARLexer( fin)
    self.mmCif = DataDictionary_mmCIF()
    mmCifParser = parser( lexer, self.mmCif, self.mmCif)
    mmCifParser.parse()
    fin.close()
    
    #
    # Set some general information
    #
    
    self.code = self.mmCif.getPdbCode()

    #mmCif.getSequenceInfo()
    #mmCif.getBondInfo()
    #mmCif.getCoordinateInfo()
  
  def convertResName(self,resName):
  
    resName = resName.strip()
    
    if redirectDict.has_key(resName):
      resName = redirectDict[resName]
      
    return resName

if __name__ == "__main__" :
  
  for pdbCode in ('1ieh','146d'):
  
    fin = open( "/Users/wim/reference/mmCif/%s.cif" % pdbCode)
    lexer = STARLexer( fin)
    mmCif = DataDictionary_mmCIF()
    mmCifParser = parser( lexer, mmCif, mmCif)
    mmCifParser.parse()

    fin.close()

    mmCif.getSequenceInfo()
    mmCif.getPdbCode()
    mmCif.getBondInfo()
    mmCif.getCoordinateInfo()
