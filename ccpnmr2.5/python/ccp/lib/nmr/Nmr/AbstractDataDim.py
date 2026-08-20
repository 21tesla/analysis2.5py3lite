"""
======================COPYRIGHT/LICENSE START==========================

AbstractDataDim.py: Utility functions for ccp.nmr.Nmr.AbstractDataDim

Copyright (C) 2005-2013 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

A copy of this license can be found in ../../../license/LGPL.license

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

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================
"""

# Additional functions for ccp.nmr.Nmr.AbstractDataDim
#
# NB All functions must have a mandatory DataSource as the first parameter
# so they can be used as AbstractDataDim methods


   
# NB change from getDataDimIsotopes
#def getDataDimIsotopes(dataDim):
def getIsotopeCodes(dataDim):
  """
  Get the shift measurement isotopes for a spectrum data dim
  
  .. describe:: Input
  
  Nmr.AbstarctDataDim (or subtypes)
  
  .. describe:: Output

  Set of Words (Nmr.ExpDimRef.isotopeCodes)
  """

  isotopes = set()
    
  for expDimRef in dataDim.expDim.expDimRefs:
    if expDimRef.measurementType in ('Shift','shift'):
      for isotopeCode in expDimRef.isotopeCodes:
        isotopes.add(isotopeCode)
 
  return isotopes
