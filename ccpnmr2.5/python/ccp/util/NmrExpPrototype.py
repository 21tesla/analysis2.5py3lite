"""
======================COPYRIGHT/LICENSE START==========================

NmrExpPrototype.py: Utility functions

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

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

# Transfer types that may connect long range, non-local atoms
#longRangeTransfers = ('through-space', 'NOESY', 'SpinDiff', 'DipDip', 'Jnonlocal')
# NBNB TBD obsolete, hardwire 'through-space' instead
longRangeTransfers = ('through-space',)


def nextFreeExpMeasurement(nmrExpPrototype, checkAttr):
  """ get ExpMeasurement in nmrExpPrototype that has fewest links to checkAttr,
  If there are several choose the one with the smallest serial.
  nmrExpPrototype is the relevant nmrExpPrototype
  checkAttr is the name of teh attribtue to check on.
  """
  expMeasurements = nmrExpPrototype.sortedExpMeasurements()
  if expMeasurements:
    # find first expMeasurement with fewest refExpDimRefs
    dd = {}
    for xx in expMeasurements:
      length = len(getattr(xx, checkAttr))
      if length not in dd:
        dd[length] = xx
    expMeasurement = dd[min(dd.keys())]
  else:
    expMeasurement = None
  #
  return expMeasurement
