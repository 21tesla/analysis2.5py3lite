"""
======================COPYRIGHT/LICENSE START==========================

DataSource.py: Utility functions for ccp.nmr.Nmr.DataSource

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

# Additional functions for ccp.nmr.Nmr.DataSource
#
# NB All functions must have a mandatory DataSource as the first parameter
# so they can be used as DataSource methods
from ccp.lib.nmr.Nmr.Experiment import getOnebondExpDimRefs, getAcqExpDim
from ccp.lib.nmr.Nmr.AbstractDataDim import getIsotopeCodes


def getDimCodes(dataSource):
  """ Get dimcode of form hx1, hx2, x1, x2, where the x's are directly bound to 
  the corresponding H. suffix '1' is given to the acquisition proton dimension.
  Dimensions not matching the specs are given code None
  """
  isotopeCodes=getIsotopeCodesList(dataSource)
  
  acqExpDim = getAcqExpDim(dataSource.experiment)
  dataDims = dataSource.sortedDataDims()
  dimCodes = [None]*dataSource.numDim
  for ii,dataDim in enumerate(dataDims):
    if isotopeCodes[ii] == '1H':
      if dataDim.expDim is acqExpDim:
        dimCodes[ii] = 'hx1'
        xCode = 'x1'
      else:
        dimCodes[ii] = 'hx2'
        xCode = 'x2'
      for tt in getOnebondDataDims(dataSource):
        if tt[0] is dataDim:
          xDim = tt[1]
          dimCodes[dataDims.index(xDim)] = xCode
          break
  #
  return dimCodes
  
  
def getXEasyDimCodes(dataSource):
  """ Get Xeasy-style dimCodes in dim order for data Source
  For use in FormatCOnverter
  """
  #
  # Get isotopecode info
  #
  isotopeCodes = getIsotopeCodesList(dataSource)
  dimCodes = getDimCodes(dataSource)
  reverse = bool( 'x2' in dimCodes and not 'x1' in dimCodes)
  found = set()
  for ii,code in enumerate(dimCodes):
    elementType = isotopeCodes[ii][-1]
    if code and (('2' in code) != reverse):
      dimCodes[ii] = elementType.lower()
    else:
      dimCodes[ii] = elementType
    
  #
  # Fix in homonuclear 2D case...
  #
  #firstDimCode = dimCodes[0]
  #if dimCodes.count(firstDimCode) == len(dimCodes):
  #  if firstDimCode == firstDimCode.lower():
  #    dimCodes[0] = firstDimCode.upper()
  #
  # and generally distinguish duplicates. Rasmus Fogh 15/6/12
  foundset = set()
  for ii,ss in enumerate(dimCodes):
    if ss in foundset:
      ss = ss.upper()
    if ss in foundset:
      ss = ss.lower()
    dimCodes[ii] = ss
    foundset.add(ss)
  #
  return dimCodes
  

# NB renamed from getSpectrumIsotopes getIsotopeCodesList
#def getSpectrumIsotopes(dataSource):
def getIsotopeCodesList(dataSource):
  """
  Give isotope code strings pertaining to the dimensions of a given NMR spectrum 
  
  .. describe:: Input
  
  Nmr.DataSource
  
  .. describe:: Output

  List of Strings (comma-joined Nmr.ExpDimRef.IsotopeCodes)
  """

  isotopes = []
  for dataDim in dataSource.sortedDataDims():
    isotopeCodes = list(getIsotopeCodes(dataDim))
    isotopes.append(','.join(sorted(isotopeCodes)) or None)
  
  return isotopes
  
       
def getOnebondDataDims(spectrum):
  """
  Get pairs of spectrum data dimensions that are connected by onebond transfers
  
  .. describe:: Input
  
  Nmr.DataSource
  
  .. describe:: Output

  List of 2-List of Nmr.DataDims
  """
  
  dataDims = []
  expDimRefs = getOnebondExpDimRefs(spectrum.experiment)

  for expDimRef0, expDimRef1 in expDimRefs:
    dataDim0 = spectrum.findFirstDataDim(expDim=expDimRef0.expDim)
    dataDim1 = spectrum.findFirstDataDim(expDim=expDimRef1.expDim)

    if dataDim0 and dataDim1:
      dataDims.append( [dataDim0,dataDim1] )

  return dataDims
