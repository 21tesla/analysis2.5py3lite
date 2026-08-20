"""
======================COPYRIGHT/LICENSE START==========================

Experiment.py: Utility functions for ccp.nmr.Nmr.Experiment

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

# Additional functions for ccp.nmr.Nmr.Experiment
#
# NB All functions must have a mandatory Experiment as the first parameter
# so they can be used as Experiment methods



def getAcqExpDim(experiment, ignorePreset=False):
  """
  ExpDim that corresponds to acquisition dimension. NB uses heuristics

  .. describe:: Input
  
   Nmr.Experiment
  
  .. describe:: Output
  
  Nmr.ExpDim
  """
  
  ll = experiment.findAllExpDims(isAcquisition=True)
  if len(ll) == 1 and not ignorePreset:
    # acquisition dimension set - return it
    result = ll.pop()
  
  else:
    # no reliable acquisition dimension set
    result = None
    
    dataSources = experiment.sortedDataSources()
    if dataSources:
      dataSource = dataSources[0]
      for ds in dataSources[1:]:
        # more than one data source. Pick one of the largest.
        if ds.numDim > dataSource.numDim:
          dataSource = ds
      
      # Take dimension with most points
      useDim = None
      currentVal = -1
      for dd in dataSource.sortedDataDims():
        if hasattr(dd, 'numPointsOrig'):
          val = dd.numPointsOrig
        else:
          val = dd.numPoints
        if val > currentVal:
          currentVal = val
          useDim = dd
      
      if useDim is not None:
        result = useDim.expDim
  
    if result is None:
      # no joy so far - just take first ExpDim
      ll = experiment.sortedExpDims()
      if ll:
        result = ll[0]
      
  #
  return result
  
  
def getOnebondExpDimRefs(experiment):
  """
  Get pairs of experiment dimensions that are connected by onebond transfers
  
  .. describe:: Input
  
  Nmr.Experiment
  
  .. describe:: Output

  List of 2-List of Nmr.ExpDimRefs
  """

  expDimRefs   = []
  expTransfers = []
  
  for expTransfer in experiment.sortedExpTransfers():
    if expTransfer.transferType in ('onebond',):
      expTransfers.append(expTransfer)
  
  for expTransfer in expTransfers:
    expDimRefs.append(expTransfer.sortedExpDimRefs())
  
  return expDimRefs
