"""
======================COPYRIGHT/LICENSE START==========================

peaksIO.py: I/O for CcpNmr Analysis tab delimited peak list file

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

import os

from memops.universal.Util import returnFloat, returnInt
from ccp.format.ccpNmr.generalIO import CcpNmrGenericFile

#####################
# Class definitions #
#####################

            
class CcpNmrPeakFile(CcpNmrGenericFile):
  
  """
  Information on file level
  """
  
  def initialize(self,assignTagSep = ' '):
  
    self.peaks = []
    
    self.hasAssignItems = True
    
    # Dictionary key is column value (after space split!).
    # Values are:
    #  - Attribute name
    #  - Extra header columns besides this one (e.g. is 2 for 'Assign F1')
    #  - Modifier function
    #
      
    self.columnInfo = {
    
      'Number':   ('rowNum',    returnInt),
      '#':        ('num',       returnInt),
      'Position': ('ppm',       returnFloat),
      'Shift':    ('ppm',       returnFloat),
      'Assign':   ('assign',    None),
      'Assign.':  ('assign',    None),
      'Height':   ('intensity', returnFloat),
      'Volume':   ('volume',    returnFloat),
      'Line':     ('width',     returnFloat),
      'Merit':    ('status',    returnFloat),
      'Details':  ('details',   None),
      'Fit':      ('heightFit', None),
      'Vol.':     ('volumeFit', None)
      
    }

  def setSpectrumInfo(self,specName,ndim):
  
    self.specNames = [specName]
    self.numDims = [ndim]

    self.dimCodes = ndim * ['']

  def read(self,verbose = False):

    if verbose:
      print "Reading %s peak list %s" % (self.format,self.name)
      
    #
    # Open and read the file
    #

    fin = open(self.name, 'rU')
    lines = fin.readlines()
    fin.close()
    
    #
    # Get header information, determine number of dimensions
    #
    
    self.headerInfo = []
    colIndex = 0
    cols = lines[0].split("\t")
    
    numDim = 0
    
    while(True):
    
      headerCol = cols[colIndex]
      
      firstColItem = headerCol.split()[0]
      
      if self.columnInfo.has_key(firstColItem):
        (attributeName,convertFunc) = self.columnInfo[firstColItem]
        
        # Track number of dimensions
        if firstColItem in ('Position','Shift'):
          numDim += 1
      
      else:
        print "  Error: unrecognized column heading %s - will be ignored" % (headerCol)
        attributeName = convertFunc = None
      
      self.headerInfo.append((attributeName,convertFunc))
      
      colIndex += 1
      
      # Exit when done
      if colIndex > len(cols) - 1:
        break
    
    #
    # Set spectrum info
    #
    
    (path,specName) = os.path.split(self.name)
    self.setSpectrumInfo(specName,numDim)
    
    #
    # Read peak information
    #
    
    headerLen = len(self.headerInfo)
    
    for line in lines[1:]:
       
      cols = line.split("\t")
      
      if len(cols) == headerLen:
        
        peak = CcpNmrPeak(self,cols)

        self.peaks.append(peak)

      else:
      
        print "  Error: Could not read peak line, has %d columns should be %d according to header." % (len(cols),headerLen)

class CcpNmrPeak:

  def __init__(self,parent,dataCols):
  
    self.parent = parent
    
    # Lists have to be set here so values can be appended to them
    self.ppm = []
    self.assign = []
    self.width = []
    
    for i in range(len(dataCols)):
      (attrName,convertFunc) = self.parent.headerInfo[i]
      
      if not attrName:
        continue
      
      value = dataCols[i].strip()
      if value == 'None':
        value = None
      
      if convertFunc and value:
        value = convertFunc(dataCols[i])
        
      if hasattr(self,attrName) and type(getattr(self,attrName)) == type([]):
        getattr(self,attrName).append(value)
      else:
        setattr(self,attrName,value)
