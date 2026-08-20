"""
======================COPYRIGHT/LICENSE START==========================

peaksIO.py: I/O for Felix peak list files

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

import os, string

# Import general functions
from memops.universal.Util import returnFloat, returnInt
from ccp.format.felix.generalIO import FelixGenericFile
from memops.universal.Io import getTopDirectory

#####################
# Class definitions #
#####################
    
class FelixPeakFile(FelixGenericFile):

  """
  Information on file level
  """
  
  def initialize(self):
  
    self.peaks = []
    
  def setSpectrumInfo(self,specName,ndim):
  
    self.specNames = [specName]
    self.numDims = [ndim]

    self.specName = specName

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading Felix peak list %s" % self.name

    fin = open(self.name, 'rU')
    
    #
    # Use file name for spectrum name... and remove any suffix
    #
    
    specName = os.path.basename(self.name)
    specName = string.join(string.split(specName,'.')[:-1])
    
    #
    # Format is not necessarily line-separated...
    # do all at once after reading xpk line
    #
    
    allColumns = string.split(fin.read())
    
    #
    # First number is dimensions (!?!)
    #
    
    numDim = 0
    i = 0
    
    while not (numDim):

      if self.patt['onlydigit'].search(allColumns[i]):

        numDim = returnInt(allColumns[i])
        self.setSpectrumInfo(specName,numDim)
        break
       
      else:
      
        i += 1        
    
    #
    # Start reading
    #
    
    allColumns = allColumns[i + 1:]
    
    while(allColumns):
    
      peaknum = allColumns.pop(0)
      
      self.peaks.append(FelixPeak(peaknum))
      
      peak = self.peaks[-1]
            
      for i in range(0,numDim):
      
        # Keep generic - need to handle this at a later stage...
        peak.column1.append(allColumns.pop(0))
        peak.column2.append(allColumns.pop(0))
        peak.column3.append(allColumns.pop(0))
        peak.column4.append(allColumns.pop(0))
        
      peak.unknown3 = returnFloat(allColumns.pop(0))
 

  def write(self, verbose = 0):

    if verbose == 1:
      print "Writing Felix peak list %s" % self.name

    #fout = open(self.name,'w')

    #numDim = len(self.dimCodes)

    #
    # Write out header
    #

    #fout.write("label dataset sw sf" + self.newline)


    #fout.close()

class FelixPeak:

  def __init__(self,num):

    self.num = returnInt(num)
      
    self.column1 = []
    self.column2 = []
    self.column3 = []
    self.column4 = []
        
    self.unknown3 = None

###################
# Main of program #
###################
    
if __name__ == "__main__":

  files = ['../reference/felix/hsqc.txt',
           '../reference/felix/ny150_h2o.txt',
           '../reference/felix/jcfiles/hhet_temp.txt',
           '../reference/felix/jcfiles/yacl_hnca.txt'
          ]

  
  for inFile in files[-2:]:
    
    peakFile = FelixPeakFile(os.path.join(getTopDirectory(),inFile))

    peakFile.read(verbose = 1)

    
    for rawPeak in peakFile.peaks:
      print rawPeak.column1, rawPeak.column2, rawPeak.column3, rawPeak.column4
