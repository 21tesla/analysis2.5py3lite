"""
======================COPYRIGHT/LICENSE START==========================

peaksIO.py: I/O for Autoassign peak list files

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
from ccp.format.autoAssign.generalIO import AutoAssignGenericFile
from memops.universal.Io import getTopDirectory

#####################
# Class definitions #
#####################
    
class AutoAssignPeakFile(AutoAssignGenericFile):

  """
  Information on file level
  """
  
  def initialize(self):
  
    self.peaks = []
    
  def setSpectrumInfo(self,specName,ndim):
  
    self.specNames = [specName]
    self.numDims = [ndim]

    self.dimCodes = ndim * ['']
    self.specName = specName
    
  def setProjectSpectrumInfo(self,specName,rootSpec,intra,sequential,throughSpace,phase,dimInfo):
    
    # These are being reset...
    self.specNames = [specName]
    self.specName = specName
    
    
    self.rootSpecName = rootSpec
    
    self.intra = intra
    self.sequential = sequential
    self.throughSpace = throughSpace
    
    self.phase = phase
    
    self.dimInfo = dimInfo

  def read(self,verbose = False):

    if verbose:
      print "Reading AutoAssign peak list %s" % self.name

    fin = open(self.name, 'rU')
    initialized = 0

    # Read the atom specs on second line
    line = fin.readline()

    # Read file
    while line:

      if self.patt['hash'].search(line) or self.patt['emptyline'].search(line):
        line = fin.readline()
        continue

      cols = string.split(line)

      #
      # Use first line to set up spectrum (should use comment line?)
      #
      
      if not initialized:
      
        initialized = 1

        #
        # Bit dangerous but AutoAssign peak list style should be fixed
        #

        specName = cols[-1]
        n = len(cols) - 3

        self.setSpectrumInfo(specName,n)
        
      if cols[0] == '*':
        break

      #
      # Line with peak info, now with validation
      #
      
      if len(cols[1:-2]) == self.numDims[0]:
        peak = AutoAssignPeak(cols[0],cols[1:-2],cols[-2],verbose=False)
        
        # Another parsing check...
        if not None in peak.ppm + [peak.intensity,peak.num]:
          self.peaks.append(peak)

      line = fin.readline()
    
    #
    # Quick check
    #
    
    fileReadOk = True
    if not self.peaks:
      fileReadOk = False
      
    return fileReadOk

  def write(self, verbose = 0):

    if verbose == 1:
      print "Writing AutoAssign peak list %s" % self.name

    fout = open(self.name,'w')

    numDim = len(self.dimCodes)
    
    #
    # Set up column lengths... not sure if this is fully correct
    #
    
    dimColLen = []
    for n in range(numDim - 1):
      dimColLen.append((9,3,3,2))
    dimColLen.append((7,2,1,2))        

    #
    # Write out header
    #

    fout.write("#       Index  ")
    for n in range(numDim):
      fout.write("%s" % (dimColLen[n][2] * ' ' + str(n+1) + 'Dim' + dimColLen[n][3] * ' '))
    fout.write("  Intensity Workbook" + self.newline)

    #
    # Write out lines
    #
    
    for peak in self.peaks:
      
      fout.write("      %5d    " % peak.num)
      
      for n in range(numDim):
        fstring = "%%%d.%df" % (dimColLen[n][0],dimColLen[n][1])
        fout.write(fstring % peak.ppm[n])
        
      fout.write("%+11.2e %-26s" % (peak.intensity,self.specName) + self.newline)    

    fout.write("*" + self.newline)
    fout.close()

class AutoAssignPeak:

  def __init__(self,num,ppm,intensity,verbose=False):

    self.num = returnInt(num,default=None,verbose=verbose)
    self.intensity = returnFloat(intensity,default=None,verbose=verbose)

    self.assign = []
    self.ppm = []

    for i in range(0,len(ppm)):
      self.ppm.append(returnFloat(ppm[i],default=None,verbose=verbose))
      self.assign.append(None)
  
###################
# Main of program #
###################
    
if __name__ == "__main__":

  files = ['../reference/autoassign/hncacb.pks',
           '../reference/autoassign/hnca.pks',
           '../reference/autoassign/hncocacb.pks',
           '../reference/autoassign/hncoca.pks',
           '../reference/autoassign/hnco.pks',
           '../reference/autoassign/hsqc.pks',
           ]

  
  for inFile in files:
    
    peakFile = AutoAssignPeakFile(os.path.join(getTopDirectory(),inFile))

    peakFile.read(verbose = 1)

    for rawPeak in peakFile.peaks:
      print rawPeak.intensity

    (path,fileName) = os.path.split(peakFile.name)
    peakFile.name = os.path.join('local',fileName)
    peakFile.write()
