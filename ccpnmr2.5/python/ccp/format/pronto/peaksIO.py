"""
======================COPYRIGHT/LICENSE START==========================

peaksIO.py: I/O for Pronto peak list files

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
from ccp.format.pronto.generalIO import ProntoGenericFile
from memops.universal.Io import getTopDirectory

#####################
# Class definitions #
#####################
    
class ProntoPeakFile(ProntoGenericFile):

  """
  Information on file level
  """
  
  def initialize(self):
  
    self.peaks = []
    self.assignSep = '.'
    
  def setSpectrumInfo(self,specName,ndim):
  
    self.specNames = [specName]
    self.numDims = [ndim]

    self.dimCodes = ndim * ['']
    self.specName = specName

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading Pronto peak list %s" % self.name

    fin = open(self.name, 'rU')

    # Read the atom specs on second line
    line = fin.readline()

    # Read file
    while line:

      if self.patt['hash'].search(line) or self.patt['emptyline'].search(line):
        line = fin.readline()
        continue

      cols = string.split(line)

      # Line with spectrum info (first one)
      if cols[0] == 'Spectrum:':
        specName = cols[-1]

      elif cols[0] == 'key':

        #
        # Bit dangerous but PRONTO peak list style should be fixed
        #

        n = cols.index('int') - 1

        self.setSpectrumInfo(specName,n)

        assign = n*['']
        ppm = n*[0]

      else:

        #
        # Line with peak info
        #

        peaknum = cols.pop(0)

        for dim in range(0,n):
          ppm[dim] = cols.pop(0)

        intensity = cols.pop(0)

        status = cols.pop(0)

        assignString = string.join(cols,' ')

        if assignString:
          self.getAssignments(assignString,assign)

        self.peaks.append(ProntoPeak(peaknum,assign,ppm,intensity,status))

      line = fin.readline()

  def getAssignments(self,assignString,assign):

    assignComps = string.split(assignString)

    nDim = len(assign)

    stringCount = 0
    assignCount = 0

    # TODO: what if an assignment missing in one dimension? How is this represented?  

    while (assignCount < nDim):

      (resLabel,resCode) = string.split(assignComps[stringCount],self.resLabelSep)
      resCode = resCode[:-1] # Get rid of trailing :

      atomNames = string.split(assignComps[stringCount+1],self.atomNameSep)

      for atomName in atomNames:
        assign[assignCount] = string.join((self.defaultMolCode,resCode,atomName),self.assignSep)
        assignCount += 1

      stringCount += 2

  def write(self, verbose = 0):

    if verbose == 1:
      print "Writing Pronto peak list %s" % self.name

    #fout = open(self.name,'w')

    #numDim = len(self.dimCodes)

    #
    # Write out header
    #

    #fout.write("label dataset sw sf" + self.newline)


    fout.close()

class ProntoPeak:

  def __init__(self,num,assign,ppm,intensity,status):

    self.num = returnInt(num)
    self.intensity = returnFloat(intensity)
    self.status = status

    self.assign = []
    self.ppm = []

    for i in range(0,len(assign)):
      self.ppm.append(returnFloat(ppm[i]))
      self.assign.append(assign[i])
  
###################
# Main of program #
###################
    
if __name__ == "__main__":

  files = ['../reference/pronto/cp40',
           '../reference/pronto/cp150'
          ]

  
  for inFile in files:
    
    peakFile = ProntoPeakFile(os.path.join(getTopDirectory(),inFile))

    peakFile.read(verbose = 1)

    
    for rawPeak in peakFile.peaks:
      print rawPeak.intensity
