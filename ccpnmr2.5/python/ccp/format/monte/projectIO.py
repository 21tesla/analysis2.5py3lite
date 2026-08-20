#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

projectIO.py: I/O for Monte project .par file

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
from ccp.format.monte.generalIO import MonteGenericFile

from ccp.format.monte.chemShiftsIO import MonteChemShiftFile
#from ccp.format.monte.sequenceIO import MonteSequenceFile

#####################
# Class definitions #
#####################

class MonteProjectFile(MonteGenericFile):

  # Information on file level

  def initialize(self):
    
    self.rootFileName = None
    self.subDirPath = None
    self.headerCols = None
    
    self.chemShiftFiles = []
    #self.sequenceFile = None

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s project file %s" % (self.format,self.name)
    

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      cols = line.split()

      if self.patt['%sComment' % self.format].search(line) or self.patt['emptyline'].search(line):
        pass

      else:
      
        if line.find('Root Filename:') > -1:
        
          line = fin.readline()
          self.rootFileName = string.strip(line)

        elif line.find('Subdirectory path:') > -1:
        
          line = fin.readline()
          self.subDirPath = string.strip(line)
        
        elif line.find('shift table column definitions') > -1 and not line.find('old solution') > -1:
        
          # Depending on fixed par file output here!!
          line = fin.readline()
          numCols = len(line.split())
          self.headerCols = [None] * numCols
          
          line = fin.readline()
          line = fin.readline()

          while not self.patt['emptyline'].search(line):
          
            cols = line.split()
            atomInfo = cols.pop(0)
            
            if cols.count('1'):
              colIndex = cols.index('1')
              self.headerCols[colIndex] = atomInfo
            
            line = fin.readline()                

      line = fin.readline()

    fin.close()
