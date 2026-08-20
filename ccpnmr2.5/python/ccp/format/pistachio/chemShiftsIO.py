"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for PISTACHIO chemical shift files

Copyright (C) 2005-2009 Wim Vranken (European Bioinformatics Institute)

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

# Import general functions
from memops.universal.Util import returnFloats, returnFloat
from memops.universal.Util import returnInt
from ccp.format.pistachio.generalIO import PistachioGenericFile

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class PistachioChemShiftFile(PistachioGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.chemShifts = []
    
    self.columnInfo = [('serial',None),
                      ('seqCode',returnInt),
                      ('resLabel',None),
                      ('atomName',None),
                      ('atomType',None),
                      ('value',returnFloat),
                      ('valueError',returnFloat),
                      ('ambCode',returnInt)]
    
  def read(self,isOutputFile = False, headerCols = None, verbose = 0):
  
    """
    Is an NMR-STAR like format, but without saveframe info - parsing this
    as-is.
    """

    if verbose == 1:
      print "Reading %s chemical shift list %s" % (self.format,self.name)

    lineErrors = []
    validLines = 0

    fin = open(self.name, 'rU')

    line = fin.readline()
    
    colLen = len(self.columnInfo)

    while line:

      if self.patt['%sComment' % self.format].search(line) or self.patt['emptyline'].search(line):

        line = fin.readline()
        continue
              
      #
      # Get the info... 
      #
      
      cols = line.split()
      
      if len(cols) == colLen:
      
        validLines += 1
        
        for i in range(0,colLen):
        
          (varName,convertFunc) = self.columnInfo[i]
          
          tvalue = cols[i]
          if convertFunc:
            tvalue = convertFunc(tvalue)
          
          # Globals can be written to, locals NOT!
          globals()[varName] = tvalue
        
        self.chemShifts.append(PistachioChemShift(value,valueError,atomName,atomType,seqCode,resLabel,self.defaultMolCode,ambCode = ambCode))

      else:
        lineErrors.append(line)

      line = fin.readline()

    fin.close()
      
    #
    # Check
    #
    
    if len(lineErrors) > min(5,validLines * 0.5):
      self.chemShifts = []
      print "  Bad %s format lines:%s" % (self.format,self.newline)
      for lineError in lineErrors:
        print lineError

  def write(self,verbose = 0):
    
    print "NOT WORKING YET"
    return

    if verbose == 1:
      print "Writing %s chemical shift list %s" % (self.format,self.name)


    fout = open(self.name,'w')

    #
    # Write out chem shifts
    #

    for chemShift in self.chemShifts:

      fout.write("%3d%s%-5s %8.3f %d" % (chemShift.seqCode,
                    	 self.seqNameSep,
                    	 chemShift.atomName,
                    	 chemShift.value,
                    	 chemShift.unknownCode))
      fout.write(self.newline)

class PistachioChemShift:

  def __init__(self,value,valueError,atomName,atomType,seqCode,resLabel,molCode,ambCode = 0):
  
    self.value = value
    self.valueError = valueError
    self.atomName = atomName
    self.atomType = atomType
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.molCode = molCode
    self.resLabel = resLabel
    
    self.ambCode = ambCode
