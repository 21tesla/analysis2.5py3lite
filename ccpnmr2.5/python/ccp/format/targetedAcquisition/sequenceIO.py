#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

sequenceIO.py: I/O for TargetedAcquisition sequence from assignment file

Copyright (C) 2011 Maxim Mayzel (Swedish NMR Centrum)

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

- contact Maxim Mayzel (maxim.mayzel@nmr.gu.se)
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

import os,string

# Import general functions
from memops.universal.Util import returnInt
from ccp.format.targetedAcquisition.generalIO import TargetedAcquisitionGenericFile
from memops.universal.Io import getTopDirectory
from ccp.format.general.formatIO import Sequence, SequenceElement
from ccp.general.Constants import code1LetterToCcpCodeDict

#####################
# Class definitions #
#####################
      
class TargetedAcquisitionSequenceFile(TargetedAcquisitionGenericFile):
  """
  Information on file level
  """
  def initialize(self):
    
    self.sequences = []
    self.sequence=''

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading TargetedAcquisition sequence file %s" % self.name
    print "Reading TargetedAcquisition sequence file %s" % self.name

    molName = None

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      sequence=self.patt[self.format + 'sequence'].search(line)
      if sequence:
        sequence=sequence.group(1)
        self.sequence=sequence
        break

      line = fin.readline()

    fin.close()

    if self.sequence:
      seqCode=1
      self.sequences.append(TargetedAcquisitionSequence())
      for code1Letter in sequence:        
        self.sequences[-1].elements.append(TargetedAcquisitionSequenceElement(self,seqCode,code1Letter))
        seqCode += 1
        

TargetedAcquisitionSequence = Sequence

class TargetedAcquisitionSequenceElement(SequenceElement):

  def setResidueCode(self,*args):
    seqCode=args[0]
    code1Letter = args[1]
    self.seqCode=seqCode
    self.code1Letter = string.upper(code1Letter)        
    self.code3Letter=code1LetterToCcpCodeDict['protein'][self.code1Letter]