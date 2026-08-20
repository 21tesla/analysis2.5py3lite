"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for TargetedAcquisition chemical shift files

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

import string

# Import general functions
from memops.universal.Util import returnFloat
from memops.universal.Util import returnInt
from ccp.format.targetedAcquisition.generalIO import TargetedAcquisitionGenericFile

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class TargetedAcquisitionChemShiftFile(TargetedAcquisitionGenericFile):
  """
  Information on file level
  """
  def initialize(self):

    self.chemShifts = []
    self.seqCodes = []

  def read(self, verbose = 0):

    if verbose == 1:
      print "Reading %s chemical shift list %s" % (self.format, self.name)

    fin = open(self.name, 'rU')

    atomCols = []

    line = fin.readline()
    while line:

      cols = line.split()

      if cols:

        seqCode1Or1LetterAndCodeSearch = self.patt[self.format + 'seqCode1Or1LetterAndCode'].search(cols[0])

        if seqCode1Or1LetterAndCodeSearch:

          figOfMeritCode = cols[1]
          figOfMerit = returnFloat(cols[2])
          if figOfMerit > 0.1:

            resCode = seqCode1Or1LetterAndCodeSearch.group(1)
            seqCode = returnInt(seqCode1Or1LetterAndCodeSearch.group(2))

            if not self.seqCodes or (seqCode, resCode) not in self.seqCodes:
              print 'New residue %i,%s' % (seqCode, resCode)
              self.seqCodes.append((seqCode, resCode))

            for colNum in range(3, len(cols)):

              value = returnFloat(cols[colNum])
              if value == 0.0:
                continue

              curSeqCode = seqCode
              curResCode = resCode
              iCode = 'i'
              atomName = atomCols[colNum - 3]
              if atomName.count('-1') > 0:
                # Check if it's already there...
                atomName = atomName[:-2]
                for i in range(len(self.chemShifts) - 1, -1, -1):
                  chemShift = self.chemShifts[i]
                  if chemShift.seqCode == (seqCode - 1) and chemShift.atomName == atomName:
                    curSeqCode = chemShift.seqCode
                    curResCode = chemShift.resLabel
                    iCode = 'i+1'
                    #if abs(chemShift.value - value) > 0.1:
                    #  print "Resetting %i %s %5.2f\twith\t%i %5.2f" % \
                    #  (curSeqCode, atomName, chemShift.value, seqCode, value)
                    #chemShift.value = (chemShift.value + value) / 2
                    chemShift.allValues['i+1'] = value
                    chemShift.valueError = abs(chemShift.value - chemShift.allValues['i+1'])
                    atomName = None
                    break

                if not atomName:
                  continue
                else:
                  # If not there, add chemShift for previous residue 
                  curSeqCode -= 1
                  curResCode = None
                  for t in self.seqCodes:
                    if t[0] == curSeqCode:
                      curResCode = t[1]
                      iCode = 'i+1'
                      break
                  if not curResCode:
                    continue

              #print "set", value, atomName, curSeqCode, figOfMerit, curResCode, self.defaultMolCode, iCode
              self.chemShifts.append(TargetedAcquisitionChemShift(value, atomName, curSeqCode, figOfMerit, curResCode, self.defaultMolCode, iCode = iCode))

        elif cols[0] == '#Res':
          for atomName in cols[3:]:
            useAtomName = atomName
            if atomName == 'N15':
              useAtomName = 'N'
            elif atomName == 'HN':
              useAtomName = 'H'
            elif 'CO' in atomName:
              useAtomName = atomName.replace('CO', 'C')
            atomCols.append(useAtomName)

        elif cols[0] == '#Seq':
          seq = cols[1]
          for i, s in enumerate(seq):
            self.seqCodes.append((i + 1, s.upper()))

      line = fin.readline()

    fin.close()

  def write(self, verbose = 0):

    print "Not relevant"

class TargetedAcquisitionChemShift:

  def __init__(self, value, atomName, seqCode, figOfMerit, resLabel, defaultMolCode, iCode = 'i'):

    self.value = returnFloat(value)
    self.valueError = 0.0
    self.atomName = atomName
    (self.seqCode, self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.molCode = defaultMolCode
    self.resLabel = resLabel
    self.figOfMerit = figOfMerit

    self.allValues = {iCode: self.value}

