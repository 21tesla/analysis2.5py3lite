"""
======================COPYRIGHT/LICENSE START==========================

chemShiftsIO.py: I/O for Sparky chemical shift export files

Copyright (C) 2005-2008 Wim Vranken (European Bioinformatics Institute)

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

import re

# Import general functions
from memops.universal.Util import returnFloat, returnInt
from ccp.format.sparky.generalIO import SparkyGenericFile

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class SparkyChemShiftFile(SparkyGenericFile):

  """
  Information on file level
  """
  
  def initialize(self):
  
    self.chemShifts = []

    self.refData = [
    
               ['residue',1,1,7,None,"Group"],
               ['atom',1,1,7,'ndim',"Atom"],
               ['atomType',1,1,5,'ndim',"Nuc"],
               ['chemShift',1,1,9,None,"Shift"],
               ['chemShiftError',1,1,8,None,"SDev"],
               ['numAssignments',1,1,7,None,"Assignments"]
              
                   ]
          
    for (sparkyCode,numHeaderCol,numDataCol,colLen,refInfo,regExpPatt) in self.refData:

      self.patt[self.format + sparkyCode] = re.compile(regExpPatt)

  def read(self,verbose = 0):

    dataColCodes = []
    dataColLengths = []
    dataColNum = []

    colNums = {}

    if verbose == 1:
      print "Reading sparky chemical shift list %s" % self.name

    fin = open(self.name, 'rU')

    #
    # Start reading, look for header line
    #

    line = fin.readline()

    while line:
      cols = line.split()

      if len(cols) == 0 or self.patt['hash'].search(line):
        pass

      else:

        #
        # First line should have column header names
        #

        numHeaderCols = 0
        numDataCols = 0

        dataColLen = 0      # No space at end of line

        for refDatum in self.refData:

          (sparkyCode,numHeaderCol,numDataCol,colLen,refInfo,regExpPatt) = refDatum

          allMatches = self.patt[self.format + sparkyCode].findall(line)

          if allMatches:

            numMatches = len(allMatches)

            #
            # Number of total columns in header and for data
            #

            numHeaderCols += numHeaderCol * numMatches
            numDataCols += numDataCol * numMatches

            #
            # Keep track of column info
            #

            for i in range(0,numMatches):
              dataColCodes.append(sparkyCode)
              dataColNum.append(numDataCol)
              dataColLengths.append(colLen)

            #
            # This one to check data lines later
            #

            dataColLen += colLen * numMatches                    

        #
        # Exit if there's column names that are not recognized: no use trying
        # to read these for Sparky (very loose format)
        #

        if numHeaderCols != len(cols):
          print "  Unrecognized column headers: can't read sparky file... "
          return

        # Set order for column headers
        colNumKeys = colNums.keys()
        colNumKeys.sort()

        for colNum in colNumKeys:
          dataColCodes.append(colNums[colNum])

        #             
        # Start reading data
        #

        line = fin.readline()

        while line:

          cols = line.split()

          if len(cols) == 0 or self.patt['hash'].search(line):
            pass

          else:

            #
            # Two methods to read file:
            #  - default one based on plain column split
            #  - second one assuming valid data line, using expected column length to parse (only works in some cases)
            #

            obligColumnFound = 0

            if len(cols) == numDataCols or len(line) == dataColLen:
              dataColCodesLen = len(dataColCodes)

            else:

              #
              # If problems, still try to read in essential bits (all until 'ndim' column)
              #

              print "  Warning: Bad column match for following Sparky line. Trying to read...\n%s" % line,

              for i in range(0,len(dataColCodes)):            

                sparkyCode = dataColCodes[i]

                for (refSparkyCode,numHeaderCol,numDataCol,colLen,refInfo,regExpPatt) in self.refData:

                  if sparkyCode == refSparkyCode:

                    if refInfo == 'ndim' and obligColumnFound == 0:
                      obligColumnFound = 1
                      break

                    elif refInfo == None and obligColumnFound == 1:
                      obligColumnFound = -1
                      break

                #
                # Obligatory columns found, other ones starting: quit
                #

                if obligColumnFound == -1:
                  dataColCodesLen = i
                  break

              if obligColumnFound != -1:
                print "  Can't read... "
                line = fin.readline()
                continue
              else:
                print "  Reading first %d columns only." % dataColCodesLen
                print

            values = {}

            colNum = 0

            for i in range(0,len(dataColCodes)):

              sparkyCode = dataColCodes[i]

              if not values.has_key(sparkyCode):
                values[sparkyCode] = []

              if len(cols) == numDataCols or obligColumnFound == -1:

                values[sparkyCode].append(" ".join(cols[colNum:colNum + dataColNum[i]]))

                if dataColNum[i] > 1:
                  colNum += dataColNum[i] - 1

              else:

                values[sparkyCode].append(line[0:dataColLengths[i]])
                line = line[dataColLengths[i]:]

              colNum += 1


            self.chemShifts.append(SparkyChemShift(self,values))

          line = fin.readline()

      line = fin.readline()

    fin.close()

  def write(self,verbose = 0):

    """
 Group   Atom  Nuc    Shift   SDev  Assignments

   1G7      C  13C  171.302  0.000      1
   1G7     CA  13C   43.514  0.000      1
   1G7     DQ  13C  115.005  0.006      2
   2G7      C  13C  170.811  0.000      1
   2G7     CA  13C   44.076  0.000      1
    """
    
    if verbose == 1:
      print "Writing %s chemical shift list %s" % (self.format,self.name)

    fout = open(self.name,'w')

    fout.write("Group          Atom      Nuc      Shift      SDev     Assignments" + self.newline)
    fout.write(self.newline)
 
    #
    # Write out chem shifts
    #

    for chemShift in self.chemShifts:
      
      assignText = "%s%s%d" % (chemShift.molCode,chemShift.resLabel,chemShift.seqCode)
      numAssignments = chemShift.numAssignments
      if not numAssignments:
        numAssignments = 1

      fout.write("%-15s%-10s%-7s%9.3f   %7.3f   %5d" % (
                       assignText,
                    	 chemShift.atomName,
                    	 chemShift.atomType,
                    	 chemShift.value,
                    	 chemShift.valueError,
                       numAssignments))  # TODO not setting numAssignments
      
      fout.write(self.newline)
          
class SparkyChemShift:

  def __init__(self,parent,values):

    self.parent = parent
    
    self.resLabel = None
    self.seqCode = None
    self.atomName = None
    self.atomType = None
    self.value = None
    self.valueError = 0.0
    self.numAssignments = 0
    self.molCode = parent.defaultMolCode
    
    self.otherSeqAndInsertCodes = []

    for valueKey in values.keys():
      self.setValue(valueKey,values[valueKey])
    
    #
    # Also set chemical shift for other seqCodes if necessary
    # WARNING: this will only work if all values set from top...
    # if values set separately will fall over!
    
    for (otherSeqCode,otherSeqInsertCode) in self.otherSeqAndInsertCodes:
      values['residue'] = [self.resLabel + str(otherSeqCode) + otherSeqInsertCode]
      self.parent.chemShifts.append(SparkyChemShift(self.parent,values))

  def setValue(self,valueKey,values):

    if len(values) == 1:
      
      value = values[0]
     
      if not self.parent.patt['emptyline'].search(value):
      
        if valueKey == 'residue':
          value = value.strip()
          self.resLabel = value[0]
          
          seqCodes = value[1:].split(',')
          (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCodes[0])
          
          if self.seqCode == None:
            self.resLabel = value
          
          # Can have multiple values - need to take care of that...
          if len(seqCodes) > 1:
            for seqCode in seqCodes[1:]:
              self.otherSeqAndInsertCodes.append(getSeqAndInsertCode(seqCode))

        elif valueKey == 'atom':
          self.atomName = value.strip()

        elif valueKey == 'atomType':
          self.atomType = value.strip()

        elif valueKey == 'chemShift':
          self.value = returnFloat(value)

        elif valueKey == 'chemShiftError':
          self.valueError = returnFloat(value)

        elif valueKey == 'numAssignments':
          self.numAssignments = returnInt(value)

    else:
    
      print "  Warning: empty or invalid data field %s for Sparky chemical shift file." % valueKey
      pass
