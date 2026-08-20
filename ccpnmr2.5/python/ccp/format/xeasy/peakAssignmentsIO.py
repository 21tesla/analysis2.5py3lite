#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

peakAssignmentsIO.py: I/O for XEasy peak assigment files

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


import os

# Import general functions
from memops.universal.Util import returnLong, returnLongs
from memops.universal.Util import returnInt, returnInts
from ccp.format.xeasy.generalIO import XEasyGenericFile
from memops.universal.Io import getTopDirectory

#####################
# Class definitions #
#####################

class XEasyPeakAssignmentFile(XEasyGenericFile):

  """
  XEasy .assign file: describes ambiguous assignments for a peak list
  """
  
  def initialize(self):
  
    self.peakAssignments = []
    self.assignmentsDim = []
  
  def setDimensions(self,ndim):
  
    self.numDim = ndim

  def read(self,verbose = 0):
  
    """
    Reader for XEasy .assign file format
    """

    if verbose == 1:
      print "Reading xeasy peak assignments file %s" % self.name

    fin = open(self.name, 'rU')

    #
    # Read number of dimensions on first line
    #

    line = fin.readline()
    cols = line.split()

    n = returnInt(cols[0])

    #
    # Create the assignment file
    #

    self.setDimensions(n)

    #
    # Read rest file
    #

    line = fin.readline()

    while line:

      cols = line.split()

      #
      # Line with peak number
      #

      if self.patt['hash'].search(line):

        peaknum = cols[1]
        assignments = []

      #
      # Line with atom serials
      #

      elif int(cols[0]) + 1 == len(cols):
        assignments.append(cols[1:])

      #
      # Line with matrix binary encoding
      #

      elif len(cols) == 4:

        #
        # First create peak assignment
        #

        peakAssignment = XEasyPeakAssignment(self,peaknum)
        self.peakAssignments.append(peakAssignment)

        #
        # Then handle binary format information
        #

        peakAssignment.setAssignCombinations(assignments,cols)

      else:

        print "Unrecognized line in XEasy assignment file:\n%s\n" % line

      line = fin.readline()

    fin.close()

  def write(self, verbose = 0):

    if verbose == 1:
      print "Writing xeasy peak assignments file %s" % self.name

    fout = open(self.name,'w')

    numDim = self.numDim

    #
    # Write out header
    #

    fout.write("%d" % numDim)
    fout.write(self.newline)

    #
    # Write out assignments
    #

    for peakAssignment in self.peakAssignments:

      #
      # Write peaknumber
      #

      fout.write("# %d" % peakAssignment.peakNum)
      fout.write(self.newline)

      binCodes = peakAssignment.getBinaryCodes()

      #
      # Write atom serials per dimension
      #

      for dim in range(0,numDim):
        if self.assignmentsDim.count(dim) == 1:
          assignIndex = self.assignmentsDim.index(dim)
          assignments = peakAssignment.assignments[assignIndex]
          numAssignments = len(assignments)
        else:
          assignments = []
          numAssignments = 0

        fout.write("%d" % numAssignments)

        for i in range(0,numAssignments):
          fout.write(" %d" % assignments[i])

        fout.write(self.newline)

      #
      # Write binary codes for matrix
      #

      line = ""

      for i in range(0,4):

        line += "%d " % binCodes[i]

      fout.write(line)
      fout.write(self.newline)

    fout.close()

class XEasyPeakAssignment:

  #
  # Do all encoding in here - none of this on top level!!!
  #

  def __init__(self,parent,peakNum):

    self.peakNum = returnInt(peakNum)
    self.parent = parent
  
    self.assignments = []
    self.assignCombinations = []

  #
  # Code to use and produce binary matrix encoding in XEasy file
  # Calculation method taken from Jens Linge's ReadXeasy.py script in the ARIA (pre-v2)
  # distributions.
  #

  def setAssignCombinations(self,assignments,binCodes):
  
    matrixDims = []
    assignmentMatrix = []
    
    #
    # Note that empty assignment dimensions will NOT be added
    # self.assignmentsDim keeps track of which dimension was used
    #
    
    for i in range(0,len(assignments)):
      
      assignmentDim = assignments[i]
          
      if len(assignmentDim) > 0:
        matrixDims.append(len(assignmentDim))
        self.assignments.append(returnInts(assignmentDim))
        
        if len(self.parent.assignmentsDim) < 2:
          self.parent.assignmentsDim.append(i)
                
    if len(matrixDims) > 2:
      print "  Error in XEasy peak assignment file: more than two dimensions have atom serials."
      
    #
    # Setup assignment matrix
    #
    
    for i in range(0,matrixDims[0]):
      
      assignmentMatrix.append(matrixDims[1] * [1])

    #
    # Decode the binary codes
    #

    binCodes = returnLongs(binCodes)
    
    for i in range(0,4):

      binCode = binCodes[i] + 1

      for exponent in range(31,-1,-1):
        if binCode == 0:
          continue
        
        #
        # Handle overflow
        #
        
        elif binCode > 0:

          binCode = - (2147483648L - binCode) - 2147483648L
          
        if binCode <= -(2L ** exponent):
          
          binCode = binCode + 2L ** exponent
          totalExponent = 32 * i + exponent
          row = totalExponent % 10
          col = totalExponent / 10

          try:
          #
          # TODO: Why is this necessary?!?! Discuss with Jens...
          #
            if len(assignments) == 2:
              assignmentMatrix[row][col] = 0
            elif len(assignments) == 3:
              assignmentMatrix[col][row] = 0
           
          except:
              print 'WARNING: problem with assignment matrix!'
              print self.peakNum, assignmentMatrix, row, col
    
    #
    # Make a list of allowed combinations...
    #
    
    for i in range(0,len(assignmentMatrix)):
      for j in range(0,len(assignmentMatrix[i])):
        if assignmentMatrix[i][j] == 1:
          self.assignCombinations.append([self.assignments[0][i],self.assignments[1][j]])
  
  def getBinaryCodes(self):
  
    #
    # Setup the assignment list for figuring out the binary codes
    #
    
    self.assignments = []
    
    for i in range(0,2):
      self.assignments.append([])
    
    for i in range(0,len(self.assignCombinations)):
      for j in range(0,2):
        if self.assignments[j].count(self.assignCombinations[i][j]) == 0:
          self.assignments[j].append(self.assignCombinations[i][j])

    
    #
    # Set up assignment matrix
    #
    
    # Hack (see setAssignCombinations) for switching col/row
    if self.parent.numDim == 2:
      assiRow = 0
      assiCol = 1
    else:
      assiRow = 1
      assiCol = 0

    assignmentMatrix = []
    
    for i in range(0,len(self.assignments[assiRow])):
      assignmentMatrix.append([])
      for j in range(0,len(self.assignments[assiCol])):
        
        if assiRow < assiCol and self.assignCombinations.count([self.assignments[0][i],self.assignments[1][j]]) != 0:
          matrixElement = 1
        elif assiRow > assiCol and self.assignCombinations.count([self.assignments[0][j],self.assignments[1][i]]) != 0:
          matrixElement = 1
        else:
          matrixElement = 0
          
        assignmentMatrix[-1].append(matrixElement)
        
    
    #
    # Get the binary code (based on calculation method Jens Linge - see higher)
    #
    
    binCodes = [-1L, -1L, -1L, -1L]
    
    for row in range(0, len(self.assignments[assiRow])):
      for col in range(0, len(self.assignments[assiCol])):
        if assignmentMatrix[row][col] == 0:  #1 means assigned
        
          findBinCode = 10 * col + row
          if findBinCode < 32:
              binCodeEl = 0
          elif findBinCode < 64:
              binCodeEl = 1
          elif findBinCode < 96:
              binCodeEl = 2
          elif findBinCode < 128:
              binCodeEl = 3
              
          exponent = findBinCode - 32 * binCodeEl
          
          #
          # Deal with overflow
          #
          
          if binCodes[binCodeEl] - 2L ** exponent > -2147483469L:
            binCodes[binCodeEl] = binCodes[binCodeEl] - 2L ** exponent
          else:
            binCodes[binCodeEl] = binCodes[binCodeEl] + 2L ** exponent
      
    return binCodes
  
###################
# Main of program #
###################

if __name__ == "__main__":  
                                                                                                            
  files = [['../reference/xeasy/2d.1.assign','local/paf1.test'],
           ['../reference/xeasy/2d.2.assign','local/paf2.test'],
           ['../reference/xeasy/n15.assign','local/paf3.test'],
           ['../reference/xeasy/3d.1.assign','local/paf4.test']
          ]

  for (inFile,outFile) in files:
    
    assignFile = XEasyPeakAssignmentFile(os.path.join(getTopDirectory(),inFile))
    
    assignFile.read(verbose = 1)

    print assignFile.assignmentsDim
    for pa in assignFile.peakAssignments:
      
      print pa.peakNum, pa.assignments, pa.assignCombinations
        
    # Writing tested on 19/03/2002 (Wim) - all OK
    # Have to switch row/col for 3d compared to 2d: hacks included in script
    assignFile.name = outFile
    assignFile.write(verbose = 1)
    
