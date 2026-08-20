#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

rdcConstraintsIO.py: I/O for PALES program rdc constraint files
 (see http://www.mpibpc.mpg.de/groups/zweckstetter/_links/software_pales.htm)

Copyright (C) 2010 Wim Vranken (European Bioinformatics Institute)
                   & Brian Smith (University of Glasgow)

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
- MSD website (http://www.ebi.ac.uk/msd/)

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
Ionides and Ernest D. Laue. The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Accepted by Proteins (2004).

===========================REFERENCE END===============================
"""

import os, string

from memops.universal.Util import returnInt, returnFloat
from memops.universal.Io import getTopDirectory

from ccp.format.pales.generalIO import PalesGenericFile
from ccp.format.pales.generalIO import PalesConstraintItem
from ccp.format.pales.generalIO import PalesConstraintMember

#####################
# Class definitions #
#####################

class PalesRdcConstraintFile(PalesGenericFile):

  def initialize(self):

    self.constraints = []

    self.constraintElements = 2
    
    self.writeKeywds = {}

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s rdc constraint list %s" % (self.format,self.name)

    rdcId = 0
    chainCode = self.defaultMolCode

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      cols = line.split()

      if len(cols) == 0 or self.patt['hash'].search(line):
        pass

      # skip sequence for now - sequence must be loaded from elsewhere
      elif cols[0] == "DATA":
        pass

      # skip VARS listing for now
      elif cols[0] == "VARS": 
        pass

      # skip format line for now
      elif cols[0] == "FORMAT":
        pass

      else:

        rdcId += 1


        self.constraints.append(PalesRdcConstraint(rdcId))
        self.constraints[-1].setRdcData(cols[6], cols[7], cols[8])
	# next line reads residue type, but I don't think the molecule, chain etc. gets created - sequence must be loaded from elsewhere
        self.constraints[-1].setAtomMembers(chainCode,cols[0],cols[1],cols[2],cols[3],cols[4],cols[5])

      line = fin.readline()

  def write(self, oneLetterSequence = "", verbose = 0):
    
    #
    # Output format is:
    #
    # DATA SEQUENCE MQIFVKTLTG KTITLEVEPS DTIENVKAKI QDKEGIPPDQ QRLIFAGKQL
    # DATA SEQUENCE EDGRTLSDYN IQKESTLHLV LRLRGG
    # 
    # VARS   RESID_I RESNAME_I ATOMNAME_I RESID_J RESNAME_J ATOMNAME_J D      DD    W
    # FORMAT %5d     %6s       %6s        %5d     %6s       %6s    %9.3f   %9.3f %.2f
    # 
    #    2    GLN      N      2    GLN     HN     -15.524     1.000 1.00
    # ...

    if verbose == 1:
      print "Writing %s rdc constraint list %s" % (self.format,self.name)

    fout = open(self.name,'w')

    self.writePalesSequence(fout, oneLetterSequence)
    
    fout.write("VARS   RESID_I RESNAME_I ATOMNAME_I RESID_J RESNAME_J ATOMNAME_J D      DD    W")
    fout.write(self.newline)
    fout.write("FORMAT %5d     %6s       %6s        %5d     %6s       %6s    %9.3f   %9.3f %.2f")
    fout.write(self.newline)
    fout.write(self.newline)
    
    for constraint in self.constraints:

      for item in constraint.items:
        for member in item.members:
          atomName = member.atomName
          seqCode = member.seqCode
	  resLabel = member.resLabel

	  fout.write("%5d%6s%6s" % (seqCode,resLabel,atomName))
        
        if constraint.error:
          error = constraint.error
        else:
          error = 0.5 # sensible minimum

        if constraint.weight:
          weight = constraint.weight
        else:
          weight = 1.00 # normally makes sense
        
        fout.write("%9.3f %9.3f %.2f" % (constraint.value,error,weight))
        fout.write(self.newline)

    fout.close()


  def writePalesSequence(self, fout, oneLetterSequence = "", verbose = 0):
  
    # one letter code sequence writer in mangeable chunks with "DATA SEQUENCE " at the start of each line
    # NB fout already opened in write mode. Should probably have some write error handler?
      
    # 50 code letters to a line. Work out how many lines
    chunkSize = 50
    (numChunks, lenLastLine) = divmod(len(oneLetterSequence),chunkSize)
    if (lenLastLine > 0):
      numChunks += 1

    for chunk in range(numChunks):
      firstEl = chunk*chunkSize
      lastEl  = ((chunk+1)*chunkSize)
      if (chunk == (numChunks - 1)):
	#in this case lastEl is set to the end of the string (i.e. the character after the last letter we want)
	lastEl = len(oneLetterSequence)
      # seq is of type ccp.api.molecule.Molecule.Molecule.seqString which is a String so slice it
      # slicing gives sequence from firstEl up to charater before lastEl
      seqChunk = oneLetterSequence[firstEl:lastEl]
      i = 0
      seqString = ""
      for seqEl in seqChunk:
	if ((i % 10) == 0):
	  # space every 10 characters
	  seqString += " "

	seqString += seqEl
	i+=1

      fout.write("DATA SEQUENCE %s" % seqString)
      fout.write(self.newline)

    fout.write(self.newline * 2)


class PalesRdcConstraint:

  def __init__(self,Id):
    
    self.Id = returnInt(Id)
    self.items = []
    
  def setRdcData(self,value,error = 0.00,weight = 1.00):
  
    self.value = returnFloat(value)
    self.error = returnFloat(error)
    self.weight = returnFloat(weight)
      
  def setAtomMembers(self,chainCode,seqCode1,resLabel1,atomName1,seqCode2,resLabel2,atomName2): 

    self.items.append(PalesConstraintItem())
     
    self.items[-1].members.append(PalesConstraintMember(chainCode,returnInt(seqCode1),string.upper(resLabel1),string.upper(atomName1)))
    self.items[-1].members.append(PalesConstraintMember(chainCode,returnInt(seqCode2),string.upper(resLabel2),string.upper(atomName2)))
 
###################
# Main of program #
###################

if __name__ == "__main__":

  files = ['../../reference/ccpNmr/pales/pales_example.pal']
  
  for file in files:
    
    file = os.path.join(getTopDirectory(), file)
    
    constraintFile = PalesRdcConstraintFile(file)

    constraintFile.read(verbose = 1)

    for constraint in constraintFile.constraints:
      print constraint.Id,

      print constraint.value, constraint.error,

      for item in constraint.items:
        for member in item.members:

          print member.seqCode, member.resLabel, member.atomName,

        print "|",

      print
    
    constraintFile.name = '../../../reference/ccpNmr/pales/rdc.testout'

    constraintFile.write(verbose = 1)

