"""
======================COPYRIGHT/LICENSE START==========================

processingScriptIO.py: NmrPipe processing script writer

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
from ccp.format.nmrPipe.generalIO import NmrPipeGenericFile

#####################
# Class definitions #
#####################

class NmrPipeProcessingScriptFile(NmrPipeGenericFile):

  def initialize(self,inputFile,outputFile,numDim):

    self.inputFile = inputFile
    self.outputFile = outputFile
    self.numDim = numDim

    self.scriptName = 'processingScriptIO.py'
    
    if self.numDim > 2:
      self.inCommand = 'xyz2pipe'
      self.outCommand = 'pipe2xyz'
      self.inArg = ' -x'
      self.outArg = ' -z -inPlace'

    else:
      self.inCommand = self.outCommand = 'nmrPipe'
      self.inArg = ''
      self.outArg = ' -ov'

    self.dataValues = []
    self.deleteImaginary = []
    for i in range(0,numDim):
      self.dataValues.append([])
      self.deleteImaginary.append(0)
        
    self.addBeginCommands = []
    self.addEndCommands = []
        
    self.dimRefs = ['x','y','z','a']
    self.dataRefs = ['NegateAlternate','Conjugate','WindowFunction','Zerofill','FourierTransform','Phase']
    
    self.dataRefCodes = ['SIGN','SIGN','SP','ZF','FT','PS']
    
  def getDataValue(self,dim,index):
    
    try:
      return self.dataValues[dim][index].value

    except:
      print "  Error: nmrPipe data reference index %d, dim %d not found." % (index,dim)

    return None

  def getDataValueRefCode(self,dim,index):
    # This one in case there's a change in the default dataRefCode...
    return self.dataValues[dim][index].dataRefCode
  
  def setDataValue(self,expDimIndex,dataRef,value,dataRefCode = None):
    dataRefIndex = self.dataRefs.index(dataRef)
    self.dataValues[expDimIndex][dataRefIndex].setValue(value)             

    # Allow the dataRefCode to change if so desired (for window functions)
    if dataRefCode:
      self.dataValues[expDimIndex][dataRefIndex].setDataRefCode(dataRefCode)

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing nmrPipe processing file %s" % self.name

    fout = open(self.name,'w')

    # First lines...
    fout.write("#!/bin/csh" + self.newline * 2)

    #
    # Write command line
    # 

    fout.write("%s -in %s%s -verb\\" % (self.inCommand,self.inputFile,self.inArg) + self.newline)

    # Get data for addBeginCommand(s) - are added as extra lines
    for addCommand in self.addBeginCommands:
      fout.write(" %s \\" % addCommand)

    #
    # Write variable bits in middle
    # 

    numFuncs = len(self.dataRefs)

    for dim in range(0,self.numDim):
      dimRef = self.dimRefs[dim]

      for i in range(0,numFuncs):
        dataRefCode = self.dataRefCodes[i]
        value = self.getDataValue(dim,i)

        dataRefCodeUpdate = self.getDataValueRefCode(dim,i)

        if dataRefCodeUpdate:
          dataRefCode = dataRefCodeUpdate

        # Only write if value not empty
        if value != None:
          # Delete imaginary
          if i == numFuncs-1 and self.deleteImaginary[dim] == 1:
            value += ' -di'

          self.writefunc(dataRefCode,value,fout)

      if dimRef == 'x':
        # Transpose after x
        fout.write("| nmrPipe -fn %-22s\\" % ('TP') + self.newline)

      if dimRef == 'y' and self.numDim > 2:
        # TODO: make func for writing in/out command and pass args to that
        # End of processing 2d for 3d:
        fout.write("| %s -out %s -y -ov" % (self.outCommand,self.outputFile) + self.newline)
        fout.write(self.newline)

        # Start processing z
        fout.write("%s -in %s -z -verb\\" % (self.inCommand,self.outputFile) + self.newline)


    # Get data for addEndCommand(s) - are added as extra lines
    for addCommand in self.addEndCommands:
      fout.write(" %s \\" % addCommand)

    #
    # Write final command line
    # 

    fout.write("| %s -verb -out %s%s" % (self.outCommand,self.outputFile,self.outArg) + self.newline)

    fout.close()

    if verbose == 1:
      print self.scriptName + " finished..."

  def writefunc(self,func,arg,fout):
    
    fout.write("| nmrPipe -fn %-22s\\" % (func + ' ' + arg) + self.newline)

class NmrPipeProcessingScript:

  def __init__(self,parent,index,value):
  
    self.index = index
    self.parent = parent
    self.setValue(value)
    self.dataRefCode = None
    # TODO: could add more info here (which dim, ...)

  def setValue(self,value):
  
    self.value = value

  def setDataRefCode(self,value):
  
    self.dataRefCode = value

###################
# Main of program #
###################
    
if __name__ == "__main__":

  print "No test available"
  
