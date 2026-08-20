"""
======================COPYRIGHT/LICENSE START==========================

processingScriptIO.py: Azara processing script file writer

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
from ccp.format.azara.generalIO import AzaraGenericFile
from ccp.format.general.Constants import indent

#####################
# Class definitions #
#####################

class AzaraProcessingScriptFile(AzaraGenericFile):

  def initialize(self,inputFile,outputFile,parFile,numDim):

    self.inputFile = inputFile
    self.outputFile = outputFile
    self.parFile = parFile
    self.numDim = numDim
    
    self.scriptName = 'processingParsIO.py'


    self.dataValues = []
    
    for i in range(0,numDim):
      self.dataValues.append([])
        
    self.addBeginCommands = []
    self.addEndCommands = []
        
    self.dimRefs = ['1','2','3','4']
    self.dataRefs = ['NegateAlternate','DataType','Conjugate','Varian','Rescale','BrukerAvance','WindowFunction','Zerofill','FourierTransform','BrukerAvancePhase','Phase','Reduce']
    
  def getDataValue(self,dim,index):    
  
    try:
      return self.dataValues[dim][index].value

    except:
      print "  Error: azara data reference index %d, dim %d not found." % (index,dim)

  def setDataValue(self,expDimIndex,dataRef,value):
    dataRefIndex = self.dataRefs.index(dataRef)
    self.dataValues[expDimIndex][dataRefIndex].setValue(value)             

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing azara processing script file %s" % self.name

    fout = open(self.name,'w')

    #
    # Write first lines...
    #

    fout.write("input %s" % self.inputFile + self.newline * 2)
    fout.write("output %s" % self.outputFile + self.newline * 2)
    fout.write("par %s" % self.parFile + self.newline * 2)

    #  
    # Get data for addBeginCommand(s) - are added as extra lines
    #

    for addCommand in self.addBeginCommands:
      fout.write("%s%s" % (addCommand,self.newline*2))

    #
    # Write variable bits in middle
    # 

    numFuncs = len(self.dataRefs)

    for dim in range(0,self.numDim):

      dimRef = self.dimRefs[dim]

      fout.write("script_com %s" % dimRef + self.newline)

      for i in range(0,numFuncs):

        value = self.getDataValue(dim,i)
        dataRef = self.dataRefs[i]

        if value != '':
          fout.write("%s%-25s! %s%s" % (indent,value,dataRef,self.newline))

      # End current script
      fout.write("end_script" + (self.newline * 2))

    fout.write(self.newline)

    #
    # Get data for addEndCommand(s) - are added as extra lines
    #

    for addCommand in self.addEndCommands:
      fout.write("%s%s" % (addCommand,self.newline*2))

    fout.close()

    if verbose == 1:
      print self.scriptName + " finished..."

class AzaraProcessingScript:

  def __init__(self,parent,index,value):
  
    self.index = index
    self.parent = parent
    self.setValue(value)
    # TODO: could add more info here (which dim, ...)

  def setValue(self,value):
  
    if value == None:
      print "Warning: no value for Azara dataRef index %s. Set to zero." % (self.parent.dataRefs[self.index])
      value = 0

    self.value = value

###################
# Main of program #
###################
    
if __name__ == "__main__":

  print "No test available"
  
