"""
======================COPYRIGHT/LICENSE START==========================

conversionParsIO.py: NmrPipe conversion parameters writer

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

class NmrPipeConversionFile(NmrPipeGenericFile):

  def initialize(self,inputFile,outputFile,numDim,command):

    self.inputFile = inputFile
    self.outputFile = outputFile
    self.numDim = numDim
    self.command = command

    self.scriptName = 'conversionParsIO.py'

    self.dataValues = []
    for i in range(0,numDim):
      self.dataValues.append([])
        
    self.addBeginCommands = []
    self.addEndCommands = []
        
    self.dimRefs = ['x','y','z','a']
    self.dataRefs = ['NumberRealPoints','NumberComplexPoints','AcquisitionMode','SpectralWidth','CarrierFrequency','CarrierPpm','AxisName']
    self.dataRefCodes = ['N','T','MODE','SW','OBS','CAR','LAB']
    self.formats = ['%15d','%15d','%15s','%15.3f','%15.3f','%15.3f','%15s']
    
  def getDataValue(self,dim,index):
    
    try:
      return self.dataValues[dim][index].value

    except:
      print "  Error: nmrPipe data reference index %d, dim %d not found." % (index,dim)

    return None

  def setDataValue(self,expDimIndex,dataRef,value):
    dataRefIndex = self.dataRefs.index(dataRef)
    self.dataValues[expDimIndex][dataRefIndex].setValue(value)             

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing nmrPipe conversion file %s" % self.name


    fout = open(self.name,'w')

    # First lines...
    fout.write("#!/bin/csh" + self.newline * 2)

    #
    # Write command line
    # 

    fout.write("%s -in %s -noaswap" % (self.command,self.inputFile))

    # Get data for addBeginCommand(s) - are added as extra lines
    for addCommand in self.addBeginCommands:
      fout.write(" %s" % addCommand)

    # Finish command line
    fout.write(" \\" + self.newline)

    #
    # Write variable bits in middle
    # 

    for i in range(0,len(self.dataRefs)):
      dataRefCode = self.dataRefCodes[i]
      outstr = ""

      for dim in range(0,self.numDim):
        dimRef = self.dimRefs[dim]
        code = dimRef + dataRefCode
        outstr += ("  -%-5s%s" % (code,self.getDataValue(dim,i)))

      fout.write(outstr + " \\" + self.newline)

    #
    # Write end of command lines     
    # 

    fout.write("  -ndim %15d  -aq2D %15s " % (self.numDim,'States'))

    # Get data for addEndCommand(s) - are added as extra lines
    for addCommand in self.addEndCommands:
      fout.write(" %s" % addCommand)

    # Finish end of command lines     
    fout.write(" \\" + self.newline)
    fout.write("  -out %s -verb -ov" % (self.outputFile) + self.newline)

    fout.close()

    if verbose == 1:
      print self.scriptName + " finished..."

class NmrPipeConversion:

  def __init__(self,parent,index,value):
  
    self.index = index
    self.parent = parent
    self.setValue(value)
    # TODO: could add more info here (which dim, ...)

  def setValue(self,value):
  
    if value == None:
      print "Warning: no value for NMRPipe dataRef index %s. Set to zero." % (self.parent.dataRefs[self.index])
      value = 0

    self.value = (self.parent.formats[self.index]) % (value)

###################
# Main of program #
###################
    
if __name__ == "__main__":

  print "No test available"
  
