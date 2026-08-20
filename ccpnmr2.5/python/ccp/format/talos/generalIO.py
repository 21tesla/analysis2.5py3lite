"""
======================COPYRIGHT/LICENSE START==========================

generalIO.py: General I/O information for Talos files

Copyright (C) 2005-2010 Wim Vranken (European Bioinformatics Institute)

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

from ccp.format.general.formatIO import FormatFile
from ccp.format.general.Constants import defaultMolCode

#####################
# Class definitions #
#####################

class TalosGenericFile(FormatFile):

  def setGeneric(self):
    
    self.format = 'talos'
    self.defaultMolCode = defaultMolCode

    self.refData = []
    self.varCodes = []
    self.formatCodes = []
    
    self.numCols = None

  def read(self, verbose = 0):
  
    if verbose == 1:
      self.printInfo('Reading')
      
    fin = open(self.name, 'rU')

    #
    # Start reading
    #

    line = fin.readline()

    while line:

      if not self.patt['emptyline'].search(line) and not self.patt[self.format + 'Remark'].search(line):

        if self.patt[self.format + 'Dataline'].search(line):

          self.handleDataLine(line)

        elif self.patt[self.format + 'Vars'].search(line):

          #
          # Get vars line, check if matches up with reference info
          #

          cols = line.split()
          colNum = 1

          for refDatum in self.refData:

            refVarCode = refDatum[0]

            # Bail out if obligatory header columns not present
            if refVarCode not in cols:
              if refDatum[4]:
                print cols, refVarCode, refDatum
                print "  Error: can't read %s file - headers don't match." % self.format
                return False
              else:
                continue
            else:
              curColNum = cols.index(refVarCode) + 1

            while curColNum > colNum + 1:
              colNum += 1
              self.varCodes.append(None)

            self.varCodes.append(refVarCode)
            colNum += 1

        elif self.patt[self.format + 'Format'].search(line):

          self.formatCodes = line.split()[1:]
          self.numCols = len(self.formatCodes)

        else:
        
          cols = line.split()

          if len(cols) == self.numCols:

            varsDict = {}

            for varCode in self.varCodes:
              varsDict[varCode] = cols.pop(0)

            self.setVarsLine(varsDict)

          else:

            print "  Warning: unparsed %s line:\n%s" % (self.format,line)

      line = fin.readline()

    fin.close()

  def write(self, use_fout = None, verbose = 0):

    if verbose == 1:
      self.printInfo('Writing')
      
    #
    # Get the data items
    #
    
    self.getVarItems()
    
    #
    # Start writing file
    #
    
    if not use_fout:
    
      fout = open(self.name,'w')

      #
      # Write out header
      #

      fout.write("REMARK File written by CcpNmrFormat converter.")
      fout.write(self.newline * 2)
    
    else:
    
      fout = use_fout   


    self.writeDataLines(fout)

    if self.varItems:
    
      fout.write(self.newline)

      fout.write("VARS  ")

      formatCodes = []

      for refDatum in self.refData:

        refVarCode = refDatum[0]
        refFormatCode = refDatum[1]

        if type(refVarCode) == type([]):

          for dimRefCode in self.refDimCodes:

            for refVarCodeSingle in refVarCode:

              outString = refVarCodeSingle.replace('n',dimRefCode)

              fout.write(" %s" % outString)

              formatCodes.append(refFormatCode)

        else:

          fout.write(" %s" % refVarCode)
          formatCodes.append(refFormatCode)

      fout.write(self.newline)

      fout.write("FORMAT")

      formatNum = 0

      for formatCode in formatCodes:

        fout.write(" " + str(formatCode))

      fout.write(self.newline * 2)

      #
      # Write out data items
      #

      for varItem in self.varItems:

        line = ""

        for refDatum in self.refData:

          refVarCode = refDatum[0]
          refFormatCode = refDatum[1]
          refAttrName = refDatum[2]
          refFunc = refDatum[3]

          value = getattr(varItem,refAttrName)

          line += (refFormatCode % refFunc(value))
          line += " "

        fout.write(line[:-1])
        fout.write(self.newline)

    if not use_fout:
      fout.close()

  def setVarsLine(self,varsDict):
    pass
    
  def printInfo(self,action):
    pass
    
  def handleDataLine(self,line):
    pass

  def writeDataLines(self,fout):
    pass

  def getVarItems(self):
    
    self.varItems = []
    
class TalosGenericDataItem:

  def __init__(self,parent):
  
    self.parent = parent
  
    for refDatum in parent.refData:

      refVarCode = refDatum[0]
      refFormatCode = refDatum[1]
      refAttrName = refDatum[2]
      refFunc = refDatum[3]
      
      setattr(self,refAttrName,None)
      
    self.initItemSpecificData()

  def setVarsDict(self,varsDict, refData = None):
  
    if not refData:
      refData = self.parent.refData
      
    for refDatum in refData:

      refVarCode = refDatum[0]
      refFormatCode = refDatum[1]
      refAttrName = refDatum[2]
      refFunc = refDatum[3]
     
      setattr(self,refAttrName,refFunc(varsDict[refVarCode]))

    self.setItemSpecificVars()
      
