#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

generalIO.py: General I/O information for Bruker files

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

from ccp.format.general.formatIO import FormatFile

from memops.universal.Util import returnFloat
from memops.universal.Util import returnFloats
from memops.universal.Util import returnInt
from memops.universal.Util import returnInts

#####################
# Class definitions #
#####################

class BrukerGenericFile(FormatFile):

  def initialize(self):
  
    pass

  def setGeneric(self):
    
    self.format = 'bruker'

class BrukerParFile(BrukerGenericFile):

  def initialize(self,parHelp = None):
  
    self.parameters = {}
    
    self.parHelp = parHelp # This is a dict: corresponds to BrukerParHelp().tags
    
    self.infoStrings = []
    self.headerDict = {}

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading Bruker acquisition parameters file %s" % self.name
    
    # For reading: based on ##
    # If value between <> or non-numeric characters (except point): is string
    # If value numeric: float if contains '.'
    #                   integer if no '.'
    #                   (data type doesn't matter for python comparisons (5 == 5.0))

    fin = open(self.name, 'rU')
    line = fin.readline()
    
    while line:
      
      #
      # Checking of empty lines not really necessary, but left in just in case
      #
      
      if not (self.patt['emptyline'].search(line)):
      
        cols = line.split()
        
        if self.patt[self.format + 'StartDoubleHashDollar'].search(cols[0]):
          
          #
          # This is a variable
          #
          
          parname = self.patt[self.format + 'StartDoubleHashDollar'].sub('',cols[0])
          parname = parname.replace('=','')

          self.parameters[parname]=BrukerPar(parname)
          #self.parameters[parname].addHelp(self.parHelp,verbose)

          # Handle value(s): is or are always in list
          valuestring = ' '.join(cols[1:])
          values = []
          valueCount = 1

          # Is an array? Indicated by (n..m) search above
          if self.patt[self.format + 'BracketsListIndicator'].search(valuestring):

            # Values can start IMMEDIATELY after right bracket! Split accordingly
            if self.patt[self.format + 'BracketCloseNumber'].search(valuestring):
              (listinfo,line) = valuestring.split(')')
            else:
              listinfo=valuestring
              line = ''

            # For info on list, get rid of brackets and split on '..', find value_len
            listinfo = self.patt['bracketEither'].sub('',listinfo)
            (listStart,listEnd) = returnInts(listinfo.split('..'))
            valueLen = listEnd - listStart + 1

            # Report error if valueMin not zero!!
            if listStart != 0:
              print "Error: Bruker array starts at non-zero for %s!" % parname

            # Start parsing list values here. If they were behind the right bracket,
            # the string containing them is the current line.
            # LISTS LIKE THIS ARE ALWAYS CONSIDERED FLOAT DATE TYPE!
            valueType = 'float'
            valueCount = 0

            while 1 == 1:
              # Always check if this line doesn't have '##'
              if self.patt[self.format + 'StartDoubleHash'].search(line):
                print "Parsing error:%s" % line + self.newline
              cols = line.split()
              values.extend(returnFloats(cols))
              valueCount += len(cols)
              if valueCount < valueLen:
                line = fin.readline()
              else:
                break

          # Is a text value (between <>)? If so, need to find matching '>',
          # since this can be multi-line.
          elif self.patt[self.format + 'SharpBracketOpen'].search(valuestring):
            while not self.patt[self.format + 'SharpBracketClose'].search(valuestring):
              line = fin.readline()
              # Always check if this line doesn't have '##'
              if self.patt[self.format + 'StartDoubleHash'].search(line):
                print "Parsing error:%s" % line
              line = line.rstrip()
              valuestring += line
            valuestring = self.patt[self.format + 'SharpBracketEither'].sub('',valuestring)
            values.append(self.patt[self.format + 'FinalReturn'].sub('',valuestring))
            valueType = 'string'

          # Is a number
          elif self.patt[self.format + 'DigitDotDigit'].search(valuestring):
            # Is a float
            if self.patt[self.format + 'Dot'].search(valuestring):
              valueType = 'float'
              values.append(returnFloat(valuestring))
            # Is an integer
            else:
              valueType = 'int'
              values.append(returnInt(valuestring))

          # Is a string (unquoted)
          else:        
            valueType = 'string'
            values.append(valuestring)

          # Add value to parameter
          self.parameters[parname].addValues(valueCount,values)
          self.parameters[parname].addType(valueType)
          self.parameters[parname].addHelp(self.parHelp,verbose)

        elif self.patt[self.format + 'StartDoubleDollar'].search(cols[0]):
          # This is a date line or the original file name with path 
          self.infoStrings.append(line[3:].strip())
        elif self.patt[self.format + 'StartDoubleHash'].search(cols[0]):
          # General information about parameter file
          (headerKeywd,headerValue) = line.strip().split("=")
          self.headerDict[headerKeywd[2:].strip()] = headerValue.strip()
        else:
          print "Error: undefined basic data line:%s" % line + self.newline

      line=fin.readline()

class BrukerPar:

  def __init__(self,name):
    # Bruker name of parameter (string)
    self.name = name

  def addValues(self,valueLen,values):
    # Add value
    if valueLen != len(values):
      print "Error: given valueLen does not match number of values (%s)" % self.name
      self.valueLen = len(values)
    else:
      # Number of values for this parameter
      self.valueLen = valueLen
    
    # Values (always make this array for ease of use)
    self.values = values
    
  def addType(self,valueType):

    self.type = valueType
      
  def addHelp(self,parHelp,verbose):
  
    if parHelp and parHelp.has_key(self.name):

      self.help = parHelp[self.name]

    else:
      if verbose > 1:
        print "Parameter %s: no help string found." % self.name
      self.help = ''

def getMatrixFilePath(inputPath, numDim):
  """
  
  Input: either a directory that contains the matrix file, 
         or a random file in that directory.
  Output:
         
  """
  if os.path.isdir(inputPath):
    dirPath = inputPath
  else:
    dirPath =os.path.dirname(inputPath)
  
  procDataFile = '%d%s' % (numDim,numDim * 'r')
  
  result = os.path.join(dirPath, procDataFile)
  if not os.path.isfile(result):
    result = inputPath
  #
  return result
