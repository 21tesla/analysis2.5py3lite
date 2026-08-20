"""
======================COPYRIGHT/LICENSE START==========================

peaksIO.py: I/O for Sparky peak list export files

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
from memops.universal.Util import returnFloat, returnFloats
from memops.universal.Util import returnInt
from ccp.format.sparky.generalIO import SparkyGenericFile
from memops.universal.Io import splitPath

#####################
# Class definitions #
#####################

      
class SparkyPeakFile(SparkyGenericFile):
  
  """
  Information on file level
  """
  
  def initialize(self):
    
    self.peaks = []
    
    self.specFiles = [None]
    self.specNumPoints = [None]
    
    self.assignSep = '-'
    
    # Only relevant for export!
    self.volumeIntensityType = ['volume']
    
    
    #
    # Field info:
    #
    # (sparkyCode,numHeaderCol,numDataCol,colLen,refInfo,regExpPatt) 
    #
    
    self.refData = [
    
           ['assignment',1,1,18,None,"Assignment"],
           ['freqDimsPpm',1,1,11,'ndim',"(w[1-5]\s*)"],
           ['freqDimsHz',1,1,11,'ndim',"(w[1-5]\s+\(Hz\)\s*)"],
           ['volume',1,2,14,None,"Volume"],
           ['dataHeight',2,1,13,None,"Data Height"],
           ['fitHeight',2,1,13,None,"Fit Height"],
           ['height',1,1,11,None,"Height"],
           ['intensity',1,1,13,None,"Intensity"],
           ['details',1,1,0,None,'Note']
              
                   ]
          
    for (sparkyCode,numHeaderCol,numDataCol,colLen,refInfo,regExpPatt) in self.refData:

      self.patt[self.format + sparkyCode] = re.compile(regExpPatt)

  def setSpectrumInfo(self,specName,ndim):
  
    self.specNames = [specName]
    self.numDims = [ndim]

  def read(self,verbose = 0):

    dataColCodes = []
    dataColLengths = []
    dataColNum = []

    ndim = 0

    if verbose == 1:
      print "Reading sparky peak list %s" % self.name

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

        dataColLen = 2      # Space and newline at end of line
        headerLine = line.strip()
        oldHeaderLine = None
        
        while (headerLine):

          for refDatum in self.refData:
  
            (sparkyCode,numHeaderCol,numDataCol,colLen,refInfo,regExpPatt) = refDatum
  
            allMatches = self.patt[self.format + sparkyCode].findall(headerLine)
  
            # Make sure that it handles the column header info from left to right!
            if allMatches and not headerLine[:headerLine.index(allMatches[0])].strip():
              numMatches = len(allMatches)
  
              #
              # Number of total columns in header and for data
              #
  
              numHeaderCols += numHeaderCol * numMatches
              numDataCols += numDataCol * numMatches
  
              #
              # Keep track of column info
              #
  
              for i in range(numMatches):
                dataColCodes.append(sparkyCode)
                dataColNum.append(numDataCol)
                dataColLengths.append(colLen)
  
              #
              # This one to check data lines later
              #
  
              dataColLen += colLen * numMatches                    
  
              #
              # Set ndim (number of dimensions)
              #
  
              if refInfo == 'ndim':
                ndim = numMatches
                
              for match in allMatches:
                headerLine = headerLine.replace(match,'').strip()

          # Exit if no matches found
          if oldHeaderLine == headerLine:
            headerLine = None
          
          oldHeaderLine = headerLine

        #
        # Exit if there's column names that are not recognized: no use trying
        # to read these for Sparky (very loose format)
        #

        if numHeaderCols != len(cols):
          print "  Unrecognized column headers: can't read sparky file... "
          print dataColCodes
          print dataColNum
          print dataColLengths

          return


        #
        # Initialize peak file
        #

        (path,specName) = splitPath(self.name)
        self.setSpectrumInfo(specName,ndim)
        
        #
        # Check if there are columns that can be empty (e.g. Note)
        #
        
        numPossibleEmptyCols = dataColLengths.count(0)

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
              
            elif numPossibleEmptyCols and (len(cols) == numDataCols - numPossibleEmptyCols):
            
              dataColCodesLen = len(dataColCodes)

            else:

              #
              # If problems, still try to read in essential bits (all until 'ndim' column)
              #

              print "  Warning: Bad column match for following Sparky line. Trying to read...\n%s" % line,

              for i in range(len(dataColCodes)):            

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

              if len(cols) == numDataCols or len(cols) == numDataCols - numPossibleEmptyCols or obligColumnFound == -1:

                values[sparkyCode].append(" ".join(cols[colNum:colNum + dataColNum[i]]))

                if dataColNum[i] > 1:
                  colNum += dataColNum[i] - 1

              else:

                values[sparkyCode].append(line[0:dataColLengths[i]])
                line = line[dataColLengths[i]:]

              colNum += 1

            self.peaks.append(SparkyPeak(self,values))          

          line = fin.readline()

      line = fin.readline()

    fin.close()



  def write(self, verbose = 0):

    if verbose == 1:
      print "Writing Sparky peak list %s" % self.name

    fout = open(self.name,'w')

    #
    # Write out header
    #

    fout.write("      Assignment")

    #
    # Write out dimension codes
    #

    for dim in range(self.numDims[0]):
      fout.write("         w%d" % (dim + 1))

    self.volumeIntensityType.sort() # Make sure Height comes first
    for volType in self.volumeIntensityType:
      fout.write("       %s" % volType.capitalize())

    fout.write(self.newline)

    #
    # Write out peaks - this is the 'Python export' format
    #  

    for peak in self.peaks:
    
      line = "%17s" % peak.getAssignmentText()

      for dim in range(self.numDims[0]):

        line += "%11.3f"  % peak.ppm[dim]

      line += ' '
      
      for volType in self.volumeIntensityType:
        volOrInt = peak.getVolOrInt(volType)
        
        if volOrInt == None:
          volOrInt = 0.0
        
        if volType == 'height':
          line += "%10.2E " % (volOrInt)
        else:  
          line += "%10.2E %2s" % (volOrInt,peak.getVolumeMethodCode())
      
      fout.write(line)
      fout.write(self.newline)

    fout.close()

class SparkyPeak:

  def __init__(self,parent,values):

    self.num = None
    self.intensityFit = None
    self.intensity = None
    self.volume = None
    self.volumeMethod = None
    self.ppm = []
    self.hz = []
    self.assign = []
    self.widths = []
    
    self.parent = parent

    for valueKey in values.keys():

      self.setValue(valueKey,values[valueKey])

  def setValue(self,valueKey,values):
    
    if len(values) == 1 and type(values[0]) == type(''):

      value = values[0]
      
      if not self.parent.patt['emptyline'].search(value):
      
        if valueKey == 'num':
          self.num = returnInt(value)

        elif valueKey == 'assignment':
          assignment = value.strip()
          self.assign = assignment.split(self.parent.assignSep)

        elif self.intensity == None and valueKey in ('dataHeight','height'):
          self.intensity = returnFloat(value)

        elif valueKey == 'fitHeight':
          self.intensityFit = returnFloat(value)

        elif valueKey == 'intensity':
          self.intensity = returnFloat(value)
   
        elif valueKey == 'volume':
          (volume,volumeMethod) = value.strip().split()
          self.volume = returnFloat(volume)

          if volumeMethod == 'ga':
            self.volumeMethod = 'Gaussian'
          elif volumeMethod == 'lo':
            self.volumeMethod = 'Lorentzian'
          elif volumeMethod == 'bx':
            self.volumeMethod = 'Sum over box'
          elif volumeMethod == 'el':
            self.volumeMethod = 'Sum over ellipse'
          elif volumeMethod == '--':
            self.volumeMethod = 'Unknown method'
          else:
            self.volumeMethod = 'Unknown method'
            print "  Error: unknown volume method %s in sparky peaks file (please define)." % volumeMethod

    elif len(values) >= 1:
    
      if valueKey == 'freqDimsPpm':
        self.ppm = returnFloats(values)

      elif valueKey == 'freqDimsHz':
        self.hz = returnFloats(values)

      elif valueKey == 'freqDimsWidth':
        self.widths = returnFloats(values[:-1])
        
    else:
    
      print "  Warning: empty data field for Sparky peak file."
      pass

  def getVolumeMethodCode(self):
  
    if self.volumeMethod == 'Gaussian':
      return 'ga'
    elif self.volumeMethod == 'Lorentzian':
      return 'lo'
    elif self.volumeMethod == 'Sum over box':
      return 'bx'
    elif self.volumeMethod == 'Sum over ellipse':
      return 'el'
    elif self.volumeMethod == 'Unknown method':
      return '--'
    else:
      return None
  
  def getVolOrInt(self,volType):
    
    volOrInt = None
    if volType == 'height':
      if self.intensity:
        volOrInt = self.intensity
      elif self.intensityFit:
        volOrInt = self.intensityFit
    else:
      volOrInt = self.volume
  
    return volOrInt
  
  def getAssignmentText(self):
  
    assignText = []

    for assign in self.assign:
    
      if not assign:
        assignText.append("?")
      else:
        assignText.append(assign)
    
    return self.parent.assignSep.join(assignText)
