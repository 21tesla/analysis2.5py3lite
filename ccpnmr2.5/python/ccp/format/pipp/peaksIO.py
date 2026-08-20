#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

peaksIO.py: I/O for Pipp peak list files

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

import os, string

# Import general functions
from memops.universal.Util import returnFloat, returnInt
from ccp.format.pipp.generalIO import PippGenericFile
from memops.universal.Io import getTopDirectory
from memops.universal.Io import splitPath

#####################
# Class definitions #
#####################

      
class PippPeakFile(PippGenericFile):
  """
  Information on file level
  """
  def initialize(self):
  
    self.peaks = []
    
    self.assignSep = '.'
    self.assignBreak = ';'
    self.assignInternalBreak = ','
    self.notAssigned = '****'
    
    self.npoints = []
    self.swhs = []
    self.o1s = []
    self.sfs = []
    self.units = []
    
    self.mapping = {}
    self.finalMapping = []

  def setSpectrumInfo(self,specName,ndim):

    self.specNames = [specName]
    self.numDims = [ndim]
    self.numDim = ndim

  def setSpectrumData(self,dataCols):
    self.npoints.append(returnInt(dataCols[0]))
    self.swhs.append(returnFloat(dataCols[1]))
    self.sfs.append(returnFloat(dataCols[2]))
    self.o1s.append(returnFloat(dataCols[3]))  # TODO: is this correct?
    self.units.append(dataCols[4])
    #TODO: unknown data in dataCols[5]

  def read(self,verbose = 0):

    dimTranslate = []

    if verbose == 1:
      print "Reading %s peak list %s" % (self.format,self.name)

    fin = open(self.name, 'rU')

    # Read, look for lines starting with DATA, FORMAT and VARS
    line = fin.readline()

    while line:
      cols = line.split()

      if len(cols) == 0 or self.patt['hash'].search(line):
        pass

      elif cols[0] == 'DATA':

        if cols[1] == 'DIMCOUNT':

          ndim = returnInt(cols[2])

          (path,specName) = splitPath(self.name)
          self.setSpectrumInfo(specName,ndim)

        elif self.patt[self.format + '?_AXIS'].search(cols[1]):
          dimCode = self.patt[self.format + '?_AXIS'].search(cols[1]).group(1)
          dimTranslate.append(dimCode)

          self.setSpectrumData(cols[2:])

      elif cols[0] == 'FORMAT':
        # Info here not really necessary as it is a space-delimited format
        pass

      elif cols[0] == 'VARS':

        defCols = {}
        assignKeywords = []
        dimChemElement = []
        
        addNDims = 0

        # Get column headers
        for i in range(1,len(cols)):
        
          colId = cols[i]

          if cols[i][:6] == 'Assign':
            if not assignKeywords:
              # Special case - hardcode for assignments
              for dim in range(ndim):
                colId = 'Assign%d' % (dim + 1)
                assignKeywords.append(colId)
                defCols[colId] = i-1 + dim
              addNDims = ndim - 1
            else:
              # Means that more than one column defined, so 'none' assign will show up...
              addNDims -= 1
              
          else:
            if len(cols[i]) <= 2:
              axisType = cols[i][-1].upper()
              chemElement = None
              #TODO: what does pipp use for 4Ds??!?
              if axisType in ['X','Y','Z']:
                colId = axisType
                if len(cols[i]) == 2:
                  chemElement = cols[i][0].upper()
                dimChemElement.append(chemElement)

            defCols[colId] = i-1

        defColsLen = len(defCols)

        # Start reading data
        line = fin.readline()

        while line:
          cols = line.split()
          localAddNDims = addNDims

          if len(cols) == 0 or self.patt['hash'].search(line):
            pass
            
          colsCopy = cols[:]
          colIndex = 0
          for i in range(len(colsCopy)):
            col = colsCopy[i]
            if col.count(self.assignInternalBreak):
              (astring1,astring2) = col.split(self.assignInternalBreak)
              (residue,atomName) = astring1.split(self.assignSep)
              
              cols[colIndex] = astring1
              cols.insert(colIndex+1,residue + self.assignSep + astring2)
              colIndex += 2
              
            elif col == self.notAssigned :
              for i in range(localAddNDims):
                cols.append(self.notAssigned)
              colIndex += (localAddNDims + 1)
              localAddNDims = 0

            else:
              colIndex += 1

          if len(cols) != defColsLen:
            print "  Warning: ignoring following line in %s peak file (columns don't match):\n%s" % (self.format,line)
            pass

          else:
            
            for peakIdName in ('PkID','Pk-ID'):
              if defCols.has_key(peakIdName):
                Id = cols[defCols[peakIdName]]
                break
              
            intensity = cols[defCols['Intensity']]

            # Get the assignment strings
            
            assign = []

            for defColKwd in assignKeywords:
              
              astring = cols[defCols[defColKwd]].strip()
              
              if astring[-1] == self.assignBreak:
                astring = astring[:-1]
              
              if astring == self.notAssigned :
                assign.append(None)
              
              else:
                assign.append(astring)
            
            for i in range(len(assign),len(dimTranslate)):
              assign.append(None)
            
            ppm = []
            dimSlice = []

            for dimCode in dimTranslate:
              ppm.append(returnFloat(cols[defCols[dimCode]]))

              if defCols.has_key('Sl.' + dimCode):
                dimSlice.append(cols[defCols['Sl.' + dimCode]])
              else:
                dimSlice.append(None)

            # Have to do some (rough) checking to determine the mapping between
            # the assignment and the columns in the file...
            
            self.determineMapping(assign,ppm)
            
            self.peaks.append(PippPeak(Id,intensity,assign,ppm,dimSlice = dimSlice))

          line = fin.readline()

      line = fin.readline()
    
    fin.close()

    #
    # Determine best mappings... gives the mapping between the ASSIGNMENT dim (column
    # position) to the PPM dim in the peak list (column value).
    #
    
    for aInd in range(0,self.numDim):
    
      if self.finalMapping == None:
        continue
        
      prefPpmInd = (-1,0)
      for ppmInd in range(0,self.numDim):
        mappingIndex = (aInd,ppmInd)
        if self.mapping.has_key(mappingIndex) and self.mapping[mappingIndex] > prefPpmInd[1]:
          prefPpmInd = (ppmInd,self.mapping[mappingIndex])
      
      if prefPpmInd[0] not in self.finalMapping:
        self.finalMapping.append(prefPpmInd[0])
      else:
        print "  Error: mapping for assignment could not be determined!"
        self.finalMapping = None
    
    #
    # Rearrange assignments for all peaks...
    #
    
    if self.finalMapping:
      
      for peak in self.peaks:      
        assign = peak.assign[:]
        for aInd in range(0,len(self.finalMapping)):
          ppmInd = self.finalMapping[aInd]
          peak.assign[ppmInd] = assign[aInd]
    
  def determineMapping(self,assign,ppm):

    for aInd in range(0,self.numDim):
      if assign[aInd]:
        (residue,atomName) = assign[aInd].split(self.assignSep)
        
        for ppmInd in range(0,self.numDim):
          value = ppm[ppmInd]
          mappingIndex = (aInd,ppmInd)
          
          if not self.mapping.has_key(mappingIndex):
            self.mapping[mappingIndex] = 0

          if atomName == 'HN' and value > 7.0 and value < 10.0:
            self.mapping[mappingIndex] += 1
          elif atomName[:2] == 'HA' and value > 4.0 and value < 5.0:
            self.mapping[mappingIndex] += 1
          elif atomName[:2] == 'HB' and value > 2.0 and value < 4.0:
            self.mapping[mappingIndex] += 1
          elif atomName[0] == 'C' and value > 70.0 and value < 100.0:
            self.mapping[mappingIndex] += 1
          elif atomName[0] == 'N' and value > 100.0 and value < 140.0:
            self.mapping[mappingIndex] += 1
          elif atomName[0] == 'C' and value > 70.0:
            self.mapping[mappingIndex] += 1

class PippPeak:

  def __init__(self,Id,intensity,assign,ppm,dimSlice = []):

    self.num = returnInt(Id)
    self.intensity = returnFloat(intensity)
    self.assign = assign[:]
    self.dimSlice = dimSlice[:]

    self.ppm = []

    for i in range(0,len(ppm)):
      self.ppm.append(returnFloat(ppm[i]))
        
###################
# Main of program #
###################

if __name__ == "__main__":  
                                                      
  files = ['../../reference/ccpNmr/michaelNewYork/050309/fhsqc_stapp.PCK',
           '../../reference/pipp/3d.1.peaks',
           '../../reference/pipp/3d.2.peaks',
           '../../reference/ccpNmr/michaelNewYork/050323/n15noehsqc.PCK',
           '../../reference/ccpNmr/michaelNewYork/050425/motafl6dipsisecryo.PCK',
           '../../reference/ccpNmr/michaelNewYork/050425/motafl6hbhaconhsecryophos.PCK',
          ]

  
  for inFile in files:
    
    peakFile = PippPeakFile(os.path.join(getTopDirectory(),inFile))

    peakFile.read(verbose = 1)
    
    for rawPeak in peakFile.peaks:
      print rawPeak.num, rawPeak.assign, rawPeak.intensity
