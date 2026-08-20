"""
======================COPYRIGHT/LICENSE START==========================

peaksIO.py: I/O for nmrView peak list file

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
from ccp.format.nmrView.generalIO import NmrViewGenericFile
from memops.universal.Io import getTopDirectory
from memops.universal.Io import splitPath

#####################
# Class definitions #
#####################

            
class NmrViewPeakFile(NmrViewGenericFile):
  
  """
  Information on file level
  """
  
  def initialize(self,assignTagSep = ' '):
  
    self.peaks = []
    self.sws = []
    self.sfs = []
    
    self.assignTagSep = assignTagSep  # This is a ' ' in newer versions...
    self.hasAssignItems = True

  def setSpectrumInfo(self,specName,ndim):
  
    self.specNames = [specName]
    self.numDims = [ndim]

    self.dimCodes = ndim * ['']

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s peak list %s" % (self.format,self.name)

    fin = open(self.name, 'rU')

    # Read the atom specs on second line
    line = fin.readline()
    line = fin.readline()
    atomlist = line.split()

    # Initialize
    # ----------

    n = len(atomlist)

    (path,specName) = splitPath(self.name)
    self.setSpectrumInfo(specName,n)

    assign = n*['']
    ppm = n*[0]
    width = n*[0]
    box = n*[0]
    o1 = n*[0]
    shape = n*['']

    # TODO: check if dimCodes, sw match experiment when using data from testProject!
    # Do this in makeTestProject??!  

    for i in range(0,n):
      self.dimCodes[i] = string.strip(atomlist[i])

    line = fin.readline()

    # Read rest file
    while line:

      # Line with peak info
      if self.patt[self.format + 'DigitSpace'].search(line):

        # Decompose line
        cols = line.split()
        peaknum = cols.pop(0)

        # This is format for v3
        if self.patt[self.format + 'CurlyBrace'].search(line):
          num_void = 2 # TODO: one of these is peak comment!

          # Assignments between {} can be separated by columns: rejoin
          # Note that will mash together any specific assignment bits - too complicated to process otherwise
          i = 0
          while i < len(cols):
            col = cols[i]

            # This is an assignment bit
            if self.patt[self.format + 'CurlyBraceStart'].search(col):

              # Join with next column if no closing bracket
              if not self.patt[self.format + 'CurlyBrace'].search(col):
                cols[i] = string.join(col + ' ' + cols.pop(i+1),'')
                continue

              # If full assignment, get rid of brackets
              else:
                cols[i] = self.patt[self.format + 'CurlyBrace'].search(col).group(1)

            i = i + 1

        # This is format for v2
        else:
          num_void = 1

        # Get info for each dimension
        assignItems = -1
        for i in range(0,n):
          (assign[i],ppm[i],width[i],box[i],shape[i]) = cols[:5]
          
          # Get rid of spaces at beginning and end, can mess up code below
          assign[i] = assign[i].strip()
          
          # Check whether assignments are all same length...
          if assignItems == -1:
            assignItems = len(assign[i].split(self.assignTagSep))
          else:
            if len(assign[i].split(self.assignTagSep)) != assignItems:
              self.hasAssignItems = False
          
          del cols[0:(5+num_void)]

        # Get remaining info
        (volume,intensity,status)=cols[0:3]

        # Only collect non-deleted peaks
        if status != -1:
          peak = NmrViewPeak(peaknum,assign,ppm,width,box,shape,volume,intensity,status)
          
          if not None in (peak.num,peak.status):
            self.peaks.append(peak)

      # Line with spectrum info
      elif self.patt[self.format + 'CurlyBrace'].search(line) or self.patt[self.format + 'NumbersNoBrace'].search(line):
        # Keeps reported sw in memory for validation purposes
        self.sws = self.getValues(line)
        line = fin.readline()
        # Keeps o1 of spectrum in MHz in memory.for ppm-Hz conversion
        # Also for validation purposes
        self.sfs = self.getValues(line)

        # TODO: check this O1 against value from API!

      line = fin.readline()


  def write(self, verbose = 0):

    if verbose == 1:
      print "Writing nmrView peak list %s" % self.name

    fout = open(self.name,'w')

    numDim = len(self.dimCodes)

    #
    # Write out header
    #

    fout.write("label dataset sw sf" + self.newline)

    #
    # Write out dimension codes
    #

    for dim in range(0,numDim):
      fout.write("%s " % self.dimCodes[dim])
    fout.write(self.newline)

    #
    # Write out spectrum name
    #

    fout.write("%s" % (self.specNames[0]) + self.newline)

    #
    # Write out spectral width and o1 information
    #

    swsline = ""
    sfsline = ""

    for dim in range(0,numDim):
      swsline += "{%.2f } " % self.sws[dim]
      sfsline += "{%8.4f } " % self.sfs[dim]

    swsline = string.rstrip(swsline)
    sfsline = string.rstrip(sfsline)

    fout.write(swsline + self.newline)
    fout.write(sfsline + self.newline)

    #
    # Write out column labels (version 3 or higher!)
    #
    # Note that there is no label for the peak number
    #

    for dim in range(0,numDim):
      for label in ['L','P','W','B','E','J','U']:
        fout.write(" %s " % (self.dimCodes[dim] + '.' + label))

    for label in ['vol','int','stat','comment','flag0']:
      fout.write(" %s " % label)

    fout.write(self.newline)

    #
    # Write out peaks
    #  
    # Format v3 or higher for now (TODO: make options for v2 later on?)
    #

    for peak in self.peaks:
      line  = "%d "      % peak.num

      for dim in range(0,numDim):

        line += " {%s} "  % peak.assign[dim]
        line += "  %5.3f "  % peak.ppm[dim]
        line += "  %5.3f "  % peak.width[dim]
        line += "  %5.3f "  % peak.box[dim]
        line += "  %s "    % peak.shape[dim]
        line += "  %5.3f "  % 0.00
        line += "  {%s}  "   % "?"


      line += "%.5f "     % peak.volume
      line += "%.5f "     % peak.intensity
      line += "%d "       % peak.status
      line += "{%s} "     % "?"
      line += "%d"        % 0
      line += self.newline

      fout.write(line)

    fout.close()


  def getValues(self,line):

    vals = []
    line = string.replace(line,'{','')
    line = string.replace(line,'}','')
    vallist = line.split()
    for i in range(0,len(vallist)):
      vals.append(returnFloat(vallist[i]))
    return vals

class NmrViewPeak:

  def __init__(self,num,assign,ppm,width,box,shape,volume,intensity,status,isReading = True):
  
    self.num = returnInt(num,default=None)
    self.volume = self.returnFloat(volume)
    self.intensity = self.returnFloat(intensity)
    self.status = returnInt(status,default=None)
    self.shape = shape[:]

    self.assign = []
    self.ppm = []
    self.width = []  # In ppm!!
    self.box = []    # In ppm!!

    for i in range(0,len(assign)):
      self.ppm.append(self.returnFloat(ppm[i]))
      self.width.append(self.returnFloat(width[i]))
      self.box.append(self.returnFloat(box[i]))
      
      # Get rid of ?, when reading
      if isReading and assign[i] == '?':
        assign[i] = None
        
      self.assign.append(assign[i])

  def returnFloat(self,value):
    
    if value == 'nan':
      value = None
    else:
      value = returnFloat(value)
    
    return value
  
