"""
======================COPYRIGHT/LICENSE START==========================

peaksIO.py: I/O for XEasy peak list files

Copyright (C) 2006-2008 Wim Vranken (European Bioinformatics Institute)

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

# Import general functions
from memops.universal.Util import returnInt, returnFloat
from memops.universal.Io import splitPath

from ccp.format.xeasy.generalIO import XEasyGenericFile

#####################
# Class definitions #
#####################

class XEasyPeakFile(XEasyGenericFile):

  def initialize(self):
  
    self.peaks = []
    
    self.cyanaFormat = False

  def setSpectrumInfo(self,specName,ndim):
  
    self.specNames = [specName]
    self.numDims = [ndim]

    self.dimCodes = ndim * ['']

  def translateIntMethod(self,source,method):
  
    translateIntMethods = [("d","Denk integration"),
                           ("r","Rectangular integration"),
	       ("e","Elliptical integration"),
	       ("m","Maximum integration"),
	       ("a","Automatic integration"),
	       ("-","Not integrated")]
    
    if source == 'ccpn':
      for (code,methodName) in translateIntMethods:
        if method == methodName:
          return code
      return 'a'        # is default
    
    elif source == 'xeasy':
      for (code,methodName) in translateIntMethods:
        if method == code:
          return methodName
      print "Unknown xeasy integration code %s" % method

    else:
      print "Unknown source for xeasy integration method translation %s" % source

    return None
    
  def read(self, verbose=0, cyanaFormat=False, cyanaTextAssignments=False):
    
    #
    # Initialisation to deal with CYANA flavours of this format
    #

    addText = ""
    self.cyanaFormat = cyanaFormat
    if cyanaFormat:
      addText += ' as CYANA format'
    self.cyanaTextAssignments = cyanaTextAssignments
    if cyanaTextAssignments:
      addText += ' with assignments as text strings'

    if verbose == 1:
      print "Reading %s peak list %s%s." % (self.format,self.name,addText)
      
    #
    # Now read file..
    #

    fin = open(self.name, 'rU')

    # Try to read number of dimensions on first line (default XEASY format)
    line = fin.readline()
    cols = line.split()

    n = returnInt(cols[-1])
    
    # If this didn't work, try something else, could be of #NnH form
    if not n:
      if line[0] == '#':
        dimCodes = line[1:].strip()
        n = len(dimCodes)
        print "  Warning: irregular format, extracting dimensions from %s string on first first line." % dimCodes
    
    if not n:
      print "  Error: number of dimensions could not be determined, aborting file read."
      return False        
    
    (path,specName) = splitPath(self.name)
    self.setSpectrumInfo(specName,n)

    assign = n*['']
    ppm = n*[0]
    
    peaknum = None
    
    # CYANA 3 files with text assignments have a column less at the end.
    if cyanaTextAssignments:
      numFinalColumns = 0
    else:
      numFinalColumns = 1
      
    
    line = fin.readline()
    # Read rest file
    while line:

      if not self.patt['hash'].search(line) and not self.patt['emptyline'].search(line):
        
        searchObj = self.patt[self.format + 'PeakInfo'].search(line)
       
        peakLine = searchObj.group(1)
        peak_el = peakLine.split()
        
        if len(peak_el) > n and (peak_el[0] != '0' or peak_el[1] != '0.000') and (not peaknum or peaknum != peak_el[0]):
        
          peaknum = peak_el.pop(0)
          
          for i in range(0,n):
            ppm[i] = peak_el.pop(0)       

          # This is a minimal format definition that came up, is handled now (Wim 18/03/2010)
          (colour,userCode,volume) = peak_el[0:3]
          del peak_el[0:3]
          
          if peak_el:
            (volumeError,intCode) = peak_el[0:2]        
            del peak_el[0:3] # deleting extra void code
          else:
            volumeError = '0.0'
            intCode = '-'

          for i in range(n):
            if peak_el:
              assign[i] = peak_el.pop(0)
            else:
              assign[i] = '0'
        
          #
          # Create the peak and check if valid
          #
          
          comment = line.replace(peakLine,'').strip()
          
          peak = XEasyPeak(peaknum,ppm,colour,userCode,volume,volumeError,intCode,assign,verbose=False,comment=comment,textAssign=self.cyanaTextAssignments)
          
          if not None in peak.ppm + peak.assign + [peak.num,peak.volume,peak.volumeError,peak.colour]:
            self.peaks.append(peak)

        else:
          
          if len(peak_el) == n:
            start = 0
          else:
            start = len(peak_el) - n - numFinalColumns
          
          for i in range(n):
            assign[i] = peak_el[i+start]
            
          self.peaks[-1].addAssignment(assign)

      elif self.patt[self.format + 'IName'].search(line):
        searchObj = self.patt[self.format + 'IName'].search(line)
        dim = returnInt(searchObj.group(1))
        dimCode = searchObj.group(2)
        self.dimCodes[dim-1] = dimCode
        
      elif self.patt[self.format + 'CyanaFormat'].search(line):
        self.cyanaFormat = True

      line = fin.readline()

    fin.close()

  def write(self, verbose=0, cyanaFormat=False, cyanaTextAssignments=False):

    if verbose == 1:
      print "Writing xeasy peak list %s" % self.name

    fout = open(self.name,'w')

    numDim = len(self.dimCodes)

    #
    # Write out header
    #

    fout.write("# Number of dimensions %d" % numDim)
    fout.write(self.newline)
    
    # Changed to use HhN for Cyana and HNHN otherwise
    dimCodes = list(self.dimCodes)
    if not cyanaFormat:
      # Added to allow reading by UNIO
      if set(dimCodes) == set(('h', 'H', 'N')):
        dimCodes[dimCodes.index('H')] = 'HN'
        dimCodes[dimCodes.index('h')] = 'H'
      elif set(dimCodes) == set(('h', 'H', 'C')):
        dimCodes[dimCodes.index('H')] = 'HC'
        dimCodes[dimCodes.index('h')] = 'H'
    
    for dim, dimCode in enumerate(dimCodes):
      fout.write("#INAME %d %s" % (dim+1,dimCode) + self.newline)
    
    if cyanaFormat:
      # Added for Cyana reading
      fout.write("#CYANAFORMAT %s" % (''.join(dimCodes)) + self.newline) 
       
    #dimCodeString = ""
    #for dim in range(0,numDim):
      #dimCode = self.dimCodes[dim]
      # RAsmus Fogh 15/6/12. Removed, as problem treated in XEasyFormat
      #if dimCode in dimCodeString:
      #  dimCode = dimCode.lower()
      

    #
    # Write out peaks
    #
    
    preAssignLen = 44 + 8 * numDim

    for peak in self.peaks:
    
      lineStart = "%4d "     % peak.num

      for dim in range(0,numDim):
        lineStart += "%7.3f " % peak.ppm[dim]

      lineStart += "%d "      % peak.colour
      lineStart += "%s "      % peak.userCode
      lineStart += "%18.3e "  % peak.volume
      lineStart += "%9.2e "   % peak.volumeError
      lineStart += "%s "      % peak.intCode
      lineStart += "%3d "     % 0
      
      line = lineStart

      for dim in range(0,numDim):
        if not cyanaTextAssignments:
          line += "%4d " % peak.assign[dim]
        else:
          line += "%-9s " % peak.assign[dim]

      if not cyanaTextAssignments:
        line += "%d"       % 0
      line += self.newline
      
      if cyanaFormat and peak.ambiguousAssign:
        for assign in peak.ambiguousAssign:
          
          if not cyanaTextAssignments:
            line += preAssignLen * ' '
          else:
            line += lineStart
            
          for dim in range(0,numDim):
            if not cyanaTextAssignments:
              line += "%4d " % assign[dim]
            else:
              line += "%-9s " % assign[dim]
          line += self.newline

      fout.write(line)

    fout.close()

#
# num  w1  w2  (w3)  colour user_defined_spec_type  volume  volume_error(%) int_method void_num  assiw1 assiw2 (assiw3) void_num

class XEasyPeak:

  def __init__(self,num,ppm,colour,userCode,volume,volumeError,intCode,assign,verbose=False,comment=None,textAssign=False):

    self.num = returnInt(num,default=None,verbose=verbose)
    self.volume = returnFloat(volume,default=None,verbose=verbose)
    self.volumeError = returnFloat(volumeError,default=None,verbose=verbose)
    self.colour = returnInt(colour,default=None,verbose=verbose)
    self.userCode = userCode
    self.intCode = intCode
    
    self.comment=comment
    self.textAssign = textAssign

    self.ppm = []
    self.assign = []
    self.ambiguousAssign = []

    for i in range(len(assign)):
      self.ppm.append(returnFloat(ppm[i],default=None,verbose=verbose))
      
      self.assign.append(self.getAssignmentValue(assign[i]))

  def addAssignment(self,assign):
    
    self.ambiguousAssign.append([])

    for i in range(len(assign)):
      self.ambiguousAssign[-1].append(self.getAssignmentValue(assign[i]))
      
  def getAssignmentValue(self,assignInitValue):

    if not self.textAssign:
      assignValue = returnInt(assignInitValue,default=None)
    else:
      assignValue = assignInitValue
  
    return assignValue
