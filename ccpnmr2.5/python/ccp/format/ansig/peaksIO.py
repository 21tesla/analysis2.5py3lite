#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

peaksIO.py: I/O for Ansig export and storage peak files

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

import os, re, string

# Import general stuff
from memops.universal.Util import returnFloat, returnInt
from memops.universal.Util import returnFloats, returnInts

from ccp.format.ansig.generalIO import AnsigGenericFile
from memops.universal.Io import getTopDirectory

# Note that this script can handle the two formats (export and storage)
# that ansig can use for writing out peak lists

#####################
# Class definitions #
#####################

#
# w1  w2  (w3) volume spectrum ? ? ? ? ? ? ? ? ?  assiw1_resn assiw2_resn (assiw3_resn) assiw1_res assiw2_res (assiw3_res) assiw1_atom assiw2_atom (assiw3_atom)

class AnsigPeakFile(AnsigGenericFile):

  def initialize(self):
  
    self.peaks = []
    self.specNames = []
    self.numDims = []

  def setSpectrumInfo(self,fileType,npeaks,ndim, specName = None):
  
    self.type = fileType
    self.npeaks = returnInt(npeaks)
    
    #
    # If ndim set for export with multiple spectra, then ndim CAN be None:
    # this happens when not clear how many dimensions spectrum has...
    #
    
    if specName:
      specName = string.strip(specName)
      self.numDims = [ndim]
      self.specNames = [specName]
      self.maxDim = ndim
      
    else:
      self.maxDim = ndim

  def setSpectrum(self,specName,ppm):
  
    specName = string.strip(specName)

    # If spectrum name already exists, see if number of dims matches
    if self.specNames.count(specName) == 1:
      specNameIndex = self.specNames.index(specName)
      numDim = self.numDims[specNameIndex]
      
      # Just to make sure not fluke peak with ppm 0.00000 as first one
      # for a particular spectrum, check again
      # TODO: this slows things down... should only do this ten times or so
      
      ndim = self.maxDim
      
      for ppmVal in ppm:
        if returnFloat(ppmVal) == 0.0:
          ndim -= 1
          
      if ndim > numDim:
        self.numDims[specNameIndex] = None
    
    # If doesn't exist, create
    else:
      self.specNames.append(specName)
      
      ndim = self.maxDim
      
      for ppmVal in ppm:
        if returnFloat(ppmVal) == 0.0:
          ndim -= 1
         
      self.numDims.append(ndim)
    
  def checkValid(self):
  
    if self.npeaks != len(self.peaks):
      print "    Warning: number of Ansig peaks reported and read does not match!"


  def read(self, verbose = False):
 
    # Format is in fortran fixed column widths!

    fin = open(self.name, 'rU')

    # Read first line
    line = fin.readline()
    fin.close()

    if line.count('export crosspeaks file'):
      self.readExported(verbose)

    elif line.count('storage crosspeaks file'):
      self.readStorage(verbose)

    else:
      print("Error, file %s is not an ANSIG v3.x export or storage crosspeak file" % self.name)
      return

    self.checkValid()

  def readExported(self, verbose = False):

    if verbose:
      print "  Reading %s as export file." % self.name

    # Code mainly taken from Rasmus' readAnsig.py script
    # Format is in fortran fixed column widths!

    fin = open(self.name, 'rU')

    # Read first line
    line = fin.readline()

    # Read second line
    line = fin.readline()

    npeaks = returnInt(line[:6])
    # Note that n is not necessarily the number of dimensions for a spectrum (could be less)!
    # This has to be determined while reading... bit of a mess
    n = returnInt(line[6:12])

    self.setSpectrumInfo('Export',npeaks,n)

    peaknum = 0

    line = fin.readline()

    # Read rest file
    while line:

      # Ignore empty lines
      if self.patt['emptyline'].search(line):
        line = fin.readline()
        continue


      assignRes = n*['']
      assignResNum = n*['']
      assignAtom = n*['']
      dimcons = n*[('','')]
      ppm = n*[0]

      peaknum += 1
      
      badLineInfo = False

      for i in range(0,n):
        ppm[i] = returnFloat(line[0:13],default = None)
        if ppm[i] == None:
          badLineInfo = True
          break
        
        line = line[13:]
        
      if badLineInfo:
        print("Error, file %s is not an ANSIG v3.x export or storage crosspeak file" % self.name)
        return
      
      volume = returnFloat(line[0:13], default = None)
      line = line[13:]

      #specName = line[0:12]
      specName = line[0:12].strip()
      line = line[12:]

      if not specName:
        deleted = 1
      else:
        deleted = 0 
      
      
      # NBNB TBD look at this!
      # this may break something, 
      # but going on with a deleted spectrum makes no sense
      if deleted:
        continue

      # Correspondences and connections
      # TODO: not implemented in class yet!
      symmcon = returnInt(line[0:6])
      line = line[6:]

      for i in range(0,n):
        x1 = returnInt(line[0:6])    
        x2 = returnInt(line[6:12])    
        line = line[12:]

        dimcons[i] = (x1,x2)

      x1 = returnInt(line[0:6])    
      x2 = returnInt(line[6:12])    
      line = line[12:]

      corrcons = (x1,x2)

      # Assignments
      for i in range(0,n):
        assignResNum[i] = line[0:4]
        line = line[4:]

      for i in range(0,n):
        assignRes[i] = line[0:4]
        line = line[4:]

      for i in range(0,n):
        assignAtom[i] = line[0:4]
        line = line[4:]

      # Determine if this is a new spectrum,
      # how many dimensions it has and if this fits with new data
      self.setSpectrum(specName,ppm)

      self.peaks.append(AnsigPeak(peaknum,deleted,ppm,volume,assignResNum,assignRes,assignAtom,n,specName))

      line = fin.readline()

    fin.close()

  def readStorage(self, verbose = False):

    if verbose:
      print "  Reading %s as storage file." % self.name

    # Code mainly taken from Rasmus' readAnsig.py script
    # Format is in fortran fixed column widths!

    fin = open(self.name, 'rU')

    # Read first line
    line = fin.readline()

    # Read second line
    line = fin.readline()

    specName = line[0:12]
    npeaks = returnInt(line[12:18])
    n = returnInt(line[18:24])

    self.setSpectrumInfo('Storage',npeaks,n,specName = specName)

    peaknum = 0
    
    line = fin.readline()

    # Check line length - different formats available
    print len(line)
    if len(line) == 131:
      formatType = 'short'
    elif len(line) == 143:
      formatType = 'long'

    # Read rest file
    while line:

      # Ignore empty lines
      if self.patt['emptyline'].search(line):
        line = fin.readline()
        continue

      assignRes = n*['']
      assignResNum = n*['']
      assignAtom = n*['']
      dimcons = n*[('','')]
      ppm = n*[0]
      otherInfo = [0,0]

      peaknum += 1

      for i in range(0,n):
        ppm[i] = returnFloat(line[0:13])
        line = line[13:]

      volume = returnFloat(line[0:13])
      line = line[13:]

      # Correspondences and connections
      # TODO: not implemented in class yet!
      symmcon = returnInt(line[0:6])
      line = line[6:]

      for i in range(0,n):
        x1 = returnInt(line[0:6])    
        x2 = returnInt(line[6:12])    
        line = line[12:]

        dimcons[i] = (x1,x2)
        
      # Check for long format, sets otherInfo, whatever may be in there...
      if formatType == 'long':
        for i in range(0,2):
          otherInfo[i] = returnInt(line[0:6])
          line = line[6:]

      # Assignments
      for i in range(0,n):
        assignResNum[i] = line[0:4]    
        line = line[4:]

      for i in range(0,n):
        assignRes[i] = line[0:4]    
        line = line[4:]

      for i in range(0,n):
        assignAtom[i] = line[0:4]    
        line = line[4:]

      self.peaks.append(AnsigPeak(peaknum,0,ppm,volume,assignResNum,assignRes,assignAtom,n,specName))

      line = fin.readline()

    fin.close()

  def write(self, verbose = False):

    if verbose:
      print "Writing nmrView peak list %s" % self.name

    # Only writing to exported format supported

    self.writeExported(verbose = verbose)

  def writeExported(self,verbose):

    fout = open(self.name,'w')

    numDim = self.maxDim

    fout.write("ANSIG v3.3 export crosspeaks file" + self.newline)

    fout.write("%6d%6d" % (self.npeaks,numDim) + self.newline)

    for peak in self.peaks:

      for i in range(0,numDim):

        if len(peak.ppm) <= i:
          value = 0.0
        else:
          value = peak.ppm[i]

        fout.write("%13.6E" % value)

      fout.write("%13.6E" % peak.volume)

      fout.write("%-12s" % peak.specName)

      fout.write("%6d" % 0)

      for i in range(0,numDim):
        fout.write("%6d%6d" % (0,0))

      fout.write("%6d%6d" % (0,0))

      for i in range(0,numDim):

        if len(peak.assignResNum) <= i:
          value = ''
        else:
          value = peak.assignResNum[i]

        fout.write("%-4s" % value)

      for i in range(0,numDim):

        if len(peak.assignRes) <= i:
          value = ''
        else:
          value = peak.assignRes[i]

        fout.write("%-4s" % value)

      for i in range(0,numDim):

        if len(peak.assignAtom) <= i:
          value = ''
        else:
          value = peak.assignAtom[i]

        fout.write("%-4s" % value)


      fout.write(self.newline)

    fout.close()

class AnsigPeak:

  def __init__(self,num,deleted,ppm,volume,assignResNum,assignRes,assignAtom,ndim,specName):

    
    self.num = returnInt(num)
    self.deleted = returnInt(deleted)
    self.volume = returnFloat(volume)
    self.specName = (string.strip(specName))[:12]

    self.ppm = returnFloats(ppm)

    self.assignRes = []
    self.assignAtom = []
    self.assignResNum = []

    for i in range(0,ndim):
      self.assignResNum.append(string.strip(assignResNum[i]))
      self.assignRes.append(string.strip(assignRes[i]))
      self.assignAtom.append(string.strip(assignAtom[i]))
