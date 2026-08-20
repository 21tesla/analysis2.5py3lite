#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

projectIO.py: I/O for Ansig project file

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
from memops.universal.Io import getTopDirectory
from memops.universal.Util import returnFloat, returnInt
from ccp.format.ansig.generalIO import AnsigGenericFile

#####################
# Class definitions #
#####################

class AnsigProjectFile(AnsigGenericFile):

  # Information on file level

  def initialize(self):
    
    (self.directory,self.fileName) = os.path.split(self.name)

    self.directories = {}
    self.files = {}
    
    self.pars = {}
    
    self.spectraFile = None

  def read(self,readSpectra = True, verbose = 0):

    if verbose == 1:
      print "Reading %s project file %s" % (self.format,self.name)
      
    #
    # Initialize
    #
    
    inControlBlock = 0
    
    #
    # Read file...
    #

    fin = open(self.name, 'rU')
    line = fin.readline()

    while line:

      # Ignore empty lines and comments
      if self.patt['emptyline'].search(line) or self.patt['exclamation'].search(line):
        line = fin.readline()
        continue
        
      cols = string.split(line)
      
      if cols[0] == 'control':
        inControlBlock = 1
      elif cols[0] == 'end_control':
        inControlBlock = 0
      elif inControlBlock:
        name = cols[0]
        value = cols[1]
        
        # Substitute directories in value
        searchObj = self.patt[self.format + 'CurlyBrace'].search(value)

        if searchObj:
          directoryName = searchObj.group(1)
          value = string.replace(value,searchObj.group(),self.directories[directoryName])
      
        # Check if directory...
        searchObj = self.patt[self.format + 'CurlyBrace'].search(name)
        
        if searchObj:
          self.directories[searchObj.group(1)] = value
        elif name not in ['binaryformat']:
          self.files[name] = os.path.join(self.directory,value)
        else:
          self.pars[name] = value
          
      line = fin.readline()

    fin.close()
    
    if readSpectra and self.files.has_key('spectra'):
      self.spectraFile = AnsigSpectraFile(self.files['spectra'],self)
      self.spectraFile.read(verbose = verbose)

    return

class AnsigSpectraFile(AnsigGenericFile):

  # Information on file level

  def initialize(self,parentFile):

    self.spectra = []
    self.ppmRange = {}
    
    self.binaryFormat = None
    if parentFile.pars.has_key('binaryformat'):
      self.binaryFormat = parentFile.pars['binaryformat']
    
    self.parentFile = parentFile

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s spectra file %s" % (self.format,self.name)
      
    #
    # Initialize
    #
    
    inSpectrum = 0
    inMatrix = 0
    
    #
    # Read file...
    #

    fin = open(self.name, 'rU')
    line = fin.readline()

    while line:

      # Ignore empty lines and comments
      if self.patt['emptyline'].search(line) or self.patt['exclamation'].search(line):
        line = fin.readline()
        continue
        
      cols = string.split(line)
      
      if cols[0] == 'spectrum':
        inSpectrum = 1
        
        specInfo = {'name': cols[2],'numDim': returnInt(cols[1][:-1]), 'dims': [], 'dimInfo': {}}
        
        if len(cols) > 3 and cols[3] == 'expnucleus':
          specInfo['dims'].append(cols[4])
          specInfo['dimInfo'][cols[4]] = {'nucleus': cols[5]}
        
      elif cols[0] == 'end_spectrum':
        inSpectrum = 0
        self.spectra.append(specInfo)
                
      elif cols[0] == 'ppm_range':
        self.ppmRange[cols[1]] = (returnFloat(cols[2]),returnFloat(cols[3]))

      elif inSpectrum:

        if cols[0] == 'expnucleus':
          specInfo['dims'].append(cols[1])
          specInfo['dimInfo'][cols[1]] = {'nucleus': cols[2]}

        elif cols[0] == 'experiment':
          specInfo['type'] = cols[1]

        elif cols[0] == 'solvent':
          specInfo['solvent'] = cols[1]
          specInfo['temperature'] = cols[3]
          specInfo['pH'] = cols[5]

        elif cols[0] == 'equivalent':
          specInfo['equivalent'] = cols[1]

        elif cols[0] == 'matrix_file':
          inMatrix = 1
          
          # Assuming it's 1H if no info was given...
          if len(specInfo['dimInfo']) != specInfo['numDim']:
            for dim in range(1,specInfo['numDim'] + 1):
              dimKey = "F%d" % dim
              if not specInfo['dimInfo'].has_key(dimKey):
                specInfo['dimInfo'][dimKey] = {'nucleus': '1H'}

        elif cols[0] == 'end_matrix_file':
          inMatrix = 0
        
        elif inMatrix and cols[0] == 'parameter_file':
          specInfo['parFile'] = self.getFullFile(cols[1])
        
        elif inMatrix and cols[0] == 'file':
          specInfo['specFile'] = self.getFullFile(cols[1])

        elif inMatrix and cols[0] == 'aliased':
          specInfo['dimInfo'][cols[1]]['aliased'] = 1

        elif inMatrix and cols[0] == 'linewidth':
          for i in range(1,len(cols),2):
            specInfo['dimInfo'][cols[i]]['linewidth'] = returnFloat(cols[i+1])
        
        elif inMatrix and cols[0] in ['noise_level','contour_base']:
          pass

        elif cols[0] in ['contour_file','scale'] or cols[0].count('{contours}') > 0:
          pass
         
        else:
          print "Warning: did not handle inside tag %s" % line
       
      else:
        print "Warning: did not handle %s" % line
        
      line = fin.readline()

    fin.close()
    
    return
  
  def getFullFile(self,fileName):
  
    # Substitute directories in value
    searchObj = self.patt[self.format + 'CurlyBrace'].search(fileName)

    if searchObj:
      directoryName = searchObj.group(1)
      fileName = string.replace(fileName,searchObj.group(),self.parentFile.directories[directoryName])
    
    return fileName
    
###################
# Main of program #
###################

if __name__ == "__main__":
                                                      
  files = ['../../reference/ccpNmr/katCam/oxime/oxime.ctr',
           '../../reference/ccpNmr/tim/forWim/g1.ctr',
           '../../reference/ccpNmr/krysEdinburgh/050121/T1_ole_ABA-1A.ctr']

  
  for inFile in files:
    
    projectFile = AnsigProjectFile(os.path.join(getTopDirectory(),inFile))

    projectFile.read(verbose = 1)
    
    print projectFile.directories
    print projectFile.files
    print projectFile.spectraFile.spectra
    print projectFile.spectraFile.binaryFormat
