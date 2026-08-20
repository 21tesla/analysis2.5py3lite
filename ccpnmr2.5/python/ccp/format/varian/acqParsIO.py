"""
======================COPYRIGHT/LICENSE START==========================

acqParsIO.py: Varian acquisition parameter file reader

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
# Import general functions
from memops.universal.Util import returnFloat
from memops.universal.Util import returnFloats
from memops.universal.Util import returnInt
from memops.universal.Util import joinUnquote
from memops.universal.Util import unquote

import os

from ccp.format.varian.generalIO import VarianGenericFile
from ccp.format.general.acqParsDefs import FreqDimReference, FidDimReference
from ccp.format.general.acqParsDefs import SpectrometerReference, GenericAcqParData

from memops.universal.Io import splitPath, joinPath

#####################
# Class definitions #
#####################

class VarianAcqParData(GenericAcqParData):  # TODO: subclass necessary?!

  def setGeneric(self):
  
    #
    # This is a hack...
    #
    
    varianGenericFile = VarianGenericFile('void')
    self.format = varianGenericFile.format

  def initialize(self,numDim,parHelp):
  
    self.numDim = numDim
    self.parHelp = parHelp
    
    (self.fileDir,void) = splitPath(self.name)
    
    self.parFiles = {}

  def get(self,verbose = 0):
    
    #############################
    # Varian, read procpar file #
    #############################
    
    self.parFiles['procpar'] = None

    files = os.listdir(self.fileDir)

    for file in files:
      # Look for procpar file
      if file == 'procpar':
      
        fileName = joinPath(self.fileDir,file)
        
        self.parFiles['procpar'] = VarianAcqParFile(fileName,parHelp = self.parHelp)
        self.parFiles['procpar'].read(verbose = verbose)

    self.process()

  def process(self):
    
    #
    # Only one file...
    #

    mainPars = self.parFiles['procpar'].parameters

    if mainPars.has_key('nt'):
      self.numScans = returnInt(mainPars['nt'].values[0])
    else:
      print "No num scans for %s" % (self.fileDir)
      self.numScans = 0

    #
    # Get frequencies for nuclei (TODO: REFLIST IN settings?)
    #
    
    self.refList = [['sfrq','tof','tn'],['dfrq','dof','dn'],['dfrq2','dof2','dn2'],['dfrq3','dof3','dn3'],['dfrq4','dof4','dn4']]
    nullValue = ''

    self.getSpecRefs(mainPars,nullValue,specTypeRef = 'console')
    
    refsNp = ['np','ni','ni2','ni3']
    refsSw = ['sw','sw1','sw2','sw3']

    self.fidDimRefs = []
    self.freqDimRefs = []

    for i in range(0,self.numDim):

      if i == 0:
        
        #
        # Different handling for acquisition dimension!
        # TODO: if data is NOT COMPLEX, take out factor 2!
        # In CCPN data model, npoints is number of COMPLEX points!
        #
        
        numPoints = mainPars[refsNp[i]].values[0] / 2

      else:

        #
        # Data always complex in indirect? Assuming so...
        # Note that ni, ni2, ... HAS A DIFFERENT DEFINITION than np!!! Is number of COMPLEX points!
        #
        
        numPoints = mainPars[refsNp[i]].values[0] 
      
      #
      # TODO: this is assumed for Varian
      #
      
      numPointsValid = numPoints
      
      #
      # Set sw and increment
      # TODO: If data is NOT COMPLEX, divide valppoint by 2!
      # (one half of spectrum is redundant)
      #

      swh = mainPars[refsSw[i]].values[0]
      valuePerPoint = 1 / swh
      
      #
      # Set FID info
      #

      self.fidDimRefs.append(FidDimReference(numPoints,numPointsValid,valuePerPoint,None))
      
      #
      # Set freq domain info
      # PROBLEM with acquisition time: Rasmus uses valppoint * (npointsvalid-1)
      # But this doesn't correspond to value from Varian (is valppoint * npointsvalid)
      # WHY? Contact Varian. Have to VERY SPECIFICALLY define these parameters!!!
      # How does this work for Bruker?
      #
      
      (numPoints,numPointsOrig,valuePerPoint) = self.getFreqDimInfo(valuePerPoint,numPointsValid)
      self.freqDimRefs.append(FreqDimReference(numPoints,numPointsOrig,valuePerPoint,None))

      #
      # Try to get phasing info (often not correct...)
      #

      refsPhase = 2 * ['']
      refsPhase[0] = ['rp','rp1','rp2','rp3']
      refsPhase[1] = ['lp','lp1','lp2','lp3']

      for po in range(0,2):
        if mainPars.has_key(refsPhase[po][i]):
          self.freqDimRefs[-1].phase[po] = mainPars[refsPhase[po][i]].values[0]
        else:
          print "  ERROR: %s has no phase for %d (dim %d)" % (self.fileDir,po,i)

class VarianAcqParFile(VarianGenericFile):

  def initialize(self,parHelp = None):
  
    self.parameters = {}
    
    self.parHelp = parHelp # This is a dict: corresponds to VarianAcqParHelp().tags

  def read(self,verbose = 0, expname = 'unknown'):

    if verbose == 1:
    
      print "Reading Varian acquisition parameters file %s" % self.name

    # For reading: Varian parameter stored format is made up of 3 OR MORE lines!
    # See VarianPar class definition for more detailed info.

    fin = open(self.name, 'rU')
    line = fin.readline()
    
    while line:
      # Checking of empty lines and hashes not really necessary, but left in just in case
      if not (self.patt['emptyline'].search(line) or self.patt['hash'].search(line)):
        cols = line.split()
        if len(cols) != 11:
          print "  ERROR: %s Varian procpar has wrong number of fields on following line:" % expname + self.newline + line
        elif self.parameters.has_key(cols[0]):
          print ("Reoccuring parameter name: %s" % cols[0]) + self.newline
        else:
          parname = cols[0]
          self.parameters[parname]=VarianAcqPar(parname,cols[1],cols[2],cols[3],cols[4],cols[5],cols[6],cols[7],cols[8],cols[9])
          self.parameters[parname].addHelp(self.parHelp,verbose)

          # Read values (LINE 2)
          line = fin.readline()
          cols = line.split()
          value_len = returnInt(cols[0])

          if self.parameters[parname].basictype == 1:
            # Case real data type
            self.parameters[parname].addvalues(expname,value_len,returnFloats(cols[1:]))
          elif self.parameters[parname].basictype == 2:
            # Case string data type
            values = []
            values.append(joinUnquote(cols[1:],'"'))

            # For multiple line string data, read other lines
            for value_num in range(1,value_len):
              # Read values
              line = fin.readline()
              cols = line.split()
              values.append(joinUnquote(cols,'"'))

            self.parameters[parname].addvalues(expname,value_len,values)

          else:
            # Case undefined data type
            print "  ERROR: %s Varian procpar has undefined basic data type. Check parameter %s." % (expname,parname)

          # Read enumerations (LAST LINE)
          line = fin.readline()
          cols = line.split()
          enum_len = returnInt(cols[0])

          if enum_len == 0:
            self.parameters[parname].addenum(expname,enum_len,'')          
          else:
            self.parameters[parname].addenum(expname,enum_len,unquote(cols[1:],'"'))          

      line=fin.readline()

class VarianAcqPar:

  def __init__(self,name,subtype,basictype,maxlim,minlim,stepsize,Ggroup,Dgroup,protection,active):
    # LINE 1
    # Varian name of parameter (string)
    self.name = name
    
    # Subtype (integer)
    #  0: undefined
    #  1: real
    #  2: string
    #  3: delay
    #  4: flag
    #  5: frequency
    #  6: pulse
    #  7: integer
    self.subtype = returnInt(subtype)
    
    # Basictype (integer)
    #  0: undefined
    #  1: real
    #  2: string
    self.basictype = returnInt(basictype)

    # Maximum limit of value parameter (real number)
    # OR IF ARRAY: index to max value in parameter parmax (in conpar file)
    # For both strings and real basic types.
    self.maxlim = returnFloat(maxlim)
    
    # Mimimum limit of value parameter (real number)
    # OR IF ARRAY: index to min value in parameter parmin (in conpar file)
    # For real basic types ONLY.
    self.minlim = returnFloat(minlim)
    
    # Step size in which parameters can be entered
    # OR IF ARRAY: index to step size in parameter parstep (in conpar file)
    # For real basic types ONLY.
    #  0: ignored
    self.stepsize = returnFloat(stepsize)
    
    # Ggroup (integer)
    #  0: ALL
    #  1: SAMPLE
    #  2: ACQUISITION
    #  3: PROCESSING
    #  4: DISPLAY
    #  5: SPIN
    self.Ggroup = returnInt(Ggroup)
    
    # Dgroup (integer)
    # Defined by specific application!
    self.Dgroup = returnInt(Dgroup)
    
    # Protection: 32 bit word made up of following bit masks:
    #  bit  value   Description
    #  0    1     	Cannot array parameter
    #  1    2     	Cannot change active/not active status
    #  2    4     	Cannot change parameter value
    #  3    8     	Causes _parameter macro to be executed
    #  4    16     	Avoids automatic redisplay
    #  5    32    	Cannot delete parameter
    #  6    64    	System parameter for spectrometer or data station
    #  7    128    	Cannot delete parameter
    #  8    256   	Cannot set array parameter
    #  9    512   	Cannot set parameter enumeral values
    #  10   1024    Cannot change parameter's group
    #  11   2048  	Cannot change protection bits
    #  12   4096  	Cannot change display group
    #  13   8192  	Take max,min,step from conpar file (parmax, parmin, parstep)
    self.protection = returnInt(protection)
    
    # Active
    #  0: NOT ACTIVE
    #  1: ACTIVE
    self.active = returnInt(active)
    
    # Intptr (not used, generally set to 64)
    #self.intptr = intptr

  def addvalues(self,expname,value_len,values):
    # LINE 2 (or group of lines starting here)
    # First field on line is number of values belonging to parameter
    # Cases:
    #  Basic type  First value  Parameters
    #    1 (real)   any number  all listed on line 2, space separated
    #    2 (string)     1       single string quoted by " listed in second field
    #    2 (string)    >1       first string quoted by " listed in second field
    #                           other quoted by " on subseq. lines, one per line
    
    if value_len != len(values):
      print "  ERROR: %s Varian procpar: given value_len does not match number of values (%s)" % (expname,self.name)
      self.value_len = len(values)
    else:
      # Number of values for this parameter
      self.value_len = value_len
    
    # Values (always make this array for ease of use)
    self.values = values
        
  def addenum(self,expname,enum_len,enums):
    # LINE 3 (or last line if group of lines)
    # Lists enumerable values of string or flag (ALWAYS STRING!)
    # First field is number of enumerations
    # All listed on one line

    # Number of enumeration
    self.enum_len = enum_len
    
    if enum_len > 0:
      if enum_len != len(enums):
        print "  ERROR: %s Varian procpar: given enum_len does not match number of enums (%s)" % (expname,self.name)
        self.enum_len = len(enums)
      # Enumeration
      self.enums = enums
      
  def addHelp(self,parHelp,verbose):
    if parHelp.has_key(self.name):
      self.help = parHelp[self.name]
    else:
      if verbose > 1:
        print "Parameter %s: no help string found." % self.name
      self.help = ''
    
