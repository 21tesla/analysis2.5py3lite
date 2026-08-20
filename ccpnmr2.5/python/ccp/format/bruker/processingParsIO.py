#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

processingParsIO.py: I/O for Bruker processed data parameter files

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
from memops.universal.Util import returnInts

from memops.general.Implementation import ApiError

import os, string

from ccp.format.general.formatIO import GenericProcessingParsFile
from ccp.format.bruker.generalIO import BrukerGenericFile, BrukerParFile, BrukerPar

from memops.universal.Io import getTopDirectory
from memops.universal.Io import splitPath, joinPath

###
##################
# Class definitions #
#####################

class BrukerProcessingParData(BrukerGenericFile,GenericProcessingParsFile):

  def initialize(self,parHelp):

    self.parHelp = parHelp
    
    (self.fileDir,void) = splitPath(self.name)

    self.scriptName = 'processingParsIO.py'
    self.parFiles = {}
    
    self.allPars = {
    
      'mainPars': [
    
        ('numDim','ndim',2,0),
        ('fileName','file',2,''),
        ('','head',2,0),
        ('','integer',1,True),
        ('','swap',1,False),
        ('BYTORDP','big_endian',1,None),
        ('','varian',-1,None),
        ('NC_proc','brukerScale',2,0),
    
       ],
       
      'dimPar': [
      
        ('','dim',2,-1)
        
       ],
    
      'dimPars': [
      
        ('SI','npts',2,None),
        ('FTSIZE','numPointsOrig',2,None),
        ('STSR','pointOffset',2,0),
        ('XDIM','block',2,None),
        ('SW_p','sw',2,1000.0),
        ('SF','sf',2,1000.0),
        ('OFFSET','refppm',2,0.0),
        ('refPoint','refpt',2,1.0),
        ('AXNUC','nuc',2,'1H'),
        ('FT_mod','dimType',2,0),
    
       ]
    }

    self.initGeneric()
  
  def get(self,verbose = 0):

    ###################################
    # Read proc(n)s and proc(n) files #
    ###################################
    
    files = set(os.listdir(self.fileDir))
    baseFileName = 'proc'      
    
    self.numDim = 1
    while True:
      i = self.numDim + 1
      if ('proc%ds' % i) in files:
        self.numDim = i
      else:
        break
    
    self.parFiles.update({baseFileName + 's': {},baseFileName:{}})

    if self.numDim >=2:
      
      for n in range(2,self.numDim + 1):

        fileName = baseFileName + str(n)
        self.parFiles.update({fileName + 's': {},fileName:{}})
        

    for file in self.parFiles.keys():
      if file in files:
        fileName = joinPath(self.fileDir,file)
        
        self.parFiles[file] = BrukerParFile(fileName,parHelp = self.parHelp)
        self.parFiles[file].read(verbose = verbose)
      
      else:
        raise ApiError(" no file called %s in %s" % (file, self.fileDir))
    
    # CHANGE - give proper error if files not found
    #for file in files:
    #
    #  if self.parFiles.has_key(file):
      
    #    #
    #    # Get proc/procs files: most information comes from here 
    #    #
    #    
    #    fileName = joinPath(self.fileDir,file)
    #    
    #    self.parFiles[file] = BrukerParFile(fileName,parHelp = self.parHelp)
    #    self.parFiles[file].read(verbose = verbose)
    
    """
    ###################################
    # Read acqu(n)s and acqu(n) files #
    ###################################
    
    fileDir = os.path.dirname(os.path.dirname(self.fileDir))
    files = set(os.listdir(fileDir))
    baseFileName = 'acqu'      
    
    parFiles = {}
    parFiles.update({baseFileName + 's': {},baseFileName:{}})

    if self.numDim >=2:
      
      for n in range(2,self.numDim + 1):

        fileName = baseFileName + str(n)
        parFiles.update({fileName + 's': {},fileName:{}})
        

    for file in parFiles.keys():
      if file in files:
        fileName = joinPath(fileDir,file)
        
        parFiles[file] = BrukerParFile(fileName,parHelp = self.parHelp)
        parFiles[file].read(verbose = verbose)
      
      else:
        raise ApiError(" no file called %s in %s" % (file, fileDir))
    
    self.parFiles.update(parFiles)
    """
    
    self.process()

  def process(self):

    #
    # ACQORDER?called 321 and 312 (cf. the chapter on 3D processing for details).
    # ???AQSEQ (for F3)
    
    #  BYTORDP  endianness of data... 0 little, 1 big?
    # FT_mod  is TPPI, ...
    #  FTSIZE  Size of real spectrum before truncation (mostly == SI)
    # MC2 (for F1,F2). Phase type of transform...
    #  OFFSET chemical shift of the first data point of the spectrum!!
    # PHC0, PHC1  phasing
    # PH_mod   mode after transform... 0 = no further processing 1 = phase correction with the parameters PHC0, PHC1 2= calculate the magnitude spectrum 3 = calculate the power (squared magnitude) spectrum
    #  SF  spectro freq.
    #  SI  Size of real spectrum (same or larger than TD, always power of 2)
    # STSR (strip start, measured in points from the beginning of the spectrum)
    # STSI (size of region in points).    
    #  SW_p Sweep width Hz
    # TDeff number of original points contributing to transform! If zero, all pts.
    # TDoff is offset (point number where transform starts) if negative, zero-fills at start...
    # WDW  window function that was applied
    #  XDIM subblock size!
    #
    
    #
    # Using current information from proc<n>s files
    #
    
    procDataFile = '%d%s' % (self.numDim,self.numDim * 'r')
    mainPars = self.parFiles['procs'].parameters
      
    for parInfo in self.allPars['mainPars']:

      parKey = parInfo[0]

      if parKey:
    
        if parKey == 'numDim':
          self.setValue(parInfo,self.numDim)
        
        elif parKey == 'fileName':

          self.setValue(parInfo,os.path.join(self.fileDir,procDataFile))
      
        else:
          value = int(mainPars[parKey].values[0])
          self.setValue(parInfo,value)
      
      # TODO: head, integer, swap?  
    
    self.initDims()
    
    for dimIndex in range(0,self.numDim):
    
      if dimIndex > 0:
        mainPars = self.parFiles['proc%ds' % (dimIndex + 1)].parameters
    
      dimPars = self.allPars['dimPars']
      for parInfo in dimPars:

        parKey = parInfo[0]
    
        if parKey:
          
          if parKey in ('SI','XDIM'):
            value = int(mainPars[parKey].values[0])
          
          elif parKey in ('FTSIZE', ):
            if mainPars.get(parKey):
              value = int(mainPars[parKey].values[0])
            else:
              value = 0
          
          elif parKey in ('STSR', ):
            if mainPars.get(parKey):
              value = int(mainPars[parKey].values[0])
              #value = max(value - 1, 0)
              # NBNB above based on RHF reading of Bruker manual.
              # It gave incorrect referencing, so the below was substituted
              value = max(value, 0)
            else:
              value = 0
            
          elif parKey in ('SW_p','SF','OFFSET'):
            value = float(mainPars[parKey].values[0])
            
          elif parKey == 'refPoint':
            value = 1.0
          
          elif parKey == 'AXNUC':
            if mainPars.get(parKey):
              value = mainPars[parKey].values[0]
            else:
              value = None

          elif parKey == 'FT_mod':
            if mainPars[parKey].values[0] == 0:
              value = self.dimTypeNotProcessed
            else:
              value = self.dimTypeFrequency

          self.setListValue(parInfo,dimIndex,value)
   
###################
# Main of program #
###################

if __name__ == "__main__":

  # Proc files in 4/4/1, 5/5/1,2,3, 11/8/1

  parlist = {}
  mol = '4'
  lab = '4'
  expset = '1'
  rawdatadir = os.path.join(getTopDirectory(),'data','den')
  
  from acqParsHelpReader import BrukerParHelp
  
  help = BrukerParHelp('help')

  for expname in ['HnCA.hc.1','HncoCA.hc.1','HncoCA.hc.2','hsqc_HN.hn.1']:
     
    indir = os.path.join(rawdatadir,mol,lab,expset,expname,'info.bruker','pdata','1')
    
    (expType,expNucl,expNum) = string.split(expname,'.')
    
    procParData = BrukerProcessingParData(os.path.join(indir,'procs'),help.tags)
    procParData.get(verbose = 1)
    
    print procParData.aPars
    print procParData.fPars
    
			
