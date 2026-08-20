"""
======================COPYRIGHT/LICENSE START==========================

projectIO.py: I/O for MARS project file

Copyright (C) 2007 Wim Vranken (European Bioinformatics Institute)

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
from ccp.format.mars.generalIO import MarsGenericFile

#####################
# Class definitions #
#####################

class MarsProjectFile(MarsGenericFile):

  # Information on file level

  def initialize(self):
  
    self.files = {
    
      'sequence': None,
      'fixConn': None,
      'fixAss': None,
      'secondary': None,
      'csTab': None,
      'pdb': None
      
      }
      
    self.cutoffs = {
    
      'CO': 0.25,
      'CA': 0.20,
      'CB': 0.50,
      'HA': 0.25
      
    }
    
    self.fragSize = 5
    self.deuterated = 0
    
    self.pdbFileInfo = {
      
      'resolution': None,     # In angstrom
      'tensor':     None,     # Alignment tensor method 0/1/2/3/4
      'nIter':      None,     # Number iterations 2/3/4
      'dObsExh':    None,     # Name of RDC table (exhaustive search)
      'dcTab':      None      # Name of RDC table (PALES format)
    }
        

  def read(self,verbose = 0):
  
    # TODO: set this up or is it not really relevant?
    
    """
    if verbose == 1:
      print "Reading %s project file %s" % (self.format,self.name)

    #
    # Read the sequence from the project file
    #
    
    self.sequenceFile = MarsSequenceFile(self.name)
    self.sequenceFile.read()
    
    #
    # Read all the peak list info and the peak lists themselves
    # Assuming that they are in the same directory!
    #

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      cols = line.split()

      if len(cols) == 0 or self.patt['hash'].search(line):
        pass

      else:
      
        if cols[0] == self.tolerancesTag:
        
          for i in range(1,len(cols),2):
            self.tolerances[cols[i]] = returnFloat(cols[i+1])

        elif cols[0] == self.spectraTag:
          
          line = fin.readline()
          
          while not self.patt['onlyDollar'].search(line):
          
            cols = line.split()
            
            phaseSearchObj = self.patt[self.format + 'specPhase'].search(line)

            if phaseSearchObj:
            
              specName = cols[0]
              rootSpec = cols[1]
              specFile = cols[2]
              intra = returnInt(cols[3])
              sequential = returnInt(cols[4])
              throughSpace = returnInt(cols[5])
              phase = phaseSearchObj.group(1)
              
              dimInfo = []

              line = fin.readline()

              while not self.patt[self.format + 'onlyCurlyBraceEnd'].search(line):

                cols = line.split()
                
                dimInfoSearchObj = self.patt[self.format + 'CurlyBrace'].search(line)
                
                if dimInfoSearchObj:
                  infoLine = dimInfoSearchObj.group(1)
                  
                  atomGroupSearchObj = self.patt[self.format + 'CurlyBrace'].search(infoLine)
                  if atomGroupSearchObj:
                    atomInfo = atomGroupSearchObj.group(1).split()
                    infoLine = infoLine.replace(atomGroupSearchObj.group(0),'')
                  else:
                    atomInfo = infoLine.split()[0]
                    infoLine = infoLine.replace(atomInfo,'')
                    atomInfo = [atomInfo]
                  
                  (lowPpm,highPpm,correction,tolerance,foldStatus) = infoLine.strip().split()
                  dimInfo.append((atomInfo,returnFloat(lowPpm),returnFloat(highPpm),correction,returnFloat(tolerance),foldStatus))
                  
                line = fin.readline()
              
              (baseName,projectFileName) = os.path.split(self.name)
              
              peakFile = MarsPeakFile(os.path.join(baseName,specFile))
              peakFile.read()
              peakFile.setProjectSpectrumInfo(specName,rootSpec,intra,sequential,throughSpace,phase,dimInfo)
              
              self.peakFiles.append(peakFile)

            else:
              line = fin.readline()

      line = fin.readline()

    fin.close()
    """

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing %s project file %s" % (self.format,self.name)

    fout = open(self.name,'w')
    
    self.writeLine(fout,"fragSize",self.fragSize)
    fout.write(self.newline)
    
    cutoffKeys = self.cutoffs.keys()
    cutoffKeys.sort()
    
    for cutoffKey in cutoffKeys:
      self.writeLine(fout,"cutoff%2s" % cutoffKey,"%4.2f" % self.cutoffs[cutoffKey])
    fout.write(self.newline)

    for fixType in ('Conn','Ass'):
      fileType = "fix%s" % fixType
      self.writeLine(fout,fileType,self.files[fileType])
    fout.write(self.newline)
    
    if self.files['pdb']:
      pdbAvailable = 1
      self.pdbFileInfo['pdbName'] = self.files['pdb']
    else:
      pdbAvailable = 0
      self.pdbFileInfo['pdbName'] = None
      
    self.writeLine(fout,"pdb",pdbAvailable)
    
    # Does the order matter here? Assuming it doesn't...
    pdbInfoKeys = self.pdbFileInfo.keys()
    pdbInfoKeys.sort()
    
    for pdbInfoKey in pdbInfoKeys:
      self.writeLine(fout,pdbInfoKey,self.pdbFileInfo[pdbInfoKey])
    fout.write(self.newline)
      
    self.writeLine(fout,"deuterated",self.deuterated)    
    fout.write(self.newline)
    
    for fileType in ('sequence','secondary','csTab'):
      self.writeLine(fout,fileType,self.files[fileType])

    fout.write(self.newline)    
    fout.close()

  def writeLine(self,fout,key,value):
  
    if value == None:
      value = "NO"
  
    fout.write("%-13s %s" % ("%s:" % key,str(value)) + self.newline)
  
