#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

projectIO.py: I/O for autoassign project file

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
from ccp.format.autoAssign.generalIO import AutoAssignGenericFile

from ccp.format.autoAssign.peaksIO import AutoAssignPeakFile
from ccp.format.autoAssign.sequenceIO import AutoAssignSequenceFile

#####################
# Class definitions #
#####################

class AutoAssignProjectFile(AutoAssignGenericFile):

  # Information on file level

  def initialize(self):

    self.peakFiles = []
    self.tolerances = {}
    self.sequenceFile = None
    
    self.toleranceAtoms = ['HN','N15','CA','CB','HA','CO']
    self.toleranceDefault = [0.02,0.3,0.5,0.5,0.1,0.5]

    self.defaultSpectra = [ ('HSQC',1,0,0),
                            ('HNCO',0,1,0),
                            ('HNcaCO',1,1,0),
                            ('HNCA',1,1,0),
                            ('HNcoCA',0,1,0),
                            ('HNCACB',1,1,0),
                            ('HNcoCACB',0,1,0),
                            ('HNHA',1,1,0),
                            ('HNcoHA',0,1,0)]

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s project file %s" % (self.format,self.name)

    #
    # Read the sequence from the project file
    #
    
    self.sequenceFile = AutoAssignSequenceFile(self.name)
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
              
              peakFile = AutoAssignPeakFile(os.path.join(baseName,specFile))
              peakFile.read()
              peakFile.setProjectSpectrumInfo(specName,rootSpec,intra,sequential,throughSpace,phase,dimInfo)
              
              self.peakFiles.append(peakFile)

            else:
              line = fin.readline()

      line = fin.readline()

    fin.close()

  def write(self,verbose = 0):

    if verbose == 1:
      print "Writing %s project file %s" % (self.format,self.name)

    fout = open(self.name,'w')
    
    molName = self.sequenceFile.sequences[0].molName

    fout.write("#" + self.newline)
    fout.write("# %s Table File for AutoAssign" % molName + self.newline)
    fout.write("#" + self.newline)

    fout.write("Protein: %s" % molName + self.newline * 2)
    
    seqString = ""
    seqElements = self.sequenceFile.sequences[0].elements
    for seqEl in seqElements:
      seqString += seqEl.code1Letter
    
    fout.write("Sequence: %d %s*" % (seqElements[0].seqCode,seqString) + self.newline * 2)
    
    toleranceString = ""
    for toleranceAtom in self.toleranceAtoms:
      toleranceString += "%s %.3f " % (toleranceAtom,self.tolerances[toleranceAtom])
    
    fout.write("Tolerances: %s" % toleranceString + self.newline * 2)
    
    fout.write("Spectra:" + self.newline)
    
    for peakFile in self.peakFiles:
      
      (specDir,specFileName) = os.path.split(peakFile.name)
      
      fout.write("%s   %s   %s   %d   %d   %d phase: {%s} {" % (peakFile.specName,
                                                                peakFile.rootSpecName,
                                                                specFileName,
                                                                peakFile.intra,
                                                                peakFile.sequential,
                                                                peakFile.throughSpace,
                                                                peakFile.phase) + self.newline)

      for (atomInfo,lowPpm,highPpm,correction,tolerance,foldStatus) in peakFile.dimInfo:
        
        if len(atomInfo) > 1:
          atomString = "{%s }" % string.join(atomInfo,' ')
        else:
          atomString = atomInfo[0]

        fout.write("{%-10s %6.2f %6.2f  %s   %.2f  %s }" % (atomString,lowPpm,highPpm,
                                                            correction,tolerance,foldStatus) + self.newline)
        
      fout.write("}" + self.newline)
        
    fout.write("$" + self.newline)    
    fout.close()
    
###################
# Main of program #
###################

if __name__ == "__main__":
                                                      
  files = [('../../reference/autoassign/fas_table.aat','local/test1.aat')]

  
  for (inFile,outFile) in files:
    
    projectFile = AutoAssignProjectFile(os.path.join(getTopDirectory(),inFile))

    projectFile.read(verbose = 1)
    
    for peakFile in projectFile.peakFiles:
      print peakFile.specNames
      print peakFile.numDims
      print peakFile.dimInfo
      print peakFile.peaks[0].ppm
      print peakFile.peaks[-1].ppm
    
    for sequence in projectFile.sequenceFile.sequences:
      for seqEl in sequence.elements:
  
        print seqEl.seqCode, seqEl.code1Letter
    
    projectFile.name = outFile
    
    projectFile.write()
    
