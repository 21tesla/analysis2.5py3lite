"""
======================COPYRIGHT/LICENSE START==========================

distanceConstraintsIO.py: I/O for XPLOR/CNS distance constraint files

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

import os

from memops.universal.Io import getTopDirectory

from memops.universal.Util import returnInt
from memops.universal.Util import returnFloat, returnFloats

from ccp.format.cns.generalIO import CnsGenericFile

#####################
# Class definitions #
#####################

class CnsDistanceConstraintFile(CnsGenericFile):

  def initialize(self):

    self.constraints = []

    self.constraintElements = 2

  def checkLinePattern(self,line):

    if self.patt[self.format + 'ChemShiftFormat'].search(line):
      print "  Can only read assigned distance restraint lists"
      return None
  
    return 1

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading cns distance constraint list %s" % self.name

    return self.readGeneric(CnsDistanceConstraint)


  def write(self,oneLinePerItem = False, verbose = 0):

    #
    # Output format is (example following, when peak info available):
    #
    # ASSI { 2818}
    #   (  segid "    " and resid 33   and name HD2%)
    #   (( segid "    " and resid 38   and name HZ  ))
    #      3.800     1.800     1.700 peak  2818 weight 1.00000E+00 volume 1.91500E+04 ppm1      0.572 ppm2      7.074
    # OR { 2818}
    #   (( segid "    " and resid 22   and name HG3 ))
    #   (( segid "    " and resid 38   and name HZ  ))
    #
    # ...

    if verbose == 1:
      self.printWriteText()

    fout = open(self.name,'w')

    for constraint in self.constraints:

      # TODO: use Id or origId or peakNum!??!?

      if hasattr(constraint,'peakNum'):
        fout.write(" ASSI {%5d}" % constraint.peakNum)

      else:
        fout.write(" assi ")

      if not oneLinePerItem:
        fout.write(self.newline)

      for i in range(len(constraint.items)):

        item = constraint.items[i] 

        for member in item.members:

          if member.chainCode != self.defaultMolCode:
            # TODO: this is not correct! Contact Jens!!
            segIdString = "segid \"%4s\" and " % member.chainCode
          else:
            segIdString = ""

          if not member.atomName or self.patt[self.format + 'MultiSign'].search(member.atomName[-1]):
            braces = 1
          else:
            braces = 2
            
          #
          # Split up cases with | in atom name...
          #
          
          if member.atomName:
            atomNameString = "and "
            if member.atomName.count("|"):
              (atomName1,atomName2) = member.atomName.split("|")
              atomNameString += "(name %-4s or name %-4s)" % (atomName1,atomName2)
            else:
              atomNameString += "name  %-4s" % member.atomName

          else:
            atomNameString = ""

          fout.write("   %-2s %sresid %-4d %s%s" % ('(' * braces,
                                                               segIdString,
		                       member.seqCode,
		                       atomNameString,
		                       ')' * braces))

          if not oneLinePerItem:
            fout.write(self.newline)


        # Write out peak and constraint info only after first item	                  
        if i == 0:
          peakString = ""

          if hasattr(constraint,'peakNum'):

            if constraint.weight == None:
              weight = 1
              
            else:
              weight = constraint.weight
              
            peakString = (" peak %5d weight  %8.5E volume  %8.5E" % (constraint.peakNum,
                                                                     weight,
			         constraint.peakVol))

            for j in range(0,len(constraint.ppms)):
              ppm = constraint.ppms[j]
              peakString += " ppm%d   %8.3f" % (j+1,ppm)
          
          # Fix by Benjamin Bardiaux, Institut Pasteur, to write out weights in absence of peaks
          elif constraint.weight != None and constraint.weight != 1.:
            peakString = " weight  %8.5E" % constraint.weight
          

          fout.write("  %9.3f %9.3f %9.3f%s" % (constraint.targetDist,
                                                constraint.minusDist,
                                                constraint.plusDist,
		                                            peakString))
          
          fout.write(self.newline)

        if i != len(constraint.items) - 1:

          if i!=0 and oneLinePerItem:
            fout.write(self.newline)

          if hasattr(constraint,'peakNum'):
            fout.write(" OR {%5d}" % constraint.peakNum)

          else:
            fout.write(" or   ")

          if not oneLinePerItem:
            fout.write(self.newline)
        
        else:
          fout.write(self.newline)
          
          if len(constraint.items) > 1 and oneLinePerItem:
            fout.write(self.newline)

    fout.close()


  def printWriteText(self):
  
    print "Writing cns distance constraint list %s" % self.name

class CnsDistanceConstraint:

  def __init__(self,Id,origId,patt = None,format = None):
    
    self.Id = returnInt(Id)
    
    if origId != None:
      self.origId = returnInt(origId)
    else:
      self.origId = None
    
    self.items = []
    self.patt = patt
    self.format = format
    self.weight = None
    
  def setDistanceData(self,targetDist,minusDist,plusDist):
  
    self.targetDist = returnFloat(targetDist)
    self.minusDist = returnFloat(minusDist)
    self.plusDist = returnFloat(plusDist)
    
  def setOtherData(self,weight,peakNum,peakVol,ppms):

    if weight != None:
      self.weight = returnFloat(weight)
    else:
      self.weight = None
      
    self.peakNum = returnInt(peakNum)
    self.peakVol = returnFloat(peakVol)
    self.ppms = returnFloats(ppms)

  def getSpecificInfo(self,assiLine): 
    
    lineRead = False

    #
    # First get distances...
    #
  
    distanceInfo = self.patt[self.format + 'RestrDistances'].search(assiLine)
  
    if distanceInfo:
      distances = assiLine[distanceInfo.start():].split()
      targetDist = distances[0]
      minusDist = distances[1]
      plusDist = distances[2]

      self.setDistanceData(targetDist,minusDist,plusDist)
      
      lineRead = True
      
    #
    # Peak and distance info line
    #
    
    if lineRead:
      
      peakInfo = self.patt[self.format + 'DistancePeakInfoLine'].search(assiLine)
     
      if peakInfo:
    
        cols = assiLine[peakInfo.start():].split()
              
        ppms = []
        weight = None
        peakVol = 0.0
        
        # Also have to account for lines with '='. Are spaces allowed around the '='? Currently NOT acceptable.
        i = 0
        while i < len(cols):
          if cols[i].count('='):
            curCols = cols[i].split('=')
            i += 1
          else:
            curCols = cols[i:i+2]
            i += 2
            
          if curCols[0] == 'peak':
            peakNum = curCols[1]
          elif curCols[0] == 'weight':
            weight = curCols[1]
          elif curCols[0] == 'volume':
            peakVol = curCols[1]
          elif curCols[0][1:3] == 'pm':   # Also takes hpm1 into account
            try:
              ppms.append(curCols[1])
            except:
              print curCols
              raise
  
        self.setOtherData(weight,peakNum,peakVol,ppms)
  
      # Fix by Benjamin Bardiaux, Institut Pasteur, to read in weights in absence of peaks
      else:
        weightString = assiLine.find('weight')
        if weightString > -1:
          weightLine = assiLine[weightString:]
          if weightLine.count('='):
            weightInfo = weightLine.split('=')
          else:
            weightInfo = weightLine.split()
          w = returnFloat(weightInfo[1])
          self.weight = w
      
    return lineRead
