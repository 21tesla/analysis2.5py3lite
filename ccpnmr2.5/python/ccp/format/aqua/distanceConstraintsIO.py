"""
======================COPYRIGHT/LICENSE START==========================

distanceConstraintsIO.py: I/O for AQUA distance constraint files

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

from ccp.format.aqua.generalIO import AquaGenericFile, AquaConstraintItem, AquaConstraintMember

from ccp.format.general.Constants import defaultLowerDist
from ccp.format.general.Constants import defaultMolCode

from memops.universal.Util import returnInt, returnFloat

#####################
# Class definitions #
#####################

class AquaDistanceConstraintFile(AquaGenericFile):

  def initialize(self):
  
    self.constraints = []
    
    self.constraintElements = 2
    
    self.comment = None
    self.originalType = None

  def read(self,verbose = 0):

    if verbose == 1:
      
      print "Reading %s distance constraint list %s" % (self.format,self.name)


    fin = open(self.name, 'rU')

    self.chainCode = self.defaultMolCode
 
    line = fin.readline()
    
    constraintDict = {}
    numConstraints = 0
    listedNumConstraints = 0
    
    while line:

      if self.patt['hash'].search(line) or self.patt['emptyline'].search(line): # TODO: is this correct for AQUA comments?
        line = fin.readline()
        continue
      
      # Check if header line
      if self.comment == None:
        self.comment = line.strip()
        line = fin.readline()
        continue
        
      # Check if count line
      countSearch = self.patt[self.format + 'Count'].search(line)
      if countSearch:
        self.originalType = countSearch.group(2).strip()
        listedNumConstraints = returnInt(countSearch.group(1))
        line = fin.readline()
        continue

      cols = string.split(line)
        
      numConstraints += 1
  
      resLabel1 = cols[1]
      seqCode1 =  cols[2]
      atomName1 = cols[3]      
      fullId1 = (seqCode1,atomName1)
      
      resLabel2 = cols[4]
      seqCode2  = cols[5]
      atomName2 = cols[6]
      fullId2 = (seqCode2,atomName2)
      
      distance  = cols[7]

      if len(cols) > 8:
        commentSearch = self.patt[self.format + 'EndComment'].search(line)
        comment = commentSearch.group(1).strip()
      else:
        comment = ""
      
      
      fullKey = (fullId1,fullId2)
      
      if constraintDict.has_key(fullKey):
        constraint = constraintDict[fullKey]

      else:
        constraint = AquaDistanceConstraint(numConstraints)
        
        constraintItem = AquaConstraintItem()
        constraint.items.append(constraintItem)
        
        constraintMember1 = AquaConstraintMember(self.chainCode,seqCode1,resLabel1,atomName1)
        constraintMember2 = AquaConstraintMember(self.chainCode,seqCode2,resLabel2,atomName2)
        
        constraint.items[-1].members.append(constraintMember1)
        constraint.items[-1].members.append(constraintMember2)
        
        self.constraints.append(constraint)
        constraintDict[fullKey] = constraint
        
      constraint.comment += comment
      
      if self.patt[self.format + 'UpperDistance'].search(cols[0]):
        constraint.setUpperDistance(distance)
        
      elif self.patt[self.format + 'LowerDistance'].search(cols[0]):
        constraint.setLowerDistance(distance)
        
      else:
        print "  Error: no distance found in %s format distance constraint file!" % self.format

      line = fin.readline()
      
    if numConstraints != listedNumConstraints:
      print "  Error: number of constraints read (%d) does not match number of constraints listed (%d)" % (numConstraints,listedNumConstraints)

    fin.close()
    
  def write(self,verbose = 0):
    
    if verbose == 1:
    
      print "Writing %s distance constraint list %s" % (self.format,self.name)
        
    #
    # Output format is (example following):
    #
    # FormatConverter AQUA distance list
    # count 7 type disgeo
    # NOEUPP  GLY    1 1HA   ALA    2 HN         3.000
    # NOEUPP  GLY    1 2HA   ALA    2 HN         3.000
    # NOEUPP  GLY    1 HN    ALA    2 HN         4.000
    # NOEUPP  GLY    1 HN    ALA    2 QB         4.000
    # NOELOW  GLY    1 HN    ALA    2 QB         2.500
    # NOELOW  GLY    1 1HA   ALA    2 HN         2.500
    # NOELOW  GLY    1 2HA   ALA    2 HN         0.500
    #
    # Lists are always unambiguous
    #

    fout = open(self.name,'w')
        
    fout.write("FormatConverter AQUA distance list" + self.newline)

    numConstraints = 0
    for constraint in self.constraints:
      if constraint.upperDist != None:
        numConstraints += 1
      elif constraint.lowerDist != None:
        numConstraints += 1  
    
    fout.write("count %d type %s" % (numConstraints,'ccpn') + self.newline)

    for constraint in self.constraints:

      firstItem = constraint.items[0]
      
      for (noeText,distance) in (('NOEUPP',constraint.upperDist),('NOELOW',constraint.lowerDist)):

        #
        # Ignore if distance is None
        #

        if distance == None:
          continue
          
        fout.write("%-8s" % noeText)

        for member in firstItem.members:
          fout.write("%-4s%4d %-5s " % (member.resLabel,member.seqCode,member.atomName))

        fout.write("%10.3f" % distance)
        
        if constraint.comment:
          fout.write(" # %s" % constraint.comment)
        
        fout.write(self.newline)

        if len(constraint.items) > 1:

          print "  Warning: AQUA can't handle ambiguous restraints - left out items:"

          for i in range(1,len(constraint.items)):
            print "    ",i,
            for j in  range(0,2):
              member = constraint.items[i].members[j]
              print ("%3d %-4s %-5s " % (member.seqCode,member.resLabel,member.atomName)),
            print

    fout.close()

class AquaDistanceConstraint:

  def __init__(self,Id):
        
    self.Id = Id
    
    self.items = []
    
    self.upperDist = None
    self.lowerDist = None
        
    self.comment = ""
    
  def setUpperDistance(self,upperDist):
  
    upperDist = returnFloat(upperDist,None)
    
    if not self.upperDist or upperDist < self.upperDist:
      self.upperDist = upperDist
    
  def setLowerDistance(self,lowerDist):
    
    lowerDist = returnFloat(lowerDist,defaultLowerDist)

    if not self.lowerDist or lowerDist > self.lowerDist:
      self.lowerDist = lowerDist
