#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

generalIO.py: General I/O information for Diana/Dyana files

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

from memops.universal.Util import returnInt, returnFloat

from ccp.format.general.formatIO import FormatFile

from ccp.format.general.Constants import defaultMolCode
from ccp.format.general.Constants import defaultLowerDist

from ccp.format.general.Util import getSeqAndInsertCode

#####################
# Class definitions #
#####################

class DyanaGenericFile(FormatFile):

  def setGeneric(self):
    
    self.setFormat()
    self.defaultMolCode = defaultMolCode

    self.linkerResidueCodes = ['PL','LL','LL2','LL5','LP']    
    self.chainCodesString = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

  def setFormat(self):
  
    self.format = 'dyana'

class DyanaGenericDistanceConstraintFile(DyanaGenericFile):

  def initialize(self,lowerDistanceFile = None):
  
    self.constraints = []
    
    self.constraintElements = 2
    
    self.lowerDistanceFile = lowerDistanceFile
    
  def readGeneric(self,ConstraintClass):
  
    self.ConstraintClass = ConstraintClass
  
    self.readGenericFile(self.name)
  
    if self.lowerDistanceFile:
    
      self.readGenericFile(self.lowerDistanceFile,upper = 0)

  def readGenericFile(self,fileName,upper = 1):

    fin = open(fileName, 'rU')

    self.restrNum = -1  # Or start at 0?

    self.chainCode = self.defaultMolCode
 
    line = fin.readline()
    
    itemList = []
    commentList = []
    
    if upper:
    
      self.constraintRefItemList = {}

    while line:

      if self.patt['hash'].search(line) or self.patt['emptyline'].search(line): # TODO: is this correct for DYANA comments?
        line = fin.readline()
        continue

      cols = line.split()

      #
      # Assuming that two columns means residue for first constr element listed separately
      #

      if len(cols) == 2:

        seqCode1 = cols[0]
        resLabel1 = cols[1]

      else:

        #
        # Assuming that nothing in beginning of line means: use previously defined seqCode1/resLabel1
        #

        if line[:9].split() == []:

          atomName1 = cols[0]
          seqCode2 = cols[1]
          resLabel2 = cols[2]
          atomName2 = cols[3]

          dist = cols[4]
          
          comments = cols[5:]

          # TODO: handle more recent/complicated DYANA
          #if len(cols) > 5:
          #  contrib = cols[5]

        #
        # Assuming this is full data line
        #

        else:

          seqCode1 = cols[0]
          resLabel1 = cols[1]
          atomName1 = cols[2]
          seqCode2 = cols[3]
          resLabel2 = cols[4]
          atomName2 = cols[5]

          dist = cols[6]
          
          comments = cols[7:]
        
        
        if seqCode1.isdigit() and seqCode2.isdigit():
          distance = returnFloat(dist,default=None,verbose=False)

          if distance == 0:

            itemList.append([(seqCode1,resLabel1,atomName1),(seqCode2,resLabel2,atomName2)])
            commentList.append(comments)

          elif distance != None:

            if itemList:

              #
              # Set up previous constraint or look for previous constraint.
              # Repeated at bottom...
              #

              self.setPreviousDistanceConstraint(itemList,upper,commentList)
              itemList = []
              commentList= []

            self.currentDist = dist
            itemList.append([(seqCode1,resLabel1,atomName1),(seqCode2,resLabel2,atomName2)])
            commentList.append(comments)
                      
      line = fin.readline()

    if itemList:
    
      self.setPreviousDistanceConstraint(itemList,upper,commentList)        

    fin.close()
  
  def setPreviousDistanceConstraint(self,itemList,upper,commentList):
  
    if upper:
          
      #
      # Always create a constraint for upper distance lists...
      #

      self.restrNum += 1
      self.constraints.append(self.ConstraintClass(self.restrNum))

      constraint = self.constraints[-1]
      constraint.setDistanceData(self.currentDist)
      
      #
      # Get constraint level comments...
      #
      
      if commentList and commentList[0]:
        for (setFuncName,commentString) in [('setPeakNum','#peak'),('setPeakNum','#PEAK'),('setQualityFactor','#QF')]:
          if commentString in commentList[0]:
            commentIndex = commentList[0].index(commentString)
            commentList[0].pop(commentIndex)
            getattr(constraint,setFuncName)(commentList[0].pop(commentIndex))            
   
      #
      # Set up new constraint item for each combination within itemList
      #

      self.constraintRefItemList[constraint] = itemList
      
      commentString = '#SUP'
      setFuncName = 'setSupport'
      
      for itemIndex in range(0,len(itemList)):

        item = itemList[itemIndex]
        
        constraint.items.append(DyanaConstraintItem())
        
        #
        # Set constraint item comments
        #

        if commentList[itemIndex]:
          if commentString in commentList[itemIndex]:
            commentIndex = commentList[itemIndex].index(commentString)
            commentList[itemIndex].pop(commentIndex)
            getattr(constraint.items[-1],setFuncName)(commentList[itemIndex].pop(commentIndex))            
        
        #
        # Set constraint item members
        #
        
        for itemMember in item:
          constraint.items[-1].members.append(DyanaConstraintMember(self.chainCode,itemMember[0],itemMember[1],itemMember[2]))

    else:
    
      #
      # Look for an existing constraint and set the lower distance...
      #
      
      for constraint in self.constraintRefItemList.keys():
        
        if itemList == self.constraintRefItemList[constraint]:
        
          constraint.setLowerDistance(self.currentDist)

  def writeGeneric(self,noAmbiguous):
  
    self.noAmbiguous = noAmbiguous

    self.writeGenericFile(self.name)
  
    if self.lowerDistanceFile:
    
      self.writeGenericFile(self.lowerDistanceFile,upper = 0)

  def writeGenericFile(self,fileName,upper = 1):
  
    #
    # Output format is (example following):
    #
    #  3 LYS  HN      2 PRO  HA      4.00
    #  3 LYS  HN      1 GLN  HA      0.00
    #
    # Second line with distance 0.00 means ambiguous combination
    # Only upper limit file currently
    #
    # If noAmbiguous, then no ambiguous restraints are written
    #

    fout = open(fileName,'w')

    for constraint in self.constraints:

      firstItem = constraint.items[0]

      if upper:
        dist = constraint.upperDist
      else:
        dist = constraint.lowerDist

      #
      # Ignore if distance is NOT None
      #
      
      if dist == None:
        continue

      for member in firstItem.members:
        fout.write("%3d %-4s %-5s " % (member.seqCode,member.resLabel,member.atomName))

      fout.write("%6.2f" % dist)
      fout.write(self.newline)

      if len(constraint.items) > 1 and self.noAmbiguous:

        print "  Warning: old DYANA can't handle ambiguous restraints - left out items:"

        for i in range(1,len(constraint.items)):
          print "    ",i,
          for j in  range(0,2):
            member = constraint.items[i].members[j]
            print ("%3d %-4s %-5s " % (member.seqCode,member.resLabel,member.atomName)),
          print

      else:

        for item in constraint.items[1:]:
          for member in item.members:
            fout.write("%3d %-4s %-5s " % (member.seqCode,member.resLabel,member.atomName))

          fout.write("%6.2f" % 0)
          fout.write(self.newline)


    fout.close()

class DyanaGenericDistanceConstraint:

  def __init__(self,Id):
    
    self.Id = returnInt(Id)
    
    self.items = []
    
    self.upperDist = None
    self.lowerDist = None
    
    self.contrib = None
    self.peakNum = None
    self.qualityFactor = None
    
  def setDistanceData(self,upperDist,lowerDist = None):
  
    self.upperDist = returnFloat(upperDist,None)    
    self.setLowerDistance(lowerDist)
    
  def setLowerDistance(self,lowerDist):
    
    if lowerDist:
      self.lowerDist = returnFloat(lowerDist,defaultLowerDist)
    else:
      self.lowerDist = defaultLowerDist
    
  def setOtherData(self,contrib):
      
    self.contrib = returnFloat(contrib)
    
  def setPeakNum(self,infoString):
    
    self.peakNum = returnInt(infoString)
  
  def setQualityFactor(self,infoString):
    
    self.qualityFactor = returnFloat(infoString)  
    
class DyanaConstraintItem:

  def __init__(self):
    
    self.members = []
    self.support = None
  
  def setSupport(self,infoString):
    
    self.support = returnFloat(infoString)
    
class DyanaConstraintMember:

  def __init__(self,chainCode,seqCode,resLabel,atomName):
    
    self.chainCode = chainCode
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.resLabel = resLabel
    self.atomName = atomName
