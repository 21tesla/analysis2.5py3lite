#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

distanceConstraintsIO.py: I/O for Concoord distance constraint files

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

from memops.universal.Io import getTopDirectory

from memops.universal.Util import returnInt, returnFloat

from ccp.format.concoord.generalIO import ConcoordGenericFile
from ccp.format.concoord.generalIO import ConcoordConstraintItem
from ccp.format.concoord.generalIO import ConcoordConstraintMember

from ccp.format.concoord.generalIO import returnSegId

#####################
# Class definitions #
#####################

class ConcoordDistanceConstraintFile(ConcoordGenericFile):

  def initialize(self):
  
    self.constraints = []
    self.constraintElements = 2
    
    self.constraintFileType = 'distance_restraints'
    
    self.addColumnNames = ['index','restrnum']
    
  def read(self,verbose = 0):

    if verbose == 1:
    
      print "Reading %s distance constraint list %s" % (self.format,self.name)
    
    fin = open(self.name, 'rU')

    self.chainCode = self.defaultMolCode
 
    #
    # First line has type
    #
    
    line = fin.readline()
    
    restraintType = string.split(line)[1]
    
    #
    # Second line has column headers
    #
    
    line = fin.readline()
    
    columnNames = string.split(line)[1:-1]
    
    line = fin.readline()
    
    constraintNum = oldConstraintNum = -1
    
    while line:

      if self.patt['emptyline'].search(line):
        line = fin.readline()
        continue

      values = string.split(line)

      constraintNum = values[columnNames.index('restrnum')]
      
      if constraintNum != oldConstraintNum:
      
        constraint = ConcoordDistanceConstraint(self,constraintNum)
        self.constraints.append(constraint)
        
        oldConstraintNum = constraintNum
        constraint.setDistances(values,columnNames)
      
      constraint.addItem(values,columnNames)
      
      line = fin.readline()
      
      
  def write(self,verbose = 0):
    
    if not self.constraints:
      return

    if verbose == 1:
 
      print "Writing %s distance constraint list %s" % (self.format,self.name)
    
    
    fout = open(self.name,'w')

    #
    # Write out distance constraints
    #
    
    fout.write("[ %s ]" % self.constraintFileType + self.newline)
    
    fout.write("[")
    
    for i in range(0,self.constraintElements):
      for columnName in self.constraints[0].itemRefs:
        fout.write(" %s" % columnName + str(i+1))
        
    for columnName in self.constraints[0].distRefs:
      fout.write(" %s" % columnName)
    
    for columnName in self.addColumnNames:
      fout.write(" %s" % columnName)

    fout.write(" ]" + self.newline)

    for constraint in self.constraints:
    
      if len(constraint.items) > 1:
        index = 1
      else:
        index = 0
        
      for item in constraint.items:
      
        fout.write(" ")

        for i in range(0,self.constraintElements):
        
          member = item.members[i]

          for columnName in constraint.itemRefs:
            (func,attrName,outFormat) = constraint.itemRefDict[columnName]
            fout.write(outFormat % getattr(member,attrName))

        for columnName in constraint.distRefs:

          (func,attrName,outFormat) = constraint.distRefDict[columnName]
          fout.write(outFormat % getattr(constraint,attrName))

        fout.write("%6d" % index)
        fout.write("%9d" % constraint.Id)      

        index += 1

        fout.write(self.newline)
    
class ConcoordDistanceConstraint:

  def __init__(self,parent,Id):
    
    self.parent = parent
    self.Id = returnInt(Id)
    
    self.items = []
    
    self.distRefs = ['lowbound','uppbound','upppscor']
    self.distRefDict = {
    
      'lowbound': [returnFloat,'lowerDist','%9.3f'],
      'uppbound': [returnFloat,'upperDist','%9.3f'],
      'upppscor': [returnFloat,'upperPseudoDist','%9.3f']
      
      }
    
    self.itemRefs = ['resid','resname','atomname','segid']
    self.itemRefDict = {
    
      'resid': [returnInt,'seqCode','%7d'],
      'resname': [None,'resLabel','%9s'],
      'atomname': [None,'atomName','%10s'],
      'segid': [returnSegId,'chainCode','%7s'],
      
      }
  
  def setDistances(self,values,columnNames):
    
    for i in range(0,len(columnNames)):
  
      columnName = columnNames[i]

      if self.distRefDict.has_key(columnName):
        
        value = values[i]
        func = self.distRefDict[columnName][0]
        
        if func:
          value = func(value)
        
        setattr(self,self.distRefDict[columnName][1],value)

  def addItem(self,values,columnNames):
  
    item = ConcoordConstraintItem()
    self.items.append(item)
    
    for i in range(1,1+self.parent.constraintElements):
    
      member = ConcoordConstraintMember()
      item.members.append(member)
      
      for j in range(0,len(columnNames)):

        columnName = columnNames[j]
        
        if columnName[-1] == str(i) and self.itemRefDict.has_key(columnName[:-1]):
         
          value = values[j]
          func = self.itemRefDict[columnName[:-1]][0]

          if func:
            value = func(value)
            
          setattr(member,self.itemRefDict[columnName[:-1]][1],value)
          
###################
# Main of program #
###################

if __name__ == "__main__":

  files = ['../reference/concoord/1d3z_concoord_pscorr_new.restr'
          ]
  
  for file in files:
    
    file = os.path.join(getTopDirectory(), file)
    
    constraintFile = ConcoordDistanceConstraintFile(file)

    constraintFile.read(verbose = 1)
    
    """
    for constraint in constraintFile.constraints:
      print constraint.Id,
      print constraint.upperDist, constraint.lowerDist, constraint.upperPseudoDist,

      #print constraint.peakNum, constraint.peakVol, constraint.ppms
      for item in constraint.items:
        for member in item.members:
          print member.chainCode,member.seqCode, member.atomName,
        print "|",
      
      print
     
    constraintFile.name = 'local/testout'

    constraintFile.write(verbose = 1)
    """
    sequence = {}
    
    for constraint in constraintFile.constraints:
      for item in constraint.items:
        for member in item.members:
          if not sequence.has_key(member.seqCode):
            sequence[member.seqCode] = member.resLabel

    seqnums = sequence.keys()

    for i in range(0,seqnums[-1]):
      
      seqnum = i +1
      
      if sequence.has_key(seqnum):
        print sequence[seqnum]
      else:
        print "ALA"
        
