#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

dihedralConstraintsIO.py: I/O for Diana/Dyana dihedral constraint files

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
from memops.universal.Io import getTopDirectory
 
from ccp.format.dyana.generalIO import DyanaGenericFile
from ccp.format.dyana.generalIO import DyanaConstraintItem
from ccp.format.dyana.generalIO import DyanaConstraintMember

from ccp.format.general.Util import getSeqAndInsertCode

# Diana/Dyana/Cyana data...
from ccp.format.cyana.cyanaLibParser import CyanaLibrary

#####################
# Class definitions #
#####################

class DyanaDihedralConstraintFile(DyanaGenericFile):

  def initialize(self):
  
    self.constraints = []
    
    self.cyanaLib = CyanaLibrary(version = self.version)
    
    self.constraintElements = 4

  def read(self,verbose = 0):

    print "Using CYANA library - courtesy of Peter Guentert."

    if verbose == 1:
      print "Reading %s dihedral constraint list %s" % (self.format,self.name)

    dihedralId = 0
    chainCode = self.defaultMolCode

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
      cols = line.split()

      if len(cols) == 0 or self.patt['hash'].search(line):
        pass

      else:

        dihedralId += 1

        keywds = {}

        if len(cols) > 5:
          keywds['energyCst'] = cols[5]

        dihedralConstraint = DyanaDihedralConstraint(dihedralId,cols[2],self.cyanaLib,verbose=False)
        if dihedralConstraint.Id != None:
          self.constraints.append(dihedralConstraint)

          self.constraints[-1].setAngleData(cols[3],cols[4],**keywds)
          if None in (self.constraints[-1].lowerAngle,self.constraints[-1].upperAngle):
            self.constraints.pop(-1)
          elif not self.constraints[-1].setAtomMembers(chainCode,cols[0],cols[1]):
            # Consistency check - remove if can't find this angle in CYANA library
            self.constraints.pop(-1)

      line = fin.readline()


  def write(self,verbose = 0):

    print "Using CYANA library - courtesy of Peter Guentert."

    #
    # Output format is (example following):
    #
    #  32 ARG+  PSI     -75.0   -15.0 1.00E-01
    #
    # Last column (weight) only listed if not 1.0
    # ...

    if verbose == 1:
      print "Writing %s dihedral constraint list %s" % (self.format,self.name)

    fout = open(self.name,'w')

    for constraint in self.constraints:

      constraint.getAngleInfo()
      
      if constraint.name:

        fout.write("%3d %-4s  %-5s%8.1f%8.1f" % (constraint.seqCode,
                                                 constraint.resLabel,
                                                 constraint.name,
                                                 constraint.lowerAngle,
                                                 constraint.upperAngle))

        if constraint.energyCst != 1.0:
          fout.write("%9.2E" % (constraint.energyCst))

        fout.write(self.newline)

      else:

        atomNames = string.join([member.atomName for member in constraint.items[0].members],'-')

        print "  Error: cannot find %s angle name for %s, atoms %s" % (self.format,constraint.resLabel,atomNames)

    fout.close()

class DyanaDihedralConstraint:

  def __init__(self,Id,name,cyanaLib,verbose=False):
    
    self.Id = returnInt(Id,default=None,verbose=verbose)
    self.items = []
    self.name = name
    self.seqCode = None
    
    self.cyanaLib = cyanaLib
    
  def setAngleData(self,lowerAngle,upperAngle,energyCst = 1.0,verbose=False):
  
    self.lowerAngle = returnFloat(lowerAngle,default=None,verbose=verbose)
    self.upperAngle = returnFloat(upperAngle,default=None,verbose=verbose)
    self.energyCst = returnFloat(energyCst,verbose=verbose)
      
  def setAtomMembers(self,chainCode,seqCode,resLabel):

    refAngle = self.cyanaLib.findAngle(resLabel,self.name)
    
    if not refAngle:
      print "  Warning: no %s angle available for %s. Ignored." % (self.name,resLabel)
      return False
    
    (self.seqCode,self.seqInsertCode) = getSeqAndInsertCode(seqCode)
    self.resLabel = resLabel
    
    self.items.append(DyanaConstraintItem())
    
    for atom in refAngle.atoms:
    
      if atom.location:
        curSeqCode = self.seqCode + atom.location
        curResLabel = ""
        
      else:
        curSeqCode = self.seqCode
        curResLabel = self.resLabel
        
      self.items[-1].members.append(DyanaConstraintMember(chainCode,curSeqCode,curResLabel,atom.name))
      
    return True
      
  def getAngleInfo(self):

    locations = []

    if len(self.items) > 1:
      print "  Warning can only handle one dihedral constraint item"

    item = self.items[0]
    locIndex = None

    if item.members[0]:

      firstSeqCode = item.members[0].seqCode

      for member in item.members:
        seqCodeDiff = member.seqCode - firstSeqCode

        locations.append(seqCodeDiff)

      if locations.count(0) >= 3:
        locIndex = locations.index(0)

      elif locations.count(1) in [2,3]:
        locIndex = locations.index(1)

      elif locations.count(-1) in [2,3]:
        locIndex = locations.index(-1)
      
    if locIndex == None:
      self.seqCode = 0
      self.resLabel = 'unknown'
    else:
      self.seqCode = item.members[locIndex].seqCode
      self.resLabel = item.members[locIndex].resLabel
 
