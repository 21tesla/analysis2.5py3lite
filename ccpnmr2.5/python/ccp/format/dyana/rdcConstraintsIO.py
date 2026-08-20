"""
======================COPYRIGHT/LICENSE START==========================

rdcConstraintsIO.py: I/O for Diana/Dyana rdc constraint files

Copyright (C) 2005-2010 Wim Vranken (European Bioinformatics Institute)

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

from memops.universal.Util import returnInt, returnFloat

from ccp.format.dyana.generalIO import DyanaGenericFile
from ccp.format.dyana.generalIO import DyanaConstraintItem
from ccp.format.dyana.generalIO import DyanaConstraintMember

# Diana/Dyana/Cyana data...
from ccp.format.cyana.cyanaLibParser import CyanaLibrary

#####################
# Class definitions #
#####################

class DyanaRdcConstraintFile(DyanaGenericFile):

  def initialize(self):

    self.constraints = []
    
    self.cyanaLib = CyanaLibrary(version = self.version)

    self.constraintElements = 2
    
    # There are two different types of file, the old ones have only one atom, new one has two atoms defined.
    self.singleAtom = True
    
    self.orientations = {}

  def read(self, verbose=False):

    print "Using CYANA library - courtesy of Peter Guentert."

    if verbose:
      print "Reading %s rdc constraint list %s" % (self.format,self.name)

    rdcId = 0
    chainCode = self.defaultMolCode
    
    #
    # Read the file
    #

    fin = open(self.name, 'rU')

    # Read, look for first line
    line = fin.readline()

    while line:
    
      cols = line.split()
        
      if line.count("Orientation") and line.count("Rhombicity"):
        # Set the file version to the most recent one
        self.singleAtom = False
        
        # Get the info for the coordinate frames
        while 1:
          coordFrameInfoLine = fin.readline()
          if self.patt['hash'].search(coordFrameInfoLine):
            break
          else:
            (orientation,magnitude,rhombicity,frameResidueNumber) = coordFrameInfoLine.split()
            self.orientations[returnInt(orientation)] = (returnFloat(magnitude),returnFloat(rhombicity),returnInt(frameResidueNumber))

      elif len(cols) == 0 or self.patt['hash'].search(line):
        pass

      else:

        rdcId += 1

        keywds = {}
        
        if self.singleAtom:
          value = cols[3]
          if len(cols) > 4:
            keywds['error'] = cols[4]
            if len(cols) > 5:
              keywds['energyCst'] = cols[5]
              
        else:
          value = cols[6]
          keywds['error'] = cols[7]
          keywds['energyCst'] = cols[8]
          keywds['orientation'] = cols[9]
          
        self.constraints.append(DyanaRdcConstraint(rdcId,self.cyanaLib))
        self.constraints[-1].setRdcData(value,**keywds)

        if self.singleAtom:
          self.constraints[-1].setAtomMembersSingle(chainCode,cols[0],cols[1],cols[2])
        else:
          self.constraints[-1].setAtomMembersDouble(chainCode,cols[:3],cols[3:6])

      line = fin.readline()

  def write(self, singleAtom=False, verbose=False):

    print "Using CYANA library - courtesy of Peter Guentert."

    #
    # Output format is (example following):
    #
    #  30 ILE  HN    4       2.0  1.4
    #
    # Last column (weight) only listed if not 1.0
    # ...

    if verbose:
      print "Writing %s rdc constraint list %s" % (self.format,self.name)
    
    #
    # Open output file
    #

    fout = open(self.name,'w')
    
    #
    # Write header
    #
    
    if not singleAtom:
    
      fout.write("#Orientation  Magnitude  Rhombicity  ORI residue number" + self.newline)
      
      orientationKeys = self.orientations.keys()
      orientationKeys.sort()
      
      for orientationKey in orientationKeys:
        (magnitude,rhombicity,frameResidueNumber) = self.orientations[orientationKey]
        fout.write("       %d       %8.3f    %8.3f          %4d" % (orientationKey,magnitude,rhombicity,frameResidueNumber))
        fout.write(self.newline)

      fout.write("#  First atom      Second atom           RDC   Error  Weight Orientation" + self.newline)

    #
    # Write constraints
    #

    for constraint in self.constraints:
        
      if singleAtom:
      
        atomNames = []
        resLabel = None

        for item in constraint.items:
          for member in item.members:
            atomNames.append(member.atomName)
            resLabel = member.resLabel
            seqCode = member.seqCode

        if not resLabel:
          print "  Error: cannot write RDC constraint (no resLabel)."

        else:
          protonName = constraint.getProtonName(resLabel,atomNames)

          if protonName:

            fout.write("%4d %-4s %-5s %-7.2f" % (seqCode,
                                                 resLabel,
	  	                                           protonName,
		                                             constraint.value))
                                                 
            if constraint.error:
              fout.write("%6.2f" % (constraint.error))

            if constraint.energyCst != 1.0:
              fout.write("%9.1f" % (constraint.energyCst))
        
            fout.write(self.newline)
      
      else:
        for item in constraint.items:
          for member in item.members:
            fout.write("  %4d  %-4s %-5s  " % (member.seqCode,
                                              member.resLabel,
	  	                                        member.atomName,
		                                          ))
                                           
        fout.write("   %8.3f %7.3f %8.3f" % (constraint.value,constraint.error,constraint.energyCst))
          
        if constraint.orientation != 0:
          fout.write("   %d" % (constraint.orientation))
        
        fout.write(self.newline)

    fout.close()

class DyanaRdcConstraint:

  def __init__(self,Id,cyanaLib):
    
    self.Id = returnInt(Id)
    self.items = []
    
    self.cyanaLib = cyanaLib
    
  def setRdcData(self,value,error=0.00, energyCst=1.0, orientation=0):
  
    self.error = returnFloat(error)
    self.value = returnFloat(value)
    self.energyCst = returnFloat(energyCst)
    
    self.orientation = returnInt(orientation)
      
  def setAtomMembersSingle(self,chainCode,seqCode,resLabel,refAtomName):
  
    refAtom = self.cyanaLib.findAtom(resLabel,refAtomName)
    
    seqCode = returnInt(seqCode)
    
    if refAtom:

      if refAtom.bondedAtomSerials.count(0) != 3:
        print "  Error: invalid single atom %s (%s %s, chain '%s'). No or multiple bonded atoms" % (refAtomName,resLabel,seqCode,chainCode)
      
      else:
      
        self.items.append(DyanaConstraintItem())
       
        self.items[-1].members.append(DyanaConstraintMember(chainCode,seqCode,resLabel,refAtomName))
  
        atomSerial = refAtom.bondedAtomSerials[0]     
        bondedAtom = self.cyanaLib.findAtomBySerial(resLabel,atomSerial)
  
        self.items[-1].members.append(DyanaConstraintMember(chainCode,seqCode,resLabel,bondedAtom.name))
    
    else:
      
      print "  Error: cannot write RDC restraint, atom %s for residue %s not found." % (refAtomName,resLabel)

  def setAtomMembersDouble(self,chainCode,residueInfo1,residueInfo2):
    
    self.items.append(DyanaConstraintItem())
    
    for residueInfo in (residueInfo1,residueInfo2):
      (seqCode,resLabel,atomName) = residueInfo

      seqCode = returnInt(seqCode)
       
      self.items[-1].members.append(DyanaConstraintMember(chainCode,seqCode,resLabel,atomName))
        
  def getProtonName(self,resLabel,atomNames):
    
    for atomName in atomNames:
    
      atom = self.cyanaLib.findAtom(resLabel,atomName)

      if atom and atom.bondedAtomSerials.count(0) == 3:
          
        return atom.name  
          
    return None
 
