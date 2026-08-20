#!/usr/bin/python

"""
======================COPYRIGHT/LICENSE START==========================

chemCompsIO.py: I/O for mol format atom/molecule information file

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
from memops.universal.Util import returnInt, returnFloat
from memops.universal.Io import getTopDirectory

# Import header reader
from ccp.format.mol.generalIO import MolGenericFile

from ccp.format.general.Constants import defaultSeqInsertCode
from ccp.format.general.Util import standardNucleusName

#####################
# Class definitions #
#####################

class MolChemCompFile(MolGenericFile):

  def initialize(self):
  
    self.chemComps = []

  def read(self,verbose = 0):

    if verbose == 1:
      print "Reading %s chemComp file %s" % (self.format,self.name)
      
    
    #
    # Initialize locals
    #
    
    chemCompRead = 0
    atomSerial = 0

    #
    # Read all information...
    #

    fin = open(self.name, 'rU')
    line = fin.readline()

    while line:
    
      #
      # Warning: atomList, stext blocks NOT treated!!
      #
      # Comments at end are also ignored.
      #

      if not chemCompRead and self.patt[self.format + 'Counts'].search(line):
        
        #
        # Counts line
        #
        
        chemCompRead = 1
        
        numAtoms = returnInt(line[0:3],verbose = 0)
        numBonds = returnInt(line[3:6],verbose = 0)
        numAtomLists = returnInt(line[6:9],verbose = 0)
        isChiral = returnInt(line[12:15],verbose = 0)
        stext = returnInt(line[15:18],verbose = 0)
        
        self.chemComps.append(MolChemComp(self,numAtoms,numBonds,isChiral))

      elif chemCompRead and self.patt[self.format + 'Atoms'].search(line):
      
        atomSerial += 1

        #
        # Atom lines
        #

        x = line[0:10]
        y = line[10:20]
        z = line[20:30]

        atomType = line[31:34]
        massDiff = line[34:36]
        charge = line[36:39]
        Hcount = line[42:45]
        stereoCare = line[45:48]
        valence = line[48:51]
        
        # NOT USED...
        stereoParity = line[39:42]
        reactionAtomAtomMap = returnInt(line[60:63],verbose = 0)
        reactionInvRetFlag = returnInt(line[63:66],verbose = 0)
        reactionChange = returnInt(line[66:69],verbose = 0)
        
        self.chemComps[-1].atoms.append(MolAtom(self.chemComps[-1],atomSerial,x,y,z,atomType,massDiff,charge,Hcount,stereoCare,valence))

      elif chemCompRead and self.patt[self.format + 'Bonds'].search(line):

        atomSerial1 = line[0:3]
        atomSerial2 = line[3:6]
        bondType = line[6:9]
        bondStereo = line[9:12]
        bondTopology = line[15:18]
        
        # NOT USED...
        reactionCenter = line[18:21]

        self.chemComps[-1].bonds.append(MolBond(self.chemComps[-1],atomSerial1,atomSerial2,bondType,bondStereo,bondTopology))

      elif chemCompRead and line[:6] == 'M  END':

        if self.chemComps[-1].numAtoms:
          print "  Warning: there are %d atoms not handled!" % self.chemComps[-1].numAtoms
          
        if self.chemComps[-1].numBonds:
          print "  Warning: there are %d bonds not handled!" % self.chemComps[-1].numBonds

        chemCompRead = 0
        atomSerial = 0

      line = fin.readline()

    fin.close()

class MolChemComp:

  def __init__(self,parent,numAtoms,numBonds,isChiral):

    self.parent = parent

    self.numAtoms = numAtoms
    self.numBonds = numBonds
    self.isChiral = isChiral

    self.atoms = []
    self.bonds = []

class MolAtom:

  def __init__(self,parent,atomSerial,x,y,z,atomType,massDiff,charge,Hcount,stereoCare,valence):
  
    self.parent = parent
    self.parent.numAtoms -= 1
    
    self.serial = atomSerial
    
    self.x = returnFloat(x,verbose = 0)
    self.y = returnFloat(y,verbose = 0)
    self.z = returnFloat(z,verbose = 0)
    
    """
    entry in periodic table or L for atom list, 
    A, Q, * for unspecified atom, and LP for 
    lone pair, or R# for Rgroup label
    """
    
    atomType = string.strip(atomType)
    
    if atomType in ['A','Q','*']:
      atomType = None
    elif atomType == 'LP':
      atomType = None
      #atomType = 'lone pair' TODO TODO can't handle this yet
    
    self.atomType = atomType

    """
    -3, -2, -1, 0, 1, 2, 3, 4 
    (0 if value beyond these limits)
    """
    
    self.massDiff = returnInt(massDiff,verbose = 0)
    
    """
    0 = uncharged or value other than these,
    1 = +3, 2 = +2, 3 = +1,  
    4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
    """
    
    chargeCode = returnInt(charge,verbose = 0)
    self.charge = [0,3,2,1,'dr',-1,-2,-3][chargeCode]    
    
    """
    1 = H0, 2 = H1, 3 = H2, 4 = H3, 5 = H4
    """
    self.Hcount = returnInt(Hcount,verbose = 0) - 1
        
        
    """
    0 = ignore stereo configuration of this 
    double bond atom, 1 = stereo 
    configuration of double bond atom 
    must match        
    """
    self.stereoCare = returnInt(stereoCare,verbose = 0)
        
    """
    0 = no marking (default) 
    (1 to 14) = (1 to 14) 15 = zero valence
    """
    self.valence = returnInt(valence,verbose = 0)

class MolBond:

  def __init__(self,parent,atomSerial1,atomSerial2,bondType,bondStereo,bondTopology):
  
    self.parent = parent
    self.parent.numBonds -= 1
        
    self.atomSerial1 = returnInt(atomSerial1,verbose = 0)
    self.atom1 = self.parent.atoms[self.atomSerial1-1]
    
    self.atomSerial2 = returnInt(atomSerial2,verbose = 0)
    self.atom2 = self.parent.atoms[self.atomSerial2-1]
    
    
    """
    1 = Single, 2 = Double, 
    3 = Triple, 4 = Aromatic, 
    5 = Single or Double, 
    6 = Single or Aromatic, 
    7 = Double or Aromatic, 8 = Any
    """

    bondTypeCode = returnInt(bondType,verbose = 0)
    self.bondType = ['void','single','double','triple','aromatic','single/double','single/aromatic','double/aromatic','unknown'][bondTypeCode]

    """        
    Single bonds: 0 = not stereo, 
    1 = Up, 4 = Either, 
    6 = Down,
    
    Double bonds: 0 = Use x-, y-, z-coords 
    from atom block to determine cis or trans, 
    3 = Cis or trans (either) double bond
    """
    
    bondStereoCode = returnInt(bondStereo,verbose = 0)
    if self.bondType == 'single':
      self.bondStereo = [None,'up',None,None,'either',None,'down'][bondStereoCode]
    elif self.bondType == 'double':
      self.bondStereo = ['coordinates',None,None,'either'][bondStereoCode]
        
    """        
    0 = Either, 1 = Ring, 2 = Chain        
    """       
    
    bondTopologyCode = returnInt(bondTopology,verbose = 0)
    self.bondTopology = ['either','ring','chain'][bondTopologyCode]
            
"""
FILE INFORMATION:

COUNTS line:

aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
Where::
aaa    = number of atoms (current max 255)*    [Generic]
bbb    = number of bonds (current max 255)*    [Generic]
lll    = number of atom lists (max 30)*     [Query]
fff    = (obsolete)
ccc    = chiral flag:  0=not chiral, 1=chiral    [Generic]
sss    = number of stext entries    [ISIS/Desktop]
xxx    = (obsolete)
rrr    = (obsolete)
ppp    = (obsolete)
iii    = (obsolete)
mmm    = number of lines of additional  properties,     [Generic]
including the M  END line. No longer 
supported, the default is set to 999.
        

ATOM lines:

xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee      

x y z    atom coordinates
aaa    atom symbol
dd    mass difference
ccc    charge
sss    atom stereo parity
hhh    hydrogen count + 1
bbb    stereo care box
vvv    valence
HHH    H0 designator
rrr    Not used
iii    Not used
mmm    atom-atom mapping 
number
nnn    inversion/retention flag
eee    exact change flag

BOND lines:

111222tttsssxxxrrrccc

111    first atom number
222    second atom number
ttt    bond type
sss    bond stereo
xxx    not used
rrr    bond topology
ccc    reacting center status

"""

###################
# Main of program #
###################

if __name__ == "__main__":  

  files = [
           '../../reference/mol/example.mol',
           '../../reference/mol/cyclic.mol',
           '../../reference/mol/mcys.mol'
          ]

  
  for inFile in files:
    
    molFile = MolChemCompFile(os.path.join(getTopDirectory(),inFile))

    molFile.read(verbose = 1)
    
    for chemComp in molFile.chemComps:
      for atom in chemComp.atoms:
        print atom.atomType,
      print
      for bond in chemComp.bonds:
        print bond.bondType
        
    #(path,baseName) = os.path.split(inFile)
    
    #molFile.name = ("local/%s" % baseName)
    
    #molFile.write(verbose = 1)
