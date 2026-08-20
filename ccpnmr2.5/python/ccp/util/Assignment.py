"""
======================COPYRIGHT/LICENSE START==========================

Assignment.py: <write function here>

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
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

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================
"""

def isChainAssigned(chain):
  """Descrn: Determines if an NMR chain has any atoms in its residues that are assigned to resonances
     Inputs: Nmr.Chain
     Output: Boolean
  """

  for residue in chain.residues:
    for atom in residue.atoms:
      if atom.atomSet:
        if atom.atomSet.resonanceSets:
          return 1
      
  return 0

def isResidueAssigned(residue):
  """Descrn: Determines if an NMR residue has any atoms that are assigned to resonances
     Inputs: Nmr.Residue
     Output: Boolean
  """

  for atom in residue.atoms:
    if atom.atomSet:
      if atom.atomSet.resonanceSets:
        return 1
  
  return 0

def isAtomAssigned(atom):
  """Descrn: Determines if an atom is assigned to any resonances
     Inputs: Nmr.Atom
     Output: Boolean
  """

  if atom.atomSet:
    if atom.atomSet.resonanceSets:
      return 1
  
  return 0

def isPeakAssigned(peak, fully=0):
  """Descrn: Determines if a peak is assigned, ether fully or partially
     Inputs: Peak, Boolean
     Output: Boolean
  """

  n = 0
  for peakDim in peak.peakDims:
    if len(peakDim.peakDimContribs) > 0:
      n +=1
      
  if n == len(peak.peakDims):
    return 1
  elif n > 0:
    if fully:
      return 0
    else:
      return 1
  else:
    return 0

def makeAtomSetsGuiName(atomSets):
  """Descrn: Give the name of a given set of atoms for GUIs
     Inputs: List of Nmr.AtomSets
     Output: Word
  """

  resonanceSet = None
  resonanceSets = []
  atomNames = []
  for atomSet in atomSets:
    for atom in atomSet.atoms:
      atomNames.append(atom.name)
      
  tryName = atomNames[0]
  end = ''
  for atomName in atomNames:
    size = len(atomName)
    for i in range(size):
      j = size - i
      subName = atomName[0:j]
      if subName == tryName[0:j]:
        tryName = tryName[0:j]
        if i > 0:
          end = '*'
        break
        
  name = tryName[:1] + tryName[1:].lower() + end
  
  return name

def makeResidueGuiName(residue):
  """Descrn: Give the name of an NMR residue for GUIs
     Inputs: Nmr.Residue
     Output: Word (case adjusted Molecule.MolResidue.ccpCode)
  """

  tlc  = residue.molResidue.ccpCode
  tlc  = tlc[:1] + tlc[1:].lower()
  return tlc
            
def makeResonanceGuiName(resonance, fullName=1):
  """Descrn: Give the name of a resonance for GUIs. Either full identifier or just the atom type name (if any)
     Inputs: Nmr.Resonance, boolean
     Output: Word
  """

  (chain, residue, name) = getResonanceAtomTuple(resonance)
  if fullName:
    #return '%s%s%s' % (chain,residue,name)
    return '%s%s' % (residue,name)
  else:
    return name
  
def getResonanceResidue(resonance):
  """Descrn: Find the residue, if any, to which a resonance is assigned via atomSets etc
     Inputs: Nmr.Resonance
     Output: Nmr.Residue or None
  """

  residue = None
  if resonance.resonanceSet:
    residue = resonance.resonanceSet.findFirstAtomSet().findFirstAtom().residue
 
  return residue

def getResonanceAtomTuple(resonance):
  """Descrn: Give a tupe of string identifiers for a resonance indicating, chain, residue and atomic assignment
     Inputs: Nmr.Resonance
     Output: Tuple of Words (Nmr.Chain.code, Nmr.Residue identifier, Nmr.Resonance identifier)
  """

  if hasattr(resonance, 'resonanceSet'):
  
    chain = ''
    res   = ''
   
    if resonance.resonanceSet:
      resonanceSet = resonance.resonanceSet
      name    = makeAtomSetsGuiName( resonanceSet.atomSets )
      if len(resonanceSet.atomSets) > 1:
        i = list(resonanceSet.resonances).index(resonance)
        name = name[:-1] + chr(ord('a')+i)
        if len(resonanceSet.atomSets[0].atoms) > 1:
          name = name + '*'
        
      residue = getResonanceResidue(resonance)
      chain   = residue.chain.code
      res     =  str(residue.seqId) + makeResidueGuiName(residue)
    elif resonance.resonanceGroup and resonance.resonanceGroup.ccpCode:
      if resonance.resonanceGroup.residue:
        chain   = resonance.resonanceGroup.residue.chain.code
        res     = str(resonance.resonanceGroup.residue.seqId) + makeResidueGuiName(resonance.resonanceGroup.residue)
        name    = '[%d]' % resonance.serial
      else:
        tlc  = resonance.resonanceGroup.ccpCode
        res  = tlc[:1] + tlc[1:].lower()
        name =  '[%d]' % (resonance.serial)
        
    elif resonance.resonanceGroup:
      res  = '{%d}' % (resonance.resonanceGroup.serial)
      name = '[%d]' % resonance.serial
    elif resonance.name:
      name = '[%s]' % resonance.name
    else:
      name = '[%d]' % resonance.serial
    return (chain, res, name)

  else:
    assignNames = resonance.assignNames
    if assignNames:
      N = len(assignNames)
      if N > 1:
        name = ''
        for i in range(N):
          name = name + assignNames[i]
          if i < (N-1):
            name = name + '|'
      else:
        name = resonance.assignNames[0]
  
  return (chain, res, name)

def deleteResonanceSet(resonanceSet):
  """Descrn: Delete a resonance set whilst still retaining old links to resonances and atom sets
             Depricated given new API: retro-fit required
     Inputs: Nmr.ResonanceSet
     Output: None
  """

  resonanceSet.oldResonances = resonanceSet.resonances
  resonanceSet.oldAtomSets   = resonanceSet.atomSets
  resonanceSet.delete()
  
  
def findAtomSetResonances(atomSets):
  """Descrn: Find resonances assigned to given atom sets
     Inputs: List of Nmr.AtomSets
     Output: List of Nmr.Resonances
     
     # NBNB currently not used (Aug 2013, Rasmus), 
     # but may be relevant where Analysis.AtomSetMapping is not appropriate.
  """

  resonanceDict = {}
  for atomSet in atomSets:
    if hasattr(atomSet, 'resonanceSets'):
      for resonanceSet in atomSet.resonanceSets:
        for resonance in resonanceSet.resonances:
          resonanceDict[resonance.serial] = resonance
  
  if len(resonanceDict) > 0:
    return resonanceDict.values()
  else:
    return None
