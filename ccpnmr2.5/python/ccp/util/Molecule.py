"""
======================COPYRIGHT/LICENSE START==========================

Molecule.py: Utility functions

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
import re
from memops.general.Implementation import ApiError
from ccp.api.molecule import MolSystem, Molecule
from ccp.api.nmr import Nmr
#from ccp.util.Assignment import deleteResonanceSet
#from ccp.general.Io import getChemComp

try:
  from memops.gui.MessageReporter import showWarning
except:
  from memops.universal.MessageReporter import showWarning

# Moved to more appropriate location.
from ccp.util.Assignment import findAtomSetResonances
from ccp.lib.MoleculeModify import makeMolecule, addMolResidues, makeLinearSequence
# The following should not be needed:
#from ccp.lib.MoleculeModify import _getLinearChemCompData as getLinearChemCompData
from ccp.lib.MoleculeModify import setMolResidueCcpCode, setMolResidueChemCompVar
from ccp.lib.MoleculeModify import nextChainCode, makeChain, renumberChainSeqCodes

# Used in Analysis v2 (only), but hardly worth copying over:
def makeChainCopy(chain):
  """Descrn: Make a duplicate of a molSystem chain, using the same underlying molecule
     Inputs: Ccp.MolSystem.Chain
     Output: Ccp.MolSystem.Chain
  """

  code = nextChainCode(chain.molSystem)
  newChain = MolSystem.Chain(chain.molSystem, code=code, molecule=chain.molecule)
  newChain.setDetails(chain.details)
  
  return newChain

# Not currently used, and ratehr strange. Rmoves ChainFragment from MOLECULE!!!
def deleteChainFragment(chainFragment):
  """Descrn: Remove a molSystem chain fragment by recreating the parent chain with fewer residues
     Inputs: Ccp.MolSystem.ChainFragment
     Output: None
  """

  # delete residues, molResidues, molResLinks and anything else?
  chain       = chainFragment.chain
  code        = chain.code
  molecule    = chain.molecule
  details     = chain.details
  molSystem   = chain.molSystem
  molResidues = [r.molResidue for r in chainFragment.residues]
  chain.delete()
  
  for molResidue in molResidues:
    for linkEnd in molResidue.molResLinkEnds:
      if linkEnd.molResLink:
        linkEnd.molResLink.delete()
    molResidue.delete() 

  chain = MolSystem.Chain(molSystem, code=code, molecule=molecule)
  chain.setDetails(details)

# Not currently used
def createCoordinatesForResidue(coordDict,residue,sourceName='euroCarbDb',chemCompCoordArchiveDir=None):
   
  """

  Input:
  
    coordDict           A dictionary (initially empty). Will eventually contain coordinates for all residues (on chain level) connected in a molecule.
    residue             The residue to create the coordinates (or coordinates for connected residues) for.
    sourceName          Source name for the chemComp coordinate system (ideal, pdb, euroCarbDb, ...)
    
  Action:
  
    Adds residues to coordDict. The coordDict[residue] then points to a dictionary with as key (chemAtomType,atomName,atomSubType), which gives the (x,y,z) coordinates.
    
  Output:
  
    None

  WARNING: This code only works if coordinates are available for the linkAtoms for a chemComp that are
           involved in defining the molecule the chain is based on!
             
  """
  
  from ccp.general.Io import getChemCompCoord
  from memops.universal.Geometry import superposeNewVectorsOnOld
  
  chemCompCoord = getChemCompCoord(residue.root,sourceName = sourceName,ccpCode = residue.ccpCode, molType = residue.molType,chemCompCoordArchiveDir=chemCompCoordArchiveDir)
  
  chemCompVar = residue.chemCompVar
  chemCompVarCoord = chemCompCoord.findFirstChemCompVarCoord(linking = chemCompVar.linking, descriptor = chemCompVar.descriptor)
  
  #
  # Set the coordinates for this residue, if required
  #
  
  if not coordDict.has_key(residue):
    coordDict[residue] = {}

    for chemAtomCoord in chemCompVarCoord.chemAtomCoords:
    
      coordDict[residue][(chemAtomCoord.chemAtom.className,chemAtomCoord.name,chemAtomCoord.subType)] = (chemAtomCoord.x,chemAtomCoord.y,chemAtomCoord.z)
  
  #
  # Now look for connected residues
  #
  
  addedResidues = []
    
  molResidue = residue.molResidue
    
  for mrle in molResidue.molResLinkEnds:
    
    #
    # Get info for this residue - always use the recalculated coordinates in coordDict!
    #
    
    linkEnd = mrle.linkEnd
    
    #
    # Get info for other end
    #
        
    mrl = mrle.molResLink    
    mrles = mrl.sortedMolResLinkEnds()
    
    otherMrle = mrles[not mrles.index(mrle)]
    
    otherMolResidue = otherMrle.molResidue
    otherResidue = residue.chain.findFirstResidue(molResidue = otherMolResidue)
    
    #
    # If this one is already done, then ignore!
    #
    
    if coordDict.has_key(otherResidue):
      continue
    else:
      addedResidues.append(otherResidue)
      
    #
    # Now get coordinate info for 'originating' residue
    #
    
    boundChemAtom = linkEnd.boundChemAtom
    boundLinkAtom = linkEnd.boundLinkAtom
    
    boundChemAtomKey = (boundChemAtom.className,boundChemAtom.name,boundChemAtom.subType)
    boundChemAtomCoords = coordDict[residue][boundChemAtomKey]
    
    boundLinkAtomKey = (boundLinkAtom.className,boundLinkAtom.name,boundLinkAtom.subType)
    boundLinkAtomCoords = coordDict[residue][boundLinkAtomKey]

    #
    # Now get coordinate info for other residue
    #
    
    otherChemCompCoord = getChemCompCoord(otherResidue.root,sourceName = sourceName,ccpCode = otherResidue.ccpCode, molType = otherResidue.molType,chemCompCoordArchiveDir=chemCompCoordArchiveDir)
     
    otherChemCompVar = otherMolResidue.chemCompVar    
    otherChemCompVarCoord = otherChemCompCoord.findFirstChemCompVarCoord(linking = otherChemCompVar.linking, descriptor = otherChemCompVar.descriptor)
    
    otherLinkEnd = otherMrle.linkEnd
    
    otherChemAtom = otherLinkEnd.boundChemAtom
    otherLinkAtom = otherLinkEnd.boundLinkAtom
    
    otherChemAtomCoord = otherChemCompVarCoord.findFirstChemAtomCoord(name = otherChemAtom.name, subType = otherChemAtom.subType)
    otherLinkAtomCoord = otherChemCompVarCoord.findFirstChemAtomCoord(name = otherLinkAtom.name, subType = otherLinkAtom.subType)
    
    newVectorCoordAtoms = [otherLinkAtomCoord,otherChemAtomCoord]
    newVectors = [(otherLinkAtomCoord.x,otherLinkAtomCoord.y,otherLinkAtomCoord.z),(otherChemAtomCoord.x,otherChemAtomCoord.y,otherChemAtomCoord.z)]
    
    for chemAtomCoord in otherChemCompVarCoord.chemAtomCoords:
      
      if not chemAtomCoord in newVectorCoordAtoms:

        newVectorCoordAtoms.append(chemAtomCoord)
        newVectors.append((chemAtomCoord.x,chemAtomCoord.y,chemAtomCoord.z))      
    
    #
    # Now get the new coordinates, and set them...
    #
            
    oldVectors = (boundChemAtomCoords,boundLinkAtomCoords)
    
    transposedVectors = superposeNewVectorsOnOld(oldVectors,newVectors)
    
    coordDict[otherResidue] = {}

    # Ignore the link atom!
    for i in range(1,len(newVectorCoordAtoms)):
      chemAtomCoord = newVectorCoordAtoms[i]
      coords = tuple(transposedVectors[i])
      coordDict[otherResidue][(chemAtomCoord.chemAtom.className,chemAtomCoord.name,chemAtomCoord.subType)] = coords
      
  #
  # Warning: can't handle cyclic molecules!
  #
  
  for addedResidue in addedResidues:
    createCoordinatesForResidue(coordDict, addedResidue, sourceName = sourceName, chemCompCoordArchiveDir=chemCompCoordArchiveDir)
  
  
# Not currently used
def createMolStructureFromChemCompCoords(chain,coordSourceName='euroCarbDb',chemCompCoordArchiveDir=None):

  """
  Input:
  
    chain               CCPN Chain object
    coordSourceName     source name for the chemComp coordinate system (ideal, pdb, euroCarbDb, ...)
    
  Action:
  
    Creates a new StructureEnsemble with one Model. This Model has coordinates that are based on the
    ChemCompCoord data (with the relevant coordSourceName)
    
  Output:
  
    The new CCPN StructureEnsemble object.
    

  WARNING: This code only works if coordinates are available for the linkAtoms for a chemComp that are
           involved in defining the molecule the chain is based on!
           
  """
  
  coordDict = {}
  
  project = chain.root
  
  #
  # Find the residue with the starting coordinates
  #
  
  refResidue = None

  #
  # For carbohydrates, start at reducing end if possible
  #
  
  if chain.molecule.molType == 'carbohydrate':
    for residue in chain.residues:
      if not 'C1' in residue.linking:
        refResidue = residue
        break
  
  #
  # If no default found, just start with first residue.
  #
  
  if not refResidue:
    refResidue = chain.findFirstResidue(seqId = 1)
    if not refResidue:
      refResidue = chain.findFirstResidue()
  
  #
  # Create the coordinates for the refResidue
  #
  # This code will automatically create the coords for all connected residues!
  #
  
  createCoordinatesForResidue(coordDict,refResidue, sourceName = coordSourceName, chemCompCoordArchiveDir=chemCompCoordArchiveDir)

  #
  # Now create the actual StructureEnsemble coordinates from the dictionary...
  #
  
  for ensembleId in range(100):
    if not project.findFirstStructureEnsemble(molSystem = chain.molSystem, ensembleId = ensembleId):
      structureEnsemble = project.newStructureEnsemble(molSystem = chain.molSystem, ensembleId = ensembleId)
      break
  
  coordChain = structureEnsemble.newChain(code = chain.code)
  
  coordinates = []
  
  for residue in coordDict.keys():
    coordResidue = coordChain.newResidue(seqCode = residue.seqCode, seqId = residue.seqId, seqInsertCode = residue.seqInsertCode)
    
    for atomKey in coordDict[residue].keys():
    
      (atomType,atomName,atomSubType) = atomKey
      
      if atomType == 'LinkAtom':
        continue
    
      coordAtom = coordResidue.newAtom(name = atomName)

      coordTriplet = coordDict[residue][atomKey]
      coordinates.extend(coordTriplet)
  
  
  model = structureEnsemble.newModel()
  model.setSubmatrixData('coordinates', coordinates)
      
  return structureEnsemble

