
"""
======================COPYRIGHT/LICENSE START==========================

SetBfactorFromShiftDiff.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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
ISOTOPE_SCALE = {'1H':100.0,'15N':20.0,'13C':10.0,}

def setBFactorFromShiftDiff(argServer):

  sdList = argServer.getMeasurementList('ShiftDifferenceList')  

  if sdList is None:
    msg = 'No Shift Difference list, cannot continue.\n'
    msg += 'Make shift difference lists by comparing peak lists'
    msg += ' or shift lists at Menu::DataAnalysis::Shift Differences'
    argServer.showWarning(msg)
    return

  structure = argServer.getStructure()

  if structure is None:
    msg = 'No structure, cannot continue.\n'
    msg += 'Load a PDB structure via Menu::Structure::Structures:[Import]'
    argServer.showWarning(msg)
    return

  question = 'Do you want to spread the B factor to all atoms in a residue?'
  if argServer.askYesNo(question):
    wholeResidue = True
  else:
    wholeResidue = False  

  atomsDict = {}
  for measurement in sdList.measurements:
    resonance = measurement.resonance
    scale = ISOTOPE_SCALE.get(resonance.isotopeCode, 10.0)
    resonanceSet = resonance.resonanceSet
    
    if resonanceSet:
      for atom in resonanceSet.findFirstAtomSet().atoms:
         atomsDict[atom] = measurement.value * scale

  for chain in structure.coordChains:
    for residue in chain.residues:
      for atom in residue.atoms:
        value = atomsDict.get(atom.atom)
      
        if value is not None:
 
          if wholeResidue:
            for atom2 in residue.atoms:
              for coord in atom2.coords:
                coord.bFactor = value
 
          else:
            for coord in atom.coords:
              coord.bFactor = value

  strucId = '%s:%d' % (structure.molSystem.code, structure.ensembleId)
  msg = 'Shift differences from list %d stored as bFactors in structure %s.\n'
  msg += 'B-factor values will be written out when you export the structure as a PDB file.'
  argServer.showInfo(msg % (sdList.serial,strucId))
