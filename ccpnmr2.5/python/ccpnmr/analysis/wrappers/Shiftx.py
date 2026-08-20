"""
======================COPYRIGHT/LICENSE START==========================

Shiftx.py: Part of the CcpNmr Analysis program

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
from os import path, mkdir, system, environ

from memops.universal.Io import getTopDirectory

rootDir = path.split(getTopDirectory())[0]

defaultDir = path.join(rootDir,'other','shiftx')

SHIFTX_DIR = environ.get('SHIFTX_DIR', defaultDir)


def testShiftxMacro(argServer):

  structure = argServer.getStructure()
  
  shiftList = shiftx(structure)


def shiftx(structure, atomType=None):

  from ccpnmr.analysis.core.StructureBasic import makePdbFromStructure

  from ccpnmr.analysis.core.MoleculeBasic import makeResidueAtomSets, DEFAULT_ISOTOPES
  
  from ccpnmr.analysis.core.AssignmentBasic import assignAtomsToRes

  memopsRoot = structure.root
  nmrProject = memopsRoot.currentNmrProject

  dataRepository = memopsRoot.findFirstRepository(name='userData')
  projPath = dataRepository.url.dataLocation

  tempDir = path.join(projPath, 'shiftx')
  if not path.exists(tempDir):
    mkdir(tempDir)

  pdbFile = path.join(tempDir, 'shiftxInput.pdb')
  outFile = path.join(tempDir, 'shiftxOutput.out')
  exeFile = path.join(SHIFTX_DIR, 'shiftx')

  chainDict = {}
  shiftList = None
  
  for model in structure.models:
    makePdbFromStructure(pdbFile, structure, model=model, useOxt=True)
  
    for coordChain in structure.sortedCoordChains():
      chain = coordChain.chain
  
      print 'CCPN SHIFTX executing for chain %s' % chain.code
  
      cmd = '%s 1%s %s %s' % (exeFile, chain.code, pdbFile, outFile)
  
      system(cmd)
  
      fileHandle = open(outFile, 'r')
      lines = fileHandle.readlines()
   
      if chain not in chainDict:
        chainDict[chain] = shiftData = {}
      else:
        shiftData = chainDict[chain]
   
      for line in lines:
        data = line.split()
   
        if not data:
            continue
  
        first = data[0]
  
        if first.upper() == 'NUM':
            atomNames = data[2:]
   
        elif first[0] == '-':
            continue
   
        elif first[0] == '*':
            continue
   
        elif first.isdigit(): # What to do with alt codes and inaccuracies? *10, 12A etc
            seqCode = int(first)
            resCode = data[1]
   
            if seqCode not in shiftData:
              shiftData[seqCode] = shiftDict = {}
	    else:
	      shiftDict = shiftData[seqCode]  
   
            for i, atomName in enumerate(atomNames):
              if atomName not in shiftDict:
                shiftDict[atomName] = []
	      
	      shiftDict[atomName].append(float(data[i+2]))
  

  for chain in chainDict:
    print 'CCPN SHIFTX filling shift list for chain %s' % chain.code
    shiftData = chainDict[chain]

    if shiftData:
      if not shiftList:
  	msCode = structure.molSystem.code
  	eId = structure.ensembleId
  	details = 'SHIFTX prediction for %s ensemble %d' % (msCode, eId)
  	shiftList = nmrProject.newShiftList(name='SHIFTX',
  					    details=details,
  					    isSimulated=True) 
      
      for residue in chain.sortedResidues():
  	shiftDict = shiftData.get(residue.seqCode)
  	
  	if not shiftDict:
  	  continue
      
  	getAtom = residue.findFirstAtom     
      
  	for atomName in shiftDict:
  	  atom = getAtom(name=atomName)
  	  if not atom:
  	    # Gly HA -> HA2, Ala HB - HB2, Val HG1 - HG12
  	    atom = getAtom(name=atomName+'2')
  	  
  	  if not atom:
  	    continue
  
  	  ppms = shiftDict[atomName]
  	  ppm = sum(ppms)/float(len(ppms))
          
  	  if ppm == 0.0: # Could actually be a real value...
  	    continue
  
  	  atomSet = atom.atomSet
  	  if not atomSet:
  	    makeResidueAtomSets(residue)
  	    atomSet = atom.atomSet
  
  	  resonance = None
  	
          for resonanceSet in atomSet.resonanceSets:
  	    resonances = list(resonanceSet.resonances)
  	    
  	    if len(resonances) == 1:
  	      resonance = resonances[0]
  	    
  	  if not resonance:
  	    isotopeCode = DEFAULT_ISOTOPES[atom.chemAtom.chemElement.symbol]
  	    resonance = nmrProject.newResonance(isotopeCode=isotopeCode)
  	    assignAtomsToRes([atomSet,], resonance)
  
  	  if shiftList.findFirstMeasurement(resonance=resonance):
            print residue.seqCode, residue.ccpCode, atom.name
          else:  
  	    shift = shiftList.newShift(value=ppm, resonance=resonance)
  
  print 'CCPN SHIFTX done'

  return shiftList
