
"""
======================COPYRIGHT/LICENSE START==========================

Mars.py: Part of the CcpNmr Analysis program

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
# Read data model, prepare temp MARS files, run MARS, interpret the output.

# REQUIREMENTS
#
# SS prediction

import os, re

reSub = re.sub

from os import path, mkdir, environ

MARS_DIR = os.environ.get('MARSHOME', '.')

INFILE_NAMES = {'config':'%s_mars.inp', 'shifts':'%s_cs.tab', 
                'fixcon':'%s_fix_con.tab', 'fixass':'%s_fix_ass.tab', 
                'sequence':'%s_fasta.tab', 'cyss':'%s_cyss.tab',
                'predcs':'%s_predcs.tab'}

ATOM_NAMES = ['H','N','C','CA','CB']

DEFAULT_OPTIONS = {'pdb':'0',
                   'resolution':'NO',
                   'pdbName':'NO',
                   'tensor':'NO',
                   'nIter':'NO',
                   'dObsExh':'NO',
                   'dcTab':'NO',}

RELIABILITY = {'H':1.0,'M':0.9,'L':0.2}


def tesMars(argServer):

  p = argServer.getProject()

  chain = p.findFirstMolSystem().findFirstChain()
  shiftList = p.currentNmrProject.findFirstMeasurementList(className='ShiftList', serial=2)

  print runMars(shiftList, chain)
  
def runMars(shiftList, chain, fragSize=5,
            cutoffCO=0.25, cutoffCA=0.2,
            cutoffCB=0.5, cutoffHA=0.25,
            cutoffN=0.5, cutoffHN=0.3,
            useConnections=True,
            useAssignment=True,
            isDeuterated=False,
            isUnfolded=False,
            pHvalue=7.0,
            corrDisulfides=False,
            syntheticShiftList=None):

  
  if not shiftList:
    print 'MARS cannot run - no shift list selected'
    return
    
  if not chain:
    print 'MARS cannot run - no chain selected'
    return
    
  marsDir, configFile = writeMarsInput(shiftList, chain, fragSize,
                              cutoffCO, cutoffCA, cutoffCB, cutoffHA,
                              cutoffN, cutoffHN, useConnections, useAssignment,
                              isDeuterated, isUnfolded, pHvalue, 
                              corrDisulfides, syntheticShiftList)
  pwd = os.getcwd() 
  
  os.chdir(marsDir) 
  
  marsExec = path.join(MARS_DIR, 'runmars')
  
  os.system('%s %s' % (marsExec, configFile))   
 
  os.chdir(pwd) 

  outFile = os.path.join(marsDir, 'assignment_AA.out')
  
  assignDict = readMarsOutput(outFile, chain)
   
  return  assignDict
   
def readMarsOutput(outFile, chain):

  assignDict = {}
  residues = chain.sortedResidues()
  nmrProject = chain.root.currentNmrProject

  file = open(outFile)
  
  results = []
  
  line = file.readline()
  while line:
    result = []
    data = line.split()
    data.pop(0)
    
    while data:
      serial =  int(data.pop(0))
      score =  RELIABILITY.get(data.pop(0)[1:-1], 0.0)
      result.append((serial, score))
  
    results.append(result)
    line = file.readline()
    
  file.close()
  
  for i, residue in enumerate(residues):
    result = results[i]
   
    if result:
      spinSystems = []
      scores = []
    
      for serial, score in result:
        spinSystem = nmrProject.findFirstResonanceGroup(serial=serial)
        
        if spinSystem:
          spinSystems.append(spinSystem)
          scores.append(score)

      if spinSystems:
        assignDict[residue] = (scores[0], spinSystems)

  return assignDict

def writeMarsInput(shiftList, chain, fragSize=5, cutoffCO=0.25, cutoffCA=0.2,
                   cutoffCB=0.4, cutoffHA=0.25, cutoffN=0.5, cutoffHN=0.3,
                   useConnections=True, useAssignment=True,
                   isDeuterated=False, isUnfolded=False,
                   pHvalue=7.0, corrDisulfides=False,
                   syntheticShiftList=None):

  from ccpnmr.analysis.core.AssignmentBasic import findConnectedSpinSystems
  from ccpnmr.analysis.core.AssignmentBasic import getSpinSystemResidues
  from ccpnmr.analysis.core.AssignmentBasic import getSpinSystemChemComps
  from ccpnmr.analysis.core.MoleculeBasic import getLinkedResidue
  from ccpnmr.analysis.wrappers.Psipred import psipredCcpn
  
  memopsRoot = chain.root

  residueNums = {}
  for i, residue in enumerate(chain.sortedResidues()):
    residueNums[residue] = i+1

  dataRepository = memopsRoot.findFirstRepository(name='userData')

  projPath = dataRepository.url.dataLocation

  prefix = '%s_%s' % (memopsRoot.name, chain.molecule.name)
  prefix = reSub('\s+',r'_',prefix)
  
  # # # # # Setup directory  # # # #
  # - use project folder - always writable

  marsDir = path.join(projPath, 'mars')
  if not path.exists(marsDir):
    mkdir(marsDir)
  
  # # # # Setup file names # # # #
  
  fileNames = {}
  for tag,text in INFILE_NAMES.items():
    fileNames[tag] = path.join(marsDir, text % prefix)
  
  # # # # Write sequence # # # #
  
  sequenceFile = open(fileNames['sequence'], 'w')
    
  seq = chain.molecule.stdSeqString.upper()
  #seq = reSub('\*','X', seq)
  seq = reSub('(\S{60})(\S)',r'\1\n\2',seq)

  sequenceFile.write('> %s\n%s\n' % (chain.molecule.name, seq))
  sequenceFile.close()
  
  # # # # Write PSIPRED secondary struct predictions # # # # 
  
  ssDict, ssFile, psipredFile = psipredCcpn(chain)
  
  # # # # Write chemical shifts # # # #
  
  shiftsFile = open(fileNames['shifts'], 'w')
  
  nmrProject = shiftList.topObject
  
  resonanceTable = {}
  usedAtomNames = set()
  
  # Extract resonances from data model
  
  # first pass - get resonanceDict, prevResidues, 
  # and decide which are real spin systems
  resonanceDicts = {}
  prevSpinSystems = {}
  rootSpinSystems = []
  for spinSystem in nmrProject.sortedResonanceGroups():
  
    # check for incompatible assignment
    assignedResidues = getSpinSystemResidues(spinSystem)
    if assignedResidues:
      assignedChains = set(x.chain for x in assignedResidues)
      if chain not in assignedChains:
        # spin system is assigned to a different chain. Skip it
        #print '### skipping', spinSystem, assignedChains
        continue
    
    resonanceDict = resonanceDicts[spinSystem] = {}
      
    for resonance in spinSystem.resonances:
      if not resonance.findFirstShift(parentList=shiftList):
        continue
    
      atomType = ','.join(resonance.assignNames)
      if atomType in ATOM_NAMES:
        resonanceDict[atomType] = resonance
    
    prevInSeq = findConnectedSpinSystems(spinSystem, -1)
    
   # print '###', 'H' in resonanceDict and 'N' in resonanceDict, spinSystem.serial, assignedResidues, prevInSeq
    
    if ((useAssignment and assignedResidues)
        or (useConnections and prevInSeq)
        or ('H' in resonanceDict and 'N' in resonanceDict)):
      # Treat as real (root) spin system
      prevSpinSystems[spinSystem] = prevInSeq
      rootSpinSystems.append(spinSystem)
  
  # second pass - write out table
  for spinSystem in rootSpinSystems:
    resonanceDict = resonanceDicts[spinSystem]
    resonanceDictPrev = {}

    columnData = []
    for atomType in ATOM_NAMES:
      columnData.append( (spinSystem, atomType) )
    
    # get non-real previous spin systems
    prevInSeq = [x for x in prevSpinSystems[spinSystem] 
                 if x not in prevSpinSystems]
    if prevInSeq:
      if len(prevInSeq) > 1:
        print 'WARNING, spin system %s has more than one i-1 pseudo-neighbour'
      prevSpinSystem = prevInSeq[0]
    
      for resonance in prevSpinSystem.resonances:
        if not resonance.findFirstShift(parentList=shiftList):
          continue
        
        atomType = ','.join(resonance.assignNames)
        if atomType in ATOM_NAMES:
          resonanceDictPrev[atomType] = resonance
          
      for atomType in ATOM_NAMES[2:]:
        columnData.append( (prevSpinSystem, atomType) )
    
    resonanceTable[spinSystem] = {}  
    for ss, atomName in columnData:
      if atomName == 'C':
        colName = 'CO'
      else:
        colName = atomName  
    
      if ss is spinSystem:
        resonance = resonanceDict.get(atomName)
      else:
        colName += '-1'
        resonance = resonanceDictPrev.get(atomName)
    
      if resonance:
        usedAtomNames.add(colName)
    
      resonanceTable[spinSystem][colName] = resonance 
    
  # Write header
  
  head = ' '*5
  colNames = list(usedAtomNames)
  colNames.sort()
  
  for name in colNames:
    head += ' %10.10s' % name
  
  shiftsFile.write('# Generated from CcpNmr Analysis\n')
  shiftsFile.write(head+'\n')
  
  # Write the text file
  
  for spinSystem in rootSpinSystems:
    line = '%5d' % spinSystem.serial
    colDict = resonanceTable[spinSystem]
    
    for colName in colNames:
      resonance = colDict.get(colName)
      
      text = '-'
      if resonance:
        shift = resonance.findFirstShift(parentList=shiftList)
      
        if shift:
          text = '%3.3f' % shift.value
        

      line += ' %10.10s' % text

    shiftsFile.write(line+'\n')
   
  shiftsFile.close()
  
  # # # # Write connections # # # #
  
  fixconFile = open(fileNames['fixcon'], 'w')
  
  if useConnections:
    for spinSystem in rootSpinSystems:
      # get real previous spin systems
      prevInSeq = [x for x in prevSpinSystems[spinSystem] 
                 if x in prevSpinSystems]
      if prevInSeq:
        if len(prevInSeq) > 1:
          print 'WARNING, spin system %s has more than one i-1 pseudo-neighbour'
        prevSpinSystem = prevInSeq[0]
      
        line = '%d\t%d\n' % (spinSystem.serial, prevSpinSystem.serial)
        fixconFile.write(line)
        
  fixconFile.close()
  
  """ S SpinSystem is real if
  1) it has a -1 neighbour
  2) if is has N and H
  3) unless a) there are N and H in the 
  """
  
  
  # # # # Write assignment # # # #
  # Rewritten RHF Jan 10 to handle more possible cases
  
  fixassFile = open(fileNames['fixass'], 'w')
  
  if useAssignment:
    for spinSystem in rootSpinSystems:
      resTypes = set()
      
      # get residue numbers to set
      ll = []
      resChemComps = set()
      for residue in getSpinSystemResidues(spinSystem):
        resChemComps.add(residue.molResidue.chemComp)
        num = residueNums.get(residue)
        if num is not None:
          ll.append(num)
      resData = [str(x) for x in sorted(ll)]
      
      # get one-letter codes to set
      chemComps = getSpinSystemChemComps(spinSystem)
      #if chemComps and not chemComps.issuperset(resChemComps):
      if chemComps and not resChemComps.issubset(chemComps):
        # there are chemComps that restrict the assignment
        resTypes = set(x.code1Letter or '*' for x in chemComps)
      else:
        resTypes = []
      
      # write out data
      if resData or resTypes:
        resData.append(''.join(resTypes))
        fixassFile .write('%d %s\n' % (spinSystem.serial, ' '.join(resData))) 
  
  fixassFile.close()
  
  # # # # Write cyssTab file (optional) # # # #
  
  cyss = [x for x in chain.sortedResidues()
          if x.ccpCode=='Cys' and 'link:SG' in x.descriptor]
  if corrDisulfides and cyss:
    cyssFile = open(fileNames['cyss'], 'w')
    for cys in cyss:
      cyssFile.write( '%s\n' % residueNums.get(cys))
    # TODO FILE writing code here
    cyssFile.close()
  else:
    fileNames['cyss'] = 'NO'
  
  # # # # Write synthetic shifts file (optional) # # # #
  
  if syntheticShiftList is None:
    fileNames['predcs'] = 'NO'
  
  else:
    predcsFile = open(fileNames['predcss'], 'w')
    # TODO FILE writing code here
    predcsFile.close()
  
  
  # # # # Write Options file # # # #
  
  configFile = open(fileNames['config'], 'w')
  
  configFile.write('fragSize:   \t%d\n' % fragSize)
  configFile.write('cutoffCO:   \t%.3f\n' % cutoffCO)
  configFile.write('cutoffCA:   \t%.3f\n' % cutoffCA)
  configFile.write('cutoffCB:   \t%.3f\n' % cutoffCB)
  configFile.write('cutoffHA:   \t%.3f\n' % cutoffHA)
  #configFile.write('cutoffN:   \t%.3f\n' % cutoffN)
  #configFile.write('cutoffHN:   \t%.3f\n' % cutoffHN)
  configFile.write('csTab:      \t%s\n' % fileNames['shifts'])
  configFile.write('fixConn:    \t%s\n' % fileNames['fixcon'])
  configFile.write('fixAss:     \t%s\n' % fileNames['fixass'])
  configFile.write('deuterated:   \t%s\n' % int(bool(isDeuterated)))
  configFile.write('sequence:   \t%s\n' % fileNames['sequence'])
  configFile.write('secondary:  \t%s\n' % psipredFile)
  #configFile.write('unfolded:   \t%s\n' % int(bool(isUnfolded)))
  #configFile.write('pHval???:   \t%.3f\n' % pHvalue)
  #
  
  for key in DEFAULT_OPTIONS:
    text = '%-12.12s' % (key+':')
    configFile.write('%s\t%s\n' % (text, DEFAULT_OPTIONS[key]))
  
  configFile.write('cyssTab:     \t%s\n' % fileNames['cyss'])
  #configFile.write('predcsTab:     \t%s\n' % fileNames['predcs'])
  
  configFile.close()
  
  # Hand back config file name
  
  return marsDir, fileNames['config']
