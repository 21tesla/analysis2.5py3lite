LICENSE = """
======================COPYRIGHT/LICENSE START==========================

PdbIo.py: Basic library for Pdb file IO

Copyright (C) 2003-2013 Rasmus Fogh, Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

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
import os

from ccp.lib import DataConvertLib
from ccp.lib.MoleculeAlign import findMatchingChains
from ccp.lib.MoleculeModify import addMolResidues, makeMolecule
from ccp.lib.MoleculeModify import makeChain, nextChainCode
from memops.universal.Util import returnInt, returnFloat

try:
  from memops.gui.MessageReporter import showOkCancel, showYesNo
  from memops.gui.MessageReporter import showWarning #, showError
except:
  from memops.universal.MessageReporter import showOkCancel, showYesNo
  from memops.universal.MessageReporter import showWarning #, showError


######################################################
#
### Structure reading
#
######################################################

def getStructureFromFile(molSystem, filePath, fileType='rough', doWarnings=True):
  """Creates a structure ensemble belonging to a molecular system 
             by loading a PDB style file.
             Option to specify the type of file loaded (proper PDB, rough PDB or CNS style)
  .. describe:: Input
  
  MolSystem.MolSystem, String (File Name), String (file type)

  .. describe:: Output

  MolStructure.StructureEmsemble
  """
  
  coordDict = {}
  
  if fileType == 'true':  
    coordDict = _makeStructureDictFromPdb(filePath, fileType='pdb')
  elif fileType == 'cns':
    coordDict = _makeStructureDictFromPdb(filePath, fileType='cns')
  elif fileType == 'rough':
    coordDict = _makeStructureDictFromRoughPdb(filePath)
  else:
    try:
      coordDict = _makeStructureDictFromPdb(filePath, fileType='pdb')
    except:
      pass
  if not coordDict:
      coordDict = _makeStructureDictFromRoughPdb(filePath)

  ensemble = None

  if coordDict:
    ensemble = makeStructureEnsemble(coordDict, molSystem, doWarnings=doWarnings)

    if ensemble:
      models = ensemble.sortedModels()
      fileRoot, fileExt = os.path.splitext(os.path.split(filePath)[1])
    
      if len(models) > 1:
        for model in models:
          model.details = 'file:%s (MODEL %d)' % (filePath, model.serial)
          model.name = '%s %d' % (fileRoot, model.serial)

      else:
        models[0].details = 'file:%s' % (filePath)
        models[0].name = fileRoot
    
  return ensemble


def getStructureFromFiles(molSystem, pdbPaths, fileType=None, doWarnings=True):
  """ Load several PDB files, taking the first model from each file.
  If only one file, load all models.
  if fileType is given, tries filretype load, then type-rough load.
  If no fileType, tries pdb, then rough.
  """
  
  if not pdbPaths:
    print ('WARNING no files to read. Skipping')
  
  elif len(pdbPaths) == 1:
    return getStructureFromFile(molSystem, pdbPaths[0], fileType, doWarnings)
  
  else:
  
    data = {}
    modelNames = []
    
    if not fileType:
      fileType = 'pdb'
    elif fileType not in ('cns', 'pdb'):
      fileType = None
    
    nSkipped = 0
    modelNum = 0
    for pdbPath in pdbPaths:
      
      # Load data
      dd = None
      if fileType is not None:
        try:
          print '### trying %s' % fileType
          dd = _makeStructureDictFromPdb(pdbPath, fileType)
        except:
          pass
   
      if not dd:
        print '### trying rough'
        dd = _makeStructureDictFromRoughPdb(pdbPath)
   
      if dd:
        #
        data[modelNum] = dd[0]
   
        # append model name
        basename = os.path.basename(pdbPath)
        indx = basename.find('.pdb')
        if indx > 0:
          modelName = basename[:indx]
        else:
          modelName, ext = os.path.splitext(basename)
        modelNames.append(modelName)
        
        modelNum += 1
        
      else:
        nSkipped += 1
        print ('Error coordinate read failed. Skipping %s'
               % pdbPath)
        
    # make ensemble
    
    if data:
      if nSkipped:
        print 'WARNING, %s files skipped' % nSkipped
      return makeStructureEnsemble(data, molSystem, modelNames=modelNames,
                                   doWarnings=doWarnings)
    
    else:
      raise Exception("ERROR no PDB data read, %s files skipped" % nSkipped)
    
    
def makeStructureEnsemble(strucDict, molSystem, modelNames=None, 
                          doWarnings=False, checkTruncation=True,
                          fixAtomNames=True):
  """Makes structure ensemble from a structure dictionary
             [model][chainCode][resNum][atomName] = (x,y,z coords)
             in a mol system. Options to supress warnings and thus
             automatically link non-identical, but similar chains
             or make a new chain if none is found,
             NB chains and resides are those present in FIRST model
                atoms are taken from all models, 
                and missing atoms are given x=y=z = 0.0; occupancy=0.0 
  .. describe:: Input
  
  Dictionary, MolSystem.MolSystem, Boolean

  .. describe:: Output

  MolStructures.StructureEnsemble
  """
  
  # Rasmus Fogh 20/6/2011: Changed to new MolStructure model
  
  # Rasmus Fogh 4/5/2013:Added automatic  name fixing: 2HB -> HB2
  
  #
  # Set up, and create structure ensemble
  #
  
  ll = sorted(strucDict.keys())
  firstModelId = ll[0]

  project = molSystem.root
  
  ensembles = molSystem.sortedStructureEnsembles()
  if ensembles:
    eId = ensembles[-1].ensembleId + 1
  else:
    eId = 1

  ensemble = project.newStructureEnsemble(molSystem=molSystem, ensembleId=eId)
  
  #models = []
  #for m in strucDict.keys():
  #  models.append((m, ensemble.newModel()))
  
  #
  # record data list - for creating and processing records in reading order
  #
  
  recordDataList = []
  recordDataDict = {}
  coordResidues = {}
  allMsResidues = {}

  usedChains = []
  failedAtoms = []
  msAtomDict = {}
  atomCheckDict = {}
  systemAtoms = set()
  chainsDict = strucDict[firstModelId]
  
  # match or make chains
  chCodes = list(chainsDict.keys())
  chCodes.sort()
  nOrigChCodes = len(chCodes)
  for iChCode,chCode in enumerate(chCodes):
    resDict = chainsDict[chCode]
    seqIds  = resDict.keys()
    seqIds.sort()

    chemComps  = []
    ccpCodes   = []
    seqIdsGood = []
    
    # get list of good chemcomps
    molTypes = set()
    resNames = []
    for seqId in seqIds:
      ll = resDict[seqId].keys()
      resName   = resDict[seqId][ll[0]]['resName']
      atomNames = [tt[0] for tt in ll]
      chemComp  = DataConvertLib.getBestChemComp(project, resName, atomNames)
      resNames.append(resName)
      
      if chemComp:
        chemComps.append(chemComp)
        ccpCodes.append(chemComp.ccpCode)
        seqIdsGood.append(seqId)
        molTypes.add(chemComp.molType)

    if not seqIdsGood:
      msg  = 'Could not find any matching CCPN ChemComps for sequence: '
      msg += ' '.join(resNames)
      showWarning('Structure Creation Failed',msg)
      continue
    
    #get matching MolSystem chains
    seqIds = seqIdsGood
    msChains, mappings = findMatchingChains(molSystem, ccpCodes,
                                            excludeChains=usedChains,
                                            doWarning=doWarnings,
                                            molTypes=molTypes)
    
    nChains = len(msChains)
    if nChains == 1:
      msChain = msChains[0]
      mapping = mappings[0]
    
    elif nChains:
      from memops.gui.DataEntry import askString
       
      codes = [c.code for c in msChains]
      msg = 'Structure chain %s matches %d \nCCPN chains' % (chCode, nChains)
      msg += ' equally well: %s\n' % ( ' '.join(codes),)
      msg += 'Which CCPN chain should be linked to?'
      
      if chCode in codes:
        default = chCode
      else:
        default = codes[0]
      
      choice = None
      while choice not in codes:
        choice = askString('Query', msg, default) or ''
        choice = choice.strip()
      
      index = codes.index(choice)                        
      msChain = msChains[index]
      mapping = mappings[index]
      
    else:
      msChain = None
      
    if nChains:
      nRes = len(mapping)
      startPair = 0
      endPair = nRes-1
      
      for i in range(nRes):
        residueA, residueB = mapping[i]
        if residueA and residueB:
          startPair = i
          break

      for endTrunc in range(nRes):
        j = endPair - endTrunc
        residueA, residueB = mapping[j]
        if residueA and residueB:
          endPair = j
          break
      
      mismatch = []
      for i in range(startPair, endPair+1):
        if None in mapping[i]:
          mismatch.append(i)
          break
      
      if doWarnings:
        if checkTruncation and ((endTrunc > 10) or (startPair > 10)):
          msg  = 'Structure truncated (Start:%s, End:%s) ' % (startPair, endTrunc)
          msg += 'compared to existing chain. Really link to chain %s?' % msChain.code
          msg += ' (Otherwise a new chain will be made)' 
          if not showYesNo('Query', msg):
            msChain = None
        
        if mismatch and msChain is not None:
          msg  = 'Imperfect match: %s non-terminal mismatches, ' % len(mismatch)
          msg += 'compared to existing chain. Really link to chain %s?' % msChain.code
          msg += ' (Otherwise a new chain will be made)' 
          if not showYesNo('Query', msg):
            msChain = None
        
      else:
      
        if checkTruncation and (startPair or endTrunc):
          print ('WARNING, Structure truncated (Start:%s, End:%s) rel. to existing chain.' 
                 % (startPair, endTrunc))
        
        if mismatch:
          print 'WARNING, Imperfect match: %s non-terminal mismatches' % len(mismatch)
        
        
    if not msChain:
      # no matching chain - make new one
      # use existing chain code from PDB file if possible
      sysChainCode = chCode.strip()
      if not sysChainCode or molSystem.findFirstChain(code=sysChainCode):
        sysChainCode = nextChainCode(molSystem)
      
      if molTypes == set(('DNA',)) and iChCode <= nOrigChCodes:
        # check for special case: 
        # double stranded DNA with a single chain code
        oneLetterCodes = [cc.code1Letter for cc in chemComps]
        if None not in oneLetterCodes:
          # All ChemComps are (variants of) Std bases. 
          # Not sure certain, but this should mostly work.
          nCodes = len(oneLetterCodes)
          halfway, remainder = divmod(nCodes, 2)
          if not remainder:
            # even number of codes
            resMap = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G' }
            for ii in range(halfway):
              if oneLetterCodes[ii] != resMap[oneLetterCodes[-1-ii]]:
                break
            else:
              # the second half is the reverse complement of the first half.
              # treat as two separate DNA chains
              # Move second half to a new chain, and remove from this chain.
              newChCode = chCode
              while newChCode in chCodes:
                newChCode=chr(ord(newChCode)+1)
              newResDict = chainsDict[newChCode] = {}
              for ii in range(halfway, nCodes):
                iix = seqIds[ii]
                newResDict[iix] = resDict[iix]
                del resDict[iix]
              
              # put both chains on end of list for later (re)processing
              chCodes.append(chCode)
              chCodes.append(newChCode)
              continue
              
      codes = (molSystem.code, sysChainCode)
      msg  = 'Structure residue sequence (chain %s) not in molecular system. ' % chCode
      msg += 'Make new molecular system %s chain (%s) for this sequence?' % codes
      if not doWarnings or (not molSystem.chains) or (doWarnings and showOkCancel('Confirm',msg)):
                
        atomNames   = resDict[seqIds[0]].keys()
 
        molType   = chemComps[0].molType
        ccpCodes0 = []
        startNum  = resDict[seqIds[0]][atomNames[0]]['seqCode']
 
        project = molSystem.root
        molecule = makeMolecule(project, molType, [])

        for i in range(len(seqIds)):
          chemComp = chemComps[i]
 
          if (chemComp.molType != molType) or (i and (seqIds[i] != seqIds[i-1]+1)):
            newMolResidues = addMolResidues(molecule, molType, ccpCodes0, 
                                            startNum=startNum)
            
            # set seqCodes and seqInsertCodes
            xSeqIds = seqIds[i-len(newMolResidues):i]
            for j,x in enumerate(xSeqIds):
              rr = newMolResidues[j]
              for dummy in resDict[x]:
                # future-safe of getting a random atomDict for seqId x
                dd = resDict[x][dummy]
                break
                        
              rr.seqCode = dd['seqCode']
              rr.seqInsertCode = dd['insertCode']
            
            ccpCodes0 = [chemComp.ccpCode,]
            molType   =  chemComp.molType
            startNum  = seqIds[i]

          else:
            ccpCodes0.append(chemComp.ccpCode)
 
        if ccpCodes0:
          addMolResidues(molecule, molType, ccpCodes0, startNum=startNum)
 
        msChain = makeChain(molSystem,molecule,code=sysChainCode)
        
        resMapping = {}
        for i, residue in enumerate(msChain.sortedResidues()):
          resMapping[i] = residue
 
        # TBD deal with HETATMs, proper Molecule name,
        # store CCPN xml as non-standard naming system?

      else:
        continue
        
    else:
      # msChain is not None
      sysChainCode = msChain.code
      resMapping = {}
      for i, residue in mapping:
        resMapping[i] = residue
    
    usedChains.append(msChain)
      
    atomNamesList = []
    msResidues = []
    for j, seqId in enumerate(seqIds):
      atomNames = [tt[0] for tt in resDict[seqId].keys()]
      if fixAtomNames:
        # change names to get right naming system.
        for ii,name in enumerate(atomNames):
          if name[0] in '123':
            atomNames[ii] = name[1:] + name[0]
      atomNamesList.append(atomNames)
      msResidue = resMapping.get(j)
      msResidues.append(msResidue)

    namingSystem = DataConvertLib.getBestNamingSystem(msResidues, atomNamesList)
    
    coordChain = ensemble.newChain(code=sysChainCode)
    #ensemble.override = True

    #
    # Done making chain - start adding residues to chain
    #
 
    for j, seqId in enumerate(seqIds):
      msResidue = msResidues[j]
 
      if not msResidue: # Structure is bigger
        continue

      resName = ccpCodes[j]
      ccpCode = msResidue.ccpCode
 
      if doWarnings:
        if resName != ccpCode:
          msg  = 'Residue names [%s,%s] don\'t match\nin'
          msg += ' loaded molecular system\nchain %s position %d'
          data = (resName,ccpCode,chCode,seqId)
          showWarning('Warning', msg % data)
          continue
 
        if msResidue is None:
          msg  = 'No equivalent molecular system residue'
          msg += '\nfor PDB chain %s residue %d'
          showWarning('Warning', msg % (chCode,seqId))
          continue
 
      coordResidue = coordChain.newResidue(seqCode=msResidue.seqCode,
                                           seqInsertCode=msResidue.seqInsertCode,
                                           seqId=msResidue.seqId)
      coordResidues[(chCode, seqId)] = coordResidue
      allMsResidues[(chCode, seqId)] = msResidue
      #
      # Now make atoms
      #
      
      for atomKey,coordDict in resDict[seqId].items():
        
        atomName, altLoc = atomKey
      
        key = '%s:%s:%s' % (atomName, msResidue, namingSystem)
        
        if key in msAtomDict:
          systemAtom = msAtomDict.get(key)
        
        else:
          systemAtom = DataConvertLib.findMatchingMolSystemAtom(atomName,
                                           msResidue, namingSystem, systemAtoms,
                                           fixAtomNames=fixAtomNames)
          msAtomDict[key] = systemAtom
 
        if (systemAtom is None) or atomCheckDict.get((systemAtom, altLoc)):
          #print '### failing', atomName
          failedAtoms.append('%s %d %s %4s' % 
                             (chCode,seqId,msResidue.ccpCode,atomName) )
          continue
            
        systemAtoms.add(systemAtom)
        #print '### adding', atomName, systemAtom.name
        atomCheckDict[(systemAtom, altLoc)] = True
        
        tt = (coordDict.get('recordId'), (chCode,seqId,atomKey), coordResidue,
              systemAtom.name, coordDict)
        recordDataList.append(tt)
        recordDataDict[(chCode,seqId,atomKey)] = tt
  
  #
  #  Finished getting data for reading order from first model
  #
  #  Now loop over other models and add new atoms to existing residues
  #  ignoring new residues and chains, to cater for varying atom presence
  #
  
  for modelInd, chainsDict in sorted(strucDict.items()[1:]):
    for chCode, resDict in sorted(chainsDict.items()):
      for seqId, atomDict in sorted(resDict.items()):
        coordResidue = coordResidues.get((chCode,seqId))
        msResidue = allMsResidues.get((chCode,seqId))
        if coordResidue is None or msResidue is None:
          continue
        for atomKey,coordDict in atomDict.items():
          tt = recordDataDict.get((chCode,seqId,atomKey))
          if tt is None:
            # new atom, make new record
            atomName, altLoc = atomKey
            key = '%s:%s:%s' % (atomName, msResidue, namingSystem)
        
            if key in msAtomDict:
              systemAtom = msAtomDict.get(key)
 
            else:
              systemAtom = DataConvertLib.findMatchingMolSystemAtom(atomName,
                                              msResidue, namingSystem, systemAtoms,
                                              fixAtomNames=fixAtomNames)
              msAtomDict[key] = systemAtom
 
            if (systemAtom is None) or atomCheckDict.get((systemAtom, altLoc)):
              failedAtoms.append('%s %d %s %4s' % 
                                 (chCode,seqId,msResidue.ccpCode,atomName) )
              continue
            print '### NBNB add extra atom:', (chCode,seqId,atomKey)
 
            systemAtoms.add(systemAtom)
            atomCheckDict[(systemAtom, altLoc)] = True
        
            tt = (coordDict.get('recordId'), (chCode,seqId,atomKey), coordResidue,
                  systemAtom.name, coordDict)
            recordDataList.append(tt)
            recordDataDict[(chCode,seqId,atomKey)] = tt
  
  #
  #  Finished getting data for reading order
  #
  
  # sort by recordId. If not set - or if duplicate recordId
  # will sort by full key.
  recordDataList.sort()
  
  #  
  nAtoms = len(recordDataList)
  
  #
  # Create atoms 
  #
  for tt in recordDataList:
    recordId, fullKey, coordResidue, atomName, coordDict = tt
    
    # make coordAtom
    coordResidue.newAtom(name=atomName, altLocationCode=fullKey[2][1])
    
  
  #
  # set data for all models
  #
  
  kk = 0
  for modelInd, chainsDict in sorted(strucDict.items()):
  
    #create model
    model = ensemble.newModel()
    if modelNames:
      model.name = modelNames[kk]
      kk += 1
  
    # set up for data reading
    coordinates = [0.0] * (3*nAtoms)  # NB deliberate this is 0.0, not None
                                      # There might be atoms missing in other models
    occupancies = [1.0] * nAtoms
    bFactors = [0.0] * nAtoms
    ii = 0
    jj = 0
    hasOccupancies = False
    hasBFactors = False
    
    for tt in recordDataList:
      recordId, fullKey, coordResidue, atomName, coordDict = tt
      
      resDict = chainsDict.get(fullKey[0])
      if resDict:
        atomDict = resDict.get(fullKey[1])
        if atomDict:
          coordDict = atomDict.get(fullKey[2])
          #try:
          #  xxx = atomDict[fullKey[2]]
          #  coordDict = xxx
          #except:
          #  print '###', atomDict.keys(), tt
          
          if coordDict:
            # get data for model
            coordinates[jj] = coordDict['x']
            jj += 1
            coordinates[jj] = coordDict['y']
            jj += 1
            coordinates[jj] = coordDict['z']
            jj += 1
            occupancy = coordDict.get('occupancy', 1.0)
            if occupancy != 1.0:
              occupancies[ii] = occupancy
              hasOccupancies = True
            bFactor = coordDict.get('bFactor', 0.0)
            if bFactor != 0.0:
              bFactors[ii] = bFactor
              hasBFactors = True
            ii += 1
          else:
            # This atom not found in this model
            # Just leave default values, set occupancy to 0
            # and increment indices
            # NBNB this leaves coordinates as 0.0 NBNBNB
            occupancies[ii] = 0.0
            hasOccupancies = True
            ii += 1
            jj += 3
    
  
    # fill data into model
    model.setSubmatrixData('coordinates', coordinates)
    if hasOccupancies:
      model.setSubmatrixData('occupancies', occupancies)
    if hasBFactors:
      model.setSubmatrixData('bFactors', bFactors)

  # final validity check
  try:
    ensemble.checkAllValid()
  except:
    ensemble.delete()
    return
 
  # reset switches
  #ensemble.override = False
  
  if failedAtoms:
    ll = sorted(set(failedAtoms))
    ll = [x for x in ll if not ('Q' in x[-5:] or '*' in x[-5:])]
    #skipAtoms = (' HT1', ' HT2', ' HT3', ' 1HT', ' 2HT', ' 3HT', ' OXT', 
    #            ' OT1', ' OT2', ' 1H', ' 2H', ' 3H')
    #ll = [x for x in failedAtoms 
    #      if 'Q' not in x 
    #      and not [x.endswith(y) for y in skipAtoms]]
    if ll:
      print ('## WARNING %s failed atoms. Unique, non-pesudo atoms are:' 
             % len(failedAtoms), ', '.join(ll))
  
  if failedAtoms and doWarnings:
    msg = 'No equivalent molecular system atoms for PDB atoms: %s'
    showWarning('Warning', msg % ( ' '.join(failedAtoms) ))

  if not ensemble.coordChains:
    ensemble.delete()
    ensemble = None

  return ensemble


def _makeStructureDictFromPdb(fileName, fileType):
  """Make a structure dictionary [model][chainCode][resNum][atomName] =
             (x,y,z coords) from a CNS or true PDB format PDB file.
  .. describe:: Input
  
  String (file name), fileType ('cns' or 'pdb')

  .. describe:: Output

  Dictionary
  """
  
  if fileType == 'cns':
    from ccp.format.cns import coordinatesIO
    coordFile = coordinatesIO.CnsCoordinateFile(fileName)
  elif fileType == 'pdb':
    from ccp.format.pdb import coordinatesIO
    coordFile = coordinatesIO.PdbCoordinateFile(fileName)
  else:
    raise Exception("Illegal fileType: %s" % fileType)
    
  coordFile.read(maxNum = 999)
  
  modelCoords = coordFile.modelCoordinates
  
  dict = {}
  for modelNum in modelCoords:
    # tracking dictionary
    recordId = 0   # used to preserve reading order
    trackDict = {}
    dict[modelNum] = {}
    lastSeqId = 0
    for coord in modelCoords[modelNum]:
      atomName  = coord.atomName
      chainId   = coord.chainId.rstrip()
      seqCode   = coord.seqCode
      insertCode   = coord.insertionCode
      segId     = coord.segId
      altLoc    = coord.altLoc
        
      if chainId == '':
        if segId:
          chainId = segId
        else:
          chainId = 'A'
      
      if dict[modelNum].get(chainId) is None:
        dict[modelNum][chainId] = {}
      
      tt1 = (chainId,seqCode,insertCode)
      seqId = trackDict.get(tt1)
      if seqId is None:
        seqId = trackDict[tt1] = lastSeqId + 1
        lastSeqId = seqId
        dict[modelNum][chainId][seqId] = {}
      
      tt2 = (atomName, altLoc)
      if dict[modelNum][chainId][seqId].get(tt2) is None:
        coordDict = {}
        dict[modelNum][chainId][seqId][tt2] = coordDict
      
      else:
        raise Exception("Duplicate record: %s %s %s %s %s" % 
                        (chainId, seqCode, insertCode, atomName, altLoc))
      
      coordDict['seqCode'] = seqCode
      coordDict['insertCode'] = insertCode
      coordDict['segId'] = segId
      
      coordDict['resName'] = coord.resName
      coordDict['x'] = coord.x
      coordDict['y'] = coord.y
      coordDict['z'] = coord.z
      
      coordDict['recordId'] = recordId
      recordId += 1
      
      if fileType == 'pdb':
        coordDict['occupancy'] = coord.occupancy
        coordDict['bFactor']   = coord.bFactor
        coordDict['atomType']  = coord.atomType
        coordDict['hetFlag']   = coord.hetFlag
       
  return dict


def _makeStructureDictFromRoughPdb(fileName):
  """Make a structure dictionary [model][chainCode][resNum][atomName] =
             (x,y,z coords) from a non-standard PDB file.
  .. describe:: Input
  
  String (file name)

  .. describe:: Output

  Dictionary
  """
  
  # Rasmus Fogh 
  # Added recordId to preserve atom writing order
  #
  # Reorganised to preserve writing order and to work with new data model
  
  residueKeys = set([])
  insertOffset = 0

  result = {}
  modelNum = 0
  fileHandle = open(fileName)
  recordId = 0
  for line in fileHandle.readlines():
    
    key = line[0:6].strip()
    if key == 'ENDMDL':
      modelNum +=1
      recordId = 0
      insertOffset = 0
      residueKeys = set([])
      
    elif key == 'ATOM':
      #serial    = returnInt(line[6:11])
      atomName  = line[12:16].strip()
      altLoc    = line[16:17]
      resName   = line[17:20].strip()
      chainId   = line[21:22].strip()
      seqCode   = returnInt(line[22:26])
      insertCode = line[26:27]
      x         = returnFloat(line[30:38])
      y         = returnFloat(line[38:46])
      z         = returnFloat(line[46:54])
      occupancy = returnFloat(line[54:60])
      bFactor   = returnFloat(line[60:66])
      segId     = line[72:76].strip()
      atomType  = line[76:78].strip()
      seqId = seqCode + insertOffset
      #charge    = line[78:80].strip()
 
      if chainId == '':
        if segId:
          chainId = segId
        else:
          chainId = 'A'
 
      if result.get(modelNum) is None:
        result[modelNum] = {}
      if result[modelNum].get(chainId) is None:
        result[modelNum][chainId] = {}
        
      resKey = (chainId, seqCode, insertCode)
      if resKey not in residueKeys:
        # Not seen this residue key before
        residueKeys.add(resKey)
        
        if result[modelNum][chainId].get(seqId) is not None:
          # SeqId is known -> seq insert
          # Residue number is one higher
          insertOffset += 1
          seqId += 1
        
        result[modelNum][chainId][seqId] = {}
      
      atomKey = (atomName, altLoc)
      if result[modelNum][chainId][seqId].get(atomKey) is None:
        coordDict = {}
        result[modelNum][chainId][seqId][atomKey] = coordDict
      
      else:
        raise Exception("Duplicate record: %s %s %s %s %s" % 
                        (chainId, seqCode, insertCode, atomName, altLoc))
      
      coordDict['resName'] = resName
      coordDict['seqCode'] = seqCode
      coordDict['insertCode'] = insertCode
      coordDict['segId'] = segId
      coordDict['x'] = x
      coordDict['y'] = y
      coordDict['z'] = z
      coordDict['occupancy'] = occupancy
      coordDict['bFactor']   = bFactor
      coordDict['atomType']  = atomType
 
      coordDict['recordId']  = recordId
      recordId += 1
      
 
  return result
  

def makePdbFromStructure(fileName,structure,model=None,useOxt=False):
  """Make a PDB file from a structure or ensemble of structures.
             If a model is passed in the PDB style file will contain only
             that model's coordinates, otherwise all models will be used. 
  .. describe:: Input
  
  String (file name), MolStructure.StructureEnsemble,
             MolStructure.Model

  .. describe:: Output

  None
  """
  """
  COLUMNS        DATA  TYPE    FIELD        DEFINITION
  -------------------------------------------------------------------------------------
   1 -  6        Record name   "ATOM  "
   7 - 11        Integer       serial       Atom  serial number.
  13 - 16        Atom          name         Atom name.
  17             Character     altLoc       Alternate location indicator.
  18 - 20        Residue name  resName      Residue name.
  22             Character     chainID      Chain identifier.
  23 - 26        Integer       resSeq       Residue sequence number.
  27             AChar         iCode        Code for insertion of residues.
  31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
  39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
  47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
  55 - 60        Real(6.2)     occupancy    Occupancy.
  61 - 66        Real(6.2)     tempFactor   Temperature  factor.
  77 - 78        LString(2)    element      Element symbol, right-justified.
  79 - 80        LString(2)    charge       Charge  on the atom.
  """

  if model:
    models = [model,]
  else:
    models = structure.sortedModels()
 
  fileHandle = open(fileName, 'w')
  lFormat = '%-80.80s\n'
  pdbFormat = '%-6.6s%5.1d %4.4s%s%3.3s %s%4.1d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  \n'
  terFormat = '%-6.6s%5.1d      %s %s%4.1d%s                                                     \n'                            
  
  fileHandle.write('REMARK 210 REMARK:\n')
  fileHandle.write('REMARK 210\n')
  fileHandle.write('REMARK 210 Generated by CcpNmr Analysis\n')
  for model in models:
    remark = ('REMARK 210 CCPN ccp.molecule.MolStructure.Model %s' 
              % model.getExpandedKey())
    fileHandle.write(lFormat  % remark)
  
  orderedAtoms = structure.orderedAtoms
  
  # map atom record name and chain code, to save time
  atomFields = {}
  chainCodes = {}
  for chain in structure.coordChains:
    if chain.chain.molecule.molType == 'nonpolymer':
      atomFields[chain] = 'HETATM'
    else:
      atomFields[chain] = 'ATOM'
    
    ch = chain.chain.pdbOneLetterCode
    if not ch.strip():
      ch = chain.chain.code[0]
    chainCodes[chain] = ch
  
  for m, model in enumerate(models):
    
    if len(models) > 1:
      line = 'MODEL     %4d' % (m+1)
      fileHandle.write(lFormat  % line)
    
    coordinates = model.coordinates
    occupancies = model.occupancies
    bFactors = model.bFactors
    
    lastChainCode = None
    offset = 0
    lineno = 0
    for ii,atom in enumerate(orderedAtoms):
      
      residue = atom.residue
      chain = residue.chain
      
      a = atomFields[chain]
      ch = chainCodes[chain]
      
      if lastChainCode is not None and lastChainCode != ch:
        # Change to a new chain. Write TER record for previous chain
        lineno += 1
        line  = terFormat % ('TER',lineno,tlc,lastChainCode,s,ins)
        fileHandle.write(line)
      lastChainCode = ch
      
      msResidue = residue.residue
      chemComp = msResidue.chemCompVar.chemComp
      tlc = chemComp.code3Letter
      s = msResidue.seqCode
      ins = msResidue.seqInsertCode
    
      n     = '%-3s' % atom.name
      if useOxt and (n=="O''"):
        n = 'OXT'
    
      alc   = atom.altLocationCode
      x = coordinates[offset]
      offset += 1
      y = coordinates[offset]
      offset += 1
      z = coordinates[offset]
      offset += 1
      o = occupancies[ii]
      b = bFactors[ii]
      e = atom.elementSymbol
      
      lineno += 1
      line  = pdbFormat % (a,lineno,n,alc,tlc,ch,s,ins,x,y,z,o,b,e)
      fileHandle.write(line)
    
    if len(models) > 1:
      fileHandle.write(lFormat  % 'ENDMDL')
  
  fileHandle.write(lFormat  % 'END')
        
  fileHandle.close()
  
