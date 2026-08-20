import os,string,sys

#
# CCPN specific...
#

from memops.api import Implementation
from memops.universal.Util import makePowerSet, drawBox

from ccp.general.Util import getOtherAtom
from ccp.general.Io import getChemComp

from ccpnmr.format.converters.Mol2Format import Mol2Format
 
from pdbe.chemComp.Constants import origMol2DataDir
from pdbe.chemComp.Io import findChemCompOrCoordFilePath
from pdbe.chemComp.Io import saveTemporaryChemCompOrCoord, consolidateTemporaryChemCompOrCoord
from pdbe.chemComp.Io import findExistingChemCompInfo, findExistingChemCompCoordInfo
from pdbe.chemComp.Constants import editChemCompCoordDataDir, testChemCompCoordDataDir
from pdbe.chemComp.Constants import editChemCompDataDir, testChemCompDataDir

def getStereoPriorities(stereoAtom,chemBondsHandled):

  priorityList = []
  priorityKeys = []
  followPriorityKeys = []
 
  for chemBond in stereoAtom.chemBonds:
    if chemBond in chemBondsHandled:
      otherAtom = None
      priorityKey = (0,0)
    else:
      otherAtom = getOtherAtom(stereoAtom,chemBond)
      priorityKey = (otherAtom.chemElement.atomNumber,bondDict[chemBond.bondType])
      chemBondsHandled.append(chemBond)
      
      if priorityKey not in priorityKeys:
        priorityKeys.append(priorityKey)
      elif priorityKey not in followPriorityKeys:
        followPriorityKeys.append(priorityKey)
    
    priorityList.append((priorityKey,otherAtom))
 
  return (priorityList,followPriorityKeys)

def findRingCarbon(chemAtom,connectedAtoms):
  
  returnConnectedAtoms = None
  
  for chemBond in chemAtom.chemBonds:
  
    otherChemAtom = getOtherAtom(chemAtom,chemBond)

    if otherChemAtom == connectedAtoms[-1]:
      returnConnectedAtoms = connectedAtoms[:]
      break
          
    elif otherChemAtom and otherChemAtom not in connectedAtoms and otherChemAtom.elementSymbol == 'C':
      
      returnConnectedAtoms = findRingCarbon(otherChemAtom,connectedAtoms[:-1] + [otherChemAtom,connectedAtoms[-1]])

      if returnConnectedAtoms:
        break
        
  return returnConnectedAtoms

def getChemAtomKey(chemAtom):

  return (chemAtom.name,chemAtom.subType)
            
def getLinkAtom(chemAtom,bondType = 'single', multi = False):

  if not multi:
    linkAtomName = "%s_1" % chemAtom.name
  else:
    linkAtomName = "%s_%d_1" % (chemAtom.name,chemAtom.subType)
  
  linkAtom = chemAtom.parent.findFirstChemAtom(className = 'LinkAtom', name = linkAtomName)
  if not linkAtom:
    linkAtom = chemAtom.parent.newLinkAtom(name = linkAtomName)
    chemAtom.parent.newChemBond(chemAtoms = [chemAtom,linkAtom], bondType = bondType)
  
  return linkAtom

  
def makeFullSugar(carboBaseName,coordSystem,baseGlycoCtCode,testMode, replace = False, saveData = True):

  carboMolType = 'carbohydrate'
  namingSystemName = 'EuroCarbDb'

  project = Implementation.MemopsRoot(name = 'chemComp')
  project.currentUserId = 'ccpnRef'
  
  #
  # Set archive dir info
  #
  
  if testMode:
    chemCompDataDir = testChemCompDataDir
    chemCompCoordDataDir= testChemCompCoordDataDir
  else:
    chemCompDataDir = editChemCompDataDir
    chemCompCoordDataDir= editChemCompCoordDataDir
  
  #
  # First import all mol2 files, pick one form as 'base' unit, adapt this one,
  # then add coords from other chemComps
  #
  # TODO: This is currently very specific!
  #
  
  #
  # TODO SET A STEREOCHEMISTRY CLASS FOR THE A/B CHEMATOMS!!
  #
 
  importDir = os.path.join(origMol2DataDir,'carbo')
  importFiles = os.listdir(os.path.join(importDir,carboBaseName))
  
  for importFile in importFiles[:]:
    if not importFile[-4:] == 'mol2':
      importFiles.pop(importFiles.index(importFile))
  
  importFiles.sort()
  
  
  mol2Format = Mol2Format(project, guiParent = None, allowPopups = False)
  molTypes = [carboMolType]
  
  chemComps = []
  rawChemComps = []
  
  for importFile in importFiles:
    
    # Should be 'a', 'b' or 'o'
    anomericCenter = importFile[0] 
    
    if anomericCenter != 'a':
      ccpCode = '%s-%s' % (anomericCenter,carboBaseName)
    else:
      ccpCode = carboBaseName

    fileName = os.path.join(importDir,carboBaseName,importFile)
    
    print "Reading mol2 file %s..." % fileName

    ccs = mol2Format.readChemComps(
                             fileName,
                             molTypes = molTypes,
                             ccpCodes = [ccpCode],
                             saveChemComp = False,
                             minimalPrompts = True,
                             makeNamingSystem = namingSystemName)
  
    chemComps.append(ccs[0])
    rawChemComps.append(mol2Format.rawChemComp)
  
  #
  # Check whether only open form available (e.g. aldi ones)
  # 
  
  if len(chemComps) == 1 and chemComps[0].ccpCode[0] == 'o':
    print "  Warning: only open form available, not creating a/b isoforms."
    hasOnlyOpenForm = True
  else: 
    hasOnlyOpenForm = False
  
  #print chemComps
  #print rawChemComps
    
  """
  for cch in project.chemCompHeads:
    print cch.molType, cch.ccpCode
    for ccv in cch.chemComp.chemCompVars:
      print ccv.descriptor
      print ccv.chemAtoms
    print
  """
  
  refChemComp = chemComps[0] # Should be the a form  
  
  #
  # Reset save location, check if file already exists
  #

  chemCompXmlFile = findChemCompOrCoordFilePath(refChemComp,testMode = testMode)
 
  # In this case, getting nothing back with replace - False means that it does exist!
  if chemCompXmlFile and not replace:
    print "  ChemComp %s, %s already exists - aborting creation." % (carboMolType,carboBaseName) 

    try:
      refChemComp = getChemComp(project, carboMolType, carboBaseName, download=False, chemCompArchiveDir = editChemCompDataDir,copyFile=False)
      
    except:
      print "WARNING: chemcomp was already loaded!"
      refChemComp = project.findFirstChemComp(molType = carboMolType, ccpCode = carboBaseName)

    return refChemComp
          
  #
  # Start creating/modifying...
  #
  
  print
  print drawBox("Creating sugar information")
  print
          
  #
  # Set the base Glyco CT code, this is always x-, except for -o only forms (aldehydes)
  #  
  # TODO: might need to hack this for substituents so know which one is which... or just do it by order? Should be fine...
  #
  
  refGlycoCtCode = "RES\n1b:%s" % baseGlycoCtCode
  
  print "Setting GlycoCT code to:\n\n%s\n" % refGlycoCtCode
  print
  
  project.override = True
  try:
    refChemComp.baseGlycoCtCode = refGlycoCtCode
    refChemComp.findFirstChemCompVar().glycoCtCode = refGlycoCtCode
  finally:
    project.override = False

  #
  # Look for a C-O-C fragment. Should be O5 for 6 rings.
  #
  # Note: ONLY works on cyclic sugars!!!
  #

  centralOAtom = None

  for chemAtom in refChemComp.sortedChemAtoms():

    # TODO: should really do this from the Var level, but...
    if chemAtom.elementSymbol == 'O' and len(chemAtom.chemBonds) == 2:
      connectedToC = True
      for chemBond in chemAtom.sortedChemBonds():
        otherChemAtom = getOtherAtom(chemAtom,chemBond)
        if otherChemAtom.elementSymbol != 'C':
          connectedToC = False
          break

      if connectedToC: 

        centralOAtom = chemAtom

        #
        # Try to generically determine the carbon atoms in the ring... trying not to depend on names
        # but if a C1 is connected to the O5, then start from there.
        #

        connectedAtoms = []

        for chemBond in centralOAtom.sortedChemBonds():
          otherChemAtom = getOtherAtom(centralOAtom,chemBond)

          if otherChemAtom.name == 'C1':
            connectedAtoms.insert(0,otherChemAtom)
          else:
            connectedAtoms.append(otherChemAtom)

        #
        # Start to loop... use recursive function to find right order
        #

        ringCarbons = findRingCarbon(connectedAtoms[0],connectedAtoms)

        if not ringCarbons:
          centralOAtom = None
          continue

        otherCarbons = []
        for searchCarbon in refChemComp.findAllChemAtoms(elementSymbol = 'C'):
          if searchCarbon not in ringCarbons:
            otherCarbons.append(searchCarbon)

        break

  if not centralOAtom and not hasOnlyOpenForm:
    raise "  Error: no central O atom found in ring!!"

  #
  # Now look for anomeric carbon and connection sites...
  #

  anomericCarbon = None
  anomericOxygen = None

  if not hasOnlyOpenForm:
    for searchCarbon in [ringCarbons[0],ringCarbons[-1]]:
      for chemBond in searchCarbon.chemBonds:
        otherChemAtom = getOtherAtom(searchCarbon,chemBond)
        if otherChemAtom and otherChemAtom.elementSymbol == 'O' and otherChemAtom != centralOAtom:
          anomericCarbon = searchCarbon 
          anomericOxygen = otherChemAtom
          break

      if anomericCarbon:
        break
      
  else:
    # Hardset these... not really anomeric but good enough
    anomericCarbon = refChemComp.findFirstChemAtom(name = 'C1')
    anomericOxygen = refChemComp.findFirstChemAtom(name = 'O1')
      

  if not anomericCarbon:
    raise "  Error: no anomeric carbon found."

  else:
    for chemBond in anomericOxygen.chemBonds:
      otherChemAtom = getOtherAtom(anomericOxygen,chemBond)
      if otherChemAtom and otherChemAtom.elementSymbol == 'H':
        anomericHydrogen = otherChemAtom 
        break

    if not anomericHydrogen:
      print "  Warning: no anomeric hydrogen found."

  
  #
  # Set the stereo information for the anomeric carbon, then create subtypes for beta and open forms
  #
  
  #anomericCarbons = {'a': None, 'b': None, 'o': None}
  
  
  stereoAtom = anomericCarbon
  bondDict = {'single': 1,  'double': 2, 'triple': 3, 'aromatic': 1.5, 'dative': 1.0, 'singleplanar': 1.5}  # TODO CHECK THIS!
  
  """
  # OK to do this because is neutral 'real' chemComp, but otherwise need to start from chemCompVar!!!
  
  #
  # TODO TODO also need to track chemBonds, but ONLY at beginning
  # then need to go back to coordinates to find out in which order the bonds appear based on the coords (when looking down
  # the main bond)!
  #
  # Might be easiest to use Martin's code (or someone else's), use the atom names to link back to here... 
  #
  totalChemBonds = chemComp.chemBonds
  chemBondsHandled = []
  
  priorityList = getStereoPriorities(stereoAtom,chemBondsHandled)
  
  priorityList.sort()
  priorityList.reverse()
  
  priorityKeys = []  
  for (priorityKey,otherAtom) in priorityList:
   priorityKeys.append(priorityKey)
  
  for priorityKey in priorityKeys:
    if priorityKeys.count(priorityKey) > 1:
      print priorityKey
  
  print priorityList

  sys.exit()
  """
  #
  # Look for binding oxygens... has to be hydroxy group connected to carbon
  #

  bindingOxygens = []
  bindingHydrogens = {}
  
  #
  # Set carbons to search for connected OH groups
  #
  
  if not hasOnlyOpenForm:
    searchCarbons = ringCarbons + otherCarbons
  else:
    searchCarbons = list(refChemComp.findAllChemAtoms(elementSymbol = 'C'))
    
    if 'aldi' in baseGlycoCtCode:
      # Don't do anything on 1 position for these - alditols
      searchCarbons.pop(searchCarbons.index(anomericCarbon))

  for searchCarbon in searchCarbons:

    validConnectedAtoms = {}
    validOHgroups = []

    for chemBond in searchCarbon.sortedChemBonds():

      otherChemAtom = getOtherAtom(searchCarbon,chemBond)

      if otherChemAtom:

        for elementSymbol in ['O','N']:

          if otherChemAtom.elementSymbol == elementSymbol and otherChemAtom != centralOAtom:

            if not validConnectedAtoms.has_key(elementSymbol):
              validConnectedAtoms[elementSymbol] = []
            validConnectedAtoms[elementSymbol].append(otherChemAtom)

            otherChemBonds = list(otherChemAtom.chemBonds)

            if elementSymbol == 'O' and len(otherChemBonds) == 2:

              otherChemBond = otherChemBonds[not otherChemBonds.index(chemBond)]

              connectedChemAtom = getOtherAtom(otherChemAtom,otherChemBond)

              if connectedChemAtom.elementSymbol == 'H':

                validOHgroups.append(otherChemAtom)


    #
    # Check if single OH (no double O, or amide, ...)
    #

    if len(validOHgroups) == 1:

      if len(validConnectedAtoms) == 1:
        # Carboxylic acid (except for ring O-C-OH!)
        if len(validConnectedAtoms['O']) == 2 and not searchCarbon == anomericCarbon:
          print "  Warning: ignoring oxygen %s - is carboxylic acid (or similar)" % validOHgroups[0].name
          validOHgroups = []
      else:
        # Amide or something similar
        print "  Warning: ignoring oxygen %s - is amide (or similar)" % validOHgroups[0].name
        validOHgroups = []

      if validOHgroups:

        bindingOxygens.append(validOHgroups[0])
        bindingHydrogens[validOHgroups[0]] = connectedChemAtom

        if not hasOnlyOpenForm and searchCarbon in otherCarbons:
          print "  Warning: setting oxygen %s as binding one (not directly connected to ring)." % (validOHgroups[0].name)
  
  
  #
  # Now create variants...
  #
  # Need to have all combinations of:
  #
  #  - anomeric a/b and free/bound
  #  - binding oxygens free/bound
  #

  bindingOxygenCombs = makePowerSet(bindingOxygens)
  origAnomericCarbon = anomericCarbon
  
  linkAtomsMapForCoordinates = {}
  
  if hasOnlyOpenForm:
    stereoTypes = ('open_1',)
  else:
    stereoTypes = ('stereo_1','stereo_2')
  

  for stereoType in stereoTypes:

    subType = int(stereoType.split('_')[1])
    anomericCarbon = refChemComp.findFirstChemAtom(name = 'C1', subType = subType)
    
    if not anomericCarbon:
      # Should only happen for stereo_2
      
      creationDict = {'name': 'C1', 'subType': subType}
      for attrName in ('elementSymbol','shortVegaType','waterExchangeable'):
        creationDict[attrName] = getattr(origAnomericCarbon,attrName)
      
      # TODO set chirality to OPPOSITE of whatever the first subtype is!
      anomericCarbon = refChemComp.newChemAtom(**creationDict)
      
      namingSystem = refChemComp.findFirstNamingSystem(name = namingSystemName)
      namingSystem.newAtomSysName(sysName = anomericCarbon.name, atomName = anomericCarbon.name, atomSubType = anomericCarbon.subType)

      # Bonds and other atoms are exactly the same!
      for chemBond in origAnomericCarbon.chemBonds:
        otherChemAtom = getOtherAtom(origAnomericCarbon,chemBond)
        refChemComp.newChemBond(chemAtoms = (anomericCarbon,otherChemAtom),bondType = chemBond.bondType, stereochem = chemBond.stereochem)
  
    for anomericBound in range(0,2):
      for i in range(0,len(bindingOxygenCombs)):
        bindingOxygens = bindingOxygenCombs[i]

        # Main list only contains real ChemAtoms!!
        currentChemAtoms = list(refChemComp.findAllChemAtoms(className = 'ChemAtom'))
        anomericCarbons = refChemComp.findAllChemAtoms(name = 'C1')
        for tempAnomericCarbon in anomericCarbons:
          if tempAnomericCarbon != anomericCarbon and tempAnomericCarbon in currentChemAtoms:
            currentChemAtoms.pop(currentChemAtoms.index(tempAnomericCarbon))

        # This is the neutral refChemComp - already exists!
        if not anomericBound and not bindingOxygens:
          continue

        # Can't be both bound on the anomeric carbon and have an oxygen link there...
        if anomericBound and len(bindingOxygens) == 1 and bindingOxygens[0] == anomericOxygen:
          continue

        # Can't have a link to 1 and something else...
        if len(bindingOxygens) > 1 and anomericOxygen in bindingOxygens:
          continue

        linkedAtomKeys = []
        linkedAtoms = {}
        linkAtoms = {}

        if anomericBound and not hasOnlyOpenForm and 'aldi' not in baseGlycoCtCode:
          # Only do C1 when relevant - for alditols it's not, always reducing end.
          currentChemAtoms.pop(currentChemAtoms.index(anomericOxygen))
          if anomericHydrogen:
            currentChemAtoms.pop(currentChemAtoms.index(anomericHydrogen))
            
          anomericCarbonKey = getChemAtomKey(anomericCarbon)
          linkedAtomKeys.append(anomericCarbonKey)

          linkedAtoms[anomericCarbonKey] = anomericCarbon

          linkAtom = getLinkAtom(anomericCarbon,multi = True)

          linkAtoms[anomericCarbonKey] = linkAtom
          currentChemAtoms.append(linkAtom)
          
          linkAtomsMapForCoordinates[(linkAtom.name,linkAtom.subType)] = anomericOxygen.name

        for bindingOxygen in bindingOxygens:

          bindingHydrogen = bindingHydrogens[bindingOxygen]

          if bindingHydrogen and bindingHydrogen in currentChemAtoms:
            currentChemAtoms.pop(currentChemAtoms.index(bindingHydrogen))
          
          bindingOxygenKey = getChemAtomKey(bindingOxygen)
          linkedAtomKeys.append(bindingOxygenKey)
          linkedAtoms[bindingOxygenKey] = bindingOxygen

          linkAtom = getLinkAtom(bindingOxygen)

          linkAtoms[bindingOxygenKey] = linkAtom
          currentChemAtoms.append(linkAtom)
          
          linkAtomsMapForCoordinates[(linkAtom.name,linkAtom.subType)] = bindingHydrogen.name

        linkedAtomKeys.sort()
        
        # Possible for C1 linkages if not handled (only reducing end)
        if not linkedAtomKeys:
          continue

        #
        # Create linkEnds
        #
        
        linkInfo = []

        for linkedAtomKey in linkedAtomKeys:
        
          if linkedAtomKey[0] == 'C1':
            linkCode = "%s_%s" % (linkedAtomKey[0],linkedAtomKey[1])
          else:
            linkCode = linkedAtomKey[0]
            
          linkInfo.append(linkCode)
        
          if not refChemComp.findFirstLinkEnd(linkCode = linkCode):
            boundChemAtom = linkedAtoms[linkedAtomKey]
            boundLinkAtom = linkAtoms[linkedAtomKey]
            linkEnd = refChemComp.newLinkEnd(linkCode = linkCode, boundChemAtom = boundChemAtom, boundLinkAtom = boundLinkAtom)

        # TODO THIS IS NOT GREAT - no way of telling what's what if atom names are messed up
        #linking = 'none'
        #descriptor = 'link:%s' % string.join(linkedAtomKeys,',')
        linking = 'link:%s' % string.join(linkInfo,',')
        
        if stereoType.count('stereo'):
          descriptor = '%s:C1' % stereoType
        else:
          descriptor = 'neutral'

        if not refChemComp.findFirstChemCompVar(linking = linking, descriptor = descriptor):
          print "  Trying %s,%s" % (linking,descriptor)
          #for ca in currentChemAtoms:
          #  if ca.className == 'LinkAtom':
          #    print "   LA:",ca.name, ca.subType
          #  else:
          #    print "   CA:",ca.name, ca.subType
          #print "   ",linkedAtomKeys
          
          # Create the stereospecific GlycoCt code
          if stereoType == "stereo_1":
            stereoCode = 'a'
          elif stereoType == "stereo_2":
            stereoCode = 'b'
          elif stereoType == 'open_1':
            stereoCode = 'o'
          
          varGlycoCtCode = "RES\n1b:%s" % (stereoCode + baseGlycoCtCode[1:])
          
          ccv = refChemComp.newChemCompVar(chemAtoms = currentChemAtoms, linking=linking, descriptor= descriptor, glycoCtCode = varGlycoCtCode, formalCharge=0, isParamagnetic=False, isAromatic=False)

  #
  # Make sure the chemElements are accessible
  #
  
  project.currentChemElementStore = project.findFirstChemElementStore()

  #
  # Reset the name and molType...
  #

  chMolType = 'carbohydrate'

  project.override = True
  try:
    refChemComp.ccpCode = carboBaseName
    refChemComp.molType = carboMolType
  finally:
    project.override = False

  # TODO SET THESE CORRECTLY? Where do I get info from for this though? Can this come from MSD? Ask Dimitris!!
  # Set the PDB/MSD name - NOTE that have to do this for the correct a/b forms! Var specific!!
  #ChemComp.ChemCompSysName(refChemComp,namingSystem = 'PDB',sysName=ccpCode,specificChemCompVars = refChemComp.sortedChemCompVars())
  #ChemComp.ChemCompSysName(refChemComp,namingSystem = 'MSD',sysName=ccpCode,specificChemCompVars = refChemComp.sortedChemCompVars())
  
  #
  # Check chemComp validatity and save
  # TODO make this all options in running script!
  #
  
  refChemComp.checkAllValid()
  
  if saveData:

    #
    # Get the original file GUID, if possible, when replacing existing file
    #
    
    if replace:

      (existingGuid,existingFile) = findExistingChemCompInfo(chemCompDataDir,refChemComp.ccpCode,refChemComp.molType)

      if existingGuid:
        project.override = True
        try:
          refChemComp.guid = existingGuid
        finally:
          project.override = False

    (tmpFilePath,existingFilePath) = saveTemporaryChemCompOrCoord(refChemComp,testMode = testMode)

    # Do a check here? Or don't bother?
    consolidateTemporaryChemCompOrCoord(refChemComp,tmpFilePath,existingFilePath,testMode=testMode,replace=replace)    

  #
  # Get the coordinates as well!
  #
      
  print "  Creating coordinates!!"
  
  chemCompCoord = project.newChemCompCoord(sourceName = coordSystem,molType = refChemComp.molType, ccpCode = refChemComp.ccpCode)
  
  for i in range(len(rawChemComps)):
  
    rawChemComp = rawChemComps[i]
    
    # Identify which one we're dealing with!!
    (dirName,baseName) = os.path.split(rawChemComp.parent.name)
    
    if baseName[0] == 'a':
      stereoDescriptor = "stereo_1:C1"
      
    elif baseName[0] == 'b':
      stereoDescriptor = "stereo_2:C1"
      
    elif baseName[0] == 'o':
      stereoDescriptor = "none"
    
    else:
      print "  Not handling type '%s' for coordinates - ignored." % baseName[0]
      continue
          
    #
    # Mark that generated by this script...
    #

    applData = Implementation.AppDataString(application = 'ccpNmr', keyword = 'origin', value = 'makeFullSugar.py')
    chemCompCoord.addApplicationData(applData)
    
    # Don't do any link atoms (yet)... could in principle use atoms that are 'missing'
    
    # TODO: should decompose descriptor here, then check...
    chemCompVars = refChemComp.findAllChemCompVars(descriptor = stereoDescriptor)
    
    chemAtomKeys = []
    
    for ccv in chemCompVars:
      for ca in ccv.sortedChemAtoms():
        caKey = (ca.name,ca.subType)
        if caKey not in chemAtomKeys:
          chemAtomKeys.append(caKey)

    #
    # Create a dictionary for the coordinates, based on the 'raw' chemComp from the mol2 file
    #
          
    chemAtomCoordDict = {}
    
    for chemAtomKey in chemAtomKeys:

      coords = None
      
      if linkAtomsMapForCoordinates.has_key(chemAtomKey):
        useChemAtomName = linkAtomsMapForCoordinates[chemAtomKey]
      else:
        useChemAtomName = chemAtomKey[0]

      for rawAtom in rawChemComp.atoms:
      
        if rawAtom.name == useChemAtomName:
          coords = (rawAtom.x, rawAtom.y, rawAtom.z)
          break
          
      if not coords:
        print "  Warning: no coordinate for %s, atom key %s." % (coordSystem,chemAtomKey)
      elif not chemAtomCoordDict.has_key(chemAtomKey):
        chemAtomCoordDict[chemAtomKey] = coords
      else:
        print "  Error: double atom key %s!" % chemAtomKey

    #
    # Set the coordinates
    #
    
    chemAtomCoords = {}

    #print [ra.name for ra in rawChemComp.atoms]

    for chemCompVar in chemCompVars:
    
      #print chemCompVar

      chemCompVarCoord = chemCompCoord.findFirstChemCompVarCoord(linking = chemCompVar.linking, descriptor = chemCompVar.descriptor)

      if not chemCompVarCoord:
        chemCompVarCoord = chemCompCoord.newChemCompVarCoord(linking = chemCompVar.linking, descriptor = chemCompVar.descriptor)
        
      #print chemCompVarCoord

      for ca in chemCompVar.sortedChemAtoms():
      
        caKey = (ca.name,ca.subType)
        
        #print "%-20s" % str(caKey),
        
        if chemAtomCoordDict.has_key(caKey):
          coords = chemAtomCoordDict[caKey]

          if coords:
          
            if chemAtomCoords.has_key(caKey):
              chemAtomCoord = chemAtomCoords[caKey]
            else:
              chemAtomCoord = chemCompCoord.newChemAtomCoord(name = caKey[0], subType = caKey[1], x = coords[0], y = coords[1], z = coords[2])
              chemAtomCoords[caKey] = chemAtomCoord
          
            if chemAtomCoord not in chemCompVarCoord.chemAtomCoords:
              chemCompVarCoord.addChemAtomCoord(chemAtomCoord)
              #print chemAtomCoord.name, chemAtomCoord.subType,
        
        #print

  chemCompCoord.checkAllValid()
  
  if saveData:

    if replace:

      (existingGuid,existingFile) = findExistingChemCompCoordInfo(chemCompCoordDataDir,coordSystem,chemCompCoord.ccpCode,chemCompCoord.molType)

      if existingGuid:
        project.override = True
        try:
          chemCompCoord.guid = existingGuid
        finally:
          project.override = False

    (tmpFilePath,existingFilePath) = saveTemporaryChemCompOrCoord(chemCompCoord,testMode = testMode)
    
    # Do a check here? Or don't bother? This is blank regeneration from reference data, so should be OK!
    consolidateTemporaryChemCompOrCoord(chemCompCoord,tmpFilePath,existingFilePath,testMode=testMode,replace=replace)    
  
  return refChemComp

if __name__ == '__main__':

  from pdbe.chemComp.Constants import editChemCompDataDir, testChemCompDataDir

  carboBaseName = sys.argv[1]
  baseGlycoCtCode = sys.argv[2]
  coordSystem = 'euroCarbDb'
   
  if '-create' in sys.argv:
    print "Warning: creating new sugar in edit/ directory!"
    testMode = False
  else:
    print "Creating in test directory!"
    testMode = True
    
  makeFullSugar(carboBaseName,coordSystem,baseGlycoCtCode,testMode,saveData = False)

