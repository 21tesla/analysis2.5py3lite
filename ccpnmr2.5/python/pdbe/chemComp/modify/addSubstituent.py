import os, re

from memops.api import Implementation
from memops.general.Io import getCcpFileString
from memops.universal.Util import drawBox
from memops.universal.Geometry import superposeNewVectorsOnOld

from ccp.general.Util import getOtherAtom, getDescriptorDict, getDescriptorFromDict, mergeDescriptorDicts
from ccp.general.Io import getChemComp, getChemCompCoord

from ccpnmr.format.converters.Mol2Format import Mol2Format

from pdbe.chemComp.Constants import origMol2DataDir
from pdbe.chemComp.Io import findChemCompOrCoordFilePath, saveTemporaryChemCompOrCoord
from pdbe.chemComp.Io import consolidateTemporaryChemCompOrCoord
from pdbe.chemComp.Io import findExistingChemCompInfo, findExistingChemCompCoordInfo

from pdbe.chemComp.modify.Constants import substituentInfo

from pdbe.chemComp.Constants import editChemCompDataDir, testChemCompDataDir
from pdbe.chemComp.Constants import editChemCompCoordDataDir, testChemCompCoordDataDir


copyDict = {

  'ChemAtom':       {'attrs': ['subType','elementSymbol','waterExchangeable','nuclGroupType']},
      # TODO: ATTRIBUTES THAT WILL HAVE TO BE RESET BY VEGA (or antechamber)... SETUP:
      # 'shortVegaType', 'chirality',

  'LinkAtom':       {'attrs': ['subType']},

  'ChemAtomSet':     {'attrs': ['subType','isEquivalent','isProchiral','distCorr']},
  
  'ChemBond':        {'attrs': ['bondType','stereochem'], 'links': ['chemAtoms']},
  
  'ChemAngle':       {'attrs': [], 'links': ['chemAtoms']},
  
  'ChemTorsion':     {'attrs': ['name'], 'links': ['chemAtoms']},
  
  'LinkEnd':         {'attrs': ['linkCode'], 'links': ['boundChemAtom','boundLinkAtom','remoteChemAtom','remoteLinkAtom']},
  # TODO here have to make sure that, if working off the CCPN directly, these SUB linkends are ignored!!
  
  'Stereochemistry': {'attrs': ['value'], 'links': ['refStereochemistry']},
  # TODO this will not work - need to find matching class beforehand

  
    }

chemCompLinkList = ['chemBonds','chemAngles','chemTorsions','linkEnds','stereochemistries']
# TODO TODO 'chemAtomSysNames': could be done, in principle, with a specific mechanism based on molType, ...

# Not treating: 'sysNames': is irrevelevant because name will change for merged chemComp.

class CreateCcpnObject:

  def __init__(self,ccpnObject):
  
    self.className = ccpnObject.className
    self.ccpnObject = ccpnObject
    
    self.creationKeywds = {}
    for attrName in copyDict[self.className]['attrs']:
      self.creationKeywds[attrName] = getattr(ccpnObject,attrName)
      
  def setLinks(self,substToBaseDict):
          
    for linkName in copyDict[self.className]['links']:
      linkedObjects = []
      useLink = True

      linkResult = getattr(self.ccpnObject,linkName)
      
      if type(linkResult) in (type(set()),type(frozenset()),type(tuple())):
        for linkedObject in linkResult:
          if not substToBaseDict.has_key(linkedObject):
            print "  Warning: aborting %s creation - %s object does not exist." % (self.className,linkedObject)
            useLink = False
            break
          else:
            linkedObjects.append(substToBaseDict[linkedObject])

      elif linkResult:
        if not substToBaseDict.has_key(linkResult):
          print "  Warning: aborting %s creation - %s object does not exist." % (self.className,linkResult)
          useLink = False
          break
        else:
          linkedObjects = substToBaseDict[linkResult]

      else:
        print "  Warning: no result for link %s, CCPN object %s!" % (linkName,self.ccpnObject)
        continue

      if not useLink:
        break    
      
      # Note: linkedObjects can be single object!
      self.creationKeywds[linkName] = linkedObjects
      
    return useLink
  
  def createNewObject(self,baseUnit,unitText='substituent'):
  
    creator = getattr(baseUnit,"new%s" % self.className)
    #print self.creationKeywds
    self.newCcpnObject = creator(**self.creationKeywds)
    
    objectKeyString = self.newObjectKeyString(self.newCcpnObject)
    print "  Added %s%s from %s unit..." % (self.className,objectKeyString,unitText)
    
    return self.newCcpnObject

  def newObjectKeyString(self,newCcpnObject):
  
    objectKeyString = ""
    
    if hasattr(newCcpnObject,'name'):
      objectKeyString = " %s" % newCcpnObject.name
    elif hasattr(newCcpnObject,'chemAtoms'):
      caStrings = []
      for chemAtom in newCcpnObject.chemAtoms:
        caStrings.append("%s,%d" % (chemAtom.name,chemAtom.subType))
      objectKeyString = " " + '-'.join(caStrings)  
  
    return objectKeyString
    
class CreateChemAtomOrSet(CreateCcpnObject):

  # Bit of a mix between ChemAtom and ChemAtomSet, could be separated but probably not worth it
  
  def setAtomSysName(self,namingSystem):
  
    if namingSystem:
      namingSystem.newAtomSysName(sysName = self.newCcpnObject.name, atomName = self.newCcpnObject.name, atomSubType = self.newCcpnObject.subType)
  
  def setName(self,substUnitIndex):
  
    self.creationKeywds['name'] = "%s_%s" % (self.ccpnObject.name,substUnitIndex)

  def setForcedName(self,atomName):
  
    self.creationKeywds['name'] = atomName
             
  def setChemAtomLinks(self,chemAtomOrSet,substToBaseDict,substRemoveAtomName,ignoreUnmapped=False):

    # Need to pass in chemAtom info for atomSet creation
    redoChemAtomSets = []
    
    if self.className == 'ChemAtomSet':
      if chemAtomOrSet.chemAtoms:
        self.creationKeywds['chemAtoms'] = []
        for chemAtom in chemAtomOrSet.chemAtoms:
          if ignoreUnmapped and not substToBaseDict.has_key(chemAtom):
            return None
          self.creationKeywds['chemAtoms'].append(substToBaseDict[chemAtom])

      else:
        # Have to create all normal chemAtomSets before these
        redoChemAtomSets.append(chemAtomOrSet)
        
    return redoChemAtomSets

  def setChemAtomSetLinks(self,substToBaseDict):
    
    self.creationKeywds['chemAtomSets'] = []
    for chemAtomSet in chemAtomSet.chemAtomSets:
      self.creationKeywds['chemAtomSets'].append(substToBaseDict[chemAtomSet])

  def checkExistence(self,baseUnit):        
            
    checkChemAtomOrSet = baseUnit.findFirstChemAtom(**self.creationKeywds)
    if not checkChemAtomOrSet:
      checkChemAtomOrSet = baseUnit.findFirstChemAtomSet(**self.creationKeywds)
      if checkChemAtomOrSet:
        print "  Error: %s %s, subType %d already exists!!!" % (self.className,checkChemAtomOrSet.name,checkChemAtomOrSet.subType)
        return True
    
    return False

  def newObjectKeyString(self,newCcpnObject):
  
    return (" %s, subType %d" % (newCcpnObject.name,newCcpnObject.subType))


def deleteVoidChemCompVar(chemCompVar,changeAttrName,atomName,baseUnit,baseUnitCcc,actionType = 'removal'):

  #
  # Initialise
  #

  attrNames = ['linking','descriptor']
  attrElements = {'descriptor': ['deprot','prot','link'], 'linking': ['link']}
  
  fixedAttrName = attrNames[not attrNames.index(changeAttrName)]
  searchDict = {fixedAttrName: getattr(chemCompVar,fixedAttrName)}

  changeAttr = getattr(chemCompVar,changeAttrName)
  changeAttrDict = getDescriptorDict(changeAttr)
  
  #
  # Get the coord var info, if present
  #
  
  chemCompVarCoord = baseUnitCcc.findFirstChemCompVarCoord(linking = chemCompVar.linking, descriptor = chemCompVar.descriptor)
  
  if chemCompVarCoord:
    addText = " and chemCompVarCoord"
  else:
    addText = ""

  #
  # Now start to check if linking/descriptor has to be changed, or the var deleted
  #
  
  for attrElement in attrElements[changeAttrName]:
    attrElementValues = [atomName]
    deleteLinkAtoms = False
  
    # Need to look for linkCode here... could be different from chemAtom.name!
    if changeAttrName == 'linking' or attrElement == 'link':
      chemAtom = chemCompVar.findFirstChemAtom(name = atomName)
      if chemAtom and chemAtom.boundLinkEnds:
        attrElementValues = [le.linkCode for le in chemAtom.boundLinkEnds]
        deleteLinkAtoms = True

    for attrElementValue in attrElementValues:
      if changeAttrDict.has_key(attrElement) and attrElementValue in changeAttrDict[attrElement]:
        changeAttrDict[attrElement].pop(changeAttrDict[attrElement].index(attrElementValue))
        changeAttr = getDescriptorFromDict(changeAttrDict)

        # Hack to set linking name correctly!
        if changeAttrName == 'linking' and changeAttr == 'neutral':
          changeAttr = 'none'

        searchDict[changeAttrName] = changeAttr

        # Only rename existing if doesn't exist yet!
        if not baseUnit.findFirstChemCompVar(**searchDict):
                   
          print "  Renaming %s for chemCompVar%s %s,%s to %s" % (changeAttrName,addText,chemCompVar.linking,chemCompVar.descriptor,changeAttr)

          # Also clean up links to any linkAtoms previously involved
          if deleteLinkAtoms:
            linkEnd = chemCompVar.findFirstLinkEnd(linkCode = attrElementValue)
            for linkAtom in (linkEnd.boundLinkAtom,linkEnd.remoteLinkAtom):
              if linkAtom:
                chemCompVar.removeChemAtom(linkAtom)
                print "    - removed linkAtom '%s' based on removed %s %s." % (linkAtom.name,attrElementValue,changeAttrName)
                
                if chemCompVarCoord:
                  cccLinkAtom = chemCompVarCoord.findFirstChemAtomCoord(name = linkAtom.name, subType = linkAtom.subType)
                  if cccLinkAtom:
                    chemCompVarCoord.removeChemAtomCoord(cccLinkAtom)               
          #
          # Now hack new info into API. Have to set a temporary 'link' for validity check, also make sure old connection is removed!
          #
          
          oldKey = (chemCompVar.linking,chemCompVar.descriptor)
          
          del(chemCompVar.parent.__dict__['chemCompVars'][oldKey])
          chemCompVar.__dict__[changeAttrName] = changeAttr
          chemCompVar.parent.__dict__['chemCompVars'][(chemCompVar.linking,chemCompVar.descriptor)] = chemCompVar
          
          if chemCompVarCoord:
            if chemCompVarCoord.parent.__dict__['chemCompVarCoords'].has_key(oldKey):
              del(chemCompVarCoord.parent.__dict__['chemCompVarCoords'][oldKey])
            chemCompVarCoord.__dict__[changeAttrName] = changeAttr
            chemCompVarCoord.parent.__dict__['chemCompVarCoords'][(chemCompVar.linking,chemCompVar.descriptor)] = chemCompVarCoord
            
        else:
          chemCompVar.delete()
          
          if chemCompVarCoord:
            chemCompVarCoord.delete()

          print "  Removed chemCompVar%s %s,%s based on %s atom %s." % (addText,chemCompVar.linking,chemCompVar.descriptor,atomName,actionType)
          

def deleteLinkEnd(chemComp,atomName,actionType = 'removal'):

  linkEnd = chemComp.findFirstLinkEnd(linkCode = atomName)
  
  if linkEnd:
    linkEnd.delete()
    linkEnd.boundLinkAtom.delete()
    print "  Removed linkEnd %s based on %s atom %s." % (atomName,atomName,actionType)

def getAtomOrigCoords(rawChemComp,chemAtomName):

  coords = None

  for rawAtom in rawChemComp.atoms:

    if rawAtom.name == chemAtomName:
      coords = (rawAtom.x, rawAtom.y, rawAtom.z)
      break

  if not coords:
    print "  Warning: no coordinate for atom name %s." % (coordSystem,chemAtomName)
    
  return coords

def copyBaseToModifiedFile(project,newBaseUnit,origBaseUnit,testMode,repository,replace = False):
  
  #
  # Reset the still empty new base unit parameters
  #
  
  newBaseUnit.__dict__['isLoaded'] = False
  newBaseUnit.__dict__['isModified'] = False
  newBaseUnit.__dict__['isStored'] = True
   
  #
  # Some pattern for identifying elements that have to be changed in the original file
  #
  
  ccpCodePatt = re.compile("ccpCode=\"([^ ]+)\"")
  namePatt    = re.compile("name=\"([^ ]+)\"")
  guidPatt    = re.compile("guid=\"([^ ]+)\"")
  
  #
  # Get the template (original) base unit file location
  #
  
  print origBaseUnit, testMode
  
  origFilePath = findChemCompOrCoordFilePath(origBaseUnit, testMode = testMode)

  #
  # Get the repository path for the new base unit
  #
  
  packageName = newBaseUnit.packageLocator.targetName
  repPath = repository.getFileLocation(packageName)
  
  if not os.path.exists(repPath):
    os.makedirs(repPath)
  
  #
  # Find out the file path for the new base unit 
  #
  
  baseFileName = "%s" % (newBaseUnit.molType)
  if newBaseUnit.className == 'ChemCompCoord':
    baseFileName = "%s+" % (getCcpFileString(newBaseUnit.sourceName)) + baseFileName
  
  newFilePath =  os.path.join(repPath,"%s+%s+" % (baseFileName,getCcpFileString(newBaseUnit.ccpCode)))
  newFilePath += "%s.xml" % newBaseUnit.guid
     
  #
  # Read in the template file and modify it for the new base unit
  #
     
  fin = open(origFilePath)
  lines = fin.readlines()
  fin.close()
  
  for lindex in range(len(lines)):
    line = lines[lindex]
    
    if line.count('<CHEM.NonStdChemComp') or line.count("<CCCO.ChemCompCoord"):
    
      ccpCodeSearch = ccpCodePatt.search(line)
      line = line.replace('ccpCode="%s"' % ccpCodeSearch.group(1), 'ccpCode="%s"' % newBaseUnit.ccpCode)
      
      # Only relevant if name set to ccpCode...
      nameSearch = namePatt.search(line)
      if nameSearch and nameSearch.group(1) == ccpCodeSearch.group(1):
        line = line.replace('name="%s"' % nameSearch.group(1), 'name="%s"' % newBaseUnit.ccpCode)
      
      guidSearch = guidPatt.search(line)
      line = line.replace('guid="%s"' % guidSearch.group(1), 'guid="%s"' % newBaseUnit.guid)
      
      lines[lindex] = line
      
      break

  #
  # If new file already exists, only replace it if replace flag set to True
  #

  if os.path.exists(newFilePath):
    if replace:
      os.remove(newFilePath)
    else:
      print "  Modified chemComp already exists - aborting..."
      return False
  
  #
  # Write out the new base unit file
  #
  
  fout = open(newFilePath,'w')
  fout.write(''.join(lines))
  fout.close()

  return True

def addSubstituentToBaseUnit(baseUnitCcpCode,
                             baseUnitMolType,
                             testMode,
                             mergeInfoList,
                             coordSystem,
                             substRemoveAtomName='SUB',
                             saveData = True,
                             replace=False,
                             namingSystemName = None,
                             resetGlycoCtCode = False):
  
  #
  # Set directories to read/write data from/to
  #

  if testMode:
    chemCompCoordDataDir = testChemCompCoordDataDir
    chemCompDataDir = testChemCompDataDir
  else:
    chemCompCoordDataDir = editChemCompCoordDataDir
    chemCompDataDir = editChemCompDataDir
    
  #
  # Now start setting up project
  #

  substIndexPatt = re.compile("(\d+)")
  substituentDir = os.path.join(origMol2DataDir,'subst')
  
  ccpCode = baseUnitCcpCode
  
  substituentList = []
  for mergeInfoDict in mergeInfoList:
    substituent = mergeInfoDict['substituent']
    if not substituent in substituentList:
      substituentList.append(substituent)
      
    ccpCode += ":%s_%s" % (mergeInfoDict['baseBindingAtomName'],substituentInfo[substituent]['shortCode'])
      
  project = Implementation.MemopsRoot(name = 'chemComp')
  project.currentUserId = 'ccpnRef'
  project.currentChemElementStore = project.findFirstChemElementStore()
    
  #
  # First import all relevant subsituent mol2 files
  # TODO could in principle have these 'pre-imported' in CCPN, in special dir...
  #      note though that the SUBST atoms have to be removed in that case...
  #
   
  mol2Format = Mol2Format(project, guiParent = None, allowPopups = False)
  
  chemComps = {}
  rawChemComps = {}
  
  for substituent in substituentList:
    
    fileName = os.path.join(substituentDir,"%s.mol2" % substituent)

    ccs = mol2Format.readChemComps(
                             fileName,
                             ccpCodes = [substituent],
                             saveChemComp = False,
                             minimalPrompts = True)
  
    chemComps[substituent] = ccs[0]
    rawChemComps[substituent] = mol2Format.rawChemComp
    
  #
  # Get the original base unit information 
  #
    
  origBaseUnit = getChemComp(project,baseUnitMolType,baseUnitCcpCode,chemCompArchiveDir = chemCompDataDir, copyFile = False)
  
  #
  # If replacing, keep GUID of original as default
  #
  
  creationKeywds = {}
  
  if replace:

    if testMode:
      dataDir = testChemCompDataDir
    else:
      dataDir = editChemCompDataDir

    (existingGuid,existingFile) = findExistingChemCompInfo(dataDir,ccpCode,baseUnitMolType)
    
    if existingGuid:
      creationKeywds['guid'] = existingGuid

  #
  # Now create the new base unit, copy the original base unit file, change the guid in the file, then load it
  #

  baseUnit = project.newNonStdChemComp(molType = baseUnitMolType, ccpCode = ccpCode, **creationKeywds)
  repository = project.findFirstRepository(name = 'userData')
  
  isNewFile = copyBaseToModifiedFile(project,baseUnit,origBaseUnit,testMode,repository,replace=replace)  
  
  if not isNewFile:
    return
   
  #print baseUnit.chemAtoms
  
  origBaseUnitCcc = getChemCompCoord(project,coordSystem,baseUnitMolType,baseUnitCcpCode, chemCompCoordArchiveDir = chemCompCoordDataDir, copyFile = False)
  if not origBaseUnitCcc:
    raise("Error: no coordinates available for %s!" % baseUnitCcpCode)
  else:
    #
    # If replacing, keep GUID of original as default
    #
    
    creationKeywds = {}

    if replace:

      if testMode:
        dataDir= testChemCompCoordDataDir
      else:
        dataDir= editChemCompCoordDataDir

      (existingGuid,existingFile) = findExistingChemCompCoordInfo(dataDir,coordSystem,ccpCode,baseUnitMolType)
      
      if existingGuid:
        creationKeywds['guid'] = existingGuid

    baseUnitCcc = project.newChemCompCoord(sourceName = coordSystem, molType = baseUnitMolType,ccpCode = ccpCode, **creationKeywds)
    copyBaseToModifiedFile(project,baseUnitCcc,origBaseUnitCcc,testMode,repository,replace=replace)
    #print baseUnitCcc.chemAtomCoords
  
  #
  # Reset the glycoCtCode for the new base unit (if relevant)
  #
  
  if resetGlycoCtCode:
  
    substGlycoCtText = ""    
    substGlycoCtInfo = {'RES': [], 'LIN': []}
    resIndex = 2
    linIndex = 1
    
    # Note: this assumes the mergeInfoList is ordered!
    
    for mergeInfoDict in mergeInfoList:
   
      baseBindingAtomName = mergeInfoDict['baseBindingAtomName']
      removeBaseAtomNames = mergeInfoDict['removeBaseAtomNames']
      substituent = mergeInfoDict['substituent']
      
      substGlycoCtInfo['RES'].append("%ds:%s" % (resIndex,substituent))
      
      parentAtomIndex = int(baseBindingAtomName[-1])
      
      if baseBindingAtomName[0] == 'O':
        parentSubstitutionType = 'o'
      else:
        parentSubstitutionType = 'd'
      
      substGlycoCtInfo['LIN'].append("%d:%d%s(%d+%d)%d%s" % (linIndex,
                                                                1,
                                                                parentSubstitutionType,
                                                                parentAtomIndex,
                                                                1,
                                                                resIndex,
                                                                'n'))
      
      resIndex +=1
      linIndex +=1  
    
    
    for tmpStr in substGlycoCtInfo['RES']:
      substGlycoCtText += "\n" + tmpStr
    substGlycoCtText += "\nLIN"
    for tmpStr in substGlycoCtInfo['LIN']:
      substGlycoCtText += "\n" + tmpStr
    
    newBaseGlycoCtCode = origBaseUnit.baseGlycoCtCode + substGlycoCtText   

    project.override = True
    try:
      baseUnit.baseGlycoCtCode = newBaseGlycoCtCode
      
      for chemCompVar in baseUnit.chemCompVars:
        chemCompVar.glycoCtCode = chemCompVar.glycoCtCode + substGlycoCtText
        
    finally:
      project.override = False
      
    print "Setting base GlycoCT code to:\n\n%s\n" % newBaseGlycoCtCode
    print

  #
  # Set naming system
  #
  
  namingSystem = None
  if namingSystemName:
    namingSystem = baseUnit.findFirstNamingSystem(name = namingSystemName)
    if not namingSystem:
      namingSystem = baseUnit.newNamingSystem(name = namingSystemName)
      print "Created new naming system %s" % namingSystemName
  
  #
  # Add substituent info, remove atoms from base unit
  #
  
  addVariants = {}
  
  for mergeInfoDict in mergeInfoList:
  
    #
    # 0. Get the info from the mergeInfoDict, print a comment
    #
   
    baseBindingAtomName = mergeInfoDict['baseBindingAtomName']
    removeBaseAtomNames = mergeInfoDict['removeBaseAtomNames']
    substituent = mergeInfoDict['substituent']
    renameSubstituentAtomNames = mergeInfoDict['renameSubstituentAtomNames']
    
    baseBindingAtoms = baseUnit.findAllChemAtoms(name = baseBindingAtomName)

    newBondType = substituentInfo[substituent]['bondType']
    newStereochem = substituentInfo[substituent]['stereochem']
    
    print
    print drawBox("Creating link between base atom %s to substituent %s" % (baseBindingAtomName,substituent))
    print
    
    #
    # 1. Set the substUnitIndex - this is the identifier that is added to the substituents
    #    when part of the base chemComp. It is taken from whichever number is part of the baseBindingAtomName
    #
    #    TODO: should be molType specific - use A,B,G,... for amino acids!!
    #
    
    substUnitSearch = substIndexPatt.search(baseBindingAtomName)
    substUnitIndex = substUnitSearch.group(1)
    
    #
    # 2. Remove relevant atoms from the base unit, all subtypes
    #    Keep track of atom directly linked to the baseBindingAtom for recalculating coordinates!
    #
    
    baseAtomCoords = {}
    
    for removeBaseAtomName in removeBaseAtomNames:
    
      removeBaseAtoms = baseUnit.findAllChemAtoms(name = removeBaseAtomName)
      
      for removeBaseAtom in removeBaseAtoms:
        
        # Search if bound to the baseBindingAtom        
        isBoundToBaseBindingAtom = False
        for chemBond in baseUnit.chemBonds:
          bondChemAtoms = list(chemBond.chemAtoms)
          if removeBaseAtom in bondChemAtoms:
            otherBondChemAtom = bondChemAtoms[not bondChemAtoms.index(removeBaseAtom)]
            if otherBondChemAtom in baseBindingAtoms:
              isBoundToBaseBindingAtom = True
              break
        
        # Track the coordinates if bound to the baseBindingAtom
        if isBoundToBaseBindingAtom:
          for baseBindingAtom in baseBindingAtoms:

            baseAtomCoords[baseBindingAtom] = {}            
            baseBindingAtomCoords = baseUnitCcc.findAllChemAtomCoords(chemAtom = baseBindingAtom)
            
            for baseBindingAtomCoord in baseBindingAtomCoords:
              baseCoord = (baseBindingAtomCoord.x,baseBindingAtomCoord.y,baseBindingAtomCoord.z)

              for chemCompVarCoord in baseBindingAtomCoord.chemCompVarCoords:
                # Get the coordinates, if any
                removeBaseAtomCoord = chemCompVarCoord.findFirstChemAtomCoord(chemAtom = removeBaseAtom)
                if removeBaseAtomCoord:
                  baseBoundCoord = (removeBaseAtomCoord.x,removeBaseAtomCoord.y,removeBaseAtomCoord.z)
                  jointCoords = (baseCoord,baseBoundCoord)
                  
                  if not baseAtomCoords[baseBindingAtom].has_key(jointCoords):
                    baseAtomCoords[baseBindingAtom][jointCoords] = []
                  baseAtomCoords[baseBindingAtom][jointCoords].append(chemCompVarCoord)
        
        # Now start deleting on chemComp and chemCompCoord levels      
        removeBaseAtomCoords = baseUnitCcc.findAllChemAtomCoords(chemAtom = removeBaseAtom)
        for removeBaseAtomCoord in removeBaseAtomCoords:
          removeBaseAtomCoord.delete()

        removeBaseAtom.delete()
        print "  Removed atom %s, subType %d from base unit..." % (removeBaseAtomName,removeBaseAtom.subType)
        
        for namingSystem in baseUnit.namingSystems:
          for asn in namingSystem.atomSysNames:
            if asn.atomName == removeBaseAtomName:
              asn.delete()
                    
        #
        # 2.1 Also rename all linkEnds and chemCompVars that have this atom in the descriptor - are now irrelevant
        #     
        
        for chemCompVar in baseUnit.chemCompVars:
          deleteVoidChemCompVar(chemCompVar,'descriptor',removeBaseAtomName,baseUnit,baseUnitCcc)
          deleteVoidChemCompVar(chemCompVar,'linking',removeBaseAtomName,baseUnit,baseUnitCcc)
          deleteLinkEnd(baseUnit,removeBaseAtomName)
                    
    #
    # 3. Add the substituent info to the base unit
    #
    # Currently this works off the mol2 file. Could add, e.g. SUB_C, SUB_O, SUB_N, depending on substituted atom... 
    # for better coordinates later on. TODO: try to implement this!!!
    #
    # TODO: look into avoiding mol2 step, just put into CCPN in temporary library, use linkAtoms as linking ones.
    # These can then also be identified by elementSymbol (O,C,N,...)
    #
    
    substUnit = chemComps[substituent]
    rawSubstChemComp = rawChemComps[substituent]
    
    #
    # 3.1 Initialise information
    #       - determine which atom to remove from substituent (will need this for coordinates further down though)
    #       - get single bond coming from this atom to identify the binding atom on the substituent side
    #       - create substToBaseDict dictionary that maps substituent objects to newly created chemComp unit objects
    #       - deal with the coordinates
    
    substRemoveAtom = substUnit.findFirstChemAtom(name = substRemoveAtomName)
    chemBond = substRemoveAtom.findFirstChemBond() # Should only have ONE single bond!
    substBindingAtom = getOtherAtom(substRemoveAtom,chemBond)
    
    substCoords = [getAtomOrigCoords(rawSubstChemComp,substRemoveAtomName),getAtomOrigCoords(rawSubstChemComp,substBindingAtom.name)]
    substCoordsBaseAtoms = [None,None]
    
    substToBaseDict = {}

    #
    # 3.2 Create chemAtoms and chemAtomSets
    #

    redoChemAtomSets = []
    substChemAtoms = []
    
    for chemAtomOrSet in substUnit.sortedChemAtoms() + substUnit.sortedChemAtomSets(): 
      
      # Ignore atoms to be substituted
      if chemAtomOrSet.name[:len(substRemoveAtomName)] == substRemoveAtomName:
        continue
      
      createChemAtomOrSet = CreateChemAtomOrSet(chemAtomOrSet)
      
      if renameSubstituentAtomNames.has_key(chemAtomOrSet.name):
        createChemAtomOrSet.setForcedName(renameSubstituentAtomNames[chemAtomOrSet.name])
      else:
        createChemAtomOrSet.setName(substUnitIndex)
      
      redoChemAtomSets.extend(createChemAtomOrSet.setChemAtomLinks(chemAtomOrSet,substToBaseDict,substRemoveAtomName)) # Only relevant for chemAtomSets!
      
      if createChemAtomOrSet.checkExistence(baseUnit):
        continue
        
      newChemAtomOrSet = createChemAtomOrSet.createNewObject(baseUnit)
      createChemAtomOrSet.setAtomSysName(namingSystem)
      
      substToBaseDict[chemAtomOrSet] = newChemAtomOrSet
      
      # Keep track of new chemAtoms, also track coordinates
      if newChemAtomOrSet.className in ('LinkAtom','ChemAtom'):
        substChemAtoms.append(newChemAtomOrSet)
        
        if chemAtomOrSet == substBindingAtom:
          baseSubstBindingAtom = newChemAtomOrSet
          substCoordsBaseAtoms[1] = newChemAtomOrSet
        else:
          substCoordsBaseAtoms.append(newChemAtomOrSet)
          substCoords.append(getAtomOrigCoords(rawSubstChemComp,chemAtomOrSet.name))
    
    
    #   
    # 3.2.1 Now set chemAtomSets that are linked to other chemAtomSets - have to do this later because
    #       otherwise could have not been created yet. Could in principle do this in above loop by
    #       organising list of chemAtomSets, but this is easier.
    #
    
    for chemAtomSet in redoChemAtomSets:
      createChemAtomOrSet = CreateChemAtomOrSet(chemAtomSet)

      if renameSubstituentAtomNames.has_key(chemAtomSet.name):
        createChemAtomOrSet.setForcedName(renameSubstituentAtomNames[chemAtomSet.name])
      else:
        createChemAtomOrSet.setName(substUnitIndex)

      createChemAtomOrSet.setChemAtomSetLinks()
      if createChemAtomOrSet.checkExistence(baseUnit):
        continue        
      newChemAtomSet = createChemAtomOrSet.createNewObject(baseUnit)
      createChemAtomOrSet.setAtomSysName(namingSystem)
      substToBaseDict[chemAtomSet] = newChemAtomSet
    
    #
    # 3.3 Now add all other objects connected to the chemComp (see global chemCompLinkList)
    #
    # TODO: this is not fully functional - have to add specific settings for some links (but probably never required, so will wait)
    #
    
    for ccLinkName in chemCompLinkList:
      for substObject in getattr(substUnit,ccLinkName):
        createCcpnObject = CreateCcpnObject(substObject)
        if not createCcpnObject.setLinks(substToBaseDict):
          continue
        newCcpnObject = createCcpnObject.createNewObject(baseUnit)
        substToBaseDict[substObject] = newCcpnObject

    """            
      # TODO: LINKS THAT REQUIRE SPECIAL TREATMENT/MERGING WITH EXISTING INFO
      # 'chemCompVars', 'applicationData',
      
      # TODO: object that require special treatment (have to be renamed, ...): chemAtomSysNames, chemCompSysNames (?)
    """
      
    #
    # 4. Connect base binding atom(s) to the new subsituent atoms
    #
        
    for baseBindingAtom in baseBindingAtoms:
    
      #
      # 4.1 First remove all chemCompVars and linkEnds that have linking vars including this atom!
      #
      
      for chemCompVar in baseBindingAtom.chemCompVars:
        deleteVoidChemCompVar(chemCompVar,'linking',baseBindingAtomName,baseUnit,baseUnitCcc,actionType = 'binding')
        deleteLinkEnd(baseUnit,baseBindingAtomName,actionType = 'binding')
              
      #
      # 4.2 Calculate the new coordinates for the substituent atoms
      #
      
      coordsForChemCompVar = {}
      
      if baseAtomCoords.has_key(baseBindingAtom):
        for jointCoords in baseAtomCoords[baseBindingAtom]:
          newSubstCoords = tuple([tuple(coord) for coord in superposeNewVectorsOnOld(jointCoords,substCoords)])
          coordsForChemCompVar[newSubstCoords] = baseAtomCoords[baseBindingAtom][jointCoords]
          
      #
      # 4.3 Now create the bond between the new substituent and the existing base unit
      #
    
      baseUnit.newChemBond(chemAtoms = (baseBindingAtom,baseSubstBindingAtom), bondType = newBondType, stereochem = newStereochem)
      
      #print "  Creating chemBond between base atom %s and subsituent atom %s" % (baseBindingAtom.name,baseSubstBindingAtom.name)
      
      for chemCompVar in baseBindingAtom.chemCompVars:
        #print "   ", baseSubstBindingAtom.name, baseSubstBindingAtom.subType,chemCompVar.linking,chemCompVar.descriptor
        for substChemAtom in substChemAtoms:
          chemCompVar.addChemAtom(substChemAtom)
        # TODO: have to RENAME the chemCompVar in case there are variants for the substituent!!
        
      #
      # 4.4 Finally set the new coordinates... ignore first one (is the substituted atom from the substituent!)
      #
      
      for newSubstCoords in coordsForChemCompVar.keys():
        for i in range(1,len(substCoordsBaseAtoms)):
          newSubstCoord = newSubstCoords[i]
          substChemAtom = substCoordsBaseAtoms[i]
          
          substChemAtomCoord = baseUnitCcc.findFirstChemAtomCoord(chemAtom = substChemAtom, x = newSubstCoord[0], y = newSubstCoord[1], z = newSubstCoord[2])          
          if not substChemAtomCoord:
            substChemAtomCoord = baseUnitCcc.newChemAtomCoord(name = substChemAtom.name, subType = substChemAtom.subType, x = newSubstCoord[0], y = newSubstCoord[1], z = newSubstCoord[2])
            
          for cccv in coordsForChemCompVar[newSubstCoords]:
            if not cccv.isDeleted:
              cccv.addChemAtomCoord(substChemAtomCoord)
        
        
      print "  Connected new substituent chemAtom %s,%s to base unit atom %s,%s, and included in relevant chemCompVars" % (
      
                        baseSubstBindingAtom.name,
                        baseSubstBindingAtom.subType,
                        baseBindingAtom.name,
                        baseBindingAtom.subType)
    
    #
    # NEXT ON LIST: make sure the chemCompVars make sense if there are any for the substituent!!!
    #
    # ALSO try this on amino acids when working for carbs - need split code and rename some variables though
    #
    
  #
  # Check validity and save
  #
  
  baseUnit.checkAllValid(complete = True)

  if baseUnitCcc:
    baseUnitCcc.checkAllValid(complete = True)
  
  if saveData:

    (filePath,existingFilePath) = saveTemporaryChemCompOrCoord(baseUnit,testMode = testMode)
    consolidateTemporaryChemCompOrCoord(baseUnit,filePath,existingFilePath,testMode = testMode,replace=replace)

    if baseUnitCcc:
      (filePath,existingFilePath) = saveTemporaryChemCompOrCoord(baseUnitCcc,testMode = testMode)
      consolidateTemporaryChemCompOrCoord(baseUnitCcc,filePath,existingFilePath,testMode = testMode,replace=replace)

if __name__ == '__main__':

  from pdbe.chemComp.Constants import testChemCompDataDir

  print "Creating in test directory!"
  chemCompSaveDir = testChemCompDataDir

  #
  # TODO TODO add GUI interface to set EuroCarbDB atom names..
  #
  # Also need GUI to select chemComp to modify, then select which substituent and where it should go...
  #
  
  if False:
    
    #
    # Carbohydrate test
    #

    baseUnitMolType = 'carbohydrate'
    baseUnitCcpCode = 'dglc-hex-1-5'

    coordSystem = 'EuroCarbDB'

    mergeInfoList = []

    # TODO: make this user interaction! Or script...

    # OMe on C1
    if False:
      mergeInfoList.append({  # This is a methoxy modification at the reducing end, directly on the C1!!

          'baseBindingAtomName': 'C1',
          'removeBaseAtomNames': ('O1','H1O'),
          'substituent': 'methoxy',
          'renameSubstituentAtomNames': {'OM': 'O1'}

          })

    # NAc on C2 (ala GlcNAc)
    if True:
      mergeInfoList.append({

          'baseBindingAtomName': 'C2',
          'removeBaseAtomNames': ('O2','H2O'),
          'substituent': 'n-acetyl',
          'renameSubstituentAtomNames': {}

          })

  if True:
    
    #
    # Protein test
    #

    baseUnitMolType = 'protein'
    baseUnitCcpCode = 'Ala'

    coordSystem = 'corina'

    mergeInfoList = []

    mergeInfoList.append({ 

        'baseBindingAtomName': 'CB',
        'removeBaseAtomNames': ('HB1'),
        'substituent': 'methyl',
        'renameSubstituentAtomNames': {'CM': 'CG', 'HM1': 'HG1', 'HM2': 'HG2', 'HM3': 'HG3'}

        })

        
  # Test data
  # mergeInfoList.append({'baseBindingAtomName': 'O3',
  #                  'removeBaseAtomNames': ('H3',),
  #                  'substituent': 'methyl'}) 

  
  addSubstituentToBaseUnit(baseUnitCcpCode,baseUnitMolType,chemCompSaveDir,mergeInfoList,coordSystem,saveData = False)
  
