from pdbe.nmrStar.IO.NmrStarHandler import NmrStarFullReaderFile, NmrStarFormat

from pdbe.nmrStar.IO.Util import getCcpnVersion, getLatestNmrStarVersion
from pdbe.nmrStar.IO.Util import getCcpn2NmrStar, getNmrStarDict

import os, re

from memops.api import Implementation
from ccp.format.general.Constants import defaultSeqInsertCode, defaultMolCode
from ccp.general.Util import findAllSysNamesByChemAtomOrSet, findChemCompSysName

from ccpnmr.format.general.Constants import assign_kw

#
# Notes on conversion project (writing data read in with NmrStarFormat):
#
#  - due to different saveframe category names, some data disappears
#
#
#
# TODO: for read have to store tags WITH name information otherwise no way of mapping!!
#       need a _Names applData setting for _tags_ and each _table_...
#
#
# Check for TODO and WARNING in this script...
#

class NmrStarExport:

  """

  Generic NMR-STAR export script - works together with mapping class
  (Ccpn_To_NmrStar) and Python dictionary derived from the original NMRSTAR
  dictionary (nmrStarDict)

  """

  setAtomsTags = ("CUSTOM_setAtoms","CUSTOM_setAtoms_individual")
 
  def __init__(self,nmrEntry, nmrStarVersion = None, forceEntryId = None):

    #
    # Get the relevant NMR-STAR dictionary and the relevant conversion class
    #
    
    ccpnVersion = getCcpnVersion()
    
    if not nmrStarVersion:
      nmrStarVersion = getLatestNmrStarVersion()
    
    print "  Using CCPN version %s, NMR-STAR version %s." % (ccpnVersion,nmrStarVersion)
    
    self.nmrEntry = nmrEntry
    self.ccpn2NmrStar = getCcpn2NmrStar(ccpnVersion,nmrStarVersion,exportClass = self)
    self.nmrStarDict = getNmrStarDict(nmrStarVersion)

    #
    # Make sure chemElementStore is set - otherwise will crash eventually
    # TODO: this can be removed once the API does this automatically!
    #
    
    if not self.nmrEntry.root.currentChemElementStore:
      self.nmrEntry.root.currentChemElementStore = self.nmrEntry.root.findFirstChemElementStore()
    
    #
    # These are dictionaries that track CCPN object to STAR ID mapping
    #
    # starKeys: dictionary with the keys for the star tags.
    #           
    #    Key: (tableName,tagName,). Multiple tagNames are possible as
    #         key, this is to handle tables/loops that correspond to
    #         multiple layers of CCPN data.
    #                    
    #    Value: [(ccpnObject,starID,parentLoopDict),...]
    #           The parentLoopDict contains information about the 'parent' object
    #
    #
    # starKeyToCcpnVar: dictionary that maps the keys from starKeys to
    #                   the corresponding CCPN class name (e.g. ConstraintItem)
    #
    # starTagToStarKey: dictionary that maps the real star tag name to the
    #                   corresponding key from starKeys.
    #
    # messages: keeps track of error and warning messages.
    #

    self.starKeys = {}
    self.starKeyToCcpnVar = {}
    self.starTagToStarKey = {}
    self.messages = {}
    
    self.titleList = []
   
    #
    # Regular expression stuff
    #
   
    self.patt = {}
    self.patt['getLocalValue'] = re.compile("^(.+)\=LOCAL$")
    self.patt['customFlag'] = re.compile("^CUSTOM")
    self.patt['insideBrackets'] = re.compile("(\(.*\))")

    #
    # This is for keeping track of CCPN objects
    #
    # Only the current relevant CCPN objects are stored in self.ccpnVar!
    #
    
    self.ccpnVar = {'nmrEntry': nmrEntry,
                    'resonances': []}               # Special case!
                    
    #
    # To keep track of application data objects
    #
    
    self.origVar = {}
    self.resonanceAuthorPdbInfo = {}
    self.coordinateInfoMap = {}   # This is specific for tagging resonances with PDB atom name information
                      
    #
    # Set up the main link (to the nmrEntry object)
    #
  
    for nmrEntryLink in self.ccpn2NmrStar.nmrEntryLinkList:
      ccpnMap = 'nmrEntry.%s' % nmrEntryLink
      self.ccpnVar[nmrEntryLink] = self.getCcpnObject(ccpnMap)
      
    #
    # Hardset the Entry ID if required
    #
    
    if forceEntryId:
      # Is handled in special way in self.starKeys - topObject so effectively None
      matchKeyName = ('Entry', 'ID')
      self.starKeys[matchKeyName] = ({None: (forceEntryId,None)},nmrEntry)
      self.starKeyToCcpnVar[matchKeyName] = ''
      
    #
    # Dictionaries to speed things up...
    #
    
    self.chemAtomSysNamesDict = {}
    
  def createFile(self,nmrStarFileName,verbose=True,profile=False,excludeSaveFrames=('distance_constraints',)):
    
    if profile:
      import hotshot
    
    self.verbose = verbose
    
    if excludeSaveFrames:
      self.excludeSaveFrames = excludeSaveFrames
    else:
      self.excludeSaveFrames = []
        
    #
    # Create the nmrStarFile object for writing
    #

    self.nmrStarFileName = nmrStarFileName
    
    self.nmrStarFile = NmrStarFullReaderFile(self.nmrStarFileName)

    #
    # Preparation stage: create resonances for NMR projects only, on individual level. This goes VIA NmrStarFormat!!!
    # Basically a hack to handle individual atom names for resonances.
    #

    self.createNmrStarFormat(initialise = True)

    #
    # First pass - setting the keys and printing error messages (only once)
    #

    if profile:
      prof = hotshot.Profile("%s.firstLoop.stats" % nmrStarFileName)
      prof.runcall(self.loopSaveFrameData,1)
      prof.close()
    else:
      self.loopSaveFrameData(1)
            
    #
    # Intermediate stage: handle the resonances. This goes VIA NmrStarFormat!!!
    #

    self.createNmrStarFormat()

    self.setWaterChainStarKeys()

    #
    # Second pass - now set the actual values
    #

    #self.tempTagInfo = {}
  
    if profile:
      prof = hotshot.Profile("%s.secondLoop.stats" % nmrStarFileName)
      prof.runcall(self.loopSaveFrameData,0)
      prof.close()
    else:
      self.loopSaveFrameData(0)
    
    #
    # Finally set additional saveframes from any application data...
    # this hangs off entry.
    #
    # Also set additional loops - go over all CCPN objects from self.entry?
    #
    # TODO: what if desired that ORIGINAL IDs are being used?! Set earlier?!
    #

    # This code was moved to loopStarLevelObjects, and is done for all save
    # frames with app. data.

    #
    # Print messages
    #
  
    self.printMessages()

  
  def writeFile(self,title = None, topComment = None, verbose = True):
  
    #
    # Writes the star file
    #
    
    keywds = {'verbose': verbose}
    
    if title:
      keywds['title'] = title

    if topComment:
      keywds['topComment'] = topComment

    self.nmrStarFile.writeGeneric(**keywds)


  ##########################
  # Second level functions #
  ##########################

  def printMessages(self):
  
    #
    # Prints out the messages
    #
    
    if self.verbose:

      messageKeys = self.messages.keys()
      messageKeys.sort()
    
      for messageKey in messageKeys:
      
        if self.messages[messageKey] > 1:
          timeString = " (%d times)" % self.messages[messageKey]
        else:
          timeString = ""
      
        print messageKey + "%s" % timeString

  def loopSaveFrameData(self,keyFlag):
  
    #
    # KeyFlag is 1 for first (source) run, 0 for second (setting) run
    #
    
    self.keyFlag = keyFlag
  
    #
    # Loop over the saveframes
    #
    
    ccpn2NmrStarSfNames = self.ccpn2NmrStar.sfDict.keys()
    
    for self.sfCatName in self.nmrStarDict.sfList:
  
      if self.sfCatName not in self.ccpn2NmrStar.sfDict.keys():
        continue
      
      if self.sfCatName in self.excludeSaveFrames:
        self.setMessage(" Excluding saveframe %s." % self.sfCatName)
        continue

      self.writeStarSfDict = self.ccpn2NmrStar.sfDict[self.sfCatName]
      self.nmrStarSfDict = self.nmrStarDict.sfDict[self.sfCatName]
      ccpn2NmrStarSfNames.pop(ccpn2NmrStarSfNames.index(self.sfCatName) )

      if self.writeStarSfDict.has_key('ccpnMap'):

        self.setSaveFrameData(self.nmrStarSfDict['name'])

      elif self.keyFlag:
      
        self.setMessage(" No ccpn mapping for %s" % self.sfCatName)
        
    if ccpn2NmrStarSfNames:
      for sfCatName in ccpn2NmrStarSfNames:
        self.setMessage(" No NMR-STAR equivalent for saveframe %s in mapping file." % sfCatName)
    
    if self.keyFlag:
    
      #
      # This handles foreign keys that should be ignored - needs to be set up so that
      # the foreign keys are only ignored outside of their saveframe (for tables)
      #
      
      self.ignoreForeignLinkStarElements = {}
      self.ignoreForeignLinks = self.ccpn2NmrStar.ignoreForeignLinks[:]
      
      for ignoreForeignLink in self.ccpn2NmrStar.ignoreForeignLinks:
        self.setMessage("  Warning: ignoring foreign link '%s' for all tags." % ignoreForeignLink)
        (starElementName,tagName) = ignoreForeignLink.split('.') # TODO get this from somewhere official!
        if not self.ignoreForeignLinkStarElements.has_key(starElementName):
          self.ignoreForeignLinkStarElements[starElementName] = []
        self.ignoreForeignLinkStarElements[starElementName].append(ignoreForeignLink)

      self.ignorePrimaryKeyStarElements = {}
      self.ignorePrimaryKeys = self.ccpn2NmrStar.ignorePrimaryKeys[:]
      
      for ignorePrimaryKey in self.ccpn2NmrStar.ignorePrimaryKeys:
        self.setMessage("  Warning: ignoring primary key '%s' for all tags." % ignorePrimaryKey)
        (starElementName,tagName) = ignorePrimaryKey.split('.') # TODO get this from somewhere official!
        if not self.ignorePrimaryKeyStarElements.has_key(starElementName):
          self.ignorePrimaryKeyStarElements[starElementName] = []
        self.ignorePrimaryKeyStarElements[starElementName].append(ignorePrimaryKey)


  def createNmrStarFormat(self,initialise = False):
  
    #
    # Create an NmrStarFormat - this is ONLY for running 
    # the .getResonanceAtomLinks() code!
    #
    # Note that because rules are different for different saveframes.
    # In chemical shift saveframes, all individual atoms have to be written out
    # In constraint saveframes, atom sets can be used
    #
  
    nmrStarFormat = NmrStarFormat(self.ccpnVar['nmrEntry'].root)
    nmrStarFormat.useOriginalResNames = False
    nmrStarFormat.molSystem = self.ccpnVar['molSystem']

    if initialise:

      resonances = []
      for measurementList in self.ccpnVar['nmrEntry'].measurementLists:
        for measurement in measurementList.measurements:
          if hasattr(measurement, "resonance"):
            resonances.append(measurement.resonance)
          elif hasattr(measurement, "resonances"):
            resonances.extend(measurement.resonances)
          
      nmrStarFormat.resonancesToLink = resonances
      nmrStarFormat.individualAtoms = True
      nmrStarFormat.individualAtomsIfNoSet = False
      nmrStarFormat.getResonanceAtomLinks()
      self.ccpnVar['nmrStarFormatIndividual'] = nmrStarFormat

    else:
      nmrStarFormat.resonancesToLink = self.ccpnVar['resonances']
      nmrStarFormat.individualAtoms = False
      nmrStarFormat.individualAtomsIfNoSet = True
      nmrStarFormat.getResonanceAtomLinks()
      self.ccpnVar['nmrStarFormatSet'] = nmrStarFormat

    # TODO: this might need to be set earlier depending on whether IDs transferred from appldata...
    self.format = nmrStarFormat.format


  def setWaterChainStarKeys(self):

    matchKeyName = ('Entity_assembly', 'ID')

    if self.starKeys.has_key(matchKeyName):

      for chain in self.starKeys[matchKeyName][0].keys():
        if chain.findFirstResidue().ccpCode == 'Hoh':
          (matchValue,parentDictInfo) = self.starKeys[matchKeyName][0][chain]
        
          for wChain in chain.molecule.chains:
            if wChain != chain:
              self.starKeys[matchKeyName][0][wChain] = (matchValue,parentDictInfo)
            
    
  #########################################################################
  # Functions that set information on saveframe/table level (starElement) #
  #########################################################################
  
  def setSaveFrameData(self,saveFrameName):
  
    #
    # Set up mapping info and IDs
    #
  
    ccpnMap = self.writeStarSfDict['ccpnMap']
    ccpnLoopInfo = self.setupCurrentIDs(ccpnMap,self.writeStarSfDict)

    #
    # Loop over all ccpn objects on saveframe level...
    #

    #print 'NAME: [%s] [%s] [%s]' % (saveFrameName, ccpnMap, ccpnLoopInfo)

    if len(ccpnLoopInfo) > 0:
      self.loopStarLevelObjects('saveFrame',saveFrameName,ccpnLoopInfo,self.writeStarSfDict,self.nmrStarSfDict)


  def setTableData(self,tableName,saveFrame = None):

    #
    # Set up mapping info and IDs
    #
  
    if self.writeStarTableDict.has_key('CONDITIONAL'):
      key = self.writeStarTableDict['CONDITIONAL']
      
      value = self.getCcpnMapValue(key[0])
      
      if value not in key[1]:
        return
        
    ccpnMap = self.writeStarTableDict['ccpnMap']
    ccpnLoopInfo = self.setupCurrentIDs(ccpnMap,self.writeStarTableDict)
    
    #
    # Set up the table at this level
    #
    
    if not self.keyFlag and saveFrame is not None:
      table = saveFrame.setupTable(tableName)
      
    else:
      table = None       

    #
    # Loop over all ccpn objects on table level...
    #

    if len(ccpnLoopInfo) > 0:
      self.loopStarLevelObjects('table',tableName,ccpnLoopInfo,self.writeStarTableDict,self.nmrStarTableDict,starElement = table)


  def loopStarLevelObjects(self,level,starElementName,ccpnLoopInfo,writeStarElementDict,nmrStarElementDict,parentCcpnLoopInfo = None,starElement = None, previousLowestCcpnLoopInfo = None):

    #
    # This loop works for BOTH saveframes and tables
    # Is recognized by level tag
    #
    # Note that previousLowestCcpnLoopInfo is necessary to set IDs that are global to a loop (e.g. atom_site), or to part of the CCPN object loop
    #

    ccpnVarKey = ccpnLoopInfo[0]
    nmrStarKey = ccpnLoopInfo[1]
    currentID = ccpnLoopInfo[2]
    
    #
    # Set unique saveframe titles... for this saveframe category at least
    #
            
    if level == 'saveFrame' and not self.keyFlag:
    
      titleDict = {}
      saveFrameTitles = {}

      for (ccpnObject,lowerCcpnLoopInfo) in ccpnLoopInfo[-1]:

        self.ccpnVar[ccpnVarKey] = ccpnObject

        title = self.searchDataTitle(ccpnObject,nmrStarElementDict)

        #
        # If no applicationData title, find out what to use
        #

        if not title:

          title = self.getCcpnMapValue(writeStarElementDict['title']) # TODO set labels as well?
          
        #
        # Keep track of titles & corresponding objects...
        #
        
        if not titleDict.has_key(title):
          titleDict[title] = []
        
        titleDict[title].append(ccpnObject)

      #
      # Now make sure all titles are unique...
      #
      
      for title in titleDict:

        ccpnObjectList = titleDict[title]
        for i in range(len(ccpnObjectList)):
          tmpTitle = None
          if len(ccpnObjectList) > 1:
            tmpTitle = "%s_%d" % (title,(i+1) )
          else:
            tmpTitle = title
        
          if not tmpTitle or tmpTitle in self.titleList:
            for j in range(1,9999):
              tmpTitle = "%s_%d" % (title,j)
              if tmpTitle not in self.titleList:
                title = tmpTitle
                break
          else:
            title = tmpTitle
        
          saveFrameTitles[ccpnObjectList[i]] = title
          
          self.titleList.append(title)

    #
    # Loop over all ccpn objects
    #

    for i in range(0,len(ccpnLoopInfo[-1]) ):
    
      (ccpnObject,lowerCcpnLoopInfo) = ccpnLoopInfo[-1][i]
      self.ccpnVar[ccpnVarKey] = ccpnObject

      if writeStarElementDict.has_key('addCcpnMap'):
        additionalMappings = writeStarElementDict['addCcpnMap']
        # Setting additional mappings - important for giant saveframes like coordinates to speed things up
        if additionalMappings.has_key(ccpnVarKey):
          for (addCcpnVarKey,addCcpnMapping) in additionalMappings[ccpnVarKey]:
            self.ccpnVar[addCcpnVarKey] = self.getCcpnMapValue(addCcpnMapping)
      
      #
      # There is an inner loop... continue
      #
      
      if lowerCcpnLoopInfo:

        if not parentCcpnLoopInfo:
          parentCcpnLoopInfo = [ccpnLoopInfo]
        else:
          parentCcpnLoopInfo.append(ccpnLoopInfo)

        previousLowestCcpnLoopInfo = self.loopStarLevelObjects(level,starElementName,lowerCcpnLoopInfo,writeStarElementDict,nmrStarElementDict,parentCcpnLoopInfo = parentCcpnLoopInfo,starElement = starElement, previousLowestCcpnLoopInfo = previousLowestCcpnLoopInfo)
        
        # Reset the previous loop info if a nmrStarKeyType given and equal to a ccpnVar name. This means the ID
        # count is reset at this level.
        if writeStarElementDict.has_key('nmrStarKeyType') and writeStarElementDict['nmrStarKeyType'] == ccpnVarKey:
          previousLowestCcpnLoopInfo = None

      else:

        #
        # Now on the lowest level - need to set everything for this one...
        #
        
        #
        # Saveframe only - saveframe object has to be created in this loop
        #

        printSFFlag = False

        if (ccpnLoopInfo[0] == '' or self.ccpnVar[ccpnLoopInfo[0] ] is not None):
          printSFFlag = True

        if nmrStarElementDict.has_key('name'):
          if nmrStarElementDict['name'] == 'Chem_shift_reference':

            #keywds = {'application': self.format,
            #          'keyword':     'chemShiftRefs',
            #          'value':       'empty'}

            #if self.ccpnVar['nmrEntry'].findFirstApplicationData(**keywds):
            #  printSFFlag = False

            printSFFlag = False

            for expt in self.ccpnVar['nmrEntry'].experiments:
              if expt.shiftReferences:
                printSFFlag = True

          elif nmrStarElementDict['name'] == 'NMR_spectrometer_list':
            printSFFlag = False

            for expt in self.ccpnVar['nmrEntry'].experiments:
              if expt.spectrometer:
                printSFFlag = True
                break

          elif nmrStarElementDict['name'] == 'Experiment_list':
            if not self.ccpnVar['nmrEntry'].experiments:
              printSFFlag = False

          elif nmrStarElementDict['name'] == 'Entity_natural_src_list':
            printSFFlag = False

            for entryMol in self.ccpnVar['nmrEntry'].entryMolecules:
              if entryMol.molecule.naturalSource:
                printSFFlag = True
                break

          elif nmrStarElementDict['name'] == 'Entity_experimental_src_list':
            printSFFlag = False

            for entryMol in self.ccpnVar['nmrEntry'].entryMolecules:
              if entryMol.experimentalSource:
                printSFFlag = True
                break

          elif nmrStarElementDict['name'] == 'Conformer_family_coord_set':
            if not self.ccpnVar[ccpnLoopInfo[0] ].structureEnsemble:
              printSFFlag = False

          elif nmrStarElementDict['name'] == 'Resonance_linker_list':
            printSFFlag = False

            for nmrLinks in ('derivedDataLists','experiments','measurementLists'):
      
              nmrLinkedObjects = getattr(self.ccpnVar['nmrEntry'], nmrLinks)
      
              for nmrLinkedObject in nmrLinkedObjects:
                if nmrLinkedObject:
                  if nmrLinkedObject.className == 'NmrProject':
                    if nmrLinkedObject.resonances:
                      printSFFlag = True
                    elif nmrLinkedObject.fixedResonances:
                      printSFFlag = True

        #print 'FLAG: [%s] [%s]' % (printSFFlag, nmrStarElementDict)

        if level == 'saveFrame':

          if not self.keyFlag and printSFFlag:

            self.currentSaveFrameTitle = saveFrameTitles[ccpnObject]

            if not self.currentSaveFrameTitle:
              self.currentSaveFrameTitle = 'void'

            self.currentSaveFrameTitle.replace(' ', '_')

            starElement = self.nmrStarFile.setupSaveFrame(self.sfCatName,self.currentSaveFrameTitle)

          else:
            starElement = None

        #
        # Relevant for both saveframe and table...
        #

        if printSFFlag:
          if self.keyFlag:
            self.setStarKeys(ccpnLoopInfo,parentCcpnLoopInfo,writeStarElementDict,nmrStarElementDict,starElementName,previousLowestCcpnLoopInfo)

          else:
            (tagNames,writeTags,conditionalTags) = self.setupStarTagInfo(writeStarElementDict,nmrStarElementDict)

            presetValues = self.setStarTagPresets(starElement,starElementName,tagNames,writeTags,conditionalTags,ccpnLoopInfo,parentCcpnLoopInfo,writeStarElementDict,nmrStarElementDict)

            if starElement:
              self.setStarTags(starElement,starElementName,tagNames,writeTags,conditionalTags,presetValues,ccpnLoopInfo,parentCcpnLoopInfo,writeStarElementDict,nmrStarElementDict)

        #
        # Saveframe only - sets the tables, handles ignoring foreign links correctly (should not happen in loops inside the saveframe itself)...
        #

        if level == 'saveFrame' and writeStarElementDict.has_key('tables'):
        
          #
          # Make sure foreign tags belonging to this SF are *not* ignored - necessary
          #
          
          if not self.keyFlag:
            self.ignoreForeignLinks = self.ccpn2NmrStar.ignoreForeignLinks[:]
            if starElementName in self.ignoreForeignLinkStarElements.keys():
              for ignoreForeignLink in self.ignoreForeignLinkStarElements[starElementName]:
                self.ignoreForeignLinks.pop(self.ignoreForeignLinks.index(ignoreForeignLink) )

            self.ignorePrimaryKeys = self.ccpn2NmrStar.ignorePrimaryKeys[:]
            if starElementName in self.ignorePrimaryKeyStarElements.keys():
              for ignorePrimaryKey in self.ignorePrimaryKeyStarElements[starElementName]:
                self.ignorePrimaryKeys.pop(self.ignorePrimaryKeys.index(ignorePrimaryKey) )

          #
          # Handle tables
          #

          for tableName in nmrStarElementDict['tableNames']:

            if tableName not in writeStarElementDict['tables'].keys():
              continue

            self.writeStarTableDict = writeStarElementDict['tables'][tableName]
            self.nmrStarTableDict = nmrStarElementDict['tables'][tableName]

            if self.writeStarTableDict.has_key('ccpnMap'):
              if not self.setTableData(tableName,saveFrame = starElement):
                continue

            elif self.keyFlag:
              self.setMessage(" No ccpn mapping for %s" % tableName)  

        #
        # Track this level for setting global IDs
        #
        
        previousLowestCcpnLoopInfo = ccpnLoopInfo

        # CJP - set Custom AppData stuff

        # If no value for writeStarElementDict['ccpnMap'], we are in the entry saveframe.

        if not self.keyFlag:
          if level == 'saveFrame':

            if not writeStarElementDict['ccpnMap']:
              self.setExtraEntryLevelTags(self.nmrEntry)
            else:
              self.setExtraAppDataTags(ccpnObject, self.sfCatName, starElement)

    return previousLowestCcpnLoopInfo


  ########################
  # Set up the star keys #
  ########################
  
  def setStarKeys(self,ccpnLoopInfo,parentCcpnLoopInfo,writeStarElementDict,nmrStarElementDict,starElementName,previousLowestCcpnLoopInfo):

    fullItemKey = ccpnLoopInfo[1]

    #
    # Note: no problem if same ID occurs twice for the same loop (in different saveframes),
    # since search on starKeys is done by CCPN object, not ID.
    #

    keyNames = nmrStarElementDict['sourcePrimaryKeys']

    if writeStarElementDict.has_key('tags'):
      writeTags = writeStarElementDict['tags']
    else:
      writeTags = {}

    for keyName in keyNames:
      
      #
      # Get the information from the nmrStar dictionary
      #
      
      (default,returnFunc,foreignKey,obligatory) = nmrStarElementDict['tags'][keyName]

      value = None

      #
      # Set matchKeyName (for matching items with a composite key)
      #
      
      currentCcpnLoopInfo = None
      
      if not fullItemKey:
      
        # No nested loop info
      
        matchKeyName = (starElementName,keyName,)
        currentCcpnLoopInfo = ccpnLoopInfo

      elif fullItemKey[-1] == keyName:
        
        # Are at the bottom (ccpnLoopInfo) level
        
        matchKeyName = (starElementName,) + fullItemKey
        currentCcpnLoopInfo = ccpnLoopInfo

        #print 'MAT KEY: [' + str(matchKeyName) + ']'
      
      else:
      
        # Determine the parent level we're at...
      
        matchKeyName = (starElementName,)
        
        for i in range(0,len(fullItemKey) ):
          
          matchKeyName += (fullItemKey[i],)
        
          if fullItemKey[i] == keyName:
            
            currentCcpnLoopInfo = parentCcpnLoopInfo[i]
            break
      
      #
      # Keep track of mapping to relevant star tag...
      #
      
      self.starTagToStarKey[starElementName + self.nmrStarFile.tagSep + keyName] = matchKeyName
      
      #
      # Set variables...
      #
      
      if not currentCcpnLoopInfo:
        self.setMessage("  Warning: %s loop has no information for tag %s." % (starElementName,keyName) )
        continue
 
      ccpnVarKey = currentCcpnLoopInfo[0]
      
      #
      # Set the current ID. If the nmrStarKeyType is set, the ID will be continued from the previous 'lowest'
      # level loop. This is, for example, to set the Atom ID in Atom_site
      #
            
      if writeStarElementDict.has_key('nmrStarKeyType'):
        if previousLowestCcpnLoopInfo and currentCcpnLoopInfo == ccpnLoopInfo and currentCcpnLoopInfo != previousLowestCcpnLoopInfo:
          currentID = previousLowestCcpnLoopInfo[2]
        else:
          currentID = currentCcpnLoopInfo[2]
      else:
        currentID = currentCcpnLoopInfo[2]
            
      #
      # Handle source keys...
      #

      if not foreignKey:

        ccpnObject = self.ccpnVar[ccpnVarKey]

        #
        # Set up starKeys if new matchKeyName
        #

        if not self.starKeys.has_key(matchKeyName):
          self.starKeys[matchKeyName] = ({},ccpnObject) # Tracking as reference object

        #
        # Check if there's a default value available
        #
        
        value = None

        if writeTags.has_key(keyName):

          ccpnTagMapping = writeTags[keyName]       
          value = self.getCcpnMapValue(ccpnTagMapping)

        #
        # If not, then start looking for a match... if no match then set everything.
        #

        (matchFound,matchValue) = self.findStarKeyMatch(matchKeyName,ccpnObject)
       
        if not value:

          #
          # If the ccpnObject present in starKeys, DO NOT increase counter!
          #

          if not matchFound:
            currentCcpnLoopInfo[2] = currentID + 1  # This is the currentID! But have to set directly...
            
          value = currentCcpnLoopInfo[2]    

        #
        # Finally set all the key and ID information...
        #
        
        if not matchFound:
          #print 'OBJ: [%s] [%s] [%s]' % (ccpnObject, value, matchKeyName)
          if fullItemKey and fullItemKey[0] != keyName:
            parentCcpnVarKey = parentCcpnLoopInfo[-1][0]
            parentDictInfo = {parentCcpnVarKey: self.ccpnVar[parentCcpnVarKey]}
          else:
            parentDictInfo = None

          self.starKeys[matchKeyName][0][ccpnObject] = (value,parentDictInfo)
          self.starKeyToCcpnVar[matchKeyName] = ccpnVarKey

        else:
          pass
          #print 'OBJ FOUND: [%s] [%s] [%s]' % (ccpnObject, value, matchKeyName)

        #print 'VAR: [%s] [%s] [%s]' % (ccpnVarKey, matchKeyName, value)    

    #
    # For custom code, sometimes have to map the relevant objects - is done here
    #

    #
    # Atoms are a special case... mapping all resonances at this stage
    #
    # NOTE/TODO: currently only useful for 'normal' resonances. When handled by individual atom other mechanism is in place!!!
    #
    
    for setAtomsTag in self.setAtomsTags:

      if writeTags.has_key(setAtomsTag):

        #
        # This creates a list of all the resonances for linking. Note that this is
        # the only place where they get collected if they are fixedResonances. The
        # link to atoms is only checked after the first stage (so after the NMR-STAR
        # keys are set.
        #

        if setAtomsTag.count('individual'):
          # Can only handle single cases!!

          objectMap = writeTags[setAtomsTag]
          measurementByIndividualAtom = self.getCcpnObject(objectMap)
          resonances = measurementByIndividualAtom.resonance
          
        else:
          ccpnMap = writeTags[setAtomsTag]       
          resonances = self.getCcpnObject(ccpnMap)

        if type(resonances) not in [type([]),type( frozenset() ), type(tuple())]:
          resonances = [resonances]

        for resonance in resonances:
          if resonance not in self.ccpnVar['resonances']:
            self.ccpnVar['resonances'].append(resonance)
                        
    #
    # Nodes for distance constraints are a special case...
    #

    if writeTags.has_key('CUSTOM_setNodes'):
      # TODO: also check code in setStarTagPresets!!


      # This is the counter for the constraint item - first one is a special case!
      if ccpnLoopInfo[2] == 1:

        #
        # Have to handle first node differently depending on whether there's multiple items...
        #

        constraintItem = self.ccpnVar['constraintItem'] 
        constraint = ccpnObject.parent
        
        # TODO: using matchKeyName here because should be last one in list?!?!

        if len(constraint.items) > 1:

          ccpnLoopInfo[2] += 1
          parentDictInfo = self.starKeys[matchKeyName][0][constraintItem][-1]
          self.starKeys[matchKeyName][0][constraintItem] = (ccpnLoopInfo[2],parentDictInfo)
          #print "changed %s value %d" % (matchKeyName,ccpnLoopInfo[2])


  def setupCurrentIDs(self,ccpnMap,ccpn2Star):
    
    #
    # Creates a fresh currentID set 
    #
    
    ccpnLoopInfo = []

    if ccpn2Star.has_key('ccpnLoop'):

      ccpnLoop = ccpn2Star['ccpnLoop']

      if type(ccpnLoop) == type([]):
      
        additionalMappings = {}
        if ccpn2Star.has_key('addCcpnMap'):
          additionalMappings = ccpn2Star['addCcpnMap']

        self.getCcpnNestedLoopValues(ccpnLoop[:],ccpnMap[:],ccpnLoopInfo,additionalMappings)

      else:
        ccpnVarKey = ccpnMap.split('.')[-1]
        ccpnLoopObjects = self.getCcpnMapValue(ccpnLoop)

        ccpnLoopInfo = [ccpnVarKey,None,0,[] ]
        
        if not ccpnLoopObjects:
          self.setMessage("No loop objects found for ccpnMap '%s', ccpnLoop '%s'!" % (ccpnMap,ccpnLoop) )
        
        else:
          for ccpnLoopObject in ccpnLoopObjects:
            ccpnLoopInfo[-1].append( (ccpnLoopObject,None) )

    else:      
      ccpnVarKey = ccpnMap.split('.')[-1]
      ccpnLoopInfo = [ccpnVarKey,None,0,[(self.getCcpnObject(ccpnMap),None)] ]
      
    return ccpnLoopInfo  


  def getCcpnNestedLoopValues(self,ccpnLoop,ccpnMap,ccpnLoopInfo,additionalMappings,nmrStarFullKey = None):
  
    #
    # Go through nested loop to get a list of all objects
    #
    # WARNING: THIS IS A RECURSIVE FUNCTION!
    #
    
    ccpnLoopObjects = self.getCcpnMapValue(ccpnLoop[0])

    if not ccpnLoopObjects:
      ccpnLoopObjects = []
    
    ccpnVarKey = ccpnMap[0][0]

    if not nmrStarFullKey:
      nmrStarFullKey = (ccpnMap[0][1],)
    else:
      nmrStarFullKey += (ccpnMap[0][1],)
         
    ccpnLoopInfo.extend([ccpnVarKey,nmrStarFullKey,0,[] ])

    for ccpnLoopObject in ccpnLoopObjects:

      if len(ccpnLoop) > 1:
        ccpnLowerLoopInfo = []
      else:
        ccpnLowerLoopInfo = None
        
      ccpnLoopInfo[-1].append( (ccpnLoopObject,ccpnLowerLoopInfo) )
        
      self.ccpnVar[ccpnVarKey] = ccpnLoopObject
      
      # Setting additional mappings - important for giant saveframes like coordinates to speed things up
      if additionalMappings.has_key(ccpnVarKey):
        for (addCcpnVarKey,addCcpnMapping) in additionalMappings[ccpnVarKey]:
          self.ccpnVar[addCcpnVarKey] = self.getCcpnMapValue(addCcpnMapping)
      
      if len(ccpnLoop) > 1:
        self.getCcpnNestedLoopValues(ccpnLoop[1:],ccpnMap[1:],ccpnLowerLoopInfo,additionalMappings,nmrStarFullKey = nmrStarFullKey[:])

  
  #########################
  # Set the star tag info #
  #########################

  def setStarTags(self,starElement,starElementName,tagNames,writeTags,conditionalTags,presetValues,ccpnLoopInfo,parentCcpnLoopInfo,writeStarElementDict,nmrStarElementDict):
    
    #
    # Loop over the tag names
    #

    for tagName in tagNames:

      #if tagName not in nmrStarElementDict['tags']:
      #  print 'TAGS: [%s]' % nmrStarElementDict['tags']

      (default,returnFunc,foreignTag,obligatory) = nmrStarElementDict['tags'][tagName]
      starTagName = starElementName + self.nmrStarFile.tagSep + tagName

      value = None
      setTag = True

      if tagName in presetValues:

        value = presetValues[tagName]
      
      elif foreignTag and foreignTag not in self.ignoreForeignLinks:

        # Only ignore the foreign key if the keyFlag 

        value = self.getForeignValue(foreignTag, tagName, writeTags, starElementName, starTagName)

        #print 'FOR: [%s] [%s] [%s] [%s] [%s] [%s]' % (value, foreignTag, tagName, writeTags, starElementName, starTagName)

      elif tagName in nmrStarElementDict['sourcePrimaryKeys'] and starTagName not in self.ignorePrimaryKeys:

        #
        # Should already be set...
        #

        matchKeyName = self.starTagToStarKey[starTagName]
        ccpnObject = self.ccpnVar[self.starKeyToCcpnVar[matchKeyName] ]
        
        (matchFound,value) = self.findStarKeyMatch(matchKeyName,ccpnObject)

        if not value:
          self.setMessage("  Error: in '%s', key %s has not been set." % (starElementName,starTagName) )

      elif default and not (writeTags and writeTags.has_key(tagName) ):

        value = default

      elif writeTags.has_key(tagName):

        ccpnTagMapping = writeTags[tagName]
        value = self.getCcpnMapValue(ccpnTagMapping)

      elif tagName in conditionalTags:

        for condition in writeTags['CONDITIONAL'].keys():

          conditionValue = self.getCcpnMapValue(condition)

          for checkConditionValue in writeTags['CONDITIONAL'][condition]:
            if conditionValue == checkConditionValue:
              conditionalWriteTags = writeTags['CONDITIONAL'][condition][conditionValue]

              if conditionalWriteTags.has_key(tagName):
                ccpnTagMapping = conditionalWriteTags[tagName]       
                value = self.getCcpnMapValue(ccpnTagMapping)
                break
                
          #
          # Handle 'DEFAULT' tag...
          #
          
          if value == None and writeTags['CONDITIONAL'][condition].has_key('DEFAULT'):
            conditionalWriteTags = writeTags['CONDITIONAL'][condition]['DEFAULT']

            if conditionalWriteTags.has_key(tagName):             
              ccpnTagMapping = conditionalWriteTags[tagName]       
              value = self.getCcpnMapValue(ccpnTagMapping)

      elif obligatory:

        value = '?'

      elif tagName == 'Sf_framecode':

        value = self.currentSaveFrameTitle

      else:

        value = None

        # TODO. If this is set, sometimes get non-rectangular matrices
        # that cannot be transposed.  Do not set it at the moment.
        #setTag = 0

      if setTag:

        #if not self.tempTagInfo.has_key(starElementName):
        #  self.tempTagInfo[starElementName] = {}

        #if not self.tempTagInfo[starElementName].has_key(starTagName):
        #  self.tempTagInfo[starElementName][starTagName] = 0

        #self.tempTagInfo[starElementName][starTagName] += 1

        if value is not None:
          value = returnFunc(value)

          #tmpValue = returnFunc(value)

          #if tmpValue is not None:
          #  value = tmpValue
          #else:
          #  tmpValue = returnFunc(str(value) )

          #  if tmpValue is not None:
          #    value = tmpValue

        starElement.setTag(self.nmrStarFile.tagStart + starTagName, value)

        #print 'SETTING: [%s] [%s] [%s]' % (starElement, starTagName, value)


  def setupStarTagInfo(self,writeStarElementDict,nmrStarElementDict):

    #
    # Order of tags the same as in the original dictionary
    #

    tagNames = nmrStarElementDict['tagNames']

    #
    # Set the 'write' tags for mapping from CCPN 
    #

    if writeStarElementDict.has_key('tags'):
      writeTags = writeStarElementDict['tags']
    else:
      writeTags = {}

    #
    # Keep track of which tags are conditional
    #

    conditionalTags = []

    if writeTags.has_key('CONDITIONAL'):
      for condition in writeTags['CONDITIONAL'].keys():
        for conditionValue in writeTags['CONDITIONAL'][condition]:
          for conditionalTag in writeTags['CONDITIONAL'][condition][conditionValue]:
            if not conditionalTag in conditionalTags:
              conditionalTags.append(conditionalTag)

    return (tagNames,writeTags,conditionalTags)  


  def setStarTagPresets(self,starElement,starElementName,tagNames,writeTags,conditionalTags,ccpnLoopInfo,parentCcpnLoopInfo,writeStarElementDict,nmrStarElementDict):

    #
    # Presetvalues are special cases...
    #

    presetValues = {}

    #
    # This for special functions...
    #

    for writeTagKey in writeTags.keys():
      searchObj = self.patt['getLocalValue'].search(writeTagKey)
      if searchObj:
        starTagValue = searchObj.group(1)
        value = writeTags[writeTagKey]
        
        presetValues[starTagValue] = value
       
    #
    # Run application data setter if required...
    #
    
    if writeTags.has_key('CUSTOM_setApplDataNames'):
      ccpnMap = writeTags['CUSTOM_setApplDataNames']
      self.setApplDataNames(ccpnMap,presetValues,nmrStarElementDict,parentCcpnLoopInfo)

    #
    # Atoms are a special case...
    #

    for setAtomsTag in self.setAtomsTags:

      if writeTags.has_key(setAtomsTag):

        if setAtomsTag.count('individual'):
          # Can only handle single cases!!
          atomHandling = 'individual'
          objectMap = writeTags[setAtomsTag]
          measurementByIndividualAtom = self.getCcpnObject(objectMap)
          ccpnMap = objectMap + '.resonance.void'

          self.setAtomPresetValues(measurementByIndividualAtom,presetValues,'',nmrStarElementDict,atomHandling)
          
        else:
          atomHandling = 'sets'
          ccpnMap = writeTags[setAtomsTag]       
          resonances = self.getCcpnObject(ccpnMap)
                    
          if type(resonances) in [type([]),type( () )]:

            # Make sure to get them ordered, if possible! Only relevant for constraintItems(?)
            # TODO Might need to add measurement stuff in here!!!!!
            if ccpnMap.count('constraintItem'):
              ccpnLoopParent = self.ccpnVar['constraintItem']
              if hasattr(ccpnLoopParent,'orderedResonances') and ccpnLoopParent.orderedResonances:
                resonances = list(ccpnLoopParent.orderedResonances)

            for i in range(0,len(resonances) ):
              starCode = '_%d' % (i+1)
              self.setAtomPresetValues(resonances[i],presetValues,starCode,nmrStarElementDict,atomHandling)

          else:
            self.setAtomPresetValues(resonances,presetValues,'',nmrStarElementDict,atomHandling)

        self.setApplDataNames(ccpnMap,presetValues,nmrStarElementDict,parentCcpnLoopInfo)
        
        break

    #
    # Nodes for distance constraints are a special case...
    #

    if writeTags.has_key('CUSTOM_setNodes'):

      ccpnMap = writeTags['CUSTOM_setNodes']       
      constraintItem = self.ccpnVar['constraintItem']
      constraint = constraintItem.parent
      constraintItems = constraint.sortedItems()

      if len(constraintItems) > 1 and constraintItem == constraintItems[0]:

        #
        # Have to handle first node differently depending on whether there's multiple items...
        #

        presetValues['Node_ID'] = 1
        presetValues['Down_node_ID'] = 2
        presetValues['Right_node_ID'] = None
        presetValues['Logic_operation'] = 'OR'

        self.setStarTags(starElement,starElementName,tagNames,writeTags,conditionalTags,presetValues,ccpnLoopInfo,parentCcpnLoopInfo,writeStarElementDict,nmrStarElementDict)

        del(presetValues['Node_ID'])
          
      presetValues['Down_node_ID'] = None
      presetValues['Logic_operation'] = None

      if constraintItem == constraintItems[-1]:
        presetValues['Right_node_ID'] = None
      else:
        presetValues['Right_node_ID'] = constraintItems.index(constraintItem) + 3
        
    return presetValues


  def setAtomPresetValues(self,infoObject,presetValues,starCode,nmrStarElementDict,atomHandling):
  
    if atomHandling == 'individual':
      nmrStarFormat = self.ccpnVar['nmrStarFormatIndividual']
      measurementByIndividualAtom = infoObject
      resonanceToAtom = measurementByIndividualAtom.resonanceToAtom
      resonance = resonanceToAtom.resonance
      
    else:
      nmrStarFormat = self.ccpnVar['nmrStarFormatSet']
      resonance = infoObject      
      
      if not nmrStarFormat.resonanceToAtoms.has_key(resonance):

        #
        # TODO: STILL HAVE TO WRITE OUT AUTHOR CODE IF AVAILABLE!!!
        #

        for tagName in ['Entity_assembly_ID','Entity_ID','Seq_ID','Comp_index_ID','Comp_ID','Atom_type']:

          self.setSimplePresetValue(presetValues,tagName + starCode,None)

        self.setResonanceAuthorPdbInfo(resonance,None,None,"")
        
        return

      #
      # Specific code to handle resonances...
      #

      resonanceToAtoms = nmrStarFormat.resonanceToAtoms[resonance]

      # TODO: how can I handle this?! Loop all combinations?!
      if len(resonanceToAtoms) > 1:

        #
        # TODO: In this case the actual resonance should be used (if resonances do enter star!)
        #

        self.setMessage("  Error: multiple mappings for resonance - using first one.")

      resonanceToAtom = resonanceToAtoms[0]

    #TODO: if unrecognized chemComp?? 'Entry_atom_ID_1': [None,returnStarString,'Atom.Entry_atom_ID',None],
    
    #
    # TODO: use 'Assembly_atom_ID' if available? Direct link to atom inside NMR-STAR, but need to 
    # have this info for all atoms in molecular system...
    #

    #
    # Tag Entity_assembly_ID
    #

    tagName = 'Entity_assembly_ID' + starCode
    chain = resonanceToAtom.chain

    (default,returnFunc,foreignKey,obligatory) = nmrStarElementDict['tags'][tagName]
    
    (matchFound,matchValue) = self.findStarKeyMatch(self.starTagToStarKey[foreignKey],chain)
    presetValues[tagName] = matchValue

    #
    # Tag Entity_ID
    #

    tagName = 'Entity_ID' + starCode
    molecule = chain.molecule

    (default,returnFunc,foreignKey,obligatory) = nmrStarElementDict['tags'][tagName]
    (matchFound,matchValue) = self.findStarKeyMatch(self.starTagToStarKey[foreignKey],molecule)
    presetValues[tagName] = matchValue

    #
    # Tag Seq_ID
    #
    
    seqId = resonanceToAtom.seqId
    self.setSimplePresetValue(presetValues,'Seq_ID' + starCode,seqId)

    #
    # Tag Comp_index_ID
    #

    tagName = 'Comp_index_ID' + starCode
    residue = chain.findFirstResidue(seqId = seqId)
    molResidue = self.ccpn2NmrStar.getActualMolResidue(residue)

    (default,returnFunc,foreignKey,obligatory) = nmrStarElementDict['tags'][tagName]
    (matchFound,matchValue) = self.findStarKeyMatch(self.starTagToStarKey[foreignKey],molResidue)
    presetValues[tagName] = matchValue

    #
    # Tag Comp_ID
    #

    self.setSimplePresetValue(presetValues,'Comp_ID' + starCode,findChemCompSysName('CIF',molResidue.chemCompVar.chemComp), default = molResidue.ccpCode)

    #
    # Tag Atom_ID
    #
    
    atomName = resonanceToAtom.atomName
    atomType = resonanceToAtom.atomType
    
    # If direct resonance-atom connection, use this info (if resonanceToAtom didn't get it)
    # TODO: if this is something that often occurs, probably need to either do something on the
    # resonanceToAtom level, or make code below more complex.
    if not atomName:
      resSet = resonanceToAtom.resonance.resonanceSet
      if resSet:
        atomSets = resSet.sortedAtomSets()
        if len(atomSets) == 1:
          atoms = atomSets[0].sortedAtoms()
          if atoms:
            if len(atoms) == 1:
              atomName = self.ccpn2NmrStar.getPdbAtomName(atoms[0])
            atomType = atoms[0].chemAtom.elementSymbol
    
      if not atomName:
        print "  Error: no single atom name found for resonance %d" % (resonanceToAtom.resonance.serial)
      if not atomType:
        print "  Error: no single atom element type found for resonance %d" % (resonanceToAtom.resonance.serial)

    self.setSimplePresetValue(presetValues,'Atom_ID' + starCode,atomName)

    #
    # Tag Atom_type
    #

    self.setSimplePresetValue(presetValues,'Atom_type' + starCode,atomType)
    
    #
    # Here pre-set information for author (nmrStar derived) and PDB identifiers, so only have to do this once!
    #
    
    self.setResonanceAuthorPdbInfo(resonance,chain,residue,atomName)
    
    """
    BELOW NOT NECESSARY ANY MORE!
    
    #
    # Check if forcing author tags... necessary for shift files that contain original PDB info (Wim 16/10/09)
    #
       
    if resonance.findFirstApplicationData(application = self.format, keyword = "forceAtomAuthorWrite", value=True):
      
      chainCode = chain.pdbOneLetterCode

      self.setSimplePresetValue(presetValues,'Auth_asym_ID' + starCode,chainCode)
      self.setSimplePresetValue(presetValues,'Auth_entity_assembly_ID' + starCode,chainCode)

      if residue.seqInsertCode != " ":
        seqCode = "%d%s" % (residue.seqCode,residue.seqInsertCode)
      else:
        seqCode = residue.seqCode

      self.setSimplePresetValue(presetValues,'Auth_seq_ID' + starCode,seqCode)

      if atomName:
          self.setSimplePresetValue(presetValues,'Auth_atom_ID' + starCode,atomName)
      
      self.setSimplePresetValue(presetValues,'Auth_comp_ID' + starCode,findChemCompSysName('CIF',molResidue.chemCompVar.chemComp), default = molResidue.ccpCode)
    """
    
  def setResonanceAuthorPdbInfo(self,resonance,chain,residue,atomName):
    
    if not self.resonanceAuthorPdbInfo.has_key(resonance):
      
      authorResonanceNames = [appData.value for appData in resonance.findAllApplicationData(application = self.format, keyword = assign_kw)] 

      #
      # Note: this will only work if a single structure ensemble is available... custom for wwPDB stuff.
      #
      # Could speed up with dictionary, but probably not that much difference!
      #
      
      pdbResonanceInfo = None
      
      if not self.coordinateInfoMap.has_key(chain):
        pdbChainCode = cChain = None
        if self.ccpnVar.has_key('structureEnsemble'):
          structureEnsemble = self.ccpnVar['structureEnsemble']
        
          if structureEnsemble:
            # Or by code?
            cChain = structureEnsemble.findFirstCoordChain(chain = chain)
            
            if cChain:
              origChainCodeApplData = cChain.findFirstApplicationData(application = 'pdb', keyword = 'originalChainCode')
            
              if origChainCodeApplData:
                pdbChainCode = origChainCodeApplData.value

        self.coordinateInfoMap[chain] = (pdbChainCode,cChain)

      else:
        (pdbChainCode,cChain) = self.coordinateInfoMap[chain]
        
      # Only do something if coordinates exist
      if pdbChainCode:
        if not self.coordinateInfoMap.has_key(residue):
          pdbResidueInfo = None
        
          # TODO or by seqId?
          cResidue = cChain.findFirstResidue(residue = residue)
          
          if cResidue:
            origSeqCodeApplData = cResidue.findFirstApplicationData(application = 'pdb', keyword = 'originalSeqCode')
            origSeqInsertCodeApplData = cResidue.findFirstApplicationData(application = 'pdb', keyword = 'originalSeqInsertCode')
            origLabelApplData = cResidue.findFirstApplicationData(application = 'pdb', keyword = 'originalResLabel')
          
            if origSeqCodeApplData:
              if origSeqInsertCodeApplData:
                pdbSeqInsertCode = origSeqInsertCodeApplData.value
              else:
                pdbSeqInsertCode = None
                
              if origLabelApplData:
                pdbResLabel = origLabelApplData.value
              else:
                pdbResLabel = None
                            
              pdbResidueInfo = (origSeqCodeApplData.value,pdbSeqInsertCode,pdbResLabel)
          
          self.coordinateInfoMap[residue] = (pdbResidueInfo,cResidue)
        
        else:
          (pdbResidueInfo,cResidue) = self.coordinateInfoMap[residue]
          
        if pdbResidueInfo:
          #origAtomName = coordAtom.findFirstApplicationData(application = formatName, keyword = 'originalName')
         
          pdbResonanceInfo = (pdbChainCode,pdbResidueInfo,atomName)
  
      #
      # Set the dictionary...
      #
    
      self.resonanceAuthorPdbInfo[resonance] = (authorResonanceNames,pdbResonanceInfo)


  #######################
  # Low level functions #
  #######################

  def setSimplePresetValue(self,presetValues,tagName,value,default = None,emptyOverwrite = True):
    
    if not value:
      value = default

    if emptyOverwrite or value or not presetValues.has_key(tagName) or not presetValues[tagName]:
      presetValues[tagName] = value

  def searchDataTitle(self,ccpnObject,nmrStarElementDict):
  
    title = None
    
    #
    # First find out if obligatory name - use that if so...
    #

    if nmrStarElementDict.has_key("saveFrameCode"):
      title = nmrStarElementDict["saveFrameCode"]

    #
    # Check whether an applicationData title exists...
    #

    elif ccpnObject and hasattr(ccpnObject, 'findFirstApplicationData'):
        
      if self.ccpn2NmrStar.mapAppDataSaveFrames.has_key(self.sfCatName):
        useSfName = self.ccpn2NmrStar.mapAppDataSaveFrames[self.sfCatName][0]
      else:
        useSfName = self.sfCatName

      applData = ccpnObject.findFirstApplicationData(application = self.format, keyword = useSfName)

      if applData:

        sfKey = applData.value
        
        keyword = useSfName +  '_title_' + sfKey
        titleAppData = ccpnObject.findFirstApplicationData(application = self.format, keyword = keyword)

        if titleAppData:
          title = titleAppData.value
          
    return title


  def setApplDataNames(self,ccpnMap,presetValues,nmrStarElementDict,parentCcpnLoopInfo):
  
    #
    # Check if anything special has to be done, or if just direct match to object
    #
    
    if ccpnMap.count('.'):
      ccpnObjectMap = '.'.join(ccpnMap.split('.')[:-1])
    else:
      ccpnObjectMap = ccpnMap
      
    ccpnObject = self.getCcpnObject(ccpnObjectMap)

    if not ccpnObject:
      return

    ccpnObjectClassName = ccpnObject.className
    usingParentObject = False
    
    #
    # This for handling one single resonance coming off a loop parent .resonances link
    #
    
    if ccpnObjectClassName == 'FixedResonance' and parentCcpnLoopInfo:
      ccpnObject = self.ccpnVar[parentCcpnLoopInfo[-1][0] ]      
      usingParentObject = True
        
    if ccpnObject and type(ccpnObject) != type(''):
      
      #
      # Resonances and fixedResonances original name information
      #
      
      if ccpnObjectClassName == 'Resonance' or hasattr(ccpnObject,'resonances'):

        # This only matches resonance info if it was ordered. For PairwiseItem stuff this is done automatically now (November 2007 - Wim)
        origAssignApplData = ccpnObject.findAllApplicationData(application = self.format, keyword = 'origAssign')
        origResLabelApplData = ccpnObject.findAllApplicationData(application = self.format, keyword = 'origResLabel')

        origDataItems = len(origAssignApplData)

        #
        # Below picks the right application data if the original residue names are not stored on the resonance itself,
        # but on the object linking to them (e.g. DistanceConstraintItem), but only 1 line is written to the NMR-STAR,
        # in the case of DistanceConstraintItem in effect a constraint item member.
        #

        if usingParentObject and origAssignApplData and origDataItems == len(ccpnObject.resonances):

          #
          # Reorganize - only select correct resonance in case of single resonance coming from .resonances link
          #
          
          resonance = self.ccpnVar['resonance']
          
          resonanceList = []
          if hasattr(ccpnObject,'orderedResonances'):
            resonanceList = list(ccpnObject.orderedResonances)
          
          if not resonanceList and hasattr(ccpnObject,'sortedResonances'):
            resonanceList = ccpnObject.sortedResonances()
            
          if not resonanceList:
            resonanceList = list(ccpnObject.resonances)

          # Warning: this will only work if this list is ordered!! See code above
          resonanceIndex = resonanceList.index(resonance)
          
          (resonanceNames,pdbResonanceInfo) = self.resonanceAuthorPdbInfo[resonance]         
          # TODO To speed up, make this dictionary so don't always have to re-search values
          #resonanceNames = [appData.value for appData in resonance.findAllApplicationData(application = self.format, keyword = assign_kw)] # TODO [remove this comment] per resonance; self.format - always nmrstar - look for res names

          # TODO: Should we put the original names as part of the measurements,
          # then the [whole] code needs to be written to handle this case;
          # Do we need the exact atom names for the imported data.

          origAssignApplDataMatch = None
          for resonanceName in resonanceNames:
            for j in range(len(origAssignApplData) ):
              if origAssignApplData[j].value == resonanceName:
                origAssignApplDataMatch = origAssignApplData.pop(j)
                #print 'RES NAME: [' + str(resonanceName) + '] [' + str(origAssignApplDataMatch) + ']'
                break

          origAssignApplData = [origAssignApplDataMatch]

          #print 'Res Ind: [' + str(ccpnMap) + ']; Ind: [' + str(resonanceIndex) + ']; Len: [' + str(len(origResLabelApplData) ) + '] [' + str(origResLabelApplData) + '] [' + str(origDataItems) + ']'

          if origResLabelApplData is not None and resonanceIndex < len(origResLabelApplData):
            origResLabelApplData = [origResLabelApplData[resonanceIndex] ]
          else:
            origResLabelApplData = None

          origDataItems = len(origAssignApplData)
       
        #
        # Special case for dihedral angles where only the angle and residue are given - need to propagate this info...
        #
        
        if not origAssignApplData and ccpnObjectClassName == 'DihedralConstraint':
          origSeqCodeApplData = ccpnObject.findFirstApplicationData(application = self.format, keyword = 'origSeqCode')
          origChainCodeApplData = ccpnObject.findFirstApplicationData(application = self.format, keyword = 'origChainCode')
          origDataItems = 4
        else:
          origSeqCodeApplData = None
        
        #
        # Make sure cleanup done - can't have None in there.
        #
        
        if origAssignApplData.count(None):
          origAssignApplData = None
        
        #
        # Now find out which values to use, whether available or not
        #  
          
        origAssignValueList = []
        pdbAssignValueList = []
        
        if origAssignApplData or origSeqCodeApplData:

          #print 'ORIG: [' + str(origAssignApplData) + ']'
        
          if ccpnObjectClassName == 'Resonance':
            resonanceList = [ccpnObject]                 # Direct resonance match
          elif origDataItems == 1:
            resonanceList = [self.ccpnVar['resonance']]  # Distance constraint item case (old saveframe)
          elif origDataItems == 2:
            resonanceList = ccpnObject.sortedResonances() # Distance constraint item case (new saveframe)
          else:
            resonanceList = list(ccpnObject.resonances)  # Dihedral constraint
            
          for i in range(len(resonanceList)):
          
            resonance = resonanceList[i]

            (resonanceNames,pdbResonanceInfo) = self.resonanceAuthorPdbInfo[resonance]         
                        
            # Because application data is now in no particular order, need to find specific match...
            origResLabelApplDataMatch = None
            if origSeqCodeApplData:
              # This is a special case for dihedrals with angle name and only residue information
              if origChainCodeApplData:
                origChainCode = origChainCodeApplData.value
              else:
                origChainCode = defaultMolCode                

              origAssignApplDataValue = "%s.%d.XXX" % (origChainCode,origSeqCodeApplData.value)
              #print origAssignApplDataValue

              if origResLabelApplData != None:
                origResLabelApplDataMatch = origResLabelApplData[0]
                
            elif len(origAssignApplData) == 1 and origAssignApplData[0]:
              origAssignApplDataValue = origAssignApplData[0].value
              if origResLabelApplData:
                origResLabelApplDataMatch = origResLabelApplData[0]
                
            else:
              origAssignApplDataValue = None
              
              # These are the original resonance names, even when data was merged, so one of these has to match the
              # resonance name that's listed as appData on the parent item (e.g. a DistanceConstraintItem)
              resonanceNames = [appData.value for appData in resonance.findAllApplicationData(application = self.format, keyword = assign_kw)] 
              #print resonanceNames
              for resonanceName in resonanceNames:
                for j in range(len(origAssignApplData) ):
                  if origAssignApplData[j].value == resonanceName:
                    origAssignApplDataValue = origAssignApplData.pop(j).value
                    if origResLabelApplData != None and j < len(origResLabelApplData):
                      origResLabelApplDataMatch = origResLabelApplData.pop(j)
                    
                    #print 'RES NAME: [' + str(resonanceName) + '] [' + str(origAssignApplDataValue) + ']'
                    break
                
                # Added Wim 19/10/2009
                if origAssignApplDataValue:
                  break
                
            if origAssignApplDataValue:
              origAssignApplDataSplit = origAssignApplDataValue.split('.')
            
            if not origAssignApplDataValue or len(origAssignApplDataSplit) != 3:
              # Could try to set values 3 and 4 into just 3.  But data types seem wrong for 1l8c
              (chainCode,seqCode,atomName) = (None,None,None)  
            else:
              (chainCode,seqCode,atomName) = origAssignApplDataSplit
              
              # Special case for dihedrals
              if atomName == 'XXX':
                atomName = None

            #print 'Res Ind2: [' + str(ccpnMap) + ']; Ind2: [' + str(i) + ']; Len2: [' + str(len(origResLabelApplData) ) + '] [' + str(origResLabelApplData) + '] [' + str(origDataItems) + ']'

            if origResLabelApplDataMatch:
              resLabel = origResLabelApplDataMatch.value
            else:
              resLabel = None

            if resLabel == 'None':
              resLabel = None

            #print 'CHAIN: [' + str(chainCode) + '] [' + str(seqCode) + '] [' + str(atomName) + '] [' + str(resLabel) + ']'
            #print

            origAssignValueList.append( (chainCode,seqCode,atomName,resLabel) )
            pdbAssignValueList.append(pdbResonanceInfo)
            
        else:
          
          defaultValues = (None,None,None,None)
          
          origDataItems = 1
          if not usingParentObject and hasattr(ccpnObject,'resonances'):
            origDataItems = len(ccpnObject.resonances)
          
          for i in range(origDataItems):
            origAssignValueList.append(defaultValues)   
            pdbAssignValueList.append(None)         

        #
        # And finally set them for NMR-STAR writing...
        #

        for i in range(len(origAssignValueList) ):

          # Only use suffix if not single resonance! Wim 19/10/2009
          if origDataItems > 1:# and ccpnObjectClassName != 'Resonance':
            starCode = '_%d' % (i+1)
          else:
            starCode = ''

          (chainCode,seqCode,atomName,resLabel) = origAssignValueList[i]

          #print 'BEFORE SET: [' + str(origAssignValueList[i]) + '] _' + starCode
          
          # Will now not overwrite if set already (Wim 19/10/2009)

          self.setSimplePresetValue(presetValues,'Auth_asym_ID' + starCode,chainCode,emptyOverwrite = False)
          #self.setSimplePresetValue(presetValues,'Auth_entity_assembly_ID' + starCode,chainCode,emptyOverwrite = False)
          self.setSimplePresetValue(presetValues,'Auth_seq_ID' + starCode,seqCode,emptyOverwrite = False)
          self.setSimplePresetValue(presetValues,'Auth_atom_ID' + starCode,atomName,emptyOverwrite = False)
          self.setSimplePresetValue(presetValues,'Auth_comp_ID' + starCode,resLabel,emptyOverwrite = False)
      
          # Now set PDB stuff, if available...
          if pdbAssignValueList[i]:
            
            pdbResonanceInfo = pdbAssignValueList[i]

            (pdbChainCode,pdbResidueInfo,pdbAtomName) = pdbResonanceInfo
            (pdbSeqCode,pdbSeqInsertCode,pdbResLabel) = pdbResidueInfo
            
            self.setSimplePresetValue(presetValues,'PDB_strand_ID' + starCode,pdbChainCode)
            self.setSimplePresetValue(presetValues,'PDB_residue_no' + starCode,pdbSeqCode)
            self.setSimplePresetValue(presetValues,'PDB_ins_code' + starCode,pdbSeqInsertCode)
            self.setSimplePresetValue(presetValues,'PDB_residue_name' + starCode,pdbResLabel)
            self.setSimplePresetValue(presetValues,'PDB_atom_name' + starCode,pdbAtomName)

        #
        # Dihedral constraint original name information has to be set at this level as well - resonance
        # author code is picking up DihedralConstraint, so need to set this in here.
        #

        if ccpnObjectClassName == 'DihedralConstraint':

          dihedralConstraint = ccpnObject

          origVar = {}

          # TODO: original residue label stuff needs to be set on resonance level...
          for origKeywd in ('origName','origResLabel'):
            origApplData = dihedralConstraint.findFirstApplicationData(application = self.format, keyword = origKeywd)
            if origApplData:
              origVar[origKeywd] = origApplData.value
            else:
              origVar[origKeywd] = None

          self.setSimplePresetValue(presetValues,'Torsion_angle_name',origVar['origName'])
          #self.setSimplePresetValue(presetValues,'Auth_seq_ID' + starCode,orig)


      #
      # Coordinate atom original name information from PDB (or NMR-STAR) files
      #
      # Fully updated for use with coordinate NMR-STAR files and remediated PDB atom naming (Wim 19/02/2008)
      #
      
      elif ccpnObjectClassName == 'Atom':

        namingSystemName = 'PDB_REMED'
                
        coordAtom = ccpnObject
        residue =   self.ccpnVar['coordResidue']
        chain =     self.ccpnVar['coordChain']

        #
        # Try and get information for both PDB and NMRSTAR.
        #
        
        for formatName in ('pdb','nmrStar'):

          if not self.origVar.has_key(formatName):
            self.origVar[formatName] = {'chain': (None,), 'residue': (None,)}
                    
          #
          # Set chain level info (speed issue)
          #
  
          if self.origVar[formatName]['chain'][0] != chain:
            origChainCodeApplData = chain.findFirstApplicationData(application = formatName, keyword = 'originalChainCode')
            if origChainCodeApplData:
              value = origChainCodeApplData.value
            elif origChainCodeApplData and origChainCodeApplData.value == defaultMolCode:
              value = None
            else:
              value = -99999
  
            self.origVar[formatName]['chain'] = (chain,value)
          
          # Continue if no data available!
          if self.origVar[formatName]['chain'][1] == -99999:
            continue  
          
  
          #
          # Set residue level info (speed issue)
          #
  
          if self.origVar[formatName]['residue'][0] != residue:
  
            origSeqCodeApplData = residue.findFirstApplicationData(application = formatName, keyword = 'originalSeqCode')
            origSeqInsertCodeApplData = residue.findFirstApplicationData(application = formatName, keyword = 'originalSeqInsertCode')
  
            origLabelApplData = residue.findFirstApplicationData(application = formatName, keyword = 'originalResLabel')
  
            seqCodeValue = None
            if origSeqCodeApplData:
              seqCodeValue = str(origSeqCodeApplData.value)
            
            seqInsertCodeValue = None
            if origSeqInsertCodeApplData:
              if origSeqInsertCodeApplData.value != defaultSeqInsertCode:
                seqInsertCodeValue = origSeqInsertCodeApplData.value
  
            labelValue = None
            if origLabelApplData:
              labelValue = origLabelApplData.value
  
            chemCompVar = residue.residue.chemCompVar
  
            if not self.chemAtomSysNamesDict.has_key(chemCompVar):
              chemAtomSysNames = findAllSysNamesByChemAtomOrSet(chemCompVar.chemComp,chemCompVar.chemAtoms,namingSystemName)
              self.chemAtomSysNamesDict[chemCompVar] = chemAtomSysNames
            else:
              chemAtomSysNames = self.chemAtomSysNamesDict[chemCompVar]
  
            self.origVar[formatName]['residue'] = (residue,seqCodeValue,seqInsertCodeValue,labelValue,chemAtomSysNames)
  
          #
          # Get atom name info, or find it from the chemAtomSysNames
          #
  
          origAtomName = coordAtom.findFirstApplicationData(application = formatName, keyword = 'originalName')
  
          if not origAtomName:
  
            for origChemAtomSysName in self.origVar[formatName]['residue'][4]:
  
              if origChemAtomSysName.atomName == coordAtom.name:
  
                origAtomName = origChemAtomSysName.sysName
                break
  
            if not origAtomName:
              origAtomName = None
  
          else:
  
            origAtomName = origAtomName.value
          
          starCode = '' # Can use this to add _1, ... tags 
          
          # This is so PDB-specific, will always set these values for this from now on...
          if formatName == 'pdb':
            #PDB_model_num set this as well?
            self.setSimplePresetValue(presetValues,'PDB_strand_ID' + starCode,self.origVar[formatName]['chain'][1])
            self.setSimplePresetValue(presetValues,'PDB_residue_no' + starCode,self.origVar[formatName]['residue'][1])
            self.setSimplePresetValue(presetValues,'PDB_ins_code' + starCode,self.origVar[formatName]['residue'][2])
            self.setSimplePresetValue(presetValues,'PDB_residue_name' + starCode,self.origVar[formatName]['residue'][3])
            self.setSimplePresetValue(presetValues,'PDB_atom_name' + starCode,origAtomName)
            
          elif formatName == 'nmrStar':  
            self.setSimplePresetValue(presetValues,'Auth_asym_ID' + starCode,self.origVar[formatName]['chain'][1])
            #self.setSimplePresetValue(presetValues,'Auth_entity_assembly_ID' + starCode,self.origVar[formatName]['chain'][1])
            
            seqCode = self.origVar[formatName]['residue'][1]
            if self.origVar[formatName]['residue'][2]:
              seqCode += self.origVar[formatName]['residue'][2] 
              
            self.setSimplePresetValue(presetValues,'Auth_seq_ID' + starCode,seqCode)
            self.setSimplePresetValue(presetValues,'Auth_comp_ID' + starCode,self.origVar[formatName]['residue'][3])
            self.setSimplePresetValue(presetValues,'Auth_atom_ID' + starCode,origAtomName)
          
          #TODO ??? 'Auth_atom_ID': None,


  def returnOrigValue(self,value):
  
    #
    # Bit dangerous but unlikely that a variable is present with the star name...
    # is mostly for setting None!
    #

    if value == 'None':
      value = None
      
    else:
      value = value

    return value


  def setMessage(self,message):
    
    #
    # Adds/counts messages
    #
  
    if self.messages.has_key(message):
      self.messages[message] += 1
    else:
      self.messages[message] = 1   


  def findStarKeyMatch(self,matchKeyName,ccpnObject):
    
    #
    # This creates small speedup
    #
    
    (starKeyDict,refCcpnObject) = self.starKeys[matchKeyName]
    
    #
    # Generic function to find a match for a key in self.starKeys
    #
  
    matchFound = 0
    matchValue = None
    
    #if not hasattr(ccpnObject,'className'):
    #  print matchKeyName, ccpnObject
    #  print len(starKeyDict)
     
    # Speedup, instead of has_key
    result = starKeyDict.get(ccpnObject)

    if result:
      (matchValue, parentDictInfo) = result
      matchFound = 1

      if parentDictInfo:
                    
        parentCcpnKey = parentDictInfo.keys()[0]
        parentCcpnObject = parentDictInfo[parentCcpnKey]
        
        if parentCcpnObject != self.ccpnVar[parentCcpnKey]:
          # TODO: What does this mean?
          matchFound = 0
    
    # Try to catch 'mismapping' errors - sped up because reference object now directly available from self.starKeys
    elif starKeyDict and type(refCcpnObject) != type(ccpnObject):
        
      self.setMessage("  Warning: for NMR-STAR key %s, trying to find %s object, but could be %s object!" % (matchKeyName,type(ccpnObject),type(refCcpnObject) ) )
         
    return (matchFound,matchValue)


  def getForeignValue(self, foreignTag, tagName, writeTags = None, starElementName = '', starTagName = ''):

    value = None

    if self.starTagToStarKey.has_key(foreignTag):
        
      foreignKey = self.starTagToStarKey[foreignTag]

      if writeTags is not None and writeTags.has_key(tagName):

        ccpnTagMapping = writeTags[tagName]

        foreignKeyCcpnObject = self.getCcpnMapValue(ccpnTagMapping)

      else:
        # TODO - hack for new foreign key types in Eldon's dictionary
        #if foreignKey[1].count('_'):
        #  patt = re.compile("^(\S+)_(\S+)$")
        #  searchObj = patt.search(foreignKey[1])
        #  foreignKey = (searchObj.group(1), searchObj.group(2) )
        foreignKeyCcpnObject = self.ccpnVar[self.starKeyToCcpnVar[foreignKey] ] # TODO: IS THIS CORRECT?

      (matchFound,matchValue) = self.findStarKeyMatch(foreignKey,foreignKeyCcpnObject)

      if matchFound:
        value = matchValue

      if not value:
        self.setMessage("  Error: in '%s', foreign tag %s does not have CCPN object." % (starElementName,starTagName) )

    elif writeTags is not None and writeTags.has_key(tagName):

      #
      # This is only necessary if there was no previous mapping available!
      #

      ccpnTagMapping = writeTags[tagName]

      value = self.getCcpnMapValue(ccpnTagMapping)

    else:
      self.setMessage("  Error: in '%s', foreign tag %s has not been set." % (starElementName,starTagName) )

    #print 'VALUE: [' + str(foreignTag) + '] [' + str(tagName) + '] [' + str(writeTags) + '] [' + str(starElementName) + '] [' + str(starTagName) + '] [' + str(value) + ']'

    return value


  #
  # Functions for decoding the mapping dictionary
  #
  
  def getAttrOrFunc(self,obj,attrOrFunc):
  
    returnObj = None
  
    if obj:
      
      # Don't bother if object is None (possible)
      
      funcSearch = self.patt['insideBrackets'].search(attrOrFunc)

      if funcSearch:
        funcArgs = funcSearch.group(1)
        funcName = attrOrFunc.replace(funcArgs,'')
        func = getattr(obj,funcName)
        returnObj = eval("func%s" % funcArgs)

      else:
        returnObj = getattr(obj,attrOrFunc)

    return returnObj
  

  def getCcpnObject(self,ccpnMap):

    findObject = None
  
    if ccpnMap:

      ccpnObjectStrings = ccpnMap.split('.')
    
      if not self.ccpnVar.has_key(ccpnObjectStrings[0]):
        raise('  Error: No definition for %s yet in ccpnVar dictionary!' % ccpnObjectStrings[0])

      else:
        startObject = self.ccpnVar[ccpnObjectStrings[0] ]

        #print 'START: [%s] [%s]' % (ccpnMap, startObject)

        if len(ccpnObjectStrings) > 1:
          try:
            findObject = reduce(lambda obj, attrOrFunc: self.getAttrOrFunc(obj,attrOrFunc),ccpnObjectStrings[1:],startObject)
          except:
            print ccpnObjectStrings, startObject
            raise
            findObject = []
        else:
          findObject = startObject

    return findObject


  def getCcpnMapValue(self,ccpnMap):

    ccpnMapValue = None

    if type(ccpnMap) == type(''):
      ccpnMapValue = self.getCcpnObject(ccpnMap)

    elif type(ccpnMap) == type( () ):
      (ccpnVarName,getFunc) = ccpnMap
      try:
        ccpnMapValue = getFunc(self.getCcpnObject(ccpnVarName) )
      except:
        print "  Cannot CCPN map value for %s, function %s" % (ccpnVarName,getFunc)
        raise
        
    elif ccpnMap == None:
    
      pass

    else:
      try:
        ccpnMapValue = ccpnMap()
      except:
        self.setMessage("  Error: ccpnLoop can only be string, function or tuple(string,function). Is %s." % str(ccpnMap) )

    return ccpnMapValue


  def setExtraEntryLevelTags(self, ccpnObject):

    for saveFrameName in self.nmrStarDict.sfList:
      self.setExtraAppDataTags(ccpnObject, saveFrameName)


  def setExtraAppDataTags(self, ccpnObject, saveFrameName, starElement = None):

    starElementFlag = False

    if starElement is None:
      starElementFlag = True

    appDataDict = self.getAppDataForSaveFrame(ccpnObject, saveFrameName)

    #print 'APP: [%s]' % appDataDict

    #appDataDict = {}

    for mainKey in appDataDict.keys():
      mainDict = self.nmrStarDict.sfDict[mainKey]

      starElementName = mainDict['name']

      for otherKey in appDataDict[mainKey].keys():
        my_temp_dict = appDataDict[mainKey][otherKey]

        if starElementFlag:
          starElement = self.nmrStarFile.setupSaveFrame(mainKey, appDataDict[mainKey][otherKey]['title'])

          starElement.setTag(self.nmrStarFile.tagStart + starElementName + self.nmrStarFile.tagSep + 'Sf_framecode',
                             appDataDict[mainKey][otherKey]['title'])

        for i in range(len(appDataDict[mainKey][otherKey]['tagNames']) ):
          tagName = appDataDict[mainKey][otherKey]['tagNames'][i]

          starTagName = starElementName + self.nmrStarFile.tagSep + tagName

          value = appDataDict[mainKey][otherKey]['tagValues'][i]

          starElement.setTag(self.nmrStarFile.tagStart + starTagName, value)

        if 'tableNames' in appDataDict[mainKey][otherKey]:

          for tableName in appDataDict[mainKey][otherKey]['tableNames']:
            starTable = starElement.setupTable(tableName)

            if not appDataDict[mainKey][otherKey].has_key(tableName):
              continue

            tableTagNames = appDataDict[mainKey][otherKey][tableName]['tableTagNames']

            tagLength = len(appDataDict[mainKey][otherKey][tableName]['tableTagValues'])

            for j in range(len(tableTagNames) ):
              starTableTagName = tableName + self.nmrStarFile.tagSep + tableTagNames[j]

              for k in range(tagLength):
                tableValue = appDataDict[mainKey][otherKey][tableName]['tableTagValues'][k][j]

                starTable.setTag(starTableTagName, tableValue)

  def getAppDataForSaveFrame(self,ccpnObject,saveFrameName):
    
    #
    # This allows mapping of a saveframe name to another one to get at application data.
    # Introduced for general_distance_constraints (Wim 29/09/09)
    #
    
    if self.ccpn2NmrStar.mapAppDataSaveFrames.has_key(saveFrameName):
      searchSaveFrameName = self.ccpn2NmrStar.mapAppDataSaveFrames[saveFrameName][0]
      mapTableNames = self.ccpn2NmrStar.mapAppDataSaveFrames[saveFrameName][1]
    else:
      searchSaveFrameName = saveFrameName
      mapTableNames = {}
      
    #
    # Now start to look for application data, put in dictionary
    #

    AppDataDict = {}
    
    if not ccpnObject or not hasattr(ccpnObject, 'findAllApplicationData'):
      return AppDataDict

    appDataList = ccpnObject.findAllApplicationData(application = self.format,keyword = searchSaveFrameName)
    
    if not appDataList:
      return AppDataDict
    else:
      self.setupAllAppDataDict(ccpnObject)
      self.curAppDataCcpnObject = ccpnObject
    
    appDataValues = [int(appData.value) for appData in appDataList]

    appDataValues.sort()

    mainKey = '%s' % saveFrameName

    if not mainKey in AppDataDict:
      AppDataDict[mainKey] = {}

    mainDict = self.nmrStarDict.sfDict[mainKey]

    maxPos, mainExtraKeys = self.getExtraKeys(mainDict)

    for appDataValue in appDataValues:

      # This is a temporary ID to track within CCPN - no real NMR-STAR value to this
      sfTempId = str(appDataValue)

      saveFrameTitle = self.getAppDataValueSfLevel(searchSaveFrameName,sfTempId,'title')

      tmpSfName = self.getAppDataValueSfLevel(searchSaveFrameName,sfTempId,'tagNames')

      saveFrameTagNames = []

      if tmpSfName:
        saveFrameTagNames  = eval(tmpSfName)

      tmpSfTag = self.getAppDataValueSfLevel(searchSaveFrameName,sfTempId,'tags')

      saveFrameTagValues = []

      if tmpSfTag:
        saveFrameTagValues = eval(tmpSfTag)

      #print '[' + sfTempId + '] [' + str(saveFrameTitle) + '] [' + str(saveFrameTagNames) + '] [' + str(saveFrameTagValues) + ']'

      saveFrameTagNames2  = [None] * maxPos
      saveFrameTagValues2 = [None] * maxPos

      for tagName in mainExtraKeys:
        pos = mainExtraKeys[tagName][0]
        foreignTag = mainExtraKeys[tagName][2]
        
        saveFrameTagNames2[pos] = tagName
        saveFrameTagValues2[pos] = None
        
        # Fixed to set foreign tags if they exist. Wim 09/01/15
        if foreignTag and foreignTag not in self.ignoreForeignLinks:
          value = self.getForeignValue(foreignTag, tagName)
          saveFrameTagValues2[pos] = value

        if saveFrameTagValues2[pos] == None:
          if tagName not in saveFrameTagNames:
            value = mainExtraKeys[tagName][1]
            saveFrameTagValues2[pos] = value
  
          else:
            tagIndex = saveFrameTagNames.index(tagName)     
            saveFrameTagValues2[pos] = saveFrameTagValues[tagIndex]

      otherKey = '%s' % sfTempId

      if not otherKey in AppDataDict[mainKey]:
        AppDataDict[mainKey][otherKey] = {}

      #print 'KEYS: [' + mainKey + '] [' + otherKey +']'

      AppDataDict[mainKey][otherKey]['title'] = saveFrameTitle
      AppDataDict[mainKey][otherKey]['tagNames'] = saveFrameTagNames2
      AppDataDict[mainKey][otherKey]['tagValues'] = saveFrameTagValues2
      
      saveFrameTableNamesAppData = self.getAppDataValueSfLevelList(searchSaveFrameName,'','tables', verbose = 0)

      #print '[' + sfTempId + '] [' + str(saveFrameTableNames) + ']'

      if not saveFrameTableNamesAppData:
        continue

      else:
        
        for saveFrameTableNameStringList in saveFrameTableNamesAppData:
          
          saveFrameTableNames = eval(saveFrameTableNameStringList)
          
          for saveFrameTableName in saveFrameTableNames:
          
            # In case of application data saveFrame name mapping, also have to convert table names!
            if mapTableNames.has_key(saveFrameTableName[1:]):
              actualSaveFrameTableName = "_" + mapTableNames[saveFrameTableName[1:]]
            else:
              actualSaveFrameTableName = saveFrameTableName
              
            #
            # Now continue with getting app data
            #
  
            if not 'tableNames' in AppDataDict[mainKey][otherKey]:
              AppDataDict[mainKey][otherKey]['tableNames'] = []
  
            AppDataDict[mainKey][otherKey]['tableNames'].append(actualSaveFrameTableName)
  
            #print 'TABLE1: [%s] [%s]' % (saveFrameTableName, sorted(self.nmrStarDict.sfDict[saveFrameName]['tables'].keys() ) )
  
            if not self.nmrStarDict.sfDict[mainKey]['tables'].has_key(actualSaveFrameTableName[1:]):
              continue
  
            mainTableDict = self.nmrStarDict.sfDict[mainKey]['tables'][actualSaveFrameTableName[1:]]
  
            tableKey = "%s_table_%s%s" % (searchSaveFrameName,sfTempId,saveFrameTableName)
  
            tmpSfTag = self.getAppDataValueTableLevel(tableKey,'tagNames')
  
            saveFrameTableTagNames = []
  
            if tmpSfTag:
              saveFrameTableTagNames = eval(tmpSfTag)
  
            maxTablePos, mainExtraTableKeys = self.getExtraKeys(mainTableDict)
  
            saveFrameTableTagNames2 = [None] * maxTablePos
  
            #print 'TABLE: [' + str(saveFrameTableTagNames) + ']'
  
            for tTagName in mainExtraTableKeys:
              pos = mainExtraTableKeys[tTagName][0]
              saveFrameTableTagNames2[pos] = tTagName
  
            if not actualSaveFrameTableName in AppDataDict[mainKey][otherKey]:
              AppDataDict[mainKey][otherKey][actualSaveFrameTableName] = {}
  
            if not 'tableTagNames' in AppDataDict[mainKey][otherKey][actualSaveFrameTableName]:
              AppDataDict[mainKey][otherKey][actualSaveFrameTableName]['tableTagNames'] = eval(str(saveFrameTableTagNames2) )
  
            tableRow = 1
  
            if not 'tableTagValues' in AppDataDict[mainKey][otherKey][actualSaveFrameTableName]:
              AppDataDict[mainKey][otherKey][actualSaveFrameTableName]['tableTagValues'] = []
  
            while(tableRow):
  
              tableRowAppData = self.getAppDataValueTableLevel(tableKey,str(tableRow) )
  
              if not tableRowAppData:
                break
  
              tableRowValues = eval(tableRowAppData)
  
              tableRowValues2 = [None] * maxTablePos
  
              for tTagName in mainExtraTableKeys:
                tpos = mainExtraTableKeys[tTagName][0]
                foreignTag = mainExtraTableKeys[tTagName][2]
                
                tableRowValues2[tpos] = None
          
                # Fixed to set foreign tags if they exist. Wim 09/01/15
                if foreignTag and foreignTag not in self.ignoreForeignLinks:
                  value = self.getForeignValue(foreignTag, tTagName)
                  tableRowValues2[tpos] = value

                if tableRowValues2[tpos] == None:
                  if tTagName not in saveFrameTableTagNames:
                    value = mainExtraTableKeys[tTagName][1]
                    tableRowValues2[tpos] = value
    
                  else:
                    tTagIndex = saveFrameTableTagNames.index(tTagName)
                    tableRowValues2[tpos] = tableRowValues[tTagIndex]
                    
              AppDataDict[mainKey][otherKey][actualSaveFrameTableName]['tableTagValues'].append(tableRowValues2)
  
              tableRow += 1

    #print AppDataDict

    return AppDataDict

  def setupAllAppDataDict(self,ccpnObject):
    
    """
    This function grabs all application data from a ccpnObject and organises it into a dictionary.
    Solely done for speedup reasons, so can avoid findFirst() functions.
    """
    
    self.allAppDataDict = {}
    
    allAppData = ccpnObject.findAllApplicationData(application = self.format)
    
    for appData in allAppData:
      keyword = appData.keyword
      value = appData.value
      if not self.allAppDataDict.has_key(keyword):
        self.allAppDataDict[keyword] = []
        
      self.allAppDataDict[keyword].append(value)

  def getAppDataValueSfLevel(self,saveFrameName,saveFrameKey,insertText, verbose = True):

    if saveFrameKey:
      saveFrameKey = "_%s" % saveFrameKey

    keyword = "%s_%s%s" % (saveFrameName,insertText,saveFrameKey)

    return self.getAppDataValue(keyword, verbose = verbose)


  def getAppDataValueSfLevelList(self,saveFrameName,saveFrameKey,insertText, verbose = True):

    if saveFrameKey:
      saveFrameKey = "_%s" % saveFrameKey

    keyword = "%s_%s%s" % (saveFrameName,insertText,saveFrameKey)

    return self.getAppDataValueList(keyword, verbose = verbose)


  def getAppDataValueTableLevel(self,tableKey,appDataId):

    keyword = "%s_%s" % (tableKey,appDataId)

    return self.getAppDataValue(keyword, verbose = False)


  def getAppDataValueTableLevelList(self,tableKey,appDataId):

    keyword = "%s_%s" % (tableKey,appDataId)

    return self.getAppDataValueList(keyword, verbose = False)


  def getAppDataValue(self,keyword, verbose = True):
 
    if self.allAppDataDict.has_key(keyword):
      # This should always be single value, if things are set up correctly at least
      value = self.allAppDataDict.pop(keyword)[0]
      
    else:
      value = None
      
      # Necessary to distinguish between local and global verbosity!
      if verbose and self.verbose:
        print "  Error: no %s appData for ccpn object %s!" % (keyword,self.curAppDataCcpnObject)

    return value


  def getAppDataValueList(self,keyword, verbose = True):

    if self.allAppDataDict.has_key(keyword):
      value = self.allAppDataDict.pop(keyword)
      
    else:
      value = []
      
      # Necessary to distinguish between local and global verbosity!
      if verbose and self.verbose:
        print "  Error: no %s appDataList for ccpn object %s!" % (keyword,self.curAppDataCcpnObject)

    return value

  def getExtraKeys(self, mainDict):

    mainExtraKeys = {}

    pos = 0

    for tagName in mainDict['tagNames']:
      mainExtraKeys[tagName] = (pos, mainDict['tags'][tagName][0], mainDict['tags'][tagName][2])
      pos += 1

    maxPos = pos

    return maxPos, mainExtraKeys
