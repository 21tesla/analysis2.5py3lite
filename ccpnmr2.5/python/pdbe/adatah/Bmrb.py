import os, re
from pdbe.adatah.Constants import archivesDataDir, bmrbUrl, bmrbRestUrl

from pdbe.adatah.Io import getTextFromHttp, getDataFromHttp, getReferenceTextFileFromHttp

#
# Constants and functions related to the BMRB
#

###################################
#                                 #
# BMRB main archive file handling #
#                                 #
###################################

bmrbDataDir = os.path.join(archivesDataDir,'bmrb')
bmrbArchiveDataDir = os.path.join(bmrbDataDir,'archive')
bmrbArchiveUrlLocation = "%s/bmrb/NMR-STAR2/" % bmrbRestUrl

bmrbReferenceDir = os.path.join(bmrbDataDir,'reference')

def getBmrbArchiveEntryList():

  """
  Read list of available BMRB archive files

  Formerly called getBmrbEntryList!
  """

  bmrbFiles = os.listdir(bmrbArchiveDataDir)
  bmrbNumbers = {}
  
  for bmrbFile in bmrbFiles:
    if bmrbFile[-4:] == '.str':
      bmrbCode = bmrbFile[3:-4]
      bmrbNumbers[int(bmrbCode)] = bmrbCode # Should this be 3:? - previous line should do that
  
  bmrbNumberList = bmrbNumbers.keys()
  bmrbNumberList.sort()
  
  bmrbCodes = []
  for bmrbNumber in bmrbNumberList:
    bmrbCodes.append(bmrbNumbers[bmrbNumber])
  
  return bmrbCodes

def getNewBmrbEntries(forceGet = False, skipCodes = None):
  
  """
  Get new BMRB entries...
  """
  
  #
  # Initialize
  #
  
  newBmrbCodeList = []

  codePattern = re.compile("\<a href\=\"bmr\d+\.str\"\>(.+)\<\/a\>")
  
  #
  # Get list of existing entries...
  #
  
  bmrbCodes = getBmrbArchiveEntryList()
  
  #
  # Get the full list of entries from the BMRB...
  #

  dataLines = getTextFromHttp(bmrbArchiveUrlLocation)
  
  #
  # Loop over entry lines
  #
  
  for dataLine in dataLines:
    if dataLine:
      
      codeSearch = codePattern.search(dataLine)
      
      if codeSearch:
        
        bmrbFile = codeSearch.group(1)
        bmrbCode = bmrbFile.replace('.str','')

        if not forceGet and bmrbCode in bmrbCodes:
          continue

        if getBmrbEntry(bmrbCode, skipCodes = skipCodes, forceGet = forceGet):
          newBmrbCodeList.append(bmrbCode)
              
  return newBmrbCodeList

def getBmrbEntry(bmrbId,skipCodes = None, forceGet = False):

  if skipCodes and bmrbId in skipCodes:
    return False
    
  bmrbFile = "%s.str" % bmrbId

  newBmrbFile = os.path.join(bmrbArchiveDataDir,bmrbFile)
  
  if not forceGet and os.path.exists(newBmrbFile):
    return True
  
  #
  # Download
  #

  print "Downloading %s..." % bmrbFile

  bmrbFileLink = '%s%s' % (bmrbArchiveUrlLocation,bmrbId)

  data = getDataFromHttp(bmrbFileLink)
  
  #
  # Exit if no data
  #
  
  if not data:
    print "  Download failed!"
    return False

  #
  # Write the new file
  #

  if not os.path.exists(bmrbArchiveDataDir):
    os.mkdir(bmrbArchiveDataDir)

  fout = open(os.path.join(newBmrbFile),'w')
  fout.write(data)
  fout.close()
  
  return True
  
def getBmrbCodeShortInfo(bmrbCode, bmrbCodeShortInfo = None, overwrite = False, saveEntry = True):

  import re
  
  strucGenPatt = re.compile("structural\s*genomics",re.I)

  #
  # If dictionary not provided, load it (if exists)
  #
  
  pickledDictFilePath = os.path.join(bmrbArchiveDataDir,'reference','bmrbCodeShortInfo.pp')  
  
  if not bmrbCodeShortInfo:
  
    from pdbe.analysis.Util import getPickledDict
    print "  Loading BMRB short info file..."
    bmrbCodeShortInfo = getPickledDict(pickledDictFilePath)

  #
  # Set the info if not yet available (and bmrbCode passed in!)
  #
  
  if bmrbCode and (not bmrbCodeShortInfo.has_key(bmrbCode) or overwrite):
  
    #
    # Get the relevant NMR-STAR info
    #

    from ccp.format.nmrStar.projectIO import NmrStarProjectFile
    from pdbe.analysis.Util import saveReferencePickledDict
    
    nmrStarFilePath = os.path.join(bmrbArchiveDataDir,"%s.str" % bmrbCode)
    
    if not os.path.exists(nmrStarFilePath):
      print "BMRB file for %s not available!" % bmrbCode
      return bmrbCodeShortInfo

    nmrStarFile = NmrStarProjectFile(nmrStarFilePath)
    nmrStarFile.read()
    
    #
    # These variables are set by code below...
    #

    origin = pH = temperature = date = shiftRef = isParamagnetic = None
    isStrucGenomics = False

    #
    # Get entry information level data - only ever one saveframe
    #
        
    saveFrame = nmrStarFile.sfs['entry_information'][0]
    # TODO THIS WILL NOT WORK WITH 3.1!
    if saveFrame.tables.has_key("_Author_ordinal"):
      origin = saveFrame.tables["_Author_ordinal"].tags['_Author_family_name'][-1]
     
    # TODO THIS WILL NOT WORK WITH 3.1!
    if saveFrame.tags.has_key("_Submission_date"):
      # Format is YEAR-MM-DD
      date = tuple(saveFrame.tags["_Submission_date"].split('-'))
      
    #
    # Get sequence info, from both shift and polymer bits, see if matches as well
    #
    # Note that sequence info will remain None if inconsistencies were found!!!
    #
    
    seqInfo = {}
    
    for sequenceFile in nmrStarFile.sequenceFiles:
      for sequence in sequenceFile.sequences:
        seqInfo[sequence.chainName] =  {'molName': sequence.molName,
                                        'molType': sequence.polymerType,
                                        'sequence': [(seqEl.seqId,seqEl.code3Letter) for seqEl in sequence.elements]}  
            
    chemShiftInfo = {}
    for chemShiftFile in nmrStarFile.chemShiftFiles:
      for chemShift in chemShiftFile.chemShifts:
        chainName = chemShift.molCode
        if not chemShiftInfo.has_key(chainName):
          chemShiftInfo[chainName] = []
          
        # Warning: resLabel can be missing for non-standard stuff
        if hasattr(chemShift,'resLabel'):
          seqEl = (chemShift.seqCode,chemShift.resLabel)
          if not seqEl in chemShiftInfo[chainName]:
            chemShiftInfo[chainName].append(seqEl)
          
    for chainName in chemShiftInfo.keys():
      if not seqInfo.has_key(chainName):
        print 'Error: chain info missing (from chemical shifts), not using for BmrbCodeShortInfo dictionary'
        seqInfo = {}
      
      else:
        for seqEl in chemShiftInfo[chainName]:
          if not seqEl in seqInfo[chainName]['sequence']:
            print 'Error: sequence element info missing (from chemical shifts), not using for BmrbCodeShortInfo dictionary'
            seqInfo = {}
            break
      
      if not seqInfo:
        break
    
    #
    # Get sample conditions information - ONLY USING FIRST SF!
    #
    
    saveFrameName = 'sample_conditions'
    if nmrStarFile.sfs.has_key(saveFrameName):

      saveFrame = nmrStarFile.sfs[saveFrameName][0]
      
      if 'pH' in saveFrame.tables["_Variable_type"].tags['_Variable_type']:
        index = saveFrame.tables["_Variable_type"].tags['_Variable_type'].index('pH')
        pH = saveFrame.tables["_Variable_type"].tags['_Variable_value'][index]

      if 'temperature' in saveFrame.tables["_Variable_type"].tags['_Variable_type']:
        index = saveFrame.tables["_Variable_type"].tags['_Variable_type'].index('temperature')
        temperature = saveFrame.tables["_Variable_type"].tags['_Variable_value'][index]
    
    #
    # Is this molSystem paramagnetic?
    #
    
    saveFrameName = 'molecular_system'
    tagKey = '_System_paramagnetic'

    if nmrStarFile.sfs.has_key(saveFrameName):
      saveFrame = nmrStarFile.sfs[saveFrameName][0]
      if saveFrame.tags.has_key(tagKey):
        value = saveFrame.tags[tagKey]

        if value == 'yes':
          isParamagnetic = True
        elif value == 'no':
          isParamagnetic = False
        else:
          isParamagnetic = value
    
    #
    # Get chemical shift referencing info - ONLY USING FIRST SF!
    #

    saveFrameName = 'chemical_shift_reference'
    tableKey = '_Mol_common_name'

    # TODO again won't work for 3.1!
    if nmrStarFile.sfs.has_key(saveFrameName):

      saveFrame = nmrStarFile.sfs[saveFrameName][0]
      shiftRef = {}
      
      if saveFrame.tables.has_key(tableKey):
 
        if saveFrame.tags.has_key("_Details"):
          details = saveFrame.tags['_Details']
        else:
          details = None
     
        table = saveFrame.tables[tableKey]

        for i in range(len(table.tags[tableKey])):

          atomType   = table.tags['_Atom_type'][i]
          refMol     = table.tags[tableKey][i]

          if table.tags.has_key('_Reference_method'):
            refMethod = table.tags['_Reference_method'][i]
          else:
            refMethod = None

          if table.tags.has_key('_Reference_type'):
            refType = table.tags['_Reference_type'][i]
          else:
            refType = None
            
          if table.tags.has_key('_Chem_shift_units'):
            unit = table.tags['_Chem_shift_units'][i]
          else:
            unit = None
            
          shiftValue = table.tags['_Chem_shift_value'][i]
          
          if table.tags.has_key('_Atom_group'):
            atomGroup = table.tags['_Atom_group'][i]
          else:
            atomGroup = None

          if table.tags.has_key('_Indirect_shift_ratio'):
            shiftRatio = table.tags['_Indirect_shift_ratio'][i]
          else:
            shiftRatio = None

          shiftRef[atomType] = (refMol,refMethod,refType,unit,shiftValue,atomGroup,shiftRatio,details)
      
    #
    # Get submission date
    #
            
    #
    # Is it structural genomics?
    #
    # 1. Search in entry and citation title
    #
    
    for (saveFrameName,tagTuple) in  (('entry_information',('_Entry_title','Title')),('entry_citation',('_Citation_title','Title'))):

      if nmrStarFile.sfs.has_key(saveFrameName):

        saveFrame = nmrStarFile.sfs[saveFrameName][0]

        # Handle 2.1.1 and 3.1
        for tagName in tagTuple:
          if saveFrame.tags.has_key(tagName):
            title = " ".join(saveFrame.tags[tagName].split("\n"))
            break

        #print saveFrameName, title
        if strucGenPatt.search(title):
          isStrucGenomics = True
          break
        
    #
    # 2. Also search in main citation keywords
    #  
 
    if not isStrucGenomics:
    
      saveFrameName = 'entry_citation'
      if nmrStarFile.sfs.has_key(saveFrameName):
      
        saveFrame = nmrStarFile.sfs[saveFrameName][0]
        
        # TODO THIS WILL NOT WORK WITH 3.1!
        if saveFrame.tables.has_key("_Keyword"):
          # Handle 2.1.1 and 3.1
          for tagName in ("_Keyword","Keyword"):
            if saveFrame.tables["_Keyword"].tags.has_key(tagName):
              for value in saveFrame.tables["_Keyword"].tags[tagName]:
                if value and strucGenPatt.search(value):
                  isStrucGenomics = True
                  break
      
    #
    # Finally set and save info...
    #

    bmrbCodeShortInfo[bmrbCode] = (isStrucGenomics,isParamagnetic,origin,pH,temperature,date,shiftRef,seqInfo)
    
    if saveEntry:
      # Forcing a write here because 'deep' structure and default saveReferencePickledDict doesn't
      # pick up on this...
      
      print "  Saving reference information for %s..." % bmrbCode
      
      saveReferencePickledDict(pickledDictFilePath,bmrbCodeShortInfo,forceWrite=True)
  
  return bmrbCodeShortInfo
  

########################################
#                                      #
# Get unique BMRB to PDB code mappings #
#                                      #
########################################

bmrbPdbMappingFile = 'adit_nmr_matched_pdb_bmrb_entry_ids.csv'
bmrbPdbMappingUrl = "%s/ftp/pub/bmrb/nmr_pdb_integrated_data/%s" % (bmrbUrl,bmrbPdbMappingFile)
bmrbPdbMappingFilePath = os.path.join(bmrbReferenceDir,bmrbPdbMappingFile)

def getBmrbPdbMappingInfo(dataFilePath=None):

  # NOTE: this gives UNIQUE mappings between PDB and BMRB ID!
  # Use getBmrbDatabaseMatches for list of all mappings.

  if not dataFilePath:
    dataFilePath = bmrbPdbMappingFilePath

  dataLines = getReferenceTextFileFromHttp(bmrbPdbMappingUrl,dataFilePath,refText = "BMRB to PDB mappings")
  
  
  bmrbPdbMappingDict = {}
  
  for dataLine in dataLines:
    if dataLine:
      (bmrbId,pdbCode) = dataLine.strip().split(',')
      
      bmrbId = int(bmrbId)
      
      if not bmrbPdbMappingDict.has_key(bmrbId):
        bmrbPdbMappingDict[bmrbId] = pdbCode.lower()
      else:
        # This should never happen... raise it?
        print "MULTIPLE MATCH FOR %d!!" % bmrbId
      
  return bmrbPdbMappingDict

##################################
#                                #
# Get all BMRB database mappings #
#                                #
##################################

bmrbDatabaseMatchFile = "dbmatch.csv"
bmrbDatabaseMatchUrl = "%s/ftp/pub/bmrb/internal_data/webdata/%s" % (bmrbUrl,bmrbDatabaseMatchFile)
bmrbDatabaseMatchFilePath = os.path.join(bmrbReferenceDir,bmrbDatabaseMatchFile)

def getBmrbDatabaseMatches(dataFilePath=None):

  """
  Read info on relation between BMRB and other database entries...
  """

  if not dataFilePath:
    dataFilePath = bmrbDatabaseMatchFilePath

  dataLines = getReferenceTextFileFromHttp(bmrbDatabaseMatchUrl,dataFilePath,refText = "BMRB database matches")

  bmrbDatabaseDict = {}

  for dataLine in dataLines:
    if dataLine:
      dataLine = dataLine.strip()
      values = dataLine.split(',')
      
      # Ignore lines with inconsistent length
      if len(values) != 11:
        continue
      
      (bmrbId,dbName,dbCode,matchData,matchScore,unknown,unknown,method,resolution,unknown,unknown) = values

      bmrbId = bmrbId.strip('"')
      if not bmrbId:
        continue        

      bmrbId = int(bmrbId)
      dbName = dbName.strip('"')
      dbCode = dbCode.strip('"')
      
      if not bmrbDatabaseDict.has_key(bmrbId):
        bmrbDatabaseDict[bmrbId] = {}
      if not bmrbDatabaseDict[bmrbId].has_key(dbName):
        bmrbDatabaseDict[bmrbId][dbName] = []     
      
      bmrbDatabaseDict[bmrbId][dbName].append(dbCode)

  return bmrbDatabaseDict

###############################
#                             #
# Get BMRB info on data types #
#                             #
###############################

def getBmrbShiftInfo():

  return getBmrbInfo(dataType='shift')

def getBmrbCouplingInfo():
  
  return getBmrbInfo(dataType='coupling')

def getBmrbRdcInfo():
  
  return getBmrbInfo(dataType='rdc')

def getBmrbT1Info():
  
  return getBmrbInfo(dataType='t1')

def getBmrbT2Info():
  
  return getBmrbInfo(dataType='t2')

def getBmrbHetNuclNoeInfo():
  
  return getBmrbInfo(dataType='hetNuclNoe')

def getBmrbOrderParamInfo():

  return getBmrbInfo(dataType='orderParam')

def getBmrbHExchangeInfo():

  return getBmrbInfo(dataType='hExchange')

def getBmrbHProtectionInfo():

  return getBmrbInfo(dataType='hProtection')

def getBmrbInfo(dataType=None):

  if not dataType:
    return None
  
  elif dataType == 'rdc':
    bmrbQueryFile = "query_1_67.html"
    valueNames = ('bmrbId','name','numRdcs','hasProtein','hasDNA','hasRNA')

  elif dataType == 'coupling':
    bmrbQueryFile = "query_1_29.html"
    valueNames = ('bmrbId','name','numCouplings','hasProtein','hasDNA','hasRNA')

  elif dataType == 't1':
    bmrbQueryFile = "query_1_33.html"
    valueNames = ('bmrbId','name','numT1s','hasProtein','hasDNA','hasRNA')

  elif dataType == 't2':
    bmrbQueryFile = "query_1_38.html"
    valueNames = ('bmrbId','name','numT2s','hasProtein','hasDNA','hasRNA')

  elif dataType == 'hetNuclNoe':
    bmrbQueryFile = "query_1_41.html"
    valueNames = ('bmrbId','name','numHetNuclNoes','hasProtein','hasDNA','hasRNA')
    
  elif dataType == 'shift':
    bmrbQueryFile = "query_1_5_4.html"
    valueNames  = ('bmrbId','name','numShifts_1H','numShifts_13C','numShifts_15N','numShifts_31P','hasProtein','hasDNA','hasRNA')
    
  elif dataType == 'orderParam':
    bmrbQueryFile = "query_1_45.html"
    valueNames  = ('bmrbId','name','numOrderParams','hasProtein','hasDNA','hasRNA')
    
  elif dataType == 'hExchange':
    bmrbQueryFile = "query_1_49.html"
    valueNames  = ('bmrbId','name','numHExchangeValues','hasProtein','hasDNA','hasRNA')

  elif dataType == 'hProtection':
    bmrbQueryFile = "query_1_53.html"
    valueNames  = ('bmrbId','name','numHProtectionValues','hasProtein','hasDNA','hasRNA')

  else:
    print "Unknown data type %s, aborting..." % dataType
    return None

  bmrbQueryUrl = "%s/search/query_grid/%s" % (bmrbUrl,bmrbQueryFile)
  bmrbQueryFilePath = os.path.join(bmrbReferenceDir,bmrbQueryFile)


  entryLinePatt = re.compile("data_library/generate_summary")
  valuePatt = re.compile("\>([^\<\>]+)\<")

  bmrbValuesDict = {}
  
  dataLines = getReferenceTextFileFromHttp(bmrbQueryUrl,bmrbQueryFilePath,refText = "BMRB %s information" % dataType)
    
  ln = 0
  
  while (ln < len(dataLines)):

    dataLine = dataLines[ln]

    if entryLinePatt.search(dataLine):
      infoDict = {}
      for valueName in valueNames:
      
        valueSearch = valuePatt.search(dataLine)
        if not valueSearch:
          print dataLine
          
        value = valueSearch.group(1)
        
        if value in ('&nbsp;','&nbsp'):
          value = False
        elif value == 'X':
          value = True
        elif valueName[:3] == 'num':
          value = int(value)
          
        infoDict[valueName] = value
        
        ln += 1
        dataLine = dataLines[ln]
        
      bmrbId = int(infoDict['bmrbId'])
      del(infoDict['bmrbId'])
      
      bmrbValuesDict[bmrbId] = infoDict

    ln += 1

  return bmrbValuesDict

############################################################
#                                                          #
# Special NMR-STAR handler for BMRB 2.1(.1) archive files  #
#                                                          #
# Also redefining NMRproject itself for specific handling  #
#                                                          #
############################################################


class NmrStarBmrbArchiveHandler:

  def readCustomNmrStarFile(self,nmrStarFileName, components = None):

    #
    # Read chemical shift data from nmrStar format
    #

    nmrStarFile = CustomNmrStarProjectFile(nmrStarFileName)

    localKeywds = {}

    if self.presets.has_key('readProject') and self.presets['readProject'].has_key('keywds'):
      localKeywds = self.presets['readProject']['keywds']

    localKeywds['components'] = ['sequence','measurements']
    if components:
      localKeywds['components'] = components

    self.formatObjectDict['NmrStar'].file = nmrStarFile

    self.formatObjectDict['NmrStar'].readProject(nmrStarFileName, entry = self.entry, molSystem = self.molSystem, version = '2.1.1', minimalPrompts = 1, chemCompPath = curChemCompRepository, **localKeywds)

    self.formatObjectDict['NmrStar'].nameMapping.isOriginalImport = True
    
    # Check if import worked - don't bother if not!
    if not self.molSystem.chains:
      return False
    else:
      return nmrStarFile

from ccp.format.nmrStar.projectIO import NmrStarProjectFile
from ccp.format.general.Constants import defaultMolCode

class CustomNmrStarProjectFile(NmrStarProjectFile):

  def readStatusCheck(self):

    returnStatus = False
    
    self.pdbCodes = []
    self.seqOffset = 0
    
    #
    # Check whether the data is monomeric and associated with a PDB file
    #
    # TODO EXTEND TO LIGANDS AND NON-PDB!!
    #
    
    if not self.sfs.has_key('ligand') and self.sfs.has_key('monomeric_polymer'):
      
      #
      # Only monomers, only one chain...
      #
      
      if len(self.sfs['monomeric_polymer']) == 1:
      
        curSf = self.sfs['monomeric_polymer'][0]

        if curSf.tables.has_key('_Database_name'):

          curTable = curSf.tables['_Database_name']

          for i in range(0,len(curTable.tags['_Database_name'])):

            dbType = curTable.tags['_Database_name'][i]

            if dbType == 'PDB':

              pdbCode = curTable.tags['_Database_accession_code'][i]

              if pdbCode not in self.pdbCodes:
                self.pdbCodes.append(pdbCode)

              returnStatus = True
    
    #
    # Check whether chemical shift data available
    #
    
    if returnStatus:
    
      if not self.sfs.has_key('assigned_chemical_shifts'):
    
        returnStatus = False     
      
    if returnStatus:
    
      for curSf in self.sfs['assigned_chemical_shifts']:
      
        #
        # Set chain code to ' ' for chemical shifts...
        #

        origName = curSf.tags['_Mol_system_component_name']
        curSf.tags['_Mol_system_component_name'] = defaultMolCode
      
      #
      # Check if renumbering of seq codes is necessary... 
      # NOTE: only for bmr < 4000 files, from then on Residue_seq_code used (???)
      #
      
      curSf = self.sfs['monomeric_polymer'][0]
      curTable = curSf.tables['_Residue_seq_code']

      if curTable.tags.has_key('_Residue_author_seq_code'):
        
        for seqId in range(0,len(curTable.tags['_Residue_author_seq_code'])):
          
          try:
            #
            # Necessary because sometimes codes are missing (?)
            #
            self.seqOffset = curTable.tags['_Residue_author_seq_code'][seqId] - 1 - seqId
            break
            
          except:
            pass

    
    return returnStatus


