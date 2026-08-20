#
# Contains constants, functions required to handle PDB files
#

import re, os

from memops.universal.Util import returnInt

from ccpnmr.format.process.matchResonToMolSys import matchCoordAtomsToMolSys

from pdbe.adatah.Constants import chemCompArchiveDataDir, archivesDataDir
from pdbe.adatah.Constants import pdbFtp, pdbFtpDir, ebiPdbFtp, atlasUrl
from pdbe.adatah.Io import getFileFromFtp, getReferenceTextFileFromHttp, getTextFromHttp

#
# Constants
#

pdbDataDir = os.path.join(archivesDataDir,'pdb')
pdbReferenceDir = os.path.join(pdbDataDir,'reference')

#
# For looking for valid PDB code...
#

# WAS PREVIOUSLY validCodeSearch!
validPdbCodeSearch = re.compile("[a-z0-9]{4}")

import os

#
# Pdb handler class
#

class PdbHandler:

  # TODO should I pass in pdb code here? Or use self.idCode directly? Problem if merging two info sets - probably
  # best to leave it as is.

  def readPdbCoordinates(self,pdbCode,chemCompArchiveDataDir = chemCompArchiveDataDir, mergeNmrStarFile = None, maxNum = 999, readPeopleCitations=True):

    from ccp.format.pdb import coordinatesIO, peopleAndCitationsIO

    #
    # Read sequence and coordinates from PDB (to get original name and such)
    #

    pdbFileName = self.getPdbFileName(pdbCode)

    # TODO TODO: Might need to reset some descriptor codes for ChemCompVars... problems with GLU, ASP.
    #            This is because the names are wrong in the PDB file - should be OK with WHATIF renamed files!

    keywds = {}
    if self.presets.has_key('readCoordinates'):

      scriptPresets = self.presets['readCoordinates']

      if scriptPresets.has_key('keywds'):

        keywds = scriptPresets['keywds']

    if not keywds.has_key('forceChainMappings'):

      # TODO TODO
      # Note that it's a bit silly to read coord file and read it again later, but
      # it's all automatic anyway so doesn't matter too much... could also 
      # pass in file if necessary. Might speed things up a bit!

      coordinateFile = coordinatesIO.PdbCoordinateFile(pdbFileName)
      coordinateFile.read(maxNum = 1)

      modelCoordKeys = coordinateFile.modelCoordinates.keys()
      modelCoordKeys.sort()

      refModel = modelCoordKeys[0]
      coords = coordinateFile.modelCoordinates[refModel]

      print "  Trying to automap coordinate atoms..."
      forceChainMappings = matchCoordAtomsToMolSys(coords,self.molSystem,test = 0)
      if forceChainMappings:
        keywds['forceChainMappings'] = forceChainMappings

    keywds['ignoreUnknownChemComps'] = True

    #
    # Special hack to merge info from NMR-STAR file in PDB file before linking coordinates.
    #

    if mergeNmrStarFile:
      keywds['mergeFunctionInfo'] = (mergeNmrStarWithPdb,mergeNmrStarFile)

    #
    # Read in PDB file coordinates, merge NMRSTAR info if necessary...
    #

    initKeywds = keywds.copy()
    self.formatObjectDict['Pdb'].readCoordinates(pdbFileName, molSystem = self.molSystem, strucGen = self.strucGen, linkAtoms = False, forceReadSequence = True, minimalPrompts = 1, chemCompPath = chemCompArchiveDataDir, maxNum = maxNum, **keywds)

    #
    # Also get people/citation info!
    # (via hack to avoid having to read file again)
    #
    
    if readPeopleCitations:

      self.formatObjectDict['Pdb'].peopleAndCitationsFile = peopleAndCitationsIO.PdbPeopleAndCitationsFile('temp')
      self.formatObjectDict['Pdb'].peopleAndCitationsFile.read(pdbFile = self.formatObjectDict['Pdb'].coordinateFile)
      self.formatObjectDict['Pdb'].file = 'temp'

      self.formatObjectDict['Pdb'].readPeopleAndCitations(None)
      self.formatObjectDict['Pdb'].file = None

    return initKeywds

  def getPdbFileName(self,pdbCode):
  
    return os.path.join(pdbDataDir, '%s.pdb' % pdbCode)

  # WAS FORMERLY returnValidPdbCcpnFile!
  # TODO Does this appear somewhere else as well!?!?
  def getValidPdbCcpnDir(self,projectDir,pdbCode,ccpnDataDir):

    validProjectDir = None

    validCode = validCodeSearch.search(pdbCode)

    if validCode and pdbCode not in ignorePdbCodes:

      pdbCodeDir = os.path.join(projectDir,pdbCode)

      if os.path.exists(pdbCodeDir):

        ccpnDir = os.path.join(pdbCodeDir,ccpnDataDir)

        if os.path.exists(ccpnDir):

          memopsDir = os.path.join(ccpnDir,"memops","Implementation")
          
          if os.path.exists(memopsDir):
          
            fileNames = os.listdir(memopsDir)
            
            for fileName in fileNames:
              if fileName[-4:] == '.xml':
                validProjectDir = ccpnDir
                break

    return validProjectDir

#
# General functions
#

# WAS FORMERLY getPdbCodeNew!
def getPdbCode(pdbCode, forceGet = False, pdbDir = None, source=''):

  if not pdbDir:
    pdbDir = pdbDataDir

  pdbCode = pdbCode.lower()
  pdbFileName = os.path.join(pdbDir,"%s.pdb" % pdbCode)
  
  if source == 'EBI':
    pdbFtpDir = ebiPdbFtpDir
    pdbFtp = ebiPdbFtp
    pdbInsertDir = os.path.join('data','structures','divided','pdb')
  else:
    pdbInsertDir = ""
    from pdbe.adatah.Constants import pdbFtpDir, pdbFtp

  if not os.path.exists(pdbFileName) or forceGet:

    #
    # Try downloading
    #

    currentFtpDir = os.path.join(pdbFtpDir,pdbInsertDir,"%s" % pdbCode[1:3])
    ftpFileName = "pdb%s.ent.gz" % pdbCode
    localFileName = os.path.join(pdbDir,"%s.gz" % pdbFileName)
        
    try:      
      getFileFromFtp(pdbFtp,currentFtpDir,ftpFileName,localFileName)
     
      if forceGet:
        os.remove(os.path.join(pdbDir,pdbFileName))

      # TODO: should be POpen object instead!!
      os.spawnlp(os.P_WAIT, 'gunzip', 'gunzip', localFileName)     

    except:
      
      print "Error: no PDB file for %s!! Directory not recognized by this computer?" % pdbCode
      return False

  return True


ebiPdbFtpDir = "pub/databases/pdb/" # removed "rcsb" to reflect changes in EBI's ftp structure
ebiPdbDerivedDataUrl = os.path.join("ftp://%s" % ebiPdbFtp,ebiPdbFtpDir,'derived_data')

xrayResolutionFile = "resolu.idx"
xrayResolutionUrl = os.path.join(ebiPdbDerivedDataUrl,'index',xrayResolutionFile)
xrayResolutionFilePath = os.path.join(pdbReferenceDir,xrayResolutionFile)

pdbEntryTypeFile = "pdb_entry_type.txt"
pdbEntryTypeUrl = os.path.join(ebiPdbDerivedDataUrl,pdbEntryTypeFile)
pdbEntryTypeFilePath = os.path.join(pdbReferenceDir,pdbEntryTypeFile)

pdbChainInfoFile = "pdb_seqres.txt"
pdbChainInfoUrl = os.path.join(ebiPdbDerivedDataUrl,"%s.gz" % pdbChainInfoFile)
pdbChainInfoFilePath = os.path.join(pdbReferenceDir,pdbChainInfoFile)

def getXrayResolution():

  """
  Read info on X-ray resolution PDB entries
  """

  resolutionDict = {}

  dataLines = getReferenceTextFileFromHttp(xrayResolutionUrl,xrayResolutionFilePath,refText = "PDB x-ray resolution information")
  
  for dataLine in dataLines:
    
    cols = dataLine.split()
    
    if len(cols) == 3 and cols[1] == ';':
      pdbCode = cols[0].lower()
      resolution = float(cols[2])
      
      if resolution != -1.0:
        resolutionDict[pdbCode] = resolution
        
  return resolutionDict

def getPdbEntryType(dataFilePath=None,updateFile=False):

  """
  Read info on experiment type of PDB entries
  """
  
  pdbEntryTypeDict = {}
  
  if not dataFilePath:
    dataFilePath = pdbEntryTypeFilePath
  
  if updateFile or not os.path.exists(dataFilePath):
    dataLines = getReferenceTextFileFromHttp(pdbEntryTypeUrl,dataFilePath,refText = "PDB entry experiment type information")
  else:
    fin = open(dataFilePath)
    dataLines = fin.readlines()
    fin.close()
  
  for dataLine in dataLines:
  
    cols = dataLine.split()
    
    if len(cols) == 3:
    
      (pdbCode,molType,expType) = cols
    
      # molType is prot, nuc or prot-nuc
      # expType is diffraction, NMR, or EM
    
      pdbEntryTypeDict[pdbCode] = (molType,expType)
      
  return pdbEntryTypeDict

def getPdbChainInfo():

  """
  Read info on chain information of PDB entries
  """
  
  pdbChainInfoDict = {}
   
  dataLines = getReferenceTextFileFromHttp(pdbChainInfoUrl,pdbChainInfoFilePath,refText = "PDB entry chain information", isGzipped = True)
  
  for dataLine in dataLines:
  
    cols = dataLine.split()
    
    if cols and cols[0][0] == '>':
      
      (pdbCode,chainCode) = cols[0][1:].split("_")
      molType = cols[1][4:]
      chainLength = returnInt(cols[2][7:])
      molName = ' '.join(cols[3:])
      
      if not pdbChainInfoDict.has_key(pdbCode):
        pdbChainInfoDict[pdbCode] = {}
      
      pdbChainInfoDict[pdbCode][chainCode] = [chainLength,molName,None]
    
    else:
    
      pdbChainInfoDict[pdbCode][chainCode][2] = cols[0]
      
  return pdbChainInfoDict

def getRelatedEntriesRcsb(pdbCode):

  """
  Get list of related PDB entries
  """
  
  relatedPdbCodes = []
  
  findClusterUrl = "http://www.rcsb.org/pdb/explore/sequenceCluster.do?structureId=%s" % pdbCode.upper()
  dataLines = getTextFromHttp(findClusterUrl)
  clusterId = None
  fullMatchPatt = re.compile("\>100\%\<\/a\>")
  clusterPatt = re.compile("cluster\=(\d+)\&amp")
  
  for dataLine in dataLines:
    if fullMatchPatt.search(dataLine):
      clusterSearch = clusterPatt.search(dataLine)
      if clusterSearch:
        clusterId = clusterSearch.group(1) 
  
  pdbCodesUrl = os.path.join("http://www.rcsb.org/pdb/explore/sequenceCluster.do?structureId=%s&entity=1&cluster=%s&seqid=100" % (pdbCode.upper(),clusterId))
  dataLines = getTextFromHttp(pdbCodesUrl)
  
  strucLinePatt = re.compile("href\=\"/pdb/explore\.do\?structureId\=([A-Z0-9]+)\"")
  
  for dataLine in dataLines:
    strucLineSearch = strucLinePatt.search(dataLine)
    if strucLineSearch:
      pdbCode = strucLineSearch.group(1)
      relatedPdbCodes.append(pdbCode.lower())

  return relatedPdbCodes  

"""
def getRelatedEntries(pdbCode):

  #Get list of related PDB entries
  
  relatedPdbCodes = []


  # TODO First check if already have a list! Update every month? Or have 'force' option?
  
  
  entryAtlasUrl = os.path.join(atlasUrl,'similarity','%s.html' % pdbCode)
  
  dataLines = getTextFromHttp(entryAtlasUrl)
  relatedEntriesText = "Related PDB entries"
  
  for dataLine in dataLines:
    if dataLine.count("Related PDB entries"):
      bits = dataLine.split("</a>")
      for bit in bits:
        if bit.count("Related PDB entries"):
          continue
        pdbCode = bit[-4:]

        if pdbCode.count(">"):
          continue
          
        relatedPdbCodes.append(pdbCode)

  return relatedPdbCodes  
"""
