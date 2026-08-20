import os, re

from pdbe.adatah.Constants import archivesDataDir, bmrbUrl
from pdbe.adatah.Io import getDataFromHttp, getReferenceTextFileFromHttp

from subprocess import Popen

#####################################
#                                   #
# NMR restraints grid file handling #
#                                   #
#####################################

nmrGridDataDir = os.path.join(archivesDataDir,'bmrb','nmrRestrGrid')
nmrGridReferenceDir = os.path.join(nmrGridDataDir,'reference')

statusConvert = {'2-parsed': 'parsed','3-converted-DOCR':'converted','4-filtered-FRED':'filtered'}

# Location of Jurgen's file with the nmrGrid information...
bmrbInfoUrl = "%s/servlet_data/viavia" % bmrbUrl

nmrGridInfoFile = "mrfile.txt"
nmrGridInfoUrl = "%s/mr_mysql_backup/%s" % (bmrbInfoUrl,nmrGridInfoFile)

# Location of servlet to get specific nmrGrid info
nmrGridServletUrl = "%s/servlets/MRGridServlet" % bmrbUrl

mrArchiveLinkPatt = re.compile("\&request_type\=archive\&")
mrLinkPatt = re.compile("MRGridServlet\?([^\>]+)\"\>")

# Joined coordinate info at the BMRB (formerly at CMBI)
#cmbiInfoUrl = "http://nmr.cmbi.ru.nl/~jd/viavia"
#cmbiGridInfoUrl = "%s/NRG/star" % cmbiInfoUrl
bmrbJoinedCoordUrl = "http://sunfish.bmrb.wisc.edu/nmrRestrGrid/"

def getBmrbNmrGridDict(dataFilePath=None, skipCodes = None):

  """
  Read info on DOCR/FRED from BMRB...
  """

  bmrbNmrGridDict = {}
  
  #1047	1017	classified	1a51	2004-12-02
  #26561	1064	parsed	1ajw	2004-12-03 
  #39888	1263	filtered	1c54	2004-12-03
  #39889	1266	converted	1c7v	2004-12-03
  # Where classified means there is a constraint file, but it's not been parsed
  # rest is obvious: converted DOCR, filtered FRED.

  if not dataFilePath:
    dataFilePath = os.path.join(nmrGridReferenceDir,nmrGridInfoFile)

  dataLines = getReferenceTextFileFromHttp(nmrGridInfoUrl,dataFilePath,refText = "NMR GRID information")
 
  for dataLine in dataLines:
    if dataLine:
      (gridId,bmrbId,statusCode,pdbCode,date) = dataLine.split()
      
      if skipCodes and pdbCode in skipCodes:
        continue
      
      if statusCode in statusConvert:
        if not bmrbNmrGridDict.has_key(pdbCode):
          bmrbNmrGridDict[pdbCode] = {}
          
        status = statusConvert[statusCode]

        if not bmrbNmrGridDict[pdbCode].has_key(status):
          webLink = '%s?pdb_id=%s&min_items=0&block_text_type=%s' % (nmrGridServletUrl,pdbCode,statusCode)
          bmrbNmrGridDict[pdbCode][status] = webLink        

  return bmrbNmrGridDict

def getBmrbNmrGridFile(pdbCode,forceGet = False, statusCode = '2-parsed'):
  
  """
  Get NMR restraints grid NMR-STAR file
  """

  pdbDir = os.path.join(nmrGridDataDir,pdbCode)

  if not os.path.exists(pdbDir):
    os.mkdir(pdbDir)
    
  nmrStarFile = "%s.str" % (statusConvert[statusCode])
  
  if not os.path.exists(os.path.join(pdbDir,nmrStarFile)) or forceGet:
  
    origDir = os.getcwd()

    # Location zipped file
    urlLocation = "%s?request_type=archive&pdb_id=%s&block_text_type=%s&file_detail=%s" % (nmrGridServletUrl,pdbCode,statusCode,statusCode)

    processBmrbNmrGridArchiveLink(urlLocation,pdbDir,nmrGridDataDir,nmrStarFile)

    os.chdir(origDir)

  return 1


def getBmrbNmrGridJointCoordinateFile(pdbCode, forceGet = False):

  """
  Get NMR restraints grid NMR-STAR file joined with coordinate info from mmCIF file
  """

  pdbDir = os.path.join(nmrGridDataDir,pdbCode)

  if not os.path.exists(pdbDir):
    os.mkdir(pdbDir)
    
  nmrStarFile = "joinedCoord.str.gz"
  
  if not os.path.exists(os.path.join(pdbDir,nmrStarFile)) or forceGet:
  
    origDir = os.getcwd()

    # Location zipped file
    urlLocation = "%s/%s/%s" % (bmrbJoinedCoordUrl,pdbCode,nmrStarFile)

    data = getDataFromHttp(urlLocation)

    localFile = os.path.join(pdbDir,nmrStarFile)
    localOut = open(localFile,'w')
    localOut.write(data)
    localOut.close()
    
    Popen(['gunzip','-f',localFile])
    
    os.chdir(origDir)

  return True


def restrFileHasIntensities(restrPdbCode):

  #
  # Load the dictionary (if exists)
  #
  
  from pdbe.analysis.Util import getPickledDict, saveReferencePickledDict
  
  pickledDictFilePath = os.path.join(nmrGridReferenceDir,'restrFileOrigIntensities.pp')  
  restrFileOrigIntensities = getPickledDict(pickledDictFilePath)

  if not restrFileOrigIntensities.has_key(restrPdbCode):
    
    #
    # Get the relevant NMR-STAR info
    #

    from ccp.format.nmrStar.distanceConstraintsIO import NmrStarFile
    
    nmrStarFilePath = os.path.join(nmrGridDataDir,restrPdbCode,'parsed.str')
    
    if not os.path.exists(nmrStarFilePath):
      print "NMR Restraints Grid file for %s not available!" % restrPdbCode
      return False

    nmrStarFile = NmrStarFile(nmrStarFilePath)
    nmrStarFile.read()

    restrFileOrigIntensities[restrPdbCode] = False
    for nmrStarConstraintFile in nmrStarFile.constraintFiles:
      for constraint in nmrStarConstraintFile.constraints:
        for node in constraint.nodes:
          if hasattr(node,'intensity'):
            if node.intensity != None:
              restrFileOrigIntensities[restrPdbCode] = True
    
    saveReferencePickledDict(pickledDictFilePath,restrFileOrigIntensities)
    
  if restrFileOrigIntensities[restrPdbCode]:
    hasOrigIntensities = True
  else:
    hasOrigIntensities = False
  
  return hasOrigIntensities

def getBmrbNmrGridCompletenessInfo(pdbCode, forceWrite = False):

  """
  Read info on completeness check for.a PDB entry at the BMRB...
  """

  fileName = 'completeness.str'
  pdbDir = os.path.join(nmrGridDataDir,pdbCode)
  
  finalFileName = os.path.join(pdbDir,fileName)
  
  if forceWrite or not os.path.exists(finalFileName):
    finalFileName = None
  
  if not finalFileName:
    urlLocation = '%s?db_username=wattos1&format=distance&pdb_id=%s&request_type=block_set&subtype=completeness&type=check' % (nmrGridServletUrl,pdbCode)

    data = getDataFromHttp(urlLocation)
    dataLines = data.split("\n")

    for dataLine in dataLines:
      if mrArchiveLinkPatt.search(dataLine):
        mrLinkSearch = mrLinkPatt.search(dataLine)
        if mrLinkSearch:
          archiveLink = "%s?%s" % (nmrGridServletUrl,mrLinkSearch.group(1))

          if not os.path.exists(pdbDir):
            os.mkdir(pdbDir)

          processBmrbNmrGridArchiveLink(archiveLink,pdbDir,nmrGridDataDir,fileName)
          
          finalFileName = os.path.join(pdbDir,fileName)
          
          break

        else:
          print "  Error: no link for completeness data for %s!" % pdbCode
          
  return finalFileName
    
def processBmrbNmrGridArchiveLink(urlLocation,saveDir,homeDir,finalFileName):

    data = getDataFromHttp(urlLocation)
    
    zipFileName = 'file.zip'
    zipFile = os.path.join(saveDir,zipFileName)
    zipOut = open(zipFile,'w')
    zipOut.write(data)
    zipOut.close()
    
    os.chdir(saveDir)
    os.spawnlp(os.P_WAIT, 'nice', 'nice', '-19', 'unzip', zipFileName)
    os.chdir(homeDir)

    os.remove(zipFile)
    
    removeFiles = ['readme.html','main.css','index.csv']
    for removeFile in removeFiles:
      removeFilePath = os.path.join(saveDir,removeFile)
      if os.path.exists(removeFilePath):
        os.remove(removeFilePath)

    files = os.listdir(saveDir)

    for fileName in files:
      if fileName[:5] == 'block' and fileName[-4:] == '.str':
        os.rename(os.path.join(saveDir,fileName),os.path.join(saveDir,finalFileName))        



"""

OBSOLETE - is there another file being generated?

from pdbe.adatah.NmrRestrGrid import bmrbInfoUrl

bmrbPdbNmrMatchFile = "score.csv"
bmrbPdbNmrMatchUrl = "%s/bmrb_pdb_match/%s" % (bmrbInfoUrl,bmrbPdbNmrMatchFile)
bmrbPdbNmrMatchFilePath = os.path.join(bmrbReferenceDir,bmrbPdbNmrMatchFile)

def getBmrbPdbNmrMatches(nmrOnly = False):

  ####
  Read info on relation between BMRB and PDB entries...
  
  Info below: the higher the score, the worse. Entries with overall 9 are ignored in final BMRB mappings.
  
  etsScore:     BMRB-PDB tracking since BMRB entry 4000
  authorScore:  Match between author names
  nmrScore:     Is this NMR structure. Only NMR structures in there, always empty
  blastScore:   Sequence match
  ligandScore:  Ligand match
  overallScore: Overall score
  ####
  
  bmrbPdbDict = {}
  
  dataLines = getReferenceTextFileFromHttp(bmrbPdbNmrMatchUrl,bmrbPdbNmrMatchFilePath,refText = "BMRB to PDB NMR matches")
    
  # First line is header  
  for dataLine in dataLines[1:]:
    if dataLine:
      (bmrbId,pdbCode,etsScore,authorScore,nmrScore,blastScore,ligandScore,overallScore) = dataLine.split(',')

      bmrbId = int(bmrbId)
      
      # Ignore some bad matches automatically...
      if '9' in (ligandScore,authorScore,blastScore,etsScore):
        continue
        
      if overallScore == '9':
        continue
        
      if not bmrbPdbDict.has_key(bmrbId):
        bmrbPdbDict[bmrbId] = []
      
      bmrbPdbDict[bmrbId].append(pdbCode)

  return bmrbPdbDict

"""

# NOT SURE IF BELOW IS STILL REQUIRED - CHECK!!!
"""
def getFixedCoordStarFile(curDir,pdbCode,forceGet = False):
    
  pdbDir = os.path.join(curDir,pdbCode)
  nmrStarFile = '%s_extra.str' % pdbCode

  if not forceGet and os.path.exists(os.path.join(pdbDir,nmrStarFile)):
    return 0

  # Location of Jurgen's link stuff...
  urlLocation = "http://tang.bmrb.wisc.edu/wattos/link/%s/%s" % (pdbCode,nmrStarFile)

  r1 = urllib.urlopen(urlLocation)
  data = r1.read()
  r1.close()

  outFileName = os.path.join(pdbDir,nmrStarFile)
  outFile = open(outFileName,'w')
  outFile.write(data)
  outFile.close()

  return 1


def getNonRedunStarFile(curDir,pdbCode,forceGet = False):
    
  pdbDir = os.path.join(curDir,pdbCode)
  nmrStarFile = 'nonredun.str'
  
  if not forceGet and os.path.exists(os.path.join(pdbDir,nmrStarFile)):
    return 0

  # Location of Jurgen's link stuff...
  urlLocation = "http://tang.bmrb.wisc.edu/wattos/link/%s/%s" % (pdbCode,nmrStarFile)

  r1 = urllib.urlopen(urlLocation)
  data = r1.read()
  r1.close()

  outFileName = os.path.join(pdbDir,nmrStarFile)
  outFile = open(outFileName,'w')
  outFile.write(data)
  outFile.close()

  return 1

"""
