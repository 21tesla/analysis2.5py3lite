import glob, os, time, shutil, sys

from memops.api import Implementation

from memops.general.Io import getCcpFileString

from ccp.general.Io import getChemCompArchiveXmlFilePath, getChemCompCoordArchiveXmlFilePath

from ccpnmr.format.general.Io import TextPipe

from pdbe.chemComp.Constants import editChemCompDataDir, testChemCompDataDir
from pdbe.chemComp.Constants import obsoleteChemCompDataDir, editLinkedChemCompDataDir

from pdbe.chemComp.Constants import editChemCompCoordDataDir, testChemCompCoordDataDir
from pdbe.chemComp.Constants import obsoleteChemCompCoordDataDir, editLinkedChemCompCoordDataDir
from pdbe.chemComp.Constants import tempChemCompDir

from pdbe.chemComp.Util import getGuidFromXmlFileName, tagAsFresh, removeFreshTag

def findChemCompOrCoordFilePath(chemCompOrCoord, testMode = True):

  """
  Gets the existing file path for the XML file for a chemComp(Coord) object.
    
  testMode        If False, will use archive/, otherwise will use test/.
  """

  if chemCompOrCoord.className == 'ChemCompCoord':

    if testMode:
      dataDir = testChemCompCoordDataDir
    else:
      dataDir = editChemCompCoordDataDir

    xmlFileName = getChemCompCoordFilePath(dataDir,chemCompOrCoord.sourceName,chemCompOrCoord.molType,chemCompOrCoord.ccpCode)

  else:

    if testMode:
      dataDir = testChemCompDataDir
    else:
      dataDir = editChemCompDataDir
      
    xmlFileName = getChemCompFilePath(dataDir,chemCompOrCoord.molType,chemCompOrCoord.ccpCode)

  return xmlFileName
  
def getNewChemCompOrCoordFilePath(chemCompOrCoord, mode = 'test'):


  """
  Gets a new file path for the XML file for a chemComp(Coord) object.
    
  mode        If 'test' will use test/
              If 'archive' will use archive/
              If 'obsolete' will use obsolete/

  """

  molType = chemCompOrCoord.molType
  ccpCode = chemCompOrCoord.ccpCode
  guid = chemCompOrCoord.guid
  
  xmlFileName = "%s+%s+%s.xml" % (molType,getCcpFileString(ccpCode),guid)
  
  if mode == 'obsolete':
    xmlFileName += ".%s" % time.strftime("%Y_%m_%d_%H_%M_%S")

  if chemCompOrCoord.className == 'ChemCompCoord':
  
    xmlFileName = "%s+%s" % (chemCompOrCoord.sourceName,xmlFileName)

    if mode == 'test':
      dataDir = testChemCompCoordDataDir
    elif mode == 'obsolete':
      dataDir = obsoleteChemCompCoordDataDir
    elif mode == 'archive':
      dataDir = editChemCompCoordDataDir
    elif mode == 'archiveLink':
      dataDir = editLinkedChemCompCoordDataDir

    xmlFileDir = getChemCompCoordArchiveXmlFilePath(dataDir,chemCompOrCoord.sourceName,molType,ccpCode)

  else:
  
    if mode == 'test':
      dataDir = testChemCompDataDir
    elif mode == 'obsolete':
      dataDir = obsoleteChemCompDataDir
    elif mode == 'archive':
      dataDir = editChemCompDataDir
    elif mode == 'archiveLink':
      dataDir = editLinkedChemCompDataDir
      
    xmlFileDir = getChemCompArchiveXmlFilePath(dataDir,molType,ccpCode)
  
  # Create if doesn't exist yet!!!
  if not os.path.exists(xmlFileDir):
    os.makedirs(xmlFileDir)

  xmlFilePath = os.path.join(xmlFileDir,xmlFileName)

  return xmlFilePath

def getChemCompFilePath(archiveDataDir,molType,ccpCode):

  fileSearchString = "%s+%s+*.xml" % (molType,getCcpFileString(ccpCode))
  searchPath = getChemCompArchiveXmlFilePath(archiveDataDir,molType,ccpCode)
    
  xmlFileNameMatches = glob.glob(os.path.join(searchPath,fileSearchString))

  xmlFileName = None
  if xmlFileNameMatches:
    xmlFileName = xmlFileNameMatches[-1]

  return xmlFileName
    
def getChemCompCoordFilePath(archiveDataDir,sourceName,molType,ccpCode):

  fileSearchString = "%s+%s+%s+*.xml" % (sourceName,molType,getCcpFileString(ccpCode))
  searchPath = getChemCompCoordArchiveXmlFilePath(archiveDataDir,sourceName,molType,ccpCode)
    
  xmlFileNameMatches = glob.glob(os.path.join(searchPath,fileSearchString))

  xmlFileName = None
  if xmlFileNameMatches:
    xmlFileName = xmlFileNameMatches[-1]

  return xmlFileName

def findExistingChemCompInfo(chemCompArchiveDataDir,ccpCode,molType):

  existingGuid = None
  existingFile = False
  
  for archiveDataDir in (editChemCompDataDir,testChemCompDataDir):
    
    xmlFileName = getChemCompFilePath(archiveDataDir,molType,ccpCode)
    
    if xmlFileName:
      if archiveDataDir == chemCompArchiveDataDir:
        existingFile = True
      existingGuid = getGuidFromXmlFileName(xmlFileName)
      break
      
  return (existingGuid,existingFile)

def findExistingChemCompCoordInfo(chemCompCoordArchiveDataDir,sourceName,ccpCode,molType):

  existingGuid = None
  existingFile = False
  
  for archiveDataDir in (editChemCompCoordDataDir,testChemCompCoordDataDir):
    
    xmlFileName = getChemCompCoordFilePath(archiveDataDir,sourceName,molType,ccpCode)
    
    if xmlFileName:
      if archiveDataDir == chemCompCoordArchiveDataDir:
        existingFile = True
      existingGuid = getGuidFromXmlFileName(xmlFileName)
      break
      
  return (existingGuid,existingFile)

def saveChemCompOrCoord(chemCompOrCoord,testMode=True,replace=False):
  
  """
  Saves the XML for a chemComp(Coord) object.
  
  testMode          If False, will save to archive/, otherwise will use test/.
  replace           If False, will not overwrite existing file.
  
  """
  
  (filePath,existingFilePath) = saveTemporaryChemCompOrCoord(chemCompOrCoord,testMode = testMode)
  consolidateTemporaryChemCompOrCoord(chemCompOrCoord,filePath,existingFilePath,testMode=testMode,replace=replace)  

def saveTemporaryChemCompOrCoord(chemCompOrCoord,testMode = True, isFresh = False):

  """
  Saves the XML for a chemComp(Coord) object to a temporary directory.
  
  testMode          If False, will use archive/, otherwise will use test/ for existingFilePath.  
  isFresh           If False, will remove tag that indicates chemComp(Coord) generated from scratch
                    If True, will set the tag.
  """
  
  #
  # Set/remove tag that indicates whether chemComp is newly generated or not
  #

  if isFresh:
    tagAsFresh(chemCompOrCoord)
  else:
    removeFreshTag(chemCompOrCoord)


  if chemCompOrCoord.className == 'ChemCompCoord':
    className = 'ChemCompCoord'
  else:
    className = 'ChemComp'
  
  tempRepositoryName = 'temporaryChemComp'
  tempRepository = chemCompOrCoord.root.findFirstRepository(name = tempRepositoryName)
  
  if not tempRepository:
    tempRepository = chemCompOrCoord.root.newRepository(name= tempRepositoryName, 
                                                        url=Implementation.Url(path=tempChemCompDir))

  curStdOut = sys.stdout
  textPipe = TextPipe([])
  sys.stdout = textPipe
  chemCompOrCoord.saveTo(tempRepository)
  sys.stdout = curStdOut
  #sys.stdout = sys.__stdout__
  
  print
  print "  *** CCPN save output ***"
  for text in textPipe.textArea:
    if text.strip():
      print "    %s" % text.strip()
  print
  
  # This was failing horribly in large runs, now fixed (Wim 2010/03/08)
  filePath = glob.glob(os.path.join(tempRepository.url.path,'ccp','molecule',className,"*%s*%s*.xml" % (chemCompOrCoord.molType,getCcpFileString(chemCompOrCoord.ccpCode))))[0]
  
  (localPath,localFileName) = os.path.split(filePath)
    
  existingFilePath = findChemCompOrCoordFilePath(chemCompOrCoord, testMode = testMode)
  
  return (filePath,existingFilePath)
  
def consolidateTemporaryChemCompOrCoord(chemCompOrCoord,filePath,existingFilePath,testMode=True,replace=False):
  
  """
  Consolidates a temporary XML file for a chemComp(Coord) object.
  
  testMode          If False, will copy existing file to obsolete/, otherwise will just overwrite.
  replace           If False, will not overwrite existing file.
  
  """

  #
  # If exists, and replace is True, move the old file and move in the new one.
  #
  
  if existingFilePath:
    
    if replace:
    
      # Only store in obsolete/ the files from the archive directory!!
      if not testMode:
        print "  Moving to obsolete/ directory: %s" %  existingFilePath
        obsoleteFilePath = getNewChemCompOrCoordFilePath(chemCompOrCoord, mode = 'obsolete')  
        shutil.move(existingFilePath,obsoleteFilePath)
       
      shutil.move(filePath,existingFilePath)
  
    else:
      print "  Warning: Not overwriting existing file %s." % existingFilePath
      os.remove(filePath)

  else:

    if testMode:
      dirMode = 'test'
    else:
      dirMode = 'archive'
  
    newFilePath = getNewChemCompOrCoordFilePath(chemCompOrCoord, mode = dirMode)
    
    print "  Creating new file %s" % newFilePath
    shutil.move(filePath,newFilePath)
    
    # Also create link in default directory, if necessary
    if dirMode == 'archive':
      linkedFilePath = getNewChemCompOrCoordFilePath(chemCompOrCoord, mode = 'archiveLink')
      if not os.path.exists(linkedFilePath) and linkedFilePath != newFilePath:
        os.symlink(newFilePath,linkedFilePath)
