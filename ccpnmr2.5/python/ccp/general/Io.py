"""
======================COPYRIGHT/LICENSE START==========================

Io.py: General I/O code for CCPN

Copyright (C) 2008 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/PDBe)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)
- PDBe website (http://www.ebi.ac.uk/pdbe/)

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

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================
"""
#
# Convenient I/O functions
#

import os, urllib, re
import pandas as pd
from memops.universal import Io as uniIo
from memops.universal import Util as uniUtil
from memops.universal.Url import fetchUrl, fetchHttpResponse
from memops.general.Implementation import ApiError
from memops.general.Io import findCcpXmlFile, getCcpFileString
from memops.api import Implementation
from memops.format.xml import XmlIO
import StringIO

from ccp.general.Constants import chemCompServer, chemCompWebPath
from ccp.general.Constants import standardResidueCcpCodes
from ccp.general.Util import setCurrentStore


def getDataPath(*args):
  
  """
  Gives location of data path. Extra args are added on as extra directories
  Result ends with the last arg (which might be either a file or a dictionary)
  """
  dataPath = uniIo.joinPath(uniIo.getTopDirectory(),'data',*args)
  return dataPath

def getChemCompArchiveDataDir():

  """
  Default directory for locally storing all chemComps
  This is now available as a CcpForge repository, see http://ccpforge.cse.rl.ac.uk/projects/ccpn-chemcomp/
  """

  return getDataPath('pdbe','chemComp','archive','ChemComp')

def getChemComp(project, molType, ccpCode, download=True, showError=None,
                partialLoad=False,chemCompArchiveDir = None, copyFile=True):
  
  """ get ChemComp corresponding to molType,ccpCode, 
  looking 1) in memory, 2) in Repositories on lookup path,
  3) in chemCompArchiveDir directory, 4) downloading from PDBe ChemComp server.
  For 3) and 4) save new ChemComp in first Repository on PAckageLocator lookup path
  Do 4) only if download==True
  
  showError is an optional GUI error handler that can be passed in.
  partialLoad controls if only the TopObject (default) or the entire file is loaded
  
  copyFile can be set to False to avoid copying file from central archive (good for testing)
  
  Optimised to avoid mass reading.
  """
  
  
  # First get it if already loaded
  # 1 May 10: below returns None if chemComp not already loaded
  # and then XmlIO.loadFromFile sets chemComp.isModified to True
  # which means that chemComp is saved
  #chemComp = project.getByNavigation(('chemComps',(molType,ccpCode)))
  chemComp = project.findFirstChemComp(molType=molType,ccpCode=ccpCode)
  
  if chemComp is None:
    # try to load it from an existing repository - avoiding mass loading
    packageName = 'ccp.molecule.ChemComp'  
    chemCompFileSearchString = "%s+%s+*.xml" % (molType,getCcpFileString(ccpCode))
    
    chemCompXmlFile = findCcpXmlFile(project, packageName, chemCompFileSearchString)
    if chemCompXmlFile:
      chemComp = XmlIO.loadFromFile(project, chemCompXmlFile, partialLoad=partialLoad)

    if chemComp is None:
      # try to get it from chemCompArchiveDir directory, if any, or to download it.
      # Use custom directory if it was passed in!!!!
      if not chemCompArchiveDir:
        chemCompPath = getChemCompArchiveDataDir()
      else:
        chemCompPath = chemCompArchiveDir
        
      ccLocator = (project.findFirstPackageLocator(targetName=packageName) or
                   project.findFirstPackageLocator(targetName='any'))
      repository = ccLocator.findFirstRepository()
      
      fileFound = getChemCompXmlFile(repository, chemCompPath, molType, ccpCode,
                                    showError=showError,copyFile=copyFile)
      if not fileFound and download:
        # Replaced by direct download from CcpForge (Wim 15/06/2010)    
        #fileFound = downloadChemCompXmlFile(repository, molType, ccpCode,
        #                                    showError=showError)
        fileFound = downloadChemCompInfoFromCcpForge(repository, molType, ccpCode,
                                                     showError=showError)        
      if fileFound:
        chemComp = XmlIO.loadFromFile(project, fileFound, 
                                      partialLoad=partialLoad)
  #
  return chemComp
  
def getChemCompArchiveXmlFilePath(chemCompPath,molType,ccpCode):
    
  chemCompXmlFilePath = uniIo.joinPath(chemCompPath,molType)
  
  if molType == 'other':
    chemCompXmlFilePath = uniIo.joinPath(chemCompXmlFilePath,ccpCode[0])
  
  return chemCompXmlFilePath
  
def getChemCompXmlFile(repository, chemCompPath, molType, ccpCode, 
                       showError=None,copyFile=True):
  
  """
  Fetch chemComp 'molType', 'ccpCode' to local repository 'repository'
  from repository defined by chemCompPath
  showError is an error display function, passed in were appropriate for GUI
  contexts
  Returns name of copied file, or None if unsuccessful
  """  
    
  fileSearchString = "%s+%s+*.xml" % (molType,getCcpFileString(ccpCode))
  fileSearchPath = getChemCompArchiveXmlFilePath(chemCompPath,molType,ccpCode)
  
  className = 'ChemComp'
  identifier = '%s.%s' % (molType,ccpCode)
  
  return getChemCompOrCoordXmlFile(repository,fileSearchString,fileSearchPath,className,identifier,showError,copyFile)

def getChemCompCoordXmlFile(repository, chemCompCoordPath, sourceName, molType, ccpCode, 
                            showError=None, copyFile=True):
  
  """
  Fetch chemCompCoord 'sourceName', 'molType', 'ccpCode' to local repository 'repository'
  from repository defined by chemCompPath
  showError is an error display function, passed in were appropriate for GUI
  contexts
  Returns name of copied file, or None if unsuccessful
  """  
    
  fileSearchString = "%s+%s+%s+*.xml" % (getCcpFileString(sourceName),molType,getCcpFileString(ccpCode))
  fileSearchPath = getChemCompCoordArchiveXmlFilePath(chemCompCoordPath,sourceName,molType,ccpCode)
  
  className = 'ChemCompCoord'
  identifier = '%s.%s.%s' % (sourceName,molType,ccpCode)
  
  return getChemCompOrCoordXmlFile(repository,fileSearchString,fileSearchPath,className,identifier,showError,copyFile)

def getChemCompOrCoordXmlFile(repository,fileSearchString,fileSearchPath,className,identifier,showError,copyFile):  
  
  result = None

  if showError is None:
    showError = uniIo.printError

  # Try to find file...
  
  import glob  
  fileNameMatches = glob.glob(uniIo.joinPath(fileSearchPath,fileSearchString))
  
  if fileNameMatches:
    if len(fileNameMatches) > 1:
      errorText = "Error: multiple matches found for %s %s - taking last one." % (className,identifier)
      showError("Multiple %s matches" % className, errorText)
    
    filePath = fileNameMatches[-1]
    (fileDir,fileName) = os.path.split(filePath)
  
    #
    # Copy file if found...
    #

    if os.path.exists(filePath):
      if copyFile:
        savePath = repository.getFileLocation('ccp.molecule.%s' % className)
        if not os.path.exists(savePath):
          os.makedirs(savePath)
        import shutil
        saveFilePath = uniIo.joinPath(savePath,fileName)
        shutil.copy(filePath,saveFilePath)
        result = saveFilePath
    
        print "  %s file %s copied to %s..." % (className,fileName,savePath)
      
      else:
        result = filePath
  
  return result
 
def getChemCompCoordArchiveXmlFilePath(chemCompPath,sourceName,molType,ccpCode):
    
  chemCompXmlFilePath = uniIo.joinPath(chemCompPath,sourceName,molType)
  
  if molType == 'other':
    chemCompXmlFilePath = uniIo.joinPath(chemCompXmlFilePath,ccpCode[0])
  
  return chemCompXmlFilePath

def getChemCompCoordArchiveDataDir():

  """
  Default directory for locally storing all chemComps
  This is now available as a CcpForge repository, see http://ccpforge.cse.rl.ac.uk/projects/ccpn-chemcomp/
  """

  return getDataPath('pdbe','chemComp','archive','ChemCompCoord')

def getChemCompCoord(project, sourceName, molType, ccpCode, download=True, showError=None,
                     partialLoad=False, chemCompCoordArchiveDir = None, copyFile=True):
  
  """ get ChemComp corresponding to molType,ccpCode, 
  looking 1) in memory, 2) in Repositories on lookup path,
  3) in allChemCompCoordPath directory, 4) downloading from PDBe ChemComp server.
  For 3) and 4) save new ChemComp in first Repository on PAckageLocator lookup path
  Do 4) only if download==True
  
  showError is an optional GUI error handler that can be passed in.
  partialLoad controls if only the TopObject (default) or the entire file is loaded
  
  copyFile can be set to False to avoid copying file from central archive (good for testing)
  
  Optimised to avoid mass reading.
  """
  
  
  # First get it if already loaded
  chemCompCoord = project.getByNavigation(('chemCompCoords',(sourceName,molType,ccpCode)))
  
  if chemCompCoord is None:
    # try to load it from an existing repository - avoiding mass loading
    packageName = 'ccp.molecule.ChemCompCoord'  
    chemCompCoordFileSearchString = "%s+%s+%s+*.xml" % (getCcpFileString(sourceName),
                                                        molType, getCcpFileString(ccpCode))
    
    chemCompCoordXmlFile = findCcpXmlFile(project, packageName, chemCompCoordFileSearchString)
    if chemCompCoordXmlFile:
      chemCompCoord = XmlIO.loadFromFile(project, chemCompCoordXmlFile, partialLoad=partialLoad)

    if chemCompCoord is None:
      # try to get it from allChemCompCoordPath directory, if any, or to download it.
      # Use custom directory if it was passed in!!!!
      if not chemCompCoordArchiveDir:
        chemCompCoordPath = getChemCompCoordArchiveDataDir()
      else:
        chemCompCoordPath = chemCompCoordArchiveDir
        
      ccLocator = (project.findFirstPackageLocator(targetName=packageName) or
                   project.findFirstPackageLocator(targetName='any'))
      repository = ccLocator.findFirstRepository()
      
      fileFound = getChemCompCoordXmlFile(repository, chemCompCoordPath, sourceName, molType, ccpCode,
                                          showError=showError,copyFile=copyFile)
      if not fileFound and download:
        # Replaced by direct download from CcpForge (Wim 15/06/2010)    
        #fileFound = downloadChemCompCoordXmlFile(repository, sourceName, molType, ccpCode,
        #                                    showError=showError)
        fileFound = downloadChemCompInfoFromCcpForge(repository, molType, ccpCode, sourceName=sourceName,
                                                     showError=showError) 
               
      if fileFound:
        chemCompCoord = XmlIO.loadFromFile(project, fileFound, 
                                           partialLoad=partialLoad)

  return chemCompCoord

def _getCcpForgeIndexUrl():
  """
  Creates the CcpForge subPaths/index Path
  """
  ccpForgeUrl = "https://raw.githubusercontent.com/VuisterLab/CcpNmr-ChemComps/master/"
  indexDir = "index"
  indexFile = "index.csv"
  ccpForgeIndexUrl = ccpForgeUrl + uniIo.joinPath(indexDir, indexFile)

  return ccpForgeUrl, indexDir, indexFile, ccpForgeIndexUrl

def getCcpForgeSubPath(molType, ccpCode, sourceName=None):
  """
  Creates the subPath info for ChemComp(Coord) downloads from CcpForge
  """
  ccpForgeUrl, indexDir, indexFile, ccpForgeIndexUrl = _getCcpForgeIndexUrl()

  if not sourceName:
    fileType = 'ChemComp'
    subPath = getChemCompArchiveXmlFilePath("", molType, ccpCode)
    fileSearchString = "%s\+%s\+" % (molType, getCcpFileString(ccpCode))

  else:
    fileType = 'ChemCompCoord'
    subPath = getChemCompCoordArchiveXmlFilePath("", sourceName, molType, ccpCode)
    fileSearchString = "%s\+%s\+%s\+" % (getCcpFileString(sourceName), molType, getCcpFileString(ccpCode))

  # need opposite
  # if fileType == 'ChemComp':
  #   (tmpMolType, tmpCcpCode, suffix) = chemCompXmlFile.split("+")
  # else:
  #   (tmpSourceName, tmpMolType, tmpCcpCode, suffix) = chemCompXmlFile.split("+")

  return (fileType, ccpForgeUrl, ccpForgeIndexUrl, fileSearchString)


def getCcpForgeUrls(molType,ccpCode,sourceName=None):
  
  """
  Creates the URL info for ChemComp(Coord) downloads from CcpForge
  """

  # ccpForge repository

  # ccpForgeUrl = "http://ccpforge.cse.rl.ac.uk/gf/project/ccpn-chemcomp/scmcvs/?action=browse&root=ccpn-chemcomp&pathrev=MAIN&path=/"
  # checkOutDir = "~checkout~"
  # archiveDir = "ccpn-chemcomp/data/pdbe/chemComp/archive/"

  # github repository

  # returns json
  # ccpForgeUrl = "http://api.github.com/repos/VuisterLab/CcpNmr-ChemComps/contents/"

  # returns .xml
  ccpForgeUrl = "http://raw.githubusercontent.com/VuisterLab/CcpNmr-ChemComps/master/"

  # e.g. fetchHttpResponse("http://raw.githubusercontent.com/VuisterLab/CcpNmr-ChemComps/master/data/pdbe/chemComp/archive/ChemComp/other/A/other+A0a+pdbe_ccpnRef_2009-03-12-11-50-15-639_00001.xml")
  # e.g data/pdbe/chemComp/archive/ChemComp/other/A/other+A0a+pdbe_ccpnRef_2009-03-12-11-50-15-639_00001.xml
  checkOutDir = ""
  archiveDir = "data/pdbe/chemComp/archive/"

  # sourceForge repository

  # ccpForgeUrl = "git://git.code.sf.net/p/ccpn-chemcomps/code/ci/master/tree/"
  # checkOutDir = ""
  # archiveDir = "data/pdbe/chemComp/archive/"
  # queryDir = "?format=raw"

  # https: // api.github.com / repos / {user} / {repo_name} / contents / {path_to_file}

  #%2Acheckout%2A%2Fccpn-chemcomp%2Fdata%2Fpdbe%2FchemComp%2Farchive%2FChemComp%2Fprotein%2Fprotein%252B004%252Bpdbe_ccpnRef_2010-09-23-14-41-20-237_00001.xml&revision=1.1
  
  if not sourceName:
    fileType = 'ChemComp'
    subPath = getChemCompArchiveXmlFilePath("",molType,ccpCode)
  
  else:
    fileType = 'ChemCompCoord'
    subPath = getChemCompCoordArchiveXmlFilePath("",sourceName,molType,ccpCode)
  
  # wb104: 15 Dec 2010: need to use joinPath because otherwise does not
  # work on Windows, but cannot stick ccpForgeUrl in joinPath() because
  # that would convert http:// to http:/
  ccpForgeDirUrl = ccpForgeUrl+uniIo.joinPath(archiveDir,fileType,subPath)
  ccpForgeDownloadUrl = ccpForgeUrl+uniIo.joinPath(checkOutDir,archiveDir,fileType,subPath)
  
  return (fileType, ccpForgeDirUrl, ccpForgeDownloadUrl)

def findCcpForgeDownloadLink(dirData,fileType,ccpCode,ccpForgeDownloadUrl):
  
  """
  Finds the relevant XML file name and download link for the ChemComp(Coord) from the repository directory data off CcpForge.
  Works by ccpCode only, assuming that info comes from right directory!
  """
  import json

  urlLocation = chemCompXmlFile = None
  dirDataDict = json.loads(dirData)

  # check that dirData is a list
  if not isinstance(dirDataDict, list):
    return (urlLocation, chemCompXmlFile)

  for entry in dirDataDict:

    # check that each entry is a dict containing the correct keys
    if entry and isinstance(entry, dict) and 'name' in entry and 'download_url' in entry:

      # get the name of the file
      chemCompXmlFile = entry['name']

      if fileType == 'ChemComp':
        (tmpMolType, tmpCcpCode, suffix) = chemCompXmlFile.split("+")
      else:
        (tmpSourceName, tmpMolType, tmpCcpCode, suffix) = chemCompXmlFile.split("+")

      if ccpCode == tmpCcpCode:
        urlLocation = entry['download_url']
        break

  return (urlLocation, chemCompXmlFile)

  # original below
  #
  # fileHtmlPatt = re.compile("a name=\"([^\"]+)\"\s+href\=\"([^\"]+)\"\s+title=\"")
  #
  # urlLocation = chemCompXmlFile = None
  #
  # for urlDirLine in dirData.split("\n"):
  #
  #   fileHtmlSearch = fileHtmlPatt.search(urlDirLine)
  #
  #   if fileHtmlSearch:
  #     chemCompXmlFile = fileHtmlSearch.group(1)
  #
  #     if fileType == 'ChemComp':
  #       (tmpMolType,tmpCcpCode,suffix) = chemCompXmlFile.split("+")
  #     else:
  #       (tmpSourceName,tmpMolType,tmpCcpCode,suffix) = chemCompXmlFile.split("+")
  #
  #     if ccpCode == tmpCcpCode:
  #
  #       urlLocation = "%s/%s&content-type=text/plain" % (ccpForgeDownloadUrl,chemCompXmlFile.replace('+','%252B'))
  #       break
  #
  # return (urlLocation, chemCompXmlFile)

def downloadChemCompInfoFromCcpForge(repository, molType, ccpCode, sourceName=None, showError=None):
  
  """
  Fetch chemComp 'molType', 'ccpCode' or, if sourceName given, the corresponding chemCompCoord, 
  to local repository 'repository'
  from chemCompServer
  showError is an error display function, passed in were appropriate for GUI
  contexts
  Returns name of copied file, or None if unsuccessful
  """
  
  if showError is None:
    showError = uniIo.printError
    
  result = None
  
  (fileType, ccpForgeDirUrl, ccpForgeDownloadUrl) = getCcpForgeUrls(molType,ccpCode,sourceName=sourceName)
  (_, ccpForgeUrl, ccpForgeIndexUrl, searchString) = getCcpForgeSubPath(molType,ccpCode,sourceName=sourceName)

  if sourceName:
    # For displaying error info
    sourceText = "%s, " % sourceName
  else:
    sourceText = ""
    
  try:

    # read the index file from the chemcomp website and then get the file location from the index
    # currently using two methods
    request = fetchHttpResponse(ccpForgeIndexUrl)
    ccIndex = request.read()
    request.close()
    csvString = StringIO.StringIO(ccIndex)
    df = pd.read_csv(csvString)

    urlLocation = None
    if df is not None:
      found = df[df['file'].str.match(searchString)]
      if found is not None and found.size:

        if found.shape[0] != 1:
          raise ValueError('Error: too many search results')

        xmlPath = found['path'].to_numpy()[0]
        xmlFile = found['file'].to_numpy()[0]

        # NOTE:ED must use '/' for url join
        urlLocation = ccpForgeUrl+'/'.join([xmlPath, xmlFile])
      else:
        return None         # nothing found, need to skip here to match original version
    else:
      raise

    try:
      if urlLocation:

        # NOTE ED - can now use  -  r2 = _fetchUrlData(urlLocation)
        #                           data = _readDataFromRequest(r2)

        r2 = fetchHttpResponse(urlLocation)
        try:
          data = r2.read()
          r2.close()

          if data and ('Not Found' in data[:min(20, len(data))] or 'Invalid request' in data[:min(20, len(data))]):

            # NOTE:ED if the server returns file containing 'not found' then don't write
            showError("File not found",
                      "File not found on server: %s" % xmlFile)

          else:
            try:
              saveChemCompPath = repository.getFileLocation('ccp.molecule.%s' % fileType)
              if not os.path.exists(saveChemCompPath):
                os.makedirs(saveChemCompPath)

              # chemCompFile = uniIo.joinPath(saveChemCompPath,chemCompXmlFile)
              chemCompFile = uniIo.joinPath(saveChemCompPath, xmlFile)
              fout = open(chemCompFile,'w')
              fout.write(data)
              fout.close()

              print ("Downloaded %s %s%s, %s from server %s, written to file %s!"
                     % (fileType,sourceText,molType,ccpCode,ccpForgeDownloadUrl,chemCompFile))
              result = chemCompFile

            except IOError, e:
              showError("Cannot write file",
                        "Cannot write %s XML file %s%s, %s: %s"
                        % (fileType,sourceText,molType,ccpCode,str(e)))
  
        except IOError, e:
          showError("Cannot read file", "Cannot read %s %s%s, %s: %s" 
                    % (fileType,sourceText,molType,ccpCode,str(e)))
      
        
      else:
        showError("Cannot find file", "Cannot find %s XML file %s%s, %s." 
                  % (fileType,sourceText,molType,ccpCode))
      
    except IOError, e:
      showError("Cannot read directory", 
                "Cannot read directory information for %s%s, %s: %s" 
                % (sourceText,molType,ccpCode,str(e)))

  except IOError, e:
    showError("No connection", 
              "Cannot connect to download server %s, or file does not exist...: %s " 
              % (ccpForgeDownloadUrl,str(e)))
  #
  return result

def getCcpCodeFromXmlFile(xmlFile):

  ccpCodePatt = re.compile("ccpCode=\"([^ ]+)\"")

  fin = open(xmlFile)
  line = fin.readline()
  
  ccpCode = None
  
  while (line):
    if line.count('<CHEM.StdChemComp') or line.count('<CHEM.NonStdChemComp') or line.count("<CCCO.ChemCompCoord"):
      ccpCodeSearch = ccpCodePatt.search(line)
      if ccpCodeSearch:
        ccpCode = ccpCodeSearch.group(1)
        break

    line = fin.readline()
      
  fin.close()
  
  return ccpCode

def getStdChemComps(project,molTypes=None):

  chemComps = []

  if not molTypes:
  
    molTypes = ['protein','RNA','DNA']

  for molType in molTypes:
  
    if standardResidueCcpCodes.has_key(molType):
  
      for ccpCode in standardResidueCcpCodes[molType]:
        chemComp = getChemComp(project, molType, ccpCode, download=False)
        if chemComp:
          chemComps.append(chemComp)

  return chemComps


def setDataSourceDataStore(dataSource, dataUrlPath, localPath, 
                           dataLocationStore=None, dataUrl=None):
  
  #
  # Get DataLocationStore
  #
  
  if not dataLocationStore:
    
    setCurrentStore(dataSource.root,'DataLocationStore')
    dataLocationStore = dataSource.root.currentDataLocationStore
  
  #
  # Get (or create) DataUrl
  #
  
  # TODO should this search function go elsewhere?
  if not dataUrl:
    for tmpDataUrl in dataLocationStore.dataUrls:
      if tmpDataUrl.url.dataLocation == dataUrlPath:
        dataUrl = tmpDataUrl
    
    if not dataUrl:
      dataUrlPath = uniIo.normalisePath(dataUrlPath)
      dataUrl = dataLocationStore.newDataUrl(url = Implementation.Url(path = dataUrlPath))
      
  #
  # Create a BlockedBinaryMatrix. TODO: could be other classes that are set up this way - rename func and make general,, pass in class?
  #  
  localPath = uniIo.normalisePath(localPath)
  blockedBinaryMatrix = dataLocationStore.newBlockedBinaryMatrix(path=localPath, 
                                                                 dataUrl=dataUrl)
  
  """
  TODO Set here as well, or do this later after returning object:

blockSizes      Int      0..*     Block sizes in dimension order  
complexStoredBy   ComplexStorage   1..1   The ordering of real and imaginary parts of hypercomplex numbers in the data matrix. See ComplexStorage type for details  
hasBlockPadding   Boolean   1..1   Are data padded to fill all blocks completely? Alternatively incomplete blocks store only the actual data.  
headerSize   Int   1..1   Header size in bytes  
isBigEndian   Boolean   1..1   Are data big-endian (alternative little-endian).  
isComplex   Boolean   0..*   Are numbers complex (if True) or real/integer (if False).  
nByte   Int   1..1   Number of bytes per number  
numPoints   Int   0..*   number of points for each matrix dimension - also defines dimensionality of matrix. The number of points is the same for real or complex data, in the sense that n complex points require 2n real numbers for storage.  
numRecords   Int   1..1   Number of matrix records in file. All other information in the object describes a single record.  
numberType   NumberType   1..1   Type of numbers held in matrix  
  
  """
  
  dataSource.dataStore = blockedBinaryMatrix
  
  return blockedBinaryMatrix


def getDataSourceFileName(dataSource):

  dataStore = dataSource.dataStore

  if not dataStore:
    return None

  return dataStore.fullPath


def setDataSourceFileName(dataSource, fileName):

  dataStore = dataSource.dataStore

  if dataStore is None:
    raise ApiError('dataStore is None')

  setDataStoreFileName(dataStore, fileName)

def setDataStoreFileName(dataStore, fileName):

  preferDataUrls=(dataStore.dataUrl,)
  (dataUrl, filePath) = getDataStoringFromFilepath(dataStore.root,
                               fullFilePath=fileName,
                               preferDataUrls=preferDataUrls,
                               dataLocationStore=dataStore.dataLocationStore)

  dataStore.dataUrl = dataUrl
  dataStore.path = filePath

def getDataStoringFromFilepath(memopsRoot, fullFilePath, preferDataUrls=None,
                               dataLocationStore=None, keepDirectories=1):
  
  # make absolute,, normalised path
  fullFilePath = uniIo.normalisePath(fullFilePath, makeAbsolute=True)
  
  dataUrl, filePath = findDataStoringFromFilepath(memopsRoot, fullFilePath, 
                                                  preferDataUrls,
                                                  dataLocationStore,
                                                  keepDirectories)
  
  if dataUrl is None:
  
    urlPath = uniIo.normalisePath((fullFilePath[:-len(filePath)]))
    dataLocationStore = memopsRoot.currentDataLocationStore
    dataUrl = dataLocationStore.newDataUrl(
                                   url=Implementation.Url(path=urlPath))
    dataUrl.name = 'auto-%s' % dataUrl.serial
  #
  return (dataUrl, filePath)

def findDataStoringFromFilepath(project, fullFilePath, preferDataUrls=None,
                               dataLocationStore=None, keepDirectories=1):
  """ Get DataUrl and relative filePath from normalised absolute filePath
  Uses heuristics to select compatible DataUrl from existing ones.
  sisterObjects is a collection of objects with a dataStore link - 
  DataUrls in use for sisterObjects are given preference in the
  heuristics.
  uses dataLocationStore or current dataLocationStore
  If no compatible DataUrl is found the routine returns dataUrl None
  and the file name plus the lowest keepDirectories directories 
  as the filePath
  """
  # NB fullFilePath *must* be absolute her for code to work properly
  # 
  if not os.path.isabs(fullFilePath):
    raise ApiError(
     "findDataStoringFromFilepath called with non-absolute file name %s"
     % fullFilePath)
  
  # get DataLocationStore
  if not dataLocationStore:
    setCurrentStore(project,'DataLocationStore')
    dataLocationStore = project.currentDataLocationStore
  
  # get DataUrl that match fullFilePath
  dataUrls = []
  for dataUrl in dataLocationStore.dataUrls:
    dirPath = uniIo.normalisePath(dataUrl.url.path)
    if fullFilePath.startswith(dirPath):
      lenPath = len(dirPath)
      ss = fullFilePath
      while len(ss) > lenPath:
        ss,junk = uniIo.splitPath(ss)
      if ss == dirPath:
        # DataUrl path matches file path
        dataUrls.append(dataUrl)
  
  # process result
  if dataUrls:
    if preferDataUrls:
      # look for DataUrls that are in use with related objects
      ll = [x for x in dataUrls if x in preferDataUrls]
      if ll:
        dataUrls = ll
        
    if len(dataUrls) == 1:
      # only one DataUrl - use it
      dataUrl = dataUrls[0]
    else:
      # use DataUrl with longest path
      ll = [(len(dataUrl.url.path),dataUrl) for dataUrl in dataUrls]
      ll.sort()
      dataUrl = ll[-1][1]
    
    # get filePath
    if ss == '/':
      ss = ''
    elif uniUtil.isWindowsOS() and ss.endswith('/') and ss[-2] == ':': # deals with case on Windows when one has drive + '/'
      ss = ss[:-1]
    else:
      ss = uniIo.joinPath(dataUrl.url.path, '') # removes file separator from end, but not if just '/', hence special cases above
    filePath = fullFilePath[len(ss)+1:]
  
  else:
    dataUrl = None
    ll = []
    ss = fullFilePath
    for dummy in range(keepDirectories + 1):
      ss,name = os.path.split(ss)
      ll.append(name)
    ll.reverse()
    filePath = uniIo.joinPath(*ll)
  
  #
  return (dataUrl, filePath)

def changeDataStoreUrl(dataStore, newPath):
  """ Change the url for this dataStore, so at the end we have
  dataStore.dataUrl.url.path = newPath.  This changes all dataUrls
  with the same old path if the old path does not exist and the
  new one does.
  """

  newPath = uniIo.normalisePath(newPath, makeAbsolute=True)
  oldDataUrl = dataStore.dataUrl
  oldUrl = oldDataUrl.url
  oldPath = oldUrl.dataLocation
  oldExists = os.path.exists(oldPath)
  if newPath != oldPath:
    dataLocationStore = dataStore.dataLocationStore
    newUrl = Implementation.Url(path=newPath)  # TBD: should use oldUrl.clone(path=newPath)

    # first check if have a dataUrl with this path
    newDataUrl = dataLocationStore.findFirstDataUrl(url=newUrl)
    if not newDataUrl:
      # if old path exists and there is more than one dataStore with
      # this dataUrl then create new one
      dataUrlStores = dataLocationStore.findAllDataStores(dataUrl=oldDataUrl)
      if oldExists and len(dataUrlStores) > 1:
        newDataUrl = dataLocationStore.newDataUrl(name=oldDataUrl.name, url=newUrl)

    # if have found or have created newDataUrl then set dataStore to point to it
    # else just change url of oldDataUrl (which could affect other dataStores)
    if newDataUrl:
      dataStore.dataUrl = newDataUrl
    else:
      oldDataUrl.url = newUrl

    # if old path does not exist and new path exists then change urls of
    # all data urls which have old path to new path (there might be none)
    if not oldExists:
      newExists = os.path.exists(newPath)
      if newExists:
        for dataUrl in dataLocationStore.dataUrls:
          if dataUrl.url == oldUrl:
            dataUrl.url = newUrl


def _fetchUrlData(urlLocation):
  """Fetch data by inserting unverified certificate into urllib2
    getProxies not tested
  """
  import urllib2
  import ssl

  context = ssl._create_unverified_context()

  proxy_support = urllib2.ProxyHandler(urllib2.getproxies())
  opener = urllib2.build_opener(proxy_support)
  urllib2.install_opener(opener)

  request = urllib2.urlopen(urlLocation, context=context, timeout=3.0)

  return request

def _readDataFromRequest(request):
  """Read string from download from request
  """
  data = StringIO.StringIO(request.read())
  request.close()
  if data and data.buf:
    return data.buf


if __name__ == '__main__':

  # NOTE:ED - testing url loading from chemcomps

  import pandas as pd

  STRINGLEN = 300

  print "\n~~~~~~~~~~~~~~~~~~~~\n>>> raw github URLLIB2 request - from Path T"
  import urllib2
  import ssl

  context = ssl._create_unverified_context()
  # context = ssl.create_default_context(purpose=ssl.Purpose.SERVER_AUTH,
  #                                      cafile=None,
  #                                      capath=None)

  proxy_support = urllib2.ProxyHandler(urllib2.getproxies())
  opener = urllib2.build_opener(proxy_support)
  urllib2.install_opener(opener)

  try:
    request = urllib2.urlopen("https://raw.githubusercontent.com/VuisterLab/CcpNmr-ChemComps/master/data/pdbe/chemComp/archive/ChemComp/other/T/other%2BTho%2Bmsd_ccpnRef_2007-12-11-10-19-31_00014.xml",
                       context=context,
                       timeout=3.0)

  except Exception as es:
    print str(es)

  else:
    data = StringIO.StringIO(request.read())
    request.close()
    if data and data.buf:
      data = data.buf
      print data[:min(STRINGLEN, len(data))]

  ccpForgeIndexUrl = 'https://raw.githubusercontent.com/VuisterLab/CcpNmr-ChemComps/master/index/index.csv'
  try:
    request = fetchHttpResponse(ccpForgeIndexUrl)
    ccIndex = request.read()
    request.close()
    csvString = StringIO.StringIO(ccIndex)
    df = pd.read_csv(csvString)
    print "fetchHttpResponse - ", df[:1]

  except Exception as es:
    print(str(es))

  try:
    request = _fetchUrlData(ccpForgeIndexUrl)
    ccIndex = _readDataFromRequest(request)
    csvString = StringIO.StringIO(ccIndex)
    df = pd.read_csv(csvString)
    print "_fetchUrlData - ", df[:1]

  except Exception as es:
    print(str(es))
