#
# Get info for CASD-NMR projects
#

import os,re

from pdbe.adatah.Io import getReferenceTextFileFromHttp, getDataFromHttp

eNmrUrl = "http://www.wenmr.eu/wenmr/"
casdNmrDataUrl = os.path.join(eNmrUrl,"casd-nmr-data-sets")

from pdbe.adatah.Constants import archivesDataDir
casdNmrDataDir = os.path.join(archivesDataDir,'casdNmr')

class CasdNmrError(StandardError):
  pass

def getCasdNmrProjectInfo(casdNmrRefFile=None):

  """
  Code to get list of CASD-NMR projects info
  """

  hrefPatt = re.compile('\<a href\=\"([^\"]+)\"')
  hrefNamePatt = re.compile('[^ ]\"\>([^\>]+)\<\/a')
  pdbCodePatt = re.compile('\>\s*([A-Za-z0-9]{4})\s*\<\/a')
  
  # This file is customisable!
  if not casdNmrRefFile:
    casdNmrRefFile = os.path.join(casdNmrDataDir,'reference','dataPage.html')
    
  
  # Get the web page...
  dataLines = getReferenceTextFileFromHttp(casdNmrDataUrl,casdNmrRefFile,refText = "CASD-NMR data", isGzipped = False)
  
  # Now get the info out...
  projectInfo = []
  for dataLine in dataLines:

    if dataLine.count('href'):
    
      if dataLine.count("assignment-software"):
        continue
       
      # Some custum hacking here - content of href lines with data not dependable enough...
      if dataLine.count("rutgers") or dataLine.count('Data for') or dataLine.count('/wenmr/files/files/'):

        hrefSearch = hrefPatt.search(dataLine)
	dataLine = dataLine.replace("<span>","")
	dataLine = dataLine.replace("</span>","")
        projectNameSearch = hrefNamePatt.search(dataLine)

        if hrefSearch:
          urlName = projectNameSearch.group(1)
          if urlName.count("/wenmr"):
            urlName = urlName.replace("/wenmr",eNmrUrl)
          projectInfo.append([urlName,hrefSearch.group(1),[]])
          
      elif dataLine.count("structureId="):
        
        pdbCodeSearch = pdbCodePatt.search(dataLine)
        
        if pdbCodeSearch:
          projectInfo[-1][-1].append(pdbCodeSearch.group(1))
        
  if not projectInfo:
    raise CasdNmrError("No files found, probably web page change!")

  return projectInfo

def getCasdNmrProjects(saveDataDir = None, forceWrite = False):

  projectInfo = getCasdNmrProjectInfo()

  # Get the files
  
  if not saveDataDir:
    saveDataDir = casdNmrDataDir

  for (projectName,dataFile,pdbCodes) in projectInfo:
  
    if dataFile[0] == '/':
      dataUrl = os.path.join(eNmrUrl,dataFile[1:])
    else:
      dataUrl = dataFile
    
    (path,fileName) = os.path.split(dataFile)
    
    localDataFilePath = os.path.join(saveDataDir,fileName)
    
    if forceWrite or not os.path.exists(localDataFilePath):
      print "  Downloading CASD-NMR project %s..." % fileName

      dataLines = getDataFromHttp(dataUrl)
      
      if os.path.exists(localDataFilePath):
        os.remove(localDataFilePath)
      
      fout = open(localDataFilePath,'w')
      fout.write(dataLines)
      fout.close()
      
