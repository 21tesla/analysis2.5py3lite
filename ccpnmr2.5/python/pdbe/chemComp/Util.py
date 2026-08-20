import os, glob, re

from memops.api import Implementation

from pdbe.chemComp.Constants import freshTag

from pdbe.chemComp.Constants import editChemCompDataDir, testChemCompDataDir
from pdbe.chemComp.Constants import editChemCompCoordDataDir, testChemCompCoordDataDir

#
# Functions
#

def initialiseChemCompScript(sysArgs,getChemCompCoordDataDir=False):

  """
  
  Initialise a chemComp script
  
  Arguments:
  
  test      Will run in test mode - i.e. write files to the test/ directories
  all       If given, will do 'all' chemComps, not just the test set (not always relevant)
  verbose   Verbose mode
  write     Write out files (if not given, won't write but will just report what is going to happen)
  
  """
  
  #
  # Set parameters from script args
  #
  
  if 'test' in sysArgs:
    testMode = True
  else:
    testMode = False
    
  if 'all' in sysArgs:
    testChemCompMode = False
  else:
    testChemCompMode = True
    
  if 'verbose' in sysArgs:
    verbose = True
  else:
    verbose = False
  
  if 'write' in sysArgs:
    writeData = True
  else:
    writeData = False
    
  #
  # Set chemComp directory to use - NOTE THIS IS NOW A UNIQUE DIRECTORY!
  #  
  
  chemCompArchiveDataDir = editChemCompDataDir

  if testMode:
    chemCompArchiveDataDir = testChemCompDataDir
  
  #
  # Get a list of ccpCodes to handle. Always get this from the archive directory!
  #

  ccpCodeList = getCcpCodeList(editChemCompDataDir, testMode = testChemCompMode)
  
  #
  # Return chemCompCoord dir if required
  #
    
  if getChemCompCoordDataDir:
    returnArchiveDataDir = editChemCompCoordDataDir
    if testMode:
      returnArchiveDataDir = testChemCompCoordDataDir
      
  else:
    returnArchiveDataDir = chemCompArchiveDataDir
 
  return (testMode,writeData,verbose,ccpCodeList,returnArchiveDataDir)


def getCcpCodeList(storageDir, testMode = False, sourceName = ""):

  """
  
  Get a list of CCP codes to handle. This is set up to work on the archive directories with their
  specific structure!!!
  
  Use sourceName to get at coordinate data files (e.g. pdb, ideal, ...)
  
  """
  
  ccpCodePatt = re.compile("ccpCode=\"([^\"]+)\"")
  
  if testMode:
  
    ccpCodeList = (('protein',('Arg',)),('DNA',('G',)),('RNA',('G',)),('other',('001',))) #,('carbohydrate',('Gal',)))
  
  else:
  
    molTypes = ['protein','DNA','RNA','other','carbohydrate']
    ccpCodeList = []


    for molType in molTypes:

      ccpCodeList.append([molType,[]])
      
      archiveStorageDir = os.path.join(storageDir,sourceName,molType)
      
      if not os.path.exists(archiveStorageDir):
        print "  Warning: directory %s does not exist!" % archiveStorageDir
        continue
      
      if sourceName:
        sourceNameText = "%s+" % sourceName
      else:
        sourceNameText = ""

      chemCompFileSearchString = "%s%s+*.xml" % (sourceNameText,molType)


      if molType == 'other':
        
        # These are organised differently
        
        subDirs = os.listdir(archiveStorageDir)
        chemCompFileNameMatches = []
        
        for subDir in subDirs:
          chemCompFileNameMatches.extend(glob.glob(os.path.join(archiveStorageDir,subDir,chemCompFileSearchString)))
      
      else:
        chemCompFileNameMatches = glob.glob(os.path.join(archiveStorageDir,chemCompFileSearchString))

      if not chemCompFileNameMatches:
        continue


      for ccpFile in chemCompFileNameMatches:

        if ccpFile[-4:] != '.xml':
          continue

        (pathName,baseName) = os.path.split(ccpFile)

        baseNameItems = baseName.split('+')
        
        if sourceName:
          molType = baseNameItems[1]
          ccpCode = baseNameItems[2]
        else:
          molType = baseNameItems[0]
          ccpCode = baseNameItems[1]
       
        # Get 'real' ccp code if filename contains underscores!
        if ccpCode.count("_"):
          ccpCode = None
          fin = open(ccpFile)
          for line in fin.readlines():
            ccpCodeSearch = ccpCodePatt.search(line)
            if ccpCodeSearch:
              ccpCode = ccpCodeSearch.group(1)
              break
              
          if not ccpCode:
            print "  Error: no ccpCode for file with name %s... ignored!" % ccpFile
          
        # Exceptions - these are handled only on CCPN side
        if molType in ['DNA','RNA'] and ccpCode == 'X':
          continue
        elif molType == 'protein' and ccpCode == 'Xxx':
          continue

        ccpCodeList[-1][-1].append(ccpCode)
      
      ccpCodeList[-1][-1].sort()
      
  return ccpCodeList

def getGuidFromXmlFileName(xmlFileName):

  components = xmlFileName[:-4].split('+')
  
  return components[-1]
  
def tagAsFresh(ccpnObject):

  appData = ccpnObject.findFirstApplicationData(keyword = freshTag)
  
  if not appData:
    appData = Implementation.AppDataBoolean(application = 'pdbeCcGeneration', keyword = freshTag, value = True)
    ccpnObject.addApplicationData(appData)
  
  appData.value = True

def removeFreshTag(ccpnObject):

  appData = ccpnObject.findFirstApplicationData(keyword = freshTag)
  
  if appData:
    ccpnObject.removeApplicationData(appData)
