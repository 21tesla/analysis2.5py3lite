import os, shutil, sys, filecmp, glob

from memops.general.Io import getCcpFileString, getTopDirectory
from pdbe.chemComp.Util import getCcpCodeList


from pdbe.chemComp.Constants import editChemCompDataDir, editChemCompCoordDataDir, standardRepChemCompDataDir

if __name__ == '__main__':

  #
  # Use this script to update the chemComp XML files in the CCPN CVS directories.
  #
  # Copies DIRECT from the editCHemCompDir and the chemCompCoordDir - CHECK IF LICENSES ARE OK!!
  #
  
  customCcpCodeList = None
  
  
  #customCcpCodeList = (
  
  #  ('DNA', ['U','I','X','A00','C00','G00','I00','T00','U00','A11','C11','G11','I11','T11','U11']),
  #  ('RNA', ['T','I','X','A00','C00','G00','I00','T00','U00','A11','C11','G11','I11','T11','U11'])
    
  #  )

  #print "CHECK FIRST"
  #sys.exit()
  
  # TODO: need decent args here that can be set by calling script!
  doFileChanges = False
  
  if not doFileChanges:
    print "  WARNING: not making any changes!"
  # TODO BELOW clearly has to come from some central Constants information!!!
  targetDirectory = os.path.join(getTopDirectory(),'..','..','stable','ccpn') #

  # WARNING WARNING: have to use specific code to get files out - are organised differently!
  #editChemCompDataDir = getChemCompArchiveDataDir()
  targetChemCompDataDir = os.path.join(targetDirectory,'data/ccp/molecule/ChemComp')

  #editChemCompCoordDataDir = getChemCompCoordArchiveDataDir()
  targetChemCompCoordDataDir = os.path.join(targetDirectory,'data/ccp/molecule/ChemCompCoord')
    
  copyList = [(editChemCompDataDir,targetChemCompDataDir,''),
              (editChemCompCoordDataDir,targetChemCompCoordDataDir,'pdb'),
              (editChemCompCoordDataDir,targetChemCompCoordDataDir,'ideal'),
              (editChemCompCoordDataDir,targetChemCompCoordDataDir,'euroCarbDb')]

  for (sourceDir,targetDir,prefix) in copyList:
      
    if not os.path.exists(targetDir):
      print "Error: target directory %s does not exist - aborting." % targetDir
      sys.exit()
      
    print "Examining source %s to target %s..." % (sourceDir,targetDir)

    if not customCcpCodeList:
      targetCcpCodeList = getCcpCodeList(targetDir,prefix = prefix)
    else:
      targetCcpCodeList = customCcpCodeList

    if prefix:
      prefixText = "%s+" % prefix
    else:
      prefixText = ""

    for (molType,ccpCodes) in targetCcpCodeList:
      
      for ccpCode in ccpCodes:          

        chemCompFileSearchString = "%s%s+%s+*.xml" % (prefixText,molType,getCcpFileString(ccpCode))
        chemCompFileNameMatches = glob.glob(os.path.join(targetDir,chemCompFileSearchString))        

        if molType == 'other':
          sourceSubDir = os.path.join(molType,ccpCode[0])
        else:
          sourceSubDir = molType
          
        if customCcpCodeList and not chemCompFileNameMatches:
          # See if new file, copy over directly in this case          
          chemCompFileNameMatches = glob.glob(os.path.join(sourceDir,sourceSubDir,chemCompFileSearchString))
    
        if chemCompFileNameMatches:
        
          chemCompFilePath = chemCompFileNameMatches[-1]
          (chemCompFileDir,chemCompFileName) = os.path.split(chemCompFilePath)

          targetFile = os.path.join(targetDir,chemCompFileName)
          sourceFile = os.path.join(sourceDir,sourceSubDir,chemCompFileName)
          
          deleteOriginalTarget = False
          # If doesn't exist, try look for a source file...
          if not os.path.exists(sourceFile):
            chemCompFileNameMatches = glob.glob(os.path.join(sourceDir,sourceSubDir,chemCompFileSearchString))
            if chemCompFileNameMatches:
              sourceFile = chemCompFileNameMatches[-1]
              (sourceFileDir,sourceFileName) = os.path.split(sourceFile)
                            
              deleteOriginalTarget = True
          
          if os.path.exists(sourceFile):
            doCopy = False
            if os.path.exists(targetFile):
            
              # TODO this does not work well any more, files are now different if re-saved even
              # if information is exactly the same. Keep for now though.
              if not filecmp.cmp(targetFile,sourceFile):
                if not deleteOriginalTarget:
                  print "Overwriting %s" % targetFile
                doCopy = True
            else:
              doCopy = True
            
            if doCopy:
            
              if deleteOriginalTarget:
                print "Using new name %s instead of %s, removing original..." % (sourceFileName,chemCompFileName)
                if doFileChanges:
                  os.remove(targetFile)
                targetFile = os.path.join(targetDir,sourceFileName)
              else:
                print "Copying %s..." % (chemCompFileName)
            
              if doFileChanges:
                shutil.copy(sourceFile,targetFile)
              
              print
            
          else:
            print " NO %s" % sourceFile
