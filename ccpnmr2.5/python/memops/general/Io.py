
"""
======================COPYRIGHT/LICENSE START==========================

Io.py: code for CCPN data model and code generation framework

Copyright (C) 2007  (CCPN Project)

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

- email: ccpn@bioc.cam.ac.uk

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

"""
#
# Convenient I/O functions
#

import gc
import os
import shutil
import time
import glob

from memops.metamodel import ImpConstants

from memops.universal.Constants import dirsep
from memops.universal.Io import normalisePath, joinPath, getTopDirectory
from memops.universal.Io import removePath, doConvertStringToFileName
from memops.universal.MessageReporter import showWarning as printWarning
from memops.universal.Util import isWindowsOS

from memops.format.xml import Util as xmlUtil

from memops.general.Implementation import ApiError

from memops.api import Implementation

def getDataDirectory():

  directory = joinPath(getTopDirectory(), 'data')

  return directory

def checkFileAtPath(path):
  """Check that file on disk ends correctly
  """
  if not os.path.exists(path):
    return False

  size = os.path.getsize(path)
  fp = open(path, 'rU')

  ss = '<!--End of Memops Data-->'
  n = len(ss)
  fp.seek(size-n-2, 0)  # -2 just to be safe in case has \r\n as line ending
  data = fp.read().strip()[-n:]
  fp.close()

  return data == ss

def checkFileIntegrity(topObject):
  """Check that topObject file on disk ends correctly
  """
  path = xmlUtil.getTopObjectPath(topObject)

  return checkFileAtPath(path)

def createTopObjectFallback(topObject):
  """
  Create backup of topObject in same directory as original file but with '.bak' appended.
  This function is not intended to be used outside this module but could be.
  """

  location = xmlUtil.getTopObjectPath(topObject)
  if not os.path.exists(location):
    return

  backupLocation = location + '.bak'
  if not checkFileAtPath(location) and checkFileAtPath(backupLocation):
    # current file no good and current backup good so do not do backup
    print 'File at location "%s" not complete so not backing up' % location
    return

  # copy rather than move because will need that much disk space in any case
  # and sometimes move fails at OS level if file with that name already exists
  directory = os.path.dirname(backupLocation)
  if not os.path.exists(directory):
    os.makedirs(directory)
  shutil.copy(location, backupLocation)

def checkRemoveIfExists(path, removeExisting = False, showYesNo = None):
  """
  If path already exists:
    If removeExisting:
      Delete the path
      Return True
    Else if showYesNo:
      Ask the user if it is ok to delete the path
      If yes, delete and return True.  If no return None.
    Else:
      Raise an IOError
  Else:
    Return False
  This function is not intended to be used outside this module but could be.
  """

  if os.path.exists(path):

    if removeExisting:
      removePath(path)
      return True

    elif showYesNo:
      if os.path.isdir(path):
        files = os.listdir(path)
        n = 5
        if len(files) > n:
          files = files[:n]
          files.append('...')
        ss = ', '.join(files)
        message = '%s already exists, remove it and all of its contents (%s) (if no, then no action taken)?' % (path, ss)
        if not showYesNo('Remove directory', message):
          return None
        else:
          removePath(path)
          return True
      else:
        message = '%s already exists (not as a directory), remove it?' % path
        if not showYesNo('Remove file', message):
          return None
        else:
          removePath(path)
          return True
    else:
      raise IOError('%s already exists' % path)

  return False

def saveProject(project, newPath = None, newProjectName = None, changeBackup = True,
                createFallback = False, removeExisting = False, showYesNo = None,
                checkValid = False, changeDataLocations = False, showWarning = None,
                revertToOriginal = False):
  """
  Save the userData for a project to a location given by newPath (the url.path
  of the userData repository) if set, or the existing location if not. 
  If userData does not exist then throws IOError.
  If newPath is not specified then it is set to oldPath.
  If newProjectName is not specified then it is set to oldProjectName if
    newPath==oldPath, otherwise it is set to basename(newPath).
  If changeBackup, then also changes backup URL path for project.
  If createFallback, then makes copy of existing modified topObjects
    files (in newPath, not oldPath) before doing save.
  If newPath != oldPath and newPath exists (either as file or as directory):
    If removeExisting:
      Delete the newPath.
    Else if showYesNo:
      Ask the user if it is ok to delete the newPath
      If yes, delete.  If no, return without saving.
    Else:
      Raise an IOError
  Elif newProjectName != oldProjectName and there exists corresponding path (file/directory):
    If removeExisting:
      Delete the path.
    Else if showYesNo:
      Ask the user if it is ok to delete the path.
      If yes, delete.  If no, return without saving.
    Else:
      Raise an IOError
  If checkValid then does checkAllValid on project
  If changeDataLocations then copy to project directory
  If revertToOriginal then after the save changes paths back to as they were originally
  If there is no exception or early return then at end userData is pointing to newPath.
  Return True if save done, False if not (unless there is an exception)
  """

  # check project valid (so don't save an obviously invalid project)
  if checkValid:
    project.checkAllValid()

  # only want to change path for userData
  userData = project.findFirstRepository(name='userData')
  if not userData:
    raise IOError('Problem: userData not found')

  oldPath = userData.url.path
  oldProjectName = project.name

  if newPath:
    # normalise newPath
    newPath = normalisePath(newPath, makeAbsolute=True)
  else:
    newPath = oldPath

  # if newProjectName isn't specified then use default
  if not newProjectName:
    if newPath == oldPath:
      newProjectName = oldProjectName
    else:
      newProjectName = os.path.basename(newPath)
      # below is because of data model limit
      newProjectName = newProjectName[:32]

  renameProject(project, newProjectName)

  # if newPath same as oldPath, check if newProjectName already exists if it's not same as oldProjectName
  if newPath == oldPath:
    if newProjectName != oldProjectName:
      location = xmlUtil.getTopObjectPath(project)
      if os.path.exists(location):
        answer = checkRemoveIfExists(location, removeExisting, showYesNo)
        if answer is None:
          project.__dict__['name'] = oldProjectName  # TBD: for now name is frozen so change this way
          return False
  else: # check instead if newPath already exists
    if os.path.exists(newPath):
      answer = checkRemoveIfExists(newPath, removeExisting, showYesNo)
      if answer is None:
        if newProjectName != oldProjectName:
          project.__dict__['name'] = oldProjectName  # TBD: for now name is frozen so change this way
        return False
      # TBD: should we be removing it?
      removePath(newPath)
    else:
      # NBNB 2008/04/03 Rasmus Fogh. Added because it otherwise fell over when
      # target path did not exist
      upDir = os.path.dirname(newPath)
      if not os.path.exists(upDir):
        os.makedirs(upDir)

    # check if any topObject activeRepository is not either of above
    refData = project.findFirstRepository(name='refData')
    genData = project.findFirstRepository(name='generalData')
    topObjects = []
    repositories = set()
    for topObject in project.topObjects:
      repository = topObject.findFirstActiveRepository()
      if repository and repository not in (userData, refData, genData):
        topObjects.append(topObject)
        repositories.add(repository)
    if topObjects:
      print('Warning, topObjects %s, in repositories %s, being left in original locations' % (topObjects, repositories))

  try:
    oldUrl = userData.url
    if changeBackup:
      # change project backup url to point to new path
      backupRepository = project.findFirstRepository(name="backup")
      if backupRepository:
        oldBackupUrl = backupRepository.url
      else:
        changeBackup = False

    # copy userData files over
    if newPath != oldPath:
      # if os.path.exists(oldPath):  # only copy if this is a directory
      if os.path.isdir(oldPath):
        # just copy everything from oldPath to newPath
        print('Copying directory %s to %s (this might take some time if there are big files)' % (oldPath, newPath))
        shutil.copytree(oldPath, newPath)

        # but need toz remove all implementation files
        implPath = joinPath(newPath, ImpConstants.modellingPackageName,
                            ImpConstants.implementationPackageName)
        #implPath = pathImplDirectory(newPath)
        removePath(implPath)

        # and need to repoint dataUrl's that were copied over
        oldPathP = oldPath + '/'
        for dataLocationStore in project.dataLocationStores:
          for dataStore in dataLocationStore.dataStores:
            oldDataPath = dataStore.fullPath
            if oldDataPath.startswith(oldPathP):
              dataUrl = dataStore.dataUrl
              oldDataUrlPath = dataUrl.url.dataLocation
              if oldDataUrlPath.startswith(oldPathP): # normally true
                newDataUrlPath = newPath + oldDataUrlPath[len(oldPath):]
                dataUrl.url = Implementation.Url(path=newDataUrlPath)
              else: # path split awkwardly between absolute and relative
                newDataUrlPath = newPath
                dataUrl.url = Implementation.Url(path=newDataUrlPath)
                dataStore.path = oldDataPath[len(oldPath):]

      # change userData url to point to new path
      userData.url = Implementation.Url(path=newPath)
      # above will set project.isModified = True

      if changeBackup:
        # change project backup repository url to point to new path
        backupRepository.url = Implementation.Url(path=newPath+'_backup')

    # change project name
    if newProjectName != oldProjectName:
      if not project.isModified: # if it isModified it will be saved below
        if createFallback:
          createTopObjectFallback(project)
        project.save()

    # create fallbacks and keep track of modified topObjects
    modifiedTopObjects = []
    if createFallback:
      for topObject in (project,)+tuple(project.topObjects):
        if not topObject.isDeleted and topObject.isModified:
          createTopObjectFallback(topObject)
          modifiedTopObjects.append(topObject)

    if changeDataLocations:
      dataLocationStores = project.sortedDataLocationStores()

      userRepository = project.findFirstRepository(name='userData')
      userPath = userRepository.url.dataLocation
      # 2010 Aug 11: remove data directory from path
      #dataPath = joinPath(userPath, 'data')
      dataPath = userPath

      # 2010 Aug 11: change name
      #dataStorePrefix = 'dataStore'
      dataStorePrefix = 'spectra'
      if os.path.exists(dataPath):
        files = [xx for xx in os.listdir(dataPath) if xx.startswith(dataStorePrefix)]
        offset = len(files)
      else:
        offset = 0

      copyingList = []
      dataUrlDict = {}
      for dataLocationStore in dataLocationStores:
        for dataStore in dataLocationStore.sortedDataStores():

          # wb104: 24 Mar 2010: below check is a complete kludge
          # we should check whether dataStore is instance of
          # NumericMatrix, etc., but those are in ccp so should
          # not be imported here
          # in any case, there is no proper way to find out if
          # a dataStore is used without explicit knowledge of class
          knt = 0
          # hicard = 1
          for attr in ('nmrDataSourceImage',):
            if hasattr(dataStore, attr):
              if getattr(dataStore, attr):
                knt += 1
          # hicard > 1
          for attr in ('externalDatas', 'nmrDataSources'):
            if hasattr(dataStore, attr):
              knt += len(getattr(dataStore, attr))
          if knt == 0:
            continue

          oldFullPath = dataStore.fullPath
          if not oldFullPath.startswith(userPath+'/'):
            # first figure out new dataUrl path
            dataUrl = dataStore.dataUrl
            oldPath = dataUrl.url.dataLocation
            if dataUrl in dataUrlDict:
              newUrl = dataUrlDict[dataUrl]
            else:
              offset += 1
              newUrlPath = '%s%d' % (dataStorePrefix, offset)
              newUrlPath = joinPath(dataPath, newUrlPath)
              newUrl = dataUrlDict[dataUrl] = Implementation.Url(path=newUrlPath)
            # then add to list to copy over if original data exists
            if os.path.exists(oldFullPath):
              newFullPath = joinPath(newUrl.dataLocation, dataStore.path)
              copyingList.append((oldFullPath, newFullPath))

      # now copy data files over
      nfilesToCopy = len(copyingList)
      for n, (oldFullPath, newFullPath) in enumerate(copyingList):
        dirName = os.path.dirname(newFullPath)
        if not os.path.exists(dirName):
          os.makedirs(dirName)
        print('Copying file %s to %s (%d of %d)' % (oldFullPath, newFullPath, n+1, nfilesToCopy))
        shutil.copy(oldFullPath, newFullPath)

      # finally change dataUrl paths
      for dataUrl in dataUrlDict:
        dataUrl.url = dataUrlDict[dataUrl]

    # save modifications
    # change way doing save in case exception is thrown
    if createFallback:
      for topObject in modifiedTopObjects:
        try:
          topObject.save()
        except:
          location = xmlUtil.getTopObjectPath(topObject)
          print('Exception working on topObject %s, file %s' % (topObject, location))
          raise
      # be safe and do below in case new modifications after
      # modifiedTopObjects has been created
      project.saveModified()
    else:
      project.saveModified()

    if not isWindowsOS():
      os.system('touch "%s"' % newPath)  # so that user can see which are most recent

    badTopObjects = []
    for topObject in modifiedTopObjects:
      if not checkFileIntegrity(topObject):
        badTopObjects.append(topObject)

    if badTopObjects:
      if showWarning:
        showWarning('Incomplete save', 'It looks like one or more files did not save completely, see console for list')
      print('It looks like one or more files did not save completely, you should check them:')
      for topObject in badTopObjects:
        print
        print('%s, path:' % topObject)
        print(xmlUtil.getTopObjectPath(topObject))
      return False
        
    if revertToOriginal: # TBD: the below does not change back dataUrl paths
      if newProjectName != oldProjectName:
        project.__dict__['name'] = oldProjectName  # TBD: for now name is frozen so change this way

      if newPath != oldPath:
        userData.url = oldUrl
        if changeBackup:
          backupRepository.url = oldBackupUrl
      
    return True

  except:
    # saveModified failed so revert to old values
    if newProjectName != oldProjectName:
      project.__dict__['name'] = oldProjectName  # TBD: for now name is frozen so change this way

    if newPath != oldPath:
      userData.url = oldUrl
      if changeBackup:
        backupRepository.url = oldBackupUrl
      try:
        removePath(newPath)
      except:
        pass

    raise

def renameProject(project, newProjectName):
  """ Rename project.
  """
  
  # change project name
  if newProjectName == project.name:
    return
  
  else:
    print '### renaming', project.name, newProjectName
  
    project.override = True # TBD: for now name is frozen so change this way
    try:
      # below constraint is not checked in setName() if override is True so repeat here
      isValid = newProjectName.isalnum()  # superfluous but faster in most cases
      if not isValid:
        for cc in newProjectName:
          if cc != '_' and not cc.isalnum():
            isValid = False
            break
        else:
          isValid = True
      if (not (isValid)):
        raise ApiError('project name must only have characters that are alphanumeric or underscore')

      # below checks for length of name as well
      project.name = newProjectName
    finally:
      project.override = False


def saveRepository(repository, newPath = None, createFallback = False,
                   removeExisting = False, showYesNo = None):
  """
  Save the repository to a location given by newPath (the url.path
  of the repository) if set, or the existing location if not. 
  For the userData repository this just calls saveProject instead.
  If createFallback, then makes copy of existing modified repository
    files (in newPath, not oldPath) before doing save.
  If newPath != oldPath and newPath exists (either as file or as directory):
    If removeExisting:
      Delete the newPath.
    Else if showYesNo:
      Ask the user if it is ok to delete the newPath
      If yes, delete.  If no, return without saving.
    Else:
      Raise an IOError
  If there is no exception or early return then at end repository is pointing to newPath.
  Return True if save done, False if not (unless there is an exception)
  """

  project = repository.root
  if repository.name == 'userData':
    saveProject(project, newPath=newPath, createFallback=createFallback,
                removeExisting=removeExisting, showYesNo=showYesNo)
    return False

  oldPath = repository.url.path

  if newPath:
    # normalise path
    newPath = normalisePath(newPath, makeAbsolute=True)
  else:
    newPath = oldPath

  if newPath != oldPath:
    answer = checkRemoveIfExists(newPath, removeExisting, showYesNo)
    if answer is None:
      return False

  # check which topObjects have this repository as activeRepositories[0]
  # and are modified but not deleted
  topObjects = set()
  for topObject in project.topObjects:
    if topObject.isModified and not topObject.isDeleted and \
           repository == topObject.findFirstActiveRepository():
      topObjects.add(topObject)

  # TBD: should this be done here?
  for topObject in topObjects:
    topObject.checkAllValid()

  # copy repository files over
  if newPath != oldPath:
    if os.path.exists(oldPath):
      # just copy everything from oldPath to newPath
      shutil.copytree(oldPath, newPath)

    # change repository url to point to new newPath
    oldUrl = repository.url
    repository.url = Implementation.Url(path=newPath)
    # above will set project.isModified = True

  # save modifications
  try:
    if createFallback:
      for topObject in topObjects:
        createTopObjectFallback(topObject)

      for topObject in topObjects:
        topObject.save()

    # the paths have changed so also need to save project
    project.save()

    return True

  except:
    # something failed so revert to old path
    if newPath != oldPath:
      repository.url = oldUrl
      try:
        removePath(newPath)
      except:
        pass

    raise

def newProject(projectName, path = None, removeExisting = False, showYesNo = None):
  """
  Create, and return, a new project using a specified path (directory).
  If path is not specified it takes the current working directory.
  The path can be either absolute or relative.
  The 'userData' repository is pointed to the path.
  The 'backup' repository is pointed to the path + '_backup'.
  If either of these paths already exist (either as files or as directories):
    If removeExisting:
      Delete the path
    Else if showYesNo:
      Ask the user if it is ok to delete the path
      If yes, delete.  If no return None.
    Else:
      Raise an IOError
  """

  # relies on knowing that repositories to move have these names, and these values for path suffix
  repositoryNameMap = {'userData': '', 'backup': '_backup'}

  if not path:
    path = os.getcwd()

  for name in repositoryNameMap.keys():
    fullPath = joinPath(path, projectName) + repositoryNameMap[name]
    answer = checkRemoveIfExists(fullPath, removeExisting, showYesNo)
    if answer is None:
      return None

  project = Implementation.MemopsRoot(name=projectName)

  for name in repositoryNameMap.keys():
    fullPath = normalisePath(joinPath(path, projectName) 
                             + repositoryNameMap[name], makeAbsolute=True)
    repository = project.findFirstRepository(name=name)
    repository.url = Implementation.Url(path=fullPath)

  return project

def loadProject(path, projectName=None, showWarning=None, askFile=None, 
                askDir=None, suppressGeneralDataDir = False):
  """
  Loads a project file and checks and deletes unwanted project repositories and
  changes the project repository path if the project has moved.  Returns the project.
  (The project repository path is effectively the userData repository.)
  showWarning (if not None) has signature showWarning(title, message)
  askFile (if not None) has signature askFile(title, message, initial_value = '')
  askDir (if not None) has signature askDir(title, message, initial_value = '')
  Throws an IOError if there is an I/O error.
  Throws an ApiError if there is an API exception.
  """
  
  from memops.format.xml import XmlIO

  if not showWarning:
    showWarning = printWarning

  path = normalisePath(path, makeAbsolute=True)

  # check if path exists and is directory
  if not os.path.exists(path):
    raise IOError('path "%s" does not exist' % path)
  if not os.path.isdir(path):
    raise IOError('path "%s" is not a directory' % path)

  projectFile = xmlUtil.getProjectFile(path, projectName)

  if projectName:
    # projectName was specified so projectFile better exist
    if not os.path.exists(projectFile):
      raise IOError('project file "%s" does not exist' % projectFile)
  else:
    # projectName was not specified so default projectFile might not exist
    if not os.path.exists(projectFile):
      projectFiles = xmlUtil.getPossibleProjectFiles(path)

      if len(projectFiles) == 0:
        raise IOError('"%s" contains no project file' % path)
      elif len(projectFiles) == 1:
        projectFile = projectFiles[0]
      elif askFile:
        projectFile = askFile('Select project file', 'Select project file', initial_value=projectFiles[0])
        if not projectFile: # cancelled
          raise IOError('Cancelled')
      else:
        raise IOError('"%s" contains %d project files, not sure which to use' % (path, len(projectFiles)))

    # TBD: should projectName be based on projectFile or on path???
    # the way it is set here do not need to change project.name
    # but if you used the path then you would need to change it
    projectName = xmlUtil.getTopObjIdFromFileName(projectFile)

  # doing the loadProject twice has a bit of an overhead, but not much

  try:
    # below assumes TopObjects exist where stated, so might fail
    project = XmlIO.loadProject(path, projectName, partialLoad=False)

  except:
    print("\nFirst loading attempt failed - compatibility problem?.")
    project = None
  
  def isGeneralDataWriteable(generalDataRepository):
    oldPath = generalDataRepository.url.path
    return (isWindowsOS() or os.access(oldPath, os.W_OK|os.X_OK))

  def isGeneralDataOk(project):
    if suppressGeneralDataDir:
      return True
    generalDataRepository = project.findFirstRepository(name='generalData')
    if generalDataRepository:
      return isGeneralDataWriteable(generalDataRepository)
    return True

  if project is not None and (not xmlUtil.areAllTopObjectsPresent(project) or
      not isGeneralDataOk(project)):
    # if not all loaded (shell) TopObjects can be found, try again
    project = None
    ###print "\nSome files unfindable - has project moved?."
  
  if project is None:
    ###print "Re-trying, skipping cached TopObjects:"
    project = XmlIO.loadProject(path, projectName, partialLoad=True)

  warningMessages = []

  # try and fix project repository path, if needed
  packageLocator = project.packageLocator
  repositories = packageLocator.repositories
  for repository in repositories:
    if repository.url.path == path:
      oldPath = path
      break
  else:
    # change first repository path to point to this one, and also change backup
    repository = repositories[0]
    oldPath = repository.url.path
    warningMessages.append('Project file has moved from\n"%s"\nto\n"%s"' % (oldPath, path))
    repository.url = Implementation.Url(path=path)
    # Necessary because activeRepositories are not set right 
    # if file names do not match:
    project.__dict__['activeRepositories'].append(repository)

  projectRepository = repository

  # check backup path
  backupRepository = project.findFirstRepository(name="backup")
  if backupRepository:
    backupUrl = backupRepository.url
    oldBackupPath = backupUrl.path
    newBackupPath = path + '_backup'
    if oldBackupPath.startswith(oldPath): # hopefully true
      if path != oldPath:
        warningMessages.append('Backup is being changed from\n"%s"\nto\n"%s"' % (oldBackupPath, newBackupPath))
        backupRepository.url = Implementation.Url(path=newBackupPath)
    elif not os.path.exists(oldBackupPath):
      warningMessages.append('Backup is being changed from\n"%s"\nto\n"%s"' % (oldBackupPath, newBackupPath))
      backupRepository.url = Implementation.Url(path=newBackupPath)

  # check if project repository is called 'userData'
  if projectRepository.name != 'userData':
    warningMessages.append('Project has non-standard repository name "%s"' % projectRepository.name)

  # repoint dataStores that are in same directory as project
  # (but only if old path does not exist and new one does)
  if path != oldPath:
    oldDirectory = os.path.dirname(oldPath)
    newDirectory = os.path.dirname(path)
    for dataLocationStore in project.dataLocationStores:
      for dataStore in dataLocationStore.dataStores:
        fullPath = dataStore.fullPath
        if not os.path.exists(fullPath):
          dataUrl = dataStore.dataUrl
          dataLocation = dataUrl.url.dataLocation
          if dataLocation == oldDirectory or \
              (dataLocation.startswith(oldDirectory) and dataLocation[len(oldDirectory)] == dirsep):
            newDataUrlPath = newDirectory + dataLocation[len(oldDirectory):]
            newPath = joinPath(newDataUrlPath, dataStore.path)
            if os.path.exists(newPath):
              warningMessages.append('DataStore %s:%s path has been changed from\n"%s"\nto\n"%s"' % (dataStore.dataLocationStore.name, dataStore.serial, oldDirectory, newDataUrlPath))
              dataUrl.url = dataUrl.url.clone(path=newDataUrlPath)

  # change refData to current one if need be
  refDataRepository = project.findFirstRepository(name='refData')
  if refDataRepository:
    oldPath = refDataRepository.url.path
    newPath = normalisePath(getDataDirectory())
    if newPath != oldPath:
      warningMessages.append('refData has been changed from\n"%s"\nto\n"%s"' % (oldPath, newPath))
      refDataRepository.url = Implementation.Url(path=newPath)
  else:
    warningMessages.append('Project has no repository with name "refData"')

  # change generalData to current one if need be
  generalDataRepository = project.findFirstRepository(name='generalData')
  if generalDataRepository and not suppressGeneralDataDir:
    oldPath = generalDataRepository.url.path
    if not os.path.exists(oldPath) or not isGeneralDataWriteable(generalDataRepository):
      newPath = normalisePath(os.path.expanduser('~/.ccpn/data'))
      if not os.path.exists(newPath):
        os.makedirs(newPath)
      warningMessages.append('generalData has been changed from\n"%s"\nto\n"%s"' % (oldPath, newPath))
      generalDataRepository.url = Implementation.Url(path=newPath)

  # check other repository paths
  for repository in project.repositories:
    if repository not in (projectRepository, refDataRepository, backupRepository, generalDataRepository):
      oldPath = repository.url.path
      if not repository.stored:
        title = 'Repository being deleted'
        msg = 'Repository "%s" with path "%s" has no packageLocators, deleting' % (repository.name, oldPath)
        repository.delete()
      elif not os.path.exists(oldPath):
        title = 'Repository path does not exist'
        if showWarning is printWarning:
          tt = ''
        else:
          tt = ' (see console for packages)'
        msg = 'Repository "%s" path "%s" does not exist%s' % (repository.name, oldPath, tt)

        print('List of packageLocators for repository "%s":' % repository.name)
        for packageLocator in repository.stored:
          print('  %s' % packageLocator.targetName)

        if askDir:
          newPath = askDir(title, msg + ': enter new path', initial_value=oldPath)
          while newPath and not os.path.exists(newPath):
            msg = 'Path "%s" does not exist' % newPath
            newPath = askDir(title, msg + ': enter new path', initial_value=newPath)
          if newPath:
            repository.url = Implementation.Url(path=normalisePath(newPath))
        else:
          warningMessages.append(msg)

  if warningMessages:
    showWarning('Warnings', '\n\n'.join(warningMessages))

  return project

def backupProject(project, dataLocationStores = None, skipRefData = True, clearOutDir = False):

  def modificationTime(path):
    return os.stat(path)[8]

  backupRepository = project.findFirstRepository(name="backup")

  if not backupRepository:
    print('Warning: no backup path set, so no backup done')
    return

  backupUrl = backupRepository.url
  backupPath = backupUrl.path

  if not dataLocationStores:
    dataLocationStores = set()

  if clearOutDir:
    removePath(backupPath)

  topObjects = tuple(project.topObjects) + (project,)

  for topObject in topObjects:
    if skipRefData:
      repository = topObject.findFirstActiveRepository(name='refData')
      if repository:
        continue

    if topObject.isModified:
      topObject.backup()
    else:
      repository = topObject.findFirstActiveRepository()
      if repository:
        origFile = xmlUtil.findTopObjectPath(repository.url.path, topObject)
        if os.path.exists(origFile):
          # problem with appending repository.name is that topObject.backup()
          # above does not do it this way, so end up with inconsistent backup
          ###backupDir = joinPath(backupPath, repository.name)
          # so use same backup pah as topObject.backup()
          backupDir = backupPath
          backupFile = xmlUtil.findTopObjectPath(backupDir, topObject)
          if not os.path.exists(backupFile) or \
              (modificationTime(backupFile) < modificationTime(origFile)):
            directory = os.path.dirname(backupFile)
            if not os.path.exists(directory):
              os.makedirs(directory)
            shutil.copy(origFile, backupFile)
        else:
          # one is stuffed
          print('Warning: could not backup %s since could not find original file "%s"' % (topObject, origFile))
      else:
        # one is stuffed
        print('Warning: could not backup %s since could not find repository' % topObject)

  dataBackupPath = joinPath(backupPath, 'data')
  for dataLocationStore in dataLocationStores:
    dataBackupDir = joinPath(dataBackupPath, dataLocationStore.name)
    for dataStore in dataLocationStore.dataStores:
      origFile = dataStore.dataUrl.url.path
      backupFile = joinPath(dataBackupDir, dataStore.path)

      if os.path.exists(origFile):
        if not os.path.exists(backupFile) or \
            (modificationTime(backupFile) < modificationTime(origFile)):
          directory = os.path.dirname(backupFile)
          if not os.path.exists(directory):
            os.makedirs(directory)
          shutil.copy(origFile, backupFile)
      else:
        print('Warning: could not backup dataStore "%s" because could not find original file "%s"' % (dataStore.name, origFile))

def modifyPackageLocators(project,repositoryName,repositoryPath,packageNames,resetPackageLocator = True,resetRepository = False):
  
  """
  Resets package locators for specified packages to specified repository.
  
  Use as, for example:
  
  modifyPackageLocators(project,'newChemComps','/mydir/data/chemComps/',('ccp.molecule.ChemComp','ccp.molecule.ChemCompCoord'))
  
  Additional args:
  
  - resetPackageLocator:  True   will reset the package locator completely, removing old info
                          False  will add the repository to the package locator.
  
  - resetRepository:      True   will reset url for the repository, even if it already exists
                          False  will not reset the url for the repository if it already exists
  
  Returns the relevant repository.
  
  """

  repository = project.findFirstRepository(name = repositoryName)
  ss = normalisePath(repositoryPath)

  if not repository:
    repository = project.newRepository(name= repositoryName, 
                                       url=Implementation.Url(path=ss))
  elif resetRepository and repository.url.path != repositoryPath:
    repository.url = Implementation.Url(path=ss)
  
  for packageName in packageNames:
    packageLocator = project.findFirstPackageLocator(targetName = packageName)
    
    if not packageLocator:
      raise ApiError("Cannot modify repository 'any' for package %s" 
                     % packageName)
    
    if resetPackageLocator:
      packageLocator.repositories = (repository,) 
    elif not repository in packageLocator.repositories:
      packageLocator.addRepository(repository)
  
  return repository
  
def packageProject(project, filePrefix=None, includeBackups=False, includeData=False, userPath=None):
  """
  Package up project userData into one gzipped tar file.
  If filePrefix is None then instead use the userData path.
  The tar file is filePrefix+".tgz".
  By default only *.xml files are packaged up.
  If includeBackups then also *.xml.bak files are included.
  If includeData then also dataStores located inside project directory are included.
  If userPath is not None then use that instead of userPath.
  """

  import tarfile

  if not userPath:
    userPath = getUserDataPath(project)
  userDir = os.path.dirname(userPath)

  if includeData:
    userPathP = userPath + '/'
    n = len(userDir) + 1
    includedDataPaths = set()
    for dataLocationStore in project.dataLocationStores:
      for dataStore in dataLocationStore.dataStores:
        fullPath = dataStore.fullPath
        if fullPath.startswith(userPathP):
          includedDataPaths.add(fullPath[n:])
        
  userPath = os.path.basename(userPath)

  def visitDir(directory):
    tarFiles = []
    files = os.listdir(directory)
    for relfile in files:
      fullfile = os.path.join(directory, relfile)
      include = False
      if os.path.isdir(fullfile):
        tarFiles.extend(visitDir(fullfile))
      elif relfile.endswith('.xml'):
        include = True
      elif includeBackups and relfile.endswith('.xml.bak'):
        include = True
      elif includeData and fullfile in includedDataPaths:
        include = True

      if include:
        print(fullfile)
        tarFiles.append(fullfile)

    if tarFiles:
      tarFiles.insert(0, directory)

    return tarFiles

  if not filePrefix:
    filePrefix = userPath

  tarFile = '%s.tgz' % filePrefix
  tarFp = tarfile.open(tarFile, 'w:gz')

  cwd = os.getcwd()
  os.chdir(userDir)
  try:
    print('Files included in tar file:')
    tarFiles = visitDir(userPath)
    for tarFile in tarFiles:
      tarFp.add(tarFile, recursive=False)
  finally:
    os.chdir(cwd)
    tarFp.close()

def logMemoryStats(project, path):
  """
  Logs stats for a project.
  It appends information to the file given by path.
  """

  path = normalisePath(path)
  fp = open(path, 'a')
  if fp:
    t = time.ctime(time.time())
    m = gc.collect()
    n = len(gc.get_objects())
    s = '%s: %d collected, %d objects' % (t, m, n)
    fp.write('%s\n' % s)
    print('logging: %s' % s)
    fp.close()

def findCcpXmlFile(project,packageName,fileSearchString):
  
  """
  Finds an XML file by a file search pattern from all available package repositories
  """

  ff = project.findFirstPackageLocator
  packageLocator = ff(targetName=packageName) or ff(targetName='any')

  xmlFileName = None
  # The repositories link is ordered!
  for repository in packageLocator.repositories:

    fileLocation = repository.getFileLocation(packageName)
    xmlFileNameMatches = glob.glob(os.path.join(fileLocation,fileSearchString))

    if xmlFileNameMatches:
      if xmlFileNameMatches[-1][-4:] == '.xml':  
        xmlFileName = xmlFileNameMatches[-1]
        break

  return xmlFileName

def getRepositoryPath(project, repositoryName):

  repository = project.findFirstRepository(name=repositoryName)
  if repository:
    path = repository.url.path
  else:
    path = None

  return path

def setRepositoryPath(project, repositoryName, path):

  from memops.api.Implementation import Url

  repository = project.findFirstRepository(name=repositoryName)
  if repository:
    if path != repository.url.path:
      # TBD: should we copy anything over from old url?
      url = Url(path=normalisePath(path))
      repository.url = url

  # TBD: should we throw an exception if repository is not found?

def getBackupPath(project):

  return getRepositoryPath(project, 'backup')

def setBackupPath(project, path):

  return setRepositoryPath(project, 'backup', path)

def getUserDataPath(project):

  return getRepositoryPath(project, 'userData')

def setUserDataPath(project, path):

  return setRepositoryPath(project, 'userData', path)

def getCcpFileString(fileString):
  """
  Changes an input string to the one used for a component of file names. 
  """
  from memops.general.Constants import defaultFileNameChar
  from memops.general.Constants import validFileNamePartChars
  
  return doConvertStringToFileName(fileString, validFileNamePartChars,
                                   defaultFileNameChar)


def findCcpnProjectDirs(inDir):
  """ Find list of (potential) CCPN project dirs within inDir
  """
  result = []
  for dirpath, dirnames, filenames in os.walk(inDir):
    if dirpath.endswith('memops') and 'Implementation' in dirnames:
      result.append(os.path.dirname(dirpath))
  #
  return result
