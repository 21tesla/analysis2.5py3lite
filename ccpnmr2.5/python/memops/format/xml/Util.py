""" Python XML I/O Utility code
"""

import os, os.path
import urllib
import string

from memops.metamodel import ImpConstants

from memops.universal import Io as uniIo

from memops.general.Implementation import ApiError

fileSuffix = ".xml"
lenFileSuffix = len(fileSuffix)
keySep = '+'

def getProjectFile(repositoryPath, projectName=None):
  """Get project file given the repositoryPath and optionally the projectName
     (if none given then determined from repositoryPath)
  """

  if not projectName:
    projectName = os.path.basename(repositoryPath)
    
  implDirectory = getImplementationDirectory(repositoryPath)
  
  return uniIo.joinPath(implDirectory, projectName + fileSuffix)

def getImplementationDirectory(repositoryPath):
  """Get implementation directory from the repositoryPath
  """

  return uniIo.joinPath(repositoryPath, ImpConstants.modellingPackageName, ImpConstants.implementationPackageName)

def getTopObjectFile(topObject):
  """Get topObject file name (not path)
  where topObject can be of class MemopsRoot or TopObject
  """

  from memops.general.Io import getCcpFileString
  
  if topObject.root is topObject:
    # This is MemopsRoot
    result = topObject.name + fileSuffix
  
  else:
    ll = [getCcpFileString(str(x)) for x in topObject.getFullKey()]
    ll.append(topObject.guid + fileSuffix)
    result = keySep.join(ll)
  
  return result[-254:]

def getTopObjectPath(topObject):
  """Get topObject (absolute) path
  where topObject can be of class MemopsRoot or TopObject
  """
  
  repositories = topObject.activeRepositories
  if repositories:
    repository = repositories[0]
  else:
    repository = topObject.packageLocator.findFirstRepository()

  repositoryPath = repository.url.getDataLocation()

  result = findTopObjectPath(repositoryPath, topObject)

  return result


def findTopObjectPath(repositoryPath, topObject):
  """Get topObject absolute file path given the repositoryPath, 
  where topObject can be of class MemopsRoot or TopObject.
  
  Will find an existing file fitting the TopObject ID.
  If none is found returns default file name
  """

  suffix = fileSuffix
  lenSuffix = lenFileSuffix
  sep = keySep
  
  if topObject.root is topObject:
    # MemopsRoot
    objId = topObject.name
  else:
    # other TopObject
    objId = topObject.guid
  
  # get default file name
  topObjectDir = uniIo.joinPath(repositoryPath, 
                                *topObject.packageName.split('.'))
  result = uniIo.joinPath(topObjectDir, getTopObjectFile(topObject))
  
  if not os.path.isfile(result):
    # default file name is not there. Look for alternative file that fits ID 
    if os.path.isdir(topObjectDir):
      for filename in os.listdir(topObjectDir):
        if filename.endswith(suffix):
          if filename.split(sep)[-1][:-lenSuffix] == objId:
            result = os.path.join(topObjectDir, filename)
            break
  
  # return whatever result we have
  return result


def areAllTopObjectsPresent(project):
  """ Input: project
  Output: Boolean - True if all loaded TopObjects exist in storage
  """
  
  # set up
  findLocator = project.findFirstPackageLocator
  anyLocator = findLocator(targetName='any')
  allLocations = {}
  
  # check for topObject presence
  result = True
  for topObject in project.topObjects:
 
    if topObject is not project and not topObject.isLoaded:
 
      # get locations
      locator = findLocator(targetName=topObject.packageName) or anyLocator
      locations = allLocations.get(locator)
      if locations is None:
        locations = [x.url.getDataLocation() for x in locator.repositories]
        allLocations[locator] = locations
 
      # check for file presence
      ll = topObject.packageName.split('.')
      ll.append(getTopObjectFile(topObject))
      for location in locations:
        if os.path.isfile(uniIo.joinPath(location, *ll)):
          # file found
          break
          
      else:
        # no file found
        result = False
        break
  
  # 
  return result

def doesRepositoryContainProject(repositoryPath, projectName=None):
  """Does repositoryPath contain project with specified projectName
     (or default projectName if not specified)?
  """

  projectFile = getProjectFile(repositoryPath, projectName)

  return os.path.exists(projectFile)

def getPossibleProjectFiles(repositoryPath):
  """Get the possible project files given the repositoryPath
  """

  if os.path.isdir(repositoryPath):
    implDirectory = getImplementationDirectory(repositoryPath)
    if os.path.isdir(implDirectory):
      files = os.listdir(implDirectory)
      files = [uniIo.joinPath(implDirectory, file) for file in files if file.endswith(fileSuffix)]

      return files

  return []

def getTopObjIdFromFileName(fileName, mustBeMultipart=None):
  """Get project name or TopObject guid from file name (relative or absolute)
  Note: TopObject ID is constrained to not need decoding
  """
  basename = os.path.basename(fileName)
  ll = basename.split(keySep)
  
  if mustBeMultipart is None:
    # no check on number of fields
    pass
    
  elif mustBeMultipart:
    # must be multi-field (normal TopObject)
    if len(ll) == 1:
      raise ApiError("TopObject fileName %s lacks field separators %s" 
                     % (fileName, keySep))
                     
  elif len(ll) != 1:
    # must be single field (Implementation)
    raise ApiError("TopObject fileName %s has field separators %s" 
                   % (fileName, keySep))
  
  
  return ll[-1][:-lenFileSuffix]
