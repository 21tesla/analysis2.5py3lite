""" Read and write models from XML

NB while ObjectDomain is in use this file must work under Jython (Python 2.1)
"""

repositoryTag = '$Name: not supported by cvs2svn $'
repositoryId  = '$Id: XmlModelIo.py,v 1.18 2010-06-25 10:10:31 rhfogh Exp $'

import types

import time
import os
import sys

import memops.universal.Io as uniIo

import memops.general.Constants as genConstants
baseDataTypeModule = genConstants.baseDataTypeModule
infinity = genConstants.infinity

import memops.general.Util as genUtil
import memops.universal.Util as uniUtil

from memops.metamodel import MetaModel
MemopsError = MetaModel.MemopsError
from memops.metamodel import ImpConstants
from memops.metamodel import TaggedValues

from memops.metamodel.ModelPortal import ModelPortal
from memops.general import TextWriter_py_2_1 
from memops.metamodel import ModelTraverse_py_2_1 

from memops.universal.ElementTree import ElementTree
from memops.universal.ElementTree import ElementInclude
XINCLUDE_FALLBACK = ElementInclude.XINCLUDE_FALLBACK
XINCLUDE_INCLUDE = ElementInclude.XINCLUDE_INCLUDE
True = not 0
False = not True

xmlTrue = 'true'
xmlFalse = 'false'

def bool2str(value):
  return value and 'true' or 'false'

def str2bool(value):
  if value in ('True', 'true', '1'):
    return True
  elif value in ('False', 'false', '0'):
    return False
  else:
    raise ValueError("String '%s' is not legal for a Boolean" % value)

def writeModel(topPackage=None, modelPortal=None, rootFileName=None, 
               rootDirName=None, releaseVersion=None, skipImplicit=None,
               **kw):
  """write XML file version of model
  """
  
  if modelPortal is None:
    if topPackage is None:
      raise MemopsError(
       "writeModel must have either modelPortal or topPackage parameter"
      )
    else:
      modelPortal = ModelPortal(topPackage)
      modelPortal.setModelFlavour('language','xmlmodel')
  
  scriptRevision = genUtil.getRepositoryInfo(repositoryId).get('revision')
  
  xmlModelGen = XmlModelGen(modelPortal=modelPortal, rootFileName=rootFileName, 
                        rootDirName=rootDirName, releaseVersion=releaseVersion,
                        scriptName='XmlModelIo', scriptRevision=scriptRevision,
                        skipImplicit=skipImplicit, **kw)
                        
  xmlModelGen.processModel()


def readModel(rootFileName=None, rootDirName=None,
              excludePackageNames=None, includePackageNames=None, 
              checkValidity=True, **kw):
  """XML model reader
  excludePackageNames is a tuple or list of package qualifiedNames
  that will be ignored together with their contents.
  includePackageNames is a tuple or list of package qualifiedNames
  that will be used together with their contents and containers, 
  with all other packages being ignored.
  NB only leaf package names should be put in includePackageNames.
  The names of containing packages will be added automatically.
  Putting e.g. 'ccp' in includePackageNames will *not* cause
  packages contained in ccp to be included.
  """
  start = time.time()
  
  scriptRevision = genUtil.getRepositoryInfo(repositoryId).get('revision')
  
  # expand included files so you also get containers
  if includePackageNames:
    newIncludes = []
    for ss in includePackageNames:
      ll = ss.split('.')
      while ll:
        ss = '.'.join(ll)
        if ss not in newIncludes:
          newIncludes.append(ss)
        ll.pop()
    includePackageNames = newIncludes
  
  xmlModelRead = XmlModelRead(rootFileName=rootFileName,
                              rootDirName=rootDirName, 
                              scriptName='XmlModelIo', 
                              scriptRevision=scriptRevision,
                              excludePackageNames=excludePackageNames,
                              releaseVersion=None,
                              includePackageNames=includePackageNames, **kw)
  
  
  topPackage = xmlModelRead.readModel()
  
  end = time.time()
  print 'Model read finished. Duration %s ' % (end-start)
  
  if checkValidity:
    start = time.time()
    topPackage.checkValid()
    end = time.time()
    print "Model validity checked. Duration %s" % (end-start)
  else:
    print "Model validity check skipped"
    
  
  #
  return topPackage


class TempHolder:
  """ class for temporary storage of parameter data
  """
  pass


class XmlModelRead(TextWriter_py_2_1.TextWriter_py_2_1):

  codeDirName = genConstants.xmlCodeDir
  classNameMapping = {}
  for clazz in MetaModel.nonAbstractClasses:
    classNameMapping[clazz.__name__] = clazz
 
  def __init__(self, **kw):
    
    """To make this work under python 2.1 we need to inline the inits
    of several other files.
    """
    
    # set input attributes
    for (tag, val) in kw.items():
      if not hasattr(self, tag):
        setattr(self, tag, val)
    
    # in-line set xml settings (in lieu of Language init)
    settings = TextWriter_py_2_1.settings['xmlmodel']
    for (tag, val) in settings.items():
      if not hasattr(self, tag):
        setattr(self, tag, val)
    
    # in-line version of TextWriter init
    for tag in TextWriter_py_2_1.mandatoryAttributes:
      if not hasattr(self, tag):
        raise MemopsError(" TextWriter lacks mandatory %s attribute" % tag)     
    
    # special parameters: optional with default values
    if self.rootFileName is None:
      self.rootFileName = ImpConstants.rootPackageDirName
    
    if self.rootDirName is None:
      self.rootDirName = uniIo.getTopDirectory()
    
    # process self.includePackageNames, including container names
    inclNames = self.includePackageNames
    if inclNames:
      for name in inclNames:
        ll = name.split('.')
        for ii in range(1,len(ll)):
          ss = '.'.join(ll[:ii])
          if ss not in inclNames:
            inclNames.append(ss)
    
    self.fp = None
    self.fileName = ''
    self.indent = 0
    self.indents = []
    self.errorMsg = ''
    self.previousLineEmpty = False # used so that do not print out two '\n' in a row
  
  def readModel(self):
    """ read model and return top package. 
    self.objMap is left as a complete guid:object map of the model
    """
    
    # set up
    self.objMap = {}
    self.delayedLoadData = []
    
    # make dummy top object and use to get appropriate file name
    dummy = MetaModel.MetaPackage(name=ImpConstants.rootPackageName,
    guid='www.ccpn.ac.uk_RasmusFogh_2006-06-21-19:13:29_00000',
     taggedValues={'packageGroup':TaggedValues.defaultPackageGroup}
    )
    if hasattr(self, 'rootFile') and self.rootFile:
      rootFile = self.rootFile
    else:
      rootFile = self.getObjFileName(dummy, addSuffix=True)
    del dummy
    
    # parse XML file and recurse over xIncluded files
    self.loadXmlFile(rootFile)
    
    # load delayed data
    self.loadDelayedData()
    
    # finalise objects
    self.finaliseObjects()

    #
    result = self.xmlRoot
    return result
    
    
  def loadXmlFile(self, absFilename):
    """ Recursively load XML file and all Xincluded XML files,
    crating Memops model objects from contents
    """
    
    # set up
    classNameMapping = self.classNameMapping
    objMap = self.objMap
    delayedLoadData = self.delayedLoadData
    currentDir = os.path.dirname(absFilename)
    
    # get tree:
    fp = open(absFilename,'rb')
    
    try:
      elemtree = ElementTree.parse(fp)
    except ImportError:
      from elementtree import SimpleXMLTreeBuilder
      fp.close()
      
      print 'WARNING, standard parser not found - trying alternative parser'
      
      fp = open(absFilename,'rb')
      elemtree = ElementTree.parse(fp,
                                   parser=SimpleXMLTreeBuilder.TreeBuilder())
    fp.close()
    
    # process elements
    # depth-first search.
    # NB - we are using our own iterator instead of ElementTrees for 
    # 1) Python 2.1 compatibility, 2) as we do not need all elements
    
    eList = [elemtree.getroot()]
    for modelElem in eList:
          
      # MetaModelElement - set up
      xmlTag = modelElem.tag
      clazz = classNameMapping.get(xmlTag)
      parameterData = clazz.parameterData
      params = {}
      _tempData = {}

      # process attributes
      for tag,val in modelElem.items():
        # process attribute
        pType = parameterData[tag].get('type')

        if tag == 'container':
          params[tag] = val

        elif pType == 'Boolean':
          if val == xmlTrue:
            params[tag] = True
          elif val == xmlFalse:
            params[tag] = False
          else:
            # NB will give error message later
            params[tag] = val

        elif pType == 'Token':
          params[tag] = val

        elif pType is None:
          _tempData [tag] = val

        else:
          if pType is types.IntType:
            # NB should be unnecessary,
            #    required for Jython 2.1 in ObjectDomain
            params[tag] = int(val)
          else:
            params[tag] = pType(val)

      # check container and make object
      contId = params.get('container')
      if contId is None:
        # RootObject or error
        if xmlTag == 'MetaPackage' and not objMap:
          currentObj = clazz(**params)
          objMap[params.get('guid')] = currentObj
          self.xmlRoot = currentObj
        else:
          raise MemopsError(
           "No container for %s in %s: attributes are:\n%s"
           % (xmlTag, absFilename, modelElem.attrib)
          )

      else:
        container = objMap.get(contId)
        if container is None:
          # container not yet created - postpone
          currentObj = TempHolder()
          currentObj.__dict__.update(params)
          currentObj._clazz = clazz
          delayedLoadData.append(currentObj)
        else:
          params['container'] = container
          currentObj = clazz(**params)
          objMap[params.get('guid')] = currentObj

      # set _tempData
      currentObj._tempData = _tempData
        
      # handle contained elements
      for elem in modelElem:
        tag = elem.tag
        # lxml also provides comments in addition to normal content
        if tag is ElementTree.Comment:
          continue
        if classNameMapping.get(tag) is None:
          # Element content - handle directly
        
          if tag == XINCLUDE_FALLBACK:
            # fallback handled as part of include
            pass
 
          elif tag == XINCLUDE_INCLUDE:
            # read XML file pointed to
 
            # check if file should be read - NB leafPackages only
            packageName = elem.text and elem.text.strip()
            if packageName:
              if (self.includePackageNames
                  and not packageName in self.includePackageNames):
                continue
              if (self.excludePackageNames
                  and packageName in self.excludePackageNames):
                continue
            
            # load file
            try:
              href = self.getHref(elem)
              href = os.path.normpath(os.path.join(currentDir, href))
            except:
              print 'Error in ', currentObj, 'in ', currentDir
              raise
            try:
              self.loadXmlFile(href)
            except IOError:
              # file not found - try fallback
              fallback = elem.find(XINCLUDE_FALLBACK)
              if fallback is None:
                raise MemopsError("Error in file %s - include file %s not found"
                                  % (absFilename, href))
              else:
                href = self.getHref(fallback, mandatory=False)
                if href:
                  # empty fallback - means we can skip
                  # NB vulnerable to parser errors, 
                  #    but they would surely be found elsewhere
                  href = os.path.normpath(os.path.join(currentDir, href))
                  self.loadXmlFile(href)
 
          elif tag == 'documentation':
            #
            setattr(currentObj, tag, elem.text)
 
          else:
            # StringDict, List, or single link
 
            pData = parameterData[tag]
            pType = pData.get('type')
 
            if pType == 'StringDict':
              # StringDict - convert and set attribtue
              dd = {}
              for ee in elem:
                dd[ee.get('tag')] = ee.text
                #if ee.get('tag') == 'python' and tag != 'typeCodes':
                  #print '\n\n\n###', tag, currentObj, 'python\n'
                  #print ee.text
              setattr(currentObj, tag, dd)
 
            elif pData.get('hicard',1) == 1:
              # single link - dereference later
              currentObj._tempData[tag] = elem.text
 
            else:
              # List - dereference later
              ll = [ee.text for ee in elem]
              currentObj._tempData[tag] = ll
        
        else:
          # Child model element - add to node list
          eList.append(elem)
        
  
  def getHref(self, elem, mandatory=True):
    """ get href attribute value. 
    NB this is a hack to handle known bugs innamespace handling for xmllib, 
    conf. elementtree.SimpleXMLTreeBuilder.
    """
    result = elem.get('href')
    if result is None:
      for tag,val in elem.items():
        if tag.split('}')[-1] == 'href':
          result = val
          break
      if result:
        print "WARNING - suspected xmllib bug. tag %s found" % tag
        
    if mandatory and not result:
      raise MemopsError("""Probable xmllib bug - 
no href attribute found for %s element
 attributes were %s""" % (elem.tag, elem.attrib))

    #
    return result
  
  
  def loadDelayedData(self):
    """create objects that have been postponed, because their container was
    not available
    """
    
    # set up
    ll = self.delayedLoadData
    objMap = self.objMap
    
    # loop over delaued data in reverse order till list is empty
    while ll:
      ii = len(ll)
      modified = False
      while ii > 0:
        ii = ii - 1
        obj = ll[ii]
        
        container = objMap.get(obj.container)
        if container:
          # remove temp info and create object
          modified = True
          del ll[ii]
          obj.container = container
          _tempData = obj._tempData
          clazz = obj._clazz
          del obj._tempData
          del obj._clazz
          newObj = clazz(**obj.__dict__)
          newObj._tempData = _tempData
          objMap[obj.guid] = newObj
          
      if not modified:
        print 'Object-Container pairs:'
        for x in ll:
          print x.name, x.guid, ' ; ', x.container
        raise MemopsError(
         "Delayed load data do not load (no container). %s objects remain"
         % len(ll)
        )
  
  def finaliseObjects(self):
    """ set remaining paramters from _tempData
    """ 
    
    objMap = self.objMap
    
    for obj in objMap.values():
      
      _tempData = obj._tempData
      if _tempData:
      
        valueTypeData = []
        parameterData = obj.parameterData
        
        while _tempData:
          # set up
          tag,val = _tempData.popitem()
          pData = parameterData[tag]
          pType = pData.get('type')
          
          # convert single value to list
          hicard =  pData.get('hicard',1)
          if hicard ==1:
           val = [val]
          
          # process values
          if pType is None:
            # valueType parameter - postpone for handling after links are set
            valueTypeData.append((tag, hicard, val))
          
          else:
          
            if pType == 'Boolean':
              # Boolean  parameter
              for ii in range(len(val)):
                vv = val[ii]
                if vv == xmlTrue:
                  val[ii] = True
                elif vv == xmlFalse:
                  val[ii] = False
                else:
                  raise MemopsError("illegal Boolean %s in %s parameter %s"
                                    % (vv, obj, tag))
                                    
            elif pData.get('isLink'):
              # Link parameter
              for ii in range(len(val)):
                try:
                  val[ii] = objMap[val[ii]] 
                except KeyError:
                  raise MemopsError(
                   "%s parameter %s - link to nonexisting object %s"
                    % (obj, tag, val[ii])
                   )
            
            else:
              # primitive data type parameter
              for ii in range(len(val)):
                try:
                  val[ii] = pType(val[ii])
                except:
                  raise MemopsError(
                   "%s parameter %s - error converting %s"
                    % (obj, tag, val[ii])
                   )
            
            # set parameter
            if hicard ==1:
              setattr(obj,tag,val[0]) 
            else:
              setattr(obj,tag,val) 
          
        for tag, hicard, val in valueTypeData:
          # valueType, type externally determined
          
          # get fromString function
          if isinstance(obj, MetaModel.MetaDataType):
            valueType = obj
          else:
            valueType = obj.valueType
          typeCode = valueType.typeCodes['python']
          if typeCode == 'Boolean':
            fromString = str2bool
          else:
            fromString = getattr(genConstants.baseDataTypeModule,
                                 typeCode).fromString
          
          # convert values
          for ii in range(len(val)):
            try:
              val[ii] = fromString(val[ii])
            except:
              raise MemopsError(
               "%s parameter %s - error converting %s"
                % (obj, tag, val[ii])
               )
            
          # set parameter
          if hicard ==1:
            setattr(obj,tag,val[0])
          else:
            setattr(obj,tag,val)
        
        
class XmlModelGen(TextWriter_py_2_1.TextWriter_py_2_1, 
                  ModelTraverse_py_2_1.ModelTraverse_py_2_1):

  codeDirName = genConstants.xmlCodeDir
  
  _xmlSpecialTags = ('name', 'guid', 'container', 'documentation')
 
  def __init__(self, **kw):
    
    """To make this work under python 2.1 we need to inline the inits
    of several other files.
    """
    
    # XmlModelGen init
    self.addModelFlavour('language', 'xmlmodel')
    for (tag, val) in kw.items():
      if not hasattr(self, tag):
        setattr(self, tag, val)
    
    # TopObject and DataRoot for us below
    ic = ImpConstants
    Impl = self.modelPortal.metaObjFromQualName('.'.join(
     [ic.modellingPackageName,ic.implementationPackageName]
    ))
    self.DataRoot = Impl.getElement(ic.dataRootName)
    self.TopObject = Impl.getElement(ic.topObjClassName)
    
    # in-line set xml settings (in lieu of Language init)
    settings = TextWriter_py_2_1.settings['xmlmodel']
    for (tag, val) in settings.items():
      if not hasattr(self, tag):
        setattr(self, tag, val)
    
    # in-line version of TextWriter init
    for tag in TextWriter_py_2_1.mandatoryAttributes:
      if not hasattr(self, tag):
        raise MemopsError(" TextWriter lacks mandatory %s attribute" % tag)     
    
    # special parameters: optional with default values
    if self.rootFileName is None:
      self.rootFileName = ImpConstants.rootPackageDirName
    if self.rootDirName is None:
      self.rootDirName = uniIo.getTopDirectory()
    
    if self.skipImplicit is None:
      self.skipImplicit = True
    
    self.fp = None
    self.fileName = ''
    self.indent = 0
    self.indents = []
    self.errorMsg = ''
    self.previousLineEmpty = False # used so that do not print out two '\n' in a row

    
    # in-line version of ModelTraverse init
    for tag in ModelTraverse_py_2_1.mandatoryAttributes:
      if not hasattr(self, tag):
        raise MemopsError(" ModelTraverse lacks mandatory %s attribute" % tag)

    # has to be done this way to allow for different initialisation orders
    if not hasattr(self, 'modelFlavours'):
      self.modelFlavours = {}
    
    # input check
    if not isinstance(self.modelPortal, ModelPortal):
      raise MemopsError("ModelTraverse input %s is not a ModelPortal"
                        % self.modelPortal)
    
    # link varNames for easier access, and check modelPortal has been processed.
    # NBNB check if needed
    if hasattr(self, 'varNames'):
      # this must have been called from ModelAdapt
      pass
    elif hasattr(self.modelPortal, 'varNames'):
      self.varNames = self.modelPortal.varNames
      self.operationData = self.modelPortal.operationData
    
    # make repositoryTag and repositoryId dictionary
    # for every qualifiedName get previous repository ID information
    # In this way you can use the previous info for newly generated models
    self.oldRepData = oldRepData = {}
    try:
      dd = kw.copy()
      del dd['scriptName']
      del dd['scriptRevision']
      del dd['releaseVersion']
      oldTop = readModel(checkValidity=False, **dd)
      oldPackages = [oldTop]
      while oldPackages:
        package = oldPackages.pop()
        pName = package.qualifiedName()
        
        # set dict entry
        ff = package._MetaModelElement__dataDict['taggedValues'].get
        x1, x2 = ff('repositoryId'), ff('repositoryTag')
        if x1 or x2:
          oldRepData[pName] = (x1, x2)
        
        ll = package.containedPackages
        if ll:
          # add contained packages to loop
          oldPackages.extend(ll)
          
        else:
          # set dict entries for contents
          for tag in ('dataTypes', 'constants', 'exceptions','dataObjTypes',
                      'classes'):
            for metaObj in getattr(package, tag):
              ff = metaObj._MetaModelElement__dataDict['taggedValues'].get
              x1, x2 = ff('repositoryId'), ff('repositoryTag')
              
              if x1 or x2:
                ss = '%s.%s' % (pName, metaObj.name)
                oldRepData[ss] = (x1,x2)
    except:
      print 'WARNING - import of XML reader failed'
      pass
    
    # set xml tag lists
    self._xmlInfo = {
     'xmlAttrs':{},
     'stringDictElements':{},
     'plainElements':{},
     'listElements':{},
    }
    
    for clazz in MetaModel.nonAbstractClasses: 
      
      xmlAttrs = self._xmlInfo['xmlAttrs'][clazz] = []
      stringDictElements = self._xmlInfo['stringDictElements'][clazz] = []
      plainElements = self._xmlInfo['plainElements'][clazz] = []
      listElements = self._xmlInfo['listElements'][clazz] = []
    
      parameterData = clazz.parameterData
 
      for tag,parData in parameterData.items():
 
        parType = parData.get('type')
 
        if tag in self._xmlSpecialTags:
          # special handling
          continue
 
        elif parData.get('setterFunc') == 'unsettable':
          # ignore
          continue
 
        elif parType == 'content':
          # ignore
          continue
 
        elif parType == 'StringDict':
          stringDictElements.append(tag)
          
        elif parData.get('hicard',1) != 1:
          listElements.append(tag)
 
        elif parData.get('isLink'):
          plainElements.append(tag)
 
        else:
          # xml attribute
          xmlAttrs.append(tag)
 
      # adapt results
      xmlAttrs.sort()
      stringDictElements.sort()
      listElements.sort()
      plainElements.sort()
      plainElements.insert(0,'documentation')
      
    
  ###########################################################################

  ###########################################################################
  ###
  ### code Looping over objects
  ###
  ###########################################################################

  ###########################################################################


  def processBranchPackage(self, package):
    
    # Implicit elements are generally not written
    if self.skipImplicit and package.isImplicit:
      return
    
    packageFile = '%s.%s' % (self.packageFile,self.fileSuffix)
    
    # clear out old model files - NBNB TBD check
    if package.container and package.container.container is None:
      # this is top level package 
      self.clearOutDir(self.getObjDirName(package), forceClearOut=True)
    
    # write package
    self.startXmlFile(package)
    self.startXmlElement(package)
    
    # handle XInclude statements for contents
    self.writeNewline()
    self.writeComment("Start packages contained in %s " % package.name)
    if package.container:
      # non-root package
      for pp in package.containedPackages:
        self.writeXmlInclude(os.path.join(pp.name, packageFile), 
                             objName=pp.qualifiedName(), isMandatory=False)
    else:
      # root package
      for pp in package.containedPackages:
        fileName = self.getObjFileName(pp, absoluteName=False, addSuffix=True)
        self.writeXmlInclude(
         os.path.join((self.getDirDepth()-1)*'../', fileName), 
         objName=pp.qualifiedName(), isMandatory=False
        )
      
    # end file
    self.writeComment("End packages contained in %s " % package.name)
    self.endXmlElement(package)
    self.endXmlFile()

  ###########################################################################

  ###########################################################################
    
  def processLeafPackage(self, package):
    """ processing actions for LeafPackage
    """
    
    # Implicit elements are generally not written
    if self.skipImplicit and package.isImplicit:
      return
 
    self.startXmlFile(package)
    self.startXmlElement(package)
    
    # handle XInclude statements for contents
    for tag in ('dataTypes', 'constants', 'exceptions','dataObjTypes', 
                'classes'):
      self.writeNewline()
      self.writeComment("Start %s contained in %s " 
                        % (tag, package.qualifiedName()))
      for obj in getattr(package, tag):
        self.writeXmlInclude("%s.%s" % (obj.name, self.fileSuffix))
      self.writeComment("End %s contained in %s " 
                        % (tag, package.qualifiedName()))
    
    self.endXmlElement(package)
    self.endXmlFile()
    
    # write content files
    for tag in ('dataTypes', 'constants', 'exceptions','dataObjTypes', 
                'classes'):
      for metaObj in getattr(package, tag):
        self.processPackageContent(metaObj)
    
  ###########################################################################

  ###########################################################################
    
  def processPackageContent(self, metaObj):
    """ processing actions package content
    """
    
    # Implicit elements are generally not written
    if self.skipImplicit and metaObj.isImplicit:
      return
 
    self.startXmlFile(metaObj)
    
    self.writeElementRecursive(metaObj)
    
    self.endXmlFile()
    
  ###########################################################################

  ###########################################################################
    
  def writeElementRecursive(self, metaObj):
    """ processing actions package content
    """
    
    # Implicit elements are generally not written
    if self.skipImplicit and metaObj.isImplicit:
      return
 
    self.writeNewline()
    self.startXmlElement(metaObj)
    
    # handle constraints
    if isinstance(metaObj,MetaModel.ConstrainedElement):
      constraints = metaObj.constraints
      if constraints:
        self.fp.write('\n')
        self.writeComment("Start constraints for %s: " % metaObj)
        for constraint in constraints:
          self.startXmlElement(constraint)
          self.endXmlElement(constraint)
        self.writeComment("End constraints for %s: " % metaObj)
    
    if isinstance(metaObj, MetaModel.MetaClass):
      supertypes = metaObj.getAllSupertypes()
    
      if (self.TopObject in supertypes and not self.DataRoot in supertypes):
        # handle special case - currentTopObject link and its operations
        #  stored on wrong side
 
        ll = [x for x in self.DataRoot._MetaModelElement__elementDict.values()
              if isinstance(x, MetaModel.MetaRole)
              and x.valueType is metaObj
              and x.otherRole is None]
        for obj in ll:
          self.writeElementRecursive(obj)
 
        ll = [x for x in self.DataRoot._MetaModelElement__elementDict.values()
              if isinstance(x, MetaModel.MetaOperation)
              and isinstance(x.target, MetaModel.MetaRole)
              and x.target.valueType is metaObj]
        for obj in ll:
          self.writeElementRecursive(obj)
    
    # handle normal cases
    ll = metaObj._MetaModelElement__elementDict.items()
    ll.sort()
    if ll:
    
      self.fp.write('\n')
 
      for junk,obj in ll:
 
        if isinstance(obj, MetaModel.MetaRole):
 
          # handle interpackage crosslinks
          otherRole = obj.otherRole
          if otherRole is None:
            if (self.DataRoot in obj.container.getAllSupertypes() 
                and self.TopObject in obj.valueType.getAllSupertypes()):
              # handle special case
              # - currentTopObject link is stored on wrong side
              continue
            else:
              pass
 
          elif not obj.canAccess(otherRole):
            # interpackage link handled from the other side. Skip.
            continue
 
          elif not otherRole.canAccess(obj):
            # interpackage link handled from this side
            # write reverse role first.
            self.writeElementRecursive(otherRole)
        
        if isinstance(obj, MetaModel.MetaOperation):
          target = obj.target
          if (self.DataRoot in obj.container.getAllSupertypes() 
              and isinstance(target, MetaModel.MetaRole)  
              and self.TopObject in target.valueType.getAllSupertypes()):
            # handle special case
            # - getter for currentTopObject link is stored on wrong side
            continue
        
        self.writeElementRecursive(obj)
 
    self.endXmlElement(metaObj)
    
  ###########################################################################

  ###########################################################################
  ###
  ### other code
  ###
  ###########################################################################
 
  ###########################################################################
 
  def startXmlFile(self, obj):
    
    # first set tagged values for repositoryId and repositoryTag
    ff = obj._MetaModelElement__dataDict['taggedValues'].get
    if not ff('repositoryId') or not ff('repositoryTag'):
      tt = self.oldRepData.get(obj.qualifiedName())
      if tt is None:
        obj.addTaggedValue('repositoryId', genConstants.emptyRepositoryId)
        obj.addTaggedValue('repositoryTag', genConstants.emptyRepositoryTag)
      else:
        obj.addTaggedValue('repositoryId', tt[0])
        obj.addTaggedValue('repositoryTag', tt[1])
    
    # do actual work
    self.openFile(self.getObjFileName(obj))
    
    self.writeOne('<?xml version="1.0"?>')
    
    self.writeMultilineComment(
     (self.getVersionString(metaobj=obj, objRevision='?', date='?') + 
     '\n' + 
     self.getCreditsString(metaobj=obj, programType='model')),
     compress=False
    )

  ###########################################################################
 
  ###########################################################################
 
  def endXmlFile(self):
    
    self.writeNewline()
    self.closeFile()

  ###########################################################################

  ###########################################################################
  
  def startXmlElement(self, metaObj):
    """ Write start of XML element for MetaObj.
    """
    
    # starting indentation level
    ind = self.indent * ' '
    indx = '\n' + ind
    
    clazz = metaObj.__class__
    
    # write start element
    # first element-independent lines
    self.fp.write('%s<%s name="%s"' % (ind, metaObj.__class__.__name__, 
                                    metaObj.name))
    self.fp.write('%s guid="%s"' % (indx, metaObj.guid))
    container = metaObj.container
    if container:
      # NB not all elements are inside the element of their container
      self.fp.write('%s container="%s"' % (indx, container.guid))
    
    # then get lines for other attributes
    strings = []
    for tag in self._xmlInfo['xmlAttrs'][clazz]:
      val = self.toXmlString(metaObj, tag)
      if val:
        strings.append(' %s="%s"' % (tag, val))
    strings = uniUtil.compactStringList(strings, maxChars=80-self.indent)
    
    # then write lines to file
    for ss in strings:
      self.fp.write(indx)
      self.fp.write(ss)
    self.fp.write('>\n')
    
    # write contents
    self.indent += self.INDENT
    for tag in self._xmlInfo['plainElements'][clazz]:
      self.writePlainXmlElement(metaObj, tag)
    
    for tag in self._xmlInfo['listElements'][clazz]:
      self.writeList(metaObj, tag)
    
    for tag in self._xmlInfo['stringDictElements'][clazz]:
      self.writeStringDict(metaObj, tag)
  
  def endXmlElement(self, metaObj):
    """ Write end of XML element for MetaObj.
    """
    
    self.indent -= self.INDENT
    self.writeOne("</%s>" % (metaObj.__class__.__name__))
    
  
  def writePlainXmlElement(self, metaObj, tag):
    """ write xmlElement for metaObj.tag
    """
    
    ss = self.toXmlString(metaObj, tag)
    if ss:
      self.writeOne('<%s>%s</%s>' % (tag, ss, tag))
    
  
  def writeStringDict(self, metaObj, tag):
    """ write StringDict element
    """
    
    value = getattr(metaObj,tag)
    
    if value:
    
      self.writeOne("<%s>" % (tag))
    
      self.indent += self.INDENT
      ind = self.indent * ' '
    
      items = value.items()
      items.sort()
      
      for key,val in items:
        # convert key
        key = key.replace("&", "&amp;")
        key = key.replace("\"", "&quot;")
        key = key.replace("<", "&lt;")
        key = key.replace(">", "&gt;")
        #convert val
        val = val.replace("&", "&amp;")
        val = val.replace("<", "&lt;")
        val = val.replace(">", "&gt;")
        # write element. NB strings can *not* be indented
        self.writeOne('<item tag="%s">%s</item>' % (key, val))
        
      self.indent -= self.INDENT
      self.writeOne("</%s>" % (tag))
  
  
  def writeList(self, metaObj, tag):
    """ write List element
    """
    
    value = getattr(metaObj, tag)
    
    if value:
    
      self.writeOne("<%s>" % (tag))
    
      self.indent += self.INDENT
      ind = self.indent * ' '
 
      for val in value:
        # write element. NB strings can *not* be indented
        self.writeOne('<item>%s</item>' 
                      % (self.toXmlString(metaObj, tag, val)))
 
      self.indent -= self.INDENT
      self.writeOne("</%s>" % (tag))
  
  
  def toXmlString(self, metaObj, tag, value=None):
    """ Convert metaObj.tag to string ready to write to XML
    """
    pData = metaObj.parameterData[tag]
    parType = pData.get('type')
    
    if value is None:
      value = getattr(metaObj, tag)
    
    if value is None:
      return  ''
    
    if parType == 'Boolean':
      return value and xmlTrue or xmlFalse
    
    elif pData.get('isLink'):
      # link
      result = value.guid
    
    elif not parType:
      # valueType, type externally determined
      # NB all legal cases will fit with this
      if isinstance(metaObj, MetaModel.MetaDataType):
        dataType = metaObj
      else:
        dataType = metaObj.valueType
      typeCode = dataType.typeCodes['python']
      if typeCode == 'Boolean':
        result = value and xmlTrue or xmlFalse
      else:
        toString = getattr(genConstants.baseDataTypeModule,typeCode).toString
        result = toString(value)
    
    else:
      # all other types
      result = str(value)
    
    # xmlify result
    result = result.replace("&", "&amp;")
    result = result.replace("\"", "&quot;")
    result = result.replace("<", "&lt;")
    result = result.replace(">", "&gt;")
    
    #
    return result
  
  def writeXmlInclude(self, filePath, objName=None, isMandatory=True):
    """write an XML include statement pointing to the file. 
    If inclusion is non-mandatory (default), the fallback will be 
    an empty string.
    """ 
    if isMandatory and not objName:
      self.writeOne(
       '<xi:include href="%s" xmlns:xi="http://www.w3.org/2001/XInclude"/>'
       % filePath
      )
    else:
      self.writeOne(
       '<xi:include href="%s" xmlns:xi="http://www.w3.org/2001/XInclude">'
       % filePath
      )
      self.indent += self.INDENT
      if objName:
        self.writeOne(objName)
      if not isMandatory:
        self.writeOne('<xi:fallback></xi:fallback>')
      self.indent -= self.INDENT
      self.writeOne('</xi:include>')
      
