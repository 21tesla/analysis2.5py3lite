""" Version for Python version >= 2.1
    
    no __init__ here,  must be in subclasses.

"""

import time
import os

from memops.metamodel import MetaModel
MemopsError = MetaModel.MemopsError
from memops.metamodel import ImpConstants

from memops.universal import Io as uniIo

from memops.general import Util as genUtil

True = not 0
False = not True

mandatoryAttributes = ('language', 'startComment', 'endComment',
 'startMultilineComment', 'endMultilineComment',
 'endStatement', 'fileSuffix', 'packageFile', 'baseDirName',
 'classesInPackage', 'releaseVersion', 'scriptName', 'scriptRevision',
 'rootFileName', 'rootDirName',
)

# settings for various contexts
settings = {
 'python': {
  # TextWriter settings :
  'language': 'python',
  'startComment': '#',
  'endComment': '',
  'startMultilineComment': '"""',
  'endMultilineComment': '"""',
  'endStatement': '',
  'fileSuffix': 'py',
  'packageFile': 'package',
  'baseDirName':'python',
  'classesInPackage':True,
   # LanguageInterface settings :
   'noneValue': 'None',
 },
 'java': {
   # TextWriter settings :
  'language': 'java',
  'startComment': '//',
  'endComment': '',
  'startMultilineComment': ('/' + 75 * '*'),
  'endMultilineComment': (75 * '*' + '/'),
  'endStatement': ';',
  'fileSuffix': 'java',
  'packageFile': '_package',
  'baseDirName':'java',
  'classesInPackage':False,
  # LanguageInterface settings :
  'noneValue': 'null',
 },
 'c': {
   # TextWriter settings :
  'language': 'c',
  'startComment': '//',
  'endComment': '',
  'startMultilineComment': ('/' + 75 * '*'),
  'endMultilineComment': (75 * '*' + '/'),
  'endStatement': ';',
  'fileSuffix': 'c',
  'packageFile': '_package',
  'baseDirName':'c',
  'classesInPackage':False,
  # LanguageInterface settings :
  'noneValue': 'NULL',
 },
 'sql': {
  'language': 'sql',
  'startComment': '--',
  'endComment': '',
  'startMultilineComment': None,
  'endMultilineComment': None,
  'endStatement': ';',
  'fileSuffix': 'sql',
  'packageFile': '',
  'baseDirName': None,
  'classesInPackage': True,
   # LanguageInterface settings :
   'noneValue': None,
 },
 'hbm': {
  'language': 'hbm',
  'startComment': '<!--',
  'endComment': '-->',
  'startMultilineComment': '<!--',
  'endMultilineComment': '-->',
  'endStatement': '',
  'fileSuffix': 'hbm.xml',
  'packageFile': '',
  'baseDirName': None,
  'classesInPackage': True,
   # LanguageInterface settings :
   'noneValue': '',
 },
 'xmlmodel': {
  # Used to write out Model.xml
  # TextWriter settings :
  'language': 'xmlmodel',
  'startComment': '<!--',
  'endComment': '-->',
  'startMultilineComment': '<!--',
  'endMultilineComment': '-->',
  'endStatement': '',
  'fileSuffix': 'xml',
  'packageFile': 'package',
  'baseDirName':'model',
  'classesInPackage':False,
   # LanguageInterface settings :
   'noneValue': '',
 },
 'html': {
  # Used to write out HTML documentation
  # TextWriter settings :
  'language': 'html',
  'startComment': '<!--',
  'endComment': '-->',
  'startMultilineComment': '<!--',
  'endMultilineComment': '-->',
  'endStatement': '',
  'fileSuffix': 'html',
  'packageFile': 'package',
  'baseDirName':None,
  'classesInPackage':False,
   # LanguageInterface settings :
   'noneValue': '',
 },
}

class TextWriter_py_2_1:

  # these are files or directories that may not be present in a directory  
  # that is supposed to be cleared out.
  clearOutGuards = ('CVS',)

  INDENT = 2

  DEBUG_OUTPUT = False

  ###########################################################################

  ###########################################################################

  # internal function
  def clearOutDir(self, topDirName, forceClearOut=False):
    """ NB DANGEROUS! clears out contents of topDirName - 
    as in 'cd topDirName; \rm -rf *'
    
    if forceClearOut directories are cleared even if they contain clearOutGuards
    The clearOutGuards are not deleted or entered into even so.
    """
    
    # initial checks
    if os.path.exists(topDirName):
      if not os.path.isdir(topDirName):
        raise MemopsError("clearOutDir input %s is not a directory" 
                          % topDirName)
    else:
      return
    
    
    stem,base =  topDirName,''
    while not base:
      stem,base =  os.path.split(stem)
    while base:
      if base in self.clearOutGuards:
        raise MemopsError("%s may not be cleared out - contains %s"
                          % (topDirName,base))
      stem,base =  os.path.split(stem)
    
    # setup
    dirs = [topDirName]
    files = []
    
    # depth-first search
    for dirName in dirs:
      ll = os.listdir(dirName)
      
      if not forceClearOut:
        for ss in self.clearOutGuards:
          if ss in ll:
            raise MemopsError("cannot clear out %s, contains file or directory %s"
                              % (dirName,ss))
      
      for name in ll:
        if name not in self.clearOutGuards:
          fullname = os.path.join(dirName,name)
          if os.path.isdir(fullname):
            dirs.append(fullname)
          else:
            files.append(fullname)
      
    # remove files
    for name in files:
      # NBNB TBD cannot be done till  directories are reorganised
      #os.remove(name)
      pass
      
    
    # remove dirs (except first one)
    dirs.reverse()
    dirs.pop()
    for name in dirs:
      # NBNB TBD cannot be done till  directories are reorganised
      #os.rmdir(name)
      pass
      
    
  ###########################################################################

  ###########################################################################

  # internal function
  def createDir(self, dirName):
    """ set up as recursive function (instead of using os.makedirs)
    as some implementations might want to do other things in each directory.
    """

    if not os.path.exists(dirName):
      self.createDir(os.path.dirname(dirName))
      try:
        os.mkdir(dirName)
      except:
        print 'Current directory =', os.getcwd()
        raise

  ###########################################################################

  ###########################################################################

  def getDirDepth(self):

    dirDepth = len(self.fileName.split('/')) - len(uniIo.getTopDirectory().split('/')) - 1

    return dirDepth

  ###########################################################################

  ###########################################################################

  def openFile(self, fileName):
    
    if self.fp:
      raise MemopsError("trying to open file '%s' while another file '%s' is still open" % (fileName, self.fileName))

    self.fileName = fileName = fileName + '.' + self.fileSuffix

    self.createDir(os.path.dirname(fileName))

    # TBD: if only have one fp what is point of having multiple indents?
    self.indents.append(self.indent)
    self.indent = 0
    self.fp = open(fileName, 'w')

  ###########################################################################

  ###########################################################################

  def closeFile(self):

    if not self.fp:
      raise MemopsError("trying to close file which is not open")

    self.fp.close()
    self.fp = None
    self.fileName = ''
    self.indent = self.indents.pop()

  ###########################################################################

  ###########################################################################

  #  convenience function to allow debug output
  def writeToFile(self, ss):

    self.fp.write(ss)

    if self.DEBUG_OUTPUT:
      indent = ' ' * self.indent
      ss = ss.strip()
      if ss and not ss.startswith(self.startComment) and not ss.startswith(self.startMultilineComment):
        import inspect
        stack = inspect.stack()
        n = 1
        while True:
          try:
            (fileName, lineNo, function) = stack[n][1:4]
          except:
            return
          fileRoot = os.path.basename(fileName)[:-3]
          if not __name__.endswith(fileRoot): # a bit dangerous this
            break
          n = n + 1
 
        for m in range(4):
          try:
            (fileName, lineNo, function) = stack[n+m][1:4]
            fileRoot = os.path.basename(fileName)[:-3]
          except:
            return
          self.fp.write('%s%s %s %d %s %s\n' % (indent, self.startComment, fileRoot, lineNo, function, self.endComment))
        self.fp.write('\n')
        
  ###########################################################################

  ###########################################################################

  def write(self, ss, compress=True):
    """ Refactored for speed Rasmus Fogh 24/11/06
    """

    if not self.fp:
      raise MemopsError("trying to write line but no file is open")

    if ss:
      
      indent = ' ' * self.indent
      
      lines = ss.splitlines()
      
      if len(lines) == 1:
      
        ss = lines[0].rstrip()
        if ss:
          ##self.fp.write(indent+ss+'\n')
          self.writeToFile(indent+ss+'\n')
          self.previousLineEmpty = False
          
        else:
          if not (compress and self.previousLineEmpty):
            self.fp.write(indent+'\n')
          self.previousLineEmpty = True
        
      else:
      
        # check first line
        ss = lines[0].rstrip()
        if compress and self.previousLineEmpty and not ss:
          lines[1] = indent + lines[1]
          del lines[0]
        else:
          lines[0] = indent + ss 
        
        # check last line
        ss = lines[-1].rstrip()
        self.previousLineEmpty = not ss
        lines[-1] = ss + '\n'
        
        # write lines
        ##self.fp.write(('\n' + indent).join(lines))
        self.writeToFile(('\n' + indent).join(lines))

  ###########################################################################

  ###########################################################################

  def writeOne(self, ss):
    """ Write one line and handle indentation - faster than 'write'
    Prepends the indentation and appends a linebreak, doing no checks or strips
    """

    if not self.fp:
      raise MemopsError("trying to write line but no file is open")

    self.previousLineEmpty = False
    self.fp.write(' ' * self.indent)
    ##self.fp.write(ss)
    self.writeToFile(ss)
    self.fp.write('\n')

  ###########################################################################

  ###########################################################################

  def writeNewline(self):
    
    if not self.previousLineEmpty:
      self.previousLineEmpty = True
      self.fp.write('\n')

  ###########################################################################

  ###########################################################################

  def writeComment(self, comment = None):

    if not comment:
      comment = ''

    lines = comment.split('\n')
    for line in lines:
      self.writeOne('%s %s%s' % (self.startComment, line, self.endComment))

  ###########################################################################

  ###########################################################################

  def writeMultilineComment(self, comment = None, compress=True):

    if not comment:
      comment = ''
    
    if self.startMultilineComment is None or self.endMultilineComment is None:
      self.writeComment(comment)
   
    else:    
      self.write(self.startMultilineComment)
      self.write(comment, compress)
      self.write(self.endMultilineComment)

  ###########################################################################

  ###########################################################################

  def writeStatement(self, s):

    self.write(s + self.endStatement)

  ###########################################################################

  ###########################################################################

  def openObjFile(self, metaObj):
    
    fileName = self.getObjFileName(metaObj)
    self.openFile(fileName)

  ###########################################################################

  ###########################################################################

  def getObjFileName(self, metaObj, absoluteName=True, addSuffix=False):
    """ Get filename for metaObj
    If absoluteName is True the file name is absolute
    otherwise it is relative to the relevant baseDir
    """
    
    # absoulte or relative path
    if absoluteName:
      pathList = [self.rootDirName, self.baseDirName]
    else:
      pathList = []
    
    if metaObj.container is None:
      # Root package only
      pathList.append(ImpConstants.modellingPackageName)
      if self.codeDirName:
        pathList.append(self.codeDirName)
      pathList.append(self.rootFileName)
      
    else:
      # any other object
      ll = metaObj.qualifiedName().split('.')
      if self.codeDirName:
        ll[1:1] = [self.codeDirName]
      pathList.extend(ll)
      
      # special handling for ModelElements that may correspond to directories:
      if (isinstance(metaObj, MetaModel.MetaPackage) 
          and (metaObj.containedPackages or not self.classesInPackage)):
        pathList.append(self.packageFile)
      
      elif not isinstance(metaObj.container, MetaModel.MetaPackage):
        raise MemopsError(" file names not implemented for objects of type %s"
                          % (metaObj.__class__.__name__,))
    
    # add suffix
    if addSuffix and self.fileSuffix:
      pathList[-1] = '%s.%s' % (pathList[-1], self.fileSuffix)
    
    # 
    return uniIo.joinPath(*pathList)

  ###########################################################################

  ###########################################################################

  def getObjDirName(self, metaObj, absoluteName=True):

    fileName = self.getObjFileName(metaObj, absoluteName=absoluteName)
    dirName = uniIo.splitPath(fileName)[0]

    return dirName

  ###########################################################################

  ###########################################################################

  def getVersionString(self, metaobj=None, scriptName=None, 
                       scriptRevision = None, objRevision=None, date=None):

    """
    Function that returns the version header that should be inserted in
    each generated file. Parameters are the metaobj, the name and revision
    of the generating script, and the date. Name and revision default to 
    self.scriptName and self.scriptRevision, if these are set.
    Date defaults to the current date - it is a parameter mainly to
    allow the suppression of date writing (for Model.py files).
    """

    unknown = '?'
    
    if date is None:
      date = time.ctime(time.time())

    dataModelVersion = self.modelPortal.dataModelVersion

    if self.releaseVersion:
      rname = self.releaseVersion.name or unknown
      releaseLine = 'Release %s version %s' % (rname, self.releaseVersion)
    else:
      releaseLine = ''

    if metaobj:
      if objRevision is None:
        objRevision = genUtil.getRepositoryRevision(metaobj)
      pkgLine = ('  from data model element %s revision %s' 
                 % (metaobj.qualifiedName(), objRevision))
    else:
      pkgLine = ''
    
    if scriptName is None:
      scriptName = self.scriptName
    if not scriptName:
      scriptName = unknown

    if scriptRevision is None:
      scriptRevision = self.scriptRevision
    if not scriptRevision:
      scriptRevision = unknown

    lines = []
    lines.append(
     '#######################################################################'
    )
    lines.append(releaseLine)
    lines.append('CCPN Data Model version %s' % dataModelVersion)
    lines.append('')
    lines.append('Autogenerated by %s revision %s on %s' % (scriptName, scriptRevision, date))
    lines.append(pkgLine)
    lines.append('')
    lines.append('#######################################################################')
    #lines.append('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    #lines.append('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')

    return '\n'.join(lines)

  ###########################################################################

  ###########################################################################

  def getCreditsString(self, metaobj, programType):

    from memops.general.license import headers
    from memops.general import Credits
   
    descriptor = metaobj.__class__.__name__
    dirDepth = self.getDirDepth()
    package = metaobj
    while not isinstance(package,MetaModel.MetaPackage):
      package = package.container
    packageGroup = package.taggedValues.get('packageGroup')
    if packageGroup is None:
      raise MemopsError("object %s package lacks packageGroup tagged value"
                        % metaobj)
    
    dd = Credits.licenseInfo[packageGroup].copy()
    dd['fileName'] = '%s.%s' % (metaobj.name, self.fileSuffix)
    dd['licenseLocation'] = dirDepth * '../' + 'license'
    dd['programType'] = programType
    dd['programFunction'] = ('%s %s for CCPN data model, %s %s' 
     % (self.language, programType, descriptor, metaobj.qualifiedName())
    )
    cr = dd.get('credits')
    if cr:
      ss = cr + '\n\n'
    else:
      ss = ''
    dd['credits'] = ss + """
This file was generated with the Memops software generation framework,
and contains original contributions embedded in the framework
"""

    return headers.getHeader(**dd)

  ###########################################################################

  ###########################################################################
