"""
======================COPYRIGHT/LICENSE START==========================

UpdateAgent.py: Part of the CcpNmr Update program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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

===========================REFERENCE END===============================

"""
import compileall
import base64
##import httplib
import os
import py_compile
import sys
import os.path as path
import urllib
import filecmp

try:
  import urllib2
except:
  pass

from shutil import copyfile
from time import ctime 

from memops.universal.Io import getTopDirectory
from memops.gui.MessageReporter import showWarning, showOkCancel
from memops.gui.DataEntry import askPassword

from memops.universal.Io import joinPath
from memops.universal.Util import isWindowsOS


UPDATE_SERVER_LOCATION = 'www2.ccpn.ac.uk'
UPDATE_DIRECTORY = 'ccpNmrUpdate2.4'
UPDATE_DATABASE_FILE = '__UpdateAgentData.db'
UPDATE_HTTP_DIR = 'updateadmin'
UPDATE_UID = 'ccpn'

# 20190322:ED updated for new server
UPDATE_SCRIPT = 'cgi-bin/update/downloadFile'
UPDATE_BASE_DIRECTORY = 'ccpNmrUpdate'
BAD_DOWNLOAD = 'Exception: '
VERSION_ROUTE= 'ccpnmr.analysis.Version'
VERSION_ATTRIBUTE = 'version'

fieldSep = '\t'
environmentFile = 'environment.txt'
SERVER_USER_ROOT = '/home'
CCPNMRCODEDIR = 'ccpnmr2.5'
CCPNMRTOPDIR = 'ccpnmr'


def getNumUninstalledUpdates(serverLocation=UPDATE_SERVER_LOCATION,
                             serverDirectory=UPDATE_DIRECTORY,
                             dataFile=UPDATE_DATABASE_FILE):

  updateAgent = UpdateAgent(serverLocation, serverDirectory,
                            dataFile, isGraphical=False)
  server = updateAgent.server
  if server:
    server.getFileUpdates()
    files = [file for file in server.fileUpdates if not file.getIsUpToDate()]
    return len(files)
  
  else:
    return 0


def areFilesIdentical(fileName1, fileName2):
  
  return filecmp.cmp(fileName1, fileName2, shallow=False)
  

class UpdateAgent:

  def __init__(self, serverLocation, serverDirectory, dataFile,
               httpDir=UPDATE_HTTP_DIR, installRoot=None, isGraphical=True,
               versionRoute=VERSION_ROUTE, admin=False,isStandAlone=False):

    module  = __import__(versionRoute, {}, {}, [VERSION_ATTRIBUTE])

    self.version         = getattr(module, VERSION_ATTRIBUTE)
    self.serverLocation  = serverLocation
    self.serverDirectory = serverDirectory
    self.dataFile        = dataFile
    self.isGraphical     = isGraphical
    self.isStandAlone    = isStandAlone

    # 20190322:ED
    self._UPDATE_DIRECTORY = UPDATE_BASE_DIRECTORY+self.version

    self.installRoot = installRoot
    if self.installRoot is None:
      self.getInstallRoot()

    self.tempDir = joinPath(self.installRoot,'python','ccpnmr','update','temp')

    if not os.path.isdir( self.tempDir ):
      os.makedirs( self.tempDir )

    # single server at the moment but could be more in future
    self.server = None
    if self.testWriteAccess():
      if self.serverLocation is not None:
        self.setServer(serverLocation, serverDirectory, dataFile, httpDir, admin=admin)
  
  def warningMessage(self, title, message):
  
    if self.isGraphical:
      showWarning(title, message)
    else:
      print 'CcpNmr UpdateAgent  - %s %s' % (title, message)
  
  def testWriteAccess(self):
  
    try:
      fileName = joinPath(self.installRoot,'__write_test__')
      file = open(fileName ,'w')
      file.close()
      os.remove(fileName)
    except:
      if not isWindowsOS():
        message = 'File updates not possible: You do not have write access to %s%s' % (self.installRoot,os.sep)
      else:
        message = 'File updates not possible: You do not have write access to %s%s' % (self.installRoot,os.sep) \
                + '\n\nTry running Analysis as administrator by right clicking on the Analysis ' \
                + 'shortcut on the desktop and selecting "Run as administrator", and then apply any updates.'
      self.warningMessage('Failure',message)
      return False

    return True
 
  def setServer(self, location, directory, dataFile, httpDir, admin=False):
  
    if self.server:
      self.server.delete()
      
    self.serverLocation  = location
    self.serverDirectory = directory
    self.server          = UpdateServer(self, location, directory, dataFile, httpDir, admin=admin)  
    

  def installNewUpdates(self):

    if self.server:
      for fileUpdate in self.server.fileUpdates:
         if fileUpdate.getIsUpToDate():
           fileUpdate.isSelected = False
         else:
           fileUpdate.isSelected = True

    self.installUpdates()


  def installUpdates(self):
  
    if self.isGraphical and not showOkCancel('Confirmation','Continue with file upgrade?'):
      return
  
    if not self.server:
      self.warningMessage('Warning','No accessible server')
      return
 
    if not self.installRoot:
      self.warningMessage('Warning','No installation path')
      return
  
    updates = self.server.getSelectedUpdates()
    
    if not updates:
      self.warningMessage('Warning','No selected updates to install')
      return
    
    needMake = False
    needMakeClean = False
    needMakeLinks = False
    for update in updates:
      update.setInstallFromCache()
      if update.fileName.endswith('.c'):
        needMake = True
      elif update.fileName.endswith('.h'):
        needMake = needMakeClean = True
      elif update.fileName == 'Makefile':
        needMake = needMakeClean = needMakeLinks = True

    #if not isWindowsOS():
    if False:  # 16 Nov 2010: do not even try to compile C code any more
      self.compileCode(needMakeClean=needMakeClean, needMake=needMake, needMakeLinks=needMakeLinks)

    if self.isGraphical and not self.isStandAlone:
      self.warningMessage('Notice','All updates complete. You must restart Analysis for changes to take effect.')

  def compileCode(self, needMakeClean=True, needMake=True, needMakeLinks=True):

    cwd = os.getcwd()
    compileDir = os.path.join(self.installRoot, 'c')
    try:
      os.chdir(compileDir)
    except Exception, e:
      print 'Error trying to cd into directory to compile C code; skipping, release will not be up-to-date:', e
      return

    cmds = []
    if needMakeClean:
      cmds.append('make clean')
    if needMake:
      cmds.append('make')
    ### replace make links code with code at bottom
    ###if needMakeLinks:
    ###  cmds.append('make links')

    try:
      if cmds:
        self.runCmds(cmds)
    finally:
      os.chdir(cwd)
    
    if needMakeLinks:
      pythonDir = os.path.join(self.installRoot, 'python')
      try:
        os.chdir(pythonDir)
      except Exception, e:
        print 'Error trying to cd into directory to make links for C code; skipping, release may not be up-to-date:', e
        return

      try:
        script = './linkSharedObjs'
        cmds = ['chmod u+x %s' % script, script]
        for directory in ('memops/c', 'ccp/c', 'ccpnmr/c'):
          os.chdir(directory)
          self.runCmds(cmds)
          os.chdir('../..')
      finally:
        os.chdir(cwd)
    
  def runCmds(self, cmds):

    cmd = ';'.join(cmds)
    os.system(cmd)

  def getInstallRoot(self):
  
    self.installRoot = getTopDirectory()

  def installLatestRelease(self):
  
    if self.server:
      msg = 'Current version is %s - ' % self.version + \
            'You must upgrade to version %s to get updates. ' % self.server.version + \
            'Continue with automatic upgrade to the latest major release? ' + \
            '(Old version will not be deleted from disk)'

      if self.isGraphical and not showOkCancel('Confirmation',msg):
        return
 
      if not self.server:
        self.warningMessage('Warning','No accessible server')
        return
 
      if not self.installRoot:
        self.warningMessage('Warning','No installation path')
        return
 
      releaseUpdate = ReleaseUpdate(self.server, self.installRoot, self.version)
      success = releaseUpdate.installRelease()

      if success:
        self.server.parent.version = self.server.version

      if success and self.isGraphical:
        
        self.server.getFileUpdates()
        if self.server.fileUpdates and showOkCancel('Confirmation','Also install any updates to new release?'):
          wasGraphical = self.isGraphical
          self.isGraphical = False
          self.installNewUpdates()
          for fileUpdate in self.server.getSelectedUpdates():
            print 'Updated %s in %s' % (fileUpdate.fileName, fileUpdate.filePath)
          self.isGraphical = wasGraphical
        
        if hasattr(self, 'parent'):
          if not self.isStandAlone:
            self.warningMessage('Notice','Release update complete. Program will exit. Restart Analysis for changes to take effect.')

          parent = self.parent
          # parent.project does not exist in some contexts (e.g. updateCheck)
          if hasattr(parent, 'project') and parent.project:
            if not parent.checkSaving():
              return

          parent.destroy()

          # Need the below code to fix the problem with Bash
          # where the commend line was getting screwed up on exit.
          if os.name == 'posix':
            os.system('stty sane')
 
          sys.exit(0)
 
        else:
          if not self.isStandAlone:
            self.warningMessage('Notice','Release update complete. You must restart CCPN software for changes to take effect.')
        

class UpdateServer:

  def __init__(self, agent, location=UPDATE_SERVER_LOCATION,
               directory=UPDATE_DIRECTORY, dataFile=UPDATE_DATABASE_FILE,
               httpDir=UPDATE_HTTP_DIR, uid=UPDATE_UID, admin=0):

    self.parent   = agent
    self.location = location
    self.dataFile = dataFile
    self.httpDir  = httpDir
    self.url      = joinPath('%s' % location, '~%s' % uid,directory)
    self.basedir  = joinPath(SERVER_USER_ROOT, '%s' % uid, 'public_html' , directory)
    self.uid      = uid
    self.identity = (location,uid,httpDir,directory)
    self.admin    = admin
    self.version  = None
    
    self.fileUpdates = []
    self.getFileUpdates()

  def setFileUpdates(self, fileUpdates=None, refresh=False, passwd=None):
  
    # Synchronise the server - This is the admin upload part
    
    if not fileUpdates:
      fileUpdates = self.fileUpdates
    
    if passwd is None:
      if self.parent.isGraphical:
        passwd = askPassword('Password Request','Enter password for %s@%s' % (self.uid,self.location))
      else:
        self.parent.warningMessage('Warning','Password must be specified when setting updates in non graphical mode')
     
    if not passwd:
      return False
    
    # Clean server files
    self.deleteFile(passwd, '__temp_*')

    added = 0
    # if self.fileUpdates:

    fileName = joinPath(self.parent.tempDir, self.dataFile)
    # file     = open(fileName, 'w')

    # 20190322:ED correct file opening
    with open(fileName, 'w') as file:
      file.write('%s\n' % self.parent.version)
      for x in self.fileUpdates:
        if (not x.isNew) and refresh and (x in fileUpdates):
          if not x.getIsUpToDate():
            x.timestamp()
            x.isNew = True

        data = fieldSep.join([x.fileName, x.filePath, x.storedAs, x.language, x.date, x.details, str(x.priority)])
        file.write('%s\n' % data)

        if x.isNew:
          copyfile(x.installedFile, x.tempFile)
          # self.uploadFile(passwd, x.tempFile,x.storedAs)

          # 20190322:ED write update files to the server
          try:
            with open(x.tempFile, 'rb') as fp:
              fileData = fp.read()
          except:
            print 'error reading file, not unicode'
            fileData = ''
          self._uploadFile('ccpn', passwd, 'https://www.ccpn.ac.uk/cgi-bin/updateadmin/uploadFile', fileData, self.identity[-1], x.storedAs)
          added += 1

    # file.close()

    # self.uploadFile(passwd, fileName,self.dataFile)

    # 20190322:ED write update file database to the server
    try:
      with open(fileName, 'rb') as fp:
        fileData = fp.read()
    except Exception, es:
      print 'error reading file, not unicode', str(es)
      fileData = ''
    self._uploadFile('ccpn', passwd, 'https://www.ccpn.ac.uk/cgi-bin/updateadmin/uploadFile', fileData, self.identity[-1], self.dataFile)

    # else:
    #   self.deleteFile(passwd, self.dataFile)
    
    self.parent.warningMessage('Notice','Server Update of %d files' % added)
    
    return True

  def uploadFile(self, passwd, localFile, serverFile):

    cc = open(localFile).read()
    # ss = joinPath(self.basedir, serverFile)

    # 20190322:ED get the update directory from the last element of identity
    ss = joinPath(self.identity[-1], serverFile)

    data = urllib.urlencode({'content': cc, 'file': ss})
    try:
      # 18 Aug 08: TEMP: TBD: remove Temp when Jenny password protects directory
      # self.callHttpScript(passwd, 'uploadFileTemp', data)

      # 20190322:ED different script
      self.callHttpScript(passwd, 'uploadFile', data)

    except Exception, e:
      self.parent.warningMessage('Server', 'Server exception: %s' % str(e))

  def deleteFile(self, passwd, serverFile):

    print 'deleteFile', serverFile


  def callHttpScript(self, passwd, script, data):

    auth = base64.encodestring(self.uid + ":" + passwd)[:-1]
    authheader = 'Basic %s' % auth
    uri = 'http://' + joinPath(self.location, 'cgi-bin', self.httpDir, script)
    req = urllib2.Request(uri)
    req.add_header("Authorization", authheader)
    req.add_data(data)
    uu = urllib2.urlopen(req)
    print uu.read()

  # 20190322:ED context manager to override quote_plus with quote - to match new code
  class _urlEncodeWithQuote(object):

    def __init__(self):
      self.should_patch = True

    def __enter__(self):
      if self.should_patch:
        self.original_handler = urllib.quote_plus
        urllib.quote_plus = urllib.quote

    def __exit__(self, *args):
      if self.should_patch:
        urllib.quote_plus = self.original_handler

  # 20190322:ED new upload method
  def _uploadFile(self, serverUser, serverPassword, serverScript, fileData, serverDbRoot, fileStoredAs):
    """New routine to upload a file to the server
    """
    import hashlib
    # import certifi
    # import urllib3

    SERVER_PASSWORD_MD5 = b'c Wo\xfc\x1e\x08\xfc\xd1C\xcb~(\x14\x8e\xdc'

    # early check on password
    m = hashlib.md5()
    m.update(serverPassword.encode('utf-8'))
    if m.digest() != SERVER_PASSWORD_MD5:
      raise Exception('incorrect password')

    auth = base64.encodestring(serverUser + ":" + serverPassword)[:-1]
    authheader = 'Basic %s' % auth

    headers = {'Content-type' : 'application/x-www-form-urlencoded;charset=UTF-8',
               'Authorization': authheader}

    # use the original quote and not quote_plus (option not available in early versions)
    with self._urlEncodeWithQuote():
      body = urllib.urlencode({'fileData': fileData, 'fileName': fileStoredAs, 'serverDbRoot': serverDbRoot}).encode('utf-8')



    # # urllib3.contrib.pyopenssl.inject_into_urllib3()       # not sure if this is needed
    # http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED',
    #                            ca_certs=certifi.where(),
    #                            timeout=urllib3.Timeout(connect=5.0, read=5.0),
    #                            retries=urllib3.Retry(1, redirect=False))
    #
    # try:
    #   response = http.request('POST', serverScript,
    #                           headers=headers,
    #                           body=body,
    #                           preload_content=False)
    #   result = response.read().decode('utf-8')


    try:
      request = urllib2.Request(serverScript,
                                headers=headers,
                                data=body
                                )
      response = urllib2.urlopen(request)
      result = response.read().decode('utf-8')


      if result.startswith(BAD_DOWNLOAD):
        print 'Error reading from server.'
      else:
        return result

    except Exception, es:
      print 'Error reading from server:', str(es)

  # 20190322:ED new download method
  def _downloadFile(self, serverScript, serverDbRoot, fileName):
    """New routine to read a file from the server
    """
    import ssl

    context = ssl._create_unverified_context()
    fileName = os.path.join(serverDbRoot, fileName)

    addr = '%s?fileName=%s' % (serverScript, fileName)
    try:
      response = urllib.urlopen(addr, context=context)

      data = response.read()  # just split for testing
      if isinstance(data, unicode):
        data = data.decode('UTF-8')
      response.close()

      if data.startswith(BAD_DOWNLOAD):
        raise Exception(data[len(BAD_DOWNLOAD):])

      return data
    except:
      return None

  def _openUrl(self, serverScript, serverDbRoot, fileName):
    """Header to open a url
    """
    # try:

    # 20190712:ED urllib3 version with where() certificate - definitely works
    # seems to fix problems on MacOS
    import ssl
    import urllib3.contrib.pyopenssl
    import certifi

    # fileName = os.path.join(serverDbRoot, fileName)
    fileName = '/'.join([serverDbRoot, fileName])

    addr = '%s?fileName=%s' % (serverScript, fileName)

    urllib3.contrib.pyopenssl.inject_into_urllib3()
    proxy_url = os.environ.get('HTTPS_PROXY', os.environ.get('HTTP_PROXY'))
    if proxy_url:
      http = urllib3.ProxyManager(proxy_url,
                                  cert_reqs='CERT_REQUIRED',
                                  ca_certs=certifi.where(),
                                  # timeout=urllib3.Timeout(connect=3.0, read=3.0),
                                  retries=urllib3.Retry(1, redirect=False))
    else:
      http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED',
                                 ca_certs=certifi.where(),
                                 # timeout=urllib3.Timeout(connect=3.0, read=3.0),
                                 retries=urllib3.Retry(1, redirect=False))

    response = http.request('GET', addr,
                            preload_content=False)
    return response

    # except Exception as es:
    #   if self.parent.isGraphical:
    #     self.parent.warningMessage('Warning', str(es))
    #   raise

  def getFileUpdates(self):
  
    # Synchronise from server
    # find available file updates on the server
    # update the local cache
    
    for fileUpdate in list(self.fileUpdates):
      fileUpdate.delete()
    
    try:
      # addr = 'http://' + joinPath(self.url,self.dataFile)
      # url  = urllib.urlopen(addr)

      # 20190322:ED use the new openUrl method
      addr = 'http://' + joinPath(UPDATE_SERVER_LOCATION, UPDATE_SCRIPT)
      url = self._openUrl(addr, UPDATE_BASE_DIRECTORY+self.parent.version, self.dataFile)

    except:
      if self.parent.isGraphical:
        self.parent.warningMessage('Warning','Cannot access update server via network')
      return

    line = None
    try:
      line = url.readline()
      version = line.split()[0]
    except:
      if self.parent.isGraphical:
        self.parent.warningMessage('Warning','Something wrong with update server: version returned "%s"' % line)
      return

    if version[0] not in '0123456789':
      self.parent.warningMessage('Warning','No updates at server location')
      self.version = None
      url.close()
      return
      
    minor_version = ''
    fields = version.split('.')
    if len(fields) == 4:
      minor_version = '.' + fields[3]
    elif len(fields) != 3:
      self.parent.warningMessage('Warning','Version at server location = %s, not of required form' % version)
      self.version = None
      url.close()
      return

    self.version = '.'.join(fields[:3])
    self.minor_version = minor_version

    if self.version != self.parent.version:
      if self.admin:
        self.parent.warningMessage('Warning','Local version %s does not match server version %s' % (self.parent.version, self.version))
     
      else:
        url.close()
        ###self.parent.warningMessage('Warning','You must upgrade to version %s (currently %s) to get updates.' % (self.version, self.parent.version))
        self.parent.installLatestRelease()
        return
     
    
    line = url.readline()
    while line:
      (fileName, filePath, storedAs, language, date, details, priority) = line.split(fieldSep)
      
      fileUpdate = FileUpdate(self, fileName, filePath, language=language, date=date, details=details, priority=int(priority), isNew=False)
      fileUpdate.getCacheFromServer()
      
      line = url.readline()
    
    url.close()
    
  def getSelectedUpdates(self):
  
    updates = []
    for update in self.fileUpdates:
      if update.isSelected:
        updates.append(update)

    return updates

  def delete(self):
  
    for update in self.fileUpdates:
      update.delete()
      
    self.parent.server = None  
    del self

  def _checkPassword(self):
    if self.parent.isGraphical:
      serverPassword = askPassword('Password Request', 'Enter password for %s@%s' % (self.uid, self.location))
      if not serverPassword:
        return False

      return serverPassword

  def _removeFileFromServer(self, fileUpdate, passwd):
    """remove a file from the server database
    """
    if passwd:
      self._uploadFile('ccpn', passwd, 'https://www.ccpn.ac.uk/cgi-bin/updateadmin/__actionFile', fileUpdate.storedAs, self.identity[-1], '')

    else:
      self.parent.warningMessage('Warning', 'Password must be specified when setting updates in non graphical mode')


class FileUpdate:

  def __init__(self, server, fileName, filePath, language='python',
               date=None, isSelected=True, details='', priority=1, isNew=True):

    self.parent     = server
    self.server     = server
    self.fileName   = fileName
    self.filePath   = filePath
    self.language   = language
    self.isSelected = isSelected
    self.isNew      = isNew
    self.details    = details
    self.priority   = priority
    self.date       = date or '%s' % ctime()

    # self.storedAs   = '__temp_'.join(filePath.split('/')) + '_' + fileName

    # 20190322:ED different separator, __temp_ will split into directory tree on the server, __sep_ will put single file in update path
    fullPath = '/'.join([filePath, fileName])
    self.storedAs   = '__sep_'.join(fullPath.split('/'))

    self.isCached   = 0

    agent = server.parent
    self.tempFile      = joinPath(agent.tempDir,    self.storedAs)
    self.installedFile = joinPath(agent.installRoot,self.filePath,self.fileName)
    self.serverFile    = joinPath(server.url,       self.storedAs)
        
    for fileUpdate in server.fileUpdates:
      if fileUpdate.fileName == fileName and fileUpdate.filePath == filePath:
        del self
        return
    
    server.fileUpdates.append(self)

  def timestamp(self):
  
    self.date = '%s' % ctime()  

  def getIsUpToDate(self):
  
    if not path.isfile(self.installedFile):
      return 0
      
    return areFilesIdentical(self.tempFile, self.installedFile)

  def getCacheFromServer(self):
    
    # copy server file to tempFile
    self.isCached = 1
    # addr = 'http://' + self.serverFile
    try:
      # url  = urllib.urlopen(addr)

      # 20190322:ED use the new openUrl method
      addr = 'http://' + joinPath(UPDATE_SERVER_LOCATION, UPDATE_SCRIPT)

      # use this if the separator is '__temp_', uploading this will put into path on the server
      # url = self.parent._openUrl(addr, UPDATE_BASE_DIRECTORY+self.parent.version, os.path.join(self.filePath, self.fileName))

      # use this if the separator is '__sep_', uploading this will put into update/basedirectory
      url = self.parent._openUrl(addr, UPDATE_BASE_DIRECTORY+self.parent.version, self.storedAs)

    except:
      self.server.parent.warningMessage('Download failed','Could not connect to %s' % addr)
      return
    
    data = url.read()

    file = open(self.tempFile, 'w')
    file.write(data)
    file.close
    url.close()

  def setInstallFromCache(self):
  
    try:
      if path.isfile(self.installedFile):
        copyfile( self.installedFile, self.installedFile+'__old' )
      print 'installing', self.installedFile
      dirname = os.path.dirname(self.installedFile)
      if not os.path.exists(dirname):
        os.makedirs(dirname)
      copyfile( self.tempFile,self.installedFile )
    except Exception, e:
      self.server.parent.warningMessage('Copy Fail','Could not update file %s: %s' % (self.installedFile, e))
      return
  
    if self.installedFile.endswith('.py'):
      py_compile.compile(self.installedFile)

  def select(self):
  
    self.isSelected = 1

  def deselect(self):
  
    self.isSelected = 0

  def delete(self):
  
    if self.isCached and path.isfile(self.tempFile):
      os.remove(self.tempFile)
    self.parent.fileUpdates.remove(self)
    del self
    
class ReleaseUpdate:

  def __init__(self, server, installRoot, currentVersion,
               releaseDir = 'temporaryReleaseDir', ccpnmrTopDir = CCPNMRTOPDIR,
               ccpnmrCodeDir = CCPNMRCODEDIR, httpServer = 'www2.ccpn.ac.uk',
               httpDir = CCPNMRTOPDIR, uid = 'ccpn'):

    self.server         = server
    self.installRoot    = installRoot
    self.currentVersion = currentVersion
    self.releaseDir     = releaseDir
    self.ccpnmrTopDir   = ccpnmrTopDir
    self.ccpnmrCodeDir  = ccpnmrCodeDir
    self.httpServer      = httpServer
    self.httpDir         = httpDir
    self.homeDir         = '~' + uid

    self.baseDir     = os.path.dirname(installRoot)
    self.installDir  = os.path.basename(installRoot)
    self.releaseFile = None

  def getLatestRelease(self):

    # assumes is in correct directory (so far)
    # fetches latest release
    # returns name of release file (not full path, just file name)
    
    ss = ''
    use_precompiled = False
    try:
      from ccpnmr.analysis.Version import platform, bits, arch, built_by
      if built_by == 'ccpn':
        if platform == 'linux':
          ss =  '_%s%s' % (platform, bits)
          use_precompiled = True
        elif platform == 'darwin':
          if arch in ('i386', 'intel'):
            arch = 'intel%s' % bits
          else:
            arch = 'ppc'
          ss =  '_%s' % arch
          use_precompiled = True
    except:
      pass

    self.use_precompiled = use_precompiled

    if use_precompiled:
      extension = 'tgz'
    else:
      extension = 'tar.gz'

    fileName = 'analysis%s%s%s.%s' % (self.server.version, self.server.minor_version, ss, extension)
    
    self.releaseFile = None

    addr = 'http://' + joinPath(self.httpServer, self.homeDir, self.httpDir, fileName)
    try:
      url  = urllib.urlopen(addr)
    except:
      self.server.parent.warningMessage('Download failed','Could not connect to %s' % addr)
      return

    data = url.read()
    url.close()

    destFile = joinPath(self.baseDir, self.releaseDir, fileName)
    try:
      fp = open(destFile, 'wb')
    except:
      self.server.parent.warningMessage('Saving release failed','Could not save release to %s' % destFile)

    fp.write(data)
    fp.close()

    self.releaseFile = fileName
  
    if self.server.parent.isGraphical:
      if showOkCancel('Query','Can CCPN log your IP address for its statistics? No other information will be taken'):
        self.logDownload()
    else:
      print 'CCPN is logging your IP address for its statistics. No other information will be taken'
      self.logDownload()
  
  def logDownload(self):
    
    addr = 'http://www2.ccpn.ac.uk/cgi-bin/update/logUserDownload.py.cgi'
    try:
      url = urllib.urlopen(addr)
      url.readline()
      url.close
    except:
      pass
    
  def installRelease(self):

    os.chdir(self.baseDir)
    if os.path.exists(self.releaseDir):
      if not os.path.isdir(self.releaseDir):
        self.server.parent.warningMessage('Failure', 'Script uses %s%s%s but that is not a directory' % (self.baseDir,os.sep, self.releaseDir))
        return False
    else:
      os.mkdir(self.releaseDir)
    os.chdir(self.releaseDir)

    self.getLatestRelease()
    
    if self.releaseFile:
      self.unpackRelease()
      if self.use_precompiled:
        # move all of ccpnmr
        self.moveAllRelease()
        self.compileAllPyCode()
        #self.runConfigScript()
      else:
        # move only ccpnmr/ccpnmrX.Y
        self.compileCCode() # What if this fails?
        self.moveCurrentRelease()
        self.moveNewRelease()
        # below used to be before moveCurrentRelease
        # but that gives misleading error messages (wrong file names)
        # so put it here instead
        self.compilePyCode()
    
    else:
      self.server.parent.warningMessage('Failure','Cannot find file for latest release - version %s' % self.server.version)
      return False

    return True

  def runCmds(self, cmds):

    cmd = ';'.join(cmds)
    os.system(cmd)

  def unpackRelease(self):

    os.chdir('%s/%s' % (self.baseDir, self.releaseDir))

    gzippedTarFile = self.releaseFile

    cmds = []
    if gzippedTarFile.endswith('.tar.gz'):
      tarFile        = gzippedTarFile[:-3]
      sourceEnvFile  = joinPath(self.baseDir, self.installDir, 'c', environmentFile)
      destEnvFile    = joinPath(self.baseDir, self.releaseDir, self.ccpnmrTopDir, self.ccpnmrCodeDir, 'c', environmentFile)
      cmds.append('gunzip %s'  % gzippedTarFile)
      cmds.append('tar xvf %s' % tarFile)
      cmds.append('gzip %s'    % tarFile)
      cmds.append('cp %s %s'   % (sourceEnvFile,destEnvFile))
    else:
      cmds.append('tar xvfz %s' % gzippedTarFile)

    self.runCmds(cmds)

  def compileCCode(self):

    os.chdir('%s/%s/%s/%s/c' % (self.baseDir, self.releaseDir, self.ccpnmrTopDir, self.ccpnmrCodeDir))

    cmds = []
    cmds.append('make')
    ### replace make links code with code at bottom
    ###cmds.append('make links')

    self.runCmds(cmds)

    # make symbolic links
    os.chdir('%s/%s/%s/%s/python' % (self.baseDir, self.releaseDir, self.ccpnmrTopDir, self.ccpnmrCodeDir))
    script = './linkSharedObjs'
    cmds = ['chmod u+x %s' % script, script]
    for directory in ('memops/c', 'ccp/c', 'ccpnmr/c', 'cambridge/c'):
      os.chdir(directory)
      self.runCmds(cmds)
      os.chdir('../..')

  def compilePyCode(self):

    # this used to be done before moving the releases around
    #directory = '%s/%s/%s/%s/python' % (self.baseDir, self.releaseDir, self.ccpnmrTopDir, self.ccpnmrCodeDir)
    os.chdir(self.baseDir)
    directory = '%s/python' % (self.installDir,)

    compileall.compile_dir(directory)

  def moveCurrentRelease(self):

    os.chdir(self.baseDir)
    renameDir = '%s_%s' % (self.installDir, self.currentVersion)
    n = 1
    while os.path.exists(renameDir):
      renameDir = '%s_%s_%d' % (self.installDir, self.currentVersion, n)
      n = n + 1

    print 'About to rename %s to %s' % (self.installDir, renameDir)
    os.rename(self.installDir, renameDir)

  def moveNewRelease(self):

    os.chdir(self.baseDir)
    dd = os.path.join(self.releaseDir, self.ccpnmrTopDir, self.ccpnmrCodeDir)
    print 'About to rename %s to %s' % (dd, self.installDir)
    os.rename(dd, self.installDir)

  def compileAllPyCode(self):

    upDir = os.path.dirname(self.baseDir)  # directory above ccpnmr
    baseDir = os.path.basename(self.baseDir)  # ccpnmr
    os.chdir(upDir)
    compileall.compile_dir(baseDir)

  def moveAllRelease(self):

    upDir = os.path.dirname(self.baseDir)
    baseDir = os.path.basename(self.baseDir)
    os.chdir(upDir)
    renameDir = '%s_%s' % (baseDir, self.currentVersion)
    n = 1
    while os.path.exists(renameDir):
      renameDir = '%s_%s_%d' % (baseDir, self.currentVersion, n)
      n = n + 1

    print 'About to rename %s to %s' % (baseDir, renameDir)
    os.rename(baseDir, renameDir)

    dd = os.path.join(renameDir, self.releaseDir, self.ccpnmrTopDir)
    print 'About to rename %s to %s' % (dd, baseDir)
    os.rename(dd, baseDir)
        
  def runConfigScript(self):

    os.chdir(self.baseDir)
    os.system('python configRelease.py')

