
"""
======================COPYRIGHT/LICENSE START==========================

Io.py: Utility code for CCPN code generation framework

Copyright (C) 2005  (CCPN Project)

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

import os
import shutil
import StringIO
import sys
import traceback
import glob

from memops.universal import Constants

def printError(title, message, parent=None):
  """ Error printing - replaces showError in non-GUI contexts
  """
  print ("%s:\n  %s" % (title, message))
  
def printExceptionAndTrace(fp = sys.stdout):

  sf = StringIO.StringIO()
  exc = sys.exc_info()
  tb = exc[2]
  traceback.print_tb(tb, file=sf)
  ss = sf.getvalue()
  sf.close()
  del sf
  del tb
  fp.write('%s: %s\n%s' % (exc[0].__name__, exc[1], ss))

def getExceptionAndTrace():

  sf = StringIO.StringIO()
  printExceptionAndTrace(fp=sf)
  ss = sf.getvalue()
  sf.close()

  return ss

# Note: makeAbsolute is Boolean so the default should be
# False but ObjectDomain does not understand True/False
def normalisePath(path, makeAbsolute=None):
  """
  Normalises the path, e.g. removes redundant .. and slashes and
  makes sure path uses '/' rather than '\\' as can happen on Windows.
  """

  dirsep = Constants.dirsep
  assert dirsep == '/', "have dirsep = '%s', not '/'" % dirsep

  # note, cannot just use os.sep below because can have window file names cropping up on unix machines
  winsep = '\\'

  if os.sep == winsep:
    path = path.replace(dirsep, winsep)
    path = os.path.normpath(path)
    path = path.replace(winsep, dirsep)
  else:
    path = path.replace(winsep, dirsep)
    path = os.path.normpath(path)

  if makeAbsolute:
    if path and path[0] != dirsep:
      path = joinPath(os.getcwd(), path)

  return path

def unnormalisePath(path):
  """
  On Unix does nothing, on Windows replaces '/' with '\\'/
  """

  if (Constants.dirsep != os.sep):
    path = path.replace(Constants.dirsep, os.sep)

  return path

def joinPath(*args):
  """
  The same as os.path.join but normalises the result.
  """

  return normalisePath(os.path.join(*args))

def splitPath(path):
  """
  The same as os.path.split but with normalisation taken into account.
  """

  (head, tail) = os.path.split(unnormalisePath(path))
  head = normalisePath(head)

  return (head, tail)

def converseSplitPath(path):
  """
  Similar to splitPath but with head being the top directory and tail being the rest.
  """

  pair = splitPath(path)
  head = None
  tail = ''
  while pair[0] and pair[0] not in ('.', '/'):
    head = pair[0]
    tail = joinPath(pair[1], tail)
    pair = splitPath(head)

  if head is None:
    (head, tail) = pair

  return (head, tail)

def getTopDirectory():

  """
  Returns the 'top' directory of the whole code and data structure.
  """
  
  func = os.path.dirname
  
  return func(getPythonDirectory())

def pathIsInDirectory(path, directory):
  """Does path point to a file that is inside directory, possibley indirectly
  NB does not check whether either exists
  """
  ss = normalisePath(path)
  directory = normalisePath(directory)
  if ss.startswith(directory):
    lendir = len(directory)
    while len(ss) > lendir:
      ss,junk = os.path.split(ss)
    return (ss == directory)
  else:
    return False
    


def setTopDirectory(newTopDirectory):

  """
  Set a new 'top' directory of the whole code and data structure.
  This is so that you can swap between different implementations.
  WARNING: This function should only be used by experts.
  """

  # There is no perfect way at getting at keys in sys.modules which
  # need to be deleted, because there are some funny keys for which
  # sys.modules[key] is None (otherwise you could use the __file__
  # attribute of the module).  So look at the directories in
  # the old python directory and see which of those have an __init__.py
  # inside and assume that those are the modules we are trying to remove
 
  newTopDirectory = normalisePath(newTopDirectory)
  oldTopDirectory = getTopDirectory()
  oldPythonDirectory = getPythonDirectory()
  # Below is a bit dangerous; could have used os.path.basename(oldPythonDirectory) instead
  # if could guarantee that python directory is only one directory below top directory
  n = len(oldTopDirectory)
  newPythonDirectory = newTopDirectory + oldPythonDirectory[n:]

  moduleNames = []
  files = os.listdir(oldPythonDirectory)
  for file in files:
    fullfile = joinPath(oldPythonDirectory, file)
    if os.path.isdir(fullfile):
      if '__init__.py' in os.listdir(fullfile):
        moduleNames.append(file)
 
  # Below does not necessarily find all paths we want because symbolic
  # links mean that sometimes the wrong names are stored in sys.paths

  paths = []
  for path in sys.path:
    path = normalisePath(path)
    if path == oldTopDirectory or path.startswith(oldTopDirectory + Constants.dirsep):
      paths.append(path)

  sys_modules = sys.modules
  sys_path = sys.path

  keys = sys_modules.keys()
  for key in keys:
    for name in moduleNames:
      if key == name or key.startswith(name + '.'):
        del sys_modules[key]
        break

  for path in paths:
    sys_path.remove(path)

  # Now that have deleted what we can from sys.modules and sys.paths, put in
  # new information in sys.paths about newTopDirectory

  sys_path.insert(0, newPythonDirectory)

def getPythonDirectory():

  """
  Returns the 'top' python directory, the one on the python path.
  """
  
  func = os.path.dirname
  
  return func(func(func(os.path.abspath(Constants.__file__))))

def removePath(path):
  """Removes path whether file or directory, taking into account whether symbolic link.
  """

  if not os.path.exists(path):
    return

  if os.path.isfile(path) or os.path.islink(path):
    os.remove(path)
  else:
    shutil.rmtree(path)

def parentDirectory(directory):
  """parent directory
  """
 
  (parent_dir, tail_dir) = splitPath(directory)
 
  if (not tail_dir): # means have trailing slash, so do it without that
    parent_dir = splitPath(parent_dir)[0]
 
  return parent_dir
 

def writePythonDir(dirname):
  """ make python directory dirname (unless it exists already)
  and write blank __init__.py into it, if none exists.
  """
 
  # make directory
  if not os.path.exists(dirname):
    os.mkdir(dirname)
 
  # write __init__ file
  initfile = joinPath(dirname,'__init__.py')
  if not os.path.exists(initfile):
    open(initfile,'w').write("pass\n")

def importModule(name):
  """ Imports module in same way as __import__ but for modules of
  form x.y returns x.y not x (as __import__ does)
  """
 
  m = __import__(name)
 
  comps = name.split('.')
  for comp in comps[1:]:
    m = getattr(m, comp)
 
  return m

def returnRealPath(fileName):

  if os.path.islink(fileName):
    fileName = os.path.realpath(fileName)
  
  return fileName

def convertStringToFileName(fileNameString):
  
  from memops.general.Constants import validCcpnFileNameChars
  from memops.general.Constants import defaultFileNameChar
  
  return doConvertStringToFileName(fileNameString, validCcpnFileNameChars,
                                   defaultFileNameChar)


def doConvertStringToFileName(fileString, validChars, defaultChar):
  
  ll = [x for x in fileString]
  for ii,char in enumerate(ll):
    if char not in validChars:
      ll[ii] = defaultChar
  #
  return ''.join(ll)

def commonSuperDirectory(*fileNames):
  """ Find lowest directory that contains all files in list 
  NB does not normalise file names.
  
  Input: a list of file names
  Output: lowest directory that contains all files. Does *not* end with a file
  """
  return os.path.dirname(os.path.commonprefix(fileNames))

def suggestFileLocations(fileNames, startDir=None):
  """ From a list of files, return a common superdirectory and a list of
  relative file names. If any of the files do not exist, search for an 
  alternative superdirectory that does contain the set of relative file names.
  Searches in either a auperdirectory of the starting/current directory,
  or in a direct subdirectory.  
  
  Input: list of file names
  Output: Superdirectory, list of relative file names.
          If no suitable location is found, superdirectory is returned as None
  """
  
  if not fileNames:
    return None, []
  
  # set up startDir
  if startDir is None:
    startDir = os.getcwd()
  startDir = normalisePath(startDir, makeAbsolute=True)
  
  # find common baseDir and paths
  files = [normalisePath(fp, makeAbsolute=True) for fp in fileNames]
  baseDir = commonSuperDirectory(*files)
  prefix = os.path.join(baseDir, '')
  lenPrefix = len(prefix)
  paths = [fp[lenPrefix:] for fp in files]
  
  if [fp for fp in files if not os.path.exists(fp)]:
    # some file not found. 
    
    # look in superdirectories
    tail = 'junk'
    baseDir = startDir
    while tail:
      for path in paths:
        fp = os.path.join(baseDir, path)
        if not os.path.exists(fp):
          # a file is not found. try new baseDir
          break
      else:
        # this one is OK - stop searching
        break
      #
      baseDir, tail = os.path.split(baseDir)
    
    else:
      # No success - try in a subdirectory (one level) of startDir
      matches = glob.glob(os.path.join(startDir, '*', paths[0]))
      for aMatch in matches:
        baseDir = normalisePath(aMatch[:-len(paths[0])])
        for path in paths:
          fp = os.path.join(baseDir, path)
          if not os.path.exists(fp):
            # a file is not found. try new baseDir
            break
        else:
          # this one is OK - stop searching
          break
      else:
        # we give up
        baseDir = None
  #
  return (baseDir, paths)
    
  
class RedirectToString:

    def __init__(self):
      self.string = ""

    def write(self, s):
      self.string += s
        
    def close(self):
      pass
