"""
======================COPYRIGHT/LICENSE START==========================

headers.py: license header generation code for CCPN framework

Copyright (C) 2005 Rasmus Fogh (CCPN Project)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/GPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
 
You should have received a copy of the GNU General Public
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
import types
import os
import sys
import re

from memops.universal import Io as uniIo

from memops.general.license.licenses import licenses
from memops.general.license.references import references
from memops.general.license.formats import formats
from memops.general.license.contacts import stdContacts

infoFileName = '_licenseInfo'

# marker lines
licenseStartLine =   "======================COPYRIGHT/LICENSE START=========================="
licenseEndLine =     "======================COPYRIGHT/LICENSE END============================"

referenceStartLine = "===========================REFERENCE START============================="
referenceEndLine =   "===========================REFERENCE END==============================="

separatorLine =      "======================================================================="

headerStartLine = licenseStartLine + '\n'
headerEndLine = referenceEndLine + '\n'

# global filtering parameters. These are always skipped
globalExcludeDirs = ['local','CVS']
#
globalExcludeFiles = ['_licenseInfo.py', '__init__.py']
#
globalExcludeFileMatches = [
 '.+\.pyc$','.+\.pyo$','.+\.jar$','.+\.class$',
 '.+\.o$','.+\.so$', '.+\.exe$', '^\..*'
]
for ii in range(len(globalExcludeFileMatches)):
  globalExcludeFileMatches[ii] = re.compile(globalExcludeFileMatches[ii])



# header and copyright holder

copyrights =  """
%(fileName)s: %(programFunction)s

Copyright (C) 2007 %(author)s (%(organization)s)
"""


# contact information

useContacts = """

To obtain more information about this %(programType)s, contact:

"""


# references

useReferences = """
If you are using this software for academic purposes, we suggest
quoting the following references:
"""
    

def setLicenses(
 directory=None, licenseDir='license', warnSkippedFiles=(not 0), mode='act'
):
  """ wrapper - set license header for 'directory' and all subdirectories
   
  licenseDir is the relative path from 'directory'
  to the directory with license texts
  """
  
  # necessary as the first directory may be a relative directory name, 
  # playing havoc with import commands
  firstLookupDir = sys.path[0]
  sys.path[0] = ''
  
  try:
    global testInfo
    testInfo = {}
 
    if directory is None:
      directory = uniIo.getTopDirectory()
    else:
      directory = os.path.abspath(directory)
 
    doSetLicenses(
     directory, licenseDir, None, 0, warnSkippedFiles, mode
    )
 
    if mode == 'test':
      print showTestInfo()
  
  finally:
    sys.path[0] = firstLookupDir


def doSetLicenses(curDir, licenseDir, infoModule, level, warnSkippedFiles, mode):
  """ set license header for a directory and all subdirectories 
  
  licenseDir is the relative path from the level 0 directory
  to the directory with license texts
  """
  print 'Checking...', curDir
  
  # set-up - handle directory changes
  olddir = os.getcwd()
  os.chdir(curDir)
  
  # get dir contents
  contents = os.listdir(curDir)
  
  # set-up - get license info
  if '%s.py' % infoFileName in contents:
    if sys.modules.has_key(infoFileName):
      del sys.modules[infoFileName]
    infoModule = __import__(infoFileName)
    setattr(
     infoModule,'moduleLocation',os.path.join(curDir,'%s.py' % infoFileName)
    )
    
  licenseLocation =  '/'.join(level*['..'] + [licenseDir])
  
  # set-up - prepare exclude patterns
  if infoModule is None:
    print 'No active _licenseInfo'
    infoRange = ()
    excludeDirs = globalExcludeDirs 
    includeDirs = ()
    excludeFiles  = 'all'
    
  else:
    infoRange = range(len(infoModule.licenseInfo))
    if infoModule.excludeFiles == 'all':
      excludeFiles = 'all'
    else:
      excludeFiles = globalExcludeFiles + infoModule.excludeFiles
    
    ll = [re.compile(x) for x in infoModule.excludeFileMatches]
    excludeFileMatches = globalExcludeFileMatches + ll
    if infoModule.includeDirs:
      includeDirs = infoModule.includeDirs
      excludeDirs = []
      if infoModule.excludeDirs:
        raise (
         "Error in %s/%s.py: \nYou cannot have both includeDirs and excludeDirs"
         % (curDir,infoFileName)
        )
    else:
      includeDirs = []
      if infoModule.excludeDirs == 'all':
        excludeDirs = 'all'
      else:
        excludeDirs = globalExcludeDirs + infoModule.excludeDirs
  # separate directory contents, filtering on exclude patterns
  dirs = []
  files = []
  for ss in contents:
    if os.path.isdir(ss):
      if excludeDirs != 'all':
        if ss not in excludeDirs:
          dirs.append(ss)
    elif os.path.isfile(ss):
      # NB this test skips symbolic links, which is what we want
      if excludeFiles != 'all':
        if ss not in excludeFiles:
          for expr in excludeFileMatches:
            if expr.match(ss):
              break
          else:
            files.append(ss)
  # reset dirs in case where includeDirs not empty
  if includeDirs:
    dirs = [x for x in includeDirs if x in dirs]
  
  files.sort()
  dirs.sort()  
  
  # Process files
  for ii in infoRange:
    dict = infoModule.licenseInfo[ii]
    oldFiles = files
    files = []
    dict['licenseLocation'] = licenseLocation
    
    infi = dict.get('includeFiles')
    if infi == 'all':
      for ss in oldFiles:
        addHeader(curDir, ss, ii, infoModule, mode)
      del dict['licenseLocation']
      break
      
    else:
      includeFileMatches = [
       re.compile(x) for x in dict['includeFileMatches']
      ]
      for ss in oldFiles:
        if ss in infi:
          addHeader(curDir, ss, ii, infoModule, mode)
        else:
          for expr in includeFileMatches:
            if expr.match(ss):
              addHeader(curDir, ss, ii, infoModule, mode)
              break
          else:
            files.append(ss)
      del dict['licenseLocation']
  
  # Process not-found files
  if files and warnSkippedFiles:
    for ss in files:
      print " No match found for %s" % os.path.join(curDir,ss)
  
  
  # Process directories:
  for ss in dirs:
    ff = os.path.join(curDir,ss)
    doSetLicenses(ff, licenseDir, infoModule, level+1, warnSkippedFiles, mode)
  
  
  # clean-up
  os.chdir(olddir)


def addHeader(directory, fName, paramIndex, infoModule, mode):
  """ set or replace license header in file fName from directory
  'directory', using parameters 'paramIndex' from 'infoModule'.
  """
  parameters = infoModule.licenseInfo[paramIndex].copy()
  
  # set-up - get format tag, fileName, and file contents
  parameters['fileName'] = fName
  
  formatTag = parameters.get('format')
  if not formatTag:
    formatTag = getFormat(fName)
  formatDict = formats[formatTag]
  commentStart = formatDict['commentStart']
  commentEnd = formatDict['commentEnd']
  
  fileName = os.path.join(directory, fName)
  lines = open(fileName).readlines()
  nlines = len(lines)
  if nlines == 0:
    return
  
  # look for pre-existing header
  try:
    headerStart = lines.index(headerStartLine)
  except:
    headerStart= None
  try:
    headerEnd = lines.index(headerEndLine) + 1
  except:
    headerEnd= None
  
  # temporaryhack
  for ii in range(len(lines)):
    if lines[ii][:len(licenseStartLine)] == licenseStartLine:
      headerStart = ii
      break
  for ii in range(len(lines)):
    if lines[ii][-len(headerEndLine):] == headerEndLine:
      headerEnd = ii + 1
      
  
  if headerStart is None and headerEnd is None:
    # temporary hack to find cases
    # where entire header is collapsed on a single line.
    for ii in range(len(lines)):
      line = lines[ii]
      if (
       line.find(licenseStartLine) == 0 and 
       line.rfind(referenceEndLine) == (len(line) - len(referenceEndLine) -1)
      ):
        headerStart = ii
        headerEnd = ii + 1
    
    
  # add headers:
  if (headerStart is not None and headerEnd is not None and 
      headerStart < headerEnd
  ):
    # there is a previous header - replace it
    
    if lines.count(headerStartLine) > 1 or lines.count(headerEndLine) > 1:
      print ("WARNING, file %s has (partly) multiple license headers"
       % fileName
      )
    
    # check for programFunction to be retained
    # programFunction is in the first non-empty line after the start,
    # the part after the first ':'line
    if not parameters['programFunction']:
      # look for programFunction in the old file
      
      for line in lines[headerStart+1:headerEnd]:
        if not line or line.isspace():
          continue
        elif line[:9] == 'Copyright':
          break
        else:
          ll = line.split(':',1)
          if len(ll) == 2:
            parameters['programFunction'] = ll[1].strip()
          break
    lines[headerStart:headerEnd] = [getHeader(**parameters)]
  
  else:
    
    if not (headerStart is None and headerEnd is None):
      # mangled headers, print warning
      print ("WARNING, file %s has partial license header"
       % fileName
      )
  
    # get header
    parameters['fileName'] = fName
    header = getHeader(**parameters)
      
    # add new header at top
    ii = -1
    while (not 0):
      ii = ii + 1
      line = lines[ii]
      
      # skip header lines and empty lines
      for expr in formatDict['ignoreStartLines']:
        if expr.match(line):
          if ii == (nlines - 1):
            return
          break
          
      else:
        # first non-skipped line. License header goes here
        
        startIndex = line.find(commentStart)
        if startIndex == -1:
          # first real line is not a comment - insert new comment
          lines[ii] = '\n'.join([commentStart, header, commentEnd, line])
          break
          
        else:
          # line is a comment start. find comment end
          
          # first look in current line
          offset =  startIndex + len(commentStart)
          ind = line[offset:].find(commentEnd)
          if ind != -1:
            ind = ind + offset
            
          # now look in following lines
          while ind == -1:
            ii = ii + 1
            if ii >= nlines:
              raise (" ERROR, file %s contains unclosed header comment"
               % fileName
              )
            line = lines[ii]
            ind = line.find(commentEnd)
          
          # found it. Insert header before comment end
          lines[ii] = '\n'.join([line[:ind], header, line[ind:]])
          break
    
    
  # final actions depend on mode
  if mode == 'act':
    # write file back out
    fp = open(fileName,'w')
    fp.writelines(lines)
    fp.close()
 
  elif mode == 'test':
    # output test information
 
    dd = testInfo.get(infoModule)
    if dd is None:
      dd = testInfo[infoModule] = {}
 
    ll = dd.get(paramIndex)
    if ll is None:
      ll = dd[paramIndex] = []
 
    ll.append((fName, formatTag, directory))
 
  else:
    raise "ERROR, illegal mode parameter %s" % `mode`

def getFormat(fName):
  """ find corrent format
  """
  for tag,dict in formats.items():
    for expr in dict['includeFileMatches']:
      if expr.match(fName):
        return tag
  #
  # default format is 'text'
  return 'text'


def showTestInfo():
  """ return formatted string with test info summary
  
  Using global testInfo dictionary for input
  """
  
  strings = []
  strapp = strings.append
  
  # reorder and sort testInfo
  ll0 = []
  for infoMod,dict in testInfo.items():
    ll0.append((infoMod.moduleLocation,infoMod,dict))
  ll0.sort()
  
  # prepare string
  for location,infoMod,dict in ll0:
    strapp(72*'#')
    strapp('#')
    strapp('# %s' % location)
    strapp('#\n')
    
    ll1 = dict.items()
    ll1.sort()
    for paramIndex, ll2 in ll1:
      strapp('# licenseInfo set %s :\n{' % paramIndex)
      ll3 = infoMod.licenseInfo[paramIndex].items()
      ll3.sort()
      for tag,val in ll3:
        strapp( "'%s' : '''%s'''," % (tag,val))
      strapp('}')
      strapp('# File, format, directory :')
      for tt in ll2:
        strapp("  %s, %s, %s" % tt)
     
    strapp('\n\n')
   
  #
  return '\n'.join(strings)

def getHeader(**parameters):
  """ function to produce header. Parameters are (all keyword) :
 
  mandatory :
  author, organization, useLicense, licenseLocation, fileName, 
 
  optional:
  programFunction, programType, stdContact, extraContact, references, credits
  """

  # defaults:
  if parameters.get('programType') is None:
    parameters['programType'] = 'Program'
    
  programFunction = parameters.get('programFunction')
  if programFunction is None:
    parameters['programFunction'] = ''
  elif len(programFunction.splitlines()) > 1:
    raise "Error in getHeader: programFunction must be a single line"
  
  if not (parameters.get('stdContact') or parameters.get('extraContact')):
    raise "Error in getHeader: either stdContact or extraContact must be given"

  # make header
  strings = []
  strapp = strings.append
  
  # license and author
  strapp(licenseStartLine)
  strapp(copyrights % parameters)
  strapp(separatorLine)
  template = licenses[parameters['useLicense']]
  strapp(template % parameters)
  strapp(licenseEndLine)
  
  # contact info
  ss = parameters.get('stdContact')
  if ss:
    strapp(stdContacts[ss])
  ss = parameters.get('extraContact')
  if ss:
    strapp(ss)
  strapp(separatorLine)
  
  # references and credits
  refs = parameters.get('references')
  if refs:
    strapp(useReferences)
  strapp(referenceStartLine)
  if refs:
    if type(refs) == types.StringType:
      refs = (refs,)
    for ref in refs:
      ss = references.get(ref)
      if ss is not None:
        strapp(ss)
      else:
        strapp(ref)
      strapp('')
  ss = parameters.get('credits')
  if ss:
    strapp(ss)
  strapp(referenceEndLine)
  
  #
  strapp('')
  return '\n'.join(strings)

if __name__ == '__main__':
  
  if sys.argv[1:]:
    directory = sys.argv[1]
  else:
    directory = None
    
  setLicenses(
   directory=directory,
   mode='test'
  )
