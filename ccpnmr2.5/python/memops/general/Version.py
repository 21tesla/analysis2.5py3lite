"""Constants for versioning
======================COPYRIGHT/LICENSE START==========================

Version.py: code for CCPN data model and code generation framework

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

import time, os

# Maps memops.general.Constants modelVersion to repository location
# Used with function getRepositoryDir below - see there for usage.
#
# 'modelVersion': (repositoryCode', (repository location list)
versionMap = {
 # history versions
 '1.1.a3': ('cvs', ('branch4',) ),
 '2.0.a0': ('cvs', ('branch_2_0_3',) ),
 '2.0.a1': ('cvs', ('stable_2_0_4',) ),
 '2.0.a2': ('cvs', ('stable_2_0_5',) ),
 '2.0.a3': ('cvs', ('stable_2_0_6',) ),
 '2.0.b1': ('cvs', ('stable_2_0_7',) ),
 '2.0.b2': ('cvs', ('stable_2_0_8',) ),
 
 '2.0.4' : ('cvs', ('stable_2_2_0',) ),
 '2.0.b3': ('svn', ('tags', 'test2.2.2_stable_A') ),
 '2.0.5' : ('svn', ('branches', 'model_2_0_5')),  # NB NOT the same as branch model2_0_5 (obsolete)
 '2.0.6' : ('svn', ('branches', 'stable_20131212')), 
 '2.1.0' : ('svn', ('tags', 'model_2_1_0') ),
 '2.1.1' : ('svn', ('tags', 'model_2_1_1') ),  # NBNB present only on rhf22 computer. FIXNE!
 #'jmci' : ('svn', ('tags', 'jmci') ),   # Temporary - JMCI python versions
 #'merge' : ('svn', ('branches', 'merge_2_1_2')),  # Temporary locatoin for 2.1.2 stable/trunk merge
 
 # Current stable/trunk versions
 's'     : ('svn', ('branches', 'stable',) ),       # stable
 't'     : ('svn', ('trunk',) ),                    # trunk
 #'marc'     : ('svn', ('branches', 'FEmarcvDijk',) ),       # Marc new xml model version
}
# Synonyms for stable/trunk, used by backwards compatibility code
versionMap['2.1.2'] = versionMap['s']
#versionMap['2.1.2'] = versionMap['t']

cvsWorkingDir = 'ccpn'

#svnWorkingDir = 'work'
svnWorkingDir = 'ccpn'

class Version:

  def __init__(self, major=0, minor=0, level='', release=0, 
               timestamp=None, name='Test'):

    if (timestamp is None):
      timestamp = time.ctime()

    self.major = int(major)
    self.minor = int(minor)
    self.level = level
    self.release = int(release)
    self.timestamp = timestamp
    self.name = name

  def __repr__(self):

    return '%s.%s.%s%s' % (self.major, self.minor, self.level, self.release)

  def __cmp__(self, other):
    
    if not isinstance(other, Version):
      return cmp(id(self), id(other))

    result = cmp( (self.major, self.minor), (other.major, other.minor) )

    if not result:
      result = cmp(self.level, other.level)

      if result:
        if '' in (self.level, other.level):
          # NB empty string compares larger than letters - special case
          result = -result
      
      else:
        result = cmp(self.release, other.release)
    
    return result
    
    
  def __hash__(self):
    return hash(('__!@#$%%memops.general.Version.Version', 
                 self.major, self.minor, self.level, self.release))
  
#versionDict = {}
#for key in versionTupleDict.keys():
#  versionDict[key] = Version(name=key, *versionTupleDict[key])
#
#def findVersion(key):
#  """Find version given key."""
#
#  return versionDict.get(key)

def getVersion(s=None, timestamp = None, name = None):
  """Get version given string representation in form '1.0.b3'.
  Amplified to handle strings of form 1.0b3, 
  for backwards compatibility. Rasmus Fogh 3May04"""
  
  if s:
    (major, minor, level, release) = parseVersionString(s)
    return Version(major, minor, level, release, timestamp=timestamp, name=name)
    
  else:
    # special case - get default version
    result = Version()
    result.timestamp = timestamp
    result.name = name
    #
    return result
  
  
def parseVersionString(s):
  """ Parse version string into (major, minor, level, release) tuple
  """
  ll = s.split('.')
  
  if len(ll) == 3:
  
    (major, minor, rest) = ll
  
  elif len(ll) == 2:
  
    (major, rest) = ll
    
    n = 0
    while (rest[n] in '0123456789'):
      n = n + 1
      
    minor = rest[:n]
    rest = rest[n:]
  
  else:
    raise ValueError("Invalid argument to getVersion: %s" % `s`)
  
  # finish processing and return result
  if rest[-1] in '0123456789':
  
    n = 0
    while (rest[n] not in '0123456789'):
      n = n + 1
    
    level = rest[:n]
    release = rest[n:]
  
  else:
    # version of form '1.0.b'. Treat as release 0
    level = rest
    release = 0
  
  return (major, minor, level, release)
  


def cmpVersionStrings(str1, str2):
  """ Compare version strings as versions
  """
  tpl1 = parseVersionString(str1)
  tpl2 = parseVersionString(str2)

  result = cmp( tpl1[:2], tpl2[:2] )

  if not result:
    ll = [tpl1[2], tpl2[2]]
    result = cmp(*ll)

    if result:
      if '' in ll:
        # NB empty string compares larger than letters - special case
        result = -result
        
    else:
      result = cmp(*[tpl1[3], tpl2[3]])
  
  return result


def getRepositoryDir(versionTag, repoTag=None):
  """ Get repository directory. Used for compatibility code generation scripts
  (CompatibilityGen) and optionally for repository navigation scripts. 
  Should *NOT* be used in released code, as it makes assumptions about the
  code repository structure.
  
  Assumes that code trees are rooted in
  for CVS:
   $CVSROOT/extraDirs/ccpn (for model only)
  for SVN:
   $CCPN_SVNROOT/work/extraDirs/ccpn
  
  extraDirs is one or more directories as given in the versionMap
  """
  
  if versionTag == '3.0.a1':
    # special case - different directory structure
    return os.path.join(os.environ.get('CCPN_MODEL'), 'trunk', 'ccpn', 'ccpnmodel', 'v_3_0_a1')
  
  progCode, extraDirs = versionMap[versionTag]
  
  if not repoTag:
    if progCode == 'cvs':
      repoTag = cvsWorkingDir
    else:
      repoTag = svnWorkingDir
  
  if progCode == 'cvs':
    # get cvs directory
    ll = [os.environ.get('CVSROOT')]
    ll.extend(extraDirs)
    ll.append(repoTag)
  
  else:
    # get svn directory

    ll = [os.environ.get('CCPN_SVNROOT'), repoTag]
    ll.extend(extraDirs)
    ll.append('ccpn')
  #
  return os.path.join(*ll)
