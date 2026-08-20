""" Python XML I/O

Data XML I/O code.

Works off XML I/O maps generated and written in memops/scripts_v2/xmlio.
See there for tag names and I/O map organisation.
For every abstract type, attribute and link there is a map dictionary that 
describes how the object should be handled. These dictionaries can be 
accessed either from the XML tag (therough loadMaps) or from the attribute 
or type name. The tag to write and the type to create or attribute to set
are kept inside the map; for version compatibility purposes the XML tag and
attribute/type name may reflect different versions.

Uses ElementTree for parsing file when reading, but writes the files directly.
The only characters escaped are '<', '>', '&', and, in attributes, '"'.

Object are created as empty shells(invalid) and attributes are set after
the fact. A number of implementation attributes (e.g. parent-child links, 
TopObjects bookkeeping) have to be set explicitly after the fact. Intra-file 
links are dereferenced after the file has been read.


There is an element or attribute for every class object, and DataObjType object.
These elements have attributes and contain  other elements  but no text. 
Out-of-package links have their own elements, that have no attributes but
contain other elements that make up, in order, the full key for the linked-to 
object. Note that the exolink tags determine the class of the element linked to,
so that a given out-of-package link may sontain objects of different classes
with different keys.
Simple data type are mostly incorporated in their containing attributes, 
but may be handled as elements with no attributes containing text only.

Attributes and links may be stored in three ways:

A) If hicard == 1 they may be put into attributes of the element that defines
the containing object. For text types this is used only for types Word and
Token, for other simple types and for intra-package links this is used whenever
possible. If the file is written with compact=False (default is compact=True)
these data are insted written according to B.

B) Intra-package links and simple type attributes that can not contain white
space are  written as an XML element per attribute containing a (white-space
separated) XML list of values. This is used when hicard != 1 (if compact=True)
or in all cases (if compact=False)

C) Attributes that may contain whitespace (most text types and DataObjTypes),
out-of-package links, and child links are written as elements without 
attributes that contain other elements. Individual out-of-package links (which 
are semantically different) are written the same way. The contents may be
class elements (for child links), DataObjType elements, SImple type elements, 
or out-of-package link elements.


Every Class element has an attribute '_ID' that gives an arbitrary object ID of
the form '_int'. Intra-package links use these IDs as object pointers.  The IDs
depend on the order in which objects are encountered while saving; it follows
that different data files cannot be compared with 'diff', because otherwise
identical dat may be using different IDs.
/home/rhf22/CCPN/cvsroots/trunk/all/python/memops/format/xml/XmlIO.py


By default all two-way intra-package links are stored at both ends, and
attributes are stored also if they are identical to the default value. 
Alternatively data can be written with the switch 'simplified=True'. In this
case attribtues are omitted if they are identical to the default value, and
one-to-many links are written only on the 'many' side (i.e. as a single list 
of IDs appearing in a single object). There is no provision for storing 
many-to-many links at one end only - this would be possible for unordered
many-to-many links, although not for ordered ones.

Inter-package links are stored on one side only - the side of the importing 
package


Attributes or links that are unset, i.e. have the value None, or are empty,
are not stored.


The code has been partially prepared for storing and loading subtrees, i.e. 
trees that are not headed by a TopObject. Some more changes would still be 
needed, though. In order todo this you would need special I/O maps, so that
out-of-tree links were written as out-of-package links.


======================COPYRIGHT/LICENSE START==========================

XmlIO.py: Generic Xml parser and I/O

Copyright (C) 2003 Rasmus Fogh (CCPN project)
 
=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license.
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

======================COPYRIGHT/LICENSE END============================

To obtain more information about this code:

- CCPN website (http://www.bio.cam.ac.uk/nmr/ccp/)

- contact Rasmus Fogh (ccpn@bioc.cam.ac.uk)

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following reference:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.
===========================REFERENCE END===============================
 
"""

import os, gc

from memops.general.Implementation import ApiError
from memops.format.xml import Util as xmlUtil
from memops.api import Implementation as ImplementationApi

implPrefix = 'IMPL'
originatorString = "CCPN Python XmlIO"

indentBySpaces = 2
indents = {0:''}

# NBNB TBD version 1.1 statt 1.0?
xmlHeader1 = """<?xml version="1.0" encoding="UTF-8"?>
"""

#"""
xmlHeader2 = """<_StorageUnit time="%s" release="%s" 
 packageGuid="%s" originator="%s">
 
"""

xmlFooter = """
</_StorageUnit>
<!--End of Memops Data-->
"""

def save(repositoryPath, topObject, mapping=None, 
         comment=None, simplified=True, compact=True, expanded=False): 
  
  fileLocation = xmlUtil.findTopObjectPath(repositoryPath, topObject)
  
  if os.path.exists(fileLocation):
    # File that fits guid already exists - overwrite it
    pass
  
  else:
    # file does not exist - make sure the directories are there.
    dirName = os.path.dirname(fileLocation)
    if not os.path.exists(dirName):
      os.makedirs(dirName)
    elif not os.path.isdir(dirName):
      raise ApiError("%s exists and is not a directory" % dirName)
  
  stream = open(fileLocation, 'w')
  try:
    saveToStream(stream, topObject, mapping, comment, simplified, compact)
  finally:
    stream.close()
    
def saveToStream(stream, topObject, mapping=None, 
         comment=None, simplified=True, compact=True, expanded=False): 
  """ wrapper function, to handle garbage collection for Xml export.
  """
   
  from memops.xml import Implementation
   
  if gc.isenabled():
    gc.disable()
    try:
      Implementation.saveToStream(stream, topObject, mapping, comment, 
                                  simplified, compact, expanded)
    finally:
      gc.enable()
  else:
    Implementation.saveToStream(stream, topObject, mapping, comment, 
                                simplified, compact, expanded)


def loadProject(repositoryPath, projectName, partialLoad=False):
  """ Load project by projectName
  """
  filePath = xmlUtil.getProjectFile(repositoryPath, projectName)
  return loadProjectFile(filePath, partialLoad=partialLoad)


def loadProjectFile(filePath, partialLoad=False):
  """ Load project by file name
  """
  topObjId = xmlUtil.getTopObjIdFromFileName(filePath, mustBeMultipart=False)
  
  filePath = os.path.abspath(filePath)
  lowestDir = os.path.dirname(filePath)
  
  project = loadFromStream(open(filePath), topObjId=topObjId, 
                           partialLoad=partialLoad)
  
  activeRepositories = project.__dict__['activeRepositories']
  if not activeRepositories:
    for repository in project.packageLocator.repositories:
      path = repository.getFileLocation('memops.Implementation')
      if path == lowestDir:
        activeRepositories.append(repository)
        break
  #
  return project


def loadTopObject(repositoryPath, topObject):
  """ Load given topObject 
  """
  filePath = xmlUtil.findTopObjectPath(repositoryPath, topObject)
  topObjId = xmlUtil.getTopObjIdFromFileName(filePath, mustBeMultipart=True)
  
  return loadFromStream(open(filePath), topObjId=topObjId, topObject=topObject)


def loadFromFile(memopsRoot, filePath, partialLoad=False):
  """ Load topObject from filePath and MemopsRoot
  """
  
  # open file
  stream = open(filePath)
   
  # load file
  try:
    topObjId = xmlUtil.getTopObjIdFromFileName(filePath, mustBeMultipart=True)
    topObject = loadFromStream(stream, topObject=memopsRoot,
                                     topObjId=topObjId,
                                     partialLoad=partialLoad)
  finally:
    stream.close()
  
  return topObject


def loadFromStream(stream, topObjId=None, topObject=None, partialLoad=False):
  """ Wrapper function, to handle garbage collection for Xml import.
  """
  from memops.xml import Implementation
  if gc.isenabled():
    gc.disable()
    try:
      result = Implementation.loadFromStream(stream, topObjId=topObjId, 
                                             topObject=topObject, 
                                             partialLoad=partialLoad)
    finally:
      gc.enable()
  else:
    result = Implementation.loadFromStream(stream, topObjId=topObjId, 
                                           topObject=topObject, 
                                           partialLoad=partialLoad)
  #
  return result
   
