"""
======================COPYRIGHT/LICENSE START==========================

Util.py: 

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

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
- PDBe website (http://www.ebi.ac.uk/pdbe/)

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

from memops.universal.Util import upperFirst

from memops.gui import Base as guiBase
from memops.gui import Util as guiUtil

from memops.editor.WebBrowser import ProjectWebBrowser

def getConcreteClasses(clazz):
 
  classes = []
 
  if (not clazz.isAbstract):
    classes.append(clazz)
 
  for subclass in clazz.subclasses:
    classes.extend(getConcreteClasses(subclass))
 
  return classes

def getKeyList(clazz, includeRoot = False):

  keyList = []

  keys = clazz.keyNames
  for key in keys:
    if (key in clazz.roleNames):
      role = clazz.elementDict[key]
      assert role.locard == role.hicard, 'role = %s, locard = %s, hicard = %s' % (role.name, role.locard, role.hicard)
      partKeyList = getKeyList(role.otherClass, includeRoot=includeRoot)
      for n in range(role.hicard):
        if (role.hicard == 1):
          start = (role,)
        else:
          start = (role, n)
        for key in partKeyList:
          keyList.append(start + key)
    else:
      attr = clazz.elementDict[key]
      keyList.append((attr,))

  if (clazz.parentRole is not None):
    role = clazz.parentRole
    start = (role,)
    clazz = role.otherClass
    if (includeRoot or clazz.parentRole is not None):
      partKeyList = getKeyList(clazz, includeRoot=includeRoot)
      for key in partKeyList:
        keyList.append(start + key)

  return keyList

def printKeyList(keyList):

  for key in keyList:
    s = []
    i = 0
    while (i < len(key)-1):
      role = key[i]
      if (role.hicard == 1):
        t = role.name
      else:
        n = key[i+1]
        i = i + 1
        t = '%s[%d]' % (role.baseName, n)
      s.append(t)
      i = i + 1
    attr= key[-1]
    t = attr.name
    s.append(t)
    print ', '.join(s)

def getObjectKey(object, keyList = None):

  if (keyList is None):
    keyList = getKeyList(object.metaclass)

  value = []
  for key in keyList:
    obj = object
    i = 0
    while (i < len(key)-1):
      role = key[i]
      obj = getattr(obj, role.name)
      if (role.hicard != 1):
        n = key[i+1]
        obj = obj[n]
        i = i + 1
      i = i + 1
    attr = key[-1]
    v = getattr(obj, attr.name)
    value.append(v)

  return value

def getAllObjects(top, metaclass):

  topClass = top.metaclass
  if (metaclass is topClass):
    return [ top ]

  clazz = metaclass
  while (clazz.parentRole.otherClass is not topClass):
    clazz = clazz.parentRole.otherClass

  func = getattr(top, 'findAll%s' % upperFirst(clazz.parentRole.otherRole.name))

  objs = func()
  allObjs = []
  for obj in objs:
    allObjs.extend(getAllObjects(obj, metaclass))

  return allObjs

def createDismissHelpButtonList(parent, *args, **kw):

  popup = guiBase.getPopup(parent)
  if (not kw.has_key('webBrowser')):
    webBrowser = ProjectWebBrowser(popup.top, popup=popup, project=popup.project)
    kw['webBrowser'] = webBrowser

  return guiUtil.createDismissHelpButtonList(parent, *args, **kw)

def createHelpButtonList(parent, *args, **kw):

  popup = guiBase.getPopup(parent)
  if (not kw.has_key('webBrowser')):
    webBrowser = ProjectWebBrowser(popup.top, popup=popup, project=popup.project)
    kw['webBrowser'] = webBrowser

  return guiUtil.createHelpButtonList(parent, *args, **kw)
