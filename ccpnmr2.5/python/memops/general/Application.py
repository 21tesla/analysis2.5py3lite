"""
======================COPYRIGHT/LICENSE START==========================

Application.py: code for CCPN data model and code generation framework

Copyright (C) 2008  (CCPN Project)

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
import types

from memops.api import Implementation as ApiImp
from memops.general import Implementation as GenImp

classMapping = {}
classMapping[types.BooleanType] = ApiImp.AppDataBoolean
classMapping[types.IntType] = ApiImp.AppDataInt
classMapping[types.FloatType] = ApiImp.AppDataFloat
classMapping[types.StringType] = ApiImp.AppDataString

allowedNotifyFuncs = set(('addAppData', 'setAppDataValue', 'setAppDataValues'))

def registerNotify(notify, classname, funcname, application = None, keyword = None):

  if application is None or keyword is None:
    assert application is None and keyword is None, 'application = %s, keyword = %s, both must be None' % (application, keyword)
    return GenImp.registerNotify(notify, classname, funcname)

  if funcname not in allowedNotifyFuncs:
    raise GenImp.ApiError('illegal funcname "%s", must be in %s' % (funcname, allowedNotifyFuncs))

  notifies = GenImp.getClassFromFullName(classname)._notifies
  notifies = notifies.setdefault((funcname, application, keyword), [])
  notifies.append(notify)

def unregisterNotify(notify, classname, funcname, application = None, keyword = None):

  if application is None or keyword is None:
    assert application is None and keyword is None, 'application = %s, keyword = %s, both must be None' % (application, keyword)
    return GenImp.unregisterNotify(notify, classname, funcname)

  if funcname not in allowedNotifyFuncs:
    raise GenImp.ApiError('illegal funcname "%s", must be in %s' % (funcname, allowedNotifyFuncs))

  try:
    notifies = GenImp.getClassFromFullName(classname)._notifies
    notifies = notifies[(funcname, application, keyword)]
    notifies.remove(notify)
  except:
    pass

def createAppData(object, application, keyword, value):

  clazz = classMapping.get(type(value))

  appData = clazz(application=application, keyword=keyword, value=value)
  object.addApplicationData(appData)

  return appData

def removeAllAppData(object, application, keyword):

  appDataList = object.findAllApplicationData(application=application, keyword=keyword)
  for appData in appDataList:
    object.removeApplicationData(appData)

  return appDataList

def doAppDataNotifies(object, funcname, application, keyword, value):

  notifies = object.__class__._notifies
  notifies = notifies.get((funcname, application, keyword))

  if notifies:
    for notify in notifies:
      notify(object)

def addAppData(object, application, keyword, value):

  createAppData(object, application, keyword, value)
  doAppDataNotifies(object, 'addAppData', application, keyword, value)

def setAppDataValue(object, application, keyword, value = None):

  removeAllAppData(object, application, keyword)

  if value is None:
    appData = None
  else:
    appData = createAppData(object, application, keyword, value)

  doAppDataNotifies(object, 'setAppDataValue', application, keyword, value)

  return appData

def setAppDataValues(object, application, keyword, values = None):

  if not values:
    values = []

  removeAllAppData(object, application, keyword)

  appDataList = []
  for value  in values:
    appData = createAppData(object, application, keyword, value)
    appDataList.append(appData)

  doAppDataNotifies(object, 'setAppDataValues', application, keyword, values)

  return appDataList

class Application:

  def __init__(self, name, alternativeNames = None):

    if not alternativeNames:
      alternativeNames = []

    self.name = name
    self.alternativeNames = alternativeNames

  # gets (first) application data value for specified keyword for object
  # if no application data exists and defaultValue is not None then creates

  def getValue(self, object, keyword, defaultValue = None, deleteAppData = False):

    #print 'getValue1', object, keyword, defaultValue, type(defaultValue)
    appData = object.findFirstApplicationData(application=self.name,
                                              keyword=keyword)
    #print 'getValue2', appData

    if not appData:
      for name in self.alternativeNames:
        appData = object.findFirstApplicationData(application=name, keyword=keyword)
        if appData:
          value = appData.value
          removeAllAppData(object, name, keyword)
          appData = createAppData(object, self.name, keyword, value)
          break

    if not appData and defaultValue is not None:
      appData = createAppData(object, self.name, keyword, defaultValue)

    if appData:
      value = appData.value
      if deleteAppData:
        removeAllAppData(object, self.name, keyword)
    else:
      value = None

    return value

  # gets (all) application data value for specified keyword for object
  # if no application data exists and defaultValues is not None then creates

  def getValues(self, object, keyword, defaultValues = None):

    appDataList = object.findAllApplicationData(application=self.name, keyword=keyword)

    if not appDataList:
      appDataList = []
      if defaultValues:
        for value in defaultValues:
          appData = createAppData(object, self.name, keyword, value)
          appDataList.append(appData)

    return [ appData.value for appData in appDataList ]

  # sets (overwrites) application data value for specified keyword for object

  def setValue(self, object, keyword, value = None):
 
    setAppDataValue(object, self.name, keyword, value)

  # adds application data value for specified keyword for object

  def addValue(self, object, keyword, value):
 
    addAppData(object, self.name, keyword, value)

  # sets (overwrites) application data values for specified keyword for object

  def setValues(self, object, keyword, values = None):
 
    setAppDataValues(object, self.name, keyword, values)

  def registerNotify(self, notify, classname, funcname, keyword = None):

    if keyword is None:
      application = None
    else:
      application = self.name

    registerNotify(notify, classname, funcname, application, keyword)

  def unregisterNotify(self, notify, classname, funcname, keyword = None):

    if keyword is None:
      application = None
    else:
      application = self.name

    unregisterNotify(notify, classname, funcname, application, keyword)
