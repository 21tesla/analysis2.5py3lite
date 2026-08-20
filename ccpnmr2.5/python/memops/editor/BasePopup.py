
"""
======================COPYRIGHT/LICENSE START==========================

BasePopup.py: <write function here>

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
import os

from memops.general.Constants import True, False
from memops.general import Application
from memops.general import Implementation

from memops.universal.Io import splitPath

import memops.gui.BasePopup

class BasePopup(memops.gui.BasePopup.BasePopup):

  def __init__(self, parent, project = None, popup_name = None,
               turnNotifiersOnOff=True, application = None, *args, **kw):

    if not hasattr(self, 'project'):
      self.project = project
    if not hasattr(self, 'application'):
      self.application = application
    self.turnNotifiersOnOff = turnNotifiersOnOff
    self.notifies = {}

    if project and popup_name:
      if not kw.has_key('location'):
        location = getPopupLocation(project, popup_name)
        if (location):
          kw['location'] = '+%s+%s' % location

    if (not kw.has_key('class_') and project and hasattr(project, 'application') and project.application):
      kw['class_'] = project.application.name

    memops.gui.BasePopup.BasePopup.__init__(self, parent, *args, **kw)

    self.initPopupLocation(popup_name, doSet=False)

  def initPopupLocation(self, popup_name, doSet = True):

    self.popup_name = popup_name

    if self.project and popup_name:
      self.bind('<Configure>', self.changedPopupLocation)
      if doSet:
        location = getPopupLocation(self.project, popup_name)
        if location:
          self.geometry('+%s+%s' % location)

  def changedPopupLocation(self, *event):

    location = self.geometry()

    # TBD: below supposedly does not give location on all platforms
    # (when called because of binding to 'Configure')
    n = location.find('+')
    xy = location[n+1:]
    (x, y) = map(int, xy.split('+'))

    setPopupLocation(self.project, self.popup_name, (x, y))

  def setTitle(self, title = ''):
 
    if not title:
      title = self.__class__.__name__
 
    if hasattr(self, 'project') and self.project:
      title = title + ' (%s)' % self.project.name
    elif hasattr(self, 'parent') and hasattr(self.parent, 'project') and self.parent.project:
      title = title + ' (%s)' % self.parent.project.name
 
    memops.gui.BasePopup.BasePopup.setTitle(self, title)

  def getTitle(self):

    title = self.title()
    n = title.rfind(' (')
    if n >= 0:
      return title[:n]
    else:
      return title

  def turnOffNotifiers(self, modifyList = False):
 
    for classname in self.notifies.keys():
      cc = self.notifies[classname]
      for key in cc.keys():
        try:
          (funcname, application, keyword) = key
        except:
          funcname = key
          application = keyword = None
        ff = cc[key]
        for notify in ff:
          self.unregisterNotify(notify, classname, funcname, application, keyword, modifyList=modifyList)
 
  def turnOnNotifiers(self):
 
    for classname in self.notifies.keys():
      cc = self.notifies[classname]
      notifies = Implementation.getClassFromFullName(classname)._notifies
      for key in cc.keys():
        try:
          (funcname, application, keyword) = key
        except:
          funcname = key
          application = keyword = None
        ff = cc[key]
        ll = notifies.get(key)
        for notify in ff:
          # do checks below in case someone does open twice in row
          if ll is None or notify not in ll:
            self.registerNotify(notify, classname, funcname, application, keyword, modifyList=False)

  def open(self):
 
    if self.turnNotifiersOnOff:
      self.turnOnNotifiers()
    memops.gui.BasePopup.BasePopup.open(self)
 
  def close(self):
 
    memops.gui.BasePopup.BasePopup.close(self)
    if self.turnNotifiersOnOff:
      self.turnOffNotifiers()
 
  def registerNotify(self, notify, classname, funcname = '', application = None, keyword = None, modifyList = True):
 
    if application is None:
      application = self.application

    if self.application:
      self.application.registerNotify(notify, classname, funcname, keyword)
    elif application:
      Application.registerNotify(notify, classname, funcname, application, keyword)
    else:
      Implementation.registerNotify(notify, classname, funcname)
 
    if modifyList:
      cc = self.notifies.get(classname)
      if cc is None:
        cc = self.notifies[classname] = {}
      if application is None:
        key = funcname
      else:
        key = (funcname, application, keyword)
      ff = cc.get(key)
      if ff is None:
        ff = cc[key] = []
      ff.append(notify)
 
  def unregisterNotify(self, notify, classname, funcname = '', application = None, keyword = None, modifyList = True):
 
    if application is None:
      application = self.application

    if modifyList:
      if application is None:
        key = funcname
      else:
        key = (funcname, application, keyword)
      try:
        cc = self.notifies[classname]
        ff = cc[key]
        ff.remove(notify)
        if not ff:
          del cc[key]
        if not cc:
          del self.notifies[classname]
      except: # already been unregistered
        pass
 
    if self.application:
      self.application.unregisterNotify(notify, classname, funcname, keyword)
    elif application:
      Application.unregisterNotify(notify, classname, funcname, application, keyword)
    else:
      Implementation.unregisterNotify(notify, classname, funcname)
 
  def destroy(self):
 
    self.turnOffNotifiers(modifyList=True)
 
    memops.gui.BasePopup.BasePopup.destroy(self)

def getPopupLocation(project, popup_name):
 
  if (project and popup_name and hasattr(project, 'application')):
    application = project.application
    if not isinstance(application, Application.Application):
      return None
 
    ### TBD v2: looks like ought to pass in topObject instead of project
    ### to hold this appData
    ###x = application.getValue(project, keyword='popupX_' + popup_name)
    x = None
    if (x is None):
      return None
 
    ###y = application.getValue(project, keyword='popupY_' + popup_name)
    y = None
    if (y is None):
      return None
 
    return (x, y)

  else:
    return None
 
def setPopupLocation(project, popup_name, location):
 
  if (project and popup_name and location and hasattr(project, 'application')):
    application = project.application
    if not isinstance(application, Application.Application):
      return

    (x, y) = location
    ### TBD v2
    """
    application.setValue(project, keyword='popupX_' + popup_name, value=x)
    application.setValue(project, keyword='popupY_' + popup_name, value=y)
"""

