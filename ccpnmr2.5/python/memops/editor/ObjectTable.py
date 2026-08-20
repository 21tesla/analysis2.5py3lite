 
"""
======================COPYRIGHT/LICENSE START==========================

ObjectTable.py: <write function here>

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

from memops.gui.ScrolledMatrix import ScrolledMatrix

from memops.editor.Util import getKeyList, getObjectKey

class ObjectTable(ScrolledMatrix):

  def __init__(self, parent, metaclass, objects = None, includeNumber = True,
               *args, **kw):

    self.metaclass = metaclass
    self.includeNumber = includeNumber

    headings = []
    if (includeNumber):
      headings.append('Number')
 
    self.keyList = keyList = getKeyList(metaclass)
    for key in keyList:
      attr = key[-1]
      s = attr.name
      if (len(key) > 1):
        if (type(key[-2]) == type(0)):
          role = key[-3]
          n = key[-2]
          t = '%s[%d]' % (role.otherClass.name, n)
        else:
          role = key[-2]
          t = role.otherClass.name
      else:
        t = metaclass.name
      headings.append('%s\n%s' % (t, s))

    ScrolledMatrix.__init__(self, parent, headingList=headings, *args, **kw)

    self.setObjects(objects)

    # do not need below for now since only have keys in ObjectTable (which are frozen)
    #Implementation.registerNotify(self.updatedObject, metaclass.qualifiedName(), '')

  """
  def destroy(self):

    #Implementation.unregisterNotify(self.updatedObject, metaclass.qualifiedName(), '')

    ScrolledMatrix.destroy(self)

  def updatedObject(self, object):

    if (object in self.objectList):
      # need to pass in objectList because otherwise currentObject set to None
      self.update(objectList=self.objectList)
  """

  def setObjects(self, objects = None):

    if (objects is None):
      objects = []

    n = 0
    textMatrix = []
    for object in objects:

      n = n + 1
      text = []
      if (self.includeNumber):
        text.append(n)

      value = getObjectKey(object, self.keyList)
      for v in value:
        text.append(v)

      textMatrix.append(text)

    try:
      self.update(objectList=objects, textMatrix=textMatrix)
    except:
      print '***', self.metaclass.name, textMatrix

if (__name__ == '__main__'):

  import Tkinter

  from memops.api import Implementation as Impl

  p = Impl.Project(name='test_project')

  root = Tkinter.Tk()
  root.grid_rowconfigure(0, weight=1)
  root.grid_columnconfigure(0, weight=1)

  table = ObjectTable(root, metaclass=p.metaclass, objects=[p])
  table.grid(sticky=Tkinter.NSEW)

  root.mainloop()
