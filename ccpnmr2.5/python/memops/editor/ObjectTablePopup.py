 
"""
======================COPYRIGHT/LICENSE START==========================

ObjectTablePopup.py: <write function here>

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
import Tkinter
 
 
from memops.general.Constants import infinity, apiSubDirs
from memops.general import Implementation
 
from memops.gui.ButtonList import ButtonList
from memops.gui.MessageReporter import showError

from memops.editor.BasePopup import BasePopup
from memops.editor.ObjectTable import ObjectTable
from memops.editor.Util import getAllObjects

class ObjectTablePopup(BasePopup):

  def __init__(self, parent, root, metaclass, onlyShow = False, *args, **kw):

    self.root = root
    self.metaclass = metaclass
    self.onlyShow = onlyShow

    BasePopup.__init__(self, parent, title='Browser for %s objects' % metaclass.name,
                       location='+50+50', *args, **kw)
 
  def body(self, master):
 
    master.grid_rowconfigure(0, weight=1)
    master.grid_columnconfigure(0, weight=1)

    self.table = ObjectTable(master, self.metaclass)
    self.table.grid(row=0, column=0, sticky=Tkinter.NSEW)

    if (self.onlyShow):
      texts = [ 'Close' ]
      commands = [ self.close ]
    else:
      texts = [ 'Ok', 'Cancel' ]
      commands = [ self.ok, self.close ]

    self.buttons = ButtonList(master, texts=texts, commands=commands,
                              direction=Tkinter.HORIZONTAL, expands=True)
    self.buttons.grid(row=1, column=0, sticky=Tkinter.EW)

    self.selected = None

    self.setObjects()

    self.doRegisters()

  def doRegisters(self):

    # do not need '' since dealing with keys
    #for func in ('', '__init__', 'delete'):
    for func in ('__init__', 'delete'):
      self.registerNotify(self.setObjects, self.metaclass.qualifiedName(), func)

  # unregisters dealt with by BasePopup

  def apply(self):

    selected = self.table.currentObject

    if (not selected):
      showError('No selection', "No object selected", parent=self)
      return False

    self.selected = selected

    return True

  def setObjects(self, *extra):

    objects = getAllObjects(self.root, self.metaclass)

    self.table.setObjects(objects)
