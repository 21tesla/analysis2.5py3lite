
"""
======================COPYRIGHT/LICENSE START==========================

ObjectsTable.py: <write function here>

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
# ObjectsTable is for displaying a selection of attributes
# of an array of objects, allowing the user to update their values
# It requires a list of attributes, how these are labelled on the
# screen (attributeHeadings, keyed on attribute), and the widget
# classes for each attribute (widgetClasses, keyed on attribute).
# The widget classes must allow creation with only parent widget
# as an argument, and also must have set and get functions.

import Tkinter


from memops.gui.ButtonList import ButtonList
from memops.gui.Frame import Frame
from memops.gui.Label import Label

class ObjectsTable(Frame):

  def __init__(self, parent, objects = None,
               attributes = None, attributeHeadings = None,
               attributeLabel = 'Attribute', widgetClasses = None,
               constructorArg = None, constructorKw = None,
               tableHeading = 'Dimension', tableHeadingOffset = 1,
               callback = None, attributeDirection = Tkinter.VERTICAL,
               includeResetUpdate = True, help_msg = '', help_url = '', *args, **kw):

    if (objects is None):
      objects = []

    if (attributes is None):
      attributes = []

    if (attributeHeadings is None):
      attributeHeadings = {}

    if (widgetClasses is None):
      widgetClasses = {}

    if (constructorArg is None):
      constructorArg = {}

    if (constructorKw is None):
      constructorKw = {}

    apply(Frame.__init__, (self, parent) + args, kw)

    self.ndim = len(objects)
    self.objects = objects
    self.attributes = attributes
    self.callback = callback

    topFrame = Frame(self)
    topFrame.grid(row=0, column=0, sticky=Tkinter.NSEW)
    bottomFrame = Frame(self)
    bottomFrame.grid(row=1, column=0, sticky=Tkinter.NSEW)

    if (attributeLabel):
      label = Label(topFrame, text=attributeLabel, relief=Tkinter.GROOVE)
      label.grid(row=0, column=0, sticky=Tkinter.EW)

    n = 0
    for i in range(self.ndim):
      text = tableHeading + ' ' + str(i + tableHeadingOffset)
      label = Label(topFrame, text=text, relief=Tkinter.GROOVE)
      if (attributeDirection == Tkinter.VERTICAL):
        r = n
        c = i + 1
      else:
        r = i + 1
        c = n
      label.grid(row=r, column=c, sticky=Tkinter.EW)

    self.widgets = {}

    for attr in self.attributes:
      n = n + 1
      text = attributeHeadings.get(attr, attr)
      label = Label(topFrame, text=text, relief=Tkinter.GROOVE)
      if (attributeDirection == Tkinter.VERTICAL):
        r = n
        c = 0
      else:
        r = 0
        c = n
      label.grid(row=r, column=c, sticky=Tkinter.EW)

      self.widgets[attr] = self.ndim * [0]
      for i in range(self.ndim):
        arg = constructorArg.get(attr, ())
        kw = constructorKw.get(attr, {})
        self.widgets[attr][i] = widgetClasses[attr](topFrame, *arg, **kw)
        if (attributeDirection == Tkinter.VERTICAL):
          r = n
          c = i + 1
        else:
          r = i + 1
          c = n
        self.widgets[attr][i].grid(row=r, column=c, sticky=Tkinter.EW)

    texts = []
    commands = []

    if (includeResetUpdate):
      texts.extend(['Reset', 'Update'])
      commands.extend([self.reset, self.update])

    if (help_url or help_msg):
      import memops.gui.HelpPopup
      if (help_url):
        help_cmd = lambda top=self.parent.top, url=help_url: memops.gui.HelpPopup.showHelpUrl(top, url)
      else:
        help_cmd = lambda top=self.parent.top, message=help_msg: memops.gui.HelpPopup.showHelpText(top, message)
      texts.append('Help')
      commands.append(help_cmd)

    if (texts):
      buttons = ButtonList(bottomFrame, texts=texts, commands=commands)
      buttons.grid(row=0, column=0, sticky=Tkinter.EW)

    self.reset()

  def reset(self):

    for attr in self.attributes:
      for i in range(self.ndim):
        value = getattr(self.objects[i], attr)
        self.widgets[attr][i].set(value)

  def update(self):

    for attr in self.attributes:
      for i in range(self.ndim):
        value = self.widgets[attr][i].get()
        setattr(self.objects[i], attr, value)

    if (self.callback):
      self.callback()
      self.reset() # in case objects changed in callback

if __name__ == '__main__':

  from memops.gui.Entry import Entry
  from memops.gui.CheckButton import CheckButton
  from memops.gui.PulldownMenu import PulldownMenu
  from memops.gui.Util import createQuitButton

  class T:

    def __init__(self, x, y, w, z):

      self.x = x
      self.y = y
      self.w = w
      self.z = z

  root = Tkinter.Tk()

  objects = [ T(8000, 600, True, 'Bruker'), T(8000, 600, True, 'Bruker'),
              T(2000, 100, False, 'Varian') ]
  attributes = [ 'x', 'y', 'w', 'z' ]
  attributeHeadings = { 'x': 'SW', 'y': 'SF', 'w': 'TPPI?', 'z': 'Manufacturer' }
  widgetClasses= { 'x': Entry, 'y': Entry, 'w': CheckButton, 'z': PulldownMenu }
  constructorKw = { 'z': { 'entries': [ 'Bruker', 'Varian' ] } }
  root.top = root
  t = ObjectsTable(root, objects=objects,
                 widgetClasses=widgetClasses, constructorKw = constructorKw,
                 attributes=attributes, attributeHeadings=attributeHeadings,
                 attributeDirection = Tkinter.HORIZONTAL, help_msg='Help message')
  t.grid(row=0, column=0, sticky=Tkinter.NSEW)

  button = createQuitButton(root)
  button.grid(row=1)

  root.mainloop()
