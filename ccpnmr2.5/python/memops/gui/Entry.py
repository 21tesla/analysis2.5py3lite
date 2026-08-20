
"""
======================COPYRIGHT/LICENSE START==========================

Entry.py: <write function here>

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
import re

import Tkinter

from memops.universal.Util import breakString

from memops.gui.Base import Base
from memops.gui.MessageReporter import showWarning

# this just makes it easier to set the entry value than does Tkinter.Entry

class Entry(Tkinter.Entry, Base):

  array_re = re.compile(',?\s*')
  separator = ', '
  joiner = '\n'

  def __init__(self, parent, text = '', returnCallback=None,
               isArray=False, maxCharsPerLine=0, valToStr=None,
               strToVal=None, grid=None, gridSpan=(1,1),
               sticky='w', docKey=None, tipText=None,
               leaveCallback=None, *args, **kw):

    if not docKey:
      docKey = text

    self.isArray = isArray
    self.maxCharsPerLine = maxCharsPerLine
    self.valToStr = valToStr
    self.strToVal = strToVal

    if (not kw.has_key('bg') and not kw.has_key('background')):
      kw['bg'] = 'white'

    Tkinter.Entry.__init__(self, parent, *args, **kw)
    Base.__init__(self, docKey=docKey, tipText=tipText)

    self.parent = parent

    self.set(text)

    #self.bind('<Enter>', self.enter)
    #self.bind('<Leave>', self.leave)
    #self.config(takefocus=0)
    
    if grid is not None:
      row, col = grid
      rowSpan, colSpan = gridSpan
      self.grid(row=row, column=col, rowspan=rowSpan,
                columnspan=colSpan, sticky=sticky)

    self.returnCallback = returnCallback
    if returnCallback:
      self.bind('<KeyPress-Return>', returnCallback)

    self.leaveCallback = leaveCallback
    if leaveCallback:
      self.bind('<Leave>', leaveCallback, '+')
    
  def destroy(self):

    if self.returnCallback:
      self.bind('<KeyPress-Return>')

    if self.leaveCallback:
      self.bind('<Leave>')

    Tkinter.Entry.destroy(self)

  def get(self):

    # TBD: str below is used to convert unicode to string
    # (in Tcl/Tk 8.5 it looks like everything is unicode)
    ###s = str(Tkinter.Entry.get(self).strip())
    s = Tkinter.Entry.get(self).strip().encode('utf-8')
    if (self.isArray):
      if (self.maxCharsPerLine):
        s = s.replace(self.joiner, self.separator)
      if (s):
        value = re.split(self.array_re, s)
      else:
        value = []
      if (self.strToVal):
        value = [self.strToVal(v) for v in value if v]
    elif (self.strToVal):
      if (s):
        try:
          value = self.strToVal(s)
        except:
          value = None
          showWarning('Error', 'Could not convert value to required type')
      else:
        value = None
    else:
      value = s

    return value

  def set(self, text):

    if (text is None):
      s = ''
    elif (type(text) == type('')):
      s = text
    elif (self.isArray):
      if (self.valToStr):
        s = self.separator.join([self.valToStr(x) for x in text])
      else:
        s = self.separator.join(text)
      if (self.maxCharsPerLine):
        s = breakString(s, self.separator, self.joiner, self.maxCharsPerLine)
    elif (self.valToStr):
      s = self.valToStr(text)
    else:
      s = text

    self.delete(0, Tkinter.END)
    self.insert(0, s)

  #def enter(self, event):

  #  self.focus()

  #def leave(self, event):

  #  self.parent.focus()

if (__name__ == '__main__'):

  from memops.gui.Button import Button

  msg = 'hello world'

  def func():

    global count

    msg = entry.get()
    entry.set(msg + '*')

  root = Tkinter.Tk()
 
  entry = Entry(root, text='hello world', tipText='Instructions')
  entry.grid()
  entry = Entry(root, text='hello world 2')
  entry.grid()
  button = Button(root, text='hit me', command=func)
  button.grid()

  root.mainloop()

