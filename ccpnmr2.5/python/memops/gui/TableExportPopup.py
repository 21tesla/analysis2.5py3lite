
"""
======================COPYRIGHT/LICENSE START==========================

TableExportPopup.py: <write function here>

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

import Tkinter


from memops.gui.BasePopup import BasePopup
from memops.gui.Button import Button
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.Util import createDismissHelpButtonList

tabFormat = 'Tab-separated'
commaFormat = 'Comma-separated'
exportFormats = (tabFormat, commaFormat)
exportJoin = { \
  tabFormat: '\t',
  commaFormat: ',',
  }


class TableExportPopup(BasePopup):

  def __init__(self, parent, headings, entries, title = '', file = '', *args, **kw):
 
    headings = [ heading.replace('\n', ' ') for heading in headings ]

    if (headings[0].lower() == 'number'):
      self.addNumber = False
    else:
      self.addNumber = True
      headings.insert(0,  'Number')

    if (not title):
      title = 'Export data'

    self.exportSelection = parent.tableExportSelection
    self.headings = headings
    self.entries = entries
    self.file = file

    kw['title'] = title
    kw['transient'] = True
    kw['modal'] = True
    BasePopup.__init__(self, parent=parent, *args, **kw)

  def body(self, master):

    master.grid_columnconfigure(2, weight=1)

    row = 0
    label = Label(master, text='Data to export:')
    label.grid(row=row, column=0, columnspan=3, sticky=Tkinter.W)

    self.check_buttons = {}
    i = 0
    for heading in self.headings:
      row = row + 1
      isSelected = self.exportSelection.get(i, True)
      self.check_buttons[i] = c = CheckButton(master, selected=isSelected, callback=self.toggleCheckButton)
      c.grid(row=row, column=1)
      label = Label(master, text=heading)
      label.grid(row=row, column=2, sticky=Tkinter.W)
      i += 1

    row = row + 1
    button = Button(master, text='File:', command=self.findFile)
    button.grid(row=row, column=0, sticky=Tkinter.W)
    self.file_entry = Entry(master, text=self.file, width=30)
    self.file_entry.grid(row=row, column=1, columnspan=2, sticky=Tkinter.EW)

    row = row + 1
    label = Label(master, text='Format:')
    label.grid(row=row, column=0, sticky=Tkinter.W)
    self.format_menu = PulldownMenu(master, entries=exportFormats)
    self.format_menu.grid(row=row, column=1, columnspan=2, sticky=Tkinter.W)

    row = row + 1
    master.grid_rowconfigure(row, weight=1)
    texts = [ 'Save' ]
    commands = [ self.ok ]
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands,
                                          dismiss_text='Cancel')
    buttons.grid(row=row, column=0, columnspan=3, sticky=Tkinter.EW)

  def toggleCheckButton(self, isSelected):
  
    for i in self.check_buttons.keys():
      self.exportSelection[i] = self.check_buttons[i].get()

  def apply(self):

    file   = self.file_entry.get()
    format = self.format_menu.getSelected()
    join   = exportJoin[format]

    n = 0
    headings = []
    for heading in self.headings:
      if self.check_buttons[n].get():
        headings.append(n)
      n += 1

    if (not headings):
      showError('Nothing selected', 'No data selected to export.', parent=self)
      return False

    try:
      fp = open(file, 'w')
    except IOError, e:
      showError('File error', str(e))
      return False

    for n in headings:
      if (n != headings[0]):
        fp.write(join)
      fp.write(self.headings[n].encode('utf-8'))
    fp.write('\n')

    m = 0
    for data in self.entries:
      if data:
        for n in headings:
          if (n != headings[0]):
            fp.write(join)
          if (self.addNumber):
            if (n == 0):
              x = m+1
            else:
              x = data[n-1]
          else:
            x = data[n]
          if (type(x) == type(float(1))):
            if (x == 0):
              t = '0'
            elif ((abs(x) > 10000) or (abs(x) < 0.001)):
              t = '%6.5e' % x
            else:
              t = '%6.5f' % x
          else:
            t = str(x)
          fp.write(t)
        fp.write('\n')
        m = m + 1

    fp.close()

    return True

  def findFile(self):
 
    format = self.format_menu.getSelected()

    if (format == tabFormat):
      file_types = [ FileType('All', ['*']),
                     FileType('Text', ['*.txt']) ]
    else:
      file_types = [ FileType('All', ['*']),
                     FileType('CSV', ['*.csv']) ]

    popup = FileSelectPopup(self, file_types=file_types)
    
    file = popup.getFile()
    if (file):
      self.file_entry.set(file)

    popup.destroy()

if __name__ == '__main__':

  root = Tkinter.Tk()

  headings = ('Name', 'Rank', 'Serial')
  entries = ( \
    ('Abc', 1, 34345),
    ('Def', 33, 7984),
    ('Ghijklm', -9998, 4456),
  )

  popup = TableExportPopup(root, headings=headings, entries=entries,
            title='Export Test Data', file='test/test1.txt')
  popup.destroy()
