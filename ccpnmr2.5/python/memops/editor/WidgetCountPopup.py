
"""
======================COPYRIGHT/LICENSE START==========================

WidgetCountPopup.py: <write function here>

Copyright (C) 2009 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

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

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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

from memops.gui.Base import getRoot, getWidgetCount
from memops.gui.Button import Button
from memops.gui.IntEntry import IntEntry
from memops.gui.Label import Label
from memops.gui.RadioButtons import RadioButtons
from memops.gui.Util import createDismissHelpButtonList

from memops.editor.BasePopup import BasePopup

class WidgetCountPopup(BasePopup):

  on_off_entries = ['on', 'off']

  def __init__(self, parent, help_msg = '', help_url = '', *args, **kw):

    self.help_msg = help_msg
    self.help_url = help_url

    BasePopup.__init__(self, parent=parent, title='Widget Count', *args, **kw)

  def body(self, master):

    self.alarm_id = None
    self.root_widget = getRoot(self)

    row = 0
    label = Label(master, text='Widget count:')
    label.grid(row=row, column=0, sticky=Tkinter.W)
    self.count_label = Label(master, text='', tipText='Number of Tkinter widget objects')
    self.count_label.grid(row=row, column=1, sticky=Tkinter.W)

    row = row + 1
    ind = 1
    label = Label(master, text='Auto count:')
    label.grid(row=row, column=0, sticky=Tkinter.W)
    self.on_off_buttons = RadioButtons(master, entries=self.on_off_entries,
                                       select_callback=self.applyAuto, selected_index=ind)
    self.on_off_buttons.grid(row=row, column=1, sticky=Tkinter.EW)
 
    row = row + 1
    label = Label(master, text='Auto frequency:')
    label.grid(row=row, column=0, sticky=Tkinter.W)
    self.freq_entry = IntEntry(master, text=60,
                               returnCallback=self.applyAuto)
    self.freq_entry.grid(row=row, column=1, sticky=Tkinter.EW)
    label = Label(master, text='seconds')
    label.grid(row=row, column=2, sticky=Tkinter.W)
 
    row = row + 1
    texts = [ 'Do Immediate Count' ]
    commands = [ self.applyManual ]
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands,
                                          help_msg=self.help_msg, help_url=self.help_url)
    buttons.grid(row=row, column=0, columnspan=3, sticky=Tkinter.EW)

    self.doCount()

  def applyAuto(self, *event):

    on_off = self.on_off_buttons.get()
    if on_off == 'on':
      self.setCountOn()
    else:
      self.setCountOff()

  def applyManual(self):

    self.doCount()

  def setCountOn(self):

    self.setCountOff()

    freq = self.freq_entry.get()
    if freq is not None and freq > 0:
      freq *= 1000  # convert from sec to msec
      self.alarm_id = self.after(freq, self.doCountAuto)

  def setCountOff(self):

    if self.alarm_id:
      self.after_cancel(self.alarm_id)
      self.alarm_id = None

  def doCountAuto(self):

    if self.on_off_buttons.get() == 'on':
      self.doCount()
      freq = self.freq_entry.get()
      if freq is not None and freq > 0:
        freq *= 1000  # convert from sec to msec
        self.alarm_id = self.after(freq, self.doCountAuto)

  def doCount(self):

    count = getWidgetCount(self.root_widget)
    self.count_label.set('%d' % count)

if __name__ == '__main__':

  root = Tkinter.Tk()
  root.withdraw()
  popup = WidgetCountPopup(root)
  root.mainloop()
