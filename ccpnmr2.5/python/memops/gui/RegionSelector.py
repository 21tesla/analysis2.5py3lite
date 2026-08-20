
"""
======================COPYRIGHT/LICENSE START==========================

RegionSelector.py: <write function here>

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

from memops.universal.Util import formatDecimals

from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.RegionScrollbar import RegionScrollbar

class RegionSelector(Frame):

  def __init__(self, parent, label = '', world_region = None, view_region = None,
               orient = Tkinter.HORIZONTAL, allow_resize = True,
               width = 20, callback = None, borderwidth = 1,
               show_text = True, text_color = '#000000', text_decimals = 2,
               units_scroll = 0.1, pages_scroll = 1.0,
               menu_entries = None, menu_callback = None,
               min_thickness = None,
               *args, **kw):
 
    self.menu_entries = menu_entries
    self.myCallback = callback

    Frame.__init__(self, parent, *args, **kw)
 
    self.text_decimals = text_decimals

    self.label = Label(self, text=label, width=4)

    self.menu = PulldownMenu(self, callback=menu_callback, entries=menu_entries)

    self.entry = FloatEntry(self, width=6, returnCallback=self.adjustScrollbar)

    self.region_scrollbar = RegionScrollbar(self, world_region=world_region,
                                view_region=view_region,
                                orient=orient, allow_resize=allow_resize,
                                width=width, callback=self.doCallback,
                                borderwidth=borderwidth, show_text=show_text,
                                text_color=text_color, text_decimals=text_decimals,
                                units_scroll=units_scroll, pages_scroll=pages_scroll,
                                min_thickness=min_thickness)

    self.gridAll()

  def gridAll(self):

    col = 0
    if (self.menu_entries and len(self.menu_entries) > 1):
      self.menu.grid(row=0, column=col, sticky=Tkinter.EW)
      col = col + 1
    else:
      self.menu.grid_forget()

    self.label.grid(row=0, column=col, sticky=Tkinter.EW)
    col = col + 1
    self.entry.grid(row=0, column=col, sticky=Tkinter.EW)
    self.grid_columnconfigure(col, weight=0)
    col = col + 1
    self.region_scrollbar.grid(row=0, column=col, sticky=Tkinter.NSEW)
    self.grid_columnconfigure(col, weight=1)
    col = col + 1
    self.grid_columnconfigure(col, weight=0)

  def setMinThickness(self, min_thickness):

    self.region_scrollbar.setMinThickness(min_thickness)

  def setMenuEntries(self, menu_entries):

    self.menu_entries = menu_entries
    self.menu.replace(menu_entries)
    self.gridAll()

  def getMenuEntry(self):

    return self.menu.getSelected()

  def adjustScrollbar(self, *event):

    try:
      x = float(self.entry.get())
    except:
      showError('Entry error', 'Need to enter float in scrollbar box')
      self.setEntry()
      return
 
    (v0, v1) = self.region_scrollbar.view_region
    d = 0.5 * (v1 - v0)
 
    self.region_scrollbar.setViewRegion(x-d, x+d, do_callback=True)

  def doCallback(self, view_region):

    if (self.myCallback):
      self.myCallback(view_region)
 
    #print 'doCallback', view_region
    self.setEntry(view_region)

  def setEntry(self, view_region = None):

    if (not view_region):
      view_region = self.region_scrollbar.view_region
    (v0, v1) = view_region
    x = 0.5 * (v0 + v1)
    s = formatDecimals(x, decimals=self.text_decimals)
    self.entry.set(s)

  def __getattr__(self, name):
 
    # dispatch everything not defined by RegionSelector to scrollbar widget
 
    try:
      return getattr(self.__dict__['region_scrollbar'], name)
    except:
      raise AttributeError, "RegionSelector instance has no attribute '%s'" % name

if __name__ == '__main__':
 
  from memops.universal.Region1D import Region1D

  root = Tkinter.Tk()
 
  wr = Region1D(20.0, -4.0)
  vr = Region1D(7.2, 6.5)
 
  r = RegionSelector(root, label='1H',
                     world_region=wr, view_region=vr)
  r.pack(fill=Tkinter.BOTH, expand=Tkinter.YES)
 
  root.mainloop()
