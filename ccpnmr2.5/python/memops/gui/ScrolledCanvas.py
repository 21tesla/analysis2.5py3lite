
"""
======================COPYRIGHT/LICENSE START==========================

ScrolledCanvas.py: <write function here>

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

from memops.gui.Canvas import Canvas
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Frame import Frame
from memops.gui.Menu import Menu
from memops.universal.Util import buttonClick, buttonRelease, buttonMotion


class ScrolledCanvas(Frame):

  def __init__(self, parent, resizeCallback = None, width=600, height=600, *args, **kw):

    self.bbox = None
    self.busy = 0
    self.initialX = None
    self.initialY = None
    self.resizeCallback = resizeCallback

    apply(Frame.__init__, (self, parent) + args, kw)
 
    self.grid_rowconfigure(0, weight=1)
    self.grid_columnconfigure(0, weight=1)

    self.menu = Menu(self, tearoff=0, include_event=True)
    self.configMenu()

    self.canvas = Canvas(self, relief='flat', borderwidth=0, width=width, height=height)
    self.canvas.configure(xscrollincrement=2, yscrollincrement=2)
    self.canvas.grid(row = 0, column = 0, sticky = Tkinter.NSEW )
    
    self.horizScrollbar = Tkinter.Scrollbar(self, bg=self.cget('bg'), command=self.canvas.xview, orient=Tkinter.HORIZONTAL, borderwidth=1)
    self.vertScrollbar  = Tkinter.Scrollbar(self, bg=self.cget('bg'), command=self.canvas.yview, orient=Tkinter.VERTICAL, borderwidth=1)
    self.canvas.configure(xscrollcommand=self.horizScrollbar.set,yscrollcommand=self.vertScrollbar.set)

    self.canvas.bind('<Configure>',        self.resizeAfter)
    self.canvas.bind('<Button-1>',         self.mouseButton1)

    # self.canvas.bind('<Button-2>',         self.mouseButton2)
    # self.canvas.bind('<Button-3>',         self.mouseButton3)
    self.canvas.bind(buttonClick(2),         self.mouseButton2)
    self.canvas.bind(buttonClick(3),         self.mouseButton3)

    self.canvas.bind('<ButtonRelease-1>',  self.mouseButtonRelease1)

    # self.canvas.bind('<ButtonRelease-2>',  self.mouseButtonRelease2)
    # self.canvas.bind('<B2-Motion>',        self.mouseScroll)
    # self.canvas.bind('<B3-Motion>',        self.doNothing)
    self.canvas.bind(buttonRelease(2),  self.mouseButtonRelease2)
    self.canvas.bind(buttonMotion(2),        self.mouseScroll)
    self.canvas.bind(buttonMotion(3),        self.doNothing)

    self.canvas.bind('<Motion>',           self.mouseEnter)
    self.canvas.bind('<Enter>',            self.mouseEnter)

  def doNothing(self, event):
  
    pass

  def mouseEnter(self, event):

    if self.menu.winfo_ismapped():
      self.removeMenu()

  def refresh(self):
  
    self.bbox = self.canvas.bbox('all')
    bbox = self.bbox
    
    if not bbox:
      self.busy = 0
      return
    
    cWidth  = int(self.canvas.cget('width'))
    cHeight = int(self.canvas.cget('height'))
    
    if cHeight > bbox[3]-bbox[1]:
      self.vertScrollbar.grid_forget()
    else:
      self.vertScrollbar.grid(row = 0, column = 1, sticky = Tkinter.NS )

    if cWidth  > bbox[2]-bbox[0]:
      self.horizScrollbar.grid_forget()
    else:
      self.horizScrollbar.grid(row = 1, column = 0, sticky = Tkinter.EW )
    
    x1 = bbox[0]
    y1 = bbox[1]
    x2 = max(bbox[0]+cWidth, bbox[2])
    y2 = max(bbox[1]+cHeight,bbox[3])
        
    self.canvas.configure(scrollregion = (x1,y1,x2,y2)  )
    self.update_idletasks()
    self.busy = 0

  def printCanvas(self, *event):
  
    fileTypes = [  FileType('PostScript', ['*.ps']), FileType('All', ['*'])]
    fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
               title = 'Print canvas to file', dismiss_text = 'Cancel',
               selected_file_must_exist = False)

    fileName = fileSelectPopup.getFile()
    
    self.bbox = bbox = self.canvas.bbox('all')
    w = bbox[2] - bbox[0]
    h = bbox[3] - bbox[1]
    self.canvas.postscript(colormode='color',file=fileName,
                           x=bbox[0], y=bbox[1], width=w+2,
                           pagewidth='21.c', height=h+2, fontmap='fontmap')

  def configMenu(self):

    items = [{ 'kind': 'command', 'label': 'Print to file', 'command' : self.printCanvas },]
      
    self.menu.setMenuItems(items)

  def removeMenu(self, *event):

    self.menu.unpost()

  def mouseButton3(self, event):

    self.menu.popupMenu(event)

  def mouseButton1(self, event):
  
    self.removeMenu()
  
    if not self.initialX:
      self.initialX = event.x
      self.initialY = event.y

  def mouseButton2(self, event):
  
    self.removeMenu()
  
    if not self.initialX:
      self.initialX = event.x
      self.initialY = event.y
      self.initialPX = self.horizScrollbar.get()[0]
      self.initialPY = self.vertScrollbar.get()[0]

  def mouseButtonRelease1(self, event):
  
    self.initialX = None
    self.initialY = None

  def mouseButtonRelease2(self, event):
  
    self.initialX = None
    self.initialY = None
    self.initialPX = None
    self.initialPY = None

  def mouseScroll(self, event):
    
    self.menu.unpost()
    bbox = self.bbox
    if not bbox:
      return
    bW = float(bbox[2] - bbox[0])
    bH = float(bbox[3] - bbox[1])
    cWidth  = int(self.canvas.cget('width'))
    cHeight = int(self.canvas.cget('height'))

    if cWidth < bW:
      dx = self.initialX - event.x
      prop = self.initialPX + ( dx/bW )
      self.canvas.xview('moveto', prop)

    if cHeight < bH:
      dy = self.initialY - event.y
      prop = self.initialPY + ( dy/bH )
      self.canvas.yview('moveto', prop)

  def resizeAfter(self, event):
  
    if self.busy:
      return    
    else:
      self.busy=1
      self.removeMenu()
      self.after_idle( lambda:self.resize(event) )

  def resize(self, event):

    if self.resizeCallback:
      self.after_idle( lambda: self.resizeCallback(event.width, event.height) )
    self.canvas.configure(width=event.width,height=event.height)
    self.after_idle( self.refresh )

