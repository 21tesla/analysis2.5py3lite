
"""
======================COPYRIGHT/LICENSE START==========================

ProgressBar.py: <write function here>

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
import time

import Tkinter

from memops.general.Constants import True, False

from memops.gui.Canvas import Canvas
from memops.gui.Label import Label
from memops.gui.BasePopup import BasePopup
 
class ProgressBar(BasePopup):

  def __init__(self, parent, text = '', progress = 0, total = 100, title=None,
               width = 200, height = 15, transient = True, *args, **kw):
    
    self.progress = float(progress)
    self.cWidth   = width
    self.cHeight  = height
    self.total    = float(total) or 1.0
    self.text     = text
    
    BasePopup.__init__(self, parent=parent, title=title or 'Progress Bar',
                       transient=transient, *args, **kw)
    
  def body(self, parent):
  
    self.bw = 2 # borderwidth
    width  = self.cWidth  + self.bw
    height = self.cHeight + self.bw
    self.label  = Label(parent,text=self.text)
    self.label.grid(row=0,column=0, sticky=Tkinter.W)
    self.percent = Label(parent,text='    %')
    self.percent.grid(row=0,column=2, sticky=Tkinter.W)
    self.canvas = Canvas(parent, background='grey70', width=width, height=height)
    self.canvas.grid(row=0,column=1, sticky=Tkinter.W)
    self.canvas.create_rectangle(0,0,width,height,outline='black',width=self.bw)
    self.bar = self.canvas.create_rectangle(self.bw,self.bw,self.bw,self.cHeight,outline='#B05848',fill='#B05848',width=self.bw)
    self.update()
    
  def update(self):
  
    p = self.get()
    width = int(self.cWidth*p)
    self.canvas.coords(self.bar,self.bw,self.bw,width,self.cHeight)
    self.percent.set( ' %3.1d' % int(100*p) + "%" )
    self.label.set(self.text)
    self.update_idletasks() 
    if self.progress == self.total:
      time.sleep(0.1)
      self.close()
    
  def increment(self, n = 1):
  
    p = self.progress + n
    self.progress = max( min(self.total, p), 0 )
    self.update()
  
  def get(self):
  
    if self.total > 0:
      return min(self.progress/self.total, 1.0)
    else:
      return 1.0
  
  def set(self, value):

    if value <=1 :
      self.progress = float(value) * self.total
    elif value <= self.total:
      self.progress = float(value)
    else:
      raise 'Cannot set to a value exceeding the total'

    self.update()

  def setText(self, text=''):
    self.text = text
    self.label.set(self.text)


if (__name__ == '__main__'):

  from memops.gui.Button import Button

  root = Tkinter.Tk()
  pb = None

  def makePB():
    pb.set(0)
    pb.open()
    for i in range(100):
      time.sleep(0.01)
      pb.increment()
      
  b = Button(root, text="Hit Me", command=makePB)
  b.pack()
  pb = ProgressBar(root, text='Increments')

  root.mainloop()
