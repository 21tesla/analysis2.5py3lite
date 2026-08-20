
"""
======================COPYRIGHT/LICENSE START==========================

ColorChooser.py: <write function here>

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
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.Color import hsbToRgb
from memops.gui.LabeledScale import LabeledScale
from memops.universal.Util import buttonClick

class ColorChooser(Frame):

  def __init__(self, parent, *args, **kw):

    apply(Frame.__init__, (self, parent) + args, kw)

    self.grid_rowconfigure(0, weight=1)
    self.grid_rowconfigure(1, weight=0)
    self.grid_rowconfigure(2, weight=0)
    self.grid_rowconfigure(3, weight=0)

    self.canvas = Canvas(self, relief='sunken', borderwidth=2, width=150)
    self.canvas.grid(row=0, column=1, columnspan=1, sticky=Tkinter.NSEW,pady=2,padx=2)

    frame = Frame(self)
    frame.grid_columnconfigure(1, weight=1)
    self.boxSize = 100
    self.numPixels = 10
    self.pixelWidth = round(self.boxSize/self.numPixels)
    self.colorBox = Canvas(frame, relief='sunken', borderwidth=2, width=self.boxSize+2, height=self.boxSize+2)
    self.colorBox.bind('<Button-1>', self.pickInColorBox)

    # self.colorBox.bind('<Button-2>', self.pickInColorBox)
    # self.colorBox.bind('<Button-3>', self.pickInColorBox)
    self.colorBox.bind(buttonClick(2), self.pickInColorBox)
    self.colorBox.bind(buttonClick(3), self.pickInColorBox)

    self.colorBox.grid(row=0, column=2, rowspan=3, sticky=Tkinter.NSEW,padx=4,pady=4)
    self.pixel = []
    self.colors = []
    self.setupColorBox()
    self.scale = Tkinter.Scale(frame, orient=Tkinter.VERTICAL, length=self.boxSize,from_=0, to=99,label='',
                               showvalue=0, command=self.refreshColorBox )
    self.scale.grid(row=0, column=3, rowspan=3, sticky=Tkinter.NS, padx=4,pady=4)
    frame.grid(row=1, column=0, columnspan=2, sticky=Tkinter.NSEW)

    labels = ('Red', 'Green', 'Blue')
    self.labeled_scale = 3*[None]

    for n in range(3):

      label = Label(frame, text=labels[n]+':', anchor=Tkinter.W)
      label.grid(row=n, column=0, sticky=Tkinter.EW)

      self.labeled_scale[n] = LabeledScale(frame, values = range(101),
                                           label_format = '%3d',
                                           set_callback = self.scaleCallback)
      self.labeled_scale[n].grid(row=n, column=1, sticky=Tkinter.EW)

  def setupColorBox(self):

    pw = self.pixelWidth
    np = self.numPixels
    for x in range(np):
       self.pixel.append( [] )
       self.colors.append( [] )
       for y in range(np):
         rgb = hsbToRgb(0,1-x/float(np),1-y/float(np))
         c = '#%02x%02x%02x' % (rgb[0]*255,rgb[1]*255,rgb[2]*255)
         p = self.colorBox.create_rectangle(2+pw*x,2+pw*y,2+(pw*x)+pw,(2+pw*y)+pw, outline=c, fill=c)
         self.pixel[x].append(p)
         self.colors[x].append(rgb)

  def pickInColorBox(self,event):
  
    x = int(event.x/self.pixelWidth)
    y = int(event.y/self.pixelWidth)
    X = min(x,self.numPixels-1)
    Y = min(y,self.numPixels-1)
    x = max(X,0)
    y = max(Y,0)
    
    rgb = self.colors[x][y]
    self.setColor(rgb[0],rgb[1],rgb[2])

  def refreshColorBox(self,H):
 
    np = self.numPixels
    H = int(H)/100.
    for x in range(np):
      for y in range(np):
        rgb = hsbToRgb(H,1-x/float(np),1-y/float(np))
        c = '#%02x%02x%02x' % (rgb[0]*255,rgb[1]*255,rgb[2]*255)
        self.colorBox.itemconfigure(self.pixel[x][y], outline=c, fill=c)
        self.colors[x][y] = rgb

  def scaleCallback(self, *index):

    c = 3*[0]
    for n in range(3):
      c[n] = int(2.559 * float(self.labeled_scale[n].getValue()))

    color = '#%02x%02x%02x' % (c[0], c[1], c[2])
    self.canvas.config(bg=color)

  def setColor(self, r, g, b):

    c = (r, g, b)
    for n in range(3):
      s = round(100.0 * c[n])
      s = min(100, s)
      self.labeled_scale[n].set(s)

  def getColor(self):

    c = 3*[0]
    for n in range(3):
      c[n] = float(self.labeled_scale[n].getValue()) / 100.0

    return tuple(c)

if __name__ == '__main__':
 
  import sys
  from memops.gui.Button import Button
 
  def quit():

    global color_chooser

    print 'color:', color_chooser.getColor()
    sys.exit(0)

  root = Tkinter.Tk()
 
  color_chooser = ColorChooser(root, width=500, height=500)
  color_chooser.pack(side=Tkinter.TOP, expand=Tkinter.YES, fill=Tkinter.BOTH)
 
  button = Button(root, text='quit', command=quit)
  button.pack(side=Tkinter.TOP)
 
  root.mainloop()
