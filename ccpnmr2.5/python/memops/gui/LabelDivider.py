
"""
======================COPYRIGHT/LICENSE START==========================

LabelFrame.py: <write function here>

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

from memops.gui.Color import hexRepr
from memops.gui.Canvas import Canvas
from memops.gui.Frame import Frame
from memops.gui.Label import Label

class LabelDivider(Frame):

  def __init__(self, parent, borderRelief='raised', text=' ', justify='left',
               font=None, sticky='ew', docKey=None, tipText=None, *args, **kw):

    Frame.__init__(self, parent, sticky=sticky, *args, **kw)

    self.borderRelief  = borderRelief
    self.bg       = self.cget('bg')
    self.justify  = justify
    self.canvas1  = Canvas(self, background=self.bg)

    if self.borderRelief  == 'sunken':
      fill1 = 'grey60'
      fill2 = 'white'
    else:
      fill1 = 'white'
      fill2 = 'grey60'
    fill3 = 'grey85'
      
    self.label = Label(self, text=text, font=font, docKey=docKey, tipText=tipText)
    self.label.grid()
    
    self.canvasL = [self.canvas1.create_rectangle(0,0,0,0,width=0,fill=fill1),
                    self.canvas1.create_rectangle(0,0,0,0,width=0,fill=fill2),
                    self.canvas1.create_rectangle(0,0,0,0,width=0,fill=fill3)]
                    
    self.canvasR = [self.canvas1.create_rectangle(0,0,0,0,width=0,fill=fill1),
                    self.canvas1.create_rectangle(0,0,0,0,width=0,fill=fill2),
                    self.canvas1.create_rectangle(0,0,0,0,width=0,fill=fill3)]    
    
    self.bind('<Configure>', self.resize)
    
    
  def setText(self, text):
  
    self.label.set(text)
    self.resize()
        
  def resize(self, event=None):
   
    if event:
      w = event.width
      h = event.height
    else:
      w = int(self.winfo_width())
      h = int(self.winfo_height())
     
    self.after_idle(lambda: self.draw(w, h))
     
  def draw(self, w, h):
    
    bw  = int(self.label.winfo_height())
    pad = bw/2
    h2  = h-(2*bw)

    textWidth = int(self.label.winfo_reqwidth())
    
    if self.justify == 'left':
      textLeft  = (3*pad) + (0.05 * w)
      textRight = min(w,textLeft + textWidth)
    elif self.justify == 'right':
      textRight = (3*pad) + (0.05 * w)
      textLeft  = w - (textWidth + textRight)
    else:
      textLeft  = (w/2) - (textWidth/2)
      textRight = textLeft + textWidth
    
    c1 = self.canvas1
    c1coords = c1.coords
    c1.place(   x=0,   y=0, anchor='nw', width=w,  height=bw)
    
    self.label.place(x=textLeft, y=0, anchor='nw')
 
    e1,e2,e3 = self.canvasL
    c1coords(e1,   pad,   pad, textLeft-pad+3, pad+3)
    c1coords(e2, pad+1, pad+1, textLeft-pad+3, pad+3)
    c1coords(e3, pad+1, pad+1, textLeft-pad+2, pad+2)

    e1,e2,e3 = self.canvasR
    c1coords(e1,   textRight+pad,   pad, w-pad+3, pad+3)
    c1coords(e2, textRight+pad+1, pad+1, w-pad+3, pad+3)
    c1coords(e3, textRight+pad+1, pad+1, w-pad+2, pad+2)
    

    
