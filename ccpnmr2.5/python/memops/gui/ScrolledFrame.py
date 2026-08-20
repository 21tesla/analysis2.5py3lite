
"""
======================COPYRIGHT/LICENSE START==========================

ScrolledFrame.py: <write function here>

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
# creates scrolled frame with attribute frame which is used
# from outside this code as the parent of the frame to be scrolled
# this only works if the size of that frame is fixed at construction
import Tkinter
from memops.gui.Frame import Frame

class ScrolledFrame(Frame):

  def __init__(self, parent, xscroll=True, yscroll=True,
               width=500, height=500, doExtraConfig=True,
               tipText=None, *args, **kw):

    Frame.__init__(self, parent, *args, **kw)

    self.grid_rowconfigure(0, weight=1)
    self.grid_columnconfigure(0, weight=1)

    self.doExtraConfig = doExtraConfig
    self.xscroll = xscroll
    self.yscroll = yscroll

    self.canvas = Tkinter.Canvas(self, borderwidth=0, bg=self.cget('bg'), highlightthickness=0, width=25, height=25)
    self.canvas.grid(row=0, column=0, sticky=Tkinter.NSEW)

    if (xscroll):
      self.xscrollbar = Tkinter.Scrollbar(self, orient=Tkinter.HORIZONTAL, borderwidth=1)
      self.xscrollbar.config(command=self.canvas.xview)
      self.canvas.config(xscrollcommand=self.xscrollbar.set)

    if (yscroll):
      self.yscrollbar = Tkinter.Scrollbar(self, orient=Tkinter.VERTICAL, borderwidth=1)
      self.yscrollbar.config(command=self.canvas.yview)
      self.canvas.config(yscrollcommand=self.yscrollbar.set)

    # reset the view (always do this if you don't use scrollbars)
    self.canvas.xview('moveto', 0)
    self.canvas.yview('moveto', 0)

    # create the frame
    self.frame = Frame(self.canvas, tipText=tipText)
    self.frame.grid(row=0, column=0, sticky=Tkinter.NSEW)
    self.frame.grid_rowconfigure(0, weight=1)
    self.frame.grid_columnconfigure(0, weight=1)
    self.waiting = 0

    # track changes to its size
    self.bind('<Configure>', self.configureCanvasAfter)

    # track changes to its contents
    self.frame.bind('<Configure>', self.configureFrame)
 
    # place the frame inside the canvas (this also
    # runs the __configure method)
    self.canvas.create_window(0, 0, window=self.frame, anchor='nw')
    #self.frameWidth  = self.frame.winfo_reqwidth()
    #self.frameHeight = self.frame.winfo_reqheight()
    self.frameWidth  = 0
    self.frameHeight = 0
    
  def configureFrame(self, event):
  
    # frame contents have changed size, but make no change to displayed canvas size
    if (event.width != self.frameWidth) or (event.height != self.frameHeight):

      # the size of the frame with widgets in the 'background'
      self.frameWidth  = event.width 
      self.frameHeight = event.height

      # the scrolled region is always the entire frame
      self.canvas.config(scrollregion="0 0 %s %s" % (self.frameWidth, self.frameHeight) )

      self.configureCanvasAfter()
 
  def configureCanvasAfter(self, event=None):

    if self.waiting:
      return
      
    else:
      self.waiting = 1
      self.after_idle(lambda: self.configureCanvas(event))

  def configureCanvas(self, event=None):
     
    d = 2 * int(self.cget('borderwidth'))
    if (event):
      # the canvas visible display size updates with the rezize event unless the axis is not scrolled
      canvasWidth  = event.width - d
      canvasHeight = event.height - d
      
    else:
      # if no resize event the canvas size is unaltered
      canvasWidth  = int(self.winfo_width()) - d
      canvasHeight = int(self.winfo_height()) - d
      
    if self.xscroll:
      if (self.frameWidth > 1):
        if (self.frameWidth > canvasWidth):
          self.xscrollbar.grid(row=1, column=0, sticky=Tkinter.EW)
        else:
          self.xscrollbar.grid_forget()
      else:
         self.xscrollbar.grid(row=1, column=0, sticky=Tkinter.EW)
    else:
      # display full width
      canvasWidth = self.frameWidth

    if self.yscroll:
      if (self.frameHeight > 1):
        if (self.frameHeight > canvasHeight):
          self.yscrollbar.grid(row=0, column=1, sticky=Tkinter.NS)
        else:
          self.yscrollbar.grid_forget()
      else:
        self.yscrollbar.grid(row=0, column=1, sticky=Tkinter.NS)        
    else:
      # display full height
      canvasHeight = self.frameHeight
 
    if (self.doExtraConfig):
      # sometimes below seems to be needed and sometimes not
      # seems to depend on what else is there in widget besides scrolled frame
      # can sometimes wildly iterate if below is done
      self.canvas.config(width=canvasWidth, height=canvasHeight)
    self.waiting = 0

if (__name__ == '__main__'):

  import Util

  def config1(event):
    print 'config1:', frame.winfo_width(), frame.winfo_height()

  def config2(event):
    print 'config2:', frame2.winfo_width(), frame2.winfo_height()

  root = Tkinter.Tk()

  frame = ScrolledFrame(root)
  #frame = ScrolledFrame(root, xscroll=false)
  # below must be pack not grid in order for expansion to work
  frame.pack(side=Tkinter.TOP, expand=Tkinter.YES, fill = Tkinter.BOTH)
  frame.bind('<Configure>', config1)

  frame2 = Frame(frame.frame)
  frame2.bind('<Configure>', config2)

  for i in range(30):
    t = 'a very long label %s' % i 
    label = Tkinter.Label(frame2, text=t)
    label.grid(sticky=Tkinter.NSEW)

  frame2.grid(sticky=Tkinter.NSEW)

  button = Util.createQuitButton(root)
  button.pack()

  root.mainloop()
