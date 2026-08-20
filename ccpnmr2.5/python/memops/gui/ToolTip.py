
"""
======================COPYRIGHT/LICENSE START==========================

ToolTip.py: <write function here>

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

FONT = "Helvetica 9"

class ToolTip(object):

  def __init__(self, parent, text='*Blank Tool Tip*', delay=1500, *args, **kw):
  
           
    self.parent = parent
    self.delay = delay
    self._text = text
    self._tag = None
    self._panel = None
    self._isAqua = parent.tk.call("tk", "windowingsystem") == 'aqua'
  
    if not text:
      return
    
    self._justifyText()
    
    parent.bind('<Enter>', self.activate, '+')
    parent.bind('<Leave>', self._mouseLeave, '+')
    parent.bind('<ButtonPress>', self._mouseLeave, '+')
    parent.bind('<Unmap>', self._mouseLeave, '+')
    parent.bind('<Destroy>', self._mouseLeave, '+')
  
  def setText(self, text):
  
    self._text = text
    self._justifyText()
    self.deactivate()

  def getText(self):
  
    return self._text
    
  text = property(getText, setText)
  
  def _justifyText(self):
    
    text = self._text.replace('\n', ' ')
  
    if len(text) > 60:
      m = len(text) // 2
      i = text[m:].find(' ')
      
      if i >= 0:
        n = m + i
        first = text[:n].rstrip()
        last = text[n:].lstrip()
        
        self._text = first + '\n' + last
      
  def _mouseLeave(self, event=None):

    self.deactivate()
    
  def activate(self, event=None):
  
    #self.deactivate()
    if self._tag:
      return # already active
    
    self._tag = self.parent.after(self.delay, self._appear)
  
  def deactivate(self):
  
    if self._tag:
      self.parent.after_cancel(self._tag)
      self._tag = None  
    
    self._disappear()
  
  def _appear(self):
  
    if self._text and not self._panel:
      parent = self.parent
      panel = self._panel = Tkinter.Toplevel(self.parent)
      panel.withdraw()
      panel.wm_overrideredirect(True)
      
      if self._isAqua:
        panel.tk.call("::tk::unsupported::MacWindowStyle",
                      "style", panel._w, "help", "none")
      
      label = Tkinter.Label(panel, text=self._text, bd=1, padx=2, pady=2,
                            relief='solid', bg='#FFFFB0', justify='left',
                            font=FONT)
      label.bind('<ButtonPress>', self._mouseLeave)
      label.pack()
      
      panel.update_idletasks()
      width = panel.winfo_reqwidth()
      height = panel.winfo_reqheight()
      screenWidth = panel.winfo_screenwidth()
      screenHeight = panel.winfo_screenheight()
      parentY = self.parent.winfo_rooty()
      
      if isinstance(parent, Tkinter.Menu):
        y = panel.winfo_pointery() + height/2
      elif isinstance(parent, Tkinter.Frame):
        y = panel.winfo_pointery() + height
      else:
        y = parentY + self.parent.winfo_height() + 3
      
      x = panel.winfo_pointerx() - (width/2)
      
      if y + height > screenHeight:
        y = parentY - height - 3
      
      x = max(0, min(x, screenWidth-width) )
      #y = max(0, min(y, screenHeight-height))

      panel.wm_geometry('+%d+%d' % (x, y))
      panel.deiconify()
      
  def _disappear(self):
  
    if self._panel:
      self._panel.destroy()
      self._panel = None
  


if __name__ == '__main__':
  
  from memops.gui.Button import Button
  
  root = Tkinter.Tk()
  button = Button(root, text='Close (hover for tip)', command=root.destroy,
                  toolTip='This will close the window')
  button.pack()

  root.mainloop()
