
"""
======================COPYRIGHT/LICENSE START==========================

BasePopup.py: <write function here>

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

from memops.gui.Base import Base
from memops.gui.Button import Button
from memops.gui.Label import Label
from memops.gui.Menu import Menu
from memops.gui.Frame import Frame

from memops.universal.Util import isWindowsOS, isMacOS

# base_popups = list of popups created using BasePopup()
# BasePopup can be subclassed or using createPopup() instead
# this never seems to have been used
# and is dangerous since popups never get garbage collected
#base_popups = []

class BasePopup(Tkinter.Toplevel, Base):

  def __init__(self, parent, title='', location='', hide=False, font=None,
               modal=False, transient=False, quitFunc=None,
               docKey=None, tipText=None, createToolTip=False, *args, **kw):

    # NOTE:ED - can only be non-transient if modal as well
    #           transient means the windows are tied to the app
    #           on MacOS his means that they cannot be raised/lowered
    if isMacOS():
      transient = modal

    self.parent = parent
    self.location = location
    self.isTransient = transient
    self.modal = modal

    if (not kw.has_key('bg') and not kw.has_key('background')):
      kw['bg'] = 'grey90'

    Tkinter.Toplevel.__init__(self, parent, *args, **kw)
    Base.__init__(self, docKey=docKey, tipText=tipText, createToolTip=createToolTip)

    if modal:
      self.config(borderwidth=5, bg='red')

    if (not location):
      x = parent.winfo_rootx() + 50
      y = parent.winfo_rooty() + 50
      location = '+%d+%d' % (x, y)

    if hasattr(parent, 'top'):
      self.top = parent.top
    else:
      self.top = self

    if (parent and transient):
      self.transient(parent)

    self.var = Tkinter.StringVar()

    self.setTitle(title)

    #self.root = self._root()

    self.geometry(location)

    if hide:
      self.withdraw()

    self.protocol('WM_DELETE_WINDOW', self.close)

    #self.bind('<Escape>', self.close)

    self.grid_rowconfigure(0, weight=1)
    self.grid_columnconfigure(0, weight=1)
    frame = Frame(self)
    frame.grid(row=0, column=0, sticky=Tkinter.NSEW)

    if quitFunc:
      self.protocol('WM_DELETE_WINDOW', quitFunc)

    self.initial_focus = self.body(frame)
    if not self.initial_focus:
      self.initial_focus = self
    self.initial_focus.focus_set()

    #base_popups.append(self)
 
    self.help_popup = None
    
    if not hasattr(self, 'font'):
      if font is None:
        if hasattr(parent, 'font'):
          font = parent.font
 
      self.font = font

    self.setFont()

    if modal:
      self.do_grab()
      self.wait_variable(self.var)
  
    if isWindowsOS():
      # Fix so that Windows popups aren't empty when opened
      self.update_idletasks()
      geometry = self.geometry()
      size = geometry.split('+')[0]
      x,y = size.split('x')
      x = str(int(x)+1)
      self.geometry('%sx%s' % (x,y))
      
  def waitCursor(self):
  
    self.configure(cursor="watch")
    self.after_idle(lambda :self.configure(cursor=""))
      
  def setFont(self, font=None):
  
    if font:
      self.font = font
    else:
      font = self.font  
  
    if font:
      childList = self.children.values()
 
      classes = [Tkinter.Button,
                 Tkinter.Label,
                 Tkinter.Menu,
                 Tkinter.Entry,
                 Tkinter.Checkbutton,
                 Tkinter.Radiobutton]
 
      for child in childList:
      
        for clazz in classes:
          if isinstance(child, clazz):
            if hasattr(child, 'font'):
              if not child.font:
                child.config( font=font)
            else:
              child.config( font=font )
            break
 
        if isinstance(child, Tkinter.Frame):
          child.font = font
 
        if hasattr(child, 'setFont'):
          child.setFont(font)
  
        if hasattr(child, 'children'):
          childList.extend( child.children.values() )

 
  def do_grab(self):

    while 1:
      try:
        self.grab_set()
        break
      except Tkinter.TclError:
        # Another application has grab.  Keep trying until grab can succeed.
        time.sleep(0.05) # 50 msec

  def open(self):

    if (self.state() == 'normal'):
      Tkinter.Toplevel.lift(self)
      self.lift()
    else:
      self.deiconify()

    if (self.modal):
      self.config(borderwidth=5, bg='red')
      self.do_grab()
      self.wait_variable(self.var)
    else:
      self.config(borderwidth=0)

  def ok(self, *event):

    if (not self.apply()):
      self.initial_focus.focus_set() # put focus back
      return
 
    self.close()
 
  def close(self, *event):
    
    if (self.state() == 'normal'):
      self.updateLocation()
      self.withdraw()

    if (self.modal):
      self.grab_release()

    self.var.set('')
    self.update_idletasks()

    if (self.parent):
      self.parent.focus_set()

  def tkraise(self, *above):

    self.lift(self, *above)

  def lift(self, *above):

    if (above and self is above):
      return

    """ this never seems to have been used
    try:

      ind = base_popups.index(self) # do this to to check exists

      if (above):
        ind = base_popups.index(above) + 1
      else:
        ind = len(base_popups)

      base_popups.remove(self)
      base_popups.insert(ind, self)

    except:

      pass # bad but not much else can do
    """

    Tkinter.Toplevel.lift(self, *above)

  def lower(self, *below):

    if (below and self is below):
      return

    """ this never seems to have been used
    try:

      ind = base_popups.index(self) # do this to to check exists

      if (below):
        ind = base_popups.index(below)
      else:
        ind = 0

      base_popups.remove(self)
      base_popups.insert(ind, self)

    except:

      pass # bad but not much else can do
    """

    Tkinter.Toplevel.lower(self, *below)
 
  def geometry(self, location=''):
    
    if location:
      if 'x' in location:
        if '+' in location:
          width, rest = location.split('x')
          height, x, y = rest.split('+')
        else:
          width, height = location.split('x')
          x = self.winfo_rootx()
          y = self.winfo_rooty()
        
      elif '+' in location:
        height, x, y = location.split('+')
        width = self.winfo_width()
        height = self.winfo_height()
      
      width, height, x, y = [int(v) for v in (width, height, x, y)]
       
      sWidth = self.winfo_screenwidth()
      sHeight = self.winfo_screenheight()
      
      width = min(sWidth, width)
      height = min(sHeight, height)
      
      if (x+width) > sWidth:
        if width == 1:
          x = sWidth // 2
        else:
          x = sWidth - width
      elif x < 0:
        x = 0
      
      if (y+height) > sHeight:
        if height == 1:
          y = sHeight // 2
        else:
          y = sHeight - height  
      elif y < 0:
        y = 0
      
      if width == height == 1:
         location = '+%d+%d' % (x, y)
      
      else:
         location = '%dx%d+%d+%d' % (width, height, x, y)
      
      # sometimes width or height < 0
      if width > 0 and height > 0:
        Tkinter.Toplevel.geometry(self, location)
       
        # guard against this routine being called before packing done
        if location != '1x1+0+0':
          self.updateLocation()
    else:
      return Tkinter.Toplevel.geometry(self)

  def updateLocation(self):

    s = self.geometry()
    n = s.index('+')
    self.location = s[n:]

  def getGeometry(self):

    g = self.geometry()
    (w, rest) = g.split('x')
    (h, x, y) = rest.split('+')
 
    return map(int, (w, h, x, y))

  def setGeometry(self, w, h, x, y):
 
    g = '%sx%s+%s+%s' % (w, h, x, y)
    self.geometry(g)

  def iconify(self):

    self.updateLocation()
    Tkinter.Toplevel.iconify(self)

    if (self.modal):
      self.grab_release()

  def deiconify(self):

    if (self.location):
      self.geometry(self.location)

    self.lift() # might happen automatically
    Tkinter.Toplevel.deiconify(self)

    if (self.modal):
      self.do_grab()
    
  def withdraw(self):

    Tkinter.Toplevel.withdraw(self)

    if (self.modal):
      self.grab_release()

  def destroy(self):

    """ this never seems to have been used
    try:
      base_popups.remove(self)
    except:
      pass
    """

    if (self.modal):
      self.grab_release()

    try:
      if (self.help_popup):
        self.help_popup.destroy()
    except:
      pass

    Tkinter.Toplevel.destroy(self)

  def setTitle(self, title = ''):

    self.title(title)
    # The documentation for iconname says:
    # Set the icon name to use when this window is iconified.
    # This method is ignored by some window managers (including Windows).
    self.iconname(title)

  def __getattr__(self, name):
 
    try:
      if self.__dict__.has_key('frame'):
        return getattr(self.__dict__['frame'], name)
    except:
      pass

    raise AttributeError, "%s instance has no attribute '%s'" % (self.__class__.__name__, name)

  def body(self, master):
    pass # this method should be overridden by subclass

  def apply(self):
    return True # this method can be overridden by subclass

def createPopup(parent, frameClass, title = '', location = '',
                hide = False, modal = False,
                popup_args = [], popup_kw = {},
                *frame_args, **frame_kw):

  #print 'frameClass =', frameClass
  #print 'title =', title
  #print 'location =', location
  #print 'hide =', hide
  #print 'popup_args =', popup_args
  #print 'popup_kw =', popup_kw
  #print 'frame_args =', frame_args
  #print 'frame_kw =', frame_kw

  popup = BasePopup(parent, title=title, location=location, hide=hide, modal=modal, *popup_args, **popup_kw)

  frame = frameClass(popup, *frame_args, **frame_kw)

  popup.grid_rowconfigure(0, weight=1)
  popup.grid_columnconfigure(0, weight=1)

  frame.grid(row=0, column=0, sticky=Tkinter.NSEW)

  popup.frame = frame

  return popup

if (__name__ == '__main__'):

  import sys
  #from Frame import Frame

  class TestFrame(Frame):

    def __init__(self, parent, *args, **kw):

      apply(Frame.__init__, (self, parent) + args, kw)

      button = Button(self, text='close', command=parent.close)
      button.pack(side=Tkinter.LEFT)

      button = Button(self, text='iconify', command=parent.iconify)
      button.pack(side=Tkinter.LEFT)

  class TestPopup(BasePopup):

    def __init__(self, root, *args, **kw):
      kw['modal'] = True
      kw['title'] = 'Test popup'
      self.result = None
      apply(BasePopup.__init__, (self, root) + args, kw)

    def body(self, master):

      self.geometry('600x400+50+50')
      label = Tkinter.Label(master, text='label 1')
      label.grid(row=0, column=0, sticky=Tkinter.EW)
      label = Tkinter.Label(master, text='label 2')
      label.grid(row=1, column=0, sticky=Tkinter.EW)
      button = Button(master, text='ok', command=self.ok)
      button.grid(row=2, column=0, sticky=Tkinter.EW)
      button = Button(master, text='cancel', command=self.close)
      button.grid(row=2, column=1, sticky=Tkinter.EW)


    def apply(self):
      self.result = 77
      return True

  popup = None

  root = Tkinter.Tk()

  def new():

    global popup, root
    popup = createPopup(root, TestFrame)

  def new2():

    global popup, root
    popup = TestPopup(root)
    print popup.result

  def lift():

    if (popup):
      popup.lift()

  def close():

    if (popup):
      popup.close()

  def open():

    if (popup):
      popup.open()

  def iconify():

    if (popup):
      popup.iconify()

  def deiconify():

    if (popup):
      popup.deiconify()

  root.protocol('WM_DELETE_WINDOW', quit)

  button = Button(root, text='new popup', command=new)
  button.pack(side=Tkinter.TOP)

  button = Button(root, text='new2 popup', command=new2)
  button.pack(side=Tkinter.TOP)

  button = Button(root, text='lift popup', command=lift)
  button.pack(side=Tkinter.TOP)

  button = Button(root, text='open popup', command=open)
  button.pack(side=Tkinter.TOP)

  button = Button(root, text='close popup', command=close)
  button.pack(side=Tkinter.TOP)

  button = Button(root, text='iconify popup', command=iconify)
  button.pack(side=Tkinter.TOP)

  button = Button(root, text='deiconify popup', command=deiconify)
  button.pack(side=Tkinter.TOP)

  button = Button(root, text='quit', command=sys.exit)
  button.pack(side=Tkinter.TOP)

  root.mainloop()
