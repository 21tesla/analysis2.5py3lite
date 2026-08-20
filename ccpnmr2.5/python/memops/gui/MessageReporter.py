
"""
======================COPYRIGHT/LICENSE START==========================

MessageReporter.py: <write function here>

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
import Tkinter, tkMessageBox


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NOTE:ED - simple message queuing for the minute
class _Messages(object):
  """Simple class to queue a set of errors that can be 'flushed' as a compound set of errors at the end
  """
  def __init__(self):
    self._initialise()

  def setParent(self, parent):
    self._parent = parent

  def _initialise(self):
    self._parent = None
    self._errorTypes = {}

  def getMessages(self, errType, title=None):
    if errType in self._errorTypes:
      err = self._errorTypes[errType]

      if title in err:
        return '\n'.join(err[title])

      elif not title:
        return '\n'.join(['\n'.join(err[tt]) for tt in err])

    return ''

  def addMessage(self, errType, title, message, parent):
    if parent:
      self.setParent(parent)

    if errType not in self._errorTypes:
      self._errorTypes[errType] = {}
    err = self._errorTypes[errType]

    if title not in err:
      err[title] = (message,)
    else:
      err[title] += (message,)

  def flushMessages(self):
    for errType, errorList in self._errorTypes.items():
      for title in errorList:
        if self._parent:
          errType(title, self.getMessages(errType, title), parent=self._parent)
        else:
          errType(title, self.getMessages(errType, title))

    self._initialise()

_queueErrors = False
_queueMessages = _Messages()

def _flushMessages():
  _queueMessages.flushMessages()
  _queueErrors = False
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def showError(title, message, parent = None):

  if not _queueErrors:
    if (parent):
      tkMessageBox.showerror(title, message, parent=parent)
    else:
      tkMessageBox.showerror(title, message)

  else:
    _queueMessages.addMessage(tkMessageBox.showerror, title, message, parent)


def showInfo(title, message, parent = None):

  if not _queueErrors:
    if (parent):
      tkMessageBox.showinfo(title, message, parent=parent)
    else:
      tkMessageBox.showinfo(title, message)

  else:
    _queueMessages.addMessage(tkMessageBox.showinfo, title, message, parent)


def showOkCancel(title, message, parent = None):

  if (parent):
    return tkMessageBox.askokcancel(title, message, parent=parent)
  else:
    return tkMessageBox.askokcancel(title, message)


def showYesNo(title, message, parent = None):

  if (parent):
    return tkMessageBox.askyesno(title, message, parent=parent)
  else:
    return tkMessageBox.askyesno(title, message)


def showWarning(title, message, parent = None):

  # if (parent):
  #   return tkMessageBox.showwarning(title, message, parent=parent)
  # else:
  #   return tkMessageBox.showwarning(title, message)

  # NOTE:ED - added simple queuing just for the minute
  if not _queueErrors:
    if (parent):
      return tkMessageBox.showwarning(title, message, parent=parent)
    else:
      return tkMessageBox.showwarning(title, message)

  else:
    _queueMessages.addMessage(tkMessageBox.showwarning, title, message, parent)


def showMulti(title, message, texts, objects=None, parent=None):

  popup = MultiChoice(title, message, texts, objects, parent)
  func = popup.get()
  popup.destroy()
  return func


class MessageReporter:

  def showError(self, title, message, parent = None, *args, **kw):
 
    showError(title, message, parent)
 
  def showInfo(self, title, message, parent = None, *args, **kw):
 
    showInfo(title, message, parent)
 
  def showWarning(self, title, message, parent = None, *args, **kw):

    showWarning(title, message, parent)

  def showOkCancel(self, title, message, parent = None, *args, **kw):
 
    return showOkCancel(title, message, parent)
 
  def showYesNo(self, title, message, parent = None, *args, **kw):
 
    return showYesNo(title, message, parent)

from memops.gui.Base import Base 
from memops.gui.Button import Button 
from memops.gui.Frame import Frame
from memops.gui.Label import Label

class MultiChoice(Tkinter.Toplevel, Base):

   def __init__(self, title, message, texts, objects=None, parent=None):
     
     if parent is None:
       parent = Tkinter.Tk()
       #parent.withdraw()
       self.root = parent
       parent.protocol('WM_DELETE_WINDOW', self.destroy)

     else:
       self.root = None

     Tkinter.Toplevel.__init__(self, parent)  
     
     self.resizable(0,0)
     self.parent = parent
     self.title(title)
     self.var = Tkinter.IntVar()
     self.objects = objects or range(len(texts)) 
    
     assert len(self.objects) == len(texts)

     x = parent.winfo_rootx() + parent.winfo_width()/2
     y = parent.winfo_rooty() + parent.winfo_height()/2
     location = '+%d+%d' % (x, y)
     
     if hasattr(parent, 'font'):
       self.font = parent.font
     else:
       self.font = None
       
     Tkinter.Toplevel.geometry(self, location)
     Tkinter.Toplevel.lift(self, parent)
     self.update_idletasks()
     self.transient(parent)

     self.protocol('WM_DELETE_WINDOW', self._null)
     self.focus_force()
     self.grab_set()
     
     self.grid_rowconfigure(0, weight=1)
     self.grid_columnconfigure(1, weight=1)
     
     label = Label(self, text=message, grid=(0,0))
     self.config(bg=label.cget('bg'))
     
     frame = Frame(self, grid=(1,0))
     
     _buttons = []
     for i, text in enumerate(texts):
       button = Button(frame, text=text,
                       command=lambda j=i: self._click(j))
       button.grid(row=i+1, column=0, padx=2, pady=2, sticky='ew')
       if objects[i] is None:
         button.config(bg='#B0FFB0')
       _buttons.append(button)
       
     self._wait()

   def _null(self):
     
     pass
   
   def get(self):
   
     i = self.var.get()
     if i < 0:
       # been destroyed
       return None
       
     self.grab_release()
     #self.destroy()
     
     return self.objects[i]
     
   def _wait(self):
   
     self.update_idletasks()
     self.wait_variable(self.var)
   
   def _click(self, index):
     
     self.var.set(index)

   def destroy(self):
     
     self.var.set(-1)
 
     Tkinter.Toplevel.destroy(self)
     
     if self.root:
       self.root.destroy()
     else:
       self.parent.focus_set()

messageReporter = MessageReporter()

if (__name__ == '__main__'):

  s = showMulti('Multiple Choice', 'Select one of the following:',
                ['Option 1','Option 2','Option 3','Cancel'],
                [1,2,3,None])
  print "Multi:", s
 
  showError('title', 'error message')
  s = showOkCancel('title', 'ok message')
  print s, type(s)
  s = showOkCancel('title', 'ok2 message')
  print s, type(s)
  s = showYesNo('title', 'yes message')
  print s, type(s)
  s = showYesNo('title', 'yes2 message')
  print s, type(s)
  showWarning('title', 'warning message')
  showInfo('title', 'info message')
  messageReporter.showError('error title', 'error message')
  messageReporter.showInfo('info title', 'info message')
  messageReporter.showWarning('warning title', 'warning message')
  print messageReporter.showOkCancel('ok cancel title', 'ok cancel message')
  print messageReporter.showYesNo('yes no title', 'yes no message')
