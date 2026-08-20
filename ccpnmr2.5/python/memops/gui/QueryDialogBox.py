
import Tkinter

from memops.gui.Base import Base
from memops.gui.Label           import Label
from memops.gui.Entry           import Entry
from memops.gui.Frame           import Frame
from memops.gui.ButtonList      import ButtonList
from memops.gui.MessageReporter import showWarning

class QueryDialogBox(Tkinter.Toplevel, Base):

  def __init__(self,parent,title='Dialog Box',prompt='Prompt',
              initialvalue=None, minvalue=None, maxvalue=None,
              position=(50,50),returnType=type('a'),hide=0 ):
  
    if parent is None:
      self.parent = Tkinter._default_root
    else:
      self.parent = parent
    
    Tkinter.Toplevel.__init__(self, parent)

    self.response = None
    self.prompt   = prompt
    self.initVal  = initialvalue
    self.minVal   = minvalue
    self.maxVal   = maxvalue
    self.hide     = hide
    
    self.returnType = returnType
    
    self.transient(parent)
    self.title(title or '')

    mainFrame = Tkinter.Frame(self)
    mainFrame.pack(padx=2, pady=2)
    
    self.drawButtons('OK','Cancel')
    self.entry = self.body(mainFrame)
    
    self.lift()     # so it is viewable
    self.grab_set() # make modal
     
    self.protocol('WM_WINDOW_DELETE', self.cancel)
    self.bind('<Return>', self.ok)
    self.bind('<Escape>', self.cancel)
    
    if parent is not None:
      dx = parent.winfo_rootx()+position[0]
      dy = parent.winfo_rooty()+position[1]
      self.geometry('+%d+%d' % (dx,dy))
    
    self.wait_window(self)

  def body(self, frame):
  
    label = Label(frame, text=self.prompt)
    
    if self.hide:
      entry = Entry(frame, text=self.initVal, bg='white', show='*')
    else:
      entry = Entry(frame, text=self.initVal, bg='white')
    
    label.grid(row=0, column=0, sticky=Tkinter.W)
    entry.grid(row=1, column=0, sticky=Tkinter.W)
    

    entry.focus_set() # for the keypresses
    return entry

  def cancel(self, *event):
    # event may come from keypress
  
    if self.parent:
      self.parent.focus_set()
      
    Tkinter.Toplevel.destroy(self)

  def ok(self, *event):
    # event may come from keypress

    self.response = self.getResponse()
    if self.response is not None:
      self.update_idletasks()
      self.destroy()
  
    else:
      self.parent.focus_set()
  
  def drawButtons(self, okText='OK', cancelText='Cancel'):
      
    texts    = [okText, cancelText]
    commands = [self.ok, self.cancel]
    
    buttonList = ButtonList(self, texts=texts, commands=commands)
    buttonList.buttons[0].config(default=Tkinter.ACTIVE)
    buttonList.pack()

  def getResponse(self):
  
    response = self.entry.get()
    if self.returnType == type(1):
      try:
        response = int(response)
      except:
        raise Exception('Value %s not an integer' % response)
      
      response = self.checkValue(response)
        
    elif self.returnType == type(1.0):
      try:
        response = float(response)
      except:
        raise Exception('Value %s not a floating point number' % response)

      response = self.checkValue(response)
    
    elif self.returnType == type('a'):
      response = response.strip()

    else:
      raise Exception('Unknown return type for query box')
    
    return response

  def checkValue(self, value):
  
    if (self.minVal is not None) and (value < self.minVal):
      showWarning('Warning', 'Value is less than minimum value (%s)' % self.minVal)
      value = None
    
    elif (self.maxVal is not None) and (value > self.maxVal):
      showWarning('Warning', 'Value is greater than maximum value (%s)' % self.maxVal)
      value = None
    
    return value

def askInteger(title, prompt, parent=None, **kw):
  kw['returnType'] = type(1)
  dialog = QueryDialogBox(parent, title, prompt, **kw)
  return dialog.response

def askFloat(title, prompt, parent=None, **kw):
  kw['returnType'] = type(1.0)
  dialog = QueryDialogBox(parent, title, prompt, **kw)
  return dialog.response

def askString(title, prompt, parent=None, **kw):
  kw['returnType'] = type('a')
  dialog = QueryDialogBox(parent, title, prompt, **kw)
  return dialog.response

def askPassword(title, prompt, parent=None, **kw):
  kw['returnType'] = type('a')
  dialog = QueryDialogBox(parent, title, prompt, hide=1, **kw)
  return dialog.response
