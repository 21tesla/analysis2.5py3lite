
"""
======================COPYRIGHT/LICENSE START==========================

MultiWidget.py: 

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
import os.path as path

from memops.gui.Frame  import Frame
from memops.gui.Label  import Label
from memops.gui.Button import Button
from memops.universal.Io import getTopDirectory
from memops.gui.ToolTip import ToolTip

blueColor = '#D0D0FF'

class MultiWidget(Frame):

  def __init__(self, parent, widgetClass, options=None, values=None,
               callback=None, minRows=3, maxRows=None, useImages=True,
               wrap=15, docKey=None, tipText=None, defaultValue=None,
               *args, **kw):

    Frame.__init__(self, parent, docKey=docKey, tipText=tipText, createToolTip=True, **kw)

    if not values:
      values = []
    if not options:
      options = [] 

    # Widget types currently: Entries, Pulldown, Checkbutton

    if useImages:
      gfxDir = path.join(getTopDirectory(),'python','memops','gui','graphics')
      self.crossImage = Tkinter.PhotoImage(file=path.join(gfxDir,'cross.gif'))
      self.tickImage  = Tkinter.PhotoImage(file=path.join(gfxDir,'tick.gif'))
      self.plusImage  = Tkinter.PhotoImage(file=path.join(gfxDir,'plus.gif'))
      self.minusImage = Tkinter.PhotoImage(file=path.join(gfxDir,'minus.gif'))
    else:
      self.crossImage = None
      self.tickImage  = None
      self.plusImage  = None
      self.minusImage = None
    
    module = widgetClass.__module__
    self.widgetType  = module.split('.')[-1]
    self.widgetClass = widgetClass
    self.callback    = callback
    self.wrap = wrap
    self.defaultValue = defaultValue
    
    self.numRows = 0
    self.values  = values
    self.options = options
    self.minRows = minRows
    self.maxRows = maxRows
    self.widgets = []
    self.labels  = []
    self.grid_columnconfigure(1, weight=1)

    self.topFrame = Frame(self)
    self.topFrame.grid(row=0, column=0, columnspan=4, sticky=Tkinter.NSEW)
    self.topFrame.grid_columnconfigure(1, weight=1)
    self.topFrame.grid_rowconfigure(0, weight=1)

    self.addButton    = None
    self.removeButton = None
    self.commitButton = None
    self.cancelButton = None
    self.allButton = None
    self.noneButton = None

    self.set(self.values, self.options)
  
    self.updateButtons()
    
    self.bind('<KeyPress-Return>', self.commit)
    self.bind('<KeyPress-Escape>', self.cancel)
    
  def get(self):
  
    values = []
    
    if self.widgetType == 'PulldownMenu':
      for i in range(self.numRows):
        value = self.widgets[i].getSelected()
        if value is not None:
          values.append( value )
          
    elif self.widgetType == 'PulldownList':
      for i in range(self.numRows):
        value = self.widgets[i].getObject()
        if value is not None:
          values.append( value )
          
    else:
      for i in range(self.numRows):
        value = self.widgets[i].get()
        if (value is not None) and (value != ''):
          values.append( value )
        elif self.defaultValue is not None:
          values.append( self.defaultValue )
        
    return values

  def set(self, values=None, options=None):


    if values is not None:
      self.values = values
    if options is not None:
      self.options = options      
    
    
    N = max( len(self.values), self.minRows)
    if self.maxRows is not None:
      N = min(N, self.maxRows)
    
    #N = max( min(len(self.values), self.maxRows), self.minRows)
    #N = len(self.values)
    if self.numRows < N:
      for i in range(self.numRows, N):
        self.addRow()

    elif self.numRows > N:
      for i in range(N, self.numRows):
        self.removeRow(doClear=True)

    for i in range(self.numRows):
      if i < len(self.values):
        value = self.values[i]
      else:
        value = None
        
      if i < len(self.options):
        option = self.options[i]
      else:
        option = None

      self.updateWidget(i, value, option)

    self.updateButtons()
      
  def setOptions(self, options):
  
    self.set(self.values, options)


  def setValues(self, values):
  
    self.set(values, self.options)


  def updateWidget(self, row, value, option):
 
    if self.labels:
      font = self.labels[0].cget('font')
    else:
      font = None
 
    widget = self.widgets[row]
    if self.widgetType in ('Entry','IntEntry','FloatEntry'):
      widget.set(value)
      label = self.labels[row]
      label.set(option or '')
      if font:
        label.config(font=font)
 
    elif self.widgetType == 'PulldownMenu':
      index = -1
      if value in self.options:
        index = self.options.index(value or self.options[0])
      widget.setup(self.options, index)

 
    elif self.widgetType == 'PulldownList':
      index = 0
      if value in self.options:
        index = self.options.index(value or self.options[0])
      widget.setup(self.options, self.options, index)

    elif self.widgetType == 'CheckButton':
      widget.set(value)
      label = self.labels[row]
      label.set(option or '')
      if font:
        label.config(font=font)

    else:
      raise 'Widget type %s not supported in MultiWidget' % self.widgetType
    

  def updateButtons(self):

    row = 1
    col = 0
    if self.widgetType != 'CheckButton':
      if (self.maxRows is None) or (self.numRows < self.maxRows):
        if not self.addButton:
          if self.plusImage:
            self.addButton = Button(self, image=self.plusImage, docKey='+', tipText='Add extra row',
                                   command=self.addRow, borderwidth=1)
          else:
            self.addButton = Button(self, text='+', command=self.addRow, tipText='Add extra row',
                                    borderwidth=1, bg=blueColor)
            
        self.addButton.grid(row=row, column=col, sticky=Tkinter.EW)
        col += 1

      elif self.addButton:
        self.addButton.grid_forget()

      if self.numRows > self.minRows:
        if not self.removeButton:    
          if self.minusImage:
            self.removeButton = Button(self, image=self.minusImage, docKey='-', tipText='Remove last row',
                                       command=self.removeRow, borderwidth=1)
          else:
            self.removeButton = Button(self, text='-', command=self.removeRow, tipText='Remove last row',
                                       borderwidth=1, bg=blueColor)
        self.removeButton.grid(row=row, column=col, sticky=Tkinter.EW)
        col += 1

      elif self.removeButton:
        self.removeButton.grid_forget()
    
    elif self.numRows > 5:
      
      if not self.allButton:
        self.allButton = Button(self, text='All', command=self._selectAll, borderwidth=1,
                                bg=blueColor, tipText='Select all options')
      
      if not self.noneButton:
        self.noneButton = Button(self, text='None', command=self._selectNone, borderwidth=1,
                                 bg=blueColor, tipText='Select no options')
      
      self.allButton.grid(row=row, column=col, sticky=Tkinter.EW)
      col += 1
      self.noneButton.grid(row=row, column=col, sticky=Tkinter.EW)
      col += 1

    elif self.allButton:
      self.allButton.grid_forget()
      self.noneButton.grid_forget()
      
    if self.callback:
      if not self.commitButton:
        if self.tickImage:
          self.commitButton = Button(self, image=self.tickImage, docKey='OK', tipText='Commit selection',
                                     command=self.commit, borderwidth=1)
        else:
          self.commitButton = Button(self, text='OK', command=self.commit, tipText='Commit selection',
                                     borderwidth=1, bg=blueColor)
 
      self.commitButton.grid(row=row, column=col, sticky=Tkinter.EW)
      self.grid_columnconfigure(col, weight=1)
      col += 1
 
      if not self.cancelButton:
        if self.crossImage:
          self.cancelButton = Button(self, image=self.crossImage, command=self.cancel, 
                                     borderwidth=1, docKey='Cancel', tipText='Abort selection')
        else:
          self.cancelButton = Button(self, text='Cancel', command=self.cancel, borderwidth=1,
                                     bg=blueColor, tipText='Abort selection')
 
      self.cancelButton.grid(row=row, column=col, sticky=Tkinter.EW)
      self.grid_columnconfigure(col, weight=1)

        
  def _selectAll(self):
  
    values = [True] * len(self.values)
    self.set(values, self.options)

  def _selectNone(self):
  
    values = [False] * len(self.values)
    self.set(values, self.options)


  def addWidget(self, i, value=None):
    
    row = i % self.wrap
    col = i // self.wrap
 
    if self.widgetType in ('Entry','IntEntry','FloatEntry'):
      option = ''
      if i < len(self.options):
        option = self.options[i]

      label = Label(self.topFrame, text=option or '', anchor='w')
      label.grid(row=row, column=2*col, sticky=Tkinter.E)
      self.labels.append(label)
      
      widget = self.widgetClass(self.topFrame, text=value, width=12)
      widget.grid(row=row, column=2*col+1, sticky=Tkinter.EW)
 
      widget.bind('<KeyPress-Return>', self.commit)
      widget.bind('<KeyPress-Escape>', self.cancel)
 

    elif self.widgetType == 'PulldownList':
      index = -1
      if value in self.options:
        index = self.options.index(value or self.options[0])
      widget = self.widgetClass(self.topFrame)
      widget.setup(self.options, self.options, index)
      widget.grid(row=row, column=2*col, columnspan=2, sticky=Tkinter.W)

    elif self.widgetType == 'PulldownMenu':
      index = -1
      if value in self.options:
        index = self.options.index(value or self.options[0])
      widget = self.widgetClass(self.topFrame)
      widget.setup(self.options, index)
      widget.grid(row=row, column=2*col, columnspan=2, sticky=Tkinter.W)

    elif self.widgetType == 'CheckButton':
      widget = self.widgetClass(self.topFrame)
      widget.set(value)
      widget.grid(row=row, column=2*col, sticky=Tkinter.W)
      
      option = ''
      if row < len(self.options):
        option = self.options[i]

      label = Label(self.topFrame, text=option or '', anchor='w')
      label.grid(row=row, column=2*col+1, sticky=Tkinter.W)
      self.labels.append(label)

    else:
      raise 'Widget type %s not supported in MultiWidget' % self.widgetType

    return widget


  def commit(self, event=None):

    values = self.get()
          
    if self.callback:
      self.callback(values)


  def cancel(self, event=None):
  
    if self.callback:
      self.callback(None)


  def addRow(self):

    if (self.maxRows is None) or (self.numRows < self.maxRows):
      i = self.numRows
      self.numRows += 1
      row = i % self.wrap
      col = i // self.wrap
      col2 = 2*col
      
      if self.numRows > len(self.widgets):
        widget = self.addWidget(i)
        self.widgets.append(widget)

      else:
        widget = self.widgets[i]
        if self.widgetType in ('Entry','IntEntry','FloatEntry'):
          label = self.labels[i]
          label.grid(row=row, column=col2, sticky=Tkinter.E)
          widget.grid(row=row, column=col2+1, sticky=Tkinter.EW)
            
        elif self.widgetType == 'CheckButton':
          widget.grid(row=row, column=col2, sticky=Tkinter.W)
          label = self.labels[i]
          label.grid(row=row, column=col2+1, sticky=Tkinter.W)
        
        else:
          widget.grid(row=row, column=col2, columnspan=2, sticky=Tkinter.W)

      if row == 0:
        self.topFrame.grid(row=0, column=0, columnspan=4, sticky=Tkinter.NSEW)


      self.updateButtons()

  def removeRow(self, doClear=0):
 
    if self.numRows > self.minRows:
      self.numRows -= 1
      row    = self.numRows
      widget = self.widgets[row]
      widget.grid_forget()
      
      if self.widgetType == 'CheckButton':
        label = self.labels[row]
        label.grid_forget()
        
      elif self.widgetType in ('Entry','IntEntry','FloatEntry'):
        label = self.labels[row]
        label.grid_forget()

      if doClear:
        self.updateWidget(row, None, '') 

      if self.numRows == 0:
        self.topFrame.grid_forget()

      self.updateButtons()

if __name__ == '__main__':

  from memops.gui.Entry        import Entry
  from memops.gui.FloatEntry   import FloatEntry
  from memops.gui.PulldownMenu import PulldownMenu
  from memops.gui.CheckButton  import CheckButton

  def fn(values):
    print values

  root = Tkinter.Tk()
  mw = MultiWidget(root, CheckButton, callback=fn, minRows=0,
                   options=['One','Two','Three'], values=[1,0,1],
                   tipText='Some advice')
  mw.grid()
  
  mw.set(options=['Easy','As','One','Two','Three'],values=[1,0,1,0,1])

  mw2 = MultiWidget(root, FloatEntry, callback=fn, minRows=3, maxRows=3, options=['1H:','13C:','15N:'], values=[1.0,0.1,0.15], useImages=False)
  mw2.grid()
  
  mw2 = MultiWidget(root, Entry, callback=fn, minRows=0, values=['Easy','As','One','Two','Three'], useImages=False)
  mw2.grid()

  mw3 = MultiWidget(root, PulldownMenu, callback=None, values=['Bjorn','Benny','Benny','Bjorn'], options=['Bjorn','Agnetha','Benny','Anni-Frid'])
  mw3.grid()
  
  root.mainloop()
