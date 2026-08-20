
"""
======================COPYRIGHT/LICENSE START==========================

PopupTemplate.py: <write function here>

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

from memops.gui.BasePopup             import BasePopup  
from memops.gui.ButtonList            import ButtonList
from memops.gui.Button                import Button
from memops.gui.CheckButton           import CheckButton
from memops.gui.DataEntry             import askString
from memops.gui.Entry                 import Entry
from memops.gui.FileSelectPopup       import FileSelectPopup
from memops.gui.Frame                 import Frame
from memops.gui.LabelFrame            import LabelFrame
from memops.gui.Label                 import Label
from memops.gui.MessageReporter       import showWarning, showYesNo
from memops.gui.ScrolledMatrix        import ScrolledMatrix
from memops.gui.Text                  import Text
from memops.gui.CheckButtons          import CheckButtons
from memops.gui.PartitionedSelector   import PartitionedSelector
from memops.gui.RadioButtons          import RadioButtons
from memops.gui.Scale                 import Scale
from memops.gui.ValueRamp             import ValueRamp
from memops.gui.ScrolledFrame         import ScrolledFrame
from memops.gui.ScrolledGraph         import ScrolledGraph
from memops.gui.ToggleLabel           import ToggleLabel
from memops.gui.ScrolledCanvas        import ScrolledCanvas

from memops.gui.PulldownMenu          import PulldownMenu

from memops.gui.FloatEntry            import FloatEntry
from memops.gui.IntEntry              import IntEntry
from memops.gui.MultiWidget           import MultiWidget
from memops.gui.ProgressBar           import ProgressBar

class Fruit:

  def __init__(self, name, color, quantity, keywords=None):
  
    if keywords is None:
      keywords = []

    self.name  = name
    self.color = color
    self.quantity = quantity
    self.keywords = keywords

class PopupTemplate(BasePopup):

  """
  **Popup Template**

  This popup still needs its documentation.
  """

  def __init__(self, parent, project=None, *args, **kw):

    self.project = project
    self.parent  = parent
    self.objects = self.getObjects()
    self.object  = None
    
    BasePopup.__init__(self, parent=parent, title='Popup Template', **kw)
                       
    self.updateObjects()

  def body(self, mainFrame):

    mainFrame.grid_columnconfigure(1, weight=1, minsize=100)
    mainFrame.config(borderwidth=5, relief='solid')

    row = 0
    label = Label(mainFrame, text="Frame (with sub-widgets):")
    label.grid(row=row, column=0, sticky=Tkinter.E)

    frame = Frame(mainFrame, relief='raised', border=2, background='#8080D0')
    # Frame expands East-West
    frame.grid(row=row, column=1, sticky=Tkinter.EW)
    # Last column expands => Widgets pusted to the West
    frame.grid_columnconfigure(3, weight=1)
    
    # Label is within the sub frame
    label = Label(frame, text='label ')
    label.grid(row=0, column=0, sticky=Tkinter.W)
    
    entry = Entry(frame, text='Entry', returnCallback=self.showWarning,
                  tipText='Hint to user')
    entry.grid(row=0, column=1, sticky=Tkinter.W)
    
    self.check = CheckButton(frame, text='Checkbutton', selected=True,
                             callback=self.updateObjects,
                             tipText='Hint to user')
    self.check.grid(row=0, column=2, sticky=Tkinter.W)
    
    # stick a button to the East wall
    button = Button(frame, text='Button', command=self.pressButton,
                    tipText='Hint to user')
    button.grid(row=0, column=3, sticky=Tkinter.E)
  
    row += 1
    label = Label(mainFrame, text="Text:")
    label.grid(row=row, column=0, sticky=Tkinter.E)
    self.textWindow = Text(mainFrame, text='Initial Text\n', width=60, height=5)
    self.textWindow.grid(row=row, column=1, sticky=Tkinter.NSEW)
    
    row += 1
    label = Label(mainFrame, text="CheckButtons:")
    label.grid(row=row, column=0, sticky=Tkinter.E)

    entries = ['Alpha','Beta','Gamma','Delta']
    selected = (False, False, True, True)
    self.checkButtons = CheckButtons(mainFrame, entries, selected=selected,
                                     select_callback=self.changedCheckButtons)
    self.checkButtons.grid(row=row, column=1, sticky=Tkinter.W)
  
    row += 1
    label = Label(mainFrame, text="PartitionedSelector:")
    label.grid(row=row, column=0, sticky=Tkinter.E)

    labels   = ['Bool','Int','Float','String']
    objects  = [type(0),type(1),type(1.0),type('a')]
    selected = [type('a')]
    self.partitionedSelector= PartitionedSelector(mainFrame, labels=labels,
                                                  objects=objects,
                                                  colors = ['red','yellow','green','#000080'],
                                                  callback=self.toggleSelector,selected=selected)
    self.partitionedSelector.grid(row=row, column=1, sticky=Tkinter.EW)

    row += 1
    label = Label(mainFrame, text="PulldownMenu")
    label.grid(row=row, column=0, sticky=Tkinter.E)
    
    entries = ['Frodo','Pipin','Merry','Sam','Bill','Gandalf','Strider','Gimli','Legolas']
    self.pulldownMenu = PulldownMenu(mainFrame, callback=self.selectPulldown,
                                     entries=entries, selected_index=2,
                                     do_initial_callback=False,
                                     tipText='Hint to user')
    self.pulldownMenu.grid(row=row, column=1, sticky=Tkinter.W)

    row += 1
    label = Label(mainFrame, text="RadioButtons in a\nScrolledFrame.frame:")
    label.grid(row=row, column=0, sticky=Tkinter.EW)
    
    frame = ScrolledFrame(mainFrame, yscroll = False, doExtraConfig = True, width=100)
    frame.grid(row=row, column=1, sticky=Tkinter.EW)
    frame.grid_columnconfigure(0, weight=1)

    self.radioButtons = RadioButtons(frame.frame, entries=entries,
                                     select_callback=self.checkRadioButtons,
                                     selected_index=1, relief='groove')
    self.radioButtons.grid(row=0, column=0, sticky=Tkinter.W)
    
    row += 1
    label = Label(mainFrame, text="LabelFrame with\nToggleLabels inside:")
    label.grid(row=row, column=0, sticky=Tkinter.E)

    labelFrame = LabelFrame(mainFrame, text='Frame Title')
    labelFrame.grid(row=row, column=1, sticky=Tkinter.NSEW)
    labelFrame.grid_rowconfigure(0, weight=1)
    labelFrame.grid_columnconfigure(3, weight=1)
    
        
    self.toggleLabel1 = ToggleLabel(labelFrame, text='ScrolledMatrix', callback=self.toggleFrame1)
    self.toggleLabel1.grid(row=0, column=0, sticky=Tkinter.W)
    self.toggleLabel1.arrowOn()

    self.toggleLabel2 = ToggleLabel(labelFrame, text='ScrolledGraph', callback=self.toggleFrame2)
    self.toggleLabel2.grid(row=0, column=1, sticky=Tkinter.W)

    self.toggleLabel3 = ToggleLabel(labelFrame, text='ScrolledCanvas', callback=self.toggleFrame3)
    self.toggleLabel3.grid(row=0, column=2, sticky=Tkinter.W)
    
    row += 1
    mainFrame.grid_rowconfigure(row, weight=1)

    label = Label(mainFrame, text="changing/shrinking frames:")
    label.grid(row=row, column=0, sticky=Tkinter.E)
    
    self.toggleRow = row
    self.toggleFrame = Frame(mainFrame)
    self.toggleFrame.grid(row=row, column=1, sticky=Tkinter.NSEW)
    self.toggleFrame.grid_rowconfigure(0, weight=1)
    self.toggleFrame.grid_columnconfigure(0, weight=1)
    
    # option 1
    
    self.intEntry = IntEntry(self, returnCallback = self.setNumber, width=8)
    
    self.multiWidget = MultiWidget(self, Entry, options=None, 
                                  values=None, callback=self.setKeywords,
                                  minRows=3, maxRows=5)

    editWidgets      = [None, None, self.intEntry,  self.multiWidget]
    editGetCallbacks = [None, None, self.getNumber, self.getKeywords]
    editSetCallbacks = [None, None, self.setNumber, self.setKeywords]
    
    tipTexts = ['Name of the whatsit','Colour of the whatsit',
                'The whatsit\'s number','The whatsit\'s keywords']
    headingList = ['Name','Colour','Number','Keywords']
    self.scrolledMatrix = ScrolledMatrix(self.toggleFrame, headingList=headingList,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets,
                                         tipTexts=tipTexts,
                                         callback=self.selectObject,
                                         multiSelect=False) 
                                         
    self.scrolledMatrix.grid(row=0, column=0, sticky=Tkinter.NSEW)

    # option 2
    self.scrolledGraph = ScrolledGraph(self.toggleFrame, width=400,
                                       height=300, symbolSize=5,
                                       symbols=['square','circle'],
                                       dataColors=['#000080','#800000'],
                                       lineWidths=[0,1] )

    self.scrolledGraph.setZoom(1.3)

    dataSet1 = [[0,0],[1,1],[2,4],[3,9],[4,16],[5,25]]
    dataSet2 = [[0,0],[1,3],[2,6],[3,9],[4,12],[5,15]]
    self.scrolledGraph.update(dataSets=[dataSet1,dataSet2],
                              xLabel = 'X axis label',
                              yLabel = 'Y axis label',
                              title  = 'Main Title')
    self.scrolledGraph.draw()

    # option 3
    self.scrolledCanvas = ScrolledCanvas(self.toggleFrame,relief = 'groove', borderwidth = 2, resizeCallback=None)
    canvas = self.scrolledCanvas.canvas
    font   = 'Helvetica 10'
    box    = canvas.create_rectangle(10,10,150,200, outline='grey', fill='grey90')
    line   = canvas.create_line(0,0,200,200,fill='#800000', width=2)
    text   = canvas.create_text(120,50, text='Text', font=font, fill='black')
    circle = canvas.create_oval(30,30,50,50,outline='#008000',fill='#404040',width=3)
     
    row += 1
    label = Label(mainFrame, text="FloatEntry:")
    label.grid(row=row, column=0, sticky=Tkinter.E)
    self.floatEntry = FloatEntry(mainFrame, text=3.14159265, tipText='Hint to user',
                                 returnCallback=self.floatEntryReturn)
    self.floatEntry.grid(row=row, column=1, sticky=Tkinter.W)
    
     
    row += 1
    label = Label(mainFrame, text="Scale:")
    label.grid(row=row, column=0, sticky=Tkinter.E)
    self.scale = Scale(mainFrame, from_=10, to=90, value=50, orient=Tkinter.HORIZONTAL)
    self.scale.grid(row=row, column=1, sticky=Tkinter.W)

    row += 1
    label = Label(mainFrame, text="Value Ramp:")
    label.grid(row=row, column=0, sticky=Tkinter.E)
    self.valueRamp = ValueRamp(mainFrame, self.valueRampCallback, speed = 1.5, delay = 50)
    self.valueRamp.grid(row=row, column=1, sticky=Tkinter.W)
  

    row += 1
    label = Label(mainFrame, text="ButtonList:")
    label.grid(row=row, column=0, sticky=Tkinter.E)
    
    texts    = ['Select File','Close','Quit']
    commands = [self.selectFile, self.close, self.quit]
    bottomButtons = ButtonList(mainFrame, texts=texts, commands=commands, expands=True) 
    bottomButtons.grid(row=row, column=1, sticky=Tkinter.EW)
  
    self.protocol('WM_DELETE_WINDOW', self.quit)

  def floatEntryReturn(self, event):
  
    value = self.floatEntry.get()
    self.textWindow.setText('%s\n' % value)

  def selectObject(self, object, row, col):
  
    self.object = object

  def getKeywords(self, object):
  
    if object :
      values = object.keywords
      self.multiWidget.set(values)
  
  def setKeywords(self, event):
  
    values = self.multiWidget.get()
    self.object.keywords = values
    self.updateObjects()

  def getNumber(self, object):

    if object :
      self.intEntry.set(object.quantity)

  def setNumber(self, event):

    value = self.intEntry.get()
    self.object.quantity = value
    self.updateObjects()

  def toggleFrame1(self, isHidden):
 
    if isHidden:
      self.scrolledMatrix.grid_forget()
      self.toggleFrame.grid_forget()
    else:
      self.scrolledGraph.grid_forget()
      self.scrolledCanvas.grid_forget()
      self.scrolledMatrix.grid(row=0, column=0, sticky=Tkinter.NSEW)
      self.toggleFrame.grid(row=self.toggleRow, column=1,sticky=Tkinter.NSEW)
      self.toggleLabel2.arrowOff()
      self.toggleLabel3.arrowOff()


  def toggleFrame2(self, isHidden):
 
    if isHidden:
      self.scrolledGraph.grid_forget()
      self.toggleFrame.grid_forget()
    else:
      self.scrolledMatrix.grid_forget()
      self.scrolledCanvas.grid_forget()
      self.scrolledGraph.grid(row=0, column=0, sticky=Tkinter.NSEW)
      self.toggleFrame.grid(row=self.toggleRow, column=1,sticky=Tkinter.NSEW)
      self.toggleLabel1.arrowOff()
      self.toggleLabel3.arrowOff()

  def toggleFrame3(self, isHidden):
 
    if isHidden:
      self.scrolledCanvas.grid_forget()
      self.toggleFrame.grid_forget()
    else:
      self.scrolledMatrix.grid_forget()
      self.scrolledGraph.grid_forget()
      self.scrolledCanvas.grid(row=0, column=0, sticky=Tkinter.NSEW)
      self.toggleFrame.grid(row=self.toggleRow, column=1,sticky=Tkinter.NSEW)
      self.toggleLabel1.arrowOff()
      self.toggleLabel2.arrowOff()

 
  def valueRampCallback(self, value):
  
    self.textWindow.setText('%s\n' % value)

  def checkRadioButtons(self, value):
  
    self.textWindow.setText('%s\n' % value)

  def selectPulldown(self, index, name):
  
    self.textWindow.setText('%d, %s\n' % (index, name))

  def toggleSelector(self, value):
  
    self.textWindow.setText('%s\n' % value)

  def changedCheckButtons(self, values):
  
    self.textWindow.setText(','.join(values) + '\n')

  def getObjects(self):
  
    objects = []
    
    objects.append( Fruit('Lemon',   '#FFFF00',1,keywords=['Bitter','Tangy'] ) )
    objects.append( Fruit('Orange',  '#FF8000',4 ) )
    objects.append( Fruit('Banana',  '#FFF000',5 ) )
    objects.append( Fruit('Pinapple','#FFD000',9 ) )
    objects.append( Fruit('Kiwi',    '#008000',12) )
    objects.append( Fruit('Lime',    '#00FF00',2 ) )
    objects.append( Fruit('Apple',   '#800000',5,keywords=['Crunchy'] ) )
    objects.append( Fruit('Pear',    '#408000',6 ) )
    objects.append( Fruit('Peach',   '#FFE0C0',2,keywords=['Sweet','Furry'] ) )
    objects.append( Fruit('Plumb',   '#800080',7 ) )
    
    return objects

  def updateObjects(self, event=None):

    textMatrix = []
    objectList = []
    colorMatrix = []
    
    for object in self.objects:
      datum = []
      datum.append( object.name )
      datum.append( None )
      datum.append( object.quantity )
      datum.append( ','.join(object.keywords) )
    
      colors = [None, object.color, None, None]
    
      textMatrix.append(datum)
      objectList.append(object)
      colorMatrix.append(colors)

    if self.check.get():
      self.scrolledMatrix.update(textMatrix=textMatrix, objectList=objectList)
    else:
      self.scrolledMatrix.update(textMatrix=textMatrix, objectList=objectList, colorMatrix=colorMatrix)

  def selectFile(self):

    fileSelectPopup = FileSelectPopup(self, title = 'Choose file', dismiss_text = 'Cancel',
                                      selected_file_must_exist = True)

    fileName = fileSelectPopup.getFile()

    self.textWindow.setText('File Selected: %s\n' % fileName)

  def showWarning(self, eventObject):
  
    self.textWindow.setText('Text Entry Return Pressed\n')
    showWarning('Warning Title','Warning Message')
    return

  def pressButton(self):
  
    self.textWindow.setText('Button Pressed\n')
    if showYesNo('Title','Prompt: Clear text window?'):
       self.textWindow.clear()
       
    return

  def quit(self):
  
    BasePopup.destroy(self)  

if __name__ == '__main__':

  root = Tkinter.Tk()
  root.withdraw()

  popup = PopupTemplate(root)
