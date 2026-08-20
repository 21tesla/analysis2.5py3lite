
# Essentially the same as the ExtendNMR GUI for now. 


import os

# added jmci; do we need this??
import Tkinter
#from memops.universal.Io import getTopDirectory

#from ccpnmr.analysis.core.ExperimentBasic import getOnebondExpDimRefs
from ccpnmr.analysis.popups.BasePopup import BasePopup

from memops.editor.Util import createDismissHelpButtonList

from memops.gui.Frame import Frame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.Text import Text
from memops.gui.LinkChart import LinkChart
from memops.gui.LinkChart import PrototypeNode
from memops.gui.Tree import Tree
from memops.gui.Tree import Node
from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList

# # # # #   B I G   T O - D O  P O I N T S   # # # # # 
#


# # # # #   L A T E R   T O - D O  P O I N T S    # # # # #
#
 
class WorkflowPopup(BasePopup):

  def __init__(self, parent, ccpnProject):

    self.parent      = parent
    self.ccpnProject = ccpnProject
    
    BasePopup.__init__(self, parent=parent, title='WMS Workflow')

                       
  def body(self, guiFrame):


    frame = WorkflowFrame(guiFrame, basePopup=self, ccpnProject=self.ccpnProject)
    frame.grid(row=0, column=0, sticky='nsew')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    
    self.geometry('650x600')
    
    self.update_idletasks()
 

class WorkflowFrame(Frame):

  def __init__(self, guiParent, basePopup):
    # Base popup required to handle notification of data model changes
    # e.g. new peak lists, so that the GUI can update to the latest
    # state
    self.basePopup = basePopup
    self.guiParent = guiParent

    self.basePopup.frameShortcuts['Workflow'] = self

    #self.registerNotify=basePopup.registerNotify
    #self.unregisterNotify=basePopup.unregisterNotify

    Frame.__init__(self, guiParent)
  
    self.grid_rowconfigure(0, weight=0)
    self.grid_rowconfigure(1, weight=0, minsize=100)
    self.grid_rowconfigure(2, weight=1)
    self.grid_rowconfigure(3, weight=0, minsize=30)
    
    self.grid_columnconfigure(0, weight=0)
    self.grid_columnconfigure(1, weight=1)
    self.grid_columnconfigure(2, weight=0)

    # build up the body.

    self.lcA = LinkChart(self, background='#8080B0', height=110)
    self.lcB = LinkChart(self, background='#8080B0')

    lcButtonOpts=['Load','Save','Clear','New Protocol']
    lcButtonCmds=[self.tmpCall,self.tmpCall,self.tmpCall,self.tmpCall]
    self.lcButtons = ButtonList(self, lcButtonOpts, lcButtonCmds)

    # hack for now. This should come from db of meta-tasks
    self.expts=['Test1','Test2','ARIA','CING','ISD']

    # needs custom version
    self.filter = FilterFrame(self, self.basePopup, text='Filter')

    # no bean udnerneath for now so mock up nodes    
    self.wfTree = Tree(self, width=35)

    wfButtonOpts=['Load','Save']
    wfButtonCmds=[self.tmpCall,self.tmpCall]
    self.wfButtons = ButtonList(self, wfButtonOpts, wfButtonCmds)
    

    self.drawFrame()

  def drawFrame(self):

    # hack the WF Tree for now. just get is looking OK
    texts=['JMCI WF 1','JMCI WF 2']
    icons=[]
    parents=[]
    callbacks=[]
    objects=['JMCI WF 1','JMCI WF 2']

    for text in texts:
      icons.append('list-add')
      callbacks.append(None)
      parents.append(None)

    self.wfTree.update(parents, objects, texts, icons, callbacks)  

    self.lcA.grid(row=1, column=1, padx=5, pady=5, sticky='new')
    self.lcB.grid(row=2, column=1, padx=5, pady=5, sticky='nsew')
    self.lcButtons.grid(row=3, column=1, padx=5, pady=5, sticky='ew')

    self.filter.grid(row=1, column=2, sticky='new')
    self.wfTree.grid(row=2, column=2, padx=5, pady=5, sticky='nsew')
    self.wfButtons.grid(row=3, column=2, padx=5, pady=5, sticky='ew')

    self.nodes = []
    x = 0
    for expt in self.expts:

      node = PrototypeNode(self.lcA, title=expt, width=38,
                           updateFunc=None, resizeable=False, shape='SQUARE',
                           coords=(50+x*90, 50), object=None)

      x += 1
    
      
    self.lcA.draw()
    self.lcB.draw()

    

  def tmpCall(self):

    pass
    
  def administerNotifiers(self, notifyFunc):

      for func in ('__init__','delete','setName'):
        notifyFunc(self.updateAllAfter, 'ccp.nmr.Nmr.Experiment', func)
        notifyFunc(self.updateAllAfter, 'ccp.nmr.Nmr.DataSource', func)

  def updateAllAfter(self, obj):

    self.after_idle(self.updateAll)

  def updateAll(self, project=None):

    if project:
      self.project = project
      self.nmrProject = project.currentNmrProject

      if not self.nmrProject:
        self.nmrProject = project.newNmrProject(name=project.name)

    if not self.project:
      return


  def quit(self):
  
    self.guiParent.parent.destroy()
    
  def destroy(self):
  
    self.administerNotifiers(self.basePopup.unregisterNotify)
    Frame.destroy(self)



class FilterFrame(LabelFrame):

  def __init__(self, guiParent, basePopup, borderRelief='raised', text=' ', justify='left',
               width=None, font=None, height=None, *args, **kw):

    # Base popup required to handle notification of data model changes
    # e.g. new peak lists, so that the GUI can update to the latest
    # state
    self.basePopup = basePopup
    self.guiParent = guiParent
    self.borderRelief = borderRelief
    self.text = text
    self.justify = justify
    self.width = width
    self.font = font
    self.height = height


    self.registerNotify=basePopup.registerNotify
    self.unregisterNotify=basePopup.unregisterNotify


    LabelFrame.__init__(self, guiParent, borderRelief, text, justify,
                        width, font, height, *args, **kw)
  
    # set up the grid

    self.grid_columnconfigure(0, weight=0, minsize=20)

    self.grid_rowconfigure(0, weight=0, minsize=5)
    self.grid_rowconfigure(1, weight=0, minsize=10)
    self.grid_rowconfigure(2, weight=0, minsize=10)
    self.grid_rowconfigure(3, weight=0, minsize=10)

    # build up the body.

    self.nameLabel = Label(self,text='Name:')
    self.nameLabel.grid(row=1, column=0, padx=2, sticky='nw')

    self.nameText = Text(self, width=15, height=1)
    self.nameText.grid(row=1, column=1, sticky='nw')

    self.userLabel = Label(self,text='User:')
    self.userLabel.grid(row=2, column=0, padx=2, sticky='nw')

    self.userText = Text(self, width=15, height=1)
    self.userText.grid(row=2, column=1, sticky='nw')

    self.clearButton = Button(self, text='Clear', command=self.clear_filter, width=7, height=1)
    self.clearButton.grid(row=3, column=0, columnspan=2, sticky='nw')

    self.filterButton = Button(self, text='Go', command=self.tmpCall, width=7, height=1)
    self.filterButton.grid(row=3, column=1, columnspan=2, sticky='ne')


  def clear_filter(self):

    self.nameText.xlear()
    self.userText.clear()
    

  def tmpCall(self, event=None):

    pass



  def administerNotifiers(self, notifyFunc):

      for func in ('__init__','delete','setName'):
        notifyFunc(self.updateAllAfter, 'ccp.nmr.Nmr.Experiment', func)
        notifyFunc(self.updateAllAfter, 'ccp.nmr.Nmr.DataSource', func)

  def updateAllAfter(self, obj):

    self.after_idle(self.updateAll)

  def updateAll(self, project=None):

      return

  def quit(self):
  
    self.guiParent.parent.destroy()
    
  def destroy(self):
  
    self.administerNotifiers(self.basePopup.unregisterNotify)
    Frame.destroy(self)




if __name__ == "__main__":

  import sys
  import Tkinter

  root = Tkinter.Tk()
  root.withdraw()
   
  if len(sys.argv) == 2:
    path = sys.argv[1]
    from ccp.gui.Io import loadProject
    ccpnProject = loadProject(root, path=path)
  else:
    ccpnProject = None

  popup = WorkflowFramePopup(root, ccpnProject=ccpnProject)
  
  root.mainloop()
