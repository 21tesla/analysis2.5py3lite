
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
from memops.gui.Tree import Tree
from memops.gui.Tree import Node
from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList
from memops.gui.ScrolledMatrix import ScrolledMatrix

# # # # #   B I G   T O - D O  P O I N T S   # # # # # 
#


# # # # #   L A T E R   T O - D O  P O I N T S    # # # # #
#
 
class ProtocolPopup(BasePopup):

  def __init__(self, parent):

    self.parent      = parent
    
    BasePopup.__init__(self, parent=parent, title='WMS Protocol')

                       
  def body(self, guiFrame):


    frame = ProtocolFrame(guiFrame, basePopup=self)
    frame.grid(row=0, column=0, sticky='nsew')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    
    self.geometry('650x600')
    
    self.update_idletasks()
 

class ProtocolFrame(Frame):

  def __init__(self, guiParent, basePopup):
    # Base popup required to handle notification of data model changes
    # e.g. new peak lists, so that the GUI can update to the latest
    # state
    self.basePopup = basePopup
    self.guiParent = guiParent

    self.basePopup.frameShortcuts['Protocol'] = self

    Frame.__init__(self, guiParent)
  
    self.grid_rowconfigure(0, weight=0)
    self.grid_rowconfigure(1, weight=0, minsize=10)
    self.grid_rowconfigure(2, weight=0)
    self.grid_rowconfigure(3, weight=1)
    self.grid_rowconfigure(4, weight=0, minsize=30)
    self.grid_rowconfigure(5, weight=0, minsize=30)
    
    self.grid_columnconfigure(0, weight=0)
    self.grid_columnconfigure(1, weight=1)
    self.grid_columnconfigure(2, weight=0)

    # build up the body.

    self.title = Label(self, font='Helvetica16')
    self.details = Label(self)

    initial_cols = ['Parameter','Type','Value']
    self.pMatrix = ScrolledMatrix(self, headingList=initial_cols, initialRows=15)
    
    editButtonOpts=['Save','Clear','New Protocol','New Expt Type']
    editButtonCmds=[self.tmpCall,self.tmpCall,self.tmpCall,self.tmpCall]
    self.editButtons = ButtonList(self, editButtonOpts, editButtonCmds)


    # needs custom version
    self.filter = FilterFrame(self, self.basePopup, text='Filter')

    # no bean udnerneath for now so mock up nodes    
    self.pTree = Tree(self, width=33)

    pButtonOpts1=['Load','Save']
    pButtonCmds1=[self.tmpCall,self.tmpCall]
    self.pButtons1 = ButtonList(self, pButtonOpts1, pButtonCmds1)

    pButtonOpts2=['New','Delete']
    pButtonCmds2=[self.tmpCall,self.tmpCall]
    self.pButtons2 = ButtonList(self, pButtonOpts2, pButtonCmds2)
    
    self.openedLinks = {'Test 1': ['ParSet 1', 'ParSet 2'],
                        'Test 2': ['Std Pars', 'Big Pars', 'Little Pars'],
                        'ARIA': ['default', 'complex'],
                        'CING': ['Std'],
                        'ISD': ['protein','DNA'] }

    self.drawFrame()

  def drawFrame(self):

    
    # hack the WF Tree for now. just get is looking OK
    texts=['Test 1','Test 2', 'ARIA', 'CING', 'ISD']
    icons=[]
    parents=[]
    callbacks=[]
    objects=['Test 1','Test 2', 'ARIA', 'CING', 'ISD']

    for text in texts:
      icons.append('list-add')
      callbacks.append(self.openLink)
      parents.append(None)

    self.pTree.update(parents, objects, texts, icons, callbacks)  

    self.title.set('Protocol: Std (ARIA)')
    self.title.grid(row=1, column=1, padx=5, sticky='w')

    self.details.set('Description: details about the exp type or protocol')
    self.details.grid(row=2, column=1, padx=5, sticky='nw')

    self.pMatrix.grid(row=3, rowspan=2, column=1, padx=5, pady=5, sticky='ew')
    self.editButtons.grid(row=5, column=1, padx=5, pady=5, sticky='ew')

    self.filter.grid(row=1, column=2, rowspan=2, sticky='new')
    self.pTree.grid(row=3, column=2, padx=5, pady=5, sticky='nsew')
    self.pButtons1.grid(row=4, column=2, padx=5, sticky='ew')
    self.pButtons2.grid(row=5, column=2, padx=5, sticky='ew')

    matrix = [];
    objs = [];

    matrix.append(['Parameter1','String','foobar'])
    objs.append('obj1');
    matrix.append(['Parameter2','Number','10'])
    objs.append('obj2');
    
    self.pMatrix.update(objectList=objs, textMatrix=matrix)
    

  def tmpCall(self):

    pass


  def openLink(self, node):

    tree = node.tree
    par = node.object

    print 'opening node "', par, '"'
    print 'got series "', self.openedLinks[par], '"'

    if len(self.openedLinks[par]) > 0:

      for pv in self.openedLinks[par]:
    
        icon='folder'
        obj = par + '::' + pv
        text = pv
        callback=None
        tree.add(node, obj, text, icon, callback)

    
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

    self.nameLabel = Label(self,text='Exp Type:')
    self.nameLabel.grid(row=1, column=0, padx=2, sticky='nw')

    self.nameText = Text(self, width=15, height=1)
    self.nameText.grid(row=1, column=1, sticky='nw')

    self.clearButton = Button(self, text='Clear', command=self.clear_filter, width=7, height=1)
    self.clearButton.grid(row=2, column=0, columnspan=2, sticky='nw')

    self.filterButton = Button(self, text='Go', command=self.tmpCall, width=7, height=1)
    self.filterButton.grid(row=2, column=1, columnspan=2, sticky='ne')


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

  popup = ProtocolFramePopup(root, ccpnProject=ccpnProject)
  
  root.mainloop()
