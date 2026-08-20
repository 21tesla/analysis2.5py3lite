
# required for WS layer
from SharedBeanService_services import *
from WSString import *

import os

# added jmci; do we need this??
import Tkinter
from memops.universal.Io import getTopDirectory

#from ccpnmr.analysis.core.ExperimentBasic import getOnebondExpDimRefs
from ccpnmr.analysis.popups.BasePopup import BasePopup

from memops.editor.Util import createDismissHelpButtonList

from memops.gui.Frame import Frame
from memops.gui.Tree import Tree
from memops.gui.Label import Label



# # # # #   B I G   T O - D O  P O I N T S   # # # # # 
#


# # # # #   L A T E R   T O - D O  P O I N T S    # # # # #
#
 
class DataExchMasterPopup(BasePopup):

  def __init__(self, parent):

    self.parent      = parent
    
    BasePopup.__init__(self, parent=parent, title='PRODECOMP')

                       
  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    frame = DataExchFrame(guiFrame, basePopup=self)
    frame.grid(row=0, column=0, sticky='nsew')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    
    self.geometry('650x600')
    
    self.update_idletasks()
 

class DataExchFrame(Frame):

  def __init__(self, guiParent, basePopup, rep1=None, rep2=None):
    
    self.basePopup = basePopup
    self.guiParent = guiParent

    # add this to shortcuts to ease navigation    
    self.basePopup.frameShortcuts['DataExch'] = self

    # FIXME

    # This is a hack just so that we can see what it would be like
    # need to set these in a different way!!

    self.repository1 = self.basePopup.repList.currentRepository
    self.repository2 = self.basePopup.repList.currentRepository

    Frame.__init__(self, guiParent)
  
    self.grid_columnconfigure(0, weight=0, minsize=10)
    self.grid_columnconfigure(1, weight=1, minsize=20)
    self.grid_columnconfigure(2, weight=0, minsize=40)
    self.grid_columnconfigure(3, weight=1, minsize=20)
    self.grid_columnconfigure(4, weight=0, minsize=10)

    self.grid_rowconfigure(0, weight=0, minsize=10)
    self.grid_rowconfigure(1, weight=0, minsize=10)
    self.grid_rowconfigure(2, weight=0, minsize=10)
    self.grid_rowconfigure(3, weight=1, minsize=10)
    self.grid_rowconfigure(4, weight=0, minsize=10)

    # build up the body. start really simple

    self.repLabel1 = Label(self, text='Repository 1')
    self.repLabel2 = Label(self, text='Repository 2')

    self.noDataWidgets = []

    self.noRep1Label = Label(self,text='No repository currently selected.')
    self.noDataWidgets.append(self.noRep1Label)
    self.noRep2Label = Label(self,text='No repository currently selected.')
    self.noDataWidgets.append(self.noRep2Label)

    # widgets for view when current repository set

    self.dataWidgets1 = []

    self.label1 = Label(self)
    self.dataWidgets1.append(self.label1)
    self.frame1 = Tree(self)
    self.dataWidgets1.append(self.frame1)
    
    self.dataWidgets2 = []

    self.label2 = Label(self)
    self.dataWidgets2.append(self.label2)
    self.frame2 = Tree(self)
    self.dataWidgets2.append(self.frame2)

    self.drawFrame()


  def drawFrame(self):

    self.repository1 = self.basePopup.repList.currentRepository
    self.repository2 = self.basePopup.repList.currentRepository

    self.repLabel1.grid(row=1, column=1, sticky='n')
    self.repLabel2.grid(row=1, column=3, sticky='n')

    if self.repository1 is None:

      for widget in self.dataWidgets1:
        widget.grid_remove()
        
      self.noRep1Label.grid(row=2, column=1, sticky='n')
      
    else :

      self.noRep1Label.grid_remove()

      self.label1.set(self.repository1.user + '@' + self.repository1.name + '  ( ' + self.repository1.connect + ' )')
      self.label1.grid(row=2,column=1, sticky='w')
      self.frame1.grid(row=3,column=1, sticky='nsew')
      self.update_tree(self.frame1)

    if self.repository2 is None:

      for widget in self.dataWidgets2:
        widget.grid_remove()
        
      self.noRep2Label.grid(row=2, column=3, sticky='n')

    else:

      self.noRep2Label.grid_remove()
      
      self.label2.set(self.repository2.user + '@' + self.repository2.name + '  ( ' + self.repository2.connect + ' )')
      self.label2.grid(row=2,column=3, sticky='w')
      self.frame2.grid(row=3,column=3, sticky='nsew')
      self.update_tree(self.frame2)


  def update_tree(self, tree):

    texts=[]
    icons=[]
    parents=[]
    callbacks=[]
    objects=[]

    # OK, this would be a good point to access the list

    # this should be associated with the repository object
    loc = SharedBeanServiceLocator()
    port = loc.getSharedBean()
    request = getList();

    # these are actually static
    request._arg0 = 'org.pimslims.applet.server.ProjectBean';
    request._arg1 = 'getList';
    request._arg2 = ''

    # get the response
    response = port.getList(request)

    # this is a hack at present. It needs to be written properly
    wsstr = WSString(response._return)
    ss = wsstr.getStruct()

    for strg in ss:
      texts.append(strg)
      icons.append('list-add')
      parents.append(None)
      objects.append(strg)
      callbacks.append(None)

    print 'UPDATE ', ss
    print 'UPDATE ', len(ss)

    if len(ss) > 0:
      print 'UPDATE: updating '
      tree.update(parents, objects, texts, icons, callbacks)
    


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


if __name__ == "__main__":

  import sys
  import Tkinter

  root = Tkinter.Tk()
  root.withdraw()
   
  popup = DataExchMasterPopup(root)
  
  root.mainloop()
