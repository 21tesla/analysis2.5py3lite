import os

# required for WS layer
from SharedBeanService_services import *
from WSString import *

# added jmci; do we need this??
import Tkinter
from memops.universal.Io import getTopDirectory

#from ccpnmr.analysis.core.ExperimentBasic import getOnebondExpDimRefs
from ccpnmr.analysis.popups.BasePopup import BasePopup

from memops.editor.Util import createDismissHelpButtonList

from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.Text import Text
from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList
from memops.gui.Canvas import Canvas
from memops.gui.Separator import Separator
from memops.gui.Frame import Frame
from memops.gui.Tree import Tree
from memops.gui.Tree import Node

from cambridge.wms.RepositoryList import RepositoryList
from cambridge.wms.Repository import Repository

# # # # #   B I G   T O - D O  P O I N T S   # # # # # 
#


# # # # #   L A T E R   T O - D O  P O I N T S    # # # # #
#
 
class FilterPopup(BasePopup):

  def __init__(self, parent):

    self.parent      = parent
    
    BasePopup.__init__(self, parent=parent, title='Repository')

                       

  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    frame = FilterFrame(guiFrame, basePopup=self)
    frame.grid(row=0, column=0, sticky='nsew')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    
    self.geometry('160x80')
    
    self.update_idletasks()
 
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


    # get a port proxy instance
    # this should probably belong to the repository directly
    
    loc = SharedBeanServiceLocator()
    self.port = loc.getSharedBean()


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

    self.nameText = Text(self, width=20, height=1)
    self.nameText.grid(row=1, column=1, sticky='nw')

    self.dateLabel = Label(self,text='Date:')
    self.dateLabel.grid(row=2, column=0, padx=2, sticky='nw')

    self.dateText = Text(self, width=20, height=1)
    self.dateText.grid(row=2, column=1, sticky='nw')

    self.userLabel = Label(self,text='User:')
    self.userLabel.grid(row=3, column=0, padx=2, sticky='nw')

    self.userText = Text(self, width=20, height=1)
    self.userText.grid(row=3, column=1, sticky='nw')

    self.clearButton = Button(self, text='Clear', command=self.clear_filter, width=10, height=1)
    self.clearButton.grid(row=4, column=0, columnspan=2, sticky='nw')

    self.filterButton = Button(self, text='Go', command=self.tmpCall, width=10, height=1)
    self.filterButton.grid(row=4, column=1, columnspan=2, sticky='ne')


  def clear_filter(self):

    self.nameText.xlear()
    self.dateText.clear()
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
   
  popup = FilterPopup(root)

  
 
  root.mainloop()
