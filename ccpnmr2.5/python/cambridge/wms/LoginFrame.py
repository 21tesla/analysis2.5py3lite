import os

# added jmci; do we need this??
import Tkinter
from memops.universal.Io import getTopDirectory

#from ccpnmr.analysis.core.ExperimentBasic import getOnebondExpDimRefs
from ccpnmr.analysis.popups.BasePopup import BasePopup

from memops.editor.Util import createDismissHelpButtonList

from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList
from memops.gui.Text import Text
from memops.gui.Frame import Frame

from cambridge.wms.RepositoryList import RepositoryList
from cambridge.wms.Repository import Repository

# # # # #   B I G   T O - D O  P O I N T S   # # # # # 
#


# # # # #   L A T E R   T O - D O  P O I N T S    # # # # #
#
 
class LoginPopup(BasePopup):

  def __init__(self, parent):

    self.parent      = parent
    
    BasePopup.__init__(self, parent=parent, title='WMS Login')

                       

  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    frame = LoginFrame(guiFrame, basePopup=self)
    frame.grid(row=0, column=0, sticky='nsew')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    
    self.geometry('280x120')
    
    self.update_idletasks()

 
class LoginFrame(Frame):

  def __init__(self, guiParent, basePopup):

    # Base popup required to handle notification of data model changes
    # e.g. new peak lists, so that the GUI can update to the latest
    # state
    self.basePopup = basePopup
    self.guiParent = guiParent

    self.registerNotify=basePopup.registerNotify
    self.unregisterNotify=basePopup.unregisterNotify


    Frame.__init__(self, guiParent)
  
    # set up the grid

    self.grid_columnconfigure(0, weight=1, minsize=10)
    self.grid_columnconfigure(1, weight=0, minsize=10)
    self.grid_columnconfigure(2, weight=0, minsize=20)
    self.grid_columnconfigure(3, weight=1, minsize=10)

    self.grid_rowconfigure(0, weight=1, minsize=5)
    self.grid_rowconfigure(1, weight=0, minsize=10)
    self.grid_rowconfigure(2, weight=0, minsize=10)
    self.grid_rowconfigure(3, weight=0, minsize=10)
    self.grid_rowconfigure(4, weight=1, minsize=5)

    # build up the body.

    # Column headers

    self.user_label = Label(self,text='Username:')
    self.user_label.grid(row=1, column=1,  padx=5, pady=5, sticky='w')

    self.user_value = Text(self, width=20, height=1, text="")
    self.user_value.grid(row=1, column=2,  padx=5, pady=5, sticky='w')


    self.pswd_label = Label(self,text='Password:')
    self.pswd_label.grid(row=2, column=1,  padx=5, pady=5, sticky='w')

    self.pswd_value = Text(self, width=20, height=1, text="")
    self.pswd_value.grid(row=2, column=2,  padx=5, pady=5, sticky='w')

    
    self.cancel_button = Button(self, width=10, height=1,
                               text="Cancel",
                               command=self.quit )
    self.cancel_button.grid(row=3, column=1,  padx=5, pady=5, sticky='e')

    
    self.login_botton = Button(self, width=10, height=1,
                               text="Login",
                               command=self.login )
    self.login_botton.grid(row=3, column=2,  padx=5, pady=5, sticky='w')


  def login(self):

    username = self.user_value.getText()

    # obviously need an authentication check here!

    self.basePopup.parent.setUsername(username)
    self.basePopup.parent.readConfig(username)
    self.basePopup.destroy()

    
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
   
  popup = LoginPopup(root)

  

 
  root.mainloop()
