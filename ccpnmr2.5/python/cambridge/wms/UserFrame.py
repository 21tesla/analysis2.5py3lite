
# Essentially the same as the ExtendNMR GUI for now. 


import os

# added jmci; do we need this??
import Tkinter
from memops.universal.Io import getTopDirectory

#from ccpnmr.analysis.core.ExperimentBasic import getOnebondExpDimRefs
from ccpnmr.analysis.popups.BasePopup import BasePopup

from memops.editor.Util import createDismissHelpButtonList

from memops.gui.Frame import Frame
from memops.gui.Label import Label


# # # # #   B I G   T O - D O  P O I N T S   # # # # # 
#


# # # # #   L A T E R   T O - D O  P O I N T S    # # # # #
#
 
class UserPopup(BasePopup):

  def __init__(self, parent):

    self.parent      = parent
    
    BasePopup.__init__(self, parent=parent, title='PRODECOMP')

                       
  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    frame = UserFrame(guiFrame, basePopup=self, ccpnProject=self.ccpnProject)
    frame.grid(row=0, column=0, sticky='nsew')
    
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    
    self.geometry('650x600')
    
    self.update_idletasks()
 

class UserFrame(Frame):

  def __init__(self, guiParent, basePopup):
    # Base popup required to handle notification of data model changes
    # e.g. new peak lists, so that the GUI can update to the latest
    # state

    print 'setting up userFrame'
    
    self.basePopup = basePopup
    self.guiParent = guiParent

    self.registerNotify=basePopup.registerNotify
    self.unregisterNotify=basePopup.unregisterNotify

    Frame.__init__(self, guiParent)
  
    self.grid_columnconfigure(0, weight=0, minsize=10)
    self.grid_columnconfigure(1, weight=1)
    self.grid_columnconfigure(2, weight=1)
    self.grid_columnconfigure(3, weight=1)
    self.grid_columnconfigure(4, weight=1)

    self.grid_rowconfigure(0, weight=0, minsize=10 )
    self.grid_rowconfigure(1, weight=0, minsize=10 )
    self.grid_rowconfigure(2, weight=0, minsize=10 )

    self.draw()

  def draw(self):

    # there must be a better way of doing this. How does the garbage
    # collector work ... do I have to formally remove widgets that
    # I no longer intend to use or it is better to adopt a strategy
    # of recyling them?

    for child in self.children.values():
      child.grid_remove()

    #
    # not clear whether this should come from an object under repList
    # or from a simple variable in root. Probably best that there is
    # a data struture that includes user and repList that is initialised
    # on Login (or built from the readConfig method)

    #
    # really need to set up a dummy structure on login that we can build
    # from if the user does not wish to log in
    
    user_label = Label(self, text="Currently Logged in as User: " + self.basePopup.repList.user.name)
    user_label.grid(row=1, column=1, columnspan=2, sticky='w')


    #
    # also need to have the settings and repositories

    rep_label = Label(self, text="Active repositories: ")
    rep_label.grid(row=2, column=1, columnspan=2, pady=10, sticky='w')

    # now rebuild. There must be a better way of doing this?
    repRow = 3

    # should maybe access through method?
    for rep in self.basePopup.repList.getRepositories():

      self.grid_rowconfigure(repRow, weight=0, minsize=10 )

      rep_name = Label(self, text=rep.name )
      rep_name.grid(row=repRow, column=1,  sticky='w')

      rep_conn = Label(self, text=rep.connect )
      rep_conn.grid(row=repRow, column=2,  sticky='w')

      if rep.type == None:
        rep_type_text = 'not specified'
      else:
        rep_type_text = rep.type
      
      rep_type = Label(self, text=rep_type_text )
      rep_type.grid(row=repRow, column=3,  sticky='w')

      repRow += 1    

    self.grid_rowconfigure(repRow, weight=1, minsize=10 )

    


    
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
   

  popup = UserPopup(root)
  
  root.mainloop()
