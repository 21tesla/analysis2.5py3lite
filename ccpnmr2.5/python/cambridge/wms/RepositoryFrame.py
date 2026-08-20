import os, re, time

# required for WS layer
from SharedBeanService_services import *
from WSString import *

# added jmci; do we need this??
import Tkinter
from memops.universal.Io import getTopDirectory

#from ccpnmr.analysis.core.ExperimentBasic import getOnebondExpDimRefs
from ccpnmr.analysis.popups.BasePopup import BasePopup

from memops.editor.Util import createDismissHelpButtonList

from memops.gui.FileSelectPopup import FileSelectPopup

from memops.gui.Label import Label
from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList
from memops.gui.Canvas import Canvas
from memops.gui.Separator import Separator
from memops.gui.Frame import Frame
from memops.gui.Tree import Tree
from memops.gui.Tree import Node
from memops.gui.ScrolledListbox import ScrolledListbox
from memops.gui.PulldownList import PulldownList
from memops.gui.RadioButtons import RadioButtons

from cambridge.wms.RepositoryList import RepositoryList
from cambridge.wms.Repository import Repository
from cambridge.wms.FilterFrame import FilterFrame

# # # # #   B I G   T O - D O  P O I N T S   # # # # # 
#


# # # # #   L A T E R   T O - D O  P O I N T S    # # # # #
#

 
class RepositoryPopup(BasePopup):

  def __init__(self, parent):

    self.parent      = parent
    
    BasePopup.__init__(self, parent=parent, title='Repository')

                       

  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    frame = RepositoryFrame(guiFrame, basePopup=self)
    frame.grid(row=0, column=0, sticky='nsew')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    
    self.geometry('200x40')
    
    self.update_idletasks()
 
class RepositoryFrame(Frame):

  def __init__(self, guiParent, basePopup):

    # Base popup required to handle notification of data model changes
    # e.g. new peak lists, so that the GUI can update to the latest
    # state
    self.basePopup = basePopup
    self.guiParent = guiParent

    # should the screen autorefresh
    self.autoRefresh = False

    # add this to shortcuts to ease navigation    
    self.basePopup.frameShortcuts['Repository'] = self

    # get a port proxy instance
    # this should probably belong to the repository directly
    # or else should belong in a dictionary in the main gui
    # layer when it can be picked up easily

    # FIXME JMCI
    # need to work out how to get a number of these!
    
    loc = SharedBeanServiceLocator()
    self.port = loc.getSharedBean()


    self.registerNotify=basePopup.registerNotify
    self.unregisterNotify=basePopup.unregisterNotify


    Frame.__init__(self, guiParent)
  
    # set up the grid

    self.grid_columnconfigure(0, weight=0, minsize=20)
    self.grid_columnconfigure(1, weight=1, minsize=10)
    self.grid_columnconfigure(2, weight=0, minsize=20)

    self.grid_rowconfigure(0, weight=0, minsize=5)
    self.grid_rowconfigure(1, weight=0, minsize=0)
    self.grid_rowconfigure(2, weight=0, minsize=10)
    self.grid_rowconfigure(3, weight=0, minsize=10)
    self.grid_rowconfigure(4, weight=0, minsize=10)
    self.grid_rowconfigure(5, weight=1, minsize=10)
    self.grid_rowconfigure(6, weight=0, minsize=10)
    self.grid_rowconfigure(7, weight=0, minsize=10)

    # widgets for view when no data

    self.noDataWidgets = []

    self.noRepLabel = Label(self,text='No repository currently selected.')
    self.noDataWidgets.append(self.noRepLabel)

    # widgets for view when current repository set

    self.dataWidgets = []

    self.repTitle = Label(self,text='Repository:', font='Helvetica16')
    self.dataWidgets.append(self.repTitle)

    #self.repLabel = Label(self,text='All Projects in Repository')
    #self.dataWidgets.append(self.repLabel)
    
    self.repTree = Tree(self, doubleCallback=self.goto_project_tab)
    self.dataWidgets.append(self.repTree)

    sel = 1
    if self.autoRefresh:
      sel = 0
    self.autoRefreshSwitchLabel = Label(self, text= 'Auto Refresh')
    self.autoRefreshSwitch = RadioButtons(self, ['on', 'off'], select_callback=self.set_refresh, selected_index = sel)

    # need to decide whether this is static or not
    r_button_texts=['Add to Basket','  Import  ','  Export  ','  Refresh ','Properties']
    r_button_cmds=[self.tmpCall,self.import_project, self.export_project, self.drawFrame,self.goto_project_tab]

    self.rep_button_list = ButtonList(self, r_button_texts, r_button_cmds)
    self.dataWidgets.append(self.rep_button_list)

    self.filterFrame = FilterFrame(self, self.basePopup, text='Filter')
    self.dataWidgets.append(self.filterFrame)

    baskets=('basket1','basket2','basket3')
    self.basketSelect = PulldownList(self, self.tmpCall, baskets)
    self.dataWidgets.append(self.basketSelect)

    basketElements=('1ay3','1ay7')
    self.basketList = ScrolledListbox(self, basketElements, 25, 18)
    self.dataWidgets.append(self.basketList)

    b_button_texts=['Remove from Basket','New Basket']
    b_button_cmds=[self.tmpCall, self.tmpCall]

    self.b_button_list = ButtonList(self, b_button_texts, b_button_cmds)
    self.dataWidgets.append(self.b_button_list)

    # draw if not automatically refreshing
    if sel == 1:
      self.drawFrame()

    # set up a loop
    self.refresh()


  def drawFrame(self):

    if self.basePopup.repList.currentRepository == None:

      for widget in self.dataWidgets:
        widget.grid_remove()
        
      self.noRepLabel.grid(row=1, column=1, sticky='n')
      
    else :

      self.repository = self.basePopup.repList.currentRepository
      
      for widget in self.noDataWidgets:
        widget.grid_remove()
        
      # Column headers
      # want to have a general admin panel; current repository, button to change, info
      # about login and so on
      
      self.repTitle.set('Repository: ' + self.repository.user + '@' + self.repository.name + '  ( ' + self.repository.connect + ' )')

      self.repTitle.grid(row=1, column=1, sticky='w')

      # self.repLabel.grid(row=2, column=1, sticky='w')

      # Repository block

      #grid and update
      self.repTree.grid(row=3, column=1, rowspan=3, sticky='nsew')
      self.updateRepTree()

      self.autoRefreshSwitchLabel.grid(row=6, column=1, padx=10, sticky='w')
      self.autoRefreshSwitch.grid(row=6, column=1, padx=10, sticky='e')

      # the buttons for functionality in repository area
      self.rep_button_list.grid(row=7, column=1, padx=3, sticky='ew')

      # filter frame
      self.filterFrame.grid(row=3, column=2, stick='nsew')

      # list of baskest. Will callback to change current basket
      self.basketSelect.grid(row=4, column=2, padx=10, stick='w')

      # lists projects per basket
      self.basketList.grid(row=5, rowspan=2, column=2, padx=10, stick='nsew')

      self.b_button_list.grid(row=7, column=2, padx=5, stick='ew')


  def openLink(self, node):

    tree = node.tree
    par = node.object

    # this should be associated with the repository object
    request = getList();

    # FIXME JMCI

    # 1. Potentially need to traverse up the node tree
    # 2. Currently there is bug in expansion/contraction
    # 3. Search needs to work off proxy method (not off
    #    generated stubs
    
    # pass the parent name for now. may need other criteria later
    h1 = {'name':  par.__str__() }
    h2 = {'project': h1}
    wsstr_in = WSString(h2)

    request._arg0 = 'org.pimslims.applet.server.ProjectVersionBean';
    request._arg1 = 'getListWithFields';
    request._arg2 = wsstr_in.str

    # get the response
    response = self.port.getList(request)

    # currently we just get a list of numbers. May want to send
    # all the fields for every object too to cut down on db calls
    wsstr_out = WSString(response._return)
    ss = wsstr_out.getStruct()

    for hm in ss:

      print 'handling ', hm, ', ', hm['status'].__str__()
      #if hm['status'].__str__() == 'AVAILABLE':
      icon='text-x-generic'
      #else:
      #  icon='emblem-readonly'
      pv = par.__str__() + '::' + hm['versionTag'].__str__()
      obj = pv
      text= hm['versionTag'].__str__()
      callback=None
      
      tree.add(node, obj, text, icon, callback)

  def set_refresh(self, text):

    print 'setting refresh ', text

    if text == 'on':
      self.autoRefresh = True
    else:
      self.autoRefresh = False


  def refresh(self):

    if self.autoRefresh:
      print 'refreshing frame ', time.time()
      self.drawFrame()
      
    self.after(5000, self.refresh)


  # FIME JMCI 

  # it would probably be a good idea to have a method that allowed
  # a user to unlock a given project version. Would need to check
  # permissions very carefully, though

  def unlock_project_version(self):

    pass



  def goto_project_tab(self, node=None):

    if node is None:
      # get the selected node
      selected_nodes = self.repTree.getSelectedNodes()
      # restrict the choice to one
      node = selected_nodes[0]

    # we really need to traverse the tree a lot more carefully

    pv = node.object.__str__()
    print 'GOTO PROJECT ', pv, ', ', node, ', ', node.object
    pat = re.compile('::')
    matcher = pat.search(pv, 0)
    name = pv[0:matcher.start()]
    versionTag = pv[matcher.end():]

    print 'GOTO PROJECT ', name, ', ', versionTag
    
    
    #versionTag = node.object.__str__()
    #name = node.parent.object.__str__()

    #print node.__dict__
    #print node.parent.__dict__

    self.basePopup.repList.currentRepository.currentProjectName = name
    self.basePopup.repList.currentRepository.currentVersionTag = versionTag

    # this is going to have to be keyed on a natural key that then
    # gets set into the criteria of a ws query

    # hack for now. build in later. Essentially, we have to handle two
    # cases; that of a specific projectVersion request and that of a
    # request on a project alone. Probably the best behaviour would be
    # for the project-only request to default to the latest trunk
    # version of the project and the most elegent way of triggering
    # this would probably be to have this as a property of the project
    # returned in the main hash from Project.getFields().

    self.basePopup.tabbedFrame.select(1)    
    if self.basePopup.frameShortcuts.has_key('Project'):
      self.basePopup.frameShortcuts['Project'].drawFrame()


  def import_project(self):

    # FIXME JMCI

    # need to preserve the current status. Take default from the
    # Wms layer 
    rootDir = self.basePopup.repList.current_import_dir

    fileSelectPopup = FileSelectPopup(self, None, rootDir)

    new_project_dir = fileSelectPopup.getDirectory()
    self.basePopup.repList.current_import_dir = new_project_dir
    
    print 'in import project with directory', new_project_dir

    idx = new_project_dir.rfind('/')
    new_project_name = new_project_dir[idx+1:]
    

    print 'in import project with project', new_project_name

    # FIXME

    # need to set the project version number somewhere. For now set to none
    # and this gets picked up and set to a default 1.1 later

    self.repository.import_project(new_project_name, None, new_project_dir)

    self.drawFrame()


  def export_project(self):

    # get the selected node
    selected_nodes = self.repTree.getSelectedNodes()

    # restrict the choice to one
    node = selected_nodes[0]

    # we really need to traverse the tree a lot more carefully
    name = node.parent.object.__str__()
    versionTag = node.object.__str__()[len(name)+2:]

    

    # need to preserve the current status. Take a default from the
    # Wms layer 
    rootDir = self.basePopup.repList.current_export_dir

    # This should be obtained from a query box. Hard code for now
    fileSelectPopup = FileSelectPopup(self, None, rootDir)

    exp_project_dir = fileSelectPopup.getDirectory()
    self.basePopup.repList.current_export_dir = exp_project_dir

    self.repository.export_project(name, versionTag, exp_project_dir)

    

  def disconnectRepository(self):

    # do we need to formally break a connection and destroy the session?
    # yes we probably should.

    # need an "are you sure" dialog box

    # need to protect
    
    self.repository.repList.repositories.remove(self.repository)
    self.grid_remove()
    for ff in self.guiParent.guiParent.parent.frames[0].children.values():
      ff.drawFrame()

  def tmpCall(self):

    return

  # This will need a separate call to the WS with an additional
  # criterion (picking up the checkout version)

  def updateRepTree(self):

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
      icons.append('folder')
      parents.append(None)
      objects.append(strg)
      callbacks.append(self.openLink)

    print 'UPDATE ', ss
    print 'UPDATE ', len(ss)

    if len(ss) > 0:
      print 'UPDATE: updating '
      self.repTree.update(parents, objects, texts, icons, callbacks)
    

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
   
  popup = RepositoryPopup(root)

  
 
  root.mainloop()
