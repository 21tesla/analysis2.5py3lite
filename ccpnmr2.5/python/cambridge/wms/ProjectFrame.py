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
from memops.gui.Button import Button
from memops.gui.Text import Text
from memops.gui.ButtonList import ButtonList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.Canvas import Canvas
from memops.gui.Separator import Separator
from memops.gui.Frame import Frame
from memops.gui.Tree import Tree
from memops.gui.Tree import Node
from memops.gui.PulldownList import PulldownList

from cambridge.wms.RepositoryList import RepositoryList
from cambridge.wms.Repository import Repository
from cambridge.wms.Project import Project
from cambridge.wms.FilterFrame import FilterFrame

# # # # #   B I G   T O - D O  P O I N T S   # # # # # 
#


# # # # #   L A T E R   T O - D O  P O I N T S    # # # # #
#


# # # # # # # # # # # ProjectFrame # # # # # # # # # # # 
#

# Widget that provides summary information about a single Project
# (well, strictly project version) and allows project-based
# methods to be fired off
 
class ProjectPopup(BasePopup):

  def __init__(self, parent):

    self.parent      = parent
    
    BasePopup.__init__(self, parent=parent, title='Repository')

                       

  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    frame = ProjectFrame(guiFrame, basePopup=self)
    frame.grid(row=0, column=0, sticky='nsew')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    
    self.geometry('200x80')
    
    self.update_idletasks()

 
class ProjectFrame(Frame):

  def __init__(self, guiParent, basePopup):

    self.basePopup = basePopup
    self.guiParent = guiParent

    self.basePopup.frameShortcuts['Project'] = self

    #self.registerNotify=basePopup.registerNotify
    #self.unregisterNotify=basePopup.unregisterNotify

    Frame.__init__(self, guiParent)
  
    # set up the grid

    self.grid_columnconfigure(0, weight=0, minsize=10)
    self.grid_columnconfigure(1, weight=0, minsize=50)
    self.grid_columnconfigure(2, weight=1, minsize=20)
    self.grid_columnconfigure(3, weight=0, minsize=40)


    self.grid_rowconfigure(0, weight=0, minsize=5)
    self.grid_rowconfigure(1, weight=0, minsize=10)
    self.grid_rowconfigure(2, weight=0, minsize=10)
    self.grid_rowconfigure(3, weight=0, minsize=10)

    # we could in theory choose to store the last visited
    # project in the user config file and load this on log in

    # if we then decide that login is mandatory then this
    # would give as a project from the start
    
    self.repository = None
    self.projectName = None
    self.versionTag = None


    # we define as much as possible once on startup and use the
    # grid manager to control the display thereafter

    # when no data
    
    self.noData_widgets = []
    
    self.null_proj_label = Label(self, text='No Project selected')
    self.noData_widgets.append(self.null_proj_label)

    # when data

    self.data_widgets = []
    
    self.proj_rep_title = Label(self, text='Project:', font='Helvetica16')
    self.data_widgets.append(self.proj_rep_title)

    self.proj_rep_label = Label(self, text='Repository:')
    self.data_widgets.append(self.proj_rep_label)
    self.proj_rep_value = Label(self, text='')
    self.data_widgets.append(self.proj_rep_value)
        
    # This is more complex than it seems. We need to have some
    # user information accoisated with this too (who has the lock)
    self.proj_status_label = Label(self, text = 'Status:')
    self.data_widgets.append(self.proj_status_label)
    self.proj_status_value = Label(self, text='')
    self.data_widgets.append(self.proj_status_value)

    self.proj_details_label = Label(self, text='Description:')
    self.proj_details = Text(self, width=75, height=6, text='')
    self.data_widgets.append(self.proj_details)

    info_modes = ['versions','tasks','raw data']
    self.info_mode_select = PulldownList(self, self.tmpCall, info_modes)
    self.data_widgets.append(self.info_mode_select)
    
    initial_cols = ['Project','Version','Created by Task', 'Created from Project version']
    self.proj_history_matrix = ScrolledMatrix(self, headingList=initial_cols, initialRows=7,
                                      doubleCallback=self.history_callback)
    self.data_widgets.append(self.proj_history_matrix)

    initial_cols = ['ID','Task','User','Status', 'Date']
    self.proj_tasks_matrix = ScrolledMatrix(self, headingList=initial_cols, initialRows=7,
                                            doubleCallback=self.history_callback)
    self.data_widgets.append(self.proj_tasks_matrix)

    initial_cols = ['name','type','size','details']
    self.proj_files_matrix = ScrolledMatrix(self, headingList=initial_cols, initialRows=7,
                                            doubleCallback=self.history_callback)
    self.data_widgets.append(self.proj_files_matrix)

    
    self.proj_export_button = Button(self, 
                                     height=1, width=14,text='Export',
                                     command=self.export_project)
    self.data_widgets.append(self.proj_export_button)

    self.proj_import_button = Button(self, 
                                     height=1, width=14,text='Import',
                                     command=self.export_project)
    self.data_widgets.append(self.proj_import_button)
   
    self.run_testtask1_button = Button(self, 
                                       height=1, width=14,text='Run TestTask1',
                                       command=self.goto_task1_tab)
    self.data_widgets.append(self.run_testtask1_button)

    # FIXME JMCI

    # insert this inside this Frame directly: not worth abstracting out

    self.filter_frame = FilterFrame(self, self.basePopup, text='Filter')
    self.data_widgets.append(self.filter_frame)


    self.drawFrame()



  # Refreshing of the Frame is meaningless without setting the full
  # connection parameters


  def drawFrame(self):

    if self.basePopup.repList:
      self.repository = self.basePopup.repList.currentRepository
      if self.repository:      
        self.projectName = self.repository.currentProjectName
        self.versionTag = self.repository.currentVersionTag
    
    print 'in drawFrame ', self.repository, self.projectName, self.versionTag

    
    if self.projectName == None:

      for widget in self.data_widgets:
        widget.grid_remove()
      
      self.null_proj_label.grid(row=1, column=1, columnspan=4, sticky='nswe')
     
    else:

    # build up the body.

      for widget in self.noData_widgets:
        widget.grid_remove()

      # OK, so we need to determine the project that we are showing
      # in this frame. To uniquely describe a project we need to
      # know the repository, the project name, and the project version
      # so all these need to be set somehow. Maybe the best approach is
      # so have a public method that takes these arguments and to force
      # every refresh to operate through this. Then the drawFrame() method
      # would be private.

      loc = SharedBeanServiceLocator()
      port = loc.getSharedBean()

      # need two calls, one for the Project and one for the version
      hm1 = {'name': self.projectName }
      wsstr_in1 = WSString(hm1)

      request1 = getFields();
      request1._arg0 = 'org.pimslims.applet.server.ProjectBean';
      request1._arg1 = 'getFields';
      request1._arg2 = wsstr_in1.str

      response1 = port.getFields(request1)

      wsstr_out1 = WSString(response1._return)
      hm_proj = wsstr_out1.getStructHM()

      # now the version
  
      hm2 = {'project' : { 'name': self.projectName },
            'versionTag': self.versionTag }

      wsstr_in2 = WSString(hm2)

      request2 = getFields();
      request2._arg0 = 'org.pimslims.applet.server.ProjectVersionBean';
      request2._arg1 = 'getFields';
      request2._arg2 = wsstr_in2.str

      response2 = port.getFields(request2)

      wsstr_out2 = WSString(response2._return)
      hm_pv = wsstr_out2.getStructHM()


      self.proj_rep_title.set('Project: ' + self.projectName + ' v ' + self.versionTag)
      self.proj_rep_title.grid(row=1, column=1, sticky='w')

      self.proj_rep_label.grid(row=2, column=1, sticky='w')
      self.proj_rep_value.set(self.repository.name)
      self.proj_rep_value.grid(row=2, column=2, sticky='w')

      status = 'N/A'
      if hm_pv['status'] != 'valIsNullPointer':
        status = hm_pv['status']
        
      self.proj_status_label.grid(row=3, column=1, sticky='w')
      self.proj_status_value.set(status)
      self.proj_status_value.grid(row=3, column=2, sticky='w')

      # also want to know:
      # * how this file was created
      # * what external file (e.g. spectrum files) it uses
      # * the log history
      # * what processes it has been used in
      # * if locked, what process, user is responsible for this

      details = 'N/A'
      if hm_proj['details'] != 'valIsNullPointer':
        details = hm_proj['details']
      
      self.proj_details_label.grid(row=4, column=1, columnspan=2, sticky='nw')
      self.proj_details.setText(details)
      self.proj_details.grid(row=4, column=2, pady=10, sticky='nw')

      self.info_mode_select.grid(row=5, column=1, sticky='w')

      # this needs to go into a separate subroutine

      # get al the version
  
      hm3 = {'project' : { 'name': self.projectName }}
      wsstr_in3 = WSString(hm3)

      request3 = getList()
      request3._arg0 = 'org.pimslims.applet.server.ProjectVersionBean';
      request3._arg1 = 'getListWithFields';
      request3._arg2 =  wsstr_in3.str

      response3 = port.getList(request3)

      wsstr_out3 = WSString(response3._return)
      vers = wsstr_out3.getStruct()

      matrix = [];
      objs = [];

      print 'PROJECT: ', vers

      for ver in vers:

        taskSerial = 'unknown'
        if (ver.has_key('createdByTask')):
          task = ver['createdByTask']
          taskSerial = task['serial']

        # FIXME JMCI

        # Really should get this from the task, of it exists

        parentVersionTag = 'unknown'
        if (ver.has_key('sourceVersion')):
          pv = ver['sourceVersion']
          pvp = pv['project']
          parentName = pvp['name']
          parentVersionTag = pvp['name'] + ' ' + pv['versionTag']

        new_row = (self.projectName, ver['versionTag'], taskSerial, parentVersionTag )
        matrix.append(new_row)
        objs.append(ver)

      self.proj_history_matrix.update(objectList=objs, textMatrix=matrix)
      self.proj_history_matrix.grid(row=6, column=1, columnspan=2, pady=10, sticky='nsew')

      self.filter_frame.grid(row=6, column=3, pady=10, sticky='n')
      

      # have a button to export. It should now be clear which version
      # of which project in which repository should be exported

      self.proj_import_button.grid(row=1, column=3, padx=3,pady=3, sticky='ew')
      self.proj_export_button.grid(row=2, column=3, padx=3,pady=3, sticky='ew')
        
      self.run_testtask1_button.grid(row=4, column=3, padx=3,pady=3, sticky='ew')

      # It is not clear precisely how methods should be triggered from
      # the project frame, as there may be >1 input to a given task with
      # a specific role. It maybe that we need to parse out a metadescription
      # of the task and use this to produce the correct buttons (maybe
      # even depending on the contents of the file ?)

      #self.task_button_texts=['Test1Task','Test2Task']
      #self.task_button_cmds=[self.runTest1Task, self.runTest2Task]
      #self.task_button_list = ButtonList(self, self.task_button_texts, self.task_button_cmds,
      #                                    None, 1, 1, True, None, Tkinter.VERTICAL)
      # 
      #self.task_button_list.grid(row=1, rowspan=3, column=4, padx=3, sticky='w')

  def history_callback (self, obj, row, col):

    if (col == 1):
      pv = obj['versionTag']
      print 'PROJECT -> PROJECT (DIRECT) ', pv
      self.repository.currentVersionTag = pv.__str__()

      self.drawFrame()
      
    elif (col == 2):
      task = obj['createdByTask']
      taskSerial = task['serial']
      
      print 'PROJECT -> TASK ', taskSerial

      self.basePopup.currentTask = taskSerial.__str__()    

      # FIXME

      # This is more complex. Need to navigate to the correct
      # panel within the Task submenu

      if self.basePopup.frameShortcuts.has_key('Task'):
        taskFrame = self.basePopup.frameShortcuts['Task']
        # really need to check on the type of task
        if taskFrame.frameShortcuts.has_key('Test1'):
          taskFrame.frameShortcuts['Test1'].drawFrame()
      # also need to select the current tabbed frame
      self.basePopup.tabbedFrame.select(2)
      

    # FIXME JMCI

    # need to look into this. Really we need both the project and
    # the version for this (in theory this could be based on a different
    # project although it is highly unlikely!)

    elif (col == 3):
      pv = obj['sourceVersion']
      parentVersion = pv['versionTag']
      
      pvp = pv['project']
      parentName = pvp['name']
      
      print 'PROJECT -> PROJECT (SOURCE) ', parentVersion

      self.repository.currentVersionTag = parentVersion.__str__()
      self.repository.currentProjectName = parentName.__str__()

      self.drawFrame()
 

  # create a new task with this project as an input

  # it is not clear how we should do this:- construct data locally
  # and then save, or save a first cut and then ammend. On the basis
  # that the serial has to be generated by the DB, we will go for the
  # former for now

  # it would be better to pass this to a processor layer that
  # determines parameters and then builds the stub.

  def goto_task1_tab(self):

    # create a new task (provisional status)

    loc = SharedBeanServiceLocator()
    port = loc.getSharedBean()

    #hm0 = { 'project': { 'name': self.projectName},
    #        'versionTag' : self.versionTag,
    #        '_handle' : 'findOrFail' }

    # FIXME

    # We have an error here (it seems to be a problem between the CCPN
    # and hibernate layers that is not picking objects as the same for
    # some reason. Therefore do not set the inputVersion for now

    hm0 = { 'project': { 'name': self.projectName},
            'versionTag' : self.versionTag}
    hm2 = { 'inputVersion' : hm0,
            'status' : 'PROVISIONAL' }

    #hm2 = { 'status' : 'PROVISIONAL' }

    wsstr_in1 = WSString(hm2)
            
    req1 = recordAndLoad();

    req1._arg0 = 'org.pimslims.applet.server.TaskBean';
    req1._arg1 = 'recordTest1';
    req1._arg2 =  wsstr_in1.str
    
    resp1 = port.recordAndLoad(req1)
    wsstr_out1 = WSString(resp1._return)
    v = wsstr_out1.getStructHM()

    self.basePopup.currentTask =   v['serial'].__str__()

    # FIXME

    # This is more complex. Need to navigate to the correct
    # panel within the Task submenu

    if self.basePopup.frameShortcuts.has_key('Task'):
      taskFrame = self.basePopup.frameShortcuts['Task']
      # really need to check on the type of task
      if taskFrame.frameShortcuts.has_key('Test1'):
        taskFrame.frameShortcuts['Test1'].drawFrame()
    # also need to select the current tabbed frame
    self.basePopup.tabbedFrame.select(2)

    
  # method for exporting repository directory to local file
  # system. This is a very basic approach, essentially first
  # listing all the files required for transfer and then
  # copying them over one by one. Note that due to limitations
  # in ZSI we cannot do much better than base64 encoding each
  # file. We could consider storing a tgz on the server and
  # exporting that. It might be a considerable improvement as
  # the compression leads to a file < 10% of the original size
  # which would cut down on traffic enormously

  def export_project(self):

    versionTag = self.versionTag
    name = self.projectName

    # need to preserve the current status. Take a default from the
    # Wms layer 
    rootDir = self.basePopup.repList.current_export_dir

    # This should be obtained from a query box. Hard code for now
    fileSelectPopup = FileSelectPopup(self, None, rootDir)

    exp_project_dir = fileSelectPopup.getDirectory()
    self.basePopup.repList.current_export_dir = exp_project_dir

    # These should come from the repository object
    loc = SharedBeanServiceLocator()
    port = loc.getSharedBean()

    req1 = getList()

    hm1 = {'projectName': name,
           'versionTag': versionTag,
           'segment': 'jionides' }
    wsstr1 = WSString(hm1)

    # these are actually static
    req1._arg0 = 'org.pimslims.applet.server.CcpnFileBean';
    req1._arg1 = 'getList';
    req1._arg2 = wsstr1.str

    # get the response
    resp1 = port.getList(req1)

    # this is a hack at present. It needs to be written properly
    wsstr = WSString(resp1._return)
    ss = wsstr.getStruct()

    print 'got array ', ss

    for strg in ss:

      print 'trying to handle file ', strg

      hm2 = {'projectName': name,
             'versionTag': versionTag,
             'segment': 'jionides',
             'fileName' : '/' + strg.__str__(),
             'random' : 'test',
             'raw' : False}
      
      print 'trying to handle file 2 ', strg
      print 'trying to handle file 2 ', hm2

      req2 = getDocWithParams()
      
      wsstr2 = WSString(hm2)

      # these are actually static
      req2._arg0 = 'org.pimslims.applet.server.CcpnFileBean';
      req2._arg1 = wsstr2.str

      # get the response
      resp2 = port.getDocWithParams(req2)

      local_file_path = exp_project_dir + '/' + name + '/' + strg

      idx = local_file_path.rfind('/')
      local_dir = local_file_path[:idx]
      print 'trying to create ', local_dir
      os.makedirs(local_dir)

      print 'writing to local area; ', local_file_path

      f = open(local_file_path, 'w')
      f.write(resp2._return)
          

  # it is likely that these methods are going to have to be highly specialised
  # as some tasks are going to have many input files and we need to know which
  # role the project we are in takes in the task.

  def runTest1Task(self):

    pass

  def runTest2Task(self):

    pass

  def view_project(self, project):

    self.project = project

    self.drawFrame()


  def tmpCall(self):

    return


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
   
  popup = ProjectPopup(root)
 
  root.mainloop()
