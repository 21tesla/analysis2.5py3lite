import os, time

# required for WS layer
from SharedBeanService_services import *
from WSString import *

# added jmci; do we need this??
import Tkinter
from memops.universal.Io import getTopDirectory

#from ccpnmr.analysis.core.ExperimentBasic import getOnebondExpDimRefs
from ccpnmr.analysis.popups.BasePopup import BasePopup

from memops.editor.Util import createDismissHelpButtonList

#from memops.gui.ScrolledFrame import ScrolledFrame
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.Text import Text
from memops.gui.Button import Button
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.RadioButtons import RadioButtons


#from cambridge.wms.FilterFrame import FilterFrame

# # # # #   B I G   T O - D O  P O I N T S   # # # # # 
#


# # # # #   L A T E R   T O - D O  P O I N T S    # # # # #
#
 
class TasksPopup(BasePopup):

  def __init__(self, parent):

    self.parent      = parent
    
 
    BasePopup.__init__(self, parent=parent, title='Repository')

                       
  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    frame = TasksFrame(guiFrame, basePopup=self)
    frame.grid(row=0, column=0, sticky='nsew')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    
    self.geometry('200x50')
    
    self.update_idletasks()


class TasksFrame(Frame):

  
  def __init__(self, guiParent, basePopup):

    self.guiParent = guiParent
    self.basePopup = basePopup

    # selection hash. This needs to be stored so that we can
    # refresh with the same criteria
    self.select = None

    # should the screen autorefresh
    self.autoRefresh = False
    
    # add this to shortcuts to ease navigation    
    self.basePopup.frameShortcuts['Tasks'] = self

    #self.registerNotify=basePopup.registerNotify
    #self.unregisterNotify=basePopup.unregisterNotify

    #self.frames = []

    # FIXME JMCI

    # This frame is specific for a single repository so we need a
    # method for identifying one from the whole set (or could have
    # subtabs)

    Frame.__init__(self, guiParent)
  
    # set up the grid

    self.grid_columnconfigure(0, weight=0, minsize=20)
    self.grid_columnconfigure(1, weight=1, minsize=30)
    self.grid_columnconfigure(2, weight=0, minsize=40)

    self.grid_rowconfigure(0, weight=0, minsize=20)
    self.grid_rowconfigure(1, weight=0, minsize=10)
    self.grid_rowconfigure(2, weight=0, minsize=0)
    self.grid_rowconfigure(3, weight=0, minsize=0)
    self.grid_rowconfigure(4, weight=0, minsize=0)
    self.grid_rowconfigure(5, weight=1, minsize=0)
    self.grid_rowconfigure(6, weight=0, minsize=10)
    self.grid_rowconfigure(7, weight=0, minsize=10)


    initial_cols = ['ID','Task','User','Status']
    self.task_matrix = ScrolledMatrix(self, headingList=initial_cols, initialRows=15,
                                      doubleCallback=self.goto_task_tab)

    # Filters for narrowing down tasks. should probably code this
    # specifically inside this package    
    #self.filter = FilterFrame(self, self.basePopup, text='Filter')
    self.filter = TaskFilterFrame(self, self.basePopup, text='Filter')

    # Frame for controlling the main client side daemon
    # seems to be a problem. add a simple button for now
    opts = ['Client','Server']
    self.daemonFrame = TabbedFrame(self, options=opts,
                                   toggleOff=False, selected=0)

    self.daemonSwitchLabel = Label(self, text= 'Task Daemon')
    self.daemonSwitch = RadioButtons(self, ['on', 'off'], select_callback=self.set_daemon, selected_index = 1)

    sel = 1
    if self.autoRefresh:
      sel = 0
    self.autoRefreshSwitchLabel = Label(self, text= 'Auto Refresh')
    self.autoRefreshSwitch = RadioButtons(self, ['on', 'off'], select_callback=self.set_refresh, selected_index = sel)

    self.label_rep = Label(self, text='', font='Helvetica16')

    # set up a loop
    self.refresh()


  def drawFrame(self):

    # build up the body.

    # This is a very simple mock up. In reality we need to be
    # able to select a number of different input files that
    # depend on the task type

    # we need to have the current repository in context
    self.repository = self.basePopup.repList.currentRepository

    if self.repository == None:

      self.label_rep.grid_remove()
      self.task_matrix.grid_remove()
      self.filter.grid_remove()
      self.daemonSwitchLabel.grid_remove()
      self.daemonSwitch.grid_remove()
      
    else:
      
      strg = self.repository.user + '@' + self.repository.name + ' (' + self.repository.connect + ')'

      self.label_rep.set('TaskManager: Repository ' + strg)
      self.label_rep.grid(row=1, column=1, columnspan=2, sticky='w')

      self.task_matrix.grid(row=2, column=1, rowspan=4, sticky='nsew')


      loc = SharedBeanServiceLocator()
      port = loc.getSharedBean()

      request1 = getList();

      request1._arg0 = 'org.pimslims.applet.server.TaskBean';
      request1._arg1 = 'getListWithFields';

      if self.select==None:
        request1._arg2 = ''
      else:
        wsstr_in = WSString(self.select)
        request1._arg2 = wsstr_in.str
        

      response1 = port.getList(request1)

      wsstr1 = WSString(response1._return)
      ss1 = wsstr1.getStruct()

      matrix = [];
      objs = [];

      for el in ss1:
        new_row = (el['serial'], 'TestTask1', 'jionides', el['status'])
        matrix.append(new_row)
        objs.append(el)

      self.task_matrix.update(objectList=objs, textMatrix=matrix)

      self.filter.grid(row=2, column=2, sticky='n')

      # FIXME JMCI

      # The tabbed frame seems to cause the GUI to hang over expansion.
      # use a simple radio button for now

      # self.daemonFrame.grid(row=3, column=2, padx=10, sticky='nsew')
      self.daemonSwitchLabel.grid(row=3, column=2, padx=10, sticky='n')
      self.daemonSwitch.grid(row=4, column=2, padx=10, sticky='n')

      self.autoRefreshSwitchLabel.grid(row=6, column=1, padx=10, sticky='w')
      self.autoRefreshSwitch.grid(row=6, column=1, padx=10, sticky='e')

  def selectTask(self):

    row = self.task_matrix.currentCell
    print 'in select task ', row


  def set_daemon(self, text):

    print 'setting daemon ', text

    if text == 'on':
      self.basePopup.taskDaemon.active = True
    else:
      self.basePopup.taskDaemon.active = False
    
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

  def goto_task_tab(self, obj, row, col):

    print obj
    self.basePopup.currentTask =  obj['serial'].__str__()   

    # This is really awkward because we actually have to go to the
    # correct tab. We have some horrible hacks here. Need a dictionary
    # that can take us straight to the correct frame based on a
    # keyword

    self.basePopup.tabbedFrame.select(2)
    #for ff in self.basePopup.tabbedFrame.frames[2].children.values():
    #  ff.tabbedFrame.select(0)

    if self.basePopup.frameShortcuts.has_key('Task'):
      taskFrame = self.basePopup.frameShortcuts['Task']
      # really need to check on the type of task
      if taskFrame.frameShortcuts.has_key('Test1'):
        taskFrame.frameShortcuts['Test1'].drawFrame()
    # also need to select the current tabbed frame
    self.basePopup.tabbedFrame.select(2)
        

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
  
    #self.administerNotifiers(self.basePopup.unregisterNotify)
    Frame.destroy(self)


class TaskFilterFrame(LabelFrame):

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


    # FIXME
    # should be an enumerated list
    self.statusLabel = Label(self,text='Type:')
    self.statusText = Text(self, width=20, height=1)

    # FIXME
    # should be an enumerated list
    self.typeLabel = Label(self,text='Type:')
    self.typeText = Text(self, width=20, height=1)

    # FIXME
    # should be an enumerated list
    self.projLabel = Label(self,text='Project:')
    self.projText = Text(self, width=20, height=1)

    # FIXME
    # should be an enumerated list
    self.userLabel = Label(self,text='User:')
    self.userText = Text(self, width=20, height=1)

    self.clearButton = Button(self, text='Clear', command=self.clear_filter, width=10, height=1)
    self.filterButton = Button(self, text='Go', command=self.refresh_filter, width=10, height=1)


    self.statusLabel.grid(row=2, column=0, padx=2, sticky='nw')
    self.statusText.grid(row=2, column=1, sticky='nw')
    self.typeLabel.grid(row=3, column=0, padx=2, sticky='nw')
    self.typeText.grid(row=3, column=1, sticky='nw')
    self.projLabel.grid(row=4, column=0, padx=2, sticky='nw')
    self.projText.grid(row=4, column=1, sticky='nw')
    self.userLabel.grid(row=5, column=0, padx=2, sticky='nw')
    self.userText.grid(row=5, column=1, sticky='nw')
    
    self.clearButton.grid(row=6, column=0, columnspan=2, sticky='nw')
    self.filterButton.grid(row=6, column=1, columnspan=2, sticky='ne')

  def clear_filter(self):

    self.statusText.clear()
    self.typeText.clear()
    self.projText.clear()
    self.userText.clear()

  def refresh_filter(self):

    status = self.statusText.getText()[:-1]
    type = self.typeText.getText()[:-1]
    proj = self.projText.getText()[:-1]
    user = self.userText.getText()[:-1]

    print 'status ',status,'.'
    print 'type ', type,'.'
    print 'proj ', proj,'.'
    print 'user ', user,'.'

    select = {}
    if status not in ["", None]:
      select['status'] = status.__str__()
    # this is going to be a little harder to work in to the Bean
    # if type not in ["", None]:
    #   select['type'] = type

    # really we need to be able to put in text here that allows more
    # than one project etc to be included. But what should the syntax be?
    if proj not in ["", None]:
      select['proj'] = proj
    if user not in ["", None]:
      select['user'] = user.__str__()

    print 'set select hash ', select

    self.guiParent.select = select
    self.guiParent.drawFrame()


if __name__ == "__main__":

  import sys
  import Tkinter

  root = Tkinter.Tk()
  root.withdraw()
   
  popup = TasksPopup(root)
  
  root.mainloop()
