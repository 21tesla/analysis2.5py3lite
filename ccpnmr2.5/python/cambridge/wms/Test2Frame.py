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

#from memops.gui.ScrolledFrame import ScrolledFrame
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.Text import Text
from memops.gui.Button import Button



# # # # #   B I G   T O - D O  P O I N T S   # # # # # 
#


# # # # #   L A T E R   T O - D O  P O I N T S    # # # # #
#
 
class Test1Popup(BasePopup):

  def __init__(self, parent):

    self.parent      = parent
    
    BasePopup.__init__(self, parent=parent, title='Repository')


  # there is a slight logical problem here. This frame refers to a
  # single task whereas the GUI as a whole contains many, many
  # tasks. How is the task represented in this frame selected?
                       
  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    frame = Test2Frame(guiFrame, basePopup=self)
    frame.grid(row=0, column=0, sticky='nsew')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    
    self.geometry('200x50')
    
    self.update_idletasks()


class Test2Frame(Frame):

  
  def __init__(self, guiParent, basePopup):

    self.guiParent = guiParent
    self.basePopup = basePopup

    #self.registerNotify=basePopup.registerNotify
    #self.unregisterNotify=basePopup.unregisterNotify

    self.task_id = None

    self.basePopup.frameShortcuts['Task'].frameShortcuts['Test2'] = self

    # NOTE: JMCI

    # We can reasonably set the repository from the current repository. It is
    # harder to know how project/projectVersion and task interract

    # TODO
    # decide whether we really need this
    self.repository = self.basePopup.repList.currentRepository
    
    Frame.__init__(self, guiParent)
  
    # set up the grid

    self.grid_columnconfigure(0, weight=0, minsize=20)
    self.grid_columnconfigure(1, weight=0, minsize=30)
    self.grid_columnconfigure(2, weight=1, minsize=20)
    self.grid_columnconfigure(3, weight=0, minsize=20)

    self.grid_rowconfigure(0, weight=0, minsize=20)
    self.grid_rowconfigure(1, weight=0, minsize=10)
    self.grid_rowconfigure(2, weight=0, minsize=10)
    self.grid_rowconfigure(3, weight=0, minsize=10)
    self.grid_rowconfigure(4, weight=0, minsize=20)


    # we define as much as possible once on startup and use the
    # grid manager to control the display thereafter

    # when no data
    
    self.noData_widgets = []
    
    self.null_task_label = Label(self, text='No Task selected')
    self.noData_widgets.append(self.null_task_label)


    # when data

    self.data_widgets = []

    self.task_rep_title = Label(self, text='Task:', font='Helvetica16')
    self.data_widgets.append(self.task_rep_title)

    self.task_rep_label = Label(self, text='Repository:')
    self.data_widgets.append(self.task_rep_label)
    self.task_rep_value = Label(self, text='')
    self.data_widgets.append(self.task_rep_value)

    self.task_status_label = Label(self, text='Status:')
    self.data_widgets.append(self.task_status_label)
    self.task_status_value = Label(self, text='')
    self.data_widgets.append(self.task_status_value)

    self.task_details_label = Label(self, text='Details:')
    self.data_widgets.append(self.task_details_label)
    self.task_details_value = Label(self, text='')
    self.data_widgets.append(self.task_details_value)

    self.task_date_label = Label(self, text='Date:')
    self.data_widgets.append(self.task_date_label)
    self.task_date_value = Label(self, text='')
    self.data_widgets.append(self.task_date_value)

    self.task_par1_label = Label(self, text='Par1:')
    self.data_widgets.append(self.task_par1_label)
    self.task_par1_value = Label(self, text='')
    self.data_widgets.append(self.task_par1_value)

    self.task_par2_label = Label(self, text='Par2:')
    self.data_widgets.append(self.task_par2_label)
    self.task_par2_value = Label(self, text='')
    self.data_widgets.append(self.task_par2_value)

    self.task_input_label = Label(self, text='Input Data:')
    self.data_widgets.append(self.task_input_label)
    self.task_input_value = Label(self, text='')
    self.data_widgets.append(self.task_input_value)

    self.task_output_label = Label(self, text='Output Data:')
    self.data_widgets.append(self.task_output_label)
    self.task_output_value = Label(self, text='')
    self.data_widgets.append(self.task_output_value)

    self.drawFrame()


  def drawFrame(self):

    # build up the body.

    # what are the rules governing when this tab should show something and
    # when it should not? In fact, there is the whole question of how the
    # two different sections (file management and task management) should
    # interract sensibly.

    # This needs to be composed. We know the task type (from the tab) so what
    # we need to know is how to handle ids. should we generate a placeholder
    # id, or would it be better not to assign an id to begin with or else
    # just use id(new)

    self.task_id = self.basePopup.currentTask

    print 'trying to render task with serial ', self.task_id 

    if self.task_id == None:

      for widget in self.data_widgets:
        widget.grid_remove()
      
      self.null_task_label.grid(row=1, column=1, columnspan=4, sticky='nswe')
     
    else:
    
      for widget in self.noData_widgets:
        widget.grid_remove()

      loc = SharedBeanServiceLocator()
      port = loc.getSharedBean()

      # need two calls, one for the Project and one for the version
      hm1 = {'serial': int(self.task_id) }
      wsstr_in1 = WSString(hm1)

      request1 = getFieldsWithMethod();
      request1._arg0 = 'org.pimslims.applet.server.TaskBean';
      request1._arg1 = 'getTest1Fields'
      request1._arg2 = wsstr_in1.str

      response1 = port.getFieldsWithMethod(request1)

      wsstr_out1 = WSString(response1._return)
      hm_tv = wsstr_out1.getStructHM()

      print 'RESULT ', hm_tv

 
      self.task_rep_title.set('Task: ' + self.task_id)
      self.task_rep_title.grid(row=1, column=1, sticky='w')

      self.repository = self.basePopup.repList.currentRepository
      self.task_rep_label.grid(row=2, column=1, sticky='w')     
      self.task_rep_value.set(self.repository.user + '@' + self.repository.name + '  ( ' + self.repository.connect + ' )')
      self.task_rep_value.grid(row=2, column=2, sticky='w')

      # FIXME JMCI

      # This actually needs to be editable

      self.task_status_label.grid(row=3, column=1, sticky='w')
      self.task_status_value.set(hm_tv['status'])
      #self.task_status_value.set('PENDING')
      self.task_status_value.grid(row=3, column=2, sticky='w')

      self.task_par1_label.grid(row=4, column=1, sticky='w')
      self.task_par1_value.set(hm_tv['par1'])
      self.task_par1_value.grid(row=4, column=2, sticky='w')

      self.task_par2_label.grid(row=5, column=1, sticky='w')
      self.task_par2_value.set(hm_tv['par2'])
      self.task_par2_value.grid(row=5, column=2, sticky='w')

      # this will depend on the status
      #self.task_date_label.grid(row=4, column=1, sticky='w')
      #self.task_date_value.set('11 FEB 2009:11:36') 
      #self.task_date_value.grid(row=4, column=2, sticky='w')

      # FIXME JMCI

      # This is not in the current model
      #self.task_details_label.grid(row=4, column=1, sticky='w')
      #self.task_details_value.set(hm_pv['details'])
      #self.task_details_value.grid(row=4, column=2, sticky='w')



      inputProject = 'N/A'
      # This needs to be populated from the fields hash
      # It is also going to be dependent on the precise task
      if hm_tv.has_key('inputVersion'):
        inputProjectName = hm_tv['inputVersion']['project']['name']
        inputVersionTag = hm_tv['inputVersion']['versionTag']
        inputProject = inputProjectName + '::' + inputVersionTag

      self.task_input_label.grid(row=6, column=1, sticky='w')
      self.task_input_value.set(inputProject)
      self.task_input_value.grid(row=6, column=2, sticky='w')


      outputProject = 'N/A'
      # This needs to be populated from the fields hash
      # It is also going to be dependent on the precise task
      if hm_tv.has_key('generatedVersion'):
        outputProjectName = hm_tv['generatedVersion']['project']['name']
        outputVersionTag = hm_tv['generatedVersion']['versionTag']
        outputProject = outputProjectName + '::' + outputVersionTag

      self.task_output_label.grid(row=7, column=1, sticky='w')
      self.task_output_value.set(outputProject)
      self.task_output_value.grid(row=7, column=2, sticky='w')

      # FIXME JMCI

      # would be sensible to have a browse/edit button here although
      # that would have to depend on the status

      #self.button_file1_browse = Button(self, text='Browse', command=self.tmpCall)
      #self.button_file1_browse.grid(row=r, column=3, sticky='w')


      #self.button_file2_browse = Button(self, text='Browse', command=self.tmpCall)
      #self.button_file2_browse.grid(row=r, column=3, sticky='w')

      # FIXME JMCI

      # again, this really depends on what the status is. Need to decide
      # the transitions
      button_submit = Button(self, text='Submit', command=self.submit)
      button_submit.grid(row=8, column=2, sticky='w')


  def submit(self):

    print 'either passing to server or running locally'

    # get a port proxy instance
    # this should probably belong to the repository directly

    pars={'id':self.value_par1.getText()[:-1], \
          'status':self.value_par2.getText()[:-1] }

    wsstr = WSString(pars)
    
    loc = SharedBeanServiceLocator()
    self.port = loc.getSharedBean()

    request = record();

    request._arg0 = 'org.pimslims.applet.server.TaskBean';
    request._arg1 = 'record';
    request._arg2 = wsstr.str

    ## # call the remote method
    response = self.port.record(request)

    # now need to redraw task-related frames
    print self.basePopup.__dict__
    if self.basePopup.frameShortcuts.has_key('Tasks'):
      self.basePopup.frameShortcuts['Tasks'].drawFrame()
    

  def tmpCall(self):

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
  
    #self.administerNotifiers(self.basePopup.unregisterNotify)
    Frame.destroy(self)


if __name__ == "__main__":

  import sys
  import Tkinter

  root = Tkinter.Tk()
  root.withdraw()
   
  popup = Test2Popup(root)
  
  root.mainloop()
