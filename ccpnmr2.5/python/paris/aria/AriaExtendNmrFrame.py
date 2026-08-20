
"""
======================COPYRIGHT/LICENSE START==========================

AriaSetup.py: <write function here>

Copyright (C) 2005 Tim Stevens and Wolfgang Rieping (University of Cambridge)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

A copy of this license can be found in ../../../license/LGPL.license

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

Rieping W., Habeck M., Bardiaux B., Bernard A., Malliavin T.E.,
Nilges M.(2007) ARIA2: automated NOE assignment and data
integration in NMR structure calculation. Bioinformatics 23:381-382
===========================REFERENCE END===============================
"""

import os, sys
import Tkinter

from memops.gui.Button          import Button
from memops.gui.ButtonList      import ButtonList, UtilityButtonList
from memops.gui.DataEntry       import askString
from memops.gui.Entry           import Entry
from memops.gui.Frame           import Frame
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FileSelect      import FileType
#from memops.gui.Label           import Label
from memops.gui.LabeledEntry	import LabeledEntry
from memops.gui.LabelFrame      import LabelFrame
from memops.gui.MessageReporter import showWarning, showOkCancel, showYesNo
from memops.gui.PulldownList    import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix

from ccp.util.NmrCalc import getRunTextParameter
from ccpnmr.analysis.core.AssignmentBasic import getShiftLists
from ccpnmr.analysis.core.ExperimentBasic import getThroughSpacePeakLists
from paris.aria.AriaRunFrame import AriaRunFrame
from paris.aria.CcpnToAriaXml import makeAriaProject
from paris.aria.CcpnToAriaXml import CNS_EXE, WORKING_DIR, TEMP_DIR

from ccpnmr.analysis.popups.BasePopup import BasePopup

ARIA = 'ARIA'
  
import os

from memops.universal.Io import getTopDirectory

CHECK_INTERVAL = 60000 # Miliseconds

class AriaFrame(Frame):

  def __init__(self, parent, application, *args, **kw):

    Frame.__init__(self, parent=parent, **kw)
    self.expandGrid(1,1)
    
    self.parent = parent
    self.application = application
    self.project = project = application.project
    self.waiting = False
    self.ariaProjectPath = None
    self.buttons = []
    self.isInstalled = False
    self.runFrame = None
    
    # # # # # # # Check ARIA2 installation # # # # # # # 
    
    self.ariaRootDir = ariaRootDir = os.environ.get('ARIA2')
    
    if not ariaRootDir:
      print 'ARIA2 environment variable not set'
      self.isAriaInstalled = False
      #label = Label(self, text='ARIA 2 environment variable not set')
      #label.grid(row=1, column=0, sticky='ew')
      #return
    
    elif not os.path.exists(ariaRootDir):
      print 'ARIA2 path does not exist'
      self.isAriaInstalled = False
      #label = Label(self, text='ARIA2 path does not exist')
      #label.grid(row=1, column=0, sticky='ew')
      #return
           
    else: 
      sys.path.append(ariaRootDir)
      self.isAriaInstalled = True

    try:
      import aria2
      
    except ImportError:
      print 'Cannot import ARIA2 modules:\nCannot run ARIA locally.'
      self.isAriaInstalled = False
      #label = Label(self, text=)
      #label.grid(row=1, column=0, sticky='ew')
      #return
    
    if self.isAriaInstalled:
      for module in aria2.PATH_MODULES:
        path = os.path.join(ariaRootDir, module)
 
        if path not in sys.path:
          sys.path.append(path)

      failedModules = []
      for module in aria2.MODULES:
        try:
          __import__(module)
 
        except ImportError, msg:
          failedModules.append(module)
          print msg
 
      if failedModules:
        text = ', '.join(failedModules)
        self.isAriaInstalled = False
        print text
        #label = Label(self, text='Cannot import modules required by ARIA 2:\n%s' % text)
        #label.grid(row=1, column=0, sticky='ew')
        #return
 
    # # # # # # # # # # # # # # # 
    
    self.isInstalled = True
    
    if project:
      nmrProject = project.currentNmrProject
      self.calcStore = project.findFirstNmrCalcStore(name=ARIA, nmrProject=nmrProject) or \
                       project.newNmrCalcStore(name=ARIA, nmrProject=nmrProject)
    
    else:
      self.calcStore = None
    
    # # # # # # # ARIA RUN SETTINGS # # # # # # # # 
    
    self.runFrame = AriaRunFrame(self, project, self.calcStore, grid=(1,0), gridSpan=(1,2))
    
    # # # # # # # # # # # # # # #  # # # # # # # # #
            
    button = Button(self, text='ARIA Project File:',bd=1,bg='#C0E0FF',
                    command=self.selectAriaProjectFile, grid=(2,0))
                    
    self.ariaProjectPathEntry = Entry(self, text=self.ariaProjectPath,
                                grid=(2,1), sticky='w')
    
    #texts = ['Launch ARIA GUI','Setup Project'] # ,'Start Run']
    #commands = [self.runAriaGui,self.setupAriaProject] # ,self.startAriaRun]
    texts = ['Launch ARIA GUI','Setup Project', 'Submit to CcpnGrid'] # ,'Start Run']
    commands = [self.runAriaGui,self.setupAriaProject,self.getCcpnGridPassword] # ,self.startAriaRun]
    buttons = ButtonList(self, texts=texts, commands=commands,
                         grid=(3,0), gridSpan=(1,2))
    
    self.buttons.append(button)
    for button0 in buttons.buttons:
      button0.config(bg='#C0E0FF')   
      self.buttons.append(button0)
    
    self.administerNotifiers(self.application.registerNotify)

  def administerNotifiers(self, notifyFunc):
  
    if self.runFrame:
      self.runFrame.administerNotifiers(notifyFunc)
    
  def checkProjectPath(self, insistCnsExe=True):
    
    if not self.isAriaInstalled:
      showWarning('Failure','No working ARIA installation', parent=self)
      return
  
    run = self.runFrame.run
    if not run:
      showWarning('Failure','No ARIA run setup', parent=self)
      return 

    if self.ariaProjectPath is None:
      showWarning('Failure','ARIA project location not set', parent=self)
      return 

    cnsExe = getRunTextParameter(run, CNS_EXE)
    if not os.path.exists(cnsExe):
      showWarning('Failure','CNS executable file %s does not exist' % cnsExe, parent=self)
      if insistCnsExe:
        return 

    elif not os.path.isfile(cnsExe):
      showWarning('Failure','CNS executable location %s is not a file' % cnsExe, parent=self)
      return 
     
    if os.path.exists(self.ariaProjectPath):
      if showYesNo('File Exists','Overwrite existing ARIA project XML file?', parent=self):
        made = self.makeAriaProject()
      else:
        made = True
          
    else:
      made = self.makeAriaProject()
    
    if made:
      return True
 
  def runAriaGui(self):
  
    if self.checkProjectPath(insistCnsExe=False):  
      import aria2
      aria2.run_gui(self.ariaProjectPath) # , self.application)
 
  def setupAriaProject(self):
      
    if self.checkProjectPath(insistCnsExe=True):
      import aria2
      aria2.setup_aria(self.ariaProjectPath, self.ariaProjectPath)
      #os.system('/usr/bin/python /home/tjs23/ccpForge/extend-nmr/aria/aria2.py -s %s' % self.ariaProjectPath)
 
  def startAriaRun(self):
  
    if self.checkProjectPath(insistCnsExe=True):
      # TBD: Should quit current CCPN project due to concurrency issues
    
      import aria2
      aria2.run_aria(self.ariaProjectPath)
      #os.system('/usr/bin/python /home/tjs23/ccpForge/extend-nmr/aria/aria2.py %s' % self.ariaProjectPath)
  
  def selectAriaProjectFile(self):
    
    file_types = [ FileType("Aria", ["*.xml"]), FileType("All", ["*"]) ]

    popup = FileSelectPopup(self, file_types, dismiss_text='Cancel')

    result = popup.getFile()
    path = popup.getFile()
    if path:
      self.ariaProjectPath = path
      self.ariaProjectPathEntry.set(path)
    
    popup.destroy()
  
  def getCcpnGridPassword(self):
  
    run = self.runFrame.run
    if not run:
      showWarning('Failure','No ARIA run setup', parent=self)
      return 
      
    data = run.findFirstInput(className='MolSystemData')
    if not data:
      msg = 'No molecular system selected for ARIA run'
      showWarning('Failure', msg, parent=self)
      return

    data1 = run.findFirstInput(className='ConstraintStoreData')
    data2 = run.findFirstInput(className='PeakListData')
    if not (data1 or data2):
      msg = 'No constraint or peak list data selected for ARIA run'
      showWarning('Failure', msg, parent=self)
      return
          
    popup = AriaAskPassword(self)    
    
  def runCcpnGrid(self, userName, password):
    
    #
    # WIM TEMPORARY HACK!
    #
    
    try:
      from ccpnmr.workflow.Aria import AriaCcpnGrid
      
    except:
      showWarning('Code not available','The code to run ARIA using CCPNGRID is not available',parent=self)
      AriaCcpnGrid = None
      
    if AriaCcpnGrid:

      import time
      
      self.ariaGrid = ariaGrid = AriaCcpnGrid(ccpnProject=self.project,useGui=False,
                                              initializeObjects=False)
      
      # TODO here need to ask for username and password to log in!
      
      tgzFileName = ariaGrid.makeCcpnProjectTgz(baseNameSuffix='_ccpnGridRun')

      # TODO here also clearly need some authentication stuff, check whether login OK...
      ariaGrid.ariaRun(userName, password, tgzFileName)
      
      # TODO should here also (depending on situation) have 'Cancel run' button, ...
      # TODO obviously has to be nicer, plus more info in status window, probably.
      self.statusPopup = CcpnGridStatusPopup(self,'Starting',ariaGrid.uniqueIdentifier)
      self.statusPopup.update()
            
      self.after(CHECK_INTERVAL, self.timedCheckStatus)

  def timedCheckStatus(self):
    
    infoDict = self.ariaGrid.getStatusPageInfo()[self.ariaGrid.uniqueIdentifier]
    status = infoDict['status']
    
    if status in ('Finished','Failed'):
      self.statusPopup.destroy()
      return
    
    else:
      statusText = "%s\n%s" % (status,infoDict['iteration'])
      self.statusPopup.setStatus(statusText)
      self.statusPopup.update()
      self.after(CHECK_INTERVAL, self.timedCheckStatus)

  def makeAriaProject(self):
    
    if not self.isAriaInstalled:
      showWarning('Failure','No working ARIA installation', parent=self)
      return

    if not self.project:
      showWarning('Failure','No active project', parent=self)
      return
    
    calcStore = self.calcStore
    
    run = self.runFrame.run
    if not run:
      showWarning('Failure','No ARIA run setup', parent=self)
      return
      
    data = run.findFirstInput(className='MolSystemData')
    if not data:
      msg = 'No molecular system selected for ARIA run'
      showWarning('Failure', msg, parent=self)
      return

    data1 = run.findFirstInput(className='ConstraintStoreData')
    data2 = run.findFirstInput(className='PeakListData')
    if not (data1 or data2):
      msg = 'No constraint or peak list data selected for ARIA run'
      showWarning('Failure', msg, parent=self)
      return
    
    self.ariaProjectPath = self.ariaProjectPathEntry.get()
    if not self.ariaProjectPath:
      showWarning('Failure','No project file specified', parent=self)
      return
    
    from paris.aria.CcpnToAriaXml import makeAriaProject
      
    makeAriaProject(self.project, self.ariaProjectPath, run=run)
    
    return True

  def updateAll(self, project=None):

    if not self.isInstalled:
      return
      
    if project is not self.project:
    
      self.project = project
      fileName = ARIA+'_'+self.project.name+'.xml'
      
      repository = project.findFirstRepository(name='userData')
      url = repository.url
      self.ariaProjectPath = os.path.join(url.path, fileName)
    
    if project:
      nmrProject = project.currentNmrProject
      self.calcStore = project.findFirstNmrCalcStore(name=ARIA, nmrProject=nmrProject) or \
                       project.newNmrCalcStore(name=ARIA, nmrProject=nmrProject)

    self.ariaProjectPathEntry.set(self.ariaProjectPath)
    self.updateRunFrame()

  def updateAfter(self, obj=None):
  
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.updateRunFrame)    

  def updateRunFrame(self):
  
    project = self.project
    nmrProject = self.project.currentNmrProject
    calcStore = project.findFirstNmrCalcStore(name=ARIA, nmrProject=nmrProject) or \
                project.newNmrCalcStore(name=ARIA, nmrProject=nmrProject)
               
    self.runFrame.update(calcStore)
    self.waiting = False

  def destroy(self):        
  
    self.administerNotifiers(self.application.unregisterNotify)
    Frame.destroy(self)


class AriaAskPassword(BasePopup):

  def __init__(self,parent):

    self.parent      = parent  # GUI parent
    self.username    = None
    self.password    = None
    
    BasePopup.__init__(self, parent=parent, title='CcpnGrid Login Details')

  def body(self, guiFrame):

    guiFrame.grid_columnconfigure(1, weight=1)
    guiFrame.grid_rowconfigure(1, weight=1)

    self.usernameEntry = LabeledEntry(guiFrame, 'username')
    self.usernameEntry.grid(row=0)

    self.passwordEntry = LabeledEntry(guiFrame, 'password', show="*")
    self.passwordEntry.grid(row=1)
    
    texts = ['Submit ARIA run']
    commands = [self.runAria]
    buttonList = UtilityButtonList(guiFrame, texts=texts, commands=commands,
                                   doHelp=False, doClone=False, expands=True,
                                   grid=(2,0))

  def runAria(self):

    warning = ''
    
    if self.usernameEntry.getEntry():
      self.username = self.usernameEntry.getEntry()
    else:
      warning += ' No username defined ' 
    
    if self.passwordEntry.getEntry():
      self.password = self.passwordEntry.getEntry()
    else:
      warning += ' No password defined '      
    
    if warning:
      showWarning('Failure', warning, parent=self) 
    
    if self.username and self.password:
      self.parent.runCcpnGrid(self.username, self.password)
      self.destroy()
      
  def close(self):
  
    BasePopup.destroy(self)
  
                
class AriaPopup(BasePopup):
  """
  """

  def __init__(self, parent, **kw):

    self.parent = parent

    BasePopup.__init__(self, parent=parent, title='ARIA Setup', **kw)
    
  def body(self, guiFrame):

    self.geometry('740x600')
  
    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    self.frame = AriaFrame(guiFrame, self.parent)
    self.frame.grid(row=0, column=0, sticky='nsew')
    self.frame.updateAll(self.project)
  
  def open(self):
  
    BasePopup.open(self)
   
   
  def close(self):
  
    BasePopup.close(self)
    
    
  def destroy(self):
  
    BasePopup.destroy(self)  

# TODO HACK WIM

from memops.gui.Label import Label
from memops.gui.Util import createDismissHelpButtonList
from memops.gui.Text import Text

import webbrowser

#from ccpnmr.format.gui.BasePopup import BasePopup

class CcpnGridStatusPopup(BasePopup):

  def __init__(self, parent, statusText, projectName):
  
    self.statusText = statusText
    self.projectName = projectName
    
    self.browser = None

    BasePopup.__init__(self, parent=parent, title='CCPN Grid progress', modal=False, transient=True)

  def body(self, master):
  
    row = 0
    columnspan = 2

    label = Label(master, text= "Project name:")
    label.grid(row=row, column=0, sticky=Tkinter.W)

    label = Label(master, text= self.projectName)
    label.grid(row=row, column=1, sticky=Tkinter.W)

    row += 1
    
    label = Label(master, text= "Run status:")
    label.grid(row=row, column=0, sticky=Tkinter.W)

    self.status = Text(master, text=self.statusText, width = 20, height = 2)
    self.status.grid(row=row, column=1, sticky=Tkinter.W)

    row += 1
    
    texts = [ 'Show web page' ]
    commands = [ self.showWebPage ]
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands) #, help_url=self.help_url
    buttons.grid(row=row, column=0, columnspan = columnspan)

  def setStatus(self,text):
   
    self.status.setText(text)
    self.update()

  def showWebPage(self):
  
    if not self.browser:
      self.browser = webbrowser.get()
      
    self.browser.open("http://webapps.ccpn.ac.uk/ccpngrid/status")

if __name__ == '__main__':

  root = Tkinter.Tk()
 
  gsp = CcpnGridStatusPopup(root,'Increments','test')
  gsp.update()
  
  #print dir(gsp)

  import time
  for i in range(5):
    gsp.setStatus(str(i))
    time.sleep(1)

  root.mainloop()
