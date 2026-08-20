import os, sys, string
import Tkinter

from memops.gui.Button          import Button
from memops.gui.ButtonList      import ButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.Entry           import Entry
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FileSelect      import FileType
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Frame           import Frame
from memops.gui.IntEntry        import IntEntry
from memops.gui.Label           import Label
from memops.gui.MessageReporter import showError, showInfo, showWarning, showYesNo
from memops.gui.MultiWidget     import MultiWidget
from memops.gui.PulldownMenu    import PulldownMenu
from memops.gui.PulldownList    import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.TabbedFrame     import TabbedFrame
from memops.gui.Util            import createDismissHelpButtonList
from memops.editor.BasePopup    import BasePopup
from memops.general.Io          import saveProject
from ccpnmr.analysis.core.ExperimentBasic import getThroughSpacePeakLists

from cambridge.isd.NmrCalcExchange import nmrCalcRunToIsd, isdToNmrCalcRun, getIsdNmrCalcStore, PROVISIONAL

redColor   = '#FF8080'
greenColor = '#80FF80'

EXIT_TIMEOUT = 2.1*60.
SAMPLER_TIMEOUT = 2.*60
CHECK_INTERVAL = 1

class SetupDataSet:

  def __init__(self, dataType, file='', format=None, key='', name='', specsObj=None):
    """Object to hold dataset info prior to and after commitment"""

    from Isd import setup

    if not format:
      format = setup.CCPN

    if key and not name:
      name = key

    self.dataType = dataType
    self.file     = file
    self.format   = format
    self.key      = key
    self.name     = name
    self.specsObj = specsObj

class IsdFrame(Frame):

  def __init__(self, parent, ccpnProject):

    ## 1) .__init__ cretes IsdFrame with empty CCPN project
    ##     (called by ExtendNmrGui.ApplicationPopup.initIsd)
    ##
    ## 2) .updateAll sets the project
    ##     (called by ExtendNmrGui.ApplicationPopup.initProject)

    self.parent = parent
    self.rowObj = None
    self.table = None
    self.cleanProject()
    self.haveIsdInstalled = False

    Frame.__init__(self, parent=parent)

    self.checkIsdInstallation()
    if not self.haveIsdInstalled:
      return

    # Ensure that the second row and first column in frame expand

    self.grid_rowconfigure(2, weight=1)
    self.grid_columnconfigure(0, weight=1, minsize=300)

    initialrows = 15

    #
    # Pulldown menu to choose ISD projects from key
    #

    self.topFrame1 = Frame(self)
    self.topFrame1.grid(row=0, column=0, sticky='ew')
    self.topFrame1.grid_columnconfigure(3, weight=1)

    label = Label(self.topFrame1, text='ISD Settings: ')
    label.grid(row=0, column=0, sticky='w')

    self.isdRunPulldown = PulldownList(self.topFrame1, callback=self.selectIsdRun)
    self.isdRunPulldown.grid(row=0, column=1, sticky='w')

    label = Label(self.topFrame1, text='   CCPN Project Name: ')
    label.grid(row=0, column=2, sticky='w')

    self.ccpnProjectNameLabel = Label(self.topFrame1, text=self.getCcpnProjectName())
    self.ccpnProjectNameLabel.grid(row=0, column=3, sticky='w')

    texts = ['New ISD Run','Delete ISD Run',]
    commands = [self.newRun, self.deleteRunSettings]
    buttons = ButtonList(self.topFrame1, texts=texts, commands=commands,
                         grid=(0,4), sticky='e')
    
    self.newButton, self.deleteButton = buttons.buttons

    self.topFrame2 = Frame(self)
    self.topFrame2.grid(row=1, column=0, sticky='ew')

    label = Label(self.topFrame2, text='Path: ')
    label.grid(row=0, column=0, sticky='w')

    self.ccpnProjectPathLabel = Label(self.topFrame2,
                                      text=self.getCcpnProjectPath())
    self.ccpnProjectPathLabel.grid(row=0, column=1, sticky='w')

    #
    # Tabbed frame
    #

    options = ['  General  ','Molecules & Structures',
               'Replica Exchange MC','  Analyses  ','Experimental Data']

    self.tabbedFrame = TabbedFrame(self, options=options, toggleOff=False, selected=0)
    self.tabbedFrame.grid(row=2, column=0, sticky='nsew')

    genFrame, molFrame, monteCarloFrame, anaFrame, dataFrame = self.tabbedFrame.frames

    # Setup generic widgets for use in tables

    self.stringEntry  = Entry(self, returnCallback=self.setValue)
    self.intEntry     = IntEntry(self, returnCallback=self.setValue,)
    self.floatEntry   = FloatEntry(self, returnCallback=self.setValue)
    self.pulldownMenu = PulldownMenu(self, callback=self.setValue, do_initial_callback=False)
    self.multiString  = MultiWidget(self, Entry, callback=self.setValue, minRows=1, values=[], useImages=False)
    self.multiCheck   = MultiWidget(self, CheckButton, callback=self.setValue, minRows=1, values=[], useImages=False)

    # A couple specific widgets for the data set table

    self.dataNameEntry      = Entry(self, returnCallback=self.setDataName)
    self.dataKeyPulldown    = PulldownMenu(self, callback=self.setDataKey, do_initial_callback=False)
    self.dataFormatPulldown = PulldownMenu(self, callback=self.setDataFormat, do_initial_callback=False)
    self.dataKeyEntry       = Entry(self, returnCallback=self.setDataKey)

    # Setup generic table headings, justification and widget getters/setters

    headingList      = ['#','Parameter','Value','Description']
    justifyList      = [None, None, 'left', 'left']
    editWidgets      = [None, None, True, None]
    editGetCallbacks = [None, None, self.getValue, None]
    editSetCallbacks = [None, None, self.setValue, None]

    #
    # General frame
    #

    genFrame.grid_columnconfigure(1, weight=1)
    genFrame.grid_rowconfigure(1, weight=1)

    self.generalMatrix = ScrolledMatrix(genFrame, headingList=headingList,
                                        justifyList=justifyList,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        multiSelect=False, initialRows=initialrows,
                                        passSelfToCallback=True,
                                        callback=self.selectRowObj)

    self.generalMatrix.grid(row=1, column=0, columnspan=2, sticky='nsew')
    self.generalMatrix.refreshFunc = self.updateGeneral

    #
    #  Mol & structure frame
    #

    molFrame.grid_columnconfigure(0, weight=1)
    molFrame.grid_rowconfigure(0, weight=1)
    self.molMatrix = ScrolledMatrix(molFrame, headingList=headingList,
                                        justifyList=justifyList,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        multiSelect=False, initialRows=initialrows,
                                        passSelfToCallback=True,
                                        callback=self.selectRowObj)

    self.molMatrix.grid(row=0, column=0, sticky='nsew')
    self.molMatrix.refreshFunc = self.updateMolStructure

    #
    #  Replica exchange MC frame
    #

    monteCarloFrame.grid_columnconfigure(0, weight=1)
    monteCarloFrame.grid_rowconfigure(0, weight=1)
    self.monteCarloMatrix = ScrolledMatrix(monteCarloFrame, headingList=headingList,
                                        justifyList=justifyList,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        multiSelect=False, initialRows=initialrows,
                                        passSelfToCallback=True,
                                        callback=self.selectRowObj)

    self.monteCarloMatrix.grid(row=0, column=0, sticky='nsew')
    self.monteCarloMatrix.refreshFunc = self.updateMonteCarlo

    #
    #  Analyses frame
    #

    anaFrame.grid_columnconfigure(0, weight=1)
    anaFrame.grid_rowconfigure(0, weight=1)
    self.analysisMatrix = ScrolledMatrix(anaFrame, headingList=headingList,
                                        justifyList=justifyList,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        multiSelect=False, initialRows=initialrows,
                                        passSelfToCallback=True,
                                        callback=self.selectRowObj)

    self.analysisMatrix.grid(row=0, column=0, sticky='nsew')
    self.analysisMatrix.refreshFunc = self.updateAnalyses

    #
    #  Experimental data frame
    #

    headingList      = ['Data Type','Name','Commit &\nEdit Details','Format','Key','Filename']
    justifyList      = [None, 'center','center','center','center','left']
    editWidgets      = [None, self.dataNameEntry, None, self.dataFormatPulldown, True, None]
    editGetCallbacks = [None, self.getDataName, self.commitDataToggle, self.getDataFormat, self.getDataKey, self.getDataFile]
    editSetCallbacks = [None, self.setDataName, None, self.setDataFormat, self.setDataKey, None]
    dataFrame.grid_columnconfigure(0, weight=1)
    dataFrame.grid_rowconfigure(0, weight=1)
    self.dataMatrix = ScrolledMatrix(dataFrame, headingList=headingList,
                                        justifyList=justifyList,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        multiSelect=False, initialRows=initialrows/2,
                                        passSelfToCallback=True,
                                        callback=self.selectExpData)

    self.dataMatrix.grid(row=0, column=0, sticky='nsew')
    self.dataMatrix.refreshFunc = self.updateExpData
    self.dataMatrix.doEditMarkExtraRules = self.dataMatrixEditRules

    texts = ['Add NOE Distance','Add Distance','Add J Coupling','Add RDC ']
    commands = [self.addNOE,self.addDistance,self.addJCoupling,self.addRDC]

    self.expDataButtons1 = ButtonList(dataFrame, texts=texts, expands=True,
                                                     commands=commands)
    self.expDataButtons1.grid(row=1, column=0, sticky='ew')

    texts = ['Add Dihedral','Add H-Bond','Add Disulfide','Remove Selected']
    commands = [self.addDihedral,self.addHBond,self.addDisulfide,self.removeSelected]

    self.expDataButtons2 = ButtonList(dataFrame, texts=texts, expands=True,
                                                     commands=commands)
    self.expDataButtons2.grid(row=2, column=0, sticky='ew')

    headingList      = ['#','Parameter','Value','Description']
    justifyList      = [None, None, 'left', 'left']
    editWidgets      = [None, None, True, None]
    editGetCallbacks = [None, None, self.getValue,   None]
    editSetCallbacks = [None, None, self.setValue,   None]
    dataFrame.grid_rowconfigure(3, weight=1)
    self.dataDetailMatrix = ScrolledMatrix(dataFrame, headingList=headingList,
                                        justifyList=justifyList,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        multiSelect=False, initialRows=initialrows/2,
                                        passSelfToCallback=True,
                                        callback=self.selectRowObj)

    self.dataDetailMatrix.grid(row=3, column=0, sticky='nsew')
    self.dataDetailMatrix.refreshFunc = self.updateExpData
    self.dataDetailMatrix.doEditMarkExtraRules = self.dataDetailMatrixEditRules

    #
    #  Bottom buttions
    #

    frame1 = Frame(self)
    frame1.grid(row=3, column=0, sticky='ew')
    frame1.grid_columnconfigure(0, weight=1)

    texts    = ['Run ISD!','Stop ISD!','Save ISD Settings']
    commands = [self.cmdRun,self.cmdStop,self.saveSettings]
    self.bottomButtons1 = ButtonList(frame1, texts=texts, commands=commands)
    self.bottomButtons1.grid(row=0, column = 0, sticky = 'ew')

    frame2 = Frame(self)
    frame2.grid(row=4, column=0, sticky='ew')
    frame2.grid_columnconfigure(1, weight=1)

    label = Label(frame2, text='Simulation:', grid=(0,0), sticky='w')

    texts    = ['Write ISD project file', 'Info','Show','Energies','Rates','Save','Report']
    commands = [self.cmdWrite, self.cmdInfo,self.cmdShow,self.cmdEnergies,self.cmdRates,
                self.cmdSave,self.cmdReport]
    self.bottomButtons2 = ButtonList(frame2, texts=texts, commands=commands)
    self.bottomButtons2.grid(row=0, column=1, sticky='ew')

    self.allWidgets = self.bottomButtons1.buttons + self.bottomButtons2.buttons + \
                      self.expDataButtons1.buttons + self.expDataButtons2.buttons + \
                      [self.isdRunPulldown,]+[self.deleteButton,self.newButton]
    self.saveButton = self.bottomButtons1.buttons[2]
    self.runButton = self.bottomButtons1.buttons[0]
    self.stopButton = self.bottomButtons1.buttons[1]
    self.writeButton = self.bottomButtons2.buttons[0]
    self.reportButtons = [self.bottomButtons2.buttons[1], self.bottomButtons2.buttons[-1]]
    self.samplerButtons = self.bottomButtons2.buttons[2:-1]
    self.datasetButtons = self.expDataButtons1.buttons + self.expDataButtons2.buttons

    self.updateAll(ccpnProject)

  def cleanProject(self):

    ## ISD project data is stores as an NmrCalc.run

    from threading import Lock, Event

    self.busylock = Lock()
    self.halted = Event()

    self.debug       = False

    ## updated in updateAll

    self.ccpnProject = None

    ## updated in updateRunPulldown

    self.run = None

    ## updated in selectIsdRun

    self.sim         = None
    self.dataSets    = []
    self.dataSet     = None

    ## updated in cmdRun

    self.simulation  = None
    self.sampler     = None

  def closeProject(self):

    ## called by ExtendNmrGui.ApplicationPopup.closeProject

    if self.busylock.locked():

      showWarning('Close CCPN project failure','ISD process is busy... Try later.')
      return False

    if self.isRunning():

      if not showYesNo('Close ISD project',
                       'ISD simulation is running. Terminate?',
                       parent=self):
        return False

      print 'DEBUG!!! cmdStop is being called from closeProject'
      if not self.cmdStop():
        return False

      self.halted.wait(timeout = EXIT_TIMEOUT)

    self.cleanProject()
    if self.haveIsdInstalled:
      self.updateRunPulldown()
        
    self.updateAll()

    return True

  def updateRunPulldown(self):

    ## Called from:
    ##   .closeProject
    ##   .saveSettings
    ##   .updateAll

    index = 0
    names = []
    runs = []
    run = self.run

    if self.ccpnProject:
      calcStore = getIsdNmrCalcStore(self.ccpnProject)
      runs  += [r for r in calcStore.sortedRuns() if r.status == PROVISIONAL]
      names += ['Run %d' % r.serial for r in runs]
      
      if runs:
        if run not in runs:
          run = runs[-1]
          
        index = runs.index(run)
          
      else:
        run = None

    self.selectIsdRun(run)
    self.isdRunPulldown.setup(names, runs, index)

  def initSettings(self):

    ## When new ISD project is created the path
    ## should be updated

    self.sim.ccpn.project_filename = self.ccpnProject.findFirstRepository(name='userData').url.path

    ## Remove deprecated (can be called only once)

    if 'n' in self.sim.replica:
      self.sim.replica.dict.pop('n')

    ## Continuing simulation:
    ##   sim.temp_path is initialised from
    ##   in sim.base_temp_path (because sim.temp_path changes
    ##   when launchSimulation, but sim.base_temp_path stays)
    ##
    ## New simulation:
    ##   sim.base_temp_path is initialised from
    ##   sim.temp_path (because sim.base_temp_path does not
    ##   exist before sim.finalise())

    if 'base_temp_path' in self.sim:
      self.sim.temp_path = self.sim.base_temp_path
    else:
      self.sim.base_temp_path = self.sim.temp_path

    ## Allow to terminate during publish (Pyro grid)
    ## because the termination is managed by the GUI

    self.sim.terminate_during_publish = True

    ## TODO: not pretty...

    if not 'gnuplot' in self.sim.analysis:
      self.sim.analysis.gnuplot = 'gnuplot'

    import Isd.gnuplot as g
    g.executable = self.sim.analysis.gnuplot
    g.old_version = g.is_old_version()

  def updateSettings(self):

    self.sim.temp_path = self.sim.base_temp_path

  def selectIsdRun(self, run):

    ## Called from:
    ##  * .updateRunPulldown
    ##     .closeProject (no project keys)
    ##     .saveSettings (adds current isd project key)
    ##     .updateAll (when new project is open)
    ##  * from the widget itself

    ## Resets the settings for the selected ISD project:
    ##   .sim, .simulation, .sampler
    ##
    ## Therefore should not continue if no .ccpnProject or
    ## simulation is running or busy

    if not self.ccpnProject:
      return

    if self.isRunning() or self.busylock.locked():
      return
    
    if not run:
      return
    
    if run is not self.run:
      self.updateNmrCalc()
      self.run = run
      self.sim = nmrCalcRunToIsd(run)
      self.dataSets = []
      self.dataSet = None

      if str(self.sim.data_sets) != '<>':
        for data_specs in self.sim.data_sets:
          dataSet = SetupDataSet(data_specs.data_type, data_specs.filename,
                                 data_specs.format, data_specs.key,
                                 data_specs.name, data_specs)
          self.dataSets.append(dataSet)

      self.simulation = None
      self.sampler = None
      self.initSettings()

      self.updateAll()

  """
  def chooseNmrProject(self, name=''):

    if not self.ccpnProject:
      return

    if not name:
      if self.ccpnProject.currentNmrProject:
        self.nmrProject = self.ccpnProject.currentNmrProject

      elif self.ccpnProject.sortedNmrProjects():
        self.nmrProject = self.ccpnProject.sortedNmrProjects()[0]

    else:
      self.nmrProject = self.ccpnReader.get_nmr_project(name)

    if not self.nmrProject:
      showWarning('Warning','Could not find NMR project within CCPN')
  """

  def saveProject(self):

    if self.ccpnProject:

      try:
        saveProject(self.ccpnProject, createFallback=True)
        print 'Successfully saved project'

      except IOError, msg:
        print 'ERROR:', msg
        showError('Saving file', str(msg), parent=self)
        return False

  def saveSettings(self):
  
    from Isd.setup import CCPN
    
    if self.sim.molecule.format == CCPN:

      if not self.sim.ccpn.export.molecular_system_name:
        showWarning('Warning','No CCPN molecular system specified (to retrieve the sequence).')

      if not self.sim.molecule.key:
        showWarning('Warning','No CCPN molecule key specified (to retrieve the sequence).')

    if not self.sim.ccpn.nmr_project_name:
      showWarning('Failure','Could not save under the CCPN project: No NmrProject specified.')
      return

    if not self.sim.ccpn.export.project.key:
      showWarning('Failure','Could not save under the CCPN project: No ISD project key specified.')
      return

    history = self.getHistory()

    if history:
      if not self.historyOk(history):
        return

    if self.ccpnProject:
      isdToNmrCalcRun(self.sim, self.ccpnProject, self.run)

      self.updateRunPulldown()
      self.saveProject()

  def newRun(self):

    if self.ccpnProject:
      calcStore = getIsdNmrCalcStore(self.ccpnProject)
      run = calcStore.newRun(status=PROVISIONAL)
      self.selectIsdRun(run) 

  def deleteRunSettings(self):

    if self.run:
      self.run.delete()

    self.updateRunPulldown()
    self.saveProject()

  def getHistory(self):

    ## returns history if it exists, None otherwise

    history_path = os.path.join(self.sim.working_path,
                                '%s_history' % self.sim.fileroot)

    if os.path.exists(history_path):

      try:
        from Isd.utils import Load
        history = Load(history_path)
        return history

      except:
        pass

    return None

  def historyOk(self, history):

    if not len(history['states']) == self.sim.replica.n_replicas:
      showWarning('Failure', 'Different number of replicas used in the previous simulation. '+\
                  'Please change the settings to match the previous simulation '+\
                  'or move the simulation results to another location.')
      return False

    data_old = history['states'][0].sub_states.keys()
    data_new = [d.name for d in self.sim.data_sets]
    data_old.sort()
    data_new.sort()

    if not data_old == data_new:
      showWarning('Failure','Different data used in the previous simulation. '+\
                  'Please change the settings to match the previous simulation '+\
                  'or move the simulation results to another location.')
      return False

    return True

  def settingsOk(self):

    ## returns False when settings are missing

    from Isd import setup

    if not self.ccpnProject or not self.sim:
      return False

    if not os.path.exists(os.path.expanduser(self.sim.cns_executable)):
      showWarning('Failure','CNS binary "%s" not found.' % self.sim.cns_executable)
      return False

    if not os.path.exists(os.path.expanduser(self.sim.replica.python)):
      showWarning('Failure','Python binary "%s" not found.' % self.sim.replica.python)
      return False

    if not os.path.exists(os.path.expanduser(self.sim.working_path)):
      if showYesNo('Warning',' Working path "%s" not found. Create path?' \
                   % self.sim.working_path):
        try:
          os.makedirs(os.path.expanduser(self.sim.working_path))
        except Exception, msg:
          print 'ERROR:', msg
          showError('Creating working path..',
                    'Could not create "%s".' % self.sim.working_path,
                    parent=self)
          return False
      else:
        return False

    if self.sim.molecule.format == setup.CCPN:

      if not self.sim.ccpn.export.molecular_system_name:
        showWarning('Failure','No CCPN molecular system specified (to retrieve the sequence). '+\
                    'Please specify molecular system name or specify correct sequence file name and format.')
        return False

      if not self.sim.molecule.key:
        showWarning('Failure','No CCPN molecule key specified (to retrieve the sequence).'+\
                    'Please specify molecular system name or specify correct sequence file name and format.')
        return False

    elif not self.sim.molecule.filename:

      showWarning('Failure','No molecule sequence file specified.')
      return False

    elif not os.path.exists(self.sim.molecule.filename):
      showWarning('Failure',': Molecule sequence file is not found.')
      return False

    if not self.dataSets:
      showWarning('Failure','No experimental data is selected for the simulation.')
      return False

    for data in self.dataSets:
      if not data.format == setup.CCPN and \
             not os.path.exists(os.path.expanduser(data.file)):
        showWarning('Failure','Data file "%s" not found. Please remove missing data set.' \
                    % os.path.expanduser(data.file))
        return False

    history = self.getHistory()

    if history:
      if not showYesNo('Warning','Previous simulation exists. If you have changed '+\
                       'simulation settings since the old run, they will be overwritten. '+\
                       'Current simulation will use conformations from the old run as an initial '+\
                       'conformation for this run. Continue?'):
        return False

      if not self.historyOk(history):
        return False

    return True

  def checkIsdInstallation(self):

    self.haveIsdInstalled = False
    isdRootDir = os.environ.get('ISD_ROOT')

    if not isdRootDir:
      label = Label(self, text='ISD_ROOT environment variable not set')
      label.grid(row=1, column=0, sticky='ew')
      return

    elif not os.path.exists(isdRootDir):
      label = Label(self, text='ISD_ROOT path does not exist')
      label.grid(row=1, column=0, sticky='ew')
      return

    modules = os.path.join(isdRootDir, 'src', 'py')
    if not modules in sys.path:
      sys.path.insert(0, os.path.join(isdRootDir, 'src', 'py'))

    try:
      from Isd import setup

    except ImportError, msg:
      Label(self, text='Cannot import ISD modules: Check installation and ISD_ROOT environment variable', grid=(1,0))
      Label(self, text='Reported error: "%s"' % msg, grid=(2,0))
      if str(msg) == 'No module named _isd':
         Label(self, text='C code may not be properly compiled', grid=(3,0))
      
      return

    self.haveIsdInstalled = True

  def nameserverOk(self):

    from Isd.PyroUtils import get_nameserver

    nshost = self.sim.replica.name_server

    try:
      ns = get_nameserver(nshost)

    except Exception, msg:
      print 'ERROR: Nameserver was not found on %s' % nshost

      showWarning('Failure', 'PyRO nameserver is not found. Please start the nameserver first '+\
                  'with pyro-ns script.')
      return False

    return True

  def isRunning(self):

    if self.sampler:
      if self.sampler.isAlive() and not self.sampler.ishalted():
        return True
    return False

  def generateReport(self):

    from Isd import report
    from cambridge.isd.isd_project_template import modify_simulation
    from Isd import setup

    try:
      simulation = self.simulation

      if not self.simulation:
        manager = setup.SetupManager(self.sim, modify_simulation)
        simulation = manager.create_simulation()

      R = report.ISDReport(simulation, superimpose=True)
      R.render()
      R.write_pdf()

    except Exception, msg:
      print 'ERROR: Exception occurred when generating a report, ', msg
      raise

    self.busylock.release()

  def __create(self):

    from cambridge.isd.isd_project_template import modify_simulation
    from Isd.generate_project_template import generate_template
    from Isd import setup

    ## write project file in working_path/isd.py
    ## (should be called before manager.create_simulation to
    ##  write the correct sim.temp_path into project file)

    project_file = os.path.join(self.sim.working_path,
                                '%s.py' % self.sim.fileroot)
    generate_template(sim = self.sim,
                      filename = project_file)

    ## create XML data files in working_path/data

    self.sim.finalize()

    manager = setup.SetupManager(self.sim, modify_simulation)
    simulation = manager.create_simulation()

    return manager, simulation

  def createSimulation(self):

    try:
      self.__create()

    except Exception, msg:
      print 'ERROR: Exception occurred when creating a simulation, ', msg
      raise

    self.busylock.release()

  def launchSimulation(self):

    ## .grid, .manager are the attributes of the frame
    ## for debugging only

    self.halted.clear()

    try:

      ## create simulation and XML data files in working_path/data

      self.manager, self.simulation = self.__create()

      ## generate samples

      print 'Starting calculation...'

      self.grid = self.manager.create_grid(self.debug)
      self.sampler = self.manager.create_sampler(self.grid, self.debug)

      initial_states = self.manager.prepare_initial_states(self.sampler)

      self.sampler.generate_sequence(self.sim.replica.n_samples, initial_states)

      print 'Calculation started.'

    except Exception, msg:
      print 'ERROR: Exception occurred when starting a simulation, ', msg
      raise

    self.busylock.release()

  def stopSimulation(self):

    cleanup(self.simulation, self.sampler)

    self.updateSettings()

    self.halted.set()

    self.busylock.release()

  def cmdInfo(self):

    if self.sampler:

      info = '''
Information on current simulation:

Total number of samples so far:    %d
Samples generated since last save: %d
Number of replicas:                %d

Results are written to %s'''

      args = (len(self.sampler.history['energies'])+1, len(self.sampler[0]),
              len(self.sampler), self.sim.working_path)

      print info % args
      showInfo( 'ISD simulation info', info % args )

    else:

      history = self.getHistory()

      if not history:
        showWarning('Failure','No simulation results in "%s".' \
                  % self.sim.working_path)
        return

      info = '''
Previous simulation results exist. Information on previous simulation:

Total number of samples so far:    %d
Number of replicas:                %d

Results are written to %s'''

      args = (len(history['energies'])+1, len(history['states']),
              self.sim.working_path)

      print info % args
      showInfo( 'ISD simulation info', info % args )

  def cmdShow(self):

    if not self.sampler:
      return

    if self.sampler.history['states']:

      self.show_index = 0
      chain = self.simulation.posterior.get_polymer()
      chain.set_torsions(self.sampler.history['states'][self.show_index].torsion_angles,1)
      chain.show()

    else:
      print 'Show: No samples available yet.'

  def cmdEnergies(self):

    if not self.sampler:
      return

    if not len(self.sampler.history['energies']):
      print 'Show: No samples available yet.'
      return

    from Isd.gnuplot import plot
    from numpy import sum, clip

    _min=-1.e100
    _max=1.e100

    self.show_index = 0

    E_replica = sum(self.sampler.history['energies'],1)
    plot(clip(E_replica[self.show_index:], _min, _max))

  def cmdRates(self, index=None):

    if not self.sampler:
      return

    from Isd.gnuplot import plot

    r = self.sampler.rates()

    if not len(r)>0:
      print 'Not enough samples yet.'
      return

    if index is None:
      plot(r)
    else:
      print 'Rate of exchange for replica %d: %.2f' % (index, r[index])

  def cmdSchedule(self):

    if not self.sampler:
      return

    from Isd.gnuplot import plot

    hs = self.sampler.heatbaths

    L = [h.replica_parameters['lambda'] for h in hs]
    Q = [h.replica_parameters['q'] for h in hs]

    print 'Lambda: \n', L
    print 'Q: \n', Q

    plot(L)
    plot(Q)

  def cmdSave(self):

    if not self.sampler:
      return

    self.sampler._dump(force=True)

  def cmdReport(self):

    if not self.getHistory():
      showWarning('Failure','No simulation results in "%s".' \
                  % self.sim.working_path)
      return

    if self.busylock.locked():
      showWarning('Report failure','ISD process is busy... Try later.')
      return False
    else:
      self.busylock.acquire()

    from threading import Thread
    t = Thread(target = self.generateReport)
    t.start()

    self.updateButtons()
    self.after(CHECK_INTERVAL, self.busylockReleased)

  def cmdWrite(self):

    if not self.settingsOk():
      return

    self.saveSettings()

    if self.busylock.locked():
      showWarning('Write ISD project failure','ISD process is busy... Try later.')
      return False
    else:
      self.busylock.acquire()

    from threading import Thread
    t = Thread(target = self.createSimulation)
    t.start()

    self.updateButtons()
    self.after(CHECK_INTERVAL, self.busylockReleased)

  def cmdRun(self):

    if not self.settingsOk():
      return

    self.saveSettings()

    if self.sim.replica.communication == 'pyro':
      if not self.nameserverOk():
        return False

    if self.busylock.locked():
      showWarning('Run ISD failure','ISD process is busy... Try later.')
      return False
    else:
      self.busylock.acquire()

    from threading import Thread
    t = Thread(target = self.launchSimulation)
    t.start()

    self.updateButtons()
    self.after(CHECK_INTERVAL, self.busylockReleased)
    self.after(CHECK_INTERVAL, self.simulationComplete)

  def cmdStop(self, force = False):

    ## called from:
    ## a) .closeProject()
    ## b) .simulationComplete with force=True
    ##    (when simulation complete is detected)
    ## c) when "Stop ISD!" button is pressed

    if self.isRunning() or force:

      if not allRemoteHostsReady(self.sampler):
        showWarning('Failure','Not all of remote hosts are initialised. '+\
                    'Please wait for several seconds and try again. '+\
                    'If the network is slow this may take a while.')
        return False

      if self.busylock.locked():
        showWarning('Stop ISD failure','ISD process is busy. '+\
                    'Please wait for several seconds and try again.')
        return False
      else:
        self.busylock.acquire()

      from threading import Thread
      t = Thread(target = self.stopSimulation)
      t.start()

      self.updateButtons()
      self.after(CHECK_INTERVAL, self.busylockReleased)

      return True

  def dataDetailMatrixEditRules(self, obj, row, col):
    # Disallow edit for immutable rows from the lower data set table

    if row < 5:
      return False
    else:
      return True

  def dataMatrixEditRules(self, obj, row, col):
    # Disallow edit of committed data sets
    from Isd import setup

    if obj.specsObj:
      return False
    elif (obj.format == setup.CCPN) and (col==5):
      return False
    else:
      return True

  def removeSelected(self):
    # Remove selected data sets from table and simulation if committed

    for datum in self.dataMatrix.currentObjects:
      self.dataSets.remove(datum)

    self.dataSet = None
    self.sim.data_sets = [ds.specsObj for ds in self.dataSets if ds.specsObj]
    self.updateExpData()

  def commitDataToggle(self, dataSet):
    from Isd import setup

    warnMsg = 'Cannot commit data set:'

    if not dataSet.key:
      showWarning('Failure','%s No key specified' % warnMsg)
      return

    if dataSet.format != setup.CCPN:
      if not file:
        showWarning('Failure','%s No file selected' % warnMsg)
        return

    if not dataSet.name:
      dataSet.name = dataSet.key

    spec = setup.setup_data(dataSet.dataType, dataSet.file, dataSet.name, dataSet.key, dataSet.format)
    dataSet.specsObj = spec

    self.sim.data_sets = [ds.specsObj for ds in self.dataSets if ds.specsObj]
    self.updateExpData()

  def getDataName(self, dataSet):

    self.dataNameEntry.set(dataSet.name)

  def setDataName(self, event): ## TODO when <None>

    name = self.dataNameEntry.get() or self.dataSet.key # key is default

    self.dataSet.name = name
    self.updateExpData()

  def getDataFormat(self, dataSet):
    from Isd import setup

    #dataType, file, format, key, name, specsObj = dataSet
    names = [setup.CCPN, setup.XML, setup.TBL]
    index = names.index(dataSet.format)
    self.dataFormatPulldown.setup(names,index)

  def setDataFormat(self, index, name=None):
    from Isd import setup

    value = self.dataFormatPulldown.getSelected()
    self.dataSet.format = value

    if value == setup.CCPN:
      self.dataSet.fileName = ''

    self.updateExpData()

  def getDataKey(self, dataSet):
    from Isd import CCPNReader
    from Isd import setup

    if dataSet.format == setup.CCPN:
      widget = self.dataKeyPulldown
      names  = []
      index  = -1

      key        = dataSet.key
      dataType   = dataSet.dataType
      nmrProject = self.ccpnProject.currentNmrProject

      if dataType == setup.NOESY:
        # These are simply distance restraints with peak intensity information (constraint.origData)
      
#         peakLists = getThroughSpacePeakLists(self.ccpnProject)

#         i = 0
#         for peakList in peakLists:
#           spectrum = peakList.dataSource
#           key0 = CCPNReader.getObjectKeyString(peakList)
#           if key0 == key:
#             index = i

#           names.append('%s - %s:%s:%d' % (key0,spectrum.experiment.name,spectrum.name,peakList.serial))
#           i += 1

        i = 0
        for constraintSet in nmrProject.sortedNmrConstraintStores():
          for constraintList in constraintSet.findAllConstraintLists(className='DistanceConstraintList'):
            key0 = CCPNReader.getObjectKeyString(constraintList)
            if key0 == key:
              index = i

            names.append('%s - %s' % (key0, constraintList.name))
            i += 1

      elif dataType == setup.DISTANCE:

        i = 0
        for constraintSet in nmrProject.sortedNmrConstraintStores():
          for constraintList in constraintSet.findAllConstraintLists(className='DistanceConstraintList'):
            key0 = CCPNReader.getObjectKeyString(constraintList)
            if key0 == key:
              index = i

            names.append('%s - %s' % (key0, constraintList.name))
            i += 1

      elif dataType == setup.JCOUPLING:

        i = 0
        for constraintSet in nmrProject.sortedNmrConstraintStores():
          for constraintList in constraintSet.findAllConstraintLists(className='JCouplingConstraintList'):
            key0 = CCPNReader.getObjectKeyString(constraintList)
            if key0 == key:
              index = i

            names.append('%s - %s' % (key0, constraintList.name))
            i += 1

      elif dataType == setup.RDC:

        i = 0
        for constraintSet in nmrProject.sortedNmrConstraintStores():
          for constraintList in constraintSet.findAllConstraintLists(className='RdcConstraintList'):
            key0 = CCPNReader.getObjectKeyString(constraintList)
            if key0 == key:
              index = i

            names.append('%s - %s' % (key0,  constraintList.name))
            i += 1

      elif dataType == setup.DIHEDRAL:

        i = 0
        for constraintSet in nmrProject.sortedNmrConstraintStores():
          for constraintList in constraintSet.findAllConstraintLists(className='DihedralConstraintList'):
            key0 = CCPNReader.getObjectKeyString(constraintList)
            if key0 == key:
              index = i

            names.append('%s - %s' % (key0, constraintList.name))
            i += 1

      elif dataType == setup.HBOND:

        i = 0
        for constraintSet in nmrProject.sortedNmrConstraintStores():
          for constraintList in constraintSet.findAllConstraintLists(className='HBondConstraintList'):
            key0 = CCPNReader.getObjectKeyString(constraintList)
            if key0 == key:
              index = i

            names.append('%s - %s' % (key0, constraintList.name))
            i += 1

      elif dataType == setup.DISULFIDE:
        pass # TBD: More tricky from CCPN, must fetch from molSystem's bonds...

      widget.setup(names,index)

    else:
      widget = self.dataKeyEntry
      widget.set(dataSet.key)

    # Bit of a hack because widgets are normally set at construction per column
    # Now setting it per row, and hence on-the-fly
    self.table.editWidget = widget

  def setDataKey(self, event, null=None):
    from Isd import setup

    dataSet = self.dataSet
    if self.dataSet.format == setup.CCPN:
      widget = self.dataKeyPulldown

      value  = widget.getSelected()

      if not value == '<None>':
        texts = value.split()
        dataSet.key  = texts[0]
        dataSet.name = texts[0]

    else:
      widget = self.dataKeyEntry
      key = widget.get()

      if key:
        dataSet.key = key

      if not dataSet.name:
        dataSet.name = key

    self.updateExpData()


  def getDataFile(self, dataSet):

    file_types = [ FileType("All", ["*"]), ]

    value = dataSet.file
    if not os.path.exists(value):
      value = None

    popup = FileSelectPopup(self, file_types, file=value, dismiss_text='Cancel')
    value = popup.getFile()
    if value:
      self.dataSet.file = value

      if not dataSet.key:
        text = os.path.split(value)[-1]
        self.dataSet.key = text.split('.')[0]

      self.updateExpData()

  def selectExpData(self, obj, row, col, table):

    self.table   = table
    self.dataSet = obj
    self.updateExpDataSet()


  def selectRowObj(self, obj, row, col, table):

    self.table = table
    self.rowObj = obj

  def addNOE(self):
    from Isd import setup
    # dataType, name, format, key, file, specsObj

    dataSet = SetupDataSet(setup.NOESY,'',setup.CCPN,'','',None)

    self.dataSets.append(dataSet)
    self.updateExpData()

  def addDistance(self):
    from Isd import setup

    dataSet = SetupDataSet(setup.DISTANCE,'',setup.CCPN,'','',None)
    self.dataSets.append(dataSet)
    self.updateExpData()

  def addJCoupling(self):
    from Isd import setup

    dataSet = SetupDataSet(setup.JCOUPLING,'',setup.CCPN,'','',None)
    self.dataSets.append(dataSet)
    self.updateExpData()

  def addRDC(self):
    from Isd import setup

    dataSet = SetupDataSet(setup.RDC,'',setup.CCPN,'','',None)
    self.dataSets.append(dataSet)
    self.updateExpData()

  def addDihedral(self):
    from Isd import setup

    dataSet = SetupDataSet(setup.DIHEDRAL,'',setup.CCPN,'','',None)
    self.dataSets.append(dataSet)
    self.updateExpData()

  def addHBond(self):
    from Isd import setup

    dataSet = SetupDataSet(setup.HBOND,'',setup.CCPN,'','',None)
    self.dataSets.append(dataSet)
    self.updateExpData()

  def addDisulfide(self):
    from Isd import setup

    dataSet = SetupDataSet(setup.DISULFIDE,'',setup.XML,'','',None)
    self.dataSets.append(dataSet)
    self.updateExpData()

  def getValue(self, rowObj):

    # get correct object, widget and get/set functions
    obj, attrName, widget, getter, setter = rowObj

    # Bit of a hack because widgets are normally set at construction per column
    # Now setting it per row, and hence on-the-fly
    self.table.editWidget = widget

    # get current value & setup widget
    getter(widget, obj, attrName)


  def setValue(self, event, null=None):
    # null is for pulldown menu callbacks passing name,index

    # get correct widget and get/set functions
    obj, attrName, widget, getter, setter = self.rowObj

    # set and check the appropriate parameter value from current edit widget

    # no setter for boolean toogles - the getter will have done all the toggling already
    # no setter for file selects - the file popup gets and sets the value and cannot be interrupted
    if setter:
      setter(widget, obj, attrName)

    self.table.refreshFunc() # Update relevant table (ScrolledMatrix)
    self.updateNmrCalc()

  def getString(self, widget, obj, attrName):

    value = getattr(obj, attrName)
    widget.set(value)

  def setString(self, widget, obj, attrName):

    value = widget.get().strip()

    if value:
      setattr(obj, attrName, value)
      self.updateNmrCalc()

  def getDirectory(self, widget, obj, attrName):
    # No widget because table placement is irrelevent and setting process not interruptable

    value = getattr(obj, attrName)
    if not os.path.exists(value):
      value = None

    popup = FileSelectPopup(self, directory=value, show_file=False)
    value = popup.getDirectory()
    if value:
      setattr(obj, attrName, value)

      self.table.refreshFunc()

  def getFile(self, widget, obj, attrName):
    # No widget because table placement is irrelevent and setting process not interruptable

    file_types = [ FileType("All", ["*"]), ]

    value = getattr(obj, attrName)
    if not os.path.exists(value):
      value = None

    popup = FileSelectPopup(self, file_types, file=value, dismiss_text='Cancel')
    value = popup.getFile()
    if value:
      setattr(obj, attrName, value)

      self.table.refreshFunc()

  def getXmlFile(self, widget, obj, attrName):
    # No widget because table placement is irrelevent and setting process not interruptable

    file_types = [ FileType("XML", ["*.xml"]), FileType("All", ["*"]) ]

    value = getattr(obj, attrName)
    if not os.path.exists(value):
      value = None

    popup = FileSelectPopup(self, file_types, file=value,dismiss_text='Cancel')
    value = popup.getFile()
    if value:
      setattr(obj, attrName, value)

      self.table.refreshFunc()

  def getInt(self, widget, obj, attrName):

    value = getattr(obj, attrName)
    widget.set(value)


  def setInt(self, widget, obj, attrName):

    value = widget.get() or 0
    setattr(obj, attrName, value)


  def setNonNegInt(self, widget, obj, attrName):

    value = widget.get() or 0

    if value >= 0:
      setattr(obj, attrName, value)
      self.updateNmrCalc()

  def setPosInt(self, widget, obj, attrName):

    value = widget.get() or 0

    if value > 0:
      setattr(obj, attrName, value)
      self.updateNmrCalc()


  def getFloat(self, widget, obj, attrName):

    value = getattr(obj, attrName)
    widget.set(value)


  def setFloat(self, widget, obj, attrName):

    value = widget.get() or 0.0
    setattr(obj, attrName, value)
    self.updateNmrCalc()

  def setFraction(self, widget, obj, attrName):

    value = min(1.0,max(0.0, widget.get() or 0.0))
    setattr(obj, attrName, value)
    self.updateNmrCalc()

  def setPosFloat(self, widget, obj, attrName):

    value = widget.get() or 0.0

    if value > 0.0:
      setattr(obj, attrName, value)
      self.updateNmrCalc()

  ######################################################

  def getCcpnProjectName(self):

    if self.ccpnProject:
      return self.ccpnProject.name

    else:
      return '<None>'

  def getCcpnProjectPath(self):

    if self.ccpnProject:
      return self.ccpnProject.findFirstRepository(name='userData').url.path

    else:
      return '<None>'

  def setNReplicas(self, widget, obj, attrName):

    value = widget.get() or 0

    if value > 0:
      #setattr(obj, , value)
      self.sim.replica.n_replicas = value
      self.sim.replica.weight_schedule.last = value
      self.sim.replica.prior_schedule.last = value
      self.updateAll()

  def getNamingSystems(self, widget, obj, attrName):
    from Isd import setup

    value = getattr(obj, attrName)
    names = [setup.IUPAC,setup.CNS]

    if value in names:
      index = names.index(value)
    else:
     index = 0

    widget.setup(names,index)

  """
  def getNmrProjectName(self, widget, obj, attrName):

    value = getattr(obj, attrName)
    names = []
    index = -1

    if self.ccpnProject:
      nmrProjects = self.ccpnProject.sortedNmrProjects()

      if nmrProjects:
        nmrProject = self.ccpnProject.findFirstNmrProject(name=value)

        names = [np.name for np in nmrProjects]

        if nmrProject:
          index = nmrProjects.index(nmrProject)
        else:
          index = 0

    widget.setup(names,index)

  def setNmrProjectName(self, widget, obj, attrName):

     value = widget.getSelected()

     if not value == '<None>':
       setattr(obj, attrName, value)
       self.chooseNmrProject(name=value)
  """

  def getMolSystemCode(self, widget, obj, attrName):

    value = getattr(obj, attrName)
    names = []
    index = -1

    if self.ccpnProject:
      molSystems = self.ccpnProject.sortedMolSystems()

      if molSystems:
        molSystem = self.ccpnProject.findFirstMolSystem(code=value)

        names = [ms.code for ms in molSystems]

        if molSystem:
          index = molSystems.index(molSystem)
        else:
          index = 0

    widget.setup(names,index)

  def getSaveStructures(self, widget, obj, attrName):
    from Isd import setup

    value = getattr(obj, attrName)
    names = [setup.REPRESENTATIVE,setup.MOST_PROBABLE]
    index = -1

    widget.setup(names,index)

  def getCcpnMoleculeKey(self, widget, obj, attrName):

    value = getattr(obj, attrName)
    names = []
    index = -1

    if self.ccpnProject:
      molSystems = self.ccpnProject.sortedMolSystems()

      if molSystems:
        code = self.sim.ccpn.export.molecular_system_name
        chains = []

        chain = None
        if code:
          molSystem = self.ccpnProject.findFirstMolSystem(code=code)

          if molSystem and molSystem.chains:
            chain = molSystem.findFirstChain()

        for molSystem0 in molSystems:
          chains.extend(molSystem0.sortedChains())

        if not chain:
          chain = chains[0]

        names = ['%s|%s' % (ch.molSystem.code,ch.code) for ch in chains]
        index = chains.index(chain)

    widget.setup(names,index)

  def setCcpnMoleculeKey(self, widget, obj, attrName):

     value = widget.getSelected()

     if not value == '<None>':

       molSystemCode, chainCode = string.split(value, '|')

       self.sim.ccpn.export.molecular_system_name = molSystemCode
       self.sim.molecule.key = value
       self.updateNmrCalc()


  def getSeqFileFormat(self, widget, obj, attrName):
    from Isd import setup

    value = getattr(obj, attrName)
    names = [setup.PDB, setup.SEQ, setup.XML, setup.CCPN]
    index = -1

    if value in names:
      index = names.index(value)
    else:
      index = 0

    widget.setup(names,index)

  def getStartConf(self, widget, obj, attrName):
    from Isd import setup

    value = getattr(obj, attrName)
    names = [setup.EXTENDED, setup.RANDOM, setup.AS_IS]
    index = -1

    if value in names:
      index = names.index(value)
    else:
      index = 0

    widget.setup(names,index)

  def getReplicaComm(self, widget, obj, attrName):

    value = getattr(obj, attrName)
    names = ['shared','pyro']
    index = -1

    if value in names:
      index = names.index(value)
    else:
      index = 0

    widget.setup(names,index)


  def setPulldownMenu(self, widget, obj, attrName):

    value = widget.getSelected()

    if not value == '<None>':
      setattr(obj, attrName, value)
      self.updateNmrCalc()

  def toggleBoolean(self, widget, obj, attrName):

    value = getattr(obj, attrName)
    setattr(obj, attrName, not value)

    self.table.refreshFunc()
    self.updateNmrCalc()

  def getMultiStrings(self, widget, obj, attrName):

    values = getattr(obj, attrName)
    widget.set(values=values, options=[])

  def getMultiOptStrings(self, widget, obj, attrName):

    dict = getattr(obj, attrName)

    options = dict.keys()
    values  = [dict.get(x) for x in options]

    widget.set(values=values, options=options)

  def textExportResidues(self, values):

    if values:

      blocks = [ [ values[0], None] ]

      x_old = values[0]

      for x in values[1:]:

        if x-x_old == 1:
          blocks[-1][1] = x
        else:
          blocks.append( [x, None] )

        x_old = x

      text_blocks = []
      for b in blocks:
        if b[1] is None:
          text_blocks.append( str(b[0]) )
        else:
          text_blocks.append( '%s-%s' % (b[0],b[1]))

      text = ','.join(text_blocks)

    else:
      text = ''

    return text

  def listExportResidues(self, text):

    if text:

      blocks = [x for x in text.split(',') if x is not '']

      values = []

      for block in blocks:

        if len(block.split('-'))>1:
          r1,r2 = [int(x) for x in block.split('-')]
          values += range(r1,r2+1)
        else:
          values.append(int(block))

    else:
      values = None

    return values

  def getExportResidues(self, widget, obj, attrName):

    values = getattr(obj, attrName) or []

    #text = ','.join([str(x) for x in values])

    text = self.textExportResidues(values)

    widget.set(text)

  def setExportResidues(self, widget, obj, attrName):

    ## takes '1-2,15,16,30-90' as input

    text = widget.get()

    #values = [int(x) for x in text.split(',') if x is not ''] or None

    try:
      values = self.listExportResidues(text)

    except Exception:
      print 'ERROR:', msg
      showError('Wrong input format',
                'Specify export residues by enumerating their numbers separated by commas and hyphens, e.g. 1,3,10-20',
                parent=self)
      return

    setattr(obj, attrName, values)
    self.updateNmrCalc()

  def setMultiStrings(self, widget, obj, attrName):

    values = widget.get()

    setattr(obj, attrName, values)

    self.updateNmrCalc()
    self.table.keyPressEscape()

  def setMultiOptStrings(self, widget, obj, attrName):

    dict = {}
    for i in range(len(widget.values)):
      value = widget.values[i]

      if value is not None:
        dict[widget.options[i]] = widget.values[i]

    self.updateNmrCalc()
    setattr(obj, attrName, dict)

    self.table.keyPressEscape()

  def textHostList(self, host_list):

    if host_list:

      hosts = []
      ns = []

      for h in host_list:
        if not h in hosts:
          hosts.append(h)
          ns.append(1)
        else:
          i = hosts.index(h)
          ns[i] += 1

      text_blocks = []

      for i in range(len(hosts)):
        if ns[i] > 1:
          text_blocks.append('%s %s' %(hosts[i], ns[i]))
        else:
          text_blocks.append( str(hosts[i]) )

      text = ', '.join(text_blocks)

    else:
      text = ''

    return text

  def listHostList(self, text):

    if text:

      blocks = [x for x in text.split(',') if x is not '']

      host_list = []

      for block in blocks:

        if len(block.split())>1:
          host, n = block.split()
          host_list += int(n)*[host.strip()]
        else:
          host_list.append(block.strip())

    else:
      host_list = None

    return host_list

  def getHostList(self, widget, obj, attrName):

    values = getattr(obj, attrName) or []

    text = self.textHostList(values)

    widget.set(text)

  def setHostList(self, widget, obj, attrName):

    ## takes 'localhost 30, localhost' as input

    text = widget.get()

    try:
      values = self.listHostList(text)

    except Exception, msg:
      print 'ERROR:', msg
      showError('Wrong input format',
                'Specify host list by their names or IP addresses '+\
                'separated by commas. Number of processes on each host '+\
                'can be specified after the name of the host. '+\
                'E.g. "localhost 10, hawk 5, localhost, duck" will create 11 processes '+\
                'on the localhost, 5 on hawk, and one process on duck',
                parent=self)
      return

    setattr(obj, attrName, values)
    self.updateNmrCalc()

  def getTempPaths(self, widget, obj, attrName):

    dict = getattr(obj, attrName)

    options = self.sim.replica.host_list
    values  = [dict.get(x) for x in options]

    widget.set(values=values, options=options)

  def simulationComplete(self):

    if self.sampler:
      if self.sampler.ishalted() and \
             (not self.busylock.locked() and not self.halted.isSet()):
        print 'Simulation complete.'
        print 'DEBUG!!! cmdStop is being called from simulationComplete'
        self.cmdStop(force = True)
        self.updateAll()
        return

    self.after(CHECK_INTERVAL, self.simulationComplete)

  def busylockReleased(self):

    if not self.busylock.locked():
      self.updateAll()
      return

    self.after(CHECK_INTERVAL, self.busylockReleased)

  def updateAll(self, project=None):

    ## called by ExtendNmrGui.ApplicationPopup.initProject

    if not self.haveIsdInstalled:
      return

    if project and (project is not self.ccpnProject):
      self.ccpnProject = project
    
    self.updateRunPulldown()
    self.updateNames()
    self.updateButtons()
    self.updateGeneral()
    self.updateMolStructure()
    self.updateMonteCarlo()
    self.updateAnalyses()
    self.updateExpData()
    
  def updateNmrCalc(self):
  
    if self.run:
      isdToNmrCalcRun(self.sim, self.ccpnProject, self.run)

  def updateNames(self):
    self.ccpnProjectNameLabel.set(self.getCcpnProjectName())
    self.ccpnProjectPathLabel.set(self.getCcpnProjectPath())

  def updateButtons(self):

    if not self.ccpnProject or self.busylock.locked():

      for b in self.allWidgets: b.disable()

    else:
      self.saveButton.enable()
      self.deleteButton.enable()
      self.newButton.enable()

      for b in self.reportButtons: b.enable()
      for b in self.datasetButtons: b.enable()

      if self.sampler:
        for b in self.samplerButtons: b.enable()
      else:
        for b in self.samplerButtons: b.disable()

      if self.isRunning():
        self.isdRunPulldown.disable()
        self.writeButton.disable()
        self.runButton.disable()
        self.stopButton.enable()

      else:
        self.writeButton.enable()
        self.isdRunPulldown.enable()
        self.runButton.enable()
        self.stopButton.disable()

  def updateGeneral(self):

    table = self.generalMatrix

    textMatrix  = []
    objectList  = []
    colorMatrix = []

    if not self.sim:
      table.update(objectList=objectList, textMatrix=textMatrix)
      return

    # Parameter Name, Parameter value, Description text
    # Object to manipulate, attribute name to manipulate, widget, getter (from obj), setter (to obj)

    textMatrix.append(['Simulation Name', self.sim.name, 'The name of the simulation'])
    objectList.append([self.sim,'name', self.stringEntry, self.getString, self.setString])

    textMatrix.append(['ISD Working Path', self.sim.working_path,'The directory in which ISD saves the results.'])
    objectList.append([self.sim,'working_path', None, self.getDirectory, None])

    textMatrix.append(['Temporary Path',self.sim.base_temp_path,'A directory accessible from all machines.'])
    objectList.append([self.sim,'base_temp_path', None,self.getDirectory, None])

    textMatrix.append(['Share Temp Path?',self.sim.shared_temp_path,'Whether temp path is accessible from all machines (cf. hosts).'])
    objectList.append([self.sim,'shared_temp_path', None, self.toggleBoolean, None])

    textMatrix.append(['CNS Executable',self.sim.cns_executable,'CNS programe executable. Required only if sequence from a SEQ file.'])
    objectList.append([self.sim,'cns_executable', None, self.getFile, None])

    textMatrix.append(['File Root Tag',self.sim.fileroot,'Start of output files, e.g. "isd" gives isd_13.pdb etc.'])
    objectList.append([self.sim,'fileroot',self.stringEntry, self.getString, self.setString])

    textMatrix.append(['Temperature',self.sim.temperature,'System temperature in Kelvin'])
    objectList.append([self.sim,'temperature', self.floatEntry, self.getFloat, self.setPosFloat])

    textMatrix.append(['Naming System',self.sim.naming_system,'Atom naming system; for sequence, PDB and data files'])
    objectList.append([self.sim,'naming_system', self.pulldownMenu, self.getNamingSystems, self.setPulldownMenu])

    #textMatrix.append(['CCPN Project File',self.sim.ccpn.project_filename,'Full path of the import/export CCPN project'])
    #objectList.append([self.sim.ccpn,'project_filename',None, self.getXmlFile, self.setFile])

    #if not self.sim.ccpn.nmr_project_name:
    #  self.sim.ccpn.nmr_project_name = self.nmrProject.name

    #textMatrix.append(['CCPN NMR Project ID',self.sim.ccpn.nmr_project_name,'Name (key) of CCPN "NmrProject" to store exported data'])
    #objectList.append([self.sim.ccpn,'nmr_project_name',
    #                   self.pulldownMenu, self.getNmrProjectName, self.setPulldownMenu])

    # TBD better defaults for below: guarantee unique
    if not self.sim.ccpn.export.project.key:
      self.sim.ccpn.export.project.key = 'defaultISD'

    #textMatrix.append(['CCPN ISD Storage Key',self.sim.ccpn.export.project.key,'Unique key of this ISD project under CCPN as an NmrCalc.run'])
    #objectList.append([self.sim.ccpn.export.project,
    #                   'key', self.stringEntry, self.getString, self.setString])

    ## NOTE: .sim.ccpn.export.project.enabled is not used anywhere
    ##       (it matters when run ISD from command line)
    ## textMatrix.append(['CCPN Save Project?',self.sim.ccpn.export.project.enabled,'Whether to export the simulation settings to CCPN'])
    ## objectList.append([self.sim.ccpn.export.project,'enabled', None, self.toggleBoolean, None])

    textMatrix.append(['CCPN Save Ensemble?',self.sim.ccpn.export.ensemble.enabled,'Whether to export the conformational ensemble to CCPN project'])
    objectList.append([self.sim.ccpn.export.ensemble,'enabled', None, self.toggleBoolean, None])

    if self.debug:
      self.sim.use_xterm = True

    textMatrix.append(['Use xterm?',self.sim.use_xterm,'For debugging only; requires X forwarding to be activated.'])
    objectList.append([self.sim,'use_xterm', None, self.toggleBoolean, None])

    textMatrix.append(['Debug?',self.debug,'For debugging only.'])
    objectList.append([self,'debug', None, self.toggleBoolean, None])

    # Add first row number column to the text maxtrix
    # Also setup blank colour matrix
    i = 0
    for datum in textMatrix:
      i += 1
      datum.insert(0,i)
      colorMatrix.append([None] * 4)

    # Add colour rules here

    from Isd import setup

    if self.sim.working_path == './':
      colorMatrix[1][2] = redColor
      colorMatrix[9][2] = redColor

    if not os.path.exists(os.path.expanduser(self.sim.cns_executable)):
      colorMatrix[4][2] = redColor

    if not self.sim.ccpn.nmr_project_name:
      colorMatrix[8][2] = redColor

    if not self.sim.ccpn.export.project.key:
      colorMatrix[9][2] = redColor

    # Finally update the ScrolledMatrix table
    table.update(objectList=objectList, textMatrix=textMatrix, colorMatrix=colorMatrix)

  def updateMolStructure(self):

    table = self.molMatrix

    textMatrix  = []
    objectList  = []
    colorMatrix = []

    if not self.sim:
      table.update(objectList=objectList, textMatrix=textMatrix)
      return

    textMatrix.append(['Structures to Save',self.sim.pdb_files.ensemble,'Set whether representative or most probable structures are saved'])
    objectList.append([self.sim.pdb_files,'ensemble', self.pulldownMenu,self.getSaveStructures,self.setPulldownMenu])

    textMatrix.append(['Saved Ensemble Size',self.sim.pdb_files.n,'Number of ensemble members that are stored as PDB file'])
    objectList.append([self.sim.pdb_files,'n',self.intEntry,self.getInt,self.setNonNegInt])

    textMatrix.append(['Exported Residues',self.textExportResidues(self.sim.pdb_files.residue_list),
                       'Residues to save in PDB files and use for quality (WhatIf etc) analysis.'])
    objectList.append([self.sim.pdb_files,'residue_list',self.stringEntry,self.getExportResidues,self.setExportResidues])

    if not self.sim.ccpn.export.molecular_system_name:
      if self.ccpnProject:
        molSystems = self.ccpnProject.sortedMolSystems()
        if len(molSystems) == 1:
          self.sim.ccpn.export.molecular_system_name = molSystems[0].code

    textMatrix.append(['CCPN Molecular System Name',self.sim.ccpn.export.molecular_system_name,'Code (key) of CCPN "MolSystem" to hold structures.'])
    objectList.append([self.sim.ccpn.export,
                       'molecular_system_name',
                       self.pulldownMenu,
                       self.getMolSystemCode,
                       self.setPulldownMenu])

    if not self.sim.molecule.key:
      if self.ccpnProject:
        molSystemCode = self.sim.ccpn.export.molecular_system_name
        molSystem = self.ccpnProject.findFirstMolSystem(code=molSystemCode)
        if molSystem:
          chains = molSystem.sortedChains()
          if len(chains) == 1:
            chainCode = chains[0].code
            self.sim.molecule.key = string.join((molSystemCode, chainCode), '|')

    textMatrix.append(['CCPN Molecule Key',self.sim.molecule.key ,'Get sequence from named molecule in CCPN project'])
    objectList.append([self.sim.molecule,'key',
                       self.pulldownMenu,
                       self.getCcpnMoleculeKey,
                       self.setCcpnMoleculeKey])

    textMatrix.append(['Sequence File',self.sim.molecule.filename,'PDB or text file containing sequence of 3-letter codes.'])
    objectList.append([self.sim.molecule,'filename',None,self.getFile,None])

    textMatrix.append(['Sequence File Format',self.sim.molecule.format,'Format of sequence file, PDB, SEQ or XML.'])
    objectList.append([self.sim.molecule,'format', self.pulldownMenu,self.getSeqFileFormat,self.setPulldownMenu])

    textMatrix.append(['Residue Start Number',self.sim.molecule.first_residue_number,'Number of the first resdiue; only set if using SEQ file'])
    objectList.append([self.sim.molecule,'first_residue_number',self.intEntry,self.getInt,self.setInt])

    textMatrix.append(['Starting Conformation',self.sim.molecule.initial_conformation,'Initial conformation used for the calculation'])
    objectList.append([self.sim.molecule,'initial_conformation', self.pulldownMenu,self.getStartConf,self.setPulldownMenu])

    textMatrix.append(['Nonbonded Term Excludes H?',self.sim.molecule.exclude_hydrogens,'Whether to exclude hydrogens in the calculation of nonbonded energy.'])
    objectList.append([self.sim.molecule,'exclude_hydrogens', None,self.toggleBoolean, None])

    # Add first row number column to the text maxtrix
    i = 0
    for datum in textMatrix:
      i += 1
      datum.insert(0,i)
      colorMatrix.append([None] * 4)

    # Add colour rules here

    from Isd import setup

    if self.sim.molecule.format == setup.CCPN:

      if not self.sim.ccpn.export.molecular_system_name:
        colorMatrix[3][2] = redColor
        colorMatrix[6][2] = redColor

      if not self.sim.molecule.key:
        colorMatrix[4][2] = redColor
        colorMatrix[6][2] = redColor

    else:
      if not os.path.exists(self.sim.molecule.filename):
        colorMatrix[5][2] = redColor
        colorMatrix[6][2] = redColor

    table.update(objectList=objectList, textMatrix=textMatrix, colorMatrix=colorMatrix)

  def updateMonteCarlo(self):

    table = self.monteCarloMatrix

    textMatrix = []
    objectList = []
    colorMatrix = []

    if not self.sim:
      table.update(objectList=objectList, textMatrix=textMatrix)
      return

    textMatrix.append(['Num. Replicas',self.sim.replica.n_replicas,'Number of replica structures in Monte Carlo simulation'])
    #objectList.append([self.sim.replica,'n_replicas',self.intEntry,self.getInt,self.setPosInt])
    objectList.append([self.sim.replica,'n_replicas',self.intEntry,self.getInt,self.setNReplicas])

    textMatrix.append(['Num. Samples',self.sim.replica.n_samples,' of replica exchange Monte Carlo samples'])
    objectList.append([self.sim.replica,'n_samples',self.intEntry,self.getInt,self.setPosInt])

    textMatrix.append(['Num. HMC Steps',self.sim.replica.hmc.steps,'Number of Hybrid Monte Carlo (HMC) steps per replica exchange step'])
    objectList.append([self.sim.replica.hmc,'steps',self.intEntry,self.getInt,self.setPosInt])

    textMatrix.append(['Num. MD Steps',self.sim.replica.hmc.md.steps,'Number of moleculear dynamics steps per HMC step'])
    objectList.append([self.sim.replica.hmc.md,'steps',self.intEntry,self.getInt,self.setPosInt])

    textMatrix.append(['HMC Time Step',self.sim.replica.hmc.stepsize,'Time step of the HMC leapfrog integration'])
    objectList.append([self.sim.replica.hmc,'stepsize',self.floatEntry,self.getFloat,self.setFloat])

    textMatrix.append(['Stepsize Adjustment',self.sim.replica.hmc.adjust_stepsize,'How many samples to adjust step size for; remaining have fixed step'])
    objectList.append([self.sim.replica.hmc,'adjust_stepsize',self.intEntry,self.getInt,self.setPosInt])

    textMatrix.append(['Communication Method',self.sim.replica.communication,
                       'Communication method between replicas. Either "shared" or "pyro"'])

    objectList.append([self.sim.replica,'communication',self.pulldownMenu,
                       self.getReplicaComm,self.setPulldownMenu])

    textMatrix.append(['Resume Override Params?',self.sim.replica.override_parameters,'After a restart, whether initial parameters override the last known ones.'])
    objectList.append([self.sim.replica,'override_parameters', None,self.toggleBoolean, None])

    textMatrix.append(['Save Interval',self.sim.replica.save_interval,'Interval in which posterior distribution states are saved.'])
    objectList.append([self.sim.replica,'save_interval',self.intEntry,self.getInt,self.setPosInt])

    textMatrix.append(['Do Full Save?',self.sim.replica.full_save,'Whether all ensembles are saved with the simulation'])
    objectList.append([self.sim.replica,'full_save', None,self.toggleBoolean, None])

    textMatrix.append(['Python Binary',self.sim.replica.python,'Name of the Python command used on remote machines'])
    objectList.append([self.sim.replica,'python',self.stringEntry, self.getString, self.setString])

    textMatrix.append(['Job Niceness',self.sim.replica.niceness,'Niceness of remote jobs: Lower value means higher scheduling priority'])
    objectList.append([self.sim.replica,'niceness',self.intEntry,self.getInt,self.setInt])

    textMatrix.append(['Backround Run?',self.sim.replica.background,'Whether to run without a terminal window (i.e. interactive Python shell)'])
    objectList.append([self.sim.replica,'background', None,self.toggleBoolean, None])

    textMatrix.append(['PyRo Name Server',self.sim.replica.name_server,'Name or IP address of the machine running the PyRo name server'])
    objectList.append([self.sim.replica,'name_server',self.stringEntry, self.getString, self.setString])

    #textMatrix.append(['Host List',self.sim.replica.host_list,'Names or IP addresses of machines used for the calculation.'])
    #objectList.append([self.sim.replica,'host_list',self.multiString,self.getMultiStrings,self.setMultiStrings])

    textMatrix.append(['Host List',self.textHostList(self.sim.replica.host_list),
                       'Names or IP addresses of machines (and the number of processes) used for the calculation.'])
    objectList.append([self.sim.replica,'host_list',self.stringEntry,self.getHostList,self.setHostList])

    textMatrix.append(['Temp Dir Paths',self.sim.replica.temp_paths or None,'If set, temp directories for each host, otherwise directory is shared'])
    objectList.append([self.sim.replica,'temp_paths',self.multiString,self.getTempPaths,self.setMultiOptStrings])

    textMatrix.append(['Weight Sched. Init',self.sim.replica.weight_schedule.initial,'Initial likelihood weight (similar to an inverse temp.); between 0 & 1'])
    objectList.append([self.sim.replica.weight_schedule,'initial',self.floatEntry,self.getFloat,self.setFraction])

    textMatrix.append(['Weight Sched. Final',self.sim.replica.weight_schedule.final,'Final likelihood weight (similar to an inverse temp.)'])
    objectList.append([self.sim.replica.weight_schedule,'final',self.floatEntry,self.getFloat,self.setFloat])

    textMatrix.append(['Weight Sched. Slope',self.sim.replica.weight_schedule.slope,'Slope of power law used to calculate intermediate weights'])
    objectList.append([self.sim.replica.weight_schedule,'slope',self.floatEntry,self.getFloat,self.setFloat])

    textMatrix.append(['Weight Sched. First',self.sim.replica.weight_schedule.first,'Replica to start weight changes at; preceding ones have initial weight'])
    objectList.append([self.sim.replica.weight_schedule,'first',self.intEntry,self.getInt,self.setPosInt])

    textMatrix.append(['Weight Sched. Last ',self.sim.replica.weight_schedule.last,'Replica to stop weight changes at; following ones have final weight'])
    #objectList.append([self.sim.replica.weight_schedule,'last',self.intEntry,self.getInt,self.setPosInt])
    objectList.append(None)

    textMatrix.append(['Prior Sched. Init',self.sim.replica.prior_schedule.initial,'Initial prior weight; between 0 & 1'])
    objectList.append([self.sim.replica.prior_schedule,'initial',self.floatEntry,self.getFloat,self.setFraction])

    textMatrix.append(['Prior Sched. Final',self.sim.replica.prior_schedule.final,'Final prior weight'])
    objectList.append([self.sim.replica.prior_schedule,'final',self.floatEntry,self.getFloat,self.setFloat])

    textMatrix.append(['Prior Sched. Slope',self.sim.replica.prior_schedule.slope,'Slope of power law used to calculate intermediate weights.'])
    objectList.append([self.sim.replica.prior_schedule,'slope',self.floatEntry,self.getFloat,self.setFloat])

    textMatrix.append(['Prior Sched. First',self.sim.replica.prior_schedule.first,'Peplica to start prior weight changes at; preceding ones have initial weight'])
    objectList.append([self.sim.replica.prior_schedule,'first',self.intEntry,self.getInt,self.setPosInt])

    textMatrix.append(['Prior Sched. Last',self.sim.replica.prior_schedule.last,'Replica to stop prior weight changes at; following ones have final weight'])
    #objectList.append([self.sim.replica.prior_schedule,'last',self.intEntry,self.getInt,self.setPosInt])
    objectList.append(None)

    # Add first row number column to the text maxtrix
    i = 0
    for datum in textMatrix:
      i += 1
      datum.insert(0,i)
      colorMatrix.append([None] * 4)

    if not os.path.exists(os.path.expanduser(self.sim.replica.python)):
      colorMatrix[10][2] = redColor

    if len(self.sim.replica.host_list) == 1 and self.sim.replica.host_list[0] == 'localhost':
      colorMatrix[14][2] = redColor

    table.update(objectList=objectList, textMatrix=textMatrix, colorMatrix=colorMatrix)

  def updateAnalyses(self):

    table = self.analysisMatrix

    textMatrix = []
    objectList = []
    colorMatrix = []

    if not self.sim:
      table.update(objectList=objectList, textMatrix=textMatrix)
      return

    textMatrix.append(['Burn-in Samples',self.sim.analysis.burnin,'Number of most recent samples analysed; prior (converging) samples are discarded'])
    objectList.append([self.sim.analysis,'burnin',self.intEntry,self.getInt,self.setPosInt])

    textMatrix.append(['Auto Report?',self.sim.analysis.report.auto,'Whether ISD creates a report whenever results are saved (cf. Save Interval).'])
    objectList.append([self.sim.analysis.report,'auto', None,self.toggleBoolean, None])

    textMatrix.append(['Keep LaTex Sources?',self.sim.analysis.report.keep_sources,'Whether to keep LaTex sources; in WORKING_PATH/analysis/data/'])
    objectList.append([self.sim.analysis.report,'keep_sources', None,self.toggleBoolean, None])

    textMatrix.append(['Structure Viewer',self.sim.analysis.pdb_viewer,'Program to display structures. Called by ISD "show" with PDB file argument'])
    objectList.append([self.sim.analysis,'pdb_viewer',self.stringEntry, self.getString, self.setString])

    textMatrix.append(['PDF to LaTex',self.sim.analysis.pdf_latex,'Command to execute PDF to LaTex converter'])
    objectList.append([self.sim.analysis,'pdf_latex',self.stringEntry, self.getString, self.setString])

    textMatrix.append(['EPS to PDF',self.sim.analysis.eps_to_pdf,'Command to execute EPS to PDF converter'])
    objectList.append([self.sim.analysis,'eps_to_pdf',self.stringEntry, self.getString, self.setString])

    textMatrix.append(['EPS to EPS',self.sim.analysis.eps_to_eps,'Command to execute eps2eps'])
    objectList.append([self.sim.analysis,'eps_to_eps',self.stringEntry, self.getString, self.setString])

    textMatrix.append(['WhatIf/WhatCheck Binary',self.sim.analysis.whatif.binary,'Command to execute the WhatIf program (or WhatCheck if below is True)'])
    objectList.append([self.sim.analysis.whatif,'binary',self.stringEntry, self.getString, self.setString])

    textMatrix.append(['Use WhatCheck?',self.sim.analysis.whatif.use_whatcheck,'Whether to run WhatCheck instead of WhatIf'])
    objectList.append([self.sim.analysis.whatif,'use_whatcheck', None,self.toggleBoolean, None])

    textMatrix.append(['Whatif Figures?',self.sim.analysis.whatif.show_traces,'Whether to generate WhatIf figures; scoring every sample'])
    objectList.append([self.sim.analysis.whatif,'show_traces', None,self.toggleBoolean, None])

    textMatrix.append(['WhatIf Quality Scores?',self.sim.analysis.whatif.enabled,'Whether to calculate quality scores with WhatIf'])
    objectList.append([self.sim.analysis.whatif,'enabled', None,self.toggleBoolean, None])

    textMatrix.append(['Prockeck Binary',self.sim.analysis.procheck.binary,'Command to execute the Procheck program'])
    objectList.append([self.sim.analysis.procheck,'binary',self.stringEntry, self.getString, self.setString])

    textMatrix.append(['Procheck Figures?',self.sim.analysis.procheck.show_traces,'Whether to generate Procheck figures; scoring every sample'])
    objectList.append([self.sim.analysis.procheck,'show_traces', None,self.toggleBoolean, None])

    textMatrix.append(['Procheck Quality Scores',self.sim.analysis.procheck.enabled,'Whether to calculate quality scores with Procheck'])
    objectList.append([self.sim.analysis.procheck,'enabled', None,self.toggleBoolean, None])

    textMatrix.append(['DSSP Binary',self.sim.analysis.dssp.binary,'Command to execute the DSSP program'])
    objectList.append([self.sim.analysis.dssp,'binary',self.stringEntry, self.getString, self.setString])

    textMatrix.append(['Use DSSP?',self.sim.analysis.dssp.enabled,'Whether to generate secondary structure assignments with DSSP'])
    objectList.append([self.sim.analysis.dssp,'enabled', None,self.toggleBoolean, None])

    textMatrix.append(['Gnuplot binary',self.sim.analysis.gnuplot,'Program to execute Gnuplot program'])
    objectList.append([self.sim.analysis,'gnuplot',self.stringEntry, self.getString, self.setString])

    #textMatrix.append(['Mas Samples',self.sim.analysis.max_samples,'For internal use; should not be modified'])
    #objectList.append([self.sim.analysis,'max_samples',self.intEntry,self.getInt,self.setPosInt

    # Add first row number column to the text maxtrix
    i = 0
    for datum in textMatrix:
      i += 1
      datum.insert(0,i)
      colorMatrix.append([None] * 4)

    if self.sim.analysis.whatif.enabled:
      if not os.path.exists(self.sim.analysis.whatif.binary):
        colorMatrix[7][2] = redColor

    if self.sim.analysis.procheck.enabled:
      if not os.path.exists(self.sim.analysis.procheck.binary):
        colorMatrix[11][2] = redColor

    if self.sim.analysis.dssp.enabled:
      if not os.path.exists(self.sim.analysis.dssp.binary):
        colorMatrix[14][2] = redColor

    table.update(objectList=objectList, textMatrix=textMatrix, colorMatrix=colorMatrix)

  def updateExpData(self):

    table = self.dataMatrix

    textMatrix = []
    objectList = []
    colorMatrix = []

    if not self.dataSets:
      table.update(objectList=objectList, textMatrix=textMatrix)
      self.updateExpDataSet()
      return

    self.dataSets.sort()

    for dataSet in self.dataSets:

      #['Data Type','Name','Commited','Format','Key','Filename']

      colors = [None] * 6

      if dataSet.specsObj:
        colors[2] = greenColor
        committed = 'Yes'
      else:
        colors[2] = redColor
        committed = 'No'

      textMatrix.append([dataSet.dataType, dataSet.name, committed, dataSet.format, dataSet.key, dataSet.file])
      objectList.append(dataSet)
      colorMatrix.append(colors)

    table.update(objectList=objectList, textMatrix=textMatrix, colorMatrix=colorMatrix)

    self.updateNmrCalc()
    self.updateExpDataSet()

  def updateExpDataSet(self):

    from Isd import setup

    table = self.dataDetailMatrix

    textMatrix = []
    objectList = []

    if not self.dataSet:
      table.update(objectList=objectList, textMatrix=textMatrix)
      return

    if self.dataSet:
      dataSet  = self.dataSet

      specs    = dataSet.specsObj
      dataType = dataSet.dataType

      textMatrix.append(['Data Type',dataType,'The type of data set to be analysed, e.g. distances or dihedrals'])
      objectList.append(None)

      textMatrix.append(['Data Name',dataSet.name,'Unique name tag for the data set'])
      objectList.append(None)

      textMatrix.append(['Data Format',dataSet.format,'The input format of the data set; XML, TBL or CCPN'])
      objectList.append(None)

      textMatrix.append(['Data Key',dataSet.key,'Unique identifier for data set; for CCPN this is the key path'])
      objectList.append(None)

      textMatrix.append(['Data File',dataSet.file,'File form which the data is loaded, if not from a CCPN project'])
      objectList.append(None)

      if specs:

        if dataType in (setup.NOESY,setup.DISTANCE):

          textMatrix.append(['Estimate Scale?',specs.theory.scale.update,'Whether to estimate the scale of the NOEs'])
          objectList.append([specs.theory.scale,'update',None,self.toggleBoolean, None])

          textMatrix.append(['Initial Error',specs.error_model.error.initial,'Initial error value (% relative error on distance scale)'])
          objectList.append([specs.error_model.error,'initial',self.floatEntry,self.getFloat,self.setFloat])

          textMatrix.append(['Estimate Error?',specs.error_model.error.update,'Whether to estimate the error of the data'])
          objectList.append([specs.error_model.error,'update',None,self.toggleBoolean, None])

        elif dataType == setup.JCOUPLING:

          textMatrix.append(['Estimate Coefficients',specs.theory.karplus_curve.update,'Whether to estimate coefficients in the Karplus relationship'])
          objectList.append([specs.theory.karplus_curve,'update',None,self.toggleBoolean, None])

          textMatrix.append(['Karplust Param. A',specs.theory.karplus_curve.A,'First Karplus coefficient'])
          objectList.append([specs.theory.karplus_curve,'A',self.floatEntry,self.getFloat,self.setFloat])

          textMatrix.append(['Karplust Param. B',specs.theory.karplus_curve.B,'Second Karplus coefficient'])
          objectList.append([specs.theory.karplus_curve,'B',self.floatEntry,self.getFloat,self.setFloat])

          textMatrix.append(['Karplust Param. C',specs.theory.karplus_curve.C,'Third Karplus coefficient'])
          objectList.append([specs.theory.karplus_curve,'C',self.floatEntry,self.getFloat,self.setFloat])

          textMatrix.append(['Initial Error',specs.error_model.error.initial,'Initial value of the error'])
          objectList.append([specs.error_model.error,'initial',self.floatEntry,self.getFloat,self.setFloat])

          textMatrix.append(['Estimate Error?',specs.error_model.error.update,'Whether to estimate the error of the data'])
          objectList.append([specs.error_model.error,'update',None,self.toggleBoolean, None])

        elif dataType == setup.RDC:

          textMatrix.append(['Estimate Tensor?',specs.theory.saupe_tensor.update,'Whether to estimate the tensor elements'])
          objectList.append([specs.theory.saupe_tensor,'update',None,self.toggleBoolean, None])

          textMatrix.append(['Tensor Element 1',specs.theory.saupe_tensor.s1,'First tensor element'])
          objectList.append([specs.theory.saupe_tensor,'s1',self.floatEntry,self.getFloat,self.setFloaat])

          textMatrix.append(['Tensor Element 2',specs.theory.saupe_tensor.s2,'Second tensor element'])
          objectList.append([specs.theory.saupe_tensor,'s2',self.floatEntry,self.getFloat,self.setFloaat])

          textMatrix.append(['Tensor Element 3',specs.theory.saupe_tensor.s3,'Third tensor element'])
          objectList.append([specs.theory.saupe_tensor,'s3',self.floatEntry,self.getFloat,self.setFloaat])

          textMatrix.append(['Tensor Element 4',specs.theory.saupe_tensor.s4,'Forth tensor element'])
          objectList.append([specs.theory.saupe_tensor,'s4',self.floatEntry,self.getFloat,self.setFloaat])

          textMatrix.append(['Tensor Element 5',specs.theory.saupe_tensor.s5,'Fifth tensor element'])
          objectList.append([specs.theory.saupe_tensor,'s5',self.floatEntry,self.getFloat,self.setFloaat])

          textMatrix.append(['Initial Error',specs.error_model.error.initial,'Initial error value'])
          objectList.append([specs.error_model.error,'initial',self.floatEntry,self.getFloat,self.setFloat])

          textMatrix.append(['Estimate Error?',specs.error_model.error.update,'Whether to estimate the error of the data'])
          objectList.append([specs.error_model.error,'update',None,self.toggleBoolean, None])

        elif dataType == setup.DIHEDRAL:

          # Below was missing in default setup
          textMatrix.append(['Initial Error',specs.error_model.error.initial,'Initial error value'])
          objectList.append([specs.error_model.error,'initial',self.floatEntry,self.getFloat,self.setFloat])

          textMatrix.append(['Estimate Error?',specs.error_model.error.update,'Whether to estimate the error of the data'])
          objectList.append([specs.error_model.error,'update',None,self.toggleBoolean, None])

        elif dataType == setup.HBOND:

          textMatrix.append(['Initial Error',specs.error_model.error.initial,'Initial error value'])
          objectList.append([specs.error_model.error,'initial',self.floatEntry,self.getFloat,self.setFloat])

          textMatrix.append(['Estimate Error?',specs.error_model.error.update,'Whether to estimate the error of the data'])
          objectList.append([specs.error_model.error,'update',None,self.toggleBoolean, None])

        elif dataType == setup.DISULFIDE:

          textMatrix.append(['Distance',specs.distance,'G-SG Distance in Angstrom; default: 2.02'])
          objectList.append([specs,'distance',self.floatEntry,self.getFloat,self.setPosFloat])

          textMatrix.append(['Initial Error',specs.error_model.error.initial,'Initial error value'])
          objectList.append([specs.error_model.error,'initial',self.floatEntry,self.getFloat,self.setFloat])

          textMatrix.append(['Estimate Error?',specs.error_model.error.update,'Whether to estimate the error of the data'])
          objectList.append([specs.error_model.error,'update',None,self.toggleBoolean, None])

    # Add first row number column to the text maxtrix
    i = 0
    for datum in textMatrix:
      i += 1
      datum.insert(0,i)

    table.update(objectList=objectList, textMatrix=textMatrix)


def allRemoteHostsReady(sampler):

    if sampler:

      service_id = sampler.heatbaths[0].SERVICE_ID
      servers = sampler.heatbaths[0].grid.servers[service_id]

      if len(servers) == len(sampler.heatbaths[0].grid.hosts):
        return True

      else:
        return False

def cleanup(simulation, sampler):

    from time import sleep
    from Isd.utils import SpinWheel

    wheel = SpinWheel()

    ## halt sampler

    sampler.halt()

    t_delta = 0.2
    t_timeout = SAMPLER_TIMEOUT

    t = 0.

    while not sampler.ishalted() and t < t_timeout:
        wheel.update('Halting sampler... ')
        sleep(t_delta)
        t += t_delta

    print

    ## halt grid

    g = sampler.heatbaths[0].grid
    g.terminate()

    t = 0.

    while not g.ishalted() and t < t_timeout:
        wheel.update('Halting grid... ')
        sleep(t_delta)
        t += t_delta

    print

    if g.ishalted():
      print 'ISD simulation is stopped!'
