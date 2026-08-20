
import subprocess
from multiprocessing import Process

import os, tempfile, glob

import Tkinter

from memops.universal import Io as uniIo

from memops.gui.DataEntry         import DataEntry
from memops.gui.FileSelect        import FileType
from memops.gui.FileSelectPopup   import FileSelectPopup
from memops.gui.TabbedFrame       import TabbedFrame
from memops.gui.PulldownList      import PulldownList
from memops.gui.Button            import Button
from memops.gui.ButtonList        import ButtonList, UtilityButtonList
from memops.gui.Entry             import Entry
from memops.gui.Label             import Label
from memops.gui.LabelDivider      import LabelDivider
from memops.gui.LabelFrame        import LabelFrame
from memops.gui.MessageReporter   import showInfo, showOkCancel, showWarning, showMulti
from memops.gui.ScrolledText      import ScrolledText
from memops.gui.ScrolledMatrix    import ScrolledMatrix 
from memops.gui.WebBrowser        import WebBrowser

from memops.gui.Frame             import Frame

from ccp.api.nmr                  import NmrCalc
from ccp.util.NmrCalc             import setRunParameter, setRunTextParameter, getRunTextParameter

from ccp.lib.StructureIo import makePdbFromStructure, getStructureFromFile

from grenoble.BlackledgeModule    import Io as moduleIo

MODULE        = 'MODULE'
MODULE_GREEN  = '#80FF80'
MODULE_RED    = '#FF8080'
MODULE_BLUE   = '#78B8F0'

INPUT       = 'input'
PROVISIONAL = 'provisional'

# NmrCalc data
MODULE_EXE              = 'MODULEExecutable'
STRUCTURE_DATA          = 'StructureEnsembleData'
RDC_CONSTRAINT_DATA     = 'RdcConstraintStoreData'
DIST_CONSTRAINT_DATA    = 'DistanceConstraintStoreData'
USER_DESCRIPTION_DATA   = 'UserDescriptionData'

# data to be sent to MODULE ?? 
DATA_TITLES = { MODULE_EXE:'MODULE Executable',
                STRUCTURE_DATA:'Structure',
                RDC_CONSTRAINT_DATA:'RDC Restraint List',
                DIST_CONSTRAINT_DATA:'Distance Restraint List',
                USER_DESCRIPTION_DATA:'User Description of MODULE Run',}


class BlackledgeModuleFrame(Frame):
  """ Frame for handling MODULE calculation.
  
  Note that the frame uses (or creates) an NmrCalcStore named 'BLACKLEDGE_MODULE'
  and linked to the current NmrProject.
  """

  def __init__(self, parent, project, closeButton=False, tempFiles=False,
               *args, **kw):

    ###########################################################################
    # INIT VARIABLES

    self.parent                   = parent
    self.project                  = project

    try:
      self.nmrProject               = (project.currentNmrProject 
                                      or project.newNmrProject(name='BLACKLEDGE_MODULE'))
    except:
      print '&&& Running MODULE popup from outside CCPN Analysis - debug only - no NmrCalc'
      self.nmrProject               = None

    if self.nmrProject:
      self.calcStore = project.findFirstNmrCalcStore(name=MODULE, nmrProject=self.nmrProject) or \
                       project.newNmrCalcStore(name=MODULE, nmrProject=self.nmrProject)

    else:
      self.calcStore = None

    self.run                          = None

    self.inputStructure               = None
    self.inputRdcConstraintList       = None
    self.inputDistanceConstraintList  = [None]
    self.inputUserDescriptionText     = None

    # path to the module executable
    modPath = subprocess.Popen(['which', 'module'], stdout=subprocess.PIPE).communicate()[0].strip()
    self.moduleExePath            = modPath or '    NB. MODULE executable not found    '

    self.waiting = False

    # for debug this could be False
    if tempFiles:
      self.useTempFiles   = True
    else:
      self.useTempFiles   = False

    # create temp files for MODULE
    if self.useTempFiles:
      self.moduleTempDir = tempfile.mkdtemp( prefix='MODULE-' )
    else:
      self.moduleTempDir = os.getcwd()

    #djo35# self.calcStore = self.resetCalcStore(name='BLACKLEDGE_MODULE')

    # END INIT OF VARIABLES
    ###########################################################################

    ###########################################################################
    # START GUI CODE

    Frame.__init__(self, parent, *args, **kw)

    self.expandGrid(0,0)

    ## Single Frame
    # frame = Frame(self, grid=(0,0))

    # or with Tabs?
    options = ['Launch','Output', 'Runs']
    tabbedFrame = TabbedFrame(self, options=options, grid=(0,0))
    frameA, frameB, frameC = tabbedFrame.frames
    self.tabbedFrame = tabbedFrame

    frameA.expandGrid(14,2)

    row = 0
    div = LabelDivider(frameA, text='MODULE Setup', grid=(row,0), gridSpan=(1,4))

    row += 1
    # allow the user to choose MODULE if either the one in PATH is incorrect or not found
    button = Button(frameA, text='Select MODULE executable:',bd=1, \
                            command=self.selectExecutable,  grid=(row,0), sticky="ew")
    self.moduleExeEntry = Entry(frameA, text=self.moduleExePath, grid=(row,1), gridSpan=(1,3), \
                            width=32, sticky="ew", bd=1)
    self.moduleExePath = self.moduleExeEntry.get()

    # separator "MODULE input"
    row += 1
    div = LabelDivider(frameA, text='MODULE input', grid=(row,0), gridSpan=(1,5))

    row += 1
    label = Label(frameA, text='Structure:',        grid=(row,1))
    self.inputStructurePulldown = PulldownList(frameA, self.changeInputStructure, \
                                                    grid=(row,2))
    # self.constraintsFileEntry.bind('<Leave>', self.updateEntryParams)

    row += 1
    label = Label(frameA, text='RDC constraints:',  grid=(row,1))
    self.inputRdcConstraintsPulldown = PulldownList(frameA, self.changeInputRdcConstraintList, \
                                                    grid=(row,2))
    #self.constraintsFileEntry.bind('<Leave>', self.updateEntryParams)

    row += 1
    label = Label(frameA, text='(Optional input)',  grid=(row,0))
    label = Label(frameA, text='Distance constraints:', \
                                                    grid=(row,1))
    self.inputDistanceConstraintsPulldown = PulldownList(frameA, self.changeInputDistanceConstraintList, \
                                                    grid=(row,2))
    #self.constraintsFileEntry.bind('<Leave>', self.updateEntryParams)

    row += 1

    subFrameDepth = 4
    subframe = LabelFrame(frameA, text='MODULE User Notes (store notes about how MODULE was run here)', \
                                                    grid=(row,0), gridSpan=(1,4))
    subframe.expandGrid(subFrameDepth,0)

    self.moduleUserText = ScrolledText( subframe )
    self.moduleUserText.grid(row=subFrameDepth, column=0, columnspan=4, sticky='nsew')

    # View Results
    row += subFrameDepth

    # row += 1
    # div = LabelDivider(frameA, text='MODULE launch', grid=(row,0), gridSpan=(1,4))

    row += 1
    button = Button(frameA, text='Run MODULE', bd=1, command=self.executeModule, \
                                                    grid=(row,0), gridSpan=(1,4), sticky="ew", bg=MODULE_GREEN)
                                                    # grid=(row,0), gridSpan=(1,2), sticky="ew", bg=MODULE_GREEN)


    ###########################################################################
    # Frame B (tab 2) Ouput
    frameB.expandGrid(4,1)
    row = 0

    subFrameDepth = 6
    subframe = LabelFrame(frameB, text='MODULE Output', \
                                                    grid=(row,0), gridSpan=(1,5))
    #subframe.grid_columnconfigure(2, weight=1)
    subframe.expandGrid(subFrameDepth,0)

    self.moduleOutputText = ScrolledText( subframe )
    self.moduleOutputText.setState( state=Tkinter.DISABLED )
    self.moduleOutputText.grid(row=subFrameDepth, column=0, columnspan=4, sticky='nsew')

    # separator "MODULE input"
    row += 1
    div = LabelDivider(frameB, text='MODULE RDC Back Values', grid=(row,0), gridSpan=(1,5))

    row += 1
    button = Button(frameB, text='Import MODULE Back Values file', bd=1, command=self.importModuleBackValues, \
                                                    grid=(row,0), gridSpan=(1,4), sticky="ew", bg=MODULE_BLUE)
                                                    # grid=(row,0), gridSpan=(2,4), sticky="ew", bg=MODULE_BLUE)

    row += 1
    self.rdcOutputTable = None
    frameB.grid_rowconfigure(row, weight=1)
    headings = ('#', 'Resonances', 'Value', 'Back Value', 'Diff.', 'Error')

    editWidgets       = [None, None, None, None, None, None]
    editGetCallbacks  = [None, None, None, None, None, None]
    editSetCallbacks  = [None, None, None, None, None, None]

    self.rdcOutputTable = ScrolledMatrix(frameB,  headingList=headings,
                                                  multiSelect=False,
                                                  editWidgets=editWidgets,
                                                  editGetCallbacks=editGetCallbacks,
                                                  editSetCallbacks=editSetCallbacks,
                                                  initialRows=4 )

    self.rdcOutputTable.grid(row=row, column=0, columnspan=4, sticky='nsew')

    row += 1
    button = Button(frameB, text='Import MODULE Structure', bd=1, command=self.importModuleStructure, \
                                                    grid=(row,0), gridSpan=(1,4), sticky="ew", bg=MODULE_BLUE)
                                                    # grid=(row,0), gridSpan=(2,4), sticky="ew", bg=MODULE_BLUE)


    ###########################################################################
    # Frame C (tab 3) NMR Calc display bits
    frameC.expandGrid(4,1)
    row = 0

    div = LabelDivider(frameC, text='Stored MODULE Runs', grid=(row,0), gridSpan=(1,5))

    # NmrCalc Run scrolled matrix
    row += 1
    self.runTable = None
    frameC.grid_rowconfigure(row, weight=1)
    headings = ('Run ID', 'notes', 'Status')

    # self.editRunNotes = DataEntry.askString('Run Notes', 'Edit notes about Run', tipText='Notes about Run', parent=self)
    # editWidgets       = [None, self.editRunNotes, None]

    editWidgets       = [None, None, None]
    editGetCallbacks  = [None, None, None]
    editSetCallbacks  = [None, None, None]

    self.runTable     = ScrolledMatrix(frameC,  headingList=headings,
                                                multiSelect=False,
                                                editWidgets=editWidgets,
                                                editGetCallbacks=editGetCallbacks,
                                                editSetCallbacks=editSetCallbacks,
                                                initialRows=4 )

    self.runTable.grid(row=row, column=0, columnspan=4, sticky='nsew')

    row += 4
    tipTexts  = ['Load Selected Run', 'Delete Selected Run']
    texts     = ['Load Selected Run', 'Delete Selected']
    commands  = [self.loadRun, self.deleteRun] 
    colours   = [MODULE_GREEN, MODULE_RED]
    self.runButtons = ButtonList(frameC, texts=texts, tipTexts=tipTexts,
                                      commands=commands, grid=(row,0), gridSpan=(1,4) )
    self.runButtons.buttons[0].config(bg=MODULE_GREEN)
    self.runButtons.buttons[1].config(bg=MODULE_RED)

    ###########################################################################
    # Keep GUI up to date

    self.updateAfter()
    self.administerNotifiers( self.parent.registerNotify )

  # END GUI CODE
  ###########################################################################

  ###########################################################################
  # FUNCTIONS CALLED FROM GUI

  ################
  # OS interaction
  def selectExecutable(self):
    """ Choose the MODULE executable appropriate for your system. """

    popup = FileSelectPopup(self, show_file=True)

    executable = popup.getFile()
    if executable:
      self.moduleExeEntry.set(executable)

    popup.destroy()

  def executeModule(self):
    """ Execute MODULE with all given parameters 
        Yeah it's a big function, bite me :-P
    """

    # check that all params have been set
    if not os.path.isfile( self.moduleExePath ):
      popup = showWarning('MODULE', 'Cannot find MODULE executable.', parent=self)
      return

    if self.inputStructure is None:
      popup = showWarning('MODULE', 'MODULE cannot run with out a structure for input.', parent=self)
      return

    if self.inputRdcConstraintList is None:
      popup = showWarning('MODULE', 'MODULE cannot run with out an RDC constraint list for input.', parent=self)
      return

    # write temp PDB file
    strucTempFile = tempfile.NamedTemporaryFile( mode='w', suffix='.pdb', dir=self.moduleTempDir, delete=False)
    strucTempFile.close()
    makePdbFromStructure( strucTempFile.name, self.inputStructure.structureEnsemble, model=self.inputStructure)

    # write temp RDC file
    rdcTempFile   = tempfile.NamedTemporaryFile( mode='w', suffix='.tab', dir=self.moduleTempDir, delete=False)

    moduleIo.writeConstraintList(rdcTempFile, self.inputRdcConstraintList)
    rdcTempFile.close()

    # command to be executed (to start MODULE)
    moduleCommands = [self.moduleExePath, strucTempFile.name, rdcTempFile.name]

    # if user has set Distance constraint list then fill add this (user optional)
    if self.inputDistanceConstraintList:
      distanceTempFile = tempfile.NamedTemporaryFile( mode='w', suffix='.tab', dir=self.moduleTempDir, delete=False)
      moduleIo.writeConstraintList(distanceTempFile, self.inputDistanceConstraintList)
      distanceTempFile.close()
      moduleCommands.append( distanceTempFile.name )

    # execute MODULE
    moduleOutput = subprocess.Popen( moduleCommands, stdout=subprocess.PIPE).communicate()[0]

    # clean out the blank lines
    for line in moduleOutput:
      if line == '':
        moduleOutput.pop( line )

    # send the output to our screen
    self.moduleOutputText.setState( state=Tkinter.NORMAL )
    self.moduleOutputText.setText( text=moduleOutput )
    self.moduleOutputText.setState( state=Tkinter.DISABLED )

    # clean up temp files (delete has been set to False)
    for tempFile in [strucTempFile.name, rdcTempFile.name]:
      if os.path.isfile( tempFile ):
        os.unlink( tempFile )
    if self.inputDistanceConstraintList:
      if os.path.isfile( distanceTempFile.name ):
        os.unlink( distanceTempFile.name )

    # create new run and store inputs
    self.newRun()
    self.updateRunInputs()

  ############################################################################
  # Outputs:

  ############################
  # Import MODULE PDB

  def findModuleExportPdbFile(self):
    """ Find the MODULE PDB export file, typically called temp%s.pdb 
      unless user has renamed it.
    """
    modulePdbFileGood = False
    def yes():    modulePdbFileGood = True
    def cancel(): modulePdbFileGood = False

    # MODULE writes pdb files to CWD no matter where you or it is
    possibleFiles = glob.glob( os.path.join( os.getcwd(), 'temp*' ) )
    # debug
    # possibleFiles = glob.glob( os.path.join( os.getcwd(), 'module', 'temp*' ) )

    if len( possibleFiles ) == 1:

      if os.path.isfile( possibleFiles[0] ):

        texts   = ['Yes', 'No, choose another PDB file', 'Cancel' ]
        objects = [ yes, self.selectModulePdbExport, cancel]
        func = showMulti('MODULE', 'Is this the PDB strcture you exported from MODULE?\n%s' % possibleFiles[0], \
                  texts=texts, objects=objects, parent=self)
        func()

        if modulePdbFileGood == True:
          modulePdbFile = possibleFiles[0]
        else:
          return None

      else:
        modulePdbFile = self.selectModulePdbExport()

    elif len( possibleFiles ) > 1:

      texts   = [ '%s' % fileName for fileName in possibleFiles ]
      objects = [ '%s' % fileName for fileName in possibleFiles ]
      modulePdbFile = showMulti('MODULE', 'Multiple possible MODULE export\nfiles found, please select one:' % possibleFiles, \
                texts=texts, objects=objects, parent=self)

      if os.path.isfile( modulePdbFile ):
        modulePdbFile = modulePdbFile

    else:
      modulePdbFile = self.selectModulePdbExport()

    return modulePdbFile

  def selectModulePdbExport(self):
    """ Choose the PDB file that was exported from MODULE. """

    file_types = [FileType("PDB structure files", ["*.pdb"]), FileType("All files", ["*"])]
    popup = FileSelectPopup(self, file_types, dismiss_text='Cancel', show_file=True)

    chosenPdbFile = popup.getFile()
    if os.path.isfile( chosenPdbFile ):
      modulePdbFile = chosenPdbFile
    else:
      warnPopup = showWarning('MODULE', 'File %s not found.' % chosenPdbFile, parent=self )
      return None

    popup.destroy()

    return modulePdbFile

  def importModuleStructure(self):
    """ Find MODULE structure file and import into Analysis (or attempt to, this will probably
      barf as MODULE pdb files are not good)
    """

    modulePdbFile = self.findModuleExportPdbFile()
    if not modulePdbFile: return

    # this stores the re-hacked MODULE pdb file and then removes it on close
    tempPdb = tempfile.NamedTemporaryFile( mode='w', suffix='.pdb', dir=self.moduleTempDir, delete=False)
    tempPdb.close()

    # if file found then add to project
    if os.path.isfile( modulePdbFile ):
      moduleIo.BlackledgeToPdbConverter( modulePdbFile, tempPdb)
      structure = getStructureFromFile(self.inputStructure.parent.molSystem, tempPdb)
      os.unlink( tempPdb )

  ############################
  # Import MODULE Back Values

  def findModuleExportBackValuesFile(self):
    """ Find the Back Values file that the User hopefully exported from MODULE """

    def yes():    moduleBvFileGood = True
    def cancel(): moduleBvFileGood = False

    # back value files are helpfully appended '*.back'
    possibleFiles  = glob.glob( os.path.join( os.getcwd(), '*.back' ) ) 
    possibleFiles += glob.glob( os.path.join( self.moduleTempDir, '*.back' ) )

    if len( possibleFiles ) == 1:

      if os.path.isfile( possibleFiles[0] ):

        texts   = ['Yes', 'No, choose another Back Values file', 'Cancel' ]
        objects = [ yes, self.selectModuleBvExport, cancel]
        func = showMulti('MODULE', 'Is this the Back Values file that you exported from MODULE?\n%s' % possibleFiles[0], \
                  texts=texts, objects=objects, parent=self)
        func()

        if moduleBvFileGood == True:
          moduleBackValueFile = possibleFiles[0]
        else:
          return None

      else:
        moduleBackValueFile = self.selectModuleBvExport()

    elif len( possibleFiles ) > 1:

      texts   = [ '%s' % fileName for fileName in possibleFiles ]
      objects = [ '%s' % fileName for fileName in possibleFiles ]
      chosenBvFile = showMulti('MODULE', 'Multiple possible MODULE Back Value\nfiles found, please select one:' % possibleFiles, \
                texts=texts, objects=objects, parent=self)

      if os.path.isfile( modulePdbFile ):
        moduleBackValueFile = chosenBvFile

    else:
      moduleBackValueFile = self.selectModuleBvExport()

    return moduleBackValueFile

  def selectModuleBvExport(self):
    """ Choose the Back Value file that was exported from MODULE. """

    file_types = [FileType("Back Value files", ["*.back"]), FileType("All files", ["*"])]
    popup = FileSelectPopup(self, file_types, dismiss_text='Cancel', show_file=True)

    chosenBvFile = popup.getFile()
    if os.path.isfile( chosenBvFile ):
      moduleBackValueFile = chosenBvFile
    else:
      warnPopup = showWarning('MODULE', 'File %s not found.' % chosenBvFile, parent=self )
      return None

    popup.destroy()

    return moduleBackValueFile

  def importModuleBackValues(self):
    """ Attempt to find the Back Values file and then import it ... """

    # find the Back Values file or ask user
    backValFile = self.findModuleExportBackValuesFile()
    if not backValFile: return

    # create a new RDC list
    if self.project:
      nmrConstraintStore = self.project.newNmrConstraintStore(nmrProject=self.project.findFirstNmrProject())

    if self.inputStructure:
      chain = self.inputStructure.parent.molSystem.findFirstChain()
    else:
      chain = None

    if nmrConstraintStore and chain:

      rdcList = moduleIo.getBackValuesListFromFile( backValFile, chain, nmrConstraintStore )

      # rdcRawData = moduleIo.getRawBackValuesFromFile( backValFile )
      # 
      # for run in runs:
      #   runText.append( [run.getSerial(), run.getDetails(), run.getStatus()] )

      self.rdcOutputTable.update(objectList=None, textMatrix=None)
        

  # end OS interaction
  ####################

  ##################################
  # get, set and update menu options

  def changeInputRun(self, run):
    self.run = run

  def changeInputStructure(self, structure):
    self.inputStructure = structure

  def changeInputRdcConstraintList(self, constraintList):
    self.inputRdcConstraintList = constraintList

  def changeInputDistanceConstraintList(self, constraintList):
    self.inputDistanceConstraintList = constraintList

  def changeInputConstraintList(self, constraintList, className):
    if className == 'RdcConstraintList':
      self.changeInputRdcConstraintList(constraintList)
    elif className == 'DistanceConstraintList':
      self.changeInputDistanceConstraintList(constraintList)

  def getInputConstraintList(self, className):
    if className == 'RdcConstraintList':
      return self.inputRdcConstraintList
    elif className == 'DistanceConstraintList':
      return self.inputDistanceConstraintList

  # end get and set
  ##################################

  ####################
  # NMR Calc runs etc.

  def loadRun(self):
    pass
  
  def deleteRun(self):
    pass

  def changeRun(self):
    pass

  def newRun(self):

    if self.calcStore:
      self.run = self.calcStore.newRun(status='provisional')

  def updateRunInputs(self, event=None):

    self.inputUserDescriptionText = self.moduleUserText.getText() or None

    # store input parameters in run
    if self.run:

      # set the path to the MODULE executable used
      setRunTextParameter(self.run, MODULE_EXE, self.moduleExePath)

      # set the input structure sent to MODULE
      sEnsemble = self.inputStructure.StructureEnsemble
      self.run.newStructureEnsembleData( sEnsemble.ensembleId, sEnsemble.molSystem.code, name=STRUCTURE_DATA, ioRole=INPUT )

      # set the input RDC constraint list
      rdcConst = self.inputRdcConstraintList
      self.run.newConstraintStoreData( rdcConst.nmrConstraintStore.serial, name=RDC_CONSTRAINT_DATA, ioRole=INPUT )

      # set the input RDC constraint list (optional)
      distConst = self.inputDistanceConstraintList
      if distConst:
        self.run.newConstraintStoreData( distConst.nmrConstraintStore.serial, name=DIST_CONSTRAINT_DATA, ioRole=INPUT )

      # set the input user data (if any)
      setRunTextParameter(self.run, USER_DESCRIPTION_DATA,  self.inputUserDescriptionText)

  def editRunNotes(self):
    pass

  # end NMR calc stuff
  ####################

  ######################
  # update and noitfiers

  def administerNotifiers(self, notifyFunc):

    for func in ['__init__', 'delete']:
      for clazz in [ 'RdcConstraintList', 'DistanceConstraintList' ]:
        notifyFunc(self.updateAfter, 'ccp.nmr.NmrConstraint.' + clazz, func)

      notifyFunc( self.updateAfter, 'ccp.nmr.NmrCalc.Run', func )

  def destroy(self):

    self.administerNotifiers(self.parent.unregisterNotify)
    Frame.destroy(self)

  def updateAll(self):

    # tab A
    self.updateModuleExecutablePath()
    self.updateInputStructure()
    self.updateInputRdcConstraintList()
    self.updateInputDistanceConstraintList()

    # tab C
    self.updateInputRuns()

    self.waiting = False

  def updateAfter(self, obj=None):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.updateAll)

  def updateInputRuns(self):

    if self.calcStore:
      runs = [r for r in self.calcStore.sortedRuns() if r.status == PROVISIONAL]

    if not runs:
      self.newRun()
      runs = [self.run]

    runText = []

    for run in runs:
      runText.append( [run.getSerial(), run.getDetails(), run.getStatus()] )

    self.runTable.update(objectList=runs, textMatrix=runText)


  def updateModuleExecutablePath(self):
    self.moduleExePath = self.moduleExeEntry.get()

  def updateInputStructure(self):

    index       = None
    # names       = ['<None>']
    # structures  = [None]
    names       = []
    structures  = []


    for molSystem in self.project.sortedMolSystems():
      for ensemble in molSystem.sortedStructureEnsembles():
        for structure in ensemble.sortedModels():
          structures.append(structure)

    if structures:

      for i, model in enumerate(structures):
        if model is None: continue
        ee = model.structureEnsemble
        name = '%s:%d:%d' % (ee.molSystem.code,ee.ensembleId, model.serial) 
        names.append(name)

      if self.inputStructure not in structures:
        self.changeInputStructure(structures[0])

      index = structures.index(self.inputStructure)

    else:
      self.inputStructure = None

    self.inputStructurePulldown.setup(names, structures, index or -1 )

  def updateInputConstraintList(self, className, obj=None):

    index           = None
    # if DistanceConstraintList then None must be an option
    # names           = ['<None>']
    # constraintLists = [None]

    if className == 'DistanceConstraintList':
      names           = ['<None>']
      constraintLists = [None]
    else:
      names           = []
      constraintLists = []

    for nmrConstraintStore in self.parent.project.sortedNmrConstraintStores():
      for constraintList in nmrConstraintStore.sortedConstraintLists():
        if not className or constraintList.className == className:
          constraintLists.append(constraintList)

    if constraintLists:

      for constList in constraintLists:
        if constList is None: continue
        store = constList.nmrConstraintStore
        name = '%d:%d' % (store.serial, constList.serial) 
        names.append(name)

      if self.getInputConstraintList(className) not in constraintLists:
        self.changeInputConstraintList(constraintLists[0], className)

      index = constraintLists.index(self.getInputConstraintList(className))

    else:
      self.changeInputConstraintList(None, className)

    if className == 'RdcConstraintList':
      self.inputRdcConstraintsPulldown.setup(names, constraintLists, index or -1)
    elif className == 'DistanceConstraintList':
      self.inputDistanceConstraintsPulldown.setup(names, constraintLists, index or -1)

  def updateInputRdcConstraintList(self):
    self.updateInputConstraintList( 'RdcConstraintList' )

  def updateInputDistanceConstraintList(self):
    self.updateInputConstraintList( 'DistanceConstraintList' )

  # end notifiers
  ##################################








