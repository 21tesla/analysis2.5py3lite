
"""
======================COPYRIGHT/LICENSE START==========================

EditCalculation.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
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

===========================REFERENCE END===============================
"""

from ccp.util.NmrCalc import getDataObjText, DATA_MISSING, getRangeString

from ccpnmr.analysis.popups.BasePopup import BasePopup

from memops.general.Util import copySubTree

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.DataEntry import askString
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.LabelDivider import LabelDivider
from memops.gui.MessageReporter import  showOkCancel, showWarning
from memops.gui.MultiWidget import MultiWidget
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.ScrolledText import ScrolledText

from memops.gui.ScrolledFrame import ScrolledFrame
from memops.gui.TabbedFrame import TabbedFrame

NICE_RED   = '#FF8080'
NICE_GREEN = '#80FF80'

STRUCTURE_DATA = 'StructureEnsembleData'
MEASUREMENT_DATA = 'MeasurementListData'
PEAK_DATA = 'PeakListData'
CONSTRAINT_DATA = 'ConstraintStoreData'
RESIDUE_DATA = 'MolResidueData'

DATA_TITLES = {STRUCTURE_DATA:'Structures',
               MEASUREMENT_DATA:'Shifts & Measurements',
               PEAK_DATA:'Peak Lists',
               CONSTRAINT_DATA:'Restraint Lists',}
               
# TO-DO
#
# More table types
#
#
#               

def testEditCalculationPopup(argServer):

  popup = EditCalculationPopup(argServer.parent)
  popup.open()
  

class EditCalculationPopup(BasePopup):
  """
  **View and Edit External NMR Calculations**
  
  This popup window is used to view the various calculations, performed by
  programs external to CCPN (like CING or ARIA), that have nonetheless been
  setup inside CCPN software to use CCPN data. 

  *This popup is somewhat of a prototype, so its functionality is neither
  complete nor throughly tested, but it is the only way to get an overview of
  certain types of data.*

  The popup window is divided into two tabs. The first "Calculation Groups" tabs
  shows the various types of software that calculations have been setup for. It
  should be noted that this system is general and where possible the user should
  consider using the dedicated graphical interfaces to run particular pieces of
  software (like those in the main Structure menu).

  The "Calculation Runs" tab shows the various jobs that have been setup, and
  possibly run, for the different kinds of software. The user can change the
  software calculation group via the top most pulldown menu, then change the
  "Run #" to look though the jobs ("runs") that have been setup. New, blank
  or duplicated, "run" specifications can be made at any time. A specification
  may also be deleted, but a run may only be edited if it has not been used by
  the calculation; a run is a record or what was done. A completed run may
  effectively be edited by making a copy first.

  The data that relates to a job/run is divided into three sub-tabs to indicate
  what the input data from the CCPN project is, what settings are used during
  the calculation, and what the output CCPN data is. The input data is further
  sub-divided into various categories of data, and the user can (if the job is
  editable) use the [Add ...] and [Remove ...] buttons to dictate which CCPN
  entities will be selected for the calculation. Not all types of data will be
  appropriate for all calculations; further incentive to use the specialised 
  program interfaces when they are available, rather than this general one.

  The "Run Settings" tab is filled with a table of all of the input parameters
  that have been associated with the currently selected run. A parameter will
  have an identifying code and perhaps a group to collect certain parameters
  together. Some parameters may just represent a simple value, but others
  may associate a value with a "Data Object"; a CCPN data model entity
  containing NMR information or a reference to external data held on disk.
  
  The "Output Summary" table lists The calculation summary or overview produced 
  by the program (if any)).
  
  The "Output Data" table lists all of the CCPN entities that were generated or
  modified by the selected calculation job (if it were run).

  It is notable that a "run" specification will persist even if some or all of
  the CCPN data entities (e.g. restraint lists) it was originally linked to have
  been deleted. In this way the user still has a record that something was run,
  and there is no restriction on the normal manipulation of CCPN data. If a CCPN
  data item has disappeared then the run will list it as "*DATA MISSING*".
  """


  def __init__(self, parent, *args, **kw):

    self.calcStore = None
    self.waiting = False

    BasePopup.__init__(self, parent, title="Other : NMR Calculations", **kw)

  def body(self, guiFrame):

    self.geometry('700x600')

    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)

    options = ['Calculation Groups','Calculation Runs']
      
    tabbedFrame = TabbedFrame(guiFrame, options=options, grid=(0,0), 
                              callback=self.update)
    self.tabbedFrame = tabbedFrame
    frameA, frameB = tabbedFrame.frames

    #
    # Calculation Groups
    #

    frameA.expandGrid(0,0)
    tipTexts = ['Row number',
                'name of calculation group (often a program name)',
                'The number of calculation run specifications stored for the group']
    headingList = ['#','Name','Num. Runs']
    editWidgets      = [None, None, None]
    editGetCallbacks = [None, None, self.editRuns]
    editSetCallbacks = [None, None, None]
    
    self.calcStoreMatrix = ScrolledMatrix(frameA, 
                                          editSetCallbacks=editSetCallbacks,
                                          editGetCallbacks=editGetCallbacks, 
                                          editWidgets=editWidgets,
                                          headingList=headingList, 
                                          callback=self.selectCalcStore,
                                          deleteFunc=self.deleteCalcStore,
                                          grid=(0,0), tipTexts=tipTexts)
     
    
    tipTexts = ['Add a new grouping (e.g. by program name) to contain calculation specifications',
                'Delete the selected calculation group (program settings) and all calculation run settings it contains']
    texts = ['Add New Calculation Group',
             'Delete Calculation Group']
    commands = [self.addCalcStore, self.deleteCalcStore] 
    self.calcStoreButtons = ButtonList(frameA, texts=texts, tipTexts=tipTexts,
                                      commands=commands, grid=(1,0))
     
    #
    #  Runs 
    #
    
    frameB.expandGrid(1,1) 
    label = Label(frameB, text='Calculation Group:', grid=(0,0))
    
    self.editRunFrame = NmrCalcRunFrame(frameB, self.project, self.calcStore,
                                       grid=(1,0), gridSpan=(1,2))
     
    tipText = 'Selects which program (or other calculation grouping) a set of calculations relate to'
    self.calcStorePulldown = PulldownList(frameB, callback=self.changedCalcStore,
                                     grid=(0,1), tipText=tipText)
    
    #
    # Main
    #
    
    buttonList = UtilityButtonList(tabbedFrame.sideFrame, grid=(0,0), sticky='e')
    
    self.administerNotifiers(self.registerNotify)
    self.editRunFrame.update()
    self.updateCalcStoresAfter()
    
  def open(self):
    
    self.update()
    BasePopup.open(self)
  
  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.NmrCalc.NmrCalcStore',
                    'ccp.nmr.NmrCalc.Run'):
        notifyFunc(self.updateCalcStoresAfter,clazz, func)
    
    self.editRunFrame.administerNotifiers(notifyFunc)
        
  def selectCalcStore(self, obj, row, col):
  
    self.calcStore = obj
  
  def addCalcStore(self):
  
    nmrProject = self.nmrProject
    name = askString('Request','Enter name for calculation group', parent=self)
    name = name and name.strip()
    if not name:
      return
    
    if nmrProject.findFirstNmrCalcStore(name=name):
      msg = 'Name "%s" already in use for another calculation group' % name
      showWarning('Failure', msg, parent=self)
      return
     
    if len(name.splitlines()) > 1:
      msg = 'Name cannot contain line breaks'
      showWarning('Failure', msg, parent=self)
      return
     
    if len(name) > 80:
      msg = 'Name must contain 80 characters or fewer'
      showWarning('Failure', msg, parent=self)
      return
  
    calcStore = self.project.newNmrCalcStore(name=name, nmrProject=nmrProject)
  
  def deleteCalcStore(self):
  
    if self.calcStore:
      msg = 'Really delete NMR Calculation group "%s"?' % self.calcStore.name
      
      if showOkCancel('Query', msg, parent=self):
        self.calcStore.delete()
        self.calcStore = None
      
     
  def editRuns(self, calcStore):
    
    self.changedCalcStore(calcStore)
    self.tabbedFrame.select(1)

  def changedCalcStore(self, calcStore):
  
    if calcStore is not self.calcStore:
      self.calcStore = calcStore
      self.updateRunsAfter()
  
  def update(self, selected=None):
    
    if selected is None:
      selected = self.tabbedFrame.selected
    
    if selected == 0:
      self.updateCalcStoresAfter()
    else:
      self.updateCalcStorePulldown()
      self.updateRunsAfter()
  
  def updateCalcStoresAfter(self, calcStore=None):
  
    self.after_idle(self.updateCalcStores)
  
  def updateCalcStores(self):
  
    textMatrix = []
    objectList = self.nmrProject.sortedNmrCalcStores()
    
    for i, calcStore in enumerate(objectList):
      datum = [i+1, calcStore.name, len(calcStore.runs)]
      textMatrix.append(datum)
     
    self.calcStoreMatrix.update(textMatrix=textMatrix, objectList=objectList)

  def updateCalcStorePulldown(self):
  
    calcStore = self.calcStore
    index = 0
    names = []
    calcStores = self.nmrProject.sortedNmrCalcStores()
    
    if calcStores:
      names = [s.name for s in calcStores]
      
      if calcStore not in calcStores:
        calcStore = calcStores[0]
      
      index = calcStores.index(calcStore)
    
    else:
      calcStore = None
      
    if calcStore is not self.calcStore:
      self.changedCalcStore(calcStore)  
    
    self.calcStorePulldown.setup(names, calcStores, index)
     
  def updateRunsAfter(self, obj=None):    
      
    if self.waiting:
      return
      
    else:
      self.waiting = True
      self.after_idle(self.updateRuns)     
  
  def updateRuns(self):
  
    self.editRunFrame.update(self.calcStore)
    self.waiting = False    

  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
        
    BasePopup.destroy(self)


class NmrCalcRunFrame(Frame):

  def __init__(self, parent, project, calcStore=None, closeButton=False,
               inputTypes=(STRUCTURE_DATA, MEASUREMENT_DATA, PEAK_DATA, CONSTRAINT_DATA),
               chainSelection=False, *args, **kw):

    self.parent = parent
    self.project = project
    self.calcStore = calcStore
    self.run = None
    self.structureData = None
    self.measurementData = None
    self.peakListData = None
    self.constraintData = None
    self.chainSelection = chainSelection
    self.inputTypes = inputTypes
    self.waiting = False
    
    Frame.__init__(self, parent, *args, **kw)
    
    self.expandGrid(0,0)
    
    options = ['Input Data','Run Settings','Output Data','Summary']
      
    tabbedFrame = TabbedFrame(self, options=options, grid=(0,0), callback=self.toggleTab)
    frameA, frameB, frameC, frameD = tabbedFrame.frames
    self.tabbedFrame = tabbedFrame
    
    label = Label(tabbedFrame.sideFrame, text='Run #:',
                  grid=(0,0), sticky='e')
                  
    tipText = 'Selects which calculation job or "run" is currently being viewed or edited'
    self.runPulldown = PulldownList(tabbedFrame.sideFrame, callback=self.changeRun,
                                    grid=(0,1), sticky='e', tipText=tipText)
    
    tipTexts = ['Make a setup for a new calculation run',
                'Make a new calculation run by copying the current one',
                'Delete the current calculation run settings']
    texts = ['New Run', 'Copy Run', 'Delete Run']
    commands = [self.newRun, self.copyRun, self.deleteRun]
    
    if closeButton:
      ButtonListClass = UtilityButtonList
    else:
      ButtonListClass = ButtonList
    
    runButtons = ButtonListClass(tabbedFrame.sideFrame, texts=texts, tipTexts=tipTexts,
                                 commands=commands, sticky='e', grid=(0,2))
    runButtons.buttons[0].config(bg='#B0FFB0')
    
    # Input
    
    frameA.expandGrid(1,0)
    frameA.expandGrid(2,0)
    frameA.expandGrid(3,0)
    frameA.expandGrid(4,0)
    
    subframe = Frame(frameA, grid=(0,0))
    subframe.grid_columnconfigure(3, weight=1)
    
    label = Label(subframe, text='Molecular System:', grid=(0,0))
    tipText = 'Selects which molecular system the calculation will be run on'
    self.molSystemPulldown = PulldownList(subframe, callback=self.changeMolSystem,
                                          grid=(0,1), tipText=tipText)
    
    if chainSelection:
      label = Label(subframe, text='Chain Selection:', grid=(0,2))
      tipText = 'Selects which particular chains to consider form the selected molecular system'
      self.chainPulldown = PulldownList(subframe, callback=self.changeChain,
                                        grid=(0,3), tipText=tipText)
        
    label = Label(subframe, text='Protocol name:', grid=(0,4))
    tipText = 'Name of calculation protocol for the current calculation'
    self.protocolEntry = Entry(subframe, text='', grid=(0,5),
                              width=32, tipText=tipText)
    self.protocolEntry.bind('<Leave>', self.changeWmsProtocol, '+')
        
    label = Label(subframe, text='Run Details:', grid=(0,6))
    tipText = 'User-editable textual comment for the current calculation'
    self.detailsEntry = Entry(subframe, text='', grid=(0,7),
                              width=32, tipText=tipText)
    self.detailsEntry.bind('<Leave>', self.changeDetails, '+')

    
    options = [DATA_TITLES[x] for x in inputTypes]
    inputTabs = TabbedFrame(frameA, options=options, grid=(1,0))
    self.inputTabs = inputTabs
    
    frameDict = {}
    for i, frame in enumerate(inputTabs.frames):
      frameDict[inputTypes[i]] = frame
    
    self.buttons = []
    
    # # # # #  S T R U C T U R E S  # # # # # 
    
    if STRUCTURE_DATA in inputTypes:
      frame = frameDict[STRUCTURE_DATA]
      frame.expandGrid(0,1)
 
      tipTexts = ['Row number',
                  'Code for the molecular system to which the structure relates',
                  'The id number of the structure ensemble within its molecular system',
                  'The numbers of the conformational models considered']
      headingList = ['#','Mol System','Ensemble','Models']
      self.structureTable = ScrolledMatrix(frame, grid=(0,0), gridSpan=(1,4),
                                           callback=self.selectModel,
                                           multiSelect=True, tipTexts=tipTexts,
                                           headingList=headingList, initialRows=2)
 
 
      tipTexts = ['Remove the selected structural models from consideration (does not delete the coordinates)',
                  'Add the selected structure ensemble as calculation input']
      texts = ['Remove Selected Models','Add Ensemble:']
      commands = [self.removeModel, self.addEnsemble]
      self.structureButtons = ButtonList(frame, texts=texts, commands=commands,
                                         grid=(1,0), tipTexts=tipTexts)
      self.buttons.extend(self.structureButtons.buttons)
      tipText = 'Selects a structure ensemble which may be added as calculation input'
      self.ensemblePulldown = PulldownList(frame, grid=(1,1), tipText=tipText)
 
      tipTexts = ['Add the selected conformational model to the calculation input',]
      texts = ['Add Model:']
      commands = [ self.addModel]
      self.modelButtons = ButtonList(frame, texts=texts, commands=commands,
                                     grid=(1,2), tipTexts=tipTexts)
      self.buttons.extend(self.modelButtons.buttons)
      tipText = 'Selects a conformational model which may be added as calculation input'
      self.modelPulldown = PulldownList(frame, grid=(1,3), tipText=tipText)
 
    # # # # #   M E A S U R E M E N T S   # # # # #
   
    if MEASUREMENT_DATA in inputTypes:
      frame = frameDict[MEASUREMENT_DATA]
      frame.expandGrid(0,0)
 
      tipTexts = ['Row number',
                  'The kind of measurement data used and calculation input',
                  'The serial number used to uniquely identify the measurement list',
                  'The textual name of the measurement list',
                  'The number of individual measurements in the list',
                  'User-specified textual comment']
      headingList = ['#','Data Type','Serial','List Name','Size','Details']
      self.measurementsTable = ScrolledMatrix(frame, grid=(0,0), gridSpan=(1,2),
                                              callback=self.selectMeasurementList,
                                              multiSelect=True, tipTexts=tipTexts,
                                              headingList=headingList, initialRows=2)
 
      tipTexts = ['Remove the selected measurement list from consideration (does not delete the data)',
                  'Add the selected measurement list as calculation input']
      texts = ['Remove Selected','Add Measurement List:']
      commands = [self.removeMeasurementList, self.addMeasurementList]
      buttons = ButtonList(frame, texts=texts, commands=commands,
                           grid=(1,0), tipTexts=tipTexts)
      self.buttons.extend(buttons.buttons)
 
      tipText = 'Selects a measurement list which may be added as calculation input'
      self.measurementListPulldown = PulldownList(frame, grid=(1,1), tipText=tipText)
   
    # # # # #   P E A K  L I S T S   # # # # # 
    
    if PEAK_DATA in inputTypes:
      frame = frameDict[PEAK_DATA]
      frame.expandGrid(0,0)
 
      tipTexts = ['Row number',
                  'The name of the experiment that contains the input peak list',
                  'The name of the spectrum that contains the input peak list',
                  'The serial number of the input peak list within its spectrum',
                  'The number of peaks in the peak list']
      headingList = ['#','Experiment','Spectrum','Serial','Peaks']
      self.peakListTable = ScrolledMatrix(frame, grid=(0,0), gridSpan=(1,2),
                                          callback=self.selectPeakListData,
                                          multiSelect=True, tipTexts=tipTexts,
                                          headingList=headingList, initialRows=2)
 

      tipTexts = ['Remove the selected peak list from consideration (does not delete the data)',
                  'Add the selected peak list as calculation input']
      texts = ['Remove Selected','Add Peak List']
      commands = [self.removePeakList, self.addPeakList]
      buttons = ButtonList(frame, texts=texts, commands=commands,
                           grid=(1,0), tipTexts=tipTexts)
      self.buttons.extend(buttons.buttons)
      
      tipText = 'Selects a peak list which may be added as calculation input'
      self.peakListPulldown = PulldownList(frame, grid=(1,1), tipText=tipText)

    # # # # #   C O N S T R A I N T  L I S T S   # # # # # 
    if CONSTRAINT_DATA in inputTypes:
      frame = frameDict[CONSTRAINT_DATA]
      frame.expandGrid(0,0)
 
      tipTexts = ['Row number',
                  'The type of restraints in the input restraint list',
                  'The restraint set that contains the input restraint list',
                  'The serial number of the input restraintlist within its restraint set',
                  'The textual name of the input restraint list',
                  'The number of restraints in the input restraint list']
      headingList = ['#','Restraint\nType','Restraint\nSet',
                     'Serial','List\nName','Restraints']
      self.constraintListTable = ScrolledMatrix(frame, tipTexts=tipTexts,
                                                callback=self.selectConstraintData,
                                                headingList=headingList, initialRows=2,
                                                multiSelect=True,
                                                grid=(0,0), gridSpan=(1,2))
 
      tipTexts = ['Remove the selected restraint list from consideration (does not delete the data)',
                  'Add the selected restraint list as calculation input']
      texts = ['Remove Selected','Add Restraint List']
      commands = [self.removeConstraintList, self.addConstraintList]
      buttons = ButtonList(frame, texts=texts, commands=commands,
                           grid=(1,0), tipTexts=tipTexts)
      self.buttons.extend(buttons.buttons)
      tipText = 'Selects a restraint list which may be added as calculation input'
      self.constraintListPulldown = PulldownList(frame, grid=(1,1), tipText=tipText)

    # # # # # # #  P A R A M S  # # # # # # # 
    
    if self.__class__ is  NmrCalcRunFrame:
      frameB.expandGrid(0,0)
      tipTexts = ['The serial number of the parameter; an internal CCPN numbering',
                  'The identifying code (one word name) of the parameter',
                  'Human readable name or decription of the parameter',
                  'The current value associated with the parameter',
                  'The data object associated with the parameter']
 
      headingList = ['#','Code','Name','Value(s)','Data Object']
      self.paramsTable = ScrolledMatrix(frameB, headingList=headingList,
                                        grid=(0,0), tipTexts=tipTexts)
   
    
   
    # # # # # # #  O U T P U T  # # # # # # # 
    
    frameC.expandGrid(0,0)
    tipTexts = ['The kind of data output by the calculation',
                'the serial number of the output data, within the calculation run',
                'The name of the output data',
                'The identity of any CPN objects that are lined to/as output data',
                'Any parameters which relate to the output data',]
                #'Any fractional weighting that was applied to the output data']
    headingList = ['Data Type','#','Name','CCPN Data','Parameters',]#'Weight']
    self.outputTable = ScrolledMatrix(frameC,headingList=headingList,
                                      grid=(0,0), tipTexts=tipTexts)
   
    
   
    # # # # # # #  O U T P U T  S U M M A R Y # # # # # # # 
    
    frameD.expandGrid(0,0)
    self.outputSummary = ScrolledText(frameD, grid=(0,0))

    #self.update()


  def administerNotifiers(self, notifyFunc):
  
    NC = 'ccp.nmr.NmrCalc.'
  
    for func in ('__init__', 'delete', 'setDetails'):
      notifyFunc(self.updateRunAfter, NC+'Run', func)

    for func in ('__init__', 'delete', 'setGroupId','setName',
                 'setFloatValue','setIntValue','setBooleanValue',
                 'setTextValue',):
      notifyFunc(self.updateRunParam, NC+'RunParameter', func)

    for className in ('ConstraintStoreData', 'DerivedListData', 'ExternalData',
                      'MeasurementListData', 'MolResidueData', 'MolSystemData',
                      'PeakListData', 'SpectrumData', 'SpinSystemData',
                      'StructureEnsembleData', 'ViolationListData'):
      for func in ('__init__', 'delete', 'setDetails', 'setIoRole',
                   'setWeight', 'setName', 'addRunParameter',
                   'removeRunParameter', 'setRunParameters'):
        notifyFunc(self.updateRunData, NC+className, func)
 
      
    notifyFunc(self.updateRunData, NC+'ConstraintStoreData', 'setConstraintStoreSerial')
    notifyFunc(self.updateRunData, NC+'ConstraintStoreData', 'setConstraintListSerials')
    notifyFunc(self.updateRunData, NC+'DerivedListData', 'setDerivedDataListSerial')
    notifyFunc(self.updateRunData, NC+'ExternalData', 'setDataStore')
    notifyFunc(self.updateRunData, NC+'MeasurementListData', 'setMeasurementListSerial')
    notifyFunc(self.updateRunData, NC+'MolResidueData', 'setMolSystemCode')
    notifyFunc(self.updateRunData, NC+'MolResidueData', 'setChainCode')
    notifyFunc(self.updateRunData, NC+'MolResidueData', 'setResidueSeqIds')
    notifyFunc(self.updateRunData, NC+'MolSystemData', 'setMolSystemCode')
    notifyFunc(self.updateRunData, NC+'MolSystemData', 'setChainCodes')
    notifyFunc(self.updateRunData, NC+'PeakListData', 'setExperimentSerial')
    notifyFunc(self.updateRunData, NC+'PeakListData', 'setDataSourceSerial')
    notifyFunc(self.updateRunData, NC+'PeakListData', 'setPeakListSerial')
    notifyFunc(self.updateRunData, NC+'SpectrumData', 'setExperimentSerial')
    notifyFunc(self.updateRunData, NC+'SpectrumData', 'setDataSourceSerial')
    notifyFunc(self.updateRunData, NC+'SpinSystemData', 'setResonanceGroupSerial')
    notifyFunc(self.updateRunData, NC+'StructureEnsembleData', 'setMolSystemCode')
    notifyFunc(self.updateRunData, NC+'StructureEnsembleData', 'setEnsembleId')
    notifyFunc(self.updateRunData, NC+'StructureEnsembleData', 'setModelSerials')
    notifyFunc(self.updateRunData, NC+'ViolationListData', 'setConstraintStoreSerial')
    notifyFunc(self.updateRunData, NC+'ViolationListData', 'setViolationListSerial ')
    
    for func in ('__init__','delete'):
      notifyFunc(self.updateMolSystemPulldown, 'ccp.molecule.MolSystem.MolSystem', func)
      notifyFunc(self.updateChainPulldown, 'ccp.molecule.MolSystem.Chain', func)

      if PEAK_DATA in self.inputTypes:
        notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.PeakList', func)

      if MEASUREMENT_DATA in self.inputTypes:
        notifyFunc(self.updateMeasurementLists, 'ccp.nmr.Nmr.ShiftList', func)
      
      if STRUCTURE_DATA in self.inputTypes:
        notifyFunc(self.updateEnsembles, 'ccp.molecule.MolStructure.StructureEnsemble', func)
        notifyFunc(self.updateModels, 'ccp.molecule.MolStructure.Model', func)
 
      if CONSTRAINT_DATA in self.inputTypes:
        notifyFunc(self.updateConstraintLists, 'ccp.nmr.NmrConstraint.NmrConstraintStore', func)
        for className in ('ChemShiftConstraintList', 'CsaConstraintList',
                          'DihedralConstraintList', 'DistanceConstraintList',
                          'HBondConstraintList', 'JCouplingConstraintList',
                          'RdcConstraintList'):
          notifyFunc(self.updateConstraintLists, 'ccp.nmr.NmrConstraint.'+className, func)
        
  
  def updateRunData(self, runData):
  
    self.updateRunAfter(runData.run)
  
  def updateRunParam(self, runParam):
  
    self.updateRunAfter(runParam.run)

  def runIsEditable(self, run=None):
  
    if not run:
      run = self.run
  
    if run:
      for datum in run.data:
        if datum.ioRole == 'output':
          return False
      
      else:
        return True
  
    return False

  def copyRun(self):
    
    # TBD: Inputs only?
  
    if self.run:
      self.configure(cursor="watch")
      run = copySubTree(self.run, self.calcStore)
      run.status = 'provisional'
 
      self.run = run
      self.updateRuns()
  
  def newRun(self):

    if self.calcStore:
      self.run = self.calcStore.newRun(status='provisional')
      

  def deleteRun(self):
  
    if self.run:
      msg = 'Really delete calculation run %d?' % self.run.serial
      
      if showOkCancel('Query', msg, parent=self):
        self.run.delete()
        self.run = None
 
  def updateRunAfter(self, run=None):
    
    if self.waiting:
      return
    
    if run and (run.nmrCalcStore is not self.calcStore):
      return
    
    self.waiting = True
    self.configure(cursor="watch")
    self.after_idle(self.update)
 
  def toggleTab(self, index):
  
    self.update()
 
  def update(self, calcStore=None):

    # Should confirm any changes where a run has generated output
    self.configure(cursor="watch")
  
    if not calcStore:
      calcStore = self.calcStore
    else:
      self.calcStore = calcStore

    if self.calcStore:
      self.project = self.calcStore.root
      self.nmrProject = self.project.currentNmrProject or \
                        self.project.findFirstNmrProject() or \
                        self.project.newNmrProject(name='default')
      
      if (not self.run) or (self.run.parent is not self.calcStore):
        self.run = self.calcStore.findFirstRun()
    
    
    buttons = self.buttons
 
    if self.runIsEditable():
      for button in buttons:
        button.enable()
    
    else:
      for button in buttons:
        button.disable()
    
    self.updateRuns()
    self.updateMolSystemPulldown()
    self.updateDetails()
    self.updateProtocol()
    
    
    if MEASUREMENT_DATA in self.inputTypes:
      self.updateMeasurementLists()
    
    if PEAK_DATA in self.inputTypes:
      self.updatePeakLists()
      
    if CONSTRAINT_DATA in self.inputTypes:
      self.updateConstraintLists()
      
    if STRUCTURE_DATA in self.inputTypes:
      self.updateEnsembles()
      self.updateModels()

    run = self.run
    
    structureData = []
    measurementData = []
    peakListData = []
    constraintData = []
    
    if run:
      for dataObj in run.sortedInputs():
        className = dataObj.className
        
        if className == STRUCTURE_DATA:
          structureData.append(dataObj)
          
        elif className == MEASUREMENT_DATA:
          measurementData.append(dataObj)

        elif className == PEAK_DATA:
          peakListData.append(dataObj)
    
        elif className == CONSTRAINT_DATA:
          constraintData.append(dataObj)
    
    # Structure table
    if STRUCTURE_DATA in self.inputTypes:
      textMatrix = []
      objectList = []
 
      for i, dataObj in enumerate(structureData):
         structure = dataObj.structureEnsemble

         if structure:
           serials = list(dataObj.modelSerials)
           
           if serials:
             rangeText = getRangeString(serials)
           else:
             rangeText = None
           
           
           datum = [i+1,
                    structure.molSystem.code,
                    structure.ensembleId,
                    rangeText]
 
         else:
           datum = [i+1, DATA_MISSING,
                    None, None]
 
         textMatrix.append(datum)
         objectList.append(dataObj)
 
      self.structureTable.update(textMatrix=textMatrix, objectList=objectList)

    # Measurement list table
    if MEASUREMENT_DATA in self.inputTypes:
      textMatrix = []
      objectList = []
    
      for i, dataObj in enumerate(measurementData):
         mList = dataObj.measurementList

         if mList:
 
           datum = [i+1, mList.className,
                    mList.serial,
                    mList.name,
                    len(mList.measurements),
                    mList.details]
 
         else:
           datum = [i+1, DATA_MISSING,
                    None, None, None, None]
 
         textMatrix.append(datum)
         objectList.append(dataObj)
 
      self.measurementsTable.update(textMatrix=textMatrix, objectList=objectList)

    # Peak List Table
    if PEAK_DATA in self.inputTypes:
      textMatrix = []
      objectList = []
 
      for i, dataObj in enumerate(peakListData):
         peakList = dataObj.peakList

         if peakList:
           serial = peakList.serial
           spectrum = peakList.dataSource
           experiment = spectrum.experiment
 
           datum = [i+1, experiment.name,
                    spectrum.name,
                    serial,
                    len(peakList.peaks)]
 
         else:
           datum = [i+1, DATA_MISSING,
                    None, None, None]
 
         textMatrix.append(datum)
         objectList.append(dataObj)
 
      self.peakListTable.update(textMatrix=textMatrix, objectList=objectList)
  
    # Constraint List Table
    if CONSTRAINT_DATA in self.inputTypes:
      textMatrix = []
      objectList = []

      for i, dataObj in enumerate(constraintData):
         cSet = dataObj.nmrConstraintStore

         if cSet:
           constraintLists = dataObj.constraintLists or cSet.sortedConstraintLists()
         else:
           constraintLists = []
 
         for constraintList in constraintLists:
           if constraintList is None: # deleted
             datum = [i+1, DATA_MISSING,
                      None, None, None, None]

           else:
             datum = [i+1, constraintList.className[:-14],
                      cSet.serial,
                      constraintList.serial,
                      constraintList.name,
                      len(constraintList.constraints)]
 
           textMatrix.append(datum)
           objectList.append(dataObj)
 
      self.constraintListTable.update(textMatrix=textMatrix, objectList=objectList)
    
    # Params table
    
    if self.__class__ is  NmrCalcRunFrame:
      textMatrix = []
      objectList = []
 
      if run:
        objectList = [x for x in run.sortedRunParameters()
                      if x.ioRole=='input']
 
        for param in objectList:
          values = [param.booleanValue, param.intValue, param.floatValue, param.textValue]
          values = [v for v in values if v is not None]
 
          if not values:
            value = None
 
          elif len(values) == 1:
            value = values[0]
 
          else:
            value = ', '.join([str(v) for v in values])
 
          datum = param.data
          if datum:
            info = getDataObjText(datum)
            obj = '%s (%s)' % (info, datum.ioRole)
 
          else:
            obj = None
 
          datum = [param.serial, param.code, param.name, value, obj]
          textMatrix.append(datum)
 
      self.paramsTable.update(textMatrix=textMatrix, objectList=objectList)
  
    # Output table
    #headingList = ['Data Type','Name','CCPN Data','Parameters','Weight']
    
    textMatrix = []
    objectList = []
    
    if run:
      outputs = run.outputs
      sortedOutputs = []
      
      for dataObj in outputs:
        dataType = dataObj.className[:-4]
        serial = dataObj.serial
        sortedOutputs.append((dataType, serial, dataObj))
      
      sortedOutputs.sort()
      
      for dataType, serial, dataObj in sortedOutputs:
        dataObjText = getDataObjText(dataObj)
      
        params = []
        for param in dataObj.runParameters:
          text = param.code
          
          if param.name:
            text += ' "%s"'% param.name
      
          text += ' Value:'
          if param.booleanValue is not None:
            text += param.booleanValue and ' True' or ' False'
            
          if param.intValue is not None:
            text += ' %d' % param.intValue

          if param.floatValue is not None:
            text += ' %f' % param.floatValue

          if param.textValue is not None:
            text += param.textValue
          
          params.append(text)
      
        paramsText = ', '.join(params)
        
        textDatum = [dataType, serial,
                     dataObj.name,
                     dataObjText,
                     paramsText,]
                     #1.0]
                     #dataObj.weight]  # NBNB TBD DataOBj no longer has weight. To fix properly NBNB
        
        objectList.append(dataObj)
        textMatrix.append(textDatum)
    
      outputPars = [x for x in run.sortedRunParameters()
                    if x.ioRole=='output' and x.data is None and x.name != 'summary']
 
      objectList.extend(outputPars)
 
 
      for param in outputPars:
        values = [param.booleanValue, param.intValue, param.floatValue, param.textValue]
        values = [v for v in values if v is not None]
 
        if not values:
          value = None
 
        elif len(values) == 1:
          value = values[0]
 
        else:
          value = ', '.join([str(v) for v in values])
 
        #datum = param.data
        #if datum:
        #  info = getDataObjText(datum)
        #  obj = '%s (%s)' % (info, datum.ioRole)
 
        #else:
        #  obj = None
 
        datum = ['Parameter', param.serial, param.name, value, ""]
        textMatrix.append(datum)
 
    self.outputTable.update(textMatrix=textMatrix, objectList=objectList)
    # Output Summary
    
    if run:
      valueObj = run.findFirstRunParameter(ioRole='output', name='summary')
      if valueObj:
        self.outputSummary.setText(valueObj.textValue)
      else:
        self.outputSummary.setText()
    else:
      self.outputSummary.setText()
    
    self.after_idle(lambda:self.configure(cursor=""))
    self.waiting = False
  
  def updateRuns(self):
  
    names = []
    index = 0
    runs = []
    run = self.run
    
    if self.calcStore:
      runs = self.calcStore.sortedRuns()
      
      if runs:
        if run not in runs:
          run = runs[-1]
          
        index = runs.index(run)
        
        names = []
        for r in runs:
          if self.runIsEditable(r):
            names.append('%d' % r.serial)
          else:
            names.append('%d (uneditable)' % r.serial)
      
      else:
        run = None
    
    else:
      run = None
      
    if run is not self.run:
      self.changeRun(run)      
  
    self.runPulldown.setup(names, runs, index) 
  
  def changeRun(self, run):
    
    if run and (run is not self.run):
      self.run = run
      self.update(run.parent)
  
  def selectModel(self, obj, row, col):
  
    self.structureData = obj
  
  def updateEnsembles(self, obj=None):
  
    names = []
    index = 0
    ensembles = []
    
    if self.run:
      for ensemble in self.project.sortedStructureEnsembles():
        eId = '%s:%d' % (ensemble.molSystem.code, ensemble.ensembleId)
        names.append(eId)
        ensembles.append(ensemble)
  
    self.ensemblePulldown.setup(names, ensembles, index)

  def addEnsemble(self):
  
    ensemble = self.ensemblePulldown.getObject()
    run = self.run
    
    if run and ensemble:
      msCode = ensemble.molSystem.code
      eId = ensemble.ensembleId
      
      if not run.findFirstData(className=STRUCTURE_DATA,
                               molSystemCode=msCode,
                               ensembleId=eId,
                               ioRole='input'):
                                                
        run.newStructureEnsembleData(molSystemCode=msCode,
                                     ensembleId=eId,
                                     ioRole='input')

  def updateModels(self, obj=None):
  
    names = []
    index = 0
    models = []
    
    run = self.run
    if run:
      for data in run.findAllInputs(className=STRUCTURE_DATA):
        ensemble = data.structureEnsemble
        
        if not ensemble:
          continue
          
        eId = ensemble.ensembleId
        for model in ensemble.sortedModels():
          serial = model.serial
        
          if serial not in data.modelSerials:
            names.append('%s:%d' % (eId, serial))
            models.append(model)
  
    self.modelPulldown.setup(names, models, index)
  
  def addModel(self):
  
    model = self.modelPulldown.getObject()
    
    if self.run and model:
      ensemble = model.structureEnsemble
      
      data = self.run.findFirstData(className=STRUCTURE_DATA,
                                    structureEnsemble=ensemble,
                                    ioRole='input')
    
      if not data:
        data = self.run.newStructureEnsembleData(structureEnsemble=ensemble,
                                                 ioRole='input')
      
      serials = set(data.modelSerials)
      serials.add(model.serial)
      serials = list(serials)
      serials.sort()
      
      data.modelSerials = serials
    
  def removeModel(self):
  
    if self.structureData:
      for dataObj in self.structureTable.currentObjects:
        dataObj.delete() 
  
  def selectConstraintData(self, obj, row, col):
  
    self.constraintData = obj
  
  def getConstraintLists(self):
  
    cLists = []
    
    nmr = self.project.currentNmrProject
    for cSet in nmr.sortedNmrConstraintStores():
      for cList in cSet.sortedConstraintLists():
        cLists.append(cList)
    
    return cLists    
  
  def updateConstraintLists(self, obj=None):
  
    names = []
    index = 0
    cLists = []
    
    run = self.run
    
    if run:          
      for cList in self.getConstraintLists():
        cSet = cList.nmrConstraintStore.serial
        data = run.findFirstData(className=CONSTRAINT_DATA,
                                 constraintStoreSerial=cSet,
                                 constraintListSerials=(cList.serial,),
                                 ioRole='input')
      
        if not data:
          cType = cList.className[:-14]
          cSet = cList.nmrConstraintStore.serial
          names.append('%s-%s:%d' % (cType, cSet, cList.serial))
          cLists.append(cList)
   
    self.constraintListPulldown.setup(names, cLists, index) 

  def addConstraintList(self):
  
    constraintList = self.constraintListPulldown.getObject()
    if self.run and constraintList:
      cSet = constraintList.nmrConstraintStore.serial
      data = self.run.findFirstData(className=CONSTRAINT_DATA,
                                    constraintStoreSerial=cSet,
                                    constraintLists=[constraintList,],
                                    ioRole='input')
      
      if not data:
        cSet = constraintList.nmrConstraintStore.serial
        
        self.run.newConstraintStoreData(constraintStoreSerial=cSet,
                                        constraintListSerials=[constraintList.serial],
                                        ioRole='input')
    
  def removeConstraintList(self):
  
    if self.constraintData:
      for dataObj in self.constraintListTable.currentObjects:    
        dataObj.delete() 
        
        
  def selectPeakListData(self, obj, row, col):
  
    self.peakListData = obj
  
  def getPeakLists(self):
  
    peakLists = []
    
    nmr = self.project.currentNmrProject
    for experiment in nmr.sortedExperiments():
      for spectrum in experiment.sortedDataSources():
        for peakList in spectrum.sortedPeakLists():
          peakLists.append(peakList)
    
    return peakLists
  
  def updatePeakLists(self, obj=None):
  
    names = []
    index = 0
    peakLists = []
    
    run = self.run
    
    if run:          
      for peakList in self.getPeakLists():
        data = run.findFirstData(className=PEAK_DATA,
                                 peakList=peakList,
                                 ioRole='input')
      
        if not data:
          spectrum = peakList.dataSource
          sName = spectrum.name
          eName = spectrum.experiment.name
          names.append('%s:%s:%d' % (eName, sName, peakList.serial))
          peakLists.append(peakList)
   
    self.peakListPulldown.setup(names, peakLists, index)
    
  def addPeakList(self):
  
    peakList = self.peakListPulldown.getObject()
    if self.run and peakList:
      data = self.run.findFirstData(className=PEAK_DATA,
                                    peakList=peakList, ioRole='input')
      
      if not data:
        spectrum = peakList.dataSource
        eSerial = spectrum.experiment.serial
        sSerial = spectrum.serial
        pSerial = peakList.serial
        
        self.run.newPeakListData(experimentSerial=eSerial,
                                 dataSourceSerial=sSerial,
                                 peakListSerial=pSerial,
                                 ioRole='input')
 
        
  def removePeakList(self):
  
    if self.peakListData:
      for dataObj in self.peakListTable.currentObjects:
        dataObj.delete() 
     
  def selectMeasurementList(self, obj, row, col):
  
    self.measurementData = obj
  
  def updateMeasurementLists(self, obj=None):
  
    names = []
    index = 0
    mLists = []
    run = self.run
    
    if run:
      nmr = self.project.currentNmrProject
      
      for mList in nmr.sortedMeasurementLists():
        data = self.run.findFirstData(className=MEASUREMENT_DATA,
                                      measurementList=mList,
                                      ioRole='input')
      
        if not data:
          mLists.append(mList)
          names.append('%s:%d' % (mList.className, mList.serial))

    self.measurementListPulldown.setup(names, mLists, index)
    
  def addMeasurementList(self):
  
    mList = self.measurementListPulldown.getObject()
    if self.run and mList:
      data = self.run.findFirstData(className=MEASUREMENT_DATA,
                                    measurementList=mList,
                                    ioRole='input')
      
      if not data:
        self.run.newMeasurementListData(measurementListSerial=mList.serial,
                                       ioRole='input')
  
  def removeMeasurementList(self):
  
    if self.measurementData:
      for dataObj in self.measurementsTable.currentObjects:
        dataObj.delete() 
  
  def changeDetails(self, event):
  
    if self.run:
      value = self.detailsEntry.get().strip() or None
      if value != self.run.details:
        self.run.details = value
  
  def changeWmsProtocol(self, event):
  
    if self.run:
      value = self.protocolEntry.get().strip() or None
      if value != self.run.wmsProtocolName:
        self.run.wmsProtocolName = value
  
  def updateDetails(self):
  
    if self.run:
      self.detailsEntry.set(self.run.details)  
  
  def updateProtocol(self):
  
    if self.run:
      self.protocolEntry.set(self.run.wmsProtocolName)  
  
  def changeMolSystem(self, molSystem):
    
    if self.run and molSystem:
      data = self.run.findFirstData(className='MolSystemData',
                                    ioRole='input')
      
      if data and (molSystem is not data.molSystem):
        data.delete()
        data = None
      
      if not data:                                       
        self.run.newMolSystemData(molSystemCode=molSystem.code,
                                  ioRole='input')
  
      self.updateChainPulldown()
  
  def changeChain(self, chainCodes):
  
    if self.run and chainCodes:
      data = self.run.findFirstData(className='MolSystemData',
                                    ioRole='input')
      
      if data:
        data.chainCodes = chainCodes
      
  
  def updateChainPulldown(self, obj=None):
  
    if not (self.chainSelection and self.run):
      return
  
    names = []
    index = 0
    chainOpts = []
    
    data = self.run.findFirstData(className='MolSystemData',
                                 ioRole='input')
   
    if data:
      molSystem = data.molSystem
      
      chains = molSystem.sortedChains()
      
      stack = [[0,],]
      for c in [c.code for c in chains]:
        for s in stack[:]:
          t = s[:]
          s.append(c)
          s[0] -= 1
          stack.append(t)

      stack.sort()
       
      chainOpts = [x[1:] for x in stack[:-1]]
      names = [','.join(x) for x in chainOpts]
    
      codes = data.chainCodes
      if not codes:
        index = 0
        
      else:
        codes = list(codes)
        codes.sort()
        name = ','.join(codes)
        
        if name in names:
          index = names.index(name)
        else:
          index = 0
          
    self.chainPulldown.setup(names, chainOpts, index)
     
  
  def updateMolSystemPulldown(self, obj=None):
    
    names = []
    index = 0
    molSystems = []
    
    run = self.run
    if run:
      data = (run.findFirstData(className='MolSystemData', ioRole='input')
              or run.findFirstData(className='MolResidueData', ioRole='input'))
      
      if data:
        molSystem = data.molSystem
        molSystems = [molSystem,]
      else:
        molSystem = None
      
      if self.runIsEditable():
        molSystems = self.project.sortedMolSystems()
      
      if molSystems:
        names = [ms.code for ms in molSystems]
        
        if molSystem not in molSystems:
          molSystem = molSystems[0]
          
        index = molSystems.index(molSystem)  
        
      else:
        molSystem = None
        
      self.changeMolSystem(molSystem)
    
    self.molSystemPulldown.setup(names, molSystems, index)
      















