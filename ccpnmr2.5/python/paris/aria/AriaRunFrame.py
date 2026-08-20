import os, sys
import Tkinter

from memops.gui.Button          import Button
from memops.gui.DataEntry       import askString
from memops.gui.Entry           import Entry
from memops.gui.Frame           import Frame
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FileSelect      import FileType
from memops.gui.Label           import Label
from memops.gui.LabelFrame      import LabelFrame
from memops.gui.MessageReporter import showWarning, showOkCancel
from memops.gui.PulldownList    import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix


from ccp.util.NmrCalc import getObjBooleanParameter, toggleObjBooleanParameter, setRunParameter, getRunTextParameter

from ccpnmr.analysis.core.ExperimentBasic import getThroughSpacePeakLists

from ccpnmr.analysis.popups.EditCalculation import NmrCalcRunFrame, PEAK_DATA, CONSTRAINT_DATA
from ccpnmr.analysis.popups.BasePopup import BasePopup

PARAM_ATTR_DICT = {type(1.0):'floatValue',
                   type(1):'intValue',
                   type(True):'booleanValue'}

FILTER_VIOL = 'FilterViol'
KEEP_ASSIGN = 'KeepAssign'
AMBIG_PROTOCOL = 'AmbigProtocol'
USE_IN_CALC = 'UseInCalc'

CNS_EXE = 'CnsExe'
FILE_PREFIX = 'FilePrefix'
WORKING_DIR = 'WorkingDir'
TEMP_DIR = 'TempDir'

class AriaRunFrame(NmrCalcRunFrame):

  def __init__(self, parent, project, calcStore, *args, **kw):

    NmrCalcRunFrame.__init__(self, parent, project, calcStore,
                             inputTypes=(PEAK_DATA, CONSTRAINT_DATA),
                             chainSelection=True, *args, **kw)
    
    self.calcStore = calcStore
    self.optPeakList = None
    self.optConstraintList = None
        
    frame = self.tabbedFrame.frames[1]
    
    headingList = ['PeakList','Filter\nViolated?',
                   'Keep\nAssignments?','Shift List']
                   
    editWidgets      = [None, None, None, None]
    editGetCallbacks = [None, self.toggleFilterViol,
                        self.toggleKeepAssign, None]
    editSetCallbacks = [None, None, None, None]
                        
    row = 0
    subFrame = Frame(frame, grid=(row,0))
    subFrame.expandGrid(0,1)

    label = Label(subFrame, text='File Name Prefix:', grid=(0,0))
    self.filePrefixEntry = Entry(subFrame, text='aria', grid=(0,1), sticky="ew")
    self.filePrefixEntry.bind('<Leave>', self.updateEntryParams)
    
    label = Label(subFrame, text='CNS Executable:', grid=(1,0))
    self.cnsExeEntry = Entry(subFrame, text='', grid=(1,1), sticky="ew")
    self.cnsExeEntry.bind('<Leave>', self.updateEntryParams)
    button = Button(subFrame, text='Select File',bd=1,bg='#C0E0FF',
                    command=self.selectCnsExe, grid=(1,2))
    
    label = Label(subFrame, text='Working Directory:', grid=(2,0))
    self.workingDirEntry = Entry(subFrame, text='', grid=(2,1), sticky="ew")
    self.workingDirEntry.bind('<Leave>', self.updateEntryParams)
    button = Button(subFrame, text='Select File',bd=1,bg='#C0E0FF',
                    command=self.selectWorkingDir, grid=(2,2))
    
    label = Label(subFrame, text='Temporary Directory:', grid=(3,0))
    self.tempDirEntry = Entry(subFrame, text='', grid=(3,1), sticky="ew")
    self.tempDirEntry.bind('<Leave>', self.updateEntryParams)
    button = Button(subFrame, text='Select File',bd=1,bg='#C0E0FF',
                    command=self.selectTempDir, grid=(3,2))

    row += 1
    frame.expandGrid(row,0)
    self.grid_rowconfigure(row, weight=1)
    self.optPeakListMatrix = ScrolledMatrix(frame, headingList=headingList,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks, 
                                         editWidgets=editWidgets,
                                         multiSelect=True, grid=(row,0),
                                         callback=self.selectOptPeakList)

    row += 1
    frame.expandGrid(row,0)
    headingList = ['Constraint List',
                   'Name',
                   'Filter\nViolated?',
                   'Ambiguous\nProtocol?',]
                   
    editWidgets      = [None, None, None, None]
    editGetCallbacks = [None, None, self.toggleFilterViol, self.toggleAmbig]
    editSetCallbacks = [None, None, None, None]
                        
    self.optConstraintMatrix = ScrolledMatrix(frame, headingList=headingList,
                                               editSetCallbacks=editSetCallbacks,
                                               editGetCallbacks=editGetCallbacks,
                                               editWidgets=editWidgets,
                                               multiSelect=True, grid=(row,0),
                                               callback=self.selectConstraintList)

    self.optConstraintMatrix.doEditMarkExtraRules = self.doEditMarkExtraRules

    self.update(calcStore)

  def update(self, calcStore=None):
  
    NmrCalcRunFrame.update(self, calcStore)
    self.updateSettings()
     
    # Need to run NmrCalcRunFrame.update before
    # this to get new self.run, self.project 
    if self.run:
      repository = self.project.findFirstRepository(name='userData')
      defaultDir = repository.url.path
    
      filePrefix = getRunTextParameter(self.run, FILE_PREFIX) or self.project.name
      cnsExe = getRunTextParameter(self.run, CNS_EXE) or '/usr/bin/cns'
      workingDir = getRunTextParameter(self.run, WORKING_DIR) or defaultDir
      tempDir = getRunTextParameter(self.run, TEMP_DIR) or defaultDir

      self.filePrefixEntry.set(filePrefix)
      self.cnsExeEntry.set(cnsExe)
      self.workingDirEntry.set(workingDir)
      self.tempDirEntry.set(tempDir)

      self.updateEntryParams()

  def updateEntryParams(self, event=None):
  
    if self.run:
      repository = self.project.findFirstRepository(name='userData')
      defaultDir = repository.url.path
 
      filePrefix = self.filePrefixEntry.get() or self.project.name
      cnsExe = self.cnsExeEntry.get() or '/usr/bin/cns'
      workingDir = self.workingDirEntry.get() or defaultDir
      tempDir = self.tempDirEntry.get() or defaultDir
 
      setRunParameter(self.run, FILE_PREFIX, filePrefix)
      setRunParameter(self.run, CNS_EXE, cnsExe)
      setRunParameter(self.run, WORKING_DIR, workingDir)
      setRunParameter(self.run, TEMP_DIR, tempDir)

  def selectCnsExe(self):
    
    fileTypes = [ FileType("All", ["*"]), FileType("EXE", ["*.exe"]) ]

    popup = FileSelectPopup(self, fileTypes)

    file = popup.getFile()

    if file:
      self.cnsExeEntry.set( file )
    
    popup.destroy()
    self.updateEntryParams()
    
  def selectWorkingDir(self):
    
    popup = FileSelectPopup(self, show_file=False)

    directory = popup.getDirectory()
    if directory:
      self.workingDirEntry.set( directory )
    
    popup.destroy()
    self.updateEntryParams()

  def selectTempDir(self):
    
    popup = FileSelectPopup(self, show_file=False)

    directory = popup.getDirectory()
    if directory:
      self.tempDirEntry.set( directory )
    
    popup.destroy()
    self.updateEntryParams()
    
  def getPeakLists(self):
    # overwrites superclass
    
    return getThroughSpacePeakLists(self.project)

  def administerNotifiers(self, notifyFunc):

    NmrCalcRunFrame.administerNotifiers(self, notifyFunc)
     
    for func in ('__init__','delete','setName'):
      #notifyFunc(self.updateSettings, 'ccp.nmr.NmrConstraint.AbstractConstraintList', func)
      notifyFunc(self.updateSettings, 'ccp.nmr.NmrConstraint.DistanceConstraintList', func)
      notifyFunc(self.updateSettings, 'ccp.nmr.NmrConstraint.DihedralConstraintList', func)
      notifyFunc(self.updateSettings, 'ccp.nmr.NmrConstraint.HBondConstraintList', func)
        
  def selectOptPeakList(self, obj, row, col):
  
    self.optPeakList = obj

  def selectConstraintList(self, obj, row, col):
  
    self.optConstraintList = obj
 
  def doEditMarkExtraRules(self, obj, row, col):
  
    if col in (2,3):
      cSet = obj.nmrConstraintStore
      if cSet:
        for cList in obj.constraintLists:
          if cList.className[:-14] != 'Distance':
            # i.e. ambig protocols and viol filtering
            # is only for dist constraints, not dihedral etc
            return False
                 
    return True
  
  def toggleUseInCalc(self, dataObj):

    toggleObjBooleanParameter(dataObj, USE_IN_CALC)

    self.updateSettings()

  def toggleFilterViol(self, dataObj):

    toggleObjBooleanParameter(dataObj, FILTER_VIOL)

    self.updateSettings()
    
  def toggleKeepAssign(self, dataObj):
  
    toggleObjBooleanParameter(dataObj, KEEP_ASSIGN)

    self.updateSettings()
  
  def toggleAmbig(self, dataObj):
    
    toggleObjBooleanParameter(dataObj, AMBIG_PROTOCOL)

    self.updateSettings()
  
  def updateSettings(self, obj=None):


    textMatrix = []
    objectList = []
    
    run = self.run
    if run:      
      peakListsData = []
      constraintData = []
      molSystemData = None
      
      for dataObj in run.sortedInputs():
        className = dataObj.className
      
        if className == 'PeakListData':
          peakList = dataObj.peakList
          if peakList:
            peakListsData.append((dataObj, peakList))
        
        elif className == 'MolSystemData':
          molSystem = dataObj.molSystem
          if molSystem:
            molSystemData = (dataObj, molSystem)
          
        elif className == 'ConstraintStoreData':
          nmrConstraintStore = dataObj.nmrConstraintStore
          if nmrConstraintStore:
            constraintLists = dataObj.constraintLists or \
                              nmrConstraintStore.sortedConstraintLists()
            
            # Chould be only one
            for constraintList in constraintLists:
              if constraintList is None:
                # Prob happens when serial no longer exists
                continue
              constraintData.append((dataObj, constraintList))
            
      
      for dataObj, peakList in peakListsData:

        spectrum   = peakList.dataSource
        experiment = spectrum.experiment

        filterViol = getObjBooleanParameter(dataObj, FILTER_VIOL)
        keepAssign = getObjBooleanParameter(dataObj, KEEP_ASSIGN)
  
        shiftList  = peakList.dataSource.experiment.shiftList
        if shiftList:
          shiftListName = '%d:%s' % (shiftList.serial,shiftList.name)
        else:
          shiftListName = None
 
        ident = '%s:%s:%d' % (experiment.name,spectrum.name,peakList.serial)
 
        textMatrix.append([ident,
                           filterViol and 'Yes' or 'No',
                           keepAssign and 'Yes' or 'No',
                           shiftListName])
        objectList.append(dataObj)
 
      self.optPeakListMatrix.update(textMatrix=textMatrix,objectList=objectList)
 
      textMatrix = []
      objectList = [] 
      
      for dataObj, constraintList in constraintData:
       
        #useInCalc = getObjBooleanParameter(dataObj, USE_IN_CALC)
        cSet = dataObj.constraintStoreSerial
        cType = constraintList.className[:-14]

        if cType == 'Distance':
        
          filterViol = getObjBooleanParameter(dataObj, FILTER_VIOL)
          ambigProtocol = getObjBooleanParameter(dataObj, AMBIG_PROTOCOL)

          ambigProtocol = ambigProtocol and 'Yes' or 'No'
          filterViol = filterViol and 'Yes' or 'No'

        else:
          ambigProtocol = None
          filterViol = None
 
        ident = '%s - %d:%d' % (cType, cSet, constraintList.serial)

        textMatrix.append([ident,
                          constraintList.name,
                          filterViol,
                          ambigProtocol])
 
        objectList.append(dataObj)
           
    self.optConstraintMatrix.update(textMatrix=textMatrix,
                                     objectList=objectList)



