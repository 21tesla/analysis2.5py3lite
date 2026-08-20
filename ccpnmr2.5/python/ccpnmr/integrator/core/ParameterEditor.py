from memops.gui.BasePopup import BasePopup
from memops.gui.ButtonList import ButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.Label import Label
from memops.gui.PulldownList import PulldownList
from memops.gui.RadioButton import RadioButton
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame
from ccpnmr.analysis.core.ExperimentBasic import getNoesyPeakLists

import Io as intIo
from ccpnmr.integrator.core import Util as intUtil

class ParameterTable(ScrolledMatrix):
  """ helper class to deal with interfaceParameters with hicard != 1
      which contain a list of objects to be selected (displayed in a 2-column table) """

  def __init__(self, guiParent, dataText, dataObjects, dataValues, **kw):

    dataValueDict = {}
    for n, dataObject in enumerate(dataObjects):
      dataValueDict[dataObject] = dataValues[n]
    self.dataValueDict = dataValueDict
    self.isUsedSet = set() # default is that none selected

    editWidgets = 2*[None]
    editSetCallbacks = 2*[None]
    editGetCallbacks = [None, self.toggleIsUsed]
    headings = [dataText, 'Is Used?']

    ScrolledMatrix.__init__(self, guiParent, headingList=headings,
                            objectList=dataObjects, editWidgets=editWidgets,
                            editGetCallbacks=editGetCallbacks,
                            editSetCallbacks=editSetCallbacks, **kw)

    self.doUpdate()

  def toggleIsUsed(self, dataObject):

    if dataObject in self.isUsedSet:
      self.isUsedSet.remove(dataObject)
    else:
      self.isUsedSet.add(dataObject)

    self.doUpdate()

  def doUpdate(self):

    objectList = []
    textMatrix = []
    for dataObject in self.objectList:
      dataValue = self.dataValueDict[dataObject]
      isUsed = 'Yes' if dataObject in self.isUsedSet else 'No'
      textMatrix.append([dataValue, isUsed])
      objectList.append(dataObject)

    self.update(objectList=objectList, textMatrix=textMatrix)

class ParameterSelector(Frame):
  """ helper class to deal with interfaceParameters with hicard != 1
      where have a selector for some object and an entry widget  """
  
  def __init__(self, guiParent, dataText, dataObjects, dataValues, helpText, selectionText, **kw):
    
    Frame.__init__(self, guiParent, **kw)
        
    label = Label(self, text=dataText+':', grid=(0,0))
    self.selectionList = PulldownList(self, objects=dataObjects, texts=dataValues, grid=(0,1))
    label = Label(self, text=helpText, grid=(1,0), gridSpan=(1,2))
    label = Label(self, text=selectionText+':', grid=(2,0))
    self.selectionEntry = Entry(self, grid=(2,1))

class ParameterEditor(BasePopup):
  """ editor for InterfaceParameters """
  
  # note that if you change what widget class a given paramType uses then the code
  # needs modifying in quite a few places
  
  def __init__(self, guiParent, run, protocolInterface, *args, **kw):
    
    self.run = run
    self.protocolInterface = protocolInterface
    self.widgetDict = {} # maps interfaceParameter --> widget
    
    BasePopup.__init__(self, guiParent, *args, **kw)
    
  def body(self, guiFrame):
    
    self.geometry('700x400')
  
    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)
   
    interfaceGroups, interfaceGroupDict = self.getInterfaceGroups()
    
    if len(interfaceGroups) > 1:
      mainFrame = TabbedFrame(guiFrame, options=interfaceGroups, grid=(0,0))
      frames = mainFrame.frames
    else:
      mainFrame = Frame(guiFrame, grid=(0, 0))
      frames = [mainFrame]
      
    for n, interfaceGroup in enumerate(interfaceGroups):
      frame = frames[n]
      
      interfaceObjects = interfaceGroupDict[interfaceGroup]
      for interfaceObject in interfaceObjects:
        self.makeInterfaceWidget(frame, interfaceObject)
    
    texts = ('Ok', 'Cancel')
    commands = (self.setupRun, self.cancelRun)
    buttonList = ButtonList(guiFrame, texts=texts, commands=commands, grid=(1,0))
    
  def makeInterfaceWidget(self, frame, interfaceObject):

    clazz = self.getInterfaceWidgetClass(interfaceObject)
    
    if clazz:
      kw = self.getInterfaceParameterDict(clazz, interfaceObject)
      widget = clazz(frame,
                     grid=(interfaceObject.row, interfaceObject.col),
                     gridSpan=(interfaceObject.rowspan, interfaceObject.colspan), **kw)
      frame.grid_rowconfigure(interfaceObject.row, weight=1)
      frame.grid_columnconfigure(interfaceObject.col, weight=1)
      if interfaceObject.className == 'InterfaceParameter':
        self.widgetDict[interfaceObject] = widget
      
  def getInterfaceWidgetClass(self, interfaceObject):
    """ return widget class corresponding to this interfaceParameter """
    
    if interfaceObject.className == 'InterfaceLabel' or not interfaceObject.isEditable:
      return Label
      
    protocolParameter = interfaceObject.protocolParameter

    if protocolParameter.container: # TBD: for now ignore these parameters
      return None
      
    paramType = protocolParameter.paramType
    
    if paramType == 'intValue':
      clazz = IntEntry
    elif paramType == 'floatValue':
      clazz = FloatEntry
    elif paramType == 'textValue':
      clazz = Entry
    elif paramType == 'booleanValue':
      clazz = CheckButton
    elif paramType in ('molResidues',):
      clazz = ParameterSelector
    elif protocolParameter.hicard == 1:
      clazz = PulldownList
    elif paramType in ('peakList',):
      clazz = ParameterTable
    else:
      clazz = ParameterSelector
      
    return clazz
      
  def getInterfaceParameterDict(self, clazz, interfaceObject):
    """ set up arguments for widget class constructor for this interfaceParameter """
    
    kw = {}
    
    if interfaceObject.className == 'InterfaceLabel':
      kw['text'] = interfaceObject.label
      return kw
      
    protocolParameter = interfaceObject.protocolParameter
    paramType = protocolParameter.paramType
    
    if paramType in ('intValue', 'floatValue', 'textValue', 'booleanValue', ):

      defaults = protocolParameter.defaultStrings
      if defaults:
        defValue = intUtil.paramTypeConverters[paramType](defaults[0])
        if defValue:
          key = 'isSelected' if paramType == 'booleanValue' else 'text'
          kw[key] = defValue
        
    else:
    
      project = self.protocolInterface.root
      nmrProject = project.currentNmrProject
      
      if paramType == 'molResidues':
        dataText = 'Select MolSystem'
        dataObjects = molSystems = project.molSystems
        dataValues = [molSystem.code for molSystem in molSystems]
        helpText = 'Residue specification style: A:;B:-4,6,9,11-13,45-;C:'
        selectionText = 'Residue selection'
        
      elif paramType == 'measurementList':
        dataObjects = nmrProject.measurementLists
        dataValues = [dataObject.serial for dataObject in dataObjects]
      
      elif paramType == 'peakList':
        dataText = 'Peak List'
        dataObjects = []
        dataValues = []
        experiments = []
        for experiment in nmrProject.experiments:
          for expTransfer in experiment.expTransfers:
            if expTransfer.transferType == 'through-space':
              experiments.append(experiment)
        for experiment in experiments:
          experimentName = experiment.name or experiment.serial
          for dataSource in experiment.dataSources:
            dataSourceName = dataSource.name or dataSource.serial
            peakLists = dataSource.peakLists
          dataObjects.extend(peakLists)
          dataValues.extend(['%s:%s:%s' % (experimentName, dataSourceName, peakList.serial) for peakList in peakLists])
            
      elif paramType == 'molSysChains':
        dataText = 'Select MolSystem'
        dataObjects = molSystems = project.molSystems
        dataValues = [molSystem.code for molSystem in molSystems]
        helpText = 'Chain code style: A,B,D'
        selectionText = 'Chain codes'
          
      elif paramType == 'structureModels':
        dataText = 'Select Structure'
        dataObjects = structureEnsembles = project.structureEnsembles
        dataValues = ['%s:%s' % (structureEnsemble.molSystem.code, structureEnsemble.ensembleId) for structureEnsemble in structureEnsembles]
        helpText = 'Model number style: -4,6,9,11-13,45-'
        selectionText = 'Model selection'
          
      elif paramType == 'constraintLists':
        dataText = 'Restraint Set'
        dataObjects = nmrConstraintStores = nmrProject.nmrConstraintStores
        dataValues = [str(nmrConstraintStore.serial) for nmrConstraintStore in nmrConstraintStores]
        helpText = 'Constraint List number style: -4,6,9,11-13,45-'
        selectionText = 'Constraint List'
      
      else:
        raise Exception("Unknown parameter type: %s for %s" %
                        (paramType, protocolParameter))
      
      # sort them by their values    
      xx = zip(dataValues, dataObjects)
      xx.sort()
      dataValues, dataObjects = zip(*xx)
      
      if clazz is ParameterTable:
        kw['dataText'] = dataText
        kw['dataObjects'] = dataObjects
        kw['dataValues'] = dataValues
      elif clazz is ParameterSelector:
        kw['dataText'] = dataText
        kw['dataObjects'] = dataObjects
        kw['dataValues'] = dataValues
        kw['helpText'] = helpText
        kw['selectionText'] = selectionText
      else: # clazz is PulldownList
        kw['objects'] = dataObjects
        kw['texts'] = dataValues
      
    return kw
    
  def getInterfaceGroups(self):
    """ determine the interface groups, and what interfaceObjects are in each group """
    
    protocolInterface = self.protocolInterface
    
    interfaceGroups = [] # use list, not set, because possibly the order is intentional
    interfaceGroupDict = {}
    
    # could loop over protocolInterface.sortedInterfaceLabels() as well but that ought to be redundant in this context
    for interfaceObject in protocolInterface.sortedInterfaceParameters() + protocolInterface.sortedInterfaceLabels():
      interfaceGroup = interfaceObject.interfaceGroup or ''
      if interfaceGroup not in interfaceGroups:
        interfaceGroups.append(interfaceGroup)
        interfaceGroupDict[interfaceGroup] = []
      interfaceGroupDict[interfaceGroup].append(interfaceObject)

    return interfaceGroups, interfaceGroupDict
    
  def setupRun(self):
    """ setup the run when the user has clicked OK button """
    
    # TBD: for now ignores the "content" of a protocol parameter
    
    run = self.run
    protocolInterface = self.protocolInterface
    
    for interfaceParameter in protocolInterface.sortedInterfaceParameters():
      nmrCalcObject = self.makeNmrCalcData(run, interfaceParameter)
      
    self.close()
       
  def cancelRun(self):
    
    self.run = None
    self.close()
     
  def makeNmrCalcData(self, run, interfaceParameter, ioRole='input', dataObject=None):
    
    from ccpnmr.integrator.core import Io as intIo
    
    widget = self.widgetDict.get(interfaceParameter)
    if not widget:
      return
    
    protocolParameter = interfaceParameter.protocolParameter
    paramType = protocolParameter.paramType
    name = protocolParameter.name

    parDD = {'name':name, 'code':protocolParameter.code, 'ioRole':ioRole}
    
    result = None
    if paramType in ('intValue', 'floatValue', 'textValue', 'booleanValue', ):
      value = widget.get()
      if value:
        parDD[paramType] = value
        # check if alredy set, and overwrite:
        if protocolParameter.hicard == 1:
          # NBNB there could also be problem with hicard != 1,
          # but we ignore that for now
          result = run.findFirstRunParameter(name=name, data=dataObject)
        if result is None:
          result = run.newRunParameter(**parDD)
        else:
          # NB will cause trouble if the parameter had data set for wrong paramType,
          # but that should never happen
          for tag, val in parDD.items():
            setattr(result, tag, val)
        
    elif paramType == 'molResidues':
      molSystem = widget.selectionList.getObject()
      if molSystem:
        ss = ':;'.join(x.code for x in molSystem.sortedChains()) + ':'
        selector = widget.selectionEntry.get()
        if selector == ss:
          result = run.newMolResidueData(chains=molSystem.sortedChains(), **parDD)
        elif selector is not None:
          residues = intIo.parseResidueExpr(molSystem, selector)
          result = run.newMolResidueData(residues=residues, **parDD)
      
    elif paramType == 'measurementList':
      measurementList = widget.getObject()
      if measurementList:
        result = run.newMeasurementListData(measurementList=measurementList, **parDD)
        
    elif paramType == 'peakList':
      peakLists = widget.isUsedSet
      for peakList in peakLists:
        result = run.newPeakListData(peakList=peakList, **parDD)
        
    elif paramType == 'molSysChains':
      molSystem = widget.selectionList.getObject()
      if molSystem:
        ss = ','.join(x.code for x in molSystem.sortedChains())
        selector = widget.selectionEntry.get()
        if selector == ss:
          result = run.newMolSystemData(chains=molSystem.sortedChains(), **parDD)
        elif selector is not None:
          chainCodes = [x.strip() for x in selector.split(',')]
          result = run.newMolSystemData(molSystemCode=molSystem.code,
                                        chainCodes=chainCodes, **parDD)
      
    elif paramType == 'structureModels':
      structureEnsemble = widget.selectionList.getObject()
      if structureEnsemble:
        maxModels = max(x.serial for x in structureEnsemble.models)
        ss = '1-%d' % maxModels
        selector = widget.selectionEntry.get()
        if selector == ss:
          result = run.newStructureEnsembleData(
                          models=structureEnsemble.sortedModels(), **parDD)
        elif selector is not None:
          modelNums = intIo.parseNumberExpr(selector, startat=1, endat=maxModels)
          result = run.newStructureEnsembleData(structureEnsemble=structureEnsemble,
                                                modelSerials=sorted(modelNums),
                                                **parDD)
      
    elif paramType == 'constraintLists':
      nmrConstraintStore = widget.selectionList.getObject()
      if nmrConstraintStore:
        ss = ','.join(str(x.serial) for x in nmrConstraintStore.sortedConstraintLists())
        selector = widget.selectionEntry.get()
        if selector == ss:
          result = run.newConstraintStoreData(
                          constraintLists=nmrConstraintStore.sortedConstraintLists(),
                          **parDD)
        else:
          maxClists = max(x.serial for x in nmrConstraintStore.constraintLists)
          clistSerials = intIo.parseNumberExpr(selector, startat=1, endat=maxClists)
          result = run.newConstraintStoreData(nmrConstraintStore=nmrConstraintStore,
                                              constraintListSerials=clistSerials,
                                              **parDD)
                                              
    else:
      raise Exception("Unknown parameter type: %s for %s" %
                      (paramType, protocolParameter))
    
    return result