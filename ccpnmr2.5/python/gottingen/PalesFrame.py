
import subprocess
import os

from memops.universal import Io as uniIo

from memops.general.Util import copySubTree

from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Frame import Frame
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.PulldownList import PulldownList
from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.Entry import Entry
from memops.gui.Label import Label
from memops.gui.LabelDivider import LabelDivider
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import  showOkCancel, showWarning
from memops.gui.ScrolledText import ScrolledText
from memops.gui.ScrolledMatrix  import ScrolledMatrix 
from memops.gui.WebBrowser  import WebBrowser

from ccp.api.nmr import NmrCalc

from gottingen import Io as palesIo

# Max length for file name shown in table (longer are truncated)
maxFileNameLength = 30

#

# Pales interface definitions

# Standard PDB input file name:
palesStrucIn = 'ref.pdb'
# Standard PDB output file name:
palesStrucOut = 'rot.pdb'
# Standard RDC input file name:
palesDcIn ='dObs.tab'
# Standard RDC output file name:
palesDcOut ='dCalc.tab'

# values are: Name of NmrCalc clas sto use, name of attribute to set, 
#            hicard of attribute

def structureEnsembleLabel(structureEnsembleData):
  modelSerials = sorted(structureEnsembleData.modelSerials)
  if not modelSerials:
    ss = 'All'
  elif len(modelSerials) == 1:
    ss = str(modelSerials[0])
  else:
    ss = str(modelSerials)
  #
  return ':'.join([structureEnsembleData.molSystemCode,
                   str(structureEnsembleData.ensembleId), ss])

def constraintListLabel(constraintStoreData):
  constraintListSerials = sorted(constraintStoreData.constraintListSerials)
  if not constraintListSerials:
    ss = 'All'
  elif len(constraintListSerials) == 1:
    ss = str(constraintListSerials[0])
  else:
    ss = str(constraintListSerials)
  
  lists = constraintStoreData.constraintLists()
  classNames = set(x.className for x in constraintStoreData.constraintLists())
  if len(lists) < len(constraintListSerials) or not classNames:
    className = 'Missing'
  elif len(classNames) == 1:
    className = classNames.pop()
  elif classNames:
     className = 'Multiple'
  
  #
  return ('%s:%s (%s)' % 
          (constraintStoreData.constraintStoreSerial, ss, className))

def measurementListLabel(measurementListData):
  
  xx = measurementListData.measurementListSerial
  target = measurementListData.measurementList
  if target is None:
    result = '%s (Missing)' % xx
  else:
    result = '%s (%s)' % (xx, target.className)
  #
  return result

NmrCalcDataTypes = {
'Int': ('RunParameter', 'intValue', 1, None),
'Float': ('RunParameter', 'floatValue', 1, None),
'String': ('RunParameter', 'textValue', 1, None),
'Boolean': ('RunParameter', 'booleanValue', 1, None),
'CommandSwitch': ('RunParameter', 'code', 1, None),
'StringEnum': ('RunParameter', 'textValue', 1, None),
'Structure': ('StructureEnsembleData', 'models', -1, structureEnsembleLabel),
'ConstraintList': ('ConstraintStoreData', 'constraintLists', -1, 
                   constraintListLabel),
'MeasurementList': ('MeasurementListData', 'measurementList', 1, 
                   measurementListLabel),
}
# Program parameters. Format is:
# guiLabel:(code, dataType, defaultValue, description)
#   code is the name used by the program. code == None is used for switches; 
#   The value is the switch text.  
#   CommandSwitch values are used for sets of mutually exclusive switches (e.g. err/noerr)
#   dataType is the data type as given in GeneralData, with special extra values

progParameters = {
'Pales Command': (None, 'CommandSwitch', 
                  (('Steric alignment', 'stPales'),
                   ('Non-protein steric alignment', 'stPalesFree'), 
                   ('Electrostatic alignment', 'elPales'),
                   ('Fit RDC', 'bestFit'), 
                   ('Add Structural Noise', 'struc'), 
                   ('Order matrix arithmetic', 'anA'), 
                   ('RDC arithmetic', 'anDc'), 
                   ('Predict from structure', 'anPdb'), 
                   ('Histogram RDC analysis', 'daHist'), 
                   ('Max likelihood RDC analysis', 'daMl'), 
                   ('RDC file conversion', 'conv'),), 
                  'Calculation type'
                 ),
'Fit Type': (None, 'CommandSwitch', 
             (('Single fit', 'single'), ('Exhaustive fit', 'exhaust')), 
             'Single/Exhaustive fit'
            ),
'Structure': ('pdb', 'Structure', palesStrucIn, 'Input Structure'),
'RDC List': ('inD', 'RdcConstraintList', palesDcIn, 'Input RDC Constraint list'),
'Error Analysis': (None, 'CommandSwitch', 
                   (('None', None), 
                    ('JackKnifing', 'jack'), 
                    ('MonteCarloRDC', 'mcDc'), 
                    ('MonteCarloStruc', 'mcStruc'), 
                   ), 
                   'Error estimation mode'
                  ),
'Number of SVD': ('map', 'NonNegInt', 0, 
                  'Number of Best-Fits in Error Mapping.'),
'Rotation Index': ('rotID', 'NonNegInt', 0, 
                   'Index of Rotation of PDB'),
'Fixed Distances': (None, 'CommandSwitch', (('Yes', 'fixedDI'), ('No', 'nofixedDI')), 
                    'Use fixed distances'),
'RDC output format': (None, 'CommandSwitch', (('Full', 'verb'), 
                      ('RDC Input', 'noverb')), 'Format for RDC output file'),
'Q calc mode': (None, 'CommandSwitch', (('Tensor', 'qDa'), ('RDC Rms', 'qRms')), 
                'Data used to calculate Q'),
'Swap z and y axes': (None, 'CommandSwitch', (('No', 'noyzInv'), ('Yes', 'yzInv')), 
                      'Exchange Ordering of y and z axis?'),
'Jackknife fraction RDC': ('dFrac', 'Float', 0.2, 'Fraction of Total DC Data.'),
'Jackknife fixed # RDC': (None, 'CommandSwitch', (('No', 'nonDfixed'), 
                          ('Yes', 'nDfixed')), 
                          'Use Fixed # of RDC in Jackknifing?'),
'RDC MC error auto': (None, 'CommandSwitch', (('Yes', 'autoMc'), ('No', 'noautoMc')),
                      'Automatically Adjust Scaling of DC Errors?'),
'RDC MC error check': (None, 'CommandSwitch', (('Yes', 'err'), ('No', 'noerr')),
                       'Skip predictions outside experimental errors'),
'RDC MC fract valid SVD': ('fracAdjustMc', 'Float', 0.7, 
                           'Required Fraction of Self-consistent Best-Fits.'),
'RDC MC fix scaling': ('kLAdjustMc', 'Float', 1.0, 
                       'Starting Scaling Factor of DC Errors'),
'Struct MC cone auto': (None, 'CommandSwitch', (('Yes', 'autoCone'), 
                       ('No', 'noautoCone')), 
                       'Automatically Adjust Structural Noise?'),
'Struct MC cone min': ('lCone', 'Float', 0.0, 
                      'Starting value of Structrual Noise.'),
'Struct MC cone incr': ('incCone', 'Float', 0.2, 
                        'Increment of Structrual Noise.'),
#'Struct MC cone fix': (None, 'Float', '', ''),
'Rdc List Out': ('outD', 'RdcList', palesDcOut, 'Dipolar Coupling Output'),
'Structure Out': ('pdbRotF', 'Structure', palesStrucOut, 'Structure Output'),
'DC Fail Output': ('outFail', 'File', 'dFail.tab', ' Mapping Fail Statistics'),
'Da Spread Output': ('outDa', 'File', 'daR.tab', 'Da & R Distribution.'),
'Ang Spread Output': ('outAng', 'File', 'orient.tab', 'Angular Deviation from Average'),
'XYZ Spread Output': ('outMap', 'File', 'world.tab', 'World Map Coordinates '),
'Tensor Output': ('outA', 'File', 'saupe.tab', 'Predicted Saupe Matrix'),
}


  
  


palesModeList = []
palesModes = {}

stdOptionals = ['Rotation Index', 'Fixed Distances', 'RDC output format', 
                'Q calc mode', 'Swap z and y axes',]

# SVD mode
tag = 'SVD'
mode = palesModes[tag] = {
 'text':'SVD fit to RDC',
 'info':'''Back-calculation of rDCs by SVD:  Losonczi et al. (1999) J. Magn. Reson., 138, 334-342.
Analysis of errors in SVD: Zweckstetter & Bax (2001) J. Biomol. NMR, 23, 127-137.''',
 'fixedpar' : ['Pales Command', 'Fit Type',],
 'fixedvalues' : ['bestFit', 'single',],
 'mandatories' :
  ['Structure', 'RDC List', 'Error Analysis','Number of SVD',],
 'optionals': stdOptionals + 
  ['Jackknife fraction RDC', 'Jackknife fixed # RDC', 'RDC MC error auto', 
    'RDC MC error check', 'RDC MC fract valid SVD', 'RDC MC fix scaling', 
    'Struct MC cone auto', 'Struct MC cone min', 'Struct MC cone incr', 
    #'Struct MC cone fix', 
  ],
 'output':['Rdc List Out', 'Structure Out'], 
 'optionalOutput':
  [('DC Fail Output', 'File', 'daR.tab'), 
   ('Da Spread Output', 'File', 'orient.tab'), 
   ('Ang Spread Output', 'File', 'world.tab'), 
   ('XYZ Spread Output', 'File', 'dFail.tab'), 
   ('Tensor Output', 'File', 'saupe.tab'), 
   
  ],
}
palesModeList.append(tag)


tag = 'template'
#  mode
mode = palesModes[tag] = {
 'mode':'', 
 'text':'',
 'info':'''''',
 'fixedpar':
  [('', '')
  ],
 'fixedvalues':
  [('', '')
  ],
 'mandatories':
  [('', '', ''), 
  ],
 'optionals':
  [('', '', ''), 
  ],
 'output':
  [('', '', ''), 
  ],
 'optionalOutput':
  [('', '', ''), 
  ],
}
palesModeList.append(tag)

def runIsEditable(run):

  if run:
    
    if run.status != 'provisional':
      return False
    
    for datum in run.data:
      if datum.ioRole == 'output':
        return False
    
    for datum in run.runParameters:
      if datum.ioRole == 'output':
        return False
    
    else:
      return True

  return False


# Names of pales modes and corrseponding pales options
#palesModes = [ 
# ('SVD', 'bestFit', 'single'),
# ('exhSVD', 'bestFit', 'exhaust'),
# ('FixS', 'bestFit'), # neither single nor exhaust, but sets saupe
# ('LQmin', 'bestFit'), # neither single nor exhaust, but sets fixed/nofixed/dadr
# ('SSIA', 'stPales'),
# ('SSIAf', 'stPalesFree'),
# ('PALES', 'elPales'),
# ('PDB', 'anPdb'),
# ('Saupe', 'anA'),
# ('DC', 'anDc'),
# ('DaHist', 'daHist'),
# ('DaMl', 'daMl'),
# ('Noise', 'struc'),
# ('Convert', 'conv'),
#             ]

from memops.gui.FloatEntry      import FloatEntry
from memops.gui.IntEntry        import IntEntry
from memops.gui.Entry           import Entry
from memops.gui.PulldownList    import PulldownList

class GenericDataMatrix(ScrolledMatrix):
  """ Generic data table
  Mofified from HaddockFrame
  """
  
  headingColor  = '#80C080'
  
  def __init__(self, parent, progParameters, initialRows=10, canEdit=True,
               headingList=['Parameter','Value','Description'],
               justifyList=['center','center', 'left'], 
               *args, **kw):
    
    self.enumerations = {}
    self.progParameters = progParameters
    self.run = None
    self.paramList = []
    
    # Standard widget
    editWidgets      = [None, canEdit, None]
    editGetCallbacks = [None, self.getValue, None]
    editSetCallbacks = [None, self.setValue, None]
    
    ScrolledMatrix.__init__(self, parent, headingList=headingList,
		            justifyList=justifyList,
		            editWidgets=editWidgets,
		            editGetCallbacks=editGetCallbacks, 
		            editSetCallbacks=editSetCallbacks,
                            initialRows=initialRows,
		            multiSelect=False,
		            passSelfToCallback=True,
		            callback=self.selectRowObj,)
    
    # Generic widgets
    self.stringEntry  = Entry(parent, returnCallback=self.setValue)
    self.intEntry     = IntEntry(parent, returnCallback=self.setValue)
    self.floatEntry   = FloatEntry(parent, returnCallback=self.setValue)
    self.pulldownList = PulldownList(parent, callback=self.setValue, 
                                     initCallback=False)
    
    self.parInfoMap = {
     'Int': (self.intEntry, self.stdGetter, self.setInt, 'Int'),
     'NonNegInt': (self.intEntry, self.stdGetter, self.setNonNegInt, 'Int'),
     'PosInt': (self.intEntry, self.stdGetter, self.setPosInt, 'Int'),
     'NegInt': (self.intEntry, self.stdGetter, self.setNegInt, 'Int'),
     'Float': (self.floatEntry, self.stdGetter, self.setFloat, 'Float'),
     'PosFloat': (self.floatEntry, self.stdGetter, self.setPosFloat, 'Float'),
     'NonNegFloat': (self.floatEntry, self.stdGetter, self.setNonNegFloat, 
                     'Float'),
     'NegFloat': (self.floatEntry, self.stdGetter, self.setNegFloat, 'Float'),
     'Fraction': (self.floatEntry, self.stdGetter, self.setFraction, 'Float'),
     'String': (self.stringEntry, self.stdGetter, self.setString, 'String'),
     'StringEnum': (self.pulldownList, self.getEnum, self.setPulldownList, 
                    'StringEnum'),
     'CommandSwitch': (self.pulldownList, self.getEnum, self.setPulldownList, 
                    'CommandSwitch'),
     'Structure': (self.pulldownList, self.getStructure, 
                   self.setDataPulldownList, 'Structure'),
     'RdcConstraintList': (self.pulldownList, self.getConstraintList, 
                           self.setDataPulldownList, 'ConstraintList'),
    }
    
  def updateGeneric(self, **kw):
  
    """
    """
    if 'run' in kw:
      self.run = kw['run']
    if 'paramList' in kw:
      self.paramList = list(kw['paramList'])
    run = self.run
    paramList = self.paramList
    
    print '### updateGeneric', run, paramList
    
    self.enumerations.clear()
    if not run or not paramList:
      self.update(objectList=[], textMatrix=[], colorMatrix=[])
      return
    
    blank = [None] * 3
 
    textMatrix = []
    objectList = []
    colorMatrix = []
    
    progParameters = self.progParameters
    
    for label in paramList:
      tt = progParameters.get(label)
      if tt:
        code, parType, default, info = tt
        result = self.setupRunObject(run, label, code, parType, default, info)
        if result is None:
          continue
        textLine, objLine = result
        textMatrix.append(textLine)
        objectList.append(objLine)
        colorMatrix.append(blank)
      
      else:
        # No parameter description - treat as heading
        colorMatrix.append([self.headingColor]*3)
        textMatrix.append([label,None,None])
        objectList.append(None)
    
    # Finally update the ScrolledMatrix table
    self.update(objectList=objectList, textMatrix=textMatrix, colorMatrix=colorMatrix)
    
  
  def selectRowObj(self, obj, row, col, table):
  
    self.rowObj = obj

  def getValue(self, rowObj):
 
    # get correct object, widget and get/set functions
    obj, attrName, widget, getter, setter, data = rowObj

    # Bit of a hack because widgets are normally set at construction per column
    # Now setting it per row, and hence on-the-fly
    self.editWidget = widget
    
    if widget is self.pulldownList:
      widget._ccpn_data = data
  
    # get current value & setup widget
    getter(widget, obj, attrName)
    
  def setValue(self, event, null=None): # null is for pulldown menu callbacks passing name,index

    # get correct widget and get/set functions
    obj, attrName, widget, getter, setter, data = self.rowObj
  
    # set and check the appropriate parameter value from current edit widget

    # no setter for boolean toogles - the getter will have done all the toggling already
    # no setter for file selects - the file popup gets and sets the value and cannot be interrupted
    if setter: 
      setter(widget, obj, attrName)

    #self.updateGeneric() # Update table
  
  def getEnum(self, widget, obj, attrName):
    
    value = getattr(obj, attrName) 
    
    names,values = widget._ccpn_data
    if value in values:
      index = values.index(value)
    else:
      raise ValueError("value %s not in allowed values %s" % 
                       (value, values))
    
    widget.setup(names, values, index)
    del widget._ccpn_data
  
  def getStructure(self, widget, obj, attrName):
    
    value = getattr(obj, attrName) 
   
    names = ['<None>']
    values = [None]
    for ensemble in self.run.memopsRoot.sortedStructureEnsembles():
      code = ensemble.molSystem.code
      ensembleId = ensemble.ensembleId
      for model in ensemble.sortedModels():
        values.append(model)
        names.append('%s:%s:%s' % (code, ensembleId, model.serial))
    
    widget.setup(names, values, index)
  
  def getConstraintList(self, widget, obj, attrName):
    
    value = getattr(obj, attrName) 
   
    names = ['<None>']
    values = [None]
    for store in self.run.nmrCalcStore.nmrProject.sortedNmrConstraintStores():
      storeSerial = store.serial
      for constraintList in store.sortedConstraintLists():
        values.append(constraintList)
        names.append('%s:%s (%s)' % 
        (storeSerial, constraintList.serial, constraintList.className))
    
    widget.setup(names, values, index)
  
  def stdGetter(self, widget, obj, attrName):

    value = getattr(obj, attrName) 
    widget.set(value)
    
  def setString(self, widget, obj, attrName):

    value = widget.get().strip()

    if value: 
      setattr(obj, attrName, value)
          
  def setPulldownList(self, widget, obj, attrName):

    value = widget.getObject()
    setattr(obj, attrName, value)
          
  def setDataPulldownList(self, widget, obj, attrName):
    
    obj, attrName, widget, getter, setter, data = self.rowObj
    value = widget.getObject()
    textLine, objLine = self.setupRunObject(self.run, *data, value=value)
    self.rowObj = objLine
    
  def toggleBoolean(self, widget, obj, attrName):

    value = getattr(obj, attrName)
    setattr(obj, attrName, not value)

    #self.updateGeneric()
    
  def getDirectory(self, widget, obj, attrName):
  
    """Is called by updateGeneral. No widget because table placement is irrelevent and setting process not interruptable"""

    value = getattr(obj, attrName) 
    if not os.path.exists(value): value = None

    popup = FileSelectPopup(self, directory=value, show_file=False)
    value = popup.getDirectory()
    if value:
      setattr(obj, attrName, value)
      #self.updateGeneric()
            
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
      #self.updateGeneric()
          
  def setInt(self, widget, obj, attrName):

    value = widget.get() or 0
    setattr(obj, attrName, value)
    
  def setNonNegInt(self, widget, obj, attrName):

    value = widget.get() or 0
    if value >= 0: 
      setattr(obj, attrName, value)
          
  def setPosInt(self, widget, obj, attrName):

    value = widget.get() or 0
    if value > 0: 
      setattr(obj, attrName, value)
          
  def setNegInt(self, widget, obj, attrName):

    value = widget.get() or 0
    if value < 0: 
      setattr(obj, attrName, value)

  def setFloat(self, widget, obj, attrName):

    value = widget.get() or 0.0
    setattr(obj, attrName, value)

  def setFraction(self, widget, obj, attrName):
  
    value = min(1.0,max(0.0, widget.get() or 0.0))
    setattr(obj, attrName, value)

  def setPosFloat(self, widget, obj, attrName):

    value = widget.get() or 0.0
    if value > 0.0: 
      setattr(obj, attrName, value)

  def setNegFloat(self, widget, obj, attrName):

    value = widget.get() or 0.0
    if value < 0.0: 
      setattr(obj, attrName, value)

  def setNonNegFloat(self, widget, obj, attrName):

    value = widget.get() or 0.0
    if value >= 0.0: 
      setattr(obj, attrName, value)



  def setupRunObject(self, run, name, code, parType, default, description, 
                     **kwArgs):
    
    
    print '### setupRunObject', run, name, code, parType, default, description
    
    names = []
    values = []
    if isinstance(default, tuple) or isinstance(default, list):
      # Enumeration: set up. Can only happen for runParameters
      names, values = zip(*default)
    
    if 'value' in kwArgs:
      hasValue = True
    else:
      hasValue = False
    value = kwArgs.get('value')
    
    info = self.parInfoMap.get(parType)
    if info is None:
      print '### WARNING parameter type not implemented: ', parType
      return
    
    entry, getter, setter, dataType = info
    className, attrName, hicard, toString = NmrCalcDataTypes[dataType]
 
    if className == 'RunParameter':
      
      data = (names, values)
      
      # RunParameter
      obj = run.findFirstRunParameter(ioRole='input', name=name)
      if obj:
        if hasValue:
          setattr(obj, attrName, value)
      else:
        obj = run.newRunParameter(name=name, code=code, ioRole='input')
        if not hasValue:
          if values:
            value = values [0]
          else:
            value = default
        setattr(obj, attrName, value)
      valueText = getattr(obj, attrName)
 
    else:
      # Data
      
      data = (name, code, parType, default, description)
      
      obj = run.findFirstData(ioRole='input', name=name)
      if obj:
        # Data object already present
        if hasValue:
          setattr(obj, attrName, value)
        valueText = toString(obj)
      
      elif value is not None:
        # we are passing a value - create a new Data object
        dd = {attrName:value}
        func = getattr(run, 'new' + className)
        obj = func(name=name, code=code, ioRole='input', **dd)
        valueText = toString(obj)
      
      else:
        valueText = '<None>'
    
    if run.status != 'provisional':
      entry = setter = None
      
    textLine = [name, valueText, description]
    objectLine = [obj, attrName, entry, getter, setter, data]
    return  (textLine, objectLine)


class PalesFrame(Frame):
  """ Frame for handling PALES calculation.
  
  Note that the frame uses (or creates) an NmrCalcStore named 'PALES'
  and linked to the current NmrProject.
  """

  def __init__(self, parent, project, closeButton=False,
               *args, **kw):

    self.parent = parent
    self.project = project
    self.nmrProject = (project.currentNmrProject 
                       or project.newNmrProject(name='PALES'))
    self.calcStore = None
    self.run = None
    self.inpStructure = None
    self.inpConstraintList = None
    self.workingDir = None
    self.waiting = False
    
    self.palesMode = None
    
    self.resetCalcStore()
    
    Frame.__init__(self, parent, *args, **kw)
    
    
    self.expandGrid(0,0)
    
    options = ['Input Data','Extra Input', 'View Results']
      
    tabbedFrame = TabbedFrame(self, options=options, grid=(0,0))
    frameA, frameX, frameB = tabbedFrame.frames
    self.tabbedFrame = tabbedFrame
    
    label = Label(tabbedFrame.sideFrame, text='Run Number:',
                  grid=(0,0), sticky='e')
                  
    tipText = 'Selects which calculation job or "run" is currently being viewed or edited'
    self.runPulldown = PulldownList(tabbedFrame.sideFrame, 
                                    callback=self.changeRun,
                                    grid=(0,1), sticky='e', tipText=tipText)
    
    tipTexts = ['Delete the current calculation run settings']
    texts = ['Delete Run']
    commands = [self.deleteRun]
    
    if closeButton:
      ButtonListClass = UtilityButtonList
    else:
      ButtonListClass = ButtonList
    
    runButtons = ButtonListClass(tabbedFrame.sideFrame, texts=texts, 
                                 tipTexts=tipTexts,
                                 commands=commands, sticky='e', grid=(0,2))
    
    
    # Input data
    
    frameA.expandGrid(2,1)
    
    row = 0
    label = Label(frameA, text='Pales mode:',
                  grid=(row,0), sticky='w')
    self.palesModePulldown = PulldownList(frameA, callback=self.changePalesMode,
                                    grid=(row,1))
    
    tipTexts = ['Make a setup for a new calculation run',
                'Make a new calculation run by copying the current one',]
    texts = ['New Run', 'Copy Run']
    commands = [self.newRun, self.copyRun]
    
    if closeButton:
      ButtonListClass = UtilityButtonList
    else:
      ButtonListClass = ButtonList
    
    runButtons = ButtonList(frameA, texts=texts, tipTexts=tipTexts,
                                 commands=commands, sticky='e', grid=(row,2))
    runButtons.buttons[0].config(bg='#B0FFB0')
                            
    row += 1  
    subframe1 = LabelFrame(frameA, text='Description', grid=(row,0), 
                           gridSpan=(1,3))
    subframe1.expandGrid(0,1)
    self.modeDescription = Label(subframe1, grid=(row,0), sticky='w')

    row += 1         
    # setup generic table headings, justification and widget getters/setters
    
    self.inputMatrix = GenericDataMatrix(frameA, progParameters)
    self.inputMatrix.grid(row=row, column=0, columnspan=3, sticky='nsew')
    
                   
    row += 1 
    label = Label(frameA, text='Comments:', grid=(row,0))
    self.detailsEntryIn = Entry(frameA, grid=(row,1),  
                              gridSpan=(1,2), sticky="ew")
    self.detailsEntryIn.bind('<Leave>', self.changeDetailsIn)
                                    
    row += 1    
    button = Button(frameA, text='Select working dir:',bd=1,
                    command=self.selectWorkingDir, grid=(row,0), sticky="ew")
    self.workingDirEntry = Entry(frameA, text='.', grid=(row,1), gridSpan=(1,2),
                         width=48, sticky="ew", bd=1)
    
    row += 1 
    button = Button(frameA, text='Execute Pales:',bd=1,
                    command=self.executePales, grid=(row,0), gridSpan=(1,3), 
                    sticky="new")
    
    
    # Extra input
    # setup generic table headings, justification and widget getters/setters
    
    frameX.expandGrid(0,0)
    
    self.extraInputMatrix = GenericDataMatrix(frameX, progParameters)
    self.extraInputMatrix.grid(row=0, column=0, sticky='nsew')
                    
    
    # View Results
    
    frameB.expandGrid(7,1)
    
    row = 0
    
      
    subframe1 = LabelFrame(frameB, text='Command Options:', grid=(row,0), 
                           gridSpan=(1,4))
    #                       gridSpan=(1,2))
    subframe1.expandGrid(0,1)
    self.palesOptionsLabel = Label(subframe1, grid=(row,0), 
                                   sticky='w')
        
    row += 1  
    div = LabelDivider(frameB, text='Data', grid=(row,0), gridSpan=(1,4))
               
    row += 1          
    
    self.outputMatrix = GenericDataMatrix(frameB, progParameters, 
                                           initialRows=4)
    self.outputMatrix.grid(row=0, column=0, sticky='nsew')
    
    self.outputMatrix.grid(row=row, column=0, columnspan=(4), sticky='nsew')
                       
    
    row += 1     
    button = Button(frameB, text='View Selected',bd=1,
                    command=self.viewPalesData, grid=(row,0), gridSpan=(1,4), 
                    sticky="ew")
                    
    row += 1 
    label = Label(frameB, text='Comments:', grid=(row,0), sticky="w")
    self.detailsEntry = Entry(frameB, grid=(row,1),  
                              gridSpan=(1,3), sticky="ew")
    self.detailsEntry.bind('<Leave>', self.changeDetails)

    
    row += 1     
    subframe2 = LabelFrame(frameB, text='Calculated Order Matrix:', grid=(row,0), 
                           gridSpan=(1,4))
    #subframe2.grid_columnconfigure(5, weight=1)
    subframe2.expandGrid(1,5)
    label = Label(subframe2, text='Daxial', grid=(0,0), sticky='ew')
    label = Label(subframe2, text='Drhombic', grid=(0,1), sticky='ew')
    label = Label(subframe2, text='Psi', grid=(0,2), sticky='ew')
    label = Label(subframe2, text='Phi', grid=(0,3), sticky='ew')
    label = Label(subframe2, text='Theta', grid=(0,4), sticky='ew')
     
    self.outputTensorLabels = ll = []
    for ii in range(5):
      label = Label(subframe2, text='<None>', grid=(1,ii), sticky='ew')
      ll.append(label)

    row += 1     
    div = LabelDivider(frameB, text='Program Output', grid=(row,0), 
                       gridSpan=(1,4), sticky='sew')
                       
    #textFrame1 = LabelFrame(frameB, text='Pales Output File', grid=(5,0), 
    #                        gridSpan=(1,6), sticky='nsew')
    #textFrame1.expandGrid(0,0)
    
    row += 1 
    self.palesOutputText = ScrolledText(frameB, xscroll=False)

    self.palesOutputText.grid(row=row, column=0, columnspan=4, sticky='nsew')
    
    self.updateAfter()
    self.administerNotifiers(self.parent.registerNotify)
  
  
  def updatePalesOptions(self):
    
    palesMode =self.palesMode
    dd = palesModes.get(palesMode)
    # NBNB TBD
  
  
  def resetCalcStore(self, calcStore=None):
    """ Reset self.calcStore if missing or deleted
    """
    
    if calcStore is None or calcStore is self.calcStore:
      nmrProject = self.nmrProject
      calcStore= self.project.findFirstNmrCalcStore(nmrProject=nmrProject,
                                                    name='PALES')
      if calcStore is None:
        calcStore = self.project.newNmrCalcStore(nmrProject=self.nmrProject,
                                                 name='PALES')
        
    if self.calcStore is not calcStore:
      self.calcStore = calcStore
      if self.run is not None:
        self.run = None
        self.updateAfter()

  def selectWorkingDir(self):
    
    popup = FileSelectPopup(self, show_file=False)

    directory = popup.getDirectory()
    if directory:
      self.workingDirEntry.set(directory)
    
    popup.destroy()
    self.updateAfter()
    

  def changeDetails(self, event):
  
    if self.run:
      value = self.detailsEntry.get().strip() or None
      if value != self.run.details:
        self.run.details = value
        
    

  def changeDetailsIn(self, event):
  
    if self.run:
      value = self.detailsEntryIn.get().strip() or None
      if value != self.run.details:
        self.run.details = value
  
  def updateDetails(self):
    if self.run:
      text = self.run.details
    else:
      text = ''
    self.detailsEntry.set(text)  
    self.detailsEntryIn.set(text) 
    
  
  def changePalesMode(self, mode):
    
    if mode and (mode is not self.palesMode):
      run = self.run
      if runIsEditable(run):
        
        palesModeObj = run.findFirstRunParameter(name='Pales Mode')
        if palesModeObj is None:
          palesModeObj = run.newRunParameter(name='Pales Mode', ioRole='input')
 
        palesModeObj.textValue = mode
        self.palesMode = mode
        self.initialiseRun()
        self.updateAfter()
  
  def updatePalesModes(self, obj=None):
  
    names = []
    index = 0
    modes = []
    
    for tag in palesModeList:
      modeInfo = palesModes[tag]
      names.append(modeInfo['text'])
      modes.append(tag)
    
    if obj is not None:
      mode = obj
    else:
      mode = self.palesMode
    
    try:
      index = modes.index(mode)
    except ValueError:
      index = 0
    
    self.palesMode = modes[index]
    
    self.palesModePulldown.setup(names, modes, index) 
    
    tag = modes[index]
    self.modeDescription.set(palesModes[tag].get('info'))
  

  def copyRun(self):
    
    # TBD: Inputs only?
  
    if self.run:
      self.configure(cursor="watch")
      run = copySubTree(self.run, self.calcStore)
      run.status = 'provisional'
      
      # remove output
      for data in run.findAllData(ioRole='output'):
        data.delete()
      for parObj in run.findAllRunParameters(ioRole='output'):
        parObj.delete()
 
      self.run = run
      self.updateAfter()
      self.after_idle(lambda:self.configure(cursor=""))
  
  def newRun(self):

    if self.calcStore:
      self.run = self.calcStore.newRun(status='provisional')
      self.initialiseRun()
      self.updateAfter()
  
  def initialiseRun(self):
    """ Create objects for known parameters
    """
    
    print '### initialiseRun', self.run and runIsEditable(self.run)
    
    run = self.run
    if run is None or not runIsEditable(run):
      return
    
    # set up
    tag = self.palesMode
    modeData = palesModes[tag]
    currentTags = set(modeData['fixedpar'] + modeData['mandatories'] +
                      modeData['optionals'])
    
    # set palesMode if missing
    self.updatePalesModes()
    
    # remove parameters no longer needed
    for obj in run.runParameters:
      if obj.name not in currentTags:
        obj.delete()
    for obj in run.data:
      if obj.name not in currentTags:
        obj.delete()
    
    self.inputMatrix.updateGeneric(run=run, paramList=(modeData['fixedpar'] + 
                                                       modeData['mandatories']))
    self.extraInputMatrix.updateGeneric(run=run, 
                                        paramList=modeData['optionals'])
    
    for ii,name in enumerate(modeData['fixedpar']):
      value = modeData['fixedvalues'][ii]
      code, parType, default, info = progParameters[name]
      self.inputMatrix.setupRunObject(run, name, code, parType, default, info, 
                                      value=value)
    
    
  def getPalesOptions(self):
    run = self.run
    if run is None:
      return
    palesMode = self.palesMode
    
    return 'NBNB TBD'
  
  def changeRun(self, run):
    print '### changeRun', self.run, run
    
    if run and (run is not self.run):
      self.run = run
      self.updateAfter()
      

  def deleteRun(self):
  
    if self.run:
      msg = 'Really delete calculation run %d?' % self.run.serial
      
      if showOkCancel('Query', msg, parent=self):
        self.run.delete()
        self.run = None
        self.updateAfter()
     
  def updateRunsAfter(self, obj=None): 
    print '### updateRunsAfter', obj,  self.waiting
      
    if self.waiting:
      return
      
    else:
      self.waiting = True
      self.after_idle(self.updateRuns)    
  
  def updateRunDataAfter(self, obj=None):
    print '### updateRunDataAfter', obj
    
    if obj is None:
      run = None
    else:
      run = obj.run
    self.updateRunAfter(run=run)
  
  def updateRunAfter(self, run=None):    
    """ update all if run is curent run, otherwise update run pulldown only
    """
    print '### updateRunAfter', run,  self.waiting
      
    if self.waiting:
      return
    
    elif run is None:
      return
    
    elif run is self.run:
      self.waiting = True
      self.after_idle(self.update)   
    
    else:
      self.waiting = True
      self.after_idle(self.updateRuns)
 
  def updateAfter(self, obj=None):
    
    if self.waiting:
      return
    
    self.waiting = True
    self.configure(cursor="watch")
    self.after_idle(self.update)
  
  def updateRuns(self, run=None):
    """ Update run pulldown only
    """
    
    print '### updateRuns', run
  
    names = []
    index = 0
    runs = []
    run = self.run
    
    if self.calcStore is None or self.calcStore.isDeleted:
      self.resetCalcStore()
    
    runs = self.calcStore.sortedRuns()
    
    if runs:
      if run not in runs:
        run = runs[-1]
        
      index = runs.index(run)
      
      names = []
      for r in runs:
        if runIsEditable(r):
          names.append('%d' % r.serial)
        else:
          names.append('%d (uneditable)' % r.serial)
    
    else:
      run = None
      
    if run is not self.run:
      self.changeRun(run)      
    
    print '###', names, runs, index
    
    self.runPulldown.setup(names, runs, index) 
    self.waiting = False    
 
  def update(self, obj=None):
    print '### update', self.run
    
    
    self.updateRuns(obj)
    
    run = self.run
    if run is None:
      palesMode = self.palesMode = None
      self.palesOptionsLabel.set('')
      for label in self.outputTensorLabels:
        label.set('<None>')
      self.palesOutputText.setText()
      
    
    else:
      
      # Pales command NBNB TBD
      #valueObj = run.findFirstRunParameter(code='command')
      #if valueObj is None:
      #  self.palesCommandLabel.set('')
      #else:
      #  self.palesCommandLabel.set(valueObj.textValue)
      
      palesModeObj = run.findFirstRunParameter(ioRole='input', 
                                                  name='Pales Mode')
      if palesModeObj is None:
        palesMode = None
      else:
        palesMode = palesModeObj.textValue
      
      
      # Output orientation matrix
      keywords =('dAxialOut', 'dRhombicOut', 'psiOut', 'phiOut', 'thetaOut')
      for ii, label in enumerate(self.outputTensorLabels):
        keyword = keywords[ii]
        valueObj = run.findFirstRunParameter(code=keyword, ioRole='output')
        if valueObj is None:
          label.set('<None>')
        else:
          label.set(valueObj.floatValue)
      
      # Pales output text
      valueObj = run.findFirstRunParameter(ioRole='output', name='Output Text')
      if valueObj:
        self.palesOutputText.setText(valueObj.textValue)
      else:
        self.palesOutputText.setText()
      
    
    self.updatePalesModes(palesMode)
    
    self.palesOptionsLabel.set(self.getPalesOptions())
    
    #self.updateInputMatrix()
    #self.updateExtraInputMatrix()
    #self.updateOutputMatrix()
    self.inputMatrix.updateGeneric()
    self.extraInputMatrix.updateGeneric()
    self.outputMatrix.updateGeneric()
    self.updateDetails()
        
    self.after_idle(lambda:self.configure(cursor=""))
    self.waiting = False


  def administerNotifiers(self, notifyFunc):
    
    notifyFunc( self.updateRunsAfter, 'ccp.nmr.NmrCalc.NmrCalcStore', 'delete')
    notifyFunc( self.updateRunsAfter, 'ccp.nmr.NmrCalc.Run', '__init__')
    notifyFunc( self.updateRunAfter, 'ccp.nmr.NmrCalc.Run', 'delete')
    notifyFunc( self.updateRunAfter, 'ccp.nmr.NmrCalc.Run', 'setDetails')
    
    
    NC = 'ccp.nmr.NmrCalc.'
    nmrCalcClasses = [
    'RunParameter', 
    'EnergyTerm', 
    'ConstraintStoreData', 
    'ViolationListData', 
    'MolSystemData', 
    'MolResidueData', 
    'SpectrumData', 
    'PeakListData', 
    'SpinSystemData', 
    'StructureEnsembleData', 
    'ExternalData', 
    'FloatMatrixData', 
    'TensorData', 
    'MeasurementListData', 
    'DerivedListData', 
    ]
    for clazz in nmrCalcClasses:
       notifyFunc(self.updateRunDataAfter, NC+clazz, '')
       
    
  def destroy(self):

    self.administerNotifiers(self.parent.unregisterNotify)
    Frame.destroy(self)

  def executePales(self):
    # NBNB TBD
    pass
    
  def viewPalesData(self):
    # NBNB TBD
    valueObj = self.outputMatrix
    
    if isinstance(valueObj, NmrCalc.ExternalData):
      dataStore = valueObj.dataStore
      if dataStore is not None:
        wb = WebBrowser(self, name='PalesData')
        wb.open(dataStore.fullPath)
      
    elif isinstance(valueObj, NmrCalc.MeasurementListData):
      mm = valueObj.measurementList
      if mm is not None:
        popup = self.parent.editMeasurementLists()
        popup.tabbedFrame.select(1)
        popup.setMeasurementList(mm)
       
    elif isinstance(valueObj, NmrCalc.StructureEnsembleData):
      models = valueObj.models
      if len(models) == 1:
        popup = self.parent.editStructures()
        popup.tabbedFrame.select(3)
        popup.changeModel(models[0])
        popup.changeTab(3)
        
    elif isinstance(valueObj, NmrCalc.ConstraintStoreData):
      constraintLists = valueObj.constraintLists
      if len(constraintLists) == 1:
        popup = self.parent.browseConstraints()
        popup.tabbedFrame.select(2)
        popup.changeRestraintList(constraintLists[0])
  
