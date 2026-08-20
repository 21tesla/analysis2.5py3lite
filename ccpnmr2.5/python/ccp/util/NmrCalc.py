"""NmrCalc related handler functions.
Mostly contain only 2-3 lines of obvious code (Rasmus: ah well).
Used (as per Oct 2012) in 
paris/aria, nijmegen/cing/CingFrame.py, 
grenoble/BlackledgeModule/BlackledgeModuleFrame.py

November 2012: Removed references to groupId
"""


DATA_MISSING = '*DATA MISSING*'
PARAM_ATTR_DICT = {type(1.0):'floatValue',
                   type(1):'intValue',
                   type(''):'textValue',
                   type(True):'booleanValue'}

def getObjBooleanParameter(dataObj, code, default=True):

  run = dataObj.run
  runParameter = dataObj.findFirstRunParameter(code=code)
  
  if not runParameter: # Not already set.
    runParameter = dataObj.run.newRunParameter(code=code)
    result = default
    runParameter.booleanValue = result
  
  elif runParameter.booleanValue is None:
    result = default
    runParameter.booleanValue = result
  
  else:
    result = runParameter.booleanValue
    
  if runParameter not in dataObj.runParameters:
    dataObj.addRunParameter(runParameter)
  
  return result

def toggleObjBooleanParameter(dataObj, code, default=True):

  run = dataObj.run
  runParameter = dataObj.findFirstRunParameter(code=code)
  
  if not runParameter: # Not already set.
    runParameter = dataObj.run.newRunParameter(code=code)
    result = default
  
  elif runParameter.booleanValue is None:
    result = default
  
  else:
    result = not runParameter.booleanValue
    
  if runParameter not in dataObj.runParameters:
    dataObj.addRunParameter(runParameter)

  runParameter.booleanValue =  result
  
  return result

def setObjRunParameter(dataObj, code, value):

  run = dataObj.run
  runParameter = dataObj.findFirstRunParameter(code=code)
  
  if not runParameter: # Not already set.
    runParameter = dataObj.run.newRunParameter(code=code)
    
  if runParameter not in dataObj.runParameters:
    dataObj.addRunParameter(runParameter)
  
  setattr(runParameter, PARAM_ATTR_DICT[type(value)], value)
  
  
def getObjRunParameter(dataObj, code):

  runParameter = dataObj.findFirstRunParameter(code=code)

  return runParameter

def deleteRunParameter(run, code, data=False):
  
  if data or data is None:
    runParameter = run.findFirstRunParameter(code=code, data=data)
  else:
    runParameter = run.findFirstRunParameter(code=code)
    
  if runParameter:
    return runParameter.delete()
      
def setRunParameter(run, code, value, data=False):

  if data or data is None:
    runParameter = (run.findFirstRunParameter(code=code, data=data) or
                    run.newRunParameter(code=code, data=data))
  
  else:
    runParameter = (run.findFirstRunParameter(code=code) or
                    run.newRunParameter(code=code))
  
  setattr(runParameter, PARAM_ATTR_DICT[type(value)], value)
 
def getRunParameter(run, code, data=False):

  if data or data is None:
    runParameter = run.findFirstRunParameter(code=code, data=data)
  else:
    runParameter = run.findFirstRunParameter(code=code)
  
  return runParameter

# Run get parameters by type
 
def getRunTextParameter(run, code, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data)
  
  if runParameter:
    return runParameter.textValue
 
def getRunBooleanParameter(run, code, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data)
  
  if runParameter:
    return runParameter.booleanValue

def getRunIntParameter(run, code, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data)
  
  if runParameter:
    return runParameter.intValue

def getRunIntParameterList(run, code, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data)
  
  if runParameter:
    data = runParameter.textValue
    if data:
      values = [int(x) for x in data.split(',')]
    else:
      values = []
  
    return values or None

def getRunTextParameterDict(run, code, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data)
  
  if runParameter:
    data = runParameter.textValue or ''
    paramDict = {}
    
    if data:
      pairs = [x.split(':') for x in data.split(',')]
    
      for key, value in pairs:
        paramDict[key] = value
  
    return paramDict

def getRunTextParameterList(run, code, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data)
  
  if runParameter:
    data = runParameter.textValue
    if data:
      values = [x for x in data.split(',')]
    else:
      values = []
  
    return values or None

def getRunFloatParameter(run, code, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data)
  
  if runParameter:
    return runParameter.floatValue

# Run set parameters by type
 
def setRunTextParameter(run, code, value, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data) or \
                 run.newRunParameter(code=code, data=data)

  runParameter.textValue = value or None
 
def setRunBooleanParameter(run, code, value, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data) or \
                 run.newRunParameter(code=code, data=data)
  
  runParameter.booleanValue = value

def setRunIntParameter(run, code, value, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data) or \
                 run.newRunParameter(code=code, data=data)
  
  runParameter.intValue = value

def setRunIntParameterList(run, code, listData, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data) or \
                 run.newRunParameter(code=code, data=data)
  
  values = [str(x) for x in listData or []]
  value = ','.join(values)
  
  runParameter.textValue = value or None

def setRunTextParameterDict(run, code, dictData, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data) or \
                 run.newRunParameter(code=code, data=data)
  
  values = ['%s:%s' % (str(key), str(value)) for key, value in dictData.items()]
  value = ','.join(values)

  runParameter.textValue = value or None
  
def setRunTextParameterList(run, code, listData, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data) or \
                 run.newRunParameter(code=code, data=data)
  
  values = [str(x) for x in listData or []]
  value = ','.join(values)
  
  runParameter.textValue = value or None


def setRunFloatParameter(run, code, value, data=None):
  """ NB default is *no* data, not *any* data"""

  runParameter = run.findFirstRunParameter(code=code, data=data) or \
                 run.newRunParameter(code=code, data=data)
                 
  runParameter.floatValue = value

# Data objects get parameters by type

def getDataTextParameter(data, code):

  runParameter = data.findFirstRunParameter(code=code)
  
  if runParameter:
    return runParameter.textValue
 
def getDataBooleanParameter(data, code):

  runParameter = data.findFirstRunParameter(code=code)
  
  if runParameter:
    return runParameter.booleanValue

def getDataIntParameter(data, code):

  runParameter = data.findFirstRunParameter(code=code)
  
  if runParameter:
    return runParameter.intValue

def getDataIntParameterList(data, code):

  runParameter = data.findFirstRunParameter(code=code)
  
  if runParameter:
    data = runParameter.textValue
    
    if data:
      values = [int(x) for x in data.split(',')]
    else:
      values = []
  
    return values or None

def getDataTextParameterDict(data, code):

  runParameter = data.findFirstRunParameter(code=code)
  
  if runParameter:
    data = runParameter.textValue or ''
    paramDict = {}
    
    if data:
      pairs = [x.split(':') for x in data.split(',')]
    
      for key, value in pairs:
        paramDict[key] = value
  
    return paramDict

def getDataTextParameterList(data, code):

  runParameter = data.findFirstRunParameter(code=code)
  
  if runParameter:
    data = runParameter.textValue
    
    if data:
      values = [x for x in data.split(',')]
    else:
      values = []
  
    return values or None

def getDataFloatParameter(data, code):

  runParameter = data.findFirstRunParameter(code=code)
  
  if runParameter:
    return runParameter.floatValue


# Data objects set parameters by type

def fetchDataRunParameter(datum, code):

  runParameter = datum.findFirstRunParameter(code=code)
  
  if not runParameter:
    runParameter = datum.run.newRunParameter(code=code, data=datum)
 
  return runParameter

def setDataTextParameter(datum, code, value):

  runParameter = fetchDataRunParameter(datum, code)
  runParameter.textValue = value or None
 
def setDataBooleanParameter(datum, code, value):

  runParameter = fetchDataRunParameter(datum, code)
  runParameter.booleanValue = value

def setDataIntParameter(datum, code, value):

  runParameter = fetchDataRunParameter(datum, code)
  runParameter.intValue = value

def setDataIntParameterList(datum, code, listData):

  runParameter = fetchDataRunParameter(datum, code)
  
  values = [str(x) for x in listData or []]
  value = ','.join(values)
  
  runParameter.textValue = value or None

def setDataTextParameterDict(datum, code, dictData):

  runParameter = fetchDataRunParameter(datum, code)
  
  values = ['%s:%s' % (str(key), str(value)) for key, value in dictData.items()]
  value = ','.join(values)

  runParameter.textValue = value or None
  
def setDataTextParameterList(datum, code, listData):

  runParameter = fetchDataRunParameter(datum, code)
  
  values = [str(x) for x in listData or []]
  value = ','.join(values)

  runParameter.textValue = value or None

def setDataFloatParameter(datum, code, value):

  runParameter = fetchDataRunParameter(datum, code)
  runParameter.floatValue = value

# Misc

def getRangeString(numbers):

  if len(numbers) == 1:
    return '%d' % (numbers[0],)

  numbers.sort()
  ranges = [[numbers[0],None]]
  
  for i, numberA in enumerate(numbers[:-1]):
    numberB = numbers[i+1]
    
    if numberB != (numberA+1):
      ranges[-1][1] = numberA
      ranges.append([numberB,None])
  
  ranges[-1][1] = numbers[-1]
  
  ranges = ['%d-%d' % tuple(r) for r in ranges]
  
  return ','.join(ranges)   

def getDataObjText(datum, missingText=DATA_MISSING):
  
  # TBD missing data info
    
  funcs = {'ConstraintStoreData':getConstraintStoreDataInfo,
           'DerivedListData':getDerivedListDataInfo,
           'FloatMatrixData':getFloatMatrixDataInfo,
           'ExternalData':getExternalDataInfo,
           'MeasurementListData':getMeasurementListDataInfo,
           'MolResidueData':getMolResidueDataInfo,
           'MolSystemData':getMolSystemDataInfo,
           'PeakListData':getPeakListDataInfo,
           'SpectrumData':getSpectrumDataInfo,
           'SpinSystemData':getSpinSystemDataInfo,
           'StructureEnsembleData':getStructureEnsembleDataInfo,
           'ViolationListData':getViolationListDataInfo}
  
  func = funcs[datum.className]

  return func(datum)

def getFloatMatrixDataInfo(datum):
  return "FloatMatrix [%s]" % datum.size

def getConstraintStoreDataInfo(datum):

  serials = list(datum.constraintListSerials)
  constraintLists = [cl for cl in datum.constraintLists if cl]
  
  if constraintLists and not serials:
    serials = [cl.serial for cl in constraintLists]
  
  if serials:
    if len(serials) == 1:
     
     if constraintLists:
       cList =  '%d - %s' % (serials[0], constraintLists[0].className[:-14])
     else:
       cList = str(serials[0])
     
     info = 'Restraint List %d:%s' % (datum.constraintStoreSerial, cList)
   
    else:
      serials.sort()
      cLists = ','.join([str(s) for s in serials])
      info = 'Restraint Lists %d:%s' % (datum.constraintStoreSerial, cLists)
 
    
  else:
    info = 'Restraint Lists %d:* all *' % (datum.constraintStoreSerial)


  return info

def getDerivedListDataInfo(datum):

  derivedDataList = datum.derivedDataList

  if derivedDataList:
    dType = derivedDataList.className[:-4]
  else:
    dType = 'Derived Data'

  info = '%s List %d' % (dType, datum.derivedDataListSerial)

  return info
  
def getExternalDataInfo(datum):

  dataStore = datum.dataStore
  
  if dataStore:
    fileName = dataStore.fullPath
  else:
    fileName = '* No file *'

  info = 'External Data; %s' % fileName

  return info
  
def getMeasurementListDataInfo(datum):

  mList = datum.measurementList

  if mList:
    dType = mList.className[:-4]
  else:
    dType = 'Measurement'

  info = '%s List %d' % (dType, datum.measurementListSerial)

  return info
  
def getMolResidueDataInfo(datum):
  
  residueSeqIds = datum.residueSeqIds
  
  if residueSeqIds:
  
    ranges = []
    start = residueSeqIds[0]
    end = start
    for i in residueSeqIds[1:]:
      if i == end+1:
        end = i
      else:
        ranges.append((start,end))
        start = end = i
 
    ranges.append((start,end))
 
    texts = []
    for start, end in ranges:
      if start == end:
        text = '%d' % (start)
      else:
        text = '%d-%d' % (start,end)
 
      texts.append(text)
  
    residues =  ','.join(texts)
    
  else:
    residues = '* all *'

  info = 'Chain %s:%s; Residues %s' % (datum.molSystemCode, datum.chainCode, residues)

  return info
  
def getMolSystemDataInfo(datum):

  chainCodes = list(datum.chainCodes)
  
  if chainCodes:
    chainCodes.sort()
    chains = ','.join(chainCodes)
    
  else:
    molSystem = datum.molSystem
    
    if molSystem:
      chains = ','.join([c.code for c in molSystem.sortedChains()])
    else:
      chains = '* all *'

  info = 'Chains %s:%s' % (datum.molSystemCode, chains)

  return info
  
def getPeakListDataInfo(datum):

  spectrum = datum.dataSource
  
  if spectrum:
    eName = spectrum.experiment.name
    sName = spectrum.name
  
  else:
    eName = str(datum.experimentSerial)
    sName = str(datum.dataSourceSerial)

  info = 'PeakList %s:%s:%d' % (eName, sName, datum.peakListSerial)

  return info
  
def getSpectrumDataInfo(datum):

  spectrum = datum.dataSource
  
  if spectrum:
    eName = spectrum.experiment.name
    sName = spectrum.name
  
  else:
    eName = str(datum.experimentSerial)
    sName = str(datum.dataSourceSerial)

  info = 'Spectrum %s:%s' % (eName, sName)  

  return info
  
def getSpinSystemDataInfo(datum):

  spinSystem = datum.resoanceGroup
  
  if spinSystem:
  
    residue = spinSystem.residue
    if residue:
      assignInfo = '%d%s' % (residue.seqCode, residue.ccpCode)
      
    elif spinSystem.ccpCode:
      assignInfo = spinSystem.ccpCode
      
    else:
      assignInfo = ''
  
  else:
    assignInfo = ''

  info = 'Spin System %d%s' % (datum.resonanceGroupSerial, assignInfo)  

  return info
  
def getStructureEnsembleDataInfo(datum):

  serials = list(datum.modelSerials)
  models = datum.models
    
  if models and not serials:
    serials = [m.serial for m in models]

  if serials:
    serials.sort()
    
    ranges = []
    start = serials[0]
    end = start
    for i in serials[1:]:
      if i == end+1:
        end = i
      else:
        ranges.append((start,end))
        start = end = i
 
    ranges.append((start,end))
 
    texts = []
    for start, end in ranges:
      if start == end:
        text = '%d' % (start)
      else:
        text = '%d-%d' % (start,end)
 
      texts.append(text)

    modelInfo =  ','.join(texts)
  
  else:
    modelInfo = '* all *'

  info = 'Structure Ensemble %s:%d; models %s' % (datum.molSystemCode, datum.ensembleId, modelInfo)

  return info
  
def getViolationListDataInfo(datum):

  info = 'Violation List %d:%d' % (datum.constraintStoreSerial, datum.violationListSerial)

  return info
  















