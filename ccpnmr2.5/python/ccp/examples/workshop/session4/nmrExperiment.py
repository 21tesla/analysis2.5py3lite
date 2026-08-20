def createExperiment(nmrProject, name, numDim, expDimRefDict):

  experiment = nmrProject.newExperiment(name=name, numDim=numDim)

  for expDim in experiment.expDims:
    dim = expDim.dim
    sf = expDimRefDict['sf'][dim-1]
    isotopeCodes = (expDimRefDict['isotopeCode'][dim-1],)
    expDim.newExpDimRef(sf=sf, isotopeCodes=isotopeCodes, unit='ppm')

  return experiment

def createDataSource(experiment, name, numDim, dataDimDict):

  dataSource = experiment.newDataSource(name=name, numDim=numDim, dataType='processed')

  for expDim in experiment.sortedExpDims():
    dim = expDim.dim
    isComplex = dataDimDict['isComplex'][dim-1]
    numPoints = dataDimDict['numPoints'][dim-1]
    numPointsOrig = dataDimDict['numPointsOrig'][dim-1]
    valuePerPoint = dataDimDict['valuePerPoint'][dim-1]
    freqDataDim = dataSource.newFreqDataDim(dim=dim, isComplex=isComplex,
                         numPoints=numPoints, numPointsOrig=numPointsOrig,
                         valuePerPoint=valuePerPoint, expDim=expDim)

    expDimRef = expDim.findFirstExpDimRef()
    refPoint = dataDimDict['refPoint'][dim-1]
    refValue = dataDimDict['refValue'][dim-1]
    freqDataDim.newDataDimRef(expDimRef=expDimRef, refPoint=refPoint, refValue=refValue)

  return dataSource

def setRefExperiment(experiment, expPrototypeSynonym, refExperimentName = None):

  root = experiment.root
  nmrExpPrototype = root.findFirstNmrExpPrototype(synonym=expPrototypeSynonym)
  if not nmrExpPrototype:
    print 'WARNING: no nmrExpPrototype with synonym %s' % expPrototypeSynonym
    return

  if refExperimentName:
    refExperiment = nmrExpPrototype.findFirstRefExperiment(name=refExperimentName)
    if not refExperiment:
      print 'WARNING: no refExperiment with name %s for nmrExpPrototype with synonym %s' % (refExperimentName, expPrototypeSynonym)
      return
  elif len(nmrExpPrototype.refExperiments) > 1:
    print 'WARNING: more than one refExperiment for nmrExpPrototype with synonym %s' % expPrototypeSynonym
    return
  else:
    refExperiment = nmrExpPrototype.findFirstRefExperiment()

  experiment.refExperiment = refExperiment

def main(root):

  # find or create NmrProject
  nmrProject = root.findFirstNmrProject()
  if not nmrProject:
    nmrProject = root.newNmrProject(name='testNmr')

  # create experiment, including expDimRefs
  dd = {}
  dd['sf'] = (800, 81)
  dd['isotopeCode'] = ('1H', '15N')
  experiment = createExperiment(nmrProject, name='testNmr', numDim=2, expDimRefDict=dd)

  # set refExperiment
  setRefExperiment(experiment, expPrototypeSynonym='15N HSQC/HMQC')

  # create dataSource, including freqDataDims and dataDimRefs
  dd = {}
  dd['isComplex'] = (False, False)
  dd['numPoints'] = (2048, 1024)
  dd['numPointsOrig'] = (2048, 1024)
  dd['valuePerPoint'] = (3.93, 0.99)
  dd['refPoint'] = (1024, 512)
  dd['refValue'] = (4.72, 117.4)
  dataSource = createDataSource(experiment, name='testSpec', numDim=2, dataDimDict=dd)

  # create peakList
  peakList = dataSource.newPeakList()

  # create dataLocationStore
  dataLocationStore = root.findFirstDataLocationStore()
  if not dataLocationStore:
    dataLocationStore = root.newDataLocationStore(name='testStore')

  # create dataStore
  from memops.api.Implementation import Url
  dd = {}
  dd['path'] = 'myHsqcData.spc'
  dd['dataUrl'] = dataLocationStore.newDataUrl(url=Url(path='/usr/myaccount/mydata'))
  dd['blockSizes'] = (256, 128)
  dd['headerSize'] = 0
  dd['isBigEndian'] = True
  dd['isComplex'] = (False, False)
  dd['nByte'] = 4
  dd['numPoints'] = [dataDim.numPoints for dataDim in dataSource.sortedDataDims()]
  dd['numberType'] = 'float'
  dataStore = dataLocationStore.newBlockedBinaryMatrix(**dd)

  dataSource.dataStore = dataStore

  root.saveModified()

if __name__ == '__main__':

  import sys

  if len(sys.argv) == 1:
    from memops.api.Implementation import MemopsRoot
    root = MemopsRoot(name='testNmr')
  elif len(sys.argv) == 2:
    from memops.general.Io import loadProject
    repositoryPath = sys.argv[1]
    root = loadProject(path=repositoryPath)

  main(root)
