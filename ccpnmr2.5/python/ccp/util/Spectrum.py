import os

from memops.universal.Io import normalisePath

from memops.universal.BlockData import determineBlockSizes

from memops.api.Implementation import Url

from memops.general.Implementation import ApiError

def createExperiment(nmrProject, name, numDim, sf, isotopeCodes):

  experiment = nmrProject.newExperiment(name=name, numDim=numDim)

  expDims = experiment.sortedExpDims()
  for n in range(numDim):
    ic = isotopeCodes[n]
    if ic:
      if type(ic) == type(''):
        ic = (ic,)
      expDims[n].newExpDimRef(sf=sf[n], unit='ppm', isotopeCodes=ic)

  return experiment

def createSpectrum(experiment, name, numPoints, sw, refppm, refpt,
                   dataStore = None, writeable = False):

  numDim = len(numPoints)

  if numDim > experiment.numDim:
    raise ApiError('numDim = %d > %d = experiment.numDim' % (numDim, experiment.numDim))

  spectrum = experiment.newDataSource(name=name, dataStore=dataStore,
                                      numDim=numDim, dataType='processed')
 
  spectrum.writeable = writeable

  expDims = experiment.sortedExpDims()
  for n in range(numDim):
    expDim = expDims[n]
    freqDataDim = spectrum.newFreqDataDim(dim=n+1, numPoints=numPoints[n],
                             isComplex=False, numPointsOrig=numPoints[n],
                             valuePerPoint=sw[n]/float(numPoints[n]), expDim=expDim)
    expDimRef = expDimRef=expDim.findFirstExpDimRef()
    if expDimRef:
      freqDataDim.newDataDimRef(refPoint=refpt[n], refValue=refppm[n], expDimRef=expDimRef)

  return spectrum

def createBlockedMatrix(dataUrl, path, numPoints, blockSizes = None,
                        isBigEndian = True, numberType = 'float', isComplex = None):

  path = normalisePath(path)

  if os.path.isabs(path):
    urlpath = normalisePath(dataUrl.url.path, makeAbsolute=True)
    if not path.startswith(urlpath):
      raise ApiError('path = %s, does not start with dataUrl path = %s' % (path, urlpath))
    if path == urlpath:
      raise ApiError('path = %s, same as dataUrl path but should be longer' % path)

    # TBD: below is a bit dangerous but should work (+1 is to remove '/')
    path = path[len(urlpath)+1:]

  if not blockSizes:
    blockSizes = determineBlockSizes(numPoints)

  if not isComplex:
    isComplex = len(numPoints) * [False]

  dataLocationStore = dataUrl.dataLocationStore

  matrix = dataLocationStore.newBlockedBinaryMatrix(dataUrl=dataUrl, path=path,
                 numPoints=numPoints, blockSizes=blockSizes, isBigEndian=isBigEndian,
                 numberType=numberType, isComplex=isComplex)

  return matrix
  
if __name__ == '__main__':

  from memops.api.Implementation import MemopsRoot

  r = MemopsRoot()
  n = r.newNmrProject(name='testNmrProject')
  e = createExperiment(n, name='testExpt', numDim=2, sf=(800, 150), isotopeCodes=('H', 'C'))
  numPoints = (512, 256)
  d = r.newDataLocationStore(name='testDLS')
  u = d.newDataUrl(name='testDataUrl', url=Url(path='/my/test/path'))
  m = createBlockedMatrix(u, 'test.spc', numPoints)
  s = createSpectrum(e, name='testSpectrum', numPoints=numPoints, sw=(8000, 2000),
                     refppm=(5, 50), refpt = (256, 128), dataStore=m)

