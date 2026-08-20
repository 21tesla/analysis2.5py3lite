import os

from memops.universal.ElementTree import ElementTree

from memops.api.Implementation import AppDataInt

from ccp.general.Io import getDataSourceFileName

def getShapeFile(spectrum):

  from memops.c import ShapeFile

  # TBD: shapeDim is way of getting 1D subset of reconstructed ND data set
  root = spectrum.root
  if hasattr(root, 'application'):
    application = root.application
    appData = spectrum.findFirstApplicationData(application=application.name, keyword='shapeDim')
    if appData:
      shapeDim = appData.value
    else:
      shapeDim = None

  fileName = getDataSourceFileName(spectrum)
  if not fileName or not os.path.exists(fileName):
    msg = 'Warning: spectrum (%s, %s): data file %s not accessible'
    print msg % (spectrum.experiment.name, spectrum.name, fileName)
    return None


  try:
    elementTree = ElementTree.parse(fileName)
    topNode = elementTree.getroot()
    ncomp = int(topNode.get('ncomponents'))
    nshapes = int(topNode.get('nshapes'))

    sizeList = []
    componentList = []
    for elem in topNode:
      if elem.tag == 'Axis':
        a = int(elem.get('a'))
        size = int(elem.get('size'))
        sizeList.append((a, size))
      elif elem.tag == 'Component':
        c = int(elem.get('c'))
        ampl = float(elem.get('ampl'))
        shapeData = []
        componentList.append((c, ampl, shapeData))
        for ee in elem:
          if ee.tag == 'Shape':
            a = int(ee.get('a'))
            size = int(ee.get('size', 0))
            offset = int(ee.get('offset', 0))
            text = ee.text.strip().split()
            values = [float(t) for t in text]
            shapeData.append((a, size, offset, values))

    if ncomp != len(componentList):
      print ('Warning: in %s, ncomp = %d, but found %d components' 
             % (fileName, ncomp, componentList))

    sizeList.sort()
    componentList.sort()

    ncomp = len(componentList)

    cmin = componentList[0][0]
    cmax = componentList[-1][0]
    if cmax - cmin + 1 != ncomp:
      print 'Warning: in %s, components are not numbered consecutively' % fileName

    aDict = {}
    points = []
    for n, (a, size) in enumerate(sizeList):
      aDict[a] = n
      points.append(size)

    if shapeDim is None:
      pts = points
    else:
      pts = points[shapeDim-1:shapeDim]

    shapeFile = ShapeFile.ShapeFile(ncomp, pts)
    for comp, (c, ampl, shapeData) in enumerate(componentList):
      shapeFile.setComponentAmplitude(comp, ampl)

      shapeData.sort()
      for (a, size, offset, values) in shapeData:
        aa = aDict[a]
        if shapeDim is None or (aa == shapeDim-1):
          if size == 0:
            size = points[aa]
          if len(values) !=size:
            raise IOError('in %s: comp = %d, a = %d, size = %d but found %d values' % (fileName, c, a, size, len(values)))
          if shapeDim is None:
            adim = 0
          else:
            adim = aa
          shapeFile.setShapeData(comp, aa, offset, values)

  except Exception, e:
    print 'Warning: in %s had exception: %s' % (fileName, e)
    shapeFile = None

  except ShapeFile.error, e:
    print 'Warning: in %s had ShapeFile error: %s' % (fileName, e)
    shapeFile = None

  return shapeFile

# get 1D shapeFile based on already created 1D spectrum and data
# number of components = len(valuesList)
# data for component c is valuesList[c]
def get1dShapeFile(spectrum, valuesList):

  from memops.c import ShapeFile

  try:
    ncomp = len(valuesList)
    points = (len(valuesList[0]),)
    shapeFile = ShapeFile.ShapeFile(ncomp, points)
    for comp in range(ncomp):
      # default is 1.0 so don't need to do below
      #shapeFile.setComponentAmplitude(comp, 1.0)
      dim = offset = 0
      shapeFile.setShapeData(comp, dim, offset, valuesList[comp])

  except Exception, e:
    print 'Warning: in %s:%s had exception: %s' % (spectrum.experiment.name, spectrum.name, e)
    shapeFile = None

  except ShapeFile.error, e:
    print 'Warning: in %s:%s had ShapeFile error: %s' % (spectrum.experiment.name, spectrum.name, e)
    shapeFile = None

  return shapeFile

# create 1D shape spectrum based on ND shape spectrum
def create1dShape(freqDataDim, applicationName):

  # can some of below be done with copySubTree??

  dim = freqDataDim.dim
  spectrum = freqDataDim.parent

  dataStore = spectrum.dataStore
  dataLocationStore = dataStore.parent
  newDataStore = dataLocationStore.newShapeMatrix(dataUrl=dataStore.dataUrl,
      path=dataStore.path, numPoints=(freqDataDim.numPoints,),
      isBigEndian=dataStore.isBigEndian, numberType=dataStore.numberType,
      isComplex=(freqDataDim.isComplex,), headerSize=dataStore.headerSize,
      numRecords=dataStore.numRecords)

  name = '%s_%d' % (spectrum.name, dim)
  experiment = spectrum.parent
  newSpectrum = experiment.newDataSource(name=name, numDim=1,
      dataStore=newDataStore, dataType=spectrum.dataType)

  appData = AppDataInt(application=applicationName, keyword='shapeDim', value=dim)
  newSpectrum.addApplicationData(appData)

  newFreqDataDim = newSpectrum.newFreqDataDim(dim=dim, numPoints=freqDataDim.numPoints,
      isComplex=freqDataDim.isComplex, numPointsOrig=freqDataDim.numPointsOrig,
      valuePerPoint=freqDataDim.valuePerPoint, expDim=freqDataDim.expDim)

  for dataDimRef in freqDataDim.dataDimRefs:
    newFreqDataDim.newDataDimRef(refPoint=dataDimRef.refPoint,
        refValue=dataDimRef.refValue, expDimRef=dataDimRef.expDimRef)

  return newSpectrum

if __name__ == '__main__':

  import sys

  if len(sys.argv) != 2:
    print 'need to specify shapeFile'
    sys.exit()

  fileName = sys.argv[1]

  shapeFile = getShapeFile(fileName)

