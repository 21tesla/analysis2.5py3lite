from memops.universal import BlockData
from ccp.format.spectra.params.AzaraParams import AzaraParams
from ccp.format.spectra.params.NmrPipeParams import NmrPipeParams, getDataFileName, getSampledValue

pipeZTemplate = '%03d'
pipeZTemplateAlt = '%04d'
pipeAZTemplate = '%02d%03d'

def convertXyzToAzara(pipeDataFileTemplate, azaraDataFile, azaraParFile = None, totalBlockSize = 4096):

  n = pipeDataFileTemplate.find(pipeAZTemplate)
  if n >= 0:
    haveAZ = True
    baseFile = getDataFileName(pipeDataFileTemplate, z=0, a=0)
  else:
    haveAZ = False
    n = pipeDataFileTemplate.find(pipeZTemplate)
    if n < 0:
      n = pipeDataFileTemplate.find(pipeZTemplateAlt)
    if n >= 0:
      haveZ = True
      baseFile = getDataFileName(pipeDataFileTemplate, z=0)
    else:
      haveZ = False
      baseFile = pipeDataFileTemplate

  baseParams = NmrPipeParams(baseFile)

  if haveAZ:
    pipeParams = []
    for a in range(baseParams.npts[3]):
      for z in range(baseParams.npts[2]):
        dataFile = getDataFileName(pipeDataFileTemplate, z=z, a=a)
        params = NmrPipeParams(dataFile)
        if params != baseParams:
          raise Exception('inconsistent header for files %s and %s' % (baseFile, dataFile))
        pipeParams.append(params)
  elif haveZ:
    pipeParams = []
    for z in range(baseParams.npts[2]):
      dataFile = getDataFileName(pipeDataFileTemplate, z=z)
      params = NmrPipeParams(dataFile)
      if params != baseParams:
        raise Exception('inconsistent header for files %s and %s' % (baseFile, dataFile))
      pipeParams.append(params)
  else:
    pipeParams = baseParams

  convertPipeToAzara(pipeParams, azaraDataFile, azaraParFile, totalBlockSize)

def convertPipeToAzara(pipeParams, azaraDataFile, azaraParFile = None, totalBlockSize = 4096):

  # first write Azara par file
  azaraParams = writeAzaraParFile(pipeParams, azaraDataFile, azaraParFile, totalBlockSize)

  # now write out Azara data file
  writeAzaraDataFile(pipeParams, azaraParams)

def writeAzaraParFile(pipeParams, azaraDataFile, azaraParFile = None, totalBlockSize = 4096):

  if not azaraParFile:
    azaraParFile = azaraDataFile + '.par'

  pointValues = []
  if type(pipeParams) == type([]):
    baseParams = pipeParams[0]
    dim = baseParams.pseudoDataDim()
    if dim is not None:  
      for params in pipeParams:
        value = getSampledValue(params.dataFile)
        pointValues.append(value)
  else:
    baseParams = pipeParams

  azaraParams = AzaraParams(parFile=azaraParFile, externalParams=baseParams)
  azaraParams.dataFile = azaraDataFile
  azaraParams.head = 0
  azaraParams.block = BlockData.determineBlockSizes(baseParams.npts, totalBlockSize)
  if pointValues:
    azaraParams.setSampledDim(dim, pointValues)
    azaraParams.nuc[dim] = '1H' # arbitrary

  azaraParams.writeParFile()

  return azaraParams

def writeAzaraDataFile(pipeParams, azaraParams):

  ndim = azaraParams.ndim
  dataFile = azaraParams.dataFile
  block = azaraParams.block
  npts = azaraParams.npts

  (blkSize, cumulBlockSize) = BlockData.cumulativeArray(block)

  nblocks = [ 1 + (npts[i] - 1)/block[i] for i in range(ndim) ]
  (nblks, cumulBlocks) = BlockData.cumulativeArray(nblocks)

  (n, cumulPoints) = BlockData.cumulativeArray(npts)
  (nblkRows, cumulBlockRow) = BlockData.cumulativeArray(block[1:])

  if type(pipeParams) == type([]): # multiple files
    multiFile = True
    baseParams = pipeParams[0]
    (n, cumulOrthogPoints) = BlockData.cumulativeArray(npts[2:])
    cumulPoints = cumulPoints[:2]
    orthogPlane = None
  else:
    multiFile = False
    baseParams = pipeParams

  fpw = open(dataFile, 'wb')

  zeroWord = '\0\0\0\0'
  fpr = None
  rowLen = block[0]
  for blk in range(nblks):
    y = BlockData.arrayOfIndex(blk, cumulBlocks)
    for row in range(nblkRows):
      z = (0,) + BlockData.arrayOfIndex(row, cumulBlockRow)
      basePoint = [ block[i]*y[i] + z[i] for i in range(ndim) ]

      if multiFile:
        plane = BlockData.indexOfArray(basePoint[2:], cumulOrthogPoints)
        offsetPoint = basePoint[:2]
      else:
        offsetPoint = basePoint

      if not fpr or (multiFile and plane != orthogPlane):
        if multiFile:
          if fpr:
            fpr.close()
          try:
            dataFile = pipeParams[plane].dataFile
          except:
            dataFile = None # will write zeros below
          orthogPlane = plane
        else:
          dataFile = pipeParams.dataFile
        if dataFile:
          fpr = open(dataFile, 'rb')

      for i in range(1, ndim):
        if (basePoint[i] >= npts[i]):
          for j in range(rowLen):
            fpw.write(zeroWord)
          break
      else:
        n = min(rowLen, npts[0]-basePoint[0])
        if n > 0:
          offset = baseParams.head + 4*BlockData.indexOfArray(offsetPoint, cumulPoints)
          fpr.seek(offset, 0)
          d = fpr.read(4*n)
          fpw.write(d)
        for j in range(rowLen-n):
          fpw.write(zeroWord)

  fpw.close()
  if fpr:
    fpr.close()

if (__name__ == '__main__'):

  import math
  import sys
  import os

  def usage(msg):

    print 'Error: %s' % msg
    print 'Arguments: <pipeDataFileTemplate> <azaraDataFile> [ <blockSize> <azaraParFile> ]'
    print '     For example: myPipeData%03d.ft myAzaraData.spc'
    print '              or: myPipeData%03d.ft myAzaraData.spc 16384'
    print '              or: myPipeData%03d.ft myAzaraData.spc 4096 myAzaraParFile.par'
    os._exit(1)

  nargs = len(sys.argv)
  if nargs < 3 or nargs > 5:
    usage('need 2 to 4 arguments')

  pipeDataFileTemplate = sys.argv[1]
  azaraDataFile = sys.argv[2]

  if nargs >= 4:
    try:
      blockSize = int(sys.argv[3])
    except:
      usage('blockSize must be integer')

    if blockSize < 1:
      usage('blockSize must be > 0')

    b = pow(2, int(math.log(blockSize)/math.log(2)))
    if b != blockSize:
      usage('blockSize must be power of 2')
  else:
    blockSize = 4096

  if nargs >= 5:
    azaraParFile = sys.argv[4]
  else:
    azaraParFile = None

  convertXyzToAzara(pipeDataFileTemplate, azaraDataFile, azaraParFile, blockSize)
