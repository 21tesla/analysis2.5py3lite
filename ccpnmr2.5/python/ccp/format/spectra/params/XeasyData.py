import array
import math

from ccp.format.spectra.params.AzaraParams import AzaraParams
from ccp.format.spectra.params.XeasyParams import XeasyParams

def _powersign(exponent):

  if exponent <= 47:
    ll = exponent - 1
    ss = 1
  else:
    ll = 95 - exponent
    ss = -1
  
  return (ll, ss)

# exponent is signed int
def xeasy8ToReal(exponent):

  (ll, ss) = _powersign(exponent)

  value = ss*math.pow(2.0, 0.5*ll)

  return value

def xeasy16ToReal(exponent, mantissa):

  (ll, ss) = _powersign(exponent)

  value = ss * ((615 + mantissa) * math.pow(2.0, 0.5*ll)) / 721

  return value

def writeAzaraParFile(xeasyParams, azaraDataFile, azaraParFile):

  azaraParams = AzaraParams(parFile=azaraParFile, externalParams=xeasyParams)
  azaraParams.dataFile = azaraDataFile
  azaraParams.head = 0
  azaraParams.nbytes = 4
  azaraParams.swap = False

  azaraParams.writeParFile()

  return azaraParams

def writeAzaraDataFile(xeasyParams, azaraParams):

  xeasyDataFile = xeasyParams.dataFile
  nbytes = xeasyParams.nbytes
  ndim = azaraParams.ndim
  azaraDataFile = azaraParams.dataFile
  block = azaraParams.block
  npts = azaraParams.npts

  nblocks = 1
  blockSize = 1
  for i in range(ndim):
    nblocks *= 1 + (npts[i] - 1)/block[i]
    blockSize *= block[i]

  fpr = open(xeasyDataFile, 'rb')
  fpw = open(azaraDataFile, 'wb')

  for i in range(nblocks):
    data = fpr.read(nbytes*blockSize)
    x = array.array('b')
    x.fromstring(data)
    z = array.array('f')
    for j in range(blockSize):
      if nbytes == 1:
        exponent = x[j]
        value = xeasy8ToReal(exponent)
      else:
        exponent = x[2*j+1]
        mantissa = x[2*j]
        if mantissa < 0:
          mantissa += 256
        value = xeasy16ToReal(exponent, mantissa)
      z.append(value)

    z.tofile(fpw)

  fpr.close()
  fpw.close()

def convertXeasyToAzara(xeasyParamFile, azaraDataFile, azaraParFile):

  xeasyParams = XeasyParams(xeasyParamFile)

  azaraParams = writeAzaraParFile(xeasyParams, azaraDataFile, azaraParFile)

  writeAzaraDataFile(xeasyParams, azaraParams)

if __name__ == '__main__':

  import sys
  import os

  def usage(msg):

    print 'Error: %s' % msg
    print 'Arguments: <XEASY param file> [<azaraDataFile> <azaraParFile> ]'
    print '     For example: HNCO.param'
    print '              or: HNCO.param HNCO.spc'
    print '              or: HNCO.param HNCO.spc HNCO.spc.par'
    # sys.exit(1)
    os._exit(1)

  nargs = len(sys.argv)
  if nargs < 2 or nargs > 4:
    usage('need 1 to 3 arguments')

  xeasyParamFile = sys.argv[1]

  if not xeasyParamFile.endswith('.param'):
    usage('first argument must be XEASY param file, so ending in ".param"')

  if nargs < 3:
    azaraDataFile = xeasyParamFile[:-5] + 'spc'
  else:
    azaraDataFile = sys.argv[2]

  if nargs < 4:
    azaraParFile = azaraDataFile + '.par'
  else:
    azaraParFile = sys.argv[3]

  convertXeasyToAzara(xeasyParamFile, azaraDataFile, azaraParFile)
