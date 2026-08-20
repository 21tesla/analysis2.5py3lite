
from memops.general.Implementation import ApiError

from ccp.format.azara.processingParsIO import AzaraProcessingParsFile

from ccp.format.spectra.params.ExternalParams import ExternalParams, fidDimType

class AzaraParams(ExternalParams):

  format = 'Azara'

  def __init__(self, parFile, **kw):

    self.parFile = parFile

    ExternalParams.__init__(self, **kw)

  # ExternalParams requires this to be defined
  def parseFile(self):

    try:
      procParData = AzaraProcessingParsFile(self.parFile)
      procParData.read()
    except IOError, e:
      raise ApiError(str(e))

    aPars = procParData.aPars
    self.setAttrs(aPars)
    
    # correct for FID data - 
    if set(self.block) == set((None,)):
      # assume that data without block information are unblocked and FID
      print 'No blocking found - setting to Unblocked FID data'
      npts = self.npts[0]
      for ii in range(self.ndim):
        self.block[ii] = npts
        self.dimType[ii] = fidDimType
        npts = 1

    if aPars.has_key('sampledValues'):
      sampledValuesList = aPars['sampledValues']
      if sampledValuesList:
        for n in range(len(sampledValuesList)):
          sampledValues = sampledValuesList[n]
          if sampledValues:
            self.setSampledDim(n, sampledValues)

  def writeParFile(self, parFile = None):

    if parFile is None:
      parFile = self.parFile

    fp = open(parFile, 'w')

    fp.write('ndim %d\n' % self.ndim)
    fp.write('file %s\n' % self.dataFile)

    if self.integer:
      fp.write('int\n')

    if self.big_endian:
      fp.write('big_endian\n')
    else:
      fp.write('little_endian\n')

    for i in range(self.ndim):
      fp.write('\n')
      fp.write('dim %d\n' % (i+1))
      fp.write('npts %d\n' % self.npts[i])
      fp.write('block %d\n' % self.block[i])
      fp.write('sw %4.3f\n' % self.sw[i])
      fp.write('sf %4.3f\n' % self.sf[i])
      fp.write('refppm %5.4f\n' % self.refppm[i])
      fp.write('refpt %5.4f\n' % self.refpt[i])
      fp.write('nuc %s\n' % self.nuc[i])
      if self.pointValues[i]:
        fp.write('params')
        for value in self.pointValues[i]:
          fp.write(' %5.4f' % value)
        fp.write('\n')

    fp.close()

if (__name__ == '__main__'):

  import sys
  if (len(sys.argv) != 2):
    print 'Error: required argument: <parFile>'
    sys.exit(1)

  parFile = sys.argv[1]
  params = AzaraParams(parFile)

  from memops.api.Implementation import MemopsRoot, Url
  from ccp.util.Spectrum import createExperiment

  r = MemopsRoot()
  l = r.newDataLocationStore(name='testDLS')
  d = l.newDataUrl(name='testDataUrl', url=Url(path='/Users/wb104/edl387'))
  n = r.newNmrProject(name='testNmrProject')
  e = createExperiment(n, name='testExpt', numDim=params.ndim, sf=params.sf, isotopeCodes=params.nuc)
  s = params.createDataSource(e, name='testSpectrum', dataUrl=d)

