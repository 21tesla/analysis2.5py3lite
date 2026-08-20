import os

from memops.general.Implementation import ApiError

from ccp.format.spectra.params import AzaraParams, BrukerParams, FactorisedParams, FelixParams, NmrPipeParams, NmrViewParams, UcsfParams, VarianParams

def openSpectrum(fileName, nmrProject=None, experiment=None, experimentName=None, spectrumName=None):

  assert nmrProject or experiment

  # try different possible formats in (some) order

  for module in (NmrPipeParams, UcsfParams, NmrViewParams, BrukerParams, VarianParams, AzaraParams, FactorisedParams, FelixParams):
    className = module.__name__.split('.')[-1]
    clazz = getattr(module, className)
    try:
      params = clazz(fileName)
      break
    except:
      pass
  else:
    if not os.path.exists(fileName):
      raise IOError('"%s" does not exist' % fileName)

    raise IOError('Could not determine file type for "%s"' % fileName)

  # create experiment if need be

  if not experiment:
    if not experimentName:
      basename = os.path.basename
      if isinstance(params, BrukerParams.BrukerParams):
        dirname = os.path.dirname
        experimentName = basename(dirname(dirname(dirname(params.dataFile))))
      else:
        experimentName = basename(params.dataFile)
        n = experimentName.rfind('.')
        if n > 0:
          experimentName = experimentName[:n]
    if not experimentName:  # unlikely but play safe
      n = 1
      while experiment.findFirstExperiment(name='Expt%d' % n):
        n += 1
      experimentName = 'Expt%d' % n
    experiment = nmrProject.newExperiment(name=experimentName, numDim=params.ndim)
    
  # create spectrum

  if not spectrumName:
    n = 1
    while experiment.findFirstDataSource(name='%d' % n):
      n += 1
    spectrumName = '%d' % n

  spectrum = params.createDataSource(experiment, name=spectrumName)

  return spectrum

if __name__ == '__main__':

  import sys
  from memops.api.Implementation import MemopsRoot

  if len(sys.argv) != 2:
    print 'Need to specify file'
    sys.exit()

  fileName = sys.argv[1]

  project = MemopsRoot(name='test')
  nmrProject = project.newNmrProject(name='test')

  spectrum = openSpectrum(fileName, nmrProject=nmrProject)
