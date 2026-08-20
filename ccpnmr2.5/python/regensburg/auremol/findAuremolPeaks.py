from ccpnmr.analysis.core.ExperimentBasic import getPrimaryDataDimRef
from ccpnmr.analysis.core.PeakBasic import setManualPeakIntensity

# internal function
def setupCPeak(peak):

  if not hasattr(peak, 'cPeakList'): # should not exist but protect against this
    peak.cPeakList = peak.peakList.cPeakList
  
  if not hasattr(peak, 'cPeak'): # should not exist but protect against this
    peak.cPeak = peak.cPeakList.addPeak()

# internal function
def setupCPeakList(peakList):

  from ccpnmr.c.PeakList import PeakList as CPeakList

  if hasattr(peakList, 'cPeakList'): # should not exist but protect against this
    return

  spectrum = peakList.dataSource
    
  # TBD: below assumes all dataDims are freqDataDims
  npoints = [ dataDim.numPoints for dataDim in spectrum.sortedDataDims() ]
   
  peakList.cPeakList = CPeakList(npoints)
  peaks = peakList.sortedPeaks()
  for peak in peaks:
    setupCPeak(peak)

# internal function
def setupPeaks(peakList, data):

  # below assumes that all positions found are in ppm and in fundamental region
  spectrum = peakList.parent
  numDim = spectrum.numDim
  dataDims = spectrum.sortedDataDims()
  dataDimRefs = [getPrimaryDataDimRef(dataDim) for dataDim in dataDims]
  for (position, height, volume, quality) in data:
    peak = peakList.newPeak()
    peakDims = peak.sortedPeakDims()
    for i in range(numDim):
      peakDim = peakDims[i]
      dataDimRef = dataDimRefs[i]
      peakDim.dataDimRef = dataDimRef
      #peakDim.position = dataDimRef.valueToPoint(position[i])
      peakDim.value = position[i]
    setManualPeakIntensity(peak, height, 'height')
    setManualPeakIntensity(peak, volume, 'volume')
    if quality >= 0: # = -1 if not calculated
      peak.figOfMerit = quality  # assumes that quality is between and 0 and 1

"""
Notes:

mode = 0: positive and negative
mode = 1: positive only
mode = 2: negative only

probability = threshold for keeping peaks (0 to 1) (e.g. 0.999)

number = keep peaks whose probability is at or above that for this peak (can be greater than number of peaks; can consider to be upper limit on number of peaks expected) (e.g. 4000)

seglevel = integration depth in respect to maximum intensity (for calculating volume) (0.001 to 1) (e.g. 0.1)

For bayes (when ready):
percentNoise = percent noise (which means what???) (0 to 100) (e.g. 30)
percentSignal = percent signal (which means what???) (0 to 100) (e.g. 30)
"""

def findAuremolPeaksThreshold(argServer = None, spectrum = None, mode = 1, useAutoThreshold = 1, threshold = 100000, seglevel = 0.1):

  from Auremol import findPeaksThreshold

  print 'Starting findAuremolPeaksThreshold'

  if not spectrum:
    spectrum = argServer.getSpectrum()

  dataStore = spectrum.dataStore
  if not dataStore:
    errMsg = 'No dataStore found for spectrum, so cannot find peaks'
    if argServer:
      argServer.messageReporter.showError('No dataStore', errMsg)
    else:
      print errMsg
    return None

  spectrumPath = dataStore.fullPath

  peakList = spectrum.newPeakList()
  if not hasattr(peakList, 'cPeakList'):
    setupCPeakList(peakList)

  cPeakList = peakList.cPeakList

  print 'About to call findPeaksThreshold'
  data = findPeaksThreshold(cPeakList, spectrumPath, mode, useAutoThreshold, threshold, seglevel)

  print 'About to call setupPeaks'
  setupPeaks(peakList, data)

  print 'Ending findAuremolPeaksThreshold with %d peaks' % len(peakList.peaks)

  return peakList

def findAuremolPeaksAdaptive(argServer = None, spectrum = None, mode = 1, number = 4000, seglevel = 0.1):

  from Auremol import findPeaksAdaptive

  if not spectrum:
    spectrum = argServer.getSpectrum()

  dataStore = spectrum.dataStore
  if not dataStore:
    errMsg = 'No dataStore found for spectrum, so cannot find peaks'
    if argServer:
      argServer.messageReporter.showError('No dataStore', errMsg)
    else:
      print errMsg
    return None

  print 'Starting findAuremolPeaksAdaptive'

  spectrumPath = dataStore.fullPath

  peakList = spectrum.newPeakList()
  if not hasattr(peakList, 'cPeakList'):
    setupCPeakList(peakList)

  cPeakList = peakList.cPeakList

  print 'About to call findPeaksAdaptive'
  data = findPeaksAdaptive(cPeakList, spectrumPath, mode, number, seglevel)

  print 'About to call setupPeaks'
  setupPeaks(peakList, data)

  print 'Ending findAuremolPeaksAdaptive with %d peaks' % len(peakList.peaks)

  return peakList

if __name__ == '__main__':

  import sys

  n = len(sys.argv)

  if n < 2:
    print 'Correct arguments: projectdir [exptName [spectrumName]]'
    sys.exit()

  projectDir = sys.argv[1]
  from memops.general.Io import loadProject

  project = loadProject(projectDir)

  nmrProject = project.findFirstNmrProject()

  expt = None
  if n > 2:
    exptName = sys.argv[2]
    expt = nmrProject.findFirstExperiment(name=exptName)
    if not expt:
      print 'Warning: experiment "%s" not found, using first one'

  if not expt:
    expt = nmrProject.findFirstExperiment()

  if not expt:
    raise Exception('no experiment found')

  spectrum = None
  if n > 3:
    spectrumName = sys.argv[3]
    spectrum = expt.findFirstDataSource(name=spectrumName)
    if not spectrum:
      print 'Warning: spectrum "%s" not found, using first one'

  if not spectrum:
    spectrum = expt.findFirstDataSource()

  if not spectrum:
    raise Exception('no spectrum found')

  findAuremolPeaksThreshold(argServer=None, spectrum=spectrum, useAutoThreshold=0, threshold=500000)

  project.saveModified()

  print 'Exiting script'

