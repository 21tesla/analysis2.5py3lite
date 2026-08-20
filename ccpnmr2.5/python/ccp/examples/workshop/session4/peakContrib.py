def main(peakList):

  dataSource = peakList.parent
  expt = dataSource.parent
  nmrProject = expt.parent

  # create Peak (this automtatically also creates peak.peakDims)
  peak = peakList.newPeak()

  # create PeakContrib
  peakContrib = peak.newPeakContrib()

  resonance = nmrProject.newResonance(isotopeCode='1H')

  for peakDim in peak.peakDims:
    peakDimContrib = peakDim.newPeakDimContrib(resonance=resonance, peakContribs=(peakContrib,))

  print 'dimension of experiment = %d, dataSource = %d' % (expt.numDim, dataSource.numDim)
  print 'peakContrib has %d peakDimContribs' % len(peakContrib.peakDimContribs)

  #peakList.root.saveModified()

if __name__ == '__main__':

  import sys

  if len(sys.argv) != 6:
    print 'need to specify: rootName, nmrProjectName, exptSerial, dataSourceSerial, peakListSerial'
    sys.exit()

  from memops.general.Io import loadProject
  repositoryPath = sys.argv[1]
  root = loadProject(path=repositoryPath)

  nmrProject = root.findFirstNmrProject(name=sys.argv[2])
  if not nmrProject:
    print 'could not find nmrProject %s' % sys.argv[2]
    sys.exit()

  expt = nmrProject.findFirstExperiment(serial=int(sys.argv[3]))
  if not expt:
    print 'could not find experiment %s' % sys.argv[3]
    sys.exit()

  dataSource = expt.findFirstDataSource(serial=int(sys.argv[4]))
  if not dataSource:
    print 'could not find dataSource %s' % sys.argv[4]
    sys.exit()

  peakList = dataSource.findFirstPeakList(serial=int(sys.argv[5]))
  if not dataSource:
    print 'could not find peakList %s' % sys.argv[5]
    sys.exit()

  main(peakList)
