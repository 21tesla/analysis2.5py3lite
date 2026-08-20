import os

from memops.universal.Io import joinPath, normalisePath

from ccp.format.ansig.AnsigParse import parseAnsigFile, AnsigParseException
from ccp.format.ansig.AnsigSequence import readSequenceFile
from ccp.format.ansig.AnsigSpectrum import readSpectraFile
from ccp.format.ansig.AnsigCrosspeak import readCrosspeakFile #, setupCrosspeaks, mergeResonances

class AnsigProject:

  def __init__(self, projectFile):

    # attributes
    self.projectFile = normalisePath(projectFile, makeAbsolute=True)
    self.sequenceFile = None
    self.spectraFile = None
    self.crosspeaksFile = None
    self.binaryformat = 'bigendian'

    # links
    self.sequence = None
    self.spectra = []
    self.crosspeaks = []
    self.resonances = set()

  def findSpectrum(self, name):

    for spectrum in self.spectra:
      if spectrum.name == name:
        return spectrum

    return None

def readProjectFile(projectFile, verbose=False):

  if verbose:
    print 'Reading Ansig project (ctr) file "%s"' % projectFile

  project = AnsigProject(projectFile)
  project.foundControl = False

  def callback(key, values, comment):
    if key == 'control':
      if project.foundControl:
        raise AnsigParseException('control found twice')
      project.foundControl = True
      controlKeywordDict = {
        'help_dir': (False, 1, (None,)),
        'html_viewer': (False, 1, (None,)),
        'dictionary': (False, 1, (None,)),
        'gr_translate': (False, 1, (None,)),
        'sequence': (False, 1, (None,)),
        'spectra': (False, 1, (None,)),
        'crosspeaks': (False, 1, (None,)),
        'backup': (False, 1, (None,)),
        'initialize': (False, 1, (None,)),
        'exit': (False, 1, (None,)),
        'binaryformat': (False, 1, (None,)),
        'font': (False, 1, (None,)),
      }
      return controlKeywordDict
    elif key in ('sequence', 'spectra', 'crosspeaks', 'binaryformat'):
      value = values[0]
      if key != 'binaryformat':
        key = key + 'File'
        if not os.path.isabs(value):
          base = os.path.dirname(project.projectFile)
          value = joinPath(base, value)
      setattr(project, key, value)

  keywordDict = {
    'control': (True, 0, ()),
  }
  directoryDict = {}
  try:
    parseAnsigFile(projectFile, callback, keywordDict,
                   directoryDict, canDefineDirectory=True)
    if not project.foundControl:
      raise IOError('no "control" statement found in file "%s"' % projectFile)
  finally:
    del project.foundControl

  if project.sequenceFile:
    readSequenceFile(project, directoryDict, verbose=verbose)

  if project.spectraFile:
    readSpectraFile(project, directoryDict, verbose=verbose)

  if project.crosspeaksFile:
    readCrosspeakFile(project, verbose=verbose)
    #setupCrosspeaks(project)
    #mergeResonances(project, verbose=verbose)

  return project

if __name__ == '__main__':

  import sys

  if len(sys.argv) != 2:
    print 'need argument: ansig project (ctr) file'
    sys.exit()

  projectFile = sys.argv[1]
  project = readProjectFile(projectFile, verbose=True)

