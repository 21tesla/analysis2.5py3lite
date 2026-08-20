from ccp.format.ansig.AnsigParse import parseAnsigFile, AnsigParseException

class AnsigResidue:

  def __init__(self, sequence, resName, resType, details=None):

    # attributes
    self.resName = resName
    self.resType = resType
    self.details = details

    # links
    self.sequence = sequence

class AnsigSequence:

  def __init__(self, project, name):

    # attributes
    self.name = name

    # links
    self.project = project
    self.residues = []

    project.sequence = self

  def addResidue(self, resName, resType, details=None):

    residue = AnsigResidue(self, resName, resType, details=details)
    self.residues.append(residue)

def readSequenceFile(project, directoryDict, verbose=False):

  if verbose:
    print 'Reading Ansig sequence file "%s"' % project.sequenceFile

  project.sequence = None

  def callback(key, values, comment):
    if key == 'sequence':
      if project.sequence:
        raise AnsigParseException('second sequence found')
      name = values[0]
      AnsigSequence(project, name)
      sequenceKeywordDict = {
        'residue': (False, 2, (None, None)),
      }
      return sequenceKeywordDict
    elif key == 'residue':
      project.sequence.addResidue(values[0], values[1], comment)

  sequenceFile = project.sequenceFile
  keywordDict = {
    'sequence': (True, 1, (None,)),
  }

  try:
    parseAnsigFile(sequenceFile, callback, keywordDict, directoryDict)
    if not project.sequence:
      raise IOError('no "sequence" statement found in file "%s"' % sequenceFile)
    if verbose:
      print 'Number of residues found = %s' % len(project.sequence.residues)
  except:
    project.sequence = None
    raise

