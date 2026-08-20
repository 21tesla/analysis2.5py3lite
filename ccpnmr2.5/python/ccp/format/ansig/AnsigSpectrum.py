import os
import re

from ccp.format.ansig.AnsigParse import parseAnsigFile, AnsigParseException, ANSIG_INT, ANSIG_FLOAT

from ccp.format.azara.processingParsIO import AzaraProcessingParsFile

ndimRe = re.compile("(\d+)D$")
dimRe = re.compile("F(\d+)$")

class AnsigSpectrum:

  def __init__(self, project, ndim, name):

    # attributes
    self.ndim = ndim
    self.name = name
    self.solvent = None
    self.temperature = None
    self.pH = None
    self.equivalent = None # the name of another spectrum which this is equivalent to
    self.scale = None
    self.linewidth = ndim*[None]
    self.transfers = []

    # attributes to do with par file
    self.haveReadParFile = False
    self.parFile = None
    self.head = None
    self.integer = False
    self.bigEndian = (project.binaryformat == 'bigendian')
    self.dataFile = None
    self.npts = ndim*[None]
    self.block = ndim*[None]
    self.sw = ndim*[None]
    self.sf = ndim*[None]
    self.refppm = ndim*[None]
    self.refpt = ndim*[None]
    self.nuc = ndim*[None]

    # links
    self.project = project
    self.crosspeaks = []

    project.spectra.append(self)

  def readParFile(self, parFile):

    self.parFile = parFile

    procParData = AzaraProcessingParsFile(parFile)
    procParData.read()

    aPars = procParData.aPars
    if aPars['ndim'] != self.ndim:
      raise IOError('spectrum "%s": ansig spd file has ndim=%d but par file "%s" has ndim=%d' % (self.name, self.ndim, parFile, aPars['ndim']))

    for key in aPars.keys():
      self.setParAttr(key, aPars[key])

    self.haveReadParFile = True

  def setParAttr(self, key, value):

    # Ignore keys that do not correspond to spectrum attributes
    # Note: so implicitly assuming here that Ansig binaryformat is more correct than what par_file says
    # because bigEndian is not a AzaraProcessingParsFile / Ansig attribute (big_endian is)

    # TBD: should change 'file' in AzaraProcessingParsFile some day
    if key == 'file':
      key = 'dataFile'
    if key in ('dataFile', 'head', 'integer'):
      if getattr(self, key) is None:
        setattr(self, key, value)
    elif key in ('npts', 'block', 'sw', 'sf', 'refppm', 'refpt', 'nuc'):
      v = getattr(self, key)
      for n in range(self.ndim):
        if v[n] is None:
          v[n] = value[n]

def readSpectraFile(project, directoryDict, verbose=False):

  if verbose:
    print 'Reading Ansig spectra file "%s"' % project.spectraFile

  project.spectra = []

  def getNdim(value):

    obj = ndimRe.match(value)
    if not obj:
      raise AnsigParseException('ndim = "%s", must be integer followed by "D" (e.g. "3D")' % value)
    ndim = int(obj.group(1))

    return ndim

  def getDim(value, ndim):

    obj = dimRe.match(value)
    if not obj:
      raise AnsigParseException('dim = "%s", must be "F" followed by integer (e.g. "F2")' % value)
    dim = int(obj.group(1))

    if dim < 1 or dim > ndim:
      raise AnsigParseException('dim = "%s", must be between "F1" and "F%d"' % (value, ndim))

    return ndim-dim+1  # Ansig dimensions are backwards

  def getFloat(value, errorMsg):

    try:
      x = float(value)
    except ValueError, e:
      raise AnsigParseException(errorMsg % value)

    return x

  def callback(key, values, comment):

    if key in ('picket_fence_bug', 'ppm_range', 'end_spectrum'):
      return

    if key == 'spectrum':
      ndim = getNdim(values[0])
      name = values[1]
     
      spectrum = AnsigSpectrum(project, ndim, name)
      spectrumKeywordDict = {
        'expnucleus': (False, (2, 5), (None, None, None, None, None)),
        'experiment': (False, (1, 3, 5), (None, None, None, None, ANSIG_FLOAT)),
        ###'symmetry': (False, ndim, ndim*(None,)), # have ( and ) around strings and not sure if white space is allowed or not
        'symmetry': (False, -1, (None,)),
        'solvent': (False, 1, (None,)),
        'temperature': (False, 1, (ANSIG_FLOAT,)),
        'pH': (False, 1, (ANSIG_FLOAT,)),
        'equivalent': (False, 1, (None,)),
        'scale': (False, 1, (ANSIG_FLOAT,)),
        'contour_file': (False, -1, (None,)),
        'matrix_file': (True, 0, ()),
      }
      return spectrumKeywordDict

    spectrum = project.spectra[-1]
    ndim = spectrum.ndim

    if key == 'matrix_file':
      rr = tuple(range(2,2*(ndim+1),2))
      matrixKeywordDict = {
        'parameter_file': (False, 1, (None,)),
        'ndim': (False, 1, (ANSIG_INT,)),
        'file': (False, 1, (None,)),
        'head': (False, 1, (ANSIG_INT,)),
        'deflate': (False, 1, (ANSIG_FLOAT,)),
        'reflate': (False, 1, (ANSIG_FLOAT,)),
        'dim': (False, 1, (ANSIG_INT,)),
        'npts': (False, 1, (ANSIG_INT,)),
        'block': (False, 1, (ANSIG_INT,)),
        'sw': (False, 1, (ANSIG_FLOAT,)),
        'sf': (False, 1, (ANSIG_FLOAT,)),
        'refppm': (False, 1, (ANSIG_FLOAT,)),
        'refpt': (False, 1, (ANSIG_FLOAT,)),
        'nuc': (False, 1, (None,)),
        'aliased': (False, 1, (None,)),
        'folded': (False, 1, (None,)),
        'linewidth': (False, rr, ndim*(None, ANSIG_FLOAT)),
        'noise_level': (False, 1, (ANSIG_FLOAT,)),
        'contour_levels': (False, -1, (ANSIG_FLOAT,)),
        'contour_base': (False, 3, (ANSIG_FLOAT, ANSIG_FLOAT, ANSIG_INT)),
      }
      spectrum.currentDim = None
      return matrixKeywordDict

    elif key == 'end_matrix_file':
      del spectrum.currentDim

    # spectrum attributes

    elif key == 'expnucleus':
      dim = getDim(values[0], ndim)
      nuc = values[1]
      spectrum.nuc[dim-1] = nuc
      if len(values) == 5:
        pass # TBD: "associate other_spectrum, other_dim" gives mapping to some other dataDim

    elif key == 'experiment':
      if len(values) == 5:
        dim1 = getDim(values[0], ndim)
        dim2 = getDim(values[1], ndim)
        expName = values[2]
        if values[3] != 'mixing_time':
          raise AnsigParseException('fourth argument = "%s", must be "mixing_time"' % values[3])
        mixing_time = getFloat(values[4], 'mixing_time = "%s", must be float')
      elif len(values) == 3:
        if ndim == 2 and values[1] == 'mixing_time':
          dim1 = 1
          dim2 = 2
          expName = values[0]
          mixing_time = getFloat(values[2], 'mixing_time = "%s", must be float')
        else:
          dim1 = getDim(values[0], ndim)
          dim2 = getDim(values[1], ndim)
          expName = values[2]
          mixing_time = None
      else: # len(values) = 1
        if ndim != 2:
          raise AnsigParseException('ndim != 2 but unspecified dims in "experiment" statement')
        dim1 = 1
        dim2 = 2
        expName = values[0]
        mixing_time = None
      spectrum.transfers.append((dim1, dim2, expName, mixing_time))

    elif key == 'symmetry':
      pass  # TBD

    elif key in ('solvent', 'temperature', 'pH', 'equivalent', 'scale'):
      setattr(spectrum, key, values[0])

    elif key == 'contour_file':
      pass  # TBD

    # matrix_file attributes

    elif key == 'parameter_file':
      parFile = values[0]
      spectrum.parFile = parFile
      if os.path.exists(parFile):
        try:
          readParFile(parFile)
        except IOError, e:
          raise AnsigParseException(str(e))

    elif key == 'ndim':
      if values[0] != ndim:
        raise AnsigParseException('spectrum "%s": ansig spd file has ndim=%d but also have set ndim=%d' % (spectrum.name, ndim, values[0]))

    elif key == 'file':
      spectrum.dataFile = values[0]

    elif key == 'head':
      spectrum.head = values[0]

    elif key in ('deflate', 'reflate'):
      pass  # TBD: not sure what to do about this (but should not be common)

    elif key == 'dim':
      spectrum.currentDim = values[0]

    elif key in ('npts', 'block', 'sw', 'sf', 'refppm', 'refpt'):
      x = getattr(spectrum, key)
      x[spectrum.currentDim] = values[0]

    elif key == 'nuc':
      if spectrum.nuc[spectrum.currentDim] == None:
        spectrum.nuc[spectrum.currentDim] = values[0]

    elif key in ('aliased', 'folded', 'noise_level', 'contour_levels', 'contour_base'):
      pass # TBD: should one do something with aliased or folded?

    elif key == 'linewidth':
      for n in range(0, len(values), 2):
        dim = getDim(values[n], ndim)
        spectrum.linewidth[dim-1] = values[n+1]

  spectraFile = project.spectraFile
  keywordDict = {
    'picket_fence_bug': (False, 0, ()),
    'ppm_range': (False, 3, (None, ANSIG_FLOAT, ANSIG_FLOAT)),
    'spectrum': (True, 2, (None, None)),
  }

  try:
    parseAnsigFile(spectraFile, callback, keywordDict, directoryDict)
    if not project.spectra:
      raise IOError('no "spectrum" statement found in file "%s"' % spectraFile)
    if verbose:
      print 'Number of spectra found = %d' % len(project.spectra)
  except:
    project.spectra = []
    raise

