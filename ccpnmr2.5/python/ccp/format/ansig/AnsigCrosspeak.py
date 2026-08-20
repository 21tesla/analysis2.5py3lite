import array

class AnsigResonance:

  def __init__(self, project):

    # links
    self.project = project
    project.resonances.add(self)
    self.crosspeakDims = set()
    self.residue = None
    self.atname = None

  def mergeWith(self, other):

    if self == other:
      return

    self.crosspeakDims.update(other.crosspeakDims)
    for (crosspeak, dim) in other.crosspeakDims:
      crosspeak.resonances[dim] = self
    self.project.resonances.remove(other)

class AnsigCrosspeak:

  def __init__(self, project, ndim, number, spectrumName, shifts, intensity, resid, resname, atname,
      symmetryConnection, prevConnection, nextConnection, prevCorrConnection, nextCorrConnection, spectrum=None):

    # attributes
    self.ndim = ndim
    self.number = number
    self.spectrumName = spectrumName
    self.shifts = shifts
    self.intensity = intensity
    self.resid = resid
    self.resname = resname
    self.atname = atname
    self.symmetryConnection = symmetryConnection
    # the next two seem to be a linked list of crosspeaks that are really the same in one dim
    # so this is a list of length ndim
    self.prevConnection = prevConnection
    self.nextConnection = nextConnection
    # the next two seem to be a linked list of crosspeaks that are really the same in all dims
    self.prevCorrConnection = prevCorrConnection
    self.nextCorrConnection = nextCorrConnection

    # links
    self.project = project
    self.spectrum = None
    self.symmetryCrosspeak = None
    self.resonances = ndim*[None]

    self.setSpectrum(spectrum)

    project.crosspeaks.append(self)

  def setSpectrum(self, spectrum):

    if self.spectrum:
      spectrum.crosspeaks.remove(self)
    self.spectrum = spectrum
    if spectrum:
      spectrum.crosspeaks.append(self)

  def setResonance(self, resonance, dim):

    self.resonances[dim] = resonance
    resonance.crosspeakDims.add((self, dim))

def readCrosspeakFile(project, verbose=False):

  def cleanString(s):
    s = s.replace('\x00', ' ')
    s = s.strip()
    return s

  if verbose:
    print 'Reading Ansig crosspeaks file "%s"' % project.crosspeaksFile

  fp = open(project.crosspeaksFile, 'rb')
  s = fp.read()
  fp.close()

  #title = s[:12]

  x = array.array('i')
  y = array.array('f')
  x.fromstring(s)
  y.fromstring(s)

  ndim = x[3]

  if ndim < 1 or ndim > 4:
    x.byteswap()
    nd = x[3]
    if nd >= 1 and nd <= 4:
      ndim = nd
      y.byteswap()
    else: # put back, but probably in trouble
      x.byteswap()

  npeaks = x[4]
  recordLength = 6*ndim + 7

  if verbose:
    deletedCount = 0
    unknownSpectrumCount = {}

  project.crosspeaks = []
  for i in range(npeaks):
    ind = (i+2) * recordLength
    shifts = y[ind:ind+ndim]
    ind += ndim
    intensity = y[ind]
    ind += 1
    spectrumName = cleanString(s[4*ind:4*(ind+3)])
    ind += 3
    connection = x[ind:ind+2*ndim+3]
    ind += 2*ndim+3
    resid = []
    for j in range(ndim):
      resid.append(cleanString(s[4*ind:4*(ind+1)]))
      ind += 1
    resname = []
    for j in range(ndim):
      resname.append(cleanString(s[4*ind:4*(ind+1)]))
      ind += 1
    atname = []
    for j in range(ndim):
      atname.append(cleanString(s[4*ind:4*(ind+1)]))
      ind += 1

    nd = ndim
    spectrum = None
    if spectrumName:
      spectrum = project.findSpectrum(spectrumName)
      if spectrum:
        nd = spectrum.ndim
      elif verbose:
        if not unknownSpectrumCount.has_key(spectrumName):
          unknownSpectrumCount[spectrumName] = 0
        unknownSpectrumCount[spectrumName] += 1
    else:
      if verbose:
        deletedCount += 1
      continue

    shifts = shifts[:nd]
    resid = resid[:nd]
    resname = resname[:nd]
    atname = atname[:nd]
    symmetryConnection = connection[0]
    prevConnection = nd * [None]
    nextConnection = nd * [None]
    prevCorrConnection = None
    nextCorrConnection = None
    for j in range(nd):
      prevConnection[j] = connection[2*j+1]
      nextConnection[j] = connection[2*j+2]
    prevCorrConnection = connection[-2]
    nextCorrConnection = connection[-1]

    # Ansig has dimensions backwards
    shifts.reverse()
    resid.reverse()
    resname.reverse()
    atname.reverse()
    prevConnection.reverse()
    nextConnection.reverse()

    AnsigCrosspeak(project, nd, i+1, spectrumName, shifts, intensity, resid, resname, atname,
      symmetryConnection, prevConnection, nextConnection, prevCorrConnection, nextCorrConnection, spectrum)

  if verbose:
    print 'Number of deleted peaks = %d (ignoring these)' % deletedCount
    print 'Number of non-deleted peaks = %d (retaining these)' % len(project.crosspeaks)
    for key in sorted(unknownSpectrumCount.keys()):
      print 'Unknown spectrum "%s" with %d crosspeaks' % (key, unknownSpectrumCount[key])

""" below was copied in but is far from ready
def setupCrosspeaks(project):

  crosspeaks = project.crosspeaks

  # setup symmetryCrosspeaks
  for crosspeak in crosspeaks:
    if crosspeak.symmetryConnection:
      crosspeak.symmetryCrosspeak = crosspeaks[crosspeak.symmetryConnection-1]

  # setup crosspeak resonances
  crosspeakDimMap = {}
  for crosspeak in crosspeaks:
    for dim in range(crosspeak.ndim):
      resid = crosspeak.resid[dim]
      atname = crosspeak.atname[dim]
      if resid and atname:
        crosspeakDims = resonanceMap.get((resid, atname))
        if not crosspeakDims:
          crosspeakDimMap[(resid, atname)] = crosspeakDims = set()
        crosspeakDims.add((crosspeak, dim))
        inverseMap[(crosspeak, dim)] = crosspeakDims
  
  project.resonances = set()
  resonanceMap = {}
  for crosspeak in crosspeaks:
    for i in range(crosspeak.ndim):
      currResonance = crosspeak.resonances[i]
      resid = crosspeak.resid[dim]
      atname = crosspeak.atname[dim]
      if resid and atname:
        key1 = (resid, atname)
        resonance1 = resonanceMap.get(key1)
      else:
        key1 = None:
        resonance1 = None:
      key2 = (crosspeak, i)
      resonance2 = crosspeak.resonances[i]
      if resonance1 and resonance2:
        resonance1.mergeWith(resonance2)
        resonance = resonance1
      elif resonance1:
        resonance = resonance1
      elif resonance2:
        resonance = resonance2
      else:
        resonance = AnsigResonance(project)
      if key1:
        resonanceMap[key1] = resonance
      crosspeak.resonances[i] = resonance

        # work in prev direction
        xpk = crosspeak
        j = i
        xpk.setResonance(resonance, j)
        while xpk.prevConnection[j]:
          xpk2 = crosspeaks[xpk.prevConnection[j]-1]
          for k in range(xpk2.ndim):
            if xpk2.nextConnection[k] == xpk.number:
              break
          else:
            # should not be here, it means linked list is broken
            # so just stop the while loop
            break
          xpk = xpk2
          j = k
          xpk.setResonance(resonance, j)

        # work in next direction
        xpk = crosspeak
        j = i
        while xpk.nextConnection[j]:
          xpk2 = crosspeaks[xpk.nextConnection[j]-1]
          for k in range(xpk2.ndim):
            if xpk2.prevConnection[k] == xpk.number:
              break
          else:
            # should not be here, it means linked list is broken
            # so just stop the while loop
            break
          xpk = xpk2
          j = k
          xpk.setResonance(resonance, j)

  # TBD: setup corr peak (i.e. same peak??) resonances
  # not really sure how ansig means these to be treated

def mergeResonances(project, verbose=False):

  if verbose:
    print 'At start of mergeResonances have %d resonances' % len(project.resonances)
    nassignedCrosspeakDims = 0

  resonanceMap = {}
  for n, crosspeak in enumerate(project.crosspeaks):
    if (n%100 == 0): print 'HERE222', n
    for dim in range(crosspeak.ndim):
      resid = crosspeak.resid[dim]
      atname = crosspeak.atname[dim]
      if resid and atname:
        if verbose:
          nassignedCrosspeakDims += 1
        resonances = resonanceMap.get((resid, atname))
        if not resonances:
          resonanceMap[(resid, atname)] = resonances = set()
        resonances.add((crosspeak, dim))

  if verbose:
    print 'Number of assigned crosspeakDims = %d' % nassignedCrosspeakDims
    print 'At end of mergeResonances have %d resonances' % len(project.resonances)
"""
