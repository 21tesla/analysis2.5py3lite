"""Functions for calculating axisCodes for NmrExpPrototypes, adn necessary utilities.
For normal cases, use only refExpFimRefCodeMap function, and get axisCodes from
refExpDimRefs and the map.

"""
#=========================================================================================
# Licence, Reference and Credits
#=========================================================================================
__copyright__ = "Copyright (C) CCPN project (www.ccpn.ac.uk) 2014 - $Date: 2014-11-26 17:28:10 +0000 (Wed, 26 Nov 2014) $"
__credits__ = "Wayne Boucher, Rasmus H Fogh, Simon Skinner, Geerten Vuister"
__license__ = ("CCPN license. See www.ccpn.ac.uk/industry/licencefees"
               "or ccpncore.memops.Credits.CcpnLicense for license text")
__reference__ = ("For publications, please use reference from www.ccpn.ac.uk/industry/licencefees"
                 " or ccpncore.memops.Credits.CcpNmrReference")

#=========================================================================================
# Last code modification:
#=========================================================================================
__author__ = "$Author: rhfogh $"
__date__ = "$Date: 2014-11-26 17:28:10 +0000 (Wed, 26 Nov 2014) $"
__version__ = "$Revision: 7915 $"

#=========================================================================================
# Start of code
#=========================================================================================

def resetAllAxisCodes(nmrProject):
  """Reset all axisCodes (ExpDimRef.name) in project to be unique, match the isotope,
  and match the standard Prototype where a prototrype is known"""

  stdCodeMap = refExpDimRefCodeMap(nmrProject.root)

  for experiment in nmrProject.sortedExperiments():

    if experiment.refExperiment is None:
      # No prototype - just use nucleus as axis code

      foundCodes = {}
      for expDim in experiment.expDims:
        for expDimRef in expDim.expDimRefs:
          isotopeCodes = expDimRef.isotopeCodes

          # Get axis code
          if expDimRef.measurementype.lower() == 'shift' and len(isotopeCodes) == 1:
            # Normal shift axis, use nucleus to set
            code = isotope2Nucleus(isotopeCodes[0])
          else:
            # Different type of axis.
            # For simplicity set to axis_1, axis_2, ... for now
            print ("WARNING, non-std axis - axisCode set to 'axis'")
            code = 'axis'
            foundCodes[code] = 0

          # add index for duplicates
          indx = foundCodes.get(code)
          if indx is None:
            foundCodes[code] = 0
          else:
            indx += 1
            foundCodes[code] = indx
            code = '%s_%s' % (code, str(indx))

          # Set the attribute
          expDimRef.name = code


    else:
      # We have a prototype - use standard axisCode map
      for expDim in experiment.expDims:
        for expDimRef in expDim.expDimRefs:
          expDimRef.name = stdCodeMap[expDimRef.refExpDimRef]



def refExpDimRefCodeMap(project):
  """get RefExpDimRef: axisCode map for all NmExpPrototypes in project"""
  result = {}

  for nxp in project.sortedNmrExpPrototypes():
    for isReversed in False, True:
      refExperiments = nxp.findAllRefExperiments(isReversed=isReversed)
      if refExperiments:
        measurementMap = _measurementCodeMap(nxp, forReversed=isReversed)
        for re in refExperiments:
          for red in re.refExpDims:
            for redr in red.refExpDimRefs:
              result[redr] = measurementMap[redr.expMeasurement]
  #
  return result


def _measurementCodeMap(nmrExpPrototype, forReversed=False):
  """get expMeasurement:axisCode map"""
  result = {}
  measurements = _orderedMeasurements(nmrExpPrototype, forReversed=forReversed)

  foundCodes = {}
  # get axisCodes per expMeasurement
  for measurement in measurements:
    code = getAxisCode(measurement)
    indx = foundCodes.get(code)
    if indx is None:
      foundCodes[code] = 0
    else:
      indx += 1
      foundCodes[code] = indx
      code = '%s_%s' % (code, str(indx))
    result[measurement] = code
  #
  return result


def _connectedShiftMeasurements(expMeasurement):
  """find expMeasurements that are of type Shift, natch a HX bond, and onebonded to the input measurement"""
  result = []

  nmrExpPrototype = expMeasurement.nmrExpPrototype

  if expMeasurement.measurementType.lower() != 'shift':
    return result

  atomSites = expMeasurement.atomSites
  if len(atomSites) != 1:
    print ("WARNING%s Shift must have single AtomSite, has: %s"
           % (expMeasurement.nmrExpPrototype.name, [x.name for x in atomSites]))
    if not atomSites:
      return result

  expSite = atomSites[0]
  expIsotope = expSite.isotopeCode
  sites = []
  if expIsotope in ('1H', '19F'):
    for et in expSite.findAllExpTransfers(transferType='onebond'):
      for site in et.atomSites:
        if site is not expSite:
          if site.isotopeCode in ('13C', '15N') and site.maxNumber != 0:
            sites.append(site)
          break

  elif expIsotope in ('13C', '15N'):
    for et in expSite.findAllExpTransfers(transferType='onebond'):
      for site in et.atomSites:
        if site is not expSite:
          if site.isotopeCode in ('1H', '19F') and site.maxNumber != 0:
            sites.append(site)
          break

  result = [x for x in nmrExpPrototype.expMeasurements if x.measurementType.lower() == 'shift'
            and x.findFirstAtomSite() in sites]
  #
  return result

def _orderedMeasurements(nmrExpPrototype, forReversed=False):
  """get ExpMeasurements in order: acquisition last, furthest from acquisition last,
  connected measurements grouped together. If forReversed, get with acquisition last
  (for reversed experiments"""

  ll = [tt[1].expMeasurement
        for tt in sorted((x.stepNumber,x)
                         for x in nmrExpPrototype.findFirstExpGraph().expSteps)]
  if not forReversed:
    ll.reverse()

  measurements = []
  for measurement in ll:
    if measurement not in measurements:
      measurements.append(measurement)
      for me in _connectedShiftMeasurements(measurement):
        if me not in measurements:
          measurements.append(me)

  # we do not search all expGraphs, so just add missing measurements to the end
  # this is a heuristic, not a rigid ordering
  for measurement in nmrExpPrototype.expMeasurements:
    if measurement not in measurements:
      measurements.append(measurement)
      for me in _connectedShiftMeasurements(measurement):
        if me not in measurements:
          measurements.append(me)
  #
  return measurements


def getAxisCode(expMeasurement):
  """Get raw expMeasurement axisCode (without number suffixes) from NmrExpPrototype.ExpMeasurement"""
  tagMapping = {
  'shift':'shift',
  'jcoupling':'J',
  'mqshift':'MQ',
  'rdc':'RDC',
  'shiftanisotropy':'ANISO',
  'troesy':'TROESY',
  'dipolarcoupling':'DIPOLAR',
  't1':'delay',
  't2':'delay',
  't1rho':'delay',
  't1zz':'delay'
  }


  em = expMeasurement
  emType = em.measurementType.lower()

  tag = tagMapping.get(emType, emType)
  if tag == 'delay':
    result = tag
  elif tag == 'shift':
    result = ''.join(sorted(getAtomSiteAxisCode(x) for x in em.atomSites))
  else:
    result = tag + ''.join(sorted(isotope2Nucleus(x.isotopeCode).lower() for x in em.atomSites))

  #
  return result


def isotope2Nucleus(isotopeCode):
  """remove integer prefix from integer+string string"""
  ii = 0
  for ii,char in enumerate(isotopeCode):
    if char not in '0123456789':
      break
  return isotopeCode[ii:]


def getAtomSiteAxisCode(atomSite):
  """Get axisCode (without number suffixes) from NmrExPrototype.AtomSite"""

  name = atomSite.name
  result = nucleus = isotope2Nucleus(atomSite.isotopeCode)

  if name == 'CA':
    result = 'CA'

  elif name == 'CO':
    result = 'CO'

  elif nucleus in 'HF':
    # H or F, check if bound to C or N
    ll = []
    for et in atomSite.findAllExpTransfers(transferType='onebond'):
      for site in et.atomSites:
        if site is not atomSite:
          if site.isotopeCode == '13C':
            ss = 'c'
          elif site.isotopeCode == '15N':
            ss = 'n'
          else:
            break
          if site.maxNumber == 0:
              ll.append('-'+ss)
          elif site.minNumber > 0:
              ll.append(ss)

          break

    # In rare cases we get both - this is nonsense, so remove them
    for tags in ('n', '-n'),('c','-c'):
      if tags[0] in ll and tags[1] in ll:
        for tag in tags:
          while tag in ll:
            ll.remove(tag)

    result += ''.join(sorted(set(ll)))

  elif nucleus in 'CN':

    ll = []
    for et in atomSite.findAllExpTransfers(transferType='onebond'):
      for site in et.atomSites:
        if site is not atomSite:
          if site.isotopeCode == '1H':
            ss = 'h'
          elif site.isotopeCode == '19F':
            ss = 'f'
          else:
            break
          if site.maxNumber == 0:
              ll.append('-'+ss)
          elif site.minNumber > 0:
              ll.append(ss)
          break

    # In rare cases we get both - this is nonsense, so remove them
    for tags in ('h', '-h'),('f','-f'):
      if tags[0] in ll and tags[1] in ll:
        for tag in tags:
          while tag in ll:
            ll.remove(tag)

    result += ''.join(sorted(set(ll)))
  #
  return result


def testExpPrototypes():
  """Test functions and make diagnostic output"""
  from ccpncore.util.Io import newProject
  project = newProject("ExpPrototypeTest", removeExisting=True)

  codeMap = refExpDimRefCodeMap(project)

  axisCodeSet = set()
  useNames = {}
  synonyms= []
  for nmrExpPrototype in project.sortedNmrExpPrototypes():
    for refExperiment in nmrExpPrototype.sortedRefExperiments():

      # check axis codes
      codes = []
      for red in refExperiment.refExpDims:
        for redr in red.refExpDimRefs:
          code = codeMap[redr]
          codes.append(code)
          axisCodeSet.add(code)
      if len(codes) != len(set(codes)):
        print ("Duplicate code in %s: %s" % (refExperiment.name, codes))

      print ("TEST %s %s" % (refExperiment.name, codes))

      # check for duplicate synonyms and collect all synonyms
      key = refExperiment.name
      usename = refExperiment.synonym or key
      ll = useNames.get(usename, [])
      ll.append(key)
      useNames[usename] = ll
      if usename != key:
        synonyms.append((usename, key))

  print ("All axisCodes: %s" % list(sorted(axisCodeSet)))

  for key, val in sorted(useNames.items()):
    if len(val) > 1:
      print ("Duplicate name: %s %s" % (key, val))

  for tt in sorted(synonyms):
    print ("SYNONYM: %s %s" % tt)

if __name__ == '__main__':
  testExpPrototypes()