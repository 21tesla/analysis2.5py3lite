LICENSE = """
======================COPYRIGHT/LICENSE START==========================

ExperimentBasic.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
import array
import math
import os
import random, re

from memops.universal.Io import joinPath
from memops.universal.Util import isBigEndian, cumulativeProductArray
from memops.universal.BlockData import determineBlockSizes, writeBlockData, cumulativeArray, arrayOfIndex

from memops.general.Implementation import ApiError

from ccp.api.general.DataLocation import ShapeMatrix
from ccp.format.general import Util as ccpGenUtil
from ccp.format.bruker import Util as brukerUtil

from ccp.general.Constants import chemShiftRefRatios, shiftRatio_13C_TMS
from ccp.general.Io import getDataStoringFromFilepath

from ccp.util.NmrExpPrototype import longRangeTransfers
from ccp.util.Spectrum import createBlockedMatrix

from ccpnmr.analysis.core.Reference import getReference
from ccpnmr.analysis.core.UnitConverter import unit_converter

from ccp.lib.nmr.Nmr.DataSource import getIsotopeCodesList as getSpectrumIsotopes
from ccp.lib.nmr.Nmr.DataSource import getOnebondDataDims
from ccp.lib.nmr.Nmr.Experiment import getAcqExpDim, getOnebondExpDimRefs
from ccp.lib.nmr.Nmr.AbstractDataDim import getIsotopeCodes as getDataDimIsotopes



filterExpr1 = re.compile('\[[^A-Z]+\]')
filterExpr2 = re.compile('H\([0123]\)')

# TBD: put several functions from Util in here

COMMON_REF_EXPTS = set(['H', 'C', 'N', 'F', 'P','H[N]', 'H[C]', 'H[N[CO]]', 
                        'H[N[co[CA]]]', 'H[N[CO[CA]]]', 'H[N[ca[CO]]]',
                        'H[N[co[{CA|ca[C]}]]]','H[N[{CA|ca[Cali]}]]',
                        'H[N]_H.through-space','H_H[N].through-space', 'H[N[CA]]', 
                        'H[C]_H.through-space','H[N[ca[HA]]]', 'H{[N]+[HA]}',
                        'H[C]_H[C].through-space','H[N]_H[C].through-space',
                        'H[N]_H[N].through-space','H[N]_H.relayed','H_H[N].relayed',
                        'HH','H_H.through-space','H_H.relayed',
                        'HC_CH.relayed','HC_cH.relayed','Hc_CH.relayed',
                        'H_H.ROESY', 'CC',
                        ])

from operator import itemgetter
from heapq import nlargest
from itertools import repeat, ifilter

# Added for use in filtering RefgExperiments. Pre-Python-2.7 version of collections.Counter
class Counter(dict):
    '''Dict subclass for counting hashable objects.  Sometimes called a bag
    or multiset.  Elements are stored as dictionary keys and their counts
    are stored as dictionary values.

    >>> Counter('zyzygy')
    Counter({'y': 3, 'z': 2, 'g': 1})

    '''

    def __init__(self, iterable=None, **kwds):
        '''Create a new, empty Counter object.  And if given, count elements
        from an input iterable.  Or, initialize the count from another mapping
        of elements to their counts.

        >>> c = Counter()                           # a new, empty counter
        >>> c = Counter('gallahad')                 # a new counter from an iterable
        >>> c = Counter({'a': 4, 'b': 2})           # a new counter from a mapping
        >>> c = Counter(a=4, b=2)                   # a new counter from keyword args

        '''
        self.update(iterable, **kwds)

    def __missing__(self, key):
        return 0

    def most_common(self, n=None):
        '''List the n most common elements and their counts from the most
        common to the least.  If n is None, then list all element counts.

        >>> Counter('abracadabra').most_common(3)
        [('a', 5), ('r', 2), ('b', 2)]

        '''
        if n is None:
            return sorted(self.iteritems(), key=itemgetter(1), reverse=True)
        return nlargest(n, self.iteritems(), key=itemgetter(1))

    def elements(self):
        '''Iterator over elements repeating each as many times as its count.

        >>> c = Counter('ABCABC')
        >>> sorted(c.elements())
        ['A', 'A', 'B', 'B', 'C', 'C']

        If an element's count has been set to zero or is a negative number,
        elements() will ignore it.

        '''
        for elem, count in self.iteritems():
            for _ in repeat(None, count):
                yield elem

    # Override dict methods where the meaning changes for Counter objects.

    @classmethod
    def fromkeys(cls, iterable, v=None):
        raise NotImplementedError(
            'Counter.fromkeys() is undefined.  Use Counter(iterable) instead.')

    def update(self, iterable=None, **kwds):
        '''Like dict.update() but add counts instead of replacing them.

        Source can be an iterable, a dictionary, or another Counter instance.

        >>> c = Counter('which')
        >>> c.update('witch')           # add elements from another iterable
        >>> d = Counter('watch')
        >>> c.update(d)                 # add elements from another counter
        >>> c['h']                      # four 'h' in which, witch, and watch
        4

        '''
        if iterable is not None:
            if hasattr(iterable, 'iteritems'):
                if self:
                    self_get = self.get
                    for elem, count in iterable.iteritems():
                        self[elem] = self_get(elem, 0) + count
                else:
                    dict.update(self, iterable) # fast path when counter is empty
            else:
                self_get = self.get
                for elem in iterable:
                    self[elem] = self_get(elem, 0) + 1
        if kwds:
            self.update(kwds)

    def copy(self):
        'Like dict.copy() but returns a Counter instance instead of a dict.'
        return Counter(self)

    def __delitem__(self, elem):
        'Like dict.__delitem__() but does not raise KeyError for missing values.'
        if elem in self:
            dict.__delitem__(self, elem)

    def __repr__(self):
        if not self:
            return '%s()' % self.__class__.__name__
        items = ', '.join(map('%r: %r'.__mod__, self.most_common()))
        return '%s({%s})' % (self.__class__.__name__, items)

    # Multiset-style mathematical operations discussed in:
    #       Knuth TAOCP Volume II section 4.6.3 exercise 19
    #       and at http://en.wikipedia.org/wiki/Multiset
    #
    # Outputs guaranteed to only include positive counts.
    #
    # To strip negative and zero counts, add-in an empty counter:
    #       c += Counter()

    def __add__(self, other):
        '''Add counts from two counters.

        >>> Counter('abbb') + Counter('bcc')
        Counter({'b': 4, 'c': 2, 'a': 1})


        '''
        if not isinstance(other, Counter):
            return NotImplemented
        result = Counter()
        for elem in set(self) | set(other):
            newcount = self[elem] + other[elem]
            if newcount > 0:
                result[elem] = newcount
        return result

    def __sub__(self, other):
        ''' Subtract count, but keep only results with positive counts.

        >>> Counter('abbbc') - Counter('bccd')
        Counter({'b': 2, 'a': 1})

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        result = Counter()
        for elem in set(self) | set(other):
            newcount = self[elem] - other[elem]
            if newcount > 0:
                result[elem] = newcount
        return result

    def __or__(self, other):
        '''Union is the maximum of value in either of the input counters.

        >>> Counter('abbb') | Counter('bcc')
        Counter({'b': 3, 'c': 2, 'a': 1})

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        _max = max
        result = Counter()
        for elem in set(self) | set(other):
            newcount = _max(self[elem], other[elem])
            if newcount > 0:
                result[elem] = newcount
        return result

    def __and__(self, other):
        ''' Intersection is the minimum of corresponding counts.

        >>> Counter('abbb') & Counter('bcc')
        Counter({'b': 1})

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        _min = min
        result = Counter()
        if len(self) < len(other):
            self, other = other, self
        for elem in ifilter(self.__contains__, other):
            newcount = _min(self[elem], other[elem])
            if newcount > 0:
                result[elem] = newcount
        return result


def selectPreferredRefExperiment(refExperiments, experiment=None):
  """
  Select preferred RefExperiment for experiment type from input list

  .. describe:: Input
  
  List of ccp.nmr.NmrExpPrototype.RefExperiment
  
  ccp.nmr.Nmr.Experiment   # Currently not used
  
  .. describe:: Output

  ccp.nmr.NmrExpPrototype.RefExperiment
  """

  if not refExperiments:
    return None
  
  # Narrow down to common ref experiments, if any
  commonExps = [x for x in refExperiments if x.name in COMMON_REF_EXPTS]
  if commonExps:
    refExperiments = commonExps
  
  # return single candidate
  if len(refExperiments) == 1:
    return refExperiments[0]
  
  # Narrow down to simplest experiments (fewest atoms in transfer pathway)
  dd = {}
  for refExp in refExperiments:
    nAtom = len(refExp.nmrExpPrototype.atomSites)
    ll = dd.get(nAtom)
    if ll:
      ll.append(refExp)
    else:
      dd[nAtom]= [refExp]
  nAtoms = sorted(dd.keys())
  refExperiments = dd[nAtoms[0]]
  
  # remove multiple refExperiments from same prototype
  foundPrototypes = {}
  for refExp in refExperiments:
    prototype = refExp.nmrExpPrototype
    lastSerial = foundPrototypes.get(prototype)
    newSerial = refExp.serial
    if lastSerial is None or lastSerial > newSerial:
      foundPrototypes[prototype] = newSerial
  #
  refExperiments = [prot.findFirstRefExperiment(serial=ser) 
                    for (prot,ser) in foundPrototypes.items()]
      
  # return single candidate
  if len(refExperiments) == 1:
    return refExperiments[0]
  
  else:
    # do not guess
    return None
  
  # select alphabetically first name, for predictability
  #dd = {}
  #for refExp in refExperiments:
  #  dd[refExp.name] = refExp
  #for name,refExp in sorted(dd.items()):
  #  break
  #
  #return refExp
  

def getRefExperimentCategories(refExperiment):
  """
  Get a set of categories (as strings) that a reference experiment
  belongs to

  .. describe:: Input
  
  ccp.nmr.NmrExpPrototype.RefExperiment
  
  .. describe:: Output

  Set of Lines
  """
  
  # Constant time?
  
  if hasattr(refExperiment, 'typeCategories'):
    return refExperiment.typeCategories
  
  cats = set()
  measuredSites = set()
  
  prototype = refExperiment.nmrExpPrototype
  shiftOnly = True
  for refExpDim in refExperiment.refExpDims:
    for refExpDimRef in refExpDim.refExpDimRefs:
      
      # Measurement type categorisation
      measurement = refExpDimRef.expMeasurement
      measurementType = measurement.measurementType
      if measurementType != 'Shift':
        shiftOnly = False
        if measurementType in ['T1', 'T2', 'T1rho', 'T1zz']:
          cats.add('Relaxation')
        elif measurementType == 'JCoupling':
          cats.add('J-resolved')
        else:
          cats.add(measurementType)
 
      # Atom site categorisation 1
      for atomSite in measurement.atomSites:
        siteName = atomSite.name
        measuredSites.add(atomSite)

        if siteName == 'CA':
          cats.add('CA Measured')

        elif siteName == 'CO':
          cats.add('CO Measured')

        elif siteName in ('HA','HB'):
          cats.add('HA/HB Measured')  
 
  name = refExperiment.name
  ndim = len(refExperiment.refExpDims)
  
  # category 'projected'
  projStr = '.%dD.{' % ndim
  if projStr in name:
    maxdim = ndim
    for rx in refExperiment.nmrExpPrototype.refExperiments:
      maxdim = max(maxdim, len(rx.refExpDims))
    if maxdim > ndim:
      cats.add('Projected')
  
  # categories by multiple paths
  expGraphs = refExperiment.nmrExpPrototype.sortedExpGraphs()
  if len(expGraphs) > 1:
    signs = set(x.peakSign for x in expGraphs)
    if len(signs) > 1:
      cats.add('Signed Peaks')
    elif not ('J-resolved' in cats):
      cats.add('Mixed Peak Type')

  
  # filtered
  if ((filterExpr1.search(name) or filterExpr2.search(name)) 
      and not ('Signed Peaks' in cats)):
    cats.add('Filtered')
  
  # transfer type categorisation
  isSimple = True
  for expGraph in prototype.expGraphs:
    for expTransfer in expGraph.expTransfers:
      transferType = expTransfer.transferType
      
      if transferType not in  ('onebond', 'Jcoupling'):
        cats.add(transferType)
        isSimple = False
  if isSimple and shiftOnly:
    if (len(prototype.expGraphs) == 1 and 'Filtered' not in cats
        and len(prototype.expMeasurements) == len(refExperiment.refExpDims)):
      cats.add(' Simple')
    else:
      cats.add('J-transfer')
        

  # Atom site categorisation 2
  for atomSite in prototype.atomSites:
    siteName = atomSite.name
    
    if siteName == 'Caro':
      cats.add('Aromatic')
    
    if siteName == 'Cmet':
      cats.add('Methyl')
    
    if siteName in ('Csugar','Cbase', 'Nbase'):
      cats.add('Nucleotide') 
    
    # non-measured CO/CA
    if atomSite not in measuredSites:
      if siteName == 'CA':
        cats.add('Through CA')

      elif siteName == 'CO':
        cats.add('Through CO') 
  
  if name in COMMON_REF_EXPTS:
    cats.add(None) # Add to top level too

  refExperiment.typeCategories = cats
  
  return cats

def getSeqAssignRefExperiments(project):
  """
  Get reference amide experiments which may be used for protein sequence
  assignment. Gives a set of amide experiment types and a subset of through
  CO specific experiment types.

  .. describe:: Input
  
  Implementation.MemopsRoot
  
  .. describe:: Output

  List of NmrExpPrototype.RefExperiments,  List of NmrExpPrototype.RefExperiments
  """

  refExps = set()
  refExpsCO = set()
  
  for prototype in project.nmrExpPrototypes:
    atomSitesN = prototype.findAllAtomSites(name='N')
    
    if not atomSitesN:
      continue
    
    if len(atomSitesN) == 2: # HNN TBC
      continue
    
    atomSitesH = prototype.findAllAtomSites(name='H')
    
    if not atomSitesH:
      continue
    
    atomSiteCO = prototype.findFirstAtomSite(name='CO')
    onebondSites = set()
    onebondSitesCO = set()
    for atomSiteN in atomSitesN:
      for atomSiteH in atomSitesH:
        sitesHN = frozenset([atomSiteH, atomSiteN])
        
        for graph in prototype.expGraphs:
          for transferType in longRangeTransfers:
            if graph.findFirstExpTransfer(transferType=transferType):
              break
          
          else:
            transfer = graph.findFirstExpTransfer(atomSites=sitesHN,
                                                  transferType='onebond')
 
            if transfer:
              onebondSites.add(sitesHN)
 
              if atomSiteCO:
                sitesNC =  frozenset([atomSiteCO, atomSiteN])
 
                if graph.findFirstExpTransfer(atomSites=sitesNC,
                                              transferType='onebond'):
                  onebondSitesCO.add(sitesHN)
    
    if not onebondSites:
      continue
      
    for refExperiment in prototype.refExperiments:
      
      refExpDims = refExperiment.refExpDims
      
      if len(refExpDims) != 3:
        continue
      
      atomSites = set()
      for refExpDim in refExpDims:
        refExpDimRef = refExpDim.sortedRefExpDimRefs()[0]
        measurement = refExpDimRef.expMeasurement
        
        if measurement.measurementType.lower() == 'shift':
          atomSites.update(measurement.atomSites)
      
      if len(atomSites) != 3:
        continue
      
      isotopes = [a.isotopeCode for a in atomSites]
      
      if isotopes.count('15N') == 2: # HNN TBC
        continue
      
      for pair in onebondSites:
        if len(atomSites & pair) == 2:
          refExps.add(refExperiment)
          break                         

      for pair in onebondSitesCO:
        if len(atomSites & pair) == 2:
          if 'N[coca]' not in refExperiment.name:
            refExpsCO.add(refExperiment)
            break                         
  
  return refExps, refExpsCO

  
def getDataDimRefFullRange(dataDimRef):
  """
  Get the full range of freq values for a data dimension reference
  taking into account spectral width and min/max unaliased freqs

  .. describe:: Input
  
  ccp.nmr.Nmr.DataDimRef
  
  .. describe:: Output

  2-List of Floats (min, max)
  """
  
  expDimRef = dataDimRef.expDimRef
  converter = unit_converter[('point', expDimRef.unit)]

  valRange = [converter(1,dataDimRef),
              converter(dataDimRef.dataDim.numPoints,dataDimRef)]
  valRange.sort()
  
  valueMin = expDimRef.minAliasedFreq # Could be 0.0
  if valueMin is None:
    valueMin = valRange[0]
  
  valueMax = expDimRef.maxAliasedFreq 
  if valueMax is None:
    valueMax = valRange[1]
  
  return [valueMin, valueMax]

def getExperimentSampledDim(experiment):
  """
  Get the first sampled data dimension for an experiment
  
  .. describe:: Input
  
  ccp.nmr.Nmr.Experiment
  
  .. describe:: Output

  List of ccp.nmr.Nmr.SampledDataDim
  """
  
  sampledDim = None
  for expDim in experiment.expDims:
    for dataDim in expDim.dataDims:
      if dataDim.className == 'SampledDataDim':
        sampledDim = dataDim
        break
    else:
      continue
    break  
  
  return sampledDim
  

def getSampledDimExperiments(nmrProject):
  """
  Get a list of experiments in an NMR project that have sammpled data dims
  
  .. describe:: Input
  
  ccp.nmr.Nmr.NmrProject
  
  .. describe:: Output

  List of ccp.nmr.Nmr.Experiments
  """
  experiments = {}

  for experiment in nmrProject.experiments:
    for spectrum in experiment.dataSources:
      for dataDim in spectrum.dataDims:
        if dataDim.className == 'SampledDataDim':
          experiments[experiment] = True
          break
  
  return experiments.keys()

def getOnebondExpDimRefs(experiment):
  """
  Get pairs of experiment dimensions that are connected by onebond transfers

  .. describe:: Input

  Nmr.Experiment

  .. describe:: Output

  List of 2-List of Nmr.ExpDimRefs
  """

  expDimRefs   = []
  expTransfers = []

  for expTransfer in experiment.sortedExpTransfers():
    if expTransfer.transferType in ('onebond',):
      expTransfers.append(expTransfer)

  for expTransfer in expTransfers:
    expDimRefs.append(expTransfer.sortedExpDimRefs())

  return expDimRefs

def getEquivalentDataDims(spectrum):
  """
  Get pairs of spectrum data dimensions that are equivalent in terms of
  isotope type
  
  .. describe:: Input
  
  Nmr.DataSource
  
  .. describe:: Output

  List of 2-List of Nmr.DataDims
  """
  
  dataDims = [dd for dd in spectrum.sortedDataDims() if dd.className == 'FreqDataDim']
  dimPairs = []
  
  N = len(dataDims)
  for i in range(N-1):
    dataDimA = dataDims[i]
    isotopeA = getDataDimIsotopes(dataDimA)
  
    for j in range(i+1,N):
      dataDimB = dataDims[j]
      isotopeB = getDataDimIsotopes(dataDimB)
   
      if isotopeB and (isotopeB == isotopeA):
        dimPairs.append([dataDimA, dataDimB])
        
  return dimPairs 
  
def getOnebondDataDims(spectrum):
  """
  Get pairs of spectrum data dimensions that are connected by onebond transfers

  .. describe:: Input

  Nmr.DataSource

  .. describe:: Output

  List of 2-List of Nmr.DataDims
  """

  dataDims = []
  expDimRefs = getOnebondExpDimRefs(spectrum.experiment)

  for expDimRef0, expDimRef1 in expDimRefs:
    dataDim0 = spectrum.findFirstDataDim(expDim=expDimRef0.expDim)
    dataDim1 = spectrum.findFirstDataDim(expDim=expDimRef1.expDim)

    if dataDim0 and dataDim1:
      dataDims.append( [dataDim0,dataDim1] )

  return dataDims

def getDiagonalBondedDataDims(spectrum):
  """
  Get pairs of spectrum data dimensions that have digaonal positions
  and are also in one bond connections.  Returns list of (dim1, dim2)
  pairs (in 3D expect at most one), where dim1 is bonded to another
  dimension, and where dim1 and dim2 have the same isotope.
  
  .. describe:: Input
  
  Nmr.DataSource
  
  .. describe:: Output

  List of 2-List of Nmr.DataDims
  """

  oneBondDimPairs = getOnebondDataDims(spectrum)
  equivDimPairs = getEquivalentDataDims(spectrum)

  dataDimPairs = []
  for dataDimPair in oneBondDimPairs:
    dataDimSet = set(dataDimPair)
    for dataDim1, dataDim2 in equivDimPairs:
      if dataDim1 in dataDimSet:
        dataDimPairs.append((dataDim1, dataDim2))
      if dataDim2 in dataDimSet:
        dataDimPairs.append((dataDim2, dataDim1))
  
  return dataDimPairs

def initExpBoundResonances(experiment):
  """
  Refresh the covalently bound status for any resonances connected
  via peaks in a given experiment.
  
  .. describe:: Input
  
  Nmr.Experiment
  
  .. describe:: Output

  None
  """
  
  from ccpnmr.analysis.core.AssignmentBasic import getBoundResonances

  resonances = {}
  for spectrum in experiment.dataSources:
    for peakList in spectrum.peakLists:
      for peak in peakList.peaks:
        for peakDim in peak.peakDims:
          for contrib in peakDim.peakDimContribs:
            resonances[contrib.resonance] = None
            
  for resonance in resonances.keys():  
    getBoundResonances(resonance, recalculate=True)


def initExpTransfers(experiment, overwrite=True):
  """
  Set up the ExpTransfers for an experiment using available refExperiment
  information. Boolean option to remove any existing transfers.
  
  .. describe:: Input
  
  Nmr.Experiment
  
  .. describe:: Output
  
  List of Nmr.ExpTransfers
  """

  if not experiment.refExperiment:
    for expTransfer in experiment.expTransfers:
      expTransfer.delete()

    initExpBoundResonances(experiment)
    return

  if experiment.expTransfers:
    if not overwrite:
      return list(experiment.expTransfers)
    else:
      for expTransfer in experiment.expTransfers:
        expTransfer.delete()
  
  visibleSites = {}
  for expDim in experiment.expDims:
    for expDimRef in expDim.expDimRefs:
      if not expDimRef.refExpDimRef:
        continue
      
      measurement = expDimRef.refExpDimRef.expMeasurement
      if measurement.measurementType in ('Shift','shift','MQShift'):
        for atomSite in measurement.atomSites:
           if atomSite not in visibleSites:
             visibleSites[atomSite] = []
           
           visibleSites[atomSite].append(expDimRef)       

  
  transferDict = {}
  for atomSite in visibleSites:
    expDimRefs = visibleSites[atomSite]
    
    for expTransfer in atomSite.expTransfers:
      atomSiteA, atomSiteB = expTransfer.atomSites
      
      if (atomSiteA in visibleSites) and (atomSiteB in visibleSites):
        if transferDict.get(expTransfer) is None:
          transferDict[expTransfer] = []
        transferDict[expTransfer].extend(expDimRefs)
   
  # Indirect transfers, e.g. Ch_hC.NOESY or H_hC.NOESY
  indirectTransfers = set()
  for expGraph in experiment.refExperiment.nmrExpPrototype.expGraphs:
    for expTransfer in expGraph.expTransfers:
      if (expTransfer not in transferDict) and \
         (expTransfer.transferType in longRangeTransfers):
        atomSiteA, atomSiteB = expTransfer.atomSites
        
        if atomSiteA not in visibleSites:
          for expTransferA in atomSiteA.expTransfers:
            if expTransferA.transferType != 'onebond':
              continue
              
            atomSites = list(expTransferA.atomSites)
            atomSites.remove(atomSiteA)
            atomSiteC = atomSites[0]
            
            if atomSiteC in visibleSites:
              atomSiteA = atomSiteC
              break
        
        if atomSiteB not in visibleSites:
          for expTransferB in atomSiteB.expTransfers:
            if expTransferB.transferType != 'onebond':
              continue
              
            atomSites = list(expTransferB.atomSites)
            atomSites.remove(atomSiteB)
            atomSiteD = atomSites[0]
            
            if atomSiteD in visibleSites:
              atomSiteB = atomSiteD
              break
        
        if (atomSiteA in visibleSites) and (atomSiteB in visibleSites):
          expDimRefsA = visibleSites[atomSiteA]
          expDimRefsB = visibleSites[atomSiteB]
          transferDict[expTransfer] = expDimRefsA + expDimRefsB
          indirectTransfers.add(expTransfer)
        
  expTransfers = []
  for refTransfer in transferDict.keys():
    expDimRefs = frozenset(transferDict[refTransfer])
    if len(expDimRefs) == 2:
      transferType = refTransfer.transferType
      expTransfer  = experiment.findFirstExpTransfer(expDimRefs=expDimRefs)
      
      if expTransfer:
        # normally this would not need setting
        # but we renamed NOESY to through-space so this catches that situation
        expTransfer.transferType = transferType
      else:
        expTransfer = experiment.newExpTransfer(transferType=transferType,
                                                expDimRefs=expDimRefs)        
      
      if refTransfer in indirectTransfers:
        isDirect = False
      else:
        isDirect = True
      
      expTransfer.isDirect = isDirect
      expTransfers.append(expTransfer)
  
  initExpBoundResonances(experiment)
  return expTransfers      

def getIndirectThroughSpaceIsotopes(experiment):
  """
  For a given experiment find the pairs of isotopes present along
  onbond transfers connected via a relayed (indirect) through space transfer.
  Returns pairs of isotopes, direct and indirect, for each experimental
  dimension. An isotope may be None if one side of the through space transfer is
  observed in experiment, e.g. H_hC.NOESY.
  
  .. describe:: Input
  
  Nmr.Experiment
  
  .. describe:: Output
  
  Dict of Nmr.ExpDimRef:(ChemElement.Isotope, ChemElement.Isotope)
  """  
  
  isotopesDict = {}
  
  expTransfer0 = experiment.findFirstExpTransfer(isDirect=False)
  refExperiment = experiment.refExperiment
  
  if expTransfer0 and refExperiment:
    expDimRefA, expDimRefB = expTransfer0.expDimRefs
    refExpDimRefA = expDimRefA.refExpDimRef
    refExpDimRefB = expDimRefB.refExpDimRef
    isotopesDict[expDimRefA] = (expDimRefA.findFirstIsotope(), None)
    isotopesDict[expDimRefB] = (expDimRefB.findFirstIsotope(), None)
 
    if refExpDimRefA and refExpDimRefB:
      
     
      # Get directly shift measured atom sites
      visibleSites = set()
      for expDim in experiment.expDims:
        for expDimRef in expDim.expDimRefs:
          if not expDimRef.refExpDimRef:
            continue
 
          measurement = expDimRef.refExpDimRef.expMeasurement
          if measurement.measurementType in ('Shift','shift','MQShift'):
            visibleSites.update(measurement.atomSites)
      
 
      # Get atom sites at or onebond from expDimRefA
      measurementA = refExpDimRefA.expMeasurement
      atomSitesA = set(measurementA.atomSites)
      for atomSiteA in measurementA.atomSites:
        for expTransferA in atomSiteA.expTransfers:
          if expTransferA.transferType == 'onebond':
            atomSitesA.update(expTransferA.atomSites)

      # Get atom sites at or onebond from expDimRefB
      measurementB = refExpDimRefB.expMeasurement
      atomSitesB = set(measurementB.atomSites)
      for atomSiteB in measurementB.atomSites:
        for expTransferB in atomSiteB.expTransfers:
          if expTransferB.transferType == 'onebond':
            atomSitesB.update(expTransferB.atomSites)
 
      # Get long range atomSite pairs
      longRange = set()
      for expGraph in refExperiment.nmrExpPrototype.expGraphs:
        for expTransfer in expGraph.expTransfers:
          if expTransfer.transferType in longRangeTransfers:
            longRange.add(expTransfer.atomSites)
 
      # Get atom site pair where at least one is not direct
      # and pair is long range
      for atomSiteA in atomSitesA:
        for atomSiteB in atomSitesB:
          if (atomSiteA not in visibleSites) or (atomSiteB not in visibleSites):
            if frozenset([atomSiteA, atomSiteB]) in longRange:
              if atomSiteA not in visibleSites:
                isotope1 = expDimRefA.findFirstIsotope()
                isotope2 = atomSiteA.isotope
                isotopesDict[expDimRefA] = (isotope1, isotope2)
 
              if atomSiteB not in visibleSites:
                isotope1 = expDimRefB.findFirstIsotope()
                isotope2 = atomSiteB.isotope
                isotopesDict[expDimRefB] = (isotope1, isotope2)

              break
        else:
          continue
        break
 
  return isotopesDict


def getThroughSpaceDataDims(dataSource):
  """
  Get the data dims of a spectrum that represent through-space
  magnetisation transfers.
  
  .. describe:: Input
  
  Nmr.DataSource
  
  .. describe:: Output
  
  2-List of Nmr.FreqDataDims
  """

  expDims = set()
  
  for expTransfer in dataSource.experiment.expTransfers:
     if expTransfer.transferType in longRangeTransfers:
       for expDimRef in expTransfer.expDimRefs:
         expDims.add(expDimRef.expDim)

  dataDims = dataSource.sortedDataDims()
  dataDims = [dd for dd in dataDims if dd.expDim in expDims]

  return dataDims

def getIndirectDataDims(dataSource):
  """
  Get the data dims of a spectrum that represent indirect 
  (e.g. relayed NOESY) magnetisation transfers.
  
  .. describe:: Input
  
  Nmr.DataSource
  
  .. describe:: Output
  
  Set of 2-List of Nmr.FreqDataDims
  """

  expDims = set()
  
  for expTransfer in dataSource.experiment.expTransfers:
     if expTransfer.transferType in longRangeTransfers:
       for expDimRef in expTransfer.expDimRefs:
         expDims.add(expDimRef.expDim)

  dataDims = dataSource.sortedDataDims()
  dataDims = [dd for dd in dataDims if dd.expDim in expDims]

  # Check for indirect transfers
  indirectDims = set()
  dataDimDict = {}
  for dataDim in dataSource.dataDims:
    for dataDimRef in dataDim.dataDimRefs:
      dataDimDict[dataDimRef.expDimRef] = dataDim

  expDimRefs = set(dataDimDict.keys())
  for expTransfer in dataSource.experiment.expTransfers:
    if expTransfer.expDimRefs.issubset(expDimRefs) and not expTransfer.isDirect:
      dataDims = [dataDimDict[edr] for edr in expTransfer.expDimRefs]
      
      indirectDims.add(tuple(dataDims))

  return indirectDims


def getThroughSpacePeakLists(project, excludeSimulated=True):
  """
  Get the peak lists that have a through-space magnetisation
  ransfer, from a given project.
  
  .. describe:: Input
  
  Implementation.MemopsRoot, Boolean
  
  .. describe:: Output
  
  List of Nmr.PeakList
  """

  experiments = []
  for experiment in project.currentNmrProject.sortedExperiments():
    if experiment.refExperiment:
      name = experiment.refExperiment.nmrExpPrototype.name
      
      for transferType in longRangeTransfers:
        if transferType in name:
          experiments.append(experiment)
          break
      else:
        for expTransfer in experiment.expTransfers:
          if expTransfer.transferType in longRangeTransfers:
            experiments.append(experiment)
            break
  
  peakLists = []
  for experiment in experiments:
    for spectrum in experiment.sortedDataSources():
      if (spectrum.dataType == 'processed') and (spectrum.numDim > 1):
        for peakList in spectrum.sortedPeakLists():
          if excludeSimulated and peakList.isSimulated:
            continue
          peakLists.append(peakList)

  return peakLists

def getNoesyPeakLists(project, guessUntyped=True):
  """
  Get the NOE peak lists from a project.
  Option to include possible peak lists that are untyped
  
  .. describe:: Input
  
  Implementation.MemopsRoot, Boolean
  
  .. describe:: Output

  List of Nmr.PeakList
  """

  experiments = []
  for experiment in project.currentNmrProject.experiments:
    if experiment.refExperiment:
      name = experiment.refExperiment.nmrExpPrototype.name
      
      if 'NOESY' in name:
        experiments.append(experiment)
      else:
        for expTransfer in experiment.expTransfers:
          if expTransfer.transferType == 'through-space':
            experiments.append(experiment)
            break
  
  peakLists = []
  for experiment in experiments:
    for spectrum in experiment.dataSources:
      if (spectrum.dataType == 'processed') and (spectrum.numDim > 1):
        for peakList in spectrum.peakLists:
          peakLists.append(peakList)

  return peakLists

def getHsqcExperiments(nmrProject):
  """
  Get the (2D) HSQC experiments for an nmrProject
  
  .. describe:: Input
  
  Nmr.NmrProject
  
  .. describe:: Output

  List of Nmr.Experiment
  """

  experiments = []

  for expt in nmrProject.sortedExperiments():
    if expt.numDim == 2 and expt.refExperiment and expt.refExperiment.nmrExpPrototype.name in ('H[N]', 'H[C]'):
      experiments.append(expt)

  return experiments
  
def getCompatibleHsqcExperiments(experiment):
  """
  Get the (2D) HSQC experiments compatible with a given experiment
  in the sense that they both have a onebond connection with same
  isotopes at both ends, and also that it has some shared molSystems
  if that is set for both experiments
  
  .. describe:: Input
  
  Nmr.Experiment
  
  .. describe:: Output

  List of (Nmr.Experiment, ExpDimRef mapping)
  """

  experiments = []
  expDimRefPairs = getOnebondExpDimRefs(experiment)

  if not expDimRefPairs:
    return experiments
    
  molSystems = set(experiment.molSystems)
  expts = getHsqcExperiments(experiment.nmrProject)
  for experment2 in expts:
    if experment2 is not experiment:
      #expDimRefMapping = {}
      molSys = set(experment2.molSystems)
      if molSys and molSystems and not molSystems.intersection(molSys):
        continue
      expDimRefPairs2 = getOnebondExpDimRefs(experment2)
      for (expDimRefA, expDimRefB) in expDimRefPairs2:
        for (expDimRefC, expDimRefD) in expDimRefPairs:
          if (expDimRefA.isotopeCodes == expDimRefC.isotopeCodes) and (expDimRefB.isotopeCodes == expDimRefD.isotopeCodes):
            #expDimRefMapping[expDimRefA] = expDimRefC
            #expDimRefMapping[expDimRefB] = expDimRefD
            experiments.append(experment2)
            break
          elif (expDimRefA.isotopeCodes == expDimRefD.isotopeCodes) and (expDimRefB.isotopeCodes == expDimRefC.isotopeCodes):
            #expDimRefMapping[expDimRefA] = expDimRefD
            #expDimRefMapping[expDimRefB] = expDimRefC
            experiments.append(experment2)
            break
        else:
          continue
        break
      #if expDimRefMapping:
        #experiments.append((experment2, expDimRefMapping))
    
  return experiments

def getCompatibleHsqcPeakLists(peakList, excludeSimulated=True):
  """
  Get the (2D) HSQC peakLists compatible with a given peakList
  in the sense that the peakList experiment have an HSQC transfer
  in it and also that the respective experiments having some shared
  molSystems if that is set for both experiments.
  
  .. describe:: Input
  
  Nmr.PeakList
  
  .. describe:: Output

  List of Nmr.PeakList
  """

  peakLists = []
  experiments = getCompatibleHsqcExperiments(peakList.dataSource.experiment)
  #for (experiment, expDimRefMapping) in experiments:
  for experiment in experiments:
    for dataSource in experiment.sortedDataSources():
      pkLists = dataSource.sortedPeakLists()
      if excludeSimulated:
        pkLists = [pkList for pkList in pkLists if not pkList.isSimulated]
      peakLists.extend(pkLists)

  return peakLists

'''
def getAcqExpDim(experiment):
  """
  ExpDim that corresponds to acquisition dimension. NB uses heuristics

  .. describe:: Input

   NmrExpPrototype.Experiment

  .. describe:: Output

  NmrExpPrototype.ExpDim
  """

  ll = experiment.findAllExpDims(isAcquisition=True)
  if len(ll) == 1:
    # acquisition dimension set - return it
    result = ll.pop()

  else:
    # no reliable acquisition dimension set
    result = None

    dataSources = experiment.sortedDataSources()
    if dataSources:
      dataSource = dataSources[0]
      for ds in dataSources[1:]:
        # more than one data source. Pick one of the largest.
        if ds.numDim > dataSource.numDim:
          dataSource = ds

      # Take dimension with most points
      useDim = None
      currentVal = -1
      for dd in dataSource.sortedDataDims():
        if hasattr(dd, 'numPointsOrig'):
          val = dd.numPointsOrig
        else:
          val = dd.numPoints
        if val > currentVal:
          currentVal = val
          useDim = dd

      if useDim is not None:
        result = useDim.expDim

    if result is None:
      # no joy so far - just take first ExpDim
      ll = experiment.sortedExpDims()
      if ll:
        result = ll[0]

  #
  return result
'''

def getAcqRefExpDimRef(refExperiment):
  """
  RefExpDimRef that corresponds to acquisition dimension

  .. describe:: Input
  
   NmrExpPrototype.RefExperiment
  
  .. describe:: Output
  
  NmrExpPrototype.RefExpDimRef
  """
  
  # get acquisition measurement
  expGraph = refExperiment.nmrExpPrototype.findFirstExpGraph()
  # even if there are several the acquisition dimension should be common.
  ll = [(expStep.stepNumber, expStep) for expStep in expGraph.expSteps]
  expSteps = [x[1] for x in sorted(ll)]
  if refExperiment.isReversed:
    acqMeasurement = expSteps[0].expMeasurement
  else:
    acqMeasurement = expSteps[-1].expMeasurement
  
  # get RefExpDimRef that fits measurement
  ll = []
  for refExpDim in refExperiment.sortedRefExpDims():
    for refExpDimRef in refExpDim.sortedRefExpDimRefs():
      if refExpDimRef.expMeasurement is acqMeasurement:
        ll.append(refExpDimRef)
  
  if len(ll) == 1:
    return ll[0]
  else:
    raise ApiError("%s has no unambiguous RefExpDimRef for acqMeasurement (%s)"
                   % (refExperiment, acqMeasurement))
    
  
  
def setRefExperiment(experiment, refExperiment):
  """
  Sets the reference experiment for an existing experiment
  and tries to map the ExpDims to RefExpDims appropriately.
  
  .. describe:: Input
  
  Nmr.Experiment, NmrExpPrototype.RefExperiment
  
  .. describe:: Output

  None
  """
  
  
  for expDim in experiment.expDims:
    if expDim.refExpDim:
      for expDimRef in expDim.expDimRefs:
        if expDimRef.refExpDimRef:
          expDimRef.setRefExpDimRef(None)
  
      expDim.setRefExpDim(None)

  experiment.setRefExperiment(refExperiment)
  if refExperiment is None:
    return

  refExpDims = refExperiment.sortedRefExpDims()
  if not refExpDims:
    # Something is wrong with the reference data
    return
  
  expDims = experiment.sortedExpDims()
  if not expDims:
    # Something is wrong with the experiment
    return
  
  acqRefExpDim = getAcqRefExpDimRef(refExperiment).refExpDim
  acqExpDim = getAcqExpDim(experiment, ignorePreset=True)
  
  if ((refExpDims.index(acqRefExpDim)*2  < len(refExpDims)) !=
      (expDims.index(acqExpDim)*2  < len(expDims))):
    # acqRefExpDim and acqExpDim are at opposite ends of their
    # respective lists. reverse refExpDims so that acquisition
    # dimensions will more likely get mapped to each other.
    refExpDims.reverse()
  
  # Rasmus 12/7/12. Must be set to None, as otherwise it is never reset below.
  # We do not want the heuristic getAcqExpDim to override everything else.
  acqExpDim = None
  
  for expDim in expDims:
    expData = []
    
    for expDimRef in expDim.expDimRefs:
      isotopes = frozenset(expDimRef.isotopeCodes)
      
      if isotopes:
        mType = expDimRef.measurementType.lower()
        expData.append((mType, isotopes))
    
    if not expData:
      continue
    
    for refExpDim in refExpDims:
      refData = [] 
      
      for refExpDimRef in refExpDim.refExpDimRefs:
        expMeasurement = refExpDimRef.expMeasurement
        isotopes = frozenset([x.isotopeCode for x in expMeasurement.atomSites])
        mType = expMeasurement.measurementType.lower()
        refData.append((mType, isotopes))
 
      if expData == refData:
        expDim.setRefExpDim(refExpDim)
        refExpDims.remove(refExpDim)
        
        if refExpDim is acqRefExpDim:
          if not acqExpDim:
            expDim.isAcquisition = True
            acqExpDim = expDim
        
        break
  
  for expDim in expDims:
    if not expDim.expDimRefs:
      continue
  
    if not expDim.refExpDim:
      expDim.setRefExpDim(refExpDims.pop(0))
    
    
    # set reference data comparison list
    refExpDimRefs = list(expDim.refExpDim.refExpDimRefs)
    refData = []
    for refExpDimRef in refExpDimRefs:
      expMeasurement = refExpDimRef.expMeasurement
      atomSites = expMeasurement.atomSites
      refData.append((frozenset(x.isotopeCode for x in atomSites),
                     expMeasurement.measurementType.lower(),
                     frozenset(x.name for x in atomSites),
                     refExpDimRef))
    
    # set experiment data comparison list
    inData = []
    for expDimRef in expDim.expDimRefs:
      inData.append((frozenset(expDimRef.isotopeCodes),
                    expDimRef.measurementType.lower(),
                    frozenset(((expDimRef.displayName or expDimRef.name),)),
                    expDimRef))
    
    # match expDimRef to refExpDimRef. comparing isotopeCodes, 
    # if equal measurementTypes, if equal name/displayname
    for end in (-1,-2,-3):
      for ii in range(len(inData)-1, -1, -1):
        for jj in range(len(refData)-1, -1, -1):
          if inData[ii][:end] == refData[jj][:end]:
            expDimRef = inData[ii][-1]
            expDimRef.setRefExpDimRef(refData[jj][-1])
            expDimRef.measurementType = refData[jj][-1].expMeasurement.measurementType
            del inData[ii]
            del refData[jj]
            break
            
  
def getPossibleRefExperiments(experiment, category=None):
  """
  Get the possible reference NMR experiments for an experiment,
  given the prototypes available within its project,
  the category and external name settings.
  Optional argument to specify the category of possible experiments
  ('use external', 'through-bond', 'through-space', 'quantification',
   'other').
   
  .. describe:: Input
  
  Nmr.Experiment, NmrExpPrototype.NmrExpPrototype.ExpCategory
  
  .. describe:: Output

  Iterator of NmrExpPrototype.RefExperiments
  """
  
  memopsRoot = experiment.root
  
  if category is None and hasattr(experiment,'category'):
    category = experiment.category
  
  extSource = None
  if hasattr(experiment,'pulProgName'):
    extName = experiment.pulProgName
    if hasattr(experiment,'pulProgType'):
      extSource = experiment.pulProgType
  else:
    extName = None
  
  if category == 'use external':
    if extSource is None or not extName:
      # no external source - return all refExperiments
      return [x for y in memopsRoot.sortedNmrExpPrototypes()
              for x in y.sortedRefExperiments()]
      
    else:
      # get refExperiments from external source
      if extSource == 'ccpn':
        result = []
        for nxp in memopsRoot.sortedNmrExpPrototypes():
          if nxp.name == extName:
            # if name is prototype name assume any refExperiment might do
            result.extend(nxp.sortedRefExperiments())
          else:
            # add refExperiments
            result.extend(nxp.findAllRefExperiments(name=extName))
 
      elif extSource == 'bruker':
        result = []
        nameMatches = brukerUtil.parsePulProgName(extName)
        for name,isReversed in nameMatches:
          expPrototype = memopsRoot.findFirstNmrExpPrototype(name=name)
          if expPrototype is None:
            # name does not match entire prototype
            # - look for matching RefExperiments
            result.extend(x for y in memopsRoot.sortedNmrExpPrototypes()
                          for x in y.sortedRefExperiments()
                          if x.name == name)
          else:
            # match to expPrototype - take RefExperiments that fit isReversed
            if isReversed is None:
              result.extend(expPrototype.sortedRefExperiments())
            else:
              result.extend(expPrototype.findAllRefExperiments(
                                                isReversed=isReversed))
                                              
      else:
        raise NotImplementedError("External expt source %s not supported" 
                                % extSource)
      #
      return result
      
  elif category is None:
    # no category - return all RefExperiments
    return [x for y in memopsRoot.sortedNmrExpPrototypes()
            for x in y.sortedRefExperiments()]
  
  else:
    # return RefExperiments by category
    return [x for y in memopsRoot.findAllNmrExpPrototypes(category=category)
            for x in y.sortedRefExperiments()]
  

def getFilteredRefExperiments(experiment, category=None):
  """
  Get the possible reference NMR experiments for an experiment,
  given the prototypes available within its project,
  the category and external name settings and the nuclei on the axes.
  Optional argument to specify the category of possible experiments
  ('use external', 'through-bond', 'through-space', 'quantification',
   'other').
              
  .. describe:: Input
  
  Nmr.Experiment, NmrExpPrototype.NmrExpPrototype.ExpCategory
  
  .. describe:: Output

  List of NmrExpPrototype.RefExperiments
  """

  shiftMeasurements = ('Shift','shift','MQShift')

  # NOTE we only compare shift and MQ shift nuclei, and treat everything else as, in effect, 'other;

  # Get isotopes per experiment dimension (list of Counters)
  isotopes = []
  for expDim in experiment.expDims:
    dimIsotopes = Counter()
    isotopes.append(dimIsotopes)
  
    for expDimRef in expDim.expDimRefs:
      isotopes0 = expDimRef.isotopeCodes
      if isotopes0 and expDimRef.measurementType in shiftMeasurements:
        dimIsotopes[' '.join(sorted(isotopes0))] += 1
  
    if not dimIsotopes:
      # Set value for empty dimensions
      # We do it this way to identify the dimension for matching (empty would always match)
      dimIsotopes[None] += 1

  refExperiments = []

  numDims = experiment.numDim
  
  for refExperiment in getPossibleRefExperiments(experiment, category):

    if len(refExperiment.refExpDims) == numDims:

      # Get isotopes per RefExperiment dimension (list of Counters)
      refIsotopes = []
      for refExpDim in refExperiment.refExpDims:
        dimIsotopes =  Counter()
        refIsotopes.append(dimIsotopes)

        for refExpDimRef in refExpDim.refExpDimRefs:
          expMeasurement = refExpDimRef.expMeasurement
          if expMeasurement.measurementType in shiftMeasurements:
            # We only count shift and MQshift dimensions
            iCodes = [ass.isotopeCode for ass in expMeasurement.atomSites]
            if expMeasurement.measurementType == 'MQShift' and len(set(iCodes)) == 1:
              # RASMUS HACK. Correctly experiments ought to have all relevant nuclei,
              # but in practice e.g. '13C','13C' DQ come in with isotopeCode just '13C'
              # So to catch these (without requiring users to set dimensions correctly)
              # we do this
              dimIsotopes[iCodes[0]] += 1
            else:
              dimIsotopes[' '.join(sorted(iCodes))] += 1

        if not dimIsotopes:
          # Set value for empty dimensions
          # We do it this way to identify the dimension for matching (empty would always match)
          dimIsotopes[None] += 1


      # remove matching elements
      for iso in isotopes:
        for riso in refIsotopes:
          if not list((iso-riso).elements()):
            # Checks if there is any element in iso that is not matched by one in riso
            refIsotopes.remove(riso)
            break

      # The lengths were the same. If all refIsotopes were removed, the two lists match
      if not refIsotopes:
        refExperiments.append(refExperiment)

  #for rx in refExperiments:
  #  print ('@~@~ matching:', rx.name)

  return refExperiments
   
def getRefExperiments(experiment):
  """
  Wrapper to getFilteredRefExperiments - set category from
  1) experiment.category (temporary) attribute (if present)
  2) experiment.refExperiment category (if present)
  and pass call on to getFilteredRefExperiments
  
  .. describe:: Input
  
  Nmr.Experiment
  
  .. describe:: Output

  List of NmrExpPrototype.RefExperiments
  """
  category = None
  if hasattr(experiment, 'category'):
    category = experiment.category
  
  if category is None and experiment.refExperiment:
    category = experiment.refExperiment.nmrExpPrototype.category
  
  return getFilteredRefExperiments(experiment, category=category)

def getSpectrumNoise(dataSource):
  """
  Get the noise level for a spectrum. If the noise level is not already set it will
  be set at an estimated value.
  
  .. describe:: Input
  
  Nmr.DataSource
  
  .. describe:: Output

  Float
  """

  noise = dataSource.noiseLevel
  if noise is None:
    noise = getNoiseEstimate(dataSource, nsamples=1000, nsubsets=100, fraction= 0.1)
    dataSource.noiseLevel = noise
  
  return noise

def getNoiseEstimate(dataSource, nsamples=1000, nsubsets=10, fraction=0.1):
  """
  Estimate the noise level for a spectrum by choosing a random nsamples points
  and finding subsets with the lowest standars devation/
  
  .. describe:: Input
  
  Nmr.DataSource, Int, Int, Float
  
  .. describe:: Output

  Float
  """

  sqrt = math.sqrt

  data = nsamples * [0]
 
  if not hasattr(dataSource, 'block_file') or not dataSource.block_file:
    #raise 'dataSource does not have attribute "block_file"'
    return 1.0  # arbitrary
 
  block_file = dataSource.block_file
  npts = [dataDim.numPoints for dataDim in dataSource.sortedDataDims()]
 
  fails = 0
  for i in range(nsamples):
    pt = [ min(n-2,int(n*random.random())) for n in npts ]
    
    try:
      d = block_file.getPointValue(pt)
      
      if d - d == 0: # Fails for inf and not-a-number
        data.append(d)
      else:
        fails += 1
        continue
      
    except:
      fails += 1
      continue
  
  if fails:
    msg = "Attempt to access %d non-existent data points in spectrum %s:%s"
    print msg % (fails, dataSource.experiment.name, dataSource.name)

  good = nsamples - fails
  if good == 0:
    return 1.0 # arbitrary
  elif good < 10:
    maxvalue = max([abs(x) for x in data])
    if maxvalue > 0:
      return 0.1 * maxvalue
    else:
      return 1.0 # arbitrary

  m = int(nsamples * fraction)
  minStd = None
  for i in range(nsubsets):
    y = getRandomSubset(data, fraction)
    n = len(data)
 
    avg = 0
    s2  = 0
    for x in data:
      avg += x
      s2  += x*x
 
    avg /= n
    s2 /= n
    s2 -= avg*avg
 
    std = sqrt(max(s2,0))
    
    if (minStd is None):
      minStd = std
    else:
      minStd = min(std, minStd)


  # multiplier a guess
  minStd = 1.1 * minStd

  if minStd == 0: # in case data file is blank
    minStd = 1.0

  minStd *= dataSource.scale

  return minStd

def calculateNoiseInBox(dataSource, boxMin, boxMax):

  try:
    import numpy
  except:
    return

  if not hasattr(dataSource, 'block_file') or not dataSource.block_file:
    return

  values = dataSource.block_file.getValues(boxMin, boxMax)
  if not values:
    return

  values = numpy.array(values)

  n = len(values)
  if n > 2:
    noise = dataSource.scale * numpy.std(values, ddof=1)

    dataSource.noiseLevel = noise
    print 'Set noise level for spectrum %s:%s to %f' % (dataSource.experiment.name, dataSource.name, noise)

def getMinMaxValues(dataSource):
  """
  Get the min and max values of a dataSource.
  Not recommended unless 1D or 2D because slow.
  
  .. describe:: Input
  
  Nmr.DataSource
  
  .. describe:: Output

  (Float, Float)
  """

  minValue = 0.0  # arbitrary
  maxValue = 1.0  # arbitrary

  if not hasattr(dataSource, 'block_file') or not dataSource.block_file:
    #raise 'dataSource does not have attribute "block_file"'
    return (minValue, maxValue)
 
  block_file = dataSource.block_file
  npts = [dataDim.numPoints for dataDim in dataSource.sortedDataDims()]
 
  (n, cumNpts) = cumulativeArray(npts)

  data = []
  fails = 0
  for i in range(n):
    pt = arrayOfIndex(i, cumNpts)
    
    try:
      d = block_file.getPointValue(pt)
      
      if d - d == 0: # Fails for inf and not-a-number
        data.append(d)
      else:
        fails += 1
        continue
      
    except:
      fails += 1
      continue
  
  if fails:
    msg = "Attempt to access %d non-existent data points in spectrum %s:%s"
    print msg % (fails, dataSource.experiment.name, dataSource.name)

  if data:
    minValue = min(data)
    maxValue = max(data)

  return (minValue, maxValue)

def getDataSlice(dataSource, position, sliceDim, unit='point', diagonalExclusion=0, normalise=True):
  """
  Get the 1D slice for the dataSource for a position which
  is fixed in all dimensions but sliceDim.
  Assumes that all dimenions are freqDataDim.
  position is of size dataSource.numDim.
  sliceDim is 1-based, not 0-based.
  The slice is returned as a Python list.
  
  .. describe:: Input
  
  Nmr.DataSource, tuple/list of Float, Int
  
  .. describe:: Output

  list of Float
  """

  from ccpnmr.analysis.core.Util import convertPosition

  if (not hasattr(dataSource, 'block_file') or not dataSource.block_file):
    raise Exception('dataSource does not have attribute "block_file"')
 
  block_file = dataSource.block_file

  ndim = dataSource.numDim
  box_min = ndim * [0]
  box_max = ndim * [0]
  
  for i in range(ndim):
    dim = i + 1
    dataDim = dataSource.findFirstDataDim(dim=dim)
    dataDimRef = getPrimaryDataDimRef(dataDim)
    
    if dim == sliceDim:
      box_min[i] = 0
      box_max[i] = dataDim.numPoints
      
      if diagonalExclusion:
        pd = int(position[i])
        if unit != 'point':
          pd = int(convertPosition(pd, dataDimRef, fromUnit=unit))
          
    else:
      p = position[i]
      if unit != 'point':
        p = convertPosition(p, dataDimRef, fromUnit=unit)
        
      p = p - 1  # points start counting from 1
      p = int(p % dataDim.numPointsOrig + dataDim.pointOffset)
      box_min[i] = p
      box_max[i] = p + 1

  values = block_file.getValues(box_min, box_max)

  if diagonalExclusion:
    nPoints = dataSource.findFirstDataDim(dim=sliceDim).numPoints
    
    for i in range(-diagonalExclusion,diagonalExclusion):
      j = pd+i
      
      if (j > 0) and (j < nPoints):
        values[j] = 0.0
        
  if normalise:
    noise  = getSpectrumNoise(dataSource)*2.0
    values = [max(noise,v)-noise for v in values]
    
    maxVal = float(max(values))
    values = [v/maxVal for v in values]

  return values
 
def calcContourLevels(baseLevel, numberLevels, levelChanger, changeMode='multiply'):
  """
  Calculate the contour levels given baseLevel, numberLevels, levelChanger, changeMode
  
  .. describe:: Input
  
  Nmr.DataSource, Float, Int, Float, Word (multiply or add)
  
  .. describe:: Output

  The contour levels as a list
  """

  levels = []
  for n in range(numberLevels):
    if (n == 0):
      level = baseLevel
    elif changeMode == 'multiply':
      level = levelChanger * levels[-1]
    else:
      level = levelChanger + levels[-1]
      
    levels.append(level)

  return levels

def changeSpectrumNContours(spectrum, delta):
  """
  Change the number of contour levels of a spectrum by a given delta
  
  .. describe:: Input
  
  Nmr.DataSource, Int
  
  .. describe:: Output

  None
  """

  from ccpnmr.analysis.core.Util import getAnalysisSpectrum

  analysisSpectrum = getAnalysisSpectrum(spectrum)

  levelChanger = analysisSpectrum.autoLevelChanger
  changeMode   = analysisSpectrum.autoLevelMode

  pos = list(analysisSpectrum.posLevels)
  neg = list(analysisSpectrum.negLevels)
  
  if pos:
    pos.sort()
    numPos = len(pos)+delta
    baseLevel = pos[0]
    
    if numPos > 0:
      pos = calcContourLevels(baseLevel, numPos, levelChanger, changeMode)
      analysisSpectrum.posLevels = pos
      analysisSpectrum.autoBaseLevel = baseLevel
      analysisSpectrum.autoNumLevels = numPos
  
  if neg:
    neg.sort()
    numNeg = len(neg)+delta
    baseLevel = neg[-1]
  
    if numNeg > 0:
      if changeMode == 'add':
        levelChanger = -levelChanger
        
      neg = calcContourLevels(baseLevel, numNeg, levelChanger, changeMode)
      analysisSpectrum.negLevels = neg
      
      if not pos:
        analysisSpectrum.autoBaseLevel = -baseLevel
        analysisSpectrum.autoNumLevels = numNeg      

  if not (pos or neg):
    numPos = delta
    
    if numPos > 0:
      baseLevel = analysisSpectrum.autoBaseLevel
      pos = calcContourLevels(baseLevel, numPos, levelChanger, changeMode)
      analysisSpectrum.posLevels = pos
      analysisSpectrum.autoNumLevels = numPos
      

def ppmDataDimBoundedRegion(region, dataDim):
  """
  Truncate a region for a data dimension at its ppm bounds.
  
  .. describe:: Input
  
  List of Floats (region), Nmr.DataDim
  
  .. describe:: Output

  List of Floats (region)
  """

  (r0, r1) = region

  if dataDim.className == "FreqDataDim":

    dataDimRef =  getPrimaryDataDimRef(dataDim)
    expDimRef = dataDimRef.expDimRef
    minAliasedFreq = expDimRef.minAliasedFreq
    maxAliasedFreq = expDimRef.maxAliasedFreq
    if (minAliasedFreq is not None):
      r0 = max(r0, minAliasedFreq)
    if (maxAliasedFreq is not None):
      r1 = min(r1, maxAliasedFreq)
    r1 = max(r0, r1)

  return (r0, r1)

def ptsDataDimBoundedRegion(region, dataDim):
  """
  Truncate a region for a data dimension at its data points bounds.
  
  .. describe:: Input
  
  List of Floats (region), Nmr.DataDim
  
  .. describe:: Output

  List of Floats (region)
  """

  (r0, r1) = region

  if dataDim.className == "FreqDataDim":

    dataDimRef =  getPrimaryDataDimRef(dataDim)
    expDimRef = dataDimRef.expDimRef
    minAliasedFreq = expDimRef.minAliasedFreq
    maxAliasedFreq = expDimRef.maxAliasedFreq
    if (maxAliasedFreq is None):
      r0 = max(r0, 0)
    if (minAliasedFreq is None):
      r1 = min(r1, dataDim.numPoints-1)
    r1 = max(r0, r1)

  return (r0, r1)


def getRandomSubset(data, fraction):
  """
  Function used by getNoiseEstimate()

  .. describe:: Input
  
  
  
  .. describe:: Output

  
  """

  n = len(data)
  m = int(n * fraction)
 
  y = []
  x = data[:]
  for i in range(m):
    # below guarantees no repeats
    k = min(n-i-1, int((n-i) * random.random()))
    # below allows repeats
    #k = min(n-1, int(n * random.random()))
    y.append(x[k])
    # below guarantees no repeats
    del x[k]
 
  return y

def getExperimentSpectra(experiment, minNumDim=1):
  """
  Gives the spectra for a given NMR experiment with a minimum number of dimensions 
  .. describe:: Input
  
  Nmr.Experiment, Integer
  
  .. describe:: Output

  List of Nmr.DataSources
  """

  if experiment and not experiment.isDeleted:
    dataSources = experiment.sortedDataSources()
    spectra = [ d for d in dataSources if not d.isDeleted and isSpectrum(d, minNumDim) ]
  else:
    spectra = []
 
  return spectra

def isSpectrum(dataSource, minNumDim=1):
  """
  Determines if a spectrum is valid for display, i.e. processed and has the min required
  number of dimensions
  
  .. describe:: Input
  
  Nmr.DataSource, Int
  
  .. describe:: Output

  Boolean
  """
 
  if (dataSource.numDim >= minNumDim) and  (dataSource.dataType == 'processed'):
    # TBD is this what we want?
    return True
    
  else:
    return False

def getSpectra(project, minNumDim=1):
  """
  Give all processed spectra in a project with a minimum number of dimensions
  
  .. describe:: Input
  
  Project, Int
  
  .. describe:: Output

  List of Nmr.DataSources
  """

  spectra = []

  for expt in project.currentNmrProject.sortedExperiments():
    spectra.extend(getExperimentSpectra(expt, minNumDim))

  return spectra

def newShiftList(project, unit='ppm'):
  """
  Make a new ShiftList for a project with an optionally specified unit (default ppm). 
  
  .. describe:: Input
  
  Project, String (Nmr.ShiftList.unit)
  
  .. describe:: Output

  Nmr.ShiftList
  """

  nmrProject = project.currentNmrProject
  shiftList = nmrProject.newShiftList(unit=unit)
  nameTemplate = 'ShiftList %d'
  
  i = 1
  name = nameTemplate % i
  while nmrProject.findFirstMeasurementList(className='ShiftList', name=name):
    i += 1
    name = nameTemplate % i
  
  shiftList.name = name
  
  return shiftList

def setExperimentShiftList(experiment, shiftList):
  """
  Set the shift list for a given experiment. Updates the shift values
  from the appropriate shift lists/
  
  .. describe:: Input
  
  Nmr.NmrExperiment, Nmr.ShiftList
  
  .. describe:: Output

  None
  """
  from ccpnmr.analysis.core.AssignmentBasic import updateResonShift, averageShiftValue

  if experiment.shiftList is shiftList:
    return

  prevShiftList = experiment.shiftList
  experiment.setShiftList(shiftList)
  
  resonanceDict = {}
  for dataSource in experiment.dataSources:
    if isSpectrum(dataSource):
      for peakList in dataSource.peakLists:
        for peak in peakList.peaks:
          for peakDim in peak.peakDims:
            for contrib in peakDim.peakDimContribs:
              resonance = contrib.resonance
              
              if resonanceDict.get(resonance.serial) is None:
                resonanceDict[resonance.serial] = True
                updateResonShift(resonance, peakDim)
                
                for shift in resonance.shifts:
                  if shift.parentList is prevShiftList:
                    averageShiftValue(shift) 
                     
  
def getSpectraByType(project, experimentType):
  """
  Give NMR spectra of a given class in a project as identified by isotope-dimension mapping
  This funtion is just a temporary kludge prior to proper implementation on spectra types   
  
  .. describe:: Input
  
  Project, String
  
  .. describe:: Output

  List of Nmr.DataSources
  """

  experimentTypes = [experimentType,]
  
  spectra = []
  experiments = project.currentNmrProject.experiments
  for experiment in experiments:
    if experiment.refExperiment:
      if experiment.refExperiment.name in experimentTypes:
        for dataSource in experiment.dataSources:
          if isSpectrum(dataSource):
            spectra.append(dataSource)
      
  return spectra
 
def getSpectrumIsotopes(dataSource):
  """
  Give isotopes pertianing to the dimensions of a given NMR spectrum 

  .. describe:: Input

  Nmr.DataSource

  .. describe:: Output

  List of Words (Nmr.ExpDimRef.IsotopeCodes)
  """

  isotopes = []
  for dataDim in dataSource.sortedDataDims():
    isotopeCodes = list(getDataDimIsotopes(dataDim))
    isotopes.append(','.join(isotopeCodes) or None)

  return isotopes

def findSpectrumDimsByIsotope(dataSource, isotopeCode):
  """
  Give dimension numbers of a spectrum that pertain to a given Isotope
  
  .. describe:: Input
  
  Nmr.DataSource, Word (Nmr.ExpDimRef.isotopeCode)
  
  .. describe:: Output

  List of Ints
  """

  dims = []
  dataDims = dataSource.sortedDataDims()
  N = len(dataDims)

  for i in range(N):
    isotopes = getDataDimIsotopes(dataDims[i])
    if isotopeCode in isotopes:
      dims.append(i)
    
  return dims 
  
def getDataDimIsotopes(dataDim):
  """
  Get the shift measurement isotopes for a spectrum data dim

  .. describe:: Input

  Nmr.FreqDataDim

  .. describe:: Output

  Set of Words (Nmr.ExpDimRef.isotopeCodes)
  """

  isotopes = set()

  for expDimRef in dataDim.expDim.expDimRefs:
    if expDimRef.measurementType in ('Shift','shift'):
      for isotopeCode in expDimRef.isotopeCodes:
        isotopes.add(isotopeCode)

  return isotopes

def getNdSpectra(project, n):
  """
  Give all NMR spectra in a given project that have n dimensions
  
  .. describe:: Input
  
  Project, Int
  
  .. describe:: Output

  List of Nmr.DataSources
  """

  spectra = []
  experiments = project.currentNmrProject.experiments
  
  for experiment in experiments:
    for dataSource in experiment.dataSources:
      if dataSource.dataType == 'processed':
        if dataSource.numDim == n:
          spectra.append(dataSource)
  
  return spectra


def getIsSpectrumBigEndian(spectrum):
  """
  Get whether spectrum's data file is big endian.
  
  .. describe:: Input
  
  Spectrum
  
  .. describe:: Output

  None if no dataStore or not BlockedBinaryMatrix, True if big endian, False if little endian
  """

  dataStore = spectrum.dataStore
  if not dataStore:
    return None

  if dataStore.className != 'BlockedBinaryMatrix':
    return None

  return dataStore.isBigEndian


def setIsSpectrumBigEndian(spectrum, isBigEndian):
  """
  Set whether spectrum's data file is big endian.
  
  .. describe:: Input
  
  Spectrum, Boolean
  
  .. describe:: Output

  None
  """

  if isShapeSpectrum(spectrum):
    return

  dataStore = spectrum.dataStore
  if not dataStore:
    return

  if dataStore.className != 'BlockedBinaryMatrix':
    return

  dataStore.isBigEndian = isBigEndian


def isSpectrumBigEndian(spectrum):
  """
  Determine (approximately) whether spectrum's data file is big endian according to actual data.
  
  .. describe:: Input
  
  Spectrum
  
  .. describe:: Output

  None if cannot determine, True if big endian, False if little endian
  """

  if isShapeSpectrum(spectrum):
    return None

  dataStore = spectrum.dataStore
  if not dataStore:
    return None

  fileName = dataStore.fullPath
  if not os.path.exists(fileName):
    return None

  numberType = dataStore.numberType
  fileHeaderSize = dataStore.headerSize
  if fileHeaderSize is None:
    fileHeaderSize = 0
    
  # TBD: below only works for BlockedBinaryMatrix
  nbytes = dataStore.nByte

  return isDataBigEndian(fileName, numberType=numberType, fileHeaderSize=fileHeaderSize, nbytes=nbytes)


def isDataBigEndian(fileName, numberType='float', fileHeaderSize=0,
                    nbytes=4, nwords=1000, threshold=1.0e9, fraction=0.01,
                    fractionNonzero=0.1):
  """
  Determine (approximately) whether data file is big endian.
  
  .. describe:: Input
  
  fileName and optionally numberType ('float' or 'int'), fileHeaderSize,
  nbytes (bytes per word), nwords (number of words to check),
  threshold (above which value considered incorrect byte ordering),
  fraction (fraction of values checked which must have incorrect byte ordering, if considered swapped)
  fractionNonzero (fraction of values checked which must be nonzero to be considered ok)
  
  .. describe:: Output

  True if big endian or if cannot determine, False if little endian
  """

  if nbytes != 4:  # TBD: not sure what else to do here right now
    print 'WARNING: isDataBigEndian() returning True for nbytes != 4 right now'
    return True

  try:
    fp = open(fileName, 'rb')
  except:
    print 'WARNING: file "%s" does not exist, isDataBigEndian() returning True' % fileName
    return True

  s = fp.read(fileHeaderSize)

  while True:
    if numberType == 'float':
      x = array.array('f')
    else:
      x = array.array('i')

    s = fp.read(nbytes*nwords)

    if not s:
      break

    x.fromstring(s)
    if not isBigEndian():
      x.byteswap()

    knt = 0
    for n in range(len(x)):
      if x[n] != 0:
        knt = knt + 1

    if knt > fractionNonzero * len(x):
      # enough values are not zero
      break

  fp.close()
  if not x:
    print 'WARNING: file "%s" seems to be mostly 0, isDataBigEndian() returning True' % fileName
    return True

  knt = 0
  nonzero = 0
  for n in range(len(x)):
    if x[n] != 0:
      nonzero = nonzero + 1
      if (abs(x[n]) > threshold):
        knt = knt + 1

  r = int(fraction * nonzero)
  if knt > r:
    return False
  else:
    return True

def isShapeSpectrum(spectrum):
  """
  Determine if spectrum is a Shape spectrum, possibly during set-up
  
  .. describe:: Input
  
  Spectrum
  
  .. describe:: Output

  Boolean
  """
  return isinstance(spectrum.dataStore, ShapeMatrix) or \
           (hasattr(spectrum, 'valuesList') and spectrum.valuesList)

def getPrimaryExpDimRef(expDim):
  """
  get expDimRef child with lowest expDimRef.serial
  
  .. describe:: Input
  
  expDim
  
  .. describe:: Output

  expDimRef or None
  """
  expDimRefs = expDim.sortedExpDimRefs()
  if expDimRefs:
    return expDimRefs[0]
  else:
    return None

def getPrimaryDataDimRef(freqDataDim):
  """
  get dataDimRef child with lowest expDimRef.serial
  
  .. describe:: Input
  
  freqDataDim
  
  .. describe:: Output

  dataDimRef or None
  """
  dataDimRefs = freqDataDim.dataDimRefs
  if dataDimRefs:
    ll = [(x.expDimRef.serial, x) for x in dataDimRefs]
    ll.sort()
    return ll[0][1]
  else:
    return None

def setAllIsotopeCodes(experiment):
  """
  set all ExpDimRef.isotopeCodes in experiment, from sf values
  Selects isotopes so that sf ratios match isotope gyromagmetic ratios
  Unknown and inconsistent isotopes efault to 1H
  
  .. describe:: Input
  
  ccp/nmr.Nmr.Experiment
  
  .. describe:: Output

  None
  """
  defaultCode = '1H'
  matchTol = 0.01
  refRatioData = [(value, key) for key,value in chemShiftRefRatios.items()]
  refRatioData.sort()
  refRatioData.reverse()
  refRatioData.append(refRatioData[0])
  
  data = []
  codesMissing = False
  for expDim in experiment.sortedExpDims():
    for expDimRef in expDim.sortedExpDimRefs():
      isotopeCodes = expDimRef.isotopeCodes
      sf = expDimRef.sf
      if isotopeCodes:
        if len(isotopeCodes) == 1:
          isotop = isotopeCodes[0]
        else:
          isotop = isotopeCodes
      else:
        if sf:
          codesMissing = True
          isotop = None
        else:
          # nothing doing without a nonzero sf value. Set to default.
          expDimRef.isotopeCodes = (defaultCode,)
          expDimRef.unit = 'ppm'
          continue
      
      if sf:
        # append data for further processing
        data.append((sf, isotop, expDimRef))
        
  
  if codesMissing:
    # try to fill in codes
    
    # sort to get largest sf first.
    data.sort()
    
    # find sf with isotope code and set ratio for calculation
    for sf, isotop, expDimRef in data:
      refRatio = chemShiftRefRatios.get(isotop)
      if refRatio:
        ratio = sf / refRatio
        break
    else:
      ratio = None
      
    if ratio:
      # use ratio to set isotopeCodes
      for sf, isotop, expDimRef in data:
        if not isotop:
          xx = sf / ratio
          for val, code in refRatioData:
            if abs(1.0 - xx/val) < matchTol:
              expDimRef.isotopeCodes = (code,)
              expDimRef.unit = 'ppm'
              break
          else:
            expDimRef.isotopeCodes = (defaultCode,)
            expDimRef.unit = 'ppm'
    
    else:
      # no usable starting codes - look for a set that matches sf.        
      
      # Note: loop tries to fit highest sf to decreasing order of refRatio,
      # starting with 1H.
      # if no refRatio matches uses last element, which 
      # sets highest sf to 1H and uses 1H for incompatible sf's.
      for refRatio, junk in refRatioData:
        ratio = sf / refRatio
        ok = True
        for sf, isotop, expDimRef in data:
          if not isotop:
            xx = sf / ratio
            for val, code in refRatioData:
              if abs(1.0 - xx/val) < matchTol:
                expDimRef.isotopeCodes = (code,)
                expDimRef.unit = 'ppm'
                break
            else:
              ok = False
              expDimRef.isotopeCodes = (defaultCode,)
              expDimRef.unit = 'ppm'
        #
        if ok:
          break
    
def isReferencingIncorrect(spectrum, fixErrors=False, printChecks=True):
  """ Check if referencing is consistent 
  and/or fitw IUPAC referencing ratios
  """
  
  sfTolerance = 1.0e-7  
  brukerTolerance = 5.0e-7  # NB this is the smallest value that will catch 
                            # the standard Bruker TMS-ref value 0.2514501008
  result = (False, '')
  
  # get data
  nucData = {}
  for dataDim in spectrum.sortedDataDims():
    if dataDim.className == 'FreqDataDim':
      for ddr in dataDim.sortedDataDimRefs():
        isotopes = ddr.expDimRef.isotopeCodes
        ll = nucData.get(isotopes)
        if ll is None:
          nucData[isotopes] = ll = []
        ll.append(ddr)
  
  # check for inconsistencies
  if printChecks:
    for isotopes, ll in sorted(nucData.items()):
      # check for consistent base frequencies
      if len(ll) > 1:
        ddr0 = ll[0]
        baseFreq = ddr0.expDimRef.baseFrequency
        for ddr in ll[1:]:
          bf2 = ddr.expDimRef.baseFrequency
          if baseFreq is None or bf2 is None:
            if baseFreq == bf2:
              continue
          
          elif abs(bf2 - baseFreq) <= sfTolerance:
            continue
            
          print ("WARNING, %s baseFrequency differs: %s:%s; %s:%s"
                 % (isotopes, ddr0, baseFreq, ddr, bf2))
  
  # check for misreferencing
  listH = nucData.get(('1H',))
  listC = nucData.get(('13C',))
  if listH and listC:
    baseFreqH = listH[0].expDimRef.baseFrequency
    baseFreqC = listC[0].expDimRef.baseFrequency
    if baseFreqH and baseFreqC:
      ratio = baseFreqC/baseFreqH
      trueRatio = chemShiftRefRatios['13C']
      ratioDelta = (ratio/trueRatio - 1.0)
      if abs(ratioDelta) > sfTolerance:
 
        text = "Carbon Referencing out by %.3f" % (1.0e6 * ratioDelta)
        if abs(1.0 - ratio / shiftRatio_13C_TMS) < brukerTolerance:
          text = text + "\nReferencing was to TMS rather than DSS"
        result = (True,text)
 
        if fixErrors:
          for ddr in listC:
            edr = ddr.expDimRef
            ddr.refValue = ddr.refValue + 1.0e6 * ratioDelta
            #sfDelta = edr.sf - edr.baseFrequency
            #newCarbonBase = edr.baseFrequency = baseFreqH * trueRatio
            #edr.sf = newCarbonBase + sfDelta
  #
  return result

def cloneExperiment(experimentOld, experimentNameNew, dataSourceNameNew=None, cloneDataFile=True):

  keys = ('date', 'details', 'nmrTubeType', 'numDim', 'numScans', 'sampleState',
          'sampleVolume', 'spinningAngle', 'spinningRate', 'volumeUnit',
          'refExperiment', 'shiftList')
  dd = dict(((x,getattr(experimentOld,x)) for x in keys))
  dd['name'] = experimentNameNew
  experimentNew = experimentOld.nmrProject.newExperiment(**dd)

  refExperiment = experimentOld.refExperiment

  expDimMap = {}
  expDimRefMap = {}
  for expDimOld in experimentOld.sortedExpDims():
    expDimNew = experimentNew.findFirstExpDim(dim=expDimOld.dim)
    expDimMap[expDimOld] = expDimNew
    keys = ('isAcquisition', 'refExpDim')
    dd = dict(((x,getattr(expDimOld,x)) for x in keys))
    for key in dd:
      setattr(expDimNew, key, dd[key])

    for expDimRefOld in expDimOld.sortedExpDimRefs():
      keys = ('baseFrequency', 'isAxisReversed', 'isFolded',
              'isotopeCodes', 'maxAliasedFreq', 'measurementType', 'minAliasedFreq',
              'name', 'nominalRefValue', 'sf', 'unit', 'variableIncrFraction',
              'refExpDimRef')
      dd = dict(((x,getattr(expDimRefOld,x)) for x in keys))
      expDimRefNew = expDimNew.newExpDimRef(**dd)
      expDimRefMap[expDimRefOld] = expDimRefNew

  for expTransferOld in experimentOld.expTransfers:
    keys = ('isDirect', 'mixingTime', 'transferSubType', 'transferType')
    dd = dict(((x,getattr(expTransferOld,x)) for x in keys))
    dd['expDimRefs'] = [expDimRefMap[expDimRefOld] for expDimRefOld in expTransferOld.sortedExpDimRefs()]
    expTransferNew = experimentNew.newExpTransfer(**dd)
        
  dataSourceOld = experimentOld.findFirstDataSource()

  keys = ('dataType', 'details', 'isSimulated', 'noiseLevel', 'numDim',
          'numShapes', 'numSparsePoints', 'recordNumber', 'scale',
          'signalLevel', 'refNmrSpectra', 'snMethod')
  if dataSourceNameNew is None:
    keys += ('name',)
  if cloneDataFile:
    keys += ('isNormalStorage', 'storageDetails',
             'compressMethod', 'dataStore', 'processMethod')
  dd = dict(((x,getattr(dataSourceOld,x)) for x in keys))
  if dataSourceNameNew is not None:
    dd['name'] = dataSourceNameNew

  dataSourceNew = experimentNew.newDataSource(**dd)

  for dataDimOld in dataSourceOld.sortedDataDims():
    keys = ('dim', 'fileDim', 'isComplex', 'numPoints', 'shapeSerial', 'unit')
    dd = dict(((x,getattr(dataDimOld,x)) for x in keys))
    expDimOld = dataDimOld.expDim
    dd['expDim'] = expDimMap[expDimOld]
    if dataDimOld.className == 'FidDataDim':
      keys = ('alternateSign', 'firstValue', 'negateImaginary',
              'numPointsValid', 'oversamplingInfo', 'phase0', 'phase1',
              'pointOffset', 'valuePerPoint')
      dd.update(((x,getattr(dataDimOld,x)) for x in keys))
      dataDimNew = dataSourceNew.newFidDataDim(**dd)
    elif dataDimOld.className == 'FreqDataDim':
      keys = ('numPointsOrig', 'phase0', 'phase1',
              'pointOffset', 'valuePerPoint', 'predictMethod')
      dd.update(((x,getattr(dataDimOld,x)) for x in keys))
      dataDimNew = dataSourceNew.newFreqDataDim(**dd)

      keys = ('localValuePerPoint', 'refPoint', 'refValue')
      dataDimRefOld = getPrimaryDataDimRef(dataDimOld)
      expDimRefOld = dataDimRefOld.expDimRef
      dd = dict(((x,getattr(dataDimRefOld,x)) for x in keys))
      dd['expDimRef'] = expDimRefMap[expDimRefOld]
      dataDimRefNew = dataDimNew.newDataDimRef(**dd)
    elif dataDimOld.className == 'SampledDataDim':
      keys = ('conditionVaried', 'details', 'pointErrors', 'pointValues')
      dd.update(((x,getattr(dataDimOld,x)) for x in keys))
      dataDimNew = dataSourceNew.newSampledDataDim(**dd)

  return experimentNew

def defaultExperiment(nmrProject, experimentName, isotopeCodes, dataSourceName=None, dataPath=None):

  numDim = len(isotopeCodes)
  experiment = nmrProject.newExperiment(numDim=numDim, name=experimentName)

  expDims = experiment.sortedExpDims()
  expDimRefs = numDim * [None]
  for n, expDim in enumerate(expDims):
    isotope = isotopeCodes[n]
    if isotope:
      sf = 500.0 * chemShiftRefRatios.get(isotope, 0.1)
      dd = { 'isotopeCodes': (isotope,),
             'sf': sf,
             'unit': 'ppm',
           }
      expDimRef = expDim.newExpDimRef(**dd)
      expDimRefs[n] = expDimRef

  if dataSourceName is None:
    dataSourceName = '1'

  dd = { 'name': dataSourceName,
         'numDim': numDim,
         'dataType': 'processed',
       }
  dataSource = experiment.newDataSource(**dd)

  numPointsArray = []
  for n, expDimRef in enumerate(expDimRefs):
    if expDimRef:
      if n == 0:
        numPoints = 1024
      else:
        numPoints = 64
      valuePerPoint = 0.1
      dd = { 'numPoints': numPoints,
             'numPointsOrig': numPoints,
             'valuePerPoint': valuePerPoint,
             'dim': n+1,
             'expDim': expDims[n],
             'isComplex': False,
           }
      dataDim = dataSource.newFreqDataDim(**dd)

      refPoint = numPoints / 2
      refValue = 4.72
      dd = { 'refPoint': refPoint,
             'refValue': refValue,
             'expDimRef': expDimRef,
           }
      dataDimRef = dataDim.newDataDimRef(**dd)
    else:
      numPoints = 10
      pointValues = 0.1 * range(numPoints)
      dd = { 'numPoints': numPoints,
             'numPointsOrig': numPoints,
             'pointValues': pointValues,
             'dim': n+1,
             'expDim': expDims[n],
             'isComplex': False,
           }
      dataDim = dataSource.newSampledDataDim(**dd)
    numPointsArray.append(numPoints)

  if dataPath:
    dataUrl, filePath = getDataStoringFromFilepath(experiment.root, dataPath)
    # blockSizes definition is pretty arbitrary
    blockSizes = len(numPointsArray) *[1]
    blockSizes[0] = numPointsArray[0]
    dataSource.dataStore = createBlockedMatrix(dataUrl, dataPath, numPointsArray, blockSizes)

  return experiment

def _calcStats(block_file, first, last):

  values = block_file.getValues(first, last)
  n = len(values)
  s = sum(values)
  s2 = sum([v*v for v in values])

  return n, s, s2

def getRegionStats(spectrum, region):

  from ccpnmr.analysis.core.PeakBasic import findDataDimRegions
  from ccpnmr.analysis.core.Util import getAnalysisSpectrum

  block_file = spectrum.block_file

  if not block_file:
    return []

  analysisSpectrum = getAnalysisSpectrum(spectrum)
  analysisProject = analysisSpectrum.analysisProject
  axisUnit = analysisProject.findFirstAxisUnit(unit='ppm')

  ndim = spectrum.numDim

  if len(region) != ndim:
    raise Exception('should have ndim (= %s) = len(region) (= %s)' % (ndim, len(region)))

  regions = ndim * [0]
  dataDims = spectrum.sortedDataDims()
  rdim = 0
  for dataDim in dataDims:
    dim = dataDim.dim - 1
    if dataDim.className == 'FreqDataDim':
      (r0, r1) = ppmDataDimBoundedRegion(region[rdim], dataDim)
      rdim += 1

      (regions[dim], min_point, max_point) = \
           findDataDimRegions(dataDim, (r0, r1), axisUnit, thickness=0, addOneToUpperPoint=False)
    else:
      npoints = dataDim.numPoints
      (r0, r1) = region[rdim]
      r0 = int(r0-1)
      r1 = int(r1)
      rdim += 1
      regions[dim] = [(r0, r1)]

  ntotal = stotal = s2total = 0
  first = ndim * [0]
  last = ndim * [0]
  (nregions, cumulative) = cumulativeProductArray([len(reg) for reg in regions])
  for r in range(nregions):
    a = arrayOfIndex(r, cumulative)
    reg = [ regions[i][a[i]] for i in range(len(a)) ]
    for dataDim in spectrum.dataDims:
      dim = dataDim.dim - 1
      first[dim] = int(reg[dim][0])
      last[dim] = int(reg[dim][1])

    n, s, s2 = _calcStats(block_file, first, last)
    ntotal += n
    stotal += s
    s2total += s2

  if ntotal > 0:
    avg = stotal / float(ntotal)
    s2total /= float(ntotal)
    std = math.sqrt(s2total - avg*avg)
  else:
    avg = std = 0

  scale = spectrum.scale

  return ntotal, scale*avg, scale*std
