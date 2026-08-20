"""
======================COPYRIGHT/LICENSE START==========================

CouplingBasic.py: Part of the CcpNmr Analysis program

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
from math import cos, sin, pi, sqrt

from ccp.api.nmr import Nmr

from memops.general import Implementation

from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes
from ccpnmr.analysis.core.WindowBasic import getDataDimAxisMapping
from ccpnmr.analysis.core.PeakBasic import getPeakDimPosition
from ccpnmr.analysis.core.Util import getAnalysisSpectrum, getAnalysisDataDim
from ccpnmr.analysis.core.MoleculeBasic import DEFAULT_ISOTOPES
from ccpnmr.analysis.core.UnitConverter import unit_converter

try:
  from memops.gui.MessageReporter import showWarning
except ImportError:
  from memops.universal.MessageReporter import showWarning

COUPLING_TYPES = ('JCoupling', 'Rdc', 'DipolarCoupling') 

EXPT_MULTIPLETS = {'H_H.DQF-COSY':'+-|-+',
                   'copy_H_H.DQF-COSY':'+-|-+',
                   'H[N].N coupled':'+|+',}

MULTIPLET_PEAK_DICT = {'+|+':((1,),(1,)),
                       '+|-':((1,),(-1,)),
                       '-|+':((-1,),(1,)),
                       '-|-':((-1,),(-1,)),
                       '++':(1,1),
                       '+-':(1,-1),
                       '-+':(-1,1),
                       '--':(-1,-1),
                       '+-|-+':((1,-1),(-1,1)),
                       '-+|+-':((-1,1),(-1,1)),
                       '+ | +':((1,0),(0,1)),
                       ' +|+ ':((0,1),(1,0)),
                       '++|++':((1,1),(1,1))}

MULTIPLET_TYPE = 'multiplet'


# TBD:

# Get obvious existing resonances for a coupling - e.g. from main assign
# Even if there are no J couplings
# Are J coupings sensible anyhow - actually tightly restricted by expt Type
# Assign jCoupling resonance atoms in Assignment Panel

def couplingNotifiers(notifyFunc):
  """
  Function to setup automatic curation when couplings ans peak
  clusters change

  .. describe:: Input
  
  Function (to issue notification)
  
  .. describe:: Output

  None
  """
  notifyFunc(updateClusterAssignments, 'ccp.nmr.Nmr.PeakCluster', '__init__')
  notifyFunc(updateClusterAssignments, 'ccp.nmr.Nmr.PeakCluster', 'setPeaks')
  notifyFunc(updateClusterAssignments, 'ccp.nmr.Nmr.PeakCluster', 'addPeak')
  notifyFunc(updateClusterAssignments, 'ccp.nmr.Nmr.PeakCluster', 'removePeak')

  for func in ('__init__','delete',):
    notifyFunc(updateContribNPeakDimAnnotation, 'ccp.nmr.Nmr.PeakDimContribN', func)

  for func in ('__init__','delete','setScalingFactor',
               'setDataDimRef', 'setPeakDimContribs',
               'addPeakDimContrib','removePeakDimContrib'):
    notifyFunc(updateComponentMeasurements, 'ccp.nmr.Nmr.PeakDimComponent', func)
  
  # PeakDimContrib notifiers are in AssignmentBasic
  # - see AssignmentBasic.makePeakDimAnnotation

def setupCouplingNotifiers():
  """
  Set up the notification calls to automatically update annotations
  couplings andpeak clusters when couplings change.
  
  .. describe:: Input
  
  None
  
  .. describe:: Output
  
  None
  """

  couplingNotifiers(Implementation.registerNotify)
  
def removeCouplingNotifiers():
  """
  Remove the notification calls to automatically update annotations
  couplings andpeak clusters when couplings change.
  
  .. describe:: Input
  
  None
  
  .. describe:: Output
  
  None
  """

  couplingNotifiers(Implementation.unregisterNotify)

def updateContribNPeakDimAnnotation(contrib):
  """
  Updates annotation string for the peak dimension of a given
  splitting contribution according to its assignment status
  
  .. describe:: Input
  
  Nmr.PeakDimContribN
  
  .. describe:: Output
  
  None
  """
  
  from ccpnmr.analysis.core.AssignmentBasic import makePeakDimAnnotation
  
  peakDim = contrib.peakDim
  for peakCluster in peakDim.peak.sortedPeakClusters():
    updateClusterAssignments(peakCluster, peakDim.peak)
    break # First only

  if not peakDim.isDeleted:
    makePeakDimAnnotation(peakDim)
  
  component = contrib.peakDimComponent
  if component:
    updateComponentMeasurements(component)
  
def updateComponentMeasurements(component):
  """
  Update the shifts /couplings of resonances assigned to the
  contributions of a particular peak dim component. Needs to be
  set because for a non-main component the bare peak dim position
  is not representative of the true shift, being a linear combination
  of ppm and coupling.
  
  .. describe:: Input
  
  Nmr.PeakDimComponent
  
  .. describe:: Output
  
  None
  """  
  
  from ccpnmr.analysis.core.AssignmentBasic import updateResonShift
  
  for contrib in component.peakDimContribs:
    expDimRef = component.dataDimRef.expDimRef
    peakDim = component.peakDim
  
    if isinstance(contrib, Nmr.PeakDimContribN):
      if expDimRef.unit == 'Hz':
        resonances = contrib.resonances
        experiment = expDimRef.expDim.experiment
        jCouplingList = getExperimentJCouplingList(experiment)
        jCoupling = jCouplingList.findFirstMeasurement(resonances=resonances)
 
        if not jCoupling:
           value = getPeakDimComponentSplitting(component)
           jCoupling = jCouplingList.newJCoupling(value=value,
                                        resonances=resonances)
         
        averageJCouplingValue(jCoupling)
    
    elif expDimRef.unit == 'ppm':
      updateResonShift(contrib.resonance, peakDim)
    
    # TBD consider what else? MQ coherence?
       
def averageJCouplingValue(jCoupling):
  """
  Recalculate a J-coupling value based upon the peakDimComponents
  to which it may be assigned via resonances. Sets any peak and
  peakDim links.
  
  .. describe:: Input
  
  ccp.nmr.Nmr.JCoupling
  
  .. describe:: Output
  
  Float
  """  
  hasApp = hasattr(jCoupling.root, 'application')

  sum1  = 0.0
  sum2  = 0.0
  N     = 0.0
  sigma2= 0.0
  sigma = 0.0
  mean  = None
  peakDims = set()
  peakDimsAdd =  peakDims.add
  
  resonances = jCoupling.resonances
  for resonance in resonances:
    for contrib in resonance.peakDimContribNs:
      component = contrib.peakDimComponent
      
      if not component:
        continue
    
      if not isinstance(contrib, Nmr.PeakDimContribN):
        continue
        
      if hasApp:
        weight = getAnalysisDataDim(component.dataDimRef.dataDim).chemShiftWeight
      else:
        weight = 1.0
  
      value = getPeakDimComponentSplitting(component)
      vw    = value * weight
      sum1 += vw
      sum2 += value * vw
      N    += weight
      peakDimsAdd(contrib.peakDim)

  if N > 0.0:
    mean  = sum1/N
    mean2 = sum2/N
    sigma2= abs(mean2 - (mean * mean))
    sigma = sqrt(sigma2)
    jCoupling.setValue( mean )
    jCoupling.setError( sigma )
    jCoupling.setPeakDims( peakDims )
    jCoupling.setPeaks( set([pd.peak for pd in peakDims]) )
    
  else:
    jCoupling.setError(0.0)
    jCoupling.setPeakDims([])
    jCoupling.setPeaks([])

  return mean

def assignPeakDimComponentResonance(component, resonance, tolerance=None):
  """
  Assign a chemical shift to a peakDim component
  (e.g. reduced dimensionality)
  
  .. describe:: Input
  
  ccp.nmr.Nmr.PeakDimComponent, ccp.nmr.Nmr.Resonance, Float
  
  .. describe:: Output
  
  ccp.nmr.Nmr.PeakDimContrib
  """
  
  # TBD version to cope with multiple resonances

  if resonance:
    dataDimRef = component.dataDimRef
    expDimRef = dataDimRef.expDimRef
    
    if resonance.isotopeCode not in expDimRef.isotopeCodes:
      return
    
    peakDimContrib = component.findFirstPeakDimContrib(resonance=resonance)
    if not peakDimContrib:
      experiment = expDimRef.expDim.experiment
      shiftList  = experiment.shiftList
      
      if not tolerance:
        tolerance = getAnalysisDataDim(dataDimRef.dataDim).assignTolerance
      
      if not shiftList:
        nmrProject = experiment.nmrProject
        shiftList  = nmrProject.findFirstMeasurementList(className='ShiftList')
 
        if not shiftList:
          shiftList = nmrProject.newShiftList(name='ShiftList 1',
                                              details='Assignment Default',
                                              unit='ppm')
 
        experiment.shiftList = shiftList
      
      shift = resonance.findFirstShift(parentList=shiftList)
      if shift:
        value = getPeakDimComponentSplitting(component)
        if abs(shift.value-value) > tolerance:
          return
      
      peakDim = component.peakDim
      peakDimContrib = peakDim.newPeakDimContrib(resonance=resonance,
                                          peakDimComponent=component)
  
  else:
    for contrib in component.peakDimContribs:
      contrib.delete()
  
    return

  return peakDimContrib
    

def getCouplingTolerance(component, dataDimRef):
  """
  Get the coupling assignment tolerance (usually Hz)
  for a given data dim.
  
  .. describe:: Input
  
  ccp.nmr.Nmr.FreqDataDim
  
  .. describe:: Output
  
  Float
  """

  freqDataDim = dataDimRef.dataDim
  delta = getAnalysisDataDim(dataDimRef.dataDim).assignTolerance
  unit = dataDimRef.expDimRef.unit
  
  if unit != 'point':
    delta /= dataDimRef.valuePerPoint
  
  tolerance = component.dataDimRef.pointToValue(delta)

  return tolerance

def assignPeakDimComponentCoupling(component, jCoupling, tolerance=None):
  """
  Assign a jCoupling to a peakDim component for a splitting.
  Any coupling value will be updated by CouplingBasic notifiers.
  JCoupling can be None to remove assignment.
  
  .. describe:: Input
  
  ccp.nmr.Nmr.PeakDimComponent, ccp.nmr.Nmr.JCoupling, Float
  
  .. describe:: Output
  
  ccp.nmr.Nmr.PeakDimContribN
  """

  if jCoupling:
    resonances = jCoupling.resonances
    isotopes = set([r.isotopeCode for r in resonances])
    peakDim = component.peakDim
    dataDimRef = component.dataDimRef
    expDimRef = dataDimRef.expDimRef
    experiment = expDimRef.expDim.experiment
    
    if isotopes != set(expDimRef.isotopeCodes):
      return

    peakDimContrib = component.findFirstPeakDimContrib(resonances=resonances)
    
    if not peakDimContrib:
      if not tolerance:
        tolerance = getCouplingTolerance(component, peakDim.dataDimRef)
        
      value = getPeakDimComponentSplitting(component)
      if abs(jCoupling.value-value) > tolerance:
        return 
        
      peakDimContrib = peakDim.newPeakDimContribN(resonances=resonances,
                                             peakDimComponent=component)
  
  else:
    for contrib in component.peakDimContribs:
      contrib.delete()
  
    peakDimContrib = None

  return peakDimContrib

def findMatchingCouplings(component, tolerance=None):
  """
  Get the J couplings that match the splitting of a 
  peakDimComponent, within the set tolerances
  
  .. describe:: Input
  
  ccp.nmr.Nmr.PeakDimComponent
  
  .. describe:: Output
  
  List of 2-Tupe of (Float (delta), ccp.nmr.Nmr.JCoupling)
  """

  peakDim = component.peakDim
  peakList = peakDim.peak.peakList
  experiment = peakList.dataSource.experiment
  value = getPeakDimComponentSplitting(component)

  if not tolerance:
    tolerance = getCouplingTolerance(component, peakDim.dataDimRef)

  jCouplingList = experiment.jCouplingList
  jCouplings = []
  jCouplingsAppend = jCouplings.append
  
  if jCouplingList:
    isotopes = set(component.dataDimRef.expDimRef.isotopeCodes)
    measurements = jCouplingList.measurements
    
    for jCoupling in measurements:
      delta = abs(jCoupling.value-value)
      if delta > tolerance:
        continue
      
      isotopesB = set([r.isotopeCode for r in jCoupling.resonances])
      if isotopes != isotopesB:
        continue
    
      jCouplingsAppend((delta, jCoupling))
  
    jCouplings.sort()
  
  return jCouplings

def getPeakDimComponentSplitting(component):
  """
  Get the splitting value (usually in Hz, but possibly in ppm)
  for a given peakDim component.
  
  .. describe:: Input
  
  ccp.nmr.Nmr.PeakDimComponent
  
  .. describe:: Output
  
  Float
  """

  peakDim = component.peakDim
  dataDimRef = component.dataDimRef
  dataDimRef2 = peakDim.dataDimRef
  
  if dataDimRef is dataDimRef2:
    return abs(peakDim.value - peakDim.realValue)
    
  else:
    conv = unit_converter[(dataDimRef2.expDimRef.unit,'point')]
    realPoint = conv(peakDim.realValue, dataDimRef2)
    deltaPoints = abs(peakDim.position-realPoint)
    # Assume ref point is zero
  
    return abs(component.scalingFactor*dataDimRef.pointToValue(deltaPoints))

def assignPrimaryClusterCoupling(cluster):
  """
  Assign the couplings of a peak cluster to any main resonance pair.
  
  .. describe:: Input
  
  ccp.nmr.Nmr.PeakCluster
  
  .. describe:: Output
  
  ccp.nmr.Nmr.JCoupling
  """
  
  from ccp.api.nmr import Nmr
  
  # Can assign some couplings automatically, e.g.
  # Homonuclear COSY etc the resonances are the same as the shift ordinates
  # Even IPAP HSCQ etc the coupling is the shift assinged resonances
  # Generally the coupled isotope will be present in the refExperiment
  # thus there is less ambiguity about which nuclei it could possibly be
  # It is only unmeasured atom sites that could cause any ambiguity
  # Accordingly only two atom sites that match the coupled isotopes unambiguously
  # means automatic assignment of the coupling once the full shift assignment is made
  # Multiple assignments would mean that there are two couplings present
  # which may not be the same
  
  # TBD: If multiply assigned ask if coupling should also
  # If not ask which, or just keep first.

  if cluster.clusterType != MULTIPLET_TYPE:
    msg = 'Peak cluster must be of multiplet type'
    showWarning('Warning', msg)
    return

  peak = cluster.findFirstPeak()
  experiment = peak.peakList.dataSource.experiment
  refExperiment = experiment.refExperiment
  
  if not refExperiment:
    msg = 'Cannot assign couplings of an untyped experiment'
    showWarning('Warning', msg)
    return
  
  refType = refExperiment.name
  if refType not in EXPT_MULTIPLETS:
    msg = 'Assignment of couplings for experiment type "%s" not implemented' % refType
    showWarning('Warning', msg)
    return
  
  peakDims = peak.sortedPeakDims()
  isotopes = set([])
  resonancesCoup = set([])
  for peakDim in peakDims:
    dataDimRef = getCouplingDataDimRef(peakDim.dataDim)
    
    if not dataDimRef:
      continue
    
    # TBD: Sort out peakContribs
      
    for contrib in peakDim.peakDimContribs:
      if not contrib.peakDimComponent:
        resonancesCoup.add(contrib.resonance)
    
    isotopes.update(dataDimRef.expDimRef.isotopeCodes)
  
  if len(resonancesCoup) > 2:
    return
  
  jCoupling = None
  for peakDim in peakDims:
    dataDimRef = getCouplingDataDimRef(peakDim.dataDim)
    expDimRef = dataDimRef.expDimRef
    
    if dataDimRef and (expDimRef.unit == 'Hz'):
      component = peakDim.findFirstPeakDimComponent(dataDimRef=dataDimRef)
      
      if component:
        splitIsotopes = set(component.dataDimRef.expDimRef.isotopeCodes)
        if splitIsotopes != isotopes:
          continue
        
        elif len(resonancesCoup) < 2:
          # not enough resonances, coupling goes 
          assignPeakDimComponentCoupling(component, None)
          
        else:
          if not jCoupling:
            experiment = expDimRef.expDim.experiment
            jCouplingList = getExperimentJCouplingList(experiment)
            jCoupling = jCouplingList.findFirstMeasurement(resonances=resonancesCoup)

            if not jCoupling:
              value = getPeakDimComponentSplitting(component)
              jCoupling = jCouplingList.newJCoupling(value=value, resonances=resonancesCoup)
 
          assignPeakDimComponentCoupling(component, jCoupling)
         
  # Spread assignments, curate annotations and shift values etc...
  updateClusterAssignments(cluster)
  
  return jCoupling

def getExperimentJCouplingList(experiment):
  """
  Get the J-Coupling list associated with an experiment, making and
  linking a new one if none exists.
  
  .. describe:: Input
  
  Nmr.Experiment
  
  .. describe:: Output
  
  Nmr.JCouplingList
  """

  jCouplingList = experiment.jCouplingList
   
  if not jCouplingList:
    nmrProject = experiment.nmrProject
    jCouplingList = nmrProject.findFirstMeasurementList(className='JCouplingList')
  
    if not jCouplingList:
      jCouplingList = nmrProject.newJCouplingList(name='default')
    
    experiment.jCouplingList = jCouplingList
 
  return jCouplingList
  
def getSpectrumMultipletPattern(spectrum):
  """
  Retrieve any multiplet pattern previously associated with a spectrum
  
  .. describe:: Input
  
  Nmr.DataSource
  
  .. describe:: Output
  
  String
  """

  analysisSpectrum = getAnalysisSpectrum(spectrum)
  return analysisSpectrum.multipletPattern


def setSpectrumMultipletPattern(spectrum, pattern):
  """
  Associate a give multiplet pattern with a spectrum for later recall.
  
  .. describe:: Input
  
  Nmr.DataSource, String
  
  .. describe:: Output
  
  None
  """
  
  analysisSpectrum = getAnalysisSpectrum(spectrum)
  analysisSpectrum.multipletPattern = pattern
 
 
def getPeakListsPeakClusters(peakLists, clusterType=MULTIPLET_TYPE):
  """
  Find any peak clusters (ordered by serial) that are represented by the input peak lists.
  
  .. describe:: Input
  
  List of Nmr.PeakLists
  
  .. describe:: Output
  
  List of Nmr.PeakClusters
  """

  clusterSet = set([])
  clusterSetAdd = clusterSet.add
  
  for peakList in peakLists:
    for peak in peakList.peaks:
      for cluster in peak.peakClusters:
        if cluster.clusterType == clusterType:
          clusterSetAdd(cluster)
  
  clusters = [(c.serial, c) for c in clusterSet]
  clusters.sort()
  clusters = [c[1] for c in clusters]
  
  return clusters

def updateClusterAssignments(cluster, refPeak=None):
  """
  Initialise a peak cluster object by enforcing consistent assignments
  across its component peaks and by setting its annotation approprately.
  
  .. describe:: Input
  
  Nmr.PeakCluster, Nmr.Peak
  
  .. describe:: Output
  
  None
  """
  
  from ccpnmr.analysis.core.AssignmentBasic import assignResToDim

  dimResonances = {}
  components = {}
  
  peaks = list(cluster.peaks)
  
  if refPeak:
    if len(refPeak.peakClusters) > 1:
      return
  
    inPeaks = [refPeak,]
    
  else:
    # Avoid mixing assignments
    inPeaks = [pk for pk in peaks if len(pk.peakClusters) < 2] 
  
  for peak in inPeaks:
    for peakDim in peak.peakDims:
      dim = peakDim.dim
      
      if dimResonances.get(dim) is None:
        dimResonances[dim] = set([])
        components[dim] = set([])
        
      for contrib in peakDim.peakDimContribs:
        
        if isinstance(contrib, Nmr.PeakDimContribN):
          component = contrib.peakDimComponent
 
          if component:
            resonances = contrib.resonances
            components[dim].add(resonances)
        
        else:
          resonance = contrib.resonance
          component = contrib.peakDimComponent
 
          if component:
            components[dim].add( (resonance,) )
 
          else:
            dimResonances[dim].add(resonance)
  
  clusterType = cluster.clusterType
  for peak in peaks:
  
    for peakDim in peak.peakDims:
      if not peakDim.dataDimRef:
        continue
        
      dim = peakDim.dim
 
      dataDimRef = None
      if clusterType == MULTIPLET_TYPE:
        dataDimRef = getCouplingDataDimRef(peakDim.dataDim)
      
       
      if dataDimRef:
        resonancesList = components.get(dim, [])
        
        if resonancesList:
          component = peakDim.findFirstPeakDimComponent(dataDimRef=dataDimRef)
 
          if not component:
            continue # Scaling factor unkown
          
          if component.dataDimRef.expDimRef.unit != 'Hz':
            continue
 
          for resonances in resonancesList:
            if len(resonances) != 2:
              continue
          
            contrib = component.findFirstPeakDimContrib(resonances=resonances)
            
            if not contrib:
              contrib = peakDim.newPeakDimContribN(resonances=resonances,
                                              peakDimComponent=component) 
            
        #else:
        # For non split dim, or non multiplet cluster: propagate assignments
        # Split dims assigned differently bacause simple position based
        # shift values are not appropriate
      
        # TBD: Symmetry requires mapping
      
      # Ensure consistent assignments, with refPeak when avail
      resonances = dimResonances[dim]
      for contrib in peakDim.peakDimContribs:
        if isinstance(contrib, Nmr.PeakDimContribN):
          continue
      
        if contrib.resonance not in resonances:
          for peakContrib in contrib.peakContribs:
            if len(peakContrib.peakDimContribs) < 2:
              peakContrib.delete()

          contrib.delete()
      
      for resonance in resonances:
        #if resonance not in existingResonancesA[peak][dim]:
        assignResToDim(peakDim, resonance)
  
  updateClusterAnnotation(cluster)

def updateClusterAnnotation(cluster):
  """
  Refresh the main peak cluster annotation (normally carries )
  
  .. describe:: Input
  
  Nmr.PeakCluster
  
  .. describe:: Output
  
  String
  """
  
  peaks = list(cluster.peaks)
  
  if peaks:
    peakDims = peaks[0].sortedPeakDims() 
    mainAnnotation = ','.join([pd.annotation or '' for pd in peakDims])
  else:
    mainAnnotation = None
  
  if cluster.annotation or mainAnnotation:
    cluster.annotation = mainAnnotation
  
  return mainAnnotation

def getCouplingDimAnnotation(cluster, dim):
  """
  Get a text string representing the assignment 
  status of couplings in a specified dimension
  represented by a peak cluster
  
  .. describe:: Input
  
  Nmr.PeakCluster, Int
  
  .. describe:: Output
  
  String
  """
  
  from ccpnmr.analysis.core.AssignmentBasic import makeResonanceGuiName
  
  annotation = ''
  
  if cluster.clusterType == MULTIPLET_TYPE:
    dimResonances = set([])
    dimResonancesAdd = dimResonances.add
    dimResonancesUpdate = dimResonances.update
    
    for peak in cluster.peaks:
      peakDim = peak.findFirstPeakDim(dim=dim)
    
      if peakDim:
        for component in peakDim.peakDimComponents:
          for contrib in component.peakDimContribs:
            if isinstance(contrib, Nmr.PeakDimContribN):
              dimResonancesUpdate(contrib.resonances)
            
            else:
              dimResonancesAdd(contrib.resonance)
  
    names = [makeResonanceGuiName(r) for r in dimResonances]
    names.sort()
    
    # TBD: Set and use component annotations
    
    annotation = 'J:' + ','.join(names)
  
  # Peak cluster annotations could carry couplings in Hz and possibly isotopes...
    
  return annotation

def getCouplingAnnotations(cluster):
  """
  Get a text string representing the assignment 
  status of couplings (in all dimensions)
  represented by a peak cluster
  
  .. describe:: Input
  
  Nmr.PeakCluster, Int
  
  .. describe:: Output
  
  String
  """

  annotations = []
  
  if cluster.clusterType == MULTIPLET_TYPE:
  
    dims = set([])
    dimsAdd = dims.add
    
    for peak in cluster.peaks:
      for peakDim in peak.peakDims:
        if peakDim.peakDimComponents:
          dimsAdd(peakDim.dim)
  
    dims = list(dims)
    dims.sort()
    
    annotations = [getCouplingDimAnnotation(cluster, dim) for dim in dims]
  
  return annotations

 
def getCoupledExpDims(experiment):
  """
  List the dimensions (ExpDims) of an experiment which carry couplings.
  
  .. describe:: Input
  
  Nmr.Experiment
  
  .. describe:: Output
  
  List of Nmr.ExpDims
  """

  expDims = []
  refExperiment = experiment.refExperiment
  
  if refExperiment:
    for expDim in experiment.expDims:
      refExpDim = expDim.refExpDim
      
      if not refExpDim: # Could be sampled dim
        continue 
    
      for refExpDimRef in refExpDim.refExpDimRefs:
        if refExpDimRef.coupledIsotopeCodes:
          expDims.append(expDim)
          break
  
  else:
    for expDim in experiment.expDims:
      for expDimRef in expDim.expDimRefs:
        if expDimRef.measurementType in COUPLING_TYPES:
          expDims.append(expDim)
          break
       
  return expDims

def getCouplingDataDimRef(dataDim):
  """
  Find the reference data dimension of the input data dimension that corresponds
  to a coupling referencing.
  
  .. describe:: Input
  
  Nmr.DataDim
  
  .. describe:: Output
  
  Nmr.DataDimRef
  """
  
  dataDimRef = None
  for dataDimRef0 in dataDim.dataDimRefs:
    if dataDimRef0.expDimRef.measurementType in COUPLING_TYPES:
      dataDimRef = dataDimRef0
      break
 
  return dataDimRef

def makeMultipletPeakCluster(peaks, multipletPattern, windowPane):
  """
  Make a peak cluster comprising of the selected peaks.
  The multiplet pattern is a recursive list of sub lists to give an object
  with the same dimensionality as the peak dims. Each element is
  either +1, 0 or -1 to specify the sign or absence of a multiplet
  peak. How to interpret a multiplet pattern is relative to input
  window pane, e.g. if spectrum is rotated. Automatically makes
  coupling expDimRefs when required.
  
  .. describe:: Input
  
  List of Nmr.Peaks, List of (Lists of)^numDim-1 Ints(-1,0,+1),
  Analysis.SpectrumWindowPane
  
  .. describe:: Output
  
  Nmr.PeakCluster
  """
  peaks      = list(peaks)
  peak0      = peaks[0]
  peakList0  = peak0.peakList
  nmrProject = peak0.topObject
  isotopes0  = getSpectrumIsotopes(peakList0.dataSource)
  shiftList  = None
  
  for peak in peaks[1:]:
    peakList = peak.peakList
    experiment = peakList.dataSource.experiment
    refExperiment = experiment.refExperiment

    if shiftList:
      if experiment.shiftList is not shiftList:
        msg = 'Experiments of selected peaks not all using the same shift list'
        showWarning('Clustering Failed', msg)
        return
    
    elif experiment.shiftList:
      shiftList = experiment.shiftList

    if not refExperiment:
        msg = 'Experiment %s does not have a reference experiment type'
        showWarning('Clustering Failed', msg % experiment.name)
        return
      
    
    if peakList is not peakList0:
      if getSpectrumIsotopes(peakList.dataSource) != isotopes0:
        msg = 'Selected peaks have different isotope dimensions'
        showWarning('Clustering Failed', msg)
        return

  # All peaks have to match - some overlap OK.
  clusters = list(getPeaksClusters(peaks, clusterType=MULTIPLET_TYPE))

  if clusters: # Use existing 
    peakCluster = clusters[0] # Arbitrary 
    
    if len(clusters) > 1:
      for cluster in clusters[1:]:
        deleteCluster(cluster) # Cleans up peakDimComponents too
  
  else:
    # Make a new one below after position checks
    peakCluster = None
  

  nDim = len(peak0.peakDims)
  dims = range(nDim)

  counts = [ 0.0 for x in dims]
  maxPos = [None for x in dims]
  minPos = [None for x in dims]
  centre = [None for x in dims]
  middle = [None for x in dims]
  
  boxWidths = [0.0 for x in dims]
  
  axisMapping = getDataDimAxisMapping(peakList0.dataSource, windowPane)
  
  for peak in peaks:
    for peakDim in peak.peakDims:
      i = peakDim.dim-1
      position   = peakDim.value 
      counts[i] += 1.0

      if (maxPos[i] is None) or (position > maxPos[i]):
        maxPos[i] = position
        
      if (minPos[i] is None) or (position < minPos[i]):
        minPos[i] = position
  
  gridSizes = []
  
  numCells = 1
  element = multipletPattern
  while type(element) in (type(()), type([])):
    size = len(element)
    gridSizes.append(size)
    element = element[0]
    numCells *= size

  for i in dims:
    mean  = (maxPos[i] + minPos[i])/2.0
    n =  float(gridSizes[i])-1.0
    
    if n:
      width = (maxPos[i] - minPos[i])/n
    else:
      width = (maxPos[i] - minPos[i])*3.0
     
    boxWidths[i] = width
    centre[i] = (mean-(width/4.0),mean+(width/4.0))
    middle[i] = mean
    
    if not counts[i]:
      showWarning('Clustering Failed',
                  'Peak dimension corrupt')
      return  # Ought never get here

  cellCoords = []
  for i in range(numCells):
    coords = []
    
    k = 1
    for size in gridSizes:
      coords.append((i/k) % size)
      k += 1      
    
    cellCoords.append(coords)

  peakCell  = {}
  peakSigns = []
  j = 0
  for coords in cellCoords:
    region = []
    multipletElement = list(multipletPattern)

    for i in dims:
      halfWidth = boxWidths[i]/2.0
      multipletElement = multipletElement[coords[i]]
      
      if axisMapping['y'].dim == i+1:
        midPt  = minPos[i] + (coords[i]*boxWidths[i])
      else:
        midPt  = maxPos[i] - (coords[i]*boxWidths[i])
        
      bounds = (midPt-halfWidth, midPt+halfWidth)
      region.append(bounds)

    for peak in peaks:
    
      if peakCell.get(peak) is None:
        for i in dims:
          peakDim = peak.findFirstPeakDim(dim=i+1)
          value = peakDim.value
   
          if (value < region[i][0]) or (value > region[i][1]):
            break
            
        else:
          peakCell[peak] = j
    
    peakSigns.append(multipletElement)
    j += 1

  centrePeaks = {}
  for size in gridSizes:
    if size % 2 == 0:
      for peak in peaks:
        for i in dims:
          peakDim = peak.findFirstPeakDim(dim=i+1)
          value = peakDim.value
          
          if (value < centre[i][0]) or (value > centre[i][1]):
            break
       
        else:
          centrePeaks[peak] = True
      
      break

  splitPeaks = [pk for pk in peaks if not centrePeaks.get(peak)]

  if len(splitPeaks) > len(peakSigns):
    patternText = '?'
    for pattern in MULTIPLET_PEAK_DICT:
      if MULTIPLET_PEAK_DICT[pattern] == multipletPattern:
        patternText = pattern
        break
    
    msg = 'More peaks selected than cluster pattern "%s"' % patternText
    msg += ' allows for (excluding any central peak).'
    showWarning('Clustering Failed', msg)
    return    
  
  if peakCluster:
    # Redefine peaks to any new selection
    if peakCluster.peaks != frozenset(peaks):
      peakCluster.setPeaks(peaks)

  else: # Make new
     peakCluster = nmrProject.newPeakCluster(peaks=peaks, clusterType=MULTIPLET_TYPE)

  assignments = {}

  for peak in splitPeaks:
    
    intensity = peak.findFirstPeakIntensity(intensityType='height')
    
    if not intensity:
      continue 
    
    dataDimRefs = []
    spectrum = peak.peakList.dataSource
    for j in dims:
      dataDim = spectrum.findFirstDataDim(dim=j+1)
    
      if dataDim.className == 'FreqDataDim':
        dataDimRef = getCouplingDataDimRef(dataDim)
        
        if not dataDimRef:
          expDim = dataDim.expDim
          
          expDimRef0 = expDim.findFirstExpDimRef(measurementType='Shift')
          if not expDimRef0:
            expDimRef0 = expDim.findFirstExpDimRef(measurementType='shift')
          
          refExpDimRef = expDimRef0.refExpDimRef
          
          if not refExpDimRef:
            msg  = 'Experiment %s is missing a reference dimension:'
            msg += ' Check experiment type'
            showWarning('Clustering Failed', msg % expDim.experiment.name)
            return
          
          isotopes = refExpDimRef.coupledIsotopeCodes
          if not isotopes:
            dataDimRefs.append(None)
            continue
            
          expDimRef = expDim.newExpDimRef(sf=1.0,
                                          isotopeCodes=isotopes,
                                          measurementType='JCoupling',
                                          isFolded=False,
                                          unit='Hz',
                                          isAxisReversed=False)
          
          dataDimRef = dataDim.newDataDimRef(expDimRef=expDimRef,
                                             refPoint=0,
                                             refValue=0.0)
                                             
        dataDimRefs.append(dataDimRef)
      
      else:
        dataDimRefs.append(None)
    
    if dataDimRefs == [None for j in dims]:
      msg = 'Reference dimensions for experiment %s have no coupled isotopes'
      showWarning('Clustering Failed', msg % spectrum.experiment.name)
    
    height = intensity.value
    i      = peakCell[peak]
    sign   = peakSigns[i]

    gridSizes.reverse()

    if sign and (sign * height > 0): # Same sign
      coord = cellCoords[i]
      
      for j in dims:
        dataDimRef = dataDimRefs[j]
        
        if not dataDimRef:
          continue
        
        size      = gridSizes[j]
        mid       = (size-1)/2.0
        scaling   = int(2*(coord[j] - mid))
        peakDim   = peak.findFirstPeakDim(dim=j+1)
        component = peakDim.findFirstPeakDimComponent(dataDimRef=dataDimRef)
        peakDim.realValue = middle[j]
        
        #print dataDimRef, height, i , sign, size, mid, scaling, peakDim

        if component:
          component.scalingFactor=scaling
                        
        else:
          component = peakDim.newPeakDimComponent(scalingFactor=scaling,
                                                  dataDimRef=dataDimRef)
  
  # Spread assignments, make C objects, update shifts etc
  updateClusterAssignments(peakCluster) 
         
  return peakCluster

def getClusterCoupling(peakCluster, dim):
  """
  Determine the coupling value in Hz for the peaks of an input peak cluster
  along the dimension corresponding to the input dimension number.
  
  .. describe:: Input
  
  Nmr.PeakCluster, Int (Nmr.DataDim.dim)
  
  .. describe:: Output
  
  Float (coupling value in Hz)
  """

  coupling = None
  
  peaks = list(peakCluster.peaks)
  peak0 = peaks[0]
  
  peakList0 = peak0.peakList
  isotopes0 = getSpectrumIsotopes(peakList0.dataSource)
  
  componentDict = {}
    
  for peak in peaks:
    peakList = peak.peakList
    
    if peakList is not peakList0:
      if getSpectrumIsotopes(peakList.dataSource) != isotopes0:
        showWarning('Clustering Failed',
                    'Selected peaks have different isotope dimensions')
        return
    
    peakDim = peak.findFirstPeakDim(dim=dim)
    dataDimRef = getCouplingDataDimRef(peakDim.dataDim)
    
    if peakDim and dataDimRef:
      component = peakDim.findFirstPeakDimComponent(dataDimRef=dataDimRef)
  
      if component:
        scale = component.scalingFactor
        
        if componentDict.get(scale) is None:
          componentDict[scale] = []
          
        componentDict[scale].append(component)
    
  scaleFactors = componentDict.keys()
  N = len(scaleFactors)
  num = 0.0
  zum = 0.0
  
  for i in range(N-1):
    scaleA = scaleFactors[i]
    componentsA = componentDict[scaleA]
    
    for j in range(i+1, N):
      scaleB = scaleFactors[j]
      componentsB = componentDict[scaleB]
  
      jFactor = abs(scaleA-scaleB)
      
      for componentA in componentsA:
        peakDimA  = componentA.peakDim
        positionA = peakDimA.position
        
        for componentB in componentsB:
          peakDimB    = componentB.peakDim
          positionB   = peakDimB.position
          deltaPoints = abs(positionA-positionB)
          value = jFactor*dataDimRef.pointToValue(deltaPoints)
           
          num += 1.0
          zum += value
  
  if num > 0.0:
    coupling = zum / num
   
  return coupling

def getPeaksClusters(peaks, clusterType=MULTIPLET_TYPE):
  """
  Find the peak clusters that may be represented by the input peaks.
  Cluster type can be set to None to get ant/all clusters
  
  .. describe:: Input
  
  List of Nmr.Peaks, Word (PeakCluster.clusterType)
  
  .. describe:: Output
  
  Nmr.PeakCluster
  """

  peaks = set(peaks)
  clusters = set([])
  clustersAdd = clusters.add
  
  for peak in peaks:
    for cluster in peak.peakClusters:
      if (clusterType is None) or (cluster.clusterType == clusterType):  
        if peaks == cluster.peaks:
          clustersAdd(cluster) 
  
  return clusters
  

def deleteCluster(cluster):
  """
  Remove the input peak cluster and any associated peakDimComponents
  that record the coupling contributions.
  
  .. describe:: Input
  
  Nmr.PeakCluster
  
  .. describe:: Output
  
  None
  """

  peaks = cluster.peaks
  cluster.setPeaks([]) # Notify on this once, not for all contribs

  for peak in peaks:
    for peakDim in peak.peakDims:
      peakDim.realValue = None
      
      for component in peakDim.peakDimComponents:
        for peakDimContrib in component.peakDimContribs:
          peakDimContrib.delete()
 
        component.delete()
 
  cluster.delete()
  
        
def getResidueJCoupling(jCouplingList, residue, atomNames):
  """
  Find the first JCoupling in a jCoupling list that corresponds to the named atoms
  in a given residue. Returns None if there is no match.
  
  .. describe:: Input
  
  Nmr.JCouplingList, MolSystem.Residue, 2-List of Words (MolSystem.Atom.names)
  
  .. describe:: Output
  
  Nmr.JCoupling
  """
  
  nmrProject = jCouplingList.topObject
  
  coupling   = None
  resonanceA, resonanceB = getCouplingAtomResonances(nmrProject, residue, atomNames, makeNew=False)
  if resonanceA and resonanceB:
    jCouplings = resonanceA.findAllJCouplings(parentList=jCouplingList)
    for jCoupling in jCouplings:
      if resonanceB in jCoupling.resonances:
        coupling = jCoupling

        break

  return coupling

def setResidueJCoupling(jCouplingList, residue, atomNames, value, error=None):
  """
  Sets the value of a JCoupling in a jCoupling list that corresponds to the named atoms.
  Makes a new JCoupling object if none can be found.
  
  .. describe:: Input
  
  Nmr.JCouplingList, MolSystem.Residue, 2-List of Words (MolSystem.Atom.names)
  
  .. describe:: Output
  
  Nmr.JCoupling
  """
  
  resonanceA, resonanceB = getCouplingAtomResonances(jCouplingList.topObject, residue, atomNames)
  coupling = None
  
  if resonanceA and resonanceB:
    jCouplings = resonanceA.findAllJCouplings(parentList=jCouplingList)
    for jCoupling in jCouplings:
      if resonanceB in jCoupling.resonances: 
        coupling = jCoupling
  
    if coupling is None:
      coupling = jCouplingList.newJCoupling(value=value, resonances=(resonanceA, resonanceB))
    else:
      coupling.setValue(value)

    if error is not None:
      coupling.setError(error)
 
  return coupling
 
 
def getCouplingAtomResonances(nmrProject, residue, atomNames, makeNew=True):
  """
  Get resonances that correspond to two named atoms in a given residue context
  - could include sequential neighbours. A prerequisite for accessign the Jcoupling
  etc for those atoms.
  
  .. describe:: Input
  
  Nmr.NmrProject, MolSystem.Residue, 2-List of Words (MolSystem.Atom.names), Boolean
  
  .. describe:: Output
  
  Nmr.Resonance, Nmr.Resonance
  """
  
  # Add FindLinkedresidue calls here so not relying on seqCode
  
  from ccpnmr.analysis.core.AssignmentBasic import assignAtomsToRes

  atomNames = tuple(atomNames)

  if 'C' in atomNames:
    if 'H' in atomNames:
      residueA = residueB = residue
    else:
      if 'C' == atomNames[0]:
        residueA = residue.chain.findFirstResidue(seqCode=residue.seqCode-1)
        residueB = residue
      else:
        residueA = residue
        residueB = residue.chain.findFirstResidue(seqCode=residue.seqCode-1)
  
  else:  
    residueA = residueB = residue
  
  if not (residueA and residueB):
    return None, None
  
  atomA = residueA.findFirstAtom(name=atomNames[0])
  atomB = residueB.findFirstAtom(name=atomNames[1])
  
  resonanceA = None
  resonanceB = None
  
  if atomA and atomB:
    atomSetA = atomA.atomSet
    atomSetB = atomB.atomSet
    
    if atomSetA and atomSetB:
      resonanceSetsA = list(atomSetA.resonanceSets)
      resonanceSetsB = list(atomSetB.resonanceSets)
      
      if not resonanceSetsA:
        if makeNew:
          isotope = DEFAULT_ISOTOPES.get(atomA.chemAtom.elementSymbol)
          resonanceA = nmrProject.newResonance(isotopeCode=isotope)
          resonanceSetA = assignAtomsToRes([atomSetA,], resonanceA)
      
      else:
        resonanceSetA = resonanceSetsA[0]
        for resonanceSet in resonanceSetsA:
          if len(resonanceSet.resonances) < len(resonanceSetA.resonances):
            resonanceSetA = resonanceSet

        resonanceA = resonanceSetA.findFirstResonance()

      if not resonanceSetsB:
        if makeNew:
          isotope = DEFAULT_ISOTOPES.get(atomB.chemAtom.elementSymbol)
          resonanceB = nmrProject.newResonance(isotopeCode=isotope)
          resonanceSetB = assignAtomsToRes([atomSetB,], resonanceB)
      
      else:
        resonanceSetB = resonanceSetsB[0]
        for resonanceSet in resonanceSetsB:
          if len(resonanceSet.resonances) < len(resonanceSetB.resonances):
            resonanceSetB = resonanceSet
        
        resonanceB = resonanceSetB.findFirstResonance()
  
  return resonanceA, resonanceB


def calculateKarplusCoupling(angle, couplingConstantTriple, angleOffset = 0):
  """
  Calculate the Karplus coupling, given angle, coupling constant and optional angle offset.
  Coupling = c0 + c1*cos(angle+angleOffset) + c2*cos(angle+angleOffset)
  
  .. describe:: Input
  
  angle in radians, (c0, c1, c2) tuple, angle offset in radians
  
  .. describe:: Output
  
  coupling
  """

  (c0, c1, c2) = couplingConstantTriple
  angle += angleOffset

  return c0 + c1*cos(angle) + c2*cos(2*angle)

def calculateKarplusCouplings(angle, couplingConstants, angleOffsets):
  """
  Calculate the Karplus couplings, given angle, coupling constants, angle offset.
  Coupling = c0 + c1*cos(angle+angleOffset) + c2*cos(angle+angleOffset)
  
  .. describe:: Input
  
  angle in radians, list of (c0, c1, c2) tuples, list of angle offsets in radians
  
  .. describe:: Output
  
  list of couplings
  """

  couplings = []
  for n in range(len(couplingConstants)):
    coupling = calculateKarplusCoupling(angle, couplingConstants[n], angleOffsets[n])
    couplings.append(coupling)

  return couplings


if __name__ == '__main__':

  #angle = pi * 134.0 / 180.0
  angle = 134.0
  couplingConstants = [[5.7,-1.2,2.5],
                       [3.1,-0.1,1.1],
                       [5.7,-1.1,2.1],
                       [3.0,-0.1,1.1],
                       [5.7,-1.2,2.5],
                       [3.1,-0.1,1.1]]
  angleOffsets = [-60, 180, 60, 120, 0, -120]
  angleOffsets = [pi*a/180.0 for a in angleOffsets]

  couplings = calculateKarplusCouplings(angle, couplingConstants, angleOffsets)
  coupling = calculateKarplusCoupling(angle, couplingConstants[0], angleOffsets[0])

  print coupling, couplings
