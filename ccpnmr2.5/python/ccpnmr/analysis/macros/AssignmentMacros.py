
"""
======================COPYRIGHT/LICENSE START==========================

AssignmentMacros.py: Part of the CcpNmr Analysis program

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
from ccpnmr.analysis.core.AssignmentBasic import assignResToDim, findSpinSystem, addSpinSystemResonance, makePeakDimAnnotation
from ccpnmr.analysis.core.ExperimentBasic import getSpectraByType, getSpectrumIsotopes
from ccpnmr.analysis.core.ExperimentBasic import getPrimaryDataDimRef
from ccpnmr.analysis.core.PeakBasic import searchPeaks, findPeaks
from ccpnmr.analysis.core.MoleculeBasic   import makeResidueAtomSets
from ccpnmr.analysis.core.UnitConverter import pnt2ppm

def refreshAssignmentStatus(argServer=None, project=None):

  assert argServer or project
  
  if argServer:
    project =  argServer.getProject()
    
  if hasattr(project, 'atomSetMappings'):
    for atomSetMapping in project.currentNmrProject.atomSetMappings:
      atomSetMapping.delete()
      
  for molSystem in project.molSystems:
    for chain in molSystem.chains:
      for residue in chain.residues:
        residue.guiAtomSets = {}
        makeResidueAtomSets(residue)

  for resonanceSet in project.currentNmrProject.resonanceSets:
    argServer.parent.initResonanceSet(resonanceSet)

def assignSpinSystemPerPeak(argServer=None, peaks=None, dims=None):
  assert argServer or peaks

  if not peaks:
    peaks = argServer.getCurrentPeaks()

  if not peaks:
    return
  
  for peak in peaks:
    spinSystem = None
    resonances = []
    for i, peakDim in enumerate(peak.sortedPeakDims()):
      if (dims is None) or (i in dims):
        if len(peakDim.peakDimContribs) == 1:
          resonance = peakDim.findFirstPeakDimContrib().resonance
          resonances.append( resonance )
    
    for resonance in resonances:
      if resonance.resonanceGroup:
        spinSystem = resonance.resonanceGroup
        
    for resonance in resonances:
      if spinSystem is None:
        spinSystem = findSpinSystem(resonance)
      addSpinSystemResonance(spinSystem, resonance)
      
def assignAllNewResonances(argServer=None, peaks=None):
  assert argServer or peaks
  if not peaks:
    peaks = argServer.getCurrentPeaks()
  # e.g. Assign all of an unassigned HSQC to different Hn and N resonances 

  for peak in peaks:
    for peakDim in peak.peakDims:
      if len(peakDim.peakDimContribs) < 1:
        assignResToDim(peakDim)
        
def splitAmbigAssign(argServer=None, peakDim=None, keepResonance = None):

  assert argServer or  peakDim
  if not peakDim:
    peakDim = argServer.getCurrentPeakDim()
    keepResonance = argServer.askYesNo('Keep resonance?')
  
  # e.g. Convert Hb# assignments to Hb2 or Hb3 for all peaks assigned to a residue
  # must be ambiguous: we must have at least two contributions
  if len(peakDim.peakDimContribs) < 2:
    return
  
  peakDimContrib = peakDim.findFirstPeakDimContrib()
  # if we don't know which contribution to keep, keep the first
  if keepResonance is None:
    keepResonance = peakDimContrib.resonance
  
  # find the resonances in the ambiguous assignment
  resonances = []
  for contrib in peakDim.peakDimContribs:
    resonance = contrib.resonance
    if resonance not in resonances:
      resonances.append( resonance )
  
  # find all the peakDims that have one of the above resonances
  peakDims =[]
  for resonance in resonances:
    contribs = resonance.peakDimContribs
    for contrib in contribs:
      peakDim = contrib.peakDim
      if peakDim not in peakDims:
        peakDim.append( peakDim )
  
  # find all the peakDims that share the same ambiguous assignment (resonances)
  matchingPeakDims = []
  for peakDim in peakDims:
    otherResonances = []
    for contrib in peakDim.peakDimContribs:
      otherResonances.append( contrib.resonance )
      
    # look for resonances that match the initial peakDim
    i = 0
    for resonance in otherResonances:
      if resonance in resonances:
        i+= 1

    # must have all the resonances that were in the initial peakDim
    if i == len(resonances):
      matchingPeakDims.append( peakDim )
      
  #from the matching peakDims delete the unwanted ambiguous contributions
  for peakDim in matchingPeakDims:
    for contrib in peakDim.contribs:
      resonance = contrib.resonance
      if resonance in resonances:
        if resonance is not keepResonance:
          contrib.delete()
  
def refreshPeakAnnotations(argServer=None, project=None):

  assert argServer or project
  if not project:
    project = argServer.getProject()
    
  for experiment in project.currentNmrProject.experiments:
    for dataSource in experiment.dataSources:
      for peakList in dataSource.peakLists:
        for peak in peakList.peaks:
          for peakDim in peak.peakDims:
            makePeakDimAnnotation(peakDim)

def test(argServer=None):

  answer = argServer.getCurrentPeak()
  print "Current peak is %s" % answer

def initialiseHSQC(argServer=None, spectrum=None):

  assert argServer or spectrum
  if not spectrum:
    spectra  = getSpectraByType(argServer.getProject(),'H[N]')
    spectra.extend( getSpectraByType(argServer.getProject(),'H[C]') )
    spectrum = argServer.getSpectrum(spectra)
  
  if spectrum:
    peakList = argServer.getPeakList(spectrum)
    peaks    = peakList.peaks
    if peaks:
      assignAllNewResonances(peaks=peaks)
      assignSpinSystemPerPeak(peaks=peaks)

  else:
    argServer.showWarning('Spectrum not found')

def initialiseHNCOorHNCOCA(argServer=None, spectrum=None):

  assert argServer or spectrum
  if not spectrum:
    spectra  = getSpectraByType(argServer.getProject(),'H[N[CO]]')
    spectra.extend( getSpectraByType(argServer.getProject(),'H[N[co[CA]]]') )
    spectrum = argServer.getSpectrum(spectra)
  
  if spectrum:
    peakList = argServer.getPeakList(spectrum)
    isotopes = getSpectrumIsotopes(spectrum)
    dims     = []
    for isotope in ('1H','15N'):
      if isotope in isotopes:
        dims.append(isotopes.index(isotope))
    
    peaks    = peakList.peaks
    if peaks and dims:
      assignAllNewResonances(peaks=peaks)
      assignSpinSystemPerPeak(peaks=peaks,dims=dims)
      
  else:
    argServer.showWarning('Spectrum not found')
    
    
def addMarksToPeaks(argServer=None, peaks=None):
  """Descrn: Adds position line markers to the selected peaks.
     Inputs: ArgumentServer, List of Nmr.Peaks
     Output: None
  """
  assert argServer or peaks

  from ccpnmr.analysis.core.MarkBasic import createPeakMark
  
  if not peaks:
    peaks = argServer.getCurrentPeaks()

  for peak in peaks:
    createPeakMark(peak)
       
def pickAssignSpecFrom2dRoot(argServer, rootPeakList=None, targetPeakList=None):
  """
  Dummy macro for backward compatibility
  """
  pickAssignSpecFromRoot(argServer, rootPeakList, targetPeakList)
       

def pickAssignSpecFromRoot(argServer, rootPeakList=None, targetPeakList=None):

  toleranceDict = {'1H':0.03,'13C':0.1,'15N':0.15}
  
  assert argServer or (rootPeakList and targetPeakList)

  if argServer:
    project = argServer.getProject()
  else:
    project = rootPeakList.root
    
  spectra = []
  for experiment in project.currentNmrProject.sortedExperiments():
    for spectrum in experiment.sortedDataSources():
      spectra.append(spectrum)
    
  if not rootPeakList:
    rootSpec = argServer.getSpectrum(spectra)
    rootPeakList = argServer.getPeakList(rootSpec)

  rootSpec = rootPeakList.dataSource
  rootIsotopes = getSpectrumIsotopes(rootSpec)
  
  tolerances = []
  for isotope in rootIsotopes:
    tolerance = argServer.askFloat('%s tolerance' % isotope, toleranceDict.get(isotope, 0.1)) or toleranceDict.get(isotope, 0.1)
    tolerances.append(tolerance)
 
  if not targetPeakList:
    targetSpectra = []
    for experiment in project.currentNmrProject.experiments:
      for spectrum in experiment.dataSources:
        isotopes = getSpectrumIsotopes(spectrum)
        n = 0
        for isotope in rootIsotopes:
          if isotope in isotopes:
            isotopes.remove(isotope)
            n += 1
            
        if n == len(rootIsotopes):
          targetSpectra.append(spectrum)
        
    targetSpec  = argServer.getSpectrum(targetSpectra)
    targetPeakList = argServer.getPeakList(targetSpec)

  targetSpec = targetPeakList.dataSource

  if not targetSpec:
    argServer.showWarning('No suitable target spectra found.')
 
  # Determine mapping of root to target dimensions according to data dim isotopes.
  #  - Only ask user if there are multiple possibilities.
  
  mapping = []
  targetIsotopes = getSpectrumIsotopes(targetSpec)
  i = 0

  for isotope in rootIsotopes:
    if targetIsotopes.count(isotope) == 1:
      mapping.append(targetIsotopes.index(isotope))      
      
    else:
      j = 0
      matches = []
      for targetIsotope in targetIsotopes:
        if targetIsotope == isotope:
          matches.append(j)
        j += 1
      
      for j in matches[:-1]:  
        if argServer.askYesNo('Match root dimension %d (%s) to target dimension %d?' % (i+1,isotope,j+1)):
          mapping.append(j)
          break

      else:
        mapping.append(matches[-1])
        
    i += 1 
        
 
  fullRegion = []
  for dataDim in targetSpec.sortedDataDims():

    if dataDim.className == 'FreqDataDim':
      dataDimRef = getPrimaryDataDimRef(dataDim)
      valueMax   = pnt2ppm(1,dataDimRef)
      valueMin   = pnt2ppm(dataDim.numPoints,dataDimRef)
      fullRegion.append([valueMin,valueMax])
    else:
      fullRegion.append([1, dataDim.numPoints])
      
  M  = len (rootPeakList.peaks)
  c  = 0
  
  foundPeaks = 0
  for peak in rootPeakList.peaks:
    
    region = list(fullRegion)
    
    for i, peakDim in enumerate(peak.sortedPeakDims()):
      if peakDim.dataDimRef:
        error = tolerances[i]
        value = peakDim.value
        region[mapping[i]] = (value-error, value+error)
        
    peaks = searchPeaks([targetPeakList,], region)   
    peaks.extend( findPeaks(targetPeakList, region) )
    
    # possible to grow the region here of no peaks are found
     
    for peak2 in peaks:
      for intens in peak.peakIntensities:
        if str(intens.value) == 'inf':
          peaks.remove(peak2)
          peak2.delete()
          break
          
    foundPeaks += len(peaks)
    c += 1
    
    for i, peakDim in enumerate(peak.sortedPeakDims()):
      for peak2 in peaks:
        peakDim2 = peak2.sortedPeakDims()[mapping[i]]
        if len(peakDim2.peakDimContribs) <1:
          for contrib in peakDim.peakDimContribs:
            assignResToDim(peakDim2, contrib.resonance, doWarning=0)
                
    print 'Root peak %d of %d' % (c,M)
    
  if argServer:
    name = '%s:%s' % (targetPeakList.dataSource.experiment.name,targetPeakList.dataSource.name)
    argServer.showInfo('Picked %d peaks in %s' % (foundPeaks,name) )
    
  return targetPeakList.peaks

def assign3dTocsyF2NewResonances(argServer,peakList=None,diagTolerance = 0.5,waterMinPpm = 4.88,waterMaxPpm = 4.94):

  assert argServer or peakList
  
  if argServer:
    diagTolerance = argServer.askFloat('Diagonal exclusion tolerance', 0.5)
    waterMinPpm   = argServer.askFloat('Minimum water exclusion ppm ', 4.88) 
    waterMaxPpm   = argServer.askFloat('Maximum water exclusion ppm ', 4.96)
  
  if not peakList:
    project  = argServer.getProject()
    spectra  = getSpectraByType(project,'H[N]_H.TOCSY')
    spectra += getSpectraByType(project,'H_H[N].TOCSY')
    if not spectra:
      argServer.showWarning('Cannot find any 3d H H N spectra')
      return
    
    if len(spectra) > 1:
      argServer.showInfo('Choose 3d TOCSY spectrum')
      spectrum = argServer.getSpectrum(spectra)
    else:
      spectrum = spectra[0]
    peakList = argServer.getPeakList(spectrum)

  spectrum = peakList.dataSource
  
  resonances = []
  for peak in peakList.peaks:
  
    peakDims = peak.sortedPeakDims()
    peakDim0 = peakDims[0]
    peakDim  = peakDims[1]
    ppm      = peakDim.value
    
    if abs( ppm - peakDim0.value) < diagTolerance:
      continue
      
    if (ppm >= waterMinPpm) and (ppm <= waterMaxPpm):
      continue  
     
    if len(peakDim.peakDimContribs) < 1:
      resonance = assignResToDim(peakDim).resonance
      resonances.append(resonance)
      contribs0 = list(peakDim0.peakDimContribs)
      if contribs0:
        resonance0 = contribs0[0].resonance
        if resonance0.resonanceGroup:
          addSpinSystemResonance(resonance0.resonanceGroup, resonance)
    else:
      resonance = peakDim.findFirstPeakDimContrib().resonance
      resonances.append(resonance)

  return resonances


def assignPeaksUnambiguous(argServer, peaks=None, shiftRanges=None, 
                           tolerances=None,aliasing=True, findAssigned=False):
  assignPeaksAutomatic(argServer, ambiguous=False, peaks=peaks, 
                       shiftRanges=shiftRanges, tolerances=tolerances, 
                       aliasing=aliasing, findAssigned=findAssigned)
  
def assignPeaksAutomatic(argServer, ambiguous=True, peaks=None, 
                         shiftRanges=None, tolerances=None, aliasing=True, 
                         findAssigned=False):
  
  from ccpnmr.analysis.core.AssignmentAdvanced import possiblePeakAssigments
  from ccpnmr.analysis.core.AssignmentBasic import assignResToDim
  
  if not peaks:
    peaks = argServer.getCurrentPeaks()
  
  for peak in peaks:
    
    resonances = possiblePeakAssigments(peak, shiftRanges=shiftRanges, 
                                        tolerances=tolerances,
                                        aliasing=aliasing, 
                                        findAssigned=findAssigned)
    
    hasMultiple = bool([ll for ll in resonances if len(ll) > 1])
    
    print ambiguous, hasMultiple, resonances, [ll for ll in resonances if len(ll) > 1]
    
    if ambiguous or not hasMultiple:
      for ii,peakDim in enumerate(peak.sortedPeakDims()):
        for resonance in resonances[ii]:
          assignResToDim(peakDim, resonance)
    
    
