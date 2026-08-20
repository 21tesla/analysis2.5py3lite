from ccpnmr.analysis.core.PeakBasic import pickPeak, setManualPeakIntensity
from ccpnmr.analysis.core.AssignmentBasic import assignResToDim

def makeDecoupledPeakList(argServer):

  peakList = argServer.getPeakList()
  spectrum = peakList.dataSource
  
  peakAssignDict = {}
  
  for peak in peakList.peaks:
    
    assignKey = []
    for peakDim in peak.sortedPeakDims():
      resonances = [c.resonance for c in peakDim.peakDimContribs]
      resonances = frozenset(resonances)
      assignKey.append(resonances)
    
    assignKey = tuple(assignKey)
    
    if assignKey not in peakAssignDict:
      peakAssignDict[assignKey] = []
      
    peakAssignDict[assignKey].append(peak)

  if not peakAssignDict:
    argServer.showWarning('No assignments detected')
    return
  
  newPeakList = spectrum.newPeakList()
  nDim = spectrum.numDim
  
  for assignKey in peakAssignDict:
    peaks = peakAssignDict[assignKey]
    position = [0.0] * nDim
    height = 0.0
    volume = 0.0
    
    for peak in peaks:
      peakDims = peak.sortedPeakDims()
      
      for i, peakDim in enumerate(peakDims):
        position[i] += peakDim.position
 
      heightIntensity = peak.findFirstPeakIntensity(intensityType='height')
      volumeIntensity = peak.findFirstPeakIntensity(intensityType='volume')
      
      if heightIntensity:
        height += heightIntensity.value
      
      if volumeIntensity:
        volume += volumeIntensity.value
      
    n = float(len(peaks))
    for i in range(nDim):
      position[i] /= n
 
    newPeak = pickPeak(newPeakList, position, doFit=False)
    newPeakDims = newPeak.sortedPeakDims()
    
    setManualPeakIntensity(newPeak, height, 'height')
    setManualPeakIntensity(newPeak, volume, 'volume')
    
    for i, peakDim in enumerate(newPeakDims):
      resonances = assignKey[i]
      
      for resonance in resonances:
        assignResToDim(peakDim, resonance, tolerance=1.0, doWarning=False)
  
  argServer.showWarning('Made new peak list with %d peaks' % len(newPeakList.peaks))
  
