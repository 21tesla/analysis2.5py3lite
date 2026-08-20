"""
CCPN Wrapper to convert CCPN input and run peak2compProcess
"""
def ccpnPeaksToInterval(peakList, iHalfRange=2.5, shiftScale=10, peakShapeModel='G',
                        halfWidthHalfHeight=2.0, intensThreshold=0.1,
                        intensityType='height'):
  """
  peakList: ccp.nmr.PeakList - input 2D HN peak list to make intervals from 
  iHalfRange: Float - Half the inteval range, i.e. initerval is +/- iHalfRange
  shiftScale: Float - Scaling factor for finding shifts relative to interval range
  peakShapeModel: Character - 'G' or 'L' - Gaussian or Lorentzian
  halfWidthHalfHeight: - Float - Peak shape width parameter
  intensThreshold: - Float - Intensity ratio to consider a peak as influencing a shape
  intensityType: - String - 'height' or 'volume'
  """
  
  from gothenburg.prodecomp.generateInterval import peak2compProcess

  # Check peak list is OK
    
  
  spectrum = peakList.dataSource
  
  if spectrum.numDim != 2:
    print 'ccpnPeaksToInterval failed'
    print 'Input peak list was not from a 2-D experiment'
    return []
  
  dataDims = spectrum.sortedDataDims()
  expDimRefs = [dd.expDim.findFirstExpDimRef() for dd in dataDims]
  
  if None in expDimRefs:
    print 'ccpnPeaksToInterval failed'
    print "Input peak list's experiment is missing dimension references"
    return []
    

  if '1H' in expDimRefs[0].isotopeCodes:
    hDim = 0
    nDim = 1
  elif '1H' in expDimRefs[1].isotopeCodes:
    hDim = 1
    nDim = 0
  else:
    print 'ccpnPeaksToInterval failed'
    print 'Input peak list had no acquisition dimension'
    return []

  # MB set manually !
  # hDim = 0
  # nDim = 1

  # Get peak data
    
  peakNums = []
  hDimPoints = []
  intensityVals = []
  hPpms = []
  nPpms = []
  
  for peak in peakList.sortedPeaks():
    intensity = peak.findFirstPeakIntensity(intensityType=intensityType)
  
    if not intensity:
      continue
  
    peakDims = peak.sortedPeakDims()
    peakDimH = peakDims[hDim]
    peakDimN = peakDims[nDim]
    
    peakNums.append(peak.serial)
    hDimPoints.append(peakDimH.position)
    intensityVals.append(intensity.value)
    hPpms.append(peakDimH.value)
    nPpms.append(peakDimN.value)
  
  
  peakData = peakNums, hDimPoints, intensityVals, hPpms, nPpms
  
    
  intervals0 = peak2compProcess(peakData, iHalfRange,
                               shiftScale, peakShapeModel,
                               halfWidthHalfHeight, intensThreshold)
  
  
  pNlst, stI, endI, comp, ppmH, ppmN = intervals0
  
  intervals = []
  for i in range(len(peakNums)):
    intervals.append([i, stI[i], endI[i], comp[i], ppmH[i], ppmN[i]])

  return intervals               
  

