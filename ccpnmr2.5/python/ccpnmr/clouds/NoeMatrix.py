
"""
======================COPYRIGHT/LICENSE START==========================

NoeMatrix.py: Part of the CcpNmr Clouds program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
def printData(data, fileName):

  print "Writing %s" % fileName
  file = open(fileName, 'w')
  file.write( getDataString(data) )
  file.write('\n')
  file.close()


def getDataString(data, level=0):

  t          = type(data)
  spacer     = '.' * level
  dataString = ''
  if t is type(list()):
    i = 0
    for datum in data:
      dataString += '%s %d:\n' % (spacer,i)
      dataString += getDataString(datum, level=level+1)
      i += 1
      
  elif t is type(tuple()):
    i = 0
    for datum in data:
      dataString += '%s %d:\n' % (spacer,i)
      dataString += getDataString(datum, level=level+1)
      i += 1

  elif t is type({}):
    keys = data.keys()
    keys.sort()
    for key in keys:
      dataString += '%s %s:\n' % (spacer,key)
      dataString += getDataString(data[key], level=level+1)
  
  elif hasattr(data, 'className'):
    if data.className == 'Peak':
      pds = data.peakDims
      anno = '-'.join( ['%s' % pd.annotation for pd in pds] )
      ppms = ','.join( ['%f' % pd.value for pd in pds])
      
      if hasattr(data, 'vol'):
        vol = 'Vol:%f' % data.vol
      else:
        vol  = ' '.join( ['%s:%f' % (pi.intensityType,pi.value) for pi in data.peakIntensities])
      
      dataString += '%s %s %s %s %s\n' % (spacer,data,anno,ppms,vol)
    else:
      dataString += '%s %s\n' % (spacer,data)
  else:
    dataString += '%s %s\n' % (spacer,data)
      
  return dataString


def pickNoePeakList(peakList, shiftList=None):

  # pick according to resonances
  # avoid water regions
  # check for symmetric peaks
  
  project = peakList.root
  
  if not shiftList:
    shiftList = project.currentNmrProject.findFirstMeasurementList(className='ShiftList')
    
  if not shiftList:
    print "No shift list"
    return  
  
  peaks = []
  
  
  return peaks


def getNoeMatrixFromPeaks(noePeaks, resonances, noesyPeakLists, mixingTimes, ratioHD=0.9, intensityType='height', analysis=None):
  """This funtion takes a list of NOE peaks and a list or resonances to create an NOE matrix for MIDGE
     NOESY peak lists are passed in to extract diagonal peaks for intensity normalisation.
     The intensity type height or volume may be selected.
     Additionally the NOE matrix is corrected for deuterium present in the sample
  """

  spectrum = noePeaks[0].peakList.dataSource
  
  # find which peak dims are 1H dims - assumes all peaks similar
  hDims = []
  i = 0
  for peakDim in noePeaks[0].sortedPeakDims():
    if '1H' in peakDim.dataDimRef.expDimRef.isotopeCodes:
      hDims.append(i)
  
    i += 1
  
  if len(hDims) != 2:
    print "Abort: Peaks must have two 1H dimensions"
    return
  
  print "Setup on-the-fly attributes"
  for peak in noePeaks:
    peakDims  = peak.sortedPeakDims()
    peak.ppm1 = peakDims[hDims[0]].value
    peak.ppm2 = peakDims[hDims[1]].value
    
    peakIntensity = peak.findFirstPeakIntensity(intensityType=intensityType)
    
    if peakIntensity:
      peak.vol = peakIntensity.value
    else:
      print "Abort: Peak %s missing %s" % (peak, intensityType) 
      return

  if len(spectrum.dataDims) == 2:
    # get a list of diagonal intensities from peak lists at three different mixing times
    diagonalPeaksList = []
    for noesyPeakList in noesyPeakLists:
      diagonalPeaksList.append( getDiagonalPeaks(noesyPeakList) )

    # should really check that the same number of diagonals is in each peak list

    # clean, normalise and symmetrise NOE peaks
 
    #printData(diagonalPeaksList, 'Check_diagonalPeaksList.txt')

    print "Get excitation profile"
    profile  = getExcitationProfile(spectrum)
    #printData(profile, 'Check_profile.txt')

    print "Correct diagonal for excitation"
    diagonal = correctDiagonalExcitation(diagonalPeaksList, profile)
    #printData(diagonal, 'Check_diagonal.txt')
    
    print "Normalise NOE peaks (%d)" % len(noePeaks)
    noePeaks = normaliseNoePeaks(noePeaks, diagonal, mixingTimes)
    #printData(noePeaks, 'Check_normaliseNoePeaks.txt')

    print "Correct NOEs for excitation"
    noePeaks = correctExcitationProfile(profile, noePeaks)
    #printData(noePeaks, 'Check_correctExcitationProfile.txt')

    print "Filter NOE peaks"
    noePeaks = filterPeaks(noePeaks, diagonalThreshold=0.04, removeUnassigned=1, minimumVolume=None, excludeRegions=None)
    #printData(noePeaks, 'Check_filterPeaks.txt')

    print "Generate symmetry weightings"
    weights  = generateSymmetryWeights(spectrum)
    #printData(weights, 'Check_SymmetryWeights.txt')

    print "Symmetrise peaks"
    noePeaks = symmetrisePeaks(noePeaks, weights, checkDuplicates=0, verbose=0)
    #printData(noePeaks, 'Check_symmetrisePeaks.txt')

  else:
    noePeaks = correct3dNoePeaks(noePeaks, analysis)

  print "Generate NOE matrix"
  N = len(resonances)

  # setup resonance names for typing: HN, H, CH3 etc
  setupResonanceNames(resonances)

  #guessImax = float(intensityMax) * math.exp(-25 * tmix/1000)  
  dict   = {}
  matrix = []
  
  i = 0
  for resonance in resonances:
    dict[resonance.serial] = i
    matrix.append( [0] * N )
    i += 1
  
  vMax = 0.0  
  # Generate the actual matrix using resonance indices  
  for peak in noePeaks:
    contribs0 = list(peak.sortedPeakDims()[0].peakDimContribs)
    contribs1 = list(peak.sortedPeakDims()[1].peakDimContribs)

    if contribs0 and contribs1: #peak is assigned in F1 F2
      resonance0 = contribs0[0].resonance
      resonance1 = contribs1[0].resonance
    
      a = dict.get(resonance0.serial)
      b = dict.get(resonance1.serial)
      if a and b:
        
        #if v > guessImax:
        #  continue
        
        v = peak.vol
        if resonance0.name == 'HN':
          v = v/ratioHD
        if resonance1.name == 'HN':
          v = v/ratioHD
       
        if v > vMax:
          vMax = v
        
        val = v/2.0
                  
        matrix[a][b] = val
        matrix[b][a] = val
 
  #for a in range(N):
  #  for b in range(N):
  #    matrix[a][b] /= vMax
 
  #printData(matrix, 'Check_matrix.txt')
  print "Done"
  return matrix


def setupResonanceNames(resonances):
  """This function sets the name attribute of resonances according to the type of atoms to which they are
     assigned. Hopefully in the future this may be set with knowledge of the experiment types only.
  """

  for resonance in resonances:
    if resonance.resonanceSet:
      atoms0 = list(resonance.resonanceSet.findFirstAtomSet().atoms)
      nAtm   = len(atoms0)
      name   = atoms0[0].name
      if (nAtm == 1) and (name == 'H'):
        resonance.name = 'HN'
      elif (nAtm == 3):
        resonance.name = 'CH3'
      elif (nAtm == 2) and (atoms0[0].residue.ccpCode in ('Tyr','Phe','Ptr')) and (name[:2] in ('HE','HD')):
        resonance.name = 'Haro'
      else:
        resonance.name = 'H'


def getDiagonalPeaks(peakList, tolerance=0.04):
  # assumes F1, F2 are both 1H dims

  diagonalPeaks = []
  for peak in peakList.peaks:
    peakDims = peak.sortedPeakDims()
    delta = abs(peakDims[0].value - peakDims[1].value)
    if delta <= tolerance:
      diagonalPeaks.append(peak)
  
  return diagonalPeaks



def correctDiagonalExcitation(diagonalPeaksList, excitationProfile, intensityType='height'):

  diagonalIntensities = []

  for peaks in diagonalPeaksList:
    intensities = []
    
    for peak in peaks:
    
      peakIntensity = peak.findFirstPeakIntensity(intensityType=intensityType)
      if not peakIntensity:
        intensities.append(None)
        continue
      else:
        intensity = peakIntensity.value  
      
      peakDims = peak.sortedPeakDims()
      ppm      = (peakDims[0].value + peakDims[1].value)/2.0
      correct  = None
      
      for ppmA, ppmB, vol in excitationProfile:
        
        if ppm <= ppmA and ppm > ppmB:
          correct = vol
          break
     
      if correct is None: 
        print 'Peak %s ppm value %f not in excitation profile ranges' % (peak, ppm)
        continue

      intensities.append( [peakDims[0].value, peakDims[1].value, intensity/correct] )

    intensities.sort()
    diagonalIntensities.append(intensities)

  return diagonalIntensities


def getExcitationProfile(spectrum, nSteps=100, isWatergate=0, weightingFactor=1.0):

  from ccpnmr.analysis.core.UnitConverter import pnt2ppm
  from math import sin

  nexcitationProfile_max = 1001
  pi = 3.1415926535

  # assumed F1 is direct dimension
  
  dataDim = None
  for dataDimB in spectrum.dataDims:
    if dataDimB.expDim.isAcquisition:
      dataDim = dataDimB
      break
      
  if not dataDim:
    sortList = [(dd.numPoints, dd) for dd in spectrum.dataDims]
    sortList.sort()
    
    aNum, aDD = sortList[0]
    bNum, bDD = sortList[-1]
    
    if aNum == bNum:
      dataDim = spectrum.sortedDataDims()[0]
    else:
      dataDim = bDD
  
  dataDimRef = dataDim.findFirstDataDimRef()
  ppm0      = pnt2ppm(0, dataDimRef)
  ppm1      = pnt2ppm(dataDim.numPointsOrig, dataDimRef)
  stepWidth = (ppm0-ppm1)/float(nSteps)
  baseFreq  = dataDimRef.expDimRef.baseFrequency
  specFreq  = dataDimRef.expDimRef.sf
  
  print  "Generating excitation profile for %s:%s (%f-%f) ppm" % (spectrum.experiment.name,spectrum.name,ppm0,ppm1)
  
  excitationProfile = []
  
  for i in range(nSteps):
    ppm = ppm0 - (i*stepWidth)
    excitationProfile.append( [ppm,ppm-stepWidth,None ] )
  
  if isWatergate and baseFreq is not None:  
    # 5/11/20 Rasmus added check for baseFreq is None
    O1  = (specFreq - baseFreq) * 1000.0
    O1p = O1/baseFreq
    d19 = raw_input("Enter d19 in seconds:")
    
    dist_nextnullHz  = 1.0/(2*d19);
    dist_nextnullppm = dist_nextnullHz/baseFreq
    
    for i in range(nSteps):
      x = pi*(excitationProfile[i][1]-O1p)/dist_nextnullppm
      excitationProfile[i][2]=abs(sin(x))
      
   
  else:
    if not isWatergate:
      print "No Watergate - equal weighting for all peaks"
    for i in range(nSteps):
      excitationProfile[i][2] = weightingFactor

  return excitationProfile



def moment(data):

  from math import sqrt
  threshold = 0.01

  n = len(data)
  if n<=1:
    print  "n must be at least 2 in moment\n"

  n = float(n)

  s = 0.0
  for x in data:
    s += x
    
  ave  = s/n
  adev = 0.0
  var  = 0.0
  skew = 0.0
  curt = 0.0
  ep   = 0.0
  
  for x in data:
    s     = x-ave
    ep   += s
    adev += abs(s)
    
    p     = s*s
    var  += p
    p    *= s
    skew += p
    p    *= s
    curt += p

  adev /= n
  var   = (var-ep*ep/n)/(n-1)
  sdev  = sqrt(var)
  
  if var:
    skew /= n*sdev*sdev*sdev
    curt  = curt/(n*var*var)-3.0

  else:
    print  "No skew/kurtosis when variance = 0 (in moment)\n"

  return (ave, adev, sdev, var, skew, curt)


def getMean(data):

  N   = float(len(data))
  sum = 0.0
  for x in data:
    sum += x

  return sum/N


def normaliseNoePeaks(noePeaks, diagonalIntensities, mixingTimes):

  from math import exp, log
  
  tmix   = [x/1000.0 for x in mixingTimes]
  sumx   = 0.0
  sumxsq = 0.0
  for t in tmix:
    sumx   += t
    sumxsq += t*t
  
  N               = len(diagonalIntensities[0])
  A0              = [None] * N
  conDiag_d       = [None] * N 
  conDiagVol      = [None] * N 
  log_conDiagVol  = [None] * N 
  sumy            = [None] * N 
  sumxy           = [None] * N
  slope           = [None] * N
  intercept       = [None] * N
  nTmix           = len(tmix)

  for i in range(N):
      
    meanDiagPpm0 = getMean( [ diagonalIntensities[j][i][0] for j in range(nTmix)] )
    meanDiagPpm1 = getMean( [ diagonalIntensities[j][i][1] for j in range(nTmix)] )
    intens       = [ diagonalIntensities[j][i][2] for j in range(nTmix)]
    logs         = [ log(x) for x in intens ] 
 
    sum = 0.0
    for j in range(nTmix):
      sum += logs[j]

    sumy[i]           = sum
    conDiag_d[i]      = [meanDiagPpm0,meanDiagPpm1] 
    conDiagVol[i]     = intens
    log_conDiagVol[i] = logs
  
  for i in range(N):
  
    sum = 0.0
    for j in range(nTmix):
      sum += tmix[j]*log_conDiagVol[i][j]
    
    sumxy[i]     = sum
    slope[i]     = (sumx*sumy[i]-nTmix*sumxy[i])/ ((sumx*sumx)-float(nTmix)*sumxsq)
    intercept[i] = (sumy[i]-slope[i]*sumx)/float(nTmix)
    A0[i]        = exp(intercept[i])
   
  rA0_ave,rA0_adev,rA0_sdev,rA0_var,rA0_skew,rA0_curt = moment(A0)

  print  "A0_ave  = %f" % rA0_ave
  print  "A0_sdev = %f" % rA0_sdev

  average_slope = 0
  for x in slope:
    average_slope += x
    
  average_slope /= float(N)

  print  "exponential constant = %f " % average_slope 

  max_cutoff = rA0_ave+rA0_sdev
  min_cutoff = rA0_ave-rA0_sdev
  
  filtave  = 0.0
  nfiltave = 0
  for i in range(N):
    if (A0[i]<max_cutoff) or (A0[i]>min_cutoff):
      filtave  += A0[i]
      nfiltave += 1
       
  if nfiltave>=1:
    filtave /= float(nfiltave)
  else:
    filtave =  rA0_ave

  for peak in noePeaks:
    peak.vol /= filtave 

  # normalised noesy peaks
  return noePeaks


def filterPeaks(peaks, diagonalThreshold=0.04, removeUnassigned=1, minimumVolume=None, excludeRegions=None):

  outPeaks = []
  for peak in peaks:
  
    if removeUnassigned:  
      n = len(peak.peakDims)
      i = 0
      for peakDim in peak.sortedPeakDims():
        if peakDim.peakDimContribs:
          i += 1
 
      if n != i:
        continue
    
    if diagonalThreshold and ( abs(peak.ppm2-peak.ppm1) < diagonalThreshold):
      continue
    
    if excludeRegions:
      exclude = 0
      for minVal, maxVal in excludeRegions:
        if (peak.ppm1 <= maxVal) and (peak.ppm1 > minVal):
          exclude = 1
          break
        
      if exclude:
        continue
    
    if minimumVolume and (peak.vol < minimumVolume):
      continue
    
    outPeaks.append(peak)

  return outPeaks

    
  
def correctExcitationProfile(excitationProfile, noePeaks, verbose=0):

  # assume on-the-fly peak.ppm1, peak.ppm2, peak.vol
  # excitationProfile is a list of upper_bound, lower_bound, multiplier

  exciteThresh = 0.5
  volThresh    = 0.0001

  N = len(noePeaks)
  M = len(excitationProfile)

  for peak in noePeaks:
    for upper, lower, factor in excitationProfile:
    
      if (peak.ppm1<=upper) and (peak.ppm1>=lower):
        if (peak.vol>=volThresh) and (factor>=exciteThresh):
          if verbose:
            print  "volume modified by the excitation profile\n"
            print  "unmodified volume = %f\n" % peak.vol
            print  "weighting factor = %f\n" % factor

          peak.vol = peak.vol/factor

          if verbose:
            print  "modified volume = %f\n" % peak.vol
                        
        break
  
  return noePeaks


def generateSymmetryWeights(spectrum, nSteps=100, isWatergate=0, weightingFactor=1.0, bleachRegions=None):

  from ccpnmr.analysis.core.UnitConverter import pnt2ppm
  from math import sin, exp
  
  Pi = 3.1415926536

  dataDim   = spectrum.sortedDataDims()[0]
  dataDimRef = dataDim.findFirstDataDimRef()
  ppm0      = pnt2ppm(0, dataDimRef)
  ppm1      = pnt2ppm(dataDim.numPointsOrig, dataDimRef)
  stepWidth = (ppm0-ppm1)/float(nSteps)
  baseFreq  = dataDimRef.expDimRef.baseFrequency
  specFreq  = dataDimRef.expDimRef.sf
  
 
  weights = []
  for i in range(nSteps):
    w = []
    p1 = ppm0-(float(i)*stepWidth)
    for j in range(nSteps):
      w.append( (p1,ppm0-(float(j)*stepWidth),1.0) )
  
    weights.append(w)
  
  print "  StepWidth", stepWidth    

  if isWatergate:
    O1  = (specFreq - baseFreq) * 1000.0
    O1p = O1/baseFreq
    d19 = raw_input("Enter d19 in seconds:")

    dist_nextnullHz=1/(2*d19);
    dist_nextnullppm=dist_nextnullHz/baseFreq
    
    for i in range(nSteps):
      wwf = abs( sin(Pi*(symmat_w1[i]-O1p)/dist_nextnullppm) )
      for j in range(nSteps):
        weights[i][j][2] *= wwf

  if bleachRegions:
    for ppm, width in bleachRegions:
      for i in range(nSteps):
        if (abs(weights[i][0][0]+threshold-ppm)<=width) or (abs(weights[i][0][0]-threshold-ppm)<=width):
          
          gaussexp = (weights[i][0][0]-ppm)/width
          gaussexp = (gaussexp * gaussexp)/-2.0
          ppm_weightingfactor = 1-exp(gaussexp)
          
          for j in range(nSteps):
            weights[i][j][2] *= ppm_weightingfactor

  return weights
  

def getPeakResonances(peak):

  assn = []
  for peakDim in peak.sortedPeakDims():
    resonances = []
    for contrib in peakDim.peakDimContribs:
      resonances.append( contrib.resonance )
    resonances.sort()
    assn.append(resonances)

  if assn == [[] for x in peak.peakDims]:
    return []

  return assn


def symmetrisePeaks(peaks, symmetryWeights, checkDuplicates=0, verbose=0):

  # set up some dictionarys for fast processing
  #  - Avoid order N^2 loops

  peakList        = peaks[0].peakList
  peakIndices     = {}
  peakResonances  = {}
  similarAssigned = {}
  
  i = 0
  for peak in peaks:
    resonances0 = getPeakResonances(peak)
    peaks1 = []
    for resonances in resonances0:
      for resonance in resonances:
        for contrib in resonance.peakDimContribs:
          peak1 = contrib.peakDim.peak
          if (peak1.peakList) is peakList and (peak1 is not peak):
            peaks1.append(peak1)

    peakIndices[peak]     = i
    peakResonances[peak]  = resonances0
    similarAssigned[peak] = peaks1


  if checkDuplicates:
    for peak0 in peaks:
      resonances0 = peakResonances[peak0]
      if resonances0:
        for peak1 in similarAssigned[peak0]: 
          resonances1 = peakResonances.get(peak1) # peak might not be in our list
          if resonances1 and (resonances1 == resonances0):
            print  "Warning: Duplicate assignments for peaks %s and %s" % (peak0,peak1)
            print resonances1, resonances0

  outPeaks = {}
  W        = len(symmetryWeights)
  nSym     = 0
  
  print "  Symmetrising", len(peaks)
  for i in range(len(peaks)-1):
    peak0 = peaks[i]
    
    if i and ( i%100 == 0):
      print "    ", i
    
    resonances0 = peakResonances[peak0]
    if not resonances0:
      continue
    
    for peak1 in similarAssigned[peak0]:
      resonances1 = peakResonances.get(peak1) # peak might not be in our list
      if not resonances1:
        continue
      
      resonances1.reverse()
      if resonances0 == resonances1:
        j = peakIndices[peak1]

        if outPeaks.get(peak0):
          if verbose:
            print "Warning: Peak %s already symmetric to another peak" % peak0
          continue
          
        if outPeaks.get(peak1):
          if verbose:
            print "Warning: Peak %s already symmetric to another peak" % peak1
          continue
        
        nSym += 1;
        
        if verbose:
          print "Found %d symmetric peaks" % nSym
          s0 = '.'.join( [ pd.annotation for pd in peak0.peakDims ] )
          s1 = '.'.join( [ pd.annotation for pd in peak1.peakDims ] )
          print "Peaks %s (%s) - %s (%s)" % (peak0, s0, peak1, s1)
        
        
        factorIJ = None
          
        for k in range(W):
          for l in range(W):
            ppm1, ppm2, weight = symmetryWeights[k][l]
            if (peak0.ppm1 <= ppm1) and (peak0.ppm2 <= ppm2):
              factorIJ = weight
              break

          else:
            continue
          break
           
        if not factorIJ:
          print  "Peak %s at %f,%f does not match symmetry weights matrix" % (peak0, peak0.ppm1, peak0.ppm2)  
          continue

        factorJI = None
        
        for k in range(W):
          for l in range(W):
            ppm1, ppm2, weight = symmetryWeights[k][l]
            if (peak1.ppm1 <= ppm1) and (peak1.ppm2 <= ppm2):
              factorJI = weight
              break
              
          else:
            continue
          break
            
        if not factorJI:
          print  "Peak %s at %f,%f does not match symmetry weights matrix" % (peak1, peak1.ppm1, peak1.ppm2)  
          continue
            
        
        # adjust shifts and volumes ccording to weighted averages
        ppm1 = (peak0.ppm1 * factorIJ) + (peak1.ppm1 * factorJI) / (factorIJ+factorJI)
        ppm2 = (peak0.ppm2 * factorIJ) + (peak1.ppm2 * factorJI) / (factorIJ+factorJI)
        vol  = (peak0.vol  * factorIJ) + (peak1.vol  * factorJI) / (factorIJ+factorJI)
        
        peak0.ppm1 = ppm1
        peak0.ppm2 = ppm2
        
        peak1.ppm1 = ppm1
        peak1.ppm2 = ppm2
        
        peak0.vol  = vol
        peak1.vol  = vol
        
        outPeaks[peak0] = 1
        outPeaks[peak1] = 1
        
        break
  
  return outPeaks.keys()
  
def correct3dNoePeaks(noePeaks, analysis):
  #print "Get excitation profile"
  #print "Correct diagonal for excitation"
  #print "Normalise NOE peaks"
  #print "Correct NOEs for excitation"
  #print "Filter NOE peaks"
  #print "Generate symmetry weightings"
  #print "Symmetrise peaks"

  factor = analysis.argumentServer.askFloat('Global Scale Factor', 1.0) or 1.0

  vMax     = 0.0
  outPeaks = []  
  resonancePeaksDict = {}
  symmResonancesDict = {}

  for peak in noePeaks:
    
    if peak.vol > vMax:
      vMax = peak.vol
    
    peakDims = peak.sortedPeakDims()[:2]
    
    if not peakDims[0].peakDimContribs:
      continue
    
    if not peakDims[1].peakDimContribs:
      continue
    
    # only fully 1H assigned
    
    serials = []
    for peakDim in peakDims:
      s = [c.resonance.serial for c in peakDim.peakDimContribs]
      s.sort()
      serials.append( tuple(s) )
      
    if serials[0] == serials[1]:
      # no diag
      continue
    
    resonancePeaksDict[tuple(serials)] = peak
    
    serials.reverse()
    symmResonancesDict[peak] = tuple(serials)
    
    outPeaks.append(peak)

  for peak in outPeaks:
    symmSerials = symmResonancesDict[peak]
    symmPeak    = resonancePeaksDict.get(symmSerials)
    
    if symmPeak:
      vMean        = (symmPeak.vol + peak.vol)/2.0
      peak.vol     = vMean
      symmPeak.vol = vMean
    
  for peak in outPeaks:
    peak.vol = peak.vol / (vMax * factor)

  return outPeaks



