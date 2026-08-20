"""
======================COPYRIGHT/LICENSE START==========================

ResonanceIdentification.py: Part of the CcpNmr Clouds program

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
import math,string
from os import system

from ccpnmr.analysis.core.AssignmentBasic import newResonance, assignResToDim, addSpinSystemResonance, \
                                            getResonanceAtomTuple,addPeakResonancesToSpinSystem, clearPeakDim
from ccpnmr.analysis.core.ExperimentBasic import getSpectraByType, findSpectrumDimsByIsotope, getNdSpectra
from ccpnmr.analysis.core.ConstraintBasic import getFixedResonance
from ccpnmr.analysis.core.PeakBasic       import getPeakDimPpm, getPeakHeight, findPeaks, searchPeaks
from ccpnmr.analysis.core.UnitConverter   import pnt2ppm, ppm2pnt
from ccpnmr.clouds.CloudBasic import code1LetterDict

from ccpnmr.analysis.core import ExperimentBasic
from ccpnmr.clouds.PseudoResonances import *

pairs = ( \
  ('HN',  'HA'  ), ('HA',  'HN'  ), ('HN',  'HA*' ), ('HA*', 'HN'  ),
  ('HN',  'HA1' ), ('HA1', 'HN'  ), ('HN',  'HA2' ), ('HA2', 'HN'  ),
  ('HA',  'HB'  ), ('HB',  'HA'  ), ('HA',  'HB1' ), ('HB1', 'HA'  ),
  ('HA',  'HB2' ), ('HB2', 'HA'  ), ('HA',  'HB*' ), ('HB*', 'HA'  ),
  ('HA*', 'HA*' ), ('HB1', 'HB2' ), ('HB2', 'HB1' ), ('HB1', 'HG1' ),
  ('HG1', 'HB1' ), ('HB1', 'HG2' ), ('HG2', 'HB1' ), ('HB2', 'HG1' ),
  ('HG1', 'HB2' ), ('HB2', 'HG2' ), ('HG2', 'HB2' ), ('HB',  'HG*' ),
  ('HG*', 'HB'  ), ('HB1', 'HG*' ), ('HG*', 'HB1' ), ('HB2', 'HG*' ),
  ('HG*', 'HB2' ), ('HB*', 'HG*' ), ('HG*', 'HB*' ), ('HB',  'HG1*'),
  ('HG1*','HB'  ), ('HD1*','HG1*'), ('HG1*','HD1*'), ('HD1*','HG11'),
  ('HG11','HD1*'), ('HD1*','HG12'), ('HG12','HD1*'), ('HD1*','HG2*'),
  ('HG2*','HD1*'), ('HB',  'HG11'), ('HG11','HB'  ), ('HB',  'HG12'),
  ('HG12','HB'  ), ('HB',  'HG2*'), ('HG2*','HB'  ), ('HB',  'HG1' ),
  ('HG1', 'HB'  ), ('HB1', 'HG'  ), ('HG',  'HB1' ), ('HB2', 'HG'  ),
  ('HG',  'HB2' ), ('HB*', 'HG'  ), ('HG',  'HB*' ), ('HB*', 'HB*' ),
  ('HB*', 'HB*' ), ('HG*', 'HD*' ), ('HG',  'HD*' ), ('HD*', 'HG'  ),
  ('HG',  'HD1*'), ('HD1*','HG'  ), ('HG',  'HD2*'), ('HD2*','HG'  ),
  ('HG1', 'HD1' ), ('HD1', 'HG1' ), ('HG1', 'HD2' ), ('HD2', 'HG1' ),
  ('HG2', 'HD1' ), ('HD1', 'HG2' ), ('HG2', 'HD2' ), ('HD2', 'HG2' ),
  ('HG1', 'HD*' ), ('HD*', 'HG1' ), ('HG2', 'HD*' ), ('HD*', 'HG2' ),
  ('HB1', 'HG*' ), ('HG*', 'HB1' ), ('HB2', 'HG*' ), ('HG*', 'HB2' ),
  ('HG*', 'HD*' ), ('HD*', 'HG*' ), ('HG1*','HD*' ), ('HD*', 'HG1*'),
  ('HG1*','HB'  ), ('HB',  'HG1*'), ('HG11','HG12'), ('HG12','HG11'),
  ('HG2*','HG11'), ('HG11','HG2*'), ('HG2*','HG12'), ('HG12','HG2*'),
  ('HG1*','HG2*'), ('HG2*','HG1*'), ('HG2*','HD*' ), ('HD*', 'HG2*'),
  ('HD*', 'HD*' ), ('HD*', 'HE*' ), ('HE*', 'HD*' ), ('HD1', 'HE*' ),
  ('HE*', 'HD1' ), ('HD2', 'HE*' ), ('HE*', 'HD2' ), ('HD1', 'HE1' ),
  ('HE1', 'HD1' ), ('HD*', 'HE'  ), ('HE',  'HD*' ), ('HD1', 'HD2' ),
  ('HD2', 'HD1' ), ('HG*', 'HD1' ), ('HD1', 'HG*' ), ('HG*', 'HD2' ),
  ('HD2', 'HG*' ), ('HD*', 'HD*' ), ('HD21','HD22'), ('HD22','HD21'),
  ('HD2*','HD2*'), ('HD2*','HD1*'), ('HD1*','HD2*'), ('HE21','HE22'),
  ('HE22','HE21'), ('HE2*','HE2*'), ('HE*', 'HZ'  ), ('HZ',  'HE*' ),
  ('HE*', 'HZ*' ), ('HZ*', 'HE*' ), ('HE3', 'HZ3' ), ('HZ3', 'HE3' ),
  ('HE*', 'HE*' ), ('HZ2', 'HH2' ), ('HH2', 'HZ2' ),
)

# for efficiency convert list to dict
pairsDict = {}
for pair in pairs:
  pairsDict[pair] = None

nameMapping = {  'H':'HN',
               'HBa':'HB1', 'HB3':'HB1', 'HBb':'HB2', 'HGa*':'HG1*', 
               'HGa':'HG1', 'HG3':'HG1', 'HGb':'HG2', 'HDa*':'HD1*', 
               'HDa':'HD1', 'HD3':'HD1', 'HDb':'HD2', 'HGb*':'HG2*',
               'HEa':'HD1', 'HE3':'HD1', 'HEb':'HE2', 'HDb*':'HD2*',
                }

olc = code1LetterDict

stdHerr  = 0.008
stdH1err = 0.015
stdH2err = 0.0293
stdNerr  = 0.0336

amideHrange = (6.0,10.5)
alphaHrange = (2.5, 6.0)
rootTwoPi = 2.506628275

"""
BACUS module
"""

def getBacusAssignments(projectName, shiftList, peakLists, is2dBacus=0):

  project = shiftList.root
  assignFileName = '%sass.out' % projectName
  peakResonances = readAssignmentFile(assignFileName, is2dBacus=is2dBacus)
  peaks = []
  for peakList in peakLists:
    peaks.extend( list(peakList.peaks) ) 

  resonances = []
  for spinSystem in project.currentNmrProject.resonanceGroups:
    hydrogens, heavies = getSpinSystemResonances(spinSystem)
    resonances.extend(hydrogens)

  for (I, r1, r2) in peakResonances:
    
    if not is2dBacus:
      i = int(I) - 1 
    else:
      i = int(I)
      
    i1, i2 = int(r1)-1, int(r2)-1
    rs     = (resonances[i1],resonances[i2])
    peak   = peaks[i]
       
    j = 0
    for peakDim in peak.peakDims:
      if '1H' in peakDim.dataDimRef.expDimRef.isotopeCodes:
        contrib = assignResToDim(peakDim, resonance=rs[j])
        if contrib is None:
          print "Assignment failed for peak %s" % peak
          #print '  ' + peak, peakDim, i, i1, i2
          print '  Position:', [p.value for p in peak.peakDims]
          print '  Resonance 1 index:%s object:%s shifts:%f' % (r1, rs[0], rs[0].findFirstShift().value)
          print '  Resonance 2 index:%s object:%s shifts:%f' % (r2, rs[1], rs[1].findFirstShift().value)
        j += 1
    
    
  return resonances


def findClosestPeak(peaks, resonances):

  shifts   = [r.findFirstShift().value for r in resonances]
  minDelta = 10.0
  bestI    = None
  bestPeak = None
  
  i = 0
  for peak in peaks:
    delta = 0.0
    j = 0
    for peakDim in peak.peakDims:
      if '1H' in peakDim.dataDimRef.expDimRef.isotopeCodes:
        v = peakDim.value - shifts[j]
        delta += v*v
        j += 1

    if delta < minDelta:
      minDelta = delta
      bestI    = i

    i += 1
  
  if bestI is not None:
    bestPeak = peaks.pop(bestI)
  else:
    print "Fail", len(peaks), shifts
    
  return bestPeak

def run3dBacus(fileRoot,shiftList,cNoesy=None,nNoesy=None):

  from ccpnmr.analysis.core.PeakBasic import copyPeakListNew
  from ccpnmr.analysis.core.AssignmentBasic import clearPeakDim

  peakLists = []

  print "Make input"
  makeBacusInput(fileRoot,cNoesy,nNoesy,shiftList)
  
  print "Run BACUS 3d"
  runBacus(fileRoot)
  
  if nNoesy:
    print "Clear 15N NOESY"
    nNoesyListNew = copyPeakListNew(nNoesy, nNoesy.dataSource)
    peakLists.append(nNoesyListNew)
    for peak in nNoesyListNew.peaks:
      for peakDim in peak.peakDims:
        clearPeakDim(peakDim)

  if cNoesy:
    print "Clear 13C NOESY"
    cNoesyListNew = copyPeakListNew(cNoesy, cNoesy.dataSource)
    peakLists.append(cNoesyListNew)
    for peak in cNoesyListNew.peaks:
      for peakDim in peak.peakDims:
        clearPeakDim(peakDim)

  print "Make assignments"
  resonances = getBacusAssignments(fileRoot, shiftList, peakLists)
 

def run2dBacus(projectName, executable, noesy,  shiftList):

  from ccpnmr.analysis.core.PeakBasic import copyPeakListNew
  from ccpnmr.analysis.core.AssignmentBasic import clearPeakDim

  project = noesy.root
  
  noesyCopy = copyPeakListNew(noesy, noesy.dataSource)
  for peak in noesyCopy.peaks:
    for peakDim in peak.peakDims:
      clearPeakDim(peakDim)
  
  peakLists = [noesyCopy,]
  ssFileName = '%sgrp.inp' % projectName
  makeSpinSystInputFile(project, ssFileName, shiftList)
  makeCosyTocyConnections(projectName)
  
  writePeakListFile(noesyCopy, 'noesy_ave.inp')

  system('%s %s' % (executable, projectName))
  system('cp assignment.out %sass.out' % projectName)

  getBacusAssignments(projectName, shiftList, peakLists, is2dBacus=1)


def makeBacusInput(projectName, cNoePeakList, nNoePeakList, shiftList):

  project = shiftList.root

  ssFileName = '%sgrp.inp' % projectName
  makeSpinSystInputFile(project, ssFileName, shiftList)

  makeCosyTocyConnections(projectName)

  freqListFileName  = '%sshi.inp' % projectName
  cosyConnFileName  = '%sclk.inp' % projectName
  tocsyConnFileName = '%stlk.inp' % projectName

  cNoeFileName = '%sccp.inp' % projectName
  writePeakListFile(cNoePeakList, cNoeFileName)
  
  cNoeDistFileName = '%scdi.inp' % projectName
  makeIntensityInputFile(cNoePeakList,cNoeDistFileName)

  nNoeFileName = '%sncp.inp' % projectName
  writePeakListFile(nNoePeakList, nNoeFileName)

  nNoeDistFileName = '%sndi.inp' % projectName
  makeIntensityInputFile(nNoePeakList,nNoeDistFileName)


def runBacus(fileRoot):

  from ccpnmr.c.Bacus import bacus
  bacus(fileRoot)


def readAssignmentFile(fileName, is2dBacus=0):

  peakResonances = []
  
  file = open(fileName)
  
  if is2dBacus:
    line = file.readline()
    while line:
      array = line.split()
      if float(array[-1]) > 0.9:
        p, r1, r2 = array[2:5]
        peakResonances.append( [int(p)-1, r2, r1] )

      line = file.readline()
  
  else:
    line = file.readline()
    while line:
      array = line.split()
      peakResonances.append( array[1:4] )
      line = file.readline()
 
    
  return peakResonances

def makeIntensityInputFile(peakList, fileName, intensityType='volume'):
 
  file = open(fileName, 'w')
 
  if not peakList:
    file.close()
    return
 
  mean = 0.0
  intensities = []
  for peak in peakList.peaks:
    intensity = peak.findFirstPeakIntensity(intensityType=intensityType)
    if intensity:
      intensities.append(intensity.value)
      mean += intensity.value
  
  N = float(len(intensities))
  mean /= N
  
  sum = 0.0
  for x in intensities:
    sum += (x-mean) * (x-mean)

  s = math.sqrt(sum/(N-1))

  k = mean*(3.0**6)

  i = 0
  for x in intensities:
    d = (k/x)**(1/6.0)
    line = '%5d %13.2f \t%2.2f\n' % (i+1, x, d)
    file.write(line)
    i += 1

  file.close()


def writePeakListFile(peakList, fileName, intensityType='height'):

  from memops.gui.DataEntry import askInteger  
  
  file = open(fileName, 'w')
  
  if peakList:

    spectrum = peakList.dataSource
    h1Dim = 0
    h2Dim = 1
    
    xDim = None
    if spectrum.numDim > 2:
      xDim = 2 # just default

      peak = peakList.findFirstPeak()
      i = 0
      for peakDim in peak.peakDims:
        isotopes = peakDim.dataDimRef.expDimRef.isotopeCodes
        if ('13C' in isotopes) or ('15N' in isotopes):
          xDim = i
          break
        i += 1
 
      text = '%s:%s:%d' % (spectrum.experiment.name, spectrum.name, peakList.serial)

      h1Dim = int(askInteger('Dimension mapping','Which dimension for peak list %s is heavy atom bound?' % text,h1Dim+1) or 1) - 1
      h2Dim = 3 - (h1Dim + xDim)

    file.write('%d\n' % len(peakList.peaks))

    i = 0
    for peak in peakList.peaks:
      peakDims = peak.peakDims
      ppms = []
      
      if xDim is not None:
        ppms.append( '%8.4f' % peakDims[xDim].value)
      
      ppms.append( '%8.4f' % peakDims[h2Dim].value)
      ppms.append( '%8.4f' % peakDims[h1Dim].value)

      peakIntensity = peak.findFirstPeakIntensity(intensityType=intensityType)
      if not peakIntensity:
        print "Missing intensity", peak
        continue
 
      line = '%4d %s     %9.3e\n' % (i+1,' '.join(ppms),peakIntensity.value)
      file.write(line)
 
      i += 1
      
  else:
    file.write('0\n')

  file.close()


def getSpinSystemResonances(spinSystem, typed=0):

  hydrogens = []
  heavies = []
  for resonance in spinSystem.resonances:
    if (not typed) or resonance.assignNames:
      if resonance.isotopeCode == '1H':
        hydrogens.append(resonance)
      else:
        heavies.append(resonance)

  return hydrogens, heavies

def makeSpinSystInputFile(project, fileName, shiftList=None):

  if not shiftList:
    shiftList = project.currentNmrProject.findFirstMeasurementList(className='shiftList')

  file = open(fileName, 'w')
  
  i = 0
  for spinSystem in project.currentNmrProject.resonanceGroups:
    tlc = spinSystem.ccpCode or '???'
    hydrogens, heavies = getSpinSystemResonances(spinSystem)
    
    seqCode = 0
    if spinSystem.residue:
      seqCode = spinSystem.residue.seqCode
    
    line = '%2d  %s  %d                      %d\n' % (i+1,olc.get(tlc,'X'),len(hydrogens),seqCode)
    file.write(line)

    for hydrogen in hydrogens:
      shift = hydrogen.findFirstShift(parentList = shiftList)
      ppm = shift.value
      
      if hydrogen.assignNames:
        assignName = hydrogen.assignNames[0]
        name = nameMapping.get(assignName) or assignName
      else:
        name = hydrogen.isotopeCode[-1] + '?'
        
      key = name[1:]
      if key[-1] == '*':
        key = key[:-1]
      elif len(key) > 2 and key[-2] in ('11','12','13','21','22','23'):
        key = key[:-1]
      elif len(key) > 1 and key[-1] in ('1','2','3'):
        key = key[:-1]
      
      name2 = 'C%s' % key
      ppm2  = 0.0
      for heavy in heavies:
        if heavy.assignNames:
          nameH = heavy.assignNames[0]
        else:
          nameH = heavy.isotopeCode[-1] + '?'
        
        if key == nameH:
          name2 = nameH
          shift2 = heavy.findFirstShift(parentList = shiftList)

          if shift2:
            ppm2 = shift2.value
            break

        elif key == nameH[1:]:
          name2 = nameH
          shift2 = heavy.findFirstShift(parentList = shiftList)
          
          if shift2:
            ppm2 = shift2.value
            break
      
      else:
        continue
      
      line = '  %6.3f  %-4s    %7.3f  %-4s\n' % (ppm, name, ppm2, name2)
      file.write(line)
    
    file.write('\n')
    i += 1
 
  file.close()
  

"""
SPAN modules
"""  
  
def makeCosyTocyConnections(name=None):

  if name is None:
    name = raw_input('input file name (5 letters): ')[:5]
  
  resonances = getResonancesFromFile(name)
  
  (cy, ty) = getCosyTocsy(resonances)
  write_hconn(resonances, cy, name+'clk')
  write_hconn(resonances, ty, name+'tlk')
  write_hconn2(resonances, cy, 'cylink')
  write_hconn2(resonances, ty, 'tylink')

  fp1 = open(name+'shi.inp', 'w')
  fp2 = open('shifts.inp', 'w')
  for i in range(len(resonances)):
    (r, hshift, xshift, hid, xid, t) = resonances[i]
    fp1.write('%3d  %6.3f  %7.3f   %-4s %1s %2d      %1s\n' % (i+1, float(hshift), float(xshift), hid, t, int(r), xid[:1]))
    fp2.write('%3d  %6.3f    %-4s %1s %2d \n' % (i+1, float(hshift), hid, t, int(r),))
  fp1.close()


def getResonancesFromFile(name):
 
  reson = []
  fp = open(name+'grp.inp')

  line = fp.readline().strip()

  while (line):
    fields = [x.strip() for x in line.split()]
    if (len(fields) < 3):
      break
    (r, t, n) = fields[:3]

    if (t[:2] in ('AR1', 'AR2')):
      t = 'R'
    elif (t == 'PH1'):
      t = 'F'
    elif (t[:2] in ('TR1', 'TR2')):
      t = 'W'
    elif (t == 'TY1'):
      t = 'Y'
    else:
      t = t[0]

    n = int(n)

    for j in range(n):
      line = fp.readline().strip()
      fields = [x.strip() for x in line.split()]
      (hshift, hid, xshift, xid) = fields[:4]
      reson.append((r, hshift, xshift, hid, xid, t))

    line = fp.readline().strip()
    line = fp.readline().strip()

  fp.close()
  return reson


def getCosyTocsy(reson):

  cy = len(reson) * [0]
  ty = len(reson) * [0]
  for i in range(len(reson)):
    cy[i] = []
    ty[i] = []
  for i in range(len(reson)):
    (r1, hshift1, xshift1, hid1, xid1, t1) = reson[i]
    for j in range(i+1, len(reson)):
      (r2, hshift2, xshift2, hid2, xid2, t2) = reson[j]
      if (r1 == r2):
        if (pairsDict.has_key((hid1, hid2))):
          cy[i].append(j)
          cy[j].append(i)
        else:
          ty[i].append(j)
          ty[j].append(i)
   
  return (cy, ty)


def write_hconn(reson, sy, name):

  fp = open(name+'.inp', 'w')
  for i in range(len(reson)):
    (r, hshift, xshift, hid, xid, t) = reson[i]
    fp.write('%3d  %6.3f  %2d\n' % (i+1, float(hshift), len(sy[i])))
    if (len(sy[i])):
      fp.write(15*' ')
      for j in range(len(sy[i])):
        fp.write('  %3d' % (sy[i][j]+1))
    fp.write('\n     \n')
  fp.close()

def write_hconn2(reson, sy, name):

  fp = open(name+'.inp', 'w')
  for i in range(len(reson)):
    (r, hshift, xshift, hid, xid, t) = reson[i]
    fp.write('%3d   %6.3f %2d\n' % (i+1, float(hshift), len(sy[i])))
    if (len(sy[i])):
      fp.write('               ')
      for j in range(len(sy[i])):
        fp.write('  %3d' % (sy[i][j]+1))
    fp.write('\n     \n')
  fp.close()


"""
CMU SPI modules 
"""

def runSpi(executable, projectName, cosy, tocsy, hsqc, hsqcTocsy, hnha, hnhb, hsqcNoesy):

  project = cosy.root
  
  makeSpiInputPeakList(cosy,      projectName + 'csy.inp')
  makeSpiInputPeakList(tocsy,     projectName + 'tsy.inp')
  makeSpiInputPeakList(hsqc,      projectName + 'hsq.inp')
  makeSpiInputPeakList(hsqcTocsy, projectName + 'hty.inp')
  makeSpiInputPeakList(hnha,      projectName + 'hna.inp')
  makeSpiInputPeakList(hnhb,      projectName + 'hnb.inp')
  makeSpiInputPeakList(hsqcNoesy, projectName + 'hny.inp')
  
  system( '%s %s' % (executable,projectName) )

  readSpiSpinSystems(project, projectName + 'grp.inp')
  
  
def makeSpiInputPeakList(peakList, fileName):

  N = peakList.dataSource.numDim
  ohmegaString = ' '.join( ['w%d' % i+1 for i in range(N)] )
  
  file = open(fileName, 'w')
  file.write(ohmegaString)
  file.write('\n\n')

  if N == 2:
    for peak in peakList.peaks:
      file.write( ' '.join( [pd.value for pd in peak.peakDims] ) )
 
  elif N == 3:
    for peak in peakList.peaks:
      positions = [pd.value for pd in peak.peakDims]
      positions[1:] = [posirions[2], positions[1]]
      file.write( ' '.join(positions) )

  file.close()


def readSpiSpinSystems(project, fileName):

  spinSystems = []
  shiftList = project.newShiftList()
  file = open(fileName, 'r')
  
  line = file.readline()
  while line:
    array = line.spit()
    if array and len(array) == 4:
      (n, ident, resType, null) = array
      line   = file.readline()
      shifts = [ float(ppm) for ppm in line.split() ]
      line   = file.readline()
      names  = line.split()
      null   = file.readline()
      
      spinSystem = project.newResonanceGroup()
      
      resonances = []
      for i in range(len(shifts)):
        isotopeCode = '1H'
        name = names[i]
        ppm  = shifts[i]
        if name[0] == 'N':
          isotopeCode = '15N'
        elif name[0] == 'C':
          isotopeCode = '13C'
      
        resonance = project.newResonance(isotopeCode=isotopeCode, name=name)
        shift = shiftList.newShift(value=ppm, resonance=resonance)
        resonances.append(resonance)
    
      spinSystem.setResonances(resonances)
      spinSystems.append(spinSystem)
    
    line = file.readLine()
  
  return spinSystems

"""
Python SPI modules (incomplete)
"""

def spi():

  # setup peak lists
  spectra = getSpectraByType(project,'HSQC')
  if len(spectra) > 1:
    argServer.showInfo('Choose HSQC spectrum')
    hsqc = argServer.getSpectrum(spectra)
  else:
    hsqc = spectra[0]
  hsqcPeaks = list( argServer.getPeakList(hsqc).peaks )

  spectra = getSpectraByType(project,'COSY')
  if len(spectra) > 1:
    argServer.showInfo('Choose COSY spectrum')
    cosy = argServer.getSpectrum(spectra)
  else:
    cosy = spectra[0]
  cosyPeaks = list( argServer.getPeakList(cosy).peaks )

  spectra = getSpectraByType(project,'TOCSY')
  if len(spectra) > 1:
    argServer.showInfo('Choose TOCSY spectrum')
    tocsy = argServer.getSpectrum(spectra)
  else:
    tocsy = spectra[0]
  tocsyPeaks = list( argServer.getPeakList(tocsy).peaks )

  spectra = getSpectraByType(project,'3dTOCSY')
  if len(spectra) > 1:
    argServer.showInfo('Choose HSQC-TOCSY spectrum')
    tocsy3d = argServer.getSpectrum(spectra)
  else:
    tocsy3d = spectra[0]
  tocsy3dPeaks = list( argServer.getPeakList(tocsy3d).peaks )

  spectra = getSpectraByType(project,'HNHA')
  if len(spectra) > 1:
    argServer.showInfo('Choose HNHA spectrum')
    hnha = argServer.getSpectrum(spectra)
  else:
    hnha = spectra[0]
  hnhaPeaks = list( argServer.getPeakList(hnha).peaks )

  spectra = getSpectraByType(project,'3dNOESY')
  if len(spectra) > 1:
    argServer.showInfo('Choose HSQC-NOESY spectrum')
    noey3d = argServer.getSpectrum(spectra)
  else:
    noesy3d = spectra[0]
  noesy3dPeakList = argServer.getPeakList(noesy3d)
  noesy3dPeaks = list(noesy3dPeakList.peaks )
  
  pseudoMolSystem = PseudoMolSystem()

  # find QN side chains
  amideSideChainPseudoSS = identifyAsnGlnAmides(pseudoMolSystem, hsqcPeaks, tocsy3dPeaks)

  # arginine side chains?

  # find TrpD1E1
  trpD1E1PseudoSS = identifyTrpD1E1(pseudoMolSystem, hsqcPeaks, cosyPeaks)

  isSideChain = {}
  for spinSystem in pseudoMolSystem.pseudoSpinSystems:
    for peak in spinSystem.getPeaks():
      isSideChain[peak] = 1
      
  # find hsqc roots
  for peak in hsqcPeaks:
    if isSideChain.get(peak) is not None:
      hsqcPeaks.remove(peak)
      
  backbonePseudoSS = identifyBackboneAmides(pseudoMolSystem, hsqcPeaks)

  # find aromatic cosyTocsy SS
  
  aromaticPseudoSS = identifyAromatic(pseudoMolSystem, cosyPeaks,tocsyPeaks)
  
  # merge TrpD1E1 with aromatic

  trpSideChainSS = mergeTrpSideChain(aromaticPseudoSS, trpD1E1PseudoSS, noesy3dPeakList, cosy)

  # join 3d TOCSY and HNHA, HNHB peaks and resoannces with backbone spin systems
  backbonePseudoSS = identify3dSideChains(pseudoMolSystem, tocsy3dPeaks, hnhaPeaks)
  
  amideAlphaSS = identify2dSideChains(pseudoMolSystem, cosyPeaks,tocsyPeaks)

  glyArgSideChainSS = identifyGlyArgSideChains(amideAlphaSS)

  link2d3d(backbonePseudoSS,amideAlphaSS)

  linkArginine(backbonePseudoSS,glyArgSideChainSS )
  # assignSpinSystems(pseudoMolSystem)

  # march ARG spin systems
  # assumes no folding
  # find arg Nd/Hd 60<15Nshift< 90
  # -> Nd resonance
  # -> Hd resonance
  # Combine with backbone
  # search 3d NOESY at He shift 
  # -> find H res <6.5 - anything in spinSystem
  # -> merge He?Ne with backbone spin system
  
  # Aliphatics
  # find PRO from tocsy/cosy spin systems
  # linking aliphatics into rooted spin systems - with Ha?
  
  # Met
  # Hbs/Hgs 2.13, 2.024 +/- 0.008
  # Read methyls cheating
  
  # Trp ring to backbone
  # Hb = 3.21 +/- 1
  
  # Link NH2 on Gln Asn
  # using NOESY 3d
  #  - finf Hb or Hy
  
  # Aromatic link
  # -> Link to Hbs 
  # Hb 3.01 +/- 0.93
  # try to ding Hd in Noesy
  
  # Hard wirted spectral width

def link2d3d(backboneSS, amideAlphaSS):

  pass
    

def linkAromatic(aromaticPseudoSS):

  for aromSS in trpSideChainSS:
    pass  
  for aromSS in aromaticPseudoSS:
    pass
    
def linkArginine(backboneSS,glyArgSideSS):

  for bb1 in backboneSS:
    amide1 = bb1.getPseudoResonByName('H')
    if amide1 and (amide1.ppm < 93) and (amide1.ppm > 60):
      for bb2 in backboneSS:
        amide2 = bb2.getPseudoResonByName('H')
        if amide2 and (amide2.ppm > 93):
          pass

def linkAliphatics():

  pass

def linkMethionine():

  pass

def linkGlnAsnSide():

  pass

def linkTryptophan():

  pass

def assignSpinSystems(pseudoMolSystem):

  pseudoResons = pseudoMolSystem.pseudoResons
  
  project = pseudoResoons[0].peakDims[0].root
  
  for pseudoSpinSyst in pseudoMolSystem.pseudoSpinSysts:
    resonances = []
    for pseudoReson in pseudoSpinSyst.pseudoResons:    
      resonance = project.newResonance( name=pseudoReson.name)
      resonances.append( resonance )
      
      for peakDim in pseudoReson.peakDims:
        assignResToDim(resonance, peakDim)
        # will make peakDim annotations and set resonance isotope code
      
    spinSystem = project.newResonanceGroup(resonances=resonances, molType='protein', name=pseudoSpinSyst.name, details='Created for CLOUDS')

def identify2dSideChains(pseudoMolSystem, cosyPeaks,tocsyPeaks):
  makePseudoSpinSystFromPeak(pseudoMolSystem, peak, name, resNames=None)

  # setup 2d side chain spin systems from peaks
  cosySS = []
  for peak in cosyPeaks:
    peakDis = peak.sortedPeakDims()
    ppm0 = getPeakDimPpm(peakDims[0])
    ppm1 = getPeakDimPpm(peakDims[1])
    if inRange(ppm0, amideHrange) and  inRange(ppm0, alphaHrange):
      pseudoSS = makePseudoSpinSystFromPeak(pseudoMolSystem, peak, '2dSide', resNames=['H','HA'])
      cosySS.append(pseudoSS)
    elif inRange(ppm1, amideHrange) and  inRange(ppm1, alphaHrange):
      pseudoSS = makePseudoSpinSystFromPeak(pseudoMolSystem, peak, '2dSide', resNames=['HA','H'])
      cosySS.append(pseudoSS)
      
  tocsySS = []
  for peak in tocsyPeaks:
    peakDis = peak.sortedPeakDims()
    ppm0 = getPeakDimPpm(peakDims[0])
    ppm1 = getPeakDimPpm(peakDims[1])
    if inRange(ppm0, amideHrange) and  inRange(ppm0, alphaHrange):
      pseudoSS = makePseudoSpinSystFromPeak(pseudoMolSystem, peak, '2dSide', resNames=['H','HA'])
      tocsySS.append(pseudoSS)
    elif inRange(ppm1, amideHrange) and  inRange(ppm1, alphaHrange):
      pseudoSS = makePseudoSpinSystFromPeak(pseudoMolSystem, peak, '2dSide', resNames=['HA','H'])
      tocsySS.append(pseudoSS)

  (commonSS,orphanedTocsySS) = link2dSpinSysts(cosySS,tocsySS)

  return commonSS

def identifyGlyArgSideChains(spinSystems):
  
  glyArgSS= []
  for ss in spinSystems:
    if len(ss.pseudoResons) > 2:
      amide = ss.getPseudoResonByName('H')
      if amide:
        c = 0
        for pseudoReson in ss.pseudoResons:
          if pseudoReson.name == 'HA':
            c+=1

        if c == 2:
          # an amide with 2x aplhas (or possibly Arg deltas)
          ss.name = 'GlyArgSide'
          glyArgSS.append(ss)
          
  return glyArgSS

def identify3dSideChains(pseudoMolSystem, tocsy3dPeaks, hnhaPeaks, pseudoSpinSysts):

  rootPseudoSS = list(pseudoSpinSysts)
  
  hnhaPseudoSS = []
  for peak in hnhaPeaks:
    pd = peak.sortedPeakDims()
    pseudoResHn   = PseudoReson(pseudoMolSystem,'H',peakDim=pd[0])
    pseudoResHa   = PseudoReson(pseudoMolSystem,'Ha',peakDim=pd[1])
    pseudoResN    = PseudoReson(pseudoMolSystem,'N',peakDim=pd[2])
    hnhaPseudoSS.append( PseudoSpinSyst(pseudoMolSystem,name='3dSide',pseudoResons=[pseudoResHn,pseudoResHa,pseudoResN]) )
    (common,uncommon,prob) = matchPseudoResons([pseudoResHn],[pseudoResHa])
    if common:
      # diagonal
      pseudoMolSystem.mergePseudoResons(pseudoResHn,pseudoResHa)
 
  # TBD repeat for hnhb if available
 
  tocsyPseudoSS = []
  for peak in tocsy3dPeaks:
    # make a pseudo spin system
    # make sure that we're not dealing with Asn Gln side chains
    pd = peak.sortedPeakDims()
    pseudoResHn   = PseudoReson(pseudoMolSystem,'H',peakDim=pd[0])
    pseudoResHx   = PseudoReson(pseudoMolSystem,'Hx',peakDim=pd[1])
    pseudoResN    = PseudoReson(pseudoMolSystem,'N',peakDim=pd[2])
    tPseudoSS     = spinSystem(name='3dSide',pseudoResons=[pseudoResHn,pseudoResN,pseudoResHx])
    tocsyPseudoSS.append( tPseudoSS )
    
    (common,uncommon,prob) = matchPseudoResons([pseudoResHn],[pseudoResHx])
    if common:
      # diagonal
      pseudoMolSystem.mergePseudoResons(pseudoResHn,pseudoResHx)
    
    for haPseudoSS in hnhaPseudoSS:
      (n,score,prob) = matchPairPseudoResons(tPseudoSS.pseudoResons,haPseudoSS.pseudoResons,[2*stdH1err,2*stdH2err,2*stdNerr])
      if n ==3:
        # if all three dims match
        for i in range(3):
          # all pseudo resonances are the same (order is preserved)
          # hnha spin syst first so these resonance names are used (i.e. Hx -> Ha)
          pseudoMolSystem.mergePseudoResons(haPseudoSS.pseudoResons[i],tPseudoSS.pseudoResons[i])
        # merge spin systems
        pseudoMolSystem.mergePseudoSpinSysts(tPseudoSS,haPseudoSS,name='3dSideHa')

  for tPseudoSS in tocsyPseudoSS:
    # check for matches to rooted spin systems
    tPseudoSS.possiblePseudoSS = []
    (pseudoResHn,pseudoResN) = tPseudoSS.resonances[0:2]
    for rootPseudoSS in rootPseudoSS:
      # match H-N root to 3d spin system
      if rootPseudoSS.name == 'AsnGlnSide':
        # 3 resonances in NH2 systems
        rootPseudoRes = rootPseudoSS.pseudoResons[0:2]
        (n,score,prob) = matchPairrootPseudoRes(rootPseudoRes,[pseudoResHn,pseudoResN],[3*stdH1err,3*stdNerr])
        if n ==2 and prob > 0.3:
          tPseudoSS.possiblePseudoSS.append((prob,rootPseudoSS,rootPseudoRes))

        rootPseudoRes = rootPseudoSS.Resonances[1:3]
        (n,score,prob) = matchPairrootPseudoRes(rootPseudoRes,[pseudoResHn,pseudoResN],[3*stdH1err,3*stdNerr])
        if n ==2 and prob > 0.3:
          tPseudoSS.possiblePseudoSS.append((prob,rootPseudoSS,rootPseudoRes))
      
      else:
        if rootPseudoSS.name == 'TrpSide':
          rootPseudoRes = [ rootPseudoSS.getPseudoResonByName('HE1'), rootPseudoSS.getPseudoResonByName('NE1') ]
        elif rootPseudoSS.name == 'backbone':
          rootPseudoRes = rootPseudoSS.pseudoResons
        else:
          continue
          
        (n,score,prob) = matchPairrootPseudoRes(rootPseudoRes,[pseudoResHn,pseudoResN],[3*stdH1err,3*stdNerr])
        if n ==2 and prob > 0.3:
          # H and N match (also includes Q,N,R,W side chains)
          tPseudoSS.possiblePseudoSS.append((prob,rootPseudoSS,rootPseudoRes))
 
      if tPseudoSS.possiblePseudoSS:
        tPseudoSS.possiblePseudoSS.sort()

  doneRoots = {}
  for tPseudoSS in tocsyPseudoSS:
    if pseudoResHn and len(tPseudoSS.possiblePseudoSS) == 1:
      # merge - tocsy with have only one matching root
      rootPseudoSS  = tPseudoSS.possiblePseudoSS[0][1]
      rootPseudoRes = tPseudoSS.possiblePseudoSS[0][2]
      if doneRoots.get(rootPseudoSS):
        print "Dplicate root match in identify3dSideChains"
        continue      
      pseudoMolSystem.mergePseudoResons(rootPseudoRes[0],tPseudoSS.resonances[0])
      pseudoMolSystem.mergePseudoResons(rootPseudoRes[1],tPseudoSS.resonances[1])
      pseudoMolSystem.mergePseudoSpinSysts(rootPseudoSS,tPseudoSS)
      tocsyPseudoSS.remove(tPseudoSS)
      doneRoots[rootPseudoSS] = 1
    
  # order in terms of most likely 
  for tPseudoSS in tocsyPseudoSS:
    # take the most probable and remove from list
    for (prob,rootPseudoSS,rootPseudoRes) in tPseudoSS.possiblePseudoSS:
      # only consider unused roots
      if doneRoots.get(rootPseudoSS) is None:
        pseudoMolSystem.mergePseudoResons(rootPseudoRes[0],tPseudoSS.resonances[0])
        pseudoMolSystem.mergePseudoResons(rootPseudoRes[1],tPseudoSS.resonances[1])
        pseudoMolSystem.mergePseudoSpinSysts(rootPseudoSS,tPseudoSS)
        doneRoots[rootPseudoSS] = 1
        break
    
def mergeTrpSideChain(aromPseudoSpinSysts,trpPseudoSpinSysts, noes3dPeakList,cosy):

  e = 2 * stdHerr
  for trpPseudoSS in trpPseudoSpinSysts:
    pseudoResHe = trpPseudoSS.getPseudoResonByName('HE1')
    pseudoResNe = trpPseudoSS.getPseudoResonByName('NE1')
    if pseudoResHe and pseudoResNe:
      for aromPseudoSS in aromPseudoSpinSysts:
        if len(aromPseudoSS.pseudoResons) == 4: # BS code!!!!!!!!!!!!!!!
          noesyPeaks = []
          for pseudoRes2 in aromPseudoSS.pseudoResons:
            n = 0
            for peakDim in pseudoRes2.peakDims:
              if peakDim.peak.peakList.dataSource is cosy:
                n +=1
            if n == 1:
              # terminal to spin system
              region = [[pseudoResHe.ppm-e,pseudoResHe.ppm+e],[pseudoRes2.ppm-e,pseudoRes2.ppm+e],[pseudoResNe.ppm-e,pseudoResNe.ppm+e]]
              peaks = searchPeaks([noesy3dPeakList], region)
              if len(peaks) == 1:
                peak = peaks[0]
                noesyPeaks.append(peak)

          if len(noesyPeaks) == 2:
            # found both sides of ring
            pseudoMolSystem.mergePseudoSpinSyst(pseudoMolSystem,trpPseudoSS,aromPseudoSS)
            # TBD assign noesy?
            break

def identifyBackboneAmides(pseudoMolSystem, hsqcPeaks):

  spinSystems  = []
  for peak in hsqcPeaks:
    pseudoSpinSyt = makePseudoSpinSystFromPeak(pseudoMolSystem, peak, 'backbone', resNames=['H','N'])
    spinSystems.append( pseudoSpinSyt )
  
  return spinSystems

def makePseudoSpinSystFromPeak(pseudoMolSystem, peak, name, resNames=None):

  res = []
  i = 0
  for peakDim in peak.peakDims:
    
    if names:
      name = resNames[i]
    else:
      name = ''
      for isotopeCode in peakDim.dataDim.expDimRef.isotopeCodes:
        name += isotopeCode
       
    res.append( PseudoReson(pseudoMolSystem,name,peakDims=[peakDim,]) )
    i += 1

  pseudoSpinSyt = PseudoSpinSyst(pseudoMolSystem,name,pseudoResons=res)
  return pseudoSpinSyt

def identifyAromatic(pseudoMolSystem, cosyPeaks, tocsyPeaks):

  aromHrange = (5.7,8.7)
  aromHrangeExtended = (4.8,8.7)

  # setup pseudo spin systems for cosy 
  aromCosy = []
  aromCosyExtended = []
  for peak in cosyPeaks:
    peakDims = peak.sortedPeakDims()
    ppm0 = getPeakDimPpm(peakDims[0])
    ppm1 = getPeakDimPpm(peakDims[1])
    if inRange(ppm0, aromHrange) and inRange(ppm1, aromHrange):
      pseudoSS = makePseudoSpinSystFromPeak(pseudoMolSystem, peak, 'cosy', resNames=['Harom','Harom'])
      pdAromCosy.append(pseudoSS)
    
    else:
      # in extended set
      if inRange(ppm0, aromHrange) and inRange(ppm1, aromHrangeExtended):
        pseudoSS = makePseudoSpinSystFromPeak(pseudoMolSystem, peak, 'cosy', resNames=['Harom','Harom'])
        aromCosyExtended.append(pseudoSS)
      if inRange(ppm1, aromHrange) and inRange(ppm0, aromHrangeExtended):
        pseudoSS = makePseudoSpinSystFromPeak(pseudoMolSystem, peak, 'cosy', resNames=['Harom','Harom'])
        aromCosyExtended.append(pseudoSS)
  
  # setup pseudo spin systems for tocsy
  aromTocsy = []
  for peak in tocsyPeaks:
    peakDims = peak.sortedPeakDims()
    ppm0 = getPeakDimPpm(peakDims[0])
    ppm1 = getPeakDimPpm(peakDims[1])
    if inRange(ppm0, aromHrange) and inRange(ppm1, aromHrange):
      pseudoSS = makePseudoSpinSystFromPeak(pseudoMolSystem, peak, 'tocsy', resNames=['Harom','Harom'])
      aromTocsy.append(pseudoSS)

    else:
      # in extended set
      if inRange(ppm0, aromHrange) and inRange(ppm1, aromHrangeExtended):
        pseudoSS = makePseudoSpinSystFromPeak(pseudoMolSystem, peak, 'tocsy', resNames=['Harom','Harom'])
        aromTocsyExtended.append(pseudoSS)

      if inRange(ppm1, aromHrange) and inRange(ppm0, aromHrangeExtended):
        pseudoSS = makePseudoSpinSystFromPeak(pseudoMolSystem, peak, 'tocsy', resNames=['Harom','Harom'])
        aromTocsyExtended.append(pseudoSS)
      
  # merge ovelapping tocsy and cosy spin systems
  for cosyPseudoSS in list(aromCosy).extend(aromCosyExtended):
    closePseudoSS = []
    for tocsyPseudoSS in list(aromTocsy).extend(aromTocsyExtended):
      (common,uncommon,prob) = matchPseudoResons(cosyPseudoSS.pseudoResons,tocsyPseudoSS.pseudoResons)
      
      if len(common) == 2:
        closePseudoSS.append([1-prob,tocsyPseudoSS,common])
        
    if closePseudoSS:
      closePseudoSS.sort()
      tocsyPseudoSS = closePseudoSS[0][1]
      common   = closePseudoSS[0][2]
      # remove external links to defunct spin syste,
      if tocsyPseudoSS in aromTocsy:
        aromTocsy.remove(tocsyPseudoSS)
      else:
        aromTocsyExtended.remove(tocsyPseudoSS)

      # spin systems are the same (peakDims transferred)
      pseudoMolSystem.mergePseudoSpinSysts(cosyPseudoSS,tocsyPseudoSS,'cosyTocsy')
      
      # resonances are the same
      combinePseudoResons(common)
      
  # test cosy and tocsy spin systems for connectivity
  (aromCosy,aromTocsy) = linkCosyTocsyPseudoSS(aromCosy,aromTocsy)

  # and again over the extended region
  aromCosy.extend(aromCosyExtended)
  aromTocsy.extend(aromTocsyExtended)
  (aromCosy,aromTocsy) = linkCosyTocsyPseudoSS(aromCosy,aromTocsy)

  allPseudoSS = aromCosy.extend( aromTocsy )
  return allPseudoSS

def link2dSpinSysts(cosyPseudoSpinSysts,tocsyPseudoSpinSysts):
  
  # combines overlapping and diagonally reflected peak derived pseudo spin systems
  
  pseudoMolSystem = cosyPseudoSpinSysts[0].parent   
  for cosyPseudoSS1 in cosyPseudoSpinSysts:
    for cosyPseudoSS2 in cosyPseudoSpinSysts:
      (commonCosy,uncommonCosy,prob) = matchPseudoResons(cosyPseudoSS1.pseudoResons,cosyPseudoSS2.pseudoResons)
      if len(commonCosy) == 2 and prob > 0.66:
        # two matching dims
        # same resonances i.e. spin system peaks are digonal reflections
        cosyPseudoSpinSysts.remove(cosyPseudoSS2)
        pseudoMolSystem.mergePseudoSpinSysts(cosyPseudoSS1,cosyPseudoSS2,'cosy')
        combinePseudoResons(commonCosy)
      
      elif len(commonCosy) == 1 and prob > 0.66:
        # one matching dim
        # two spin systems share a resonance
        for tocsyPseudoSS in tocyPseudoSpinSysts:
          (common,null,prob) = matchPseudoResons(uncommonCosy,tocsyPseudoSS.pseudoResons)
          
          if len(common) == 2 and prob > 0.5:
            # cosy spin systems have a confirming tocsy spin system: same spin system
            combinedPseudoSS = pseudoMolSystem.mergePseudoSpinSysts(cosyPseudoSS1,cosyPseudoSS2,'cosyLinked')
            cosyPseudoSpinSysts.remove(cosyPseudoSS2)
            combinePseudoResons(commonCosy)
            
            # tocsy joins the two, uncommon cosy resonances
            tocsyPseudoSpinSysts.remove(tocsyPseudoSS)
            combinedPseudoSS = pseudoMolSystem.mergePseudoSpinSysts(combinedPseudoSS,tocsyPseudoSS,'cosyTocsyLinked')
            combinePseudoResons(common)
                        
  return (cosyPseudoSpinSysts,tocsyPseudoSpinSysts)

def combinePseudoResons(groups):

  for (pseudoRes1,pseudoRes2) in groups:
    pseudoMolSystem = pseudoRes1.parent
    pseudoMolSystem.mergePseudoResons(pseudoRes1,pseudoRes2)

def matchPpm(ppmA,ppmB):

  diff = abs(ppmA-ppmB)
  if diff < stdHerr:
    return 1
    
  return 0

def matchPairPseudoResons(pseudoRes1s,pseudoRes2s,errors):

  N = len(pseudoRes1s)
  
  n      = 0
  prob   = 1
  score  = 0 
  common = []
  for i in range(N):
    ppmA = pseudoRes1s[i].ppm
    ppmB = pseudoRes2s[i].ppm
    diff = abs(ppmA-ppmB)
    if diff < 2*errors[i]:
      e = diff/errors[i]
      score += e*e
      prob  *= math.exp(-0.5*e*e)
      n +=1
  
  return (n,math.sqrt(score),prob)

def matchPseudoResons(pseudoRes1s,pseudoRes2s):

  n = 0
  prob = 0
  matchDict = {}
  for pseudoRes1 in pseudoRes1s:
    matchDict[pseudoRes1] = [pseudoRes1]
    ppmA = pseudoRes1.ppm
    for pseudoRes2 in pseudoRes2s:
      matchDict[pseudoRes2] = [pseudoRes2]
      ppmB = pseudoRes2.ppm
      diff = abs(ppmA-ppmB)
      if diff < 2*stdHerr:
        n += 1
        prob += math.exp( -0.5 * (diff/stdHerr) * (diff/stdHerr) )
                                    
        matchDict[pseudoRes1].extend(matchDict[pseudoRes2])
        matchDict[pseudoRes2].extend(matchDict[pseudoRes1])

  if n > 0:
    prob /= n

  # remove any duplicate groups
  commonDict = {}
  for group in matchDict.values():
    commonDict[tuple(group)] = 1
  common = commonDict.keys()
 
  uncommon = []
  for pseudoRes1 in pseudoRes1s:
    if len(matchDict[pseudoRes1]) < 2:
      uncommon.append(pseudoRes1)

  for pseudoRes2 in pseudoRes2s:
    if len(matchDict[pseudoRes2]) < 2:
      uncommon.append(pseudoRes2)
 
  return (common, uncommon, prob)  

def identifyAsnGlnAmides(pseudoMolSystem, hsqcPeaks, tocsy3dPeaks):

  spinSystems = []
  
  # tolerances - eventually these will be deteremined by linewidth
  tol1H = 0.04
  tol15N = 0.15
  pairTol15N = 0.05
  
  # setup peak lists
  hPeaks = hsqcPeaks
  tPeaks = tocsy3dPeaks
  tocsy3dPeakList = tPeaks[0].peakList
  tNdim  = tocsy3dPeakList.dataSource.numDim
  
  # get the dimension numbers for heteroatom, direct and indirect proton 
  h15Ndim = findSpectrumDimsByIsotope(hPeaks[0].peakList.dataSource,'15N')[0]
  h1Hdim  = findSpectrumDimsByIsotope(hPeaks[0].peakList.dataSource,'1H')[0]
  t15Ndim = findSpectrumDimsByIsotope(tocsy3dPeakList.dataSource,'15N')[0]
  t1H1dim = findSpectrumDimsByIsotope(tocsy3dPeakList.dataSource,'1H')[0]
  t1H2dim = findSpectrumDimsByIsotope(tocsy3dPeakList.dataSource,'1H')[1]
  
  # select 15N shift range for HSQC peaks
  shifts15N = []
  for peak in hPeaks:
    peak.annotation = ''
    peakDim  = peak.peakDims[h15Ndim]
    shift15N = getPeakDimPpm(peakDim)
    
    # remove peaks with inappropriate 15N shifts
    if shift15N > 120.0:
      hPeaks.remove(peak)
    elif shift15N < 100.0:
      hPeaks.remove(peak)
    else:
      shifts15N.append( (shift15N, peak) )
     
  # sort HSQC with good 15N shifts 
  shifts15N.sort()
  
  # match 15N shifts in HSQC pairs
  # - compare each 15N shift of ech peak with each other in range
  close15NPairs = []
  amideQNPairs  = []
  i  = 0
  for (shift, peak) in shifts15N[:-1]:
    for (shift2, peak2) in shifts15N[i:]:
      if abs(shift2 - shift) < pairTol15N:
        close15NPairs.append( (peak,peak2) )
      else:
        break      
    i += 1

  # for each matching HSQC pair look for tocsy return peaks
  for pair in close15NPairs:
    (p1, p2) = pair
    
    # get 15N shift
    shift15N1 = getPeakDimPpm(p1.peakDims[h15Ndim])
    shift15N2 = getPeakDimPpm(p2.peakDims[h15Ndim])
    shift15N = (shift15N1 + shift15N2) /2

    # get amide 1H shift
    shift1H1 = getPeakDimPpm(p1.peakDims[h1Hdim])
    shift1H2 = getPeakDimPpm(p2.peakDims[h1Hdim])
    
    # define first return peak region
    region1 = tNdim * [None]
    region1[t1H1dim] = (shift1H1-tol1H,shift1H1+tol1H)
    region1[t1H2dim] = (shift1H2-tol1H,shift1H2+tol1H)
    region1[t15Ndim] = (shift15N-tol15N,shift15N+tol15N)
    
    # define second return peak region
    region2 = tNdim * [None]
    region2[t1H1dim] = (shift1H2-tol1H,shift1H2+tol1H)
    region2[t1H2dim] = (shift1H1-tol1H,shift1H1+tol1H)
    region2[t15Ndim] = (shift15N-tol15N,shift15N+tol15N)
    
    # search for return peaks in the regions of tocsy
    returnPeak1 = searchPeaks([tocsy3dPeakList], region1)
    returnPeak2 = searchPeaks([tocsy3dPeakList], region2)
    
    # if both return peaks are found mark HQSC and tocsy peaks as QN amide
    if returnPeak1 and returnPeak2:
      amideQNPairs.append( pair )
      # build initially from HSQC
      pseudoResH1 = PseudoReson(pseudoMolSystem,'Hde',peakDims=[p1.peakDims[h1Hdim] ])
      pseudoResH2 = PseudoReson(pseudoMolSystem,'Hde',peakDims=[p2.peakDims[h1Hdim] ])
      pseudoResN  = PseudoReson(pseudoMolSystem,'Nde',peakDims=[p1.peakDims[h15Ndim]])
      nqPseudoSS  = PseudoSpinSyst(pseudoMolSystem, name='GlnAsnSide',pseudoResons=[pseudoResH1,pseudoResH2,pseudoResN] )
            
      # don't bother deal with all of TOCSY_HSQC later
      #for p3 in returnPeak1:
      #  peakDims = p3.peakDims
      #  pseudoResH1.addPeakDim( peakDims[t1H1dim] )
      #  pseudoResH2.addPeakDim( peakDims[t1H2dim] )
      #  pseudoResN.addPeakDim( peakDims[t15Ndim] )

      #for p4 in returnPeak2:
      #  peakDims = p4.peakDims
      #  pseudoResH1.addPeakDim( peakDims[t1H2dim] )
      #  pseudoResH2.addPeakDim( peakDims[t1H1dim] )
      #  pseudoResN.addPeakDim( peakDims[t15Ndim] )
        
      spinSystems.append(nqPseudoSS)
        
  return spinSystems

def identifyTrpD1E1(pseudoMolSystem, hsqcPeaks, cosyPeaks):

  trpHD1range = (   6.5,  8.0 )
  trpHE1range = (   8.9, 12.0 )
  trpNE1range = ( 127.0,170.0 )

  spinSystems = []
  
  possibleHsqc = []
  for peak in hsqcPeaks:
    (ppm0,ppm1) = getPeakPpms(peak)

    if inRange(ppm0,trpHE1range):
      if inRange(ppm1,trpNE1range):
        possibleHsqc.append(peak)

  possibleCosy = []
  for peak in cosyPeaks:
    (ppm0,ppm1) = getPeakPpms(peak)
    
    if inRange(ppm0,trpHE1range) and inRange(ppm1,trpHD1range):
      possibleCosy.append(peak)
    elif inRange(ppm1,trpHE1range) and inRange(ppm0,trpHD1range):
      possibleCosy.append(peak)

  for hsqcPeak in possibleHsqc:
    (ppm0,ppm1) = getPeakPpms(hsqcPeak)
    
    possibleMatchesA = []
    possibleMatchesB = []
    for cosyPeak in possibleCosy:
      (ppm2,ppm3) = getPeakPpms(cosyPeak)
      if ppm2 < ppm3:
        diff = abs(ppm0-ppm3)
        if diff < 3*stdHerr:
          possibleMatchesA.append([diff,cosyPeak])
      else:
        diff = abs(ppm0-ppm2)
        if diff < 3*stdHerr:
          possibleMatchesB.append([diff,cosyPeak])
          
    possibleMatchesA.sort()
    possibleMatchesB.sort()
    if possibleMatchesA or possibleMatchesB:
      
      peakDims = hsqcPeak.sortedPeakDims()
      pseudoResHe = PseudoReson(pseudoMolSystem,'HE1',peakDims=peakDims[0])
      pseudoResNe = PseudoReson(pseudoMolSystem,'NE1',peakDims=peakDims[1])
      pseudoResHd = PseudoReson(pseudoMolSystem,'HD1')
           
      if possibleMatchesA:
        cosyPeak = possibleMatchesA[0][1]
        pseudoResHd.addPeakDim(cosyPeak.peakDims[0])

      if possibleMatchesB:
        cosyPeak = possibleMatchesB[0][1]
        pseudoResHd.addPeakDim(cosyPeak.peakDims[1])
      
      trpPseudoSS = PseudoSpinSyst(pseudoMolSystem,name='TrpSide',pseudoResons=[pseudoResHe,pseudoResNe,pseudoResHd])
      
      spinSystems.append(trpPseudoSS)

  return spinSystems

def assignPeakResonances(peak, resonances=None):

  if not resonances:
    resonances = []

  i = 0
  for peakDim in peak.peakDims:
    if len(resoances <= i):
      resonances.append(None)

    if resonances[i]:
      assignResToDim(peakDim, resonances[i])
    else:
      if peakDim.peakDimContribs:
        resonance = peakDim.findFirstPeakDimContrib().resonance
      else:
        contrib = assignResToDim(peakDim)
        resonance = contrib.resonance
      resonances[i] = resonance
      
    i +=1
    
  return resonances

def getPeakPpms(peak):

  position = []
  for peakDim in peak.peakDims:
    position.append( peakDim.value )
    
  return position

def inRange(value,valRange):

  if value >= valRange[0]:
    if value <= valRange[1]:
      return 1
  return 0
    
def makeNoeAdcs(resonances, spec, constraintHead, diagExclusion=0.4, water=4.92,
                allowedAtomTypes= ['H','HA','HA1','HA2']):

  hTol = 0.02
  allowedAtomTypes= None
  from ccpnmr.analysis.core.ExperimentBasic import getNoiseEstimate
  print "Make DCL"
  distConstraintList  = constraintHead.newDistanceConstraintList()

  shiftList = spec.experiment.shiftList
  getValue = spec.block_file.getValue
  analysisSpec = spec.analysisSpectrum
  noise = 0.5 * min(analysisSpec.posLevels) * analysisSpec.analysisProject.globalContourScale
  print "Estimated Noise", noise
  
  atomTypes = []
  residues = []
  resonances0 = []
  ppms = []
  for i, resonance in enumerate(resonances):
    
    shift = resonance.findFirstShift(parentList=shiftList)
    if not shift:
      continue
    
    atomType = None
    residue = None
    if resonance.resonanceSet:
      atom = resonance.resonanceSet.findFirstAtomSet().findFirstAtom()
      atomType = atom.name
      residue = atom.residue
    
    atomTypes.append(atomType)
    residues.append(residue)
    resonances0.append(resonance)
    ppms.append(shift.value)
    
  resonances = resonances0
  N = len(resonances)
  dataDims = spec.sortedDataDims()
  dataDimRef0   = ExperimentBasic.getPrimaryDataDimRef(dataDims[0])
  dataDimRef1   = ExperimentBasic.getPrimaryDataDimRef(dataDims[1])
  for i in range(N-1):
    if not allowedAtomTypes or (atomTypes[i] in allowedAtomTypes):
      #if 1:
      
      residue0 = residues[i]
      fixedResonanceI = getFixedResonance(constraintHead,resonances[i])
      ppm0 = ppms[i]
      pos0a = ppm2pnt(ppm0-hTol,dataDimRef0) - 1
      pos0b = ppm2pnt(ppm0+hTol,dataDimRef0) - 1
      pos0c = ppm2pnt(ppm0,dataDimRef0) - 1
      
      if pos0c < 0:
        continue
      if pos0c > dataDims[0].numPoints:
        continue
      
      for j in range(i+1,N):
        if not allowedAtomTypes or (atomTypes[j] in allowedAtomTypes):
          #if 1:
          residue1 = residues[j]
          fixedResonanceJ = getFixedResonance(constraintHead,resonances[j])
          ppm1 = ppms[j]
          
          if abs(ppm1-ppm0) < diagExclusion:
            continue
          
          if residue1 and residue0:
            if (residue1.chain is residue0.chain) \
                and (abs(residue1.seqCode-residue0.seqCode) < 3):
              continue
          
          pos1a = ppm2pnt(ppm1-hTol,dataDimRef1) - 1
          pos1b = ppm2pnt(ppm1+hTol,dataDimRef1) - 1
          pos1c = ppm2pnt(ppm1,dataDimRef1) - 1
          
          if pos1c < 0:
            continue
          if pos1c > dataDims[1].numPoints:
            continue
      
          height1 = getValue( (pos0a,pos1a) )
          height2 = getValue( (pos0a,pos1b) )
          height3 = getValue( (pos0b,pos1a) )
          height4 = getValue( (pos0b,pos1b) )
          height5 = getValue( (pos0c,pos1c) )
          
          for h in (height1, height2, height3 ,height4, height5):
            if abs(h) > noise:
              break
          
          else:
            constraint = distConstraintList.newDistanceConstraint(weight=1, targetValue=75, upperLimit=100, lowerLimit=5.5, error=5.0)
            item = constraint.newDistanceConstraintItem(resonances=[fixedResonanceI,fixedResonanceJ])
        
  return distConstraintList

def getCloudsResonanceList(argServer, hsqcPeakList=None, tocsy3dPeakList=None, noesy2dPeakList=None ):

  #Noesy3dPeakList    = getNamedSpectrum(argServer, '3dNOESY')

  if not hsqcPeakList:
    hsqcPeakList       = getHsqcPeakList(argServer) 

  #carbonHsqcPeakList = getCarbonHsqcPeakList(argServer) 

  if not tocsy3dPeakList:
    tocsy3dPeakList    = getNamedSpectrum(argServer, '3dTOCSY')

  if not noesy2dPeakList:
    noesy2dPeakList    = getNamedSpectrum(argServer, '2dNOESY')

  amides             = findAmideResonances(argServer,hsqcPeakList)
  #methyls            = findMethylResonances(argServer,carbonHsqcPeakList)
  aromatic           = []
  
  intensityFactors = []
  amideDict = {}
  checkDict = {}
  for i in range(len(amides)):
    amideDict[amides[i].serial] = i
    intensityFactors.append(1.0)

  otherNonAromatic = []
  
  pickAssignSpecFrom2dRoot(argServer, rootPeakList=hsqcPeakList, targetPeakList=tocsy3dPeakList)
  #tocsyResonances = assign3dTocsyF2NewResonances(argServer, peakList=tocsy3dPeakList)
  #for resonance in tocsyResonances:
  #  if amideDict.get(resonance.serial) is None:
  #    otherNonAromatic.append(resonance)
  for peak in tocsy3dPeakList.peaks:
    peakDims = peak.sortedPeakDims()
  
    if len(peakDims[0].peakDimContribs) > 0:
      if len(peakDims[1].peakDimContribs) > 0:
        resonance = peakDims[1].peakDimContribs[0].resonance
        if resonance not in otherNonAromatic:
          if amideDict.get(resonance.serial) is None:
            otherNonAromatic.append(resonance)
  
  #noesyPeaks = pickAssignSpecFrom2dRoot(argServer, rootPeakList=hsqcPeakList, targetPeakList=Noesy3dPeakList)
  noesyPeaks = list(noesy2dPeakList.peaks)
  for peak in noesyPeaks:
    peakDims = peak.sortedPeakDims()
  
    if len(peakDims[0].peakDimContribs) > 0:
      resonance0 = peakDims[0].peakDimContribs[0].resonance
      if len(peakDims[1].peakDimContribs) > 0:
        resonance = peakDims[1].peakDimContribs[0].resonance
        checkDict[resonance.serial] = 1
        if resonance not in otherNonAromatic:
          if amideDict.get(resonance.serial) is None:
            otherNonAromatic.append(resonance)
          
  for resonance in amides:
    resonance.name = 'HN'
    if not checkDict.get(resonance.serial):
      print "Discarded resonance", getResonanceAtomTuple(resonance)
      amides.remove(resonance)

  for resonance in otherNonAromatic:
    resonance.name = 'H'
  
  maxVolume = 0
  for peak in hsqcPeakList.peaks:
    peakDims = sortedPeakDims()
    i = amideDict.get(peakDims[0].findFirstPeakDimContrib().resonance.serial)
    if i:
      volume = peak.findFirstPeakIntensity(intensityType='volume').value
      intensityFactors[i] = volume
      if volume > maxVolume:
        maxVolume = volume
 
  for i in range(len(amides)):
    intensityFactors[i] = maxVolume/intensityFactors[i]
  
  resonances = []
  resonances.extend(amides)
  #resonances.extend(methyls)
  #resonances.extend(aromatic)
  resonances.extend(otherNonAromatic)
  
  return (resonances,noesyPeaks,intensityFactors)

def pickAssignSpecFrom2dRoot(argServer=None, rootPeakList=None, targetPeakList=None):
  """
  Only copes with 2d roots at the moment
  """
  
  # errors in ppm
  XError = 0.4
  HError = 0.06
  
  assert argServer or (rootPeakList and targetPeakList)

  if argServer:
    project = argServer.getProject()
  else:
    project = rootPeakList.root
    
  if not rootPeakList:
    spectra2d = getNdSpectra(project, 2)
    rootSpec  = argServer.getSpectrum(spectra2d)
    rootPeakList = argServer.getPeakList(rootSpec)

  rootSpec = rootPeakList.dataSource

  if not targetPeakList:
    spectra2d = getNdSpectra(project, 2)
    rootSpec  = argServer.getSpectrum(spectra2d)
    rootPeakList = argServer.getPeakList(rootSpec)

  HgenDims = findSpectrumDimsByIsotope(rootSpec, '1H')
  rootHDim = HgenDims[0]
  rootXDim = findSpectrumDimsByIsotope(rootSpec, '15N')[0]
  XIsotope = '15N'
  if not rootXDim:
    rootXDim   = findSpectrumDimsByIsotope(rootSpec, '13C')[0]
    XIsotope = '13C'
    
  if len(HgenDims) != (len(rootSpec.dataDims)-1):
    argServer.showWarning('Spectrum must have one and only one X isotope (%s) dimension' % XIsotope)
        
  if not targetPeakList:
    spectra3d = getNdSpectra(project, 3)
    for spec in spectra3d:
      targetXDims = findSpectrumDimsByIsotope(rootSpec, XIsotope)
      if len(targetXDims) != 1:
        spectra3d.remove(spec)
        
    targetSpec  = argServer.getSpectrum(spectra3d)
    targetPeakList = argServer.getPeakList(targetSpec)

  targetSpec = targetPeakList.dataSource

  if not targetSpec:
    argServer.showWarning('No suitable search target spectra')
 
  targetXDim = findSpectrumDimsByIsotope(targetSpec, XIsotope)[0]
  hydroDims = findSpectrumDimsByIsotope(targetSpec, '1H')
  targetHDim  = hydroDims[0]
  # Assumes the root Hgen dimension is the first 1H dim
  targetH2Dim = hydroDims[1]
  origPoints = []
  for dataDim in targetSpec.dataDims:
    origPoints.append(dataDim.numPoints)

  #print "ROOT SPEC", rootSpec.experiment.name, rootSpec.name, rootHDim, rootXDim
  #print "TARGET SPEC", targetSpec.experiment.name, targetSpec.name, targetHDim, targetXDim
  N = targetSpec.numDim

  M = len (rootPeakList.peaks)
  c = 0
  
  dataDimRef = ExperimentBasic.getPrimaryDataDimRef(targetSpec.dataDims[targetH2Dim])
  H2PositionMax =  pnt2ppm(1,dataDimRef)
  H2PositionMin =  pnt2ppm(targetSpec.dataDims[targetH2Dim].numPoints,dataDimRef)
  H2Range = (H2PositionMin, H2PositionMax)
  
  foundPeaks = 0
  for peak in rootPeakList.peaks:
    
    region = N * [None]
    peakDims = peak.sortedPeakDims()
    XPosition = pnt2ppm(peakDims[rootXDim].position,peakDims[rootXDim].dataDimRef)
    HPosition = pnt2ppm(peakDims[rootHDim].position,peakDims[rootHDim].dataDimRef)
    XRange  = (XPosition - XError,XPosition + XError)
    HRange  = (HPosition - HError,HPosition + HError)
    
    region[targetXDim]  = XRange
    region[targetHDim]  = HRange
    region[targetH2Dim] = H2Range
    
    peaks = searchPeaks([targetPeakList,], region)   
    #peaks.extend( findPeaks(targetPeakList, region) )
    for peak2 in peaks:
      for intens in peak.peakIntensities:
        if str(intens.value) == 'inf':
          peaks.remove(peak2)
          peak2.delete()
          break
    foundPeaks += len(peaks)

    c +=1
    
    contribs = list(peakDims[rootHDim].peakDimContribs)
    if contribs:
      resonance = contribs[0].resonance
      for peak2 in peaks:
        peakDim = peak2.peakDims[targetHDim]
        if len(peakDim.peakDimContribs) <1:
          assignResToDim(peakDim, resonance)

    contribs = list(peakDims[rootXDim].peakDimContribs)
    if contribs:
      resonance = contribs[0].resonance
      for peak2 in peaks:
        peakDim = peak2.peakDims[targetXDim]
        if len(peakDim.peakDimContribs) <1:
          assignResToDim(peakDim, resonance)
    
    print 'Spin system %d of %d' % (c,M)
    
  if argServer:
    name = '%s:%s' % (targetPeakList.dataSource.experiment.name,targetPeakList.dataSource.name)
    argServer.showInfo('Picked %d peaks in %s' % (foundPeaks,name) )
    
  return targetPeakList.peaks

def assign3dTocsyF2NewResonances(argServer,peakList=None,diagTolerance = 0.5,waterMinPpm = 4.88,waterMaxPpm = 4.94):

  assert argServer or spectrum
  if argServer:
    diagTolerance = argServer.askFloat('Diagonal exclusion tolerance', 0.5)
    waterMinPpm   = argServer.askFloat('Minumum water exclusion ppm ', 4.88) 
    waterMaxPpm   = argServer.askFloat('Maximum water exclusion ppm ', 4.96)
  
  if not peakList:
    project = argServer.getProject()
    spectra = getSpectraByType(project,'3dTOCSY')
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
    peaKDims = peak.sortedPeakDims()
    peakDim0 = peakDims[0]
    peakDim  = peakDims[1]
    ppm      = getPeakDimPpm(peakDim)
    
    if abs( ppm - getPeakDimPpm(peakDim0)) < diagTolerance:
      continue
      
    if (ppm >= waterMinPpm) and (ppm <= waterMaxPpm):
      continue  
     
    if len(peakDim.peakDimContribs) < 1:
      resonance = assignResToDim(peakDim).resonance
      resonances.append(resonance)
      contribs0 = peakDim0.peakDimContribs
      if contribs0:
        resonance0 = contribs0[0].resonance
        if resonance0.resonanceGroup:
          addSpinSystemResonance(resonance0.resonanceGroup, resonance)
    else:
      resonance = peakDim.peakDimContribs[0].resonance
      resonances.append(resonance)

  return resonances

def getHsqcPeakList(argServer):

  project = argServer.getProject()
  spectra = getSpectraByType(project,'HSQC')
  if not spectra:
    argServer.showWarning('Cannot find any nitrogen HSQC spectra')
    return 
  
  if len(spectra) > 1:
    argServer.showInfo('Pick HSQC spectrum')
    spec = argServer.getSpectrum(spectra)  
  else:
    spec = spectra[0]

  if len(spec.peakLists) > 1:
    argServer.showInfo('Pick HSQC peak list')
    peakList = argServer.getPeakList(spec)
  else:
    peakList = spec.peakLists[0]

  return peakList

def getNamedSpectrum(argServer, specType):

  project = argServer.getProject()
  spectra = getSpectraByType(project,specType)
  if not spectra:
    argServer.showWarning('Cannot find any nitrogen %s spectra' % specType)
    return 
  
  if len(spectra) > 1:
    argServer.showInfo('Pick % spectrum' % specType)
    spec = argServer.getSpectrum(spectra)  
  else:
    spec = spectra[0]

  if len(spec.peakLists) > 1:
    argServer.showInfo('Pick peak list')
    peakList = argServer.getPeakList(spec)
  else:
    peakList = spec.peakLists[0]

  return peakList

def getCarbonHsqcPeakList(argServer):

  project = argServer.getProject()
  spectra = getSpectraByType(project,'CHSQC')
  if not spectra:
    argServer.showWarning('Cannot find any nitrogen Carbon HSQC spectra')
    return 
  
  if len(spectra) > 1:
    argServer.showInfo('Pick Carbon HSQC spectrum')
    spec = argServer.getSpectrum(spectra)  
  else:
    spec = spectra[0]

  if len(spec.peakLists) > 1:
    argServer.showInfo('Pick Carbon HSQC peak list')
    peakList = argServer.getPeakList(spec)
  else:
    peakList = spec.peakLists[0]

  return peakList

def findMethylResonances(argServer=None, peakList=None):

  assert argServer or project

  if not peakList:
    peakList = getCarbonHsqcPeakList(argServer)
    if not peakList:
      return
  
  cSpec = peakList.dataSource
  dimH = findSpectrumDimsByIsotope(cSpec, '1H')[0]
  dimC = findSpectrumDimsByIsotope(cSpec, '13C')[0]
  peaks = []
  for peak in peakList.peaks:
    peakDims = peak.peakDims
    ppmH = getPeakDimPpm(peakDims[dimH])
    ppmC = getPeakDimPpm(peakDims[dimC])
    if ppmH < 1.8:
      if ppmC < 22.0:
        peaks.append([getPeakHeight(peak),peak])
      elif ppmH < 1.2:
        if ppmC < 28.0:
          peaks.append([getPeakHeight(peak),peak])
  
  peaks.sort()
  peaks = [x[1] for x in peaks]
  
  N = getNumMethyls(argServer)
  
  if N < len(peaks):
    peaks = peaks[-N:]
  
  resonances = []
  for peak in peaks:
    peak.annotation = '*'
    peakDim = peak.sortedPeakDims()[0]
    if len(peakDim.peakDimContribs) < 1:
      resonance = assignResToDim(peakDim).resonance
    else:
      resonance = peakDim.findFirstPeakDimContrib().resonance
    resonance.nuclGroupType = 'CH3'
    resonance.name = 'CH3'
    resonances.append(resonance)

  return resonances

def findAmideResonances(argServer=None, peakList=None):
 
  assert argServer or peakList
      
  if not peakList:
    peakList = getHsqcPeakList(argServer)
    if not peakList:
      return
  
  nSpec = peakList.dataSource

  dimH = findSpectrumDimsByIsotope(nSpec, '1H')[0]
  resonances = []
  for peak in peakList.peaks:
    peakDim = peak.peakDims[dimH]
    if len(peakDim.peakDimContribs) < 1:
      resonance = assignResToDim(peakDim).resonance
    else:
      resonance = peakDim.findFirstPeakDimContrib().resonance
    resonance.nuclGroupType = 'HN'
    resonance.name = 'HN'
    resonances.append(resonance)

  return resonances

def constrainSpinSystems(distanceConstraintList):

  constraintHead = distanceConstraintList.nmrConstraintStore
  resonances = []
  resDict = {}
  resDict2 = {}
  project = distanceConstraintList.root
  for resonance in project.currentNmrProject.resonances:
    resDict[resonance.serial] = resonance
    
  for constraint in distanceConstraintList.constraints:
    for item in constraint.items:
      for fixedResonance in item.resonances:
        resonance = resDict.get(fixedResonance.resonanceSerial)
        if resonance is not None:
          resDict[resonance.serial] = None
          resonances.append(resonance)
          resDict2[resonance] = 1

  spinSystemDict = {} 
  for resonance in resonances:
    if resonance.resonanceGroup:
      spinSystemDict[resonance.resonanceGroup] =1
      
  spinSystems = spinSystemDict.keys()
  
  for spinSystem in spinSystems:
    resonances = spinSystem.resonances
    n = len(resonances)
    for i in range(n-1):
      resonanceI = resonances[i]
      if resDict2.get(resonanceI):
        for j in range(i+1,n):
          resonanceJ = resonances[j]
          if resDict2.get(resonanceJ):
            fixedResonance0 = getFixedResonance(constraintHead,resonanceI)
            fixedResonance1 = getFixedResonance(constraintHead,resonanceJ)
            constraint = distanceConstraintList.newDistanceConstraint(weight=1, targetValue=10, upperLimit=11.5, lowerLimit=1.5, error=10)
            item = constraint.newDistanceConstraintItem(resonances=[fixedResonance0,fixedResonance1])

  return distanceConstraintList

def getNumMethyls(argServer, chain=None):

  assert argServer or chain

  if not chain:
    chain = argServer.getChain()

  numMethyls = 0
  for residue in chain.residues:
    chemCompLoc = residue.molResidue.chemCompLoc
    for chemAtomSet in chemCompLoc.chemAtomSets:
      if not chemAtomSet.chemAtomSet:
        if len(chemAtomSet.chemAtoms) == 3:
          if chemAtomSet.findFirstChemAtom().elementSymbol == 'H':
            numMethyls += 1
           
  print 'Number of methyls is %s' % (numMethyls)
  return numMethyls
