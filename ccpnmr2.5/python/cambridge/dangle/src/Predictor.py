"""
============================COPYRIGHT START=============================

Predictor.py: Part of the DANGLE package (release v1.1.1)

DANGLE: Dihedral ANgles from Global Likelihood Estimates 
Copyright (C) 2009 Nicole Cheung, Tim Stevens, Bill Broadhurst (University of Cambridge)

========================================================================


If you make use of the software or any documents in this package, 
please give credit by citing this package, its authors and references
in the literature with the same authors.

We would appreciate hearing of any problems you may encounter, 
but the programs, the documents and any files created by the programs 
are provided WITHOUT ANY WARRANTY and without even the implied warranty of
CORRECTNESS, MERCHANTABILITY or FITNESS FOR A PARTICULAR OR GENERAL USE.

THE RESPONSIBILITY FOR ANY ADVERSE CONSEQUENCES FROM THE USE OF PROGRAMS OR
DOCUMENTS OR ANY FILE OR FILES CREATED BY USE OF THE PROGRAMS OR DOCUMENTS
LIES SOLELY WITH THE USERS OF THE PROGRAMS OR DOCUMENTS OR FILE OR FILES AND
NOT WITH AUTHORS OF THE PROGRAMS OR DOCUMENTS.

You are not permitted to use any pieces of DANGLE in other programs, make
modifications to DANGLE, or make what a lawyer would call a "derived work" 
in any other way without the consent from either author.
 

============================COPYRIGHT END===============================

for further information, please contact the authors:

- Nicole Cheung   : msc51@cam.ac.uk

- Tim Stevens     : tjs23@cam.ac.uk

- Bill Broadhurst : r.w.broadhurst@bioc.cam.ac.uk

========================================================================

If you are using this software for academic purposes, we suggest
quoting the following reference:

===========================REFERENCE START==============================

Cheung MS, Maguire ML, Stevens TJ, Broadhurst RW. DANGLE: A Bayesian 
inferential method for predicting protein backbone dihedral angles and
secondary structure. J Magn Reson. 202(2010):223-233.

===========================REFERENCE END================================
"""


import sys, os
from math import log, e, sqrt, radians, degrees, atan2, cos, sin, exp
from os.path import isfile

BG = -1
DEGREE_BIN = 10
SEA_LEVEL = 1.0/65536

class Predictor:
  
  def __init__(self, protein, topMatches, reference, writePgm=True):
    
    self.query      = protein
    self.topMatches = topMatches
    self.reference  = reference
    self.writePgm   = writePgm
    self.gleScores  = []
    self.predictions = []
    
    self.readScattergram()
    
    
    
  def predictPhiPsiFromDatabaseMatches(self, progressBar=None):
    
    filename = os.path.join(self.reference.outDir,'danglePred.txt')
    fopen = open(filename,'w')
    fopen.write('!  ResNum  ResName  NumOfIsland  PhiExp  PhiUpper  PhiLower  PsiExp  PsiUpper  PsiLower  Omega  SS\n')
    fopen.write('! =================================================================================================\n')
    
    query = self.query
    res0  = query.res0
    qSequence = query.sequence
    chemShiftsDict = query.chemShiftsDict
    topMatches = self.topMatches

    cns = self.reference.cns
    if cns:
      filename = os.path.join(self.reference.outDir,'danglePred.tbl')
      fcns = open(filename,'w')
      
    self.gleScores  = {}
    self.predictions = {}
    gleScores = self.gleScores
    predictions = self.predictions
    
    if progressBar:
      progressBar.setText('Generating GLE diagrams')
      progressBar.set(0)
      progressBar.open()
      progressBar.parent.update_idletasks()
  
    for resNum in topMatches:
    
      if progressBar:
        progressBar.increment()
      
       
      resName = qSequence[resNum-res0]
      
      prediction = [resNum, resName, None, None, None,
                    None, None, None, None, None, None]
      
      if (resNum+1-res0 < len(qSequence)):
        nextResName = qSequence[resNum+1-res0]
      else:
        nextResName = None
      
      if (resNum-1-res0 >= 0):
        prevResName = qSequence[resNum-1-res0]
      else:
        prevResName = None
      
      aaClass = 'GEN'
      if (resName == 'G'):
        aaClass = 'GLY'
      elif (resName == 'P'):
        aaClass = 'PRO'
      elif (nextResName == 'P'):
        aaClass = 'PRE'
        
      ss = self.findSS(topMatches[resNum])
      
      phiPsiList = [(match[1], match[2]) for match in topMatches[resNum]]

      queryScat = self.getQueryScatterSmoothed(phiPsiList)
      
      scoreList = self.makeScorogram(resNum, queryScat, len(phiPsiList), aaClass)
      gleScores[resNum] = scoreList
      
      islands = self.islanding(scoreList)
      
      if (len(islands) > 0) and (len(islands) <= self.reference.rejectThresh):
      
        majorIsland = self.findMajorIsland(islands, scoreList)
        
        # cis/trans predition for PRO
        omega = 180   # default
        if (resName == 'P')and(chemShiftsDict.get(resNum)):
          shiftCb = chemShiftsDict[resNum].get('CB')
          shiftCg = chemShiftsDict[resNum].get('CG')
          if (shiftCb is not None)and(shiftCg is not None):
            shiftDiff = abs(shiftCb - shiftCg)
            if (shiftDiff > 9):
              omega = 0
            else:
              if (prevResName in ['W','Y','F','G','C']) or (nextResName in ['W','Y','F','G','C']):
                if shiftDiff >= 8:
                  omega = 0
                  
        # primary island for phi & psi
        phiExp, psiExp = self.findPredictionInIsland(scoreList, majorIsland)
        phiUpper, phiLower, psiUpper, psiLower = self.findUpperLowerOfIsland(majorIsland, phiExp, psiExp)
 
        fopen.write('%4d%3s%3d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8s\n' % (resNum,resName,len(islands),phiExp,phiUpper,phiLower,psiExp,psiUpper,psiLower,omega, ss))
        
        prediction = [resNum, resName, len(islands), phiExp, phiUpper,
                      phiLower, psiExp, psiUpper, psiLower, omega, ss]
        
        if cns:
        
          # phi
          upLimit = abs(phiUpper-phiExp)
          if upLimit > 180:
            upLimit = 360 - upLimit
          lowLimit = abs(phiExp-phiLower)
          if lowLimit > 180:
            lowLimit = 360 - lowLimit
          limit = max(upLimit,lowLimit)
          fcns.write(' ASSIGN  (resid %-4s and name C   ) (resid %-4s and name N   )\n         (resid %-4s and name CA  ) (resid %-4s and name C   )  1.00 %3.2f %3.2f 2\n\n' % (resNum-1,resNum,resNum,resNum,phiExp,limit))
        
          # psi
          upLimit = abs(psiUpper-psiExp)
          if upLimit > 180:
            upLimit = 360 - upLimit
          lowLimit = abs(psiExp-psiLower)
          if lowLimit > 180:
            lowLimit = 360 - lowLimit
          limit = max(upLimit,lowLimit)
          fcns.write(' ASSIGN  (resid %-4s and name N   ) (resid %-4s and name CA  )\n         (resid %-4s and name C   ) (resid %-4s and name N   )  1.00 %3.2f %3.2f 2\n\n' % (resNum,resNum,resNum,resNum+1,psiExp,limit))
          
      else:
        fopen.write('%4d%3s%3d%8s%8s%8s%8s%8s%8s%8s%8s\n' % (resNum,resName,len(islands),'None','None','None','None','None','None','None','None'))
        
      predictions[resNum] = prediction
        
    fopen.close()
    
    if progressBar:
      progressBar.destroy()
    
    if self.reference.cns:
      fcns.close()
      
    return self.predictions
    
    
  def makeScorogram(self, resNum, query, n, aaClass):
    
    gList = self.getGScores(query, n)
    scoreList = self.gScoreToBayesianProb(gList, aaClass)
    
    if (not self.reference.angleOnly) and self.writePgm:
      if self.reference.ppm:
        filename = os.path.join(self.reference.outDir, 'Res_%d.ppm' % resNum)
        self.writePpm(scoreList, filename)
      else:
        filename = os.path.join(self.reference.outDir, 'Res_%d.pgm' % resNum)  # default is pgm
        self.writePgm(scoreList, filename)
    
    return scoreList
  
  
  
  def writePgm(self, scoreList, filename):

    fpgm = open(filename,'w')
    fwrite = fpgm.write
    fwrite('P2\n')
    fwrite('#DANGLE GLE PGM FILE\n')
    fwrite('360 360\n')
    
    maxVal = max(scoreList)
    
    fwrite('65535\n')

    i = 0
    while i < 1295:
      for p in range(10):
        for q in range(36):
          for r in range(10):
            if (scoreList[i+q] == BG):
              fwrite('0\n')
            else:
              frac = scoreList[i+q] / maxVal
              fwrite('%d\n' % int(frac*65535))
              
      i += 36

    fpgm.close()
    
  
  def writePpm(self, scoreList, filename):
    
    fppm = open(filename, 'w')
    fwrite = fppm.write
    fwrite('P3\n')
    fwrite('#DANGLE GLE PPM FILE\n')
    fwrite('360 360\n')
    fwrite('255\n')
    
    maxVal = max(scoreList)
    base   = 1 / 65535.0
    
    i = 0
    while i < 1295:
      for p in range(10):
        for q in range(36):
          for r in range(10):
            if (scoreList[i+q] == BG):
              fwrite('255 255 255\n')
            else:
              frac = scoreList[i+q] / maxVal
              if (frac < base):
                fwrite('255 255 255\n')
              else:
                fwrite('%d 0 0\n' % ((1-i)*255))   
              
      i += 36

    fppm.close()
    
 
      
  def getGScores(self, query, n):
    # Get G-scores by comparing query scatter with each database scattergram
    
    scattergram = self.scattergram
    n = float(n)   # sample size
    gList = []
    gListAppend = gList.append
    
    for y in range(36):
      for x in range(36):

        i = x
        j = 36 - 1 - y
          
        if scattergram[i][j] is None:
          gListAppend(BG)
          continue

        plot = scattergram[i][j]

        s = 0.0

        for p in range(36):
          for q in range(36):

            if query[p][q] == 0:
              continue
          
            obs = query[p][q]
            exp = plot[p][q] * n
            s += obs * log((obs/exp))

        s *= 2
        gListAppend(s)
    
    return gList
    

    
  def gScoreToBayesianProb(self, scoreList, aaClass):
    # Convert G-score to bayesian probability
    
    angleProbDict = self.reference.angleProbLookup[aaClass]

    total = 0.0
    
    for g in scoreList:
      if (g != BG):
        total += exp(-g)

    phi = 0
    psi = 35
    newScoreList = []
    
    for g in scoreList:
      if (g == BG):
        newScoreList.append(BG)
      else:
        new_g = exp(-g) / total
        angProb = angleProbDict[phi][psi]
        newScoreList.append(new_g * angProb)
      
      phi += 1
      if (phi > 35):
        phi = 0
        psi -= 1

    return newScoreList
  
  
  def readScattergram(self):
    
    SCATTERGRAM = os.path.join(self.reference.scatDir, 'Plot_%d_%d.int')

    self.scattergram = {}
    
    for i in range(36):

      self.scattergram[i] = {}
      
      for j in range(36):
          
        filename = SCATTERGRAM % (i,j)
        if not isfile(filename):
          self.scattergram[i][j] = None
          continue 
          
        plot = {}
        
        fopen = open(filename, 'r')
        
        for line in fopen.readlines():
          if (line[0] == '\0')or(line[0] == '#'):
            continue

          array = line.split()
          x     = int(array[0])
          y     = int(array[1])
          num   = float(array[2])
          
          if (not plot.get(x)):
            plot[x] = {}
          plot[x][y] = num

        self.scattergram[i][j] = plot

        fopen.close()
        
       
        
  def getQueryScatterSmoothed(self, phiPsiList):

    query = {}
    
    for i in range(36):
      query[i] = [0] * 36

    for (phi0, psi0) in phiPsiList:

      phi0 = (phi0+180)/10.0
      if (phi0 >= 36):
        phi0 = 0 
      phi_bin = int(phi0)

      psi0 = (psi0+180)/10.0
      if (psi0 >= 36):
        psi0 = 0
      psi_bin = int(psi0)

      if (phi0 >= phi_bin+0.5):
        if (psi0 >= psi_bin+0.5):

          # data point in upper right quarter of the bin
          x = 1 - (phi0-0.5 - phi_bin)
          y = 1 - (psi0-0.5 - psi_bin)
          next_phi_bin = phi_bin+1
          if (next_phi_bin >= 36):
            next_phi_bin = 0
          next_psi_bin = psi_bin+1
          if (next_psi_bin >= 36):
            next_psi_bin = 0
          query[phi_bin][psi_bin]           +=  x*y
          query[phi_bin][next_psi_bin]      += (1-y)*x
          query[next_phi_bin][psi_bin]      += (1-x)*y
          query[next_phi_bin][next_psi_bin] += (1-x)*(1-y)
          
        else:

          # data point in lower right quarter of the bin
          x = 1 - (phi0-0.5 - phi_bin)
          y = psi0+0.5 - psi_bin
          next_phi_bin = phi_bin+1
          if (next_phi_bin >= 36):
            next_phi_bin = 0
          prev_psi_bin = psi_bin-1
          if (prev_psi_bin < 0):
            prev_psi_bin = 35
          query[phi_bin][psi_bin]           +=  x*y
          query[phi_bin][prev_psi_bin]      += (1-y)*x
          query[next_phi_bin][psi_bin]      += (1-x)*y
          query[next_phi_bin][prev_psi_bin] += (1-x)*(1-y)
          
      else:
        if (psi0 >= psi_bin+0.5):

          # data point in upper left quarter of the bin
          x = phi0+0.5 - phi_bin
          y = 1 - (psi0-0.5 - psi_bin)
          prev_phi_bin = phi_bin-1
          if (prev_phi_bin < 0):
            prev_phi_bin = 35
          next_psi_bin = psi_bin+1
          if (next_psi_bin >= 36):
            next_psi_bin = 0
          query[phi_bin][psi_bin]           +=  x*y
          query[phi_bin][next_psi_bin]      += (1-y)*x
          query[prev_phi_bin][psi_bin]      += (1-x)*y
          query[prev_phi_bin][next_psi_bin] += (1-x)*(1-y)
          
        else:

          # data point in lower left quarter of the bin
          x = phi0+0.5 - phi_bin
          y = psi0+0.5 - psi_bin
          prev_phi_bin = phi_bin-1
          if (prev_phi_bin < 0):
            prev_phi_bin = 35
          prev_psi_bin = psi_bin-1
          if (prev_psi_bin < 0):
            prev_psi_bin = 35
          query[phi_bin][psi_bin]           +=  x*y
          query[phi_bin][prev_psi_bin]      += (1-y)*x
          query[prev_phi_bin][psi_bin]      += (1-x)*y
          query[prev_phi_bin][prev_psi_bin] += (1-x)*(1-y)

    return query
  
  
    
  def islanding(self, scoreList):
  
    maxVal = max(scoreList)
    k = 65535.0 / maxVal

    islands  = []
    aboveSea = {}

    # 1. Every bin above sea-level is a separate island by its own
    
    i = 0
    j = 35
    
    for fx in scoreList:

      if (fx != BG)and(int(fx*k) > 0):
        islands.append(set([(i,j),]))
        aboveSea[(i,j)] = True

      i += 1
      
      if (i > 35):
        i = 0
        j -= 1

    if (len(islands) == 0):
      return islands


    # 2. Traverse each island, join when necessary

    for j in range(35,-1,-1):
      if (j == 35):
        jUp = 0
      else:
        jUp = j+1   
      
      if (j == 0):
        jDown = 35
      else:
        jDown = j-1 
      
      for i in range(36):
        
        if not aboveSea.get((i,j)):
          continue

        # 2.1 identify neighbouring islands

        islandCenter = None
        islandUp     = None
        islandDown   = None
        islandLeft   = None
        islandRight  = None
        islandUL     = None
        islandUR     = None
        islandDL     = None
        islandDR     = None
      
        if (i == 35):
          iRight = 0
        else:
          iRight = i+1

        if (i == 0):
          iLeft = 35
        else:
          iLeft = i-1
          
        up    = (i,jUp)
        down  = (i,jDown)
        left  = (iLeft, j)
        right = (iRight,j)
        upL   = (iLeft, jUp)
        upR   = (iRight,jUp)
        downL = (iLeft, jDown)
        downR = (iRight,jDown)
          
        
        for island in islands:

          # the 9 islands are mutually exclusive
          if (i,j) in island:
            islandCenter = island
          elif up in island:
            islandUp = island
          elif down in island:
            islandDown = island
          elif left in island:
            islandLeft = island
          elif right in island:
            islandRight = island
          elif upL in island:
            islandUL = island
          elif upR in island:
            islandUR = island
          elif downL in island:
            islandDL = island
          elif downR in island:
            islandDR = island
        

        # 2.2 join neighbouring islands
        #     remove individual islands that have been joined
        #     add the joined island back in
    
        biggerIsland = islandCenter
        islands.remove(islandCenter)

        if (islandUp is not None):
          biggerIsland = biggerIsland.union(islandUp)
          islands.remove(islandUp)
        if (islandDown is not None):
          biggerIsland = biggerIsland.union(islandDown)
          islands.remove(islandDown)
        if (islandLeft is not None):
          biggerIsland = biggerIsland.union(islandLeft)
          islands.remove(islandLeft)
        if (islandRight is not None):
          biggerIsland = biggerIsland.union(islandRight)
          islands.remove(islandRight)
        if (islandUL is not None):
          biggerIsland = biggerIsland.union(islandUL)
          islands.remove(islandUL)
        if (islandUR is not None):
          biggerIsland = biggerIsland.union(islandUR)
          islands.remove(islandUR)
        if (islandDL is not None):
          biggerIsland = biggerIsland.union(islandDL)
          islands.remove(islandDL)
        if (islandDR is not None):
          biggerIsland = biggerIsland.union(islandDR)
          islands.remove(islandDR)

        islands.append(biggerIsland)

    return islands


  
  
  def findMajorIsland(self, islands, scoreList):
    
    numOfIslands = len(islands)

    if numOfIslands == 1:
      return islands[0]
    
    height = max(scoreList)
    index = scoreList.index(height)
    maxPhi = index % 36
    maxPsi = 35 - (index / 36)

    for island in islands:
      if (maxPhi,maxPsi) in island:
        return island
        


  def findPredictionInIsland(self, scoreList, island):
    
    newScoreList = self.makeScoreListFromIsland(scoreList, island=island)
    
    phi1D, psi1D = self.collapseInto1D(newScoreList)
    phiCos = 0.0
    phiSin = 0.0
    psiCos = 0.0
    psiSin = 0.0
    
    for i in range(36):
      
      ang = radians(i * 10 + 5 - 180)
      cosAng = cos(ang)
      sinAng = sin(ang)
      
      freq = phi1D[i]
      phiCos += cosAng * freq
      phiSin += sinAng * freq
      
      freq = psi1D[i]
      psiCos += cosAng * freq
      psiSin += sinAng * freq

    phiExp = degrees(atan2(phiSin, phiCos))
    psiExp = degrees(atan2(psiSin, psiCos))
    
    return (phiExp, psiExp)
  
  
  def findUpperLowerOfIsland(self, island, phi_p, psi_p):

    phi = int((phi_p+180) / 10)
    psi = int((psi_p+180) / 10)

    phi_upper = phi
    phi_lower = phi
    psi_upper = psi
    psi_lower = psi

    phiSet = set([phi0 for (phi0,psi0) in island])
    psiSet = set([psi0 for (phi0,psi0) in island])
    
    for i in range(1,19):
      phi_upper += 1
      if (phi_upper > 35):
         phi_upper = 0
      if (phi_upper not in phiSet):
         break
  
    for i in range(1,18):
      phi_lower -= 1
      if (phi_lower < 0):
         phi_lower = 35
      if (phi_lower not in phiSet):
         break


    for i in range(1,19):
      psi_upper += 1
      if (psi_upper > 35):
         psi_upper = 0
      if (psi_upper not in psiSet):
         break
      
    for i in range(1,18):
      psi_lower -= 1
      if (psi_lower < 0):
         psi_lower = 35
      if (psi_lower not in psiSet):
         break
   
    phi_upper = (phi_upper*10-180+10)
    phi_lower = (phi_lower*10-180)
    psi_upper = (psi_upper*10-180+10)
    psi_lower = (psi_lower*10-180)
    
    if phi_upper == 180:
      phi_upper = 179.99
    if psi_upper == 180:
      psi_upper = 179.99

    return (phi_upper, phi_lower, psi_upper, psi_lower)
    
          
  
  def makeScoreListFromIsland(self, scoreList, island=None):
    # if island is not None, make scoreList using that island only
    # otherwise, make scoreList using all islands above sea-level
    
    newScoreList0 = []

    i = 0
    j = 35
    fx_sum = 0.

    for fx in scoreList:

      if island is not None:
        
        if ((i,j) in island):
          newScoreList0.append(fx)
          fx_sum += fx
        else:
          newScoreList0.append(BG)
          
      else:
        
        if (fx != BG) and (fx >= SEA_LEVEL):
          newScoreList0.append(fx)
          fx_sum += fx
        else:
          newScoreList0.append(BG)

      i += 1
      
      if (i > 35):
        i = 0
        j -= 1

    newScoreList1 = []
    
    for fx in newScoreList0:
      if (fx != BG):
        newScoreList1.append(fx/fx_sum)
      else:
        newScoreList1.append(BG)

    return newScoreList1
  
  
  
  def collapseInto1D(self, scoreList):

    phi1D = [0] * 36
    psi1D = [0] * 36

    phi = 0
    psi = 35

    for fx in scoreList:

      if (fx != BG):
        phi1D[phi] += fx
        psi1D[psi] += fx
        
      phi += 1
  
      if (phi > 35):
        phi = 0
        psi -= 1

    return (phi1D, psi1D)
    
    
  
  def findSS(self, top10Matches):
  
    
    ssList = [match[3] for match in top10Matches]
  
    nTotal = float(len(ssList))
    
    nH = 0
    nE = 0
    nC = 0
    
    for ss in ssList:
      if (ss == 'H'):
        nH += 1
      elif (ss == 'E'):
        nE += 1
      else:
        nC += 1

    if (nH/nTotal >= 0.6):
      return 'H'
    if (nE/nTotal >= 0.6):
      return 'E'
    return 'C'

  
  
  
