"""
============================COPYRIGHT START=============================

dangle.py: Part of the DANGLE package (release v1.1.1)

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

DANGLE_CITE = """
Please cite the following reference for DANGLE:

Cheung MS, Maguire ML, Stevens TJ, Broadhurst RW. DANGLE: A Bayesian 
inferential method for predicting protein backbone dihedral angles and
secondary structure. J Magn Reson. 202(2010):223-233.
"""

import os, sys, cPickle
from Reference import Reference
from Protein   import Protein
from Predictor import Predictor

OUTDIR = 'DanglePred'

class Dangle:
  
  def __init__(self, location, inputFile, outputDir=None, cns=False, reject=None,
               angleOnly=False, ppm=False, progressBar=None, writePgm=True):
    
    self.input       = inputFile
    self.progressBar = progressBar
    
    print 'DANGLE (version 1.1)'
    print DANGLE_CITE
  
    # 1. read config file for location of reference information  
    self.reference = Reference(os.path.dirname(location))
    self.reference.outDir  = outputDir or OUTDIR
    if not os.path.isdir(self.reference.outDir):
      os.makedirs(self.reference.outDir)
      
    self.reference.cns       = cns
    self.reference.ppm       = ppm
    self.reference.angleOnly = angleOnly
    
    if (reject is not None):
      self.reference.rejectThresh = reject
    
    # 2. read shifts of query protein (input) and calculate secondary shifts
    self.query = Protein(self.reference)
    self.query.readShiftsFromXml(inputFile)
    
    # 3. compare with DB
    print 'STEP1: Shift search'
    self.topMatches = self.compareWithShiftDB()

    # 4. make preditions from scorograms
    print 'STEP2: GLE generation'
    self.predictor   = Predictor(self.query, self.topMatches, self.reference, writePgm)
    self.predictions = self.predictor.predictPhiPsiFromDatabaseMatches(progressBar=self.progressBar)
   
    

  def compareWithShiftDB(self):
    
    
    # Add timings
    
            
    topMatches = {}
    query      = self.query
    ref        = self.reference
    database   = ref.database
    overhang   = ref.overhang
    
    numMatches   = ref.numMatches
    qSequence    = query.sequence
    kValues      = ref.kValues
    qSequenceLen = len(qSequence)
    qRes0        = query.res0
    rAtomList    = [(an, kValues[an], ref.missingShift[an], ref.shiftDiffLimit[an]) for an in ref.atomList]
    
    homologyLookup = ref.homologyLookup
    secondaryShiftsDict = query.secondaryShiftsDict
    offsetRange = [[i, 0] for i in range(-overhang,overhang+1)]
    offsetRange[0][1] = -1  
    offsetRange[-1][1] = 1  
    # i.e. offsetRange = [[0,-1],[1,0],[2,0],[3,0],[4,1]]
    
    if self.progressBar:
      self.progressBar.total = qSequenceLen - 2*overhang
      self.progressBar.setText('Comparing Shifts With Database')
      self.progressBar.set(0)
    
    for i in range(overhang,qSequenceLen-overhang):  
      
      if self.progressBar:
        self.progressBar.increment()
      
      resNum = i + qRes0
      topMatches[resNum] = resTopMatches = []
      thresh = 99999
            
      unknownWin = qSequence[i-overhang:i+overhang+1]
      
      # Could extract pairs
      for knownWin in database:
        
        databaseKnownWin = database[knownWin]        
        hScore = 0.0

        for offset, pos in offsetRange:
          winIndex = offset+overhang
          homologyScore = homologyLookup[knownWin[winIndex]][unknownWin[winIndex]]
          hScore += homologyScore * homologyScore * kValues['Homology'][pos]

          if hScore > thresh:
            break
 
        if hScore > thresh:
          continue        
        
        for shiftSets, (phi, psi, ss) in databaseKnownWin:
          score  = hScore
          aIndex = 0
        
          for atomName, kValuesAtomName, missingShiftAtomName, atomShiftDiffLimit in rAtomList:
      
            for offset, pos in offsetRange:
              iOffset = i+offset
              qSeqI   = qSequence[iOffset]
          
              if (atomName=='CB') and (qSeqI=='G'):
                continue
                
              if (atomName=='N') and (qSeqI=='P'):
                continue
                
              if (atomName=='C') and (iOffset+1 < qSequenceLen) and (qSequence[iOffset+1]=='P'):
                continue
              
              cs2_known = shiftSets[offset+overhang][aIndex]
            
              if cs2_known is None:
                score += missingShiftAtomName * kValuesAtomName[pos]
            
              else:
                ssDict = secondaryShiftsDict.get(resNum+offset)
              
                if ssDict:
                  cs2_unknown = ssDict.get(atomName)
 
                  if cs2_unknown is not None:
                    val = min(abs(cs2_known-cs2_unknown), atomShiftDiffLimit)
                    score += val * val * kValuesAtomName[pos]
 
                  else:
                    score += missingShiftAtomName * kValuesAtomName[pos]
 
                else:
                  score += missingShiftAtomName * kValuesAtomName[pos]
 
              if score > thresh:
                break
           
            aIndex += 1
               
            if score > thresh:
              break
          
          if score < thresh:
            # improve here - len(topMatches) is limited to numMatches
          
            resTopMatches.append((score, phi, psi, ss, knownWin, unknownWin))
            
            if len(resTopMatches) > numMatches:
              resTopMatches.sort()
              resTopMatches.pop()
              thresh = resTopMatches[-1][0]
      
    #if self.progressBar:
    #  self.progressBar.destroy()
                        
    return topMatches


  

if (__name__ == '__main__'):
  
  if (len(sys.argv) < 2)or('-help' in sys.argv):
    print 'Usage     : python dangle.py inputShiftFile [-dir outputDir] [-cns] [-reject numOfIsland] [-angleOnly] [-ppm]'
    print '-cns      : make a cns .tbl file for dihedral angle constraints'
    print '-reject x : make no prediction for GLE with > x islands'
    print '-angleOnly: only the prediction table (no GLE PGM files) will be generated'
    print '-ppm      : write GLE in colored PPM format instead of grayscaled PGM format.'
    sys.exit(0)
  
  location  = os.path.dirname(os.path.abspath(sys.argv[0]))  # directory of dangle.py
  inputFile = sys.argv[1]
  outputDir = None
  cnsFlag   = False
  reject    = None
  angleOnly = False
  ppm       = False
  if ('-dir' in sys.argv):
    index = sys.argv.index('-dir')
    if (index+1 < len(sys.argv)):
      outputDir = sys.argv[index+1]
  if ('-cns' in sys.argv):
    cnsFlag = True
  if ('-reject' in sys.argv):
    index = sys.argv.index('-reject')
    if (index+1 < len(sys.argv)):
      reject = int(sys.argv[index+1])
      if (reject < 1):
        print 'Error: Rejection threshold must be larger than 0.'
        sys.exit(0)
  if ('-angleOnly' in sys.argv):
    angleOnly = True
  if ('-ppm' in sys.argv):
    ppm = True
      
  d = Dangle(location, inputFile, outputDir=outputDir, cns=cnsFlag, reject=reject, angleOnly=angleOnly, ppm=ppm)
  
  
  
