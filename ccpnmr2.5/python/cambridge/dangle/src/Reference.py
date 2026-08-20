"""
============================COPYRIGHT START=============================

Reference.py: Part of the DANGLE package (release v1.1.1)

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
from os.path import isfile


code3to1 = {
    'ALA':'A', 'GLY':'G', 'LEU':'L', 'VAL':'V', 'TRP':'W',
    'CYS':'C', 'ILE':'I', 'MET':'M', 'PHE':'F', 'PRO':'P',
    'ARG':'R', 'ASP':'D', 'GLU':'E', 'HIS':'H', 'LYS':'K',
    'ASN':'N', 'GLN':'Q', 'SER':'S', 'THR':'T', 'TYR':'Y'
    }
    
code1to3 = {
    'A':'ALA', 'G':'GLY', 'L':'LEU', 'V':'VAL', 'W':'TRP',
    'C':'CYS', 'I':'ILE', 'M':'MET', 'F':'PHE', 'P':'PRO',
    'R':'ARG', 'D':'ASP', 'E':'GLU', 'H':'HIS', 'K':'LYS',
    'N':'ASN', 'Q':'GLN', 'S':'SER', 'T':'THR', 'Y':'TYR'
    }


class Reference:
  
  def __init__(self, location):
    
    self.location                  = location
    self.outDir                    = None
    self.cns                       = False
    self.ppm                       = False
    self.homologyLookup            = None
    self.randomCoilShiftLookup     = None
    self.randomCoilShiftCorrLookup = None
    self.angleProbLookup           = {}
    self.database                  = {}
    self.rejectThresh              = 1296   # default: no rejection
    self.angleOnly                 = False
    
    self.winSize      = 5
    self.numMatches   = 10
    self.overhang     = self.winSize / 2
    self.atomList     = ['HA','CA','CB','C','N']   # the order is important
    
    self.shiftDiffLimit = {'N' : 20.781, 
                           'CA': 6.182,
                           'C' : 4.063,
                           'CB': 4.006,
                           'HA': 0.259   }

    self.missingShift   = {'N' : 20.7936,
                           'CA': 6.2001,
                           'C' : 4.0804,
                           'CB': 4.0401,
                           'HA': 0.2601  }
  
    self.readConfigFile()
    self.readKValues()

    
  def readConfigFile(self):
    
    config  = os.path.join(self.location,'config')
    self.scatDir = os.path.join(self.location,'scattergrams')
    
    try:
      fopen = open(config,'r')
    except:
      print 'Error: Cannot open config file: %s' % config
      #sys.exit(0)
      
    for line in fopen.readlines():
      if (line == '\n')or(line[0] == '#'):
        continue
      array = line.split()
      fileName = os.path.join(self.location,array[1])
      
      if (array[0] == 'HOMOLOGY'):
        self.homologyLookup = self.readHomologyTable(fileName)
      elif (array[0] == 'REF_SHIFTS'):
        self.randomCoilShiftLookup = self.readRandomCoilShift(fileName)
      elif (array[0] == 'SHIFT_CORR'):
        self.randomCoilShiftCorrLookup = self.readRandomCoilShiftCorr(fileName)
      elif (array[0] == 'GEN_ANGLE_PROB'):
        self.angleProbLookup['GEN'] = self.readAngleProb(fileName)
      elif (array[0] == 'GLY_ANGLE_PROB'):
        self.angleProbLookup['GLY'] = self.readAngleProb(fileName)
      elif (array[0] == 'PRO_ANGLE_PROB'):
        self.angleProbLookup['PRO'] = self.readAngleProb(fileName)
      elif (array[0] == 'PRE_ANGLE_PROB'):
        self.angleProbLookup['PRE'] = self.readAngleProb(fileName)
      elif (array[0] == 'DATABASE'):
        self.database = self.readDatabase(fileName)
            
    fopen.close()
      

  def readHomologyTable(self, filename):
  
    if not isfile(filename):
      print 'Error: Homology file %s does not exist.' % filename
      #sys.exit(0)
  
    try:
      fopen = open(filename,'r')
    except:
      print 'Error: Cannot open %s.' % filename
      #sys.exit(0)
  
    dict = None
    aaList = []
  
    for line in fopen.readlines():
      if (line == '\n')or(line[0] == '#'):
        continue
      array = line.split()
      if dict is None:
        dict = {}
        for item in array:
          dict[item] = {}
          aaList.append(item)
      else:
        aa1 = array[0]
        for i in range(1,len(array)):
          aa2 = aaList[i-1]
          dict[aa1][aa2] = float(array[i])
      
    fopen.close()
    return dict
  
  
  def readRandomCoilShift(self, filename):
    
    if not isfile(filename):
      print 'Error: Reference random coil shift file %s does not exist.' % filename
      #sys.exit(0)
  
    try:
      fopen = open(filename,'r')
    except:
      print 'Error: Cannot open %s.' % filename
      #sys.exit(0)
      
    dict = {}
    for line in fopen.readlines():
      if (line == '\n')or(line[0] == '#'):
        continue
      
      array = line.split()
      tlc, atom, element, value = array
      
      if (tlc == 'PRO_t'):
        tlc = 'PRO'
      tlc = code3to1.get(tlc)
      
      if tlc is None:
        continue
      
      if dict.get(tlc) is None:
        dict[tlc] = {}
  
      dict[tlc][atom] = float(value)

    fopen.close()
    
    if dict['G'].get('HA2'):
      dict['G']['HA'] = dict['G']['HA2']
  
    return dict
  
  
   
  def readRandomCoilShiftCorr(self, filename):
    
    if not isfile(filename):
      print 'Error: Random coil shift correction file %s does not exist.' % filename
      #sys.exit(0)
  
    try:
      fopen = open(filename,'r')
    except:
      print 'Error: Cannot open %s.' % filename
      #sys.exit(0)
      
    dict = {}
    for line in fopen.readlines():
      if (line == '\n')or(line[0] == '#'):
        continue
      array = line.split()
      
      tlc, atom, element, valA, valB, valC, valD = array
      
      tlc = code3to1.get(tlc)
      if dict.get(tlc) is None:
        dict[tlc] = {}
  
      dict[tlc][atom] = {}
      dict[tlc][atom][-2] = float(valA)
      dict[tlc][atom][-1] = float(valB)
      dict[tlc][atom][1]  = float(valC)
      dict[tlc][atom][2]  = float(valD)

    if dict['G'].get('HA2'):
      dict['G']['HA'] = dict['G']['HA2']
    
    fopen.close()
    return dict



  def readAngleProb(self, filename):
    
    if not isfile(filename):
      print 'Error: Angle probability file %s does not exist.' % filename
      #sys.exit(0)
  
    try:
      fopen = open(filename,'r')
    except:
      print 'Error: Cannot open %s.' % filename
      #sys.exit(0)
    
    dict = {}
    
    for line in fopen.readlines():
      if (line == '\n')or(line[0] == '#'):
        continue
      array = line.split()
      
      i = int(array[0])
      j = int(array[1])
      p = float(array[2])
      
      if dict.get(i) is None:
        dict[i] = {}
      
      dict[i][j] = p
    
    fopen.close()
    return dict
  
  
  def readKValues(self):
    
    data = """
    Homology   0.7390   1.4780   0.7390
    N          0.1596   0.1752   0.1972
    HA        14.6650  17.5390  15.2510
    C          1.1455   1.2051   1.0422
    CA         0.7213   0.9857   0.7178
    CB         0.7624   0.9092   0.6990
    """
    
    dict = {}
    lines = data.split('\n')
    for line in lines:
      array = line.split()
      if array:
        type, k0, k1, k2 = array
      
        if dict.get(type) is None:
          dict[type] = {}
  
        dict[type][-1] = float(k0)
        dict[type][0]  = float(k1)
        dict[type][1]  = float(k2)
        
    self.kValues = dict
    
  

  def getShiftCorrection(self, resName, atomName, offset):
    
    if self.randomCoilShiftCorrLookup.get(resName):
      if self.randomCoilShiftCorrLookup[resName].get(atomName):
        if self.randomCoilShiftCorrLookup[resName][atomName].get(offset):
          return self.randomCoilShiftCorrLookup[resName][atomName][offset]
    
    return 0.0
  
  
  def readDatabase(self, dbFile):
    
    dict = {}
    fopen = open(dbFile,'r')
    for line in fopen.readlines():
      if (line == '\n')or(line[0] == '#'):
        continue
      winSeq = (line[0:10]).replace(' ','')
      if (dict.get(winSeq) is None):
        dict[winSeq] = []
        
      """
      bmrbId = line[10:20]
      bmrbResNum = line[20:30]
      pdbId = line[30:40]
      chainCode = line[40:50]
      pdbResNum = line[50:60]
      """
    
      shifts = []
      index = 60
      for i in range(self.winSize):
        shiftSet = []
        for atomName in self.atomList:   
          shift = (line[index:index+10]).replace(' ','')
          if (shift == 'None'):
            shift = None
          else:
            shift = float(shift)
          index += 10
          shiftSet.append(shift)
        shifts.append(shiftSet)
      
      phi = float((line[index:index+10]))
      index += 10
      psi = float((line[index:index+10]))
      index += 10
      ss  = line[index+4]
      
      dict[winSeq].append([shifts, (phi,psi,ss)])      
    
    fopen.close()
    return dict
  

      
      
