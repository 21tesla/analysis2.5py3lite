"""
============================COPYRIGHT START=============================

Protein.py: Part of the DANGLE package (release v1.1.1)

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


import sys, string
from Reference import Reference


class Protein:

  def __init__(self, ref):

    self.reference  = ref
    
    self.res0                = None
    self.sequence            = None
    self.chain               = None
    self.chemShiftsDict      = {}
    self.secondaryShiftsDict = None
    self.anglesDict          = None
    
    
  def readShiftsFromXml(self, xml):
    
    self.xmlFile = xml
    
    try:
      fopen = open(xml,'r')
    except:
      print 'Error: Cannot open input xml file %s.' % xml
      #sys.exit(0)
      
    lines = fopen.readlines()
    for i in range(len(lines)):
      
      line = lines[i]
      if (line == '\n')or(line[0] == '#'):
        continue
    
      if (line.find('<res_0>') != -1):
        i, self.res0 = self.getString(i, lines, '<res_0>', '</res_0>')
        self.res0 = int(self.res0)
        if (self.res0 < 0):
          print 'Error: Res0 cannot be negative. Please check your chemical shift file.'
          #sys.exit(0)
      elif (line.find('<chain>') != -1):
        i, self.chain = self.getString(i, lines, '<chain>', '</chain>')
      elif (line.find('<seq_1>') != -1):
        i, self.sequence = self.getString(i, lines, '<seq_1>', '</seq_1>')
        self.sequence = ''.join(self.sequence.split())
        self.sequence = self.sequence.upper()
      elif (line.find('<cs_data>') != -1):
        i, shiftsDict = self.getShifts(i+1, lines)
        self.chemShiftsDict = shiftsDict
    
    fopen.close()
    
    self.calculateSecondaryShifts(self.chemShiftsDict)
    
    
  def getString(self, i, lines, startTag, endTag):
    
    line = lines[i]
    startIndex = line.index(startTag) + len(startTag)
    string = line[startIndex:]
    
    endIndex = string.find(endTag)
    while (endIndex == -1):  # not found
      i += 1
      line = lines[i]
      string += line
      endIndex = string.find(endTag)
    
    string = string[:endIndex]    
    string = ''.join(string.split('\n'))
    return (i+1, string)
      
      
  def getShifts(self, i, lines):
    
    dict = {}
    
    line = lines[i]    
    while (line.find('</cs_data>') == -1):   
      if not ((line == '\n')or(line[0] == '#')or(line.find('<!--') != -1)):       
        array = line.split()
        resNum, resName, atomName, shift = array
        if len(resName) == 1:
          resName = self.reference.code1to3[resName]
        if (atomName == 'CO'):
          atomName = 'C'
        if atomName in ['C','N','HA','CA','CB','CG']:  
          if not dict.get(int(resNum)):
            dict[int(resNum)] = {} 
          dict[int(resNum)][atomName] = float(shift)
        if (resName == 'GLY')and(atomName in ['HA2','HA3']):
          if not dict.get(int(resNum)):
            dict[int(resNum)] = {} 
          if not dict[int(resNum)].get('HA'):
            dict[int(resNum)]['HA'] = float(shift)
          else:
            dict[int(resNum)]['HA'] = (float(shift)+dict[int(resNum)]['HA'])/2.0
      
      i += 1
      if (i >= len(lines)):
        break
      line = lines[i]
    
    return (i+1, dict)
  
  
  def calculateSecondaryShifts(self, shiftsDict):
      
    self.secondaryShiftsDict = {}
    ref = self.reference
    
    for resNum in shiftsDict:
      self.secondaryShiftsDict[resNum] = {}
      resName = self.getResidueOneLetterCode(resNum)
      if resName is None:
        continue
      for atomName in shiftsDict[resNum]:
        if (ref.randomCoilShiftLookup.get(resName)):
          randomCoilShift = ref.randomCoilShiftLookup[resName].get(atomName)
          if randomCoilShift is not None:
            correction = 0.0
            for offset in [-2,-1,1,2]:
              resName2 = self.getResidueOneLetterCode(resNum+offset)
              if resName2 is not None:
                offset2 = -1*offset
                correction += ref.getShiftCorrection(resName2, atomName, offset2)

            cs2 = shiftsDict[resNum][atomName] - (randomCoilShift + correction)
            self.secondaryShiftsDict[resNum][atomName] = cs2
        
        
      
  def getResidueOneLetterCode(self, resNum):

    index = resNum-self.res0
    if (index < 0)or(index >= len(self.sequence)):
      return None

    return self.sequence[resNum-self.res0]

