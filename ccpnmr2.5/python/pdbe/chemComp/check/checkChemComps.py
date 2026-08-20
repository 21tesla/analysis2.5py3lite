"""

Script that checks bond types, ... of chemComps

"""

import os, re

# General stuff
from memops.universal.Util import drawBox

from memops.api import Implementation

from ccp.general.Util import getOtherAtom
from ccp.general.Io import getChemComp

from pdbe.chemComp.Util import initialiseChemCompScript

def checkChemComp(chemComp,verbose=False):

  if verbose:
    printAtomsBonds(chemComp)

  checkAtomBinding(chemComp,verbose=verbose)
  
def printAtomsBonds(chemComp):

  print drawBox("ChemBond information",indent = "  ")

  for chemBond in chemComp.sortedChemBonds():
    chemAtomNames = ["%s (%d)" % (chemAtom.name,chemAtom.subType) for chemAtom in chemBond.chemAtoms]
    print "    %-12s-%-12s: %s" % (chemAtomNames[0],chemAtomNames[1],chemBond.bondType)
    
  # Print CCV info
  print
  print drawBox("ChemCompVar atom information",indent = "  ")
  for ccv in chemComp.sortedChemCompVars():
  
    print "  %s, %s" % (ccv.linking,ccv.descriptor)
    
    chemAtoms = ccv.sortedChemAtoms()
    print "    %s" % ', '.join(["%s (%d)" % (chemAtom.name,chemAtom.subType) for chemAtom in chemAtoms])
    
    otherChemAtoms = []
    for chemAtom in chemComp.sortedChemAtoms():
      if not chemAtom in chemAtoms:
        otherChemAtoms.append(chemAtom)
        
    print "    NOT INCLUDED: %s" % ', '.join(["%s (%d)" % (chemAtom.name,chemAtom.subType) for chemAtom in otherChemAtoms])
    print
    
  print
  
def checkAtomBinding(chemComp,verbose=False):

  for ccv in chemComp.sortedChemCompVars():
    
    ccvAtoms = ccv.chemAtoms
    ccvBonds = ccv.chemBonds
    
    #
    # Set up reference info
    #
    
    ccvAtomBindings = {}
    ccvAtomNames = []
    ccvAtomDict = {}
    
    for ccvAtom in ccvAtoms:
      ccvAtomBindings[ccvAtom] = []
      ccvAtomNames.append(ccvAtom.name)
      ccvAtomDict[ccvAtom.name] = ccvAtom
      for bond in ccvAtom.chemBonds:
        if bond in ccvBonds:
          ccvAtomBindings[ccvAtom].append(bond.bondType)
    
    #
    # Sort the atom names for consistent output
    #
    
    ccvAtomNames.sort()
    
    #
    # Check basic bindings...
    #
    
    errorList = []
    
    for ccvAtomName in ccvAtomNames:
      
      ccvAtom = ccvAtomDict[ccvAtomName]
      
      errorText = None
      
      if not ccvAtomBindings[ccvAtom]:
        if len(ccvAtoms) != 1:
          errorText =  "    ChemAtom %s has no bonds." % ccvAtom.name
    
      elif ccvAtom.elementSymbol == 'H':
        if ccvAtomBindings[ccvAtom] != ['single']:
          errorText =  "    %s: Proton has bonds %s." % (ccvAtom.name,ccvAtomBindings[ccvAtom])
      
      else:
      
        numBond = 0
        for bondType in ccvAtomBindings[ccvAtom]:
          if bondType == 'single':
            numBond += 1
          elif bondType == 'double':
            numBond += 2
          elif bondType == 'triple':
            numBond += 3
          elif bondType in ['aromatic','singleplanar']:
            if ccvAtom.elementSymbol == 'C' and bondType == 'aromatic' and ccvAtomBindings[ccvAtom].count('aromatic') > 2 and numBond == 3:
              numBond += 1
            elif ccvAtom.elementSymbol == 'O' and bondType == 'aromatic' and ccvAtomBindings[ccvAtom].count('aromatic') > 1 and numBond == 1.5:
              numBond += 1
            else:
              numBond += 1.5
          elif bondType == 'dative':
            numBond = 0
          else:
            print " ERROR: NOT USING BONDTYPE %s" % bondType
        
        if ccvAtom.elementSymbol == 'C' and numBond != 4:
          errorText =  "    %-4s: Carbon has %.1f bonds %s." % (ccvAtom.name,numBond,ccvAtomBindings[ccvAtom])
      
        elif ccvAtom.elementSymbol == 'N' and not (3 <= numBond <= 4):
          errorText =  "    %-4s: Nitrogen has %.1f bonds %s." % (ccvAtom.name,numBond,ccvAtomBindings[ccvAtom])

        elif ccvAtom.elementSymbol == 'O' and not (1.5 <= numBond <= 2.5):
          errorText =  "    %-4s: Oxygen has %.1f bonds %s." % (ccvAtom.name,numBond,ccvAtomBindings[ccvAtom])
      
      if errorText:
        errorList.append(errorText)
        
    if errorList:

      print "  %s, %s" % (ccv.linking,ccv.descriptor)
      
      for errorText in errorList:
        print errorText
      
      print
       
###################
# Main of program #
###################

if __name__ == "__main__":  

  import sys
  
  (testMode,writeData,verbose,ccpCodeList,chemCompDataDir) = initialiseChemCompScript(sys.argv)
    
  for (molType,ccpCodes) in ccpCodeList:
  
    if not ccpCodes:
      print "NO %s" % molType
      continue
    
    for ccpCode in ccpCodes[:1]:

      project = Implementation.MemopsRoot(name = 'tempData')

      print drawBox("ChemComp %s, %s" % (molType, ccpCode))

      chemComp = getChemComp(project,molType,ccpCode,download=False,chemCompArchiveDir='lala',copyFile=False)

      if chemComp:
        checkChemComp(chemComp,verbose = verbose)
      else:
        print "  ERROR: not available!"
        print
