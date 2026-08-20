import sys

from pdbe.chemComp.Constants import editChemCompDataDir, testChemCompDataDir

from ccp.general.Io import getChemComp
from ccp.general.Util import createNewDescriptors, getDescriptorDict

from memops.api import Implementation

if __name__ == '__main__':

  #
  # This script adds linkings to chemComps
  #
  # Use as:
  #
  # python addLinkingToChemComp.py edit|test molType ccpCode link:heavyAtomName@remove:protonAtomName
  #

  argList = ['edit','test']
  if not len(sys.argv) == 5 or sys.argv[1] not in argList:
    print "Need one of %s as first argument to run script,\nthen give molType, ccpCode and new linking/descriptor to add with atom to remove (e.g. addLinkingToChemComp.py edit protein ASN link:ND2@remove:HD2)" % str(argList)
    sys.exit()
    
  if sys.argv[1] == 'test':
    chemCompDir = testChemCompDataDir
  elif sys.argv[1] == 'edit':
    chemCompDir = editChemCompDataDir
  
  molType = sys.argv[2]
  ccpCode = sys.argv[3]
  newLinking = sys.argv[4]
  removeAtomNames = []
  renameAtoms = {}
  
  #
  # Also allows renaming of atoms: E.g. HD21/HD22: remove:HD22.rename:HD21,HD2
  #
    
  if '@' in newLinking:
    (newLinking,changeAtomText) = newLinking.split('@')
    changeAtomTextEls = changeAtomText.split('.')
    for changeAtomTextEl in changeAtomTextEls:
      for allowedEl in ('remove:','rename:'):
        if changeAtomTextEl.count(allowedEl):
          infoList = changeAtomTextEl.replace(allowedEl,'').split(',')
          if allowedEl == 'remove:':
            removeAtomNames.extend(infoList)
          elif allowedEl == 'rename:':
            renameAtoms[infoList[0]] = infoList[1]
  
  newDict = getDescriptorDict(newLinking)
  
  if not newDict:
    print "  Exiting: cannot handle '%s' as descriptor/linking code." % newLinking
    sys.exit()
  
  #
  # Load the chemComp...
  #

  tempProject = Implementation.MemopsRoot(name = 'tempData')
  
  chemComp = getChemComp(tempProject,molType,ccpCode,chemCompArchiveDir=chemCompDir)
    
  if not chemComp:
    print "  Exiting: unknown chemComp '%s','%s'" % (molType,ccpCode)
    sys.exit()

  createNewDescriptors(chemComp,newDict,removeAtomNames,renameAtoms)
