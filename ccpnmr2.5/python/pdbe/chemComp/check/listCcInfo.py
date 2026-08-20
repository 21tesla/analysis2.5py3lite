"""

List information for a chemComp.

"""
from memops.api import Implementation
from pdbe.chemComp.Constants import editChemCompDataDir, testChemCompDataDir
from ccp.general.Io import getChemComp

if __name__ == '__main__':

  testMode = True

  molType = 'carbohydrate'
  ccpCode = 'dgal-hex-1-5:C1_OMe:C2_NAc'
  
  if testMode == True:
    chemCompArchiveDir = testChemCompDataDir
    
  else:
    chemCompArchiveDir = editChemCompDataDir
  #refDir = '/ebi/msd/nmrqual/chemComps/test'
  #refDir = '/ebi/msd/nmrqual/msdWeb/chemCompXml.v1.0.1'

  project = Implementation.MemopsRoot(name = 'test')
  chemComp = getChemComp(project,molType,ccpCode,chemCompArchiveDir = chemCompArchiveDir, copyFile = False)
  
  chemAtomOrSets = []
  i = -1
  
  for tempChemAtomOrSet in chemComp.sortedChemAtoms() + chemComp.sortedChemAtomSets():
    for i in range(0,len(chemAtomOrSets)):
      if tempChemAtomOrSet.name < chemAtomOrSets[i].name:
        chemAtomOrSets.insert(i,tempChemAtomOrSet)
        break

    if i == -1 or i == len(chemAtomOrSets)-1:
      chemAtomOrSets.append(tempChemAtomOrSet)
  
  
  namingSystemList = chemComp.sortedNamingSystems()
  
  for chemAtomOrSet in chemAtomOrSets:

    print chemAtomOrSet.name
    
    chemAtomSysNames = []
    i = -1
    
    for namingSystem in namingSystemList:
      for tempCasn in namingSystem.findAllAtomSysNames(atomName = chemAtomOrSet.name, atomSubType = chemAtomOrSet.subType):
        print "   %15s %s" % (namingSystem.name, chemAtomSysName.sysName)
