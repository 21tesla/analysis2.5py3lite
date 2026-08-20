import os
from memops.universal.Io import getTopDirectory, joinPath

from ccp.general.Io import getDataPath

from pdbe.general.Io import pdbeDataDir

#
# Tag for newly created chemComp(Coord)s - will be removed if modified afterwards!!!
#

freshTag = 'objectGeneratedFromScratch'

#
# Directory containing reference information for substituents
#

origMol2DataDir = getDataPath('ccp','mol2')

#
# Work chemComp directory, in all/
#

try:
  from localConstants import workChemCompDir
except:
  # Have to work on the actual repository directory for editing!
  workChemCompDir = joinPath(pdbeDataDir,'chemComp')

refDataDir =      os.path.join(workChemCompDir,'refData')

tempChemCompDir = os.path.join(workChemCompDir,'tmp')

testChemCompDir = os.path.join(workChemCompDir,'test')
testChemCompDataDir = os.path.join(testChemCompDir,'ChemComp')
testChemCompCoordDataDir = os.path.join(testChemCompDir,'ChemCompCoord')

obsoleteChemCompDir = os.path.join(workChemCompDir,'obsolete')
obsoleteChemCompDataDir = os.path.join(obsoleteChemCompDir,'ChemComp')
obsoleteChemCompCoordDataDir = os.path.join(obsoleteChemCompDir,'ChemCompCoord')

#
# ChemComp directory for editing, in ccpn-chemcomp repository directly
#

try:
  from localConstants import editChemCompDir
except:
  editChemCompDir = joinPath(getTopDirectory(),'..','ccpn-chemcomp','data','pdbe','chemComp','archive')

editChemCompDataDir = joinPath(editChemCompDir,'ChemComp')
editChemCompCoordDataDir = joinPath(editChemCompDir,'ChemCompCoord')

editLinkedChemCompDir = os.path.join(workChemCompDir,'archive')
editLinkedChemCompDataDir = os.path.join(editLinkedChemCompDir,'ChemComp')
editLinkedChemCompCoordDataDir = os.path.join(editLinkedChemCompDir,'ChemCompCoord')


#
# Repository for standard residues, in ccpn sourceforge repository directly
#

try:
  from localConstants import standardRepChemCompDataDir
except:
  standardRepChemCompDataDir = joinPath(getTopDirectory(),'..','ccpn','data')
