
from ccp.general.Io import getChemComp
from ccp.general.Util import MakeAtomSets
from memops.api import Implementation
from pdbe.chemComp.Constants import editChemCompDataDir
from pdbe.chemComp.Util import getCcpCodeList

###################
# Main of program #
###################

if __name__ == "__main__":

  makeAtomSets = MakeAtomSets()

  ccpCodeList = getCcpCodeList(editChemCompDataDir)

  for (molType,ccpCodes) in ccpCodeList:

    for ccpCode in ccpCodes[:5]:

      project = Implementation.MemopsRoot(name = 'tempData')

      chemComp = getChemComp(project,molType,ccpCode,copyFile=False,download=False)
      chemComp.isModifiable = True

      print("\nChemComp code %s, molType %s" % (chemComp.ccpCode,chemComp.molType))

      makeAtomSets.setChemComps(chemComp)

      #if chemComp.isModified:
      #  chemComp.save()
