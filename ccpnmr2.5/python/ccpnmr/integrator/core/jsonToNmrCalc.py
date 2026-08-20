import sys, json, traceback

from memops.general import Io as genIo
from ccpnmr.integrator.core import Util as intUtil


if __name__ == '__main__':
  """ Make NmrCalc.Run from json file.
  """
  
  # get input arguments
  if len(sys.argv) == 3:
    projectDir, jsonFile  = sys.argv[1:3]
    nmrProjectName = None
  
  else:
    raise Exception("Usage: projectToJson projectDir jsonFilePath")
  
  try:
    memopsRoot = genIo.loadProject(projectDir, suppressGeneralDataDir=True)
    jsonObject = json.load(open(jsonFile))
    tmpFile = jsonFile + '_tmp'
    json.dump(jsonObject, open(tmpFile, 'w'),  sort_keys=True, indent=2)
    nmrCalcId = intUtil.makeNmrCalc(memopsRoot, jsonObject)
    memopsRoot.saveModified()
    jsonObject['CCPN.nmrCalcId'] = nmrCalcId
    fp = open(jsonFile, 'w')
    json.dump(jsonObject, fp,  sort_keys=True, indent=2)
    fp.close()
  except:
    traceback.print_exc()
    raise
