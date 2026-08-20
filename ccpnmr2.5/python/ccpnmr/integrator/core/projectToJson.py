import os, sys, json

from memops.general import Io as genIo
from ccpnmr.integrator.core import Util as intUtil


if __name__ == '__main__':
  """ Create project summary JSON file
  """
  
  # get input arguments
  if len(sys.argv) == 3:
    projectDir, jsonFile  = sys.argv[1:3]
    nmrProjectName = None
    
  elif len(sys.argv) == 4:
    projectDir, jsonFile, nmrProjectName = sys.argv[1:4]
  
  else:
    raise Exception("Usage: projectToJson projectDir jsonFilePath [nmrProjectName optional]")
  
  memopsRoot = genIo.loadProject(projectDir, suppressGeneralDataDir=True)
  jsonObject = intUtil.wmsProjectSummary(memopsRoot, 
                                         nmrProjectName=nmrProjectName)
  
  # If jsonFile is relative this puts it relative to project repository
  # If it is absoolute it goes where it should
  jsonFile = os.path.join(
              memopsRoot.findFirstActiveRepository().url.dataLocation,
              jsonFile
             )
  fp = open(jsonFile, 'w')
  json.dump(jsonObject, fp,  sort_keys=True, indent=2)
  fp.close()
