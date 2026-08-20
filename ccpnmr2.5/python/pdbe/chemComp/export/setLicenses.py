from memops.scripts.license.headers import setLicenses

from pdbe.chemComp.Constants import editChemCompDataDir

import sys

if __name__ == '__main__':

  if sys.argv[1:]:
    mode = sys.argv[1]
  else:
    mode = 'test'

  setLicenses(directory = editChemCompDataDir,mode = mode)

