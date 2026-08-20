import os
from memops.universal.Io import getTopDirectory

# TODO: make these local-settable as well?

try:
  from localConstants import dataDir
except:
  dataDir = os.path.join(getTopDirectory(),'data')

try:
  from localConstants import tmpDataDir
except:
  tmpDataDir = os.path.join(dataDir,'tmp')

try:
  from localConstants import archivesDataDir
except:
  archivesDataDir = os.path.join(dataDir,'archives')

generalReferenceDir = os.path.join(archivesDataDir,'reference')
alignMatrixDir = os.path.join(dataDir,'alignMatrices')

try:
  from localConstants import customArchivesDataDir
except:
  customArchivesDataDir = os.path.join(archivesDataDir,'custom')
  
# TODO THIS ONE MIGHT GO!
try:
  from localConstants import archivesCcpnDataDir
except:
  archivesCcpnDataDir = os.path.join(dataDir,'archivesCcpn')

try:
  from localConstants import pythonCommand
except:
  pythonCommand = "python2.5"

try:
  from localConstants import numCpu
except:
  numCpu = 1

try:
  from localConstants import chemCompArchiveDataDir
except:
  from ccp.general.Io import getChemCompArchiveDataDir
  chemCompArchiveDataDir = getChemCompArchiveDataDir()

#
# Archive sites
#

bmrbUrl = "http://www.bmrb.wisc.edu"
bmrbRestUrl = "http://rest.bmrb.wisc.edu"

try:
  from localConstants import pdbFtp
except:
  pdbFtp = "ftp.wwpdb.org" 

try:
  from localConstants import pdbFtpDir
except:
  pdbFtpDir = "pub/pdb/data/structures/divided/pdb/"
  
try:
  from localConstants import ebiPdbFtp
except:
  ebiPdbFtp = "ftp.ebi.ac.uk"

try:
  from localConstants import atlasUrl
except:
  atlasUrl = "http://www.ebi.ac.uk/pdbe/pdbelite/atlas/"
