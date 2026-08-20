from ccp.general.Io import getDataPath

# This maybe not so good? Or only relevant for large scale testing, not for people to work on own
# workflow
try:
  from localConstants import importDataDir
except:
  importDataDir=getDataPath('ccpnmr','workflow')

try:
  from localConstants import cnsExec
except:
  cnsExec = 'cns_solve'


programList = ['Aria','Cing','Cns','Fc','Haddock','Isd']

