try:
  from pdbe.general.localConstants import pdbeDataDir
except:
  from ccp.general.Io import getDataPath
  pdbeDataDir = getDataPath('pdbe')
