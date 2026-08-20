importVersionSep = '_'

unknownMapping = 'unknown_mapping'

try:
  from localConstants import mappingImportLocation
except:
  mappingImportLocation = 'pdbe.nmrStar.IO'