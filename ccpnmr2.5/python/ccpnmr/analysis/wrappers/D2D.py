import os
import tempfile

from memops.universal.Url import uploadFile

from ccpnmr.format.general.Conversion import FormatConversion

RES_NUM_KEY = 'num'
SS_KEY = 'SS'
SEC_STRUC_KEYS = ('Helix', 'Beta', 'Coil', 'PPII')
SEC_STRUC_TIPS = tuple(['Probability residue is %s' % xx for xx in SEC_STRUC_KEYS])

UNRELIABLE = '*'

SERVER = 'http://www-vendruscolo.ch.cam.ac.uk/d2D'
SCRIPT = 'd2D.cc536.php'
URL = os.path.join(SERVER, SCRIPT)

def runD2D(chain, shiftList):

  fileName = tempfile.mkstemp()[1]
  predictionDict = {}

  try:

    # FIRST: run Format Converter to export shifts to ShiftY file
    chains = [chain]
    fc = FormatConversion(ccpnProject=chain.root)
    fc.exportFile('shifts', 'shifty', fileName, addKeywords={'measurementList': shiftList, 'chains': chains})

    # SECOND: run D2D on remote server
    fields = {'MAX_FILE_SIZE': 7000000, 'ph': 'No' }
    fileKey = 'userfile'
    response = uploadFile(URL, fileKey, fileName, fields)
    response = response.split('\n')

    # THIRD: parse results and return
    ss = 'File failed to upload correctly'
    tt = '#num'
    vv = '#DONE!'
    fields = None
    for line in response:
      if ss in line:
        raise Exception('D2D: %s' % ss)
      if line.startswith(tt):
        fields = line[1:].split()
      elif line.startswith(vv):
        break
      elif fields:
        values = line.split()
        if len(values) == len(fields):
          residue = secStrucCode = None
          probabilityDict = {}
          for n, value in enumerate(values):
            field = fields[n]
            if field in SEC_STRUC_KEYS:
              probabilityDict[field] = float(value)
            elif field == RES_NUM_KEY:
              seqCode = int(value)
              residue = chain.findFirstResidue(seqCode=seqCode)
            elif field == SS_KEY:
              isReliable = not value.endswith(UNRELIABLE)
              if not isReliable:
                value = value[:-1]
              secStrucCode = value
          predictionDict[residue] = (secStrucCode, isReliable, probabilityDict)

  finally:
    try: 
      # sometimes, at least on Windows, get following error:
      # The process cannot access the file because it is being used by another process
      os.remove(fileName)
    except:
      pass

  return predictionDict

