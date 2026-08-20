import os

from memops.universal.Url import fetchUrl

from ccp.general.Constants import ccpCodeToCode1LetterDict
proteinDict = ccpCodeToCode1LetterDict['protein']

SERVER = 'http://www-vendruscolo.ch.cam.ac.uk/'
SCRIPTS = ['camcoil%s.ggt23.php' % ss for ss in ('', '_lfp')]

COIL_SCRIPT = 'Random coil shifts'
LFP_SCRIPT = 'Protein loops shifts'
SCRIPT_TEXTS = [COIL_SCRIPT, LFP_SCRIPT]
SCRIPT_PHS = ['2.0', '6.1']

# mapping from CamCoil atom names to CCPN atom names
# the C0 is a kludge because of a typo in the CamCoil LFP output
atomNameMapping = {'Ca': 'CA', 'Cb': 'CB', 'CO': 'C', 'HN': 'H', 'Ha': 'HA', 'N': 'N', 'C0': 'C'}

def runCamCoil(chain, pH=None, script=COIL_SCRIPT):

  if script not in SCRIPT_TEXTS:
    script = SCRIPT_TEXTS[0]

  if pH is None:
    pH = ''
  elif script==LFP_SCRIPT:
    pH = SCRIPT_PHS[0]
  elif pH not in SCRIPT_PHS:
    pH = SCRIPT_PHS[0]

  scriptNumber = SCRIPT_TEXTS.index(script)

  atomNames = []
  atomShiftDict = {}

  # FIRST: extract 1-letter codes
  residues = []
  sequence = []
  for residue in chain.sortedResidues():
    if residue.molType == 'protein':
      ccpCode = residue.ccpCode
      if ccpCode in proteinDict:
        residues.append(residue)
        sequence.append(proteinDict[ccpCode])

  if not sequence:
    raise Exception('No valid protein residue found')

  #print ''.join(sequence)

  # SECOND: run CamCoil on remote server
  ss = 'No Valid Sequence Uploaded!'
  values = {'sequence': ''.join(sequence), 'pH': pH}
  url = os.path.join(SERVER, SCRIPTS[scriptNumber])
  response = fetchUrl(url, values)
  if ss in response:
    raise Exception('No valid sequence uploaded')

  response = response.split('\n')

  # THIRD: parse results and return
  tt = 'Pos. Res.'
  vv = '</font>'
  resNum = 0
  for line in response:
    if line.startswith(tt):
      atomNames = [atomNameMapping.get(atomName, atomName) for atomName in line.split()[2:]]
    elif line.startswith(vv):
      break
    elif atomNames:
      values = line.split()
      if int(values[0]) != resNum+1:
        raise Exception('Inconsistent residue number in output, expecting %s, got %s' % (resNum+1, int(values[0])))
 
      values = values[2:]
      if len(values) == len(atomNames):
        residue = residues[resNum]
        for n, value in enumerate(values):
          try:
            value = float(value)
          except:  # value = '-----' or similar
            value = None
          atomName = atomNames[n]
          if residue.linking == 'start' and atomName == 'H':
            atomName = 'H1'
          atom = residue.findFirstAtom(name=atomName)
          if atom:
            atomShiftDict[atom] = value
      resNum += 1
      if resNum >= len(residues):
        break

  return atomShiftDict

if __name__ == '__main__':

  import sys

  if len(sys.argv) != 2:
    print 'need to specify project directory'
    sys.exit()

  projectDir = sys.argv[1]

  from memops.general.Io import loadProject

  project = loadProject(projectDir)
  chain = project.findFirstMolSystem().findFirstChain()

  atomShiftDict = runCamCoil(chain)

  for residue in chain.sortedResidues():
    for atom in residue.sortedAtoms():
      value = atomShiftDict.get(atom)
      if value is not None:
        print '%s\t%s\t%s\t%s' % (residue.seqCode, residue.ccpCode, atom.name, value)

