"""
======================COPYRIGHT/LICENSE START==========================

Psipred.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
# Wrapper for PSIPRED V2.6 by David Jones
# Predicts protein backbone secondary structure
# from its primary sequence

import os, re

osPath = os.path
reSub = re.sub

PSIPRED_DIR = os.environ.get('PSIPRED_DIR', '.')

def testCcpnPsipred(argServer):

  p = argServer.getProject()

  chain = p.findFirstMolSystem().findFirstChain()

  ssDict = psipredCcpn(chain)
  
  print ssDict
  
def psipredCcpn(chain, smooth=1, alphaBias=1.0, betaBias=1.0):
  """Descrn: Get the secondary structure predictions for the residues
             of a chain. Option to adjust degree of smoothing.
             Options to bias alpha helix or beta sheet content.
     Inputs: MolSystem.Chain, Int, Float, Float
     Output: Dict of MolSystem.Residue:Char (DSSP code C/H/E),
             Line, Line (output files)
  """
  
  memopsRoot = chain.root

  dataRepository = memopsRoot.findFirstRepository(name='userData')

  projPath = dataRepository.url.dataLocation

  prefix = '%s_%s' % (memopsRoot.name, chain.molecule.name)
  prefix = reSub('\s+',r'_',prefix)
 
  psipredDir = osPath.join(projPath, 'psipred')
  if not osPath.exists(psipredDir):
    os.mkdir(psipredDir)
    
  seqFileName = osPath.join(psipredDir, '%s_seq.fasta' % prefix)
  _writeChainSeq(chain, seqFileName)
 
  execDir = osPath.join(PSIPRED_DIR, 'bin')
  dataDir = osPath.join(PSIPRED_DIR, 'data')

  ssFile, prettyFile = psipred(seqFileName, psipredDir, prefix,
                               smooth, alphaBias, betaBias,
                               execDir, dataDir)

  ssCodes = _readSecStrucCodes(ssFile)

  ssCodeDict = {}
  for i, residue in enumerate(chain.sortedResidues()):
    ssCodeDict[residue] = ssCodes[i]
    
  return ssCodeDict, ssFile, prettyFile

def _readSecStrucCodes(fileName):
  """ Read PSIPRED output secondary structure codes """

  ssCodes = []
  file = open(fileName)
  
  line = file.readline()
  while line:
    data = line.split()
    
    if data and (data[0][0] != '#'):
      ssCodes.append(data[2])
  
    line = file.readline()
  
  file.close()
  
  return ssCodes

def _writeChainSeq(chain, fileName):
  """ Write a chain's 1-letter sequence to a FASTA file """

  seq = chain.molecule.stdSeqString.upper()
  seq = reSub('\*','X', seq)
  _writeSeq(chain.molecule.name, seq, fileName)


def _writeSeq(name, seq, fileName):
  """ Wite a single sequence FASTA style file """
  
  sequenceFile = open(fileName, 'w')
  seq = reSub('\s+','', seq)
  seq = reSub('(\S{60})(\S)',r'\1\n\2',seq)
  sequenceFile.write('> %s\n%s\n' % (name, seq))
  sequenceFile.close()


def psipred(seqFile, psipredDir, prefix,
            smooth=1, alphaBias=1.0, betaBias=1.0,
            execDir='./bin', dataDir='./data'):
  """ Python wrapper to run PSIPRED  V2.6 by David Jones """

  matrixFile = osPath.join(psipredDir, '%s.mtx' % prefix)
  horizFile = osPath.join(psipredDir, '%s.horiz' % prefix)
  ssFile = osPath.join(psipredDir, '%s.ss' % prefix)
  ssP2file = osPath.join(psipredDir, '%s.ss2' % prefix)

  seq2mtx = osPath.join(execDir, 'seq2mtx')

  os.system('%s %s > %s' % (seq2mtx, seqFile, matrixFile))

  print "PSIPRED: Predicting secondary structure based on single sequence ..."

  print "PSIPRED: Pass1 ..."
  
  psipred = osPath.join(execDir, 'psipred')
 
  dat1 = osPath.join(dataDir, 'weights_s.dat')
  if osPath.exists(dat1):
    dat2 = osPath.join(dataDir, 'weights_s.dat2')
    dat3 = osPath.join(dataDir, 'weights_s.dat3')
  else:
    # looks like file name changed
    dat1 = osPath.join(dataDir, 'weights.dat')
    dat2 = osPath.join(dataDir, 'weights.dat2')
    dat3 = osPath.join(dataDir, 'weights.dat3')
 
  data = (psipred, matrixFile, dat1, dat2, dat3, ssFile)
  os.system('%s %s %s %s %s > %s' % data)

  print "PSIPRED: Pass2 ..."
  
  datP2 = osPath.join(dataDir, 'weights_p2.dat')
  psipass2 = osPath.join(execDir, 'psipass2')
  
  data = (psipass2, datP2, smooth, alphaBias,
          betaBias, ssP2file, ssFile, horizFile)
  os.system('%s %s %d %f %f %s %s > %s' % data)

  os.unlink(matrixFile)
  #os.unlink('error.log')

  # more ...

  print 'PSIPRED: Done'

  return ssP2file, horizFile
