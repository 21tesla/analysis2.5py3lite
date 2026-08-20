LICENSE = """
======================COPYRIGHT/LICENSE START==========================

MoleculeAlign.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2013 Wayne Boucher, REasmus Fogh, Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

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
import re

from ccp.lib.MoleculeQuery import DNA_MOLTYPE, RNA_MOLTYPE
from ccp.lib.MoleculeQuery import PROTEIN_MOLTYPE, DNARNA_MOLTYPE

BLOSUM62={'A':{'A': 4,'R':-1,'N':-2,'D':-2,'C': 0,'Q':-1,'E':-1,'G': 0,'H':-2,'I':-1,'L':-1,
               'K':-1,'M':-1,'F':-2,'P':-1,'S': 1,'T': 0,'W':-3,'Y':-2,'V': 0,'X':0},
          'R':{'A':-1,'R': 5,'N': 0,'D':-2,'C':-3,'Q': 1,'E': 0,'G':-2,'H': 0,'I':-3,'L':-2,
               'K': 2,'M':-1,'F':-3,'P':-2,'S':-1,'T':-1,'W':-3,'Y':-2,'V':-3,'X':0},
          'N':{'A':-2,'R': 0,'N': 6,'D': 1,'C':-3,'Q': 0,'E': 0,'G': 0,'H': 1,'I':-3,'L':-3,
               'K': 0,'M':-2,'F':-3,'P':-2,'S': 1,'T': 0,'W':-4,'Y':-2,'V':-3,'X':0},
          'D':{'A':-2,'R':-2,'N': 1,'D': 6,'C':-3,'Q': 0,'E': 2,'G':-1,'H':-1,'I':-3,'L':-4,
               'K':-1,'M':-3,'F':-3,'P':-1,'S': 0,'T':-1,'W':-4,'Y':-3,'V':-3,'X':0},
          'C':{'A': 0,'R':-3,'N':-3,'D':-3,'C': 9,'Q':-3,'E':-4,'G':-3,'H':-3,'I':-1,'L':-1,
               'K':-3,'M':-1,'F':-2,'P':-3,'S':-1,'T':-1,'W':-2,'Y':-2,'V':-1,'X':0},
          'Q':{'A':-1,'R': 1,'N': 0,'D': 0,'C':-3,'Q': 5,'E': 2,'G':-2,'H': 0,'I':-3,'L':-2,
               'K': 1,'M': 0,'F':-3,'P':-1,'S': 0,'T':-1,'W':-2,'Y':-1,'V':-2,'X':0},
          'E':{'A':-1,'R': 0,'N': 0,'D': 2,'C':-4,'Q': 2,'E': 5,'G':-2,'H': 0,'I':-3,'L':-3,
               'K': 1,'M':-2,'F':-3,'P':-1,'S': 0,'T':-1,'W':-3,'Y':-2,'V':-2,'X':0},
          'G':{'A': 0,'R':-2,'N': 0,'D':-1,'C':-3,'Q':-2,'E':-2,'G': 6,'H':-2,'I':-4,'L':-4,
               'K':-2,'M':-3,'F':-3,'P':-2,'S': 0,'T':-2,'W':-2,'Y':-3,'V':-3,'X':0},
          'H':{'A':-2,'R': 0,'N': 1,'D':-1,'C':-3,'Q': 0,'E': 0,'G':-2,'H': 8,'I':-3,'L':-3,
               'K':-1,'M':-2,'F':-1,'P':-2,'S':-1,'T':-2,'W':-2,'Y': 2,'V':-3,'X':0},
          'I':{'A':-1,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-3,'E':-3,'G':-4,'H':-3,'I': 4,'L': 2,
               'K':-3,'M': 1,'F': 0,'P':-3,'S':-2,'T':-1,'W':-3,'Y':-1,'V': 3,'X':0},
          'L':{'A':-1,'R':-2,'N':-3,'D':-4,'C':-1,'Q':-2,'E':-3,'G':-4,'H':-3,'I': 2,'L': 4,
               'K':-2,'M': 2,'F': 0,'P':-3,'S':-2,'T':-1,'W':-2,'Y':-1,'V': 1,'X':0},
          'K':{'A':-1,'R': 2,'N': 0,'D':-1,'C':-3,'Q': 1,'E': 1,'G':-2,'H':-1,'I':-3,'L':-2,
               'K': 5,'M':-1,'F':-3,'P':-1,'S': 0,'T':-1,'W':-3,'Y':-2,'V':-2,'X':0},
          'M':{'A':-1,'R':-1,'N':-2,'D':-3,'C':-1,'Q': 0,'E':-2,'G':-3,'H':-2,'I': 1,'L': 2,
               'K':-1,'M': 5,'F': 0,'P':-2,'S':-1,'T':-1,'W':-1,'Y':-1,'V': 1,'X':0},
          'F':{'A':-2,'R':-3,'N':-3,'D':-3,'C':-2,'Q':-3,'E':-3,'G':-3,'H':-1,'I': 0,'L': 0,
               'K':-3,'M': 0,'F': 6,'P':-4,'S':-2,'T':-2,'W': 1,'Y': 3,'V':-1,'X':0},
          'P':{'A':-1,'R':-2,'N':-2,'D':-1,'C':-3,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-3,'L':-3,
               'K':-1,'M':-2,'F':-4,'P': 7,'S':-1,'T':-1,'W':-4,'Y':-3,'V':-2,'X':0},
          'S':{'A': 1,'R':-1,'N': 1,'D': 0,'C':-1,'Q': 0,'E': 0,'G': 0,'H':-1,'I':-2,'L':-2,
               'K': 0,'M':-1,'F':-2,'P':-1,'S': 4,'T': 1,'W':-3,'Y':-2,'V':-2,'X':0},
          'T':{'A': 0,'R':-1,'N': 0,'D':-1,'C':-1,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-1,'L':-1,
               'K':-1,'M':-1,'F':-2,'P':-1,'S': 1,'T': 5,'W':-2,'Y':-2,'V': 0,'X':0},
          'W':{'A':-3,'R':-3,'N':-4,'D':-4,'C':-2,'Q':-2,'E':-3,'G':-2,'H':-2,'I':-3,'L':-2,
               'K':-3,'M':-1,'F': 1,'P':-4,'S':-3,'T':-2,'W':11,'Y': 2,'V':-3,'X':0},
          'Y':{'A':-2,'R':-2,'N':-2,'D':-3,'C':-2,'Q':-1,'E':-2,'G':-3,'H': 2,'I':-1,'L':-1,
               'K':-2,'M':-1,'F': 3,'P':-3,'S':-2,'T':-2,'W': 2,'Y': 7,'V':-1,'X':0},
          'V':{'A': 0,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-2,'E':-2,'G':-3,'H':-3,'I': 3,'L': 1,
               'K':-2,'M': 1,'F':-1,'P':-2,'S':-2,'T': 0,'W':-3,'Y':-1,'V': 4,'X':0},
          'X':{'A': 0,'R': 0,'N': 0,'D': 0,'C': 0,'Q': 0,'E': 0,'G': 0,'H': 0,'I': 0,'L': 0,
               'K': 0,'M': 0,'F': 0,'P': 0,'S': 0,'T': 0,'W': 0,'Y': 0,'V': 0,'X':0}
          }
 
 
NAMAT={'A':{'A': 2,'C':-5,'G':-5,'T':-5,'U':-5,'R': 1,'Y':-5,'N':0,'X':0},
       'C':{'A':-5,'C': 2,'G':-5,'T':-5,'U':-5,'R':-5,'Y': 1,'N':0,'X':0},
       'G':{'A':-5,'C':-5,'G': 2,'T':-5,'U':-5,'R': 1,'Y':-5,'N':0,'X':0},
       'T':{'A':-5,'C':-5,'G':-5,'T': 2,'U': 2,'R':-5,'Y': 1,'N':0,'X':0},
       'U':{'A':-5,'C':-5,'G':-5,'T': 2,'U': 2,'R':-5,'Y': 1,'N':0,'X':0},
       'R':{'A': 1,'C':-5,'G': 1,'T':-5,'U':-5,'R': 1,'Y':-5,'N':0,'X':0},
       'Y':{'A':-5,'C': 1,'G':-5,'T': 1,'U': 1,'R':-5,'Y': 1,'N':0,'X':0},
       'N':{'A': 0,'C': 0,'G': 0,'T': 0,'U': 0,'R': 0,'Y': 0,'N':0,'X':0},
       'X':{'A': 0,'C': 0,'G': 0,'T': 0,'U': 0,'R': 0,'Y': 0,'N':0,'X':0}}


def findMatchingChains(molSystem, ccpCodes, excludeChains=None, molTypes=None, doWarning=True):
  """Find the mol system chains that best match the input ccpCodes
             (like three letter codes). Useful for trying to match structures
             to existing molecular data. Optional argument to specify which
             chains cannot be matched.
  .. describe:: Input
  
  MolSystem.MolSystem, List of Strings (MolSystem.MolResidue.ccpCodes),
             List of MolSystem.Chains

  .. describe:: Output
  
  Tuple of (MolSystem.Chain, Int (index of first matching ccpCode),
             Int (first matching MolSystem.Residue.seqId))
  """

  chains0 = []
  for chain in molSystem.sortedChains():
    if excludeChains and (chain in excludeChains):
      continue
    
    #sequence = []
    if not chain.residues:
      continue
    
    chains0.append(chain)

  
  bestMapping = None
  bestChain   = None

  scoreList = []
  for chain in chains0:
    mapping, score = getSequenceResidueMapping(chain, ccpCodes)
   
    if score:
      scoreList.append((score, chain, mapping))
  
  chains = []
  mappings = [] 
  if scoreList:
    scoreList.sort()
    scoreList.reverse()
    bestScore, bestChain, bestMapping = scoreList[0]
    lenSeq = len(ccpCodes)
    
    # Find perfect stretches
    stretch = set([])
    i = 0
    for pair in bestMapping:
      if None in pair:
        stretch.add(i)
        i = 0
      else:
        i += 1
    
    stretch.add(i)
    
    if (bestScore/float(len(bestChain.residues)) > 2.0) or \
       ((max(stretch) == lenSeq) and (lenSeq > 20)):
      # No really bad alignments
    
      chains = [bestChain, ]
      mappings = [bestMapping, ]
      for score, chain, mapping in scoreList[1:]:
        if score == bestScore:
          chains.append(chain)
          mappings.append(mapping)   

  return chains, mappings


   
def getSequenceResidueMapping(chain, sequence):
  """Match a sequence of residue ccpCodes to an existing chain and calculate the alignment score.
  .. describe:: Input
  
  MolSystem.Chain, List of Words (MolSystem.Residue.ccpCodes)

  .. describe:: Output
  
  List of 2-List of [Int or None, MolSystem.Residue or None], Float (score)
  """
  from ccp.general.Constants import ccpCodeToCode1LetterDict
  
  molType = chain.molecule.molType

  if ccpCodeToCode1LetterDict.get(molType) is None:
    letterDict = {}
  else:
    letterDict = ccpCodeToCode1LetterDict[molType]
  
  seq1 = chain.molecule.stdSeqString
  seq2 = ''.join([letterDict.get(ccpCode, 'X') for ccpCode in sequence])
  
  if molType in (DNA_MOLTYPE, RNA_MOLTYPE, DNARNA_MOLTYPE):
    seqA, seqB, score = _sequenceAlign(seq1,seq2,NAMAT)
 
  else:
    seqA, seqB, score = _sequenceAlign(seq1,seq2,BLOSUM62)
 
  sortedResidues = chain.sortedResidues()

  x = 0
  y = 0
  mapping = []
  for i in range(len(seqA)):
    mapping.append([None,None])
    if seqA[i] != '-':
      mapping[i][1] = sortedResidues[x] 
      x += 1
      
    if seqB[i] != '-':
      mapping[i][0] = y 
      y += 1

  return mapping, score
  
  

def getChainResidueMapping(chainA, chainB):
  """Find the corresponding pairs of residues in two sequence similar chains.
             Matches independent chain fragments, which can be of different molTypes
  .. describe:: Input
  
  MolSystem.Chain, MolSystem.Chain

  .. describe:: Output
  
  List of List of [MolSystem.Residue or None]
  """
  
  fragmentsA = []
  for fragment in chainA.chainFragments:
    molType  = fragment.findFirstResidue().molResidue.molType
    seq      = ''
    residues = fragment.residues
    for residue in residues:
      seq += residue.chemCompVar.chemComp.code1Letter or 'X'

    fragmentsA.append( (residues, molType, seq, fragment) )
  

  fragmentsB = []
  for fragment in chainB.chainFragments:
    molType  = fragment.findFirstResidue().molResidue.molType
    seq      = ''
    residues = fragment.residues
    for residue in residues:
      seq += residue.chemCompVar.chemComp.code1Letter or 'X'
  
    fragmentsB.append( (residues, molType, seq, fragment) )
    
  
  scores = []  
  for residuesA, molTypeA, seqA, fragmentA in fragmentsA:
    for residuesB, molTypeB, seqB, fragmentB in fragmentsB:
      if molTypeB != molTypeA:
        continue
      
      if molType == PROTEIN_MOLTYPE:
        seq1, seq2, score = _sequenceAlign(seqA,seqB,BLOSUM62)

      elif molType in (DNA_MOLTYPE, RNA_MOLTYPE):
        seq1, seq2, score = _sequenceAlign(seqA,seqB,NAMAT)
      
      else:
        seq1, seq2 = seqA, seqB
        
        while len(seq1) < len(seq2):
          seq1 += '-'
        while len(seq2) < len(seq1):
          seq2 += '-'
        
        n1 = float(len(seqA))
        n2 = float(len(seqB))
        score = min(n1,n2)/max(n1,n2)
        
        if residuesA[0].seqCode == residuesB[0].seqCode:
          score += 1.0
  
      scores.append((score, seq1, seq2, residuesA, residuesB, fragmentA, fragmentB))
  
  scores.sort()
  
  done = {}
  totalScore = 0.0  
  mapping = []
  while scores:
    score, seq1, seq2, residuesA, residuesB, fragmentA, fragmentB = scores.pop()
    
    if done.get(fragmentA):
      continue

    if done.get(fragmentB):
      continue

    done[fragmentA] = True
    done[fragmentB] = True
    totalScore += score
    
    x = 0
    y = 0
    j = len(mapping)
    for i in range(len(seq1)):
      k = i+j
      mapping.append([None,None])
      residueA = None
      residueB = None
      
      if seq1[i] != '-':
        residueA = residuesA[x]
        if done.get(residueA):
          residueA = None
        else:
          mapping[k][0] = residueA
          
        x += 1
 
      if seq2[i] != '-':
        residueB = residuesB[y]
        if done.get(residueB):
          residueB = None
        else:
          mapping[k][1] = residueB
 
        y += 1
      
      if residueA and residueB:
        done[residueA] = True
        done[residueB] = True

  return mapping, totalScore
  
  
      
def _sequenceAlign(seq1,seq2,matrix,inspen=2,extpen=1):
  """Aligns two sequence strings (of one letter codes)
              - results are gapped with '-'; first tries an exact substring match,
             if that failts it then uses dynamicAlign
  .. describe:: Input
  
  String, String, Dict of Dicts (homology score matrix), Int, Int, Int, Int

  .. describe:: Output
  
  String, String (aligned & gapped sequence strings)
  """

  def _substringAlign(s1, s2, n1, n2, n):
    
    if n1 < n2:
      score = -inspen - extpen*(n2-n1-1)
    else:
      score  = 0
    for x in s1:
      score += matrix[x][x]

    sA = n*'-' + s1 + (n2-n1-n)*'-'
    sB = s2

    return (sA, sB, score)

  nX = len(seq1)
  nY = len(seq2)

  if nX <= nY:
    n = seq2.find(seq1)
    if n >= 0:
      return _substringAlign(seq1, seq2, nX, nY, n)
  else:
    n = seq1.find(seq2)
    if n >= 0:
      seqB, seqA, score = _substringAlign(seq2, seq1, nY, nX, n)
      return seqA, seqB, score

  return _dynamicAlign(seq1, seq2, matrix, inspen, extpen)

def _dynamicAlign(seq1,seq2,matrix,inspen=2,extpen=1):
  """Aligns two sequence strings (of one letter codes) - results are gapped with '-'
  .. describe:: Input
  
  String, String, Dict of Dicts (homology score matrix), Int, Int, Int, Int

  .. describe:: Output
  
  String, String (aligned & gapped sequence strings)
  """

  seq1 = re.sub('-|\s','',seq1)
  seq2 = re.sub('-|\s','',seq2)

  seq1 = re.sub('\*','X',seq1)
  seq2 = re.sub('\*','X',seq2)

  maxScore = 0
          
  nX    = len(seq1)+1
  nY    = len(seq2)+1

  M     = [[0] * nY]
  R     = [[2] * nY]
  route = 0

  R[0][0] = 0
  for x in range(1,nX):
    M.append([0,])
    R.append([1,])
    for y in range(1,nY):
      p1 = inspen
      p2 = inspen
      if route == 1:
        p1 = extpen
      elif route == 2:
        p2 = extpen

      score = matrix[seq1[x-1]][seq2[y-1]]
      paths = [M[x-1][y-1]+score,
               M[x-1][y]-p1,
               M[x][y-1]-p2]

      best  = max(paths)
      route = paths.index(best)

      M[x].append( best )
      R[x].append( route )

      if best >= maxScore:
        maxScore = best

  best = M[-1][-1] - 1
  for x0 in range(1, nX):
    m = M[x0][-1]
    if m > best:
      best = m
      x = x0
      y = nY-1
  for y0 in range(1, nY):
    m = M[-1][y0]
    if m > best:
      best = m
      x = nX-1
      y = y0

  route = R[x][y]
  dx = nX - x - 1
  dy = nY - y - 1
  if dx > 0:
    seqA = seq1[-dx:]
    seqB = dx * '-'
  elif dy > 0:
    dy = nY - y - 1
    seqA = dy * '-'
    seqB = seq2[-dy:]
  else:
    seqA = ''
    seqB = ''

  while  x > 0 or y > 0 :
    if route == 0:
      seqA = seq1[x-1] + seqA
      seqB = seq2[y-1] + seqB
      x   -= 1
      y   -= 1

    elif route == 1:
      seqB = '-' + seqB
      seqA = seq1[x-1] + seqA
      x   -= 1

    elif route == 2:
      seqA = '-' + seqA
      seqB = seq2[y-1] + seqB
      y   -= 1

    route = R[x][y]

  return seqA, seqB, maxScore


##########################################################################
#
#  Rasmus Fogh New alignment functions, May 2013:

def matchBySeqCode(molSystem, residToChemComp, excludeChains=None):
  """ return chain, {resId:msResidue.seqId} dictionary, 
  Tries to use seqCode as a mapping.
  input molSystem: ccp./molecule.MolSystem
        residToChemComp {Int resId: ccp.molecule.ChemComp} dictionary
  fitting to chain not in excludeChains
  Rasmus Fogh May 2013. Works for any combination of ChemComps, 
  also for incompletely specified sequences.
  
  Experimental - used in dataIo.DataMapper
  """
  
  # set up
  inputSeq = sorted(residToChemComp.items())
  
  for chain in molSystem.sortedChains():
    if not (excludeChains and chain in excludeChains):
      msResidues = chain.sortedResidues()
      seqCodeIndex = dict((x.seqCode,ii) for ii,x in enumerate(msResidues))
      ccSequence = [x.molResidue.chemComp for x in msResidues]
      
      if len(seqCodeIndex) == len(msResidues):
        # All seqCodes unique. Try mapping based on seqCode
        for resId,chemComp in inputSeq:
          if chemComp is not ccSequence[seqCodeIndex[resId]]:
            break
        else:
          # seqCode mapping matches. Use it.
          return (chain, dict((seqCode,msResidues[indx].seqId)
                               for seqCode,indx in seqCodeIndex.items()))
  #
  return None, {}

    
def matchByOffset(molSystem, residToChemComp, excludeChains=None):
  """ return chain, {resId:msResidue.seqId} dictionary, 
  Aligns to seqId assuming a gap-free offset.
  input molSystem: ccp./molecule.MolSystem
        residToChemComp {Int resId: ccp.molecule.ChemComp} dictionary
  fitting to chain not in excludeChains
  NBNB TODO now assumes resids fit within single chain
  Rasmus Fogh May 2013. Works for any combination of ChemComps, 
  also for incompletely specified sequences.
  
  Experimental - used in dataIo.DataMapper
  """
  
  # set up
  inputSeq = sorted(residToChemComp.items())
  minId = inputSeq[0][0]
  maxId = inputSeq[-1][0]
  
  for chain in molSystem.sortedChains():
    if not (excludeChains and chain in excludeChains):
      ccSequence = dict((x.seqId, x.molResidue.chemComp) 
                        for x in chain.sortedResidues())
                      
      for offset in range(1-minId, len(ccSequence) - maxId):
        for resId,chemComp in inputSeq:
          if chemComp is not ccSequence[resId+offset]:
            break
        else:
          # found a match
          break
      else:
        # no match, try next chain
        continue
      #
      return (chain, dict((seqId-offset, seqId) for seqId in ccSequence))
  #
  return None, {}


def matchResidueNumbers(molSystem, resIds, excludeChains=None):
  """ Map resIds to a chain in MolSystem assuming either seqCode or seqId 
  numbering matches
  molSystem: ccp.molecule.MolSystem to link to
  resIds: collection of Int residue ID
  excludeChains: chains that should not be used.
  Rasmus Fogh May 2013. For otherwise hopeless case where there are no resNames.
  
  Experimental - used in dataIo.DataMapper
  """
  
  residSet = set(resIds)
  
  for chain in molSystem.sortedChains():
    if not (excludeChains and chain in excludeChains):
      msResidues = chain.sortedResidues()
      
      if residSet.issubset(x.seqCode for x in msResidues):
        # assume resIds correspond to seqCodes
        return (chain, dict((x.seqCode,x.seqId) for x in msResidues))
        
      elif residSet.issubset(x.seqId for x in msResidues):
        # assume resIds corresond to seqIds
        return (chain, dict((x.seqId,x.seqId) for x in msResidues))
  #
  return (None, {})
  
  
def matchSequences(molSystem, resIdToChemComp, chains=None):
  """ Wrapper routine for TJS findMatchingChains
  NB can handle gaps, but is designed for pretty complete sequences
  mostly with standard residues.
  
  Experimental - used in dataIo.DataMapper
  """
  
  # Initialize sequence from 1 to max resId + 20.
  # That should get a map including missing residues and terminal deletions, 
  # as far as possible, even with lots of residues missing.
  # The alignment should bring us back to sense anyway
  # Of course with sparse residues present and gaps in the match 
  # the results may still be imperfect.
  maxResId = max(resIdToChemComp.keys())
  ccpCodes = ['???'] * (maxResId + 20)
  for ii, chemComp in resIdToChemComp.items():
    ccpCodes[ii] = chemComp.ccpCode
  
  # get the mapping
  chains, mappings = findMatchingChains(molSystem, ccpCodes,
                                                      excludeChains=chains)
  if mappings:
    chain = chains[0]
    mapping = mappings[0]
    seqMap = dict((tt[0], tt[1].seqId) for tt in mappping.items() 
                  if None not in tt)
    return chain, seqMap
  #
  return None, {}
  
