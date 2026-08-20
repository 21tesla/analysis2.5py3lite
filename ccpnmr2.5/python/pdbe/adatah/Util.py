import sys

from subprocess import Popen
from time import sleep

from pdbe.adatah.Constants import pythonCommand, numCpu

def runConversionJobs(pdbCodes,addOptions,scriptName,useNumCpu = None, isPython = True, sleepTime = 3):
    
  startCode = 'start'

  if not '-noGui' in addOptions:
    addOptions.append('-noGui')

  if '-reverse' in addOptions:
    addOptions.pop(addOptions.index('-reverse'))
    pdbCodes.reverse()

  currentPdbCode = None
  currentProcesses = {startCode: None}
  currentIndex = -1
  endPdbCode = pdbCodes[-1]
  
  if not useNumCpu:
    useNumCpu = numCpu
  
  scriptArgs = []
  if isPython:
    scriptArgs.append(pythonCommand)
  scriptArgs.append(scriptName)

  outputHandle = sys.__stdout__

  while (currentProcesses):

    if startCode in currentProcesses.keys():
      del(currentProcesses[startCode])

    if len(currentProcesses.keys()) < useNumCpu:

      tempIndex = currentIndex + 1
      for i in range(currentIndex,currentIndex + useNumCpu - len(currentProcesses.keys())):
        # Don't start a job if it's at the end!
        if currentPdbCode != endPdbCode:

          pdbCode = pdbCodes[tempIndex]

          #
          # TODO: ensure that stdout and stdin OK for this job!! Might want to reroute!!
          #

          process = Popen(['nice', '-19'] + scriptArgs + [pdbCode] + addOptions)

          currentProcesses[pdbCode] = process
          currentPdbCode = pdbCode

          outputHandle.write("\n*** Job %s started ***\n\n" % pdbCode)

          tempIndex += 1

          # Allow jobs to start one by one... nicer for output
          sleep(sleepTime)

      currentIndex = tempIndex - 1

    sleep(sleepTime)

    #
    # Check if finished
    #

    for pdbCode in currentProcesses.keys():

      # Finished...
      if currentProcesses[pdbCode].poll() != None:
        del(currentProcesses[pdbCode])
        outputHandle.write("\n *** Job %s finished ***\n" % pdbCode) 
        if not currentProcesses:
          currentProcesses = {startCode: None}

    outputHandle.flush()

#
# Sequence related
#

class AlignNeedlemanWunsch:

  def __init__(self,seq1,seq2,output = None,verbosity=1,compareToUncertain=False):

    """
    Code taken and adapted from http://snippets.dzone.com/posts/show/3387

    The Needleman-Wunsch algorithm preforms a global alignment 
    off two sequences (of length n and m) 
    For a given similarity matrix s
    (containig the penalties for character match-mismatch) 
    and a LINEAR gap penalty the algorithm is guaranteed 
    to find the alignment with highest score (in O(nm) time).   
    The algorithm is outlined through comments to the source.
    """
    
    self.output = output
    if not output:
      self.output = sys.stdout

    if seq1.count('X') == len(seq1) or seq2.count("X") == len(seq2):
      # NBNB TBD remove seq output again?
      self.output.write('Error: one sequence only consists of X, aborting.|%s|%s|\n'
      % (seq1, seq2))
      self.cols = self.rows = 0
      self.seq1 = self.seq2 = self.aseq1 = self.aseq2 = ""
      self.align = None
 
    else:
      self.verbosity = verbosity
      if verbosity:
        self.output.write('Calculating alignment...\n')
  
      self.rows=len(seq1)+1
      self.cols=len(seq2)+1
      
      self.seq1 = seq1
      self.seq2 = seq2
  
      try:
        #use fast numerical arrays if we can
        from numpy import zeros
        self.align=zeros((self.rows,self.cols),float)
      except ImportError:
        #use a list if we have to
        self.align=[]
        for i in range(self.rows):
          self.align+=[[0.]*self.cols]
  
      #################################################
      ##              Needleman-Wunsch               ##
      #################################################
  
      self.match=1.
      self.halfMatch = 0.5
      self.mismatch=-1.
      
      # Not sure if this is going to be OK overall...
      if not compareToUncertain:
        self.gap=-1.
      else:
        self.gap=-0.1
  
      # Set up and modify similarity matrix (ADAPTED WIM!)
      self.similarity = {}
      chars = "ABCDEFGHIJKLMNOPQRSTUVWYZ"
      
      if not compareToUncertain:
        chars += 'X'
      
      for char1 in chars:
        for char2 in chars:
          if char1 == char2:
            self.similarity[char1 + char2] = self.match
            
            if compareToUncertain:
              self.similarity[char1 + char2.lower()] = self.halfMatch
              self.similarity[char1.lower() + char2] = self.halfMatch
                        
          else:
            self.similarity[char1 + char2] = self.mismatch
            
            if compareToUncertain:
              self.similarity[char1 + char2.lower()] = self.mismatch
              self.similarity[char1.lower() + char2] = self.mismatch
      
      # No penalty for X?
      if compareToUncertain:
        for char1 in chars:
          self.similarity[char1 + 'X'] = 0.
          self.similarity['X' + char1] = 0.
          self.similarity[char1.lower() + 'X'] = 0.
          self.similarity['X' + char1.lower()] = 0.
  
        self.similarity['XX'] = 0.5
        self.similarity['tV'] = 0. # Reset, t is possibly a Val if a CG1 is missing somehow.
        self.similarity['Vt'] = 0. # Reset, t is possibly a Val if a CG1 is missing somehow.
        
      else:
        self.similarity['XX'] = 0.
        
      self.similarity['CM'] = 0.
      self.similarity['MC'] = 0.
  
      # Original matrix
      #s={
      #'AA':     match,'AG':mismatch,'AC':mismatch,'AT':mismatch,\
      #'GA':mismatch,'GG':     match,'GC':mismatch,'GT':mismatch,\
      #'CA':mismatch,'CG':mismatch,'CC':     match,'CT':mismatch,\
      #'TA':mismatch,'TG':mismatch,'TC':mismatch,'TT':     match,\
      #}
  
      for i in range(self.rows):
        self.align[i][0] = 0
      for j in range(self.cols):
        self.align[0][j] = 0
      for i in range(1,self.rows):
        for j in range(1,self.cols):
          # Dynamic programing -- aka. divide and conquer:
          # Since gap penalties are linear in gap size
          # the score of an alignmet of length l only depends on the   
          # the l-th characters in the alignment (match - mismatch - gap)
          # and the score of the one shorter (l-1) alignment,
          # i.e. we can calculate how to extend an arbritary alignment
          # soley based on the previous score value.  
          choice1 = self.align[i-1][j-1] + self.similarity[(seq1[i-1] + seq2[j-1])]
          choice2 = self.align[i-1][j] + self.gap
          choice3 = self.align[i][j-1] + self.gap
          
          #if i == j: # Or whatever alignment you're looking for...
          #  print i, j, choice1, choice2, choice3, (seq1[i-1],seq2[j-1])
          
          self.align[i][j] = max(choice1, choice2, choice3)
  
  
      aseq1 = ''
      aseq2 = ''
      #We reconstruct the alignment into aseq1 and aseq2, 
      i = len(seq1)
      j = len(seq2)
      while i>0 and j>0:
        if i%10==0 and verbosity:
          self.output.write('.')
  
        #by preforming a traceback of how the matrix was filled out above,
        #i.e. we find a shortest path from a[n,m] to a[0,0]
        score = self.align[i][j]
        score_diag = self.align[i-1][j-1]
        score_up = self.align[i][j-1]
        score_left = self.align[i-1][j]
        
        #print i, j, i-j, score, score_diag, score_up, score_left
        
        if score == score_diag + self.similarity[seq1[i-1] + seq2[j-1]]:
          aseq1 = seq1[i-1] + aseq1
          aseq2 = seq2[j-1] + aseq2
          i -= 1
          j -= 1
        elif score == score_left + self.gap:
          aseq1 = seq1[i-1] + aseq1
          aseq2 = '_' + aseq2
          i -= 1
        elif score == score_up + self.gap:
          aseq1 = '_' + aseq1
          aseq2 = seq2[j-1] + aseq2
          j -= 1
        else:
          #should never get here..
          print 'ERROR'
          i=0
          j=0
          aseq1='ERROR';aseq2='ERROR';seq1='ERROR';seq2='ERROR'
      while i>0:
        #If we hit j==0 before i==0 we keep going in i.
        aseq1 = seq1[i-1] + aseq1
        aseq2 = '_' + aseq2
        i -= 1                
  
      while j>0:
        #If we hit i==0 before i==0 we keep going in j. 
        aseq1 = '_' + aseq1
        aseq2 = seq2[j-1] + aseq2
        j -= 1
          
      self.aseq1 = aseq1
      self.aseq2 = aseq2
      
      if verbosity:
        self.output.write("\n")
    
  #################################################
  #################################################
  ##              Full backtrack                 ##
  #################################################

  #To reconstruct all alinghments is somewhat tedious..
  def make_graph(self):
  #the simpilest way is to make a graph of the possible constructions of the values in a 
    graph={}
    for i in range(1,self.cols)[::-1]:
      graph[(i,0)] = [(i-1,0)]
      graph[(0,i)] = [(0,i-1)]
      for j in range(1,self.cols)[::-1]:
        graph[(i,j)]=[]
        score = self.align[i][j]
        score_diag = self.align[i-1][j-1]
        score_up = self.align[i][j-1]
        score_left = self.align[i-1][j]
        if score == score_diag + self.similarity[self.seq1[i-1] + self.seq2[j-1]]:
          graph[(i,j)] += [(i-1,j-1)]
        if score == score_left + self.gap:
          graph[(i,j)] += [(i-1,j)]
        if score == score_up + self.gap:
          graph[(i,j)] += [(i,j-1)]

    self.graph = graph

  def searchPaths(self):
  
    self.tracks = self.find_all_paths((self.cols-1,self.rows-1),(0,0))

  def find_all_paths(self, start, end, path=[], depth = 0):
  
    if depth == 30:
      if self.verbosity:
        print "  ERROR: recursive loop problem. Aborting at depth 30."
      return []
      
  #and then to recursivly find all paths 
  #from bottom right to top left..
    path = path + [start]
  #    print start
    if start == end:
      return [path]
    if not self.graph.has_key(start):
      return []
    paths = []
    for node in self.graph[start]:
      if node not in path:
        newpaths = self.find_all_paths(node, end, path,depth = depth + 1)
        for newpath in newpaths:
          paths.append(newpath)    
          
    return paths
  
  def printScores(self,verbosity = 2):
  
    baseqs1=[]
    baseqs2=[]
    for track in self.tracks:
    #using these we can reconstruct all optimal alig.-s 
      baseq1 = ''
      baseq2 = ''
      last_step=(self.cols-1,self.rows-1)
      for step in track:
        i,j=last_step
        if i==step[0]:
          baseq1 = '_' + baseq1
          baseq2 = self.seq2[j-1] + baseq2
        elif j==step[1]:
          baseq1 = self.seq1[i-1] + baseq1
          baseq2 = '_' + baseq2
        else:
          baseq1 = self.seq1[i-1] + baseq1
          baseq2 = self.seq2[j-1] + baseq2

        last_step=step
      baseqs1+=[baseq1]
      baseqs2+=[baseq2]
    #################################################
    
    if verbosity > 1:

      print ''
      print '# Using:  match='+repr(self.match)+'; mismatch='+repr(self.mismatch)+'; gap='+repr(self.gap)              
      print self.seq1
      print self.seq2
      print '# We get e.g.:'
      print self.aseq1
      print self.aseq2
      print ''
      
    gaps=0
    mms=0
    ms=0
    for i in range(len(self.aseq1)):
      if self.aseq1[i]==self.aseq2[i]:
        self.aseq1=self.aseq1[:i]+'='+self.aseq1[i+1:]
        self.aseq2=self.aseq2[:i]+'='+self.aseq2[i+1:]
        ms+=1
      else:
        if self.aseq1[i]=='_' or self.aseq2[i]=='_':
          gaps+=1
        else:
          mms+=1

    if verbosity > 0:
      print self.aseq1
      print self.aseq2
      print ''
      print ms,' matches; ',mms,' mismatches; ',gaps,' gaps.' 
      print '# With a score of'
      print self.align[self.rows-2][self.cols-2],'/',min(len(self.seq1),len(self.seq2))

      print 'Optimal alignment is ',len(self.tracks),' times degenerate:'
      print ''
      for i in range(len(self.tracks)):
        print i+1,'.'
        print baseqs1[i]
        print baseqs2[i]

  def getBestMatchInfo(self,verbosity = 0):
  
    self.make_graph()
    self.searchPaths()
    
    baseqs1=[]
    baseqs2=[]
    for track in self.tracks:
    #using these we can reconstruct all optimal alig.-s 
      baseq1 = ''
      baseq2 = ''
      last_step=(self.cols-1,self.rows-1)
      for step in track:
        i,j=last_step
        if i==step[0]:
          baseq1 = '_' + baseq1
          baseq2 = self.seq2[j-1] + baseq2
        elif j==step[1]:
          baseq1 = self.seq1[i-1] + baseq1
          baseq2 = '_' + baseq2
        else:
          baseq1 = self.seq1[i-1] + baseq1
          baseq2 = self.seq2[j-1] + baseq2

        last_step=step
      baseqs1+=[baseq1]
      baseqs2+=[baseq2]
    #################################################
      
    gaps=0
    mismatches=0
    matches=0
    for i in range(len(self.aseq1)):
      if self.aseq1[i]==self.aseq2[i]:
        self.aseq1=self.aseq1[:i]+'='+self.aseq1[i+1:]
        self.aseq2=self.aseq2[:i]+'='+self.aseq2[i+1:]
        matches+=1
      else:
        if self.aseq1[i]=='_' or self.aseq2[i]=='_':
          gaps+=1
        else:
          mismatches+=1

    if verbosity > 0:
      print self.aseq1
      print self.aseq2
      print ''
    
    
    if self.align is not None:
      score = self.align[self.rows-2][self.cols-2]
      total = min(len(self.seq1),len(self.seq2))
      degeneracy = len(self.tracks)
    else:
      score = total = degeneracy = 0
    
    return (score, total, degeneracy, matches, mismatches, gaps,(self.aseq1,self.aseq2))

def getAlignmentInfo(seqAlignInfo):

  """
  Takes the sequence results from the Needleman-Wunsch algorithm (last value returned by getBestMatchInfo() method)
  and return a list of seqDifferences that can be used for forceChainMappings.
  
  WARNING: is *swapping* the output, PDB first, then BMRB! Need to sort this out!!
  
  """
          
  #
  # The PDB sequence is the reference one used to create the molecular system! Always starts at 1
  #

  seqDifferences = []

  pdbIndex = bmrbIndex = 1

  (bmrbAlignInfo,pdbAlignInfo) = seqAlignInfo
  
  totalAlignLen = len(bmrbAlignInfo)

  for i in range(totalAlignLen):
  
    bmrbStatus = bmrbAlignInfo[i]
    pdbStatus = pdbAlignInfo[i]

    if bmrbStatus == '_':
      seqDifferences.append((pdbIndex,pdbStatus,bmrbIndex,None))
      pdbIndex += 1
    elif pdbStatus == '_':
      seqDifferences.append((pdbIndex,None,bmrbIndex,bmrbStatus))
      bmrbIndex += 1

    else:
      if pdbStatus != '=':
        seqDifferences.append((pdbIndex,pdbStatus,bmrbIndex,bmrbStatus))

      pdbIndex +=1
      bmrbIndex += 1

  return seqDifferences

#
# CCPN related
#

from memops.api import Implementation

def copyOriginalEntryInfoToNewEntry(originalEntry,newEntry):

  #
  # Copy over general information
  #

  molSystem = originalEntry.molSystem
  measurementLists = originalEntry.measurementLists

  newEntry.molSystem = molSystem
  newEntry.measurementLists = measurementLists

  #
  # Copy all application data from the old entry...
  #

  for applDatum in originalEntry.applicationData:
    applDataClass = getattr(Implementation,applDatum.className)
    newEntry.addApplicationData(applDataClass(value=applDatum.value, application=applDatum.application, keyword=applDatum.keyword))

  return newEntry

#
# 
# duplicateResonances.py: Duplicates a set of FixedResonances (e.g. for other chain in homodimer)
#
# TODO IS THIS THE RIGHT PLACE FOR THIS CODE?!?!
#

#
# General API imports...
#

from ccp.api.nmr import NmrConstraint
from ccp.general.Copy import copyAttributeInfo, copyResonanceInfo
from ccpnmr.format.general.Util import getNameInfo, getResName
from ccpnmr.format.general.Constants import assign_kw
from ccpnmr.format.general.Util import updateResonanceNamesDict, getResNameText

def getMappedResonances(constraintOrItem,resonanceMappingDict,newChainCode):

  newResonances = []
  oldToNewResonanceDict = {}

  for resonance in constraintOrItem.resonances:
    newResonance = resonance
    if resonanceMappingDict.has_key(resonance):
      if resonanceMappingDict[resonance].has_key(newChainCode):
        newResonance = resonanceMappingDict[resonance][newChainCode]
        oldToNewResonanceDict[resonance] = newResonance

    newResonances.append(newResonance)

  return (newResonances, oldToNewResonanceDict)

def duplicateResonances(nmrConstraintStore,format,mappingDict):
  
  #
  # Check which resonances have to be copied, and copy...
  #
  
  resonanceMappingDict = {}
  
  copyResonances = []
  resonanceDict = {}

  copyConstraints = []
  constrLinks = ['dihedralConstraints','chemShiftConstraints']

  copyConstraintItems = {}
  copyConstraintItemsList = []
  
  #
  # First determine which resonances are eligible for copying, and create a dictionary with resNames...
  #
  
  for resonance in nmrConstraintStore.sortedFixedResonances():
         
    applData = resonance.findFirstApplicationData(application = format, keyword = assign_kw)

    if applData:

      resName = applData.value
      
      resonanceDict[resName] = resonance

      (chain,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resName)
              
      if mappingDict.has_key(chain):
        
        copyResonances.append(resonance)
        
  #
  # Then determine which constraint/items have to be copied (all resonances have to have the same chainCode!!)
  #
  
  for resonance in copyResonances:
    
    #
    # Constraints
    #

    for constrLink in constrLinks:
      for constr in getattr(resonance,constrLink):
        copyConstraint = 1
        for otherRes in constr.resonances:
          if otherRes not in copyResonances:
            copyConstraint = 0
            break
        
        if copyConstraint and constr not in copyConstraints:
          j = 0
          constrKey = (constr.parent.serial,constr.serial)
          for j in range(0,):
            if constrKey < (copyConstraints[j].parent.serial,copyConstraints[j].serial):
              copyConstraints.insert(j,constr)
              break
          listLen = len(copyConstraints)
          if j == (listLen - 1) or listLen == 0:
            copyConstraints.append(constr)
            

    #
    # ConstraintItems
    #

    for constrItem in resonance.sortedPairwiseConstraintItems():

      copyConstraintItem = 1
      for otherRes in constrItem.resonances:
        if otherRes not in copyResonances:
          copyConstraintItem = 0
          break

      if copyConstraintItem:
        if not copyConstraintItems.has_key(constrItem.constraint):
          copyConstraintItems[constrItem.constraint] = []
          j = 0

          constrKey = (constrItem.constraint.parent.serial,constrItem.constraint.serial)

          for j in range(0,len(copyConstraintItemsList)):
            if constrKey < (copyConstraintItemsList[j].parent.serial,copyConstraintItemsList[j].serial):
              copyConstraintItemsList.insert(j,constrItem.constraint)
              break
          listLen = len(copyConstraintItemsList)
          #print constrKey, j, listLen
          if j == (listLen - 1) or listLen == 0:
            copyConstraintItemsList.append(constrItem.constraint)

        if constrItem not in copyConstraintItems[constrItem.constraint]:
          copyConstraintItems[constrItem.constraint].append(constrItem)

  #
  # Now copy the resonances or find an existing resonance...
  #
  
  for resonance in copyResonances:
         
    applData = resonance.findFirstApplicationData(application = format, keyword = assign_kw)
    applDataClass = applData.__class__

    resName = applData.value

    (chain,seqCode,spinSystemId,seqInsertCode,atomName) = getNameInfo(resName)
      
    resonanceMappingDict[resonance] = {}
    
    resonance.removeApplicationData(applData)
    
    newValue = getResName(mappingDict[chain][0],seqCode,atomName,seqInsertCode = seqInsertCode)
    resonance.addApplicationData(applDataClass(application = format, keyword = assign_kw, value = newValue))

    for i in range(1,len(mappingDict[chain])):

      newChainCode = mappingDict[chain][i]
      newResName = getResName(newChainCode,seqCode,atomName,seqInsertCode = seqInsertCode)
      
      if resonanceDict.has_key(newResName):
        newResonance = resonanceDict[newResName]
        print "  Using existing resonance %s..." % newResName

      else:

        print "  Creating new resonance %s..." % newResName

        newResonance = nmrConstraintStore.newFixedResonance(isotopeCode = resonance.isotopeCode)
        resonanceDict[newResName] = newResonance

        copyResonanceInfo(resonance,newResonance,toResName = newResName)

        applData = newResonance.findFirstApplicationData(application = format, keyword = assign_kw)
        applDataClass = applData.__class__

        newResonance.removeApplicationData(applData)
        newResonance.addApplicationData(applDataClass(application = format, keyword = assign_kw, value = newResName))

        #
        # Constraints and constraint items have to be deleted and recreated later...
        #
        
        for constrLink in constrLinks:
          for constr in getattr(newResonance,constrLink):
            constr.delete()

        for constrItem in newResonance.pairwiseConstraintItems:
          constrItem.delete()

      resonanceMappingDict[resonance][newChainCode] = newResonance

  #
  # Loop over the 'new' chain codes...
  #
  
  print
  print "########################## "
  print "# Duplicating resonances # "
  print "########################## "
  print
  
  for i in range(1,len(mappingDict[chain])):

    newChainCode = mappingDict[chain][i]
  
    #
    # Copy the relevant dihedral constraints...
    #

    for constraint in copyConstraints:
    
      print "  Copying constraint %d.%d" % (constraint.parent.serial,constraint.serial)
      print "    Resonances %s" % ', '.join([getResNameText(res) for res in constraint.resonances])

      (newResonances, oldToNewResonanceDict) = getMappedResonances(constraint,resonanceMappingDict,newChainCode)
      
      print "            to %s" % ', '.join([getResNameText(res) for res in newResonances])

      newConstr = constraint.__class__(constraint.parent,resonances = newResonances)
      copyAttributeInfo(constraint,newConstr)

      for constrItem in constraint.sortedItems():

        newConstrItem = constrItem.__class__(newConstr, upperLimit = 1.0, lowerLimit = 0.0)
        copyAttributeInfo(constrItem,newConstrItem)
      
      print

    #
    # ..and the constraints from the relevant constraintItems...
    #

    for constraint in copyConstraintItemsList:
   
      print "  Copying constraint %d.%d" % (constraint.parent.serial,constraint.serial)

      newConstr = constraint.__class__(constraint.parent)
      copyAttributeInfo(constraint,newConstr)
      
      if hasattr(constraint,'weight'):
        newConstr.weight = constraint.weight

      for constrItem in constraint.sortedItems():
        print "    Copying item with resonances %s" %  ', '.join([getResNameText(res) for res in constrItem.resonances])
      
        if constrItem in copyConstraintItems[constraint]:

          (newResonances, oldToNewResonanceDict) = getMappedResonances(constrItem,resonanceMappingDict,newChainCode)
        
        else:
        
          newResonances = constrItem.resonances
          oldToNewResonanceDict = {}

      
        print "                              to %s" % ', '.join([getResNameText(res) for res in newResonances])

        newConstrItem = constrItem.__class__(newConstr,resonances = newResonances)

        if constrItem.className in ('DistanceConstraintItem', 'HBondConstraintItem', 'JCouplingConstraintItem', 'RdcConstraintItem'):
          firstResonance = constrItem.firstResonance
          if oldToNewResonanceDict:
            firstResonance = oldToNewResonanceDict[firstResonance] 

          newConstrItem.firstResonance = firstResonance
 
        copyAttributeInfo(constrItem,newConstrItem)
      
      print

  print "###################### "
  print "# End of duplication # "
  print "###################### "
  print
