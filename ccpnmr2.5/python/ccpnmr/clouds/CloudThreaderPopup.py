"""
======================COPYRIGHT/LICENSE START==========================

CloudThreaderPopup.py: Part of the CcpNmr Clouds program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

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
A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
import Tkinter, re

from math import sqrt

from os import listdir

from memops.general import Implementation

from memops.editor.BasePopup    import BasePopup

from memops.gui.CheckButton     import CheckButton
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.IntEntry        import IntEntry
from memops.gui.Entry           import Entry
from memops.gui.Label           import Label
from memops.gui.PulldownMenu    import PulldownMenu
from memops.gui.Util            import createDismissHelpButtonList
from memops.gui.ScrolledGraph   import ScrolledGraph
from memops.gui.ProgressBar     import ProgressBar

#from ccpnmr.clouds.CloudThreader import cloudThreader, scoreAssignment, fitSideChains
from ccpnmr.clouds.CloudThreader import makeStructureDictFromXml
from ccpnmr.clouds.CloudBasic   import getFileNamesFromPattern, getCloudsFromFile
from ccpnmr.clouds.CloudBasic   import code1LetterDict

from ccpnmr.analysis.core.AssignmentBasic import getResonanceAtomTuple, assignSpinSystemResidue

from memops.gui.MessageReporter import showWarning, showYesNo

from ccpnmr.analysis.core.SpinSystemTyping import getSpinSystemScore, getSpinSystemShifts

def cloudThreaderMacro(argServer):

  popup = CloudThreaderPopup(argServer.parent)
  popup.open()
  
globalIndex = 0

sideChainAtoms = ['HA','HA1','HA2','HB','HB1','HB2','HB3',
                  'HG','HG2','HG3','HG11','HG12','HG13',
                  'HG21','HG22','HG23','HD1','HD2','HD3',
                  'HD11','HD12','HD13','HD21','HD22','HD23',
                  'HE','HE1','HE2','HE3','HE21','HE22','HZ',
                  'HZ1','HZ2','HZ3','HH','HH2']


class CloudThreaderPopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.guiParent = parent
    self.project   = parent.getProject()
    self.molSystem = None
    self.chain     = None
    self.assignment = None
    self.scores     = []

    BasePopup.__init__(self, parent, title="Cloud Threader", **kw)
  
  def body(self, guiFrame):

    guiFrame.grid_columnconfigure(3, weight=1)
    
    row = 0
    label = Label(guiFrame, text='Molecular system: ')
    label.grid(row=row, column=0, sticky=Tkinter.NW)
    self.molSysPulldown = PulldownMenu(guiFrame, self.changeMolSystem, selected_index=-1, do_initial_callback=0)
    self.molSysPulldown.grid(row=row, column=1, sticky=Tkinter.NW)

    label = Label(guiFrame, text='File name: ')
    label.grid(row=row, column=2, sticky=Tkinter.NW)
    self.filenameEntry = Entry(guiFrame,text='nexus15.pdb_\d+')
    self.filenameEntry.grid(row=row, column=3, sticky=Tkinter.NW)


    row += 1
    label = Label(guiFrame, text='Chain: ')
    label.grid(row=row, column=0, sticky=Tkinter.NW)
    self.chainPulldown = PulldownMenu(guiFrame, self.changeChain, selected_index=-1, do_initial_callback=0)
    self.chainPulldown.grid(row=row, column=1, sticky=Tkinter.NW)

    label = Label(guiFrame, text='Search steps: ')
    label.grid(row=row, column=2, sticky=Tkinter.NW)
    self.numStepsEntry = IntEntry(guiFrame,text=5000)
    self.numStepsEntry.grid(row=row, column=3, sticky=Tkinter.NW)

    row += 1

    label = Label(guiFrame, text='Keep existing assignments?')
    label.grid(row=row, column=0, sticky=Tkinter.NW)
    self.keepAssignSelect = CheckButton(guiFrame, callback=None)
    self.keepAssignSelect.grid(row=row, column=1, sticky=Tkinter.NW)
    self.keepAssignSelect.set(0)

    label = Label(guiFrame, text='Assignment Threshold: ')
    label.grid(row=row, column=2, sticky=Tkinter.NW)
    self.thresholdEntry = FloatEntry(guiFrame,text=-4.9)
    self.thresholdEntry.grid(row=row, column=3, sticky=Tkinter.NW)
    row += 1

    label = Label(guiFrame, text='Global score: ')
    label.grid(row=row, column=0, sticky=Tkinter.NW)
    self.globalScoreLabel = Label(guiFrame, text='')
    self.globalScoreLabel.grid(row=row, column=1, sticky=Tkinter.NW)

    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    self.graph = ScrolledGraph(guiFrame, width=300, height=200, graphType='scatter', dataColors=['#008000','#F00000'])
    self.graph.grid(row=row, column=0, columnspan=4, sticky = Tkinter.NSEW)

    row += 1
    texts    = ['Run','Assign!','Score Current','Fit Side Chains']
    commands = [self.run, self.assignSpinSystems,
                self.scoreCurrentAssignments, self.assignSideChains]
    bottomButtons = createDismissHelpButtonList(guiFrame,texts=texts,commands=commands,expands=0,help_url=None)
    bottomButtons.grid(row=row, column=0, columnspan=4, sticky=Tkinter.EW)
    self.assignButton = bottomButtons.buttons[1]

    for func in ('__init__','delete'):
      Implementation.registerNotify(self.updateMolSystems, 'ccp.molecule.MolSystem.MolSystem', func)
      Implementation.registerNotify(self.updateChains, 'ccp.molecule.MolSystem.Chain', func)
    
    self.updateMolSystems()
    self.updateChains()

  def assignSideChains(self):

    if self.assignment and self.scores:
      print "Fitting side chains"
      pattern = self.filenameEntry.get()
      files     = getFileNamesFromPattern(pattern , '.')
      clouds    = getCloudsFromFile(files, self.chain.root)
      threshold = self.thresholdEntry.get() or -4.0
      fitSideChains(self.scores, self.assignment, clouds, threshold)

  def scoreCurrentAssignments(self):
  
    if self.chain:
      pattern = self.filenameEntry.get()
      files   = getFileNamesFromPattern(pattern , '.')
      clouds  = getCloudsFromFile(files, self.chain.root)
      score, self.scores, self.assignment = scoreAssignment(self.chain, clouds, spinSystems=None, graph=self.graph)
      self.update()

  def update(self):
  
    if self.assignment and self.scores:
      self.assignButton.enable()
    else:
      self.assignButton.disable()  

  def run(self):
  
    #import profile, pstats
    nSteps   = self.numStepsEntry.get()
    if not nSteps:
      self.scoreCurrentAssignments()

    elif self.chain:
      pattern  = self.filenameEntry.get()
      pgb      = ProgressBar(self, text='Searching', total=nSteps)
      files    = getFileNamesFromPattern(pattern , '.')
      preserve = self.keepAssignSelect.get()
      if not files:
        pgb.destroy()
        return

      clouds  = getCloudsFromFile(files, self.chain.root)
      score, self.scores, self.assignment = fastCloudThreader(self.chain, clouds, nSteps=nSteps, graph=self.graph,
                                                              preserveAssign=preserve, progressBar=pgb)

      pgb.destroy()
      self.globalScoreLabel.set(str(score))
      self.update()

      """
      self.guiParent.arguments = (self.chain, clouds,nSteps,self.graph,0,pgb)
      self.guiParent.func      = self.cloudThread
      profile.run('top.func(top.arguments)','cloudThreadOut.txt')
      def cloudThread(self, args):
 
        cloudThreader(args[0], args[1], nSteps=args[2], graph=args[3],
                      preserveAssign=args[4], progressBar=args[5])"""

  def assignSpinSystems(self):
   
    if self.assignment and self.scores:
      if showYesNo('Query','Are you sure?'):
        threshold = self.thresholdEntry.get() or -4.0
        i = 0
 
        for residue in self.assignment.keys():
          #if self.scores[residue] > threshold:
          spinSystem = self.assignment[residue]
          spinSystem.name = None
          #assignSpinSystemResidue(spinSystem,residue=None)
 
        for residue in self.assignment.keys():
          if self.scores[residue] > threshold:
            i += 1
            spinSystem = self.assignment[residue]
            name = '%d%s' % (residue.seqCode,residue.ccpCode)
            spinSystem.name = name
            print 'Assign %d%s - {%d) %s' % (residue.seqCode,residue.ccpCode,spinSystem.serial,name)
            #assignSpinSystemResidue(spinSystem,residue=residue)
      
      showWarning('Done','%d residues assigned' % i)
      
  def getMolSystems(self):
  
    names = []
    for molSystem in self.project.molSystems:
      if molSystem.chains:
        names.append( '%s' % (molSystem.code) )
    return names


  def changeMolSystem(self, i, name):
  
    self.molSystem = self.project.findFirstMolSystem(code=name)


  def updateMolSystems(self, *opt):
  
    names = self.getMolSystems()
    if names:
      if not self.molSystem:
        self.molSystem = self.project.findFirstMolSystem(code=names[0])
      self.molSysPulldown.setup(names, names.index(self.molSystem.code))


  def getChains(self):
  
    chains = []
    if self.molSystem:
      for chain in self.molSystem.chains:
        chains.append( [chain.code, chain] )
	
    return chains


  def changeChain(self, i, name=None):
    
    if not name:
      i = self.chainPulldown.selected_index
    
    chains = self.getChains()
    if chains:
      self.chain = chains[i][1]
    
    
  def updateChains(self, *chain):
  
    chains = self.getChains()
 
    if chains:
      names = [x[0] for x in chains]
      if (not self.chain) or (self.chain.code not in names):
        self.chain = chains[0][1]
      self.chainPulldown.setup(names, names.index(self.chain.code) )

    self.update()

  def destroy(self):

    for func in ('__init__','delete'):
      Implementation.unregisterNotify(self.updateMolSystems, 'ccp.molecule.MolSystem.MolSystem', func)
      Implementation.unregisterNotify(self.updateChains, 'ccp.molecule.MolSystem.Chain', func)

    BasePopup.destroy(self)

import cPickle
from os.path import exists, isfile, isdir
from os import listdir, path
from math import sqrt, log, exp
from random import randint, random
from memops.gui.MessageReporter import showWarning

from ccpnmr.analysis.core.SpinSystemTyping import getSpinSystemTypes
from ccpnmr.analysis.core.AssignmentBasic  import makeResonanceGuiName

#from ccpnmr.analysis.StructureBasic import makeStructureDictFromRoughPdb
olc = code1LetterDict

logDict = {}
nmrRefStoreDict = {}
chemAtomShiftDict = {}
global J
J = 0
global inter
global intra
inter = None

def scoreAssignment(chain, clouds, spinSystems=None, graph=None):


  global inter
  project  = chain.root
  if not spinSystems:
    spinSystems = list(project.currentNmrProject.resonanceGroups)
  
  assignment0 = {}
  for ss in spinSystems:
    ss.resonanceCoords = {}
    if ss.residue:
      assignment0[ss.residue] = ss
    
  residueScores = {}
  if inter is None:
    intra, inter = getDistDistributions()
  
  for residue in chain.residues:
    ss = assignment0.get(residue)
    if ss:
      ll = getResidueLikelihood(residue, ss, assignment0, clouds, inter, {}) or 0.0
    else:
      ll = -10.0
    residueScores[residue] = ll
  
  if graph:
    updateGraph(graph, chain.residues, assignment0, 0)
    
  return 1.0, residueScores, assignment0

def cloudThreader(chain, clouds, spinSystems=None, nSteps=3000, graph=None, preserveAssign=0, progressBar=None):

  project  = chain.root
  residues = list(chain.residues)

  if not spinSystems:
    spinSystems = list(project.currentNmrProject.resonanceGroups)

  for residue in residues:
    residue.spinSystemScoreDict = {}

  for ss in spinSystems:
    ss.resonanceCoords = {}
    
  score, scores, assignment = searchPosterior(project, residues, spinSystems,
                                              clouds, nSteps, graph=graph,
                                              preserveAssign=preserveAssign,
                                              progressBar=progressBar )

  return score, scores, assignment


def searchPosterior(project, residues, spinSystems, clouds, known=None, nSteps=100, inter=None, graph=None, verbose=0, preserveAssign=0, progressBar=None):

  chain = residues[0].chain
  for spinSystem in spinSystems:
    for resonance in spinSystem.resonances:
      if not resonance.resonanceSet:
        pass # resonance.delete()
      else:
        resonance.name == 'r'

  sequence = ''
  for residue in chain.residues:
    sequence += olc.get(residue.ccpCode, 'X')

  for residue in residues:
    residue.atomDict = {}
    for atom in residue.atoms:
      residue.atomDict[atom.name] = atom

  print "Generating intial assignment"
  
  assignment0 = getInitialAssignment(residues, spinSystems, known or {}, graph=graph)
  assignmentB = assignment0
    
  intra = {}
  if not inter:
    print "Reading distance distributions"
    intra, inter = getDistDistributions()

  p0 = generatePosterior(residues, assignment0, clouds, intra, inter, known, 0, preserveAssign)
  p  = p0
  pB = p0
  
  if progressBar:
    progressBar.parent.update_idletasks()
    progressBar.total = 100
    progressBar.label.set('Searching cloud')
    progressBar.set(0)
    progressBar.open()

  q = max(1, nSteps/100)
  
  if graph:
    updateGraph(graph, chain.residues, assignment0, 0)

  print "Generate posterior"
  for i in range(nSteps):
    t = i/float(nSteps) * 20

    #assignment = swapSpinSystem(assignment0, clouds, inter, known, 1, i)  
    assignment = swapSpinSystem2(residues, assignment0)  
    p = generatePosterior(residues, assignment, clouds, intra, inter, known, i, preserveAssign)

    success = 0
    if p > p0:
      success = 1
    
    else:
      delta = (p0 - p)/p0
      delta *= t + 50
      r = random()
     
      if delta and delta < r:
	success = 1
         
    if success: 
      p0 = p
      print p0
      assignment0 = assignment
      
      if p > pB:
        pB = p
        assignmentB = assignment
        if graph:
          updateGraph(graph, chain.residues, assignment0, i)
 
      
      if verbose:
        found = 0
        assign = ''
        typed  = ''
        foundCodes = ''
        scores = ''
        for residue in chain.residues:
          ss = assignment0[residue]
          if hasattr(residue, 'likelihood'):
            ll = residue.likelihood
          else:
            ll = 0.0
          score = int(ll+10)
          if score > 9:
            scores += '*'
          else:
            scores += str(score)
          if ss.residue:
            if ss.residue is residue:
              found +=1
              foundCodes += olc.get(residue.ccpCode, 'X')
              #if ll == -10.0:
              #  print residue.seqCode, residue.ccpCode
              #  getResidueLikelihood(residue, ss, assignment0, cloud, inter, known, v=1)
              
            else:
              foundCodes += '-'

            if ss.residue.ccpCode == residue.ccpCode:
              typed += olc.get(residue.ccpCode, 'X')
            else:
              typed += '-'

            assign += olc.get(ss.residue.ccpCode, 'X')
          else:
            assign += '-'
            foundCodes += '-'
            typed += '-'

        print i, p, p0
        print 'Found: %d' % found
        print scores
        print foundCodes
        print typed
        print sequence
        print assign
    
    if progressBar and (i % q == 0):
      progressBar.increment()
    
  residueScores = {}
  for residue in chain.residues:
    ss = assignment0.get(residue)
    if ss:
      ll = getResidueLikelihood(residue, ss, assignmentB, clouds, inter, known) or 0.0
      residueScores[residue] = ll
   
  return pB, residueScores, assignmentB

def updateGraph(graph, residues, assignment, step):

  title  = 'Clouds Threader Fitting Scores - step %d' % step
  xLabel = 'Residue'
  yLabel = 'Score'

  dataSets = [[], []]
  for residue in residues:
    ss = assignment.get(residue)
    if ss is not None:
      residue2 = ss.residue
      if residue2 is None:
        check = 1
      elif residue2 is residue:
        check = 0
      else:
        check = 1
 
      if hasattr(residue, 'likelihood'):
        ll = residue.likelihood
      else:
        ll = -10.0
      dataSets[check].append(  [residue.seqCode, ll] )
      #dataSet2.append( [residue.seqCode, ll] )
 

  graph.update(dataSets=dataSets,xLabel=xLabel,yLabel=yLabel,title=title)
  graph.parent.update_idletasks()

def printResidueScores(residues, assignment, clouds, interDistrib, known):

  for residue in residues:
    spinSystem = assignment[residue]
    l = getResidueLikelihood(residue, spinSystem, assignment, clouds, interDistrib, known) or 0.0
    s = getSpinSystemScore(residue.ccpCode, spinSystem) or 0.0
    print "%d %s %f %f" % (residue.seqCode,residue.ccpCode, s, l)
 
 
def generatePosterior(residues, assignment, clouds, intra, inter, known, i, preserveAssign=0):

  
  prior = generatePrior(residues, assignment)
  likelihood = generateLikelihood(residues, assignment, clouds, intra, inter, known, preserveAssign)

  p = prior * likelihood 
  # no normalisation (yet)
  
  return p


def generatePrior(residues, assignment):

  # how well do the shifts of the resonances in a spin system match a residue
  # could the shift possibly be in the spin system
  # no conditional dependance in the intial instance

  score = 0.0
  N = 0

  for residue in residues:
    spinSystem = assignment[residue]
    score2 = getSpinSystemScore(residue.ccpCode, spinSystem)
    
    N += 1
    score += score2

  if N > 0:
    return exp(score/len(residues))

  else:
    return 0.0


def getResidueLikelihood(residue1, spinSystem1, assignment, clouds, interDistrib, known=None, v=0, floor=0.0025):

  if not known:
    known = {}

  chain   = residue1.chain
  amide1  = getResonanceCoords(clouds, spinSystem1, 'H') or getResonanceCoords(clouds, spinSystem1, 'HD2')
  alpha1  = getResonanceCoords(clouds, spinSystem1, 'HA') or getResonanceCoords(clouds, spinSystem1, 'HA2')
  betas1  = getResonanceCoords(clouds, spinSystem1, 'HB')
  coords1 = [amide1,alpha1,betas1]
  coords2 = None 
     
  if not hasattr(chain, 'lookupSeqID'):
    chain.lookupSeqID = {}
    for residue in chain.residues:
      chain.lookupSeqID[residue.seqId] = residue
  lookupSeqID = chain.lookupSeqID
  
  r2 = []
  
  N = 0
  p = 0.0
  badPrev = 0
  #for delta in (-2,-1,1,2):
  for delta in (-3,-2,-1,1,2,3):
    #for delta in (-1,1):
    p2 = 0.0
    N2 = 0
    residue2 = lookupSeqID.get(residue1.seqId+delta)
    if residue2 is None:
      continue
  
    spinSystem2 = assignment.get(residue2)
    if spinSystem2 is None:
      spinSystem2 = known.get(residue2)
      if spinSystem2 is None:
        continue
  
    r2.append( '%d%s' % (residue2.seqCode, residue2.ccpCode))
  
    dist = []
    
    amide2  = getResonanceCoords(clouds, spinSystem2, 'H') or getResonanceCoords(clouds, spinSystem2, 'HD2')
    alpha2  = getResonanceCoords(clouds, spinSystem2, 'HA') or getResonanceCoords(clouds, spinSystem1, 'HA2')
    betas2  = getResonanceCoords(clouds, spinSystem2, 'HB')
    coords2 = [amide2,alpha2,betas2]
    
    for coord1 in coords1:
      for coord2 in coords2:
        if coord1 and coord2:
          dist.append(getEnsembleCoordsDist(coord1,coord2))
        else:
          dist.append(None)

    for j in range(9):
      if dist[j] is not None:
        key = str(round(dist[j],1))  
        q   = interDistrib[delta][j].get(key, 0.0)
          
        if q:
          q =  max(q, floor)
        else:
          q = floor

 	p2 += log(q)
 	N2 += 1
                    
    if N2 and p2 > -9.9:
      p2 /= float(N2)
      N2 = 1
      
    N += N2
    p += p2
  
  
  if N == 0:
    if v:
      print "barf"
      print r2
      print coords1
      print coords2
      
    
    residue1.likelihood = -10.0
    return -10.0
  
  out = p/float(N)  
  
  if v:
    print "ok", out
    print r2
    print coords1
    print coords2
  
    
  
  residue1.likelihood = out 
  return out


def Log(q):
  
  key = str(q)
  if logDict.get(key) is not None:
    return logDict[key]
  else:
    v = log(q)
    logDict[key] = v
    return v


def getSpinSystemPairLikelihood(spinSystem1, spinSystem2, clouds, interDistrib, v=0, floor=0.0025):

  if not hasattr(spinSystem1, 'pairLikelihood'):
    spinSystem1.pairLikelihood = {}

  val = spinSystem1.pairLikelihood.get(spinSystem2)
  if val is not None:
    return val

  amide1  = getResonanceCoords(clouds, spinSystem1, 'H') or getResonanceCoords(clouds, spinSystem1, 'HD2')
  alpha1  = getResonanceCoords(clouds, spinSystem1, 'HA') or getResonanceCoords(clouds, spinSystem1, 'HA2')
  betas1  = getResonanceCoords(clouds, spinSystem1, 'HB')
  coords1 = [amide1,alpha1,betas1]

  amide2  = getResonanceCoords(clouds, spinSystem2, 'H') or getResonanceCoords(clouds, spinSystem2, 'HD2')
  alpha2  = getResonanceCoords(clouds, spinSystem2, 'HA') or getResonanceCoords(clouds, spinSystem1, 'HA2')
  betas2  = getResonanceCoords(clouds, spinSystem2, 'HB')
  coords2 = [amide2,alpha2,betas2]
   
  dist = []
  for coord1 in coords1:
    for coord2 in coords2:
      if coord1 and coord2:
        dist.append(getEnsembleCoordsDist(coord1,coord2))
      else:
        dist.append(None)

  N = 0.0
  p = 0.0
  for j in range(9):
    if dist[j] is not None:
      key = str(round(dist[j],1))  
      q  = interDistrib[1][j].get(key, 0.0) # +1 only
        
      if q:
        q = max(q, floor)
      else:
        q = floor

      p += log(q)
      N += 1.0
                  
  if N:
    p /= N
  else:
    p = -10.0
    
  spinSystem1.pairLikelihood[spinSystem2] = p
  
  return p


def Log(q):
  
  key = str(q)
  if logDict.get(key) is not None:
    return logDict[key]
  else:
    v = log(q)
    logDict[key] = v
    return v


def generateLikelihood(residues, assignment, clouds, intraDistrib, interDistrib, known, preserveAssign=0):

  p = 0.0
  N = 0
  for residue1 in residues:
      
    spinSystem1 = assignment[residue1]
    if preserveAssign and spinSystem1.residue:
      if spinSystem1.residue is residue1:
        p2 = 0.0
      else:
        p2 = -10.0
    else:
      p2 = getResidueLikelihood(residue1, spinSystem1, assignment, clouds, interDistrib, known)
    
    p += p2
    N += 1
    
  if N == 0:
    return 0.0
 
  return exp(p/len(residues))

def swapSpinSystem2(residues, assignment, n=1):

  #nexus = residues[0].root.nexus


  global globalIndex

  newAssignment = {}

  for residue in residues:
    newAssignment[residue] = assignment[residue]
  
  N = len(residues)
  
  
  # Random swap
  for i in range(n):
    globalIndex += 1
    if globalIndex > (N-1):
      globalIndex = 0
  
    j = globalIndex
    k = j
    while k == j:
      k = randint(0,N-1)

    ss1 = newAssignment[residues[j]]
    ss2 = newAssignment[residues[k]]
    newAssignment[residues[j]] = ss2
    newAssignment[residues[k]] = ss1

  return newAssignment

def swapSpinSystem(assignment, clouds, interDistrib, known, n, iteration):
  
  
  threshold = -5.9
  threshold2 = -4.3
  newAssignment = {}
  residues = assignment.keys()
  for residue in residues:
    newAssignment[residue] = assignment[residue]
  
  N = len(residues)
  
  for i in range(n):
    j = randint(0,N-1)
    spinSystem1 = newAssignment[residues[j]]
    
    h = 0
    while h<25 and getSpinSystemScore(residues[j].ccpCode,spinSystem1) > threshold and getResidueLikelihood(residues[j],spinSystem1, newAssignment, clouds, interDistrib, known) > threshold2:
      j = randint(0,N-1)
      spinSystem1 = newAssignment[residues[j]]
      h += 1

    k = j
 
    h = 0
    while (k == j) and h<25:
      k = randint(0,N-1)
      
      spinSystem2 = newAssignment[residues[k]]
      if k != j:
        if getSpinSystemScore(residues[k].ccpCode,spinSystem1) < threshold:
          if getSpinSystemScore(residues[j].ccpCode,spinSystem2) < threshold:
            k = j
        
      
        if (k != j) and (iteration <500):
          s0  = getResidueLikelihood(residues[j],spinSystem1, newAssignment, clouds, interDistrib, known)
          s0 += getResidueLikelihood(residues[k],spinSystem2, newAssignment, clouds, interDistrib, known)
          s1  = getResidueLikelihood(residues[k],spinSystem1, newAssignment, clouds, interDistrib, known)
          s1 += getResidueLikelihood(residues[j],spinSystem2, newAssignment, clouds, interDistrib, known)
          if s1 <= s0:
            k = j
      
      h += 1
    
    #spinSystem2 = newAssignment[residues[k]]
    newAssignment[residues[j]] = spinSystem2
    newAssignment[residues[k]] = spinSystem1

  return newAssignment

def getInitialAssignment(residues, spinSystems, init, graph=None, progressBar=None):

  print "Gen initial assignment"
   
  shiftList = residues[0].chain.root.findFirstNmrMeasurementList(className='ShiftList')
  #typeScores, bestGuess = getSpinSystemTypes(residues, spinSystems, False, ('1H','13C'), shiftList=shiftList, graph=graph, progressBar=progressBar)
  
  ssDict   = {}
  for ss in spinSystems:
    if ss.residue:
      ssDict[ss.residue] = ss
  
  sss = list(spinSystems)
  assignment = {}
  
  for residue in residues:
    if sss:
      if init.get(residue):
        ss = init[residue]
        assignment[residue] = ss
        if ss in sss:
          sss.remove(ss)

      else:
        assignment[residue] = sss.pop(randint(0, len(sss)-1))

    else:
      assignment[residue] = spinSystems[0]
      
  return assignment


def fitSideChains(scores, assignment, clouds, threshold, shiftList=None):

  # assignment[residue] = spinSystem

  # Get prior form ato chemical shift match
  
  # Search with likelihood from clouds
  
  print 'Calculating distance distributions'
  
  intraDistrib = readDistribution('intraDistribs001.txt')
  
  print 'Monte Carlo...'

  # Get only the good assignments
  assign = {}
  for residue in assignment.keys():
    score = scores.get(residue)
    if score is not None:
      if score > threshold:
        assign[residue] = assignment[residue]

  residues = assign.keys()
  project  = residues[0].root
  
  if shiftList is None:
    shiftList = project.currentNmrProject.findFirstMeasurementList(className='ShiftList')
  
  for residue in residues:
  
    #if residue.ccpCode != 'LYS':
    #  continue
  
    spinSystem = assign[residue]
    shifts = []
    for resonance in spinSystem.resonances:
      resonance.name = None
      shift = resonance.findFirstShift(parentList=shiftList)
      if shift:
        shifts.append(shift)
    
    if shifts:
      distrib = intraDistrib.get(residue.ccpCode) or {}
      fitSideChain(shifts, residue, clouds, distrib)

def fitSideChain(shifts, residue, clouds, distrib, useAssignNames=False):
  
  ccpCode    = residue.ccpCode
  atomScores = {}
  check = {}
  project = residue.root

  from ccpnmr.analysis.core.SpinSystemTyping import getAtomNames, getAtomShiftScore

  atomSets = {}
  shifts2  = []
  for shift in shifts:
    
    """scores[shift] = None
    
    name = None
    if shift.resonance.assignNames:
      name = shift.resonance.assignNames
      
    elif shift.resonance.name:
      name = shift.resonance.name  
      
    name = shift.resonance.name
    
    if name == 'HN':
      name = 'H'"""
 
    if shift.resonance.isotopeCode != '1H':
      continue
 
    shifts2.append(shift)
    done = {}
    atomScores[shift] = {}
    check[shift] = {}

    for atom in residue.atoms:
      if atom.chemAtom.elementSymbol != 'H':
        continue
        
      if not atom.chemAtom.waterExchangeable:
        atomSet = atom.atomSet
        if atomSet and (done.get(atomSet) is None):
          boundAtom = None
          for atom2 in list(atom.chemAtom.chemBonds)[0].chemAtoms:
            if atom2 is not atom.chemAtom:
              boundAtom = atom2
        
          boundShift = None
          if shift.resonance.covalentlyBound:
            boundResonance = list(shift.resonance.covalentlyBound)[0]
            boundShift = boundResonance.findFirstShift(parentList=shift.parentList)
        
          atomSetName       = atomSet.name
          atomName          = atom.name
          done[atomSet]     = True
          atomSets[atomSet] = True
 
          if atom.chemAtom.elementSymbol == shift.resonance.isotope.chemElement.symbol:
            score = getAtomShiftScore(project, atomSetName, ccpCode, shift)
            if not score:
              score = getAtomShiftScore(project, atomName, ccpCode, shift)
            
            if score and boundAtom and boundShift:
              score2 = getAtomShiftScore(project, boundAtom.name, ccpCode, boundShift)
              score  = score2#sqrt(((score*score) + 4*(score2*score2)) /5.0)
 
              check[shift][atomSetName] = boundAtom.name, boundShift.value, score2
              
            if score:
                
              atomScores[shift][atomSetName] = log(score)
            else:
              atomScores[shift][atomSetName] = -8.0
              check[shift][atomSetName] = None
 
          else:
            atomScores[shift][atomSetName] = None
            check[shift][atomSetName] = None
  
  shifts = shifts2
  atomNames0 = [(ass.findFirstAtom().name, ass.name) for ass in atomSets.keys()] # getAtomNames(residue)
    
  best    = -20000.0
  bestM   = None
  passes  = 0
  J       = 0  
  i       = 0
  
  nameMap = {}
  for atomName, atomSetName in atomNames0:
    nameMap[atomName] = atomSetName
  
  print "Fitting residue %d %s" % (residue.seqCode, residue.ccpCode)
  #print atomNames0
  
  sz = len(atomNames0)
  while passes < (sz*sz*2):
  
    i += 1
    if bestM is None:
      bestM = {}
      shifts1 = list(shifts)
      
      for atomName, atomSetName in atomNames0:
        if shifts1:
          j = randint(0,len(shifts1)-1)
          bestM[atomName] = shifts1.pop(j)
        else:
          bestM[atomName] = None
      
      mapping = bestM.copy()
      
    else:
      mapping = bestM.copy()
      
      if len(atomNames0) < 2:
        break
      
      j = J
      J += 1
      if J >= len(atomNames0):
        J = 0
      
      j = randint(0, len(atomNames0)-1)  
      
      k = j
      while k == j:
        k = randint(0, len(atomNames0)-1)    

      swap = mapping[atomNames0[j][0]]
      mapping[atomNames0[j][0]] = mapping[atomNames0[k][0]]
      mapping[atomNames0[k][0]] = swap
       
    prior = 0.0
    for atomName, atomSetName in atomNames0:
      shift = mapping[atomName]
      if shift:
        s = atomScores[shift][atomSetName]
        if s:
          prior += s
        else:
          prior += -50.0
    
    
    likelihood = getAtomCloudScores(mapping, ccpCode, clouds, distrib)
          
    score = likelihood + prior 
    #score = prior 
    #score = likelihood

    if score > best:
      #print i, score, likelihood, prior
      #print i, prior
      passes = 0
      best   = score
      bestM  = mapping
      
      #print i, '%.3f' % likelihood,'%.3f' % prior,'%.3f' % score, 
      #print ' '.join(['%s:%.3f' % (key,mapping[key].value) for key in mapping.keys() if mapping.get(key)])
    else:
      passes += 1
 
 
  if bestM:
    for atomName in sideChainAtoms:
     #for atomName in bestM.keys():
     
     shift = bestM.get(atomName)
     if shift:
       s = getAtomCloudScores2(bestM, ccpCode, clouds, distrib, atomName) or -10.0
       print '%-5.5s S: %.3f D: %.3f' % (atomName, atomScores[shift][nameMap[atomName]], s)

  # do some provisional assignment
  if bestM:
    for atomName in sideChainAtoms:
      #for atomName in bestM.keys():
      
      shift = bestM.get(atomName)
      if shift:
        resonance = shift.resonance
        resonance.name = None
      
        cloudScore = getAtomCloudScores2(bestM, ccpCode, clouds, distrib, atomName) or -10.0
        typeScore  = atomScores[shift][nameMap[atomName]]
        #resonance.setName(atomName)
        #print atomName, 
        #continue
        print ' ',
        if atomName in ('HA','HA1','HA2'):
          resonance.setName(atomName)
          print atomName, 
        
        else:
          if (cloudScore + typeScore) > -7.0:
            resonance.setName(atomName)
            print atomName, 
          else:
            print ''
            break
              
        print ''
         
        """
        atomNames4 = []
        if resonance.resonanceSet:
          for atomSet in resonance.resonanceSet.atomSets:
            for atom4 in atomSet.atoms:
              atomNames4.append(atom4.name) 
        if atomNames4 and (atomName not in atomNames4):
          #if residue.seqCode == 8:
          #print i, best
          print "******ERR*******", atomName, shift.value
          print 'Debug %5.5s %3.3f - %s' % (atomName, shift.value, makeResonanceGuiName(resonance))
          for atom4 in residue.atoms:
            if atom4.chemAtom.elementSymbol == 'H':
              if not atom4.chemAtom.waterExchangeable:
                print '     %3.3f %5.5s %.3f' % (shift.value, atom4.atomSet.name, atomScores[shift][atom4.atomSet.name] or -1000.0),
                print check[shift].get(atom4.atomSet.name)
        """  
  return bestM

def getAtomCloudScores(mapping, ccpCode, clouds, distrib, atomNames=None):

  score = 0.0
  if not atomNames:
    atomNames = mapping.keys()
  
  
  count = 0.0
  
  for atomName1 in mapping.keys():
    shift1  = mapping[atomName1]
    if shift1:
      ensemble1 = []
      for cloud in clouds:
        coords = []
        coord = getResonanceCoord(cloud, shift1.resonance)
        if coord:
          coords.append(coord)

        if coords:
          coord = averageCoord(coords)
          ensemble1.append(coord)

      for atomName2 in atomNames:
        if atomName1 != atomName2:
          shift2 = mapping[atomName2]
          if shift2:
            ensemble2 = []
            for cloud in clouds:
              coords = []
              coord = getResonanceCoord(cloud, shift2.resonance)
              if coord:
                coords.append(coord)

              if coords:
                coord = averageCoord(coords)
                ensemble2.append(coord)
 
            if ensemble1 and ensemble2:
              dist = getEnsembleCoordsDist(ensemble1,ensemble2)
              d2   = dist
              key  = str(round(dist,1))
              val  = distrib[atomName1][atomName2].get(key)
              while (not val) and d2 > 0.0:
                d2 -= 0.1
                key = str(round(d2,1))
                val = distrib[atomName1][atomName2].get(key)
                
              if val:
                score += log(val)
                count += 1.0
              """else:
                print '  %3.3s' % atomName1, '%3.3s' % atomName2,
                print makeResonanceGuiName(shift1.resonance),
                print makeResonanceGuiName(shift2.resonance),
                print '%.3f' % log(val or 1.0), '%.3f' % dist
                kk = distrib[atomName1][atomName2].keys()
                kk.sort
                print kk
              
                #score += -10.0"""

  if count:
    score /= count
  else:
    score = -10.0

  #print ' '.join(['%s:%.3f' % (key,mapping[key].value) for key in mapping.keys()]), '%.3f' % score
  
  return score

def getAtomCloudScores2(mapping, ccpCode, clouds, distrib, atomName, atomNames=None):

  score = 0.0
  if not atomNames:
    atomNames = mapping.keys()
  
  
  count = 0.0
  
  for atomName1 in [atomName,]:
    shift1  = mapping[atomName1]
    if shift1:
      ensemble1 = []
      for cloud in clouds:
        coords = []
        coord = getResonanceCoord(cloud, shift1.resonance)
        if coord:
          coords.append(coord)

        if coords:
          coord = averageCoord(coords)
          ensemble1.append(coord)

      for atomName2 in atomNames:
        if atomName1 != atomName2:
          shift2 = mapping[atomName2]
          if shift2:
            ensemble2 = []
            for cloud in clouds:
              coords = []
              coord = getResonanceCoord(cloud, shift2.resonance)
              if coord:
                coords.append(coord)

              if coords:
                coord = averageCoord(coords)
                ensemble2.append(coord)
 
            if ensemble1 and ensemble2:
              dist = getEnsembleCoordsDist(ensemble1,ensemble2)
              d2   = dist
              key  = str(round(dist,1))
              val  = distrib[atomName1][atomName2].get(key)
              while (not val) and d2 > 0.0:
                d2 -= 0.1
                key = str(round(d2,1))
                val = distrib[atomName1][atomName2].get(key)
                
              if val:
                score += log(val)
                count += 1.0
              """else:
                print '  %3.3s' % atomName1, '%3.3s' % atomName2,
                print makeResonanceGuiName(shift1.resonance),
                print makeResonanceGuiName(shift2.resonance),
                print '%.3f' % log(val or 1.0), '%.3f' % dist
                kk = distrib[atomName1][atomName2].keys()
                kk.sort
                print kk
              
                #score += -10.0"""

  if count:
    score /= count

  #print ' '.join(['%s:%.3f' % (key,mapping[key].value) for key in mapping.keys()]), '%.3f' % score
  
  return score

def getDistDistributions(dataDirName='/home/tjs234/ccpn/recoord/', intraFileName='intraDistribs001.txt', interFileName='interDistribs001.txt'):

  if not( exists(intraFileName) and isfile(intraFileName) and exists(interFileName) and isfile(interFileName)):
    calcDistDistributions(dataDirName, intraFileName, interFileName)

  inter = readDistribution(interFileName)
  intra = readDistribution(intraFileName)

  return intra, inter

def calcDistDistributions(dirName, intraFileName='intraDistribs001.txt', interFileName='interDistribs001.txt'):

  intra = {}
  inter = {}
  
  files = listdir(dirName)

  for fileName in files:
    path = '%s%s/%s.xml' % (dirName, fileName, fileName)
    path2 = '%s%s/' % (dirName, fileName)
    if not exists(path):
      continue
  
    dict =  makeStructureDictFromXml(path2, path)
    #model = dict.keys()[0]

    for chain in dict.keys():
      chainDict = dict[chain]
      for residue1 in chainDict.keys():
        ccpCode = chainDict[residue1]['resName']
        atoms   = chainDict[residue1]['atoms'].keys()
        if intra.get(ccpCode) is None:
          intra[ccpCode] = {}
          for atom1 in atoms:
            intra[ccpCode][atom1] = {}
            for atom2 in atoms:
              intra[ccpCode][atom1][atom2] = {}
          
        amide1 = getAmideCoord(chainDict[residue1]['atoms'])
        alpha1 = getAlphaCoord(chainDict[residue1]['atoms'])
        betas1 = getBetasCoord(chainDict[residue1]['atoms'])
        coords1 = [amide1,alpha1,betas1]

        if not (amide1 and alpha1):
          continue

        for residue2 in chainDict.keys():
          delta = residue2 - residue1
          if delta == 0:
            continue
          
          amide2 = getAmideCoord(chainDict[residue2]['atoms'])
          alpha2 = getAlphaCoord(chainDict[residue2]['atoms'])
          betas2 = getBetasCoord(chainDict[residue2]['atoms'])
          coords2 = [amide2,alpha2,betas2]

          if not (amide1 and alpha2):
            continue
            if coord1 and coord2:
              dist = getCoordsDist(coord1, coord2)
              key  = str(round(dist, 1))

          if inter.get(delta) is None:
            inter[delta] = [{} for x in range(9)]
          
          index = 0
          for coord1 in coords1:
            for coord2 in coords2:
              if coord1 and coord2:
                dist = getCoordsDist(coord1, coord2)
                key  = str(round(dist, 1))
                inter[delta][index][key] = inter[delta][index].get(key, 0) + 1
              index += 1   
          
        for atom1 in atoms:
          coord1 = chainDict[residue1]['atoms'][atom1]
          for atom2 in atoms:
            if atom1 == atom2:
              continue
            
            coord2 = chainDict[residue1]['atoms'][atom2]
            dist = getCoordsDist(coord1, coord2)
            key = str(round(dist, 1))
            if not intra[ccpCode].get(atom1):
              intra[ccpCode][atom1] = {}
              
            if not intra[ccpCode][atom1].get(atom2):
              intra[ccpCode][atom1][atom2] = {}
              
            intra[ccpCode][atom1][atom2][key] = intra[ccpCode][atom1][atom2].get(key, 0) + 1

  # normalise inter residue dists
  for delta in inter.keys():
    for index in range(9):
      distr = inter[delta][index]
      N = 0.0
      for key in distr.keys():
        N += float(distr[key])
      
      for key in distr.keys():
        distr[key] /= N
        
      inter[delta][index] = distr
      

  # normalise intra residue dists
  for ccpCode in intra.keys():
    for atom1 in intra[ccpCode].keys():
      for atom2 in intra[ccpCode][atom1].keys():
        distr = intra[ccpCode][atom1][atom2]
        N = 0.0
        for key in distr.keys():
          N += float(distr[key])
        
        for key in distr.keys():
          distr[key] /= N
          
        intra[ccpCode][atom1][atom2] = distr

  intraFile = open(intraFileName, 'w')
  cPickle.dump(intra, intraFile)
  interFile = open(interFileName, 'w')
  cPickle.dump(inter, interFile)
  
  intraFile.close()
  interFile.close()
  
def getCoordsDist(coord1, coord2):

  dx = coord1[0] - coord2[0]
  dy = coord1[1] - coord2[1]
  dz = coord1[2] - coord2[2]
  dist = sqrt((dx*dx) + (dy*dy) + (dz*dz))

  return dist

def readDistribution(fileName):

  
  file = open(fileName, 'r')
  dict = cPickle.load(file)

  return dict

def getAlphaCoord(dict):

  coord = None

  if dict.get('HA'):
    coord = dict['HA']
  else:
    coords = []
    for atom in dict.keys():
      if atom[:2] == 'HA':
        coords.append(dict[atom])

    if coords:
      coord = averageCoord(coords)
  
  return coord


def getAmideCoord(dict):

  coord = None

  if dict.get('H'):
    coord = dict['H']
  elif dict.get('HN'):
    coord = dict['HN']
  else:
    coords = []
    for atom in dict.keys():
      if atom[:2] == 'HD':
        coords.append(dict[atom])

    if coords:
      coord = averageCoord(coords)
  
  return coord


def getBetasCoord(dict):    
  
  coord = None

  if dict.get('HB'):
    coord = dict['HB']
  else:
    coords = []
    for atom in dict.keys():
      if atom[:2] == 'HB':
        coords.append(dict[atom])

    if coords:
      coord = averageCoord(coords)

  return coord


def averageCoord(coordList):

  sum = [0.0, 0.0, 0.0]
  for coord in coordList:
    for i in range(3):
      sum[i] += coord[i]

  N = float(len(coordList))
  for i in range(3):
    sum[i] /= N

  return sum


def getEnsembleCoordsDist(ensemble1, ensemble2):

  n    = len(ensemble1) 
  dist = 0.0

  for j in range(n):
    coords1 = ensemble1[j]
    coords2 = ensemble2[j]
    sum = 0.0
    for i in range(3):
      diff = coords2[i]-coords1[i]
      sum += diff * diff
    
    dist += sqrt(sum)
  dist /= float(n)
    
  return dist


def getResonanceCoords(clouds, spinSystem, atomType):

  # find spin system in cloud
  # find atom type within spin system
  # get coords for each cloud of an FOC

  if not hasattr(spinSystem, 'resonanceCoords'):
    spinSystem.resonanceCoords = {}
    
  if spinSystem.resonanceCoords.get(atomType):
    return spinSystem.resonanceCoords[atomType]

  resonances = []
  for resonance in spinSystem.resonances:
    if resonance.name == atomType:
      resonances.append(resonance)
    elif atomType in resonance.assignNames:
      resonances.append(resonance)
    elif resonance.name and resonance.name[:2] == atomType:
      resonances.append(resonance)

  ensemble = []
  for cloud in clouds:
    coords = []
    for resonance in resonances:
      coord = getResonanceCoord(cloud, resonance)
      if coord:
        coords.append(coord)

    if coords:
      coord = averageCoord(coords)
      ensemble.append(coord)

  spinSystem.resonanceCoords[atomType] = ensemble

  return ensemble
  
def getResonanceCoord(cloud, resonance):

  # not all spinSystem resonances are represented in a cloud
  return cloud.get(resonance)

def fastCloudThreader(chain, clouds, spinSystems=None, nSteps=5000, graph=None, preserveAssign=False, progressBar=None):

  print "Loading distance distributions"
  intra, inter = getDistDistributions()
  
  project   = chain.root
  shiftList = project.currentNmrProject.findFirstMeasurementList(className='ShiftList')
  residues  = list(chain.residues)

  if not spinSystems:
    spinSystems = list(project.currentNmrProject.resonanceGroups)

  for residue in residues:
    residue.spinSystemScoreDict = {}

  print "Initialisation"
  residueDict = {}
  for ss in spinSystems:
    ss.shifts = getSpinSystemShifts(ss, shiftList, ['1H','13C','15N'] )
    ss.resonanceCoords = {}
    ss.codeScoreDict   = {} # Clear cache
    ss.pairLikelihood  = {}
    if ss.residue:
      residueDict[ss.residue] = ss

  for i in range(len(residues)-1):
    r1  = residues[i]
    r2  = residues[i+1]
    ss1 = residueDict.get(r1)
    ss2 = residueDict.get(r2)
    
    if ss1 and ss2:
      s1  = getSpinSystemPairLikelihood(ss1, ss2, clouds, inter)
      s2  = getSpinSystemScore(r1.ccpCode, ss1)
      dat = (r1.seqCode,r1.ccpCode,r2.seqCode,r2.ccpCode,s1,s2)
      print '%d%s %d%s %.3f %.3f' % dat


  print "Peptide search"
  dict, hits   = peptideSearch(8, residues, spinSystems, clouds, inter, exclude=None)

  dict2, hits2 = peptideSearch(7, residues, spinSystems, clouds, inter, exclude=dict)


  for obj in dict2.keys(): # From potential second round - shorter windows
    dict[obj] = dict2[obj]
 
  for residue in hits2.keys(): # From potential second round - shorter windows
    hits[residue] = hits2[residue]

  dict2, hits2 = peptideSearch(6, residues, spinSystems, clouds, inter, exclude=dict)

  for obj in dict2.keys(): # From potential second round - shorter windows
    dict[obj] = dict2[obj]
 
  for residue in hits2.keys(): # From potential second round - shorter windows
    hits[residue] = hits2[residue]

  dict2, hits2 = peptideSearch(5, residues, spinSystems, clouds, inter, exclude=dict)

  for obj in dict2.keys(): # From potential second round - shorter windows
    dict[obj] = dict2[obj]
 
  for residue in hits2.keys(): # From potential second round - shorter windows
    hits[residue] = hits2[residue]


  print "Finding unique mappings"
  assignment = {}
  matchDict = {}
  nhits = 0
  for residue in residues:
    if not dict.get(residue):
      continue

    if hits.get(residue):
    
      match = None
      ss = dict[residue][0]
      if len(dict[ss]) < 2:
        if hits[residue][ss] > 0:
          match = ss
      
      else:
        for ss in hits[residue].keys():
          if hits[residue][ss] > 2:
            match = ss
            break
       
      if match:   
        assignment[residue] = match
        matchDict[match] = residue
        gb = "?????"
        if match.residue:
          if match.residue is residue:
            gb = 'Good!'
            nhits += 1
          else:
            gb = '*BAD*'
            
        print "Match %d%s - %d [%s]" % (residue.seqCode,residue.ccpCode,match.serial,gb)

  print "HITS: ", nhits

  #return

  residues = [] # Can also eliminate the good ones!
  for residue in chain.residues:
    if residue.ccpCode != 'PRO':
      ss = assignment.get(residue)
      if ss is None:
        residues.append(residue)

  spinSystems2 = []
  for spinSystem in spinSystems:
    if matchDict.get(spinSystem) is None:
      spinSystems2.append(spinSystem)

  nSteps0 = len(residues)
  nSteps0 *= nSteps0 * 5

  print "Searching posterior 1"
  score1, scores, assignment3 = searchPosterior(project, residues, spinSystems2,
                                               clouds, nSteps=nSteps0, graph=graph, known=assignment,
                                               preserveAssign=preserveAssign, inter=inter,
                                               progressBar=progressBar )

  for residue in assignment3.keys():
    assignment[residue] = assignment3[residue]
  
  score2, scores, assignment4 = searchPosterior(project, chain.residues, spinSystems,
                                               clouds, nSteps=100, graph=graph, known=assignment, 
                                               preserveAssign=preserveAssign, inter=inter,
                                               progressBar=progressBar )
  #return score, scores, assignment4
  
  spinSystems2 = []
  residues = []
  for residue in assignment4.keys():
    score0 = scores[residue]
    if score0 < -4.9:
      ss = assignment4.get(residue)
      if ss:
        spinSystems2.append(ss)
      residues.append(residue)

  if residues and spinSystems2:
    nSteps0 = len(residues)
    nSteps0 *= nSteps0 * 5
 
    print "Searching posterior 2"
    score3, scores, assignment4 = searchPosterior(project, residues, spinSystems2,
                                                 clouds, nSteps=nSteps0, graph=graph, known=assignment,
                                                 preserveAssign=preserveAssign, inter=inter,
                                                 progressBar=progressBar )
  for residue in assignment4.keys():
    assignment[residue] = assignment4[residue]

     
  print "Searching posterior 3"
  score, scores, assignment5 = searchPosterior(project, chain.residues, spinSystems,
                                               clouds, nSteps=nSteps, graph=graph, known=assignment, 
                                               preserveAssign=preserveAssign, inter=inter,
                                               progressBar=progressBar )


  return score, scores, assignment5

def followRoute(peptide, spinSystems, ensemble, inter, route=None,
               typeThreshold=-5.8, scoreThreshold=-5.0, totalScore=0.0):
  
  if route is None:
    route = []
  
  index = len(route)
  
  #print index, totalScore, [ss.serial for ss in route]
  
  if index == len(peptide):
    return totalScore, route
  
  bestScore = None
  bestRoute = []
  for spinSystem in spinSystems:
    if spinSystem in route:
      continue
  
  
    typeScore = getSpinSystemScore(peptide[index].ccpCode, spinSystem)
    if typeScore > typeThreshold:
      #print spinSystem.serial, 
      
      score = 0.0
      if index > 0:
        spinSystem0 = route[-1]
        score = getSpinSystemPairLikelihood(spinSystem0, spinSystem, ensemble, inter)
        if score < scoreThreshold:
          continue
        
      route2 = list(route)
      route2.append(spinSystem)
      totalScore += typeScore + score
      
      partScore, partRoute = followRoute(peptide, spinSystems, ensemble, inter, route2,
                                         typeThreshold, scoreThreshold, totalScore )
  
      if len(partRoute) < len(peptide):
        continue
  
      if (bestScore is None) or (partScore > bestScore):
        bestScore = partScore
        bestRoute = partRoute
  
  """if (index == 0) and (bestScore is None):
    if len(peptide) > 5:
      peptide = peptide[:-1]
      bestScore, bestRoute = followRoute(peptide, spinSystems, ensemble, inter, None,
                                         typeThreshold, scoreThreshold, totalScore )"""
  
  return bestScore, bestRoute  
      
def peptideSearch(win, residues, spinSystems, clouds, inter, exclude=None):
  if not exclude:
    exclude = {}

  i = 0
  hits = {}
  partners = {}
  while i+win < len(residues):

    if exclude.get(residues[i]):
      i += 1
      continue

    peptide  = residues[i:i+win]
    """ccpCodes = [r.ccpCode for r in peptide]
    if 'PRO' in ccpCodes:
      i += 1
      continue"""

    score, route = followRoute(peptide, spinSystems, clouds, inter)

    pp = ' '.join(['%3d%s' % (r.seqCode,r.ccpCode) for r in peptide])
    if route:
      rr = [ss.residue for ss in route]
      dd = ' '.join(['%d%s' % (r.seqCode,r.ccpCode) for r in rr if r])
      print "Best for", pp, '\n        ', dd
      #print "Match for", pp
      for j in range(len(route)):
        residue    = peptide[j]
        spinSystem = route[j]

        if exclude.get(residue) is not None:
          continue

        if hits.get(residue) is None:
          hits[residue] = {}

        hits[residue][spinSystem] = hits[residue].get(spinSystem, 0) + 1

        if partners.get(residue) is None:
          partners[residue] = []

        if partners.get(spinSystem) is None:
          partners[spinSystem] = []

        if residue not in partners[spinSystem]:
          partners[spinSystem].append(residue)

        if spinSystem not in partners[residue]:
          partners[residue].append(spinSystem)
    else:
      #print '%d Missing ' % i, pp
      i += 1
      continue

    i += 1

  return partners, hits

if __name__ == '__main__':
  calcDistDistributions('/home/tjs234/ccpn/recoord/')
  
  

