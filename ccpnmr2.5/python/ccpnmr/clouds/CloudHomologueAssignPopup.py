"""
======================COPYRIGHT/LICENSE START==========================

CloudHomologueAssignPopup.py: Part of the CcpNmr Clouds program

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

from os import listdir

from memops.general import Implementation

from memops.editor.BasePopup    import BasePopup

from memops.gui.FloatEntry      import FloatEntry
from memops.gui.IntEntry        import IntEntry
from memops.gui.Entry           import Entry
from memops.gui.Label           import Label
from memops.gui.PulldownMenu    import PulldownMenu
from memops.gui.Util            import createDismissHelpButtonList
from memops.gui.ScrolledGraph   import ScrolledGraph
from memops.gui.ProgressBar     import ProgressBar

from ccpnmr.clouds.CloudHomologueAssign import cloudHomologueAssign
from ccpnmr.clouds.CloudBasic           import getFileNamesFromPattern, getCloudsFromFile

from ccpnmr.analysis.core.AssignmentBasic import getResonanceAtomTuple, assignSpinSystemResidue

from memops.gui.MessageReporter import showWarning

class CloudHomologueAssignPopup(BasePopup):

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

    label = Label(guiFrame, text='Clouds files: ')
    label.grid(row=row, column=2, sticky=Tkinter.NW)
    self.filenameEntry = Entry(guiFrame,text='perfect00.pdb')
    self.filenameEntry.grid(row=row, column=3, sticky=Tkinter.NW)


    row += 1
    label = Label(guiFrame, text='Chain: ')
    label.grid(row=row, column=0, sticky=Tkinter.NW)
    self.chainPulldown = PulldownMenu(guiFrame, self.changeChain, selected_index=-1, do_initial_callback=0)
    self.chainPulldown.grid(row=row, column=1, sticky=Tkinter.NW)

    label = Label(guiFrame, text='Thread steps: ')
    label.grid(row=row, column=2, sticky=Tkinter.NW)
    self.numStepsEntry = IntEntry(guiFrame,text=3000)
    self.numStepsEntry.grid(row=row, column=3, sticky=Tkinter.NW)
    row += 1

    label = Label(guiFrame, text='Homologue PDB file: ')
    label.grid(row=row, column=0, sticky=Tkinter.NW)
    self.pdbEntry = Entry(guiFrame,text='')
    self.pdbEntry.grid(row=row, column=1, sticky=Tkinter.NW)

    label = Label(guiFrame, text='Dist. Threshold: ')
    label.grid(row=row, column=2, sticky=Tkinter.NW)
    self.distEntry = FloatEntry(guiFrame,text=3.0)
    self.distEntry.grid(row=row, column=3, sticky=Tkinter.NW)

    row += 1

    label = Label(guiFrame, text='Global score: ')
    label.grid(row=row, column=0, sticky=Tkinter.NW)
    self.globalScoreLabel = Label(guiFrame, text='')
    self.globalScoreLabel.grid(row=row, column=1, sticky=Tkinter.NW)

    label = Label(guiFrame, text='Assignment Threshold: ')
    label.grid(row=row, column=2, sticky=Tkinter.NW)
    self.thresholdEntry = FloatEntry(guiFrame,text=-4.5)
    self.thresholdEntry.grid(row=row, column=3, sticky=Tkinter.NW)

    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    self.graph = ScrolledGraph(guiFrame, width=300, height=200)
    self.graph.grid(row=row, column=0, columnspan=4, sticky = Tkinter.NSEW)

    row += 1
    texts    = ['Run','Assign!']
    commands = [self.run, self.assignSpinSystems]
    bottomButtons = createDismissHelpButtonList(guiFrame,texts=texts,commands=commands,expands=0,help_url=None)
    bottomButtons.grid(row=row, column=0, columnspan=4, sticky=Tkinter.EW)
    self.assignButton = bottomButtons.buttons[1]

    for func in ('__init__','delete'):
      Implementation.registerNotify(self.updateMolSystems, 'ccp.molecule.MolSystem.MolSystem', func)
      Implementation.registerNotify(self.updateChains, 'ccp.molecule.MolSystem.Chain', func)
    
    self.updateMolSystems()
    self.updateChains()

  def update(self):
  
    if self.assignment and self.scores:
      self.assignButton.enable()
    else:
      self.assignButton.disable()  

  def run(self):
  
    if self.chain:
      pattern = self.filenameEntry.get()
      nSteps  = self.numStepsEntry.get() or 4000
      pdbFile = self.pdbEntry.get()
      dist    =  self.distEntry.get() or 3.0
      pgb     = ProgressBar(self, text='Searching', total=nSteps)
      files   = getFileNamesFromPattern(pattern , '.')
      if not files:
        return
      clouds  = getCloudsFromFile(files, self.chain.root)
       
      score, self.scores, self.assignment = cloudHomologueAssign(self.chain, clouds, pdbFile, dist, nSteps, self.graph, pgb)
 
      pgb.destroy()
      self.globalScoreLabel.set(str(score))
      self.update()

  def assignSpinSystems(self):
   
    if self.assignment and self.scores:
      if showWarning('Query','Are you sure?'):
        threshold = self.thresholdEntry.get() or -4.0
        i = 0
 
        for residue in self.assignment.keys():
          if self.scores[residue] > threshold:
            spinSystem = self.assignment[residue]
            assignSpinSystemResidue(spinSystem,residue=None)
 
        for residue in self.assignment.keys():
          if self.scores[residue] > threshold:
            i += 1
            spinSystem = self.assignment[residue]
            assignSpinSystemResidue(spinSystem,residue=residue)
      
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
