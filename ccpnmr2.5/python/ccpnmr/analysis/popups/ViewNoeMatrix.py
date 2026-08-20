"""
======================COPYRIGHT/LICENSE START==========================

ViewNoeMatrix.py: Part of the CcpNmr Analysis program

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

import Tkinter


from memops.general import Implementation

from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.Label import Label
from memops.gui.MessageReporter import showWarning
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes
from ccpnmr.analysis.frames.NoeMatrix import NoeMatrix

# TBD More colour options?

class ViewNoeMatrix(BasePopup):
  """
  **Display a Density Matrix of Residue-Residue Contact Information**

  This popup window is based around the "Interaction Matrix" in the first tab
  and is designed to give a graphical representation of the residue to residue
  interactions that are present in NMR data. At present this informations comes
  from peak lists; to show assignment connectivities (e.g. NOE) and restraint
  lists; to show distance and H-bond restraints.

  To operate this system the user opens the "Peak Lists" and "Restraint Lists"
  tabs and toggles the "Use?" columns to enable or disable peak and restraint
  lists in the analysis. After pressing the [Draw] button at the top, the
  interaction matrix in the first tab is redrawn to display the connectivity
  information form the selected sources.

  The main interaction matrix is colour coded so that darker squares represent 
  a greater number of observed interactions. The coordinates of each square
  represents the intersection between two residues from the axes. The residues
  that are used on the X- and Y-axes many be specified via the "Residues" tab,
  both in terms of sequence range and molecular chain. To determine the precise
  residues that correspond to a given interaction the user can hover the mouse
  cursor over a small square in the matrix and the identities of the residues
  are displayed on the two axes.

  **Tips**

  If the interaction matrix display is too large or too small the size of
  the chart can be adjusted via the [+] and [-] buttons.

  A PostScript file of the density matrix can be saved via the right mouse menu
  of the matrix tab.
  """

  def __init__(self, parent, *args, **kw):
 
    self.guiParent = parent
    self.peakLists = []
    self.constraintLists = []
    self.constraintSet = None
    BasePopup.__init__(self, parent=parent, title="Chart : Residue Interaction Matrix", **kw)
 
  def body(self, guiFrame):

    
    self.geometry('600x700')

    self.noeMatrix = None
    
    self.xMol = None
    self.yMol = None

    guiFrame.expandGrid(0,0)
    
    tipTexts = ['The main colour density matrix that displays residue-residue interaction strength',
                'Selects which peak lists to show residue-residue interactions for',
                'Selects which distance restraint lists to show residue-residue interactions for',
                'Specifies which molecular chains and residue ranges to consider on the matrix axes']
    options = ['Interaction Matrix','Peak Lists','Restraint Lists','Residues']
    tabbedFrame = TabbedFrame(guiFrame, options=options,
                              grid=(0,0), tipTexts=tipTexts)
    frameA, frameB, frameC, frameD = tabbedFrame.frames
    

    
    #
    # Matrix
    #
    
    frameA.expandGrid(0,0)
    
    self.noeMatrix = NoeMatrix(frameA, borderwidth=1, relief='flat',
                               background='darkGrey', labelAxes=False)
    self.noeMatrix.grid(row=0, column=0, sticky='nsew')
    self.noeMatrix.updateAfter()
    
    #
    # Peak Lists
    #
    
    frameB.expandGrid(0,0)
    
    tipTexts = ['The experiment:spectrum name of the peak list that may be considered',
                'The serial number of the peak lists within its spectrum',
                'Sets whether or not the peak list will be used as a source of residue interaction information']
    headingList = ['Spectrum','PeakList','Use?']
    editWidgets      = [None, None, None]
    editGetCallbacks = [None, None, self.togglePeakList]
    editSetCallbacks = [None, None, None]
    self.peakMatrix = ScrolledMatrix(frameB,headingList=headingList,
                                     editSetCallbacks=editSetCallbacks,
                                     editGetCallbacks=editGetCallbacks,
                                     editWidgets=editWidgets,
                                     multiSelect=False, grid=(0,0),
                                     tipTexts=tipTexts)

    #
    # Restraints
    #

    frameC.expandGrid(1,1)
    
    label = Label(frameC, text='Restraint Set: ', grid=(0,0))
    tipText = 'Selects which set of restraints to select restraint connectivity information from'
    self.constraintSetPulldown = PulldownList(frameC, self.changeConstraintSet,
                                              grid=(0,1), tipText=tipText)
    
    tipTexts = ['The serial number, within the restraint set, and name of the restraint list',
                'Whether the restraint list is a Distance restraint list or an H-Bond restraint list',
                'Sets whether or not the restraint list will be used as a source of residue interaction information']
    headingList = ['List','Type','Use?']
    editWidgets      = [None, None, None]
    editGetCallbacks = [None, None, self.toggleConstraintList]
    editSetCallbacks = [None, None, None]
    self.constraintMatrix = ScrolledMatrix(frameC,headingList=headingList,
                                           editSetCallbacks=editSetCallbacks,
                                           editGetCallbacks=editGetCallbacks,
                                           editWidgets=editWidgets,
                                           multiSelect=False, grid=(1,0),
                                           tipTexts=tipTexts, gridSpan=(1,2))
 
    #
    # Residues
    #

    frameD.expandGrid(2,5)

    label = Label(frameD,text='X axis: ', grid=(0,0))
    label = Label(frameD,text='Y axis: ', grid=(1,0))
    
    tipText = 'Selects which molecular chain to use as the sequence along the horizontal axis of the interaction matrix'
    self.xMolPulldown = PulldownList(frameD, callback=self.changeMolX,
                                     grid=(0,1), tipText=tipText)
    tipText = 'Selects which molecular chain to use as the sequence along the vertical axis of the interaction matrix'
    self.yMolPulldown = PulldownList(frameD, callback=self.changeMolY, 
                                     grid=(1,1), tipText=tipText)
    
    tipText = 'Sets the number of the first residue that appears at the start (left) of the horizontal matrix axis'
    self.xEntryStart = IntEntry(frameD,text='',returnCallback=self.setMolRanges, 
                                width=5, grid=(0,2), tipText=tipText)
    tipText = 'Sets the number of the first residue that appears at the start (bottom) of the vertical matrix axis'
    self.yEntryStart = IntEntry(frameD,text='',returnCallback=self.setMolRanges, 
                                width=5, grid=(1,2), tipText=tipText)
    
    label = Label(frameD,text=' to ', grid=(0,3))
    label = Label(frameD,text=' to ', grid=(1,3))
    
    tipText = 'Sets the number of the last residue that appears at the end (right) of the horizontal matrix axis'
    self.xEntryStop  = IntEntry(frameD,text='',returnCallback=self.setMolRanges, 
                                width=5, grid=(0,4), tipText=tipText)
    tipText = 'Sets the number of the last residue that appears at the end (top) of the vertical matrix axis'
    self.yEntryStop  = IntEntry(frameD,text='',returnCallback=self.setMolRanges,  
                                width=5, grid=(1,4), tipText=tipText)

    #
    # Main
    #

    tipTexts = ['Using the selected peak lists, restraint lists and residue ranges draw a colour density residue interaction matrix',
                'Zoom in on the matrix view; make the boxes larger',
                'Zoom out on the matrix view; make the boxes smaller']
    commands = [self.updateNoes,self.zoomIn,self.zoomOut]
    texts    = ['Draw','[+]','[-]']
    bottomButtons = UtilityButtonList(tabbedFrame.sideFrame, commands=commands,
                                      texts=texts, helpUrl=self.help_url,
                                      grid=(0,0), sticky = 'e', tipTexts=tipTexts)

    self.updateConstraintSets()
    self.updateMols()
    self.update()
    #self.updateNoes()
 
    self.administerNotifiers(self.registerNotify)
    
  def administerNotifiers(self, notifyFunc):
 
    for func in ('__init__', 'delete', 'setName',):
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.DataSource', func)
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.Experiment', func)
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.PeakList', func)

    for func in ('__init__', 'delete'):
      notifyFunc(self.updateMols, 'ccp.molecule.MolSystem.MolSystem',  func)
      notifyFunc(self.updateMols, 'ccp.molecule.MolSystem.Chain',      func)

    for func in ('__init__', 'delete'):
      notifyFunc(self.updateConstraintSets,'ccp.nmr.NmrConstraint.NmrConstraintStore', func)

    for func in ('__init__', 'delete', 'setName'):
      for clazz in ('ccp.nmr.NmrConstraint.DistanceConstraintList',
                    'ccp.nmr.NmrConstraint.HBondConstraintList'):
        notifyFunc(self.updateConstraintLists,clazz,func)
        
    
  def open(self):
  
    self.updateConstraintSets()
    self.updateMols()
    self.update()
    BasePopup.open(self)
    
  def update(self):
    
    self.updatePeakLists()
    self.updateConstraintLists()
  
  def togglePeakList(self, peakList):
  
    if peakList in self.peakLists:
      self.peakLists.remove(peakList)
    else:
      self.peakLists.append(peakList)

    self.updatePeakLists()
  
  def toggleConstraintList(self, constraintList):
  
    if constraintList in self.constraintLists:
      self.constraintLists.remove(constraintList)
    else:
      self.constraintLists.append(constraintList)

    self.updateConstraintLists()

  def getMolKeys(self):
      
    names = []
    for molSystem in self.project.sortedMolSystems():
      for chain in molSystem.sortedChains():
        if chain.residues:
          name = '%s:%s' % (molSystem.code, chain.code)
          names.append( (name,chain) )
        
    return names

  def updateMols(self):

    molKeys = self.getMolKeys()
    names = [x[0] for x in molKeys]
    chains = [x[1] for x in molKeys]

    indexX = 0
    indexY = 0
    
    if chains:
      if self.xMol not in chains:
        self.xMol = chains[0]
      
      if self.yMol not in chains:
        self.yMol = chains[0]

      indexX = chains.index(self.xMol)
      indexY = chains.index(self.yMol)
      
    self.xMolPulldown.setup(names, chains, indexX)
    self.yMolPulldown.setup(names, chains, indexY)

  def changeMolX(self, chain):
  
    self.xMol = chain
    residues = chain.sortedResidues()
    self.xEntryStart.set(residues[0].seqCode)
    self.xEntryStop.set(residues[-1].seqCode)
    self.setMolRanges()
       
  def changeMolY(self, chain):

    self.yMol = chain
    residues = chain.sortedResidues()
    self.yEntryStart.set(residues[0].seqCode)
    self.yEntryStop.set(residues[-1].seqCode)
    self.setMolRanges()

  def getPeakLists(self):
    # tbd noesy only
    
    peakLists = []
    for experiment in self.project.currentNmrProject.sortedExperiments():
      for spectrum in experiment.sortedDataSources():
        if (spectrum.dataType == 'processed') and (spectrum.numDim > 1):
          #isotopes = getSpectrumIsotopes(spectrum)
          #if isotopes.count('1H') > 1:
          #  spectra.append(spectrum)
          for peakList in spectrum.peakLists:
            peakLists.append(peakList)

    return peakLists
  
  def zoomIn(self):
  
    z = self.noeMatrix.zoom
    self.noeMatrix.setZoom(z*1.4)
  
  def zoomOut(self):
  
    z = self.noeMatrix.zoom
    self.noeMatrix.setZoom(z*0.7)

  def changeConstraintSet(self, constraintSet):    

    if constraintSet is not self.constraintSet: 
      self.constraintSet = constraintSet
      self.updateConstraintLists()   

  def updateConstraintSets(self, *opt):
    
    index = 0
    constraintSet = self.constraintSet
    constraintSets = self.nmrProject.sortedNmrConstraintStores()
    constraintSetNames = ['%d' % cs.serial for cs in constraintSets]
    
    if constraintSets:

      if constraintSet not in constraintSets:
        constraintSet = constraintSets[0]
          
      index = constraintSets.index(constraintSet)
 
    else:
      constraintSet = None
    
    self.constraintSetPulldown.setup(constraintSetNames, constraintSets, index)

    if self.constraintSet is not constraintSet:
      self.constraintSet = constraintSet
      self.updateConstraintLists()

  def getConstraintLists(self):
  
    constraintLists = []
    if self.constraintSet:
      for constraintList in self.constraintSet.constraintLists:
        if constraintList.className in ('DistanceConstraintList','HBondConstraintList'):
          constraintLists.append( constraintList )
    
    return constraintLists    
    
  def updateConstraintLists(self, constraintList=None):
  
    if constraintList and (constraintList.parent is not self.constraintSet):
      return
  
    textMatrix = []
    colorMatrix = []
    objectList = []
    
    if self.constraintSet:
      objectList = self.getConstraintLists()
      
      for constraintList in objectList:
        name = '%d:%s' % (constraintList.serial, constraintList.name)
        lType = constraintList.className[:-14]
    
        if constraintList in self.constraintLists:
          isUsed = 'Yes'
          colors = [None,None,'#C0FFC0']
 
        else:
          isUsed = 'No'
          colors = [None,None,None]
 
        datum = [name,lType,isUsed]
        textMatrix.append(datum)
        colorMatrix.append(colors)
      
    self.constraintMatrix.update(objectList=objectList,
                                 textMatrix=textMatrix,
                                 colorMatrix=colorMatrix)
  
  def updatePeakLists(self, *opt):
  
    textMatrix = []
    colorMatrix = []
    objectList = self.getPeakLists()
    for peakList in objectList:
      
      spectrum = peakList.dataSource
      experiment = spectrum.experiment
      
      if peakList in self.peakLists:
        isUsed = 'Yes' 
        colors = [None,None,'#C0FFC0']
     
      else:
        isUsed = 'No'
        colors = [None,None,None]
      
      datum = ['%s:%s' % (experiment.name,spectrum.name),peakList.serial,isUsed]
      textMatrix.append(datum)
      colorMatrix.append(colors)
      
    self.peakMatrix.update(objectList=objectList,textMatrix=textMatrix,colorMatrix=colorMatrix)
  
  def setMolRanges(self, *opt):
  
    xStart = self.xEntryStart.get() or None
    
    xStop = self.xEntryStop.get() or None
    
    yStart = self.yEntryStart.get() or None
    
    yStop = self.yEntryStop.get() or None

    #if xStop < xStart:
    #  (xStart,xStop) = (xStop,xStart)

    #if yStop < yStart:
    #  (yStart,yStop) = (yStop,yStart)
          
    if self.xMol and self.yMol:
      xChain = self.xMol
      yChain = self.yMol
      xResidues = xChain.sortedResidues()
      yResidues = yChain.sortedResidues()
      xFirst = xResidues[0].seqCode
      xLast  = xResidues[-1].seqCode
      yFirst = yResidues[0].seqCode
      yLast  = yResidues[-1].seqCode
      
      if xStart is None:
        xStart = xFirst
      if yStart is None:
        yStart = yFirst
      if xStop is None:
        xStop = xLast
      if yStop is None:
        yStop = yLast
      
      xStart = max(xStart, xFirst)
      yStart = max(yStart, yFirst)
      xStop  = min(xStop, xLast)
      yStop  = min(yStop, yLast)
      
    self.xEntryStart.set(xStart)
    self.yEntryStart.set(yStart)
    self.xEntryStop.set(xStop)
    self.yEntryStop.set(yStop)
    
    return (xStart,xStop,yStart,yStop)
 
  def updateNoes(self):
  
    #if not self.peakLists:
    #  return
  
    if not (self.noeMatrix and self.xMol and self.yMol):
      return  

    xChain = self.xMol
    yChain = self.yMol
    
    (xStart,xStop,yStart,yStop) = self.setMolRanges()
    
    xResidues = []
    for residue in xChain.sortedResidues():
      if (residue.seqCode >= xStart) and (residue.seqCode <= xStop):
        xResidues.append(residue)
    
    yResidues = []
    for residue in yChain.sortedResidues():
      if (residue.seqCode >= yStart) and (residue.seqCode <= yStop):
        yResidues.append(residue)
    
    self.noeMatrix.peakLists = self.peakLists
    self.noeMatrix.restraintLists = self.constraintLists
    self.noeMatrix.xResidues = xResidues
    self.noeMatrix.yResidues = yResidues
    self.noeMatrix.updateAfter()
      
def destroy(self):

    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self) 
    
