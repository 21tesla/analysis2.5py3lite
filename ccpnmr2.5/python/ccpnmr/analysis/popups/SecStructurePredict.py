"""
======================COPYRIGHT/LICENSE START==========================

SecStructurePredict.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2012 Wayne Boucher and Tim Stevens (University of Cambridge)

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

from memops.gui.Label               import Label
from memops.gui.LabelFrame          import LabelFrame
from memops.gui.ScrolledMatrix      import ScrolledMatrix
from memops.gui.PulldownList        import PulldownList
from memops.gui.MessageReporter     import showError, showInfo

from memops.editor.Util             import createDismissHelpButtonList

from ccpnmr.analysis.popups.BasePopup       import BasePopup

from ccpnmr.analysis.core.AssignmentBasic   import getShiftLists

from ccpnmr.analysis.wrappers.D2D           import SEC_STRUC_KEYS, SEC_STRUC_TIPS, runD2D

class SecStructurePredictPopup(BasePopup):
  """
  **Predict Protein Secondary Structure**
  
  This popup window is designed to allow the prediction of secondary structure
  for a protein chain given chemical shifts, using the (external) program D2D.

  The Options to select are the Chain and the Shift List, for which the
  prediction is then made.

  The Secondary Structure Predictions table lists the residues in the chain.
  For each residue the residue number, residue type and current secondary
  structure set for that residue is given.  The remaining columns are for
  the predictions made by D2D, and of course are only filled in once D2D
  is run.  The predicted secondary structure is listed first, followed by
  the probability of that residue being Helix, Beta, Coil or PPII (the
  predicted secondary structure will be specified by the maximum of these).
  
  To run the prediction click on the "Run D2D Prediction!" button.  This
  does not store this information in the project.  To do that you have to
  click on the "Commit Predicted Secondary Structure" button.

  **Caveats & Tips**

  The predicted secondary structure cell is coloured red if the prediction
  is unreliable.  Unreliable predictions are not stored with the "Commit"
  button but all reliable ones are.  If you need to edit the secondary
  structure for a residue then use the Secondary Structure Chart:

  .. _Secondary Structure Chart: SecStructureGraphPopup.html

  **References**

  The D2D programme:

  http://www-vendruscolo.ch.cam.ac.uk/d2D/index.php

  *C. Camilloni, A. De Simone, W. Vranken and M. Vendruscolo.
  Determination of Secondary Structure Populations in Disordered States of Proteins using NMR Chemical Shifts.
  Biochemistry 2012, 51: 2224-2231
  """

  def __init__(self, parent, *args, **kw):

    self.chain     = None
    self.shiftList = None
    self.predictionDict = {}
     
    BasePopup.__init__(self, parent=parent, title='Structure : Predict Secondary Structure')

  def body(self, guiFrame):

    self.geometry('700x500')
   
    guiFrame.expandGrid(1,0)
    
    row = 0
    
    # TOP LEFT FRAME
    
    frame = LabelFrame(guiFrame, text='Options')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.columnconfigure(5, weight=1)
        
    label = Label(frame, text='Chain')
    label.grid(row=0, column=0, sticky='w')
    self.chainPulldown = PulldownList(frame, callback=self.changeChain,
                                      tipText='Choose the molecular system chain to make predictions for')
    self.chainPulldown.grid(row=0, column=1, sticky='w')
    
    label = Label(frame, text='Shift List')
    label.grid(row=0, column=2, sticky='w')
    self.shiftListPulldown = PulldownList(frame, callback=self.changeShiftList,
                                          tipText='Select the shift list to take input chemical shifts from')
    self.shiftListPulldown.grid(row=0, column=3, sticky='w')
    
    row += 1
    
    # BOTTOM LEFT FRAME
    
    frame = LabelFrame(guiFrame, text='Secondary Structure Predictions')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)
    
    tipTexts = ('Residue number in chain',
                'Residue type code',
                'Current stored secondary structure code',
                'Predicted secondary structure code') + SEC_STRUC_TIPS

    headingList = ('Res\nNum','Res\nType','Current\nSS', 'Predicted\nSS') + SEC_STRUC_KEYS
                          
    n = len(headingList)
    editWidgets        = n * [None]
    editGetCallbacks   = n * [None]
    editSetCallbacks   = n * [None]
                          
    self.predictionMatrix = ScrolledMatrix(frame, 
                                           headingList=headingList, 
                                           tipTexts=tipTexts,
                                           editWidgets=editWidgets,
                                           editGetCallbacks=editGetCallbacks,
                                           editSetCallbacks=editSetCallbacks)
    self.predictionMatrix.grid(row=0, column=0, sticky='nsew')
    
    row += 1

    tipTexts = ['Run the D2D method to predict secondary structure',
                'Store the secondary structure predictions in the CCPN project']
                
    texts = ['Run D2D Prediction!','Commit Predicted\nSecondary Structure']
    commands = [self.runD2D, self.storeSecondaryStructure]
    self.buttonList = createDismissHelpButtonList(guiFrame, texts=texts, commands=commands,
                                                 help_url=self.help_url,
                                                 expands=True, tipTexts=tipTexts)
    self.buttonList.grid(row=row, column=0, columnspan=2, sticky='ew')
    
    self.update()
 
    self.notify(self.registerNotify)

  def destroy(self):
      
    self.notify(self.unregisterNotify)
    BasePopup.destroy(self)
 
  def notify(self, notifyfunc):
     
    for func in ('__init__', 'delete'):
      notifyfunc(self.updateChainPulldown, 'ccp.molecule.MolSystem.Chain', func)
      
    for func in ('__init__', 'delete', 'setName'):
      notifyfunc(self.updateShiftListPulldown, 'ccp.nmr.Nmr.ShiftList', func)
      
    for func in ('setSecStrucCode',):
      notifyfunc(self.updatePredictionMatrixAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)

  def update(self):
    
    self.updateShiftListPulldown()
    self.updateChainPulldown()
    self.updatePredictionMatrixAfter()
  
  def runD2D(self):
    
    chain     = self.chain
    shiftList = self.shiftList
    
    if (not chain) or (not shiftList):
      showError('Cannot Run D2D', 'Please specify a chain and a shift list.', parent=self)
      return
      
    self.predictionDict[chain] = runD2D(chain, shiftList)

    self.updatePredictionMatrix()
   
  def storeSecondaryStructure(self):
  
    if not self.chain:
      return
    
    getSpinSystem = self.nmrProject.findFirstResonanceGroup
    newSpinSystem = self.nmrProject.newResonanceGroup
    
    predDict = self.predictionDict.get(self.chain, {})

    n = 0
    for residue in predDict:
      (ssCode, isReliable, probabilityDict) = predDict[residue]
      
      if isReliable:
        spinSystem = getSpinSystem(residue=residue)
      
        if not spinSystem:
          spinSystem = newSpinSystem(residue=residue, ccpCode=residue.ccpCode)
      
        spinSystem.secStrucCode = ssCode
        n += 1
      
    showInfo('Info', 'Stored secondary structure types for %d residues.' % n, parent=self)
   
  def updatePredictionMatrixAfter(self, index=None, text=None):
    
    self.after_idle(self.updatePredictionMatrix)
    
  def updatePredictionMatrix(self):

    objectList  = []
    textMatrix  = []
    colorMatrix = []

    n = len(SEC_STRUC_KEYS)
    chain = self.chain
    if chain:
      predDict = self.predictionDict.get(chain, {})

      getSpinSystem = self.nmrProject.findFirstResonanceGroup
      for residue in chain.sortedResidues():
        spinSystem = getSpinSystem(residue=residue)
        currentSsCode = spinSystem and spinSystem.secStrucCode

        data = [residue.seqCode, residue.ccpCode, currentSsCode]
        colors = 3*[None]
        if residue in predDict:
          (ssCode, isReliable, probabilityDict) = predDict[residue]
          data.append(ssCode)
          if isReliable:
            colors.append(None)
          else:
            colors.append('#FF3333')
          for key in SEC_STRUC_KEYS:
            data.append(probabilityDict[key])
          colors.extend(n*[None])
        else:
          data.extend((n+1)*[None])
          colors.extend((n+1)*[None])

        textMatrix.append(data)
        objectList.append(residue)
        colorMatrix.append(colors)
      
    self.predictionMatrix.update(textMatrix=textMatrix,
                                 objectList=objectList,
                                 colorMatrix=colorMatrix)
    
  def changeChain(self, chain):
    
    if chain is not self.chain:
      self.chain = chain
      self.updatePredictionMatrixAfter()
    
  def changeShiftList(self, shiftList):
  
    if shiftList is not self.shiftList:
      self.shiftList = shiftList
      self.updatePredictionMatrixAfter()
    
  def updateChainPulldown(self, obj=None):
    
    index = 0
    names = []
    chains = []
    chain = self.chain
    
    for molSystem in self.project.molSystems:
      msCode = molSystem.code

      for chainA in molSystem.chains:
        residues = chainA.residues

        if not residues:
          continue

        for residue in residues:
          # Must have at least one protein residue
          if residue.molType == 'protein':
            names.append('%s:%s' % (msCode, chainA.code))
            chains.append(chainA)
            break
 
    if chains:
      if chain not in chains:
        chain = chains[0]
        
      index = chains.index(chain)
        
    else:
      chain = None
    
    if chain is not self.chain:
      self.chain = chain
      self.updatePredictionMatrixAfter()
        
    self.chainPulldown.setup(names, chains, index)
    
  def updateShiftListPulldown(self, obj=None):
    
    index = 0
    names = []
    shiftLists = getShiftLists(self.nmrProject)
        
    if shiftLists:
      if self.shiftList not in shiftLists:
        self.shiftList = shiftLists[0]
        
      index = shiftLists.index(self.shiftList)
      names = ['%s:%d' % (sl.name,sl.serial) for sl in shiftLists]  
    
    else:
      self.shiftList = None
      
    self.shiftListPulldown.setup(names, shiftLists, index)
    
