"""
======================COPYRIGHT/LICENSE START==========================

SequenceShiftPredict.py: Part of the CcpNmr Analysis program

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

from ccpnmr.analysis.core.AssignmentBasic   import getShiftLists, findFirstAtomShiftInShiftList

from ccpnmr.analysis.wrappers.CamCoil       import runCamCoil, SCRIPT_TEXTS, SCRIPT_PHS, LFP_SCRIPT

class SequenceShiftPredictPopup(BasePopup):
  """
  **Predict Protein Shifts from Sequence**
  
  This popup window is designed to allow the prediction of chemical shifts
  for a protein chain from the sequence (so with no structural information),
  using the (external) program CamCoil.

  The Options to select are the Chain for which the prediction is made, and
  the prediction type and the pH used for the prediction, and also the Shift
  List, which is not used for the prediction but is used for the comparison
  with the prediction.

  CamCoil has two variations, one for the prediction of random coil chemical
  shifts and one for prediction of protein loops chemical shifts.

  The Chemical Shift Predictions table lists the atoms in the chain.
  For each atom the data listed is the residue number, residue type, atom
  name, first shift found for that atom in the chosen shiftList, chemical
  shift predicted by CamCoil, and the difference between the actual shift
  and the predicted shift (if both exist).
  
  To run the prediction click on the "Run CamCoil Prediction!" button.  This
  does not store any predicted shifts in the project.

  **Caveats & Tips**

  **References**

  The CamCoil programme:

  http://www-vendruscolo.ch.cam.ac.uk/camcoil.php

  *A. De Simone, A. Cavalli, S-T. D. Hsu, W. Vranken and M. Vendruscolo
  Accurate random coil chemical shifts from an analysis of loop regions in native states of proteins.
  J. Am. Chem. Soc. 131(45):16332-3
  """

  def __init__(self, parent, *args, **kw):

    self.chain = None
    self.shiftList = None
    self.predictionDict = {}
     
    BasePopup.__init__(self, parent=parent, title='Data Analysis : Predict Shifts from Sequence')

  def body(self, guiFrame):

    self.geometry('700x500')
   
    guiFrame.expandGrid(1,0)
    
    row = 0
    
    # TOP LEFT FRAME
    
    frame = LabelFrame(guiFrame, text='Options')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.columnconfigure(7, weight=1)
        
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

    label = Label(frame, text='Type')
    label.grid(row=0, column=4, sticky='w')
    self.scriptPulldown = PulldownList(frame, texts=SCRIPT_TEXTS,
                                       callback=self.changeScript,
                                       tipText='Select the algorithm script for this chain')
    self.scriptPulldown.grid(row=0, column=5, sticky='w')
    
    self.pHLabel = Label(frame, text='pH')
    self.pHLabel.grid(row=0, column=6, sticky='w')
    self.pHPulldown = PulldownList(frame, texts=SCRIPT_PHS,
                        tipText='Select the pH to make the prediction for')
    self.pHPulldown.grid(row=0, column=7, sticky='w')
    
    row += 1
    
    # BOTTOM LEFT FRAME
    
    frame = LabelFrame(guiFrame, text='Chemical Shift Predictions')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)
    
    tipTexts = ['Residue number in chain',
                'Residue type code',
                'Atom name',
                'Actual shift (first one it finds for atom in chosen shiftList)',
                'CamCoil predicted shift',
                'Predicted - Actual']

    headingList        = ['Res\nNum','Res\nType','Atom\nName',
                          'Actual\nShift','Predicted\nShift','Difference']
                          
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

    tipTexts = ['Run the CamCoil method to predict chemical shifts from sequence']
                
    texts = ['Run CamCoil Prediction!']
    commands = [self.runCamCoil]
    self.buttonList = createDismissHelpButtonList(guiFrame, texts=texts, commands=commands,
                                                 help_url=self.help_url,
                                                 expands=True, tipTexts=tipTexts)
    self.buttonList.grid(row=row, column=0)
    
    self.update()
 
    self.notify(self.registerNotify)

  def destroy(self):
      
    self.notify(self.unregisterNotify)
    BasePopup.destroy(self)
 
  def notify(self, notifyfunc):
     
    for func in ('__init__', 'delete'):
      notifyfunc(self.updateChainPulldown, 'ccp.molecule.MolSystem.Chain', func)
      
    for func in ('setValue',):
      notifyfunc(self.updatePredictionMatrixAfter, 'ccp.nmr.Nmr.Shift', func)
      
  def update(self):
    
    self.updateShiftListPulldown()
    self.updateChainPulldown()
    self.updatePredictionMatrixAfter()
  
  def runCamCoil(self):
    
    chain = self.chain
    script = self.scriptPulldown.getText()
    if script == LFP_SCRIPT:
      pH = ''
    else:
      pH = self.pHPulldown.getText()
    
    if not chain:
      showError('Cannot Run CamCoil', 'Please specify a chain.', parent=self)
      return
      
    self.predictionDict[chain] = runCamCoil(chain, pH=pH, script=script)

    self.updatePredictionMatrix()
   
  def updatePredictionMatrixAfter(self, index=None, text=None):
    
    self.after_idle(self.updatePredictionMatrix)
    
  def updatePredictionMatrix(self):

    objectList  = []
    textMatrix  = []

    chain = self.chain
    shiftList = self.shiftList
    if chain:
      atomShiftDict = self.predictionDict.get(chain, {})

      for residue in chain.sortedResidues():
        for atom in residue.sortedAtoms():
          currShift = shiftList and findFirstAtomShiftInShiftList(atom, shiftList)
          value = currShift and currShift.value
          predShift = atomShiftDict.get(atom)
          if currShift and predShift:
            delta = predShift - value
          else:
            delta = None
          data = [residue.seqCode, residue.ccpCode, atom.name, value, predShift, delta]

          textMatrix.append(data)
          objectList.append(atom)
      
    self.predictionMatrix.update(textMatrix=textMatrix,
                                 objectList=objectList)
    
  def changeChain(self, chain):
    
    if chain is not self.chain:
      self.chain = chain
      self.updatePredictionMatrixAfter()
    
  def changeShiftList(self, shiftList):

    if shiftList is not self.shiftList:
      self.shiftList = shiftList
      self.updatePredictionMatrixAfter()

  def changeScript(self, script):
  
    if script == LFP_SCRIPT:
      self.pHLabel.grid_forget()
      self.pHPulldown.grid_forget()
    else:
      self.pHLabel.grid(row=0, column=6, sticky='w')
      self.pHPulldown.grid(row=0, column=7, sticky='w')

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
 
