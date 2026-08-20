"""
======================COPYRIGHT/LICENSE START==========================

SpinSystemTyping.py: Part of the CcpNmr Analysis program

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
from memops.gui.ButtonList   import ButtonList, UtilityButtonList
from memops.gui.CheckButton  import CheckButton
from memops.gui.FloatEntry   import FloatEntry
from memops.gui.Frame        import Frame
from memops.gui.IntEntry     import IntEntry
from memops.gui.LabelFrame   import LabelFrame
from memops.gui.Label        import Label
from memops.gui.ProgressBar  import ProgressBar
from memops.gui.PulldownList import PulldownList


from memops.gui.ScrolledMatrix      import ScrolledMatrix
from memops.gui.PartitionedSelector import PartitionedSelector
from memops.gui.ScrolledGraph       import ScrolledGraph
from memops.gui.MessageReporter     import showOkCancel

from ccpnmr.analysis.popups.BasePopup        import BasePopup
from ccpnmr.analysis.core.ChemicalShiftBasic import getShiftsChainProbabilities, lookupAtomProbability
from ccpnmr.analysis.core.SpinSystemTyping import getSpinSystemTypes
from ccpnmr.analysis.core.AssignmentBasic  import getResonanceName, assignSpinSystemType, assignSpinSystemResidue, \
                                                  deassignResonance, getShiftLists, removeSpinSystemResonance, \
                                                  assignResonanceType
from ccpnmr.analysis.core.MoleculeBasic import DEFAULT_ISOTOPES, getResidueCode

COLOR_DICT = {'1H':'#a0d0a0','15N':'#a0a0d0','13C':'#d0d0a0'}

class SpinSystemTypeScoresPopup(BasePopup):
  """
  **Predict Residue Type for a Spin System of Resonances**
  
  This tool aims to predict the residue type of a spin system based upon the
  chemical shifts of the resonances that it contains. The general principle is
  that different kinds of atoms in different kinds of residues have different
  observed distributions of chemical shifts. This system uses chemical shift
  distributions from the RefDB database, or otherwise from the BMRB where data
  is not available in RefDB. The observed chemical shifts of a spin system are
  compared to the per-atom distributions for each residue type and the residue
  types with the best matches are deemed to be more likely.

  This system can work with various levels of information, although the more
  information the better. Naturally, the more chemical shifts you have in a spin
  system then the better the prediction of type, and 13C resonances are more
  distinctive than 1H on the whole. Also, setting the atom type of a resonance
  can have a big influence on the type of residue predicted, for example knowing
  that a 13C resonance at 63 ppm is of type CB points  very strongly toward the
  residue being a serine. Atom type information can come from two sources: from
  a specific type assignment made by the user (via this popup or elsewhere) or
  by virtue of assignment in an experimental dimension that detects a 
  restricted class of atom - e.g. 13C resonances in an HNCA experiment, assuming
  their shift matches, are of CA type as far as this prediction is concerned.
  Resonances that do not have a known atom type are compared with all of the
  unallocated types to find the combination that is most likely.

  The residue type prediction is based on the list of resonances displayed in the
  upper table. Here the user can see the chemical shifts (from the selected
  shift list) and any specific atom type setting. The user may set the atom type
  for any of the resonances, which would normally be done to reduce prediction
  ambiguity, by double-clicking in the "Atom Type" column.

  The lower table shows a ranked list of the probable residue types. All
  probability scores are normalised and represented as a percentage of the total
  of all scores, considering residue types in the selected chain. The type of a
  spin system may be set by clicking on a row of the lower table (hopefully a
  unique and high-scoring option) and  then selecting [Assign Spin System Type].
  If the user attempts to change the type of a spin system that is currently
  assigned to a specific residue then there is an opportunity to back out of the
  assignment, but otherwise any sequence specific information will be removed.

  **Caveats & Tips**

  It is assumed that the spectra from which the chemical shifts are derived are
  fairly well referenced. 

  A type prediction will always be given, no matter how few resonances are
  present in a spin system. This system says which of the available types are
  most likely, *not how reliable* the prediction is; the latter depends largely
  on the amount of information present. The user should not for example make a
  judgement based only on amide resonances. Reliability scores will be added in
  the future.

  Rouge resonances in a spin system often adversely affect the prediction, if
  something is not genuinely in the spin system it should be removed.

  The system will never predict the residue type to be something that does not
  appear in the selected molecular chain. Thus, make sure the chain selection is
  appropriate for your prediction.

  **Reference**
  
  The residue type prediction method is not published independently but is very
  similar to the Bayesian method presented in: *Marin A, Malliavin TE, Nicolas P,
  Delsuc MA. From NMR chemical shifts to amino acid types: investigation of the
  predictive power carried by nuclei. J Biomol NMR. 2004 Sep;30(1):47-60.*

  One major difference however is that probabilities for resonances not being
  observed are not used. The CCPN prediction method is not only for complete
  spin systems and may be used at any time during the assignment process; here
  missing resonances are mostly due to the current assignment state and not such
  a useful indicator of residue type. """

  def __init__(self, parent, spinSystem=None, chain=None, *args, **kw):
  
    self.spinSystem = spinSystem
    self.shiftList  = None
    self.resonance  = None
    self.isotopes   = ('1H','13C','15N')
    self.chain      = chain
    self.ccpCode    = None
    self.waiting    = False
    self.atomTypes  = {}
  
    self.project   = parent.project
    self.guiParent = parent
  
    BasePopup.__init__(self, parent, title="Spin System Type Scores", **kw)


  def body(self, guiFrame):
  
    guiFrame.grid_columnconfigure(3, weight=1)
    
    row = 0
    label = Label(guiFrame, text='Spin System: ', grid=(row,0))
    tipText = 'Indicates which spin system the residue type prediction is done for'
    self.spinSystemLabel = Label(guiFrame, text='Serial:   Assignment:',
                                 grid=(row,1), gridSpan=(1,3), tipText=tipText)

    row += 1
    label = Label(guiFrame, text='Shift List: ', grid=(row,0))
    tipText = 'Selects which shift list is the source of chemical shift information to make the residue type prediction'
    self.shiftListPulldown = PulldownList(guiFrame, tipText=tipText,
                                          callback=self.setShiftList,
                                          grid=(row,1))

    label = Label(guiFrame, text='Chain: ', grid=(row,2))
    tipText = 'Selects which molecular chain the prediction is for; sets prior probabilities for the various residue types'
    self.chainPulldown = PulldownList(guiFrame, self.changeChain,
                                      grid=(row,3), tipText=tipText)

    row += 1
    labelFrame = LabelFrame(guiFrame, text='Resonances', grid=(row,0), gridSpan=(1,4))
    labelFrame.expandGrid(0,0)
    
    self.atomTypePulldown = PulldownList(self, callback=self.setAtomType)
    
    editWidgets = [ None, None, None, None, self.atomTypePulldown ]
    editGetCallbacks = [ None, None, None, None, self.getAtomType]
    editSetCallbacks = [ None, None, None, None, self.setAtomType]
    
    tipTexts = ['The nuclear isotope type of the resonance within the current spin system',
                'The assignment annotation for the spin system resonance within the current spin system',
                'The chemical shift of the resonance in the stated shift list',
                'The weighted standard deviation of the resonance chemical shift',
                'The current atom type of the resonance; when set this helps refine residue type prediction']
    headingList = ['Isotope','Name','Shift\nValue','Shift\nError','Atom\nType']
    self.resonanceMatrix = ScrolledMatrix(labelFrame,
                                          editWidgets=editWidgets, multiSelect=False,
                                          editGetCallbacks=editGetCallbacks,
                                          editSetCallbacks=editSetCallbacks,
                                          headingList=headingList,
                                          callback=self.selectResonance,
                                          grid=(0,0), tipTexts=tipTexts)

    tipTexts = ['Remove the selected resonance from the current spin system',
                'Remove residue type information from the current spin system',
                'Show a table of information for the  selected resonance, including a list of all peak dimension positions',
                'Show a table of the peaks to which the selected resonance is assigned']
    texts = ['Remove From\nSpin System', 'Deassign\nResidue Type',
             'Resonance\nInfo', 'Show\nPeaks']
    commands = [self.removeResonance, self.deassignType,
                self.showResonanceInfo, self.showPeaks]
    buttonList = ButtonList(labelFrame, texts=texts, commands=commands,
                            grid=(1,0), tipTexts=tipTexts)
    self.resButtons = buttonList.buttons

    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    labelFrame = LabelFrame(guiFrame, text='Type Scores', grid=(row,0), gridSpan=(1,4))
    labelFrame.expandGrid(0,0)
    
    tipTexts = ['The ranking of the residue type possibility for the current spin system',
                'The CCPN residue code for the type',
                'The estimated percentage probability of the spin system being the residue type']
    headingList = ['Rank','Ccp Code','% Probability']
    self.scoresMatrix = ScrolledMatrix(labelFrame,
                                       headingList=headingList,
                                       callback=self.selectCcpCode,
                                       grid=(0,0), tipTexts=tipTexts)
 
    row += 1
    tipTexts = ['Assign the residue type of the current spin system to the kind selected in the lower table',]
    texts    = ['Assign Spin System Type']
    commands = [self.assign]
    bottomButtons = UtilityButtonList(guiFrame, texts=texts, commands=commands,
                                      helpUrl=self.help_url, grid=(row,0),
                                      gridSpan=(1,4), tipTexts=tipTexts)
    self.assignButton = bottomButtons.buttons[0]

    self.updateShiftLists()
    self.updateChains()
    self.getChainAtomTypes()
    self.update()
  
    self.curateNotifiers(self.registerNotify)
  
  def curateNotifiers(self, notifyFunc):

    for func in ('addResonance','removeResonance',
                 'setResonances','delete','setName'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)

    for func in ('setResonanceSet','addAssignName',
                 'removeAssignName','setAssignNames'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.Resonance', func)

    for func in ('__init__','delete'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.ResonanceSet', func)
      notifyFunc(self.updateChains, 'ccp.molecule.MolSystem.Chain', func)
      notifyFunc(self.updateShiftLists, 'ccp.nmr.Nmr.ShiftList', func)
   
    for func in ('__init__','setValue'):
      notifyFunc(self.updateShiftAfter, 'ccp.nmr.Nmr.Shift', func)
  
  def getChainAtomTypes(self):
  
    doneResType = {}
    
    atomTypes = atomTypes = {}
    for isotope in DEFAULT_ISOTOPES.values():
      atomTypes[isotope] = set()
    
    if self.chain:
      for residue in self.chain.residues:
        molResidue = residue.molResidue
        ccpCode = molResidue.ccpCode
        molType = molResidue.molType
        key = '%s:%s:%s' % (ccpCode, molResidue.linking, molResidue.descriptor)
        
        if doneResType.get(key):
          continue
          
        doneResType[key] = True
          
        for atom in residue.atoms:
          chemAtom = atom.chemAtom
          element = chemAtom.elementSymbol
          isotope = DEFAULT_ISOTOPES.get(element)
          
          if not isotope:
            continue
            
          atomTypes[isotope].add((ccpCode, atom.name, molType))
   
    self.atomTypes = atomTypes
  
  def getAtomType(self, resonance):
  
    index = 0
    atomNames = set(['<None>',])
    
    assignNames = resonance.assignNames
    if assignNames:
      for atomName in assignNames:
        atomNames.add(atomName)
      
      if len(assignNames) > 1:
        orig = ','.join(assignNames)
        atomNames.add(orig)
    
    shift = resonance.findFirstShift(parentList=self.shiftList)
    
    if shift and self.chain:
    
      project = self.project
      atomTypes = self.atomTypes.get(resonance.isotopeCode, [])
 
      for ccpCode, atomName, molType in atomTypes:
         prob = lookupAtomProbability(project, ccpCode, atomName,
                                      shift.value, molType=molType)
  
         if prob >= 0.001:
           atomNames.add(atomName)
 
    atomNames = list(atomNames)
    atomNames.sort()
 
    if resonance.assignNames:
      orig = ','.join(assignNames)
      index = atomNames.index(orig)
  
    atomNameObjs = atomNames[:]
    atomNameObjs[0] = None
  
    self.atomTypePulldown.setup(atomNames, atomNameObjs, index)
  
  def setAtomType(self, obj):
  
    atomNameStr = self.atomTypePulldown.getObject()
  
    if self.resonance:
      if atomNameStr:
        atomNames = atomNameStr.split(',')
        assignResonanceType(self.resonance, assignNames=atomNames)
      else:
        assignResonanceType(self.resonance, assignNames=None)
  
  def removeResonance(self):
  
    if self.resonance and self.spinSystem and (self.resonance in self.spinSystem.resonances):
      if showOkCancel('Confirm','Really remove resonance from spin system?', parent=self):
        self.spinSystem.codeScoreDict = {}
        deassignResonance(self.resonance, clearAssignNames=False)
        removeSpinSystemResonance(self.spinSystem, self.resonance)

  def showResonanceInfo(self):
  
    if self.resonance:
      self.guiParent.browseResonanceInfo(self.resonance)

  def showPeaks(self):
  
    if self.resonance:
      peaksDict = {}
      for contrib in self.resonance.peakDimContribs:
        peaksDict[contrib.peakDim.peak] = 1

      peaks = peaksDict.keys()
      if len(peaks) > 0:
        self.guiParent.viewPeaks(peaks)

  def deassignType(self):
  
   if self.spinSystem:
     residue = self.spinSystem.residue
     
     if residue:
       resText = '%d%s' % (residue.seqCode, residue.ccpCode)
       msg = 'Spin system assigned to %s. Continue and deassign residue?' 
       if showOkCancel('Warning', msg % resText, parent=self):
         assignSpinSystemResidue(self.spinSystem, None)
         assignSpinSystemType(self.spinSystem, None)     
  
     else:
       assignSpinSystemType(self.spinSystem,None)     

  def assign(self):
  
    if self.spinSystem and self.ccpCode:
       if self.spinSystem.residue and (self.spinSystem.residue.ccpCode != self.ccpCode):
         resText = '%d%s' % (self.spinSystem.residue.seqCode, self.spinSystem.residue.ccpCode)
         msg = 'Spin system is already assigned to %s. Continue?'
         if showOkCancel('Warning', msg % resText, parent=self):
           assignSpinSystemResidue(self.spinSystem,residue=None)
           
         else:
           return
    
       if self.spinSystem.ccpCode != self.ccpCode:
         assignSpinSystemType(self.spinSystem,self.ccpCode,'protein')
         self.update()


  def getChains(self):
  
  
    chains = []
    if self.project:
      for molSystem in self.project.sortedMolSystems():
        for chain in molSystem.sortedChains():
          if chain.molecule.molType in ('protein',None):
	    text = '%s:%s' % (molSystem.code, chain.code)
            chains.append( [text, chain] )
	
    return chains


  def changeChain(self, chain):
    
    if self.chain is not chain:
      self.chain = chain
      self.getChainAtomTypes()
      self.updateAfter()
         
    
  def updateChains(self, *chain):
  
    data = self.getChains()
    names = [x[0] for x in data]
    chains = [x[1] for x in data]
    chain = self.chain
    index = 0
 
    if chains:
      if chain not in chains:
        chain = chains[0]
     
      index = chains.index(chain)

    if chain is not self.chain:
      self.chain = chain
      self.getChainAtomTypes()
      self.updateAfter()
    
    self.chainPulldown.setup(names, chains, index)


  def updateShiftLists(self, *opt):
   
    shiftLists = getShiftLists(self.nmrProject)
    shiftList = self.shiftList
    names = ['%s [%d]' % (x.name or '<No name>', x.serial) for x in shiftLists]
    index = 0
    
    if names:
      if shiftList not in shiftLists:
        shiftList = shiftLists[0]
      
      index = shiftLists.index(shiftList)	
    
    if shiftList is not self.shiftList:
      self.shiftList = shiftList 
      self.updateAfter()
    
    self.shiftListPulldown.setup(names, shiftLists, index)
  
  
  def setShiftList(self, shiftList):
  
    if self.shiftList is not shiftList:
      self.shiftList = shiftList
      self.updateAfter()
 
  def updateButtons(self):
  
    if self.resonance:
      self.resButtons[0].enable()
      self.resButtons[2].enable()
      self.resButtons[3].enable()
    
    else:
      self.resButtons[0].disable()
      self.resButtons[2].disable()
      self.resButtons[3].disable()
    
    if self.spinSystem and (self.spinSystem.residue or self.spinSystem.ccpCode):
      self.resButtons[1].enable()
    else:
      self.resButtons[1].disable()
    
    if self.ccpCode and self.spinSystem:
      self.assignButton.enable()
    else:
      self.assignButton.disable()
 
  def selectResonance(self, resonance, row, col):
  
    self.resonance = resonance
    self.updateButtons()
  
  def selectCcpCode(self, ccpCode, row, col):
  
    self.ccpCode = ccpCode
    self.assignButton.enable()

  def clearSpinSystemCache(self):
  
    if self.spinSystem:
      self.spinSystem.sstTypes = []
      self.spinSystem.ssScore = None
      self.spinSystem.codeScoreDict = {}
      
 
  def updateShiftAfter(self, shift):
  
    if shift.parentList is not self.shiftList:
      return
      
    resonance = shift.resonance
    
    if resonance.resonanceGroup is not self.spinSystem:
      return
    
    self.clearSpinSystemCache()
    
    self.updateAfter()
       
 
  def updateAfter(self, obj=None):

    if obj and obj is self.spinSystem:
      self.clearSpinSystemCache()
      if obj.isDeleted:
        self.spinSystem = None
    
    if self.spinSystem:
      if obj is not None:
        if obj.className == 'ResonanceSet':
          for resonance in obj.resonances:
            if resonance.resonanceGroup is self.spinSystem:
              break
 
          else:
            return
 
        elif obj.className == 'Resonance':
          if obj.resonanceGroup is not self.spinSystem:
            return
      
      self.clearSpinSystemCache()
        
    if self.waiting:
      return
      
    else:
      self.waiting = True
      self.after_idle(self.update)

  def update(self, spinSystem=None, chain=None, shiftList=None):
  
    if spinSystem is not None:
      if spinSystem is not self.spinSystem:
        self.ccpCode = None
        self.spinSystem = spinSystem
        
      if chain:
        self.chain = chain
      else:
        self.chain = None
        
      self.updateChains()
      self.getChainAtomTypes()
    
    if shiftList:
      self.shiftList = shiftList
      self.updateShiftLists()   
 
    if self.resonance:
      if not self.spinSystem:
        self.resonance = None
        
      elif self.resonance.resonanceGroup is not self.spinSystem:
        self.resonance = None
    
    
    self.updateButtons()
  
    textMatrix = []
    objectList = []
    if self.spinSystem:
      
      if self.spinSystem.residue:
        resText = 'Assignment: %d%s' % (self.spinSystem.residue.seqCode, self.spinSystem.residue.ccpCode)
      elif self.spinSystem.ccpCode:
        resText = 'Type: %s' % self.spinSystem.ccpCode
      elif self.spinSystem.name:
        resText = 'Name: %s' % self.spinSystem.name
      else:
        resText = 'Unassigned'
    
      self.spinSystemLabel.set('Serial: %d  %s' % (self.spinSystem.serial, resText))
    
      for resonance in self.spinSystem.resonances:
        shift = resonance.findFirstShift(parentList=self.shiftList)
        if shift:
          datum = [resonance.isotopeCode,
                   getResonanceName(resonance),
                   shift.value,
                   shift.error,
                   '/'.join(resonance.assignNames)]

 
          objectList.append(resonance)
          textMatrix.append(datum)
 
    self.resonanceMatrix.update(textMatrix=textMatrix, objectList=objectList)

    textMatrix = []
    objectList = []
    colorMatrix = []
    if self.spinSystem and self.chain and self.shiftList:
      shifts = []
      for resonance in self.spinSystem.resonances:
        if resonance.isotopeCode in self.isotopes:
          shift = resonance.findFirstShift(parentList=self.shiftList)
          if shift:
            shifts.append(shift)
    
      scores = getShiftsChainProbabilities(shifts, self.chain)
      total = sum(scores.values())
      
      scoreList = []
      
      if total:
        ccpCodes  = self.getCcpCodes(self.chain)
        baseLevel = 100.0/len(ccpCodes)
        for ccpCode in ccpCodes:
          scoreList.append( (100.0*scores[ccpCode]/total, ccpCode) )
      
      scoreList.sort()
      scoreList.reverse()
      
      i = 0
      for score, ccpCode in scoreList:
        if not score:
          continue
      
        i += 1
        datum = [i, ccpCode,score]
        
        if score >= min(100.0,5*baseLevel):
          color = '#80ff80'
        elif score > 2*baseLevel:
          color = '#ffff80'
        elif score > baseLevel:
          color = '#ffc080'
        else:
          color = '#ff8080'
        
        colors = [color, color, color]
        objectList.append(ccpCode)
        textMatrix.append(datum)
        colorMatrix.append(colors)
 
    self.scoresMatrix.update(textMatrix=textMatrix, colorMatrix=colorMatrix, objectList=objectList)
    self.waiting = False

  def getCcpCodes(self, chain):
  
    ccpDict = {}
    for residue in chain.residues:
      ccpCode = residue.ccpCode
      
      #if (ccpCode == 'Cys') and (residue.descriptor == 'link:SG'):
      #  ccpCode = 'Cyss'
      
      ccpDict[ccpCode] = True

    ccpCodes = ccpDict.keys()
    ccpCodes.sort()
    
    return ccpCodes

  def destroy(self):

    self.curateNotifiers(self.unregisterNotify)
  
    BasePopup.destroy(self)

class SpinSystemTypingPopup(BasePopup):
  """
  **Predict Which Types of Residue Spin Systems Represent**
  
  This popup window uses chemical shift information, obtained from the
  resonances of a spin system group, to predict which kind of residue a spin
  system could be. Naturally, the more resonances/shifts there are in a spin
  system the better the prediction will be. Predictions are either made for a
  single spin system in isolation, or for a whole chain; by shuffling the
  residue to spin system mapping to give the optimum arrangement. Prediction for
  individual spin systems, as accessed by [Show Individual Classification], is
  covered in the `Spin System Type Scores`_ popup and will not be discussed
  here.

  .. _`Spin System Type Scores`: ../popups/SpinSystemTypeScoresPopup.html
  
  **Finding the Optimal Arrangement of Residue Types**
  
  This system attempts to find the best match of spin system to residue type by
  performing a Monte Carlo search to arrange each group of chemical shifts
  amongst the residue types found in the chain. The main principle is that the
  chain's residues dictate how many spin systems of a given type may be found,
  such that for example if there is only one Threonine residue then only one
  spin system may be predicted to be of Threonine type. This result comes
  naturally from shuffling the spin systems, with their chemical shifts, amongst
  the residue slots (disregarding sequence position). 

  The prediction is made my selecting the chain and shift list to use, then
  which kinds of isotope to consider, by toggling the relevant buttons, choosing
  various Monte Carlo search options and finally selecting [Run Typing]. The
  prediction may take some time to run, depending upon the number of residues
  and spin systems that are being matched, but gives a graphical output of the
  progress. If the final prediction looks good the "Highest Scoring Mappings"
  display may be closed and [Assign Types] may be used to set the residue types
  for all of the spin systems in the main table that match only a single type
  and have a score above the assignment threshold value. Spin systems that
  already have a type or full residue assignment will not be affected.

  The default Monte Carlo search options ought to be appropriate for a small
  protein (100 residues) with 1H and 13C chemical shift information. Increasing
  the number of search steps may help if the search does not converge; still
  swaps between optimal assignments toward the end of the search. Increasing the
  ensemble size (how many test mappings are optimised at the same time) may help
  if the prediction gets stuck in local minima; different runs predict different
  arrangement, but larger ensembles require more search steps to converge. Where
  a spin system has multiple residue types predicted these are the types that
  come out of the final ensemble for that set of shifts, i.e. at this point the
  ensemble of solutions differ.

  Overall, it should be noted that if a human being cannot readily predict the
  probably types of a spin system from its shifts alone, then this search tool
  cannot be expected to do a good job; it is merely an optimiser to address the
  problem of shuffling spin systems within a chain.

  The scores are currently unnormalised log-odds values and are not especially
  meaningful in the human sense, other than higher is better (closer to zero
  for negative values). This issue will be addressed in the future.
  Spin systems without a unique type prediction will not get a final score.
  """


  def __init__(self, parent, *args, **kw):

    self.guiParent = parent
    self.chain     = None
    self.waiting   = False
    self.shiftList = None
    self.spinSystem  = None
    self.preserveTypes = 0
    self.progressBar = None
    self.scoresPopup = None
    self.threshold   = -20.0
    
    BasePopup.__init__(self, parent, title="Assignment : Spin System Typing", **kw)

  def body(self, guiFrame):

    guiFrame.grid_columnconfigure(3, weight=1)

    self.progressBar = TypingEnsemblePopup(self,total=100)
    self.progressBar.close()
    
    row = 0
    label = Label(guiFrame, text=' Chain: ', grid=(row,0))
    tipText = 'Selects which molecular chain the spin system residue types will be predicted for; determines which range of types are available'
    self.chainPulldown = PulldownList(guiFrame, self.changeChain,
                                      grid=(row,1), tipText=tipText)

    tipText = 'Selects which shift list will be used as the source of chemical shift information to make the residue type predictions'
    label = Label(guiFrame, text='Shift List: ', grid=(row,2))
    self.shiftListPulldown = PulldownList(guiFrame, callback=self.setShiftList,
                                          grid=(row,3), tipText=tipText)

    utilButtons = UtilityButtonList(guiFrame, helpUrl=self.help_url)
    utilButtons.grid(row=row, column=4, sticky='w')

    row += 1
    frame = LabelFrame(guiFrame, text='Options', grid=(row,0), gridSpan=(1,5))
    frame.grid_columnconfigure(3, weight=1)

    frow = 0
    label = Label(frame, text='Keep existing types?',
                  grid=(frow,0), sticky='e')
    tipText = 'Whether any existing residue type information should be preserved, when predicting the type of others'
    self.preserveTypesSelect = CheckButton(frame, grid=(frow,1), selected=False, 
                                           callback=self.selectPreserveTypes,
                                           tipText=tipText)
  

    label = Label(frame, text='Assignment threshold: ',
                  grid=(frow,2), sticky='e')
    tipText = 'The lower limit for the predicted residue type to be set with "Assign Types"; needs to be adjusted according to result statistics and amount of shift data'
    self.thresholdEntry = FloatEntry(frame, text=self.threshold,
                                     width=8, grid=(frow,3), tipText=tipText)

    frow += 1
    label = Label(frame, text='Ensemble size: ', grid=(frow,0), sticky='e')
    tipText = 'The number of best scoring residue type mappings, from the Monte Carlo search, to use un the prediction'
    self.ensembleEntry = IntEntry(frame,text=20,width=4,
                                  grid=(frow,1), tipText=tipText)

    label = Label(frame, text='Num Search Steps: ', grid=(frow,2), sticky='e')
    tipText = 'The number of iterative steps that will be used in the Monte Carlo search of best spin system to residue type mappings'
    self.stepsEntry = IntEntry(frame, text=100000, width=8,
                               tipText=tipText, grid=(frow,3))

    frow += 1
    label = Label(frame, text='Isotope shifts to consider:',
                  grid=(frow,0), gridSpan=(1,4))
    
    frow += 1
    self.isotopes = ['1H','13C']
    isos   = ['1H','13C','15N']
    colors = [COLOR_DICT[x] for x in isos] 
    tipText = 'Selects which kinds of resonances, in terms of isotope, the residue type predictions will be made with'
    self.isotopeCheckButtons = PartitionedSelector(frame, labels=isos,
                                                   objects=isos, colors=colors,
                                                   callback=self.toggleIsotope,
                                                   selected=self.isotopes,
                                                   grid=(frow,0),
                                                   gridSpan=(1,4), tipText=tipText)
        
    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    labelFrame = LabelFrame(guiFrame, text='Spin Systems', grid=(row,0), gridSpan=(1,5))
    labelFrame.expandGrid(0,0)
    
    tipTexts = ['The spin system serial number',
                'The residue to which the spin system may currently be assigned',
                'Set whether to include a particular spin system in the type predictions',
                'The spin system to residue type match score for a prediction; higher (less negative) is better',
                'The predicted types of residue that the spin system may be',
                'The chemical shifts in the spin system that will be used in the analysis']
    headingList = ['#','Residue','Use?','Score','Types','Shifts']
    justifyList = ['center','center','center','center','center','left']
    editWidgets      = [None, None, None, None, None, None]
    editGetCallbacks = [None, None, self.toggleInclude, None, None, None]
    editSetCallbacks = [None, None, None, None, None, None]
    self.scrolledMatrix = ScrolledMatrix(labelFrame, headingList=headingList,
                                         justifyList=justifyList,
 					 editSetCallbacks=editSetCallbacks,
                                         editWidgets=editWidgets,
 					 editGetCallbacks=editGetCallbacks,
                                         callback=self.selectCell,
                                         grid=(0,0), tipTexts=tipTexts)

    row += 1
    tipTexts = ['Execute the Monte Carlo search that will make the residue type predictions for the spin systems',
                'Assign the residue type of spin systems with a unique type prediction and prediction score above the stated threshold',
                'Show a residue type prediction for the selected spin system alone; only considers that spin system of shifts, not how all spin systems fit to the sequence',
                'Show a table of peaks that are assigned to the resonances of the selected spin system']
    texts    = ['Run\nTyping','Assign\nTypes',
                'Show Individual\nClassification',
                'Show\nPeaks']
    commands = [self.run, self.assign,
                self.individualScore,
                self.showPeaks]
    bottomButtons = ButtonList(guiFrame, texts=texts, commands=commands,
                               grid=(row,0), gridSpan=(1,5), tipTexts=tipTexts)
    
    self.runButton    = bottomButtons.buttons[0]
    self.assignButton = bottomButtons.buttons[1]
    self.scoreButton  = bottomButtons.buttons[2]
    self.peaksButton  = bottomButtons.buttons[2]
    self.runButton.config(bg='#B0FFB0')
    
    for func in ('__init__','delete'):
      self.registerNotify(self.updateChains, 'ccp.molecule.MolSystem.Chain', func)
      self.registerNotify(self.updateShiftLists, 'ccp.nmr.Nmr.ShiftList', func)
      
    for func in ('__init__','delete','setCcpCode',
                 'setResidue','addResonance', 'setName',
                 'removeResonance','setResonances'):
       self.registerNotify(self.updateSpinSystemsAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)

    self.updateChains()
    self.updateShiftLists()
    self.updateSpinSystems()

  def individualScore(self):
  
    if self.scoresPopup:
      self.scoresPopup.open()
      self.scoresPopup.update(self.spinSystem, self.chain)
    
    else:
      self.scoresPopup = SpinSystemTypeScoresPopup(self.guiParent, spinSystem=self.spinSystem, chain=self.chain)

  def showPeaks(self):
  
    if self.spinSystem:
      peaks = []
      for resonance in self.spinSystem.resonances:
        for contrib in resonance.peakDimContribs:
          peaks.append(contrib.peakDim.peak)

      if len(peaks) > 0:
        self.guiParent.viewPeaks(peaks)
 

  def assign(self):

    for spinSystem in self.scrolledMatrix.objectList:
      if spinSystem.sstSelected == 'Yes':
	if (len(spinSystem.sstTypes) == 1) and (spinSystem.ssScore > self.threshold):
          ccpCode = spinSystem.sstTypes[0]
          
          if spinSystem.residue and (spinSystem.residue.ccpCode != ccpCode):
            assignSpinSystemResidue(spinSystem,residue=None)
          
          if spinSystem.ccpCode != ccpCode:
            assignSpinSystemType(spinSystem,ccpCode,'protein')
         
    self.updateSpinSystems()


  def toggleIsotope(self, isotope):
  
    if isotope in self.isotopes:
      self.isotopes.remove(isotope)
      
    else:
      self.isotopes.append(isotope)
    
    self.updateSpinSystemsAfter()

  def selectPreserveTypes(self, boolean):
  
    self.preserveTypes = boolean

  def toggleInclude(self, *obj):
  
    if self.spinSystem:
      if self.spinSystem.sstSelected == 'Yes':
        self.spinSystem.sstSelected = 'No'
      else:
        self.spinSystem.sstSelected = 'Yes'	
    
      self.updateSpinSystemsAfter()

  def selectCell(self, spinSystem, row, col):
  
    self.spinSystem = spinSystem
    self.updateButtons()


  def updateShiftLists(self, *opt):
   
    shiftLists = getShiftLists(self.nmrProject)
    names = ['Shift List %d' % x.serial for x in shiftLists]
    index = 0
    
    if shiftLists:
      if self.shiftList not in shiftLists:
        self.shiftList = shiftLists[0]
        
      index = shiftLists.index(self.shiftList)	
    
    self.shiftListPulldown.setup(names, shiftLists, index)
  
  def setShiftList(self, shiftList):
  
    self.shiftList = shiftList
    self.updateSpinSystemsAfter()

  def updateButtons(self):
  
    if self.chain and self.scrolledMatrix.objectList:
      self.runButton.enable()
      self.assignButton.enable()
  
    else:
      self.runButton.enable()
      self.assignButton.enable()
      
    if self.spinSystem:
      self.scoreButton.enable()
      self.peaksButton.enable()

    else:
      self.scoreButton.disable()
      self.peaksButton.disable()


  def updateSpinSystemsAfter(self, spinSystem=None):

    if spinSystem:
      spinSystem.sstTypes      = []
      spinSystem.ssScore       = None
      spinSystem.codeScoreDict = {}

    if self.waiting:
      return
      
    else:
      self.waiting = True
      self.after_idle(self.updateSpinSystems)

  def updateSpinSystems(self):

    textMatrix = []
    objectList = []
    if self.project:
      for spinSystem in self.nmrProject.resonanceGroups:
        if not spinSystem.resonances:
          continue
      
        if self.chain:
          if spinSystem.residue and (spinSystem.residue.chain is not self.chain):
            continue
        
          if spinSystem.chains and (self.chain not in spinSystem.chains):
            continue
        
	if hasattr(spinSystem, 'sstSelected'):
	  includeText = spinSystem.sstSelected
      
        else:
	  spinSystem.sstSelected = 'Yes'
	  includeText = 'Yes'

	if not hasattr(spinSystem, 'sstTypes'):
	  spinSystem.sstTypes = []

	if not hasattr(spinSystem, 'ssScore'):
	  spinSystem.ssScore = None
      
        if spinSystem.ssScore:
          scoreText = '%.2f' % spinSystem.ssScore
        else:
          scoreText = None
      
        typesText = ' '.join(spinSystem.sstTypes)
      
        residueText = ''
	
	if spinSystem.residue:
          resCode = getResidueCode(spinSystem.residue)
	  residueText = '%d%s' % (spinSystem.residue.seqCode,resCode)
	
        elif spinSystem.residueProbs:
          resTexts = []
          resSeqs = []
          resCodes = set()
 
          for residueProb in spinSystem.residueProbs:
            if not residueProb.weight:
              continue
 
            residue = residueProb.possibility
            seq = residue.seqCode
            resCode = getResidueCode(residue)
            resText = '%d?%s' % (seq, resCode)

            resTexts.append(resText)
            resSeqs.append('%d?' % seq)
            resCodes.add(resCode)
 
          if len(resCodes) == 1:
            residueText = '/'.join(resSeqs) + resCodes.pop()
          else:
            residueText = '/'.join(resTexts)
          
	elif spinSystem.ccpCode:
	  getResidueCode(spinSystem) 
	
	shifts = []
	if self.shiftList:
	  for resonance in spinSystem.resonances:
            if resonance.isotopeCode in self.isotopes:
	      shift = resonance.findFirstShift(parentList = self.shiftList)
              if shift:
	        shifts.append('%.2f' % shift.value)
	    
	shifts.sort()
	
	shiftsText = ' '.join(shifts)
	
	data = []
	data.append(spinSystem.serial)
	data.append(residueText)
	data.append(includeText)
	data.append(scoreText)
	data.append(typesText)
	data.append(shiftsText)
	
	objectList.append(spinSystem)
        textMatrix.append(data)
        

    self.scrolledMatrix.update(textMatrix=textMatrix, objectList=objectList)
    self.updateButtons()
    self.waiting = False

  def getChains(self):
  
    chains = []
    if self.project:
      for molSystem in self.project.sortedMolSystems():
        for chain in molSystem.sortedChains():
          if chain.molecule.molType in ('protein',None):
            # None moltype may be mixed, including protein component
	    text = '%s:%s' % (molSystem.code, chain.code)
            chains.append( [text, chain] )
	
    return chains


  def changeChain(self, chain):
    
    self.chain = chain
    self.updateSpinSystemsAfter()
    
  def updateChains(self, *chain):
  
    data = self.getChains()
    names = [x[0] for x in  data]
    chains = [x[1] for x in data]
    index = 0
 
    if chains:
      if self.chain not in chains:
       if self.spinSystem:
         if self.spinSystem.residue:
           self.chain = self.spinSystem.residue.chain
         
       else:
         self.chain = chains[0]
      
      index = chains.index(self.chain)
    
    self.chainPulldown.setup(names, chains, index)

    self.updateButtons()
  
  def run(self):
    
    spinSystems = []
    for spinSystem in self.scrolledMatrix.objectList:
      if spinSystem.sstSelected == 'Yes':
        spinSystems.append(spinSystem)
    
    if self.chain and spinSystems:
      if self.progressBar:
        self.progressBar.destroy()
      self.progressBar = TypingEnsemblePopup(self,total=100)
      residues = self.chain.sortedResidues()
      numBest  = self.ensembleEntry.get() or 20
      numSteps = max(100, self.stepsEntry.get() or 100000)
      graph    = self.progressBar.graph
      typeScores, cc0 = getSpinSystemTypes(residues, spinSystems, self.preserveTypes, isotopes=self.isotopes,
                                           shiftList=self.shiftList, numSteps=numSteps, numBest=numBest,
                                           graph=graph, progressBar=self.progressBar)
      threshold = self.thresholdEntry.get()
      
      for ss in typeScores.keys():
        ss.sstTypes = []
        ss.ssScore  = None
        for ccpCode in typeScores[ss].keys():
	  if ccpCode and typeScores[ss][ccpCode] > threshold:
	    if ccpCode not in ss.sstTypes:
	      ss.sstTypes.append(ccpCode)
	
        if len(ss.sstTypes) == 1:
          ss.ssScore = typeScores[ss][ss.sstTypes[0]]
        
      self.updateSpinSystemsAfter()

  def destroy(self):

    for func in ('__init__','delete'):
      self.unregisterNotify(self.updateChains, 'ccp.molecule.MolSystem.Chain', func)
      self.unregisterNotify(self.updateShiftLists, 'ccp.nmr.Nmr.ShiftList', func)
      
    for func in ('__init__','delete','setCcpCode',
                 'setResidue','setName','addResonance',
                 'removeResonance','setResonances'):
       self.unregisterNotify(self.updateSpinSystemsAfter, 'ccp.nmr.Nmr.ResonanceGroup', func)

    if self.scoresPopup:
      self.scoresPopup.destroy()

    BasePopup.destroy(self)

class TypingEnsemblePopup(ProgressBar):

  def __init__(self, parent, ensembleSize=20, *args, **kw):

    self.guiParent = parent
    self.ensembleSize = ensembleSize or 20
    self.labels = []

    ProgressBar.__init__(self, parent, title="Highest Scoring Mappings", text=' Progress', **kw)

  def body(self, guiFrame):

    ProgressBar.body(self, guiFrame)
    guiFrame.expandGrid(2,3)
    
    self.stepLabel = Label(guiFrame, text='Best step:  ', grid=(0,3))

    
    frame = LabelFrame(guiFrame, text='Best mappings', grid=(1,0), gridSpan=(1,4))

    row = 0
    for i in range(self.ensembleSize):
      label = Label(frame, text='', pady=0, font='Courier 10',
                    borderwidth=0, grid=(row,0), sticky='ew')
      self.labels.append(label)
      row +=1

    guiFrame.grid_rowconfigure(2, weight=1)
    self.graph = ScrolledGraph(guiFrame, width=450, height=300,
                               graphType='scatter', title='Typing Scores',
                               xLabel='Spin System', yLabel='Score',
                               grid=(2,0), gridSpan=(1,4))


    self.buttonList = ButtonList(guiFrame, texts=['Close',], commands=[self.done],
                                 grid=(3,0), gridSpan=(1,4))
    self.buttonList.buttons[0].disable() 

  def updateSequences(self, num, data):

    labels = self.labels
    N = len(labels)
    self.stepLabel.set('Best scoring step: %d' % num)

    for i, datum in enumerate(data):
      if i < N:
        labels[i].set(datum)

    self.update_idletasks() 

  def update(self):
  
    if self.progress == self.total:
      width = int(self.cWidth)
      self.canvas.coords(self.bar,self.bw,self.bw,width,self.cHeight)
      self.percent.set( ' %3.1d %%' % 100)
      self.buttonList.buttons[0].enable() 
  
    else:
      ProgressBar.update(self)
   
  def done(self):
  
    self.progress = self.total
    ProgressBar.update(self)
    

