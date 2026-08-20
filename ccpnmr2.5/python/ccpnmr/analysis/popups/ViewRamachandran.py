
"""
======================COPYRIGHT/LICENSE START==========================

ViewRamachandran.py: Part of the CcpNmr Analysis program

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

from ccpnmr.analysis.popups.BasePopup      import BasePopup
from ccpnmr.analysis.core.StructureBasic import getResiduePhiPsi
from ccpnmr.analysis.core.MoleculeBasic  import getResidueCode
from ccpnmr.analysis.core.Util import getHueSortedColorSchemes

from ccp.gui.ViewRamachandranFrame   import ViewRamachandranFrame

from memops.gui.ButtonList      import ButtonList, UtilityButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.Color           import grey, hexToHsb, standardColors
from memops.gui.Frame           import Frame
from memops.gui.Label           import Label
from memops.gui.PulldownList    import PulldownList
from memops.gui.ScrolledMatrix  import ScrolledMatrix

def testRamachandranPopup(argServer):
  
  popup = ViewRamachandranPopup(argServer.parent)
  popup.open()


LABEL_MODES = ('None','Disallowed','All')
DEFAULT_SCHEME = 'rgb'
REGION_TEXT = 'Region Statistics: Core;%5.1f%%, Allowed;%5.1f%%,  Disallowed;%5.1f%%'

class ViewRamachandranPopup(BasePopup):
  """
  **Display Protein Backbone Phi & Psi Angles**
   
  This graphical display allow the user to display phi and psi protein backbone
  dihedral angles on a Ramachandran plot that indicates the likelihood (database
  abundance) of those angles. This can be used as a form of structure quality
  control to detect residues in a calculated three dimensional structure that
  are distorted away from regular protein-like conformations. Although a few
  atypical angles may truly occur in a given protein structure the presence of
  many unusual angles and abnormal distributions of angles over a structure
  ensemble with many models may indicate a poor quality structure.

  With this popup window the user selects a molecular system and then a
  structure ensemble that relates to that system. Th e user can choose to display
  dihedral angle information for all models (conformations) in the structure ensemble
  or just one model, by changing the "Model" pulldown menu. The user can control which 
  residues are considered by selecting either a particular type of residue via
  "Ccp Code" or a specific residue in the sequence. Other options control how the
  phi & psi angle information is presented on screen. The "Labels" options
  control which angle spots on the Ramachandran chart have residue sequence number and
  type displayed; by default only the "disallowed" angles in unusual (white) regions of the
  plot are labelled. The "Spot Size" dictates how large the angle markers are and
  the "Colours" is a colour scheme to differentiate the different models
  from within the selected ensemble.
  
  The [Previous Residue] and [Next Residue] buttons allow the user to quickly
  scan though all of the residues in the structure, to check for
  unusual/disallowed backbone angles or unlikely distributions of within the
  ensemble. In the Ramachandran matrix plot the wite areas represent phi & psi
  angle combinations that are very unusual and thus usually indicative of a poor
  structure at that point. The grey and red areas represent the more common, and
  thus more expected, angle combinations. The more red the colour of a square
  the greater the likelihood of the phi & psi angles. It should be noted that
  when displaying angle points for only a single residue or residue type, the
  Ramachandran chart in the background is at  residue-specific version, i.e. it
  shows angle likelihoods for only that one kind of residue. This is particularly
  important for Pro and Gly residues that have notably different distributions
  (and hence angle expectations) to other residues.
 """

  def __init__(self, parent, *args, **kw):

    self.molSystem   = None
    self.structure   = None
    self.model       = None
    self.waiting     = False
    self.colorScheme = None
    self.ccpCode     = None # Which Rama background to use
    self.residue     = None  
    self.labelMode   = LABEL_MODES[1]

    BasePopup.__init__(self, parent=parent, title='Chart : Ramachandran Plot', **kw)

  def open(self):
  
    BasePopup.open(self)
    self.updateAfter()


  def body(self, guiFrame):

    self.geometry('700x700')

    self.update_idletasks()

    guiFrame.grid_columnconfigure(0, weight=1)
 
    row = 0
    
    frame = Frame(guiFrame, relief='raised', bd=1, grid=(row,0), sticky='ew')
    frame.grid_columnconfigure(8, weight=1)
    
    label = Label(frame, text='MolSystem:', grid=(0,0))
    tipText = 'Selects which molecular system to display data for; from which a structure is selected'
    self.molSystemPulldown = PulldownList(frame, callback=self.setMolSystem,
                                          grid=(0,1), tipText=tipText)
    
    label = Label(frame, text=' Structure:', grid=(0,2))
    tipText = 'Selects which structure ensemble to display phi/psi backbone angle data for'
    self.structurePulldown = PulldownList(frame, callback=self.setStructure,
                                          grid=(0,3), tipText=tipText)
    
    label = Label(frame, text='    Model:', grid=(0,4))
    tipText = 'Selects which conformational model(s), from the structure ensemble, to display data for'
    self.modelPulldown = PulldownList(frame, callback=self.setModel,
                                      grid=(0,5), tipText=tipText)
        

    label = Label(frame, text=' Labels:', grid=(0,6))
    tipText = 'Sets which phi/psi points carry residue labels, e.g. according to whether values are "disallowed" (very uncommon)'
    self.labelPulldown = PulldownList(frame, texts=LABEL_MODES, index=1,
                                      callback=self.updateLabels, grid=(0,7),
                                      sticky='e', tipText=tipText)

    utilButtons = UtilityButtonList(frame, helpUrl=self.help_url, grid=(0,9))

    
    label = Label(frame, text='Ccp Code:', grid=(1,0))
    tipText = 'Allows the phi/psi points to be restricted to only those from a particular residue type'
    self.ccpCodePulldown = PulldownList(frame, callback=self.setCcpCode,
                                        grid=(1,1), tipText=tipText)
    
    label = Label(frame, text=' Residue:', grid=(1,2))
    tipText = 'Allows the display of phi/psi points from only a single residue'
    self.residuePulldown = PulldownList(frame, callback=self.setResidue,
                                        grid=(1,3), tipText=tipText)
    
    label = Label(frame, text=' Spot Size:', grid=(1,4))
    sizes = [1,2,3,4,5,6,7,8,9,10]
    texts = [str(s) for s in sizes]
    tipText = 'Sets how large to display the circles that indicate the phi/psi points'
    self.spotSizePulldown = PulldownList(frame, texts=texts, objects=sizes,
                                         callback=self.updateAfter, index=1,
                                         grid=(1,5), tipText=tipText)
       
    label = Label(frame, text=' Colours:', grid=(1,6)) 
    tipText = 'Selects which colour scheme to use to distinguish phi/psi points from different conformational models'
    self.schemePulldown = PulldownList(frame, callback=self.selectColorScheme,
                                       grid=(1,7), tipText=tipText)
 
    row +=1
    tipText = 'The percentage of residues found in the different regions of the Ramachandran plot, according to PROCHECK cetegories' 
    self.regionLabel = Label(guiFrame, text=REGION_TEXT % (0.0,0.0,0.0),
                             grid=(row, 0), tipText=tipText)

    row +=1
    self.colorRow = row
    self.colorFrame = Frame(guiFrame, sticky='ew')
    self.colorFrame.grid_columnconfigure(0, weight=1)
    self.colorLabels = []

    tipText = 'The colors of the ensembles in the plot'
    label = Label(self.colorFrame, text='Colors: ', grid=(0,0), sticky='w', tipText=tipText)

    row +=1

    guiFrame.grid_rowconfigure(row, weight=1)
    self.plot = ViewRamachandranFrame(guiFrame, relief='sunken',
                                      bgColor=self.cget('bg'))
    self.plot.grid(row=row, column=0, sticky='nsew')

    row +=1
    tipTexts = ['Show phi/spi spots for the previous residue in the sequence; when there is no current residue, starts fron the first in sequence there is data for',
                'Show phi/spi spots for the next residue in the sequence; when there is no current residue, starts fron the first in sequence there is data for']
    texts    = ['Previous Residue','Next Residue']
    commands = [self.prevResidue,self.nextResidue]
    self.bottomButtons = ButtonList(guiFrame, commands=commands,
                                    texts=texts, grid=(row,0),
                                    tipTexts=tipTexts)

    self.updateColorSchemes()
    self.updateMolSystems()
    self.updateAfter()
    
    self.administerNotifiers(self.registerNotify)
     
  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete',):
      notifyFunc(self.updateStructuresAfter, 'ccp.molecule.MolStructure.StructureEnsemble', func)

    for func in ('__init__', 'delete','setColors'):
      notifyFunc(self.updateColorSchemes, 'ccpnmr.AnalysisProfile.ColorScheme', func)

  def stepResidue(self, step=1):
  
    structure = self.structure
    if structure:
      if structure == 'All':
        structure = self.getStructures()[0]
      if not self.residue:
        chain = structure.sortedCoordChains()[0]
        residue = chain.sortedResidues()[1]
 
      else:
        residue = self.residue
        chain =  self.residue.chain
        residues = chain.sortedResidues()[1:]
 
        index = residues.index(residue)
        index = (index+step) % len(residues)
        residue = residues[index]
    
      self.residue = residue
      self.updateResidues()
      self.updateAfter()

  def prevResidue(self):
  
    self.stepResidue(-1)
  
  def nextResidue(self):
  
    self.stepResidue(1)

  def setResidue(self,residue):
  
    if residue is not self.residue:
      if residue:
        self.ccpCode = None
        self.ccpCodePulldown.setIndex(0)
    
      self.residue = residue
      self.updateAfter()

  def setCcpCode(self, ccpCode):
  
    if ccpCode != self.ccpCode:
      if ccpCode:
        self.residue = None
        self.residuePulldown.setIndex(0)
        
      self.ccpCode = ccpCode
      self.updateAfter()

  def updateLabels(self, selection):

    if selection is not self.labelMode:
      self.labelMode = selection
      self.updateAfter()

  def updateStructuresAfter(self, structure):

    if structure.molSystem is self.molSystem:
      self.updateStructures()

  def selectColorScheme(self, scheme):
  
    if scheme is not self.colorScheme:
      self.colorScheme = scheme
      
      if self.structure:
        self.updateAfter()
    

  def updateColorSchemes(self, scheme=None):
  
    schemes = getHueSortedColorSchemes(self.analysisProfile)
    names = [s.name for s in schemes]
    colors = [list(s.colors) for s in schemes]
    index = 0
    
    scheme = self.colorScheme
    
    if schemes:
      if scheme not in schemes:
        scheme = self.analysisProfile.findFirstColorScheme(name=DEFAULT_SCHEME)
        if not scheme:
          scheme = schemes[0]
        
      index = schemes.index(scheme)  
    
    else:
      scheme = None
      
    if scheme is not self.colorScheme:
      self.colorScheme = scheme
      
      if self.structure:
        self.updateAfter()
      
    self.schemePulldown.setup(names, schemes, index, colors=colors)

  def updateResidues(self):
  
    resNames = ['<All>',]
    ccpNames = ['<All>',]
    ccpCodes = [None,]
    residues = [None,]
    resCats  = [None,]
    indexR = 0
    indexC = 0
   
    structure = self.structure
    if structure:
      if structure == 'All':
        structure = self.getStructures()[0]
      models = structure.sortedModels()
      models.append(None)
      
      ccpCodes0 = set()
      for chain in structure.sortedCoordChains():
        chainCode = chain.code
      
        for residue in chain.sortedResidues()[1:]:
          sysResidue = residue.residue
          
          seqCode = residue.seqCode
          ccpCode = sysResidue.ccpCode
          ccpCodes0.add(ccpCode)
          resName = '%d %s' % (seqCode, ccpCode)
          n = 10*int(seqCode/10)
          resCat = '%s %d-%d' % (chainCode,n,n+10)
          
          resNames.append(resName)
          residues.append(residue)
          resCats.append(resCat)
      
      ccpCodes0 = list(ccpCodes0)
      ccpCodes0.sort()
      
      ccpCodes += ccpCodes0
      ccpNames += ccpCodes0
    
    doUpdate = False
    if self.residue not in residues:
      if self.residue and len(residues) > 1:
        chain = residues[1].chain
        self.residue = chain.findFirstResidue(residue=self.residue.residue)
      else: 
        self.residue = None
      doUpdate = True
    
    indexR = residues.index(self.residue)
      
    if self.ccpCode not in ccpCodes:
      self.ccpCode = None
      doUpdate = True
    
    indexC = ccpCodes.index(self.ccpCode)    
      
    self.residuePulldown.setup(resNames, residues, indexR, categories=resCats)
    self.ccpCodePulldown.setup(ccpNames, ccpCodes, indexC)
    
    if doUpdate:
      self.updateAfter() 
       
  def updateModels(self, model=None):
  
    structure = self.structure
    # looks like model is always None
    # possibly need to worry about below otherwise if structure = 'All'
    if model and (model.structureEnsemble is not structure):
      return
    
    self.updateResidues()
    
    models = []
    names = []
    index = 0
    model = self.model
    
    if structure:
      if structure != 'All':
        models = self.structure.sortedModels()
      models.append(None)
     
    if models:
      if model not in models:
        model = models[0]
    
      names = ['%d' % m.serial for m in models[:-1]]
      names.append('<All>')
      index = models.index(model)
    
    else:
      model = None

    self.modelPulldown.setup(names, models, index)
    
    if self.model is not model:
      self.model = model
      self.updateAfter()   
 
  def setModel(self, model):
  
    if model is not self.model:
      self.model = model
      self.updateAfter()
      
  
  def setStructure(self, structure):

    if structure is not self.structure:
      self.structure = structure
      self.model = None
      self.updateModels()
      self.updateAfter()

  def getStructures(self, molSystem=None):

    molSystem = molSystem or self.molSystem

    structures = []
    if molSystem:
      for structure in molSystem.sortedStructureEnsembles():
        for chain in structure.coordChains:
          for residue in chain.residues:
            if residue.residue.molResidue.molType == 'protein':
              structures.append(structure)
              break
              
          else:
            continue
            
          break  
 
    return structures

  def updateStructures(self, structure=None):
 
    if structure and (structure.molSystem is not self.molSystem):
      return    
 
    names = []
    index = 0
    structures = self.getStructures()
    
    if structures:
      names = [str(x.ensembleId) for x in structures]
      names.append('<All>')
      structures.append('All')
      if self.structure in structures:
        structure = self.structure
      else:
        structure = structures[0]
 
      index = structures.index(structure)
 
    else:
      structure = None

    self.structurePulldown.setup(names, structures, index)

    if structure is not self.structure:
      self.structure = structure
      self.model = None
      self.updateModels()

  def setMolSystem(self, molSystem):

    if molSystem is not self.molSystem:
      self.molSystem = molSystem
      self.structure = None
      self.model = None
      self.updateStructures()


  def getMolSystems(self):
 
    molSystems = []
    if self.project:
      for molSystem in self.project.sortedMolSystems():
        if self.getStructures(molSystem):
          molSystems.append(molSystem)

    return molSystems


  def updateMolSystems(self, *object):
 
    names = []
    index = 0
    
    molSystems = self.getMolSystems()
    if molSystems:
      if self.molSystem in molSystems:
        molSystem = self.molSystem
      else:
        molSystem = molSystems[0]
 
      names = [x.code for x in molSystems]
      index = molSystems.index(molSystem)
 
    else:
      molSystem = None
 
    self.molSystemPulldown.setup(names, molSystems, index)

    if molSystem is not self.molSystem:
      self.molSystem = molSystem
      self.structure = None
      self.model = None
      self.updateStructures()
 
 
  def updateAfter(self, object=None):

    if self.waiting:
      return
 
    self.waiting = True
    self.after_idle(self.update)

  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

  def update(self):
  
    self.plot.setAminoAcid(self.ccpCode)
    self.updatePhiPsi()
    self.waiting = False
    
  def updatePhiPsi(self):
    # labels='outliers', ''
    
    if not self.structure:
      return
    
    model = self.model
    ccpCode = self.ccpCode
    residueSel = self.residue
    labelMode = self.labelMode
    getValue = self.plot.getIntensityValue
    
    ccpCode = self.ccpCode
    if self.residue:
      ccpCode = self.residue.residue.ccpCode
    self.plot.setAminoAcid(ccpCode)
    
    phiPsiAccept = []
    plotObjects  = []
    resLabels    = []
    colors = []
    
    if self.colorScheme:
      scheme = list(self.colorScheme.colors)
    else:
      scheme = ['#800000','#008000','#000080']
      
    nCols = len(scheme)

    nCore = 0
    nAllowed = 0
    nDisallowed = 0
    
    colorFrame = self.colorFrame
    colorLabels = self.colorLabels
    structure = self.structure
    if structure == 'All':
      structures = self.getStructures()
      nstructures = len(structures)
      ncolorLabels = len(colorLabels)
      n = nstructures - ncolorLabels
      if n > 0:
        for i in range(n):
          label = Label(colorFrame, grid=(0,i+ncolorLabels+1))
          colorLabels.append(label)
      for i, structure0 in enumerate(structures):
        text = 'Structure %d' % structure0.ensembleId
        label = colorLabels[i]
        label.set(text)
        color = scheme[i % nCols]
        label.config(bg=color)
        if color == '#000000':
          label.config(fg='#FFFFFF')
        else:
          label.config(fg='#000000')
      colorFrame.grid(row=self.colorRow, column=0)
    else:
      structures = [structure]
      colorFrame.grid_forget()

    for structure0 in structures:
      if model:
        models = [model,]
      else:
        models = list(structure0.models)
        
      nModels = len(models)
    
      for chain in structure0.coordChains:
        for residue in chain.residues:
          sysResidue = residue.residue
          sysCode = getResidueCode(sysResidue)
          resLabel = '%d%s' % (sysResidue.seqCode,sysCode)
 
          if sysResidue.molResidue.molType != 'protein':
            continue

          if residue and residueSel and (residue.residue is not residueSel.residue):
            continue
 
          if ccpCode and (sysCode != ccpCode):
            continue
       
          for model0 in models:
            phi, psi = getResiduePhiPsi(residue, model=model0)
          
            if None in (phi,psi):
              continue
          
            value = getValue(phi,psi)
          
            if nModels == 1:
              resLabels.append(resLabel)
            else:
              resLabels.append( '%s:%d' % (resLabel, model0.serial) )
          
            doLabel = False
          
            if value < 6.564e-5:
              if labelMode == LABEL_MODES[1]:
                doLabel = True
              nDisallowed += 1
            elif value < 0.000821:
              nAllowed += 1
            else:
              nCore += 1

            if labelMode == LABEL_MODES[0]:
              doLabel = False
            
            elif labelMode == LABEL_MODES[2]:
              doLabel = True
              
            if structure == 'All':
              ind = structures.index(structure0)
            else:
              ind = model0.serial - 1
            colors.append(scheme[ind % nCols])
            plotObjects.append((residue, model0))
            phiPsiAccept.append((phi,psi,doLabel))
              
    spotSize = self.spotSizePulldown.getObject()
    
    nRes = 0.01*float(nDisallowed+nAllowed+nCore)
    if nRes:
      self.regionLabel.set(REGION_TEXT % (nCore/nRes,nAllowed/nRes,nDisallowed/nRes))
    else:
      self.regionLabel.set(REGION_TEXT % (0.0, 0.0, 0.0))
      
    self.plot.cirRadius = spotSize
    self.plot.updateObjects(phiPsiAccept, plotObjects, resLabels, colors)
