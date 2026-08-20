
"""
======================COPYRIGHT/LICENSE START==========================

ViewStructure.py: Part of the CcpNmr Analysis program

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

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.Label import Label
from memops.gui.Frame import Frame
from memops.gui.PulldownList import PulldownList
from memops.gui.MessageReporter import showWarning

from ccp.gui.ViewStructureFrame import ViewStructureFrame

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.StructureBasic import getAtomSetsDistance, getAtomSetsDihedral, alignStructures
from ccpnmr.analysis.core.ExperimentBasic import getThroughSpacePeakLists
from ccpnmr.analysis.core.PrintBasic import getPrintOption, setPrintOption



distanceMethods = {'NOE sum': 'noe',
                   'Min dist':'min',}

class ViewStructurePopup(BasePopup):
  """
  **A Simple Graphical Display for Macromolecule 3D Coordinates**

  CcpNmr Analysis contains a simple 3D structure viewing module which is used to
  display NMR derived information on macromolecular coordinates. For example
  this viewer can be used to display alternative possibilities for NOE
  assignments as dashed lines that connect different parts of a molecular
  structure. The structural model may be moved and rotated by various means
  listed below. Also, specific atoms may be selected and de-selected in the
  display by left clicking.

  The structures that may be displayed with this system are loaded in to the
  CCPN project via the main Structures_ popup window. A structure is chosen for
  display by selecting from the "MolSystem", "Ensemble" and "Model" pulldown
  menus, i.e. at present only one conformational model is displayed at a time.

  The "Peak List" selection is used in combination with the [Show Peaks] button
  at the bottom, which brings up a table listing all of the peaks, within the
  selected peak list, that relate to that atoms chosen (left mouse click) in the
  structural view.

  Much of the NMR-derived information that is presented in the graphical display
  will be controlled via separate popups, for example the [Show On Structure]
  button of the `Assignment Panel`_ or [Show Selected On Structure] in the
  `Restraints and Violations`_ popup. Nonetheless, some data can be added to the
  display via the viewer directly; ensemble RMSDs and other validation parameters
  can be superimposed as coloured spheres of various sizes.

  **View Controls**
  
  To move and rotate the three-dimensional coordinate display the following
  keyboard controls may be used:
  
  * Rotate: Arrow keys
  
  * Zoom: Page Up & Page Down keys

  * Translate: Arrow keys + Control key

  Or alternatively the following mouse controls:
  
  * Rotate: Middle button click & drag
  
  * Zoom: Mouse wheel or middle button click + Shift key & drag up/down

  * Translate: Middle button click & drag + Control key

  Also an options menu appears when the right mouse button is clicked
  and the left mouse button is used to select and de-select atoms in
  the current model view.

  .. _Structures: EditStructuresPopup.html
  .. _`Assignment Panel`: EditAssignmentPopup.html
  .. _`Restraints and Violations`: BrowseConstraintsPopup.html

  """

  def __init__(self, parent, *args, **kw):
   
    self.molSystem   = None
    self.structure   = None
    self.model       = None
    self.distMethod  = 'noe'
    self.waiting     = False
    self.connections = []
    self.selectedAtoms = set([])
    self.selectedResidues = set([])
    self.guiParent = parent
    self.peakList = None
    BasePopup.__init__(self, parent=parent, title="Structure : Structure Viewer", **kw)


  def body(self, guiFrame):
  
  
    guiFrame.grid_columnconfigure(0, weight=1)
 
    row = 0
    
    frame = Frame(guiFrame, grid=(row,0), sticky='ew') 
    frame.grid_columnconfigure(6, weight=1)
    
    label = Label(frame, text='MolSystem:', grid=(0,0))
    tipText = 'Selects which molecular system to select a structure for'
    self.molSystemPulldown = PulldownList(frame, callback=self.setMolSystem,
                                          grid=(0,1), tipText=tipText)
    
    label = Label(frame, text=' Ensemble:', grid=(0,2))
    tipText = 'Selects which structure ensemble to display, for the specified molecular system'
    self.structurePulldown = PulldownList(frame, callback=self.setStructure, 
                                          grid=(0,3), tipText=tipText)
    
    label = Label(frame, text=' Model:', grid=(0,4))
    tipText = 'Selects which conformational model of the selected structure/ensemble to display'
    self.modelPulldown = PulldownList(frame, callback=self.setModel, 
                                      grid=(0,5), tipText=tipText)
  
    label = Label(frame, text='Peak List:', grid=(0,7), sticky='e') 
    tipText = 'When using the "Show Peak" option, sets which peak list is used to display atom connectivities'
    self.peakListPulldown = PulldownList(frame, callback=self.setPeakList, 
                                         grid=(0,8), sticky='e', tipText=tipText)

    label = Label(frame, text=' Dist Method:', grid=(0,9), sticky='e') 
    tipText = 'Where the distances between sets of atoms are displayed, sets whether to use the NOE equivalent (sum r^-6 intensities) or minimum distance'
    self.distMethodPulldown = PulldownList(frame,
                                           callback=self.setDistMethod,
                                           texts=distanceMethods.keys(),
                                           grid=(0,10), sticky='e', tipText=tipText)
    
    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)

    analysisProject = self.analysisProject
    getOption = lambda key, defaultValue: getPrintOption(analysisProject, key, defaultValue)
    setOption = lambda key, value: setPrintOption(analysisProject, key, value)
    self.structFrame = ViewStructureFrame(guiFrame, project=self.project,
                                          radiiScale=0.0, bondWidth=1,
                                          atomCallback=self.selectAtom,
                                          getPrintOption=getOption,
                                          setPrintOption=setOption,
                                          grid=(row,0))
   
    row += 1
    
    frame = Frame(guiFrame)
    frame.grid(row=row, column=0, sticky='ew') 
    frame.grid_columnconfigure(2, weight=1)
 
    tipTexts = ['Remove all highlights and connections from the structure display',
                'For an ensemble, calculate the per-atom coordinate root mean square deviations and adjust atom size and colours accordingly',
                'Display the selected structural parameters on the structure, adjusting atom labels, size and colours accordingly']
    texts    = ['Reset', 'RMSDs', 'Display Params:',]
    commands = [self.clearConnections,
                self.displayAtomRmsds,
                self.displayStrucParams,]
    
    self.paramButtons = ButtonList(frame, texts=texts, commands=commands,
                                   grid=(0,0), tipTexts=tipTexts)
    
    tipText = 'Selects which structural parameters, from those calculated, to display on the structure'
    self.strucParamPulldown = PulldownList(frame, grid=(0,1), tipText=tipText)   
    
    tipTexts = ['In the stated peak list, display peak assignment connectivities between the highlighted atoms (left click to select atoms in the display)',]
    texts    = ['Show Peaks']
    commands = [self.showPeaks,]
    self.bottomButtons = UtilityButtonList(frame, texts=texts, commands=commands,
                                           helpUrl=self.help_url, grid=(0,3),
                                           tipTexts=tipTexts)

    self.updateMolSystems()
    self.updatePeakLists()
    self.administerNotifiers(self.registerNotify)
    
  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete', 'setName'):
      notifyFunc(self.updateMolSystems, 'ccp.molecule.MolSystem.MolSystem', func)

    for func in ('__init__', 'delete', 'setName'):
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.DataSource', func)
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.Experiment', func)
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.PeakList', func)

    for func in ('delete', '__init__'):
      notifyFunc(self.updateStructures, 'ccp.molecule.MolStructure.StructureEnsemble', func)

    for func in ('delete', '__init__'):
      notifyFunc(self.updateModels, 'ccp.molecule.MolStructure.Model', func)
  
  def updatePeakLists(self, obj=None):
  
    index = 0
    names = []
    peakLists = getThroughSpacePeakLists(self.project)
    
    if peakLists:
      if self.peakList not in peakLists:
        self.peakList = peakLists[0]
      
      index = peakLists.index(self.peakList)
      names = ['%s:%s:%d' % (pl.dataSource.experiment.name, \
               pl.dataSource.name, pl.serial) for pl in peakLists ]
      
      peakLists.append(None)
      names.append('<All Avail>')
    
    else:
      self.peakList = None
    
    self.peakListPulldown.setup(names, peakLists, index)
  
  def setPeakList(self, peakList):
  
    if self.peakList is not peakList:
      self.peakList = peakList
  
  def showPeaks(self):

    if not self.selectedAtoms:
      msg = 'No atoms selected'
      showWarning('Warning', msg, parent=self)
      return
    
    if self.peakList:
      peakLists = [self.peakList,]
    else:
      peakLists = set(getThroughSpacePeakLists(self.project))
    
    resonances = []
    for coordAtom in self.selectedAtoms:
      atom = coordAtom.atom
      atomSet = atom.atomSet
      
      if not atomSet:
        continue
        
      for resonanceSet in atomSet.resonanceSets:
        for resonance in resonanceSet.resonances:
          resonances.append(resonance)
    
    peaks = set([])
    for resonance in resonances:
      for contrib in resonance.peakDimContribs:
        peak = contrib.peakDim.peak
        if peak.peakList in peakLists:
          peaks.add(peak)
   
    self.connections = []
    self.structFrame.clearConnections()
    self.structFrame.clearHighlights()
   
    if peaks:
      for peak in peaks:
        self.showPeakConnection(peak)
      self.guiParent.viewPeaks(list(peaks))
      
    self.updateAfter()

  def displayAtomParamsList(self, atomParams, size=3.0):
    
    self.connections = []
    self.structFrame.clearConnections()
  
    configAtom = self.structFrame.highlightAtom

    if self.structure and atomParams:
      atomDict = {}
      
      for chain in self.structure.coordChains:
        for residue in chain.residues:
          for atom in residue.atoms:
            atomDict[atom.atom] = atom
    
      atomParams.sort()
      minVal = atomParams[0][0]
      maxVal = atomParams[-1][0] 
      scale = maxVal-minVal
    
      if not scale:
        return
      
      for value, atoms in atomParams:
        frac = (value-minVal)/scale
        rgb = [frac, frac, 1-frac]
        label = None

        for atom in atoms:
          coordAtom = atomDict.get(atom)
          if not coordAtom:
            continue
            
          if not label:
            residue = atom.residue
            label = '%d%s%s %.3f' % (residue.seqCode, residue.ccpCode,
                                     atom.name, value) 
          else:
            label = atom.name

          configAtom(coordAtom, rgb, frac*size, label=label)        
    

  
  def displayAtomRmsds(self, scale=5.0):

    self.connections = []
    self.structFrame.clearConnections()
  
    configAtom = self.structFrame.highlightAtom
  
    if self.structure:
      nModels = len(self.structure.models)
      if nModels < 2:
        msg = 'Cannot calculate atom RMSDs:\n'
        if nModels == 1:
          msg += 'Only one model present'
        else:
          msg += 'No models'
          
        showWarning('Failure', msg)
        return
    
      #structures, error, structureRmsds, atomRmsdDict
      data = alignStructures([self.structure,])
      atomRmsdDict = data[3]
    
      for atom, rmsd in atomRmsdDict.items():
      
        frac = min(scale, rmsd)/scale
      
        if atom.name == 'CA':
          residue = atom.residue
          label = '%d%s' % (residue.seqCode, residue.residue.ccpCode) 
        else:
          label = ''
      
        rgb = (frac, frac, 1-frac)
      
        configAtom(atom, rgb, rmsd/scale, label=label)        
    
      self.updateParams()
    
  def displayStrucParams(self):
  
    paramKey = self.strucParamPulldown.getObject()
    
    if self.structure and paramKey:
      context, keyword, validStore = paramKey
      self.displayResidueParams(validStore, context, keyword)
      
  def displayResidueParams(self, validStore, context, keyword):
  
    self.connections = []
    self.structFrame.clearConnections()
    self.structFrame.clearHighlights()
  
    configResidue = self.structFrame.customResidueStyle
  
    if self.structure and (validStore.structureEnsemble is self.structure):
      scores = []
      residues = []
      unscoredResidues = []
      
      if (context=='CING') and (keyword=='ROGscore'):
        for chain in self.structure.coordChains:
          for residue in chain.residues:
            validObj = residue.findFirstResidueValidation(context=context, keyword=keyword)
 
            if validObj:
              value = validObj.textValue
              if value:
                scores.append(value)
                residues.append(residue)
            else:
              unscoredResidues.append(residue)
      else:
        for chain in self.structure.coordChains:
          for residue in chain.residues:
            validObj = residue.findFirstResidueValidation(context=context, keyword=keyword)
 
            if validObj:
              value = validObj.floatValue
              if value is not None:
                scores.append(value)
                residues.append(residue)
            else:
              unscoredResidues.append(residue)
      
      if scores:
        if (context=='CING') and (keyword=='ROGscore'):
          for i, residue in enumerate(residues):
            score = scores[i]
            if not score:
              continue
            
            if score == 'red':
              rgb = (0.9, 0.0, 0.0)
              size = 1.0
            elif score == 'orange':
              rgb = (0.9, 0.7, 0.0)
              size = 0.5
            else:
              rgb = (0.0, 0.5, 0.0)
              size = 0.1
           
            configResidue(residue, label='', color=rgb, size=size)
      
        else:
          upper = max(scores)
          lower = min([s for s in scores if s is not None])
          delta = upper-lower
 
          if delta:
            for i, residue in enumerate(residues):
              score = scores[i]
              frac = 0.9 * (score-lower)/delta
              if context == 'RPF':
                frac = 1-frac
                
              rgb = (frac, frac, 1-frac)
              size = frac
              configResidue(residue, label='%.2f' % score, color=rgb, size=size)
            
        for residue in unscoredResidues:
          configResidue(residue, label='', color=(0.5, 0.5, 0.5), size=0.1)
              
  def selectAtom(self, atom, hilightSize=0.4):
  
    if atom in self.selectedAtoms:
      self.selectedAtoms.remove(atom)
      self.structFrame.clearAtomHighlight(atom, atomSize=hilightSize)
      
    else:
      self.selectedAtoms.add(atom)
      self.structFrame.highlightResidue(atom, color=(0.0,1.0,0.0), atomSize=hilightSize)
    
    selectedResidues = set([a.residue for a in self.selectedAtoms])
    
    for residue in self.selectedResidues:
      if residue not in selectedResidues:
        self.structFrame.clearResidueHighlight(residue)
    
    self.selectedResidues = selectedResidues
    
    #c = atom.findFirstCoord()
    #self.structFrame.clearHighlights()  
    self.structFrame.drawStructure()

  def clearConnections(self):
    
    for atom in self.selectedAtoms:
      self.structFrame.clearAtomHighlight(atom, atomSize=0.4)
    
    for residue in self.selectedResidues:
      self.structFrame.clearResidueHighlight(residue)
    
    self.selectedAtoms = set([])
    self.selectedResidues = set([])
    self.connections = []
    self.structFrame.clearConnections()
    self.structFrame.clearHighlights()
    self.updateAfter()
  
  def showPeakConnection(self, peak):
  
    peakContribs = peak.peakContribs
    dimAtomSets = []
    dimIsotopes = []
    
    for peakDim in peak.sortedPeakDims():
      isotope = None
      
      for contrib in peakDim.peakDimContribs:
        resonance = contrib.resonance
        resonanceSet = resonance.resonanceSet
        
        if resonanceSet:
          atomSets = resonanceSet.atomSets
          isotope = resonance.isotopeCode
          dimAtomSets.append( (peakDim, contrib.peakContribs, atomSets, isotope) )
      
      dimIsotopes.append(isotope)
      
    M = len(dimAtomSets)
    atomSetsPairs = set()
    
    if dimIsotopes.count('1H') > 1:
      hOnly = True
    else:
      hOnly = False  
    
    for i in range(M-1):
      peakDimI, peakContribsI, atomSetsI, isotopeI = dimAtomSets[i]
    
      if hOnly and (isotopeI != '1H'):
        continue
    
      for j in range(i+1,M):
        peakDimJ, peakContribsJ, atomSetsJ, isotopeJ = dimAtomSets[j]
        
        if peakDimJ is peakDimI:
          continue
    
        if isotopeI != isotopeJ:
          continue
    
        if hOnly and (isotopeJ != '1H'):
          continue
        
        if atomSetsI == atomSetsJ:
          continue
    
        if peakContribs:
          for peakContrib in peakContribsJ:
            if peakContrib in peakContribsI:
              atomSetsPairs.add( frozenset([atomSetsI, atomSetsJ]) )
              break
        
        else:
          atomSetsPairs.add( frozenset([atomSetsI, atomSetsJ]) )

    if len(atomSetsPairs):
      for atomSetsI, atomSetsJ in atomSetsPairs:
        value = getAtomSetsDistance(atomSetsI, atomSetsJ,
                                    self.structure,
                                    self.model,
                                    method=self.distMethod)

        color = self.getNoeColor(value)
        self.showAtomSetsConnection(atomSetsI, atomSetsJ,
                                    value, color=color)

    # If the AtomSets could not be paired just highlight the respective atoms.
    else:
      for i in range(M):
        peakDimI, peakContribsI, atomSetsI, isotopeI = dimAtomSets[i]
        self.highlightAtoms(atomSetsI)
  
  def showResonancesConnection(self, resonance1,resonance2):

    if resonance1.resonanceSet and resonance2.resonanceSet:
    
      atomSets1 = resonance1.resonanceSet.atomSets
      atomSets2 = resonance2.resonanceSet.atomSets
      value = getAtomSetsDistance(atomSets1,
                                  atomSets2,
                                  self.structure,
                                  self.model,
                                  method=self.distMethod)
                                                                    
      color = self.getNoeColor(value)
      self.showAtomSetsConnection(atomSets1, atomSets2, value, color=color)

  def getNoeColor(self, value):
  
    color = (0.0,1.0,0.0)
    if value > 8.0:
      color = (1.0,0.0,0.0)
    elif value > 5.5:
      color = (1.0,1.0,0.0)
    
    return color

  def showResonancesDihedral(self, resonances):
  
    atomSets = []
    for resonance in resonances:
      if not resonance.resonanceSet:
        return
      atomSets.append(list(resonance.resonanceSet.atomSets))
    
    angle = getAtomSetsDihedral(atomSets, self.structure)
    color = (0.0,1.0,1.0) # Cyan - so it doesn't look like NOEs 

    self.showAtomSetsConnection( atomSets[0], atomSets[1], color=color )
    self.showAtomSetsConnection( atomSets[1], atomSets[2], color=color )
    self.showAtomSetsConnection( atomSets[2], atomSets[3], color=color )
    self.showAtomSetsConnection( atomSets[3], atomSets[0], value=angle, color=color )

  def showConstraintConnections(self, constraint):
  
    if not hasattr(constraint, 'items'):
      print "Structure display not supported for constraint type %s" % constraint.className
    
    if constraint.className == 'DihedralConstraint':
      self.showResonancesDihedral(constraint.resonances)
    
    else: 
      for item in constraint.items:
        resonances = list(item.resonances)
        self.showResonancesConnection(resonances[0],resonances[1])
 
  
  def showAtomSetsConnection(self, atomSets1, atomSets2, value=None, color=(1.0,1.0,1.0)):
  
  
    atomSet1 = list(atomSets1)[0]
    atomSet2 = list(atomSets2)[0]
 
    molSystem = atomSet1.findFirstAtom().residue.chain.molSystem
    if self.molSystem is not molSystem:
      self.molSystem = molSystem
      self.updateMolSystems()
      # will find a structure if there's a valid one
  
    if self.structure:
      residue1 = atomSet1.findFirstAtom().residue
      residue2 = atomSet2.findFirstAtom().residue
      chain1   = residue1.chain
      chain2   = residue2.chain
      
      coordChain1 = self.structure.findFirstCoordChain(code=chain1.code)
      coordChain2 = self.structure.findFirstCoordChain(code=chain2.code)
      
      if not (coordChain1 and coordChain2):
        print "No coord chain found"
        return
      
      coordRes1 = coordChain1.findFirstResidue(seqId=residue1.seqId)
      coordRes2 = coordChain2.findFirstResidue(seqId=residue2.seqId)
  
      if not (coordRes1 and coordRes2):
        print "No coord res found"
        return
      
      coordAtoms1 = []
      for atomSet in atomSets1:
        for atom in atomSet.atoms:
          coordAtom = coordRes1.findFirstAtom(name=atom.name)
          if coordAtom:
            coordAtoms1.append(coordAtom)

      coordAtoms2 = []
      for atomSet in atomSets2:
        for atom in atomSet.atoms:
          coordAtom = coordRes2.findFirstAtom(name=atom.name)
          if coordAtom:
            coordAtoms2.append(coordAtom)
     
      for atom1 in coordAtoms1:
        self.structFrame.highlightAtom(atom1)

      for atom2 in coordAtoms2:
        self.structFrame.highlightAtom(atom2)

      if coordAtoms1 and coordAtoms2:
        if value is None:
          cBond = self.structFrame.drawConnection(coordAtoms1,coordAtoms2,color=color)
 
        else:
          cBond = self.structFrame.drawConnection(coordAtoms1,coordAtoms2,color=color,label='%.3f' % value)
       
        self.connections.append([atomSets1, atomSets2, cBond, color, value])

    else:
      print "Display attempted for atoms without structure"
  
    self.updateAfter()
    
  def highlightAtoms(self, atomSetsIn):
    atomSets = list(atomSetsIn)
    
    molSystem = atomSets[0].findFirstAtom().residue.chain.molSystem
    if self.molSystem is not molSystem:
      self.molSystem = molSystem
      self.updateMolSystems()
      # will find a structure if there's a valid one
  
    if self.structure:
      residue = atomSets[0].findFirstAtom().residue
      chain   = residue.chain
      
      coordChain = self.structure.findFirstCoordChain(code=chain.code)
      
      if not (coordChain):
        print "No coord chain found"
        return
      
      coordRes = coordChain.findFirstResidue(seqId=residue.seqId)
  
      if not (coordRes):
        print "No coord res found"
        return
      
      coordAtoms = []
      for atomSet in atomSets:
        for atom in atomSet.atoms:
          coordAtom = coordRes.findFirstAtom(name=atom.name)
          if coordAtom:
            coordAtoms.append(coordAtom)

      for atom in coordAtoms:
        self.structFrame.highlightAtom(atom)
        
    else:
      print "Display attempted for atoms without structure"
  
    self.updateAfter()
  
  def setDistMethod(self, name):
  
    distMethod = distanceMethods.get(name, 'noe')
    if distMethod != self.distMethod:
      self.distMethod = distMethod
      
      for i, data in enumerate(self.connections):
        atomSets1, atomSets2, cBond, oldColor, oldVal = data
        value = getAtomSetsDistance(atomSets1,
                                    atomSets2,
                                    self.structure,
                                    self.model,
                                    method=self.distMethod)
                                    
        color = self.getNoeColor(value)
        cBond.setColor(color)
        cBond.setAnnotation('%.3f' % value)
        
        self.connections[i][3] = color
        self.connections[i][4] = value

      self.updateAfter()

  def updateModels(self, model=None):
  
    self.updateParams()
  
    if model and (model.structureEnsemble is not self.structure):
      return
  
    names = []
    models = []
    index = 0
    model = self.model
    
    if self.structure:
      models = self.structure.sortedModels()
      
    if models:
      if model not in models:
        model = models[0]
    
      names = ['%d' % m.serial for m in models]
      index = models.index(model)
    
    else:
      model = None

    self.modelPulldown.setup(names, models, index)
    
    if self.model is not model:
      self.setModel(model)
 
  def setModel(self, model):
  
    if model is not self.model:
      self.model = model
      
      oldConn = self.connections[:]
      drawConn = self.showAtomSetsConnection
      
      self.connections = []
      self.structFrame.clearConnections()
      self.update()
      
      for atomSets1, atomSets2, cBond, color, value in oldConn:
        value = getAtomSetsDistance(atomSets1,
                                    atomSets2,
                                    self.structure,
                                    self.model,
                                    method=self.distMethod)
                                    
        drawConn(atomSets1, atomSets2, value, color)
      
  
  def setStructure(self, structure):

    if structure is not self.structure:
      self.structure = structure
      self.model = None
      self.updateModels()


  def getStructures(self):

    structures = []
    if self.molSystem:
      for structure in self.molSystem.sortedStructureEnsembles():
        structures.append(structure)

    return structures

  def updateStructures(self, structure=None):
 
    if structure:
      if self.molSystem and (structure.molSystem is not self.molSystem):
        return
       
      self.molSystem = structure.molSystem
      self.updateMolSystems()
 
    names = []
    index = 0
    structures = self.getStructures()
    
    if structures:
      if self.structure in structures:
        structure = self.structure
      else:
        structure = structures[0]
 
      names = [str(x.ensembleId) for x in structures]
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
        if molSystem.structureEnsembles:
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
 
 
  def update(self, structure=None):
    
    # Main display
     
    if structure:
      self.structure = structure
      self.molSystem = structure.molSystem
      self.updateModels()

    if self.model:
      if self.structFrame.model is not self.model:
        self.structFrame.update(self.model)
        self.structFrame.highlightBackbone()
              
      self.structFrame.drawStructure()
      self.updateStructures()
      self.updateMolSystems()
  
    self.updateParams()
  
    self.waiting = False
  
  
  def updateParams(self):
    # Parameter pulldown
    
    structure = self.structure
      
    names = []
    index = 0
    validKeys = []
    
    if structure:
      validKeys = set()

      for validStore in structure.validationStores:
        validObjs = validStore.findAllValidationResults(className='ResidueValidation')
        
        for validObj in validObjs:
          context = validObj.context
          key = (validObj.context, validObj.keyword, validStore)
          validKeys.add(key)
      
      validKeys = list(validKeys)
      validKeys.sort()
      names = ['%s:%s' % (c,k) for c,k,s in validKeys]
    
    if names:
      index = min(self.strucParamPulldown.index, len(names)-1)
 
    self.strucParamPulldown.setup(names, validKeys, index)
  
    self.waiting = False

  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)
  
  
  
  
