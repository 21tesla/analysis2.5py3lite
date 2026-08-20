
"""
======================COPYRIGHT/LICENSE START==========================

EditStructures.py: Part of the CcpNmr Analysis program

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
from math import sqrt

from memops.universal import Geometry

from memops.general import Implementation

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Frame import Frame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showOkCancel, showYesNo, showWarning, showError
from memops.gui.PulldownList import PulldownList
from memops.gui.RadioButtons import RadioButtons
from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccp.util.Validation import getModelValidation, storeEnsembleValidation

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.StructureBasic import getResiduePhiPsi, makeEnsemble, alignStructures, compareEnsembles
from ccpnmr.analysis.core.StructureBasic import getStructureFromFile, makePdbFromStructure
from ccpnmr.analysis.core.StructureBasic import copyModelToEnsemble, makeEmptyEnsembleCopy
from ccpnmr.analysis.core.ValidationBasic import getEnsembleValidationStore, ANALYSIS_RMSD_CONTEXT
from ccpnmr.analysis.core.ValidationBasic import RMSD_KEYWORDS

# TBD one-chain duplexes for DNA

GRAPH_COLORS = (('#0000A0','#A08000'),
                ('#008080','#A00000'),
                ('#800080','#80A000'),
                ('#008000','#A00000'))

GRAPH_COLORS2 = ('#0000A0','#A08000', '#008080','#A00000', '#800080','#80A000', '#008000','#A00000')

COMPARE_OPTIONS = ('All atoms', 'Backbone only')

class EditStructuresPopup(BasePopup):
  """
  **Manage Coordinate Structure Ensembles and Models**
  
  This system is used to list all of the three-dimensional coordinate
  information in the CCPN project. This data is organised as structure ensembles
  that may contain several models, each of which represents a different
  three-dimensional conformation (albeit only slightly so). All structures in
  CCPN are describes as ensembles, even if there is only one conformation
  present. A structure ensemble can be thought of as being equivalent to the
  ATOM records in a PDB file.

  The various tabs list structure ensemble data in tabular form, if a graphical
  representation of the coordinates is required the user can click on [Viewer] at
  the upper right of the popup window, after selecting an ensemble in any of the
  panels. 

  **Ensembles**

  The first tab contains a table that lists the structure ensembles in the
  current CCPN project. All structures are associated with a specific molecular
  system; a group of molecular chains that go together in the same complex or
  NMR sample. Hence, the "Molecular System" pulldown menu allows the user to
  switch between structures that correspond to different molecular groupings,
  although often there will be only one molecular system present.

  Structure ensembles may be imported from PDB and PDB-like files via the
  [Import] function. The associated [Export] will write the selected ensemble to
  disk in a PDB format; albeit only the ATOM records and using IUPAC atom
  nomenclature. For a more formal PDB file export the FormatConverter_ should be
  used.

  The conformational models of a structure ensemble may be split apart into
  separate entities via [Split Models]; making structures with one model each.
  Likewise several structure, or ensembles, can be combined into a larger
  ensemble with [Merge Into Ensemble], assuming that all the structures
  represent the same set of atoms.

  The [Superpose Coords] function will move and rotate the relative relative
  orientations of the models within an ensemble so that they superpose (spatially
  align) in an optimal way. This allows for ensemble analyses of structural
  variation, including the calculation of root mean square deviation (RMSD)
  values, which gives a measure of coordinate spread. Once the models of an
  ensemble are superposed the RMSD values will appear in the relevant tables.

  **Structure Models**

  The second tab lists all of the conformational models that are present within
  a structure ensemble (identified by molecular system code and ensemble number)
  selected in the "Ensemble" pulldown menu. The table lists all of the models
  with various textual annotations (that can be changes by the user) and RMSD
  values of each model from the overall ensemble. It should be noted that RMSD
  values will only be displayed if the models of the ensemble have been
  superposed, which can be achieved in this tab using [Superpose & Calculate
  RMSDs].

  The lower graph represents the heavy atom backbone and all atom RMSD values
  for residues across the whole ensemble. Thus the user is able to see any areas
  of regional variation as far as the molecular sequence is concerned.

  **Residues**

  The "Residues" tab lists per-residue information for a structure ensemble
  selected in the upper left pulldown menu. Additionally the user can choose to
  display only certain kinds of residue by changing the "Ccp Code" pulldown
  menu. The main table lists the residues with any RMSD values (if they are
  calculated during model superposition), backbone dihedral angles and various
  other structural parameters. Many of the tables columns can be added or
  removed by checking the various options above the table.

  The "Validation Parameter" selection at the bottom right of the panel relates
  to structure-derived values that may be used to asses the quality of a
  structure. As standard these will represent various kinds of RMSD value,
  although other validation parameters may appear from time to time. The chosen
  set of parameters may be deleted at any time, causing any corresponding column
  of the table to disappear. The [Display Params] function will use the
  per-residue values to colour a three-dimensional display of the structure,
  thus helping to visualise the location of any structural problems.

  **Coords**

  The last tab lists all of the coordinate records for a selected
  conformational model in a selected structure ensemble. These are equivalent
  to the ATOM records on a MODEL of a PDB file. The table is for display only
  and cannot currently be edited.

  **Model Superposition & RMSD Calculation**

  The way that coordinate superposition is done in Analysis aims to eliminate
  any need to make any residue or atom selections. The superposition from the
  inbuilt method should be very good indeed, and aims to give a solution fairly
  close to global more probabilistic methods. It is however not appropriate for
  aligning more minor regions of structure that may be connected via a flexible
  linker to a major globular part; the superposition will converge on the major
  structured region.

  The underlying basis of the method is singular value decomposition to
  calculate an optimum rotation for two coordinate sets with co-located
  centroids. Each atom is weighted, both for the SVD and centroid calculation.
  There are multiple rounds of superimposition to get the ensemble from pairwise
  model comparison and to better refine the weights.

  Initially the weights come from the atomic masses, but in the later stages of
  the superposition the weights come from the atomic RMSDs calculated in the
  earlier round. Accordingly, the dissimilar parts have proportionately little
  influence on the final ensemble (and thus you don't need to choose the
  'structured' residues). Precisely, the weights are adjusted to be one over the
  mean squared deviation (so for an atom RMSD of 0.1 A its weight is 100, but
  for an RMSD of 3.0 A its weight is 0.111.) and then in the last stage the
  wighting is  exp(-(RMSD/k)^2), which is even more strict (with k=0.8, 0.1 A
  gives weight=0.98, 3.0 A gives weight=0.00000078),

  The superimposition procedure used in Analysis is as follows:

   * Translate coordinates to co-locate (weighted) centroids

   * SVD superimpose all models to first model; rotation
   
   * Calculate mean and closest model to mean

   * Re-calculate weights based upon RMSD

   * SVD superimpose again with new weights to closest to mean

   * Repeat once again using RMSD-based weights rather than atomic weights



  .. _FormatConverter: FormatConverter.html

  """


  def __init__(self, parent, *args, **kw):

    self.guiParent = parent
    self.structGen = 'Any'
    self.molSystem = False  # why False and not None??
    self.molSystem2 = None
    self.structure = None
    self.residue = None
    self.ccpCode = None
    self.model = None
    self.coordAtom = None
    self.compareRmsdsDict = {}
    self.compareRmsdsDict[True] = {}  # backbone only rmsd cache
    self.compareRmsdsDict[False] = {}  # all atoms rmsd cache
    self.waiting = False
    self.waitingCoords = False
   
    BasePopup.__init__(self, parent=parent, title="Structure : Structures", **kw)

  def body(self, guiFrame):
    
    self.geometry("600x600")
    
    guiFrame.expandGrid(0,0)

    tipTexts = ['A table of the structure ensembles, according to the molecular system they represent',
                'A table representing the conformational models in an ensemble and a graph of per-residue RMSDs',
                'A table of residues found in a structure and several per-residue parameters',
                'A table of the coordinate locations of atoms within a conformational model',
                'A table and graph of the comparison of structure ensembles']
    options = ['Ensembles','Structure Models','Residues','Coords', 'Compare Ensembles']
    tabbedFrame = TabbedFrame(guiFrame, options=options, tipTexts=tipTexts,
                              callback=self.changeTab, grid=(0,0))
    structFrame, modelFrame, resFrame, coordsFrame, compareFrame = tabbedFrame.frames
    self.tabbedFrame = tabbedFrame 
    

    # Ensembles
    label = Label(structFrame, text='Molecular System:', grid=(0,0))
    tipText = 'Selects which molecular system (grouping of chains) to show structure ensembles for'
    self.systemPulldown = PulldownList(structFrame, callback=self.changeMolSystem,
                                       grid=(0,1), tipText=tipText)

    structFrame.expandGrid(1,1)
    
    tipTexts = ['The identification number of the structure ensemble for the molecular system it represents',
                'Indicates which chain codes the coordinates of the structure cover',
                'The number of residues the coordinates of the structure cover',
                'The number of conformational models in the structure ensemble',
                'The calculated all-atom root mean square deviation of the whole ensemble, in units of Angstroms',
                'The name of any structure generation run/job that generated the structure ensemble',
                'The number of violation lists that have been calculated for the structure ensemble',
                ]
    colHeadings = ['#', 'Chains', 'Residues', 'Models', 'RMSD',
                   'Generation\nRun', 'Violation\nLists']
    self.structMatrix = ScrolledMatrix(structFrame, initialRows=7, tipTexts=tipTexts,
                                       headingList=colHeadings, callback=self.selectStructure,
                                       multiSelect=True, deleteFunc=self.deleteStructure,
                                       grid=(1,0), gridSpan=(1,2))

    tipTexts = ['Show a table of atom coordinates for a model from the selected ensemble',
                'Import a structure from file into the CCPN project',
                'Export the selected structure ensemble as a pseudo-PDB format file; no header, ATOM coordinate records only with IUPAC naming',
                'Delete the selected structure ensemble from the project; any source PDB file remains',
                'Merge the selected structures into a single milti-model ensemble; only possible of atom records are identical',
                'Split the selected structure ensemble into separate structure records, each containing only one conformational model',
                'Perform a two-stage RMSD-weighted superposition of the conformational models in an ensemble; generates RMSD values']
    texts    = ['Show\nCoords','Import','Export','Delete',
                'Merge Into\nEnsemble','Split\nModels','Superpose\nCoords']
    commands = [self.showCoords, self.importStructure,
                self.exportStructure, self.deleteStructure,
                self.mergeModels, self.splitModels, self.alignEnsemble]
    self.structButtons = ButtonList(structFrame,texts=texts, commands=commands,
                                    grid=(2,0), gridSpan=(1,2), tipTexts=tipTexts)
    
    # Models
    
    modelFrame.expandGrid(1,0)
    modelFrame.grid_rowconfigure(2, minsize=350)
    
    frame = Frame(modelFrame, grid=(0,0))
    frame.expandGrid(None,2)

    label = Label(frame, text='Ensemble:', grid=(0,0))
    tipText = 'Selects which structure ensemble to show the conformational models and RMSD graph for'
    self.modelsEnsemblePulldown = PulldownList(frame, grid=(0,1), tipText=tipText,
                                               callback=self.changeEnsemble)
   
    commands = [self.alignEnsembleM,
                self.deleteModel]
    tipTexts = ['Calculate all-atom and backbone atom coordinate root mean square deviations for the selected ensemble; fills the graph',
                'Delete the currently selected structural model from the ensemble']
    texts = ['Superpose & Calculate RMSDs',
             'Delete Model']
    buttons = ButtonList(frame, grid=(0,3), commands=commands,
                         texts=texts, tipTexts=tipTexts)
    
    self.modelNameEntry = Entry(self, returnCallback=self.setModelName, width=12)
    self.modelDetailsEntry = Entry(self, returnCallback=self.setModelDetails, width=12)
                                
    editWidgets      = [None, self.modelNameEntry,
                        None, self.modelDetailsEntry]
    editGetCallbacks = [None, self.getModelName,
                        None, self.getModelDetails]
    editSetCallbacks = [None, self.setModelName,
                        None, self.setModelDetails]
                        
    tipTexts = ['The number of the conformational model within the structure ensemble',
                'A short name for the conformational model, if needed',
                'The calculated all-atom root mean square deviation of the model coordinates relative to the whole ensemble',
                'A user-provided verbose textual comment about the conformational model']
    headingList = ['Model #','Name','RMSD','Details']
    self.modelMatrix = ScrolledMatrix(modelFrame,headingList=headingList,
                                      editWidgets=editWidgets,
                                      editGetCallbacks=editGetCallbacks,
                                      editSetCallbacks=editSetCallbacks,
                                      callback=self.selectModel,
                                      grid=(1,0), tipTexts=tipTexts)
    
    
    self.rmsdGraph = ScrolledGraph(modelFrame, symbolSize=2, reverseX=False, width=550,
                                   height=220, title='Residue RMSDs',
                                   xLabel='Residue', yLabel='RMSD', grid=(2,0))

    # Residues
    resFrame.expandGrid(2,2)    
    
    frame = Frame(resFrame)
    frame.grid(row=0, column=0, columnspan=3, sticky='ew')
    frame.grid_columnconfigure(4, weight=1)

    label = Label(frame, text='Ensemble:', grid=(0,0))
    tipText = 'Selects which structure ensemble to show residue records for'
    self.resEnsemblePulldown = PulldownList(frame, callback=self.changeEnsemble,
                                            grid=(0,1), tipText=tipText)
    
    label = Label(frame, text='Ccp Code:', grid=(0,2))
    tipText = 'Allows the residue display to be restricted to only a specific kind of residue, if required'
    self.ccpCodePulldown = PulldownList(frame, callback=self.changeCcpCode,
                                        grid=(0,3), tipText=tipText)
    
    frame = Frame(resFrame, grid=(1,0), gridSpan=(1,3), sticky='ew')
    frame.grid_columnconfigure(8, weight=1)
    
    label = Label(frame, text=' Angles', grid=(0,0))
    tipText = 'Whether to show backbone dihedral angles in the table'
    self.angleCheck = CheckButton(frame, selected=True, tipText=tipText, 
                                  callback=self.updateResidues, grid=(0,1))
    
    label = Label(frame, text='RMSDs', grid=(0,2))
    tipText = 'Whether to show per-residue coordinate RMSD values in the table; calculated within the ensemble'
    self.rmsdCheck = CheckButton(frame, selected=True, tipText=tipText,
                                 callback=self.updateResidues, grid=(0,3))
        
    label = Label(frame, text='Rama.\nparams', grid=(0,4))
    tipText = 'Whether to show the Ramachandran plot regional classification using backbone angles'
    self.ramaCheck = CheckButton(frame, selected=True, tipText=tipText,
                                 callback=self.updateResidues, grid=(0,5))
        
    label = Label(frame, text='Validation\nresults', grid=(0,6))
    tipText = 'Whether to show per-residue structure validation results (of various types where calculated) in the table'
    self.validCheck = CheckButton(frame, selected=False, tipText=tipText,
                                  callback=self.updateResidues, grid=(0,7))
        
    label = Label(frame, text='Physical\nparams')
    #label.grid(row=0, column=6, sticky='w')
    self.physCheck = CheckButton(frame, selected=False,
                                 callback=self.updateResidues)
    #self.physCheck.grid(row=0, column=7, sticky='w')

    tipTexts = ['The code of the molecular chain to which the residue belongs',
                'The sequence number of the residue in its chain',
                'The CCPN code for identifying the kind of residue']
    self.residueTipTexts = tipTexts  
    self.resColHeadings = ['Chain','Seq\nCode','Ccp\nCode']
    self.residueMatrix = ScrolledMatrix(resFrame,headingList=self.resColHeadings,
                                        callback=self.selectResidue,
                                        grid=(2,0), gridSpan=(1,3), tipTexts=tipTexts)
 
    tipTexts = ['Superpose the conformational models of the selected structure and calculate per-residue and overall RMSD values',
                'Open a graphical structure display and label/highlight the selected residue',
                'Show the selected set of validation parameters as colors on a graphical structure display',
                'Delete the selected set of structure validation parameters']
    texts    = ['Calculate\nRMSDs',
                'View Residue',
                'Display\nParams',
                'Delete\nParams']
    commands = [self.alignEnsembleR,
                self.viewResidue,
                self.displayStrucParams,
                self.deleteStrucParams]
                
    self.resButtons = ButtonList(resFrame,texts=texts, tipTexts=tipTexts,
                                 commands=commands, grid=(3,0))
                                 
    label = Label(resFrame, text='Validation\nParameter:', grid=(3,1))
    
    tipText = 'Selects a structure validation parameter to consider; on a graphical structure display of for deletion'
    self.strucParamPulldown = PulldownList(resFrame, grid=(3,2), tipText=tipText)      
   
    # Coords
    coordsFrame.expandGrid(1,4)
    
    label = Label(coordsFrame, text='Ensemble', grid=(0,0))
    tipText = 'Selects which structure ensemble to select a conformational model from'
    self.ensemblePulldown = PulldownList(coordsFrame, callback=self.changeEnsemble,
                                         grid=(0,1), tipText=tipText)
    
    label = Label(coordsFrame, text='Model', grid=(0,2))
    tipText = 'Selects which conformational model within the selected ensemble to display coordinates for'
    self.modelPulldown = PulldownList(coordsFrame, callback=self.changeModel,
                                      grid=(0,3), tipText=tipText)

    tipTexts = ['Number of atom within the structure',
                'Name of the atom with coordinates',
                'The CCPN code of the residue the atom is in',
                'The code of the molecular chain the atom is in',
                'The number of the residue in the sequence; may differ from the molecular system numbering',
                'The atom location along the X-axis, in units of Angstroms',
                'The atom location along the Y-axis, in units of Angstroms',
                'The atom location along the Z-axis, in units of Angstroms',
                'The fractional occupancy of this coordinate position by the atom',
                'The X-ray crystallographic B-factor',
                'An alternative location code for coordinate record']
    colHeadings = ['#', 'Atom', 'Residue', 'Chain',
                   'SeqId', 'x', 'y', 'z', 'Occupancy',
                   'B Factor', 'Locn Code']
    self.coordsMatrix = ScrolledMatrix(coordsFrame, headingList=colHeadings,
                                       callback=self.selectCoordCell,
                                       grid=(1,0), gridSpan=(1,5), tipTexts=tipTexts)

    # Compare Ensembles

    compareFrame.expandGrid(1,0)
    compareFrame.grid_rowconfigure(2, minsize=50)
    
    frame = Frame(compareFrame, grid=(0,0))
    frame.expandGrid(None,2)

    label = Label(frame, text='Molecular System:', grid=(0,0))
    tipText = 'Selects which molecular system (grouping of chains) to compare structure ensembles for'
    self.systemPulldown2 = PulldownList(frame, callback=self.changeMolSystem2,
                                       grid=(0,1), tipText=tipText)

    tipTexts = ('Compare all atoms', 'Compare only backbone atoms')
    self.compareRadio = RadioButtons(frame, entries=COMPARE_OPTIONS, select_callback=self.changeComparison,
                                     tipTexts=tipTexts, grid=(0,3))

    commands = [self.compareEnsembles]
    tipTexts = ['Compare the ensembles; fills the table and graph with RMSD values']
    texts = ['Compare Ensembles']
    buttons = ButtonList(frame, grid=(0,4), commands=commands,
                         texts=texts, tipTexts=tipTexts)
    
    tipTexts = ['The ensembles']
    headingList = ['Ensembles']
    self.compareMatrix = ScrolledMatrix(compareFrame, headingList=headingList,
                                        grid=(1,0), tipTexts=tipTexts)
    
    self.compareGraph = ScrolledGraph(compareFrame, symbolSize=2, reverseX=False, width=550,
                                      height=220, title='Residue comparison RMSDs',
                                      xLabel='Residue', yLabel='RMSD', grid=(2,0))
       
    # Main
 
    tipTexts = ['Open a graphical display of the current structure',]
    texts = [ 'Viewer' ]
    commands = [ self.viewStruct ]
    self.bottomButtons = UtilityButtonList(tabbedFrame.sideFrame, helpUrl=self.help_url,
                                           commands=commands, texts=texts, grid=(0,0),
                                           sticky = 'e', tipTexts=tipTexts)
    
    self.updateMolSystems()
    self.updateMolSystems2()
    self.updateAfter() 

    self.administerNotifiers(self.registerNotify)

  def updateRmsdGraph(self):
  
    if self.structure:
      data = (self.structure.molSystem.code,self.structure.ensembleId)
      self.rmsdGraph.title = '%s:%d Residue RMSDs' % data
      self.rmsdGraph.zoom = 1.0

      validStore = getEnsembleValidationStore(self.structure,
                                              ANALYSIS_RMSD_CONTEXT,
                                              RMSD_KEYWORDS)
      
      dataNames = []
      dataColors = []
      dataSets = []
      
      chains = self.structure.sortedCoordChains()
      maxColor = len(GRAPH_COLORS)
      
      for i, chain in enumerate(chains):
        chainData = [[], []]
      
        for residue in chain.sortedResidues():
          findValid = residue.findFirstResidueValidation
      
          for j in (0,1):
            
            validObj = findValid(validationStore=validStore,
                                 context=ANALYSIS_RMSD_CONTEXT,
                                 keyword=RMSD_KEYWORDS[j])
 
            if validObj:
              chainData[j].append( (residue.seqCode, validObj.floatValue) )
        
        
        dataSets.extend(chainData)
        dataColors.extend(GRAPH_COLORS[i % maxColor])
        
        if len(chains) > 1:
          dataNames.append('Chain %s backbone' % chain.code)
          dataNames.append('Chain %s all atom' % chain.code)
        else:
          dataNames.append('Backbone')
          dataNames.append('All atom')
      
      self.rmsdGraph.update(dataSets=dataSets,
                            dataColors=dataColors,
                            dataNames=dataNames)  

  def updateModels(self):
  
    objectList = []
    textMatrix = []

    if self.structure:
      validStore = getEnsembleValidationStore(self.structure,
                                              ANALYSIS_RMSD_CONTEXT,
                                              RMSD_KEYWORDS)
            
      objectList = self.structure.sortedModels()
      for model in objectList:
        rmsdValid = getModelValidation(validStore, model,
                                       ANALYSIS_RMSD_CONTEXT,
                                       RMSD_KEYWORDS[1])

        if rmsdValid:
          rmsd = rmsdValid.floatValue
        else:
          rmsd = None  
        
        datum = [model.serial, model.name, rmsd, model.details]
        textMatrix.append(datum)
      

    self.modelMatrix.update(textMatrix=textMatrix,
                            objectList=objectList)

  def getModelName(self, model):

    if self.model:
      self.modelNameEntry.set(self.model.name)
  
  def setModelName(self, *whatever):
  
    if self.model:
      value = self.modelNameEntry.get().strip() or None
      
      if value and len(value) > 80:
        msg = 'Name too long, truncated at 80 characters.\n'
        msg += 'Consider using the details field instead'
        showWarning('Warning', msg, parent=self)
        value = value[:80]
      
      self.model.setName(value)

  def getModelDetails(self, model):

    if self.model:
      self.modelDetailsEntry.set(self.model.name)
  
  def setModelDetails(self, *whatever):
  
    if self.model:
      value = self.modelDetailsEntry.get().strip() or None
      self.model.setDetails(value)
  
  def selectModel(self, obj, row, col):
  
    self.model = obj
  
  def deleteModel(self):
  
    if self.model:
      ensemble = self.model.structureEnsemble
      data = (self.model.serial, ensemble.molSystem.code, ensemble.ensembleId)
      msg = 'Really delete model %s from ensemble %s:%d?' % data
    
      if showOkCancel('Confirm', msg, parent=self):
        self.model.delete()
        self.model = None
    
        validStore = getEnsembleValidationStore(ensemble,
                                                ANALYSIS_RMSD_CONTEXT,
                                                RMSD_KEYWORDS)

        if validStore:
          self.alignEnsembleM(ensemble)
    
  def displayStrucParams(self):
  
    paramKey = self.strucParamPulldown.getObject()
    
    if self.structure and paramKey:
      context, keyword, validStore = paramKey
      self.parent.viewStructure(self.structure)
      popup = self.parent.popups['view_structure']
      popup.displayResidueParams(validStore, context, keyword)
      
  def deleteStrucParams(self):
  
    paramKey = self.strucParamPulldown.getObject()
    
    if self.structure and paramKey:
      context, keyword, validStore = paramKey
      
      for chain in self.structure.coordChains:
        for residue in chain.residues:
          getValidObj = residue.findFirstResidueValidation    
          residueValidation = getValidObj(context=context, keyword=keyword)
    
          if residueValidation:
            residueValidation.delete()
      
      self.updateResidues()
         
  def changeTab(self, index):

    if index == 0:
      self.updateAfter()
    
    elif index == 1:
      self.updateEnsemblePulldowns()
      self.updateModels()
      self.updateRmsdGraph()
      
    elif index == 2:
      self.updateCcpCodePulldown()
      self.updateEnsemblePulldowns()
      self.updateResidues()
      
    elif index == 3:
      self.updateEnsemblePulldowns()
      self.updateCoordModels()
    
    elif index == 4:
      self.updateEnsemblePulldowns()
      self.updateCompareTable()
      self.updateCompareGraph()

  def selectResidue(self, residue, row, col):
  
    self.residue = residue
    self.updateButtons()


  def updateResidues(self, obj=None):
  
    objectList = []
    textMatrix = []
    colorMatrix = []
    headingList = self.resColHeadings[:]
    tipTexts = self.residueTipTexts[:]
  
    doAngles = self.angleCheck.get()
    doRmsd = self.rmsdCheck.get()
    doRama = self.ramaCheck.get()
    #doPhys = self.physCheck.get()
    doValid = self.validCheck.get()
      
    ccpCode = self.ccpCode

    if doRmsd:
      headings = ['Backbone\nRMSD','All atom\nRMSD',]
      #             u'C\u03B1\nRMSD',u'C\u03B2\nRMSD',
      #            'Hn\nRMSD','O\nRMSD',]
      headingList.extend(headings)
      tipTexts.append('The coordinate root mean square deviation for the residue within the ensemble, considering only backbone atoms')
      tipTexts.append('The coordinate root mean square deviation for the residue within the ensemble, considering all atoms')
    
    if doAngles:
      headingList.extend([u'\u03A6\nAngle',u'\u03A8\nAngle'])
      tipTexts.append('The phi protein backbone dihedral angle of the residue (C-N-CA-C)')
      tipTexts.append('The psi protein backbone dihedral angle of the residue (N-CA-C-N)')
      
    if doRama:
      headingList.extend(['Rama.\nRegion',]) # 'Rama.\nProb'])
      tipTexts.append('The region of a Ramachandran plot protein backbone dihedral angles lie in')
      
    validKeys = []
    if self.structure:
      
      if doValid:
        validKeys = set()
        
        for validStore in self.structure.validationStores:
          validObjs = validStore.findAllValidationResults(className='ResidueValidation')
          
          for validObj in validObjs:
            context = validObj.context
            
            if 'CcpNmrAnalysis' in context:
              continue 
          
            key = (context, validObj.keyword, validStore)
            validKeys.add(key)
            
        validKeys = list(validKeys)
        validKeys.sort()
        headingList.extend(['%s\n%s' % (c,k) for c,k,s in validKeys])
        tipFormat = 'Structure validation results for keyword "%s" in the context of "%s"'
        tipTexts.extend([tipFormat % (k,c) for c,k,s in validKeys])
       
      rmsdStore = getEnsembleValidationStore(self.structure,
                                             ANALYSIS_RMSD_CONTEXT)
      
      nCols = len(headingList)
        
      for chain in self.structure.sortedCoordChains():
        chainCode = chain.code
        
        for residue in chain.sortedResidues():
          sysResidue = residue.residue
          if ccpCode and ccpCode != sysResidue.ccpCode:
            continue

          datum = [chainCode, residue.seqCode, sysResidue.ccpCode]
          colors = [None] * nCols
          
          findValid = residue.findFirstResidueValidation
          
          if doRmsd:
            for keyword in RMSD_KEYWORDS[:2]:
              validObj = findValid(validationStore=rmsdStore,
                                   context=ANALYSIS_RMSD_CONTEXT,
                                   keyword=keyword)
 
              if validObj:
                datum.append(validObj.floatValue)
              else:
                datum.append(None)
          
          phiPsi = None
          if doAngles:
            phiPsi = getResiduePhiPsi(residue)
            datum.extend(phiPsi)
          
          if doRama:
            if not phiPsi:
              phiPsi = getResiduePhiPsi(residue)
          
            prob = None
            category = None
            
            phi, psi = phiPsi
            if (phi is not None) and phi > 0:
              category = u'+\u03A6'
            
            elif psi is not None:
              if -110 < psi < 40:
                category = u'\u03B1'
              else:
                category = u'\u03B2'
                
            datum.extend([category,])# prob])
          
          if doValid:
            for context, keyword, validStore in validKeys:
              validObj = findValid(validationStore=validStore,
                                   context=context,
                                   keyword=keyword)
 
              if validObj.floatValue is not None:
                datum.append(validObj.floatValue)
                
              elif validObj.intValue is not None:
                datum.append(validObj.intValue)
                
              elif validObj.textValue is not None:
                text = validObj.textValue
                datum.append(text)
                col = len(datum) - 1
                if text == 'red':
                  colors[col] = '#FF8080'
                
                elif text == 'orange':
                  colors[col] = '#E0D080'
                
                elif text == 'green':
                  colors[col] = '#B0FFB0'
                
              else:
                datum.append(None)
              
            
          objectList.append(residue)
          textMatrix.append(datum)
          colorMatrix.append(colors)
      
      if doRmsd and rmsdStore.findFirstValidationResult(className='ResidueValidation',
                                                        context=ANALYSIS_RMSD_CONTEXT):
          for keyword in RMSD_KEYWORDS[:2]:
            validKeys.append((ANALYSIS_RMSD_CONTEXT, keyword, rmsdStore))
      
    self.residueMatrix.update(textMatrix=textMatrix,
                              objectList=objectList,
                              headingList=headingList,
                              colorMatrix=colorMatrix,
                              tipTexts=tipTexts)
     
    names = []
    index = 0
    objects = []
    
    if validKeys:
      for validKey in validKeys:
        context, keyword, validStore = validKey
        
        names.append('%s:%s' % (context,keyword))
        objects.append(validKey)
      
      index = min(self.strucParamPulldown.index, len(names))
                              
    self.strucParamPulldown.setup(names, objects, index)
    
    
  def updateEnsemblePulldowns(self):

    texts = []
    index = 0
    ensemble = self.structure
    
    ensembles = []
    empty = []
    for molSystem in self.project.molSystems:
      msCode = molSystem.code
      
      for e1 in molSystem.structureEnsembles:
        if e1.models:
          ensembles.append([msCode, e1.ensembleId, e1])
        else:
          empty.append(e1)
     
    if empty:
      msg = 'One or more completely blank structures are present '
      msg + 'in this project and will be removed.'
      showWarning('Warning', msg, parent=self)
      
      for e1 in empty:
        e1.delete()
    
    ensembles.sort()
    ensembles = [x[2] for x in ensembles]
    
    if ensembles:
      texts = ['%s:%d' % (e.molSystem.code, e.ensembleId) for e in ensembles]
      
      if ensemble not in ensembles:
        ensemble = ensembles[-1]
 
      index = ensembles.index(ensemble)
    
    else:
      ensemble = None
      
    if self.structure is not ensemble:
      self.changeEnsemble(ensemble)
 
    self.resEnsemblePulldown.setup(texts,ensembles,index)
    self.ensemblePulldown.setup(texts,ensembles,index)
    self.modelsEnsemblePulldown.setup(texts,ensembles,index)

    self.updateCcpCodePulldown()

  def updateCcpCodePulldown(self):
  
    texts = ['<All>',]
    ccpCodes = [None, ]
    index = 0
    
    ccpCodes0 = set()
    if self.structure:
      for chain in self.structure.sortedCoordChains():
        for residue in chain.sortedResidues():
          sysResidue = residue.residue
          ccpCodes0.add(sysResidue.ccpCode)
    
    ccpCodes0 = list(ccpCodes0)
    ccpCodes0.sort()
    
    texts.extend(ccpCodes0)
    ccpCodes.extend(ccpCodes0)
    
    if self.ccpCode not in ccpCodes:
      self.ccpCode = None
    
    index = ccpCodes.index(self.ccpCode)
  
    self.ccpCodePulldown.setup(texts, ccpCodes, index)

  def changeCcpCode(self, ccpCode):
  
    if ccpCode is not self.ccpCode:
      self.ccpCode = ccpCode
      self.updateResidues()

  def changeEnsemble(self, ensemble):
      
    if ensemble is not self.structure:

      if self.model and (ensemble is not self.model.structureEnsemble):
        self.model = ensemble.sortedModels()[0]
      else:
        self.model = None

      self.structure = ensemble
      index = self.tabbedFrame.selected
      
      if index == 1:
        self.updateModels()
        self.updateRmsdGraph()
         
      elif index == 2:    
        self.updateResidues()
        
      elif index == 3:    
        self.updateCoordModels()
        self.updateCoords()
     
      self.updateButtons()
        
  def changeModel(self, model):
  
    if model is not self.model:
      self.model = model
      self.updateCoords()

  def updateCoordModels(self):

    texts = []
    index = 0
    model = self.model
    
    models = []
    if self.structure:
      models = self.structure.sortedModels()

      if models:
        texts = ['%d' % m.serial for m in models]
        
        if model not in models:
          model = models[0]
      
        index = models.index(model)
      
      else:
        models = None
           
    else:
      models = None
      
    if self.model is not model:
      self.model = model
      self.updateCoords()  

    self.modelPulldown.setup(texts,models,index)

  def administerNotifiers(self, notifyFunc):

    for func in ('__init__', 'delete', 'setName', 'setDetails'):
       for clazz in ('ccp.molecule.MolStructure.Model',):
        notifyFunc(self.updateModelsAfter, clazz, func)

    for func in ('__init__', 'delete','setX','setY','setZ',
                 'setOccupancy','setBFactor','setaltLocationCode'):
       for clazz in ('ccp.molecule.MolStructure.Coord',):
        notifyFunc(self.updateCoordsAfter, clazz, func)

    for func in ('__init__', 'delete','addMolStructure',
                 'removeMolStructure','setMolStructures','setDetail'):
     for clazz in ('ccp.nmr.Nmr.StructureGeneration',
                   'ccp.nmr.NmrConstraint.ViolationList'):
        notifyFunc(self.updateAfter, clazz, func)
        
    for func in ('__init__', 'delete','setName'):
      notifyFunc(self.updateMolSystemsAfter, 'ccp.molecule.MolSystem.MolSystem', func)

    # __init__ removed due to unloading during structure chain creation
    for func in ('__init__','delete','setAtomNamingSystem','setResNamingSystem'):
      notifyFunc(self.updateAfter, 'ccp.molecule.MolStructure.StructureEnsemble', func)

  def open(self):
  
    self.changeTab(self.tabbedFrame.selected)
    
    BasePopup.open(self)

  def mergeModels(self):

    ensembles = self.structMatrix.currentObjects
    
    if ensembles and self.molSystem:
      models = [x for y in ensembles for x in y.sortedModels()]
      if models:
        makeEnsemble(models, ensemble=models[0].structureEnsemble, 
                     replaceStructures=True)
    
  
  def splitModels(self):
  
    for ensemble in self.structMatrix.currentObjects:
      for model in ensemble.sortedModels():
        newEnsemble = makeEmptyEnsembleCopy(ensemble)
        copyModelToEnsemble(newEnsemble, model, identicalEnsembles=True)
        
      ensemble.delete()
 
  def showCoords(self):
  
    self.updateButtons()

    if self.structure:
      if self.model and (self.model.structureEnsemble is not self.structure):
        self.model = self.structure.findFirstModel()
    
    self.updateCoordModels()
    self.tabbedFrame.select(3)  
  
  def changeMolSystem(self, molSystem):
  
    if molSystem is not self.molSystem:
      self.molSystem = molSystem
      if self.structure and (self.molSystem is not molSystem):
        self.structure = None
    
    self.updateAfter()

  def changeMolSystem2(self, molSystem2):
  
    if molSystem2 is not self.molSystem2:
      self.molSystem2 = molSystem2
      self.updateCompareTable()

  def updateMolSystemsAfter(self, molSystem=None):

    self.after_idle(self.updateMolSystems)
    self.after_idle(self.updateMolSystems2)

  def updateMolSystems(self):
  
    molSystems = [None,] + self.project.sortedMolSystems()
    texts = ['<New>']
    index = 0
    
    if len(molSystems) > 1:
      texts += [ms.code for ms in molSystems[1:]]
      if self.molSystem not in molSystems:
        self.molSystem = molSystems[1]
      
      index = molSystems.index(self.molSystem)
      
    else:
      self.molSystem = None

    self.systemPulldown.setup(texts,molSystems,index)  
    
  def updateMolSystems2(self):
  
    molSystems = self.project.sortedMolSystems()
    index = 0
    
    texts = [ms.code for ms in molSystems]
    if molSystems and self.molSystem2 not in molSystems:
      self.molSystem2 = molSystems[0]
      index = molSystems.index(self.molSystem2)
      self.updateCompareTable()
    else:
      self.molSystem2 = None
      self.updateCompareTable()

    self.systemPulldown2.setup(texts,molSystems,index)  

  def importStructure(self, fileNames=None):
  
    if not fileNames:
      from os.path import dirname
 
      dataRepository = self.project.findFirstRepository(name='userData')
      projPath = dirname(dataRepository.url.path)
 
      fileTypes = [  FileType('PDB', ['*.pdb']), FileType('All', ['*'])]
      fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
                 title = 'Import PDB file', dismiss_text = 'Cancel',
                 selected_file_must_exist = True, multiSelect=True,
                 directory=projPath)

      fileNames = fileSelectPopup.file_select.getFiles()
    
    ensemble = None
    
    for file in fileNames:
      if not file:
        continue
      
      if self.molSystem is None:
        codeNum = 1
        while self.project.findFirstMolSystem(code='MS%d' % codeNum):
          codeNum += 1
      
        self.molSystem = self.project.newMolSystem(code='MS%d' % codeNum,name='Default',
                                                   details='Default made for structures')
       
      try:
        ensemble = getStructureFromFile(self.molSystem, file)
      except Exception, e:
        showError('Structure from file', str(e), parent=self)
    
    if ensemble:
      self.structure = ensemble
        
    self.updateAfter()

  def alignEnsembleM(self, structure=None):
    
    self.alignEnsemble(self.structure)
    self.updateModels()
    self.updateRmsdGraph()

  def alignEnsembleR(self, structure=None):
    
    self.alignEnsemble(self.structure)
    self.updateResidues()
    
  def alignEnsemble(self, structure=None):
  
    if not structure:
      structure = self.structure
  
    if structure:
      models = structure.models
      
      if len(models) == 1: 
        msg = 'Cannot superpose coodinates in a structure with only one model'
        showWarning('Failure', msg, parent=self)
        return

      msg = 'Really superpose %d models of structure ensemble %s? (The change will persist)' % (len(models),structure.ensembleId)
      
      if showOkCancel('Confirm',msg, parent=self):
        structures, error, structureRmsds, atomRmsdDict = alignStructures([structure,])
        
        if structureRmsds:
          rmsd = sqrt(sum([x*x for x in structureRmsds])/len(structureRmsds))
          rmsdStore = getEnsembleValidationStore(structure, ANALYSIS_RMSD_CONTEXT, RMSD_KEYWORDS)
          storeEnsembleValidation(rmsdStore, ANALYSIS_RMSD_CONTEXT, 'all', structure, rmsd)
    
      self.updateAfter()

  def exportStructure(self):
        
    fileTypes = [  FileType('PDB', ['*.pdb']), FileType('All', ['*'])]
    fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
               title = 'Export PDB file', dismiss_text = 'Cancel',
               selected_file_must_exist = False)

    file = fileSelectPopup.getFile() 
    if file :
      makePdbFromStructure(file,self.structure,model=None) # All models
          
          
  def deleteStructure(self, *event):

    structures = self.structMatrix.currentObjects
    msg = 'Really remove %d structure ensembles?' % len(structures)
    
    if structures and showOkCancel('Confirm',msg, parent=self):
    
      for structure in structures[:]:
        if self.model is structure:
          self.model = None
          
        structure.delete()
 
      self.structure = None
      
 
  def selectStructure(self, structure, row, col):
  
    self.changeEnsemble(structure)
    self.updateButtons()


  def selectCoordCell(self, coordAtom, row, col):
  
    self.coordAtom = coordAtom

    self.updateButtons()

      
  def updateAfter(self, object=None):

    if self.tabbedFrame.selected != 0:
      return
  
    if self.waiting:
      return
      
    self.waiting = True
    self.after_idle(self.update)


  def updateCoordsAfter(self, coord):
    
    if self.tabbedFrame.selected != 3:
      return
    
    if coord.model is not self.model:
      return
  
    elif self.waitingCoords:
      return
   
    else:
      self.waitingCoords = True
      self.after_idle(self.updateCoords)

  def updateModelsAfter(self, model):
  
    index = self.tabbedFrame.selected
    
    if index == 0:
      self.updateAfter()
  
    elif index == 1:
      if model.structureEnsemble is self.structure:
        self.after_idle(self.updateModels)
        self.after_idle(self.updateRmsdGraph)

    elif index == 2:
      if model.structureEnsemble is self.structure:
        self.after_idle(self.updateResidues)

    elif index == 3:
      if model.structureEnsemble is self.structure:
        self.updateCoordModels()
      if model is self.model:
        self.after_idle(self.updateCoordsAfter)
        
  
  def destroy(self):
    
    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)
    
  def updateCoords(self, *opt):  
    """
    NBNB Displays seqCode but not seqInsertCode 
    """
    
    model = self.model
    
    objectList  = []
    textMatrix  = []
    
    if model:
    
      coordinates = model.coordinates
      occupancies = model.occupancies
      bFactors = model.bFactors
      
      offset = 0
      for ii,atom in enumerate(model.structureEnsemble.orderedAtoms):
      
        coordResidue = atom.residue
        chainCode = coordResidue.chain.code
        seqCode = coordResidue.seqCode
        
        residue = coordResidue.residue
        if residue:
          ccpCode = residue.ccpCode
        else:
          # Should be unnecessary but does not hurt
          ccpCode = '???'  
        
        datum   = [ii, atom.name, ccpCode,
                   chainCode, seqCode,
                   coordinates[offset], coordinates[offset+1], 
                   coordinates[offset+2],
                   occupancies[ii], bFactors[ii],
                   atom.altLocationCode]
        
        offset += 3
        
        objectList.append(atom)         
        textMatrix.append(datum)
          
    self.coordsMatrix.update(objectList=objectList, textMatrix=textMatrix)
    self.waitingCoords = False
    
  def update(self):
  
    ensembles  = []
    if self.molSystem:
      ensembles = self.molSystem.sortedStructureEnsembles()

    objectList  = []
    textMatrix  = []
    for structure in ensembles:
    
      codes = [c.code for c in structure.coordChains]
      codes.sort()

      numResidues = 0
      for chain in structure.coordChains:
        numResidues += len(chain.residues)
            
      strucGenStr = None
      sg = structure.structureGeneration
      if sg:
        strucGenStr = '%d:%s' % (sg.serial, sg.name) 

      models = structure.models
      violLists = set()
      for model in models:
        violLists.update(model.violationLists)

      validObj = structure.findFirstEnsembleValidation(context=ANALYSIS_RMSD_CONTEXT, keyword='all')
      if validObj:
        rmsd = validObj.floatValue
      else:
        rmsd = None

      datum = [structure.ensembleId,
               ' '.join(codes),
               numResidues,
               len(models),
               rmsd,
               strucGenStr,
               len(violLists)]
      
      textMatrix.append(datum)
      objectList.append(structure)
            
    self.updateButtons()
    self.structMatrix.update(objectList=objectList, textMatrix=textMatrix)

    self.waiting = False
 
  def updateButtons(self):
    
    if len(self.structMatrix.currentObjects) > 1:
      multiple = True
    else:
      multiple = False  
    
    buttons = self.structButtons.buttons
  
    if self.structure:
      for n in (0, 2, 3, 5, 6):
        buttons[n].enable()
    else:
      for n in (0, 2, 3, 5, 6):
        buttons[n].disable()
    
    if self.structure:
      self.bottomButtons.buttons[0].enable()
    else:   
      self.bottomButtons.buttons[0].disable()
    
    if multiple:
      buttons[4].enable()
    else:
      buttons[4].disable()

  def viewResidue(self):
  
    if self.residue:
      molType = self.residue.residue.molType
      atom = self.residue.findFirstAtom(name='CA') or \
             self.residue.findFirstAtom(name='C1') or \
             self.residue.findFirstAtom()

      self.viewStruct()
    
      popup = self.parent.popups.get('view_structure')

      structFrame = popup.structFrame
      structFrame.highlightResidue(atom)
      structFrame.focusOnAtom(atom)
 
  def viewStruct(self):
    
    index = self.tabbedFrame.selected
    
    if index == 3 and self.model:
      self.parent.viewStructure(self.model.structureEnsemble)
    elif self.structure:
      self.parent.viewStructure(self.structure)
      
  def compareEnsembles(self):

    if self.molSystem2:
      compareBackboneOnly = (self.compareRadio.get() == COMPARE_OPTIONS[1])
      rmsdsDict = self.compareRmsdsDict[compareBackboneOnly]
      structures = self.molSystem2.sortedStructureEnsembles()
      for structure1 in structures:
        for structure2 in structures:
          if structure2 is not structure1:
            key = frozenset((structure1, structure2))
            rmsdsDict[key] = compareEnsembles(structure1, structure2, compareBackboneOnly)

      self.updateCompareTable()
      self.updateCompareGraph()

  def changeComparison(self, entry):

    self.updateCompareTable()
    self.updateCompareGraph()

  def updateCompareTable(self):

    textMatrix  = []
    tipTexts = ['The ensembles']
    headingList = ['Ensemble Number']

    if self.molSystem2:
      compareBackboneOnly = (self.compareRadio.get() == COMPARE_OPTIONS[1])
      rmsdsDict = self.compareRmsdsDict[compareBackboneOnly]

      structures = self.molSystem2.sortedStructureEnsembles()
      for structure1 in structures:
        ensembleId = structure1.ensembleId
        tipTexts.append('Ensemble %d RMSDs' % ensembleId)
        headingList.append('RMSD for Ensemble %d' % ensembleId)
        datum = ['%d' % ensembleId]
        for structure2 in structures:
          if structure2 is structure1:
            rmsd = None
          else:
            key = frozenset((structure1, structure2))
            if key in rmsdsDict:
              rmsd, residueRmsdDict = rmsdsDict[key]
            else:
              rmsd = None
          datum.append(rmsd)
        textMatrix.append(datum)
    else:
      structures = []

    self.compareMatrix.update(headingList=headingList, objectList=structures, textMatrix=textMatrix, tipTexts=tipTexts)

  def updateCompareGraph(self):
  
    if self.molSystem2:
      dataNames = []
      dataColors = []
      dataSets = []
      
      self.compareGraph.title = '%s Residue RMSDs' % self.molSystem2.code
      self.compareGraph.zoom = 1.0
      maxColor = len(GRAPH_COLORS2)
      knt = 0

      compareBackboneOnly = (self.compareRadio.get() == COMPARE_OPTIONS[1])
      rmsdsDict = self.compareRmsdsDict[compareBackboneOnly]

      structures = self.molSystem2.sortedStructureEnsembles()
      for i, structure1 in enumerate(structures):
        chains = structure1.sortedCoordChains()
        for structure2 in structures[i+1:]:
          key = frozenset((structure1, structure2))
          if key in rmsdsDict:
            rmsd, residueRmsdDict = rmsdsDict[key]

            for chain in chains:
              chainData = []
      
              for residue in chain.sortedResidues():
                res = residue.residue
                value = residueRmsdDict.get(residue.residue, 0.0)
            
                chainData.append((residue.seqCode, value))
        
              dataSets.append(chainData)
              dataColors.append(GRAPH_COLORS2[knt % maxColor])
              knt += 1
        
              if len(chains) > 1:
                dataNames.append('%d-%d (chain %s)' % (structure1.ensembleId, structure2.ensembleId, chain.code))
              else:
                dataNames.append('%d-%d' % (structure1.ensembleId, structure2.ensembleId))
      
      self.compareGraph.update(dataSets=dataSets,
                               dataColors=dataColors,
                               dataNames=dataNames)  

