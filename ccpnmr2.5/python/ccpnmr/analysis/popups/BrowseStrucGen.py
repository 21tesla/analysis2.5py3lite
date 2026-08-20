"""
======================COPYRIGHT/LICENSE START==========================

BrowseStrucGen.py: Part of the CcpNmr Analysis program

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

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.Entry import Entry
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showOkCancel
from memops.gui.ScrolledMatrix import ScrolledMatrix

from ccpnmr.analysis.popups.BasePopup import BasePopup

# Removed - not used in file/ 24/4/08 Rasmus Fogh
#from ccpnmr.analysis.core.StructureBasic import makeStructures, makeStructureDictFromCnsPdb, makePdbFromStructure

class BrowseStrucGenPopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.guiParent  = parent
    self.strucGen   = None
    self.violList   = None
    self.constrList = None
    self.waiting       = 0

    BasePopup.__init__(self, parent=parent, title="Structure Generation Runs", **kw)
    
    #self.geometry("+150+150")

  def body(self, guiFrame):
    
    guiFrame.grid_columnconfigure(0, weight=1)

    row = 0
    strucGenFrame = LabelFrame(guiFrame, text='Generation Runs')
    strucGenFrame.grid(row = row, column = 0, columnspan=1, sticky='nsew')
    strucGenFrame.grid_columnconfigure(0, weight=1)
    strucGenFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(row, weight=0)
    
    #self.editDetailsEntry = Entry(self,text='',returnCallback = self.setDetails, width=12)
    #editWidgets      = [None, None, None, None, None, self.editDetailsEntry]
    #editGetCallbacks = [None, None, None, None, None, self.getDetails]
    #editSetCallbacks = [None, None, None, None, None, self.setDetails]
    colHeadings = ['#','Constraint\nLists','Violation\nLists','Structures','Fixed\nAtom Sets','Fixes\nResonance Sets','Chain\nStates','Database\nEntries','Resonance\nClouds']
    editWidgets      = [None, None, None, None, None, None, None, None, None]
    editGetCallbacks = [None, None, None, None, None, None, None, None, None]
    editSetCallbacks = [None, None, None, None, None, None, None, None, None]
    self.structGenMatrix = ScrolledMatrix(strucGenFrame, editSetCallbacks=editSetCallbacks, editGetCallbacks=editGetCallbacks, editWidgets=editWidgets,initialRows=3,initialCols=6, headingList=colHeadings, callback=self.selectStructGenCell, objectList=[], textMatrix=[[],])
    self.structGenMatrix.grid(row = 0, column = 0, sticky='nsew')

    texts    = ['View Structures','Delete']
    commands = [self.viewStructures,self.deleteStrucGen]
    self.structGenButtons = ButtonList(strucGenFrame,texts=texts, expands=True, commands=commands)
    self.structGenButtons.grid(row=1, column=0, sticky='ew')
    
    row += 1
    constrFrame = LabelFrame(guiFrame, text='Constraint Lists')
    constrFrame.grid(row=row, column=0, columnspan=1, sticky='nsew')
    constrFrame.grid_columnconfigure(0, weight=1)
    constrFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(row, weight=1)
    
    colHeadings = ['#','Type','Name','Constraints','Experiments','Details','Unit']
    editWidgets      = [None, None, None, None, None, None, None]
    editGetCallbacks = [None, None, None, None, None, None, None]
    editSetCallbacks = [None, None, None, None, None, None, None]
    self.constrListMatrix = ScrolledMatrix(constrFrame, editSetCallbacks=editSetCallbacks, editGetCallbacks=editGetCallbacks, editWidgets=editWidgets,initialRows=10, headingList=colHeadings, callback=self.selectConstrListCell, objectList=[], textMatrix=[[],])
    self.constrListMatrix.grid(row = 0, column = 0, sticky='nsew')

    texts    = ['View Constraints','Create List','Delete List']
    commands = [self.viewConstraints,self.createConstraints,self.deleteConstraints]
    self.constrListButtons = ButtonList(constrFrame,texts=texts, expands=True, commands=commands)
    self.constrListButtons.grid(row=1, column=0, sticky='ew')
    self.constrListButtons.buttons[1].disable()
    
    row += 1
    violFrame = LabelFrame(guiFrame, text='Violation Lists')
    violFrame.grid(row=row, column=0, columnspan=1, sticky='nsew')
    violFrame.grid_columnconfigure(0, weight=1)
    violFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(row, weight=1)
    
    colHeadings = ['#','Violations','Structures','Details',]
    editWidgets      = [None, None, None, None]
    editGetCallbacks = [None, None, None, None]
    editSetCallbacks = [None, None, None, None]
    self.violListMatrix = ScrolledMatrix(violFrame, editSetCallbacks=editSetCallbacks, editGetCallbacks=editGetCallbacks, editWidgets=editWidgets,initialRows=10, headingList=colHeadings, callback=self.selectViolListCell, objectList=[], textMatrix=[[],])
    self.violListMatrix.grid(row = 0, column = 0, sticky='nsew')

    texts    = ['View Violations','Delete List']
    commands = [self.viewViolations,self.deleteViolations]
    self.violListButtons = ButtonList(violFrame,texts=texts, expands=True, commands=commands)
    self.violListButtons.grid(row=1, column=0, sticky='ew')

    row += 1
    self.bottomButtons = UtilityButtonList(guiFrame, helpUrl=self.help_url)
    self.bottomButtons.grid(row = row, column = 0, columnspan=1, sticky = 'ew')
    self.update()

    for func in ('__init__', 'delete','setName','setDetails','setUnit','setExperiments','addExperiment','removeExperiment'):
      for clazz in ('ccp.nmr.Nmr.ChemShiftConstraintList','ccp.nmr.Nmr.DihedralConstraintList','ccp.nmr.Nmr.DistanceConstraintList',
                    'ccp.nmr.Nmr.HBondConstraintList','ccp.nmr.Nmr.JCouplingConstraintList','ccp.nmr.Nmr.RdcConstraintList'):
        self.registerNotify(self.updateAfter, clazz, func)
    for func in ('__init__', 'delete',):
      for clazz in ('ccp.nmr.Nmr.ChemShiftConstraint','ccp.nmr.Nmr.DihedralConstraint','ccp.nmr.Nmr.DistanceConstraint',
                    'ccp.nmr.Nmr.HBondConstraint','ccp.nmr.Nmr.JCouplingConstraint','ccp.nmr.Nmr.RdcConstraint'):
        self.registerNotify(self.updateAfter, clazz, func)

    for func in ('__init__', 'delete','setChainStates','addChainState','removeChainState',
                 'addEntry','removeEntry','setResStructures','addResStructure','setEntries',
                 'removeResStructure','setStructures','addStructure','removeStructure'):
        self.registerNotify(self.updateAfter, 'ccp.nmr.Nmr.StructureGeneration', func)

    for func in ('__init__', 'delete','setDetails'):   
      for clazz in ('ccp.nmr.Nmr.ViolationList',):
        self.registerNotify(self.updateAfter, clazz, func)

    for func in ('__init__', 'delete',):
      for clazz in ('ccp.nmr.Nmr.Violation',):
        self.registerNotify(self.updateAfter, clazz, func)

  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)
    
  def deleteViolations(self):

    if self.violList and showOkCancel('Confirm','Really delete violation list?', parent=self):
      self.violList.delete()
      self.violList = None
      self.violListButtons.buttons[0].disable()
      self.violListButtons.buttons[1].disable()

  def deleteConstraints(self):

    if self.constrList and showOkCancel('Confirm','Really delete constraint list?', parent=self):
      self.constrList.delete()
      self.constrList = None
      self.constrListButtons.buttons[0].disable()
      self.constrListButtons.buttons[2].disable()

  def createConstraints(self):

    pass

  def viewViolations(self):
    
    self.guiParent.browseViolations()
    popup = self.guiParent.popups['browse_violations']
    popup.structGen  = self.strucGen
    popup.violList   = self.violList
    popup.updateAfter()

  def viewConstraints(self):
  
    self.guiParent.browseConstraints()
    popup = self.guiParent.popups['browse_constraints']
    popup.structGen  = self.strucGen
    popup.constrList = self.constrList
    popup.updateAfter()
      
  def viewStructures(self):
    
    self.guiParent.editStructures()
    popup = self.guiParent.popups['edit_structures']
    popup.structGen = self.strucGen
    popup.updateAfter()
    
  def deleteStrucGen(self):

    if self.strucGen and showOkCancel('Confirm','Really delete structure generation run?', parent=self):
      self.strucGen.delete()
      self.strucGen   = None
      self.constrList = None
      self.violist    = None
      self.structGenButtons.buttons[0].disable()
      self.structGenButtons.buttons[1].disable()

  def selectStructGenCell(self, strucGen, row, col):
  
    self.strucGen = strucGen
    if strucGen:
      self.structGenButtons.buttons[1].enable()
      if len(strucGen.molStructures) > 0:
        self.structGenButtons.buttons[0].enable()
      if self.constrList and (self.constrList.structureGeneration is not self.strucGen):
        self.constrList = None
      self.updateConstrLists(self.strucGen)
      
  def selectConstrListCell(self, constrList, row, col):
  
    self.constrList = constrList
    if constrList:
      self.constrListButtons.buttons[0].enable()
      self.constrListButtons.buttons[2].enable()

  def selectViolListCell(self, violList, row, col):
  
    self.violList = violList
    if violList:
      self.violListButtons.buttons[0].enable()
      self.violListButtons.buttons[1].enable()
      
  def updateAfter(self, object=None):
  
    if self.waiting:
      return
    
    strucGen = None
    name = object.className
    if name == 'StructureGeneration':
      self.waiting = True
      self.after_idle(self.update)
      return
      
    elif name == 'ViolationList':
      strucGen =  object.structureGeneration
    elif name == 'Violation':
      strucGen =  object.violationList.structureGeneration
    elif name[-14:] == 'ConstraintList':
      strucGen =  object.structureGeneration
    elif name[-10:] == 'Constraint':
      strucGen =  object.parentList.structureGeneration

    if (object is None) or (strucGen is self.strucGen) :
      self.waiting = True
      self.after_idle(self.update)
       
  def destroy(self):

    for func in ('__init__', 'delete','setName','setDetails','setUnit','setExperiments','addExperiment','removeExperiment'):
      for clazz in ('ccp.nmr.Nmr.ChemShiftConstraintList','ccp.nmr.Nmr.DihedralConstraintList','ccp.nmr.Nmr.DistanceConstraintList',
                    'ccp.nmr.Nmr.HBondConstraintList','ccp.nmr.Nmr.JCouplingConstraintList','ccp.nmr.Nmr.RdcConstraintList'):
        self.unregisterNotify(self.updateAfter, clazz, func)
    for func in ('__init__', 'delete',):
      for clazz in ('ccp.nmr.Nmr.ChemShiftConstraint','ccp.nmr.Nmr.DihedralConstraint','ccp.nmr.Nmr.DistanceConstraint',
                    'ccp.nmr.Nmr.HBondConstraint','ccp.nmr.Nmr.JCouplingConstraint','ccp.nmr.Nmr.RdcConstraint'):
        self.unregisterNotify(self.updateAfter, clazz, func)

    for func in ('__init__', 'delete','setChainStates','addChainState','removeChainState',
                 'addEntry','removeEntry','setResStructures','addResStructure','setEntries',
                 'removeResStructure','setStructures','addStructure','removeStructure'):
        self.unregisterNotify(self.updateAfter, 'ccp.nmr.Nmr.StructureGeneration', func)

    for func in ('__init__', 'delete','setDetails'):   
      for clazz in ('ccp.nmr.Nmr.ViolationList',):
        self.unregisterNotify(self.updateAfter, clazz, func)

    for func in ('__init__', 'delete',):   
      for clazz in ('ccp.nmr.Nmr.Violation',):
        self.unregisterNotify(self.updateAfter, clazz, func)

    BasePopup.destroy(self)
    
  def updateConstrLists(self, *opt):  
    
    strucGen = self.strucGen
    
    objectList  = []
    textMatrix  = []
    
    if strucGen:
      for constraintList in strucGen.constraintLists:
        objectList.append(constraintList)
    else:
      textMatrix.append([])
      
    for constrList in objectList:
      datum   = []
      expText = ''
      for e in constrList.experiments:
        if expText:
          expText += ' '
        expText += e.name
      datum.append(constrList.serial)
      datum.append(constrList.className[0:-14])
      datum.append(constrList.name)
      datum.append(len(constrList.constraints))
      datum.append(expText)
      datum.append(constrList.details)
      datum.append(constrList.unit)
      textMatrix.append(datum)
    
    if self.constrList:
      self.constrListButtons.buttons[0].enable()
      self.constrListButtons.buttons[2].enable()
    else:
      self.constrListButtons.buttons[0].disable()
      self.constrListButtons.buttons[2].disable()
     
    self.constrListMatrix.update(objectList=objectList, textMatrix=textMatrix)
        
  def updateViolLists(self, *opt):  
    
    strucGen = self.strucGen
    
    objectList  = []
    textMatrix  = []
    
    if strucGen:
      for violationList in strucGen.violationLists:
        objectList.append(violationList)
    else:
      textMatrix.append([])
      
    for violationList in objectList:
      datum   = []
      datum.append(violationList.serial)
      datum.append(len(violationList.violations))
      datum.append(len(violationList.molStructures))
      datum.append(violationList.details)
      textMatrix.append(datum)
    
    if self.violList:
      self.violListButtons.buttons[0].enable()
      self.violListButtons.buttons[1].enable()
    else:
      self.violListButtons.buttons[0].disable()
      self.violListButtons.buttons[1].disable()
    
    self.violListMatrix.update(objectList=objectList, textMatrix=textMatrix)
    
  def update(self):
  
    objectList  = []
    textMatrix  = []
    
    project = self.project

    for strucGen in self.nmrProject.structureGenerations:
      objectList.append(strucGen)
    
    if not objectList:
      textMatrix.append([])
   
    for strucGen in objectList:
      datum   = []
      datum.append(strucGen.serial)
      datum.append(len(strucGen.constraintLists))
      datum.append(len(strucGen.violationLists))
      datum.append(len(strucGen.molStructures))
      datum.append(len(strucGen.fixedAtomSets))
      datum.append(len(strucGen.fixedResonanceSets))
      datum.append(len(strucGen.chainStates))
      datum.append(len(strucGen.entries))
      datum.append(len(strucGen.resStructures))
      textMatrix.append(datum)
          
    if not self.strucGen:
      self.structGenButtons.buttons[0].disable()
      self.structGenButtons.buttons[1].disable()
      
    self.structGenMatrix.update(objectList=objectList, textMatrix=textMatrix)
    self.updateConstrLists()
    self.updateViolLists()
    
    self.waiting = False
 
