
"""
======================COPYRIGHT/LICENSE START==========================

SetupStructureCalcFrame.py: Part of the CcpNmr Analysis program

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

import os

from memops.gui.ButtonList import ButtonList
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Frame import Frame
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.MessageReporter import showWarning
from memops.gui.DataEntry import askString

from ccpnmr.analysis.core.AssignmentBasic import getShiftLists
from ccpnmr.analysis.core.ExperimentBasic import getThroughSpacePeakLists

from ccpnmr.analysis.popups.BasePopup import BasePopup

KEY_GENERAL = 'general'
KEY_CCPN = 'ccpn'
KEY_MOLECULAR_SYSTEM = 'mol_system'
KEY_NOES = 'noe_data'
KEY_DIHEDRALS = 'dihedral_data'
KEY_DISTANCES_AMBIG = 'ambig_distance_data'
KEY_DISTANCES_UNAMBIG = 'unambig_distance_data'
KEY_JCOUPLINGS = 'jcoupling_data'
KEY_ID = 'id'

class testPopup(BasePopup):

  def ___init__(self, parent, project):
  
    self.project = project
    BasePopup.__init__(self, parent)
    
  def body(self, guiFrame):
  
    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)

    frame = SetupStructureCalcFrame(guiFrame, self)
    frame.grid(row=0, column=0, sticky='nsew')
    

def testSetupStructureCalcFrame(argServer):

  popup = testPopup(argServer.parent)
  popup.open()


# MolSystem
# Starting structure
# Constraint set
# ConstraintLists
# PeakList:ShiftList

def getObjectKeyString(object, delimiter='|'):
  """Descrn: Make an object identifier string for a CCPN data model object
     Inputs: CCPN data model object
     Output: String
  """

  keys = object.getExpandedKey()
  
  for i in range(len(keys)):
    key = keys[i]
    
    keyType = type(key)
    if keyType is type([]):
      keys[i] = delimiter.join([str(k) for k in key])
    elif keyType is not type(''):
      keys[i] = str(key)
  
  return delimiter.join(keys)
  
class SetupStructureCalcFrame(Frame):

  def __init__(self, parent, application, *args, **kw):

    self.parent = parent
    self.molSystem = None
    self.constraintSet = None
    self.application = application
    self.project = application.project
    self.waiting = False
    self.tableObj = None

    Frame.__init__(self, parent=parent, **kw)
    self.grid_columnconfigure(3, weight=1)
 
    label = Label(self, text='Molecular system: ')
    label.grid(row=0, column=0, sticky='nw')
    self.molSysPulldown = PulldownMenu(self, self.changeMolSystem,
                                       tipText='Select which molecular system to use in the calculation')
    self.molSysPulldown.grid(row=0, column=1, sticky='nw')

    label = Label(self,text='Constraint Set: ')
    label.grid(row=0, column=2, sticky='nw')
    self.constraintSetPulldown = PulldownMenu(self, self.changeConstraintSet,
                                              tipText='Select which set of restraints to use in the calculation')
    self.constraintSetPulldown.grid(row=0, column=3, sticky='nw')

    self.shiftListPulldown = PulldownMenu(parent, callback=self.setShiftList,
                                          do_initial_callback=False,
                                          tipText='select the shift list to use in the calculation')
    
    tipTexts=['The type of ccp object used in the calculation',
              'The identifier for the input data',
              'A toggle whether to use the data',
              'Whether to filter violated restraints',
              'Whether to keep CCPN assignments (otherwise new restraint assignments will be used)',
              'Whether the restraints are overly anbiguous and need refining',
              'The shift list to use for shift-peak matching']
    headingList = ['Input Class','Id','Use?','Filter\nViolated',
                   'Keep\nAssignments?','Ambiguous\nProtocol','Shift List']
    editWidgets      = [None, None, None, None, None, None, self.shiftListPulldown]
    editGetCallbacks = [None, None, self.toggleUse, self.toggleFilterViol,
                        self.toggleKeepAssign, self.toggleAmbig, self.getShiftList]
    editSetCallbacks = [None, None, None, None, None, None, self.setShiftList]
                        
    self.grid_rowconfigure(1, weight=1)
    self.scrolledMatrix = ScrolledMatrix(self, headingList=headingList,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks, 
                                         editWidgets=editWidgets,
                                         multiSelect=True,
                                         tipTexts=tipTexts,
                                         callback=self.selectCell)

    self.scrolledMatrix.grid(row=1, column=0, columnspan=4, sticky='nsew')
    self.scrolledMatrix.doEditMarkExtraRules = self.doEditMarkExtraRules
    
    #texts    = ['','','']
                
    #commands = [None,None,None]
                
    #self.bottomButtons = ButtonList(self, texts=texts, expands=True, commands=commands)
    #self.bottomButtons.grid(row=2, column=0, columnspan=4, sticky='ew')
    
    self.updateMolSystems()
    self.updateConstraintSets()

    notify = self.application.registerNotify
    for func in ('__init__','delete'):
      notify(self.updateAfter, 'ccp.nmr.Nmr.PeakList', func)
      notify(self.updateConstraintSets, 'ccp.nmr.NmrConstraint.NmrConstraintStore', func)
      notify(self.updateMolSystems, 'ccp.molecule.MolSystem.MolSystem', func)
      notify(self.updateMolSystems, 'ccp.molecule.MolSystem.Chain', func)
 
    for func in ('__init__','delete','setName'):
      notify(self.updateAfter, 'ccp.nmr.Nmr.Experiment', func)
      notify(self.updateAfter, 'ccp.nmr.Nmr.DataSource', func)


  def doEditMarkExtraRules(self, obj, row, col):
  
    className = obj.className
    if className == 'PeakList':
      if col == 5:
        return False
    
    elif className == 'DistanceConstraintList':
      if col in (4,6):
        return False
        
    else:
      if col in (4,5,6):
        return False
                 
    return True
  
  def toggleUse(self, obj):

    bool = not obj.setupStructureCalcFrameDict['use']
    obj.setupStructureCalcFrameDict['use'] = bool

    self.updateAfter()
  
  def toggleFilterViol(self, obj):

    bool = not obj.setupStructureCalcFrameDict['filter']
    obj.setupStructureCalcFrameDict['filter'] = bool

    self.updateAfter()

  def toggleKeepAssign(self, obj):
  
    bool = not obj.setupStructureCalcFrameDict['preserve']
    obj.setupStructureCalcFrameDict['preserve'] = bool

    self.updateAfter()
  
  def toggleAmbig(self, obj):
  
    bool = not obj.setupStructureCalcFrameDict['ambig']
    obj.setupStructureCalcFrameDict['ambig'] = bool

    self.updateAfter()
  
  def selectCell(self, obj, row, col):
  
    self.tableObj = obj
  
  def getShiftList(self, peakList):

    names = []
    index = -1
    shiftLists = getShiftLists(self.project.currentNmrProject)
    
    if shiftLists:
      shiftList = peakList.dataSource.experiment.shiftList
      
      if shiftList not in shiftLists:
        shiftList = shiftLists[0]
      
      names = [sl.name or str(sl.serial) for sl in shiftLists]
      index = shiftLists.index(shiftList)
        
    print "getShiftList", names, index

    self.shiftListPulldown.setup(names,index)

  def setShiftList(self, index, name=None):
  
    print "setShiftList", index, name
  
    if not name:
      index = self.shiftListPulldown.getSelectedIndex()
      
    shiftList = getShiftLists(self.project.currentNmrProject)[index]
    
    self.tableObj.setupStructureCalcFrameDict['shiftList'] = shiftList

    self.updateAfter()
      
  def updateMolSystems(self, *opt):
  
    names = []
    index = -1
  
    molSystems = self.getMolSystems()
    if molSystems:
      names = [ms.code for ms in molSystems]
    
      if self.molSystem not in molSystems:
        self.molSystem = molSystems[0]
 
      index = molSystems.index(self.molSystem)
    
    self.molSysPulldown.setup(names, index)
    
  def getMolSystems(self):
  
    molSystems = []
    
    if self.project:
      for molSystem in self.project.molSystems:
        if molSystem.chains:
          molSystems.append(molSystem)
           
    return molSystems 

  def getTableObjects(self):
  
    objs = []
    
    if self.project:
      objs = getThroughSpacePeakLists(self.project)
    
      if self.constraintSet:
        objs.extend(list(self.constraintSet.constraintLists))
    
    return objs

  def updateConstraintSets(self, obj=None):
  
    names = []
    index = -1
    
    if self.project and self.project.currentNmrProject:
      constraintSets = self.project.currentNmrProject.sortedNmrConstraintStores()
     
      if constraintSets:
        if self.constraintSet not in constraintSets:
          self.constraintSet = constraintSets[0]
 
        names = [str(cs.serial) for cs in constraintSets]
        index = constraintSets.index(self.constraintSet)

    names.append('<None>')
 
    self.constraintSetPulldown.setup(names, index)
    
  def changeMolSystem(self, i, name):
  
    self.molSystem = self.project.findFirstMolSystem(code=name)
    self.updateAfter()

  def changeConstraintSet(self, i, name):
 
    if name == '<None>':
      self.constraintSet = None
    else:
      self.constraintSet = self.project.currentNmrProject.sortedNmrConstraintStores()[i]
    
    self.updateAfter()
  
  def update(self):

    textMatrix = []
    objectList = []
    #headingList = ['Input Class','Use?','Filter\nViolated',
    #               'Keep\nAssignments?','Ambigous\nProtocol','Shift List']
    
    for obj in self.getTableObjects():
      
      if not hasattr(obj,'setupStructureCalcFrameDict'):
        obj.setupStructureCalcFrameDict = dict = {}
        dict['use'] = True
        dict['filter'] = True
        dict['preserve'] = False
        dict['ambig'] = False
        dict['shiftList'] = None
        new = True
      else:
        dict = obj.setupStructureCalcFrameDict
        new = False
      
      className = obj.className
      
      if className == 'PeakList':
        spectrum   = obj.dataSource
        experiment = spectrum.experiment
        shiftList  = dict['shiftList'] or experiment.shiftList
        
        if shiftList:
          shiftListName = shiftList.name
        else:
          shiftListName = '<None>'
        
        if new and (self.molSystem not in experiment.molSystems):
          dict['use'] = False
      
        ident = '%s:%s:%d' % (experiment.name,spectrum.name,obj.serial)
        used = dict['use'] and 'Yes' or 'No'
        filtered = dict['filter'] and 'Yes' or 'No'
        preserved = dict['preserve'] and 'Yes' or 'No'
        ambig = None
      
      elif className == 'DistanceConstraintList':
        shiftListName = None
        
        if new:
          dict['filter'] = True
          dict['ambig'] = False
          
        ident = '%s:%d' % (obj.name or '', obj.serial)
        used = dict['use'] and 'Yes' or 'No'
        filtered = dict['filter'] and 'Yes' or 'No'
        preserved = None
        ambig = dict['ambig'] and 'Yes' or 'No'
      
      else:
        shiftListName = None

        if new:
          dict['filter'] = False
          dict['ambig'] = False
          
        ident = '%s:%d' % (obj.name or '', obj.serial)
        used = dict['use'] and 'Yes' or 'No'
        filtered = None
        preserved = None
        ambig = None
      
      
      datum = [className,ident,used,filtered,preserved,ambig,shiftListName]
 
      textMatrix.append(datum)
      objectList.append(obj)
        
    
    self.waiting = False
    self.scrolledMatrix.update(textMatrix=textMatrix,objectList=objectList)

  def getAriaData(self):
    
    if not self.project:
      showWarning('Failure','No active project', parent=self)
      return
    
    if not self.molSystem:
      showWarning('Failure','No molecular system selected', parent=self)
      return
      
    if len(self.molSystem.chains) > 1:
      code = self.molSystem.findFirstChain().code
      code = askString('Multiple chains','Enter chain code:',code,parent=self) 
      chain = self.molSystem.findFirstChain(code=code)
      
    else:
      chain = self.molSystem.findFirstChain()
    
    peakData   = []
    dihedrals  = []
    ambigs     = []
    unambigs   = []
    jCouplings = []
    
    for obj in self.scrolledMatrix.objectList:
      opts = obj.setupStructureCalcFrameDict
      
      if not opts['use']:
        continue
        
      className = obj.className
      key = getObjectKeyString(obj)
      
      if className == 'PeakList':
        shiftList = opts['shiftList']
        
        if not shiftList:
          continue
      
        datum = {'shifts': {KEY_ID: getObjectKeyString(shiftList)},
                 'peaks':  {KEY_ID: key},'use_assignments': opts['preserve'],
                 'trust_assigned_peaks': opts['filter']}
      
        peakData.append(datum)
      
      elif className == 'DistanceConstraintList':
        if opts['ambig']:
          ambigs.append({KEY_ID: key, 'filter_contributions':opts['filter']})
          
        else:
          unambigs.append({KEY_ID: key, 'filter_contributions':opts['filter']})
        
      elif className == 'JCouplingConstraintList':
        jCouplings.append({KEY_ID: key})
      
      elif className == 'DihedralConstraintList':
        dihedrals.append({KEY_ID: key})   

    projFileName = os.path.join(self.project.url.path, self.project.path)
    dict = {KEY_GENERAL: {'project_name': self.project.name},
         KEY_CCPN: {'filename': projFileName},
         KEY_MOLECULAR_SYSTEM: {KEY_ID: getObjectKeyString(chain) },
         KEY_NOES: peakData,         
         KEY_DIHEDRALS: dihedrals,
         KEY_DISTANCES_AMBIG: ambigs,
         KEY_DISTANCES_UNAMBIG: unambigs,
         KEY_JCOUPLINGS: jCouplings}

    return dict

  def updateAll(self):
    
    self.project = self.application.project
  
    self.updateMolSystems()
    self.updateConstraintSets()
    self.updateAfter()

  def updateAfter(self, obj=None):
  
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)    

  def destroy(self):        

    notify = self.application.unregisterNotify
    for func in ('__init__','delete'):
      notify(self.updateAfter, 'ccp.nmr.Nmr.PeakList', func)
      notify(self.updateConstraintSets, 'ccp.nmr.NmrConstraint.NmrConstraintStore', func)
      notify(self.updateMolSystems, 'ccp.molecule.MolSystem.MolSystem', func)
      notify(self.updateMolSystems, 'ccp.molecule.MolSystem.Chain', func)
 
    for func in ('__init__','delete','setName'):
      notify(self.updateAfter, 'ccp.nmr.Nmr.Experiment', func)
      notify(self.updateAfter, 'ccp.nmr.Nmr.DataSource', func)
  
    Frame.destroy(self)
    

