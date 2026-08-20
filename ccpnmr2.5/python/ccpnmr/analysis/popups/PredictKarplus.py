
"""
======================COPYRIGHT/LICENSE START==========================

PredictKarplus.py: Part of the CcpNmr Analysis program

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
from ccpnmr.analysis.popups.BasePopup       import BasePopup
from ccpnmr.analysis.core.CouplingBasic   import getResidueJCoupling, setResidueJCoupling
from ccpnmr.analysis.core.MoleculeBasic   import getResidueCode

from memops.gui.ButtonList import UtilityButtonList
from memops.gui.FloatEntry      import FloatEntry
from memops.gui.Label           import Label
from memops.gui.LabelFrame      import LabelFrame
from memops.gui.MessageReporter import showOkCancel
from memops.gui.PulldownMenu    import PulldownMenu
from memops.gui.ScrolledMatrix  import ScrolledMatrix

couplingAtoms = [('H','HA'),
                 ('H','C'),
                 ('H','CB'),
                 ('C','HA'),
                 ('C','C'),
                 ('C','CB')]

def testPredictKarplusPopup(argServer):

  popup = PredictKarplusPopup(argServer.parent)
  popup.open()
  

class PredictKarplusPopup(BasePopup):

  def __init__(self, parent, *args, **kw):
  
    self.waiting    = False
    self.molSystem  = None
    self.residue    = None
    self.jCouplings = []
    self.atomNames  = []
    self.measureJCouplingList = None
    self.predictJCouplingList = None
    
    self.coeffAtomNames = []
    self.coefficients = {}
    self.coefficient = None
    for i in range(3):
      self.coefficients[i] = {}
  
    BasePopup.__init__(self, parent=parent, title='Predict Karplus Coefficients from J Coupling', **kw)

    self.updateJCouplingLists()
    self.updateMolSystems()
    self.updateAfter()

  def body(self, guiFrame):
  
    guiFrame.grid_columnconfigure(0, weight=1, minsize=500)
  
    row = 0
     
    frame = LabelFrame(guiFrame, text='General Options')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(1, weight=1)
   
    label = Label(frame, text='Molecular System:')
    label.grid(row=0,column=0,sticky='w')
    self.molSystemPulldown = PulldownMenu(frame, self.changeMolSystem, do_initial_callback=False)
    self.molSystemPulldown.grid(row=0, column=1, sticky='w')
   
    label = Label(frame, text='Measured J coupling list:')
    label.grid(row=0,column=2,sticky='w')
    self.measureCouplingListPulldown = PulldownMenu(frame, self.changeMeasureCouplingList)
    self.measureCouplingListPulldown.grid(row=0, column=3, sticky='w')
   
    label = Label(frame, text='Predicted J coupling list:')
    label.grid(row=0,column=4,sticky='w')
    self.predictCouplingListPulldown = PulldownMenu(frame, self.changePredictCouplingList)
    self.predictCouplingListPulldown.grid(row=0, column=5, sticky='w')
    
    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    
    frame = LabelFrame(guiFrame, text='Angles & 3J Couplings')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)
    
    headingList = ['Residue','Phi','Chi^2']
    editWidgets      = [None, None, None]
    editGetCallbacks = [None, None, None]
    editSetCallbacks = [None, None, None]
    
    
    self.jCouplingEntry = FloatEntry(self,width=8,returnCallback=lambda event: self.setJCoupling())
    i = 0
    for names in couplingAtoms:
      headingList.extend(['Pred\nJ(%s,%s)' % names,'Expt\nJ(%s,%s)' % names])
      editWidgets.extend([None, self.jCouplingEntry])
      editGetCallbacks.extend([None, lambda obj, i=i:self.getJCoupling(obj, i)]) 
      editSetCallbacks.extend([None, lambda event, i=i:self.setJCoupling(i)]) 
      i += 1
      

    self.couplingMatrix = ScrolledMatrix(frame,
                                         headingList=headingList,
                                         editWidgets=editWidgets,
                                         editGetCallbacks=editGetCallbacks,
                                         editSetCallbacks=editSetCallbacks,
                                         callback=self.selectJCoupling,
                                         multiSelect=False)
    self.couplingMatrix.grid(row=0, column=0, sticky='nsew')
    self.couplingMatrix.doEditMarkExtraRules = self.doEditMarkExtraRules
   
    row += 1

    frame = LabelFrame(guiFrame, text='Karplus Coefficients')
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)

    headingList = ['Phase']
    editWidgets      = [None,]
    editGetCallbacks = [None,]
    editSetCallbacks = [None,]
    
    self.coefficientEntry = FloatEntry(self,width=8,returnCallback=lambda event: self.setCoefficient())
    for names in couplingAtoms:
      headingList.append('J(%s,%s)' % names)
      editWidgets.append(self.coefficientEntry)
      editGetCallbacks.append(lambda obj, n=names:self.getCoefficient(obj, n)) 
      editSetCallbacks.append(lambda event, n=names:self.setCoefficient(n)) 
   
    headingList = ['Phase','J(H,HA)','J(H,C)','J(H,CB)','J(C,HA)','J(C,C)','J(C,CB)']
    self.coefficientMatrix = ScrolledMatrix(frame,
                                          headingList=headingList,
                                          editWidgets=editWidgets,
                                          editGetCallbacks=editGetCallbacks,
                                          editSetCallbacks=editSetCallbacks,
                                          callback=self.selectCoefficient,
                                          maxRows=3,
                                          multiSelect=False)
    self.coefficientMatrix.grid(row=0, column=0, sticky='nsew')
  
    row += 1
    texts    = []
    commands = []
    self.bottomButtons = UtilityButtonList(guiFrame, commands=commands, texts=texts,
                                           helpUrl=self.help_url)
    self.bottomButtons.grid(row=row, column=0, sticky = 'ew')


    for func in ('__init__', 'delete','setValue','setError'):
      for clazz in ('ccp.nmr.Nmr.JCoupling',):
        self.registerNotify(self.updateAfter, clazz, func)

    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.JCouplingList',):
        self.registerNotify(self.updateJCouplingLists, clazz, func)


  def getJCouplingLists(self, isSimulated=False):
  
     jCouplingLists = self.nmrProject.findAllMeasurementLists(className='JCouplingList',isSimulated=isSimulated)
   
     return jCouplingLists

  def updateJCouplingLists(self, *obj):
  
    jCouplingLists = self.getJCouplingLists(isSimulated=False)
    names = ['%d:%s' % (jcl.serial, jcl.name or '') for jcl in jCouplingLists]
    names.append('<New>')
    index = -1
    
    if jCouplingLists:
      if self.measureJCouplingList in jCouplingLists:
        index = jCouplingLists.index(self.measureJCouplingList)

      else:
        self.measureJCouplingList = jCouplingLists[0]
        index = 0   
    else:
      self.measureJCouplingList = None

    self.measureCouplingListPulldown.setup(names, index)
  
    jCouplingLists = self.getJCouplingLists(isSimulated=True)
    names = ['%d:%s' % (jcl.serial, jcl.name) for jcl in jCouplingLists]
    names.append('<New>')
    index = -1
    
    if jCouplingLists:
      if self.predictJCouplingList in jCouplingLists:
        index = jCouplingLists.index(self.predictJCouplingList)

      else:
        self.predictJCouplingList = jCouplingLists[0]
        index = 0   
    else:
      self.predictJCouplingList = None
    
    self.predictCouplingListPulldown.setup(names, index)
    
  def changeMeasureCouplingList(self, index, name):

    if name == '<New>':
      self.measureJCouplingList = None
      
    else:
      jCouplingLists = self.getJCouplingLists()
      jCouplingList = jCouplingLists[index]

      if jCouplingList is not self.measureJCouplingList:
        self.measureJCouplingList = jCouplingList

        if self.predictJCouplingList is self.measureJCouplingList:
          for jCouplingList2 in jCouplingLists:
            if self.measureJCouplingList is not jCouplingList2:
              self.predictJCouplingList = jCouplingList2
              break
              
          else:
            self.predictJCouplingList = None

    self.updateAfter()


  def changePredictCouplingList(self, index, name):
  
    if name == '<New>':
      self.predictJCouplingList = None
      self.updateAfter()
      
    else:
      jCouplingLists = self.getJCouplingLists()
      jCouplingList  = jCouplingLists[index]
 
      if jCouplingList is not self.predictJCouplingList:
        self.predictJCouplingList = jCouplingList
        
        if self.predictJCouplingList is self.measureJCouplingList:
          for jCouplingList2 in jCouplingLists:
            if self.predictJCouplingList is not jCouplingList2:
              self.measureJCouplingList = jCouplingList2
              break
              
          else:
            self.measureJCouplingList = None
        
    self.updateAfter()


  def doEditMarkExtraRules(self, obj, row, col):
  
    i = (col-3)/2
  
    if i > -1:
      atoms = couplingAtoms[i]
      residue = obj[0]
    
      for atomName in atoms:
        if not residue.findFirstAtom(name=atomName):
          return False
  
    return True
  
  def open(self):
  
    self.updateMolSystems()
  
    BasePopup.open(self)
  
  def selectJCoupling(self, object, row, col):
  
    self.residue = object[0]
    self.jCouplings = object[1:]
  
 
  def selectCoefficient(self, object, row, col):
  
    self.coefficient = object
    
  
  def getCoefficient(self, object, atomNames):
  
    self.coeffAtomNames = atomNames
  
    coefficient = self.coefficients[object].get(atomNames)
      
    self.coefficientEntry.set(coefficient)
  
  
  def setCoefficient(self, atomNames=None):
  
    if not atomNames:
      atomNames = self.coeffAtomNames
  
    value = self.coefficientEntry.get()
    if atomNames and (self.coefficient is not None):
      self.coefficients[self.coefficient][atomNames] = value
      #print self.coefficient, atomNames,  value
      
    self.updateCoefficients()  
  
  def getJCoupling(self, object, index):
  
    atomNames = couplingAtoms[index]
    self.atomNames = atomNames
  
    coupling = object[index+1]
    
    if not coupling:
      residue  = object[0]
      coupling = self.getResidueJCoupling(self.measureJCouplingList, residue, atomNames)
    
    if coupling:
      value = coupling.value
    else: 
      value = None
      
    self.jCouplingEntry.set(value)
  
  def setJCoupling(self, index=None):

    if not index:
      atomNames = self.atomNames
      index = couplingAtoms.index(atomNames)
    else:
      atomNames = couplingAtoms[index] 
  
    value = self.jCouplingEntry.get()
    if self.residue and atomNames:
      coupling = self.jCouplings[index]
      
      if value is None:
        if coupling:
          if showOkCancel('Confirm',
                          'Really remove J coupling?', parent=self):
                         
            coupling.delete()
       
      else:
        if not coupling:
          self.setResidueJCoupling(self.measureJCouplingList, self.residue,
                                   atomNames, value, isSimulated=False)
 
        else:
          coupling.setValue(value)
          
  def getResidueJCoupling(self, jCouplingList, residue, atomNames):
  
    if jCouplingList:
      return getResidueJCoupling(jCouplingList, residue, atomNames)
        
  
  def setResidueJCoupling(self, jCouplingList, residue, atomNames, value, isSimulated=False):
  
    if jCouplingList is None:
      jCouplingList = self.nmrProject.newJCouplingList(isSimulated=isSimulated)
      
      if isSimulated:
        self.predictJCouplingList = jCouplingList
      else:
        self.measureJCouplingList = jCouplingList
     
    setResidueJCoupling(jCouplingList, self.residue, atomNames, value)
   
   
   
  def updateMolSystems(self, *obj):

    index = -1
    names = []
    
    molSystems = self.project.sortedMolSystems()
    if molSystems:
      
      if self.molSystem not in molSystems:
        self.molSystem = molSystems[0]
      
      names = [ms.code for ms in molSystems]
      index = molSystems.index(self.molSystem) 
    
    
    self.molSystemPulldown.setup(names, index)

  def changeMolSystem(self, index, name):

    molSystems = self.project.sortedMolSystems()
   
    if molSystems:
      molSystem = molSystems[index]
    else:
      molSystem = None
      
    if molSystem is not self.molSystem:
      self.molSystem = molSystem
      self.updateAfter()
        

  def updateAfter(self, object=None):
  
    current = False
    if object and (object.className == 'JCoupling'):
      for resonance in object.resonances:
        resonanceSet = resonance.resonanceSet
        
        if resonanceSet:
          molSystem = resonanceSet.findFirstAtomSet().findFirstAtom().topObject
          
          if molSystem is self.molSystem:
            current = True
            break
        
    else:
      current = True
  
    if self.waiting or (not current):
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
  

  def update(self):
  
    textMatrix = []
    objectList = []
    
    if self.molSystem:
      chains =  self.molSystem.sortedChains()
      
      if len(chains) > 1:
        doChains = False
      else:
        doChains = True
      
      for chain in chains:
        if doChains:
          chainCode = chain.code
        else:
          chainCode = ''  
      
        for residue in chain.sortedResidues():
          name  = '%s%d%s' % (chainCode, residue.seqCode, getResidueCode(residue))
          phi   = 0.0
          chiSq = 0.0
          datum = [name,
                   phi, chiSq]
          
          object = [residue, ]          
          for atomNames in couplingAtoms:
            couplingM = self.getResidueJCoupling(self.measureJCouplingList, residue, atomNames)
            couplingP = self.getResidueJCoupling(self.predictJCouplingList, residue, atomNames)
            
            pred = None
            expt = None
            
            if couplingM:
              expt = couplingM.value

            if couplingP:
              pred = couplingP.value
             
            datum.append(pred)
            datum.append(expt)
            object.append(couplingM)
            
          objectList.append(object)
          textMatrix.append(datum)
    
    self.couplingMatrix.update(textMatrix=textMatrix, objectList=objectList) 
    
    self.updateCoefficients()
    
    self.waiting = False

  def updateCoefficients(self):
    
    textMatrix = []
    objectList = []
    for i in range(3):
      name  = 'C%d' % i
      datum = [name,]
      
      for atomNames in couplingAtoms:
        datum.append( self.coefficients[i].get(atomNames) )
      
      textMatrix.append(datum)
      objectList.append(i)
    
    self.coefficientMatrix.update(textMatrix=textMatrix, objectList=objectList) 

  def destroy(self):

    for func in ('__init__', 'delete','setValue','setError'):
      for clazz in ('ccp.nmr.Nmr.JCoupling',):
        self.unregisterNotify(self.updateAfter, clazz, func)

    for func in ('__init__', 'delete'):
      for clazz in ('ccp.nmr.Nmr.JCouplingList',):
        self.unregisterNotify(self.updateJCouplingLists, clazz, func)
  
    BasePopup.destroy(self)  
