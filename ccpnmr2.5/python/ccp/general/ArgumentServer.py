
"""
======================COPYRIGHT/LICENSE START==========================

ArgumentServer.py: Part of the CcpNmr Analysis program

Copyright (C) 2005 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in ../../../license/CCPN.license.

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
import re

from memops.universal.Io import importModule
from memops.general import Implementation 
ApiError = Implementation.ApiError
from memops.metamodel import ImpConstants
from ccp.general.SelectObject import getDisplayAttrNames

userResidueCodesDict = {'DNA':{'A':'Ade','T':'Thy','G':'Gua','C':'Cyt','U':'Ura'},
                        'RNA':{'A':'Ade','T':'Thy','G':'Gua','C':'Cyt','U':'Ura','I':'Ino'},
                        'protein':{},
                        'carbohydrate':{}
                        }

class ArgumentServer:

  def __init__(self, parent, inGui=None):
  
    if inGui:
      self.popup           = None
      self.dataEntryDialog = importModule('memops.gui.DataEntry')
      self.messageReporter = importModule('memops.gui.MessageReporter')
      try:
        self.chooseObjectDialog = self.objectSelectionPopup
        self.dataModelPopup  = importModule('ccp.general.SelectObject')
      except ImportError:
        self.chooseObjectDialog = self.pickObjectByAttribute
    else:
      self.chooseObjectDialog = self.pickObjectByAttribute
      self.dataEntryDialog = importModule('memops.universal.DataEntry')
      self.messageReporter = importModule('memops.universal.MessageReporter')
  
    self.inGui   = inGui
    self.parent  = parent

  def destroy(self):

    if self.inGui:
      if self.popup:
        self.popup.destroy()
        self.popup = None

    self.dataModelPopup  = None
    self.dataEntryDialog = None
    self.messageReporter = None

  def getFile(self, file_types = None, dismiss_text='Cancel',
              selected_file_must_exist = False):

    if self.inGui:
      module = importModule('memops.gui.FileSelectPopup')
      popup = module.FileSelectPopup(self.parent, file_types=file_types,
                                     dismiss_text=dismiss_text,
                                     selected_file_must_exist=selected_file_must_exist)
      fileName = popup.getFile()
      popup.destroy()
    else:
      fileName = raw_input('Enter file name: ').strip()

    return fileName

  def getDirectory(self, dismiss_text='Cancel', selected_file_must_exist=False):

    if self.inGui:
      module = importModule('memops.gui.FileSelectPopup')
      popup = module.FileSelectPopup(self.parent,
                                     dismiss_text=dismiss_text,
                                     selected_file_must_exist=selected_file_must_exist)
      directory = popup.getDirectory(needFileSelected=False)
      popup.destroy()
    else:
      directory = raw_input('Enter directory: ').strip()

    return directory

  def getProgressBar(self, text='', progress=0, total=100,
                     width=200, height=15, transient=True):
  
    if self.inGui:
      from memops.gui.ProgressBar import ProgressBar
      return ProgressBar(self.parent, text=text, progress=progress, total=total,
                         width=width, height=height, transient=transient)

  def chooseObject(self, objects, key=None, labels=None, objectName=None,
                   refClassName=None):
  
    if len(objects) == 0:
      return None
    elif len(objects) == 1:
      return objects[0]
    else:
      return self.chooseObjectDialog(objects, key=key, labels=labels, 
                                     objectName=objectName, 
                                     refClassName=refClassName)

  def objectSelectionPopup(self,objects,key=None,labels=None,objectName=None,
                           refClassName=None):
  
    self.parent.argumentVar = Tkinter.IntVar()    
    self.parent.argumentObject = None   
    if self.popup:
      self.popup.update(objects, key=key, labels=labels, objectName=objectName,
                        refClassName=refClassName)
      self.popup.open()
    else:
      self.popup = self.dataModelPopup.SelectObjectPopup(self.parent,objects,
                                                         key=key,
                                                         labels=labels,
                                                         objectName=objectName,
                                                         refClassName=refClassName)
    self.popup.wait_variable(self.parent.argumentVar)
    obj = self.parent.argumentObject

    del self.parent.argumentVar
    del self.parent.argumentObject

    return obj

  def askString(self, message, default):
  
    return self.dataEntryDialog.askString('String Input',message,default,parent=self.parent)

  def askInteger(self, message, default):
  
    return self.dataEntryDialog.askInteger('Integer Input',message,default,parent=self.parent)

  def askFloat(self, message, default):
  
    return self.dataEntryDialog.askFloat('Float Input',message,default,parent=self.parent)
  
  def askYesNo(self, message):
  
    return self.messageReporter.showYesNo('Question', message,parent=self.parent)

  def showWarning(self, message):
  
    return self.messageReporter.showWarning('Warning', message,parent=self.parent)

  def showInfo(self, message):
  
    return self.messageReporter.showInfo('Info', message,parent=self.parent)
 
  def getSpectrum(self,spectra=None):

    project = self.getProject()
    labels = {}
    if spectra is None:
      spectra = []
      for experiment in project.currentNmrProject.sortedExperiments():
        for dataSource in experiment.sortedDataSources():
          spectra.append(dataSource)
      if len(spectra) < 1:
        self.messageReporter.showWarning('Warning','No spectra in current project')
        return

    nameList = []
    for spec in spectra:
      nameList.append(spec.experiment.name+':'+spec.name)
      
    labels['Name'] = nameList

    return self.chooseObject(spectra, key='Name', labels=labels,
                             objectName='Spectrum')

  def getPeakList(self,spectrum=None):
  
    peakList = None
    if spectrum is None:
      spectrum = self.getSpectrum()

    if not spectrum:
      return
            
    N = len(spectrum.peakLists)
    if N == 1:
      peakList = spectrum.findFirstPeakList()
    else:
      peakList = self.chooseObject(spectrum.sortedPeakLists())
      
    return peakList
  
  def getPeak(self, spectrum=None):
  
    if spectrum is None:
      spectrum = self.getSpectrum()
    
    if not spectrum:
      return
      
    peakList = self.getPeakList(spectrum=spectrum)
    return self.chooseObject(peakList.sortedPeaks())

  def getProject(self):
    
    parent = self.parent
    if hasattr(parent,'project'):
      project = parent.project
      if not project:
        self.messageReporter.showWarning('Warning','No project active')
        return
    else:
      self.messageReporter.showWarning('Warning',
                                       'No project defined - code error')
      return
    #
    return project
    
  def getResonance(self, resonances=None):
  
    labels = {}
    project = self.getProject()
    if resonances is None:
      resonances = project.currentNmrProject.sortedResonances()
  
    annotations = []
    for resonance in resonances:
      annotations.append(self.getResonanceGuiName(resonance))
      
    labels['Annotation'] = annotations
      
    return self.chooseObject(resonances,labels=labels)

  def getViolationList(self, constraintSet=None):
  
    project = self.getProject()
    if constraintSet is None:
      constraintSet = self.chooseObject(project.sortedNmrConstraintStores())

    if constraintSet:
      violationList = self.chooseObject(constraintSet.sortedViolationLists())
      return violationList

    return None

  def getConstraintSet(self):
  
    nmrProject = self.getProject().currentNmrProject
    constraintSet = self.chooseObject(nmrProject.sortedNmrConstraintStores(), 
                                      objectName='Constraint Set')

    return constraintSet

  def getConstraintList(self, constraintSet=None):
  
    if constraintSet is None:
      constraintSet = self.getConstraintSet()

    if constraintSet:
      constraintList = self.chooseObject(constraintSet.sortedConstraintLists(),
       refClassName='ccp.nmr.NmrConstraint.AbstractConstraintList'
      )
      return constraintList

    return None

  def getConstraint(self, constraintSet=None, constraintList=None):
  
    def getResonanceGuiName(resonance):
      if not hasattr(resonance, 'guiName'):
        resonance.guiName = self.getResonanceGuiName(resonance)
      return resonance.guiName

  
    if constraintList is None:
      constraintList = self.getConstraintList(constraintSet=constraintSet)
    
    if constraintList:
      constrType  = constraintList.className[:-14]
      constraints = constraintList.sortedConstraints()
    
      resNames = []
      for constraint in constraints:
        texts = []
        items = constraint.sortedItems()
        
        if constrType == 'Dihedral':
          resonanceNames = [getResonanceGuiName(r) for r in constraint.resonances]
          resonanceText  = '%s-%s--%s-%s' % tuple(resonanceNames)

        elif constrType == 'ChemShift':
          resonanceText = getResonanceGuiName(constraint.resonance)
        else:
          resonanceNames = [getResonanceGuiName(r) for r in items[0].resonances]
          resonanceText  = '%s-%s' % tuple(resonanceNames)
        texts.append(resonanceText)
        
        text = str(constraint.serial) + ': ' + ','.join(texts) 
        if len(items) > 1:
          text += '...'
          
        resNames.append(text) 
 
      labels = {}
      labels['Name'] = resNames

      #return self.chooseObject(chains,key='Name',labels=labels)
      return self.chooseObject(constraints,key='serial',labels=labels)

    return None

  def getStructure(self):
  
    project = self.getProject()
    
    structures = []
    for molSystem in project.molSystems:
      for structure in molSystem.sortedStructureEnsembles():
        structures.append(structure)

    return self.chooseObject(structures)

  def getStructureGeneration(self):
  
    nmrProject = self.getProject().currentNmrProject
    return self.chooseObject(nmrProject.sortedStructureGenerations())
    
  def getExperiment(self):
  
    project = self.getProject()
    return  self.chooseObject(project.currentNmrProject.sortedExperiments())
    
  def getShiftList(self):
         
    return self.getMeasurementList('ShiftList')

  def getMeasurementList(self, className=None):
  
    project = self.getProject()
    mLists = []
    for measurementList in project.currentNmrProject.sortedMeasurementLists():
      if (not className) or (measurementList.className == className):
        mLists.append(measurementList)
         
    return self.chooseObject(mLists)
    
  def getSpinSystem(self):

    nmrProject = self.getProject().currentNmrProject
    spinSystems = nmrProject.sortedResonanceGroups()
    labels = {}
    ll = labels['residue'] = ['']*len(spinSystems)
    for ii,xx in enumerate(spinSystems):
      res = xx.residue
      if res is not None:
        ll[ii] = ':'.join(str(x) for x in res.getFullKey())
    
    return self.chooseObject(nmrProject.sortedResonanceGroups(), 
                             objectName='SpinSystem', labels=labels)
    # put residue
  def getMolSystem(self):
  
    project = self.getProject()
    molSystems = project.sortedMolSystems()

    return self.chooseObject(molSystems,key='code')
 
  def getChain(self, molSystem = None):
  
    project = self.getProject()
    if molSystem is not None:
      molSystems = [molSystem,]
    else:
      molSystems = project.sortedMolSystems()
    
    chains = []
    names  = []
    for molSystem in molSystems:
      for chain in molSystem.chains:
        chains.append(chain)
        names.append('%s:%s' % (molSystem.code, chain.code))

    labels = {}
    labels['Name'] = names  

    return self.chooseObject(chains,key='Name',labels=labels)
    
  def getResidue(self, chain = None):
  
    if chain is None:
      chain = self.getChain()
    
    labels = {}
    tlcs = []  
    for residue in chain.residues:
      tlcs.append(self.getResidueCode(residue))
      
    labels['Three Letter Code'] = tlcs  
      
    return self.chooseObject(chain.sortedResidues(),labels=labels)
    
  def getAtomSet(self, residue = None):
  
    project = self.getProject()
    if not residue:
      residue = self.getResidue()
    
    atomSets = []
    for atom in residue.atoms:
      atomSet = atom.atomSet
      if atomSet:
        if atomSet not in atomSets:
          atomSets.append(atomSet)
    
    if not atomSets:
      atomSets = project.currentNmrProject.sortedAtomSets()
      
    return self.chooseObject(atomSets)

  def getAtom(self, residue = None):
  
    if not residue:
      residue = self.getResidue()
      
    atoms = residue.sortedAtoms()
    return self.chooseObject(atoms)

  def getIsotope(self):
  
    labels  = {}
    project = self.getProject()
    isotopes = []
    chemElementNames = []
    isotopeNames = []
    for chemElement in project.currentChemElementStore.chemElements:
      for isotope in chemElement.isotopes:
        isotopes.append(isotope)
        chemElementNames.append(chemElement.symbol)
        isotopeNames.append(str(isotope.massNumber) + chemElement.symbol)
    labels['element'] = chemElementNames 
    labels['isotope'] = isotopeNames
    
    return self.chooseObject(isotopes,key='isotope', labels=labels)
    
  def pickObjectByAttribute(self, objects, key=None, labels=None, 
                            objectName=None, refClassName=None):
    """ 9/1/09 Rasmus Fogh
      changed __dict__ access to getattr
      changed explicit loops to ','.join(alist)
      changed to work with refClass and external call to getDisplayAttrNames
    """
  
    assert objects
    
    if len(objects) == 1:
      return objects[0]
      
    if refClassName:
      apiClass = Implementation.getClassFromFullName(refClassName)
      refClass = apiClass._metaclass
    else:
      apiClass = objects[0].__class__
      refClass = objects[0].metaclass
      refClassName = refClass.qualifiedName()

    for obj in objects:
      if not isinstance(obj, apiClass):
        raise ApiError("Object %s is not instance of %s" % (obj, refClassName))
      
    if not objectName:
      objectName = refClass.name 
    
    if labels:
      for label,ll in labels.items():
        if len(ll) != len(objects):
          msg = 'Length of %s label list must be the same as the number of objects' % label
          raise Exception(msg)
      labels = labels.copy()
    else:
      labels = {}
    
    # get attribute names for display
    attrNames, keyNames = getDisplayAttrNames(refClass, labels, minimal=True)
  
    # set up for key parameter and fullkey
    # fullKey will be used as key if none othe is present,
    # and will be put at the end otherwise
    if key is None:
      if keyNames:
        attrNames.append('fullKey')
      else:
        keyNames = ['fullKey']
        # key in class not available, no key passed in. Use fullKekey = 'fullKey'
        attrNames.insert(0,'fullKey')
    else:
      try:
        attrNames.remove(key)
      except ValueError:
        raise ApiError("key %s not in attributes or labels: %s"
                        % (key, attrNames))
      keyNames = [key]
      attrNames.insert(0,key)
      attrNames.append('fullKey')
    
    # add 'fullKey to labels
    labels['fullKey'] = [str(getattr(obj, 'getFullKey')()) 
           for obj in objects]

    # set standard texts
    keyText = ','.join(keyNames)
    attrText = ','.join(attrNames)
    text = 'Available %s objects [%s]:\n' % (objectName, attrText)
    
    # set default text
    ll = []
    for keyName in keyNames:
      xx = labels.get(keyName)
      if xx:
        val = xx[0]
      else:
        try:
          val = getattr(objects[0],keyName)
          if isinstance(val, frozenset):
            val = list(val)
            val.sort()
        except AttributeError:
          self.messageReporter.showError('Error',
                                         'Object %s lacks key attribute %s' 
                                         % (obj, keyName))
      ll.append(val)
    defaultText = ','.join(str(x) for x in ll)
    
    # set objectTexts list
    objectTexts = []
    for ii,obj in enumerate(objects):
      ll = []
      for attrName in attrNames:
        xx = labels.get(attrName)
        if xx:
          val = xx[ii]
        else:
          try:
            val = getattr(obj,attrName)
            if isinstance(val, frozenset):
              val = list(val)
              val.sort()
          except AttributeError:
            self.messageReporter.showError('Error',
                                           'Object %s lacks attribute %s' 
                                           % (obj, attrName))
        ll.append(val)
      objectTexts.append(ll)
    objectTexts.sort()
    
    # truncate individual fields - NB must be done after sorting
    for ll in objectTexts:
      for ii,val in enumerate(ll):
        val = str(val)
        if len(val) > 32:
          val = val[:32] + '...'
        ll[ii] = val
        
    sortedObjectTexts = [','.join(str(x) for x in ll) for ll in objectTexts]
    
    i = 0
    for objectText in sortedObjectTexts:
      i += len(objectText)
      if i>80:
        text = text + '\n'
        i = 0
      text = text + ' ' + objectText + ';'

    text = text + '\n\nChoose %s by %s\n' % (objectName,keyText)
    
    # get selection input
    try:
      stringInput = self.dataEntryDialog.askString('Choose %s' % objectName,text,defaultText)
    except Exception:
      return
    inputVals = re.split(',',stringInput)
    
    print 'inputVals', inputVals
    print 'labels', labels
    print 'keyNames', keyNames
    
    # select object
    for ii,obj in enumerate(objects):
      for jj,keyName in enumerate(keyNames):
        inputVal = inputVals[jj]
        xx = labels.get(keyName)
        if xx:
          if xx[ii] != inputVal:
            # mismatch
            break
        elif not hasattr(obj, keyName):
          # mismatch
          break
        
        else:
          val = getattr(obj, keyName)
          if isinstance(val, frozenset):
            val = list(val)
            val.sort()
            if inputVal != str(val):
              # mismatch
              break
      else:
        # we have a match
        return obj
    
        
    self.messageReporter.showWarning('Warning','No %s with matching %s found' % (objectName,keyText))
    return
  
  def getResonanceGuiName(self, resonance, fullName=True):
    """ default resonance name
    """
    
    if fullName:
      if hasattr(resonance, 'guiName'):
        return resonance.guiName
    
    if hasattr(resonance, 'label'):
      return resonance.label
    else:  
      return '[%d]' % resonance.serial
  
  def getResidueCode(self, obj):
    """generic version of getResidueCode
    """
  
    ccpCode = obj.ccpCode
    molType = obj.molType
 
    if molType is None:
      molType = 'protein'
 
    ccpCodeDict = userResidueCodesDict.get(molType, {})
 
    return ccpCodeDict.get(ccpCode, ccpCode)
        
