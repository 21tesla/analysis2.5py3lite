"""
======================COPYRIGHT/LICENSE START==========================

EditResStructures.py: Part of the CcpNmr Analysis program

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
import re
import os

from math import sqrt

from memops.general import Implementation

from memops.universal.Io import joinPath, splitPath

from memops.gui.ButtonList import UtilityButtonList
from memops.gui.DataEntry import askString
from memops.gui.FloatEntry import FloatEntry
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Label import Label
from memops.gui.MessageReporter import showOkCancel
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.Spacer import Spacer

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.StructureBasic import getAtomSetCoords

class EditResStructuresPopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.guiParent      = parent
    self.structure      = None
    self.constraintSet = None
    self.cloud          = None
    self.cloudRmsdDict  = {}
    self.strucRmsdDict  = {}
    self.waiting        = False
    
    BasePopup.__init__(self, parent=parent, title="Resonance Cloud Structures",**kw)
    
  def body(self, guiFrame):
    
    row = 0
    guiFrame.grid_columnconfigure(1, weight=1)
    guiFrame.grid_rowconfigure(0, weight=0)
    guiFrame.grid_rowconfigure(1, weight=1)
    
    self.generationLabel    = Label(guiFrame, text = 'Structure Generation:')
    
    constraintSets     = []
    constraintSetNames = []
    index = -1
    for constraintSet in self.project.nmrConstraintStores:
      index += 1
      constraintSets.append(constraintSet)
      constraintSetNames.append(str(constraintSet.serial))
      self.constraintSet = constraintSet
     
    self.constrSetPulldown = PulldownMenu(guiFrame, self.changeConstraintSet, constraintSetNames, selected_index=index, do_initial_callback=False)

    self.generationLabel.grid   (row=row, column=0, columnspan=1, sticky = 'e')
    self.constrSetPulldown.grid(row=row, column=1, columnspan=1, sticky = 'w')
 
    strucLabel = Label(guiFrame, text='Comparison structure')
    strucLabel.grid(row=row, column=2, sticky='e')
    self.strucPulldown = PulldownMenu(guiFrame, entries=self.getStructures(), callback=self.setStructure, selected_index=0, do_initial_callback=False)
    self.strucPulldown.grid(row=row, column=3, sticky='w')
    
    sdTolLabel = Label(guiFrame, text='Tolerance (SDs):')
    sdTolLabel.grid(row=row, column=4, sticky='e')
    self.sdToleranceEntry = FloatEntry(guiFrame, text=2.0, width=6)
    self.sdToleranceEntry.grid(row=row, column=5, stick='w') 
    
    row += 1
    colHeadings = ['#','File name','RMSD to mean','RMSD to structure']
    self.scrolledMatrix = ScrolledMatrix(guiFrame,initialRows=10, headingList=colHeadings, callback=self.selectCell, objectList=[], textMatrix=[[],])
    self.scrolledMatrix.grid(row = row, column = 0, columnspan=6, sticky='nsew')
    
    row += 1
    texts    = ['Calc\nRMSD','Make Cloud\nfrom structure','Delete','Delete\nbad']
    commands = [self.calcRmsd, self.makeStrucCloud,self.deleteCloud,self.filterClouds]
    self.bottomButtons = UtilityButtonList(guiFrame, texts=texts, expands=False,
                                           commands=commands, helpUrl=self.help_url)
    self.bottomButtons.grid(row = row, column = 0, columnspan=6, sticky = 'nsew')
    self.update()

    for func in ('__init__', 'delete'):
      self.registerNotify(self.updateStructureGen, 'ccp.nmr.Nmr.StructureGeneration', func)

  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)
    

  def getStructures(self):
  
    names = ['<None>',]
    for molSystem in self.project.sortedMolSystems():
      for structure in molSystem.sortedStructureEnsembles():
        names.append('%s:%d' % (molSystem.name,structure.ensembleId) )
        
    return names
  
  def setStructure(self, index, name=None):

    if index < 1:
      self.structure = None
    else:
      structures = []
      for molSystem in self.project.sortedMolSystems():
        for structure in molSystem.sortedStructureEnsembles():
          structures.append( structure )
 
      self.structure = structures[index-1]

  def filterClouds(self):
  
    if self.constraintSet:
      sdTolerance = self.sdToleranceEntry.get() or 2.0
      keptClouds  = []
      clouds = self.guiParent.application.getValues(self.constraintSet, 'clouds')

      meanGroupRmsds = []
      for cloud in clouds:
        rmsd = self.cloudRmsdDict.get(cloud)
        if rmsd is not None:
          meanGroupRmsds.append(rmsd)
        
      meanRmsd = 0.0
      N = 0
      for rmsd in meanGroupRmsds:
        meanRmsd += rmsd
        N += 1
 
      if N > 0:
        meanRmsd /= float(N)
 
      sd = 0.0
      for rmsd in meanGroupRmsds:
        sd += (rmsd-meanRmsd) * (rmsd-meanRmsd)
 
      if N > 0:
        sd /=float(N-1)
 
      sd = sqrt(sd)
 
      print meanRmsd, '+/-', sd
 
      n = 0
      for cloud in clouds:
        rmsd = self.cloudRmsdDict.get(cloud)
        if rmsd is None:
          keptClouds.append(cloud)
        elif abs(rmsd - meanRmsd) > (sdTolerance*sd):
          print 'Cloud %s is bad' % (cloud)
        else:
          keptClouds.append(cloud)
            
      self.guiParent.application.setValues(self.constraintSet, 'clouds', values=keptClouds)
      self.updateAfter()
 
  def makeStrucCloud(self):

    if self.structure:
      serials = self.guiParent.application.getValues(self.constraintSet, 'cloudsResonances')
      pdbFileName = 'CloudForStructure.pdb'
      #from ccpnmr.clouds.AtomCoordList import AtomCoordList
      from ccpnmr.c.AtomCoordList import AtomCoordList
      atomCoordList = AtomCoordList()
      resDict = {}
      hmass = 25
      
      print "L1", len(serials)
      
      for resonance in self.nmrProject.resonances:
        resDict[resonance.serial] = resonance
 
      print "L2", len(resDict)
      
      resonances = []
      for serial in serials:
        if resDict.get(serial) is not None:
          resonances.append( resDict[serial] )

      print "L3", len(resonances)
      
      C = 0
      for resonance in resonances:
        resonanceSet = resonance.resonanceSet
        if resonanceSet:
          i = resonanceSet.sortedResonances().index(resonance)
          atomSet = resonance.resonanceSet.sortedAtomSets()[i]
          coords = getAtomSetCoords(atomSet, self.structure)
          coord = coords[0]
          atomCoordList.add(hmass, coord.x, coord.y, coord.z)
          
          C +=1
      print "L4", len(atomCoordList)
    
      from ccpnmr.clouds.FilterClouds import writeTypedPdbCloud
      writeTypedPdbCloud(atomCoordList, pdbFileName, resonances)
      
      clouds = self.guiParent.application.getValues(self.constraintSet, 'clouds')
      clouds.append(pdbFileName)
      
      self.guiParent.application.setValues(self.constraintSet, 'clouds', values=clouds)
      self.updateAfter()

  def calcRmsd(self):
    
    if self.constraintSet:
      clouds = self.guiParent.application.getValues(self.constraintSet, 'clouds')
      from ccpnmr.clouds.FilterClouds import filterClouds
      rmsds  = filterClouds(clouds) 
      n = len(clouds)
      
      for i in range(n):
        cloud = clouds[i]
        rmsd  = rmsds[i] 
        self.cloudRmsdDict[cloud] = rmsd
        
      self.updateAfter()

  def changeConstraintSet(self, i, name):
  
    project = self.project
    if project.nmrConstraintStores:
      constraintSet = project.nmrConstraintStores[i]
    else:
      constraintSet = None
    
    if constraintSet is not self.constraintSet: 
      self.constraintSet  = constraintSet
      self.cloud = None
    
    self.updateAfter()
          
  def deleteCloud(self):

    if self.constraintSet and self.cloud and showOkCancel('Confirm','Really delete resonance cloud?', parent=self):
      clouds = self.guiParent.application.getValues(self.constraintSet, 'clouds')
      if clouds:
        clouds.remove(self.cloud)
        self.cloud = None
        self.guiParent.application.setValues(self.constraintSet, 'clouds', values=clouds)
        self.updateAfter()
      
  def selectCell(self, cloud, row, col):
    
    self.cloud = cloud
    self.bottomButtons.buttons[1].enable()
     
  def updateAfter(self, *opt):
  
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)

  def getConstraintSetNames(self):
  
    names = []
    constraintSets = self.project.nmrConstraintStores
    for set in constraintSets:
      names.append('%d' % set.serial)
      
    return names

  def updateStructureGen(self, *opt):
   
    project    = self.project
    constraintSets = self.project.sortedNmrConstraintStores
    
    if constraintSets:
      constraintSetNames = self.getConstraintSetNames()
      
      # set defaults
      if self.constraintSet not in constraintSets:
        self.constraintSet = constraintSets[0]
        self.cloud = None
                 
      i =  constraintSets.index(self.constraintSet)
      self.constrSetPulldown.setup(constraintSetNames,i)

    else:
      self.constraintSet  = None
      self.cloud = None
      self.constrSetPulldown.setup([],-1)
    
  def destroy(self):

    for func in ('__init__', 'delete'):
      self.unregisterNotify(self.updateStructureGen, 'ccp.nmr.Nmr.StructureGeneration', func)

    BasePopup.destroy(self)
    
  def update(self):
       
    objectList  = []
    textMatrix  = []
    if self.constraintSet:
      clouds = self.guiParent.application.getValues(self.constraintSet, 'clouds')
      if clouds:
        objectList = list(clouds)
     
    i = 0
    for cloud in objectList:
      i += 1
      datum = []
      datum.append(i)
      datum.append(cloud)
      datum.append(self.cloudRmsdDict.get(cloud) or '-')
      datum.append(self.strucRmsdDict.get(cloud) or '-')      
      textMatrix.append( datum )
            
    if not self.cloud: 
      self.bottomButtons.buttons[1].disable()

    self.scrolledMatrix.update(objectList=objectList, textMatrix=textMatrix)
    self.waiting = False
 
