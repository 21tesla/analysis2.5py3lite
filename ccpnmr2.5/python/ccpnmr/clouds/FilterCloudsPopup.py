"""
======================COPYRIGHT/LICENSE START==========================

FilterCloudsPopup.py: Part of the CcpNmr Clouds program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
A. Grishaev and M. Llinas (2002).
CLOUDS, a protocol for deriving a molecular proton density via NMR.
Proc Natl Acad Sci USA. 99, 6707-6712.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""
import Tkinter
import re
import os

from math import sqrt

from memops.general import Implementation

from memops.gui.DataEntry import askString
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Label import Label
from memops.gui.MessageReporter import showOkCancel, showWarning
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.Spacer import Spacer
from memops.gui.Util import createDismissHelpButtonList
from memops.universal.Io import joinPath, splitPath
from memops.universal.Geometry import inverseMatrix

from ccpnmr.clouds.CloudBasic import writeTypedPdbCloud, getCloudsFromFile, getFileNamesFromPattern
from ccpnmr.clouds.FilterClouds import centerCoords, alignToMeanCloud, matrixVecMultiply, filterClouds, alignClouds, alignCloudsToRef, alignCoordinates, getMeanCoords

#from ccpnmr.analysis.ArgumentServer import ArgumentServer
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.StructureBasic import getAtomSetCoords

atomTypeList = [None, ['H'],['H','HA','HA1','HA2'],['H','HA','HA1','HA2','HB','HB1','HB2']]

class FilterCloudsPopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.guiParent  = parent
    self.structure  = None
    self.name       = None
    self.clouds     = []
    self.rmsds      = []
    self.names      = []
    self.atomTypes  = None
    self.waiting = 0
    
    BasePopup.__init__(self, parent=parent, title="Filter Clouds",**kw)
    
  def body(self, guiFrame):
    
    row = 0
    guiFrame.grid_columnconfigure(3, weight=1)

    label  = Label(guiFrame, text = 'Cloud file names:')
    label.grid   (row=row, column=0, sticky = Tkinter.W)
    self.fileNameEntry = Entry(guiFrame, text='testHistone\d+.pdb', returnCallback=self.loadClouds)
    self.fileNameEntry.grid(row=row, column=1, sticky = Tkinter.W)
 
    strucLabel = Label(guiFrame, text='Comparison structure')
    strucLabel.grid(row=row, column=2, sticky=Tkinter.W)
    self.strucPulldown = PulldownMenu(guiFrame, entries=self.getStructures(), callback=self.setStructure, selected_index=0, do_initial_callback=0)
    self.strucPulldown.grid(row=row, column=3, sticky=Tkinter.W)
    
    row += 1
    sdTolLabel = Label(guiFrame, text='Tolerance (SDs):')
    sdTolLabel.grid(row=row, column=0, sticky=Tkinter.W)
    self.sdToleranceEntry = FloatEntry(guiFrame, text=2.0, width=6)
    self.sdToleranceEntry.grid(row=row, column=1, stick=Tkinter.W) 
    
    atomTypes = ['All','H','H HA','H HA HB']
    label = Label(guiFrame, text='RMSD Atom Types:')
    label.grid(row=row, column=2, sticky=Tkinter.W)
    self.atomsPulldown = PulldownMenu(guiFrame, entries=atomTypes, callback=self.setAtomTypes, selected_index=0, do_initial_callback=0)
    self.atomsPulldown.grid(row=row, column=3, sticky=Tkinter.W)
    
    
    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    colHeadings = ['#','File name','RMSD to mean']
    self.scrolledMatrix = ScrolledMatrix(guiFrame,initialRows=10, headingList=colHeadings,
                                         callback=self.selectCell, objectList=[],
                                         textMatrix=[[],], multiSelect=1)
    self.scrolledMatrix.grid(row = row, column = 0, columnspan=4, sticky=Tkinter.NSEW)
    
    row += 1
    texts    = ['Load\nClouds','Align\nClouds','Calc\nRMSD','Make Cloud\nfrom structure','Remove','Remove\nbad']
    commands = [self.loadClouds,self.alignClouds,self.calcRmsd, self.makeStrucCloud,self.deleteClouds,self.filterClouds]
    self.bottomButtons = createDismissHelpButtonList(guiFrame, texts=texts, expands=1,
                           commands=commands, help_url=self.help_url)
    self.bottomButtons.grid(row = row, column = 0, columnspan=4, sticky = Tkinter.NSEW)
    self.update()

  def alignClouds(self):
  
    pattern     = self.fileNameEntry.get()
    self.names  = getFileNamesFromPattern(pattern, '.')
    self.clouds = getCloudsFromFile(self.names, self.guiParent.project)
    alignClouds(self.clouds, self.names)

  def loadClouds(self):
  
    pattern     = self.fileNameEntry.get()
    self.names  = getFileNamesFromPattern(pattern, '.')
    self.clouds = getCloudsFromFile(self.names, self.guiParent.project)
    self.name   = None
    self.rmsds  = [None for x in range(len(self.clouds))]
    self.updateAfter()

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
      for molSystem in self.project.molSystems:
        for structure in molSystem.structures:
          structures.append( structure )
 
      self.structure = structures[index-1]
      
    self.updateButtons()

  def setAtomTypes(self, index, name=None):

    self.atomTypes = atomTypeList[index]

  def filterClouds(self):
  
    if self.clouds:
      sdTolerance = self.sdToleranceEntry.get() or 2.0
      keptClouds  = []
        
      meanRmsd = 0.0
      N = 0
      for rmsd in self.rmsds:
        meanRmsd += rmsd or 0.0
        N += 1
 
      if N > 0:
        meanRmsd /= float(N)
 
      sd = 0.0
      for r in self.rmsds:
        rmsd = r or 0.0
        sd += (rmsd-meanRmsd) * (rmsd-meanRmsd)
 
      if N > 0:
        sd /=float(N-1)
 
      sd = sqrt(sd)
 
      print meanRmsd, '+/-', sd
 
      for i in range(len(self.clouds),0,-1):
        rmsd = self.rmsds[i]
        if abs(rmsd - meanRmsd) > (sdTolerance*sd):
          self.rmsds.pop(i)
          self.names.pop(i)
          self.clouds.pop(i)
          #print 'Cloud %s is bad' % (cloud)

      self.updateAfter()
 
  def makeStrucCloud(self):

    if self.structure and self.clouds:
      pdbFileName   = 'CloudForStructure.pdb'
      atomCoordList = []
      atomCoordList0 = []
      resDict = {}
      hmass = 25
      resonances  = self.clouds[0].keys()
      resonances2 = []
      
      C = 0
      for resonance in resonances:
        if resonance == 'rmsd':
          continue
      
        resonanceSet = resonance.resonanceSet
        if resonanceSet:
          i = list(resonanceSet.resonances).index(resonance)
          atomSet = list(resonance.resonanceSet.atomSets)[i]
          coords = getAtomSetCoords(atomSet, self.structure)
          coord = coords[0]
          atomCoordList.append(  [coord.x, coord.y, coord.z] )
          atomCoordList0.append( [coord.x, coord.y, coord.z] )
          resonances2.append(resonance)
          
          C +=1
    
      print len(atomCoordList)
      print len(resonances), len(resonances2)

      print "Generating Mean"
      cloudsList = []
      for cloud in self.clouds:
        orderCloud = []
        for resonance in resonances2:
          x,y,z = cloud.get(resonance) or (0.0,0.0,0.0)
        
          orderCloud.append([-x,-y,-z])
        cloudsList.append(orderCloud)
      (meanCloud,cloudsList) = alignToMeanCloud(cloudsList)

      weights = [1.0 for x in atomCoordList]
      centerCoords(atomCoordList)
      print "init cen", getMeanCoords(atomCoordList)
      print "mean cen", getMeanCoords(meanCloud)
      
      print "Print aligning struct clouds to mean", len(meanCloud), len(atomCoordList), len(weights)
      atomCoordsList, error, rotMat = alignCoordinates(meanCloud,atomCoordList,weights)    

      print "  Rotation", rotMat
      writeTypedPdbCloud(atomCoordList, pdbFileName, resonances2)
      
      print "Getting centres"
      oldCentre = getMeanCoords(atomCoordList0)
      newCentre = getMeanCoords(atomCoordList)
      delta     = [newCentre[i]-oldCentre[i] for i in range(len(oldCentre))]
      
      print "  New centre", newCentre
      print "  Old centre", oldCentre
      print "  Delta", delta
      
      #inverseRot = inverseMatrix(rotMat)

      model = self.structure.findFirstModel()
      coordinates = model.coordinates
      offset = 0
      iis = (0,1,2)
      for atom in self.structure.orderedAtoms:
        next = offset + 3
        coords = [coordinates[offset+ii] + delta[ii] for ii in iis]
        coords  = matrixVecMultiply(rotMat, coords)
        coordinates[offset:next] = coords
        offset = next
      model.setSubmatrixData('coordinates',  coordinates) 
      
      clouds = getCloudsFromFile([pdbFileName,], self.structure.root)
      self.clouds.append(clouds[0])
      self.rmsds.append(None)
      self.names.append(pdbFileName)
      
      self.updateAfter()

  def calcRmsd(self):
    
    if self.clouds:
    
      if len(self.scrolledMatrix.currentObjects) < 2:
        clouds = self.clouds  
      else:
        clouds = []
        for name in self.scrolledMatrix.currentObjects:
          clouds.append(self.clouds[self.names.index(name)])
    
      self.rmsds  = filterClouds(clouds, atomTypes=self.atomTypes) 
      self.updateAfter()

          
  def deleteClouds(self):

    if self.names and self.name and showOkCancel('Confirm','Really remove selected clouds?'):
      indices = []
      for name in self.scrolledMatrix.currentObjects:
        i = self.names.index(name)
        indices.append(i)
      indices.sort()
      indices.reverse()

      for i in indices:
        self.clouds.pop(i)
        self.rmsds.pop(i)
        self.names.pop(i)
        
      self.name = None
      self.updateAfter()
         
  def selectCell(self, name, row, col):
    
    self.name = name
    self.updateButtons()
     
  def updateAfter(self, *opt):
  
    if self.waiting:
      return
    else:
      self.waiting = 1
      self.after_idle(self.update)
    
  def destroy(self):

    BasePopup.destroy(self)
  
  def updateButtons(self):
  
    if self.names:
      self.bottomButtons.buttons[1].enable()
      self.bottomButtons.buttons[2].enable()
      self.bottomButtons.buttons[5].enable()
    else:
      self.bottomButtons.buttons[1].disable()
      self.bottomButtons.buttons[2].enable()
      self.bottomButtons.buttons[5].disable()
            
    if self.name:
      self.bottomButtons.buttons[4].enable()
    else:
      self.bottomButtons.buttons[4].disable()

    if self.structure and self.clouds:
      self.bottomButtons.buttons[3].enable()
    else:
      self.bottomButtons.buttons[3].disable()
    
  def update(self):
       
    textMatrix = []
    objectList = self.names
    self.updateButtons()
     
    i = 0
    for name in objectList:
      datum = []
      datum.append(i+1)
      datum.append(name)
      datum.append(self.rmsds[i])
      textMatrix.append( datum )
      i += 1
            
    if not objectList:
      textMatrix = [[],]
       
    self.scrolledMatrix.update(objectList=objectList, textMatrix=textMatrix)
    self.waiting = 0
 
