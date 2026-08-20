
"""
======================COPYRIGHT/LICENSE START==========================

ResonanceFrame.py: Part of the CcpNmr Analysis program

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

from memops.general import Implementation

from memops.gui.Frame import Frame
from memops.gui.ScrolledMatrix import ScrolledMatrix

from ccp.api.nmr import Nmr
from ccpnmr.api import Analysis

from ccp.util.NmrExpPrototype import longRangeTransfers

from ccpnmr.analysis.core.AssignmentBasic import makeResonanceGuiName, newResonance, findMatchingPeakDimShifts, getResonanceMolSystem
from ccpnmr.analysis.core.UnitConverter import unit_converter, pnt2ppm, ppm2pnt
from ccpnmr.analysis.core.PeakBasic import findNumAliasing
from ccpnmr.analysis.core.StructureBasic import getAtomSetsDistance
from ccpnmr.analysis.core.Util import getAnalysisDataDim
from ccpnmr.analysis.core.CouplingBasic import findMatchingCouplings

TIP_TEXTS = ['Row number',
             'Resonance name',
             'Difference between peak dimension position and chemical shift',
             'Chemical shift value of assignment',
             'Standard deviation of chemical shift',
             'Resonance-resonance distance in structure']


TIP_TEXTS2 = ['Row number',
              'Names of coupled resonances',
              'Difference between peak splitting and chemical shift',
              'Coupling value of peak splitting',
              'Standard deviation of chemical shift']

class ResonanceFrame(Frame):

  def __init__(self, parent, guiParent, *args, **kw):
   
    Frame.__init__(self, parent, *args, **kw)
   
    self.guiParent = guiParent
    self.aliasing = 0
    self.peakDim = None
    self.contrib = None
    self.component = None
    self.structure = None
    self.resonances = []
    self.jCouplings = []
    
    headingList = ['#','Name','Delta','Shift','SD','Dist']
    self.scrolledMatrix = ScrolledMatrix(self, initialRows=5, headingList=headingList,
                                         tipTexts=TIP_TEXTS,
                                         callback=self.setCurrentResonance, highlightType=None)

    row = 0
    self.scrolledMatrix.grid(row=row, column=0, columnspan=2, sticky='nsew', padx=1)
    self.grid_rowconfigure(row, weight=1)

    self.grid_columnconfigure(0, weight=1)
    self.grid_columnconfigure(1, weight=1)
    
    self.update(None,None)


  def update(self, contrib, peakDim, aliasing=False, structure=None,
             limitedResonances=None, molSystems=None, component=None,
             noDisplay=False, doubleTol=False, resetScrollbars=True,
             showMolSystem=False ):
  
    self.aliasing = aliasing
    self.contrib = contrib
    self.resonances = []
    self.jCouplings = []
    self.structure = structure
    self.component = component
    atomSets1 = []
          
    if peakDim:
      self.peakDim = peakDim
    elif contrib:
      self.peakDim = None
      peakDim = contrib.peakDim
    else:
      self.scrolledMatrix.update(objectList=[], textMatrix=[[],])
      return

    if component and (component.dataDimRef.expDimRef.unit != 'ppm'):
      headingList = ['#','Names','Delta','Coupling','SD']
      tipTexts = TIP_TEXTS2
      isCoupling = True
    else:
      tipTexts = TIP_TEXTS
      headingList = ['#','Name','Delta','Shift','SD','Dist']
      isCoupling = False

    if showMolSystem:
      tipTexts = tipTexts[:]
      tipText = 'Mol System of resonance'
      if isCoupling:
        tipText += 's'
      tipTexts.append(tipText)
      headingList.append('MS')
      
    dataDimRef = peakDim.dataDimRef
    expDimRef = dataDimRef.expDimRef
    dataDim = dataDimRef.dataDim
    textMatrix  = []
    colorMatrix = []
    
    if isCoupling:

      ppm = getAnalysisDataDim(peakDim.dataDim).assignTolerance
      points = peakDim.dataDimRef.valueToPoint(ppm)
      tolerance = component.dataDimRef.pointToValue(points)

      if doubleTol:
        tolerance *= 2
      
      # get couplings that might match the splitting
      for delta, coupling in findMatchingCouplings(component, tolerance):
        resonanceA, resonanceB = coupling.resonances
        nameA = makeResonanceGuiName(resonanceA)
        nameB = makeResonanceGuiName(resonanceB)
        
        self.jCouplings.append(coupling)

        colors = [None,None,None,None,None]
        if coupling.error >= tolerance:
          colors[4] = '#d0a0a0'
 
        textMatrix.append( ['%d,%d' % (resonanceA.serial, resonanceB.serial),
                            '%s-%s' % (nameA, nameB),
                            '%5.3f' % round(delta, 3),
                            coupling.value,
                            '%5.3f' % round(coupling.error, 3)] )
        colorMatrix.append( colors )

        if showMolSystem:
          molSystemA = getResonanceMolSystem(resonanceA) or ''
          molSystemB = getResonanceMolSystem(resonanceB) or ''
          if molSystemA == molSystemB:
            texts.append(molSystemA and molSystemA.code)
          else:
            texts.append('%s,%s' % (molSystemA and molSystemA.code, molSystemB and molSystemB.code))
          colors.append(None)

      
      objectList = self.jCouplings
      
    else:
      # set up atomSets for distance calculations where possible
      if self.structure:
        for expTransfer in expDimRef.expTransfers:
          if expTransfer.transferType in longRangeTransfers:
            expDimRefs = list(expTransfer.expDimRefs)
            expDimRefs.remove(expDimRef)
 
            for peakDim1 in peakDim.peak.sortedPeakDims():
              if peakDim1.dataDimRef.expDimRef is expDimRefs[0]:
                for contrib1 in peakDim1.peakDimContribs:
                  resonance1 = contrib1.resonance
 
                  if resonance1 and resonance1.resonanceSet:
                    for atomSet in resonance1.resonanceSet.atomSets:
                      atomSets1.append(atomSet)
      
      if component:
        # E.g. MQ or reduced dimensionality
        dataDimRef2 = component.dataDimRef
        
        if dataDimRef is dataDimRef2:
          ppm = abs(peakDim.value - peakDim.realValue)
        
        else:
          realPoint = ppm2pnt(peakDim.realValue, dataDimRef)
          deltaPoints = abs(peakDim.position-realPoint)
          ppm = abs(component.scalingFactor*dataDimRef2.pointToValue(deltaPoints))
 
        peakUnit = dataDimRef2.expDimRef.unit or 'point'
      
      else:
        ppm = peakDim.realValue
        peakUnit = expDimRef.unit or 'point'
        
      shiftList = peakDim.peak.peakList.dataSource.experiment.shiftList
      shiftUnit = shiftList.unit
  
      tolerance = getAnalysisDataDim(dataDim).assignTolerance
      if doubleTol:
        tolerance *= 2
        
        
      shifts = self.givePossAssignments(peakDim, tolerance, ppm)
      
      if len(shifts)>0: # resonance matches
        temp = []
        for shift in shifts:
          resonance = shift.resonance
          if limitedResonances is not None:
            if resonance not in limitedResonances:
              continue
 
          if molSystems:
            molSystem = getResonanceMolSystem(resonance)
            if molSystem and (molSystem not in molSystems):
              continue
            
            # Kludge until chain states
            if (not molSystem) and resonance.resonanceGroup and resonance.resonanceGroup.chains:
              molSystems2 = set([ch.molSystem for ch in resonance.resonanceGroup.chains])
              
              if not molSystems2.intersection(molSystems):
                continue
 
          if peakUnit != shiftUnit:
            shiftValue = unit_converter[(shiftUnit,peakUnit)](shift.value,dataDimRef)
          else:
            shiftValue = shift.value
 
          points = unit_converter[(peakUnit,'point')](shiftValue, dataDimRef)
 
          numAliasing = findNumAliasing(dataDimRef, points)
          deltaAliasing = numAliasing - peakDim.numAliasing
 
          if deltaAliasing:
            points -= dataDim.numPointsOrig * deltaAliasing
            shiftValue = unit_converter[('point',peakUnit)](points, dataDimRef)
            shift.isDiffAliasing = 1
          else:
            shift.isDiffAliasing = 0
 
          delta = abs(ppm-shiftValue)
          temp.append( (delta, shift) )
          shift.delta = delta
 
        temp.sort()
        shifts = [ x[1] for x in temp ]
 
        textMatrix  = []
        colorMatrix = []
        for shift in shifts:
          resonance = shift.resonance
          if limitedResonances is not None:
            if resonance not in limitedResonances:
              continue
 
          # wb104, 13 Apr 2016, added below (could check this instead in if statement
          # but the problem is that shift.delta also needs to exist further down
          # and also both attributes get deleted at bottom of the loop)

          if not hasattr(shift, 'isDiffAliasing'):
            continue

          if shift.isDiffAliasing:
            colors = [None,None,'#d0d0a0',None,None,None]
          else:
            colors = [None,None,None,None,None,None]
 
          name = makeResonanceGuiName(resonance)
          self.resonances.append(resonance)

          if shift.error >= tolerance:
            colors[4] = '#d0a0a0'
 
          if self.structure and atomSets1:
            if resonance.resonanceSet:
              atomSets2 = list(resonance.resonanceSet.atomSets)
              dist = getAtomSetsDistance(atomSets1, atomSets2,
                                         self.structure, method='noe')
            else:
              dist = None
          else:
            dist = None
 
          texts = [resonance.serial,
                   name,
                   '%5.3f' % round(shift.delta, 3),
                   shift.value,
                   '%5.3f' % round(shift.error, 3),
                   dist or ' '] 

          if showMolSystem:
            molSystem = getResonanceMolSystem(resonance)
            texts.append(molSystem and molSystem.code)
            colors.append(None)

          textMatrix.append(texts)
          colorMatrix.append( colors )
 
          del shift.isDiffAliasing
          del shift.delta
      
      objectList = self.resonances
    
    
    if not noDisplay:            
      self.scrolledMatrix.update(headingList=headingList,
                                 tipTexts=tipTexts,
                                 objectList=objectList,
                                 textMatrix=textMatrix,
                                 colorMatrix=colorMatrix,
                                 resetScrollbars=resetScrollbars)
    
  def givePossAssignments(self, peakDim, tolerance, ppm):

   return findMatchingPeakDimShifts(peakDim,aliasing=self.aliasing,
                                    tolerance=tolerance,ppm=ppm)
                       
  def setCurrentResonance(self,object, row, col):
  
    # will leave atom set mappings in here for suggestions based upon chemical shift
    if object is None:
      project = self.peakDim.root
      object  = newResonance(project)
      
    self.guiParent.peakDim = self.peakDim
    self.guiParent.assign(object, self.component)
    self.scrolledMatrix.currentObject = None
       
    
    
    
    
    
    
    

    
    
    
