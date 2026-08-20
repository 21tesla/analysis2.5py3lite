"""
======================COPYRIGHT/LICENSE START==========================

NoeMatrix.py: Part of the CcpNmr Analysis program

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
from memops.gui.ScrolledDensityMatrix import ScrolledDensityMatrix
from ccpnmr.analysis.core.AssignmentBasic import getResidueResonances
from ccpnmr.analysis.core.MoleculeBasic import makeGuiName, getResidueCode

class NoeMatrix(ScrolledDensityMatrix):

  def __init__(self, parent, peakLists=None, restraintLists=None, *args,**kw):
  
    ScrolledDensityMatrix.__init__(self, parent, labelAxes=False)
    
    self.peakLists = peakLists or []
    self.restraintLists = restraintLists or []
    self.doAtoms   = 0
    self.displayedAtoms = ['H','HA']
    self.xResidues = []
    self.yResidues = []
    self.waiting = False 
    
  def setXrange(self, residues):
  
    self.xResidues = residues
    self.updateAfter()
  
  def setYrange(self, residues):
      
    self.yResidues = residues
    self.updateAfter()
      
  def setPeakLists(self, peakLists):
  
    self.peakLists = peakLists
    self.updateAfter()
      
  def getLabels(self, residues):
  
    texts = []
    for residue in residues:
      ccpCode = getResidueCode(residue)
      molType = residue.molType
      text = '%d %s' % (residue.seqCode, ccpCode)
         
      if self.doAtoms:
        for atom in residue.atoms:
          if atom.name in self.displayedAtoms:
            texts.append( '%s %s' % (text,makeGuiName(atom.name, atom.chemAtom.elementSymbol, molType)) )
            
      else:
        texts.append(text)

    return texts
 
  def getNumResidueContacts(self):
  
    contactDict = {}
    xResDict = {}
    yResDict = {}
    
    allResidues = self.xResidues + self.yResidues
    for residueA in allResidues:
      contactDict[residueA] = {}
      
      for residueB in allResidues:
        contactDict[residueA][residueB] = 0
  
    xPeaks = {}
    for xResidue in self.xResidues:
      xResDict[xResidue] = True
      xResonances = getResidueResonances(xResidue)
      xPeaks[xResidue] = {}
      for resonance in xResonances:
        for contrib in resonance.peakDimContribs:
          peak = contrib.peakDim.peak
          if peak.peakList in self.peakLists:
            xPeaks[xResidue][peak] = contrib.peakDim

            
    for yResidue in self.yResidues:
      yResDict[yResidue] = True
      yResonances = getResidueResonances(yResidue)
      yPeaks = {}
      for resonance in yResonances:
        for contrib in resonance.peakDimContribs:
          peak = contrib.peakDim.peak
          if peak.peakList in self.peakLists:
            yPeaks[peak] = contrib.peakDim
      
      for xResidue in self.xResidues:
        for peak in xPeaks[xResidue]:
          if yPeaks.get(peak) is not None:
            if yPeaks[peak] is not xPeaks[xResidue][peak]:
              if xResidue is not yResidue:
                contactDict[xResidue][yResidue] += 1
                contactDict[yResidue][xResidue] += 1

                
    for constraintList in self.restraintLists:
      for constraint in constraintList.constraints:
        items  = constraint.items
        weight = 1.0/len(items)
      
        for item in items:
          resonances = list(item.resonances)
          
          if None not in resonances:
            resonanceSetA = resonances[0].resonanceSet
            resonanceSetB = resonances[1].resonanceSet
            
            if resonanceSetA and resonanceSetB:
              residueA = resonanceSetA.findFirstAtomSet().findFirstAtom().residue
              residueB = resonanceSetB.findFirstAtomSet().findFirstAtom().residue
              if residueA is residueB:
                continue
              
              if xResDict.get(residueA) and yResDict.get(residueB):
                contactDict[residueA][residueB] += weight
             
              if xResDict.get(residueB) and yResDict.get(residueA):
                contactDict[residueB][residueA] += weight
         
    return contactDict
  
  def updateAssignment(self, peakDimContrib):
  
    # another peak has been assigned
    if peakDimContrib.peak.peakList in self.peakLists:
      resonanceSet = peakDimContrib.resonance.resonanceSet
      if resonanceSet:
        residue = resonanceSet.findFirstAtomSet().findFirstAtom().residue
        if (residue in self.xResidues) or (residue in self.yResidues):
          self.updateAfter()
  
  def updateAfter(self, *opt):
    
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
    
  def update(self):
   
    
    self.matrix = []
    if self.doAtoms:
      xPeaks = {}
      xAtoms = []
      for xResidue in self.xResidues:
        for xAtom in xResidue.atoms:
          if xAtom.chemAtom.elementSymbol == 'H':
            if xAtom.name in self.displayedAtoms:
              xAtomSet = xAtom.atomSet
              if xAtomSet and xAtomSet.resonanceSets:
                xAtoms.append(xAtom)
                xPeaks[xAtom] = {}
                for resonanceSet in xAtomSet.resonanceSets:
                  for resonance in resonanceSet.resonances:
                    for contrib in resonance.peakDimContribs:
                      xPeaks[xAtom][contrib.peakDim.peak] = contrib.peakDim
                         
      yPeaks = {}
      yAtoms = []
      for yResidue in self.yResidues:
        for yAtom in yResidue.atoms:
          if yAtom.chemAtom.elementSymbol == 'H':
            if yAtom.name in self.displayedAtoms:
              yAtomSet = yAtom.atomSet
              if yAtomSet and yAtomSet.resonanceSets:
                yAtoms.append(yAtom)
                yPeaks[yAtom] = {}
                for resonanceSet in yAtomSet.resonanceSets:
                  for resonance in resonanceSet.resonances:
                    for contrib in resonance.peakDimContribs:
                      yPeaks[yAtom][contrib.peakDim.peak] = contrib.peakDim
            
      for xAtom in xAtoms:
        row = []
        for yAtom in yAtoms:
          n = 0
          for yPeak in yPeaks[yAtom]:
            if xPeaks[xAtom].get(yPeak):
              if xPeaks[xAtom][yPeak] is not yPeaks[yAtom][yPeak]:
                n += 1
          row.append(n)
            
        row.reverse()
        self.matrix.append(row)
                  
    else:
      numContacts = self.getNumResidueContacts()
      blankRow = [0] * len(self.yResidues)
      for xResidue in self.xResidues:
        row = list(blankRow)
        for i, yResidue in enumerate(self.yResidues):
          row[i] = numContacts[xResidue][yResidue]
          
        row.reverse()
        self.matrix.append(row)

    if not self.matrix:
      self.matrix += [[0] * 50 for i in range(50)]  
      
      for i in range(50):
        for j in range(50):
           self.matrix[i][j] = (i/50.0) * (j/50.0)

    #self.matrix.reverse()
    self.xLabels = self.getLabels(self.xResidues)
    self.yLabels = self.getLabels(self.yResidues)
    self.yLabels.reverse()
    self.setMatrixMaxVal()
    self.draw()
    
    self.waiting = False

    
