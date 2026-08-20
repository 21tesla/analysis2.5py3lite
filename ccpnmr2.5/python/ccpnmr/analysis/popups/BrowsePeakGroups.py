
"""
======================COPYRIGHT/LICENSE START==========================

BrowsePeakGroups.py: Part of the CcpNmr Analysis program

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
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showOkCancel
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.MarkBasic import createRuler
from ccpnmr.analysis.core.ExperimentBasic import getSpectrumIsotopes
from ccpnmr.analysis.core.WindowBasic import findOrthogonalWindows, isSpectrumInWindowPane, displayStrips, getPeakDimAxisMapping

axisNames = ['x','y','z1','z2','z3','z4']

class BrowsePeakGroupsPopup(BasePopup):
  """
  **Display Groups of Peaks Linked by a Root Assignment**
  
  This popup window is used to display the results of operations in Analysis that
  generate grouped peaks, usually linked by virtue of a common 'root' assignment.
  For example a selected set of peaks may be macthed to other peaks with
  similar positions in a spectrum dimension, where each group corresponds
  to a separate C-H or N-H position. For example the right-click
  window menu option "Locate Peaks::Match multiple peaks" generates
  such groupings to find rows or columns of peaks that share chemical shifts.
  
  Whatever operation generated the groups, they are listed in ranked order in the
  table; best scoring at the top. The number of peaks in the group and any common
  assigment and position is also listed. The idea is that the groups are candidates
  for some searh or macthing process that the user triggered. The best scoring
  groups may be used to control the spectrum window display, thus displaying any
  peak matches, by using [Display Groups in Strips] after selecting  an approprate
  window. This system is often used for NOESY experiments where resonances that are
  close in space share common sets of peak positions; they are close to the same
  set of other resonances.

  """
  
  def __init__(self, parent, *args, **kw):
    
    self.groups         = []
    self.variableDim    = 1
    self.windowPane     = None
    self.rulers         = []
    self.peakLists      = []
    self.orthogonalDict = {}
    self.parent         = parent
    self.searchPeaks    = []

    BasePopup.__init__(self, parent, title= "Peak Groups", borderwidth=5, **kw)
    
  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)
  
    row = 0 
    frame = LabelFrame(guiFrame, text='Matched Peak Groups')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)
    frame.grid(row=row, sticky='nsew')
    
    headingList, tipTexts = self.getHeadingList()
    self.scrolledMatrix = ScrolledMatrix(frame, initialRows=15, multiSelect=True, 
                                        headingList=headingList, tipTexts=tipTexts,
                                        grid=(0,0), gridSpan=(1,3))
 
    tipTexts = ['Remove the selected peak groups from the table',
                'Show 1D positional ruler lines for the selected groups',
                'Remove any ruler lines previously added for peak group',
                'Display selected peak groups within strips of selected window']
    texts = ['Remove\nGroups','Show\nRulers',
             'Delete\nRulers','Display Groups\nIn Strips']
    commands = [self.removeGroups,self.addRulers,
                self.removeRulers,self.stripGroups]
    self.buttons = ButtonList(frame, texts=texts, commands=commands, tipTexts=tipTexts)
    self.buttons.grid(row=1, column=0, sticky='ew')
    
    tipText = 'Selects the spectrum window in which to display strips & ruler lines'
    label = Label(frame, text='Target window:', grid=(1,1))
    self.windowPulldown = PulldownList(frame, callback=None,
                                       grid=(1,2), tipText=tipText)
    
    row+= 1
    bottomButtons = UtilityButtonList(guiFrame, helpUrl=self.help_url, expands=True, doClone=False)
    bottomButtons.grid(row = row, sticky = 'ew')

    self.update()

    for func in ('__init__', 'delete', 'setName'):
      self.registerNotify(self.updateWindowsAfter, 'ccpnmr.Analysis.SpectrumWindow', func)

  def open(self):
  
    self.update()
    BasePopup.open(self)
    
  def updateWindowsAfter(self, opt=None):
 
    self.after_idle(self.updateWindows)

  def updateWindows(self):
    
    
    if not self.groups:
      self.windowPulldown.clear()
      return

    selected = self.windowPulldown.getText()
    groups   = self.scrolledMatrix.currentObjects or self.groups
    self.orthogonalDict = {}
    windowPane = self.windowPane
    if windowPane and groups:

      meanVarPos = 0.0
      for peak in groups[0][2]:
        meanVarPos += peak.sortedPeakDims()[self.variableDim].value
      meanVarPos /= float(len(groups[0][2]))
     
      xyz = []
      for score, coords0, peaks in groups:
      
        peak       = peaks[0]
        peakDims   = peak.sortedPeakDims()
        dimMapping = getPeakDimAxisMapping(peak, self.windowPane)

        coords     = list(coords0)
        coords.insert(self.variableDim, meanVarPos)
        
        posDict = {}
        for i, peakDim in enumerate(peakDims):
          posDict[peakDim] = coords[i] 

        pos = []
        for axisPanel in windowPane.sortedAxisPanels():
          peakDim = dimMapping[axisPanel.label]
          pos.append( posDict[peakDim] )
          
        xyz.append(pos)
      
      windowZPlanes = findOrthogonalWindows(windowPane, xyz,
                                            excludeQuery=False)
      for windowPosition in windowZPlanes:
        if windowPosition:
          key, pane, xyz = windowPosition
          self.orthogonalDict[key] = (pane, xyz)
 
    xyPlaneKeys = self.orthogonalDict.keys()
    xyPlaneKeys.sort()
    if xyPlaneKeys:
      index = min(self.windowPulldown.index, len(xyPlaneKeys)-1)

      if selected in xyPlaneKeys:
        index = xyPlaneKeys.index(selected)
    
      self.windowPulldown.setup(xyPlaneKeys, xyPlaneKeys, index)

     
  def removeGroups(self):
  
    groups = self.scrolledMatrix.currentObjects
    
    if groups:
            
      newGroups = []
      for group in self.groups:
        if group not in groups:
          newGroups.append(group)
      
      self.update(groups=newGroups)
  
  def stripGroups(self):

    self.updateWindows()
    name = self.windowPulldown.getText()
    groups = self.scrolledMatrix.currentObjects
    
    if groups and name and (name != '<None>'):
      selected = self.windowPulldown.getText()
      windowPane, positions = self.orthogonalDict[selected]
      window = windowPane.spectrumWindow

      N = len(positions)
      msg = '%d positions selected. Really make %d strips in window %s'
      if N > 10 and not showOkCancel('Warning', msg % (N,N,window.name ), parent=self):
        return

      displayStrips(self.parent, positions, orthoPositions=None,
                    spectrum=None, windowPane=windowPane)
      # often fails first time...
      self.update_idletasks()
      displayStrips(self.parent, positions, orthoPositions=None,
                    spectrum=None, windowPane=windowPane)

  
  def removeRulers(self):
  
    for ruler in self.rulers:
      if not ruler.isDeleted:
        ruler.delete()
    
    self.rulers = []
  
  def addRulers(self):
     
    groups = self.scrolledMatrix.currentObjects
    if groups and self.windowPane:
      
      positions  = []
      panelTypes = []
      for peak in self.searchPeaks:
        peakDim    = peak.sortedPeakDims()[self.variableDim]
        dimMapping = getPeakDimAxisMapping(peak, self.windowPane)
        for axisPanel in self.windowPane.axisPanels:
          if dimMapping[axisPanel.label] is peakDim:
            positions.append(peakDim.value)
            panelTypes.append(axisPanel.panelType)
            break

      color = '#808080'
      for i in range(len(positions)):
        ruler = createRuler(positions[i], panelTypes[i],
                            dashLength=3, gapLength=1,
                            color=color, remove=False)
        self.rulers.append(ruler)


  def update(self, groups=None, variableDim=None, 
             peakLists=None, searchPeaks=None):
  
    self.removeRulers()

    if groups is not None:
      self.groups = groups
    
    if variableDim is not None:
      self.variableDim = variableDim

    if peakLists:
      self.peakLists = peakLists
      spectrum = self.peakLists[0].dataSource
      for spectrumWindow in self.analysisProject.spectrumWindows:
        for windowPane in spectrumWindow.spectrumWindowPanes:
          if isSpectrumInWindowPane(windowPane, spectrum):
            if len(windowPane.axisPanels) == len(spectrum.dataDims):
              self.windowPane = windowPane
              break

    self.searchPeaks = searchPeaks or []

    objectList = []
    textMatrix = []
    for i, group in enumerate(self.groups):
      
      (score,positions,peaks) = group
      datum = ['%d' % (i+1),
               '%f' % score,
               '%d' % len(peaks)]
      
      annotations = []
      for j in range(len(peaks[0].peakDims)):
        if j != self.variableDim:
          annotation = ''
          
          for peak in peaks:
            peakDims = peak.sortedPeakDims()
            peakDim = peakDims[j]
            
            if peakDim.annotation:
              if annotation and (annotation != peakDim.annotation):
                annotation = '*'
              else:
                annotation = peakDim.annotation
          annotations.append(annotation)
     
      datum.append(' '.join(annotations))
      
      for position in positions:
        datum.append('%f' % position)
     
      objectList.append(group)
      textMatrix.append(datum)

    headingList, tipTexts = self.getHeadingList()
    self.scrolledMatrix.update(textMatrix=textMatrix,
                               objectList=objectList,
                               headingList=headingList,
                               tipTexts=tipTexts)
    self.updateWindows()

  def getHeadingList(self):
  
    tipTexts = ['The ranking of the peak group, comparing its score with others',
                'The match score of the peak group',
                'Number of peaks within the matched peak group',
                'The root (usually amide) assignment, common to peaks in a group']

    headingList = ['Rank','Score','Num\npeaks','Root\nAssignment']
    
    if self.groups:
      n = len( self.groups[0][2][0].peakDims )
      
      for i in range(n):
        if i != self.variableDim:
          headingList.append( 'F%d position' % (i+1) )
          tipTexts.append('Average location of the group in the F%d dimension' % (i+1))
          
    else:
      headingList.extend(['F1 position','F2 position'])
      tipTexts.append('Average location of the group in the F1 dimension')
      tipTexts.append('Average location of the group in the F2 dimension')
    
    return headingList, tipTexts
    
  def destroy(self):
  
    for func in ('__init__', 'delete', 'setName'):
      self.unregisterNotify(self.updateWindowsAfter, 'ccpnmr.Analysis.SpectrumWindow', func)
    
    BasePopup.destroy(self)
