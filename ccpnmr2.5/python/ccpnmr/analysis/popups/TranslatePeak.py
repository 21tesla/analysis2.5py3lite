
"""
======================COPYRIGHT/LICENSE START==========================

TranslatePeak.py: Part of the CcpNmr Analysis program

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
from memops.general import Implementation

from memops.gui.ButtonList     import UtilityButtonList
from memops.gui.CheckButton    import CheckButton
from memops.gui.Label          import Label
from memops.gui.LabelFrame     import LabelFrame
from memops.gui.PulldownList   import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.PeakBasic import translateSpectrumUsingPeaks

class TranslatePeakPopup(BasePopup):

  def __init__(self, parent, fixedPeak, movePeak, *args, **kw):

    self.movePeakDim = None
    self.fixedPeak = fixedPeak
    self.movePeak = movePeak
    self.setupDimMapping()

    BasePopup.__init__(self, parent=parent, title='Reference By Peak',
                       modal=True, *args, **kw)

  def body(self, guiFrame):

    fixedSpectrum = self.fixedPeak.peakList.dataSource
    moveSpectrum  = self.movePeak.peakList.dataSource
    
    guiFrame.grid_columnconfigure(0, weight=1)

    row = 0
    frame = LabelFrame(guiFrame, text='Working Spectra',
                       grid=(row,0), sticky='ew')
    frame.grid_columnconfigure(1, weight=1)
    
    s = '%s:%-20s' % (fixedSpectrum.experiment.name, fixedSpectrum.name)
    
    Label(frame, text='Fixed:', grid=(0,0), sticky='e')
    tipText = 'The experiment:spectrum name of the reference spectrum, which will not be moved'
    Label(frame, text=s, grid=(0,1), tipText=tipText)

    s = '%s:%-20s' % (moveSpectrum.experiment.name, moveSpectrum.name)

    Label(frame, text='Moving:', grid=(0,2), sticky='e')
    tipText = 'The experiment:spectrum name of the spectrum which will be re-referenced using peak positions'
    Label(frame, text=s, grid=(0,3), tipText=tipText)
    

    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)
    frame = LabelFrame(guiFrame, text='Peak Position Mapping', grid=(row,0))
    frame.expandGrid(0,0)
    
    tipText = 'Whether the more than one moving spectrum dimension may be mapped to the same fixed spectrum dimension; typically not used but can be useful for re-referencing with a diagonal peak'
    self.fixedDimsButton = CheckButton(frame, text='Allow repeated fixed dim',
                                       selected=False, grid=(0,0),
                                       callback=self.changedFixedDimsButton,
                                       tipText=tipText)

    self.dimPulldown = PulldownList(self, callback=self.setDimMapping)
    
    editWidgets      = [None, None, self.dimPulldown,   None, None, None]
    editGetCallbacks = [None, None, self.getDimMapping, None, None, None]
    editSetCallbacks = [None, None, self.setDimMapping, None, None, None]
     
    tipTexts = ['The dimension number of the spectrum being re-referenced',
                'The isotope type that corresponds to both the fixed and re-reference spectrum dimensions',
                'The dimension number of the fixed spectrum, used as a reference',
                'The position of the moving peak in the relevant dimension, before any changes (for the spectrum being re-referenced)',
                'The position of the fixed peak in the relevant dimension; potential new position of the moving peak after re-referencing changes',
                'The chemical shift difference between the two peak positions on one dimension']
    headingList = ['Moving\nDim','Isotope','Fixed\nDim','Original\nPosition',
                   'New\nPosition',u'\u0394']
    self.scrolledMatrix = ScrolledMatrix(frame, headingList=headingList, 
                                         callback=self.selectPeakDim,
                                         editWidgets=editWidgets, maxRows=5,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         grid=(1,0), tipTexts=tipTexts)

    row += 1
    tipTexts = ['Use the stated mapping of dimensions to re-reference the moving spectrum so that the selected peaks lie at the same position',]
    texts        = ['Commit']
    commands     = [self.ok]
    buttons = UtilityButtonList(guiFrame, texts=texts, commands=commands, 
                                closeText='Cancel', helpUrl=self.help_url,
                                grid=(row,0), tipTexts=tipTexts)

    self.update()

  def getDataDimIsotopes(self, dataDim):
  
    isotopes = []
    
    if dataDim.className == 'FreqDataDim':
      for dataDimRef in dataDim.dataDimRefs:
        for isotopeCode in dataDimRef.expDimRef.isotopeCodes:
          isotopes.append(isotopeCode)
        
    return isotopes    

  def setupDimMapping(self):

    self.dimMapping = {}
    usedDataDims = set()
    for movePeakDim in self.movePeak.sortedPeakDims():
      fixedDataDims = self.getMatchingDataDims(movePeakDim)
      fixedDataDims = [dataDim for dataDim in fixedDataDims if dataDim not in usedDataDims]
      if fixedDataDims:
        moveDataDim = movePeakDim.dataDim
        fixedDataDim = fixedDataDims[0]
        self.dimMapping[moveDataDim.dim] = fixedDataDim.dim
        usedDataDims.add(fixedDataDim)

  def getMatchingDataDims(self, movePeakDim):

    moveIsotopes = self.getDataDimIsotopes(movePeakDim.dataDim)
    fixedDataDims = []

    for fixedPeakDim in self.fixedPeak.sortedPeakDims():
      fixedDataDim  = fixedPeakDim.dataDim
      isotopes = self.getDataDimIsotopes(fixedDataDim)

      if isotopes == moveIsotopes:
        fixedDataDims.append(fixedDataDim)

    return fixedDataDims

  def getDimMapping(self, movePeakDim):
    
    fixedDataDims = [None,]
    fixedDataDims.extend(self.getMatchingDataDims(movePeakDim))
    possibleDims = [('%d' % fixedDataDim.dim if fixedDataDim else 'None') for fixedDataDim in fixedDataDims]
    
    mappedDim = self.dimMapping.get(movePeakDim.dim)
    if mappedDim:
      try:
        # SF bug 3137169 had problem with string not being in list
        # not sure why that might be happening but play safe
        index = possibleDims.index(str(mappedDim))
      except:
        index = 0
    else:
      index = 0

    self.dimPulldown.setup(possibleDims, fixedDataDims, index)
    
  def setDimMapping(self, fixedDataDim):
    
    if fixedDataDim and fixedDataDim.className == 'PeakDim':
      # some kind of bug in ScrolledMatrix which means that this being called with table object
      # instead of menu object when you click outside of menu to stop editing
      return

    # fixedDataDim is the new suggested dataDim which self.movePeakDim.dataDim gets mapped to
    movePeakDim = self.movePeakDim
    moveDim = movePeakDim.dim
    if fixedDataDim is None:
      fixedDimNew = None
    else:
      fixedDimNew = fixedDataDim.dim
      fixedDimOld = self.dimMapping.get(moveDim)

    if fixedDimNew:
      if not self.fixedDimsButton.get():
        # do not allow repeated fixed dims, so need to swap things around
        # old: movePeakDim --> fixedDimOld, ??? --> fixedDimNew
        # new: movePeakDim --> fixedDimNew, ??? --> fixedDimOld
        # do not have to do anything if fixedDimNew does not appear (i.e. ??? non-existent)
 
        for moveDim2 in self.dimMapping:
          if self.dimMapping[moveDim2] == fixedDimNew:
            self.dimMapping[moveDim2] = fixedDimOld
            break
      self.dimMapping[moveDim] = fixedDimNew
 
    elif moveDim in self.dimMapping:
      del self.dimMapping[moveDim]
    
    self.update()
        
  def changedFixedDimsButton(self, selected):

    if not selected:
      keys = self.dimMapping.keys()
      keys.sort()
      usedSet = set()
      for moveDim in keys:
        fixedDim = self.dimMapping[moveDim]
        if fixedDim in usedSet:
          del self.dimMapping[moveDim]
        else:
          usedSet.add(fixedDim)
      
      self.update()

  def selectPeakDim(self, obj, row, col):
  
    self.movePeakDim = obj

  def update(self):
  
    textMatrix = []
    objectList = []
    
    for peakDim in self.movePeak.sortedPeakDims():
      dataDim = peakDim.dataDim
      
      if dataDim.className != 'FreqDataDim':
        continue
    
      dim         = peakDim.dim
      isotopes    = self.getDataDimIsotopes(dataDim)
      selectedDim = self.dimMapping.get(dim)
 
      refValue = delta = None
      if selectedDim:
        peakDim0 = self.fixedPeak.findFirstPeakDim(dim=selectedDim)
        
        if peakDim0:
          refValue = '%4.3f' %  peakDim0.value
          delta    = '%4.3f' % (peakDim0.value - peakDim.value)
      else:
        selectedDim = 'None'

      datum = [dim,','.join(isotopes),selectedDim,
               peakDim.value,refValue,delta]
               
      textMatrix.append(datum)
      objectList.append(peakDim)
  
    self.scrolledMatrix.update(objectList=objectList, textMatrix=textMatrix)

  def apply(self):

    translateSpectrumUsingPeaks(self.fixedPeak, self.movePeak, self.dimMapping)

    self.update()

    return True
