
"""
======================COPYRIGHT/LICENSE START==========================

EditPeakAliasing.py: Part of the CcpNmr Analysis program

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

from memops.gui.ButtonList import UtilityButtonList
from memops.gui.IntEntry import IntEntry
from memops.gui.LabelDivider import LabelDivider
from memops.gui.ScrolledMatrix import ScrolledMatrix

from ccpnmr.analysis.core.AssignmentBasic import aliasedPeakDimPosition
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.PeakBasic import setPeakDimNumAliasing
from ccpnmr.analysis.core.UnitConverter import unit_converter

class EditPeakAliasingPopup(BasePopup):
  """
  **Move Aliased Peaks to Their Underlying Resonance Position**
  
  This popup window is used to move aliased 'ghost' peaks to their real
  underlying resonance locations by adding or subtracting a whole number of
  spectrum widths to the position in one or more peak dimensions. This is used
  when a resonance, that causes a peak, lies outside the normal recorded bounds
  of the spectrum but the peak nonetheless still appears within the spectrum, as
  an aliased signal that re-appears as if wrapped back onto the opposite side of
  the spectrum.

  The minimum and maximum aliased frequency values for a spectrum dimension, as
  edited in the "Referencing" table of the main Spectra_ popup, may be set
  extend the contour display beyond the normal sweep width (ppm range) of the
  spectrum and thus cover the real ppm position of peaks that have been
  unaliased. In such instances the contours are extended by tiling; one or more
  copies of the contours (not mirror image) are made and placed sequentially
  next to the normal, fundamental region. If a peak is unaliased to a position
  outside the displayed spectrum limits then the contour display will naturally
  be extended to cover the new peak position; all peaks will be visible within
  the spectrum by default. However, the user may at any time reset the minimum
  and maximum aliased frequency for a spectrum display (see the Spectra_ popup);
  deleting all values will reset bounds to the original sweep with, but any
  values may be chosen, within reason.

  A peak can be moved to the unaliased position of its underlying resonances by
  editing the "Num. Aliasing" column of this popup; double-clicking and typing in
  the appropriate number. When this number is changed for a peak dimension the peak
  will be instantly moved to its new location. An aliasing value of zero means that
  the peak lies within the sweep width of the spectrum. For a ppm scale having a
  positive aliasing value will *reduce* the ppm value, placing the peak a number of
  spectrum widths above or to the right of the spectrum bounds. Likewise a negative
  aliasing value will place a peak at a higher ppm value; below or to the left.
  Often aliasing values will be 1 or -1, where a peak has just fallen off the edge
  of a spectrum. For example a glycine amide nitrogen really at 100 ppm may be just
  outside the top of a 15N HSQC and be wrapped back into the bottom to appear as a
  peak at around 135 ppm, which means that the aliasing value should be set to 1,
  moving the peaks position up by a sweep with of 35 ppm, from 135 ppm to 100 ppm.
  The sign of the aliasing number may seem to be backwards, but it is perhaps the
  ppm scale that is 'backwards'. More complex aliasing is often seen in 3D 13C
  spectra where the 13C axis can be folded (wrapped) to reduce the sweep width that
  needs to be recorded, but still avoiding peak overlap because of the way shifts
  correlate. For example a 13C HSQC NOESY may be recorded with a sweep with that
  covers the CA, CB range 25-75 ppm but aliased methyl carbons below 25 ppm and
  aromatic carbons between 110 ppm and 140 ppm will be present; the methyls will
  have aliasing values of 1, and the aromatics -1 or -2.

  It should be noted that picking peaks in the tiled copies of a contour display,
  i.e. outside the sweep width, will automatically set the aliasing value for the
  peak to reflect the displayed chemical shift value. Thus, the user does not need
  to explicitly unalias the peak position. Any peaks that are moved by virtue of
  being unaliased will have their contribution to the chemical shifts, of any
  assigned resonances, adjusted automatically. Chemical shift values are always
  calculated using the underlying resonance positions, not the apparent peak
  position. Also, if it is known that many peaks share the same aliasing values,
  i.e. are in the same sweep width tile, then the user can propagate the aliasing
  value from one peak to many others in a single step via the right-click window
  menu; "Peak::Unaliasing propagate".

  .. _Spectra: EditSpectrumPopup.html
  
  """
  def __init__(self, parent, peak=None, *args, **kw):

    self.peak      = peak
    self.peakDim   = None
    self.guiParent = parent

    BasePopup.__init__(self, parent=parent, title="Edit Peak Aliasing", **kw)

  def body(self, guiFrame):
    
    self.geometry("500x250")

    self.numAliasingEntry = IntEntry(self,text='', returnCallback = self.setNumAliasing, width=4)
  
    guiFrame.expandGrid(1,0)
    
    
    div = LabelDivider(guiFrame, text='Peak Dimension Positions', grid=(0,0))
    utilButtons = UtilityButtonList(guiFrame, doClone=False, 
                                    closeCmd=self.close,
                                    helpUrl=self.help_url, grid=(0,1))


    tipTexts = ['The peak/spectrum dimension number',
                'The kind of isotope measured in the dimension',
                'The position of the peak in this dimension, in units of ppm',
                'The frequency position of the peak in this dimension, in units of Hz',
                'The data point position (in the spectrum matrix) of the peak in this dimension',
                'Sets the number of spectrum sweep withs to add to the peak dimension position to locate it at its real ppm value. Note an aliasing of "1" moves a peak to a lower ppm',
                'The assignment annotation for the peak dimension']
                
    headingList = ['Dimension','Isotope','ppm','Hz','Points','Num.\nAliasing','Annotation']
    editWidgets      = [None, None, None, None, None, self.numAliasingEntry, None]
    editGetCallbacks = [None, None, None, None, None, self.getNumAliasing,   None]
    editSetCallbacks = [None, None, None, None, None, self.setNumAliasing,   None]
    self.scrolledMatrix = ScrolledMatrix(guiFrame, tipTexts=tipTexts,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets,
                                         initialCols=5,
                                         initialRows=3,
                                         headingList=headingList,
                                         callback=self.selectCell,
                                         grid=(1,0), gridSpan=(1,2))


    for func in ('__init__','delete','setAnnotation',
                 'setNumAliasing','setPosition'):
      self.registerNotify(self.updateAfter, 'ccp.nmr.Nmr.PeakDim', func)
    
    self.waiting = False
    self.updateAfter(self.peak)

  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)

  def setNumAliasing(self, event):
  
    value = self.numAliasingEntry.get()
    if (value is not None) and self.peakDim:
      setPeakDimNumAliasing(self.peakDim, value)
      self.updateAfter()
  
  def getNumAliasing(self, peakDim):
  
    if peakDim:
      self.numAliasingEntry.set(peakDim.numAliasing)

  def selectCell(self, object, row, col):
  
    self.peakDim = object

  def updateAfter(self, object=None):
  
    if object:
      if object.className == 'Peak':
        self.peak = object
      elif object.peak is not self.peak:
        # object is peakDim & function was called by notifier
        # return if not my peak
        return
      
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)

  def update(self):
  
    objectList  = []
    textMatrix  = []
    colorMatrix = []
    colors = [None] * 7
    colors[5] = '#B0FFB0'
    
    if self.peak:
      for peakDim in self.peak.sortedPeakDims():
        dataDimRef = peakDim.dataDimRef
        if dataDimRef:
          objectList.append(peakDim)
    else:
      textMatrix.append([])

    for peakDim in objectList:
      dataDimRef = peakDim.dataDimRef
      expDimRef = dataDimRef.expDimRef
      position = aliasedPeakDimPosition(peakDim)
      datum = [peakDim.dim,
               '/'.join(expDimRef.isotopeCodes),
               unit_converter[('point','ppm')](position,dataDimRef),
               unit_converter[('point','Hz') ](position,dataDimRef),
               position,
               peakDim.numAliasing,
               peakDim.annotation]
 
      textMatrix.append(datum)
      colorMatrix.append(colors)

    self.scrolledMatrix.update(objectList=objectList,
                               textMatrix=textMatrix,
                               colorMatrix=colorMatrix)
    self.waiting = False
  
  def destroy(self):

    for func in ('__init__','delete','setAnnotation',
                 'setNumAliasing','setPosition'):
      self.unregisterNotify(self.updateAfter, 'ccp.nmr.Nmr.PeakDim', func)

    BasePopup.destroy(self)

