
"""
======================COPYRIGHT/LICENSE START==========================

EditMarks.py: Part of the CcpNmr Analysis program

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
from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Label import Label
from memops.gui.LabelDivider import LabelDivider
from memops.gui.PulldownList import PulldownList
from memops.gui import Color

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.MarkBasic import setNonPeakMarkColor, setRulerColor
from ccpnmr.analysis.core.Util import  getHueSortedColorSchemes

class EditMarksPopup(BasePopup):
  
  """
  **Setup Spectrum Marker and Ruler Line Options**
  
  This popup window is used to specify how many marker and ruler lines may
  be added by the user to spectrum windows, and what colours these will be.
  
  A marker line is a multi-dimensional cross mark that is typically added to
  the spectrum display with the "m" key, using the mouse cursor location as
  the center point. A mark of this kind gives a cross that is displayed on a
  spectrum window (in the plane of the screen), although the same lines will
  be visible in other spectra that show the relevant frequency/ppm range. A
  mark will also place lines in any orthogonal (z) positions of the spectrum
  view; these will not be visible in the window to which the mark was first
  added but will be seen in any windows that show orthogonal view views of
  the marked point.

  A ruler is a single line that is added to the spectrum window displays,
  typically using the "h" key for a horizontal line and the "v" key for a
  vertical line. A ruler only places a single line in the plane of the
  screen at the mouse cursor location. A ruler may be visible in more than
  one spectrum window, but only if the view is in the right range and for
  axes that share the same "Panel Type"; a sub-division of axis types that
  are otherwise specified in terms of isotope. The panel type for a window
  axis is set via the main Windows_ table.

  This system allows the user to specify the maximum number of marks and
  rulers that may be displayed on the screen at any one time. It should be
  noted that this setting only applies to marks and rulers that are directly
  added by the user, others may be added by various tools throughout
  Analysis. When adding marks or rulers to a spectrum window, if the limit
  is reached then adding another will remove the oldest mark/ruler.
  Naturally, if the limit is set to one, then the mark or ruler is always
  replaced. By decreasing the maximum limit the oldest will be removed,
  although the usual way of removing (all) marks and rulers is via the "n"
  key after clicking in a spectrum window.

  Marks and rulers that are added by the user may use a colour scheme, such
  that each subsequent addition uses the next colour from the scheme.

  .. _Windows: EditWindowPopup.html
  
  """
  
  maxAllowedMarks = 30
  maxAllowedRulers = 30

  def __init__(self, parent, *args, **kw):

    BasePopup.__init__(self, parent=parent, title='Window : Marks and Rulers', **kw)

  def body(self, guiParent):

    self.geometry('320x120')
    analysisProject = self.analysisProject
    guiParent.expandGrid(5,3)
    guiParent.expandGrid(5,4)

    row = 0
    div = LabelDivider(guiParent, text='Multi-dimensional Marks',
                       grid=(row,0), gridSpan=(1,4))
 
    buttons = UtilityButtonList(guiParent, doClone=False, helpUrl=self.help_url)
    buttons.grid(row=row, column=4, sticky='w')
    
    row += 1

    label = Label(guiParent, text='Maximum marks: ', grid=(row,0))
 
    values = range(1, self.maxAllowedMarks+1)
    entries = [ str(x) for x in  values]
    value = analysisProject.maxMarks
    tipText = 'Sets the maximum number of multi-dimensional cross-marks that the user can add to spectrum window displays'
    self.marks_menu = PulldownList(guiParent, callback=self.setMaxMarks, grid=(row,1),
                                   texts=entries, objects=values,
                                   index=value-1, tipText=tipText)

    label = Label(guiParent, text=' Mark colour: ', grid=(row,2))
 
    tipText = 'Sets the line colour of the multi-dimensional cross-marks, excepting those that go through peak centers'
    self.marksSchemePulldown = PulldownList(guiParent, callback=self.setMarksColor,
                                            grid=(row,3), tipText=tipText)

    row += 1
    label = Label(guiParent, text='(Non-peak marks only)', grid=(row,2), gridSpan=(1,2))
 
    row += 1
    div = LabelDivider(guiParent, text='1-dimensional Rulers',
                       grid=(row,0), gridSpan=(1,5))
    
    row += 1
    label = Label(guiParent, text='Maximum rulers: ', grid=(row,0))
 
    values = range(1, self.maxAllowedRulers+1)
    entries = [ str(x) for x in  values]
    value = analysisProject.maxRulers
    tipText = 'Sets the maximum number of 1-dimensional ruler lines that the user can add to spectrum window displays'
    self.rulers_menu = PulldownList(guiParent, callback=self.setMaxRulers, grid=(row,1),
                                    texts=entries, objects=values, index=value-1, tipText=tipText)
 
    label = Label(guiParent, text=' Ruler colour: ', grid=(row,2))
 
    tipText = 'Sets the line colour of the 1-dimensional ruler lines'
    self.rulersSchemePulldown = PulldownList(guiParent, callback=self.setRulersColor,
                                             grid=(row,3), tipText=tipText)

    self.updateSchemePulldowns()

  def updateSchemePulldowns(self):

    profile = self.analysisProfile
    
    indexR = 0
    indexM = 0

    schemes = getHueSortedColorSchemes(profile)
    colors = [list(s.colors) for s in schemes]
    names = [s.name for s in schemes]
    
    if profile.marksColor not in schemes:
      profile.marksColor = schemes[0]
    
    if profile.rulersColor not in schemes:
      profile.rulersColor = schemes[0]
    
    indexM = schemes.index(profile.marksColor)
    indexR = schemes.index(profile.rulersColor)
    
    self.marksSchemePulldown.setup(names, schemes, indexM, colors=colors)
    self.rulersSchemePulldown.setup(names, schemes, indexR, colors=colors)

  def setMaxMarks(self, value):
 
    self.analysisProject.maxMarks = value
 
  def setMarksColor(self, scheme):

    profile = self.analysisProfile
    profile.marksColor = scheme

    setNonPeakMarkColor(self.project)
 
  def setMaxRulers(self, value):
 
    self.analysisProject.maxRulers = value
 
  def setRulersColor(self, scheme):

    profile = self.analysisProfile
    profile.rulersColor = scheme

    setRulerColor(self.project)
