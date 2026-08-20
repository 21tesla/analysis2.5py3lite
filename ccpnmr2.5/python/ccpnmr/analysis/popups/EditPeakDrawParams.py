
"""
======================COPYRIGHT/LICENSE START==========================

EditPeakDrawParams.py: Part of the CcpNmr Analysis program

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


from memops.gui.Button import Button
from memops.gui.ButtonList import UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.LabeledEntry import LabeledEntry
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showError
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.PeakBasic import refreshPeakAnnotations
from ccpnmr.analysis.core.PeakDrawUtil import peak_draw_methods, getCPeakDrawMethod

class EditPeakDrawParamsPopup(BasePopup):
  """
  **Control How Peak Symbols and Annotations Are Drawn**
  
  This popup window controls the display of picked peak positions within
  spectrum windows and also how their annotations are displayed in
  many of the tables throughout the rest of Analysis.

  **Peak Annotations**

  The "Annotation Style" tab, as the name suggests, controls which bits of
  information go to form the assignment annotation that is used to label peaks
  in spectra and when listing their identities in tables. The user can toggle
  various options on or off to dictate what should be included. With no options
  selected at all peaks will show only resonance number information like
  "[123],[45]", irrespective of any atomic assignments. Switching the atom
  assignments on will substitute the atom and residue name for  the resonance
  names. The spin system toggle is used to control the display of spin system
  numbers like the curly bracketed numbers of "{89}[123],[45]" or a residue name
  if the spin system is assigned, although the atom assignment option will also
  give any residue name. The chain and molecular system toggles are used when
  there are several sequences in the assignment that need to be distinguished, so
  that you can identify "A25LsyH,N" separately to "B25LysH,N" or even
  "MS2:A25LysH,N". The details toggle allows the user to contribute to peak
  annotations; such details are set in the relevant column of the peak tables
  (e.g. `Selected Peaks`_)or the "Peak::Set details" option of the right-click
  window menu. The last "Minimal Annotations" option is special because it gives
  a quick replacement for spectrum annotations that overrides all other settings,
  for example "25LysH,N" becomes "25K", which is useful to reduce clutter when
  printing out an HSQC.

  It should be noted that, with the exception of the "Minimal Annotations"
  option, any changes to the annotation options for existing peak assignments
  will not be visible until [Update Full Annotations] is pressed, which can take
  a while to update if the number of peaks in the project is large. Any new
  assignments will observe the current settings.

  **Draw Size**
  
  The second tab controls how the peak cross symbols, which indicate the pick
  position, are drawn in the spectrum display. Here, the user has four choices:
  the default is to use a fixed number of screen pixels, so that all symbols are
  the same size, irrespective of zoom level; uniform size in ppm (with the 
  appropriate values for each isotope set in the lower table that appears)
  allows the peak symbol maintain a constant size relative to the contours, and
  thus differ in on-screen size according to the contour zoom level; and two
  choices for a "scaled" size where the peak symbol is larger or smaller
  depending on the intensity of the peak, either relative to its peak list or an
  absolute globally applied value.

  **Figure of Merit**

  The last "Merit Symbols" tab allows the user to add textual markers to a peak
  depending on the figure-of-merit value, which is an indication of reliability
  (although the precise context is at the discretion of the user). The merit
  value for a peak is set via the peak tables or the "Peak::Set merit" option of
  the right-click window menu. Aside from display purposes the merit value is
  used at a few points throughout Analysis to optionally exclude peaks. Typically
  a peak with a merit value of 0.0 may be ignored during certain analyses, e.g.
  for making NOE derived distance restraints.

  **Caveats & Tips**

  When assigning peaks via the `Assignment Panel`_ the resonance annotation will
  always display spin system information, given that this grouping is critical
  to assignment, but the chain and atom assignment will reflect the peak
  display, and thus these details may be missing if unset in this Draw
  Parameters popup.


  .. _`Assignment Panel`: EditAssignmentPopup.html
  .. _`Selected Peaks`: SelectedPeaksPopup.html
  """
  def __init__(self, parent, *args, **kw):
 
    BasePopup.__init__(self, parent=parent, title='Peak : Draw Parameters', **kw)

  def body(self, guiParent):
 
    self.geometry('380x250')

    analysisProject = self.analysisProject

    self.smallFont    = '-schumacher-clean-medium-r-normal--14-140-75-75-c-70-iso646.1991-irv'
    guiParent.grid_columnconfigure(0, weight=1)
    guiParent.grid_rowconfigure(0, weight=1)
    
    row = 0
    
    tipTexts = ['Options specifying how peak annotations are drawn in spectrum windows and in tables etc.',
                'Options to specify how large the peak crosses/symbols are drawn in spectrum windows',
                'Enables the user to specify special symbols to annotate peaks according to quality']
    options=['Annotation Style','Draw Size','Merit Symbols']
    tabbedFrame = TabbedFrame(guiParent, options=options,
                              grid=(row,0), tipTexts=tipTexts)
    
    frame1, frame2, frame3 = tabbedFrame.frames
    frame1.expandGrid(7,1)
    frame2.expandGrid(2,0)
    frame3.expandGrid(3,1)
    
    self.drawFrame = frame2
    
 
    # Annotations

    tipText = 'Whether to show unassigned spin system numbers like "{123}" in peak annotations; does not affect the display of full residue assignments'
    self.spinSystAnnoSelect = CheckButton(frame1, callback=self.setSpinSystAnno, grid=(0,0),
                                          selected=analysisProject.doSpinSystemAnnotations,
                                          tipText=tipText)
    spinSystAnnoLabel = Label(frame1, text='Spin System Info', grid=(0,1))

    tipText = 'Whether to show assigned atom names in peak annotations rather than resonance numbers, e.g. "21LeuH,N" or "21Leu[55],[56]"'
    self.resonanceAnnoSelect = CheckButton(frame1, callback=self.setResonanceAnno, grid=(1,0),
                                           selected=analysisProject.doAssignmentAnnotations,
                                           tipText=tipText)
    resonanceAnnoLabel = Label(frame1, text='Atom Assignment', grid=(1,1))

    tipText = 'Whether to show the molecular chain code letter in peak annotations, e.g. "21LeuH,N" or "A21LeuH,N"'
    self.chainAnnoSelect = CheckButton(frame1, callback=self.setChainAnno, grid=(2,0),
                                       selected=analysisProject.doChainAnnotations,
                                       tipText=tipText)
    chainAnnoLabel = Label(frame1, text='Chain Assignment', grid=(2,1))

    tipText = 'Whether to show the molecular system code name in peak annotations, e.g. "MS1:21LeuH,N" or "21LeuH,N"'
    self.molSysAnnoSelect = CheckButton(frame1, callback=self.setMolSysAnno, grid=(3,0),
                                       selected=analysisProject.doMolSysAnnotations,
                                       tipText=tipText)
    molSysAnnoLabel = Label(frame1, text='Molecular System', grid=(3,1))

    tipText = 'Whether to show a symbol indicating the quality of a peak; these symbols must be set in the "Merit Symbols" tab'
    self.meritAnnoSelect = CheckButton(frame1, callback=self.setMeritAnno, grid=(4,0),
                                       selected=analysisProject.doMeritAnnotations,
                                       tipText=tipText)
    meritAnnoLabel = Label(frame1, text='Merit Symbol', grid=(4,1))

    tipText = 'Whether to show the details text for a peak in spectrum windows (set in the peak tables)'
    self.detailAnnoSelect = CheckButton(frame1, callback=self.setDetailAnno, grid=(5,0),
                                        selected=analysisProject.doDetailAnnotations,
                                        tipText=tipText)
    detailAnnoLabel = Label(frame1, text='Details', grid=(5,1))

    tipText = 'Whether to ignore above settings, and use a short peak annotation using only residue numbers and one-letter codes, e.g. "21L"'
    self.simpleAnnoSelect = CheckButton(frame1, callback=self.setSimpleAnnotations,
                                        grid=(6,0), tipText=tipText,
                                        selected=analysisProject.doMinimalAnnotations)
    simpleAnnoLabel = Label(frame1, text='Minimal Annotations (overriding option)', grid=(6,1))

    tipText = 'Manually cause an update of all peak annotations; may take some time but recommended to see the immediate effect of changes'
    self.updateFullAnnoButton = Button(frame1, text='Update Full Annotations',
                                       command=self.updateAnnotations, grid=(7,0),
                                       gridSpan=(1,2), sticky='nsew', tipText=tipText)

    # Size
    
    tipText = 'The width of the displayed peak cross/symbol in screen pixels; the number of pixels from the centre point'
    self.pixel_entry = LabeledEntry(frame2, label='Number of pixels',
                                    entry=analysisProject.peakPixelSize,
                                    returnCallback=self.updatePixelSize,
                                    entry_width=8, tipText=tipText)

    tipText = 'Specifies which peak height corresponds to the specific per-isotope ppm size; actual peak symbol sizes are scaled, according to relative height (and volume)'
    self.intensity_entry = LabeledEntry(frame2, label='Height scale',
                                        entry=analysisProject.peakIntensityScale,
                                        returnCallback=self.updateIntensityScale,
                                        entry_width=8, tipText=tipText)
                                        
    tipText = 'Specifies which peak volume corresponds to the specific per-isotope ppm size; actual peak symbol sizes are scaled, according to relative volume (and height)'
    self.volume_entry = LabeledEntry(frame2, label='Volume scale',
                                     entry=analysisProject.peakVolumeScale,
                                     returnCallback=self.updateVolumeScale,
                                     entry_width=8, tipText=tipText)

    tipTexts = ['The kind of nuclear isotope, as found on the axis if a spectrum window',
                'The ppm width of the peak cross/symbol along the specific kind of isotope axis']
    headings = ['Isotope', 'Peak Size (ppm)']
    self.peakSizeEntry = FloatEntry(self, returnCallback=self.setPeakSize, width=10)
    editWidgets = [None, self.peakSizeEntry]
    editGetCallbacks = [None, self.getPeakSize]
    editSetCallbacks = [None, self.setPeakSize]
    self.peak_size_table = ScrolledMatrix(frame2, headingList=headings,
                                          initialRows=5, editWidgets=editWidgets,
                                          editGetCallbacks=editGetCallbacks,
                                          editSetCallbacks=editSetCallbacks,
                                          tipTexts=tipTexts)

    self.selected_method = None
    self.method_widgets = {}
    self.method_widgets[peak_draw_methods[0]] = [self.pixel_entry]
    self.method_widgets[peak_draw_methods[1]] = [self.peak_size_table]
    self.method_widgets[peak_draw_methods[2]] = [self.intensity_entry, self.volume_entry, self.peak_size_table]
    self.method_widgets[peak_draw_methods[3]] = [self.peak_size_table]
    self.method_widgets[peak_draw_methods[4]] = []
    self.method_widgets[peak_draw_methods[5]] = []

    selected_index = peak_draw_methods.index(analysisProject.peakDrawMethod)
    tipText = 'Selects which mode to use for determining peak cross/symbol sizes; a fixed pixel/ppm value or scaled by intensity for all peaks or just one peak list'
    self.method_menu = PulldownList(frame2, texts=peak_draw_methods, gridSpan=(1,2), 
                                    grid=(0,0), tipText=tipText,
                                    callback=self.setMethod, index=selected_index)
    
    # Merit symbol
    
    label = Label(frame3, text='Good merit (>0.66)', grid=(0,0))
    tipText = 'Although typically blank, sets a symbol to display for good quality peaks (merit value > 0.66), if the relevant option is selected in the "Annotation Style" tab'
    self.meritGoodEntry = Entry(frame3, text=analysisProject.meritAnnotationGood,
                                width=8, grid=(0,1), sticky='ew', tipText=tipText,
                                returnCallback=self.setMeritAnnotation)

    label = Label(frame3, text='Medium merit', grid=(1,0))
    tipText = 'Sets a symbol to display for mediocre quality peaks (merit value between 0.34 & 0.66), if the relevant option is selected in the "Annotation Style" tab'
    self.meritUglyEntry = Entry(frame3, text=analysisProject.meritAnnotationMediocre,
                                width=8, grid=(1,1), sticky='ew', tipText=tipText,
                                returnCallback=self.setMeritAnnotation)

    label = Label(frame3, text='Poor merit (<0.34)', grid=(2,0))
    tipText = 'Sets a symbol to display for poor quality peaks (merit value < 0.34), if the relevant option is selected in the "Annotation Style" tab'
    self.meritBadEntry  = Entry(frame3, text=analysisProject.meritAnnotationBad,
                                width=8, grid=(2,1), sticky='ew', tipText=tipText,
                                returnCallback=self.setMeritAnnotation)

    tipText = 'Commit the changes to the merit symbols; this will not automatically update the display - see the "Annotation Style" tab'
    setButton = Button(frame3, text='Set symbols', command=self.setMeritAnnotation,
                       grid=(3,0), gridSpan=(1,2), tipText=tipText, sticky='nsew')

    buttons = UtilityButtonList(tabbedFrame.sideFrame, grid=(0,0), 
                                helpUrl=self.help_url, sticky='e')

    self.updatePeakSizeTable()
    self.updateMethod()

    for func in ('__init__', 'delete', ''):                                           
      self.registerNotify(self.updatePeakSizeTable, 'ccpnmr.Analysis.AxisType', func)

    self.registerNotify(self.changedDoSpinSystemAnnotations, 'ccpnmr.Analysis.AnalysisProject', 'setDoSpinSystemAnnotations')
    self.registerNotify(self.changedDoAssignmentAnnotations, 'ccpnmr.Analysis.AnalysisProject', 'setDoAssignmentAnnotations')
    self.registerNotify(self.changedDoChainAnnotations, 'ccpnmr.Analysis.AnalysisProject', 'setDoChainAnnotations')
    self.registerNotify(self.changedPeakDrawMethod, 'ccpnmr.Analysis.AnalysisProject', 'setPeakDrawMethod')
    self.registerNotify(self.changedPeakPixelSize, 'ccpnmr.Analysis.AnalysisProject', 'setPeakPixelSize')

  def open(self):
  
    self.updatePeakSizeTable()
    self.updateMethod()
    BasePopup.open(self)
    

  def close(self):
  
    self.setMeritAnnotation()
    BasePopup.close(self)

  def destroy(self):
 
    for func in ('__init__', 'delete', ''):                                           
      self.unregisterNotify(self.updatePeakSizeTable, 'ccpnmr.Analysis.AxisType', func)

    self.unregisterNotify(self.changedDoSpinSystemAnnotations, 'ccpnmr.Analysis.AnalysisProject', 'setDoSpinSystemAnnotations')
    self.unregisterNotify(self.changedDoAssignmentAnnotations, 'ccpnmr.Analysis.AnalysisProject', 'setDoAssignmentAnnotations')
    self.unregisterNotify(self.changedDoChainAnnotations, 'ccpnmr.Analysis.AnalysisProject', 'setDoChainAnnotations')
    self.unregisterNotify(self.changedPeakDrawMethod, 'ccpnmr.Analysis.AnalysisProject', 'setPeakDrawMethod')
    self.unregisterNotify(self.changedPeakPixelSize, 'ccpnmr.Analysis.AnalysisProject', 'setPeakPixelSize')

    BasePopup.destroy(self)

  def setMeritAnnotation(self, *opt):
    
    analysisProject = self.analysisProject
    analysisProject.meritAnnotationGood     = self.meritGoodEntry.get() or None
    analysisProject.meritAnnotationMediocre = self.meritUglyEntry.get() or None
    analysisProject.meritAnnotationBad      = self.meritBadEntry.get() or None
   
  def updateAnnotations(self):
  
    self.setMeritAnnotation()
    refreshPeakAnnotations(self.project)
  
  def setSimpleAnnotations(self, trueOrFalse):
  
    self.analysisProject.doMinimalAnnotations = trueOrFalse

  def setSpinSystAnno(self, trueOrFalse):
  
    self.analysisProject.doSpinSystemAnnotations = trueOrFalse

  def setResonanceAnno(self, trueOrFalse):
  
    self.analysisProject.doAssignmentAnnotations = trueOrFalse
    
  def setMolSysAnno(self, trueOrFalse):
  
    self.analysisProject.doMolSysAnnotations = trueOrFalse
    
  def setChainAnno(self, trueOrFalse):
  
    self.analysisProject.doChainAnnotations = trueOrFalse
    
  def setMeritAnno(self, trueOrFalse):

    self.analysisProject.doMeritAnnotations = trueOrFalse

  def setDetailAnno(self, trueOrFalse):   

    self.analysisProject.doDetailAnnotations = trueOrFalse
    
  def changedDoSpinSystemAnnotations(self, analysisProject):

    self.spinSystAnnoSelect.set(analysisProject.doSpinSystemAnnotations)

  def changedDoAssignmentAnnotations(self, analysisProject):

    self.resonanceAnnoSelect.set(analysisProject.doAssignmentAnnotations)

  def changedDoChainAnnotations(self, analysisProject):

    self.chainAnnoSelect.set(analysisProject.doChainAnnotations)

  def changedPeakDrawMethod(self, analysisProject):
 
    self.updateMethod()

  def changedPeakPixelSize(self, analysisProject):
 
    self.pixel_entry.setEntry(analysisProject.peakPixelSize)

  def setMethod(self, text):
  
    self.analysisProject.peakDrawMethod = text

  def updateMethod(self):

    selected = self.analysisProject.peakDrawMethod
    if selected == self.selected_method:
      return

    self.unsetMethod()

    widgets = self.method_widgets[selected]

    row = 1
    for widget in widgets:
      if isinstance(widget, ScrolledMatrix):
        widget.grid(row=row, column=0, sticky='nsew')
        self.drawFrame.grid_rowconfigure(row, weight=1)
      else:
        widget.grid(row=row, column=0, sticky='new')
        self.drawFrame.grid_rowconfigure(row, weight=0)
  
      row += 1

    self.selected_method = selected

  def unsetMethod(self):

    selected = self.selected_method
    if selected is None:
      return

    widgets = self.method_widgets[selected]

    for widget in widgets:
      widget.grid_forget()

    self.selected_method = None

  def updatePeakSizeTable(self, *extra):
 
    axisTypes = self.parent.getAxisTypes()
 
    textMatrix = []
    n = 0
    for axisType in axisTypes:
      if (len(axisType.isotopeCodes) == 1):
        n = n + 1
        text = []
        text.append(', '.join(axisType.isotopeCodes))
        text.append(axisType.peakSize)
        textMatrix.append(text)
 
    self.peak_size_table.update(objectList=axisTypes, textMatrix=textMatrix)

  def getPeakSize(self, axisType):
 
    self.peakSizeEntry.set(axisType.peakSize)
 
  def setPeakSize(self, *extra):
 
    axisType = self.getAxisType()
    try:
      axisType.peakSize = self.peakSizeEntry.get()
    except Implementation.ApiError, e:
      showError('Setting peak size', e.error_msg, parent=self)

  def getAxisType(self):

    return self.peak_size_table.currentObject

  def updatePixelSize(self, *extra):

    try:
      size = int(self.pixel_entry.getEntry())
      self.analysisProject.peakPixelSize = size
    except:
      pass

  def updateIntensityScale(self, *extra):

    try:
      scale = float(self.intensity_entry.getEntry())
      self.analysisProject.peakIntensityScale = scale
    except:
      pass

  def updateVolumeScale(self, *extra):

    try:
      scale = float(self.volume_entry.getEntry())
      self.analysisProject.peakVolumeScale = scale
    except:
      pass

