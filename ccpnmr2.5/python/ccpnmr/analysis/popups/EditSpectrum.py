
"""
======================COPYRIGHT/LICENSE START==========================

EditSpectrum.py: Part of the CcpNmr Analysis program

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
import os
import memops.gui.Color as Color

from memops.general import Implementation

from memops.universal.Util import formatFloat, isMacOS

from memops.gui.Button import Button
from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.Frame import Frame
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FloatEntry import FloatEntry
from memops.gui.FontList import FontList
from memops.gui.IntEntry import IntEntry
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.LabelDivider import LabelDivider
from memops.gui.MessageReporter import showYesNo, showError, showWarning, showOkCancel
from memops.gui.MultiWidget import MultiWidget
from memops.gui.PulldownList import PulldownList
from memops.gui.RadioButtons import RadioButtons
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccp.general.Io import getDataSourceFileName, setDataSourceFileName

from ccp.gui.DataLocationFrame import DataLocationFrame

from ccpnmr.analysis.core.ExperimentBasic import getNoiseEstimate, isDataBigEndian, isShapeSpectrum
from ccpnmr.analysis.core.ExperimentBasic import getPrimaryDataDimRef, setAllIsotopeCodes
from ccpnmr.analysis.core.MoleculeBasic import STANDARD_ISOTOPES
from ccpnmr.analysis.core.Reference import getReference, shiftDataDimRef
from ccpnmr.analysis.core.Util import shift_key_modifier, ctrl_key_modifier
from ccpnmr.analysis.core.Util import spectrum_shortcuts, getHueSortedColorSchemes
from ccpnmr.analysis.popups.BasePopup import BasePopup

# TBD 
# spec colour scheme backward compat re-point
#
# Test
#


from ccp.api.nmr import Nmr

from ccp.general.Constants import chemShiftRefRatios

standardHeadingList = ('Dim', 'Isotope', 'Spectrometer\nfrequency (MHz)',
                       'Spectral\nwidth (ppm)',  'Spectral\nwidth (Hz)', 'Reference\nppm', 'Reference\npoint',
                       'Orig. number\nof points', 'Point\noffset',
                       'Minimum aliased\nfrequency (ppm)',
                       'Maximum aliased\nfrequency (ppm)')
standardTipTexts = ('Spectrum dimension or sub dimension',
  'Isotope being recorded',
  'Carrier frequency for the particular isotope, in MHz',
  'Actual spectral width, so possibly after region truncated, in ppm',
  'Actual spectral width, so possibly after region truncated, in Hz',
  'Ppm value at chosen reference point',
  'Point value where dimension has chosen reference ppm', 
  'Original number of points, so after Fourier transform and before any truncation of the region',
  'Offset of region relative to original region after Fourier transform, so 0 if the region has not been truncated',
  'Minimum frequency, in ppm, of where any signal is expected, if not set then uses minimum frequency of fundamental region',
  'Maximum frequency, in ppm, of where any signal is expected, if not set then uses maximum frequency of fundamental region',
)

reducedDimHeadingList = standardHeadingList[:1] + ('Sub Dim\nName',) + standardHeadingList[1:2] + ('Scaling\nFactors',) + standardHeadingList[2:]

reducedTipTexts = standardTipTexts[:1] + ('Name of sub-dimension',) + standardTipTexts[1:2] + ('Scaling factors',) + standardTipTexts[2:]

number_bytes_entries = range(1, 9)
number_type_entries = ('int', 'float')

shortcut_entries = [str(None)] + spectrum_shortcuts
if isMacOS():
  shortcut_entries = shortcut_entries + ['Alt-%s' % x for x in spectrum_shortcuts]
  shortcut_entries = shortcut_entries + ['Shift-%s' % x for x in spectrum_shortcuts]
else:
  shortcut_entries = shortcut_entries + ['Ctrl-%s' % x for x in spectrum_shortcuts]

class EditSpectrumPopup(BasePopup):

  """
  **Curate Spectrum Parameters**
  
  This popup window allows the user to view and edit information that is
  associated with spectra. In this regard a spectrum is a record of the NMR data
  that was *obtained*. Parameters that relate to what was *done experimentally*
  to get the data are specified elsewhere; experiment level information is set
  via the main Experiments_  popup window. Many spectrum level parameters in
  Analysis relate to how the spectrum data files are interpreted and displayed
  graphically, i.e. as contours. Data that derive from the analysis of spectra,
  like peak lists and assignments, are managed elsewhere in specialised
  displays.

  The layout of this popup window involves several tab panels that sub-divide
  the spectrum settings into related groups. The functionality within each
  tabbed sub-division is largely independent, although a spectrum entity
  selected in any of the tables or pulldown menus is subject to the [Delete
  Spectrum] function.

  **Main "Spectra" tab**
  
  This table presents an overview of all of the spectra that have been loaded
  into the CCPN project, what their dimensions are and how they relate to
  experiments and peak lists. Some of the parameters may be edited by
  double-clicking on a cell in the table, thus for example the user may change
  the name and set the active peak list (the one that is used for peak peaking
  and selection). In the main however, this table usually serves simply to
  illustrate what is present, which experiment each spectrum belongs to and
  what its NMR dimensions are.

  **Display Options**
  
  This table is used to control how spectra are visualised within the spectrum
  window displays of Analysis. Specifically, the user can: change the colours
  of the displayed contours and 1D slice traces; get access to the contour
  level settings, via the `Spectrum Contour Levels`_ tool; access the `Spectrum
  Contour Files`_ setup; set keyboard toggles; and administer other graphical
  elements like bounding box lines, peak fonts and peak pointer lines. It
  should be noted that the colours for peak text and peak symbols are set for
  individual peak lists, so that different peak lists within the same spectrum
  may be distinguished. Accordingly these parameters are set via the main 
  `Peak Lists`_ option.

  The "Rank" value of a spectrum dictates the order in which spectra are drawn
  on screen and presented in the toggle options at the top of spectrum windows.
  The rank may be set on an individual basis or following the order of spectra
  in the table using [Set Ranks From Order]; this is especially useful if the
  table is first sorted in a particular way (by clicking on a table heading).

  **Referencing**
  
  This table is used to control how the points in a spectrum data file are
  interpreted in terms of NMR frequency, on each of their axes. Usually, this
  means to know how positions in the spectrum data, stored as an array of
  intensity data points, relate to ppm units and Hz frequency values of the
  spectrometer. Such referencing occurs by specifying how the width of a
  spectrum dimension in data points relates to its width in frequency units,
  defining the scaling factor, and how the two scales are anchored together; by
  giving a reference data point and its corresponding a frequency/ppm value.

  The spectrum referencing values are normally extracted from the spectrum
  header or parameter file when the spectrum was first loaded. Often these
  values do not need to be altered. However, the user may still adjust all of
  the parameters after load time, with the sole exception of the dimensions
  isotope type. Changes to spectrum referencing may be made after loading to
  correct for mistakes in the initial data interpretation and to tweak spectra
  so that they give a better overlay of signals. Note that the table columns
  'Spectrometer Frequency', 'Minimum Aliased Frequency', 'Maximum aliased
  frequency', and 'Sub Dim Name' (normally hidden) apply to the experiment, and
  will are shared between all spectra in the same experiment.

  Peaks that are picked within spectra have positions defined in terms of
  spectrum data points, i.e. normally indicating an extremum in intensity. The
  referencing values of a spectrum's dimensions are used to calculate what the
  position of a peak is in ppm or Hz units. Accordingly, if spectrum
  referencing values are changed then peak ppm locations will naturally move.
  If peaks should not move with the spectrum data when referencing is changed,
  to instead stay at the same ppm values, then the user can select "ppm" in the
  "reference changes keep constant peak:" options.

  Some less common experiments have more than one type of signal (shift) on 
  each axis, and need a reference line to describe each experiment type. One
  obvious case is projection experiments, where signals in one dimension are
  linear combinations of frequencies for different nuclei. Another would be the
  15N/13C HSQC-NOESY, where 15N and 13C are measured in parallel on one axis.
  The special options for this situation are turned on by selecting  'Use
  reduced dimensionality options'. The individual reference lines are called
  'Sub-dimensions'; they are shown explicitly when there is more than one. They
  can be added, removed, and edited individually. The Sub-dimensions set for the
  first Spectrum serve to (re)define those for the Experiment; once there is
  more than one Spectrum in an Experiment it is not possible to add extra
  Sub-dimensions. Projection spectra can only be read correctly from some
  formats (Bruker and USF3) and even these cases may require editing. The
  special 'Sub Dim Name' column gives a name that is used by some external
  programs (e.g. Prodecomp) to identify a given signal type across different
  experiments. The special 'Scaling Factors' column gives the scaling
  coefficient that applies to a particular signal in this spectrum. For instance
  scaling factors  CO:1; CA:0; CB:-2 would denote a signal of the form wN-2wCB.
  There can be more than one scaling factor where different types of peaks are
  possible, e.g. wN OR wC. The system has only been tested for 2D projection
  spectra.


  The minimum and maximum aliased frequency values for a spectrum dimension, as
  displayed in the last columns of the table, may be set to cut off the
  spectrum contour display at a particular value or extend the contour display.
  Contour displays that are extended beyond the normal sweep width (ppm range)
  of the spectrum do so by tiling; a copy of the contours (not mirror image) is
  made and placed next to the normal, fundamental region. This mechanism is
  used to cope with aliased resonances that appear (as 'ghost' peak images)
  within the main spectrum width in places that are not their true ppm values;
  the resonance is really outside the sweep width of the spectrum. Making tiled
  copies of the contours allows the spectrum to extend and cover the real ppm
  values of aliased peaks. Aliased peaks may be moved to their true ppm
  locations or picked at that location in the first instance by using the
  appropriate tile region.

  **Tolerances**
  
  This section is used to control the weightings and tolerances that apply to
  individual spectrum dimensions during the resonance assignment process. There
  are two kinds of value that can be set. One is the "Assignment Tolerance",
  which states how close (typically in ppm) chemical shift values of resonances
  match to peak peak locations; this is an upper limit such that resonances
  with shifts outside the tolerance are not considered as assignment
  possibilities. The second thing is the "Shift Weighting" that specifies how
  much influence resonance assignments have on the averaged chemical shift
  values. Both of these parameters affect individual spectrum dimensions
  because each axis has is its own resonance line shape and data resolution.

  The normal chemical shift calculation within Analysis calculates values as
  the weighted mean over all of the peak positions to which a resonance is
  assigned. If a dimension has a shift weighting of 2.0 but others have 1.0
  then any resonance assignment in the former will have twice as much influence
  on the shift mean, as if it were assigned to two peak dimensions rather than
  one. By default all weights are set to 1.0, but the user can increase this if
  a spectrum dimension is particularly precise (narrow peaks) or decrease it if
  peaks are broad or if resonance assignments are highly ambiguous; for
  speculative NOE assignments the user usually doesn't want the shift values to
  change according to those assignments. 

  The assignment tolerance and shift weighting values are set for the
  dimensions of a selected spectrum in the upper table. The lower "Spectrum
  Overview" table shows these values across all spectra in the project, so that
  the user can see the overall situation and range of values. Although
  tolerances are not set in the lower table on an individual basis, the
  assignment tolerances and shift weights can be spread from one spectrum (the
  last selected) to any others selected in the table (using left-click + <Crtl>
  or <Shift), where the isotopes of dimensions match.

  **File Details**
  
  This tab is used to indicate the location of a spectrum's binary data file 
  and how the binary data within is read. In essence the binary values of the
  file need to be interpreted as spectrum intensity values, but in order to do
  this the system has to know how the numbers are stored. This involves knowing
  how numbers are represented per se; how many bytes a number uses, whether
  numbers are integers or floating point and what orientation the bits are in
  (big or little endian). Also, most spectrum data formats involve binary
  blocks; chunks of data organised into a grid across all the dimensions. The
  benefit of using blocks comes when looking at large data sets, e.g. using 3D
  and higher dimensionality data, because it means that when only looking at a
  sub-region of a spectrum, like a single plane of a 3D, only part of the
  spectrum data needs to be read; most of the spectrum data that isn't
  displayed does not need to be loaded. Naturally, the system needs to know
  what the sizes of the data block along each dimension are. These are
  displayed in the lower table and may be adjusted if there are any mistakes,
  which may have been inherited from the spectrum processing.

  The "File name" field shows the present location of the binary data. If this
  file is missing the spectrum will still be represented in Analysis, but the
  operations which rely of having intensity data, like picking peaks, cannot be
  performed. It is notable that the binary data file for a spectrum can be
  changed at any time, either by selecting a different file name (typing
  directly or using [Choose File]) or by directly replacing the binary data
  inside the file. If the binary data or any of the parameters that influence
  its interpretation has changed, then the user can click on [Refresh Spectrum]
  to see the effect if the changes.

  **Data Locations**
  
  This last table is used to show the locations of all of the spectrum data
  files used by the whole project. These locations need to be changed if the
  spectrum data is moved in disk, otherwise the system will not know where to
  look. Each data location is broken down into two components; an "absolute"
  directory path and a "relative" path, which contains the file name. Both of
  these components may be changed.

  The absolute part of the location is from the root of the file system (the
  all-encompassing directory) and is often shared by several spectrum files; it
  is typical to store spectra in one repository. For convenience, changing the
  absolute part of the path for one spectrum will usually change the value for
  all the other spectra that use the same directory. Also, if multiple files
  are moved to the same directory then [Propagate Absolute Path] can be used to
  state that several files are in the same location.

  The relative part of the spectrum data path is individual to the spectrum and
  contains the name of the data file and any extra directory location, relative
  to the absolute part of the path. Which parts of the spectrum file path is
  considered to be absolute or relative may be changed using the "Shift
  Directory..." buttons. For example an absolute : relative split of
  "/data/spectra/" : "protein/123.spec" can be changed to
  "/data/spectra/protein" : "123.spec", so that several spectra can use the
  absolute part.

  **Caveats & Tips**
   
  The spectrum chemical shift weighting used for synthetic peak lists is set to
  0.0001 of the the regular weighting for the spectrum dimension set via this
  popup.

  The fonts used for the printing of spectra (e.g. via PostScript of PDF) may
  be set separately to the main spectrum font specification, which is used
  within Analysis spectrum displays. 

  .. _Experiments: EditExperimentPopup.html
  .. _`Spectrum Contour Levels`: EditContourLevelsPopup.html
  .. _`Spectrum Contour Files`: EditContourFilesPopup.html
  .. _`Peak Lists`: EditPeakListsPopup.html
  """

  def __init__(self, parent, spectrum=None, modal=False, useReducedDim=False, 
               *args, **kw):
    
    self.spectrum = spectrum
    self.analysisSpec = spectrum
    self.tolSpectrum = spectrum
    self.refSpectrum = spectrum
    self.fileSpectrum = spectrum
    self.waiting  = False
    self.waitingF = False
    self.waitingA = False
    self.waitingT = False
    self.waitingR = False
    self.analysisDataDim = None
    self.fileDim = None
    self.isModal = modal
    self.dataDimRef = None
    self.useReducedDim = useReducedDim
   
    if modal:
      self.cancelled = False
      title = 'Verify Spectra'
    else:
      title = 'Experiment : Spectra'
    
    BasePopup.__init__(self, parent=parent, modal=modal, title=title, **kw)

  def cancel(self):

    self.cancelled = True
    self.close()
  
  def close(self):
    
    if not (hasattr(self, 'cancelled') and self.cancelled):
      if self.isModal and len(self.spectrum.experiment.dataSources) == 1:
        # Isotopes are editable - make sure they are set:
        spectrum = self.refSpectrum
        if spectrum:
          setAllIsotopeCodes(spectrum.experiment)
    #
    BasePopup.close(self)
   
  def body(self, guiFrame):
    
    self.geometry('780x500')
  
    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    if self.isModal:
      options = ['Verify Referencing','Verify File Details']
      tipTexts = ['Information for interpreting data points in terms of PPM & Hz',
                  'Location of spectrum data files and how data is stored in blocks & points',]
      self.tabbedFrame = TabbedFrame(guiFrame, options=options, tipTexts=tipTexts)
      frameC, frameE = self.tabbedFrame.frames
     
    else:
      options = ['Spectra','Display Options','Referencing',
                 'Tolerances','File Details', 'Data Locations']
      tipTexts = ['Basic information about spectra, their experiments, peak lists & dimensions',
                  'How spectra are displayed on screen in terms of colours & contours etc.',
                  'Information for interpreting data points in terms of PPM & Hz',
                  'Per-dimension assignment tolerances & chemical shift calculation contributions',
                  'Location of spectrum data files and how data is stored in blocks & points',
                  'A table of all spectrum data locations for en masse manipulation']
      self.tabbedFrame = TabbedFrame(guiFrame, options=options, tipTexts=tipTexts,
                                     callback=self.updateSpectrumButtons)
      frameA, frameB, frameC, frameD, frameE, frameF = self.tabbedFrame.frames
    
    self.tabbedFrame.grid(row=0, column=0, sticky='nsew')
    # Main spectrum table
    
    if not self.isModal:
    
      frameA.grid_rowconfigure(0, weight=1)
      frameA.grid_columnconfigure(0, weight=1)
 
      self.expNameEntry  = Entry(self, width=20,
                                 returnCallback=self.setExperimentName)
      self.specNameEntry = Entry(self, width=20,
                                 returnCallback=self.setSpectrumName)
      self.peakListPulldown = PulldownList(self, callback=self.setActivePeakList)
      self.scaleEntry = FloatEntry(self, returnCallback=self.setScale, width=10)
      self.noiseEntry = FloatEntry(self, returnCallback=self.setNoise, width=15)
 
      headingList = ['#', 'Experiment', 'Spectrum', 'Dimensions', 'Active\nPeak List',
                     'Num Peak\nLists', 'Scale', 'Noise\nLevel', 'Data Type',
                     'Simulated?','Data\nFile']
 
      tipTexts = [ 'Row number',
                   'The name of the experiment to which the spectrum belongs',
                   'A short name for the spectrum, often used with experiment name for unique identification "experiment:spectrum". e.g. "HSQC:115"',
                   'Number of dimensions/axes in the spectrum (can be fewer than in the experiment) and their isotope types',
                   'The active peak list, for example the peak list into which picked peaks are put',
                   'Number of separate peak lists within the spectrum',
                   'Spectrum scale relative to stored data points, e.g. for drawing',
                   'Estimated intensity level of spectrum noise',
                   'Whether fully processed, FID or part-processed',
                   'Whether the spectrum is simulated',
                   'Location of the data file relative to the data directory']

      editWidgets = [ None,
                      self.expNameEntry, self.specNameEntry,
                      None, self.peakListPulldown,
                      None, self.scaleEntry,
                      self.noiseEntry, None, None, None ]
 
      editGetCallbacks = [ None,
                           self.getExperimentName, self.getSpectrumName,
                           None, self.getActivePeakList,
                           None, self.getScale,
                           self.getNoise, None, None, None]
 
      editSetCallbacks = [ None,
                           self.setExperimentName, self.setSpectrumName,
                           None, self.setActivePeakList,
                           None, self.setScale,
                           self.setNoise, None, None, None]
 
      justifyList = ['center'] * len(headingList)
      justifyList[3] = 'left'
 
      self.spectrumTable = ScrolledMatrix(frameA,headingList=headingList,
                                          callback=self.selectSpectrum,
                                          editWidgets=editWidgets, multiSelect=True,
                                          editGetCallbacks=editGetCallbacks,
                                          editSetCallbacks=editSetCallbacks,
                                          justifyList=justifyList,
                                          tipTexts=tipTexts,
                                          deleteFunc=self.deleteSpectra)
      self.spectrumTable.grid(row=0, column=0, sticky='nsew')

      texts = [ 'Delete', 'File\nDetails', 'Referencing', 'Tolerances']
      commands = [ self.deleteSpectra, self.editFileDetails,
                   self.editReferencing, self.editTolerances]
 
      # Display options
 
      frameB.grid_rowconfigure(0, weight=1)
      frameB.grid_columnconfigure(0, weight=1)
 
      headingList = ['#','Spectrum', 'Positive\nColours', 'Negative\nColours',
                     'Slice\nColour', 'Box\nVisible?', 'Shortcut', 'Rank',
                     'Positive\nContours','Negative\nContours','Use Contour\nFiles',
                     'Peak\nFont', 'Peak\nPointer?']

      tipTexts = [ 'Row number',
                   'Identifying name Experiment:Spectrum',
                   'Colour scheme for positive contours (can be a single colour)',
                   'Colour scheme for negative contours (can be a single colour)',
                   'Colour of 1D slice and trace lines',
                   'Whether fundamental region (sweep width bounding) box should be drawn',
                   'Shortcut key for toggling spectrum on/off',
                   'Determines order spectra drawn in (lowest appears on top)',
                   'Positive contour levels, relative to global scale',
                   'Negative contour levels, relative to global scale',
                   'Whether pre-calculated contours should be used',
                   'The typeface used for peak annotation in spectrum window displays',
                   'Whether the peak pointer (a triangle/line from annotation to centre), if it exists, should be drawn' ]
 
      mode = self.analysisProfile.graphicsHandler

      self.rankEntry = IntEntry(self, returnCallback=self.setRank, width=10)
      self.shortcutPulldown = PulldownList(self, callback=self.setShortcut, texts=shortcut_entries)
      self.fontMenu  = FontList(self, callback=self.setSpectrumFont, mode=mode)
      self.posColorPulldown = PulldownList(self, callback=self.setPosColor)
      self.negColorPulldown = PulldownList(self, callback=self.setNegColor)
      self.sliceColorPulldown = PulldownList(self, callback=self.setSliceColor)
 
      editWidgets = [None, None, self.posColorPulldown, self.negColorPulldown,
                     self.sliceColorPulldown, None, self.shortcutPulldown, self.rankEntry,
                     None, None, None, self.fontMenu, None ]
 
      editGetCallbacks = [None, None, self.getPosColor, self.getNegColor,
                          self.getSliceColor, self.toggleBoxVisible,
                          self.getShortcut, self.getRank,
                          self.getContourLevels, self.getContourLevels,
                          self.toggleTryPrecalcContours, self.getSpectrumFont,
                          self.togglePointerVisible]
 
      editSetCallbacks = [None, None, self.setPosColor, self.setNegColor,
                          self.setSliceColor, None,
                          self.setShortcut, self.setRank,
                          None, None,
                          None, self.setSpectrumFont,
                          None]
 
      self.displayTable = ScrolledMatrix(frameB,headingList=headingList,
                                         callback=self.selectAnalysisSpec,
                                         editWidgets=editWidgets, multiSelect=True,
                                         editGetCallbacks=editGetCallbacks,
                                         editSetCallbacks=editSetCallbacks,
                                         tipTexts=tipTexts)
      self.displayTable.grid(row=0, column=0, sticky='nsew')
      self.displayTable.doEditMarkExtraRules = self.doEditMarkExtraRules
 
      texts = ['Contour\nLevels', 'Contour\nFiles',
               'Propagate\nContours', 'Set Ranks\nFrom Order']
      tipTexts = ['Open popup for detailed setting of contour levels',
                  'Open popup for administering pre-calculated contour files',
                  'Propagate contour levels from row selected last to others',
                  'Set spectrum rank based on current order in table (e.g. after sorting on a column)']
      commands = [self.editLevels, self.editContourFiles,
                  self.propagateContours, self.orderRanks]
      self.displayButtons = ButtonList(frameB, texts=texts, tipTexts=tipTexts, expands=True, commands=commands)
      self.displayButtons.grid(row=1, column=0, sticky='ew')
    
    # Referencing
    
    frameC.grid_columnconfigure(4, weight=1)

    row = 0
    tipText = 'Select the experiment:spectrum to adjust'
    label = Label(frameC, text='Spectrum: ')
    label.grid(row=row, column=0, sticky='w')
    self.refSpecPulldown = PulldownList(frameC, callback=self.selectRefSpectrum,
                                        tipText=tipText,
                                        grid=(row,1))
                                        
    label = Label(frameC, text=' Use reduced\ndimensionality options:')
    label.grid(row=row, column=2, sticky='w')
    self.reducedDimSelect = CheckButton(frameC, callback=self.toggleReducedDim,
                                        tipText='Used for reduced dimensionality or ' + \
                                        'projected experiments when you have sub-dimensions')
    self.reducedDimSelect.grid(row=row, column=3, sticky='w')
    self.reducedDimSelect.set(self.useReducedDim)
    
    if not self.isModal:
      frame = Frame(frameC)
      frame.grid(row=row, column=4, sticky='nsew')
      label = Label(frame, text='Reference changes\nkeep constant peak: ')
      label.grid(row=0, column=0, sticky='ew')
      entries = [ 'point', 'ppm' ]
      tipTexts = [ 'If you change the referencing then peak position in points stays the same',
                   'If you change the referencing then peak position in ppm stays the same' ]
      self.peakRefSelect = RadioButtons(frame, entries=entries, direction='horizontal', tipTexts=tipTexts)
      self.peakRefSelect.grid(row=0, column=1, sticky='ew')

    row = row + 1
    frameC.grid_rowconfigure(row, weight=1)

    # only allow editing of isotope if no existing spectrum for that experiment
    if self.isModal and len(self.spectrum.experiment.dataSources) == 1:
      self.getIsotope = self.getIsotopeFunc
      self.setIsotope = self.setIsotopeFunc
      self.isotopeMenu = PulldownList(self, texts=STANDARD_ISOTOPES, callback=self.setIsotope)
  
    else:
      self.isotopeMenu = None
      self.getIsotope = None
      self.setIsotope = None
      
    self.sfEntry          = FloatEntry(self, returnCallback=self.setSf, width=10, formatPlaces=8)
    self.swEntry          = FloatEntry(self, returnCallback=self.setSw, width=12, formatPlaces=8)
    self.refppmEntry      = FloatEntry(self, returnCallback=self.setRefppm, width=10, formatPlaces=8)
    self.refptEntry       = FloatEntry(self, returnCallback=self.setRefpt, width=12, formatPlaces=8)
    self.origNptsEntry    = IntEntry(self, returnCallback=self.setOrigNpts, width=6)
    self.pointOffsetEntry = IntEntry(self, returnCallback=self.setPointOffset, width=5)
    self.minFreqEntry     = FloatEntry(self, returnCallback=self.setMinFreq, width=10, formatPlaces=6)
    self.maxFreqEntry     = FloatEntry(self, returnCallback=self.setMaxFreq, width=10, formatPlaces=6)
    self.scaleFactorsWidget = MultiWidget(self, FloatEntry, callback=self.setDimScaling,
                                          minRows=0, useImages=False) 
    self.subDimNameEntry = Entry(self, returnCallback=self.setSubDimName, width=10)

    headingList = standardHeadingList
    tipTexts = standardTipTexts
    self.dataDimRefTable = ScrolledMatrix(frameC, headingList=headingList,
                                          callback=self.selectDataDimRef,
                                          tipTexts=tipTexts)
    self.dataDimRefTable.grid(row=row, column=0, columnspan=5, sticky='nsew')

    row += 1
    
    texts = ['Add Sub-dimension Copy','Remove Sub-dimension']
    tipTexts = ['Add a new sub dimension (e.g. reduced dimensionality component) by copying the selected sub-dimension',
                'Remove selected sub-dimension']
    commands = [self.addSubDim, self.removeSubDim]
      
    self.refButtons = ButtonList(frameC, texts=texts, tipTexts=tipTexts,
                                 commands=commands, expands=True)
    self.refButtons.grid(row=row, column=0, columnspan=5, sticky='ew')
    self.addSubDimButton, self.removeSubDimButton = self.refButtons.buttons[:2]

    # Tolerances
   
    if not self.isModal:
   
      frameD.expandGrid(3,1)
 
      label = Label(frameD, text='Spectrum: ', grid=(0,0))
      self.tolSpecPulldown = PulldownList(frameD, callback=self.selectTolSpectrum, grid=(0,1))

      headingList = ('Dimension', 'Isotope', 'Assignment\nTolerance','Shift\nWeighting')
      tipTexts = ('Spectrum dimension',
                  'Isotope recorded in this dimension',
                  'When assigning, the maximum difference between a peak dimension position and a known chemical shift value',
                  'Contribution, in this dimension, of peak positions to the averaged, weighed chemical shift values')
      self.assignTolEntry = FloatEntry(self, returnCallback=self.setAssignmentTol, width=6)
      self.shiftWeightEntry = FloatEntry(self, returnCallback=self.setShiftWeight, width=6)
      editWidgets = [ None, None, self.assignTolEntry, self.shiftWeightEntry]
      editGetCallbacks = [ None, None, self.getAssignmentTol, self.getShiftWeight ]
      editSetCallbacks = [ None, None, self.setAssignmentTol, self.setShiftWeight ]
 
      self.tolTable = ScrolledMatrix(frameD, headingList=headingList, tipTexts=tipTexts,
                                     callback=self.selectDataDim,
                                     editWidgets=editWidgets,
                                     editGetCallbacks=editGetCallbacks,
                                     editSetCallbacks=editSetCallbacks,
                                     grid=(1,0), gridSpan=(1,2))
    
      div = LabelDivider(frameD, text='Spectrum Overview', grid=(2,0), gridSpan=(1,2))
      
      headingList = ['#','Spectrum','Dim Assign\nTolerances','Dim Shift\nWeights']
      tipTexts = ['Row number',
                  'Experiment name : Spectrum name',
                  'Peak vs shift assignment tolerances in the various dimensions, as set in the Spectrum table',
                  'Chemical shift weighting contribution, in the various dimensions, as set in the Spectrum table']
      self.tolTableAll = ScrolledMatrix(frameD, headingList=headingList, tipTexts=tipTexts,
                                        callback=self.selectTolSpectrum,
                                        multiSelect=True,
                                        grid=(3,0), gridSpan=(1,2))
      
      texts = ['Propagate Assignment Tolerances','Propagate Shift Weights']
      tipTexts = ['Propagate peak vs shift assignment tolerance from last selected spectrum to other selected spectra in Spectrum Overview table',
                  'Propagate chemical shift weighting contribution from last selected spectrum to other selected spectra in Spectrum Overview table']
      commands = [self.propagateAssignTols, self.propagateShiftWeights]
      buttons = ButtonList(frameD, texts=texts, tipTexts=tipTexts, commands=commands, grid=(4,0), gridSpan=(1,2))
    
    
    # File details

    frameE.expandGrid(2,1)
    
    tipText = 'Select the experiment:spectrum to operate on'
    label = Label(frameE, text='Spectrum: ', tipText=tipText)
    label.grid(row=0, column=0, sticky='w')
    self.fileSpecPulldown = PulldownList(frameE, callback=self.selectFileSpectrum, tipText=tipText)
    self.fileSpecPulldown.grid(row=0, column=1, sticky='w')
    
    fileFrame = LabelFrame(frameE, text="Data File Details")
    fileFrame.grid(row=1,column=0, columnspan=2, sticky='ew')
    fileFrame.grid_columnconfigure(1, weight=1)
    
    tipText = 'Full path to the file containing the spectrum data (typically binary)'
    label = Label(fileFrame, text='File name:', tipText=tipText)
    label.grid(row=0, column=0, sticky='w')
    self.specFileEntry = Entry(fileFrame, width=32, tipText=tipText)
    self.specFileEntry.grid(row=0, column=1, sticky='ew')
 
    self.dataOptFrame = Frame(fileFrame)
    self.dataOptFrame.grid(row=1, column=0, columnspan=2, sticky='ew')
    self.dataOptFrame.grid_columnconfigure(7, weight=1)
    
    tipText = 'Whether binary numbers in data are big endian (e.g. Mac PowerPC, SGI)\n' + \
              'or little endian (e.g. Intel), which depends on how and where data was collected and processed'
    label = Label(self.dataOptFrame, text='Is data big endian?:', tipText=tipText)
    label.grid(row=0, column=0, sticky='w')
    self.endianSelect = CheckButton(self.dataOptFrame, tipText=tipText)
    self.endianSelect.grid(row=0, column=1, sticky='w')

    tipText = 'Bytes per data point (normally 4)'
    label = Label(self.dataOptFrame, text='Bytes per word:', tipText=tipText)
    label.grid(row=0, column=2, sticky='w')
    self.numBytesPulldown = PulldownList(self.dataOptFrame, texts=map(str, number_bytes_entries), tipText=tipText)
    self.numBytesPulldown.grid(row=0, column=3, sticky='w')

    tipText = 'Whether data is floating point or integer (normally floating point)'
    label = Label(self.dataOptFrame, text='Data number type:', tipText=tipText)
    label.grid(row=0, column=4, sticky='w')
    self.numTypePulldown = PulldownList(self.dataOptFrame, texts=number_type_entries, tipText=tipText)
    self.numTypePulldown.grid(row=0, column=5, sticky='w')

    tipText = 'The header size in bytes, which depends on the file format'
    label = Label(self.dataOptFrame, text='File header size:', tipText=tipText)
    label.grid(row=0, column=6, sticky='w')
    self.headerSizeEntry = IntEntry(self.dataOptFrame, width=6, tipText=tipText)
    self.headerSizeEntry.grid(row=0, column=7, sticky='w')
 
    texts = [ 'Choose File','Refresh Spectrum' ]
    tipTexts = [ 'Select a new spectrum data file using a dialog, rather than typing it above',
                 'Save the information as chosen and refresh spectrum display' ]
    commands = [ self.browseSpecFile, self.commitFileDetails ]
    self.fileDetailsButtons = ButtonList(fileFrame, texts=texts, expands=True, commands=commands, tipTexts=tipTexts)
    self.fileDetailsButtons.grid(row=2, column=0, columnspan=2, sticky='ew')
    
    headingList = ('Dimension', 'Number\nPoints', 'Block\nSize')
    tipTexts = ('Dimension in data file',
                'Number of data points in this dimension',
                'Size of blocks into which data points are partitioned for efficient reading')
    self.npointsEntry = IntEntry(self, returnCallback=self.setNpoints, width=8)
    self.blockSizeEntry = IntEntry(self, returnCallback=self.setBlockSize, width=5)
    editWidgets = [ None, self.npointsEntry, self.blockSizeEntry ]
    editGetCallbacks = [ None, self.getNpoints, self.getBlockSize ]
    editSetCallbacks = [ None, self.setNpoints, self.setBlockSize ]
    self.fileDimTable = ScrolledMatrix(frameE, headingList=headingList,
                                       tipTexts=tipTexts, initialRows=5,
                                       callback=self.selectFileDim, 
                                       editWidgets=editWidgets,
                                       editGetCallbacks=editGetCallbacks,
                                       editSetCallbacks=editSetCallbacks)
    self.fileDimTable.grid(row=2, column=0, columnspan=2, sticky='nsew')

    # Data Locations
   
    if not self.isModal:

      frameF.expandGrid(0,0)
      self.dataLocationFrame = DataLocationFrame(frameF, self.project, grid=(0,0))
      
    # Main popup
 
    if self.isModal:
      texts = ['Commit',]
      tipTexts = ['Use edited version of data']
      commands = [self.ok,]
      closeText = 'Cancel'
      closeCmd = self.cancel
    else:
      texts = ['Delete Spectrum']
      tipTexts = ['Delete selected spectrum or spectra']
      commands = [self.deleteSpectra]
      closeText = None
      closeCmd = None   
      
    self.bottomButtons = UtilityButtonList(self.tabbedFrame.sideFrame,
                                           texts=texts,
                                           tipTexts=tipTexts,
                                           commands=commands,
                                           closeText=closeText,
                                           closeCmd=closeCmd,
                                           helpUrl=self.help_url)
                                                    
    self.bottomButtons.grid(row=0, column=0, sticky='e')
    self.bottomButtons.buttons[0].config(bg='#E0D0C0')
    
    self.updateAfter()
    self.administerNotifiers(self.registerNotify)

  def updateAfter(self, spectrum=None):

    if self.waiting:
      return
      
    self.waiting = True  
    self.after_idle(self.updateAll)  

  def updateAll(self):
  
    if not self.isModal:
      self.updateSpectra()
      self.updateAnalysisSpectra()
      self.updateTolSpectra()
      self.updateTolerances()
      
    self.updateRefSpectra()
    self.updateReferencing()
    self.updateFileSpectra()
    self.updateFileDetails()

    self.waiting = False


  def administerNotifiers(self, notifyFunc):
  
    #notifyFunc(self.updateFileDetailsAfter, 'ccp.general.DataLocation.NumericMatrix', '')
    #notifyFunc(self.updateAfter, 'ccp.general.DataLocation.NumericMatrix', 'setPath')
    notifyFunc(self.updateFileDetailsAfter, 'ccp.nmr.Nmr.DataSource', 'setDataStore')
    
    for clazz in ('ccp.nmr.Nmr.FreqDataDim',
                   'ccp.nmr.Nmr.DataDimRef',
                   'ccp.nmr.Nmr.ExpDimRef',
                   'ccp.nmr.Nmr.DimensionScaling'):
      notifyFunc(self.updateReferencingAfter, clazz, '__init__')
      notifyFunc(self.updateReferencingAfter, clazz, '')
      notifyFunc(self.updateReferencingAfter, clazz, 'delete')

    for clazz in ('FreqDataDim', 'FidDataDim', 'SampledDataDim'):
      notifyFunc(self.updateDataDim, 'ccp.nmr.Nmr.%s' % clazz, '__init__')
      notifyFunc(self.updateDataDim, 'ccp.nmr.Nmr.%s' % clazz, '')
    
    notifyFunc(self.updateDataStore, 'ccp.general.DataLocation.BlockedBinaryMatrix', '')
    
    if not self.isModal:

      for func in ('__init__', 'delete', 'setName', 'setScale',
                   'setNoiseLevel', 'setActivePeakList', 'setDataStore'):
        notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.DataSource', func)

      for func in ('__init__', 'delete'):
        notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.PeakList', func)

      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.Experiment', 'setName')

      # TBD Any notifiers from app data?

      for func in ('__init__', 'delete'):
        notifyFunc(self.updateAfter, 'ccpnmr.Analysis.StoredContour', func)

      for func in ('setChemShiftWeight','setAssignTolerance'):
        notifyFunc(self.updateTolerancesAfter, 'ccpnmr.Analysis.AnalysisDataDim', func)

      for func in ('setIsotopeCodes',):
        notifyFunc(self.updateTolerancesAfter, 'ccp.nmr.Nmr.ExpDimRef', func)

      for func in ('setNegLevels','setPosLevels',
                   'setRank','setShortcut',
                   'setSliceColor','setUseBoundingBox',
                   'setUsePeakArrow','setUsePrecalculated',
                   'setNegColors','setPosColors', 'setFont'):
        notifyFunc(self.updateAnalysisSpectraAfter, 'ccpnmr.Analysis.AnalysisSpectrum', func)
     
      self.dataLocationFrame.administerNotifiers(notifyFunc)
 
  def updateDataStore(self, dataStore):
  
    if self.waiting:
      return
  
    if dataStore is not self.fileSpectrum.dataStore:
      return
  
    self.updateFileDetails(self.fileSpectrum)
  
  def updateDataDim(self, dataDim):
  
    if self.waiting:
      return
  
    self.updateFileDetails(dataDim.dataSource)

  
  def open(self):
  
    self.updateAfter()
    BasePopup.open(self)

  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)

  def getSpectra(self):
  
    spectra = []
    
    if self.nmrProject:
      #for experiment in self.nmrProject.sortedExperiments():
      #  spectra += experiment.sortedDataSources()
      experiments = sorted(self.nmrProject.experiments, key=lambda experiment: experiment.name)
      for experiment in experiments:
        spectra.extend(sorted(experiment.dataSources, key=lambda dataSource: dataSource.name))
    
    return spectra

  def getAnalysisSpectra(self):
  
    spectra = []
    
    if self.nmrProject:
      for spectrum in self.getSpectra():
        analysisSpec = spectrum.analysisSpectrum
        
        if analysisSpec:
          spectra.append(analysisSpec)
    
    return spectra    

  # Main spectrum functions
  
  def updateSpectra(self):
    
    textMatrix = []
    spectra = self.getSpectra()
    
    for i, spec in enumerate(spectra):
    
      activePl = None
      if spec.activePeakList:
        activePl = spec.activePeakList.serial
    
      path = None
      if spec.dataStore:
        path = spec.dataStore.path
      
      isotopes = []
      for dataDim in spec.sortedDataDims():
        isotopeStr = '-'
        if dataDim.className == 'FreqDataDim':
          dataDimRef = getPrimaryDataDimRef(dataDim)
          if dataDimRef is not None:
            isotopeStr = '/'.join(dataDimRef.expDimRef.isotopeCodes)
    
        isotopes.append(isotopeStr)
    
      dimText = '%s: %s' % (spec.numDim, ','.join(isotopes)) 
      datum = [i+1,
               spec.experiment.name,
               spec.name,
               dimText,
               activePl,
               len(spec.peakLists),
               spec.scale,
               self.getSpectrumNoise(spec),
               spec.dataType,
               spec.isSimulated and 'Yes' or 'No',
               path]
               
      textMatrix.append(datum)
    
    self.spectrumTable.update(textMatrix=textMatrix,
                              objectList=spectra)
  
  def updateSpectrumButtons(self, index=None):
    
    if index is None:
      index = self.tabbedFrame.selected
    
    selectedObjDict = {0:self.spectrum,
                       1:self.analysisSpec,
                       2:self.tolSpectrum,
                       3:self.refSpectrum,
                       4:self.fileSpectrum,
                       5:True} # Too much of a pain to delve
    
    if selectedObjDict[index]:
      self.bottomButtons.buttons[0].enable()
    else:
      self.bottomButtons.buttons[0].disable()
  
  def selectSpectrum(self, obj, row, col):
  
    self.spectrum = obj
    self.updateSpectrumButtons()
  
  def getExperimentName(self, spectrum):

    width = max(20, len(spectrum.experiment.name))
    self.expNameEntry.config(width=width)
    self.expNameEntry.set(spectrum.experiment.name)

  def setExperimentName(self, *extra):
    
    name = self.expNameEntry.get()
    experiment = self.spectrum.experiment
    
    if experiment.name != name:
      if self.nmrProject.findFirstExperiment(name=name):
        showWarning('Failure','Name %s already in use' % name, parent=self)
        return
 
      if len(name) > 80:
        showWarning('Failure','Name too long', parent=self)
        return
      
      if len(name.splitlines()) > 1:
        showWarning('Failure','Name cannot contain line breaks', parent=self)
        return
        
      experiment.name = name
  
  def getSpectrumName(self, spectrum):

    width = max(20, len(spectrum.name))
    self.specNameEntry.config(width=width)
    self.specNameEntry.set(spectrum.name)

  def setSpectrumName(self, *extra):

    name = self.specNameEntry.get()
    spectrum = self.spectrum
    
    if spectrum.name != name:
    
      expt = spectrum.experiment
      if expt.findFirstDataSource(name=name):
        msg = 'Name %s already in use in experiment %s' % (name, expt.name)
        showWarning('Failure', msg, parent=self)
        return
 
      if len(name) > 80:
        showWarning('Failure','Name too long', parent=self)
        return
      
      if len(name.splitlines()) > 1:
        showWarning('Failure','Name cannot contain line breaks', parent=self)

      spectrum.name = name
        
  def setActivePeakList(self, null):

    peakList = self.peakListPulldown.getObject()
    self.spectrum.activePeakList = peakList
  
  def getActivePeakList(self, spectrum):

    peakList  = spectrum.activePeakList

    index = 0
    peakLists = spectrum.sortedPeakLists()
    names = ['%d' % pl.serial for pl in peakLists]

    if peakList in peakLists:
      index = peakLists.index(peakList)
    

    self.peakListPulldown.setup(names, peakLists, index)

  def getSpectrumNoise(self, spectrum):

    noise = spectrum.noiseLevel
    if noise is None:
      noise = getNoiseEstimate(spectrum)
      spectrum.noiseLevel = noise
    
    return noise
    
  def getNoise(self, spectrum):

    noise = self.getSpectrumNoise(spectrum)

    self.noiseEntry.set(noise)

  def setNoise(self, *extra):

    spectrum = self.spectrum
    noise    = self.noiseEntry.get()
    if (noise is not None) and spectrum:
      spectrum.noiseLevel = noise
      self.updateAfter()

  def getScale(self, spectrum):

    self.scaleEntry.set(spectrum.scale)

  def setScale(self, *extra):

    spectrum = self.spectrum
    scale = self.scaleEntry.get() or 0.0
    
    if scale <= 0.0:
      showError('Scale', 'Scale must be a positive value', parent=self)
      return
    
    spectrum.scale = scale

  def editReferencing(self):

    if self.refSpectrum is not self.spectrum:
      self.refSpectrum = self.spectrum
      self.updateReferencing()

    self.tabbedFrame.select(2)

  def editTolerances(self):

    if self.tolSpectrum is not self.spectrum:
      self.tolSpectrum = self.spectrum
      self.updateTolerances()

    self.tabbedFrame.select(3)
    
  def editFileDetails(self):

    if self.fileSpectrum is not self.spectrum:
      self.fileSpectrum = self.spectrum
      self.updateFileDetails()
      
    self.tabbedFrame.select(4) 

  def deleteSpectra(self, *event):

    index = self.tabbedFrame.selected
    
    spectra = []
    
    if index == 0:
      # need [:] otherwise deletion loop does not delete all spectra
      spectra = self.spectrumTable.currentObjects[:]
    
    elif index == 1:
      spectra = [ass.dataSource for ass in self.displayTable.currentObjects]
    
    elif index == 2:
      spectra = [self.refSpectrum,]
    
    elif index == 3:
      spectra = [self.tolSpectrum,]
    
    elif index == 4:
      spectra = [self.fileSpectrum,]
    
    elif index == 5:
      spectra = set()
      dataStores = self.dataLocationFrame.dataStoreTable.currentObjects
      for dataStore in dataStores:
        spectra.update(dataStore.nmrDataSources)
      spectra = list(spectra)  
    
    if not spectra:
      return
      
    if len(spectra) ==1:
      msg = 'Are you sure you want to delete spectrum "%s"?' % spectra[0].name
    else:
      names = ','.join(['%s:%s' % (s.experiment.name, s.name) for s in spectra])
      msg = 'Are you sure you want to delete %d spectra [ %s ]?' % (len(spectra), names)
          
    if showYesNo('Confirm', msg, parent=self):
      dataStores = set()
      
      for spec in spectra:
        if spec.isDeleted:
          data = (spec.experiment.name, spec.name)
          print 'WARNING: Experiment %s spectrum %s already deleted...' % data
          continue
        
        dataStores.add(spec.dataStore)
        spec.delete()
        
       
      for dataStore in dataStores:
        if dataStore and not dataStore.nmrDataSources:
          dataStore.delete()
        
  # Display functions      

  def propagateContours(self):
  
    analysisSpectra = self.displayTable.currentObjects
    
    if self.analysisSpec and len(analysisSpectra) > 1:
      spectrum = self.analysisSpec.dataSource
      eName = spectrum.experiment.name
      sName = spectrum.name
      msg = 'Propagate contour levels from %s:%s to selection?' % (eName, sName)
    
      if showYesNo('Confirm', msg, parent=self):
        posLevels = self.analysisSpec.posLevels
        negLevels = self.analysisSpec.negLevels

        for analysisSpectrum in analysisSpectra:
          if analysisSpectrum is not self.analysisSpec:
            analysisSpectrum.posLevels = posLevels
            analysisSpectrum.negLevels = negLevels

  def formatLevels(self, levels):

    if (len(levels) <= 3):
      t = ','.join(map(formatFloat, levels))
    else:
      t = formatFloat(levels[0]) + (',..[%d]..,' % (len(levels)-2)) + formatFloat(levels[-1])

    return t

  def doEditMarkExtraRules(self, analysisSpec, row, col):

    if self.displayTable.headingList[col] != 'Try Precalc\nContours':
      return True

    dataSource = analysisSpec.dataSource
    if hasattr(dataSource, 'block_file') \
       and dataSource.block_file \
       and analysisSpec.storedContours:
      return True
      
    else:
      return False
     

  def updateAnalysisSpectraAfter(self, obj=None):

    if self.waitingA:
      return
  
    else:
      self.waitingA = True
      self.after_idle(self.updateAnalysisSpectra)

  def updateAnalysisSpectra(self):
     
    self.updateAnalysisSpecButtons() 
     
    textMatrix = []
    colorMatrix = []
    analysisSpectra = self.getAnalysisSpectra()
    getColorScheme = self.analysisProfile.findFirstColorScheme
    
    for i, spec in enumerate(analysisSpectra):
      dataSource = spec.dataSource
      expt = dataSource.experiment
      
      name = '%s:%s' % (expt.name, dataSource.name)

      posColors = spec.posColors
      negColors = spec.negColors
      
      scheme = getColorScheme(colors=(spec.sliceColor,))
      if scheme:
        sliceColor = scheme.name
      else:
        sliceColor = spec.sliceColor
     
      
      font = spec.font
      
      scheme = getColorScheme(colors=posColors)
      if scheme:
        posColor = scheme.name
      else:
        posColor = 'unnamed'  
        
      scheme = getColorScheme(colors=negColors)
      if scheme:
        negColor = scheme.name  
      else:
        negColor = 'unnamed'  

      pos = self.formatLevels(spec.posLevels)
      neg = self.formatLevels(spec.negLevels)

      
      if hasattr(dataSource, 'block_file') and dataSource.block_file:
        if spec.storedContours:
          usePrecalculated = spec.usePrecalculated

        else:
          usePrecalculated = False
          
      elif spec.storedContours:
        usePrecalculated = True
        
      else:
        usePrecalculated = spec.usePrecalculated

      datum = [i+1,
               name,
               posColor,
               negColor,
               sliceColor,
               spec.useBoundingBox and 'Yes' or 'No',
               spec.shortcut,
               spec.rank,
               pos,
               neg,
               usePrecalculated and 'Yes' or 'No',
               font,
               spec.usePeakArrow and 'Yes' or 'No']
      
      colors = [None] * 13
      
      colors[4] = spec.sliceColor
               
      textMatrix.append(datum)
      colorMatrix.append(colors)
    
    self.displayTable.update(textMatrix=textMatrix,
                             objectList=analysisSpectra,
                             colorMatrix=colorMatrix)

    self.waitingA = False

  def updateAnalysisSpecButtons(self):
  
    buttons = self.displayButtons.buttons
    
    if self.analysisSpec:
      for button in buttons:
        button.enable()
    
    else:
      for button in buttons:
        button.disable()
  
  def selectAnalysisSpec(self, obj, row, col):
  
    self.analysisSpec = obj
    self.updateAnalysisSpecButtons()
    self.updateSpectrumButtons()
    
  def editLevels(self):

    if self.analysisSpec:
      spectrum = self.analysisSpec.dataSource
      self.parent.editContourLevels(spectrum)

  def editContourFiles(self):

    self.parent.editContourFiles()

  def getContourLevels(self, analysisSpec):

    if self.analysisSpec:
      spectrum = analysisSpec.dataSource
      self.parent.editContourLevels(spectrum)

  def toggleBoxVisible(self, analysisSpec):

    analysisSpec.useBoundingBox = not analysisSpec.useBoundingBox

  def togglePointerVisible(self, analysisSpec):

    analysisSpec.usePeakArrow = not analysisSpec.usePeakArrow

  def toggleTryPrecalcContours(self, analysisSpec):
 
    analysisSpec.usePrecalculated = not analysisSpec.usePrecalculated

  def getSliceColor(self, analysisSpec):

    cObjs  = Color.standardColors
    colors = [c.hex for c in cObjs]
    names  = [c.name for c in cObjs]
    index  = 0
    color  = analysisSpec.sliceColor
    
    if color in colors:
      index = colors.index(color)
    
    self.sliceColorPulldown.setup(names, colors, index, colors)

  def setSliceColor(self, color):

    self.analysisSpec.sliceColor = self.sliceColorPulldown.getObject()
    
  def getShortcut(self, analysisSpec):

    shortcut = analysisSpec.shortcut
    if shortcut is None:
      shortcut = shortcut_entries[0]
      
    self.shortcutPulldown.setSelected(shortcut)

  def setShortcut(self, *extra):

    self.analysisSpec.shortcut = self.shortcutPulldown.getText()

  def getSpectrumFont(self, analysisSpec):

    font = analysisSpec.font
    self.fontMenu.set(font)

  def setSpectrumFont(self, *extra):

    self.analysisSpec.font = self.fontMenu.getText()

  def orderRanks(self):
    
    specs = self.displayTable.objectList
    
    for i, spec in enumerate(specs):
      spec.rank = i 
  
  def getRank(self, analysisSpec):

    order = analysisSpec.rank
    self.rankEntry.set(order)

  def setRank(self, *extra):
    
    rank = self.rankEntry.get() or 0
    
    replaceSpec = self.analysisProject.findFirstAnalysisSpectrum(rank=rank)
    self.analysisSpec.rank = rank
    
    if replaceSpec and (replaceSpec is not self.analysisSpec):
      sortSpectra = [(asp.rank, asp) for asp in \
                     self.analysisProject.analysisSpectra \
                     if asp is not self.analysisSpec]
      sortSpectra.sort()
      
      analysisSpectra = [x[1] for x in sortSpectra]
      ranks = [x[0] for x in sortSpectra]
    
      index = analysisSpectra.index(replaceSpec)
      
      analysisSpectra.insert(index, self.analysisSpec)
      ranks.insert(index, rank)
      
      for i, analysisSpectrum in enumerate(analysisSpectra[:-1]):
        nextSpec = analysisSpectra[i+1]
        
        if nextSpec.rank <= analysisSpectrum.rank:
          nextSpec.rank = analysisSpectrum.rank+1 
    

  def getPosColor(self, analysisSpec):

    schemes = getHueSortedColorSchemes(self.analysisProfile)
    names = [s.name for s in schemes]
    colors = [list(s.colors) for s in schemes]
    index = 0
    
    try:
      posColors = list(analysisSpec.posColors)
      
      if posColors in colors:
        index = colors.index(posColors)
      
    except Implementation.ApiError:
      analysisSpec.posColors = ['#808080',]
      print 'Warning %s missing positive color scheme' % analysisSpec
      
    self.posColorPulldown.setup(names, schemes, index, colors)

  def setPosColor(self, *extra):

    self.analysisSpec.posColors = self.posColorPulldown.getObject().colors
    
  def getNegColor(self, analysisSpec):

    schemes = getHueSortedColorSchemes(self.analysisProfile)
    names = [s.name for s in schemes]
    colors = [list(s.colors) for s in schemes]
    index = 0
    
    try:
      negColors = list(analysisSpec.negColors)
      
      if negColors in colors:
        index = colors.index(negColors)
      
    except Exception:
      analysisSpec.negColors = ['#808080',]
      print 'Warning %s missing negative color scheme' % analysisSpec
      
    self.negColorPulldown.setup(names, schemes, index, colors)

  def setNegColor(self, scheme):

    self.analysisSpec.negColors = self.negColorPulldown.getObject().colors
    
  # Tolerance functions
  
  def selectDataDim(self, obj, row, col):
  
    self.analysisDataDim = obj
        
  def updateTolSpectra(self):  
  
    analysisSpectra = self.getAnalysisSpectra()
    spectra = [ass.dataSource for ass in analysisSpectra]
    names = ['%s:%s' % (s.experiment.name, s.name) for s in spectra]
    index = 0
    spec = self.tolSpectrum
    
    if spectra:
      if spec not in spectra:
        spec = spectra[0]
    
      index = spectra.index(spec)
    
    else:
      spec = None
    
    if spec is not self.tolSpectrum:
      self.tolSpectrum = spec
      self.updateTolerances()
  
    self.tolSpecPulldown.setup(names, spectra, index)
  
  def selectTolSpectrum(self, spectrum, row=None, col=None):
  
    if self.tolSpectrum is not spectrum:
      self.tolSpectrum = spectrum
      self.updateTolSpectra()  
      self.updateTolerances()
      self.updateSpectrumButtons()  
      
  def updateTolerancesAfter(self, obj=None):
  
    if self.waitingT:
      return
    
    self.waitingT = True
    self.after_idle(self.updateTolSpectra)  
    self.after_idle(self.updateTolerances)

  def updateTolerances(self):
  
    textMatrix = []
    analysisDataDims = []
    
    if self.tolSpectrum:
      analysisDataDims = self.tolSpectrum.analysisSpectrum.analysisDataDims
      analysisDataDims = [(add.dataDim.dim, add) for add in analysisDataDims]
      analysisDataDims.sort()
      analysisDataDims = [x[1] for x in analysisDataDims]
    
      for add in analysisDataDims:
        dataDim = add.dataDim
        if dataDim.className == 'SampledDataDim':
          continue

        dataDimRef = getPrimaryDataDimRef(dataDim)
        isotopes = '/'.join(dataDimRef.expDimRef.isotopeCodes)
        datum = [dataDim.dim,
                 isotopes,
                 add.assignTolerance,
                 add.chemShiftWeight]

        textMatrix.append(datum)
 

    self.tolTable.update(textMatrix=textMatrix,
                         objectList=analysisDataDims)   
    
    textMatrix = []
    objectList = []
    
    if self.nmrProject:
      i = 1
      for experiment in self.nmrProject.sortedExperiments():
        eName = experiment.name
        
        for spectrum in experiment.sortedDataSources():
          name = '%s:%s' % (eName, spectrum.name)
          
          dataDims = [dd for dd in spectrum.sortedDataDims() if dd.analysisDataDim]
          analysisDataDims = [dd.analysisDataDim for dd in dataDims]
          
          tols = []
          weights = []
          
          for j, dataDim in enumerate(dataDims):
            if dataDim.className == 'SampledDataDim':
              continue
          
            analysisDataDim = analysisDataDims[j]
            dataDimRef = getPrimaryDataDimRef(dataDim)
            isotopes = '/'.join(dataDimRef.expDimRef.isotopeCodes)            
            
            tols.append('%s:%.2f' % (isotopes,analysisDataDim.assignTolerance))
            weights.append('%s:%.2f' % (isotopes,analysisDataDim.chemShiftWeight))
          
          datum = [i,
                   name,
                   ' '.join(tols),
                   ' '.join(weights)]
          
          textMatrix.append(datum)
          objectList.append(spectrum)
          i += 1
    
    self.tolTableAll.update(textMatrix=textMatrix,
                            objectList=objectList)   
    
    if self.tolTableAll.currentObject is not self.tolSpectrum:
      self.tolTableAll.currentObject = self.tolSpectrum
      #self.tolTableAll.currentObjects = []
      self.tolTableAll.clearHighlights()
      
      if self.tolSpectrum:
        self.tolTableAll.hilightObject(self.tolSpectrum)
    
    self.waitingT = False
                               
  def getShiftWeight(self, analysisDataDim):

    self.shiftWeightEntry.set(analysisDataDim.chemShiftWeight)

  def setShiftWeight(self, *extra):

    analysisDataDim = self.analysisDataDim
    
    value = self.shiftWeightEntry.get()
    if value is None:
      value = analysisDataDim.chemShiftWeight
      
    value = min(100.0, max(value,0.0))
    
    analysisDataDim.chemShiftWeight = value

  def getAssignmentTol(self, analysisDataDim):

    self.assignTolEntry.set(analysisDataDim.assignTolerance)

  def setAssignmentTol(self, *extra):

    analysisDataDim = self.analysisDataDim
    value = self.assignTolEntry.get()
    
    if value is None:
      value = analysisDataDim.assignTolerance
      
    analysisDataDim.assignTolerance = value 
  
  def propagateAssignTols(self):
  
    primary = self.tolTableAll.currentObject
    
    if not primary:
      return
    
    others = self.tolTableAll.currentObjects
    if not others or (primary not in others):
      return

    others.remove(primary)
    
    if not others:
      return
    
    specName = '%s:%s' % (primary.experiment.name, primary.name)
      
    msg = 'Propagate assignment tolerances '
    msg += 'from spectrum %s to %d others?' % (specName, len(others))
  
    if showOkCancel('Confirm', msg, parent=self):
      
      isotopeDict = {}
      dimDict = {}
      for dataDim in primary.sortedDataDims():
        if dataDim.className == 'SampledDataDim':
          continue
          
        add = dataDim.analysisDataDim
        if not add:
          continue
        
        dataDimRef = getPrimaryDataDimRef(dataDim)
        isotopes = '/'.join(dataDimRef.expDimRef.isotopeCodes)
        
        if not isotopeDict.has_key(isotopes):
          isotopeDict[isotopes] = add.assignTolerance
    
        dimDict[dataDim.dim] = (isotopes, add.assignTolerance)
    
      for spectrum in others:
        for dataDim in spectrum.dataDims:
          if dataDim.className == 'SampledDataDim':
            continue
 
          add = dataDim.analysisDataDim
          if not add:
            continue
          
          dataDimRef = getPrimaryDataDimRef(dataDim)
          isotopes = '/'.join(dataDimRef.expDimRef.isotopeCodes)
          
          isotopesA, value = dimDict.get(dataDim.dim, (None, None)) 
          
          if isotopesA == isotopes:
            add.assignTolerance = value
          
          else:
            value = isotopeDict.get(isotopes)
            if value is not None:
              add.assignTolerance = value
  
  def propagateShiftWeights(self):
  
    primary = self.tolTableAll.currentObject
    
    if not primary:
      return
    
    others = self.tolTableAll.currentObjects
    if not others or (primary not in others):
      return

    others.remove(primary)
    
    specName = '%s:%s' % (primary.experiment.name, primary.name)
      
    msg = 'Propagate chemical shift weights '
    msg += 'from spectrum %s to %d others?' % (specName, len(others))
    
    if showOkCancel('Confirm', msg, parent=self):
      
      isotopeDict = {}
      dimDict = {}
      for dataDim in primary.sortedDataDims():
        if dataDim.className == 'SampledDataDim':
          continue
          
        add = dataDim.analysisDataDim
        if not add:
          continue
        
        dataDimRef = getPrimaryDataDimRef(dataDim)
        isotopes = '/'.join(dataDimRef.expDimRef.isotopeCodes)
        
        if not isotopeDict.has_key(isotopes):
          isotopeDict[isotopes] = add.chemShiftWeight
    
        dimDict[dataDim.dim] = (isotopes, add.chemShiftWeight)
    
      for spectrum in others:
        for dataDim in spectrum.dataDims:
          if dataDim.className == 'SampledDataDim':
            continue
 
          add = dataDim.analysisDataDim
          if not add:
            continue
          
          dataDimRef = getPrimaryDataDimRef(dataDim)
          isotopes = '/'.join(dataDimRef.expDimRef.isotopeCodes)
          
          isotopesA, value = dimDict.get(dataDim.dim, (None, None)) 
          
          if isotopesA == isotopes:
            add.chemShiftWeight = value
          
          else:
            value = isotopeDict.get(isotopes)
            if value is not None:
              add.chemShiftWeight = value
  
  # Referencing functions
  
  def updateReferencingAfter(self, dim=None):

    if self.waitingR:
      return
    
    # Filter on dims parents
    
    self.waitingR = True
    self.updateRefSpectra()
    self.updateReferencing()

  def updateRefSpectra(self):
  
    index = 0
    spectra = self.getSpectra()
    names = ['%s:%s' % (s.experiment.name,s.name) for s in spectra]
    spec = self.refSpectrum
    
    if spectra:
      if spec not in spectra:
        spec = spectra[0]
      
      index = spectra.index(spec)
    
    else:
      spec = None
      
    if spec is not self.refSpectrum:
      self.refSpectrum = spec
      self.updateReferencing()  
  
    self.refSpecPulldown.setup(names, spectra, index)
    
  def updateReferencing(self, spectrum=None):

    self.updateRefButtons()
  
    if spectrum and (spectrum is not self.refSpectrum):
      return

    if self.refSpectrum:
      dataDims = self.refSpectrum.sortedDataDims()
    else:
      dataDims = []

    textMatrix = []
    objectList = []
    for dataDim in dataDims:
      dim = dataDim.dim
    
      if isinstance(dataDim, Nmr.FreqDataDim):
        
        # TBD: more general - FIDs ?
        
        ll = [(x.expDimRef.serial, x) for x in dataDim.dataDimRefs]
        ll.sort()
        dataDimRefs = [x[1] for x in ll]

        if len(dataDimRefs) > 1:
          multiRef = True
        else:
          multiRef = False
        
        for dataDimRef in dataDimRefs:
          expDimRef  = dataDimRef.expDimRef
          isotopes = '/'.join(expDimRef.isotopeCodes)
 
          if multiRef:
            label = '%d (%d)' % (dim, expDimRef.serial)
          else:
            label = dim
          
          # get sw etc. NB cannot use DataDimRef.valuePerPoint 
          # as this breaks when sf not set
          swhz = ((dataDimRef.localValuePerPoint or dataDim.valuePerPoint)
                  * dataDim.numPoints)
          sf = expDimRef.sf
          if sf:
            swppm = swhz / sf
          else:
            swppm = None
 
          if self.useReducedDim:
          
            subDimName = expDimRef.displayName
            if not subDimName:
              if expDimRef.refExpDimRef:
                measurement = expDimRef.refExpDimRef.expMeasurement
                subDimName = ','.join([ass.name for ass in measurement.atomSites])
                expDimRef.displayName = subDimName
            
            dimensionScaling = dataDim.findFirstDimensionScaling(expDimRef=expDimRef)
            if dimensionScaling:
              scaleFactors = ','.join(['%.3f' % v for v in dimensionScaling.scalingFactors])
            else:  
              scaleFactors = None
            datum = [label,
                     subDimName,
                     isotopes,
                     scaleFactors,
                     sf,
                     swppm,
                     swhz,
                     dataDimRef.refValue,
                     dataDimRef.refPoint,
                     dataDim.numPointsOrig,
                     dataDim.pointOffset,
                     expDimRef.minAliasedFreq,
                     expDimRef.maxAliasedFreq]
          
          else:

            swppm = dataDimRef.valuePerPoint * dataDim.numPoints
            datum = [label,
                     isotopes,
                     sf,
                     swppm,
                     swhz,
                     dataDimRef.refValue,
                     dataDimRef.refPoint,
                     dataDim.numPointsOrig,
                     dataDim.pointOffset,
                     expDimRef.minAliasedFreq,
                     expDimRef.maxAliasedFreq]
 
 
          textMatrix.append(datum)
          objectList.append(dataDimRef)
 
    if self.useReducedDim:
      headingList = reducedDimHeadingList
      editWidgets = [ None, self.subDimNameEntry,
                      self.isotopeMenu, self.scaleFactorsWidget, self.sfEntry, None,
                      self.swEntry, self.refppmEntry, self.refptEntry, self.origNptsEntry,
                      self.pointOffsetEntry, self.minFreqEntry, self.maxFreqEntry ]
      editGetCallbacks = [ None, self.getSubDimName,
                           self.getIsotope, self.getDimScaling, self.getSf, None,
                           self.getSw, self.getRefppm, self.getRefpt, self.getOrigNpts,
                           self.getPointOffset, self.getMinFreq, self.getMaxFreq ]
      editSetCallbacks = [ None, self.setSubDimName,
                           self.setIsotope, None, self.setSf, None, self.setSw,
                           self.setRefppm, self.setRefpt, self.setOrigNpts,
                           self.setPointOffset, self.setMinFreq, self.setMaxFreq ]
      tipTexts = reducedTipTexts
      
    else:
      headingList = standardHeadingList
      editWidgets = [ None, self.isotopeMenu, self.sfEntry, None, self.swEntry,
                      self.refppmEntry, self.refptEntry, self.origNptsEntry,
                      self.pointOffsetEntry, self.minFreqEntry, self.maxFreqEntry ]
      editGetCallbacks = [ None, self.getIsotope, self.getSf, None, self.getSw,
                           self.getRefppm, self.getRefpt, self.getOrigNpts,
                           self.getPointOffset, self.getMinFreq, self.getMaxFreq ]
      editSetCallbacks = [ None, self.setIsotope, self.setSf, None, self.setSw,
                           self.setRefppm, self.setRefpt, self.setOrigNpts,
                           self.setPointOffset, self.setMinFreq, self.setMaxFreq ]
      tipTexts = standardTipTexts
 
    self.dataDimRefTable.update(objectList=objectList,
                                textMatrix=textMatrix,
                                tipTexts=tipTexts,
                                headingList=headingList,
                                editWidgets=editWidgets,
                                editGetCallbacks=editGetCallbacks,
                                editSetCallbacks=editSetCallbacks)
    self.waitingR = False

  def guessIsotopeCode(self, sf, ratio):

    for key in chemShiftRefRatios.keys():
      r = sf / chemShiftRefRatios[key]
      if abs((r-ratio)/ratio) < 0.01:
        return key

    return None
    
  def selectRefSpectrum(self, spectrum):
  
    if self.refSpectrum is not spectrum:
      self.refSpectrum = spectrum
      self.updateReferencing()
      self.updateSpectrumButtons()  
      
  def selectDataDimRef(self, obj, row, col):
  
    self.dataDimRef = obj
    self.updateRefButtons()
   
  def updateRefButtons(self):
    
    ddr = self.dataDimRef
    
    if ddr and self.useReducedDim:
      nDdr = len(ddr.dataDim.dataDimRefs)
      if ((self.isModal and len(self.spectrum.experiment.dataSources) == 1) 
          or (nDdr < len(ddr.expDimRef.expDim.expDimRefs))):
        self.addSubDimButton.enable()
      
      if nDdr > 1:
        self.removeSubDimButton.enable()
      else:
        self.removeSubDimButton.disable()
        
    else:
      self.addSubDimButton.disable()
      self.removeSubDimButton.disable()
  
  def addSubDim(self):
    """Add sub-dimension that copies selected dataDimRef in so far as possible
    """
  
    if self.dataDimRef:
      dataDimRef = self.dataDimRef
      self.useReducedDim = True
      self.reducedDimSelect.set(self.useReducedDim)
  
      dataDim = dataDimRef.dataDim
      edr = dataDimRef.expDimRef
      expDim = edr.expDim
      
      # no existing spectrum for that experiment
      # we are making new experiment 
      if self.isModal and len(self.spectrum.experiment.dataSources) == 1:
        useExpDimRef  = expDim.newExpDimRef(sf=edr.sf, unit=edr.unit, 
                                            isotopeCodes=edr.isotopeCodes,
                                            measurementType=edr.measurementType,
                                            baseFrequency=edr.baseFrequency)
        prevDataDimRef = dataDimRef
      else:
        # working in existing experiment
        inUse = set(x.expDimRef for x in dataDim.dataDimRefs)
        useExpDimRef = None
        for newExpDimRef in expDim.sortedExpDimRefs():
          if newExpDimRef not in inUse:
            if newExpDimRef.isotopeCodes == edr.isotopeCodes:
              useExpDimRef = newExpDimRef
              prevDataDimRef = dataDimRef
              break
            elif useExpDimRef is None:
              useExpDimRef = newExpDimRef
            
        else:
          if useExpDimRef is None:
            return
          else:
            prevDataDimRef = useExpDimRef.findFirstDataDimRef()
      
      # make new DataDimRef and set previous referencing, if any
      newDataDimRef = dataDim.newDataDimRef(expDimRef=useExpDimRef)
      if prevDataDimRef:
        for tag in ('refPoint', 'refValue', 'localValuePerPoint'):
        #for tag in ('refPoint', 'refValue'):
          setattr(newDataDimRef, tag, getattr(prevDataDimRef, tag))
  
  def removeSubDim(self):
    
    dataDimRef = self.dataDimRef
    if dataDimRef:
      expDimRef = dataDimRef.expDimRef
      dataDimRef.delete()
      if self.isModal and len(self.spectrum.experiment.dataSources) == 1:
        expDimRef.delete()
   
  def toggleReducedDim(self, *null):
  
    self.useReducedDim = not self.useReducedDim
    self.updateReferencing()  


  def getDimScaling(self, dataDimRef):
  
    dataDim   = dataDimRef.dataDim
    expDimRef = dataDimRef.expDimRef
  
    dimensionScaling = dataDim.findFirstDimensionScaling(expDimRef=expDimRef)

    if dimensionScaling:
      values = dimensionScaling.scalingFactors
    else:  
      values = None
     
    self.scaleFactorsWidget.set(values, options=None)
  
  def setDimScaling(self, *event):

    values = self.scaleFactorsWidget.get() 
    self.dataDimRefTable.keyPressEscape()

    dataDim   = self.dataDimRef.dataDim
    expDimRef = self.dataDimRef.expDimRef
    dimensionScaling = dataDim.findFirstDimensionScaling(expDimRef=expDimRef)
    
    if dimensionScaling and values:
      dimensionScaling.scalingFactors = values
      
    elif dimensionScaling:
      dimensionScaling.delete()
    
    elif values:
      dimensionScaling = dataDim.newDimensionScaling(expDimRef=expDimRef,
                                                     scalingFactors=values)
                                                     

  def getSubDimName(self, dataDimRef):
  
    self.subDimNameEntry.set(dataDimRef.expDimRef.displayName)
  
  def setSubDimName(self, *event):

    if self.dataDimRef:
      name = self.subDimNameEntry.get() or None
      self.dataDimRef.expDimRef.displayName = name

  def getIsotopeFunc(self, dataDimRef):
    isotopeCodes = dataDimRef.expDimRef.isotopeCodes
    if isotopeCodes:
      self.isotopeMenu.set(isotopeCodes[0])
 
  def setIsotopeFunc(self, *extra):

    isotopeCode = self.isotopeMenu.getObject()
    try:
      oldIsotopeCodes = self.dataDimRef.expDimRef.isotopeCodes
      self.dataDimRef.expDimRef.isotopeCodes = (isotopeCode,)
      
    except Implementation.ApiError, e:
      showError('Setting isotope codes', e.error_msg, parent=self)
    
    # reset SF if changing isotope
    if oldIsotopeCodes != (isotopeCode,):
      ratio = None
      newGr = chemShiftRefRatios.get(isotopeCode)
      if newGr:
        for dataDimRef in self.dataDimRefTable.objectList:
          if dataDimRef is not self.dataDimRef:
            expDimRef = dataDimRef.expDimRef
            sf = expDimRef.sf
            gr = chemShiftRefRatios.get('/'.join(expDimRef.isotopeCodes))
            
            if sf and gr:
              ratio = sf/gr
              break

        if ratio is not None:
          self.dataDimRef.expDimRef.sf = ratio * newGr
    
    
  def getSf(self, dataDimRef):

    self.sfEntry.set(dataDimRef.expDimRef.sf)

  def setSf(self, *extra):

    dataDimRef = self.dataDimRef
    
    if not self.isModal and self.peakRefSelect.getIndex():
      oldReference = getReference(dataDimRef)
      shiftPeaks = True
    else:
      shiftPeaks = False

    try:
      dataDimRef.expDimRef.sf = self.sfEntry.get()
    except Implementation.ApiError, e:
      showError('Setting spectrometer frequency', e.error_msg, parent=self)

    if shiftPeaks:
      try:
        shiftDataDimRef(dataDimRef, oldReference)
      except Implementation.ApiError, e:
        showError('Shifting peaks', e.error_msg, parent=self)

  def getSw(self, dataDimRef):

    dataDim = dataDimRef.dataDim
    swhz = ((dataDimRef.localValuePerPoint or dataDim.valuePerPoint)
            * dataDim.numPoints)
    self.swEntry.set(swhz)
    #dnp = dataDim.numPoints
    #self.swEntry.set("%s - %s" % (dataDimRef.localValuePerPoint, dataDim.valuePerPoint))

  def setSw(self, *extra):

    dataDimRef = self.dataDimRef
    dataDim = dataDimRef.dataDim

    if not self.isModal and self.peakRefSelect.getIndex():
      oldReference = getReference(dataDimRef)
      shiftPeaks = True
    else:
      shiftPeaks = False
    
    try:
      hzPerPoint = self.swEntry.get() / dataDim.numPoints 
      ll = [(x.expDimRef.serial, x) for x in dataDim.dataDimRefs]
      dataDimRefs = [x[1] for x in sorted(ll)]
      
      if dataDimRef is dataDimRefs[0]:
        dataDimRef.localValuePerPoint = None
        dataDim.valuePerPoint = hzPerPoint
      else:
        dataDimRef.localValuePerPoint = hzPerPoint
        
    except Implementation.ApiError, e:
      showError('Setting value per point', e.error_msg, parent=self)
    
    if shiftPeaks:
      try:
        shiftDataDimRef(dataDimRef, oldReference)
      except Implementation.ApiError, e:
        showError('Shifting peaks', e.error_msg, parent=self)

  def getRefppm(self, dataDimRef):

    self.refppmEntry.set(dataDimRef.refValue)

  def setRefppm(self, *extra):

    dataDimRef = self.dataDimRef
    dataDim = dataDimRef.dataDim

    if not self.isModal and self.peakRefSelect.getIndex():
      oldReference = getReference(dataDimRef)
      shiftPeaks = True
    else:
      shiftPeaks = False

    try:
      dataDimRef.refValue = self.refppmEntry.get()
    except Implementation.ApiError, e:
      showError('Setting reference value', e.error_msg, parent=self)

    if shiftPeaks:
      try:
        shiftDataDimRef(dataDimRef, oldReference)
      except Implementation.ApiError, e:
        showError('Shifting peaks', e.error_msg, parent=self)

  def getRefpt(self, dataDimRef):

    self.refptEntry.set(dataDimRef.refPoint)

  def setRefpt(self, *extra):

    dataDimRef = self.dataDimRef
    dataDim = dataDimRef.dataDim
    
    if not self.isModal and self.peakRefSelect.getIndex():
      oldReference = getReference(dataDimRef)
      shiftPeaks = True
    else:
      shiftPeaks = False

    try:
      dataDimRef.refPoint = self.refptEntry.get()
    except Implementation.ApiError, e:
      showError('Setting reference point', e.error_msg, parent=self)

    if shiftPeaks:
      try:
        shiftDataDimRef(dataDimRef, oldReference)
      except Implementation.ApiError, e:
        showError('Shifting peaks', e.error_msg, parent=self)

  def getOrigNpts(self, dataDimRef):

    dataDim = dataDimRef.dataDim
    self.origNptsEntry.set(dataDim.numPointsOrig)

  def setOrigNpts(self, *extra):

    dataDimRef = self.dataDimRef
    dataDim = dataDimRef.dataDim
    
    if not self.isModal and self.peakRefSelect.getIndex():
      oldReference = getReference(dataDimRef)
      shiftPeaks = True
    else:
      shiftPeaks = False

    try:
      dataDim.numPointsOrig = self.origNptsEntry.get()
    except Implementation.ApiError, e:
      showError('Setting original number of points', e.error_msg, parent=self)

    if shiftPeaks:
      try:
        shiftDataDimRef(dataDimRef, oldReference)
      except Implementation.ApiError, e:
        showError('Shifting peaks', e.error_msg, parent=self)

  def getPointOffset(self, dataDimRef):

    dataDim = dataDimRef.dataDim
    self.pointOffsetEntry.set(dataDim.pointOffset)

  def setPointOffset(self, *extra):

    dataDimRef = self.dataDimRef
    dataDim = dataDimRef.dataDim
    
    if not self.isModal and self.peakRefSelect.getIndex():
      oldReference = getReference(dataDimRef)
      shiftPeaks = True
    else:
      shiftPeaks = False

    try:
      dataDim.pointOffset = self.pointOffsetEntry.get()
    except Implementation.ApiError, e:
      showError('Setting point offset', e.error_msg, parent=self)

    if shiftPeaks:
      try:
        shiftDataDimRef(dataDimRef, oldReference)
      except Implementation.ApiError, e:
        showError('Shifting peaks', e.error_msg, parent=self)

  def getMinFreq(self, dataDimRef):

    self.minFreqEntry.set(dataDimRef.expDimRef.minAliasedFreq)

  def setMinFreq(self, *extra):

    dataDimRef = self.dataDimRef
    try:
      dataDimRef.expDimRef.minAliasedFreq = self.minFreqEntry.get()
      
    except Implementation.ApiError, e:
      showError('Setting min aliased frequency', e.error_msg, parent=self)

  def getMaxFreq(self, dataDimRef):

    self.maxFreqEntry.set(dataDimRef.expDimRef.maxAliasedFreq)

  def setMaxFreq(self, *extra):

    dataDimRef = self.dataDimRef
    try:
      dataDimRef.expDimRef.maxAliasedFreq = self.maxFreqEntry.get()
      
    except Implementation.ApiError, e:
      showError('Setting max aliased frequency', e.error_msg, parent=self)

  # below is only called if isModal, in which case ok function calls this
  def apply(self):

    if not self.commitFileDetails():
      return False

    spectrum = self.refSpectrum
    
    if spectrum:
      ratio = None
      dataDimRefs = self.dataDimRefTable.objectList
      
      for dataDimRef in dataDimRefs:
        dataDim = dataDimRef.dataDim
        expDimRef = dataDimRef.expDimRef
        isotopeCode = '/'.join(expDimRef.isotopeCodes)
        
        dimId = expDimRef.displayName
        if not dimId:
          if len(dataDim.dataDimRefs) > 1:
            dimId = '%s (%s)' % (dataDim.dim, expDimRef.serial)
          else:
            dimId = '%d' % dataDim.dim
        
        known = chemShiftRefRatios.has_key(isotopeCode)
        if ratio:
          possibleIsotopeCode = self.guessIsotopeCode(expDimRef.sf, ratio)
          
          if known and not expDimRef.sf:
            # sf is not set, offer to set it
            newSf = ratio * chemShiftRefRatios[isotopeCode]
            msg = ('For %s in dim %s, Spectrometer Frequency ought to be %s. Reset it?' 
                   % (isotopeCode, dimId, newSf))
            if showYesNo('Reset SF?', msg, parent=self):
              expDimRef.sf = newSf
              return False
          
          if (known or possibleIsotopeCode) and (isotopeCode != possibleIsotopeCode):
            msg = ('Given sf, the isotope (%s) looks incorrect for dim %s' 
                   % (isotopeCode, dimId))
            if possibleIsotopeCode:
              msg = msg + ' (perhaps it should be %s)' % possibleIsotopeCode
            msg = msg + ': continue as is?'
            
            if not showYesNo('Isotope ok?', msg, parent=self):
              return False
              
        elif known:
          ratio = expDimRef.sf / chemShiftRefRatios[isotopeCode]

    return True
    
  # File details functions
  
  def selectFileDim(self, obj, row, col):
  
    self.fileDim = obj

  def updateFileDetailsAfter(self, spec):

    if self.waitingF:
      return
    
    self.waitingF = True
    
    if spec.className == 'DataSource':
      self.updateFileSpectra()
      self.updateFileDetails(spec)
      
    elif self.fileSpectrum in spec.nmrDataSources:
      self.updateFileSpectra()
      self.updateFileDetails()
  
  def updateFileSpectra(self):
  
    index = 0
    spectra = self.getSpectra()
    names = ['%s:%s' % (s.experiment.name,s.name) for s in spectra]
    spec = self.fileSpectrum
    
    if spectra:
      if spec not in spectra:
        spec = spectra[0]
      
      index = spectra.index(spec)
    
    else:
      spec = None
      
    if spec is not self.fileSpectrum:
      self.fileSpectrum = spec
      self.updateFileDetails()  
  
    self.fileSpecPulldown.setup(names, spectra, index)
  
  def updateFileDetails(self, spectrum=None):
    
    if spectrum and (spectrum is not self.fileSpectrum):
      return
    
    textMatrix = []
    spectrum = self.fileSpectrum

    isShape = spectrum and isShapeSpectrum(spectrum)

    # Remove unwanted options for factorised spectra
    
    if spectrum and isShape:
      self.dataOptFrame.grid_forget()
    
    else:
      self.dataOptFrame.grid(row=1, column=0, columnspan=2, sticky='ew')
    
    if spectrum:
      # Update table
      
      dataDims = spectrum.sortedDataDims()
      
      # get block sizes
      dataStore = spectrum.dataStore
      if isShape:
        blockSizes = ['']*len(dataDims)
      else:
        blockSizes = None
        if dataStore:
          blockSizes = dataStore.blockSizes
        
        if not blockSizes:
          blockSizes = [0]*len(dataDims)
      
      # update matrix contents
      for ii,dataDim in enumerate(dataDims):
        datum = [dataDim.dim,
                 dataDim.numPoints,
                 blockSizes[ii]]
                 
        textMatrix.append(datum)
        
    else:
      dataDims = []

    self.fileDimTable.update(objectList=dataDims,
                        textMatrix=textMatrix)

    # Set properties

    if spectrum:
      self.specFileEntry.set(getDataSourceFileName(spectrum))
    else:
      self.specFileEntry.set('')

    if spectrum and spectrum.dataStore and not isShape:
      dataStore = spectrum.dataStore
      
      self.endianSelect.set(dataStore.isBigEndian)
      
      nByte = dataStore.nByte
      
      if nByte in number_bytes_entries:
        self.numBytesPulldown.setSelected(str(nByte))
        
      else: # unlikely
        entries = map(str, number_bytes_entries+[nByte])
        self.numBytesPulldown.setup(entries, entries, len(self.numBytesPulldown.objects))

      if not dataStore.numberType:
        dataStore.numberType = 'float'
        
      if dataStore.numberType in number_type_entries:
        self.numTypePulldown.setSelected(dataStore.numberType)
        
      else:
        entries = list(number_type_entries)+[dataStore.numberType]
        self.numTypePulldown.setup(entries, entries, len(self.numTypePulldown.objects))

      self.headerSizeEntry.set(dataStore.headerSize)

    else:

      self.endianSelect.set(True)
      self.numBytesPulldown.setSelected(str(4))
      self.numTypePulldown.setSelected('float')
      self.headerSizeEntry.set(0)
   
    if spectrum and not isShape:
      self.fileDetailsButtons.buttons[0].enable()
    else:
      self.fileDetailsButtons.buttons[0].disable()
   
    self.waitingF = False
 
  def selectFileSpectrum(self, spectrum):
  
    if self.fileSpectrum is not spectrum:
      self.fileSpectrum = spectrum
      self.updateFileDetails()
      self.updateSpectrumButtons()  
  
  def browseSpecFile(self):

    fileName = self.specFileEntry.get()
    if fileName:
      directory = os.path.dirname(fileName)
      if not os.path.exists(directory):
        directory = os.getcwd()
    else:
      directory = os.getcwd()

    popup = FileSelectPopup(self, directory=directory)
    fileName = popup.getFile()
    popup.destroy()
    
    if fileName:
      self.specFileEntry.set(fileName)

  def getNpoints(self, dataDim):

    self.npointsEntry.set(dataDim.numPoints)

  def setNpoints(self, *extra):

    dataDim = self.fileDim
    
    npoints = self.npointsEntry.get()
    try:
      dataDim.numPoints = npoints
    except Implementation.ApiError, e:
      showError('Setting number of points', e.error_msg, parent=self)

  def getBlockSize(self, dataDim):

    dataStore = dataDim.dataSource.dataStore

    if dataStore and (dataStore.className == 'BlockedBinaryMatrix'):
      i = dataDim.dim-1
      blockSizes = dataStore.blockSizes
      try:
        blockSize = blockSizes[i]
      except: # blockSizes might not have been set
        blockSize = 1
      self.blockSizeEntry.set(blockSize)

  def setBlockSize(self, *extra):

    blockSize = self.blockSizeEntry.get()
    dataDim = self.fileDim
    
    dataSource = dataDim.dataSource
    dataStore = dataSource.dataStore

    if dataStore and (dataStore.className == 'BlockedBinaryMatrix'):
      blockSizes = list(dataStore.blockSizes)
      if not blockSizes:
        blockSizes = dataSource.numDim * [1]
      i = dataDim.dim-1
      blockSizes[i] = blockSize
      
      try:
        dataStore.blockSizes = blockSizes
      except Implementation.ApiError, e:
        showError('Setting block size', e.error_msg, parent=self)

  def checkFileExists(self):

    fileName = self.specFileEntry.get()

    return os.path.exists(fileName) and os.path.isfile(fileName)

  def checkEndianess(self):

    fileName = self.specFileEntry.get()
    numberType = self.numTypePulldown.getObject()
    fileHeaderSize = self.headerSizeEntry.get()
    if (fileHeaderSize is None):
      fileHeaderSize = 0
    nByte = int(self.numBytesPulldown.getObject())

    isBigEndian = self.endianSelect.get()

    test = isDataBigEndian(fileName, numberType=numberType,
                           fileHeaderSize=fileHeaderSize,
                           nbytes=nByte)
                           
    if isBigEndian == test:
      return True
    else:
      return False

  def commitFileDetails(self):

    spectrum = self.fileSpectrum

    if not spectrum:
      return True

    isShape = isShapeSpectrum(spectrum)

    if not isShape:
      isBigEndian = self.endianSelect.get()

    if self.isModal:
      if self.checkFileExists():
        if not isShape and not self.checkEndianess():
          if isBigEndian:
            s = 'little'
            t = 'big'
          else:
            s = 'big'
            t = 'little'
          
          msg = 'The data looks %s endian, not %s endian, should it be set to %s endian?'  
          if showYesNo('Data endianess', msg % (s, t, s), parent=self):
            isBigEndian = not isBigEndian
            # no point doing below because quitting dialog in any case
            #if (isBigEndian):
            #  self.endianSelect.set(True)
            #else:
            #  self.endianSelect.set(False)
          else:
            msg = 'Data endianess looks wrong, cancel operation (otherwise continues)?'
            if showYesNo('Data endianess', msg, parent=self):
              return False
              
      else:
        fileName = self.specFileEntry.get()
        if fileName:
          msg  = 'The chosen data file does not exist (or exists and is not a file). '
          msg += 'Go back to the verification dialog (otherwise continues)?'
          if showYesNo('File does not exist', msg, parent=self):
            return False
          
    else:
      msg = 'Update spectrum "' + spectrum.name + '" file details with chosen values?'
      if not showYesNo('Update file details', msg, parent=self):
        return False

    fileName = self.specFileEntry.get()
    if fileName:
      if not isShape:
        nByte = int(self.numBytesPulldown.getObject())
        numberType = self.numTypePulldown.getObject()

        try:
          fileHeaderSize = self.headerSizeEntry.get()
          assert (fileHeaderSize is not None) and (fileHeaderSize >= 0)
        except:
          msg = 'File header size must be non-negative integer'
          showError('Illegal file header size', msg, parent=self)
          return False

      # NBNB TBD system for adding new files need changing
      setDataSourceFileName(spectrum, fileName)

      if not isShape:
        dataStore = spectrum.dataStore
        dataStore.isBigEndian = isBigEndian
        dataStore.nByte = nByte
        dataStore.numberType = numberType
        dataStore.headerSize = fileHeaderSize

    if not self.isModal:
      self.parent.updatedSpectrumFileDetails(spectrum)
      self.parent.drawWindows()

    return True
