
"""
======================COPYRIGHT/LICENSE START==========================

CalcRates.py: Part of the CcpNmr Analysis program

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

from math import log, sqrt

from ccpnmr.analysis.popups.EditFitGraph import EditFitGraphPopup
from ccpnmr.analysis.core.DataAnalysisBasic import DataFitting
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.DataAnalysisBasic import getNmrExpSeriesSampleConditions, getPeakSeriesIntensities
from ccpnmr.analysis.core.DataAnalysisBasic import matchSeriesPeaks, makeRatesList, getPeakSampledDimIntensities
from ccpnmr.analysis.core.DataAnalysisBasic import matchSampledExperimentPeaks, getFitMethodInfo, getFitErrorMethods
from ccpnmr.analysis.core.DataAnalysisBasic import getFitErrorMethod, setFitErrorMethod
from ccpnmr.analysis.core.ExperimentBasic import getSpectrumNoise, getSampledDimExperiments, getExperimentSampledDim

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.LabelDivider import LabelDivider
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showWarning
from memops.gui.ProgressBar import ProgressBar
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.core.Util import getAnalysisDataDim

from ccp.api.nmr import Nmr

GRADIENT_STRENGTH = 'gradient strength'
ALLOWED_CONDITION_TYPES = ('delay time', 'Time','num delays',
                           'mixing time', 'pulsing frequency', GRADIENT_STRENGTH)

class RowObject:
  def __init__(self, dataDimName, dataDim, tol):
    self.dataDimName = dataDimName
    self.dataDim = dataDim
    self.tol = tol

class CalcRatesPopup(BasePopup):
  """
  **Follow Peak Intensity Series for Rate Calculations etc**
  
  The purpose of this system is to expedite the extraction of peak intensity data
  from series of related experiments and subsequently fit a selected function to the
  data to extract parameters such as relaxation rates or NOE decay. The general idea
  is that the user sets up an "NMR series" that contains an array of experiments, or
  a single experiment stacked into separate planes (e.g. a pseudo 3D), where each
  experiment or plane in the series represents a different value for some time or
  frequency measurement being investigated. Examples of this  include T1, T2 & T1rho
  delay times for estimating relaxation rates and NOE mixing times. It should be noted
  that this system is specifically designed for the extraction of parameters derived
  from spectrum peak intensities where the position of related peaks across a series
  *do not move significantly*. For following peak movements and chemical shift
  changes the `Follow Shift Changes`_ tool may be used.

  The layout of the popup window is split into two tabs to reduce clutter. The
  first tab allows the user to setup and adjust all of the options used to
  perform the peak intensity extraction and function fitting. The second tab is
  used to actually perform the operations and display the results on a table of
  peak groups. Each group corresponds to a series of peaks, usually with a
  common resonance assignment, with one peak for each time value (i.e. per
  experiment or plane).

  The general idea is that the user selects a reference peak list, which will
  determine the location and assignment identities of the peaks being analysed.
  For proteins this reference is often a 15N HSQC peak list, in which case the
  analysis  operates on peak groups that correspond to amides of individual
  residues. The reference need not be part of the series of analysed
  experiments, but naturally peak locations should match. Also, the user chooses
  an NMR series that has been setup elsewhere with the relevant experiments or
  planes and their corresponding time parameters. The [Edit NMR Series] button
  will open the `NMR Series` popup window to create and adjust such series. The
  "Fitting Function" option is adjusted to say what kind of curve should be
  fitted to the peak intensity data (often this is a two-parameter exponential).
  The "Rate Type" and "Coherence Type" options are only used if results are
  stored in the CCPN project as a measurement list and do not affect the initial
  analysis. The Error Method determines how the errors in the parameters of the
  fitted function (e.g. error in relaxation rate) will be calculated.

  The "covariance" error method can be used if the measurement errors are
  normally distributed (which is often a reasonable assumption).  For each
  parameter the error (standard deviation) estimate is the square root of
  (the chi squared value times the covariance matrix diagonal term for that
  parameter).
  Reference: section 15.6, "Confidence Limits on Estimated Model Parameters"
  in Numerical Recipes, second edition.

  The "bootstrap" error method uses repeated sampling to provide an estimate
  of the error.  If there are N (x, y) points to be fit then each sampling
  takes N of those (x, y), but with replacement allowed, so some of the (x, y)
  might be repeated and some might be left out.  For each sampling the best
  fit is calculated and that determines the parameters for this specific
  sampling, which in turn allows an estimate of the error (standard deviation)
  over all samplings.  Analysis uses 1000 samples.
  Reference: "Bootstrap Methods for Standard Errors, Confidence Intervals and
  Other Measures of Statistical Accuracy", B. Efron and R. Tibshirani,
  Statistical Science, 1986, Vol. 1, No. 1, 54-77.

  The "jiggling" error method uses repeated sampling but here the (x, y) are
  both sampled from a normal distribution with mean the actual value and
  standard deviation the estimated data errors.  There is no real scientific
  basis for this estimate, so probably best avoided.

  The peak picking and search tolerance sections control how peaks are grouped
  together so that their intensities may be analysed. The basic process is that
  each reference peak position is used to locate a corresponding peak in each
  plane/experiment of the series. How exactly this is done depends on which
  options are checked. If a peak, at the right point of the NMR series, has the
  same assignment as the reference peak then that peak is used in preference to
  any others, irrespective of location. If a peak with a matching assignment
  cannot be found, the position of the reference peak is then used to locate the
  peaks for its group. When looking for peaks based on location the system
  checks to see if there are any existing, picked peaks in the series that are
  close to the reference (within the search tolerances). If no existing peaks
  are found for a point in the series then, should the "Pick new peaks" option
  be set, an attempt is made to pick a new peak extremum within the stated
  tolerances. If there is no extremum to pick, then the system may still add a
  "non-maxima" peak at exactly the reference position; useful where peak
  intensities dip into the noise levels, but are still helpful in a function
  fit. Also, having the "Assign groups?" option set means that after the first
  peak grouping, peaks will be linked via assignment and subsequent peak searches
  are not generally required.

  The peak grouping and function fitting is performed using the [Group & Fit
  Peaks] function. After the initial grouping the intensity curve fitting may be
  redone with via one of the "Re-fit" buttons; this useful if the fitting
  function is changed. When the curve fitting is done the parameter results from
  the fit, e.g. the "A" and "B" from an "A exp(-Bx)" equation, are immediately
  made available from the results table. Also, where relevant, any time constant
  values (one divided by the rate) are also presented. In the "Peak Groups &
  Analysis" table the user can see the fit results and analyse or adjust the
  peak groups. It is commonplace to look through all the intensity curves for
  each of the peak groups by using [Show Fit Graph]; here the user can check how
  well the curve-fit worked and whether any adjustments (e.g. in peak picking)
  need to be made or groups removed. See the `Fit Graph`_ documentation for
  details about how the resultant popup window operates. The "Y" value of the
  curves naturally come from the selected type of peak intensity and the "X"
  values come from those that were entered for the experimental points/planes in
  the NMR series. When the results have been checked, if the data is of a kind
  that corresponds to the "Rate Type" in the settings, then the user may save
  the values of the time constant, like T1 or T2, in a measurement list within
  the CCPN project using the "Make List" function. Alternatively the results may
  be used by directly exporting the fitted parameters from the table.

  The user may fit the peak intensity data outside of Analysis by exporting the
  values for the individual intensities (and any fitted parameters) using the
  [Export Data] button at the bottom right. This will produce an aligned,
  whitespace-separated plain text file that aims to be easy to analyse with
  external programs and scripts. If the fitting functions that are available in
  CCPN are not required, or cause problems, the fitting function may be set to
  "<None>", which means that the peaks are still grouped and that their
  intensities are available for export.

  **Caveats & Tips**
  
  Each peak group need not contain the same number of peaks if data is missing.
  
  A subset of peaks in a series may be analysed by reducing the number of peaks
  in the reference peak list. For example the user could make a copy of an HSQC
  peak list and then remove and peak locations that are not required in the
  analysis, e.g. for side chain NH2 peaks or severely overlapped peaks.

  Any peak picking done by the system uses the same spectrum peak finding
  parameters as is normally used in Analysis. Such parameters may be adjusted
  via the `Peak Finding`_ popup.

  If there are problems with grouping peaks together the user may assign all
  peaks that ought to go in the same group to the same resonances, thus
  connecting peaks together.   

  As with all analyses based upon peak intensity, the user should be cautious in
  regions of spectra where peaks are severely overlapped. In such cases using
  peak "height" rather than "volume" integral may help to a degree.
  
  The user should be cautious of manually copying or importing peaks into an
  intensity series from other spectra. Potentially both the peak position and
  intensity might not reflect the real data in the series. Accordingly, for
  copied peaks the user should recalculate intensities and check that the peak
  positions are at extrema, re-centring (<Ctrl> + <p>) as required.
  Alternatively, the automated series picking is normally adequate if the
  reference peak list matches reasonably and the picking PPM tolerances are
  appropriate.
  
  The NMR series that will be considered by this system are currently limited to
  the following types: "delay time", "time","num delays", "mixing time" and 
  "pulsing frequency".

  .. _`Follow Shift Changes`: FollowShiftChangesPopup.html
  .. _`NMR Series`: EditExperimentSeriesPopup.html
  .. _`Peak Finding`: EditPeakFindParamsPopup.html
  .. _`Fit Graph`: EditFitGraphPopup.html
  
  """

  def __init__(self, parent, *args, **kw):

    self.fitGraphPopup  = None
    self.waiting        = False
    self.expSeries      = None
    self.peakGroup      = None
    self.peakGroups     = []
    self.tolerances     = []
    self.peakList       = None
    self.fitFunc        = 3
    self.pickPeaks      = True
    self.pickNonMaxima  = True
    self.skipZeroMerit  = True
    self.doAssignGroups = True
    self.seriesType     = None
    self.gyroMagneticRatios = []
    self.guiParent = parent
    
    BasePopup.__init__(self, parent=parent, title='Data Analysis : Follow Intensity Changes', **kw)

  def body(self, guiFrame):
 
    self.geometry('700x600')
 
    self.rateType = self.getRateTypes()[0]
    self.intensityType = self.getIntensityTypes()[0]
    self.tolEntry = FloatEntry(self, returnCallback=self.setTolPpm, width=5)
    self.noise = 1.0
 
    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)
    
    tipTexts = ['Allows the user to setup how peak intensities are followed and which experiments to analyse',
                'The peaks that have been grouped by location, and the various equation parameters estimated for each']
    options = ['Settings','Peak Groups & Analysis']
      
    tabbedFrame = TabbedFrame(guiFrame, options=options, tipText=tipTexts)
    tabbedFrame.grid(row=0, column=0, sticky='nsew')
    self.tabbedFrame = tabbedFrame
    frameA, frameB = tabbedFrame.frames

    #
    # Settings
    #

    frameA.grid_columnconfigure(0, weight=1)

    row = 0
    div = LabelDivider(frameA, text='Experiment Series')
    div.grid(row=row, column=0, sticky='ew')

    row += 1
    frame = Frame(frameA)
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(3, weight=1)

    tipText = 'Selects the peak list used to confer assignments and groupings to analysed peaks '
    label = Label(frame, text='Reference Peak List', grid=(0,0))
    self.peakListPulldown = PulldownList(frame, self.changePeakList, tipText=tipText, grid=(0,1))
    
    tipText = 'The NMR series that carries details of which experiments to follow and their sampled conditions (e.g. delay times)'
    label = Label(frame, text='  NMR Experiment Series:', grid=(0,2))
    self.expSeriesPulldown = PulldownList(frame, self.changeExpSeries, tipText=tipText, grid=(0,3))
   
    Label(frame, text='', grid=(1,0))
    
    row += 1
    div = LabelDivider(frameA, text='Data Fitting')
    div.grid(row=row, column=0, sticky='ew')

    row += 1
    frame = Frame(frameA)
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(5, weight=1)

    tipText = 'Selects the type of graph to fit to the peak intensities within each grouping'
    label = Label(frame, text='Fitting Function:', grid=(0,0))
    self.fitFuncPulldown  =  PulldownList(frame, self.changeFitFunc, tipText=tipText, grid=(0,1))

    tipText = 'Whether to fit to peak heights or volume integrals'
    label = Label(frame, text='  Intensity Type:', grid=(0,2))
    self.intensityTypePulldown = PulldownList(frame, self.changeIntensityType, tipText=tipText, grid=(0,3))

    tipText = 'Selects which method is used to estimate errors in the graph fitting'
    label = Label(frame, text='  Error Method:', grid=(0,4))
    fitErrorMethods = getFitErrorMethods()
    fitErrorMethod = getFitErrorMethod(self.analysisProject)
    ind = fitErrorMethods.index(fitErrorMethod)
    self.fitErrorPulldown = PulldownList(frame, callback=self.setFitErrorMethod,
                                         texts=fitErrorMethods, index=ind,
                                         tipText=tipText, grid=(0,5))
 
    
    Label(frame, text='', grid=(1,0))
    row += 1
    div = LabelDivider(frameA, text='Relaxation Options')
    div.grid(row=row, column=0, sticky='ew')

    row += 1
    frame = Frame(frameA)
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(3, weight=1)

    tipText = 'The type of rate experiment performed; important for storage in CCPN project'
    label = Label(frame, text='Rate Type:', grid=(0,0))
    self.rateTypePulldown  =  PulldownList(frame, self.changeRateType, tipText=tipText, grid=(0,1))

    tipText = 'Sets which sub-type of rate experiment was performed'
    self.seriesTypeLabel = Label(frame, text='  Series Type:', grid=(0,2))
    self.seriesTypePulldown = PulldownList(frame, self.changeSeriesType, tipText=tipText, grid=(0,3))
 
    tipText = 'The total experimental delay time, for analyses that require such information'
    self.totalDelayLabel = Label(frame, text='  Total delay time (msec):')
    self.totalDelayEntry = FloatEntry(frame, text='', width=8, tipText=tipText)
    #self.totalDelayLabel.grid(row=3, column=2,sticky='e')
    #self.totalDelayEntry.grid(row=3, column=3, sticky='w')
    Label(frame, text='', grid=(1,0))
 
    row += 1
    div = LabelDivider(frameA, text='DOSY Parameters')
    div.grid(row=row, column=0, sticky='ew')
  
    row += 1
    frame = Frame(frameA)
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(3, weight=1)
  
    tipText = 'The total diffusion time in seconds between NMR gradients (from pulse sequence)'
    label = Label(frame, text='Total diffusion time (s):', grid=(0,0))
    self.diffTimeEntry = FloatEntry(frame, text=None, tipText=tipText, grid=(0,1), width=8)
  
    tipText = 'The duration of the gradient in the pulse sequence'
    label = Label(frame, text='  Gradient length (s):', grid=(0,2))
    self.gradLengthEntry = FloatEntry(frame, text=None, tipText=tipText, grid=(0,3), width=8)
  
    tipText = 'The maximum gradient strength (i.e. at 100%); used if relative gradient strengths have been entered'
    label = Label(frame, text='Full gradient strength (G/cm):', grid=(1,0))
    self.gradStrengthEntry = FloatEntry(frame, text=None, tipText=tipText, grid=(1,1), width=8)
  
    tipText = 'The short time interval between bipolar gradients, if bipolar gradients were used'
    label = Label(frame, text='  Bipolar separation (s):', grid=(1,2))
    self.bipolarSepEntry = FloatEntry(frame, text=None, tipText=tipText, grid=(1,3), width=8)
    
    Label(frame, text='', grid=(2,0))
 
    row +=1
    div = LabelDivider(frameA, text='Peak Picking')
    div.grid(row=row, column=0, sticky='ew')

    row +=1
    frame = Frame(frameA)
    frame.grid(row=row, column=0, sticky='nsew')
    frame.grid_columnconfigure(7, weight=1)
    frame.grid_rowconfigure(1, weight=1)
    
    tipText = 'Whether to pick new peaks, according to the reference peaks, should none of the right assignment be found'
    label = Label(frame, text='Pick new peaks?')
    label.grid(row=0, column=0,sticky='w')
    self.nonMaximaSelect = CheckButton(frame, callback=self.setPickPeaks, tipText=tipText)
    self.nonMaximaSelect.grid(row=0, column=1, sticky='w')
    self.nonMaximaSelect.set(True)

    tipText = 'Whether to pick new peaks at exactly the reference position, if no extremum exists'
    label = Label(frame, text='Pick non-maxima peaks?\n(At reference position)')
    label.grid(row=0, column=2,sticky='w')
    self.nonMaximaSelect = CheckButton(frame, callback=self.setNonMaxima, tipText=tipText)
    self.nonMaximaSelect.grid(row=0, column=3, sticky='w')
    self.nonMaximaSelect.set(True)

    tipText = 'Whether to assign newly picked peaks, using the reference assignment for the group'
    label = Label(frame, text='Assign groups?')
    label.grid(row=0, column=4,sticky='w')
    self.assignSelect = CheckButton(frame, callback=self.setAssign, tipText=tipText)
    self.assignSelect.grid(row=0, column=5, sticky='w')
    self.assignSelect.set(True)

    tipText = 'Whether to exclude peaks from analysis if they have figure-of-merit set to zero'
    label = Label(frame, text='Skip zero merit peaks?', grid=(0,6))
    self.assignSelect = CheckButton(frame, callback=self.setSkipMerit,
                                    grid=(0,7), selected=True, tipText=tipText)
    
    Label(frame, text='', grid=(1,0))
    
    row +=1
    div = LabelDivider(frameA, text='Peak Search Tolerances')
    div.grid(row=row, column=0, sticky='ew')

    row +=1
    frame = Frame(frameA)
    frame.grid(row=row, column=0, sticky='ew')
    Label(frame, text='Use noise threshold for finding maxima?', grid=(0,0))
    tipText = 'Whether to search for peaks above a special noise threshold, rather than the regular pick level.'
    self.noiseSelect = CheckButton(frame, tipText=tipText, grid=(0,1))
    self.noiseSelect.set(False)

    Label(frame, text='Threshold intensity: ', grid=(0,2))
    tipText = 'When searching for peaks, the minimum spectrum intensity level to consider. '
    tipText += 'Below this level peaks may still be picked at the exact reference position.'
    self.noiseEntry = FloatEntry(frame, text=None, tipText=tipText, grid=(0,3), width=16)

    row +=1
    frameA.grid_rowconfigure(row, weight=1)
    tipTexts = ['Spectrum dimension for tolerance setting',
                'Maximum ppm distance to group or pick peaks, relative to reference position']
    headingList = ['Dimension','Tolerance']
    editWidgets = [None, self.tolEntry]
    editGetCallbacks = [None, self.getTolPpm]
    editSetCallbacks = [None, self.setTolPpm]
    self.toleranceMatrix = ScrolledMatrix(frameA, headingList=headingList,
                                          editWidgets=editWidgets, tipTexts=tipTexts,
                                          editSetCallbacks=editSetCallbacks,
                                          editGetCallbacks=editGetCallbacks)
    self.toleranceMatrix.grid(row=row, column=0, sticky='nsew')
    
    #
    # Settings
    #

    frameB.grid_columnconfigure(0, weight=1)
    frameB.grid_rowconfigure(0, weight=1)

    
    headingList, tipTexts = self.getHeadingList(2)
    self.scrolledMatrix = ScrolledMatrix(frameB, headingList=headingList,
                                         callback=self.selectGroup,
                                         multiSelect=True, tipTexts=tipTexts,
                                         deleteFunc=self.removeGroup)
    self.scrolledMatrix.grid(row=0, column=0, sticky='nsew')

    tipTexts = ['Remove the selected peak groups from the table',
                'Redo the intensity graph fitting for the selected groups',
                'Show a graph of peak intensities and the fitted function',
                'Show a table of peaks within the selected groups',
                'Save the results as a data list in the CCPN project, using selected type options']
    texts = ['Remove\nSelected Groups', 'Re-fit\nSelected',
             'Show Fit\nGraph', 'Show\nPeaks',
             'Make Rates List']
    commands = [self.removeGroup, self.recalcRates,
                self.showFunctionFit, self.showPeaks,
                self.makeRatesList]
    bottomButtons = ButtonList(frameB, expands=True, tipTexts=tipTexts,
                                commands=commands, texts=texts)
    bottomButtons.grid(row=1, column=0, sticky='ew')

    self.removeButton  = bottomButtons.buttons[0]
    self.recalcButton  = bottomButtons.buttons[1]
    self.editFitButton = bottomButtons.buttons[2]
    self.peaksButton   = bottomButtons.buttons[3]
    self.makeMeasurementsButton = bottomButtons.buttons[4]
    
    #
    # Main
    #
    
    buttons = UtilityButtonList(tabbedFrame.sideFrame, helpUrl=self.help_url)
    buttons.grid(row=0, column=0, sticky='e')
    
    tipTexts = ['Search for peaks near to reference positions to create peak groups',
                'Perform graph fitting for the intensities of peaks in all groups',
                'Edit the selected NMR experiment series, e.g. to setup delay times',
                'Export the results as a flat text file']
    texts = ['Group & Fit Peaks','Re-fit All Groups',
             'Edit NMR Series','Export Data']
    commands = [self.matchPeaks, self.calcAllRates,
                self.editExpSeries, self.exportData]
    buttonList = ButtonList(guiFrame,  commands=commands,
                            texts=texts, tipTexts=tipTexts)
    buttonList.grid(row=1, column=0, sticky='ew')
    buttonList.buttons[0].config(bg='#B0FFB0')

    self.updateIntensityTypes()
    self.updateSeriesTypes()
    self.updateRateTypes()
    self.updateExpSeries()
    self.updatePeakLists()
    self.updateTolerances()
    self.updateFitFuncs()
    self.update()
    self.administerNotifiers(self.registerNotify)
  
  def administerNotifiers(self, notifyFunc):

    for func in ('__init__','delete','setName'):
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.Experiment', func)
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.DataSource', func)
    
    notifyFunc(self.updateNoiseLevel, 'ccp.nmr.Nmr.DataSource', 'setNoiseLevel')

    for func in ('__init__','delete'):
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.PeakList', func)
    
    notifyFunc(self.updatePeaks, 'ccp.nmr.Nmr.Peak', 'delete')

    for func in ('__init__','delete','setVariableUnit',
                 'setVariableDescriptor','setVariableName','setName',
                 'setExperiment','addExperiment','removeExperiment'):
      notifyFunc(self.updateExpSeries, 'ccp.nmr.Nmr.NmrExpSeries', func)

    for func in ('setConditionVaried',): 
      notifyFunc(self.updateExpSeries,'ccp.nmr.Nmr.SampledDataDim', func)


  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

  def exportData(self):
  
    from memops.gui.FileSelect import FileType
    from memops.gui.FileSelectPopup import FileSelectPopup
    from memops.gui.MessageReporter import showMulti

    texts = ['Fixed format', 'Tab separated', 'Comma separated', 'Cancel']
    objects = ['', '\t', ',', None]
    sep = showMulti('Export format', 'Select export format:', texts, objects, self)
    if sep is None:
      return

    lines = []
    if self.peakGroups and self.expSeries:
      expSeries = self.expSeries
      
      sampledDim = None
      if expSeries.className == 'Experiment':
        dataDim    = getExperimentSampledDim(expSeries)
        sampledDim = dataDim.dim
        values = dataDim.pointValues
     
      else:
        values = []
        conditionDict = getNmrExpSeriesSampleConditions(expSeries)
        for conditions in conditionDict.values():
          for condition in conditions:
            values.append([condition.value, condition.sampleConditionSet.experiments])
          
        values.sort()
        
        experiments = [x[1] for x in values]
        values = [x[0] for x in values]
    
      if sep == '':
        line = '!' + ' ' * 20
        for value in values:
          line += '%13.13s' % ('%.3f' % value)
        line += '%15.15s' % (' Fit Error')
        line += ' Fit Parameters'
      else:
        for peakGroup in self.peakGroups:
          refPk  = peakGroup.refObject
          ndim = len(refPk.peakDims)
          break
        else:
          ndim = 1
        fields = ['Dim %d' % (n+1) for n in range(ndim)]
        fields.extend(['%.3f' % value for value in values])
        fields.append('Fit Error')
        fields.append('Fit Parameters')
        line = sep.join(fields)
      lines.append(line)
    
      nVal = len(values)
    
      for peakGroup in self.peakGroups:
        refPk  = peakGroup.refObject
        peaks  = peakGroup.objects
        if peakGroup.fitError is None:
          fitErr = ''
        else:
          fitErr = '%.5f' % peakGroup.fitError
        params = ['%.5f' % x for x in peakGroup.parameters]
      
        iValues = ['-'] * nVal
        maxValue = 0
        texts = [pd.annotation or '-' for pd in refPk.sortedPeakDims() if pd.dataDimRef]
        if sep == '':
          text = ' '.join(texts)
          line = '%22.22s' % text
        else:
          fields = texts
        
        if sampledDim is None:
          for peak in peaks:
            intensity = peak.findFirstPeakIntensity(intensityType=self.intensityType)
            
            if intensity:
              i = 0
              for experiments0 in experiments:
                if peak.peakList.dataSource.experiment in experiments0:
                  iValues[i] = intensity.value
                  maxValue = max(maxValue, abs(intensity.value))
                  break
                  
                i += 1
        
        else:
          for peak in peaks:
            peakDim   = peak.findFirstPeakDim(dim=sampledDim)
            index     = int(peakDim.position)-1
            intensity = peak.findFirstPeakIntensity(intensityType=self.intensityType)
 
            if intensity:
              iValues[index] = intensity.value
              maxValue = max(maxValue, abs(intensity.value))
                
        for (i, value) in enumerate(iValues):
          if value != '-':
            if maxValue < 10000:
              iValues[i] = '%.5f' % value
            else:
              iValues[i] = '%.2f' % value

        if sep == '':
          line += ' '.join(['%12.12s' % x for x in iValues])
          line += ' %14.14s ' % fitErr
          line += ' '.join(['%14.14s' % x for x in params])
        else:
          fields.extend(iValues)
          fields.append(fitErr)
          fields.extend(params)
          line = sep.join(fields)
        lines.append(line)
        
      
    if lines:
      fileTypes = [  FileType('All', ['*'])]
      fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
                 title = 'Export Rates Data', dismiss_text = 'Cancel',
                 selected_file_must_exist = False)

      fileName = fileSelectPopup.getFile()
      if fileName:
        file = open(fileName, 'w')
        
        for line in lines:
          file.write(line)
          file.write('\n')
        
        file.close()      

  def getTolerances(self):
  
    if not self.tolerances:
      tolerances = []
 
      if self.peakList:
        for dataDim in self.peakList.dataSource.sortedDataDims():
          if dataDim.className == 'SampledDataDim':
            tolerances.append( ['Sampled',dataDim,None])
            continue
        
          expDimRefs = dataDim.expDim.expDimRefs
          if expDimRefs:
            name = self.getDataDimName(dataDim)
            isotopes = []
            for expDimRef in expDimRefs:
              isotopes += expDimRef.isotopeCodes
              
            tol = getAnalysisDataDim(dataDim).assignTolerance
            if '1H' in isotopes:
              tolerances.append( [ name,dataDim,tol ] )
            elif '13C' in isotopes:
              tolerances.append( [ name,dataDim,tol ] )
            elif '15N' in isotopes:
              tolerances.append( [ name,dataDim,tol ] )
            else:
              tolerances.append( [ name,dataDim,tol ] )
      self.tolerances = [RowObject(name,dataDim,tol) for (name,dataDim,tol) in tolerances]

    return self.tolerances

  def updateTolerances(self):

    textMatrix = []
    objectList = []
    for obj in self.getTolerances():
      if obj.tol is None:
        continue
        
      datum = []
      datum.append(obj.dataDimName)
      datum.append(obj.tol)
      textMatrix.append(datum)
      objectList.append(obj)
    
    self.toleranceMatrix.update(objectList=objectList,textMatrix=textMatrix)
 
  def setTolPpm(self, event):
    
    obj = self.toleranceMatrix.currentObject
    if obj:
      # Actually setting the self.tolerances list
      obj.tol = self.tolEntry.get() 

    self.updateTolerances()
  
  def getTolPpm(self, obj):
    
    if obj:
      self.tolEntry.set(obj.tol)

  def setPickPeaks(self, boolean):
  
    self.pickPeaks = boolean

  def setNonMaxima(self, boolean):
  
    self.pickNonMaxima = boolean

  def setAssign(self, boolean):
  
    self.doAssignGroups = boolean

  def setSkipMerit(self, boolean):
  
    self.skipZeroMerit = boolean

  def setFitErrorMethod(self, value):

    # TBD: registerNotify to notice when this set externally
    setFitErrorMethod(self.analysisProject, value)

  def removeGroup(self, *event):
  
    peakGroups = self.scrolledMatrix.currentObjects
    
    if peakGroups:
      for peakGroup in peakGroups:
        self.peakGroups.remove(peakGroup)
        peakGroup.delete()
      
      self.peakGroup = None
      self.updateAfter()
    
  def showFunctionFit(self, noOpen=0):
    
    dataFitting = self.peakGroup
  
    if dataFitting:
      peaks    = dataFitting.objects
      refPeak  = dataFitting.refObject
      function = dataFitting.fitFunction
      title    = '.'.join([pd.annotation or '?' for pd in refPeak.sortedPeakDims() if pd.dataDimRef])
      
    else:
      title    = ''
      
    if self.fitGraphPopup:
      self.fitGraphPopup.update(dataFitting=dataFitting,
                                xLabel=self.rateType,
                                yLabel=self.intensityType,
                                graphTitle=title )
      if not noOpen:
        self.fitGraphPopup.open()
    else:
      self.fitGraphPopup = EditFitGraphPopup(self, dataFitting,
                                             self.getPeakSeriesData, self.updateGroup,
                                             self.rateType, self.intensityType,
                                             nextSetFunction=self.showNextGroup,
                                             prevSetFunction=self.showPrevGroup,
                                             graphTitle=title)
      self.guiParent.popups['fit_graph'] = self.fitGraphPopup

  def showNextGroup(self):
  
    peakGroups = self.scrolledMatrix.objectList
    
    if peakGroups:
      if self.peakGroup and (self.peakGroup in peakGroups):
        i = peakGroups.index(self.peakGroup)
        i += 1
        if i >= len(peakGroups):
          i = 0
        self.peakGroup = peakGroups[i]
        self.scrolledMatrix.selectNthObject(i)
        self.showFunctionFit()

  def showPrevGroup(self):
  
    peakGroups = self.scrolledMatrix.objectList
    
    if peakGroups:
      if self.peakGroup and (self.peakGroup in peakGroups):
        i = peakGroups.index(self.peakGroup)
        i -= 1
        if i < 0:
          i = len(peakGroups)-1
        self.peakGroup = peakGroups[i]
        self.scrolledMatrix.selectNthObject(i)
        self.showFunctionFit()

  def updateGroup(self, dataFitting):

    if self.peakGroup:
      if dataFitting.parameters:
        self.addDerivedParams(dataFitting) 
                  
      self.updateAfter()

  def getPeakSeriesData(self, dataFitting):
  
    peaks   = dataFitting.objects
    xData   = []
    yData   = []
    xWidths = []
    yWidths = []
    errData = None
    if self.expSeries:
      expSeries = self.expSeries
      
      useTotalDelay = False
      if expSeries.className == 'Experiment':
        (xData,yData,xWidths,yWidths, refIntensity) = getPeakSampledDimIntensities(peaks, expSeries, self.intensityType)
      
        dataDim   = getExperimentSampledDim(expSeries)
        delayTime = (self.totalDelayEntry.get() or 0.0)/ 1000.0
        
        if (refIntensity is not None) and (refIntensity > 0) and delayTime:
          if dataDim:
            if dataDim.conditionVaried == 'num delays':
              xData = [delayTime/(8.0*x) for x in xData]

            if dataDim.conditionVaried in ('num delays', 'pulsing frequency'):
              yData = [log(refIntensity/abs(y))/delayTime for y in yData]
              spectrum = dataDim.dataSource
              if spectrum.noiseLevel:
                scale = spectrum.noiseLevel / (delayTime * refIntensity)
                errData = [scale*sqrt(1+(y/refIntensity)**2) for y in yData]

      else:
        # wb104: 24 Mar 2010: sampleConditionType added to get around
        # default 'delay time' limitation in getPeakSeriesIntensities
        sampleConditionType = self.getExpSeriesCondition(expSeries)
        (xData,yData,xWidths,yWidths) = getPeakSeriesIntensities(peaks, expSeries, self.intensityType, sampleConditionType)
    
    seriesType, unit = self.getExpSeriesType(self.expSeries)
    if seriesType == GRADIENT_STRENGTH and unit == '%':
      # wb104: 10 Mar 2014: convert to G/cm from %
      gradStren = 0.01*(self.gradStrengthEntry.get() or 100.0)
      xData = [gradStren*x for x in xData]
    
    dataFitting.dataX   = xData  
    dataFitting.dataY   = yData  
    dataFitting.dataErr = errData  
    dataFitting.errorsX = xWidths
    dataFitting.errorsY = yWidths
    
    return dataFitting

  def calcAllRates(self):

    if not (self.peakGroups and self.expSeries):
      return
    
    if not self.fitFunc:
      showWarning('Warning','No fitting function selected', parent=self)
      for dataFitting in self.peakGroups:
        dataFitting.resetFit()
        self.updateAfter()
      return
    
    self.tabbedFrame.select(1)
    
    progressBar = ProgressBar(self, text="Calculating rates ",total = len(self.peakGroups))
    failedPeaks = []

    for dataFitting in self.peakGroups:
    
      progressBar.increment()
      peaks = dataFitting.objects
      peaks = [p for p in peaks if p.figOfMerit]
      if len(peaks) < 2:
        continue
      
      dataFitting.objects = peaks
      dataFitting = self.getPeakSeriesData(dataFitting)
      dataFitting.noiseLevel   = self.noise
      dataFitting.fitFunction  = self.fitFunc
      dataFitting.fitErrorFunction = getFitErrorMethod(self.analysisProject)
      try:
        dataFitting.fit()
        ok = True
      except:
        failedPeaks.append(dataFitting.refObject)
        ok = False

      if ok:
        self.addDerivedParams(dataFitting)
    
    progressBar.destroy()
    
    if failedPeaks:
      self.warnFailedPeakGroups(failedPeaks)
      
    self.updateAfter()

  def warnFailedPeakGroups(self, failedPeaks):

    numFail   = len(failedPeaks)
    positions = []
    
    for peakF in failedPeaks:
      values = ['%.3f' % pd.value for pd in peakF.sortedPeakDims() if pd.dataDimRef]
      positions.append('(' + ','.join(values) + ')')
    
    failText = ', '.join(positions)
    showWarning('Fit failure','Fit failed for %d groups at positions %s' % (numFail,failText), parent=self)

  def addDerivedParams(self, dataFitting):
    
    seriesType, unit = self.getExpSeriesType(self.expSeries)
    
    if seriesType == GRADIENT_STRENGTH:
      if (dataFitting.fitFunction == 10) and self.gyroMagneticRatios:
        gyro = self.gyroMagneticRatios[0]
        q = dataFitting.parameters[1]
        qError = dataFitting.parameterErrors[1]
 
        diffTime = self.diffTimeEntry.get()
        gradLen = self.gradLengthEntry.get()
        gradStren = self.gradStrengthEntry.get() or 100.0
        bipolSep = self.bipolarSepEntry.get() or 0.0
 
        if diffTime and gradLen:
          # wb104: 10 Mar 2014: below if doesn't belong because have converted xData values before the fit
          #if unit == '%':
          #  q *= gradStren * 1e2 # Gauss per cm to Tesla per metre / 100
          #  qError *= gradStren * 1e2
          #else:
          q *= 1e4 # Gauss per cm to Tesla per metre
          qError *= 1e4
        
          f = diffTime - (gradLen/3.0) - (bipolSep/2.0)

          scale = gyro * gyro * gradLen * gradLen * f

          dc = q/scale
          dcError = qError/scale

          dataFitting.derivedData = dc
          dataFitting.derivedDataErrors = dcError
        
        else:

          dataFitting.derivedData = None
          dataFitting.derivedDataErrors = None
    
    else:
      if dataFitting.fitFunction in (2,3,4,8):
        rate = dataFitting.parameters[1]
        rateError = dataFitting.parameterErrors[1]
        c = None
        if rate:
          c    = 1.0/float(rate)
          cMin = 1.0/float(rate+rateError)
          if rateError != rate:
            cMax = 1.0/float(rate-rateError)
          else:
            cMax = c + c - cMin
 
        cError = None
        if c is not None:
          cError = (cMax-cMin)/2.0
 
        dataFitting.derivedData = c
        dataFitting.derivedDataErrors = cError
 
  def recalcRates(self):

    peakGroups = self.scrolledMatrix.currentObjects
    if peakGroups and self.expSeries:
      if not self.fitFunc:
        showWarning('Warning','No fitting function selected', parent=self)
        for dataFitting in peakGroups:
          dataFitting.resetFit()
          self.updateAfter()
        return
      
      progressBar = ProgressBar(self, text="Recalculating rates ",total=len(self.peakGroups))
      
      failedPeaks = []  
      for dataFitting in peakGroups:
        progressBar.increment()
        peaks = dataFitting.objects
        peaks = [p for p in peaks if p.figOfMerit]
        if len(peaks) < 2:
          continue
        
        dataFitting.objects = peaks
        dataFitting = self.getPeakSeriesData(dataFitting)
        dataFitting.noiseLevel   = self.noise
        dataFitting.fitFunction  = self.fitFunc
        dataFitting.fitErrorMethod = getFitErrorMethod(self.analysisProject)
        
        try:
          dataFitting.fit()
          ok = True
        except:
          failedPeaks.append(dataFitting.refObject)
          ok = False

        if ok:
          self.addDerivedParams(dataFitting)         
      
      progressBar.destroy()
      
      if failedPeaks:
        self.warnFailedPeakGroups(failedPeaks)  
              
      self.updateAfter()

  def matchPeaks(self):
  
    if not (self.peakList and self.expSeries):
      return 
    
    self.tabbedFrame.select(1)
    
    tolerances  = [ obj.tol for obj in self.tolerances ]
    
    progressBar = ProgressBar(self, text="Grouping peaks ")
    
    if self.noiseSelect.get():
      noiseThreshold = self.noiseEntry.get()
    else:
      noiseThreshold = None
    
    expSeries = self.expSeries
    if expSeries.className == 'Experiment':
      groups = matchSampledExperimentPeaks(self.peakList ,expSeries, tolerances,
                                           self.pickPeaks, self.pickNonMaxima,
                                           self.doAssignGroups, progressBar=progressBar,
                                           noiseThreshold=noiseThreshold)

    else:
      groups = matchSeriesPeaks(self.peakList, expSeries, tolerances,
                                self.pickPeaks, self.pickNonMaxima,
                                self.doAssignGroups, progressBar=progressBar,
                                noiseThreshold=noiseThreshold)
 
    progressBar.destroy()

    self.clearPeakGroups()
    
    failedPeaks = []
    
    i = 0
    progressBar = ProgressBar(self, text="Fitting data ",total = len(groups))
    for peaks in groups:
      progressBar.increment()
      
      if self.skipZeroMerit:
        peaks2 = [p for p in peaks[1:] if p.figOfMerit]
      else:
        peaks2 = [p for p in peaks[1:]]
      
      if len(peaks2) < 2:
        continue
    
      dataFitting = DataFitting(refObject=peaks[0],
                                objects=peaks2,
                                fitFunction=self.fitFunc,
                                fitErrorFunction=getFitErrorMethod(self.analysisProject),
                                noiseLevel=self.noise,
                                guiParent=self)
      dataFitting = self.getPeakSeriesData(dataFitting)
      
      try:
        dataFitting.fit()
        ok = True
      except:
        failedPeaks.append(peaks[0])
        ok = False

      if ok or self.fitFunc == 0:
        self.addDerivedParams(dataFitting)     
        self.peakGroups.append(dataFitting)
      
      i += 1
    
    progressBar.destroy()
    
    if failedPeaks:
      self.warnFailedPeakGroups(failedPeaks)
    
    self.updateAfter()

  def getIntensityTypes(self):
  
    intensityTypes = ['height','volume']
    
    return intensityTypes

  def changeIntensityType(self, intensityType):
  
    self.intensityType = intensityType

  def updateIntensityTypes(self, *opt):
  
    names = self.getIntensityTypes()
    if self.intensityType not in names:
      self.intensityType = names[0]
    
    index = names.index(self.intensityType)
    
    self.intensityTypePulldown.setup(names, names, index)

  def getSeriesTypes(self):

    def getEnumeration(className):
      return list(self.project._metaclass.metaObjFromQualName('ccp.nmr.Nmr.'+className).enumeration)

    if self.rateType == 'T1':
      return getEnumeration('T1CoherenceType')
      
    elif self.rateType == 'T1rho':
      return getEnumeration('T2CoherenceType')
      
    elif self.rateType == 'T2':
      return getEnumeration('T2CoherenceType')
      
    elif self.rateType == 'H/D exchange':
      return []
      
    elif self.rateType == 'H/D protection':
      return ['binding', 'folding', 'intrinsic']
      
  def changeSeriesType(self, seriesType):

    self.seriesType = seriesType


  def updateSeriesTypes(self, *opt):
   
    if self.rateType[0] == 'T':
      self.seriesTypeLabel.set('  Coherence Type')

    elif self.rateType[0] == 'H':
      self.seriesTypeLabel.set('  Protection Type')

    else:
      self.seriesTypeLabel.set('  Series Type')
    
    index = 0  
    names = self.getSeriesTypes()
    if names:
      if (self.seriesType) or (self.seriesType not in names):
        self.seriesType = names[0]
    
      index = names.index(self.seriesType)
    
    self.seriesTypePulldown.setup(names, names, index)


  def getFitFuncs(self):
  
    fitFunctions = [x[0] for x in getFitMethodInfo()]
    fitFunctions = ['<None>'] + fitFunctions

    return fitFunctions

  def changeFitFunc(self, funcIndex):
  
    self.fitFunc = funcIndex
  
  def updateFitFuncs(self, *opt):
  
    names = self.getFitFuncs()
    if self.fitFunc > len(names):
      self.fitFunc = 1
      
    indices = range(len(names))
    
    self.fitFuncPulldown.setup(names, indices, self.fitFunc)

  def getRateTypes(self):
    
    rateTypes = ['T1','T2','T1rho']
    
    return rateTypes
 
  def changeRateType(self, rateType):

    self.rateType = rateType
    self.updateSeriesTypes()
    self.updateButtons()
     
  def updateRateTypes(self, *rateTypes):
    
    names = self.getRateTypes()
    if names:
      if (not self.rateType) or ( self.rateType not in names ):
        self.rateType = names[0]
    
    
      self.rateTypePulldown.setup(names, names, names.index(self.rateType))

  def getVariableDescriptor(self, rateType):
    
    descriptor = '%s interval' % (rateType)
    
    return descriptor
    
  def makeRatesList(self):
  
    if self.peakGroups and self.expSeries:
    
      peakGroups = []
      values    = []
      errors    = []
      for dataFitting in self.peakGroups:
        peakGroups.append( dataFitting.objects )
        values.append( dataFitting.derivedData )
        errors.append( dataFitting.derivedDataErrors )
     
      measurementList = makeRatesList(self.rateType, self.expSeries,
                                      self.seriesType, peakGroups,
                                      values, errors=errors)
  
      if measurementList:
        self.parent.editMeasurements(measurementList=measurementList)
      
  def showPeaks(self):
    
    peakGroups = self.scrolledMatrix.currentObjects
    
    if peakGroups:
      peaks = []
      for peakGroup in peakGroups:
        peaks.extend(peakGroup.objects)
      
      self.guiParent.viewPeaks(peaks)
    
  def selectGroup(self, peakGroup, row, col):
  
    if peakGroup:
      self.peakGroup = peakGroup
      self.updateButtons()
 
  def editExpSeries(self):
  
    self.guiParent.editExpSeries()
    
    
  def getExpSeries(self):
  
    expSeries = getSampledDimExperiments(self.nmrProject) + self.nmrProject.sortedNmrExpSeries()

    """wb104: 08 Mar 2010: changed to allow everything, Tim not
       sure this is good idea, so leave below code here for now
    expSeries = []
    for experiment in getSampledDimExperiments(self.nmrProject):
      dataDim = getExperimentSampledDim(experiment)
      if dataDim.conditionVaried in ALLOWED_CONDITION_TYPES:
        expSeries.append(experiment)
    
    for expSeries0 in self.nmrProject.nmrExpSeries:
      for conditionName in ALLOWED_CONDITION_TYPES:
        if conditionName in expSeries0.conditionNames:
          expSeries.append(expSeries0)
          break
"""
     
    return expSeries

  
  def getExpSeriesName(self, series):
  
    if series.className == 'Experiment':
      dataDim = getExperimentSampledDim(series)
      return '%s:%s' % (series.name,dataDim.conditionVaried)
      
    else:
      conditionType = ','.join([n for n in series.conditionNames])
      return '%s:%s:%s' % (series.serial,series.name,conditionType)
   
  def getExpSeriesType(self, series):
  
    if series.className == 'Experiment':
      dataDim = getExperimentSampledDim(series)
      
      return dataDim.conditionVaried, dataDim.unit
      
    else:
      conditionTypes = series.conditionNames
      conditionType = ','.join([n for n in conditionTypes])
      
      unit = None
      
      experiment = series.findFirstExperiment()
      if experiment:
        scSet = experiment.sampleConditionSet
        
        if scSet:
          for condition in conditionTypes:
             sampleCondition = scSet.findFirstSampleCondition(condition=condition)
 
             if sampleCondition:
               unit = sampleCondition.unit
               break
         
      return conditionType, unit
 
  def getExpSeriesCondition(self, series):
  
    if series.className == 'Experiment':
      dataDim = getExperimentSampledDim(series)
      return dataDim.conditionVaried
      
    else:
      if series.conditionNames:
        conditionType = series.conditionNames[0]
      else:
        conditionType = ''
      return conditionType
  
  def getExpSeriesDataSource(self, expSeries):
      
    dataSource = None
    
    if expSeries:
      if expSeries.className == 'Experiment':
        dataDim = getExperimentSampledDim(expSeries)
        dataSource = dataDim.dataSource
 
      else:
        experiment = expSeries.findFirstExperiment()
 
        if experiment:
          dataSource = experiment.findFirstDataSource()
    
    return dataSource
  
  def changeExpSeries(self, expSeries):
  
    if self.expSeries is not expSeries:
      self.expSeries  = expSeries
      self.clearPeakGroups()
      self.updatePeakLists()
      
      dataSource = self.getExpSeriesDataSource(expSeries)
      if dataSource:  
        self.noiseEntry.set(dataSource.noiseLevel or None)

      self.updateAfter()
      
  def updateNoiseLevel(self, *dataSource):
    
    if dataSource is self.getExpSeriesDataSource(self.expSeries):
      self.updateAfter()
     
  def updateExpSeries(self, *expSeries0):
    
    expSeries = self.expSeries
    expSeriesList = self.getExpSeries()
    names = [ self.getExpSeriesName(es) for es in expSeriesList ]
    index = 0

    if expSeriesList:
      if expSeries not in expSeriesList:
        expSeries = expSeriesList[0]
    
      index = expSeriesList.index(expSeries)
    
    else:
      expSeries = None
        
    if expSeries is not self.expSeries:
      self.changeExpSeries(expSeries)
        
    self.expSeriesPulldown.setup(names, expSeriesList, index)

  def clearPeakGroups(self):
  
    for peakGroup in self.peakGroups:
      peakGroup.delete()
      
      self.peakGroups = []

  def getDataDimName(self,dataDim):
  
    isotopes = set()
    if dataDim.expDim.expDimRefs:
      for expDimRef in dataDim.expDim.expDimRefs:
        for isotope in expDimRef.isotopeCodes:
          isotopes.add(isotope)
  
    isotopes = list(isotopes)
    isotopes.sort()
    isotopeStr = ','.join(isotopes)
  
    name = 'F%d %s' % (dataDim.dim, isotopeStr)
    return name
    
  def getPeakListName(self, peakList):
  
    spectrum = peakList.dataSource
    return '%s:%s:%d' % (spectrum.experiment.name,spectrum.name,peakList.serial)
    
  def getMinNumDim(self):

    minNumDim = 2
    expSeries = self.expSeries
    if expSeries:
      if isinstance(expSeries, Nmr.Experiment):
        for dataSource in expSeries.dataSources:
          # TBD: below assumes exactly one SampledDataDim
          minNumDim = min(minNumDim, dataSource.numDim-1)
      else:
        for experiment in expSeries.experiments:
          for dataSource in experiment.dataSources:
            minNumDim = min(minNumDim, dataSource.numDim)

    return minNumDim

  def getPeakLists(self):
  
    minNumDim = self.getMinNumDim()
    names = []
    for experiment in self.nmrProject.sortedExperiments():
      for spectrum in experiment.sortedDataSources():
        if (spectrum.dataType == 'processed') and (spectrum.numDim >= minNumDim):
          for peakList in spectrum.sortedPeakLists():
            names.append( [self.getPeakListName(peakList),peakList] )
          
    return names
  
  def changePeakList(self, peakList):
  
    if peakList is not self.peakList:
      self.tolerances = []
      self.peakList = peakList
      self.gyroMagneticRatios = []
      
      for expDim in peakList.dataSource.experiment.sortedExpDims():
        expDimRefs = expDim.sortedExpDimRefs()
        
        if expDimRefs:
          isotope = expDimRefs[0].isotopes[0]
          self.gyroMagneticRatios.append(isotope.gyroMagneticRatio)
      
    self.updateTolerances()
    self.noise = getSpectrumNoise(self.peakList.dataSource)
    if self.peakGroups:
      pass
      # do anything here?
  
  def updatePeaks(self, peak):
  
    if peak.isDeleted:
    
      i = 0
      for group in self.peakGroups:
        peaks = group.objects
        
        if peak in peaks:
          peaks.remove(peak)
      
          if len(peaks) < 2:
            self.peakGroups.pop(i)
            group.delete()

          self.updateAfter()
          break
        
        i += 1

  
  def updatePeakLists(self, *peakList):
  
    data = self.getPeakLists()
    if not data:
      return
    
    peakList = self.peakList  
    names = [x[0] for x in data]
    peakLists = [x[1] for x in data]
    
    if peakLists:
      if peakList not in peakLists:
        peakList = peakLists[0]
    
      index = peakLists.index(peakList)
    
    else:
      peakList = None
      
    if self.peakList is not peakList:  
      self.changePeakList(peakList)
      
    self.peakListPulldown.setup(names, peakLists, index)
      
      
  def getHeadingList(self, numDim):
  
    if self.expSeries:
      seriesType, seriesUnit = self.getExpSeriesType(self.expSeries)
    else:
      seriesType = None  
  
    headingList = ['#',]
    tipTexts = ['Peak group row number',]
    
    for i in range(numDim):
      headingList.append( 'Assign\nF%d' % (i+1) )
      tipTexts.append('Assignment in F%d dimension, from reference peak list, for group' % (i+1) )
    
    if seriesType == GRADIENT_STRENGTH:
       headingList.extend(['Diffusion\nCoefficient','DC\nError'])
       tipTexts.extend(['Diffucsion coeffient in m^2/s; calculated using the DOSY parameters and fitted decay',
                        'Estimated error in the diffucsion coeffient'])
    
    else:
       headingList.extend(['Time\nConstant','TC\nError'])
       tipTexts.extend(['Time constant, e.g. T1 or T2; one over the decay rate for exponential fits',
                        'Estimated error in the time constant'])
    
    
    headingList.extend(['Fit\nError','Num\nPeaks','Function'])
    tipTexts.extend(['Error in the fit of the peak intensities to the fitting function',
                     'Number of peaks in groups',
                     'Function used to fit to peak intensities'])
        
    ordA = ord('A')
    numParams = getFitMethodInfo()[self.fitFunc-1][1]
    for i in range(numParams):
      headingList.append('Fit\nParam %s' % (chr(ordA+i), ))
      tipTexts.append('Parameter "%s" from the graph fitting equation' % (chr(ordA+i), ) ) 
      
    for i in range(numParams):
      headingList.append('Param\nError %s' % (chr(ordA+i), ))
      tipTexts.append('Estimated error in parameter "%s", according to selected error method' % (chr(ordA+i), ) ) 

    return headingList, tipTexts
  
  
  def updateButtons(self):
  
    if self.peakGroup:
      self.removeButton.enable()
      self.recalcButton.enable()
      self.editFitButton.enable()
      self.peaksButton.enable()
    else:
      self.removeButton.disable()
      self.recalcButton.disable()
      self.editFitButton.disable()
      self.peaksButton.disable()
 
    if self.peakGroups:
      self.makeMeasurementsButton.enable()
    else:
      self.makeMeasurementsButton.disable()
      
    if self.rateType:
      self.makeMeasurementsButton.config(text='Make\n%s List' % self.rateType)
    else:
      self.makeMeasurementsButton.config(text='Make\nRates List')
 
  def updateAfter(self, *object):
    
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
    
  def update(self):
   
    self.updateButtons()
    
    # wb104: 27 Mar 2014: below causes exception when switch between expSeries 
    # (that exception only happens if the functionFit is open)
    # this means that some graphs might need to be updated manually
    #if self.fitGraphPopup:
    #  self.showFunctionFit(noOpen=True)

    nDim = 2
    if self.peakList:
      dataDims = self.peakList.dataSource.findAllDataDims(className='FreqDataDim')
      nDim = len(dataDims)
  
    refPlanes = 0
    useTotalDelay = False
    if self.expSeries and (self.expSeries.className == 'Experiment'):
      dataDim = getExperimentSampledDim(self.expSeries)
      refPlane = getAnalysisDataDim(dataDim).refSamplePlane
      if refPlane is not None:
        refPlanes = 1
      
      if dataDim.conditionVaried in ('num delays','pulsing frequency'):
        useTotalDelay = True
      
    if useTotalDelay:  
      self.totalDelayLabel.grid(row=0, column=4, sticky='w')
      self.totalDelayEntry.grid(row=0, column=5, sticky='w')
    else:   
      self.totalDelayLabel.grid_forget()
      self.totalDelayEntry.grid_forget()

    dataSource = self.getExpSeriesDataSource(self.expSeries)
    if dataSource:  
      noiseLevel = self.noiseEntry.get()
      if noiseLevel is not None:
        dataSource.noiseLevel = noiseLevel
    
    nParams     = getFitMethodInfo()[self.fitFunc-1][1]
    headingList, tipTexts = self.getHeadingList(nDim)
    objectList  = self.peakGroups
    textMatrix  = []
    
    functionNames = [None, ] + [x[0] for x in getFitMethodInfo()]
    
    i = 0
    for dataFitting in objectList:
      i += 1
      datum = []
      datum.append(i)
      
      peak = dataFitting.refObject
      peakDims = peak.sortedPeakDims()
      
      for j in range(nDim):
        datum.append(peakDims[j].annotation)
      
      datum.append(dataFitting.derivedData)
      datum.append(dataFitting.derivedDataErrors)
      datum.append(dataFitting.fitError)
      datum.append(len(dataFitting.objects)-refPlanes)
      datum.append(functionNames[dataFitting.fitFunction])

      params = dataFitting.parameters
      errors = dataFitting.parameterErrors
       
      if dataFitting.fitFunction == 10:
        # wb104: 10 Mar 2014: dangerous so take copy to play safe
        params = params[:]
        params[1] = sqrt(params[1])
        errors = errors[:]
        errors[1] = sqrt(errors[1]) # Only temporary!
      
      for j in range(nParams):
        if j < len(params):
          datum.append(params[j])
        else:
          datum.append(None)

      for j in range(nParams):
        if j < len(errors):
          datum.append(errors[j])
        else:
          datum.append(None)
     
      textMatrix.append( datum )
   
    self.scrolledMatrix.update(textMatrix=textMatrix, objectList=objectList,
                               headingList=headingList, tipTexts=tipTexts)
    
    self.waiting = False

