
"""
======================COPYRIGHT/LICENSE START==========================

FollowShiftChanges.py: Part of the CcpNmr Analysis program

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

from math import sqrt

from ccpnmr.analysis.popups.EditFitGraph import EditFitGraphPopup
from ccpnmr.analysis.popups.BasePopup import BasePopup

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.CheckButton import CheckButton
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.LabelDivider import LabelDivider
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showWarning, showYesNo
from memops.gui.ProgressBar import ProgressBar
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.core.AssignmentBasic   import propagatePeakAssignments, getResonanceName
from ccpnmr.analysis.core.DataAnalysisBasic import functionFitData, getNmrExpSeriesSampleConditions
from ccpnmr.analysis.core.DataAnalysisBasic import makeDataList, getFitMethodInfo, getFitErrorMethods
from ccpnmr.analysis.core.DataAnalysisBasic import getFitErrorMethod, setFitErrorMethod
from ccpnmr.analysis.core.DataAnalysisBasic import DataFitting
from ccpnmr.analysis.core.ExperimentBasic   import newShiftList, getSpectrumIsotopes, getSpectrumNoise
from ccpnmr.analysis.core.PeakBasic         import findSameAssignmentPeaks, findClosePeaks, arePeaksAssignedSame
#from ccpnmr.analysis.core.MoleculeBasic     import getResidueCode

MAX_STEP_SIZES   = {'1H':0.05,'15N':0.50,'13C':0.50,'13C,15N':0.50}
ISOTOPE_WEIGHTS = {'1H':1.00,'15N':0.15,'13C':0.10,'13C,15N':0.10}

CONC_UNITS = ('M','mM','uM','nM','pM','fM')
CONC_FACTORS = {'M':1.0,
                'mM':0.001,
                'uM':1e-6,
                'nM':1e-9,
                'pM':1e-12,
                'fM':1e-15,}

# TBD 
# Generic parameter bounds
# Fix/remove gradient
# Cope with disappearing peaks
# Macro to export shift deltas per concentration

def FollowTitrationMacro(argServer):

  popup = FollowShiftChangesPopup(argServer.parent)
  popup.open()


def exportSeriesShifts(experimentSeries, sep=''):

  from ccpnmr.analysis.core.AssignmentBasic import getResonanceResidue, getResonanceName

  shiftData = {}
  sampleConditions = {}
  conditionName = experimentSeries.conditionNames[0]
  for experiment in experimentSeries.experiments:

    sampleConditionSet = experiment.sampleConditionSet
    if not sampleConditionSet:
      continue 

    sampleCondition = sampleConditionSet.findFirstSampleCondition(condition=conditionName)
    if not sampleCondition:
      continue

    sampleConditions[sampleCondition] = True

    for spectrum in experiment.dataSources:
      numDim = spectrum.numDim
      peakList = spectrum.activePeakList

      for peak in peakList.peaks:
        key = []
      
        spinSystem = None
        for peakDim in peak.sortedPeakDims():
          for contrib in peakDim.peakDimContribs:
            resonance = contrib.resonance
            spinSystemB = resonance.resonanceGroup
            
            if (spinSystem is None) or spinSystemB:
              spinSystem = spinSystemB
              key.append(resonance)
              
              break

         
        if key:
          key = tuple(key)
          if shiftData.get(key) is None:
            shiftData[key] = {}
            
          shiftData[key][sampleCondition] = ['%.3f' % peakDim.value for peakDim in peak.sortedPeakDims()]       

  conditions = [(sc.value, sc) for sc in sampleConditions.keys()]
  conditions.sort()
  sampleConditions = [x[1] for x in conditions]
  
  keys = []
  maxResonances = 0
  for resonances in shiftData.keys():
    spinSystems = set([r.resonanceGroup for r in resonances if r.resonanceGroup])
    
    if spinSystems:
      ss = spinSystems.pop()
    else:
      ss = None
      
    resNames = []
    for resonance in resonances:
      if resonance.resonanceSet:
        resNames.append(getResonanceName(resonance))
      else:
        resNames.append('%s' % resonance.serial)
    
    maxResonances = max(maxResonances, len(resonances))
    if ss is None:
      keys.append(['?', ss, resNames, resonances])
    elif ss.residue:
      keys.append(['%5.5d' % ss.residue.seqCode, ss, resNames, resonances])
    else:
      keys.append(['{%5.5d}' % ss.serial, ss, resNames, resonances])
  
  keys.sort()

  lines = []


  if sep == '':
    line = '!' + ' ' * 11
    line += ' ' * (5*maxResonances)
    
    for value, sc in conditions:
      line += '%15.15s' % ('%.3f' % value)
  else:
    fields = ['SysCode', 'SeqCode', 'CcpCode']
    fields.extend(['Resonance %d' % (n+1) for n in range(maxResonances)])
    for value, sc in conditions:
      v = '%.3f' % value
      fields.extend(maxResonances*[v])
    line = sep.join(fields)

  lines.append(line)

  for name, ss, resNames, resonances in keys:
    
    seqCode = '-'
    ccpCode = '-'
    
    if ss:
      sysCode = '%d' % ss.serial 
      residue = ss.residue
      if residue:
        seqCode = '%d' % residue.seqCode
        ccpCode = residue.ccpCode # getResidueCode(residue)
    
    else:
      sysCode = 'None'
    
    ppms = []
    for condition in sampleConditions:
      values = shiftData[resonances].get(condition)
      if sep == '':
        # TBD: this assumes that all numDim are the same
        if values is None:
          ppms.append(' '.join(['%7.7s' % '-' for v in range(numDim)]))
        else:
          ppms.append(' '.join(['%7.7s' % v for v in values]))
      else:
        # TBD: should be using maxResonances somehow
        if values is None:
          # TBD: this assumes that all numDim are the same
          ppms.extend(['%7.7s' % '-' for v in range(numDim)])
        else:
          ppms.extend(['%7.7s' % v for v in values])

    if sep == '':
      line = ''
      for text in (sysCode, seqCode, ccpCode):
        line += '%5.5s' % text
    
      for resName in resNames:
        line += '%5.5s' % resName
    
      for ppm in ppms:
        line += ' '
        line += ppm
    else:
      fields = [sysCode, seqCode, ccpCode]
      fields.extend(resNames)
      fields.extend(ppms)
      line = sep.join(fields)
      
    lines.append(line)
   
  return lines
  
class FollowShiftChangesPopup(BasePopup):
  """
  **Follow Chemical Shift Changes During Titrations**
  
  This system is designed to efficiently extract changing peak position data,
  which for example may occur during a titration experiment, and then fit the
  chemical shift changes to an equation curve for the extraction of parameters
  that relate to the peak movements. This analysis may be used to measure such
  things as dissociation constants (e.g. Kd) and temperature coefficients; where
  the position of one or more peaks in a group  of related spectra changes according
  to some experimental condition. Many kinds of experimental condition may be
  varied, but series with changing in concentration (e.g. ligand), temperature
  or pH are commonplace. This system is closely related to the `Follow Intensity
  Changes`_ tool, but here the peak grouping and function fitting is for peak
  locations that do move, rather than stay in the same place and change
  intensity.

  The general idea is that the user sets up an "NMR series" that contains an
  array of experiments where each experiment is one point in the series and
  represents a different value for a parameter (like concentration) being
  investigated. Based initially on the reference peak positions, trajectories
  where picked peak positions move in the different experiments are tracked by
  finding the peak groups that best fits the stated function which relates
  chemical shift to series parameter. It should be noted that NMR series that
  are comprised of experiment planes stacked into a higher-dimensionality
  "pseudo-nD" experiment cannot be used in this analysis. The reason for this
  relates to the way that experiments link to chemical shift lists in the CCPN
  data model; there is no mechanism to record a changing chemical shift within
  a single experiment.

  The layout of the popup window is split into two tabs to reduce clutter. The
  first tab allows the user to setup and adjust all of the options used to
  follow the peak positions as they form "trajectories" and do the function
  fitting. The second tab is used to display the results on a table of peak
  groups, where each group corresponds to a series of peaks which have common
  assignments, with one peak for each experiment.

  The general idea is that the user selects a reference peak list, which will
  give assignment identities to the peak groups being analysed. Typically the
  reference peak list will be from an experiment *in the NMR series*, so that 
  the positions of the peaks, as they move in the different experiments,
  cross or end at the reference position. For proteins this reference will often
  be a 15N HSQC peak list, in which case the analysis operates on peak groups
  that correspond to amides of individual residues. 

  The user chooses an NMR series that has been setup in the `NMR Series`_ popup,
  accessible using the [Edit NMR Series] button. For the analysis to proceed
  properly the selected NMR series will need to contain all of the experiments
  that form part of the titration/analysis and the values of the condition being
  studied (e.g. concentration) must be correctly set. The "Data List Type" and
  corresponding unit indicate what kind of experimental condition/parameter, as
  dictated by the NMR series, will be fitted to chemical shift distance. The
  "Fitting Function" option is adjusted to say what kind of curve should be
  fitted to the peak intensity data; the linear "Ax + B" is common for
  temperature coefficients while "A(B+x-sqrt((B+x)^2-4x))" is used for
  protein-ligand systems in fast exchange, and "A((B+4x-sqrt((B+4x)^2-(4x)^2))/4x - C)"
  is used for monomer-dimer systems in fast exchange. The Error Method
  determines how the errors in the parameters of the fitted function (e.g. error
  in Kd rate) will be calculated. 
  
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

  The peak groups that are analysed for the fit of experimental parameter to
  chemical shift distance may be manually specified by giving a common
  assignment to peaks that derive from the same resonances. Alternatively, an
  automatic method is used to find peaks which are not assigned. This automation
  follows the positions of peaks in their trajectories, choosing the best
  combination of peaks that a) roughly follows a straight line and b) fits the
  selected function equation, in terms of the expected chemical shift distance
  for the experiment. Having the "Assign groups?" option set means that after
  the first peak grouping, peaks will be linked via assignment and subsequent
  peak searches are not generally required. Peak positions may be tracked in one
  or more spectrum dimensions, according to the "Followed Dimensions" selection.
  When multiple dimensions are used, chemical shift difference for dissimilar
  isotopes are combined using the "Shift Weighting" values. The "Max Step Size"
  value are important for the automated peak grouping, given that they limit 
  which peaks are considered when going from one experiment to the next. When
  step sizes are too large the grouping calculation can take a long time. When
  step sizes are too small peaks will be missed and grouping may fail.

  The peak grouping and function fitting is performed using the [Group & Fit
  Peaks] function. After the initial grouping the curve fitting may be
  redone via one of the "Re-fit" buttons; this useful if the fitting
  function is changed. When the curve fitting is done the parameter results from
  the fit, e.g. the "A" and "B" from an "Ax + B" equation, are immediately
  made available from the results table. Also, where relevant, any Kd values are also presented; this requires that the
  binding site concentration was specified. In the "Peak Groups & Analysis"
  table the user can see the fit results and analyse or adjust the peak groups.
  It is commonplace to look through all the intensity curves for each of the
  peak groups by using [Show Fit Graph]; here the user can check how well the
  curve-fit worked and whether any adjustments (e.g. in peak picking) need to be
  made or groups removed. See the `Fit Graph`_ documentation for details about
  how the resultant popup window operates. The "Y" value of the curves come
  from the (isotope weighted) chemical shift distance for each peak of the group
  *along the trajectory from its start point* and the "X" values come from those
  that were entered for the experimental points/planes in the NMR series. When
  the results have been checked, they may be used by directly  exporting the
  fitted parameters from the table.

  **Caveats & Tips**
  
  Each peak group need not contain the same number of peaks if data is missing.
  
  Peaks must be picked in all of the analysed experiments beforehand for this
  system to function.

  Choosing an assigned reference peak list that is postioned in the centre of the
  moving peak trajectories will give the quickest and most reliable peak
  groupings; the trajectory search radius is minimised.

  If there are problems with grouping peaks together the user may assign all
  peaks that ought to go in the same group to the same resonances, for example
  using the "propagate" assignment option, thus connecting peaks together.   

  It is expected that each experiment of the analysed series, because peaks move
  position significantly, will be linked to a separate chemical shift record; so
  that there is a shift value for each condition point. If the experiments of a
  series do no not use separate shift list the user will be propted to set this
  up.

  For groups where the peaks do not move significantly, between experiments in
  the series, a curve will still be fitted to the trajectory. This is because
  the value fitted is a weighted chemical shift distance, from one position to
  the next in the trajectory, and distances will always be positive, including
  when points double-black (within the "Shift Error" tolerance).

  If the peak grouping is taking a significant amount of time, consider reducing
  the "Max Step Size" values; but still leave vales large enough to jump from
  one experiment to the next.

  The fit of the equation curve to the chemical shift changes is naturally
  limited to how many experimental points there are in the series and how well
  spread they are. For example when measuring ligand Kd values, where possible,
  if is best to have some experiments at low concentration and some at high
  concentration, near saturation.

  A subset of peaks in a series may be analysed by reducing the number of peaks
  in the reference peak list. For example the user could make a copy of an HSQC
  peak list and then remove and peak locations that are not required in the
  analysis, e.g. for side chain NH2 peaks or for resonance which don't move
  significantly enough for analysis.
  
  .. _`Follow Intensity Changes`: CalcRatesPopup.html
  .. _`NMR Series`: EditExperimentSeriesPopup.html
  .. _`Fit Graph`: EditFitGraphPopup.html
  """

  def __init__(self, parent, *args, **kw):

    self.waiting        = False
    self.guiParent      = parent
    self.fitFunc        = 1
    self.peakList       = None
    self.doAssignGroups = True
    self.fitGraphPopup  = None
    self.expSeries      = None
    self.peakGroup      = None
    self.peakGroups     = []
    self.dataList       = None
    self.noise          = 0.1
    self.dataDims       = []
   
    BasePopup.__init__(self, parent=parent, title='Data Analysis : Follow Shift Changes', **kw)

    #self.updateAfter()
 
  def body(self, guiFrame):

    self.geometry('700x400')
  
    guiFrame.grid_columnconfigure(0, weight=1)
    guiFrame.grid_rowconfigure(0, weight=1)
   
    tipTexts = ['Allows the user to setup how resonance trajectories are followed and which experiments to analyse',
                'The peaks that have been grouped into resonance trajectories, and the various equation parameters estimated for each']
    options = ['Settings','Peak Groups & Analysis']
      
    tabbedFrame = TabbedFrame(guiFrame, options=options, grid=(0,0), tipTexts=tipTexts)
    self.tabbedFrame = tabbedFrame
    frameA, frameB = tabbedFrame.frames

    #
    # Settings
    #
        
    frameA.expandGrid(6,0)
    
    row = 0
    div = LabelDivider(frameA, text='Experiment Setup', grid=(row,0))
    
    row +=1
    frame = Frame(frameA, grid=(row,0))
    frame.grid_columnconfigure(3, weight=1)

    label = Label(frame, text='Reference Peak List:', grid=(0,0), sticky='e')
    tipText = 'Selects which peak list is the source of assignments and representative positions for peak group trajectories (need not be at the end of the trajectory)'
    self.peakListPulldown = PulldownList(frame, self.changePeakList,
                                         grid=(0,1), tipText=tipText)

    label = Label(frame, text=' NMR Experiment Series:', grid=(1,0), sticky='e')
    tipText = 'Selects which series of NMR experiments to perform the peak-following and equation fitting analysis for'
    self.expSeriesPulldown = PulldownList(frame, self.changeExpSeries,
                                          grid=(1,1), tipText=tipText)

    label = Label(frame, text='  Data List Type:', grid=(0,2), sticky='e')
    tipText = 'Displays the kind of values described by the NMR series being fitted, e.g. temperature, ligand concentration...'
    self.dataListTypeLabel = Label(frame, text='', grid=(0,3), tipText=tipText)

    label = Label(frame, text='Data List Unit:', grid=(1,2), sticky='e')
    tipText = 'Displays the units of measurement for the values described by the NMR series being fitted'
    self.dataListUnitLabel = Label(frame, text='', grid=(1,3), tipText=tipText)
    
    Label(frame, text='', grid=(2,0))

    row +=1
    div = LabelDivider(frameA, text='Function Fitting')
    div.grid(row=row,column=0, sticky='ew')

    row +=1
    frame = Frame(frameA, grid=(row,0))
    frame.grid_columnconfigure(5, weight=1)

    label = Label(frame, text='Fitting Function:',
                  grid=(0,0), sticky='e')
    tipText = 'Selects which kind of parameterised equation to use in the fitting of chemical shift distance to NMR series value'
    self.fitFuncPulldown  =  PulldownList(frame, self.changeFitFunc,
                                          grid=(0,1), gridSpan=(1,2), tipText=tipText)

    label = Label(frame, text='Followed Dimensions:',
                  grid=(1,0), sticky='e')
    tipText = 'Selects which dimensions of the spectra to follow chemical shifts for; can use different isotope types given the stated weightings'
    self.dataDimsPulldown = PulldownList(frame, self.changeDataDim,
                                         grid=(1,1), tipText=tipText)

    label = Label(frame, text='Error method:', grid=(2,0), sticky='e')
    fitErrorMethods = getFitErrorMethods()
    fitErrorMethod = getFitErrorMethod(self.analysisProject)
    ind = fitErrorMethods.index(fitErrorMethod)
    tipText = 'Selects how errors in the fit of the selected equation to the chemical shift data are estimated'
    self.fitErrorPulldown = PulldownList(frame, callback=self.setFitErrorMethod,
                                         texts=fitErrorMethods, index=ind, grid=(2,1),
                                         tipText=tipText)

    label = Label(frame, text='Assign groups?',
                  grid=(3,2), sticky='e')
    tipText = 'Whether to assign all peaks in the same group to the same resonances; where possible inherited from the reference peak list'
    self.assignSelect = CheckButton(frame, selected=-True,
                                    grid=(3,3), tipText=tipText)

    label = Label(frame, text='  Ignore Zero Merit Peaks?',
                  grid=(4,2), sticky='e')
    tipText = 'When grouping peaks along resonance trajectories, whether to ignore peaks with zero figure-of-merit value'
    self.meritSelect = CheckButton(frame, selected=True,
                                   grid=(4,3), tipText=tipText)

    label = Label(frame, text='  Shift Error (points):',
                  grid=(3,0), sticky='e')
    tipText = 'The amount of grace, in data points, for grouping resonance locations that do not move significantly; allows trajectories to backtrack a little'
    self.shiftErrorEntry = FloatEntry(frame, text=4.0, width=8,
                                      grid=(3,1), tipText=tipText)

    self.proteinConcLabel = Label(frame, text='Binding site concentration:',) # No grid
    tipText = 'For fitting equation that allow calculation of dissociation constants, Kd, the known concentration of the site being reversibly bound to'
    self.proteinConcEntry = FloatEntry(frame, text=None, width=8, tipText=tipText,
                                       returnCallback=self.updateAfter) # No grid
    self.proteinConcEntry.bind('<Leave>', self.updateAfter, '+')
    
    self.unitLabel = Label(frame, text=' unit:',) # No grid
    tipText = 'The units that the binding site concentration is given in'
    self.unitPulldown = PulldownList(frame, texts=CONC_UNITS, tipText=tipText) # No grid

    label = Label(frame, text='Min F2/F1 Grad:', grid=(4,0), sticky='e')
    tipText = 'The lower limit of the chemical shift trajectory gradient; difference in second dimension over difference in first dimension'
    self.minGradientEntry = FloatEntry(frame, text=None, width=8,
                                      grid=(4,1), tipText=tipText)

    label = Label(frame, text='Max F2/F1 Grad:', grid=(5,0), sticky='e')
    tipText = 'The upper limit of the chemical shift trajectory gradient; difference in second dimension over difference in first dimension'
    self.maxGradientEntry = FloatEntry(frame, text=None, width=8,
                                       grid=(5,1), tipText=tipText)
     
    Label(frame, text='', grid=(6,0))
   
    row +=1
    div = LabelDivider(frameA, text='Isotope Parameters')
    div.grid(row=row,column=0, sticky='ew')
    
    row +=1
    frame = Frame(frameA, grid=(row,0))
    frame.expandGrid(1,6)

    label = Label(frame, text='Shift Weighting: ', grid=(0,0))
     
    label = Label(frame, text=' 1H', grid=(0,1))
    tipText = 'The scaling factor for 1H dimensions, used to give equivalency to ppm distances for different kinds of isotope'
    self.tolHydrogenEntry = FloatEntry(frame, text=ISOTOPE_WEIGHTS['1H'],
                                       width=8, grid=(0,2), tipText=tipText)

    label = Label(frame, text=' 15N', grid=(0,3))
    tipText = 'The scaling factor for 15N dimensions, used to give equivalency to ppm distances for different kinds of isotope; default is relative to a 1H weight of 1.0'
    self.tolNitrogenEntry = FloatEntry(frame, text=ISOTOPE_WEIGHTS['15N'],
                                       width=8, grid=(0,4), tipText=tipText)

    label = Label(frame, text=' 13C', grid=(0,5))
    tipText = 'The scaling factor for 13C dimensions, used to give equivalency to ppm distances for different kinds of isotope; default is relative to a 1H weight of 1.0'
    self.tolCarbonEntry = FloatEntry(frame, text=ISOTOPE_WEIGHTS['13C'],
                                     width=8, grid=(0,6), tipText=tipText)
     
    label = Label(frame, text='Max Step Size: ', grid=(1,0))

    label = Label(frame, text=' 1H', grid=(1,1))
    tipText = 'The 1H ppm difference limit, within which each subsequent point along a resonance trajectory may be found; to make peak groups'
    self.stepHydrogenEntry = FloatEntry(frame, text=MAX_STEP_SIZES['1H'],
                                        width=8, grid=(1,2), tipText=tipText)

    label = Label(frame, text=' 15N', grid=(1,3))
    tipText = 'The 15N ppm difference limit, within which each subsequent point along a resonance trajectory may be found; to make peak groups'
    self.stepNitrogenEntry = FloatEntry(frame, text=MAX_STEP_SIZES['15N'],
                                        width=8, grid=(1,4), tipText=tipText)

    label = Label(frame, text=' 13C', grid=(1,5))
    tipText = 'The 13C ppm difference limit, within which each subsequent point along a resonance trajectory may be found; to make peak groups'
    self.stepCarbonEntry = FloatEntry(frame, text=MAX_STEP_SIZES['13C'],
                                      width=8, grid=(1,6), tipText=tipText)

    #
    # Peak Groups
    #
        
    frameB.expandGrid(0,0)
    
    headingList, tipTexts = self.getHeadingList(3)
    self.scrolledMatrix = ScrolledMatrix(frameB,
                                         headingList=headingList,
                                         callback=self.selectGroup,
                                         multiSelect=True,
                                         deleteFunc=self.removeGroup, 
                                         grid=(0,0), tipTexts=tipTexts)

    row +=1
    tipTexts = ['',
                '',
                '',
                '']
    texts = ['Remove\nSelected Groups','Re-fit\nSelected',
             'Show Fit\nGraph','Show\nPeaks']
    commands = [self.removeGroup, self.calculateFit,
                self.showFunctionFit ,self.showPeaks]
                
    self.bottomButtons = ButtonList(frameB,  commands=commands, texts=texts,
                                    grid=(row,0), tipTexts=tipTexts)

    #
    # Main
    #
        
    tipTexts = ['',
                '',
                '',
                '']
    texts = ['Group & Fit Peaks', 'Re-fit All Groups',
             'Edit NMR Series', 'Export Shifts']
    commands = [self.groupPeaks, self.calculateAllFits,
                self.editExpSeries, self.exportShifts]
    self.buttonList = ButtonList(guiFrame, grid=(1,0), tipTexts=tipTexts,
                                 commands=commands, texts=texts)
 
 
    buttons = UtilityButtonList(tabbedFrame.sideFrame, helpUrl=self.help_url,
                                grid=(0,0), sticky='e')

    self.updateExpSeries()
    self.updatePeakLists()
    self.updateDataDims()
    self.updateFitFuncs()
    self.updateAfter()

    self.administerNotifiers(self.registerNotify)
  
  def administerNotifiers(self, notifyFunc):

    for func in ('__init__','delete','setName'):
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.Experiment', func)
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.DataSource', func)

    for func in ('__init__','delete'):
      notifyFunc(self.updatePeakLists, 'ccp.nmr.Nmr.PeakList', func)
    
    notifyFunc(self.updatePeaks, 'ccp.nmr.Nmr.Peak', 'delete')
      
    for func in ('__init__','delete','setVariableUnit',
                 'setVariableDescriptor','setVariableName','setName',
                 'setExperiment','addExperiment','removeExperiment'):
      notifyFunc(self.updateExpSeries, 'ccp.nmr.Nmr.NmrExpSeries', func)

    for func in ('setConditionVaried',): 
      notifyFunc(self.updateExpSeries,'ccp.nmr.Nmr.SampledDataDim', func)

  def open(self):
  
    self.updateExpSeries()
    self.updatePeakLists()
    self.updateDataDims()
    self.updateFitFuncs()
    self.updateAfter()
    BasePopup.open(self)

  def close(self):
  
    BasePopup.close(self)

  def setFitErrorMethod(self, method):

    setFitErrorMethod(self.analysisProject, method)

  def exportShifts(self):
  
    from memops.gui.FileSelect import FileType
    from memops.gui.FileSelectPopup import FileSelectPopup
    from memops.gui.MessageReporter import showMulti

    texts = ['Fixed format', 'Tab separated', 'Comma separated', 'Cancel']
    objects = ['', '\t', ',', None]
    sep = showMulti('Export format', 'Select export format:', texts, objects, self)
    if sep is None:
      return

    if self.expSeries:
      lines = exportSeriesShifts(self.expSeries, sep)
    
    
    if lines:
      fileTypes = [  FileType('All', ['*'])]
      fileSelectPopup = FileSelectPopup(self, file_types = fileTypes,
                 title = 'Export Series Shifts', dismiss_text = 'Cancel',
                 selected_file_must_exist = False)

      fileName = fileSelectPopup.getFile()
      if fileName:
        file = open(fileName, 'w')
        
        for line in lines:
          file.write(line)
          file.write('\n')
        
        file.close()  

  def calculateAllFits(self):
  
    self.tabbedFrame.select(1)
  
    self.updateIsotopeWeights()
    self.updateMaxStepSizes()
  
    if self.peakGroups and self.expSeries:
      if not self.fitFunc:
        showWarning('Warning','No fitting function selected', parent=self)
        return
  
      failedGroups = []
      progressBar = ProgressBar(self, text="Calculating Fit ",total = len(self.peakGroups))
      for group in self.peakGroups:

        try:
          self.fitPeakGroup(group)
          progressBar.increment()
        except:
          failedGroups.append(group.refObject)
  
      progressBar.destroy()
    
      if failedGroups:
        self.warnFailedPeakGroups(failedGroups)
    
    self.updateAfter()
   
   
  def calculateFit(self):
  
    self.updateIsotopeWeights()
    self.updateMaxStepSizes()
  
    if self.peakGroups and self.expSeries:
      if not self.fitFunc:
        showWarning('Warning','No fitting function selected', parent=self)
        return
  
      failedGroups = []
      progressBar = ProgressBar(self, text="Calculating Fit ",total = len(self.peakGroups))
      for group in self.scrolledMatrix.currentObjects:

        try:
          self.fitPeakGroup(group)
          progressBar.increment()
        except:
          failedGroups.append(group.refObject)
 
      progressBar.destroy()
      
      if failedGroups:
        self.warnFailedPeakGroups(failedGroups)
    
    self.updateAfter()

  def warnFailedPeakGroups(self, failedPeaks):

    numFail   = len(failedPeaks)
    positions = []
    
    for peakF in failedPeaks:
      values = ['%.3f' % pd.value for pd in peakF.sortedPeakDims() if pd.dataDimRef]
      positions.append('(' + ','.join(values) + ')')
    
    failText = ', '.join(positions)
    showWarning('Fit failure','Fit failed for %d groups at reference positions %s' % (numFail,failText), parent=self)

  def fitPeakGroup(self, dataFitting, func=None):
  
    if func is None:
      func = self.fitFunc
    
    peaks = dataFitting.objects
    if (len(peaks) < 2) or ((func > 1) and (len(peaks) < 3)):      
      dataFitting.resetFit()
      return
      
    dataFitting = self.getPeakSeriesData(dataFitting)
    dataFitting.fitFunction  = func
    dataFitting.fitErrorFunction = getFitErrorMethod(self.analysisProject)
    dataFitting.noiseLevel   = self.noise

    dataFitting.fit()

  def updateIsotopeWeights(self):

    ISOTOPE_WEIGHTS['1H']  = self.tolHydrogenEntry.get() or 1.0
    ISOTOPE_WEIGHTS['15N'] = self.tolNitrogenEntry.get() or 0.15
    ISOTOPE_WEIGHTS['13C'] = self.tolCarbonEntry.get() or 0.1
    ISOTOPE_WEIGHTS['13C,15N'] = ISOTOPE_WEIGHTS['13C']
    ISOTOPE_WEIGHTS['15N,13C'] = ISOTOPE_WEIGHTS['13C']

  def updateMaxStepSizes(self):

    MAX_STEP_SIZES['1H']  = self.stepHydrogenEntry.get() or 1.0
    MAX_STEP_SIZES['15N'] = self.stepNitrogenEntry.get() or 0.15
    MAX_STEP_SIZES['13C'] = self.stepCarbonEntry.get() or 0.1
    MAX_STEP_SIZES['13C,15N'] = MAX_STEP_SIZES['13C']
    MAX_STEP_SIZES['15N,13C'] = MAX_STEP_SIZES['13C']
  
  def getBindingSiteConc(self):
  
    unitFactor = CONC_FACTORS[self.unitPulldown.getText()]
    proteinConc = self.proteinConcEntry.get() or 1.0
    proteinConc *= unitFactor

    return proteinConc
  
  def getPeakSeriesData(self, dataFitting):

    peaks = dataFitting.objects 
    pointError = self.shiftErrorEntry.get() or 0.0
  
    xData   = []
    yData   = []
    yDelta = 0.0
    xWidths = [] # From expt series
    yWidths = [] # Based upon the digital resolution at the moment

    proteinConc = self.getBindingSiteConc()
    fitFunc = self.fitFunc

    if self.peakList and self.dataDims and self.expSeries:
    
      yData   = [0.0,]
      dimNums = [dataDim.dim for dataDim in self.dataDims]
      sampleConditions = self.getSampleConditions()
      experimentDict = {}
      for sampleCondition in sampleConditions:
        for experiment in sampleCondition.sampleConditionSet.experiments:
          experimentDict[experiment] = sampleCondition
     
      locations = []
      for peak in peaks:
        experiment = peak.peakList.dataSource.experiment
        sampleCondition = experimentDict.get(experiment)
 
        if sampleCondition:
          xValue = sampleCondition.value
          
          if fitFunc in (7,9):
            xValue *= CONC_FACTORS.get(sampleCondition.unit, 1.0)
          
          if (fitFunc == 9) and (sampleCondition.unit in CONC_FACTORS): # TBD: Check this is the only one
            xValue /= proteinConc
            
          xData.append(xValue)
          xWidths.append(sampleCondition.error)

        else:
          xData.append(0.0)
          xWidths.append(0.0)
          
        location = []
        yWidth   = 0.0
        for peakDim in peak.sortedPeakDims():
          dataDimRef = peakDim.dataDimRef
          i = peakDim.dim
          if i in dimNums:
            isotope = ','.join(dataDimRef.expDimRef.isotopeCodes)
            weight  = ISOTOPE_WEIGHTS.get(isotope, 1.0)
            error   = weight * dataDimRef.valuePerPoint * pointError
            yWidth += error * error
            location.append((peakDim.value, weight))
        
        yWidths.append(sqrt(yWidth))
        locations.append(location)
      
      for i in range(1,len(locations)):
        
        location1 = locations[i]
        location0 = locations[i-1]
        shiftDist = 0.0
        
        for j in range(len(location1)):
          value1, weight1 = location1[j]
          value0, weight0 = location0[j]
         
          diff = (value1 - value0) * weight1
          shiftDist += diff * diff
      
        yData.append( yData[-1] + sqrt(shiftDist) )

      if len(locations) > 1:

        location1 = locations[-1]
        location0 = locations[0]
        shiftDist = 0.0
        
        for j in range(len(location1)):
          value1, weight1 = location1[j]
          value0, weight0 = location0[j]
         
          diff = (value1 - value0) * weight1
          shiftDist += diff * diff
      
        yDelta = sqrt(shiftDist)

    if not dataFitting:
      dataFitting = DataFitting(guiParent=self)
      
    dataFitting.dataX = xData
    dataFitting.dataY = yData
    dataFitting.deltaY = yDelta
    dataFitting.errorsX = xWidths
    dataFitting.errorsY = yWidths
     
    return dataFitting

  def getSampleConditions(self):
    
    sampleConditions = []
    
    if self.expSeries:
      conditionName = self.expSeries.conditionNames[0]
      sampleConditions = getNmrExpSeriesSampleConditions(self.expSeries).get(conditionName, [])
     
    return sampleConditions
    
  def groupPeaks(self):

    self.tabbedFrame.select(1)

    self.updateMaxStepSizes()
    self.minGrad = self.minGradientEntry.get()
    self.maxGrad = self.maxGradientEntry.get()
    
    # Score using regularity of spacing - optional switch?
    
    noise        = self.noise
    fitFunc      = self.fitFunc
    fitErrorFunc = getFitErrorMethod(self.analysisProject)
    pointError   = self.shiftErrorEntry.get() or 0.0

    if not self.fitFunc:
      showWarning('Warning','No fitting function selected', parent=self)
      return
    
    def getAlignScore(peaks):
      dataFitting = DataFitting(objects=peaks, guiParent=self)
      dataFitting = self.getPeakSeriesData(dataFitting)
      dataFitting.noiseLevel   = noise
      dataFitting.fitFunction  = fitFunc
      dataFitting.fitErrorFunction = fitErrorFunc
      try:
        dataFitting = functionFitData(dataFitting)
      except Exception, e:
        ###showWarning('Function fit failure', e, parent=self)
        print e
        return 1.0e20  # arbitrary large number

      fitError = dataFitting.fitError
      dataFitting.delete()
      
      locations = []
      for peak0 in peaks:
        location = [peakDim.value for peakDim in peak0.sortedPeakDims()]
        locations.append(location)
      
      N = len(locations)-2
      if N > 0:
        M = len(locations[0])
        meanCos = 0.0

        for i in range(N):
          vector1 = []
          vector2 = []
          sumSq1  = 0.0
          sumSq2  = 0.0
          product = 0.0

          for j in range(M):
            coord1   = locations[i+1][j] - locations[i][j]
            coord2   = locations[i+2][j] - locations[i][j]
            sumSq1  += coord1 * coord1
            sumSq2  += coord2 * coord2
            product += coord1 * coord2
      
          if sumSq1 and sumSq2:
            meanCos += product/( sqrt(sumSq1) * sqrt(sumSq2) )
      
        if meanCos:
          fitError /= meanCos * meanCos
      
      else:
        fitError = 1.0
        
      return fitError
                  
    def getRoutes(refPeak, peakSets, tolerances, routes, weights, route=None, i=0, dists=None):

      
      if not route:
        route = []
      
      if not dists:
        dists = []
                        
      N = len(peakSets)
      if i >= N:
        routes.append(route)
        return
        
      peaks = peakSets[i]   

      # if peaks is empty but the next one is not...
      # carry on but adjust the valid distance 

      if (not peaks) and (i < N):
        route2 = list(route)
        dists2 = list(dists)
        getRoutes(refPeak, peakSets, tolerances, routes, weights, route2, i+1, dists2)
      

      for peak in peaks:
        
        peakDims = peak.sortedPeakDims()
          
        distE  = None
        tolE   = None
        deltas = []
        route2 = list(route)
        route2.append(peak)
        dists2 = list(dists)

        if len(route2) > 1:
          ppms1 = [pd.value for pd in route2[ 0].sortedPeakDims()]
          ppms2 = [pd.value for pd in route2[-2].sortedPeakDims()]
          ppms3 = [pd.value for pd in route2[-1].sortedPeakDims()]
        
          distE = 0.0
          tolE = 0.0
          for j in range(len(ppms1)):
            weight = weights[j]
            deltaE = ppms1[j] - ppms3[j]
            deltaN = ppms2[j] - ppms3[j]
            distE += deltaE * deltaE * weight * weight
            deltas.append(abs(deltaN))
            errorE = peakDims[j].dataDimRef.valuePerPoint * pointError * weight
            tolE += errorE * errorE
        
          distE = sqrt(distE)
          tolE  = sqrt(tolE)
          
          if (self.maxGrad is not None) or (self.minGrad is not None):
            if deltas[0]:
              grad = (deltas[1]*weights[1])/(deltas[0]*weights[0])
              if (self.maxGrad is not None) and (grad > self.maxGrad):
                continue
                
              if (self.minGrad is not None) and (grad < self.minGrad):
                continue
 
            else:
              if (self.maxGrad is not None):
                continue
        
        if distE is not None:
          if dists and ( (dists[-1] - distE) > tolE ): # backward check
            if (len(peaks) != 1) or ( not arePeaksAssignedSame(peak, refPeak) ):
              continue
        
          dists2.append(distE)

        for j in range(len(deltas)):
          valuePerPoint = peakDims[j].dataDimRef.valuePerPoint
          if (deltas[j] - tolerances[j]) > valuePerPoint * pointError:
            if (len(peaks) != 1) or ( not arePeaksAssignedSame(peak, refPeak) ):
              route2 = list(route)
              dists2 = list(dists)
              getRoutes(refPeak, peakSets, tolerances, routes, weights, route2, i+1, dists2)
              break

        else:
          getRoutes(refPeak, peakSets, tolerances, routes, weights, route2, i+1, dists2)
     
     
      return routes


    if not (self.expSeries and self.peakList):
      return
  
    if self.assignSelect.get():
      self.checkShiftLists()
  
    self.clearPeakGroups()

    sampleConditions = [(sc.value, sc) for sc in self.getSampleConditions()]
    sampleConditions.sort()
    
    experimentDict = {}
    for value, sampleCondition in sampleConditions:
      for experiment in sampleCondition.sampleConditionSet.experiments:
        experimentDict[experiment] = sampleCondition
    
    experiments = self.expSeries.experiments
    refIsotopes = getSpectrumIsotopes(self.peakList.dataSource)

    peakLists   = {}   
    for experiment in experiments:
      sampleCondition = experimentDict.get(experiment)
      if not sampleCondition:
        continue
    
      for spectrum in experiment.dataSources:
        if (spectrum.dataType == 'processed') and ( getSpectrumIsotopes(spectrum) == refIsotopes):
          peakList = spectrum.activePeakList
          if peakList:
            peakLists[sampleCondition] = peakList
            break
    
    sortedPeakLists = []
    for value, sampleCondition in sampleConditions:
      peakList = peakLists.get(sampleCondition)
      if peakList:
        sortedPeakLists.append(peakList)
    
    refIndex = sortedPeakLists.index(self.peakList) 
    
    searchWidths = []
    for isotope in refIsotopes:
      delta = MAX_STEP_SIZES.get(isotope, 0.1)
      searchWidths.append(delta)
    
    def isPeakAssigned(peak):
      for peakDim in peak.peakDims:
        if peakDim.peakDimContribs:
          return True
    
      return False
        
    # Find candidate peaks for following a titration from a reference peak
    mapping = {}
    candidates = {}
    for peak in self.peakList.peaks:      
      peaks0 = []     
      i = 0
      for peakList in sortedPeakLists:
        if peakList is self.peakList:
          peaks2 = [peak,]
        
        else:  
          peaks2 = findSameAssignmentPeaks(peak, peakList) # Could get more than one
          peaks4 = []

          if not peaks2:
            intervals = abs(refIndex-i)
            tolerances = [x*intervals for x in searchWidths]
            peaks3 = findClosePeaks(peak, peakList, tolerances=tolerances, pickNewPeaks=False)
            peaks2 = []

            for peak3 in peaks3:
              if not isPeakAssigned(peak3):
                peaks2.append(peak3)
      
        peaks0.append(peaks2)
        i += 1
                  
      candidates[peak] = peaks0
           
    progressBar = ProgressBar(self, text="Grouping Peaks ",total = len(self.peakList.peaks))
    
    grouping = {}
    peakGroups = {}
    for peak in self.peakList.peaks:      
      peaks      = []
      weights    = []
      tolerances = []
      peakGroups[peak] = {}
      
      for peakDim in peak.sortedPeakDims():
        isotope = ','.join(peakDim.dataDimRef.expDimRef.isotopeCodes)
        tolerances.append(MAX_STEP_SIZES.get(isotope))
        weights.append(ISOTOPE_WEIGHTS.get(isotope))
     
      
      routes = getRoutes(peak, candidates[peak],tolerances,[],weights)
      
            
      if routes:
        bestRoute = routes[0]
        bestScore = getAlignScore(bestRoute)
        
        for route in routes:
          score = getAlignScore(route)
                         
          if score < bestScore:
            bestScore = score
            bestRoute = route
      
        peakGroups[peak] = bestRoute
        
        for peak1 in bestRoute:
          if grouping.get(peak1) is None:
            grouping[peak1] = []
          
          grouping[peak1].append([bestScore, peak])
          
      
      progressBar.increment()

    progressBar.destroy()
    
    for peak1 in grouping.keys():
      matches = grouping[peak1]
      matches.sort()
      grouping[peak1] = matches[0][1]

    ignorePoorMerit = self.meritSelect.get()
    for peak in self.peakList.peaks:
      peaks1 = peakGroups[peak]
      peaks2 = []
      
      for peak1 in peaks1:
        if grouping[peak1] is peak:
          peaks2.append(peak1)
      
      if ignorePoorMerit:
        peaks2 = [pk for pk in peaks2 if pk.figOfMerit]
        
      peakGroups[peak] = peaks2
      dataFitting = DataFitting(fitFunction=self.fitFunc, refObject=peak, objects=peaks2, guiParent=self)
      self.peakGroups.append(dataFitting)
    
    if self.assignSelect.get():
      for peak in self.peakList.peaks:
        propagatePeakAssignments(peakGroups[peak], refPeak=peak, cleanNonRef=True)
      
    self.calculateAllFits()

    

  def editExpSeries(self):
  
    self.guiParent.editExpSeries()


  def checkShiftLists(self):
  
    dict = {}
    
    experiments = self.expSeries.experiments
    for experiment in experiments:
      shiftList = experiment.shiftList
      if shiftList:
        if not dict.get(shiftList):
          dict[shiftList] = []
      
        dict[shiftList].append(experiment)
    
    refExperiment = self.peakList.dataSource.experiment
    
    if len(experiments) != len(dict.keys()):
      msg  = 'Experiments in series do not have individual shift lists. '
      msg += 'Link each experiment to its own shift list?'
      if showYesNo('Query', msg, parent=self):
        
        for shiftList in dict.keys():
          experiments1 = dict[shiftList]
          if len(experiments1) > 1:
 
            if refExperiment in experiments1:
              experiment1 = refExperiment
            else:
              experiment1 = experiments1[0]
 
            for experiment2 in experiments1:
              if experiment2 is not experiment1:
                experiment2.shiftList = newShiftList(refExperiment.root, unit=shiftList.unit)

        for experiment in experiments:
          if not experiment.shiftList:
            experiment.shiftList = newShiftList(refExperiment.root)
          

  def getDataListInfo(self):
  
    name = ''
    unit = ''
    
    if self.expSeries:
      name = '%s series' % (','.join(self.expSeries.conditionNames))
      unit = 'None'
      dict = getNmrExpSeriesSampleConditions(self.expSeries)
      
      for conditionName in self.expSeries.conditionNames:
        for sampleCondition in dict.get(conditionName, []):
          unit = sampleCondition.unit
          break
          
        else:
          continue
        break  
  
    return name, unit
   

  def updateDataListInfo(self):
  
    name, unit = self.getDataListInfo()
  
    self.dataListTypeLabel.set(name)
    self.dataListUnitLabel.set(unit)
      
        
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
      
  def removeGroup(self, *event):
  
    peakGroups = self.scrolledMatrix.currentObjects
    
    if peakGroups:
      for peakGroup in peakGroups:
        self.peakGroups.remove(peakGroup)
        peakGroup.delete()
      
      self.peakGroup = None
      self.updateAfter()
  
  def getPeakLists(self):

    peakLists = []
    
    if self.expSeries:
      experiments = self.expSeries.experiments
    else:
      experiments = self.nmrProject.experiments
    
    for experiment in experiments:
      for spectrum in experiment.dataSources:
        if (spectrum.dataType == 'processed') and (spectrum.numDim > 1):
          if spectrum.findFirstDataDim(className='SampledDataDim'):
            continue
            
          for peakList in spectrum.peakLists:
            peakLists.append(peakList)
  
    return peakLists
  
  
  def changePeakList(self, peakList):
      
    if peakList is not self.peakList:
      self.peakList = peakList
      self.updateDataDims()
      self.updateButtons()
  
      
  def getFitFuncs(self):
  
    fitFunctions = [x[0] for x in getFitMethodInfo()]
    fitFunctions = ['<None>'] + fitFunctions
    
    return fitFunctions

  def changeFitFunc(self, funIndex):
  
    self.fitFunc = funIndex
    
    if self.fitFunc in (7,9):
      self.proteinConcLabel.grid(row=5, column=2, sticky='w')
      self.proteinConcEntry.grid(row=5, column=3, sticky='w')
      self.unitLabel.grid(row=5, column=4, sticky='w')
      self.unitPulldown.grid(row=5, column=5, sticky='w')
   
    else:
      self.proteinConcLabel.grid_forget()
      self.proteinConcEntry.grid_forget()
      self.unitLabel.grid_forget()
      self.unitPulldown.grid_forget()
  
  def updateFitFuncs(self, *opt):
  
    names = self.getFitFuncs()
    if self.fitFunc > len(names):
      self.changeFitFunc(1)
    
    indices = range(len(names))
    
    self.fitFuncPulldown.setup(names, indices, self.fitFunc)
  
  
  def getHeadingList(self, numDim):
  
    fitFunc = self.fitFunc
  
    headingList = ['#',]
    tipTexts = ['Number of the peak group',]
    
    for i in range(numDim):
      headingList.append( 'Assign\nF%d' % (i+1) )
      tipTexts.append('The assignment of the reference peak for the group in the F%d dimension' % (i+1))
    
    headingList.extend(['Traj\nDist','Shift\nDist','Fit\nError','Num\nPeaks','Fitted\nFunction'])
    tipTexts.extend(['The isotope-weighted chemical shift path length; following each pair of points in the trajectory of the peak group',
                     'The isotope-weighted chemical shift distance between the first and last point in the peak group',
                     'The error in the goodness fo fit of the selected graph to the observed chemical shift data',
                     'The number of peaks contained in the analysis group',
                     'The kind of function used in the graph fitting, to extract parameters, for this group',])
    
    if fitFunc in (7,9):
       headingList.append('Kd')
       headingList.append('Kd\nerror')
       tipTexts.append('The estimated dissociation constant for binding, if using certain types of fit equation')
       tipTexts.append('The estimated error in the dissociation constant for binding')
    
    numParams = getFitMethodInfo()[fitFunc-1][1]
    for i in range(numParams):
      letter = chr(ord('A') + i)
      headingList.append('Fit\nParam %s' % letter)
      tipTexts.append('The estimated value of parameter "%s", obtained by fitting the selected equation to the chemical shift data' % letter)
      
    for i in range(numParams):
      letter = chr(ord('A') + i)
      headingList.append('Param\nError %s' % letter)
      tipTexts.append('The estimated error, using the selected method, in fit parameter "%s"' % letter)

    return headingList, tipTexts

  def showFunctionFit(self, noOpen=False):
  
    if self.peakGroup:
      peaks = self.peakGroup.objects
      peak  = self.peakGroup.refObject
      func  = self.peakGroup.fitFunction
      
      if not peaks:
        return
      
      title  = '.'.join([pd.annotation or '?' for pd in peak.sortedPeakDims()])
      xLabel = ','.join(self.expSeries.conditionNames)
      yLabel = 'Shift Distance'
 
      if self.fitGraphPopup:
        self.fitGraphPopup.update(dataFitting=self.peakGroup,
                                  xLabel=xLabel, yLabel=yLabel,
                                  graphTitle=title, force=True)
        if not noOpen:
          self.fitGraphPopup.open()
      else:
        self.fitGraphPopup = EditFitGraphPopup(self, self.peakGroup,
                                               self.getPeakSeriesData, self.updateGroup,
                                               xLabel, yLabel,
                                               nextSetFunction=self.showNextGroup,
                                               prevSetFunction=self.showPrevGroup,
                                               graphTitle=title)
 
                                             
  def updateGroup(self, dataFitting):

    if self.peakGroup:
      self.updateAfter()
      
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

  def updateAfter(self, *obj):
    
    if self.waiting:
      return
      
    else:
      self.waiting = True
      self.after_idle(self.update)

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

  def getPeakListName(self, peakList):
  
    spectrum = peakList.dataSource
    return '%s:%s:%d' % (spectrum.experiment.name,spectrum.name,peakList.serial)
 
  
  def updatePeakLists(self, *peakList):
  
    index = 0
    names = []
    peakLists = self.getPeakLists()
    
    if peakLists:
      names = [self.getPeakListName(pl) for pl in peakLists]
      
      if self.peakList not in peakLists:
        self.peakList = peakLists[0]
  
      index = peakLists.index(self.peakList)
      
      self.updateDataDims()
      
    else:
      self.peakList = None
      
    self.peakListPulldown.setup(names, peakLists, index)


  def updateButtons(self):
  
   midButtons = self.buttonList.buttons
   botButtons = self.bottomButtons.buttons

   midButtons[2].enable() # Cam always create new if none exists.

   if self.peakList and self.expSeries:
     midButtons[0].enable()
     
     if self.peakGroups:
       midButtons[1].enable()
       midButtons[3].enable()
     else:
       midButtons[1].disable()
       midButtons[3].disable()
       
     if self.peakGroup and self.peakGroups:
       for i in range(4):
         botButtons[i].enable()
 
     else:
       for i in range(4):
         botButtons[i].disable()
       
   else:
     midButtons[0].disable()
     midButtons[1].disable()
     midButtons[3].disable()
     for i in range(4):
       botButtons[i].disable()  
  

  def update(self):
  
    self.updateButtons()
    
    # wb104: 26 Mar 2014: below causes exception when switch between expSeries 
    # (that exception only happens if the functionFit is open)
    # this means that some graphs might need to be updated manually
    #if allowShowFunctionFit and self.fitGraphPopup:
    #  self.showFunctionFit(noOpen=True)
    
    fitFunc = self.fitFunc
    nParams = getFitMethodInfo()[fitFunc-1][1]
    
    numDim = 2
    if self.peakList:
      numDim = self.peakList.dataSource.numDim
   
    headingList, tipTexts  = self.getHeadingList(numDim)
    objectList   = self.peakGroups
    textMatrix   = []
    fitFunctions = self.getFitFuncs()
    proteinConc = self.getBindingSiteConc()
    
    i = 0
    for dataFitting in objectList:
        
      i += 1

      peak = dataFitting.refObject
      peakDims = peak.sortedPeakDims()
      
      datum = []
      datum.append(i)
      for j in range(numDim):
        datum.append(peakDims[j].annotation)
      
      numPeaks = len(dataFitting.objects)
      
      yData = dataFitting.dataY
      if hasattr(dataFitting, 'deltaY'):
        yDelta = dataFitting.deltaY
      else:
        yDelta = None

      if yData:
        dY = yData[-1]
      else:
        dY = None  
      
      datum.append(dY)
      datum.append(yDelta)
      datum.append(dataFitting.fitError)
      datum.append(len(dataFitting.objects))
      datum.append(fitFunctions[dataFitting.fitFunction])
      
      if fitFunc in (7,9):
        if (proteinConc is None) or not dataFitting.parameters:
          datum.append(None)
          datum.append(None)
        else:
          b = dataFitting.parameters[1]
          kd = (b-1) * proteinConc
          datum.append(kd)
          datum.append(dataFitting.parameterErrors[1] * proteinConc)
      
      params = dataFitting.parameters
      for j in range(nParams):
        if j < len(params):
          datum.append(params[j])
        else:
          datum.append(None)

      errors = dataFitting.parameterErrors
      for j in range(nParams):
        if j < len(errors):
          datum.append(errors[j])
        else:
          datum.append(None)

      textMatrix.append( datum )
   
    self.scrolledMatrix.update(textMatrix=textMatrix,
                               objectList=objectList,
                               headingList=headingList,
                               tipTexts=tipTexts)
  
    self.waiting = False
      

  def makeDataList(self):

    if self.peakGroups and self.expSeries:
    
      name, unit = self.getDataListInfo()
      
      peakGroups    = []
      values        = []
      errors        = []
      peakDimGroups = []
      shiftGroups   = []
      for dataFitting in self.peakGroups:
        peaks = dataFitting.objects
        peak = peaks[0]
       
        shiftList = peak.peakList.dataSource.experiment.shiftList or \
                    peak.topObject.findFirstMeasurementList(className='ShiftList')

        peakDims = []
        for peakB in peaks:
          for peakDim in peakB.sortedPeakDims():
            if peakDim.dataDimRef.dataDim in self.dataDims:
              peakDims.append(peakDim)
        
        shifts = []
        for peakDim in peak.sortedPeakDims():
          if peakDim.dataDimRef.dataDim in self.dataDims:
            for contrib in peakDim.peakDimContribs:
              shift = contrib.resonance.findFirstShift(parentList=shiftList)
              if shift:
                shifts.append(shift)
        
        values.append(dataFitting.derivedData)
        errors.append(dataFitting.derivedDataErrors)
        peakGroups.append(dataFitting.objects)
        shiftGroups.append(shifts)
        peakDimGroups.append(peakDims)
        
      dataList = makeDataList(name, unit, values, errors=errors,
                              peaks=peakGroups, peakDims=peakDimGroups, measurements=shiftGroups)

      #print dataList
      
      self.parent.editDerivedData(dataList=dataList)

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

  def getDataDimOptName(self,dataDims):

    return ','.join([self.getDataDimName(dataDim) for dataDim in dataDims])
  
  def getDataDimOpts(self):
  
    opts = []
    
    if self.peakList:
      dataDims = self.peakList.dataSource.sortedDataDims()
      
      opts.append(dataDims)
              
      if len(dataDims) > 2:
        for dataDim in dataDims:
          dataDims2 = list(dataDims)
          dataDims2.remove(dataDim)
          opts.append(dataDims2)
        
      for dataDim in dataDims:
        opts.append([dataDim,])
     
    return opts
    
  def changeDataDim(self, dataDims):
    
    if dataDims != self.dataDims:
      self.dataDims = dataDims
      if self.peakGroups:
        self.calculateAllFits()
  
  
  def updateDataDims(self, *obj):
    
    dataDimOpts = self.getDataDimOpts()
    dataDims    = self.dataDims
    names       = [self.getDataDimOptName(opt) for opt in dataDimOpts]
    index       = -1

    if dataDimOpts:
      if dataDims not in dataDimOpts:
        dataDims = dataDimOpts[0]
      
      index = dataDimOpts.index(dataDims)
    
    if dataDims != self.dataDims:
      self.dataDims = dataDims
      if self.peakGroups:
        self.calculateAllFits()
      
    self.dataDimsPulldown.setup(names, dataDimOpts, index)  
      
  def getExpSeries(self):
  
    return self.nmrProject.sortedNmrExpSeries()

    # Could add option to filter out some series
    
  def getExpSeriesName(self, series):
  
    data = [str(series.serial),','.join([n for n in series.conditionNames])]
    if series.name:
      data.append(series.name)
    return ':'.join(data)
  
  def clearPeakGroups(self):
  
    for peakGroup in self.peakGroups:
      peakGroup.delete()
      
      self.peakGroups = []
  
  def changeExpSeries(self, expSeries):
    
    if self.expSeries is not expSeries:
      self.expSeries = expSeries
      
      sampleConditions = self.getSampleConditions()
      if sampleConditions and (sampleConditions[0].unit in CONC_FACTORS) and (self.unitPulldown.getText() not in CONC_FACTORS):
        self.unitPulldown.set(sampleConditions[0].unit)
      
      self.clearPeakGroups()
      self.updatePeakLists()
      self.updateDataListInfo()
      self.updateAfter()

      
  def updateExpSeries(self, *obj):
    
    
    expSeries     = self.expSeries
    expSeriesList = self.getExpSeries()
    names         = [ self.getExpSeriesName(es) for es in expSeriesList ]
    index         = -1

    if expSeriesList:
      if expSeries not in expSeriesList:
        expSeries = expSeriesList[0]
      
      index = expSeriesList.index(expSeries)
        
    if expSeries != self.expSeries:
      self.expSeries  = expSeries
      
      sampleConditions = self.getSampleConditions()
      if sampleConditions and (sampleConditions[0].unit in CONC_FACTORS):
        self.unitPulldown.set(sampleConditions[0].unit)

      self.clearPeakGroups()
      self.updateDataListInfo()
      self.updateAfter()
    
    self.updatePeakLists()    
    self.expSeriesPulldown.setup(names, expSeriesList, index)

  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
    
    BasePopup.destroy(self)
