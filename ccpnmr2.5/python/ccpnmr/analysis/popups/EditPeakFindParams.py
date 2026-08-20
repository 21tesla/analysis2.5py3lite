
"""
======================COPYRIGHT/LICENSE START==========================

EditPeakFindParams.py: Part of the CcpNmr Analysis program

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

from ccp.api.nmr import Nmr

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError, showWarning
from memops.gui.PulldownList import PulldownList
from memops.gui.RadioButtons import RadioButtons
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.Util import getNmrIsotopes, getIsotopeExclusion, setIsotopeExclusion
from ccpnmr.analysis.core import PeakFindParams
from ccpnmr.analysis.core import PeakBasic
from ccpnmr.analysis.core import ExperimentBasic
from ccpnmr.analysis.core import UnitConverter
 
class RegionCondition:

  def __init__(self, condition, regionMin=None, regionMax=None):

    self.condition = condition
    self.regionMin = regionMin
    self.regionMax = regionMax

class EditPeakFindParamsPopup(BasePopup):
  """
  ** Peak Settings and Non-Interactive Peak Finding **
  
  The purpose of this dialog is to allow the user to select settings for
  finding and integrating peaks, and also to be able to find peaks in an
  arbitrary region that is specified in a table rather than via a spectrum
  window.
  
  ** Find Parameters tab **

  This can be used to specify how peak finding works.

  First of all, you can search for just positive peaks, just negative
  peaks or both, and the default is that it is just positive peaks.
  However, this is further filtered by what the contour levels are.
  If there are no positive contour levels for a given spectrum then
  positive peaks are not found even if this dialog says they can be,
  and similarly if there are no negative contour levels for a given
  spectrum then negative peaks are not found even if this dialog says
  they can be.

  The peak finding algorithm looks for local extrema (maximum for
  positive peaks and minima for negative peaks).  But on a grid there
  are various ways to define what you mean by an extremum.  Suppose
  you are trying to determine if point p is a maximum (similar
  considerations apply for minimum).  You would want the intensity
  at all nearby points to be less than or equal to the intensity at p.
  You can just check points that are just +- one point from p in each
  dimension, or you can also check "diagonal" points.  For
  example, if you are looking at point p = (x, y) in 2D, then the
  former would mean checking the four points (x-1, y), (x+1, y)
  (x, y-1) and (x, y+1), whereas for the latter you would also have
  to check (x-1, y-1), (x-1, y+1), (x+1, y-1) and (x+1, y+1).  In
  N dimensions the "diagonal" method involves checking 3^N-1 points
  whereas the "non-diagonal" method involves checking only 2N points.
  In general the "non-diagonal" method is probably the one to use,
  and it is the default.

  Peaks are only found above (for positive peaks) or below (for negative
  peaks) some threshold.  By default this is determined by the contour level
  for the spectrum.  For positive peaks the threshold is the minimum
  positive contour level, and for negative peaks the threshold is the
  maximum negative contour level.  However these levels can be scaled up
  (or down) using the "Scale relative to contour levels" option (default
  value 1).  For example, if you have drawn the contour levels low to
  show a bit of noise, but do not want the noise picked as peaks, then
  you could select a scale of 2 (or whatever) to increase the threshold.

  The "Exclusion buffer around peaks" is so that in crowded regions you
  do not get too many peaks near one location.  By default the exclusion
  buffer is 1 point in each dimension, but this can be increased to make
  the algorithm find fewer peaks.

  By default the peak finding only looks at the orthogonal region that
  is displayed in the given window where peak finding is taking place.
  Sometimes it looks like a peak should be found because in x, y you
  can see an extremum, but unless it is also an extremum in the orthogonal
  dimensions it is not picked.  You can widen out the points being
  examined in the orthogonal dimensions by using the "Extra thickness in
  orthogonal dims" option, which is specified in points.

  The "Minimum drop factor" is by what factor the intensity needs to drop
  from its extreme value for there to be considered to be a peak.  This
  could help remove sinc wiggle peaks, for example.  The default is that
  the drop factor is 0, which in effect means that there is no condition.

  The "Volume method" is what is used to estimate the volume of peaks that
  are found.  The default is "box sum", which just looks at a fixed size
  box around the peak centre and sums the intensities in that.  The size
  of the box is set in the table in the Spectrum Widths tab.  The
  "truncated box sum" is the same as "box sum" except that the summing
  stops in a given direction when (if) the intensities start increasing.
  The "parabolic" fit fits a quadratic equation in each dimension to the
  intensity at the peak centre and ad +- 1 points and then uses the
  equivalent Gaussian fit to estimate the volume.

  ** Spectrum Widths **

  This can be used to specify minimum linewidths (in Hz) for there to be
  considered a peak to exist in the peak finding algorithm.  It is also
  where the Boxwidth for each dimension in each spectrum is specified.

  ** Diagonal Exclusions **

  This can be used to exclude peaks from being found in regions near
  the diagonal (so in homonuclear experiments).  The exclusion region
  is specified in ppm and is independent of spectrum.

  ** Region Peak Find **

  This can be used to find peaks non-interactively (so not having to
  control shift drag inside a spectrum window).  The region being
  analysed is specified in the table.  There are two types of conditions
  that can be specified, "include" for regions that should be included
  and "exclude" for regions that should be excluded.  The regions are
  specified in ppm.

  The "Whole Region" button will set the selected row in the table to be
  the entire fundamental region of the spectrum.

  The "Add Region" button adds an extra row to the table, and the "Delete
  Region" button removes the selected row.

  The "Adjust Params" button goes to the Find Parameters tab.

  The "Find Peaks!" button does the peak finding.

"""

  def __init__(self, parent, *args, **kw):
 
    self.spectrum = None

    BasePopup.__init__(self, parent=parent, title='Peak : Peak Finding', **kw)

  def body(self, guiFrame):
 
    self.geometry('600x350')

    guiFrame.expandGrid(0,0)
    
    tipTexts = ['',
                '',
                '',
                '']
    options = ['Find Parameters','Spectrum Widths',
               'Diagonal Exclusions','Region Peak Find']
    tabbedFrame = TabbedFrame(guiFrame, options=options, grid=(0,0))
    
    frameA, frameB, frameC, frameD = tabbedFrame.frames
    self.tabbedFrame = tabbedFrame
 
    # Find Params
 
    frameA.expandGrid(2,0)

    row = 0
    label = LabelFrame(frameA, text='Extrema to search for:',
                       grid=(row,0), gridSpan=(1,2))
    label.expandGrid(0,1)
    
    entries = ['positive and negative', 'positive only', 'negative only']
    tipTexts = ['Sets whether peak picking within spectra find intensity maxima, minima or both maxima and minima',]
    self.extrema_buttons = RadioButtons(label, entries=entries,
                                        select_callback=self.apply,
                                        direction='horizontal',
                                        grid=(0,0), tipTexts=tipTexts)
    
    row += 1
    label = LabelFrame(frameA, text='Nearby points to check:',
                       grid=(row,0), gridSpan=(1,2))
    label.expandGrid(None,1)

    entries = ['+-1 in at most one dim', '+-1 allowed in any dim']
    tipTexts = ['Sets how permissive the peak picking in when searching for intensity extrema; by adding extra points to the selected search region',]
    self.adjacent_buttons = RadioButtons(label, entries=entries,
                                         select_callback=self.apply,
                                         direction='horizontal',
                                         grid=(0,0), tipTexts=tipTexts)

    row += 1
    labelFrame = LabelFrame(frameA, text='Other parameters:',
                            grid=(row,0), gridSpan=(1,2))
    labelFrame.expandGrid(5,2)
    
    frow = 0
    label = Label(labelFrame, text='Scale relative to contour levels:',
                   grid=(frow,0), sticky='e')
    tipText = 'Threshold above which peaks are picked, relative to the lowest displayed contour; 1.0 means picking exactly what is visible'
    self.scale_entry = FloatEntry(labelFrame, grid=(frow,1), tipText=tipText,
                                  returnCallback=self.apply, width=10)
    self.scale_entry.bind('<Leave>', self.apply, '+')

    frow += 1
    label = Label(labelFrame, text='Exclusion buffer around peaks (in points):',
                  grid=(frow,0), sticky='e')
    tipText = 'The size of the no-pick region, in data points, around existing picked peaks; eliminates duplicate picking'
    self.buffer_entry = IntEntry(labelFrame, returnCallback=self.apply,
                                 grid=(frow,1), width=10, tipText=tipText)
    self.buffer_entry.bind('<Leave>', self.apply, '+')

    frow += 1
    label = Label(labelFrame, text='Extra thickness in orthogonal dims (in points):',
                  grid=(frow,0), sticky='e')
    tipText = 'Sets whether to consider any additional planes (Z dimension) when calculating peak volume integrals'
    self.thickness_entry = IntEntry(labelFrame, returnCallback=self.apply,
                                    width=10, grid=(frow,1), tipText=tipText)
    self.thickness_entry.bind('<Leave>', self.apply, '+')

    frow += 1
    label = Label(labelFrame, text='Minimum drop factor (0.0-1.0):',
                  grid=(frow,0), sticky='e')
    tipText = ''
    self.drop_entry = FloatEntry(labelFrame, returnCallback=self.apply,
                                 width=10, grid=(frow,1), tipText=tipText)
    self.drop_entry.bind('<Leave>', self.apply, '+')

    frow += 1
    label = Label(labelFrame, text='Volume method:',
                  grid=(frow,0), sticky='e')
    tipText = 'Selects which method to use to calculate peak volume integrals when peaks are picked; box sizes are specified in "Spectrum Widths"'
    self.method_menu = PulldownList(labelFrame, texts=PeakBasic.PEAK_VOLUME_METHODS,
                                    grid=(frow,1), callback=self.apply, tipText=tipText)
    
    # Spectrum widths
    
    frameB.expandGrid(1,1)

    label = Label(frameB, text='Spectrum: ')
    label.grid(row=0, column=0, sticky='e')

    tipText = 'The spectrum which determines the widths being shown'
    self.expt_spectrum = PulldownList(frameB, tipText=tipText,
                                      callback=self.setSpectrumProperties)
    self.expt_spectrum.grid(row=0, column=1, sticky='w')
 
    self.editLinewidthEntry = FloatEntry(self,text='',returnCallback = self.setLinewidth, width=10)
    self.editBoxwidthEntry = FloatEntry(self,text='',returnCallback = self.setBoxwidth, width=10)
    tipTexts = ['The number of the spectrum dimension to which the settings apply',
                'The nuclear isotope measures in the spectrum dimension',
                'The smallest value for the linewidth of a peak for it to be picked',
                'The size of the spectrum region to perform the volume integral over']
    headingList = ['Dimension','Isotope','Minimum Linewidth (Hz)','Boxwidth']
    editSetCallbacks = [None, None, self.setLinewidth, self.setBoxwidth]
    editGetCallbacks = [None, None, self.getLinewidth, self.getBoxwidth]
    editWidgets =      [None, None, self.editLinewidthEntry, self.editBoxwidthEntry]
    self.spectrumMatrix = ScrolledMatrix(frameB, initialRows=6,
                                         editSetCallbacks=editSetCallbacks,
                                         editGetCallbacks=editGetCallbacks,
                                         editWidgets=editWidgets,
                                         headingList=headingList,
                                         callback=self.selectCell,
                                         tipTexts=tipTexts)
    self.spectrumMatrix.grid(row=1, column=0, columnspan=2, sticky='nsew')

    # Diagonal Exclusions

    frameC.expandGrid(0,0)
 
    tipTexts = ['The isotope as measures on the axis of a spectrum window',
                'The distance from the homonuclear diagonal line within which no peak picking can occur']
    self.exclusionEntry = FloatEntry(self,text='',returnCallback = self.setExclusion, width=10)
    headingList = ['Isotope','Diagonal Exclusion (ppm)']
    editSetCallbacks = [None, self.setExclusion]
    editGetCallbacks = [None, self.getExclusion]
    editWidgets =      [None, self.exclusionEntry]
    self.isotopeMatrix = ScrolledMatrix(frameC,
                                        editSetCallbacks=editSetCallbacks,
                                        editGetCallbacks=editGetCallbacks,
                                        editWidgets=editWidgets,
                                        headingList=headingList,
                                        grid=(0,0), tipTexts=tipTexts)

    # Region peak find
    
    self.regionFindPeakList = None
    self.regionCondition = None
    self.regionConditions = []
    self.regionCol = 1

    
    row = 0

    label = Label(frameD, text='Peak List: ', grid=(0,0))
    tipText = 'Selects which peak list to perform region-wide peak picking for'
    self.regionPeakListPulldown = PulldownList(frameD, callback=self.changeRegionPeakList,
                                               grid=(0,1), tipText=tipText)
 
    row += 1
    frameD.expandGrid(row,1)

    self.regionEntry = FloatEntry(self, text='', returnCallback=self.setRegion, width=10)
    self.conditionMenu = PulldownList(self, texts=('include', 'exclude'),
                                      callback=self.setCondition)

    tipTexts = ['Whether to include or exclude the states region from region-wide peak picking',]
    headingList = ['Condition']
    editSetCallbacks = [None]
    editGetCallbacks = [None]
    editWidgets = [self.conditionMenu]
    self.regionFindMatrix = ScrolledMatrix(frameD, headingList=headingList,
                                         callback=self.selectRegionCell,
                                         editWidgets=editWidgets,
                                         editGetCallbacks=editGetCallbacks,
                                         editSetCallbacks=editSetCallbacks,
                                         grid=(row,0), gridSpan=(1,2))

    row += 1
    tipTexts = ['Sets the currently selected region row to cover the whole spectrum',
                'Add a new region row, which may them be set for exclusion or inclusion when peak picking large areas',
                'Remove the selected region specification',
                'Go to the panel for setting the parameters that control how peaks extrema are picked',
                'Using the stated regions and parameters, perform region-wide peak picking']
    texts = [ 'Whole Region', 'Add Region', 'Delete Region',
              'Adjust Params', 'Find Peaks!' ]
    commands = [self.wholeRegion, self.addCondition, self.deleteCondition,
                self.adjustParams, self.regionFindPeaks]
                
    buttons = ButtonList(frameD, texts=texts,
                         commands=commands, grid=(row,0),
                         gridSpan=(1,2), tipTexts=tipTexts)
    buttons.buttons[4].config(bg='#B0FFB0')

    utilButtons = UtilityButtonList(tabbedFrame.sideFrame, grid=(0,0),
                                    helpUrl=self.help_url, sticky='e')
    
    self.dataDim = None
    self.setParamsEntries()
    self.updateSpectrum()
    self.setIsotopeProperties()
    self.updateRegionPeakLists()
    
    self.administerNotifiers(self.registerNotify)

  def administerNotifiers(self, notifyFunc):
    
    # Many more needed here, esp on the AnalysisProject prams
        
    for func in ('__init__', 'delete','setName'):
      notifyFunc(self.updateRegionPeakLists, 'ccp.nmr.Nmr.DataSource', func)
      notifyFunc(self.updateRegionPeakLists, 'ccp.nmr.Nmr.Experiment', func)

    for func in ('__init__', 'delete'):
      notifyFunc(self.updateRegionPeakLists, 'ccp.nmr.Nmr.PeakList', func)

    for clazz in ('Experiment', 'DataSource'):
      for func in ('__init__', 'delete', 'setName'):
        notifyFunc(self.updateSpectrumTable, 'ccp.nmr.Nmr.%s' % clazz, func)

    for func in ('setPeakFindBoxWidth', 'setPeakFindMinLineWidth'):
      notifyFunc(self.updateSpectrumTable, 'ccpnmr.Analysis.AnalysisDataDim', func)

  def destroy(self):
  
    self.administerNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)

  def updateSpectrum(self, spectrum=None):

    if not spectrum:
      spectrum = self.spectrum

    spectra = self.parent.getSpectra()
    if spectra:
      if spectrum not in spectra:
        spectrum = spectra[0]
      index = spectra.index(spectrum)
      names = ['%s:%s' % (x.experiment.name, x.name) for x in spectra]
    else:
      index = 0
      names = []

    self.expt_spectrum.setup(names, spectra, index)

    self.setSpectrumProperties(spectrum)

  def updateNotifier(self, *extra):

    self.updateSpectrum()

  def setLinewidth(self,*event):
  
    value = self.editLinewidthEntry.get()
    if value is not None:
      PeakFindParams.setPeakFindMinLinewidth(self.dataDim, value)
      self.setSpectrumProperties(self.dataDim.dataSource)
  
  def getLinewidth(self,dataDim):
  
    if dataDim:
      self.editLinewidthEntry.set(PeakFindParams.getPeakFindMinLinewidth(self.dataDim))

  def setBoxwidth(self,*event):
  
    value = self.editBoxwidthEntry.get()
    if value is not None:
      PeakFindParams.setPeakFindBoxwidth(self.dataDim, value)
      self.setSpectrumProperties(self.dataDim.dataSource)
  
  def getBoxwidth(self,dataDim):
  
    if dataDim:
      self.editBoxwidthEntry.set(PeakFindParams.getPeakFindBoxwidth(self.dataDim))

  def selectCell(self, object, row, col):

    self.dataDim = object

  def setExclusion(self, *extra):

    isotope = self.isotopeMatrix.currentObject
    if not isotope:
      return

    value = self.exclusionEntry.get()
    if value is not None:
      setIsotopeExclusion(isotope, value)
      self.setIsotopeProperties()

  def getExclusion(self, isotope):

    value = getIsotopeExclusion(isotope)
    self.exclusionEntry.set(value)

  def setParamsEntries(self):

    project = self.project
    params = PeakFindParams.getPeakFindParams(project)

    self.scale_entry.set(params['scale'])
    self.buffer_entry.set(params['buffer'])
    self.thickness_entry.set(params['thickness'])
    self.drop_entry.set(params['drop'])
    volumeMethod = params['volumeMethod']
    if volumeMethod == 'parabolic fit':
      volumeMethod = PeakBasic.PEAK_VOLUME_METHODS[0]
    self.method_menu.set(params['volumeMethod'])

    if (params['nonadjacent']):
      n = 1
    else:
      n = 0
    self.adjacent_buttons.setIndex(n)

    have_high = params['haveHigh']
    have_low = params['haveLow']
    if (have_high and have_low):
      n = 0
    elif (have_high):
      n = 1
    else:
      n = 2
    self.extrema_buttons.setIndex(n)

  def apply(self, *extra):

    params = {}
    params['scale'] = self.scale_entry.get()
    params['buffer'] = self.buffer_entry.get()
    params['thickness'] = self.thickness_entry.get()
    params['drop'] = self.drop_entry.get()
    params['volumeMethod'] = self.method_menu.getText()

    n = self.adjacent_buttons.getIndex()
    if (n == 0):
      nonadjacent = False
    else:
      nonadjacent = True
    params['nonadjacent'] = nonadjacent
    
    n = self.extrema_buttons.getIndex()
    if (n == 0):
      have_high = True
      have_low = True
    elif (n == 1):
      have_high = True
      have_low = False
    elif (n == 2):
      have_high = False
      have_low = True
    params['haveHigh'] = have_high
    params['haveLow']  = have_low

    project = self.project
    try:
      PeakFindParams.setPeakFindParams(project, params)
    except Implementation.ApiError, e:
      showError('Parameter error', e.error_msg, parent=self)

  def reset(self):
 
    self.setParamsEntries()

    spectrum = self.spectrum
    if spectrum:
      self.dataDim = None
      self.spectrumMatrix.update(objectList=None, textMatrix=None)
 
  def setSpectrumProperties(self, spectrum):
 
    if spectrum is self.spectrum:
      return

    self.spectrum = spectrum
    self.updateSpectrumTable()

  def updateSpectrumTable(self, *extra):

    spectrum = self.spectrum

    textMatrix = []
    dataDims   = []
    if spectrum:
      dataDims = spectrum.sortedDataDims()
 
    objectList = []
    for i in range(len(dataDims)):
      dataDim = dataDims[i]
      if isinstance(dataDim, Nmr.FreqDataDim):
        objectList.append(dataDim)
        textMatrix.append( [dataDim.dim, ExperimentBasic.getPrimaryDataDimRef(dataDim).expDimRef.isotopeCodes[0],
                          PeakFindParams.getPeakFindMinLinewidth(dataDim), PeakFindParams.getPeakFindBoxwidth(dataDim)] )
    
    self.spectrumMatrix.update(objectList=objectList, textMatrix=textMatrix)    

  def setIsotopeProperties(self):

    project  = self.project
    isotopes = getNmrIsotopes(project)

    textMatrix = []
    for isotope in isotopes:
      name = '%d%s' % (isotope.massNumber, isotope.chemElement.symbol)
      exclusion = getIsotopeExclusion(isotope)
      textMatrix.append([name, exclusion])
    
    self.isotopeMatrix.update(objectList=isotopes, textMatrix=textMatrix)    
  # Region peak find functions

  def updateRegionPeakLists(self, *extra):
  
    index = 0
    names = []
    peakLists = []
    peakList = self.regionFindPeakList
    
    for experiment in self.nmrProject.sortedExperiments():
      eName = experiment.name
      
      for spectrum in experiment.sortedDataSources():
        sName = spectrum.name
        
        if spectrum.dataType == 'processed':
          for pl in spectrum.sortedPeakLists():
            names.append('%s:%s %d' % (eName, sName, pl.serial))
            peakLists.append(pl)
    
    if peakLists:
      if peakList not in peakLists:
        peakList = peakLists[0]
      
      index = peakLists.index(peakList) 
    
    if peakList is not self.regionFindPeakList:
      self.changeRegionPeakList(peakList)
      
    self.regionPeakListPulldown.setup(names, peakLists, index)

  def wholeRegion(self):
 
    if not self.regionFindPeakList:
      return

    if not self.regionCondition:
      msg = 'No condition selected'
      showWarning('Failure', msg, parent=self)
      return

    spectrum = self.regionFindPeakList.dataSource
    
    regionMin = []
    regionMax = []
    
    for dataDim in spectrum.sortedDataDims():
      (rMin, rMax) = self.getWholeRegion(dataDim)
      regionMin.append(rMin)
      regionMax.append(rMax)

    self.regionCondition.regionMin = regionMin
    self.regionCondition.regionMax = regionMax
    
    self.updateRegionFindTable()
 
  def getDimMin(self, dataDim):

    if dataDim.className == 'SampledDataDim':
      r = 1.0
    else:
      converter = UnitConverter.pnt2ppm
      dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
      r = converter(float(dataDim.numPoints), dataDimRef)

    return r

  def getDimMax(self, dataDim):

    if dataDim.className == 'SampledDataDim':
      r = float(dataDim.numPoints)
    else:
      converter = UnitConverter.pnt2ppm
      dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
      r = converter(1.0, dataDimRef)

    return r

  def getWholeRegion(self, dataDim):

    rMin = self.getDimMin(dataDim)
    rMax = self.getDimMax(dataDim)

    return (rMin, rMax)

  def addCondition(self):

    if not self.regionFindPeakList:
      return

    spectrum = self.regionFindPeakList.dataSource

    ndim = spectrum.numDim
    regionMin = 2*ndim*[None]
    regionMax = 2*ndim*[None]
    
    condition = RegionCondition('exclude', regionMin, regionMax)
    self.regionConditions.append(condition)

    self.updateRegionFindTable()
    self.regionFindMatrix.selectObject(condition)

  def deleteCondition(self):

    if self.regionCondition:
      self.regionConditions.remove(self.regionCondition)
      self.regionCondition = None
      
      self.updateRegionFindTable()

      objs = self.regionFindMatrix.objectList
      if objs:
        self.regionFindMatrix.selectObject(objs[-1])

  def adjustParams(self):

    self.tabbedFrame.select(0)


  def checkRegion(self, dataDims, rMin, rMax):

    for i, dataDim in enumerate(dataDims):

      pMax = rMax[i]
      pMin = rMin[i]
      
      absMin = self.getDimMin(dataDim)
      absMax = self.getDimMax(dataDim)

      if pMin is None:
        pMin = self.getDimMin(dataDim)

      if pMax is None:
        pMax = self.getDimMax(dataDim)

      pMax = min(absMax, pMax)
      pMin = max(absMin, pMin)
      
      if pMin > pMax:
        pMin, pMax = pMax, pMin

      rMax[i] = pMax
      rMin[i] = pMin

    return (rMin, rMax)

  def convertToPoints(self, dataDims, region):

    pointsRegion = []
    for i, dataDim in enumerate(dataDims):
      rMin, rMax = region[i]
    
      if dataDim.className != 'SampledDataDim':
        converter  = UnitConverter.ppm2pnt
        dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)
        rMin = converter(rMin, dataDimRef)
        rMax = converter(rMax, dataDimRef)
        rMin = max(1, rMin)
        rMax = min(dataDim.numPoints, rMax)

      rMin = rMin - 1
      rMax = rMax - 1

      if rMin > rMax:
        rMax, rMin = rMin, rMax

      pointsRegion.append((rMin, rMax))

    return pointsRegion

  def regionFindPeaks(self):

    peakList = self.regionFindPeakList
    if not peakList:
      return

    spectrum = peakList.dataSource
    dataDims = spectrum.sortedDataDims()
    nDims = len(dataDims)
    includedRegions = []
    excludedRegions = []

    for regionCondition in self.regionConditions:
      
      instruction = regionCondition.condition
      rMin = regionCondition.regionMin
      rMax = regionCondition.regionMax
      rMin, rMax = self.checkRegion(dataDims, rMin, rMax)
      
      region = []
      for i in range(nDims):
        r = (rMin[i], rMax[i])
        region.append(r)

      if instruction == 'include':
        includedRegions.append(region)

      else:
        # this region supposed to be in points
        region = self.convertToPoints(dataDims, region)
        excludedRegions.append(region)

    for region in includedRegions:
      PeakBasic.findPeaks(peakList, region, excludedRegions=excludedRegions)

  def selectRegionCell(self, obj, row, col):

    self.regionCondition = obj
    self.regionCol = col

  def getCondition(self, obj):
    
    self.conditionMenu.set(obj.condition)

  def setCondition(self, obj):
    
    self.regionCondition.condition = self.conditionMenu.getText()
    self.updateRegionFindTable()

  def getRegionMin(self, obj, dim):
    
    self.regionEntry.set(self.regionCondition.regionMin[dim])

  def setRegionMin(self, obj, dim):
    
    rMin = self.regionCondition.regionMin
    rMax = self.regionCondition.regionMax
    rMin[dim] = self.regionEntry.get()

    peakList = self.regionFindPeakList
    dataDims = peakList.dataSource.sortedDataDims()
    rMin, rMax = self.checkRegion(dataDims, rMin, rMax)
    self.regionCondition.regionMin = rMin
    self.regionCondition.regionMax = rMax
    self.updateRegionFindTable()

  def getRegionMax(self, obj, dim):
    
    self.regionEntry.set(self.regionCondition.regionMax[dim])

  def setRegionMax(self, obj, dim):
 
    rMin = self.regionCondition.regionMin
    rMax = self.regionCondition.regionMax
    rMax[dim] = self.regionEntry.get()

    peakList = self.regionFindPeakList
    dataDims = peakList.dataSource.sortedDataDims()
    rMin, rMax = self.checkRegion(dataDims, rMin, rMax)
    self.regionCondition.regionMin = rMin
    self.regionCondition.regionMax = rMax
    self.updateRegionFindTable()

  def setRegion(self, *event):

    if self.regionCondition and self.regionCol is not None:
      col = self.regionCol - 1 # -1 because of condition
      dim = col / 2
      if col % 2: # max
        self.setRegionMax(self.regionCondition, dim)
      else: # min
        self.setRegionMin(self.regionCondition, dim)

  def changeRegionPeakList(self, peakList):
 
    if peakList is self.regionFindPeakList:
      return

    if peakList and self.regionFindPeakList:
      spectrum1 = peakList.dataSource
      spectrum2 = self.regionFindPeakList.dataSource
      self.regionFindPeakList = peakList

      if spectrum1 is spectrum2:
        return

      if spectrum1.numDim == spectrum2.numDim:
        for n in range(spectrum1.numDim):
          dim = n + 1
          dataDim1 = spectrum1.findFirstDataDim(dim=dim)
          dataDim2 = spectrum2.findFirstDataDim(dim=dim)
          if dataDim1.className != dataDim2.className:
            break
          
          if (dataDim1.className == 'FreqDataDim') and (dataDim2.className == 'FreqDataDim'):
            if ExperimentBasic.getPrimaryDataDimRef(dataDim1).expDimRef.isotopeCodes != ExperimentBasic.getPrimaryDataDimRef(dataDim2).expDimRef.isotopeCodes:
              break
        else:
          return # just use what is there already as sensible default
    else:
      self.regionFindPeakList = peakList

    if peakList:
      spectrum = peakList.dataSource
      ndim = spectrum.numDim
      dataDims = spectrum.sortedDataDims()
    else:
      spectrum = None
      ndim = 0
      dataDims = []
 
    tipTexts = ['Whether to include or exclude the states region from region-wide peak picking',]
    headingList = ['Condition']
    textRow = ['include']
    regionMin = []
    regionMax = []
    editWidgets = [self.conditionMenu] + 2*ndim*[self.regionEntry]
    editGetCallbacks = [self.getCondition]
    editSetCallbacks = [self.setCondition]
    for dataDim in dataDims:
      dim = dataDim.dim
      headingList.extend(['Dim %d Min' % dim, 'Dim %d Max' % dim])
      tipTexts.append('Lower value bound of peak picking inclusion/exclusion region for spectrum dimension %s' % dim)
      tipTexts.append('Upper value bound of peak picking inclusion/exclusion region for spectrum dimension %s' % dim)
      (rMin, rMax) = self.getWholeRegion(dataDim)
      textRow.append(rMin)
      textRow.append(rMax)
      regionMin.append(rMin)
      regionMax.append(rMax)

      i = dim-1
      editGetCallbacks.append(lambda row, i=i: self.getRegionMin(row, i))
      editGetCallbacks.append(lambda row, i=i: self.getRegionMax(row, i))
      editSetCallbacks.append(lambda row, i=i: self.setRegionMin(row, i))
      editSetCallbacks.append(lambda row, i=i: self.setRegionMax(row, i))

    condition = RegionCondition('include', regionMin, regionMax)
    objectList = [ condition ]
    textMatrix = [ textRow ]
    self.regionConditions = [condition]
    self.regionFindMatrix.update(objectList=objectList,
                                 textMatrix=textMatrix,
                                 headingList=headingList,
                                 tipTexts=tipTexts,
                                 editSetCallbacks=editSetCallbacks,
                                 editGetCallbacks=editGetCallbacks,
                                 editWidgets=editWidgets)

  def updateRegionFindTable(self):

    peakList = self.regionFindPeakList

    if peakList:
      ndim = peakList.dataSource.numDim
    else:
      ndim = 0
 
    textMatrix = []
    objectList = self.regionConditions
    objectList.sort()
    for regionCondition in objectList:
      instruction = regionCondition.condition
      rMin = regionCondition.regionMin
      rMax = regionCondition.regionMax

      textRow = [instruction, ]
      
      for i in range(ndim):
        textRow.append(rMin[i])
        textRow.append(rMax[i])
      
      textMatrix.append(textRow)

    self.regionFindMatrix.update(objectList=objectList, textMatrix=textMatrix)
