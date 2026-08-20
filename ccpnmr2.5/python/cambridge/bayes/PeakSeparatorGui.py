#!/usr/bin/env python
# encoding: utf-8
"""
PeakSeparatorGui.py

Created by Daniel O'Donovan on 2008-09-22.
Copyright (c) 2008 University of Cambridge. All rights reserved.

"""

import sys, os

from memops.gui.ButtonList          import ButtonList, UtilityButtonList
from memops.gui.FloatEntry          import FloatEntry
from memops.gui.Frame               import Frame
from memops.gui.IntEntry            import IntEntry
from memops.gui.Label               import Label
from memops.gui.LabelDivider        import LabelDivider
from memops.gui.LabelFrame          import LabelFrame
from memops.gui.Menu                import Menu
from memops.gui.MessageReporter     import showError
from memops.gui.PulldownList        import PulldownList
from memops.gui.RadioButtons        import RadioButtons
from memops.gui.ScrolledMatrix      import ScrolledMatrix
from memops.gui.TabbedFrame         import TabbedFrame

from ccp.api.nmr.Nmr                        import FreqDataDim

from ccpnmr.analysis.core.ExperimentBasic   import getPrimaryDataDimRef
from ccpnmr.analysis.core.WindowBasic       import getSpectrumViews
from ccpnmr.analysis.core.UnitConverter     import pnt2hz, hz2pnt
from ccpnmr.analysis.popups.BasePopup       import BasePopup

from PeakSeparatorParams            import PeakSeparatorParams
from PeakSeparatorPeakList          import getPeakListParams
from PeakSeparatorRegion            import getRegionParams
from PeakSeparator                  import SeparatePeakRoutine, SeparatePeaksInPeakList

class PeakSeparatorGui(BasePopup):
  """
  **Separate Merged Peaks Using Peak Models**

  The Peak Separator code uses a Markov Chain Monte Carlo search which, using
  idealised peak shapes, attempts to deconvolve overlapped peak regions into 
  their separate constituent peaks.
  
  This routine is also suitable for accurately fitting model shapes to single
  peaks in order to calculate precise intensities.
  
  **Options Peak Separator Parameters**
  *Min. Number of peaks* is by default set to one, it is not possible to set 
  this to a value less than one.
  *Max. Number of peaks* is by default set to one, increasing this value allows
  the search routine to fit more models. The best fit may be found with fewer than
  the maximum number models. Higher numbers slow the routine, and setting this
  value to 0 allows the routine to (effectively) fit unlimited peaks.
  *Only pick positive peaks*. If you are not interested in negative peaks, removing
  the possibility of fitting negative peaks can reduce search time.
  *Peak Model* fits the spectra with either a Gaussian peak model or a Lorentzian
  peak model.

  **Options Region**
  *Peak List* choose which peak list newly picked peaks should be added to. Peaks
  picked using this method will have their details appended with 'PeakSepartor' 
  so you know where they came from.
  *Region Table* shows which area of the current spectrum is about to be searched.
  *Add Region*. Once an area of spectra has been highlighted clicking this button
  will pass it's details on to the Peak Separator.
  *Reset All* will reset all search parameters.
  *Separate Peaks* will run the Peak Separator code with your current settings. This
  may take a few minutes to run, depending on the size of the spectral region being
  searched, the number of peaks being fitted and the speed of your machine. Please
  wait while this completes.
  
  After a successful Peak Separation run, the found peaks will be added to the 
  selected peak list. These peaks intensties (volume) have been found using the
  peak model selected.

  **Advanced Settings Tab**
  *Rate* affects the speed of the Markov Chain Monte Carlo routine. A smaller value
  results in longer execution, but possibly higher quality results. The default 
  setting is deemed sensible for the majority of runs.
  *Line Width* offers a finer degree of control over maximum and minimum peak widths
  for each dimension. The default values are *very* stupid and could do with 
  re-checking for each experiment.
  *Re-Pick Entire Peak List* if you would like to use the Peak Separator to repick
  *every* peak in your peak list, try this option - but note that this may take
  a very long time!

  """

  def __init__(self, parent, programName='Peak Separator', **kw ):

    self.parent         = parent
    self.programName    = programName
    self.versionInfo    = 'Version 0.2'
    self.help_url       = 'http://www.ccpn.ac.uk/'

    self.window         = None
    self.waiting        = False
    self.rootWindow     = None

    # just used for display - PeakSeparator will not see this
    self._minSigmaHz     = None
    self._maxSigmaHz     = None

    self.customSigma     = False
    self.rePickPeakList  = False

    self._sampleStartPpm = None
    self._sampleEndPpm   = None

    try:
      self.project         = parent.project
    except:
      pass

    self.params = PeakSeparatorParams()

    BasePopup.__init__(self, parent=parent, title=programName, 
                        location='+100+100', **kw)

    if not self.analysisProject:
      print '&&& init: No analysis project found ...'
    try:
      if parent.argumentServer:
        self.argServer = parent.argumentServer
      else:
        print '&&& init: No argument server found...'
    except:
      print '&&& init: Test'

  ###########################################################################

  def body(self, guiFrame):

    self.geometry('450x500')

    guiFrame.grid_rowconfigure(0,weight=1)
    guiFrame.grid_columnconfigure(0,weight=1)

    options = ['Peak Separator', 'Advanced Settings']

    tabbedFrame = TabbedFrame(guiFrame, options=options)
    tabbedFrame.grid(row=0, column=0, sticky='nsew')

    buttons = UtilityButtonList(tabbedFrame.sideFrame,
                                helpUrl=self.help_url)
    buttons.grid(row=0, column=0, sticky='e')

    self.tabbedFrame = tabbedFrame
    frameA, frameB = tabbedFrame.frames

    #
    # FrameA : Main Settings
    #

    frameA.grid_columnconfigure(1, weight=1)
    row = 0 # Label row
    
    row += 1
    div = LabelDivider(frameA, text='Peak Separator Parameters')
    div.grid(row=row, column=0, columnspan=2, sticky='ew')

    row += 1
    label = Label(frameA, text='Min. number of peaks:')
    label.grid(row=row, column=0, sticky='w')
    self.minPeaksEntry = IntEntry(frameA, returnCallback=self.applyChange, width=10, \
          tipText='Minimum number of peaks to find (must be > 0)')
    self.minPeaksEntry.grid(row=row, column=1, sticky='n')
    self.minPeaksEntry.bind('<Leave>', self.applyChange, '+')

    row += 1
    label = Label(frameA, text='Max. number of peaks:')
    label.grid(row=row, column=0, sticky='w')
    self.maxPeaksEntry = IntEntry(frameA, returnCallback=self.applyChange, width=10, \
          tipText='Maximum number of peaks to find (0 is unlimited - not recommended)')
    self.maxPeaksEntry.grid(row=row, column=1, sticky='n')
    self.maxPeaksEntry.bind('<Leave>', self.applyChange, '+')

    row += 1
    label = Label(frameA, text='Only pick positive peaks:')
    label.grid(row=row, column=0, sticky='w')
    entries = ['False', 'True']
    self.posPeaksButtons = RadioButtons(frameA, entries=entries,
                                    select_callback=self.applyChange,
                                    direction='horizontal', 
                                    tipTexts=['Search for both positive and negative intensity peaks',
                                              'Limit search to only positive peaks'])
    self.posPeaksButtons.grid(row=row, column=1, sticky='n')

    row += 1
    label = Label(frameA, text='Peak Model:')
    label.grid(row=row, column=0, sticky='w')
    ### G/L Mixture works, but volume calculation involves Gamma function
    # entries = ['Gaussian', 'Lorentzian', 'G/L Mixture']
    entries = ['Gaussian', 'Lorentzian']
    self.shapeButtons = RadioButtons(frameA, entries=entries,
                                        select_callback=self.applyChange,
                                        direction='horizontal',
                                        tipTexts=['Choose a Gaussian model peak shape to fit to peaks',
                                                  'Choose a Lorentzian model peak shape to fit to peaks'])
    self.shapeButtons.grid(row=row, column=1, sticky='n')

    row += 1
    div = LabelDivider(frameA, text='Region', tipText='Region that search will limit itself to')
    div.grid(row=row, column=0, columnspan=2, sticky='ew')

    row += 1
    label = Label(frameA, text='Peak List:')
    label.grid(row=row, column=0, sticky='nw')
    self.peakListPulldown = PulldownList(frameA, callback=self.setManuallyPickPeakList, 
                                        tipText='Select which peak list new peaks are to be added to')
    self.peakListPulldown.grid(row=row, column=1, sticky='nw')

    # tricky scrolled matrix
    row += 1
    self.regionTable = None
    frameA.grid_rowconfigure(row, weight=1)
    headings = ('dim.', 'start (ppm)', 'end (ppm)', 'actual size')

    self.editDimEntry   = IntEntry(self,   returnCallback=self.applyChange, width=5, tipText='Dimension number')
    self.editStartEntry = FloatEntry(self, returnCallback=self.applyChange, width=5, tipText='Search area lower bound')
    self.editEndEntry   = FloatEntry(self, returnCallback=self.applyChange, width=5, tipText='Search area upper bound')

    editWidgets = [self.editDimEntry, self.editStartEntry, self.editEndEntry, None]

    editGetCallbacks = [None, None, None, None]
    editSetCallbacks = [None, None, None, None]

    self.regionTable = ScrolledMatrix(frameA,   headingList=headings,
                                                multiSelect=False,
                                                editWidgets=editWidgets,
                                                editGetCallbacks=editGetCallbacks,
                                                editSetCallbacks=editSetCallbacks,
                                                initialRows=5)

    self.regionTable.grid(row=row, column=0, columnspan=2, sticky='nsew')

    # Run Button
    row += 1
    texts = ['Add Region' ]
    commands = [ self.updateFromRegion ]
    self.addResetButtons = ButtonList(frameA, texts=texts, commands=commands, 
          tipTexts=['Add selected specrtral region'])
    self.addResetButtons.grid(row=row, column=0, columnspan=2, sticky='ew')

    row += 1
    texts = [ 'Separate Peaks' ]
    commands = [ self.runPeakSeparator ]
    self.runButton = ButtonList(frameA, texts=texts, commands=commands, expands=True, tipTexts=['Run peak search now'])
    self.runButton.grid(row=row, column=0, columnspan=2, sticky='nsew')

    #
    # FrameB : Further Settings
    #

    frameB.grid_columnconfigure(0, weight=1)

    row = 0

    div = LabelDivider(frameB, text='Rate:')
    div.grid(row=row,column=0,columnspan=2,sticky='ew')
    row += 1

    label = Label(frameB, text='Rate of MCMC step size change')
    label.grid(row=row, column=0, columnspan=1, sticky='w')

    self.rateEntry = FloatEntry(frameB, returnCallback=self.applyChange, width=10, \
          tipText='Rate effects speed of run, smaller values take longer but may produce better results')
    self.rateEntry.grid(row=row, column=1, sticky='n')
    self.rateEntry.bind('<Leave>', self.applyChange, '+')
    self.rateEntry.set(self.params.rate)

    # tricky scrolled matrix for line width
    row += 2
    div = LabelDivider(frameB, text='Line Width (Hz):')
    div.grid(row=row,column=0,columnspan=2,sticky='ew')

    row += 1
    label = Label(frameB, text="Descr.")
    label.grid(row=row, rowspan=2, column=0, sticky='w')

    row += 1
    self.lineWidthTable = None
    frameB.grid_rowconfigure(row, weight=1)
    lineWidthHeadings = ('dim.', 'min. σ (Hz)', 'max. σ (Hz)')

    self.editMinSigmaEntry  = FloatEntry(self, returnCallback=self.applyChange, width=5, tipText='Minimum line width (Hz)')
    self.editMaxSigmaEntry  = FloatEntry(self, returnCallback=self.applyChange, width=5, tipText='Maximum line width (Hz)')

    # self.editDimEntry is also from regionTable
    initialWidthRows = 4

    editLineWidthWidgets      = [None, self.editMinSigmaEntry,  self.editMaxSigmaEntry]
    editLineWidthGetCallbacks = [None, self.getSigmaMin,        self.getSigmaMax]
    editLineWidthSetCallbacks = [None, self.setSigmaMin,        self.setSigmaMax]

    self.lineWidthTable = ScrolledMatrix(frameB,  headingList=lineWidthHeadings,
                                                  multiSelect=False,
                                                  editWidgets=editLineWidthWidgets,
                                                  editGetCallbacks=editLineWidthGetCallbacks,
                                                  editSetCallbacks=editLineWidthSetCallbacks,
                                                  initialRows=initialWidthRows)

    self.lineWidthTable.grid(row=row, column=0, columnspan=2, sticky='nsew')

    # option to 'repick' exisiting peak list
    row += initialWidthRows
    div = LabelDivider(frameB, text='(optional - repick entire peak list)')
    div.grid(row=row,column=0,columnspan=2,sticky='ew')
    row += 1

    self.repickListPulldown = PulldownList(frameB, callback=self.setRePickPeakList, 
                                                          tipText='Select which peak list to repick (new peaks will be put into a new peak list)')
    self.repickListPulldown.grid(row=row, column=0, sticky='nw')

    texts = [ 'Repick Peak List' ]
    commands = [ self.runRepickPeaks ]
    self.runButton = ButtonList(frameB, texts=texts, commands=commands, expands=True, tipTexts=['Repick selected peak list into a new peak list.'])
    self.runButton.grid(row=row, column=1, columnspan=1, sticky='nsew')

    row += 1
    div = LabelDivider(frameB)
    row += 1
    texts = [ 'Separate Peaks' ]
    commands = [ self.runPeakSeparator ]
    self.runButton = ButtonList(frameB, texts=texts, commands=commands, expands=True, tipTexts=['Run peak search now'])
    self.runButton.grid(row=row, column=0, columnspan=2, sticky='nsew')

    self.setWidgetEntries()

    self.administerNotifiers( self.registerNotify )

  def administerNotifiers(self, notifyFunc):

    for func in ('__init__', 'delete'):
      notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.PeakList', func)

    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.Experiment', 'setName')
    notifyFunc(self.updateAfter, 'ccp.nmr.Nmr.DataSource', 'setName')

  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)

  ###########################################################################
  # update parameters from PS Region

  def updateFromRegion(self):

    if not self.params.peakList:
      print '&&& update from region: Need a peak list'
      return

    if (self.argServer.parent.currentRegion) == None:
      showError('No Region','Please select a peak region to be separated')
      return

    self.rePickPeakList = False

    getRegionParams(self.params, argServer=self.argServer)

    if not self.customSigma: self.initSigmaParams()

    self.setWidgetEntries()

  ###########################################################################
  # update parameters from PS PeakList
  
  def updateFromPeakList(self):

    if not self.params.peakList:
      print '&&& update from peakList: Need a peak list'
      return

    getPeakListParams(self.params)

    if not self.customSigma: self.initSigmaParams()

    self.setWidgetEntries()


  ###########################################################################
  # Run the C library!

  def runPeakSeparator(self):
    """ run the peak separator """

    # hack for Macs - focus isn't always lost on mouse move 
    # so bind event not always called. Shouldn't affect other OS.
    self.applyChange()

    if not self.params.peakList:
      print '&&& Peak list not yet set'
    else:
      # SeparatePeakRoutine(self.params, self.params.peakList, routine='pymc' )
      SeparatePeakRoutine(self.params, self.params.peakList, routine='bayesys' )

  def runRepickPeaks( self ):
    """ Run the Peak Separator on entire chosen peak list """
    # hack for Macs - focus isn't always lost on mouse move 
    # so bind event not always called. Shouldn't affect other OS.
    self.applyChange()

    if not self.params.peakList:
      print '&&& Peak list not yet set'
    else:
      SeparatePeaksInPeakList( self.params )

  ###########################################################################

  def setWidgetEntries(self):

    ### Page One widgets
    self.minPeaksEntry.set(self.params.minAtoms)
    self.maxPeaksEntry.set(self.params.maxAtoms)

    if self.params.positivePeaks == 1:
      self.posPeaksButtons.set('True') # only pick pos peaks
    else:
      self.posPeaksButtons.set('False')

    # do something fancy if different shapes for each dim!
    n = self.params.peakShape - 3       # shape is only 3, 4, (5)
    self.shapeButtons.setIndex(n)

    if self.project is not None:
      self.updatePeakListList()
    self.updateSpectrumWindow()

    if self.params.sampleStart and self.params.peakList:

      if not self.rePickPeakList:
        objectList = []
        textMatrix = []

        if len(self.params.samplePpmStart) != self.params.Ndim: return

        for i in range( self.params.Ndim ):
          dim_entry = []
          dim_entry.append( '%2d' % (i + 1) )
          dim_entry.append( '%7.3f' % self.params.samplePpmStart[i] )
          dim_entry.append( '%7.3f' % self.params.samplePpmEnd[i] )
          dim_entry.append( '%3d'   % self.params.sampleSize[i] )
          textMatrix.append( dim_entry )

        self.regionTable.update( textMatrix=textMatrix, objectList=objectList )

    ### Page Two widgets
    self.rateEntry.set(self.params.rate)

    if self.params.peakList and self.params.Ndim:

      textMatrix = []
      objectList = []

      for i in range( self.params.Ndim ):
        if self.params.isFreqDim[i]:
          dim_entry = []
          objectList.append( i )
          dim_entry.append( '%2d' % (i + 1) )
          dim_entry.append( '%7.3f' % self._minSigmaHz[i] )
          dim_entry.append( '%7.3f' % self._maxSigmaHz[i] )
          textMatrix.append( dim_entry )

      self.lineWidthTable.update( textMatrix=textMatrix, objectList=objectList )

  def applyChange(self, *event):
    """ Upon change, add settings to params """

    # Page One apply changes
    self.params.minAtoms = self.minPeaksEntry.get()
    self.params.maxAtoms = self.maxPeaksEntry.get()

    if self.posPeaksButtons.get() == 'True': # asked only pick pos peaks
      self.params.positivePeaks = 1
    else:
      self.params.positivePeaks = 0

    # do something fancy if different shapes for each dim!
    n = self.shapeButtons.getIndex() # shape is only 3, 4, (5)
    self.params.peakShape = n + 3

    # Page Two apply changes
    self.params.rate = float( self.rateEntry.get() )

    self.updateSigmaParams()

  ###########################################################################
  # Peak list functions provide PeakSeparator some inherited params

  def getPeakListList(self):
    """ given a spectrum, get list of peak lists """
    project = self.project

    peakLists = []
    for experiment in self.nmrProject.experiments:
      for spectrum in experiment.dataSources:
        for peakList in spectrum.peakLists:
          peakLists.append(['%s:%s:%d' % (experiment.name,spectrum.name,peakList.serial), peakList])
    peakLists.sort()
    return peakLists

  def updatePeakListList(self):
    """ set the peaklist list in the pulldown menu """
    peakListData = self.getPeakListList()

    index    = -1
    names    = []
    peakList = self.params.peakList

    if peakListData:
      names     = [x[0] for x in peakListData]
      peakLists = [x[1] for x in peakListData]

      if peakList not in peakLists:
        peakList = peakLists[0]

      index = peakLists.index(peakList)

    else:
      peakList = None
      peakLists = []

    if peakList is not self.params.peakList:
      self.params.peakList = peakList

    self.peakListPulldown.setup(names, peakLists, index)
    self.repickListPulldown.setup(names, peakLists, index)

  def setRePickPeakList(self, peakList):
    """ Set the peak list to be repicked (and hit a Flag) """
    self.rePickPeakList = True
    self.setPeakList(peakList)

  def setManuallyPickPeakList(self, peakList):
    """ Set the peak list to add new peaks to (and hit a Flag) """
    self.rePickPeakList = False
    self.setPeakList(peakList)

  def setPeakList(self, peakList):
    """ Sets the Peak List """
    if peakList is not self.params.peakList:
      self.params.peakList = peakList
      # # interrogate the peak list and get all the usefull parameters out
      self.updateFromPeakList()
      self.updateSpectrumWindow()
      self.setWidgetEntries()

  ###########################################################################
  # TBD I suspect this is for matching region with peak list, but may be obsolete now

  def getSpectrumWindowList(self):
    """ get list of windows which spectrum could be in """
    windows = {}
    if self.params.peakList:
      views = getSpectrumViews(self.params.peakList.dataSource)
      for view in views:
        windows[view.spectrumWindowPane.spectrumWindow] = None

    return [[w.name, w] for w in windows.keys()]

  def updateSpectrumWindow(self):
    """ update the spectrum window """
    windowData = self.getSpectrumWindowList()

    index  = -1
    names  = []
    window = self.rootWindow

    if windowData:
      names   = [x[0] for x in windowData]
      windows = [x[1] for x in windowData]

      if window not in windows:
        window = windows[0]

      index = windows.index(window)

    else:
      window = None
      windows = []

    if window is not self.rootWindow:
      self.rootWindow = window

  ###########################################################################
  # get and set sigma stuff
  def setSigmaMin( self, dim ):

    value = self.editMinSigmaEntry.get()
    self._minSigmaHz[dim] = value

    # dont go and re-write users settings
    self.customSigma = True

    # make sure changes are in params object
    self.updateSigmaParams( dim )
    self.setWidgetEntries()

  def getSigmaMin( self, dim ):

    if dim is not None:
      self.editMinSigmaEntry.set( self._minSigmaHz[dim] )

  def setSigmaMax( self, dim ):

    value = self.editMaxSigmaEntry.get()
    self._maxSigmaHz[dim] = value

    # dont go and re-write users settings
    self.customSigma = True

    # make sure changes are in params object
    self.updateSigmaParams( dim )
    self.setWidgetEntries()

  def getSigmaMax( self, dim ):

    if dim is not None:
      self.editMaxSigmaEntry.set( self._maxSigmaHz[dim] )

  def updateSigmaParams( self, dim=None ):
    """ updateSigmaParams Just updates the parameters (params obj) for sigma values. 
        If dim is None, do this for each dim
    """

    dataDimRefs = self.params.dataDimRefs

    if not dataDimRefs: return

    if not self.params.minSigma or len(self.params.minSigma) != self.params.Ndim:
      self.params.minSigma = [0.] * self.params.Ndim

    if not self.params.maxSigma or len(self.params.maxSigma) != self.params.Ndim:
      self.params.maxSigma = [0.] * self.params.Ndim

    def updateSigmaParam( dim, dataDimRefs ):
      """ Convert and update sigma for dim """

      if self.params.isFreqDim[dim]:
        # note factor of two!
        self.params.minSigma[dim] = self.rHz2pnt( self._minSigmaHz[dim], dataDimRefs[dim] ) / 2.
        self.params.maxSigma[dim] = self.rHz2pnt( self._maxSigmaHz[dim], dataDimRefs[dim] ) / 2.
      else:
        self.params.minSigma[dim] = 1.0
        self.params.maxSigma[dim] = 1.0

    if dim:
      updateSigmaParam( dim, dataDimRefs )
    else:
      for dim in range( self.params.Ndim ):
        updateSigmaParam( dim, dataDimRefs )

  # utility functions for sigma values
  def pnt2rHz( self, point, dataDimRef ):
    """ Point to relative Hz frequency relative to frequency at Zeroeth point
        Necessary when (for example) looking for width of peak in Hz
    """
    assert point, dataDimRef

    sigmaBase   = pnt2hz(     0, dataDimRef )
    sigmaHz     = pnt2hz( point, dataDimRef )

    return abs( sigmaHz - sigmaBase )

  def rHz2pnt( self, freq, dataDimRef ):
    """ Relative Hz to point frequency relative to frequency at Zeroeth point
        Necessary when (for example) looking for width of peak in Hz
    """
    assert freq, dataDimRef

    sigmaBase   = hz2pnt(    0, dataDimRef )
    sigmaPoint  = hz2pnt( freq, dataDimRef )

    return abs( sigmaPoint - sigmaBase )

  def initSigmaParams(self):
    """ Set some initial default values for sigma """

    self._minSigmaHz = []
    self._maxSigmaHz = []

    if self.params.Ndim:
      for dim in range( self.params.Ndim ):
        self._minSigmaHz.append(  6. )
        self._maxSigmaHz.append( 28. )

  ###########################################################################

  def updateAll(self):

    self.updateSpectrumWindow()
    self.updatePeakListList()

    self.waiting = False

  def updateAfter(self, obj=None):

    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.updateAll)


def launchPeakSeparator():

  import Tkinter

  global top

  root = Tkinter.Tk()
  root.withdraw()
  top  = PeakSeparatorGui(root)

  top.update_idletasks()

  # root.mainloop()

if __name__ == '__main__':

  launchPeakSeparator()
