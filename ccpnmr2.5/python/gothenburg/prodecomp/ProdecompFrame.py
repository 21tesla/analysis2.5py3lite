# NBNB TBD get rid of self.shapeParams and use sself.dataDimRefDixt instead ??

import os

from ccpnmr.analysis.core.ExperimentBasic import getOnebondExpDimRefs
from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.WindowBasic import createSpectrumWindow
from ccpnmr.analysis.core.Util import getAnalysisSpectrum

#from ccpnmr.api.Analysis import ApiError
from memops.api.Implementation import ApiError

from ccp.general.Constants import chemShiftRefRatios

from memops.editor.Util import createDismissHelpButtonList

from memops.gui.Color import hsbToRgb, hexRepr
from memops.gui.Frame import Frame
from memops.gui.LabelFrame import LabelFrame
from memops.gui.TabbedFrame import TabbedFrame
from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.ButtonList import ButtonList
from memops.gui.Label import Label
from memops.gui.MessageReporter import showWarning, showInfo
from memops.gui.Entry import Entry
from memops.gui.IntEntry import IntEntry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.PulldownList import PulldownList
from memops.gui.LabelDivider import LabelDivider
from memops.gui.RadioButtons import RadioButtons

from gothenburg.prodecomp import Projection
from gothenburg.prodecomp import CcpnProdecomp
from gothenburg.prodecomp.PeaksToInterval import ccpnPeaksToInterval
shapeMatrixHeadings = ['Shape Name','Nucleus', 'Sweep Width\n(ppm)',
             'Carrier Freq\nOffset (ppm)', 'Left edge\nOffset (ppm)',
             'Comment']

# # # # #   B I G   T O - D O  P O I N T S   # # # # #
#
# Store all interval data - or write out XML as we go.
#
# Notifers - init/del spectra



# # # # #   L A T E R   T O - D O  P O I N T S  # # # # #
#
# Output of high-dim CCPN spectra
# - coordinate with single decomposition XML file
#
# 1D shaped viewer using spectrum windows
# - Can assign etc
#
# Make the automated interval generation smarter

def prodecompCcpnMacro(argServer):

  project = argServer.getProject()
  popup = ProdecompPopup(argServer.parent, project)
  popup.open()

class ProdecompPopup(BasePopup):

  def __init__(self, parent, ccpnProject):

    self.parent    = parent
    self.analysis  = parent #NB so that it works in both ExtendNmr and Analysis
    self.ccpnProject = ccpnProject

    BasePopup.__init__(self, parent=parent, title='PRODECOMP')


  def body(self, guiFrame):

    guiFrame.grid_rowconfigure(0, weight=1)
    guiFrame.grid_columnconfigure(0, weight=1)

    frame = ProdecompFrame(guiFrame, basePopup=self, ccpnProject=self.ccpnProject)
    frame.grid(row=0, column=0, sticky='nsew')
    frame.grid_rowconfigure(0, weight=1)
    frame.grid_columnconfigure(0, weight=1)

    #self.geometry('650x600')
    width = int(self.winfo_screenwidth()*0.666)
    height = int(self.winfo_screenheight())
    self.geometry('%dx%d+%d+%d' % (width,height,0,0))

    print frame.printOutDocString

    self.update_idletasks()


#def getPpmRange(sw, tfo, pointsRange, numPoints):
#
#  bVal = tfo + sw/2.0
#  fac = (-sw/numPoints)
#  return [x*fac+bVal for x in pointsRange]


#def getPpmRange(refppm, refpt, ppmPerPoint, pointsRange):
#  
#  return [refppm + (x - refpt)* ppmPerPoint for x in pointsRange]


class ProdecompFrame(Frame):

  # doc string to print out when starting
  printOutDocString = CcpnProdecomp.prodecomp.publicDocumentation

  def __init__(self, guiParent, basePopup, ccpnProject=None):
    # Base popup required to handle notification of data model changes
    # e.g. new peak lists, so that the GUI can update to the latest
    # state

    self.basePopup = basePopup
    self.guiParent = guiParent
    self.project = ccpnProject

    if ccpnProject:
      self.nmrProject = ccpnProject.currentNmrProject
    else:
      self.nmrProject = None

    self.spectra = []
    self.output  = None
    self.graphs  = []
    self.peakList = None
    self.refType = None
    self.shapeParams = {}
    self.dataDimRefDict = {}
    self.shapeNames = []
    self.acqName = 'HN'
    self.defsMatrix = []
    self.intervalId = 0
    #self.interval = None
    self.decompositions = []

    # NB values in this dictionary contain numpy objects! CAREFUL!
    self.prodecompOutput = {}

    Frame.__init__(self, guiParent)

    self.grid_rowconfigure(0, weight=1)
    self.grid_columnconfigure(0, weight=1)

    options = ['Input Data Sources','Exp Parameters','Decomposition','Output & Results']
    self.tabbedFrame = TabbedFrame(self, options=options)
    self.tabbedFrame.grid(row=0,column=0,sticky='nsew')
    frameA, frameB, frameC, frameD = self.tabbedFrame.frames

    frameA.grid_columnconfigure(0, weight=1)
    frameA.grid_rowconfigure(0, weight=1)
    frameB.grid_columnconfigure(0, weight=1)
    frameB.grid_rowconfigure(1, weight=1)
    frameC.grid_columnconfigure(0, weight=1)
    frameC.grid_rowconfigure(3, weight=1)
    frameD.grid_columnconfigure(0, weight=1)
    frameD.grid_rowconfigure(1, weight=1)

    #
    # Input Data Sources
    #

    headingList = ['Spectrum', 'Active?', 'Acq Nucl', 'nPoints']
    self.inputDataMatrix = ScrolledMatrix(frameA, multiSelect=True,
                        callback=self._selectSpectrum,
                        headingList=headingList)
    self.inputDataMatrix.grid(row=0,column=0,sticky='nsew')

    texts = ['Activate\nSelected',
         'Inactivate\nSelected',
         'Remove\nSelected']
    commands = [self.activateSelected,
          self.inactivateSelected,
          self.removeSelected]
    buttonList = ButtonList(frameA,texts=texts,
                commands=commands,
                expands=True)
    buttonList.grid(row=1,column=0,sticky='ew')

    #
    # Exp Parameters
    #

    frame = LabelFrame(frameB, text='Input Parameters')
    frame.grid(row=0,column=0,sticky='ew')
    frame.grid_columnconfigure(3, weight=1)
    frame.grid_rowconfigure(1, weight=1)

    label = Label(frame, text='Number of iterations')
    label.grid(row=0,column=0,sticky='w')
    self.numIterEntry = IntEntry(frame, text=2, width=6)
    self.numIterEntry.grid(row=0,column=1,sticky='w')

    label = Label(frame, text='Regularisation factor')
    label.grid(row=0,column=2,sticky='w')
    self.regFacEntry = FloatEntry(frame, text=0.1, width=8)
    self.regFacEntry.grid(row=0,column=3,sticky='w')

    #
    #self.sweepEntry   = FloatEntry(self, text=None, returnCallback=self.setSweep, width=8)
    #self.tOffsetEntry = FloatEntry(self, text=None, returnCallback=self.setTransOffset, width=8)


    #editWidgets    = [None, None, self.sweepEntry,
    #          self.tOffsetEntry, None,]
    #editGetCallbacks = [None, None, self.getSweep,
    #          self.getTransOffset,  None,]
    #editSetCallbacks = [None, None, self.setSweep,
    #          self.setTransOffset,  None,]
    editWidgets    = [None, None, None, None, None, None,]
    editGetCallbacks = [None, None, None, None, None, None,]
    editSetCallbacks = [None, None, None, None, None, None,]

    frame = LabelFrame(frameB, text='Shape Referencing')
    frame.grid(row=1,column=0,sticky='nsew')
    frame.grid_columnconfigure(0, weight=1)
    frame.grid_rowconfigure(0, weight=1)

    headingList = shapeMatrixHeadings
    self.shapeParamMatrix = ScrolledMatrix(frame, multiSelect=True,
                         callback=self.selectRefType,
                         editWidgets=editWidgets,
                         editGetCallbacks=editGetCallbacks,
                         editSetCallbacks=editSetCallbacks,
                         headingList=headingList)
    self.shapeParamMatrix.grid(row=0,column=0,sticky='nsew')

    #
    # Calculation
    #

    div = LabelDivider(frameC, text='Interval Generation Parameters')
    div.grid(row=0,column=0,sticky='ew')

    frame = Frame(frameC)
    frame.grid(row=1,column=0,sticky='ew')

    #

    label = Label(frame, text='Guide Peak List:')
    label.grid(row=0,column=0,sticky='e')

    self.peakListPulldown = PulldownList(frame, callback=self.selectPeakList)
    self.peakListPulldown.grid(row=0,column=1,sticky='ew')

    label = Label(frame, text='Interval Width (points):')
    label.grid(row=0,column=2,sticky='w')

    self.iWidthEntry = FloatEntry(frame, text=5, width=8)
    self.iWidthEntry.grid(row=0,column=3,sticky='w')

    label = Label(frame, text='Shift Width Scale:')
    label.grid(row=0,column=4,sticky='w')

    self.sScaleEntry = FloatEntry(frame, text=10, width=8)
    self.sScaleEntry.grid(row=0,column=5,sticky='w')

    #

    label = Label(frame, text='Width @ Half-height:')
    label.grid(row=1,column=0,sticky='w')

    self.hwhhEntry = FloatEntry(frame, text=2, width=8)
    self.hwhhEntry.grid(row=1,column=1,sticky='w')

    label = Label(frame, text='Intensity Threshold:')
    label.grid(row=1,column=2,sticky='w')

    self.iThreshEntry = FloatEntry(frame, text=0.1, width=8)
    self.iThreshEntry.grid(row=1,column=3,sticky='w')

    label = Label(frame, text='Intensity Type:')
    label.grid(row=1,column=4,sticky='w')

    self.iTypePulldown = PulldownList(frame, texts=['height','volume'])
    self.iTypePulldown.grid(row=1,column=5,sticky='ew')

    #

    label = Label(frame, text='Peak Shape Model:')
    label.grid(row=2,column=0,sticky='w')

    self.peakModelSelect = RadioButtons(frame, entries=['Gaussian','Lorentzian'])
    self.peakModelSelect.grid(row=2, column=1, columnspan=5, sticky='w')

    #


    div = LabelDivider(frameC, text='Interval Table')
    div.grid(row=2,column=0,sticky='ew')

    self.startEntry = IntEntry(self, text=None, returnCallback=self.setIntervalStart, width=8)
    self.endEntry   = IntEntry(self, text=None, returnCallback=self.setIntervalEnd,   width=8)
    self.compEntry  = IntEntry(self, text=None, returnCallback=self.setNumComps, width=8)

    editWidgets    = [None, None, self.startEntry, self.endEntry, self.compEntry]
    editGetCallbacks = [None, None, self.getIntervalStart, self.getIntervalEnd, self.getNumComps]
    editSetCallbacks = [None, None, self.setIntervalStart, self.setIntervalEnd, self.setNumComps]

    headingList = ['#','Complete?','Interval\nStart','Interval\nEnd','Num Components']
    self.intervalMatrix = ScrolledMatrix(frameC, multiSelect=True,
                                         editWidgets=editWidgets,
                                         editGetCallbacks=editGetCallbacks,
                                         editSetCallbacks=editSetCallbacks,
                                         callback=self._selectInterval,
                                         headingList=headingList)

    self.intervalMatrix.grid(row=3, column=0,sticky='nsew')

    texts = ['Add\nInterval','Remove\nInterval','Generate\nIntervals']
    commands = [self.addInterval, self.removeInterval, self.generateIntervals]
    buttonList = ButtonList(frameC,texts=texts,
                            commands=commands,
                            expands=True)
    buttonList.grid(row=4,column=0,sticky='ew')

    texts = ['Plot\n1D Spectra','Plot\nGraphs','Run\nSelected','Run\nAll']
    commands = [self.plotWindows, self.plot, self.runInterval,self.runProdecompAll]
    buttonList = ButtonList(frameC,texts=texts,
                            commands=commands,
                            expands=True)
    buttonList.grid(row=5,column=0,sticky='ew')

    
    #
    # Output & Results
    #

    frame = Frame(frameD)
    frame.grid(row=0, column=0,sticky='ew')
    frame.grid_columnconfigure(2, weight=1)
    frame.grid_rowconfigure(1, weight=1)

    label = Label(frame, text='Displayed Interval: ')
    label.grid(row=0, column=0,sticky='w')
    self.displayIntervalPulldown = PulldownList(frame, callback=self.selectPulldownInterval)
    self.displayIntervalPulldown.grid(row=0, column=1, sticky='w')

    label = Label(frame, text='Output XML file (.usf3 format): ')
    label.grid(row=0, column=2,sticky='e')
    self.fileNameEntry = Entry(frame, text='ProdecompOut.usf3', width=32)
    self.fileNameEntry.grid(row=0, column=3,sticky='ew')
    self.graphFrame = None

    texts = ['Write Out XML',]
    commands = [self.saveXml,]
    buttonList = ButtonList(frame,texts=texts,
                            commands=commands,
                            expands=True)
    buttonList.grid(row=0,column=4, sticky='w')

    #
    if isinstance(self.parent, BasePopup):

      bottomButtons = createDismissHelpButtonList(self.tabbedFrame.sideFrame,
                                                  texts=texts, commands=commands,
                                                  expands=True)
      bottomButtons.grid(row=0,column=0,sticky='e')

    self.updateAll()
    self.administerNotifiers(self.basePopup.registerNotify)

  def administerNotifiers(self, notifyFunc):

    for func in ('__init__','delete','setName'):
      notifyFunc(self.updateAllAfter, 'ccp.nmr.Nmr.Experiment', func)
      notifyFunc(self.updateAllAfter, 'ccp.nmr.Nmr.DataSource', func)
      notifyFunc(self.updatePeakListsAfter, 'ccp.nmr.Nmr.PeakList', func)

  def updateAllAfter(self, obj):

    self.after_idle(self.updateAll)

  def updatePeakListsAfter(self, peakList):

    self.after_idle(self.updatePeakLists)

  def open(self):

    self.updateAll()

    Frame.open(self)

  def updateAll(self, project=None):

    if project:
      self.project = project
      self.nmrProject = project.currentNmrProject

      if not self.nmrProject:
        self.nmrProject = project.newNmrProject(name=project.name)

    if not self.project:
      return

    self.update()
    self.updatePeakLists()
    self.updateShapeParams()


  def selectRefType(self, obj, row, col):

    self.refType = obj

  #def getSweep(self, refType):
  #
  #  self.sweepEntry.set(refType[1])

  #def getTransOffset(self, refType):
  #
  #  self.tOffsetEntry.set(refType[2])

  #def setSweep(self, *event):
  #
  #  self.refType[1] = self.sweepEntry.get() or 0.0
  #  self.shapeParams[self.refType[0]] = self.refType[1:]
  #  self.updateShapeParams()

  #def setTransOffset(self, *event):
  #
  #  self.refType[2] = self.tOffsetEntry.get() or 0.0
  #  self.shapeParams[self.refType[0]] = self.refType[1:]
  #  self.updateShapeParams()

  def selectPulldownInterval(self, interval):

    #self.interval = interval
    self.intervalMatrix.selectObject(interval)
    self.plot()

  def selectPeakList(self, obj):

    self.peakList = obj

  def generateIntervals(self):

    iWidth  = self.iWidthEntry.get() or 1.0
    sScale  = self.sScaleEntry.get() or 1.0
    hwhh  = self.hwhhEntry.get() or 1.0
    iThresh = self.iThreshEntry.get() or 0.1
    pModel  = self.peakModelSelect.get()[0] # 'G' or 'L'
    iType   = self.iTypePulldown.getText()

    if self.peakList:

      intervals = ccpnPeaksToInterval(self.peakList, iHalfRange=iWidth/2.0,
                                      shiftScale=sScale, peakShapeModel=pModel,
                                      halfWidthHalfHeight=hwhh,
                                      intensThreshold=iThresh,
                                      intensityType=iType)

      for num, start, end, nComp, ppmH, ppmN in intervals: #@UnusedVariable
        # The False is the 'Done' toggle
        self.decompositions.append([start, end, nComp, False])

      self.updateIntervals()


  def updatePeakLists(self):

    index = 0
    peakLists = []
    texts = []

    if self.project:
      nmrProject = self.project.currentNmrProject

      for experiment in nmrProject.experiments:

        eName = experiment.name

        # get list of ExpDimRef pairs connected by onebond transfers
        expDimRefList = getOnebondExpDimRefs(experiment)

        for dataSource in experiment.dataSources:
          dataDims = dataSource.sortedDataDims()
          
          if (dataSource.dataType == 'processed' and dataSource.numDim == 2
              and len(dataDims) == 2):   #necessary with some unfinished specs

            expDimRefs0 = set(x.expDimRef for x in dataDims[0].dataDimRefs)
            expDimRefs1 = set(x.expDimRef for x in dataDims[1].dataDimRefs)

            # check for onebond transfer in use
            for xdr0,xdr1 in expDimRefList:
              if (xdr0 in expDimRefs0 and xdr1 in expDimRefs1 or
                xdr1 in expDimRefs0 and xdr0 in expDimRefs1):
                break
            else:
              # datasource not used
              continue

            # dataSource used - look for peakLists
            sName = dataSource.name

            for peakList in dataSource.peakLists:
              if peakList.peaks:
                peakLists.append(peakList)
                texts.append('%s:%s:%d' % (eName,sName,peakList.serial))

    if peakLists:
      if self.peakList not in peakLists:
        self.peakList = peakLists[0]

      index = peakLists.index(self.peakList)

    self.peakListPulldown.setup(texts=texts, objects=peakLists, index=index)

  def runInterval(self):

    #self.runProdecomp([self.interval,])
    self.runProdecomp(self.intervalMatrix.currentObjects)

  def getIntervalStart(self, interval):

    self.startEntry.set(interval[0])

  def getIntervalEnd(self, interval):

    self.endEntry.set(interval[1])

  def getNumComps(self, interval):

    self.compEntry.set(interval[2])

  def setIntervalStart(self, *event):
    
    interval = self.intervalMatrix.currentObject
    if interval is None:
      return

    value = self.startEntry.get() or interval[0]

    if interval[1] is not None:
      end = interval[1]-1
      value = min(end, value)

    interval[0] = max(0, value)
    self.updateIntervals()

  def setIntervalEnd(self, *event):
    
    interval = self.intervalMatrix.currentObject
    if interval is None:
      return

    value = self.endEntry.get() or interval[1]

    if interval[0] is not None:
      start = interval[0]+1
      value = max(value, start)

    interval[1] = value
    self.updateIntervals()

  def setNumComps(self, *event):
    
    interval = self.intervalMatrix.currentObject
    if interval is None:
      return

    self.updateIntervals()
    value = self.compEntry.get() or interval[2]
    interval[2] = max(1, value)
    self.updateIntervals()

  def updateIntervals(self):

    textMatrix  = []
    objectList  = []
    colorMatrix = []

    texts = []
    intervals = []

    for i, decomposition in enumerate(self.decompositions):

      color = [None,] * 5
      if decomposition[3]:
        complete = 'Yes'
        color[1] = '#489A48'
        texts.append('%d' % (i+1,))
        intervals.append(decomposition)

      else:
        complete = 'No'


      datum  = [i+1,complete] + decomposition[:3]

      textMatrix.append(datum)
      objectList.append(decomposition)
      colorMatrix.append(color)

    self.intervalMatrix.update(textMatrix=textMatrix,
                               objectList=objectList,
                               colorMatrix=colorMatrix)

    interval = self.displayIntervalPulldown.getObject()
    if interval in intervals:
      index = intervals.index(interval)
    else:
      index = 0

    self.displayIntervalPulldown.setup(texts, intervals, index)
    #self.interval = self.displayIntervalPulldown.getObject()
    if intervals:
      interval = intervals[index]
    else:
      interval = None
    self.intervalMatrix.selectObject(interval)

  def addInterval(self):

    self.decompositions.append([None, None ,1,None])
    self.updateIntervals()

  def removeInterval(self):
    
    interval = self.intervalMatrix.currentObject
    if interval:
      self.decompositions.remove(interval)
      self.intervalMatrix.keyPressEscape()
      self.updateIntervals()

  def _selectInterval(self, obj, row, col):

    self.displayIntervalPulldown.set(obj)


  def updateShapeParams(self):

    swTolerance = 0.001
    ###tfoTolerance = 0.001
    tfoTolerance = 0.1

    textMatrix  = []
    objectList  = []
    colorMatrix = []

    headingList = shapeMatrixHeadings

    nColumns = len(headingList)
    dataSources = [spec for spec in self.spectra if spec.isProdecompActive]

    if dataSources and self.shapeNames:

      self.numPoints = [None, None]
      self.numPointsOrig = [None, None]

      acqString = self.acqName + ' (acq)'
      tags = [self.acqName] + self.shapeNames
      dataDimRefDict = self.dataDimRefDict

      # first get dataDimRefs list for acquisition dimension
      activeExps = set(spec.experiment for spec in dataSources)
      expDimRefs = set()
      sf = 0.0
      for xp in activeExps:
        expDim = xp.findFirstExpDim(isAcquisition=True)
        for expDimRef in expDim.sortedExpDimRefs():
          if expDimRef.measurementType in ('shift', 'Shift'):
            expDimRefs.add(expDimRef)
            if not sf:
              sf = expDimRef.sf
      #
      self.acqSf = sf

      ll = []
      sizes = set()
      origsizes = set()
      for dataSource in dataSources:
        dataDim = dataSource.sortedDataDims()[0]
        sizes.add(dataDim.numPoints)
        if hasattr(dataDim, 'numPointsOrig'):
          origsizes.add(dataDim.numPointsOrig)
        for dataDimRef in dataDim.sortedDataDimRefs():
          if dataDimRef.expDimRef in expDimRefs:
            ll.append(dataDimRef)
      #
      dataDimRefDict[self.acqName] = ll
      if len(sizes) == 1:
        self.numPoints[0] = sizes.pop()
      if len(origsizes) == 1:
        self.numPointsOrig[0] = origsizes.pop()

      # now get dataDimRefs lists for other shapes
      for shapeName in self.shapeNames:
        dataDimRefDict[shapeName] = []

      sizes = set()
      origsizes = set()
      for dataSource in dataSources:
        dataDim = dataSource.sortedDataDims()[1]
        sizes.add(dataDim.numPoints)
        if hasattr(dataDim, 'numPointsOrig'):
          origsizes.add(dataDim.numPointsOrig)
        for dataDimRef in dataDim.sortedDataDimRefs():

          shapeName = dataDimRef.expDimRef.displayName
          ll = dataDimRefDict.get(shapeName)
          if ll is not None:
            ll.append(dataDimRef)
      if len(sizes) == 1:
        self.numPoints[1] = sizes.pop()
      if len(origsizes) == 1:
        self.numPointsOrig[1] = origsizes.pop()

      if None in self.numPointsOrig:
        print ('WARNING, selected spectra were not the same original shape')
      if None in self.numPoints:
        print ('WARNING, selected spectra were not the same shape')
      else:
        for tag in tags:

          # set textMatrix
          if textMatrix:
            # normal case
            datum = [tag]
          else:
            # first time, special case
            datum = [acqString]
          textMatrix.append(datum)

          bad = False
          commList = []
          dataDimRefs = dataDimRefDict[tag]
          isotopeCodes = set()
          llsw = []
          lltfo = []
          llstartppm = []
          for dataDimRef in dataDimRefs:
            isotopeCodes.add(dataDimRef.expDimRef.isotopeCodes)
            llsw.append(dataDimRef.spectralWidth)
            ddim = dataDimRef.dataDim
            lltfo.append(dataDimRef.pointToValue(ddim.numPointsOrig/2
                                                 - ddim.pointOffset + 1))
            llstartppm.append(dataDimRef.pointToValue(1.0))
          if len(isotopeCodes) == 1:
            datum.append(isotopeCodes.pop()[0])
          else:
            bad = True
            commList.append("Isotopes")
            datum.append(str(tuple(isotopeCodes)))

          maxval = max(llsw)
          minval = min(llsw)
          avval = sum(llsw)/len(llsw)
          if abs(maxval - minval) > swTolerance:
            bad = True
            commList.append("SW %.3f - %.3f" % (maxval, minval,))
          datum.append(avval)

          maxval = max(lltfo)
          minval = min(lltfo)
          avval = sum(lltfo)/len(lltfo)
          if abs(maxval - minval) > tfoTolerance:
            bad = True
            commList.append("Offset %.3f - %.3f" % (maxval, minval,))
          datum.append(avval)

          # set objectList
          ll = self.shapeParams[tag] = datum[-3:]
          objectList.append(ll)

          avval = sum(llstartppm)/len(llstartppm)
          datum.append(avval)

          # set colors
          if bad:
            datum.append("Bad " + ', '.join(commList))
            colorMatrix.append(nColumns * ['#FF0000'])
          else:
            datum.append(None)
            colorMatrix.append(nColumns * [None])

    self.shapeParamMatrix.update(textMatrix=textMatrix,
                                 objectList=objectList,
                                 headingList=headingList,
                                 colorMatrix=colorMatrix)

  def quit(self):

    self.guiParent.parent.destroy()

  def destroy(self):

    self.administerNotifiers(self.basePopup.unregisterNotify)
    Frame.destroy(self)

  def runProdecompAll(self):

    self.runProdecomp(self.decompositions)

  def runProdecomp(self, intervals):

    inpDataSources = self.spectra
    regFact = self.regFacEntry.get() or 0.0
    nIter = self.numIterEntry.get() or 0

    errs = []

    dataSources = [spec for spec in inpDataSources if spec.isProdecompActive]

    if not inpDataSources:
      errs.append('No input data sources selected')

    elif not dataSources:
      errs.append('No active data sources')

    if regFact <= 0:
      errs.append('Regularisation factor must be grater than zero')

    if nIter < 1:
      errs.append('Must have at least one iteration')

    for decomposition in intervals:
      start, end, comp, null = decomposition #@UnusedVariable

      if start < 1:
        errs.append('Interval start must be positive integer')

      if end < 1:
        errs.append('Interval end must be positive integer')

      if comp < 1:
        errs.append('Component must be positive integer')

      if errs:
        msg = '\n'.join(errs)
        showWarning('Could not run PRODECOMP',
              'The following problems were encountered:\n\n' + msg)
        return

      # WARNING - fdir and f are numpy objects
      fdir, f = CcpnProdecomp.runProdecomp(dataSources, self.defsMatrix,
                                           (start,end), (comp,comp),
                                           regFact, nIter)


      intervalId = self.intervalId + 1
      self.intervalId = intervalId
      self.prodecompOutput[intervalId] = (fdir, f, self.acqName,
                                          self.shapeNames)
      decomposition[3] = intervalId

    showInfo("Info","Decomposition is done!", parent=self)

    self.updateIntervals()
    self.selectPulldownInterval(decomposition)
    #self.plot()

  def activateSelected(self):

    for spec in self.inputDataMatrix.currentObjects:
      spec.isProdecompActive = True

    self.update()
    self.updateShapeParams()

  def inactivateSelected(self):

    for spec in self.inputDataMatrix.currentObjects:
      spec.isProdecompActive = False

    self.update()
    self.updateShapeParams()

  def removeSelected(self):

    remove = self.inputDataMatrix.currentObjects

    self.spectra = [spec for spec in self.spectra if spec not in remove]

    self.update()

  def _selectSpectrum(self, obj, row, col):

    pass


  def update(self):
    """ NB this assumes 2D spectra, first dim is acquisition,
    only dimensionScalings with a single factor and displayName set.
    Inappropriate cases are filtered out in
    Projection functions
    """

    textMatrix  = []
    objectList  = []
    colorMatrix = []
    headingList = ['Spectrum','Active?', 'nPoints']

    spectra = Projection.getProjectionSpectra(self.nmrProject)
    self.spectra = spectra

    allShapeNames, allDefsMatrix = Projection.getProjectionData(spectra)

    headingList = headingList + allShapeNames

    ff = Projection.formatScalingFactor
    for i, spec in enumerate(spectra):

      dataDims = spec.sortedDataDims()

      if not hasattr(spec, 'isProdecompActive'):
        spec.isProdecompActive = True

      name = '%s:%s' % (spec.experiment.name, spec.name)
      nPoints = '%s:%s' % (dataDims[0].numPoints, dataDims[1].numPoints)
      datum  = [name,
                spec.isProdecompActive and 'Yes' or 'No',
                nPoints]

      datum.extend(ff(x) for x in allDefsMatrix[i])

      colors = [None] * len(datum)
      if spec.isProdecompActive:
        colors[1] = '#B0FFB0'

      textMatrix.append(datum)
      objectList.append(spec)
      colorMatrix.append(colors)



    # now set up active specs and generate calculation data

    # get projection spectra, active spectra, and active experiments
    # and do first pass
    activeExps = set()
    activeSpecs = []
    for spec in spectra:
      if spec.isProdecompActive:
        activeExps.add(spec.experiment)
        activeSpecs.append(spec)

    # find and set shapeNames and definitions matrix
    shapeNames, defsMatrix = Projection.getProjectionData(activeSpecs)
    self.shapeNames = shapeNames
    self.defsMatrix = defsMatrix

    # find and set acquisition dimension name
    acqName = '?'
    acqNames = set()
    for xp in activeExps:
      expDim = xp.findFirstExpDim(isAcquisition=True)
      for expDimRef in expDim.sortedExpDimRefs():
        if expDimRef.measurementType in ('Shift', 'shift'):
          sname = expDimRef.displayName
          if sname:
            acqNames.add(sname)
          else:
            try:
              acqName = expDimRef.isotopeCodes[0]
            except IndexError:
              raise ApiError(" No isotopeCodes found for %s" % expDimRef)

    if acqNames:
      if len(acqNames)!= 1:
        print ('WARNING, inconsistent displayNames for acquisition dimension:',
             tuple(acqNames))
      acqName = acqNames.pop()

    if acqName in shapeNames:
      print ('WARNING, acquisition display name % duplicates other shape:',
           acqNames)

    self.acqName = acqName

    # update at end
    self.inputDataMatrix.update(textMatrix=textMatrix,
                                objectList=objectList,
                                headingList=headingList,
                                colorMatrix=colorMatrix)

  def getDataColors(self, n):

    colors = []
    s = 1.0
    v = 0.5
    for i in range(int(n)):
      h = i/float(n)
      r,g,b = hsbToRgb(h,s,v)
      colors.append(hexRepr(r,g,b))

    return colors

  def findMatchingWindowSpectrum(self, windowPane):

    for view in windowPane.spectrumWindowViews:
      analysisSpectrum = view.analysisSpectrum
      spectrum = analysisSpectrum.dataSource
      if spectrum.name == windowPane.name:
        break
    else:
      spectrum = None

    return spectrum

  def plotWindows(self):

    analysis = self.basePopup.analysis
    if not analysis:
      msg = "CcpNmr Analysis not available to show spectra"
      showWarning('Failure', msg, parent=self)
      return

    interval = self.intervalMatrix.currentObject
    if not (interval and interval[3]):#  i.e. Not calculated
      return

    intervalA, intervalB, numComps, intervalId = interval
    output = self.prodecompOutput[intervalId]
    outfdir, outf, acqName, shapeNames = output #@UnusedVariable
    isotope, sw, tfo = self.shapeParams[self.acqName]

    # Acquisition dim first
    isotopes = [isotope,] + [self.shapeParams[label][0] for label in shapeNames]

    numDim = len(isotopes)
    analysisProject = self.project.currentAnalysisProject
    nmrProject = analysisProject.nmrProject

    # try and find an existing experiment, otherwise create one
    experiments = nmrProject.findAllExperiments(name='Prodecomp', numDim=numDim)
    for experiment in experiments:
      expDims = experiment.sortedExpDims()
      for (n, expDim) in enumerate(expDims):
        expDimRef = expDim.findFirstExpDimRef()
        isotopeCodes = expDimRef.isotopeCodes
        if len(isotopeCodes) == 1 and isotopeCodes[0] != isotopes[n]:
          break # this experiment no good
      else:
        break # this experiment good
    else:
      # need to create an experiment since no matching one found
      experiment = nmrProject.newExperiment(name='Prodecomp', numDim=numDim)
      expDims = experiment.sortedExpDims()
      for (n, expDim) in enumerate(expDims):
        isotope = isotopes[n]
        if n == 0:
          sfRatio = chemShiftRefRatios.get(isotope,0.0)
          if sfRatio:
            sfRatio = self.acqSf/sfRatio
          sf = self.acqSf
        else:
          sf = chemShiftRefRatios.get(isotope,0.0)*sfRatio
        expDim.newExpDimRef(sf=sf, isotopeCodes=(isotope,))

    valueAxis = analysisProject.findFirstAxisType(name='value')

    """
    axisTypes = []
    windows = []
    for name in isotopes:
      # In a hurry, I assume this exists
      axisTypes.append(analysisProject.findFirstAxisType(name=name))
    
    windows1d = []
    for window in analysisProject.spectrumWindows:
      for windowFrame in window.spectrumWindowFrame:
        axisPanels = window.axisPanels
        if len(axisPanels) != 2:
          continue

        if not window.findFirstAxisPanel(axisType=valueAxis):
          continue

        windows1d.append(windowFrame)

    for i, name in enumerate(isotopes):
      windowName = '%s_%d' % (name, i+1)
      for window in windows1d:
        if (window.name == windowName) and window.findFirstAxisPanel(axisType=axisTypes[i]):
          windows.append(window)
          windows1d.remove(window)
          break

      else:
        project = analysisProject.root
        window = createSpectrumWindow(project, windowName, [(axisTypes[i],valueAxis),],
                                      spectrum=None, background=None, width=None,
                                      height=None, ncols=1, nrows=1, regions=None)
        windows.append(window)
    """

    window = analysisProject.findFirstSpectrumWindow(name='Prodecomp')
    if not window:
      axisTypePairs = []
      for name in isotopes:
        axisType = analysisProject.findFirstAxisType(name=name)
	axisTypePairs.append( (axisType, valueAxis) ) 
	 
      window = createSpectrumWindow(analysisProject.root, 'Prodecomp',
                                    axisTypePairs, spectrum=None,
				    ncols=1, nrows=1, regions=None)   

    self.update_idletasks()
    for i, windowPane in enumerate(window.sortedSpectrumWindowPanes()):
      windowPane.name = '%s_%d' % (isotopes[i], i+1)
    
    # geometry
    
    """width = int(self.winfo_screenwidth())/2

    nWin = len(windows)
    if nWin < 6:
      heightJ = int(self.winfo_screenheight())/nWin
      for i, window in enumerate(windows):
        popup = analysis.getWindowPopup(window.name, doOpen=True)
        x = width
        y = i*heightJ
        popup.geometry('%dx%d+%d+%d' % (width,heightJ,x,y))

    else:
      j = nWin/2
      k = nWin-j
      heightJ = int(self.winfo_screenheight())/j
      heightK = int(self.winfo_screenheight())/k

      for i in range(j):
        window = windows[i]
        popup = analysis.getWindowPopup(window.name, doOpen=True)
        x = 0
        y = i*heightJ
        popup.geometry('%dx%d+%d+%d' % (width,heightJ,x,y))

      for i in range(k):
        window = windows[i+j]
        popup = analysis.getWindowPopup(window.name, doOpen=True)
        x = width
        y = i*heightK
        popup.geometry('%dx%d+%d+%d' % (width,heightK,x,y))
    """
    
    width = int(self.winfo_screenwidth())/3
    height = int(self.winfo_screenheight())
    popup = analysis.getWindowPopup(window.name, doOpen=True)
    popup.geometry('%dx%d+%d+%d' % (width,height,0,0))
    
    

    # scrollbars
    """for window in windows:
      for label in ('x', 'y'):
        axisPanel = window.findFirstAxisPanel(label=label)
        axisPanel.isVisible = False
    """
    
    for windowPane in window.sortedSpectrumWindowPanes():
      for label in ('x', 'y'):
        axisPanel = windowPane.findFirstAxisPanel(label=label)
        axisPanel.isVisible = False
      
     
    for spec in self.spectra:
      if spec.isProdecompActive:
        break
    
    dataDims = spec.sortedDataDims()
    numPointsDir = dataDims[0].numPoints - 1
    numPointsShape = dataDims[1].numPoints - 1

    # Direct dim

    # Find H window for direct dim
    windowPanes = window.sortedSpectrumWindowPanes()
    windowPane = windowPanes[0]

    dataDimRef = self.dataDimRefDict[self.acqName][0]
    xPpmPoints = [dataDimRef.pointToValue(x) 
                  for x in range(intervalA,intervalB+1)]
    #xPpmPoints =  getPpmRange(sw,tfo, range(intervalA,intervalB+1),
    #                          numPointsDir )


    colors = self.getDataColors(numComps)

    # Find new or existing spectrum
    spectrum = self.findMatchingWindowSpectrum(windowPane)

    valuesList = []
    for i in range(numComps):
      yValuePoints = [outfdir[j][i] for j in range(len(xPpmPoints))]
      valuesList.append(yValuePoints)

    sw *= float(intervalB-intervalA+1) / numPointsDir
    self.plot1dWindowShape(windowPane, experiment, 1, spectrum,
                           sw, colors, xPpmPoints, valuesList)

    # Indirect dims

    pointsRange = range(numPointsShape)
    for k, label in enumerate(shapeNames):

      isotope, sw, tfo = self.shapeParams[label]
      dataDimRef = self.dataDimRefDict[label][0]

      # Find window for this isotope
      windowPane = windowPanes[k+1]
      spectrum = self.findMatchingWindowSpectrum(windowPane)

      if sw is None:
        xPpmPoints = pointsRange
      else:
        xPpmPoints = [dataDimRef.pointToValue(x) for x in pointsRange]
        #xPpmPoints = getPpmRange(sw, tfo, pointsRange, numPointsShape)

      valuesList = []
      for i in range(numComps):
        yValuePoints = [outf[j][k][i] for j in range(len(xPpmPoints))]
        valuesList.append(yValuePoints)

      self.plot1dWindowShape(windowPane, experiment, k+2, spectrum,
                             sw, colors, xPpmPoints, valuesList)

    # make sure only what should be visible is visible
    for windowPane in windowPanes:
      for view in windowPane.sortedSpectrumWindowViews():
        analysisSpectrum = view.analysisSpectrum
        spectrum = analysisSpectrum.dataSource
        if spectrum.name != windowPane.name:
          view.isInToolbar = view.isNegVisible = view.isPosVisible = view.isSliceVisible = False

  def plot1dWindowShape(self, windowPane, experiment, dim,
                        spectrum, sw, colors, xPpms, valuesList):

    expDim = experiment.findFirstExpDim(dim=dim)
    numPoints = len(valuesList[0])
    if sw is None:
      sw = numPoints
    else:
      # Prodecomp has sw in Hz, not ppm
      sw *= expDim.findFirstExpDimRef().sf

    # make a new spectrum if needed
    if not spectrum:
      spectrum = experiment.newDataSource(name=windowPane.name, numDim=1, dataType='processed')
      dataDim = spectrum.newFreqDataDim(dim=1, numPoints=numPoints, numPointsOrig=numPoints,
                                        isComplex=False, valuePerPoint=sw/numPoints, expDim=expDim)
      dataDim.newDataDimRef(refPoint=1, refValue=xPpms[0], expDimRef=expDim.findFirstExpDimRef())
    else:
      dataDim = spectrum.findFirstDataDim()
      dataDim.numPoints = numPoints
      dataDim.numPointsOrig = numPoints
      dataDim.valuePerPoint = sw/numPoints
      dataDim.findFirstDataDimRef().refValue = xPpms[0]

    minValue = maxValue = None
    for values in valuesList:
      if minValue is None:
        minValue = min(values)
        maxValue = max(values)
      else:
        minValue = min(minValue, min(values))
        maxValue = max(maxValue, max(values))
    delta = maxValue - minValue
    if delta > 0:
      delta *= 0.025
      minValue -= delta
      maxValue += delta
      windowPane.sliceRange = (minValue, maxValue)

    axisPanel = windowPane.findFirstAxisPanel(label='x')
    axisRegion = axisPanel.findFirstAxisRegion()
    axisRegion.region = (min(xPpms), max(xPpms))

    spectrum.valuesList = valuesList

    analysisSpectrum = getAnalysisSpectrum(spectrum)
    analysisSpectrum.components = range(len(valuesList)) # i.e. draw all components
    analysisSpectrum.sliceColors = colors

    analysis = self.basePopup.analysis
    analysis.finishInitSpectrum(spectrum)
    analysis.drawWindows()


  def plot(self):
    """ NBNB TOXIC CODE!
    Relies on numpy objects outfdir and outf
    """
    
    allPaneWidth = int(self.winfo_screenwidth()*0.1)
    allPaneHeight = 20

    self.updateIntervals()
    
    interval = self.displayIntervalPulldown.getObject()
    if not interval:
      return
    
    
    #for interval in self.intervalMatrix.currentObjects:
    #  if interval[3]:
    self.tabbedFrame.select(3)

    intervalA, intervalB, numComps, intervalId = interval
    output = self.prodecompOutput[intervalId]
    outfdir, outf, acqName, shapeNames = output

    isotope, sw, tfo = self.shapeParams[self.acqName]

    for spec in self.spectra:
      if spec.isProdecompActive:
        break

    dataDims = spec.sortedDataDims()
    numPointsDir = dataDims[0].numPoints - 1
    numPointsShape = dataDims[1].numPoints - 1
    
    dataDimRef = self.dataDimRefDict[self.acqName][0]
    ppmxfd = [dataDimRef.pointToValue(x) 
              for x in range(intervalA,intervalB+1)]
    #ppmxfd =  getPpmRange(sw,tfo,
    #                     range(intervalA,intervalB+1),
    #                     numPointsDir )
    
    dataSets = []
    names = []
    for i in range(numComps):
      dataSet = [(x, outfdir[j][i]) for j, x in enumerate(ppmxfd)]
      dataSets.append(dataSet)
      names.append('Component %d' % (i+1) )

    while self.graphs:
      graph = self.graphs.pop()
      del graph


    if self.graphFrame:
      index = self.graphFrame.selected
      for frame in self.graphFrame.frames:
        frame.grid_forget()
      self.graphFrame.grid_forget()
      del self.graphFrame
    else:
      index = None

    options = ['ALL', 'Direct Dim',] + shapeNames
    self.graphFrame = TabbedFrame(self.tabbedFrame.frames[3], options=options)
    self.graphFrame.grid(row=1, column=0, sticky='nsew')
    frames =  self.graphFrame.frames

    frames[0].grid_columnconfigure(1, weight=10)
    frames[0].grid_columnconfigure(0, weight=1)
    for frame in frames[1:]:
      frame.grid_columnconfigure(0, weight=1)
      frame.grid_rowconfigure(0, weight=1)

    if (index is not None) and (index < len(frames)):
      self.graphFrame.select(index)

    # N: a: -41.106/Np; b = 7024.99999999873/60.810663 + 41.1064623173443/2 = 136.076
    
    graphParams = {'xLabels':None, 'symbolSize':1, 
                   'xTicks':True, 'yTicks':False, 'graphType':'line',
                   'reverseX':True, 'dataColors':None, 'lineWidths':None, 
                  }
    
    title = 'Direct Dimension Components'
    graph = ScrolledGraph(frames[1], dataSets=dataSets,
                          width=300, height=200, title=title,
                          xLabel='[ppm]', zoom=1.5, showCoords=True, 
                          yLabel=acqName, dataNames=names, **graphParams)
    
    graph.grid(row=0, column=0, sticky='nsew')
    self.graphs.append(graph)

    # NB New Sep 2010 - All graphs tab
    graph2 = ScrolledGraph(frames[0], dataSets=dataSets,
                           width=allPaneWidth, height=allPaneHeight, 
                           showCoords=False, 
                           xLabel='[ppm]',  zoom=2.4,
                           xGrid=True, yGrid=False,
                           yLabel=acqName, **graphParams)
    graph2.grid(row=0, column=0, sticky='nsew')
    nextGrid0 = 1
    nextGrid1 = 0

    pointsRange = range(numPointsShape)
    for k, label in enumerate(shapeNames):

      isotope, sw, tfo = self.shapeParams[label] #@UnusedVariable
      dataDimRef = self.dataDimRefDict[label][0]

      if sw is None:
        ppmxf = pointsRange
      else:
        
        ppmxf = [dataDimRef.pointToValue(x) for x in pointsRange]
        #ppmxf =  getPpmRange(sw, tfo, pointsRange, numPointsShape)

      dataSets = []
      names = []
      for i in range(numComps):
        # NB New Sep 2010 Gothenburg. Scale shapes
        data = [outf[j][k][i] for j in pointsRange]
        for ii, x in enumerate(data):
          # get rid of NaN
          if x !=x:
            data[ii] = 0
        scale = sum(data)
        if scale == 0:
          scale = 1.0
        dataSet = [(x, data[j]/scale) for j, x in enumerate(ppmxf)]
        #dataSet = [(x, outf[j][k][i]) for j, x in enumerate(ppmxf)]
        
        dataSets.append(dataSet)
        names.append('Component %d' % (i+1) )

      title = 'Indirect Dimension Components for %s' % label
      graph = ScrolledGraph(frames[k+2], dataSets=dataSets,
                            width=300, height=200, title=title,
                            xLabel='[ppm]',  zoom=1.5, showCoords=True, 
                            yLabel=label, dataNames=names, **graphParams)
 
                          
      graph.grid(row=0, column=0,sticky='nsew')

      self.graphs.append(graph)
    
      # NB New Sep 2010 - All graphs tab
      if isotope == '13C':
        column = 1
        row = nextGrid1
        nextGrid1 += 1
        #dataNames = names
        dataNames = None
      else:
        column = 0
        row = nextGrid0
        nextGrid0 += 1
        dataNames = None
      graph2 = ScrolledGraph(frames[0], dataSets=dataSets, zoom=2.4,
                             width=allPaneWidth, height=allPaneHeight, 
                             yLabel=label, dataNames=dataNames, 
                             xLabel='[ppm]', showCoords=False, 
                             xGrid=True, yGrid=False,
                             **graphParams)
      graph2.grid(row=row, column=column, sticky='nsew')
    frames[0].grid_rowconfigure(max(nextGrid0,nextGrid1)-1 , weight=1)
      

  def saveXml(self):
    """ WARNING! TOXIC CODE!
    Relies on numpy objects outfdir and outf
    """

    # default values for isReconstructable, isResolved.
    defaultReconstructable = False
    defaultResolved = False

    # XML OUTPUT of PRODECOMP
    #Create the minidom document
    import xml.dom.minidom as minidom


    if not self.decompositions:
      return

    fileName = self.fileNameEntry.get()

    if fileName[-5:] != '.usf3':
      fileName += '.usf3'

    inpDataSources = self.spectra
    spec = inpDataSources[0]

    # find relevant
    intervals = [x for x in self.decompositions if x[3]]
    allComp = sum(x[2] for x in intervals)

    if intervals and allComp:


      dataSources = [spec for spec in inpDataSources if spec.isProdecompActive]
      aSet = set(spec.experiment for spec in dataSources)
      if len(aSet) > 1:
        ll = [(xp.serial, xp) for xp in aSet]
        experiments = [x[1] for x in sorted(ll)]
      else:
        experiments = list(aSet)

      refExperiment= getRefExperiment(experiments, self.acqName,
                       self.shapeNames)

      #aSet = set(xp.refExperiment for xp in experiments)
      #if len(aSet) == 1:
      #  refExperiment = aSet.pop()
      #else:
      #  refExperiment = None


      # start document
      run0 = intervals[0]
      doc = minidom.Document()

      # create Decomposition element
      intervalId = run0[3]
      output = self.prodecompOutput[intervalId]
      outfdir, outf, acqName, shapeNames = output
      decomposition = doc.createElement("Decomposition")
      decomposition.setAttribute("n", "0")
      decomposition.setAttribute("name", experiments[0].name)
      decomposition.setAttribute("nshapes", str(len(shapeNames)+1))
      decomposition.setAttribute("ncomponents", str(allComp))
      decomposition.setAttribute("nregions", str(len(intervals)))
      decomposition.setAttribute("nprojsets", str(len(experiments)))
      decomposition.setAttribute("reconstructable", 
                                 (defaultReconstructable and 'true' or 'false'))
      decomposition.setAttribute("resolved", 
                                 (defaultResolved and 'true' or 'false'))
      if refExperiment is not None:
        decomposition.setAttribute("refexperiment", refExperiment.name)

      decomposition.setAttribute("method","PRODECOMP")
      doc.appendChild(decomposition)

      com = doc.createComment("PRODECOMP Version 3.0 (07 Nov 2008 05:24:36)")
      decomposition.appendChild(com)

      # create Axis elements
      for ii, tag in enumerate([acqName]+shapeNames):
        isotope, sw, tfo = self.shapeParams[tag]
        dataDimRef = self.dataDimRefDict[tag][0]
        Axis = doc.createElement("Axis")
        Axis.setAttribute("a",str(ii))
        Axis.setAttribute("name",tag)
        Axis.setAttribute("nucleus", isotope)
        Axis.setAttribute("domain","freq")
        Axis.setAttribute("type","real")
        if ii == 0:
          sfRatio = chemShiftRefRatios.get(isotope,0.0)
          if sfRatio:
            sfRatio = self.acqSf/sfRatio
          Axis.setAttribute("size",str(self.numPoints[0]))
          Axis.setAttribute("origsize",str(self.numPointsOrig[0]))
          Axis.setAttribute("sfo",str(self.acqSf))
        else:
          Axis.setAttribute("size",str(self.numPoints[1]))
          Axis.setAttribute("origsize",str(self.numPointsOrig[1]))
          Axis.setAttribute("sfo",
                    str(chemShiftRefRatios.get(isotope,0.0)*sfRatio))
        Axis.setAttribute("swppm",str(sw))
        Axis.setAttribute("carppm",str(tfo))
        #Axis.setAttribute("startppm",str(tfo+sw/2.0)) # NBNB incorect
        Axis.setAttribute("startppm",str(dataDimRef.pointToValue(1.0)))
        decomposition.appendChild(Axis)

      # create ProjSet element and contents
      for jj, experiment in enumerate(experiments):

        specs = [x for x in dataSources if x.experiment is experiment]

        # get dimensions
        dimAxes = []
        for expDim in experiment.sortedExpDims()[1:]:
          for expDimRef in expDim.sortedExpDimRefs():
            if expDimRef.measurementType in ('shift', 'Shift'):
              tag = expDimRef.displayName
              if tag in self.shapeNames:
                dimAxes.append(self.shapeNames.index(tag)+1)
        dimAxes.sort()

        projset = doc.createElement("Projset")
        projset.setAttribute("s",str(jj))
        projset.setAttribute("name",experiment.name)
        projset.setAttribute("ndims",str(len(dimAxes)+1))
        projset.setAttribute("nproj",str(len(specs)))
        refExperiment = experiment.refExperiment
        if refExperiment is not None:
          projset.setAttribute("refexperiment", refExperiment.name)
        decomposition.appendChild(projset)

        # projdim 0
        projdim = doc.createElement("Projdim")
        projdim.setAttribute("d","0")
        projdim.setAttribute("axes","0")
        projset.appendChild(projdim)

        # projdim 1

        projdim = doc.createElement("Projdim")
        projdim.setAttribute("d","1")
        projdim.setAttribute("axes"," ".join(str(x) for x in dimAxes))
        projset.appendChild(projdim)

        # write projections
        for ii,spec in enumerate(specs):
          index = dataSources.index(spec)
          ll1 = self.defsMatrix[index]
          # NB all projections have factor 1 for the first(acquisition) dimension
          ll2 = ["+1"]
          for factor in ll1:
            if factor == 0:
              ll2.append(' 0')
            elif int(factor) == factor:
              ll2.append('%+d' % int(factor))
            else:
              ll2.append('%.3f' % factor)
          projection = doc.createElement("Projection")
          projection.setAttribute("p",str(ii))
          projection.setAttribute("name",spec.name)
          projection.setAttribute("factors",' '.join(ll2))
          projset.appendChild(projection)

      #Write out regions
      for index, interval in enumerate(intervals):
        intervalA, intervalB, nComp, intervalId = interval #@UnusedVariable
        region = doc.createElement("Region")
        region.setAttribute("r",str(index))
        region.setAttribute("ncomp",str(nComp))
        decomposition.appendChild(region)

      #Write out components
      cIndex = -1 # set component index to start at zero
      for index, interval in enumerate(intervals):

        intervalA, intervalB, nComp, intervalId = interval #@UnusedVariable
        output = self.prodecompOutput[intervalId]
        outfdir, outf, acqName, shapeNames = output

        for i in range(nComp):
          cIndex += 1
          component = doc.createElement("Component")
          decomposition.appendChild(component)
          component.setAttribute("c",str(cIndex))
          component.setAttribute("regionid",str(index))
          component.setAttribute("status",'raw')
          #
          #
          # Below is to be put into a higher-level container for all shapes
          # i.e. Decomposition data params
          #
          #component.setAttribute("data", defdat)

          # calculate the amplitude of a component
          ampl = max(outfdir[:,i])
          af = [max(outf[:,s,i]) for s in range(len(shapeNames))]
          for ca in range(len(af)):
            ampl = ampl*af[ca]
          component.setAttribute("ampl",str(ampl))

          # DIRECT shape
          shape = doc.createElement("Shape")
          shape.setAttribute("a",str(0))
          shape.setAttribute("size",str(outfdir.shape[0]))
          shape.setAttribute("offset",str(intervalA))
          for x in outfdir[:,i]:
            #shape.appendChild(doc.createTextNode("\n"))
            shape.appendChild(doc.createTextNode(str(x)))
          peaks = doc.createElement("Peaks")
          peaks.setAttribute("list","pk,posppm,intensity")
          #shape.appendChild(doc.createTextNode("\n"))
          shape.appendChild(peaks)
          #shape.appendChild(doc.createTextNode("\n"))

          #component.appendChild(doc.createTextNode("\n"))
          component.appendChild(shape)

          # INDIRECT shapes
          for j in range(len(shapeNames)):
            shape = doc.createElement("Shape")
            shape.setAttribute("a",str(j+1))
            for x in outf[:,j,i]:
              #shape.appendChild(doc.createTextNode("\n"))
              shape.appendChild(doc.createTextNode(str(x)))
            peaks = doc.createElement("Peaks")
            peaks.setAttribute("list","pk,posppm,intensity")
            #shape.appendChild(doc.createTextNode("\n"))
            shape.appendChild(peaks)
            #shape.appendChild(doc.createTextNode("\n"))

            #component.appendChild(doc.createTextNode("\n"))
            component.appendChild(shape)

          #component.appendChild(doc.createTextNode("\n"))
          #decomposition.appendChild(doc.createTextNode("\n"))

    #mwaydec.appendChild(doc.createTextNode("\n"))

    # NBNB now ovrwrites previous file every time
    #FL = open(fileName,'w+')
    FL = open(fileName,'w')
    doc.writexml(FL, addindent="  ", newl='\n')
    FL.close()

    showInfo("Info","Shapes are saved to \n"+fileName)


def getRefExperiment(experiments, acqName, shapeNames):
  """ Get refExperiment that describes an nD experiment for the
  result of a PRODECOMP decomposition starting from experiments with
  acquisition displayName acqName and shape names shapeNames in active use
  """

  nameSet = set(shapeNames)
  nameSet.add(acqName)

  aSet = set(xp.refExperiment for xp in experiments)
  if len(aSet) == 1:
    refExperiment = aSet.pop()
  else:
    refExperiment = None

  if refExperiment is None:
    return None

  else:
    expMeasurements = set()
    for xp in experiments:
      for xpDim in xp.sortedExpDims():
        for xpDimRef in xpDim.sortedExpDimRefs():

          refExpDimRef = xpDimRef.refExpDimRef
          if refExpDimRef is None:
            print 'WARNING, %s has no refExpDimRef' % xpDimRef
            return None

          elif xpDimRef.displayName in nameSet:
            expMeasurements.add(refExpDimRef.expMeasurement)


  isReversed = refExperiment.isReversed
  ll = []
  for rxp in refExperiment.nmrExpPrototype.refExperiments:
    if rxp.isReversed == isReversed:
      aSet = set()
      for rxd in rxp.refExpDims:
        for rxdr in rxd.refExpDimRefs:
          aSet.add(rxdr.expMeasurement)
      if aSet == expMeasurements:
        ll.append((len(rxp.refExpDims), rxp))
  if ll:
    ll.sort()
    result = ll[-1][1]
  else:
    result = None
  #
  return result

if __name__ == "__main__":

  import sys
  import Tkinter

  root = Tkinter.Tk()
  root.withdraw()

  if len(sys.argv) == 2:
    path = sys.argv[1]
    from ccp.gui.Io import loadProject
    ccpnProject = loadProject(root, path=path)
  else:
    ccpnProject = None

  popup = ProdecompPopup(root, ccpnProject=ccpnProject)

  root.mainloop()
