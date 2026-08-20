
"""
======================COPYRIGHT/LICENSE START==========================

EditFitGraph.py: Part of the CcpNmr Analysis program

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
from memops.gui.CheckButton import CheckButton
from memops.gui.Frame import Frame
from memops.gui.Label import Label
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledGraph import ScrolledGraph
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.MessageReporter import showYesNo

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.DataAnalysisBasic import getFitMethodInfo
from ccpnmr.analysis.core.MarkBasic import createPeakMark
from ccpnmr.analysis.core.WindowBasic import getActiveWindows, isSpectrumInWindowPane, getWindowPaneName

from math import  sqrt

METHOD_NAMES = ['<None>',] + [x[0] for x in getFitMethodInfo()]

class EditFitGraphPopup(BasePopup):
  """
  **Analyse Function Curve Fitting to Peak Series Data**
  
  This popup is used to display the fit of a curve, of the displayed equation
  type, to data that has been extracted for a group of spectrum peaks from an
  NMR series. Precisely which kind of data is being fitted depends on the parent
  tool that this popup window was launched from. For example, for the `Follow
  Intensity Changes`_ tool the displayed graph is of a time or frequency value
  on the "X" axis (e.g. T1) verses peak intensity. For the `Follow Shift
  Changes`_ system the plot is for the parameters of the experiment titration on
  the "X" axis (e.g concentration) verses chemical shift distance.

  The upper graph shows the peak data that was used; the blue points, and the
  points of the fitted theoretical curve; red points. The lower table lists the
  data points for all of the peaks in the group, to which the data relates.
  There will be one peak for each row of the table, and hence point on the "X"
  axis. For each data point in the table the corresponding peak from which the
  data was extracted may be located with the "Follow in window" option, using
  the stated spectrum window and making marker lines as desired. Alternatively,
  the peak may be viewed in a table with [Show Peak]

  In general operation this system is used to check how well the selected
  function curve fits the peak data. The data that is displayed comes from the
  analysis tool that launched this popup. Any outliers may be removed using
  [Remove Point] (only with good reason) or the peaks themselves may be
  corrected if something has gone awry. Given that the display is launched from
  the selection of a specific group of peaks from a popup like `Follow Shift
  Changes`_  or `Follow Intensity Changes`_ and there are usually several peak
  groups that need to be analysed, the [Previous Set] and [Next Set] buttons can
  be used to quickly jump to the next peak group and see the next curve in the
  analysis results.

  .. _`Follow Intensity Changes`: CalcRatesPopup.html
  .. _`Follow Shift Changes`: FollowShiftChangesPopup.html

  """
  def __init__(self, parent, dataFitting, getXYfunction, updateFunction,
               xLabel, yLabel, showObjectFunction=None,
               nextSetFunction=None, prevSetFunction=None, 
               graphTitle='', **kw):

    self.guiParent       = parent
    self.getXYfunction   = getXYfunction
    self.updateFunction  = updateFunction
    self.nextSetFunction = nextSetFunction
    self.prevSetFunction = prevSetFunction
    self.dataFitting     = dataFitting
    self.object = None
    self.xLabel = xLabel
    self.yLabel = yLabel
    self.x = []
    self.y = []
    self.yFit   = []
    self.params = ()
    self.chiSq  = 0
    self.method  = dataFitting.fitFunction
    self.waiting = False
    self.noiseLevel = dataFitting.noiseLevel
    self.graphTitle = graphTitle
    self.windowPane = None
    self.mark = None
    
    BasePopup.__init__(self, parent=parent, title='Fit Graph', **kw)
    parent.protocol("WM_DELETE_WINDOW", self.close)

  def body(self, guiFrame):

    guiFrame.grid_columnconfigure(0, weight = 1)
    row = 0

    self.scrolledGraph = ScrolledGraph(guiFrame, width=400, height=300,
                                       symbolSize=5, symbols=['square','circle'],
                                       dataColors=['#000080','#800000'],
                                       lineWidths=[0,1], grid=(row,0))

    #self.scrolledGraph.setZoom(0.7)

    row += 1
    frame = Frame(guiFrame, grid=(row,0), sticky='ew')
    
    label = Label(frame, text='Fitting Function:', grid=(0,0))
    
    tipText = 'Selects which form of function to fit to the experimental data'
    self.methodPulldown  =  PulldownList(frame, self.changeMethod,
                                         grid=(0,1), tipText=tipText)

    row += 1
    frame = Frame(guiFrame, grid=(row,0), sticky='ew')
    
    tipText = 'The effective equation of the final fitted graph, incorporating all parameters'
    self.equationLabel = Label(frame, text='Equation:',
                               grid=(0,0),tipText=tipText)
    
    tipText = 'The error in the fit of the selected parameterised function to the experimental data'
    self.errorLabel =    Label(frame, text='Fit Error:',
                               grid=(0,1), tipText=tipText)
    
    row += 1
    frame = Frame(guiFrame, grid=(row,0), sticky='ew')
    
    label = Label(frame, text='Include x origin?:', grid=(0,0))
    
    tipText = 'Whether to include the x=0 point in the drawing'
    self.xOriginSelect = CheckButton(frame, callback=self.draw,
                                     grid=(0,1), selected=False, tipText=tipText)
    
    label = Label(frame, text='Include y origin?:', grid=(0,2))
    
    tipText = 'Whether to include the y=0 point in the drawing'
    self.yOriginSelect = CheckButton(frame, callback=self.draw,
                                     grid=(0,3), selected=False, tipText=tipText)

    label = Label(frame, text='Include y error?:', grid=(0,4))
    
    tipText = 'Whether to include the y error bars in the drawing (if these exist)'
    self.yErrorSelect = CheckButton(frame, callback=self.draw,
                                    grid=(0,5), selected=False, tipText=tipText)

    row += 1
    frame = Frame(guiFrame, grid=(row,0), sticky='ew')
    
    label = Label(frame, text='Navigation Window:', grid=(0,0))
    
    tipText = 'Selects which spectrum window will be used for navigating to peak positions'
    self.windowPanePulldown = PulldownList(frame, self.changeWindow, 
                                           grid=(0,1), tipText=tipText)
    
    label = Label(frame, text='Follow in window?:', grid=(0,2))
    
    tipText = 'Whether to navigate to the position of the reference peak (for the group), in the selected window'
    self.followSelect = CheckButton(frame, callback=self.windowPaneNavigate,
                                    grid=(0,3), selected=False, tipText=tipText)
    
    label = Label(frame, text='Mark Ref Peak?:', grid=(0,4))
    
    tipText = 'Whether to put a multi-dimensional cross-mark through the reference peak position, so it can be identified in spectra'
    self.markSelect = CheckButton(frame, callback=None, tipText=tipText,
                                  grid=(0,5), selected=False)


    row += 1
    guiFrame.grid_rowconfigure(row, weight = 1)
    tipTexts = ['The number of the data point, in order of increasing X-axis value',
                'For each point, the value of the parameter which is varied in the NMR series, e.g. T1, temperature, concentration etc.',
                'For each point, the experimental value being fitted, e.g. peak intensity of chemical shift distance',
                'The value of the best-fit function at the X-axis location',
                'The difference between the experimental (Y-axis) value and the fitted value',
                'The error in the experimental (Y-axis) value']
    headingList = ['Point','x','y','Fitted y',u'\u0394','y error']
    self.scrolledMatrix = ScrolledMatrix(guiFrame, headingList=headingList,
                                         callback=self.selectObject,
                                         tipTexts=tipTexts, grid=(row,0))

    row += 1
    tipTexts = ['Remove the selected data point, optionally removing the underlying peak',
                'Show a table of spectrum peaks that correspond to the selected data point']
    texts = ['Remove Point','Show Peak']
    commands = [self.removePoint, self.showObject]
    
    if self.prevSetFunction:
      texts.append('Previous Set')
      tipTexts.append('Move to the previous set of fitted values; the next group of peaks, often corresponding to a different residue or resonance')
      commands.append(self.prevSet)
    
    if self.nextSetFunction:
      tipTexts.append('Move to the next set of fitted values; the next group of peaks, often corresponding to a different residue or resonance')
      texts .append('Next Set')
      commands.append(self.nextSet)
      
    bottomButtons = UtilityButtonList(guiFrame, texts=texts, commands=commands,
                                      helpUrl=self.help_url, tipTexts=tipTexts,
                                      grid=(row,0), doClone=False)
    self.removeButton = bottomButtons.buttons[0]

    for func in ('__init__', 'delete', 'setName'):
      self.registerNotify(self.updateWindows, 'ccpnmr.Analysis.SpectrumWindow', func)

    self.update()

  def destroy(self):
  
    for func in ('__init__', 'delete', 'setName'):
      self.unregisterNotify(self.updateWindows, 'ccpnmr.Analysis.SpectrumWindow', func)
  
    BasePopup.destroy(self)

  def nextSet(self):
  
    if self.nextSetFunction:
      self.nextSetFunction()
      self.windowPaneNavigate()
  
  def prevSet(self):

    if self.prevSetFunction:
      self.prevSetFunction()
      self.windowPaneNavigate()

  def windowPaneNavigate(self, event=None):
  
    if self.followSelect.get() and self.dataFitting:
      peaks = self.dataFitting.objects
      if self.windowPane and peaks:
        peak = peaks[int(len(peaks)/2)]
        
        windowFrame = self.windowPane.getWindowFrame()
        windowFrame.gotoPeak(peak)
          
        if self.markSelect.get():
          if self.mark and not self.mark.isDeleted:
            self.mark.delete()
            
          self.mark = createPeakMark(self.dataFitting.refObject)
          
 
  def getWindows(self):
  
    peaks = self.scrolledMatrix.objectList
    
    windowPanes = []
    if peaks:
      peak = peaks[0]
      project    = peak.root
      spectrum   = peak.peakList.dataSource
      tryWindows = getActiveWindows(project)
      for window in tryWindows:
        for windowPane in window.sortedSpectrumWindowPanes():
          if isSpectrumInWindowPane(windowPane, spectrum):
            windowPanes.append( windowPane )
    
    return windowPanes    
  
  def changeWindow(self,  windowPane):
    
    if windowPane is not self.windowPane:
      self.windowPane = windowPane
      self.windowPaneNavigate()
  
  def updateWindows(self, window=None):
  
    windowPanes = self.getWindows()
    index = -1
    names = [getWindowPaneName(wp) for wp in windowPanes]
    windowPane = self.windowPane
    
    if windowPanes:
      if windowPane not in windowPanes:
        windowPane = windowPanes[0]
      
      index = windowPanes.index(windowPane)
          
    if windowPane is not self.windowPane:
      self.windowPane = windowPane
      self.windowPaneNavigate()
    
    self.windowPanePulldown.setup(names, windowPanes, index)

  def showObject(self):
  
    
    if self.object and (self.object.className == 'Peak'):
      peaks = [self.object]
      self.guiParent.guiParent.viewPeaks(peaks)

  def close(self):
  
    BasePopup.close(self)


  def changeMethod(self, index):
  
    if self.dataFitting:
      self.method = index
      self.dataFitting.fitFunction = self.method

    self.updateAfter()
  
  def updateMethods(self, *opt):
  
    n = len(METHOD_NAMES)
    if self.method >= n:
      self.method = 1
    
    self.methodPulldown.setup(METHOD_NAMES, range(n), self.method)

  def selectObject(self, object, row, col):
   
    if object:
      self.object = object
      self.updateButtons()

  def removePoint(self):
  
    if self.object and (len(self.dataFitting.objects) > 2 ):
      self.dataFitting.objects.remove(self.object)
      if showYesNo('Query','Delete the corresponding %s?' % (self.object.className), parent=self):
        self.object.delete()
      self.updateAfter()
  
  def draw(self, *junk):
    
    title = '%s Function fit %s' % (self.graphTitle, METHOD_NAMES[self.method])
    dataSet1 = []
    dataSet2 = []
    useErr = self.yErrorSelect.isSelected()
    if not useErr:
      useErr = None
    for i in range(len(self.scrolledMatrix.textMatrix)):
      row = self.scrolledMatrix.textMatrix[i]
      x  = row[1]
      y  = row[2]
      y2 = row[3]
      err = useErr and row[5]
      dataSet1.append( [x,y, err] )
      dataSet2.append( [x,y2] )
          
    dataSet1.sort()
    dataSet2.sort()
    
    if dataSet1 and dataSet2:
      drawOriginX = self.xOriginSelect.isSelected()
      drawOriginY = self.yOriginSelect.isSelected()
      self.scrolledGraph.update(dataSets=[dataSet1,dataSet2],
                                xLabel=self.xLabel,
                                yLabel=self.yLabel,
                                title=title,
                                drawOriginX=drawOriginX,
                                drawOriginY=drawOriginY)
      self.scrolledGraph.draw()
      
  def updateAfter(self, *object):
    
    if self.waiting:
      return
    else:
      self.waiting = True
      self.after_idle(self.update)
    
  def update(self, dataFitting=None, xLabel=None, yLabel=None, graphTitle=None, force=False):
           
    if (not force) and dataFitting and (self.dataFitting is dataFitting):
      if (self.method, self.xLabel, self.yLabel, self.noiseLevel) == \
         (self.dataFitting.fitFunction, xLabel, yLabel, self.dataFitting.noiseLevel):
        self.waiting = False
        return
    
    if dataFitting:
      self.dataFitting = dataFitting
 
    self.xLabel     = xLabel     or self.xLabel
    self.yLabel     = yLabel     or self.yLabel
    self.graphTitle = graphTitle or self.graphTitle
   
    if not self.dataFitting:
      self.waiting = False
      return
    else:
      dataFitting = self.dataFitting
      dataFitting = self.getXYfunction(dataFitting)
      if dataFitting.fitFunction is not None:
        self.method = dataFitting.fitFunction
      self.noiseLevel = dataFitting.noiseLevel  or self.noiseLevel
    self.updateMethods()

    isFitted = dataFitting.fit()
    
    textMatrix  = []
    objectList  = []
    
    if isFitted and self.method:
      methodInfo = getFitMethodInfo()[self.method-1]
      textFormat, indices = methodInfo[2]
      textParams   = [dataFitting.parameters[i] or 0 for i in indices]
      equationText = textFormat % tuple(textParams)
      errorText    = '%4f' % dataFitting.fitError
        
      
      x    = dataFitting.dataX
      y    = dataFitting.dataY
      yFit = dataFitting.fittedY
      N    = len(x)
      err = hasattr(dataFitting, 'dataErr') and dataFitting.dataErr
      if not err:
        err = None
      
      # wb104: 10 Mar 2014: not sure why the below was here, it isn't needed
      #if self.method == 10:
      #  x = [sqrt(v) for v in x]

      data = [ (x[i],y[i],yFit[i],y[i]-yFit[i],dataFitting.objects[i], err and err[i]) for i in range(N) ]
      data.sort()
      for i in range(N):
        xi, yi, yFiti, deltai, object, erri = data[i]
        textMatrix.append( [i+1, xi, yi, yFiti, deltai, erri] )
        objectList.append(object)
        
    else:
      equationText = 'Equation: <No Fit>'
      errorText    = '<None>'
      x = dataFitting.dataX
      y = dataFitting.dataY

      for i,x in enumerate(x):
        textMatrix.append( [i+1, x, y[i], None, None, None] )
        objectList.append(dataFitting.objects[i])

    self.equationLabel.set('Equation: %s' % equationText)
    self.errorLabel.set('Fit Error: %s' % errorText)
    self.scrolledMatrix.update(textMatrix=textMatrix, objectList=objectList)

    self.updateWindows()
    self.updateButtons()

    self.draw()

    if self.updateFunction:
      self.updateFunction(dataFitting)

    self.waiting = False
  
  def updateButtons(self):
  
    if self.object:
      self.removeButton.enable()
    else:
      self.removeButton.disable()
