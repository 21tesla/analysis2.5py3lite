
"""
======================COPYRIGHT/LICENSE START==========================

NewWindow.py: Part of the CcpNmr Analysis program

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

from memops.gui.ButtonList import ButtonList, UtilityButtonList
from memops.gui.Entry import Entry
from memops.gui.ButtonList import ButtonList
from memops.gui.Label import Label
from memops.gui.LabelFrame import LabelFrame
from memops.gui.Frame import Frame
from memops.gui.LabelDivider import LabelDivider
from memops.gui.MessageReporter import showError
from memops.gui.PulldownList import PulldownList
from memops.gui.ScrolledMatrix import ScrolledMatrix

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core import WindowBasic
from ccpnmr.analysis.core import Util


# Below is not used
# DEFAULT_COLOR = 'white'

AXIS_LABELS = ['x', 'y'] + ['z%d' % (i+1) for i in range(4)]

STRIP_NUMS = range(1, 6)

class NewWindowPopup(BasePopup):
  """
  **Create New Windows to Display Spectra**
  
  This tool is used to make new windows for the graphical display of spectra,
  which will usually be as contours. It is notable that some spectrum windows
  will be made automatically when spectra are loaded if there is no existing
  appropriate window to display a spectrum. However, making new spectrum windows
  allows the user to specialise different windows for different tasks and gives
  complete freedom as to which types of axis go in which direction. For example
  the user may wish to make a new window so that a spectrum can be viewed from
  an orthogonal, rotated aspect.

  A new spectrum window is made by first considering whether it is similar to
  any existing windows. If so, then the user selects the appropriate template
  window to base the new one upon. The user then chooses a name for the window
  via the "New window name" field, although the name may be changed after the
  window is created. Usually the user does not need to consider the "Strips"
  section, but if required the new window can be created with starting strips
  and orthogonal sub-divisions (although these are not permanent). After setting
  the required axes and spectrum visibility, as described below, the user clicks
  on [Create Window!] to actually make the new spectrum display window.

  **Axes**

  The number and type of axes for the new window are chosen using the pulldown
  menus in the "Axes" section. The idea is that the user chooses which NMR
  isotopes should appear on the X, Y, & Z axes. Naturally, X and Y axes must
  always be set to something, to represent the plane of the screen, but the Z
  axes are optional. If not required, a Z axis may be set to "None" to indicate
  that it will not be used. Up to four Z axes may be specified, labelled as
  "z1", "z2" etc., and these represent extra dimensions orthogonal to the plane
  of the  screen, which are often conceptualised as depth axes.
  
  It should be noted that the Y axis type may be set to "value", which refers to
  a spectrum intensity axis, rather than an NMR isotope axis. Setting the Y axis
  to "value" and only having the X axis set to an isotope is used to create
  windows that can show 1D spectra. Such value axes can also be used for 2D and
  higher dimensionality spectra, to show the data as an intensity graph (a
  "slice") rather than as contours.

  **Spectra**

  The lower "Viewed Spectra" section lists all of the spectra within the project
  that may be shown by a window with the selected axes. All spectra with
  isotopes in their data dimensions that match the isotope types of the window
  axes can potentially be displayed. This system also allows for displayed 
  spectra to have fewer dimensions than the axes has windows, as long as at
  least the X and Y axes are present. For example a 2D H-N spectrum can be
  shown in a 3D H-N-H window, but not a 3D H-H-N.

  For spectra that have more than one data dimension of the same isotope, then
  which data dimension goes with which window axis is not always known to
  Analysis. Where there is ambiguity, this system will simply map the spectrum
  data dimensions in order to the next matching window axis. If this mapping
  turns out to be wrong, then it may be changed at any time via the main
  _Windows settings; toggling the "Dim. Mapping" of the "Spectrum & Peak List
  Mappings" tab.
  
  For the spectra listed in the lower table, which may be placed in the new
  window, the user has control over whether the spectra actually will appear.
  Firstly the user can change the "Visible?" column, either via a double-click
  or by using the appropriate lower buttons. By default spectra are set as not
  being visible in new windows, and the user toggles the ones that should be
  seen to "Yes". This basic spectrum visibility can readily be changed by the
  toggle buttons that appear in the "toolbar" at the top of the spectrum
  display, so the "Visible?" setting here is only about what initially appears. 

  The "In Toolbar?" setting of a spectrum is a way of allowing the user to state
  that a spectrum should never appear in the window, and not even allow it to be
  toggled on later via the toolbar at the top of the windows. This is a way of
  reducing clutter, and allows certain windows to be used for particular subsets
  of spectra. For example the user may wish to put the spectra for a temperature
  series in one window, but not in other windows used for resonance assignment
  where they would get in the way. The "In Toolbar" setting can be changed
  after a window has been made, but only via the main Windows_ settings popup.

  .. _Windows: EditWindowPopup.html

  """
  def __init__(self, parent, *args, **kw):

    self.visibleSpectra = parent.visibleSpectra
    self.toolbarSpectra = parent.toolbarSpectra
    self.waiting = False
    self.window = None

    BasePopup.__init__(self, parent=parent, title='Window : New Window', **kw)

  def body(self, guiFrame):

    guiFrame.grid_columnconfigure(2, weight=1)

    row = 0
    label = Label(guiFrame, text='Template window: ', grid=(row,0))
    tipText = 'Selects which window to use as the basis for making a new spectrum window; sets the axis types accordingly'
    self.window_list = PulldownList(guiFrame, grid=(row,1),
                                    callback=self.setAxisTypes, tipText=tipText)
    frame = LabelFrame(guiFrame, text='Strips',
                       grid=(row,2), gridSpan=(2,1))
    buttons = UtilityButtonList(guiFrame, doClone=False,
                                helpUrl=self.help_url, grid=(row,3))


    row += 1
    label = Label(guiFrame, text='New window name: ', grid=(row,0))
    tipText = 'A short name to identify the spectrum window, which will appear in the graphical interface'
    self.nameEntry = Entry(guiFrame, width=16,
                           grid=(row,1), tipText=tipText)

    
    row += 1

    label = Label(frame, text='Columns: ', grid=(0,0))
    tipText = 'The number of vertical strips/dividers to initially make in the spectrum window'
    self.cols_menu = PulldownList(frame, objects=STRIP_NUMS, grid=(0,1),
                                  texts=[str(x) for x in STRIP_NUMS],
                                  tipText=tipText)

    label = Label(frame, text='Rows: ', grid=(0,2))
    tipText = 'The number of horizontal strips/dividers to initially make in the spectrum window'
    self.rows_menu = PulldownList(frame, objects=STRIP_NUMS, grid=(0,3),
                                  texts=[str(x) for x in STRIP_NUMS],
                                  tipText=tipText)
    row += 1
    div = LabelDivider(guiFrame, text='Axes',
                       grid=(row,0), gridSpan=(1,4))

    row += 1
    self.axis_lists = {}
    frame = Frame(guiFrame, grid=(row,0), gridSpan=(1,4))   

    col = 0
    self.axisTypes = {}
    self.axisTypesIncludeNone = {}
    for label in AXIS_LABELS:
      self.axisTypes[label] = None
      w = Label(frame, text=' ' + label)
      w.grid(row=0, column=col, sticky='w')
      col += 1

      if label in ('x', 'y'):
        includeNone = False
        tipText = 'Sets the kind of measurement (typically ppm for a given isotope) that will be used along the window %s axis' % label
      else:
        includeNone = True
        tipText = 'Where required, sets the kind of measurement (typically ppm for a given isotope) that will be used along the window %s axis' % label
      self.axisTypesIncludeNone[label] = includeNone
        
      getAxisTypes = lambda label=label: self.getAxisTypes(label)
      callback = lambda axisType, label=label: self.changedAxisType(label, axisType)
      self.axis_lists[label] = PulldownList(frame, callback=callback, tipText=tipText)
      self.axis_lists[label].grid(row=0, column=col, sticky='w')
      col +=1

    frame.grid_columnconfigure(col, weight=1)

    row += 1
    div = LabelDivider(guiFrame, text='Viewed Spectra',
                       grid=(row,0), gridSpan=(1,4))

    row += 1
    guiFrame.grid_rowconfigure(row, weight=1)

    editWidgets      = [None,None,None,None]
    editGetCallbacks = [None,self.toggleVisible,
                        self.toggleToolbar,None]
    editSetCallbacks = [None,None,None,None]
    tipTexts = ['The "experiment:spectrum" name for the spectrum that may be viewed in the new window, given the axis selections',
                'Sets whether the spectrum contours will be visible in the new window',
                'Sets whether the spectrum appears at all in the window; if not in the toolbar it cannot be displayed',
                'The number of peak lists the spectrum contains']
    headingList = ['Spectrum','Visible?',
                   'In Toolbar?','Peak Lists']

    self.scrolledMatrix = ScrolledMatrix(guiFrame,
                          headingList=headingList,
			  editWidgets=editWidgets,
                          editGetCallbacks=editGetCallbacks,
                          editSetCallbacks=editSetCallbacks,
			  multiSelect=True, grid=(row,0),
                          gridSpan=(1,4), tipTexts=tipTexts)

    row += 1
    tipTexts = ['Creates a new spectrum window with the specified parameters',
                'Sets the contours of the selected spectra to be visible when the new window is made',
                'Sets the contours of the selected spectra to not be displayed when the new window is made',
                'Sets the selected spectra as absent from the window toolbar, and thus not displayable at all']
    texts    = ['Create Window!','Selected\nVisible',
                'Selected\nNot Visible','Selected\nNot In Toolbar']
    commands = [self.ok,
                self.setSelectedDisplay,
                self.setSelectedHide,
                self.setSelectedAbsent]
    buttonList = ButtonList(guiFrame, texts=texts, grid=(row,0),
                            commands=commands, gridSpan=(1,4),
                            tipTexts=tipTexts)
    buttonList.buttons[0].config(bg='#B0FFB0')

    self.updateAxisTypes()
    self.updateWindow()
    self.updateWindowName()

    self.administerNotifiers(self.registerNotify)
    self.updateAfter()

  def administerNotifiers(self, notifyFunc):
    
    self.registerNotify(self.updateWindowName, 'ccpnmr.Analysis.SpectrumWindow', '__init__')
    
    for clazz in ('ccp.nmr.Nmr.Experiment','ccp.nmr.Nmr.DataSource'):
      for func in ('__init__','delete','setName'):
         notifyFunc(self.updateAfter, clazz, func)

    for func in ('__init__','delete'):
      notifyFunc(self.updateAfter,'ccp.nmr.Nmr.PeakList', func)

    for func in ('__init__', 'delete', 'setName', 'addSpectrumWindowGroup',
                 'removeSpectrumWindowGroup', 'setSpectrumWindowGroups'):
      notifyFunc(self.updateWindow, 'ccpnmr.Analysis.SpectrumWindow', func)

    for func in ('addSpectrumWindow', 'removeSpectrumWindow', 'setSpectrumWindows'):
      notifyFunc(self.updateWindow, 'ccpnmr.Analysis.SpectrumWindowGroup', func)

    for func in ('__init__', 'delete', 'setName'):
      notifyFunc(self.updateAxisTypes, 'ccpnmr.Analysis.AxisType', func)

# Set visible contours, on commit, according to selection
#   Get hold of spectrumWindowView ASAP

  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)
    
    BasePopup.destroy(self)

  def getSpectra(self):
  
    spectrumIsotopes = {}
  
    spectra = []
    for experiment in self.nmrProject.sortedExperiments():
      name = experiment.name
      
      for spectrum in experiment.sortedDataSources():
        spectrumIsotopes[spectrum] = []
        spectra.append(['%s:%s' % (name, spectrum.name), spectrum])
        
        for dataDim in spectrum.dataDims:
          dimTypes = []
          
          if dataDim.className != 'SampledDataDim':
            for dataDimRef in dataDim.dataDimRefs:
              expDimRef = dataDimRef.expDimRef
              isotopes = set()
            
              for isotopeCode in expDimRef.isotopeCodes:
                 isotopes.add(isotopeCode)
          
              dimTypes.append((expDimRef.measurementType.lower(), isotopes))
          else:
            dimTypes.append('sampled')
          
          spectrumIsotopes[spectrum].append(dimTypes)
    
    axisIsotopes = {}
    for label in AXIS_LABELS:
      if label not in ('x', 'y'):
        if self.axis_lists[label].getSelectedIndex() == 0:
          continue
          
      axisType = self.axis_lists[label].getObject()
      axisIsotopes[label] = (axisType.measurementType.lower(), set(axisType.isotopeCodes))

    spectraSel = []
    axes = axisIsotopes.keys()
    axes.sort()
    for name, spectrum in spectra:
      dimIsotopes = spectrumIsotopes[spectrum]
      
      for label in axes:
        mType, selected = axisIsotopes[label]
        if label == 'y' and mType == 'none':
          continue  # value axis
        
        for i, dimTypes in enumerate(dimIsotopes):
          for dimType in dimTypes:
            if dimType == 'sampled':
              if label != 'z1':
                continue
              axisType = self.axis_lists[label].getObject()
              if axisType.name != 'sampled':
                continue
              dimIsotopes.pop(i)
              break
            else:
              measurementType, isotopes = dimType
              if (mType == measurementType) and (selected <= isotopes):
                dimIsotopes.pop(i)
                break
          else:
            continue
          break       

        else:
          if label in ('x', 'y'):
            break
          
      else:
        
        if not dimIsotopes:
           spectraSel.append([name, spectrum])
 
    return spectraSel

  def setSelectedAbsent(self):
  
    for spectrum in self.scrolledMatrix.currentObjects:
      self.visibleSpectra[spectrum] = False     
      self.toolbarSpectra[spectrum] = False

    self.updateAfter()

  def setSelectedDisplay(self):

    for spectrum in self.scrolledMatrix.currentObjects:
      self.visibleSpectra[spectrum] = True       
      self.toolbarSpectra[spectrum] = True

    self.updateAfter()

  def setSelectedHide(self):

    for spectrum in self.scrolledMatrix.currentObjects:
      self.visibleSpectra[spectrum] = False       
      self.toolbarSpectra[spectrum] = True
      
    self.updateAfter()

  def toggleToolbar(self, spectrum):
  
    boolean = not self.toolbarSpectra.get(spectrum, True)
    self.toolbarSpectra[spectrum] = boolean
    
    if boolean is False:
      self.visibleSpectra[spectrum] = False       
    
    self.updateAfter()

  def toggleVisible(self, spectrum):
  
    boolean = not self.visibleSpectra.get(spectrum, False)
    self.visibleSpectra[spectrum] = boolean
    
    if boolean:
      if not self.toolbarSpectra.get(spectrum, True):
        self.toolbarSpectra[spectrum]  = True
    
    self.updateAfter()

  def updateAfter(self, object=None):
  
    if self.waiting:
      return

    else:
      self.waiting = True
      self.after_idle(self.update)
    
  def update(self):

    for spectrum in self.visibleSpectra.keys():
      if spectrum.isDeleted:
        del self.visibleSpectra[spectrum]

    textMatrix = []
    objectList = []
    colorMatrix = []
    
    for name, spectrum in self.getSpectra():
      colours = [None,None,None,None]
    
      if self.visibleSpectra.get(spectrum):
        colours[0] = '#60F060'
        isVisible = 'Yes'
        self.visibleSpectra[spectrum] = True # do not need this but play safe in case above if changed
      else:
        isVisible = 'No'  
        self.visibleSpectra[spectrum] = False

      if self.toolbarSpectra.get(spectrum, True):
        inToolbar = 'Yes'
        self.toolbarSpectra[spectrum] = True
      else:
        colours[0] = '#600000'
        inToolbar = 'No'  
        self.toolbarSpectra[spectrum] = False

      datum = [name,
               isVisible,
               inToolbar,
               ','.join(['%d' % pl.serial for pl in spectrum.peakLists])]
              
      textMatrix.append(datum)
      objectList.append(spectrum)
      colorMatrix.append(colours)
	 
    self.scrolledMatrix.update(objectList=objectList,
                               textMatrix=textMatrix,
			       colorMatrix=colorMatrix)

    self.waiting = False

  def updateAxisTypes(self, *extra):

    for label in AXIS_LABELS:
      names = []
      objects = []
      includeNone = self.axisTypesIncludeNone[label]
      
      if includeNone:
        names.append('None')
        objects.append(None)
      axisType = self.axisTypes[label]
      axisTypes = self.getAxisTypes(label)
      objects.extend(axisTypes)
      
      if axisTypes:
        if axisType not in objects:
          axisType = objects[0]
        index = objects.index(axisType)
        names.extend([x.name for x in axisTypes])

      else:
        index = 0

      self.axis_lists[label].setup(names, objects, index)
      self.changedAxisType(label, axisType)

  def changedAxisType(self, label, axisType):

    if axisType is not self.axisTypes[label]:
      self.axisTypes[label] = axisType
    
    self.updateAfter()

  def updateWindow(self, *extra):

    window = self.window
    windows = self.parent.getWindows()
    if windows:
      if window not in windows:
        window = windows[0]
      index = windows.index(window)
      names = [x.name for x in windows]
    else:
      index = 0
      names = []

    self.window_list.setup(names, windows, index)
    self.setAxisTypes(window)

  def updateWindowName(self, *extra):

    self.nameEntry.set(WindowBasic.defaultWindowName(self.project))

  def getAxisTypes(self, label):

    axisTypes = self.parent.getAxisTypes()

    if label == 'z1':
      axisTypes = [ axisType for axisType in axisTypes if axisType.name == 'sampled' or not axisType.isSampled ]
    else:
      axisTypes = [ axisType for axisType in axisTypes if not axisType.isSampled ]
    
    if label != 'y':
      axisTypes = [ at for at in axisTypes if at.name != WindowBasic.VALUE_AXIS_NAME ]

    return axisTypes

  def setAxisTypes(self, window):

    project = self.project
    if not project:
      return

    if window is self.window:
      return

    self.window = window
    if window:
      windowPanes = window.sortedSpectrumWindowPanes()
      if windowPanes:
        # might be empty because notifier called before windowPanes set up
        windowPane = windowPanes[0]
    
        for label in AXIS_LABELS:
          axisPanel = windowPane.findFirstAxisPanel(label=label)
        
          if axisPanel and axisPanel.panelType \
              and axisPanel.panelType.axisType:
            self.axis_lists[label].setSelected(axisPanel.panelType.axisType)
        
          elif label not in ('x', 'y'):
            self.axis_lists[label].setIndex(0)
      
      self.updateAfter()

  def apply(self):

    project = self.project
    if not project:
      return False

    name = self.nameEntry.get().strip()
    if not name:
      showError('No name', 'Need to enter name', parent=self)
      return False

    names = [window.name for window in self.analysisProject.spectrumWindows ]

    if (name in names):
      showError('Repeated name', 'Name already used', parent=self)
      return False

    axisTypes = []
    for label in AXIS_LABELS:
      if label not in ('x', 'y'):
        if self.axis_lists[label].getSelectedIndex() == 0:
          continue
          
      axisType = self.axis_lists[label].getObject()
      axisTypes.append(axisType)

    ncols = self.cols_menu.getObject()
    nrows = self.rows_menu.getObject()

    window = WindowBasic.createSpectrumWindow(project, name, [axisTypes,],
                                              ncols=ncols, nrows=nrows)

    return True
