
"""
======================COPYRIGHT/LICENSE START==========================

EditAxisPanel.py: Part of the CcpNmr Analysis program

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
from memops.gui.FloatEntry import FloatEntry
from memops.gui.IntEntry import IntEntry
from memops.gui.ItemSelectPopup import ItemSelectPopup
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError, showYesNo
from memops.gui.ScrolledMatrix import ScrolledMatrix
from memops.gui.TabbedFrame import TabbedFrame

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.popups.CreateAxisUnit import CreateAxisUnitPopup
from ccpnmr.analysis.popups.CreateAxisType import CreateAxisTypePopup
from ccpnmr.analysis.popups.CreatePanelType import CreatePanelTypePopup
from ccpnmr.analysis.core.Util import default_axis_units_table, default_axis_types_table

class EditAxisPanelPopup(BasePopup):
  """
  **Make And Edit Spectrum Window Axis Types**
  
  The purpose of this popup window is to allow the user to create new types of
  axis for the display of spectra within windows. This requirement is relatively
  rare and under most circumstances the standard, default axis descriptions
  will not need to be changed. The popup is divided into three tabs which display different aspects of the
  spectrum display axis definitions.
  
  The first "Axis Units" tab lists the different units of measurement that may
  be used in an axis. Although new kinds of axis may be created, only the ppm,
  Hz and "arbitrary" units are currently used by Analysis, i.e. spectrum windows
  only display standard NMR isotope axes and intensity for a "value" Y axis.

  The "Axis Types" tab allows the user to specify new kinds of NMR data axis,
  for example to be able to work with spectra that have an isotope type that is
  not part of the normal Analysis setup. Clicking [Create] will prompt the user
  to specify the characteristics of the axis, and some of these can be adjusted
  after the initial definition. Under normal circumstances the only features
  that the user is likely to edit are the "Region", "Decimal places" and "Peak
  Size" columns. The region describes the upper and lower bounds of the axis
  vales, and dictates the maximum extent of a spectrum window display when fully
  zoomed out; a spectrum window cannot display outside the bounds of its axes,
  however these bounds may be changed.

  The "Panel Types" tab allows the user to create sub-categories within a given
  kind of axis. These are not used for very much, but are used with the 1D ruler
  lines, so that a line is only drawn on axes that share the same panel type.
  Accordingly, a ruler line drawn on an amide 1H axis (panel type "H1") will not
  show on the indirect 1H axis (panel type "H2"). The user can create extra
  panel types if this kind of distinction is required, although "H1" and "H2"
  types are present by default. The panel type of a spectrum window axis is set
  via the "Windows & Axes" section of the main Windows_ popup.

  .. _Windows: EditWindowPopup.html

  """
 
  def __init__(self, parent, *args, **kw):
 
    BasePopup.__init__(self, parent=parent, title='Window : Axes', **kw)

  def body(self, guiFrame):
 
    self.geometry('600x350')
 
    guiFrame.expandGrid(0,0)
   
    options = ['Axis Units','Axis Types','Panel Types']
    tabbedFrame = TabbedFrame(guiFrame, options=options, grid=(0,0))
    frameA, frameB, frameC = tabbedFrame.frames
    self.tabbedFrame = tabbedFrame
    
    #
    # Units
    #
    
    frameA.expandGrid(0,0)
    
    tipTexts = ['Row number',
                'Short text name for the unit of measurement, for graphical display',
                'Whether the axis values decrease left to right & bottom to top. For example "ppm" does, but most other units do not']
    headings = ('#', 'Name', 'Is backwards?')
    self.axis_unit_table = ScrolledMatrix(frameA, headingList=headings,
                                          callback=self.selectAxisUnit,
                                          deleteFunc=self.deleteAxisUnit,
                                          tipTexts=tipTexts, grid=(0,0))

    tipTexts = ['Create a new specification of a unit of measurement',
                'Delete the selected measurement unit']
    texts = [ 'Create', 'Delete' ]
    commands = [ self.createAxisUnit, self.deleteAxisUnit ]
    self.axis_unit_buttons = ButtonList(frameA, texts=texts, tipTexts=tipTexts,
                                        grid=(1,0), commands=commands)
    
    #
    # Types
    #
    
    frameB.expandGrid(0,0)

    tipTexts = ['Row number',
                'Name of window axis type, for graphical interface etc.',
                'Which isotopes the axis definition covers',
                'What kind of physical property is measured along the window axis',
                'Whether the axis represents discretely sampled values or a continuum of values (albeit fixed to a data grid)',
                'The upper and lower bounds for numerical values allowed on the axis',
                'The number of decimal places used to round axis values in graphical displays',
                'The relative scale for the peak symbol (i.e the "X" shape) size compared to other axes',
                'Units of measurement allowed for this kind of axis']
                
    headingList = ('#', 'Name', 'Isotope\ncodes', 'Measurement\nType',
                   'Dim\nSampled?', 'Region', 'Decimal\nplaces',
                   'Peak\nSize', 'Allowed\nUnits')
                
    self.regionEntry = FloatEntry(self, isArray=True, returnCallback=self.setRegion, width=12)
    self.decimalEntry = IntEntry(self, returnCallback=self.setDecimal, width=5)
    self.peakSizeEntry = FloatEntry(self, returnCallback=self.setPeakSize, width=5)
    
    editWidgets = [ None, None, None, None, None, self.regionEntry,
                    self.decimalEntry, self.peakSizeEntry, None ]
    editGetCallbacks = [ None, None, None, None, None, self.getRegion,
                         self.getDecimal, self.getPeakSize, self.addUnit ]
    editSetCallbacks = [ None, None, None, None, None, self.setRegion,
                         self.setDecimal, self.setPeakSize, None ]
                         
    self.axisTypeMatrix = ScrolledMatrix(frameB, tipTexts=tipTexts,
                                         headingList=headingList,
                                         initialRows=5, grid=(0,0),
                                         callback=self.selectAxisType,
                                         editWidgets=editWidgets,
                                         editGetCallbacks=editGetCallbacks,
                                         editSetCallbacks=editSetCallbacks,
                                         deleteFunc=self.deleteAxisType)
                                         
    tipTexts = ['Create a new kind of axis to use in spectrum windows',
                'Delete the selected window axis type specification']
    texts = [ 'Create', 'Delete' ]
    commands = [ self.createAxisType, self.deleteAxisType ]
    self.axis_type_buttons = ButtonList(frameB, texts=texts, tipTexts=tipTexts,
                                        commands=commands, grid=(1,0))
    
    #
    # Types
    #
    
    frameC.expandGrid(0,0)

    tipTexts = ['Row number',
                'Name of panel type specification',
                'Which kind of axis the panel is a subtype of']
    headings = ('#', 'Name', 'AxisType')
    self.panelTypeMatrix = ScrolledMatrix(frameC, tipTexts=tipTexts,
                                           headingList=headings,
                                           initialRows=5,
                                           callback=self.selectPanelType,
                                           deleteFunc=self.deletePanelType)
                                           
    self.panelTypeMatrix.grid(row=0, column=0, sticky='nsew')

    tipTexts = ['Add a new panel type specification (a subtype of a given kind of axis)',
                'Delete the selected panel type specification']
    texts = [ 'Create', 'Delete' ]
    commands = [ self.createPanelType, self.deletePanelType ]
    self.panel_type_buttons = ButtonList(frameC, texts=texts, grid=(1,0),
                                         commands=commands, tipTexts=tipTexts)

    #
    # Main
    #

    buttons = UtilityButtonList(tabbedFrame.sideFrame,
                                helpUrl=self.help_url,
                                grid=(0,0), sticky='e')


    self.updateAxisUnitTable()
    self.updateAxisTypeTable()
    self.updatePanelTypeTable()
    self.selectAxisUnit()
    self.selectAxisType()
    self.selectPanelType()
    self.administerNotifiers(self.registerNotify)

  def administerNotifiers(self, notifyFunc):
  
    for func in ('__init__', 'delete', ''):
      notifyFunc(self.updateAxisUnitTable, 'ccpnmr.Analysis.AxisUnit', func)
      notifyFunc(self.updateAxisTypeTable, 'ccpnmr.Analysis.AxisType', func)
      notifyFunc(self.updatePanelTypeTable, 'ccpnmr.Analysis.PanelType', func)  


  def open(self):
  
    self.updateAxisUnitTable()
    self.updateAxisTypeTable()
    self.updatePanelTypeTable()
    self.selectAxisUnit()
    self.selectAxisType()
    self.selectPanelType()
    BasePopup.open(self)
    
    
  def destroy(self):
 
    self.administerNotifiers(self.unregisterNotify)
    BasePopup.destroy(self)


  def getAxisPanels(self):

    axisPanels = []

    for window in self.analysisProject.sortedSpectrumWindows():
      for windowPane in window.sortedSpectrumWindowPanes():
        axisPanels.extend(windowPane.sortedAxisPanels())

    return axisPanels

  def updateAxisUnitTable(self, *extra):

    axisUnits = self.parent.getAxisUnits()

    textMatrix = []
    n = 0
    for axisUnit in axisUnits:
      n = n + 1
     
      if (axisUnit.isBackwards):
        isBackwards = 'Yes'
      else:
        isBackwards = 'No'

      textMatrix.append([n, axisUnit.unit, isBackwards])
 
    self.axis_unit_table.update(objectList=axisUnits, textMatrix=textMatrix)

  def selectAxisUnit(self, axisUnit = None, *extra):

    axisUnits = []
    for axisType in self.analysisProject.axisTypes:
      axisUnits.extend(axisType.axisUnits)

    if (not axisUnit or \
        axisUnit in axisUnits or \
        axisUnit in [ axisPanel.axisUnit for axisPanel in self.getAxisPanels() ] or \
        axisUnit.unit in [ entry[0] for entry in default_axis_units_table ]):
      state = 'disabled'
    else:
      state = 'normal'
 
    self.axis_unit_buttons.buttons[1].config(state=state)

  def createAxisUnit(self):

    popup = CreateAxisUnitPopup(self, self.project)
    popup.destroy()

  def deleteAxisUnit(self, *event):

    axisUnit = self.getAxisUnit()
    if (axisUnit):
      if (showYesNo('Delete axis unit',
          'Are you sure you want to delete axis unit "%s"?' % axisUnit.unit, self)):
        axisUnit.delete()

  def getAxisUnit(self):

    return self.axis_unit_table.currentObject

  def updateAxisTypeTable(self, *extra):

    axisTypes = self.parent.getAxisTypes()

    textMatrix = []
    n = 0
    for axisType in axisTypes:
      n = n + 1
      if (axisType.isSampled):
        isSampled = 'Yes'
      else:
        isSampled = 'No'
        
      text = [n,
              axisType.name,
              ', '.join(axisType.isotopeCodes),
              axisType.measurementType,
              isSampled,
              '%2.1f, %2.1f' % axisType.region,
              axisType.numDecimals,
              axisType.peakSize,
              ', '.join([axisUnit.unit for axisUnit in axisType.axisUnits])]
 
      textMatrix.append(text)
 
    self.axisTypeMatrix.update(objectList=axisTypes, textMatrix=textMatrix)

  def selectAxisType(self, axisType = None, *extra):

    markDims = []
    for mark in self.analysisProject.marks:
      markDims.extend(mark.markDims)

    if (not axisType or \
        axisType in [ axisPanel.axisType for axisPanel in self.getAxisPanels() ] or \
        axisType in [ markDim.axisType for markDim in markDims ] or \
        axisType.name in [ entry[0] for entry in default_axis_types_table ]):
      state = 'disabled'
    else:
      state = 'normal'
 
    self.axis_type_buttons.buttons[1].config(state=state)

  def getRegion(self, axisType):

    self.regionEntry.set(axisType.region)

  def setRegion(self, *extra):

    axisType = self.getAxisType()
    try:
      axisType.region = self.regionEntry.get()
    except Implementation.ApiError, e:
      showError('Setting region', e.error_msg, parent=self)

  def getDecimal(self, axisType):

    self.decimalEntry.set(axisType.numDecimals)

  def setDecimal(self, *extra):

    axisType = self.getAxisType()
    try:
      axisType.numDecimals = self.decimalEntry.get()
    except Implementation.ApiError, e:
      showError('Setting number of decimals', e.error_msg, parent=self)

  def getPeakSize(self, axisType):

    self.peakSizeEntry.set(axisType.peakSize)

  def setPeakSize(self, *extra):

    axisType = self.getAxisType()
    try:
      axisType.peakSize = self.peakSizeEntry.get()
    except Implementation.ApiError, e:
      showError('Setting peak size', e.error_msg, parent=self)

  def addUnit(self, axisType):

    axisUnits = self.analysisProject.sortedAxisUnits()
    
    for axisUnit in axisType.axisUnits:
      axisUnits.remove(axisUnit)
    popup = ItemSelectPopup(self,
              label = 'Axis unit: ',
              entries=[axisUnit.unit for axisUnit in axisUnits],
              message='Select axis unit to add to axis type "%s":' % axisType.name,
              select_text='Add')
    name = popup.item
    popup.destroy()
    if name:
      axisUnit = self.analysisProject.findFirstAxisUnit(unit=name)
      axisType.addAxisUnit(axisUnit)

  # needed by CreateAxisTypePopup
  def getAxisUnits(self, *args, **kw):

    return self.parent.getAxisUnits(*args, **kw)

  # needed by CreateAxisTypePopup
  def getMeasurementTypes(self, *args, **kw):

    return self.parent.getMeasurementTypes(*args, **kw)

  def createAxisType(self):

    popup = CreateAxisTypePopup(self, self.project)
    popup.destroy()

  def deleteAxisType(self, *event):

    axisType = self.getAxisType()
    if axisType:
      if (showYesNo('Delete axis type',
          'Are you sure you want to delete axis type "%s"?' % axisType.name, self)):
        axisType.delete()

  def getAxisType(self):

    # when measurementMenu created calls this before self.axisTypeMatrix
    # exists so protect against this
    try:
      return self.axisTypeMatrix.currentObject
    except:
      return None

  def updatePanelTypeTable(self, *extra):

    panelTypes = self.parent.getPanelTypes()

    textMatrix = []
    n = 0
    for panelType in panelTypes:
      n = n + 1
      text = []
      text.append(n)
      text.append(panelType.name)
      text.append(panelType.axisType.name)
      textMatrix.append(text)
 
    self.panelTypeMatrix.update(objectList=panelTypes, textMatrix=textMatrix)

  def selectPanelType(self, panelType = None, *extra):

    if (not panelType or \
        panelType in [ axisPanel.panelType for axisPanel in self.getAxisPanels() ] or \
        panelType in [ ruler.panelType for ruler in self.analysisProject.rulers ]):
      state = 'disabled'
    else:
      state = 'normal'
 
    self.panel_type_buttons.buttons[1].config(state=state)

  # needed by CreatePanelTypePopup
  def getAxisTypes(self, *args, **kw):

    return self.parent.getAxisTypes(*args, **kw)

  def createPanelType(self):

    popup = CreatePanelTypePopup(self, self.project)
    popup.destroy()

  def deletePanelType(self, *event):

    panelType = self.getPanelType()
    if (panelType):
      if (showYesNo('Delete panel type',
          'Are you sure you want to delete panel type "%s"?' % panelType.name, self)):
        panelType.delete()

  def getPanelType(self):

    return self.panelTypeMatrix.currentObject
