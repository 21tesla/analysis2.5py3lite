
"""
======================COPYRIGHT/LICENSE START==========================

GetCenterCoords.py: Part of the CcpNmr Analysis program

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

from memops.gui.ButtonList import UtilityButtonList
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError
from memops.gui.PulldownList import PulldownList

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.WindowBasic import getWindowPaneName
from ccpnmr.analysis.core import ExperimentBasic
from ccpnmr.analysis.core import UnitConverter

class GetCenterCoordsPopup(BasePopup):

  def __init__(self, parent):

    BasePopup.__init__(self, parent=parent, title='Center coordinates', modal=True)

  def body(self, master):

    self.center = None
    self.region_number = (0, 0)
    self.previousView = None

    windowPane = self.parent.activeWindowPane
    xregions = windowPane.findFirstAxisPanel(label='x').sortedAxisRegions()
    yregions = windowPane.findFirstAxisPanel(label='y').sortedAxisRegions()

    master.grid_columnconfigure(2, weight=1)
    master.grid_rowconfigure(0, weight=1)

    row = 0
    label = Label(master, text='Enter new xy center for window %s' % getWindowPaneName(windowPane))
    label.grid(row=row, column=0, columnspan=3, sticky='w')
    row = row + 1

    tipText = 'Sets the location for a new window centre position along the horizontal axis, using the below stated units'
    self.x_center_entry = FloatEntry(master, tipText=tipText)

    n = len(xregions)
    if (n > 1):
      label = Label(master, text='x strip number:')
      label.grid(row=row, column=0, sticky='w')
      tipText = 'Selects which strip, parallel to X axis, will be centered at the stated coordinates'
      indices = range(n)
      texts = ['%d' % (x+1) for x in indices]
      self.x_region_menu = PulldownList(master, objects=indices, texts=texts,
                                        callback=self.setXCenter, tipText=tipText)
      self.x_region_menu.grid(row=row, column=1, sticky='ew')
      row = row + 1
    else:
      self.setXCenter(0)

    label = Label(master, text='x center position:')
    label.grid(row=row, column=0, sticky='w')
    self.x_center_entry.grid(row=row, column=1, sticky='ew')
    row = row + 1

    tipText = 'Sets the location for a new window centre position along the vertical axis, using the below stated units'
    self.y_center_entry = FloatEntry(master, tipText=tipText)

    n = len(yregions)
    if (n > 1):
      label = Label(master, text='y strip number:')
      label.grid(row=row, column=0, sticky='w')
      tipText = 'Selects which strip, parallel to Y axis, will be centered at the stated coordinates'
      indices = range(n)
      texts = ['%d' % (x+1) for x in indices]
      self.y_region_menu = PulldownList(master, objects=indices, texts=texts,
                                        callback=self.setYCenter, tipText=tipText)
      self.y_region_menu.grid(row=row, column=1, sticky='ew')
      row = row + 1
    else:
      self.setYCenter(0)

    label = Label(master, text='y center position:')
    label.grid(row=row, column=0, sticky='w')
    self.y_center_entry.grid(row=row, column=1, sticky='ew')
    row = row + 1

    label = Label(master, text='Display unit:')
    label.grid(row=row, column=0, sticky='w')
    unit = self.parent.getPositionUnit()
    units = self.parent.getPositionUnits()
    
    tipText = 'Whether to display center coordinates and cursor location (top toolbar) in ppm or Hz units; using the stated spectrum for Hz referencing'
    self.unit_list = PulldownList(master, texts=units,
                                  callback=self.changedUnit, tipText=tipText)
    self.unit_list.set(unit)
    self.unit_list.grid(row=row, column=1, sticky='w')
    row = row + 1

    self.view_label = Label(master, text='  using spectrum:')
    views = self.parent.getPositionViews()
    texts = []
    for view in views:
      spectrum = view.analysisSpectrum.dataSource
      expt = spectrum.experiment
      text = '%s:%s' % (expt.name, spectrum.name)
      texts.append(text)

    tipText = 'Selects which spectrum to get spectrometer frequency from when using the Hz unit display'
    self.view_list = PulldownList(master, texts=texts, objects=views,
                                  tipText=tipText, callback=self.changedSpectrum)
    self.view_row = row
    row = row + 1

    texts = commands = []
    tipTexts = ['Move the spectrum window display so that its center, or that of the stated strip, is at the specified position',]
    texts = [ ' Commit ',]
    commands = [ self.ok,]
    buttons = UtilityButtonList(master, texts=texts, commands=commands, doClone=False,
                                helpUrl=self.help_url, tipTexts=tipTexts)
    buttons.grid(row=row, column=0, columnspan=2, sticky='ew')

    self.changedUnit(unit)

  def changedUnit(self, unit):

    if unit == 'Hz':
      self.previousView = None
      view = self.parent.getPositionView()
      if view:
        self.view_list.set(view)
        self.changedSpectrum(view)
      self.view_label.grid(row=self.view_row, column=0, sticky='ew')
      self.view_list.grid(row=self.view_row, column=1, sticky='ew')
    else:
      self.view_label.grid_forget()
      self.view_list.grid_forget()
      view = self.previousView
      if view:
        cx = self.x_center_entry.get()
        cy = self.y_center_entry.get()

        for (c, label, entry) in ((cx, 'x', self.x_center_entry), (cy, 'y', self.y_center_entry)):
          if c is not None:
            # first convert to points
            dataDimRef = self.getDataDimRef(view, label)
            if dataDimRef:
              c = UnitConverter.hz2pnt(c, dataDimRef)
              # now convert to new ppm
              c = UnitConverter.pnt2ppm(c, dataDimRef)
              entry.set(c)
      self.previousView = None

  def changedSpectrum(self, view):

    cx = self.x_center_entry.get()
    cy = self.y_center_entry.get()

    for (c, label, entry) in ((cx, 'x', self.x_center_entry), (cy, 'y', self.y_center_entry)):
      if c is not None:
        # first convert to ppm using old view
        if self.previousView:
          # c is in Hz
          oldDataDimRef = self.getDataDimRef(self.previousView, label)
          if oldDataDimRef:
            c = UnitConverter.hz2pnt(c, oldDataDimRef)
            c = UnitConverter.pnt2ppm(c, oldDataDimRef)
        # now convert to new Hz
        dataDimRef = self.getDataDimRef(view, label)
        if dataDimRef:
          c = UnitConverter.ppm2pnt(c, dataDimRef)
          c = UnitConverter.pnt2hz(c, dataDimRef)
          entry.set(c)

    self.previousView = view

  def setXCenter(self, selected):

    self.setCenter(selected, 0, self.x_center_entry, 'x')

  def setYCenter(self, selected):

    self.setCenter(selected, 1, self.y_center_entry, 'y')

  def setCenter(self, selected, n, entry, label):

    axisPanel = self.parent.activeWindowPane.sortedAxisPanels()[n]
    region = axisPanel.sortedAxisRegions()[selected].region
    c = 0.5 * (region[0] + region[1])
    if hasattr(self, 'unit_list'): # first time called this is not yet set up
      unit = self.unit_list.getText()
      if unit == 'Hz':
        # TBD: c always in ppm, so have to convert here...
        view = self.view_list.getObject()
        dataDimRef = self.getDataDimRef(view, label)
        if dataDimRef:
          c = UnitConverter.ppm2pnt(c, dataDimRef)
          c = UnitConverter.pnt2hz(c, dataDimRef)
    entry.set(c)

  def getDataDimRef(self, view, label):

    axisMapping = view.findFirstAxisMapping(label=label)
    if not axisMapping:
      return None

    dataDim = axisMapping.analysisDataDim.dataDim
    dataDimRef = ExperimentBasic.getPrimaryDataDimRef(dataDim)

    return dataDimRef

  def apply(self):

    cx = self.x_center_entry.get()
    cy = self.y_center_entry.get()

    if (cx is None or cy is None):
      showError('Blank entry', 'Entries must be floating point numbers.', parent=self)
      return False

    cs = [cx, cy]
    unit = self.unit_list.getText()
    if unit == 'Hz':
      for (c, label, n) in ((cx, 'x', 0), (cy, 'y', 1)):
        # TBD: need c always in ppm, so have to convert here...
        view = self.view_list.getObject()
        dataDimRef = self.getDataDimRef(view, label)
        if dataDimRef:
          c = UnitConverter.hz2pnt(c, dataDimRef)
          c = UnitConverter.pnt2ppm(c, dataDimRef)
          cs[n] = c
    self.center = tuple(cs)

    windowPane = self.parent.activeWindowPane

    xregions = windowPane.findFirstAxisPanel(label='x').sortedAxisRegions()
    if (len(xregions) == 1):
      rx = 0
    else:
      rx = self.x_region_menu.getObject()

    yregions = windowPane.findFirstAxisPanel(label='y').sortedAxisRegions()
    if (len(yregions) == 1):
      ry = 0
    else:
      ry = self.y_region_menu.getObject()

    self.region_number = (rx, ry)

    unit = self.unit_list.getText()
    view = self.view_list.getObject()
    self.parent.setPositionUnit(unit)
    self.parent.setPositionView(view)

    return True
