
"""
======================COPYRIGHT/LICENSE START==========================

CreateAxisType.py: Part of the CcpNmr Analysis program

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
from memops.gui.CheckButtons import CheckButtons
from memops.gui.BooleanPulldownMenu import BooleanPulldownMenu
from memops.gui.Label import Label
from memops.gui.Entry import Entry
from memops.gui.FloatEntry import FloatEntry
from memops.gui.IntEntry import IntEntry
from memops.gui.MessageReporter import showError
from memops.gui.PulldownList import PulldownList
from memops.gui.Separator import Separator

from ccpnmr.analysis.popups.BasePopup import BasePopup

from ccpnmr.api import Analysis

class CreateAxisTypePopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.measurementType = None

    BasePopup.__init__(self, parent=parent,
                       title='Create axis type', modal=True, **kw)

  def body(self, master):

    master.grid_columnconfigure(1, weight=1)

    row = 0
    label = Label(master, text='Axis name: ', grid=(row, 0))
    tipText = 'Short text name for new type of axis e.g. "17O"'
    self.name_entry = Entry(master, width=15, grid=(row, 1), tipText=tipText)

    row += 1
    label = Label(master, text='Axis region: ', grid=(row, 0))
    tipText = 'Comma separated values for the upper and lower bound of the axis allowed range of values'
    self.region_entry = FloatEntry(master, text=[0.0, 1.0], isArray=True,
                                   width=15, grid=(row, 1), tipText=tipText)

    row += 1
    label = Label(master, text='Measurement type:', grid=(row, 0))
    tipText = 'The physical entity that is being measured along the axis'
    self.measurement_list = PulldownList(master, tipText=tipText)
    self.measurement_list.grid(row=row, column=1, sticky='w')

    row += 1
    label = Label(master, text='Dimension is sampled: ', grid=(row, 0))
    tipText = 'Whether the axis is discretely sampled or a continuous range (albeit on a grid)'
    self.sampled_popup = BooleanPulldownMenu(master, grid=(row, 1), tipText=tipText)

    row += 1
    label = Label(master, text='Decimal places: ', grid=(row, 0))
    tipText = 'The number of decimal places that the axis values are rounded to for display purposes'
    self.decimals_entry = IntEntry(master, text=0, width=15,
                                   grid=(row, 1), tipText=tipText)

    row += 1
    label = Label(master, text='Peak size: ', grid=(row, 0))
    tipText = 'The relative scale for the peak symbol (i.e the "X" shape) size compared to other axes'
    self.peak_size_entry = FloatEntry(master, text=1.0, width=15,
                                      grid=(row, 1), tipText=tipText)

    row += 1
    label = Label(master, text='Allowed axis units:', grid=(row, 0))
    tipTexts = ['Units of measurement allowed for this kind of axis',]
    units = [au.unit for au in self.parent.getAxisUnits()]
    selected = [True] * len(units)
    self.units_list = CheckButtons(master, units, selected=selected,
                                   direction='vertical',
                                   grid=(row, 1), tipTexts=tipTexts)

    row += 1
    tipTexts = ['Make a new axis specification of the selected type and close this popup']
    texts = [ 'Create' ]
    commands = [ self.ok ]
    buttons = UtilityButtonList(master, texts=texts, commands=commands, doClone=False,
                                closeText='Cancel', helpUrl=self.help_url, grid=(row, 0),
                                gridSpan=(1,2), tipTexts=tipTexts)

    master.grid_rowconfigure(row, weight=1)

    self.update()

  def update(self, *extra):

    measurementType = self.measurementType
    measurementTypes = self.parent.getMeasurementTypes()
    if measurementTypes:
      if measurementType not in measurementTypes:
        self.measurementType = measurementType = measurementTypes[0]
      index = measurementTypes.index(measurementType)
    else:
      index = 0
      self.measurementType = None

    self.measurement_list.setup(measurementTypes, None, index)

  def apply(self):

    name = self.name_entry.get()
    if (not name):
      showError('No name', 'Need to enter name', self)
      return False

    names = [ axisType.name for axisType in self.analysisProject.axisTypes ]
    if (name in names):
      showError('Repeated name', 'Name already used', self)
      return False

    region = self.region_entry.get()
    if ((region is None) or (len(region) != 2)):
      showError('Region error', 'Region must be float array of length two', self)
      return False

    if (region[0] >= region[1]):
      showError('Region error', 'Region must have first number < second number', self)
      return False

    measurementType = self.measurement_list.getText()
    isSampled = self.sampled_popup.getSelected()

    numDecimals = self.decimals_entry.get()
    if ((numDecimals is None) or (numDecimals < 0)):
      showError('Decimals error', 'Number of decimal places must be >= 0', self)
      return False

    peakSize = self.peak_size_entry.get()
    if ((peakSize is None) or (peakSize <= 0)):
      showError('Peak size error', 'Peak size must be > 0', self)
      return False

    selected = self.units_list.getSelected()
    allUnits = self.parent.getAxisUnits()
    axisUnits = [au for au in allUnits if au.unit in selected]
    
    self.analysisProject.newAxisType(name=name, region=region,
                      isSampled=isSampled, axisUnits=axisUnits,
                      numDecimals=numDecimals, peakSize=peakSize,
                      measurementType=measurementType)
 
    return True
