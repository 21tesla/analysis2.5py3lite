
"""
======================COPYRIGHT/LICENSE START==========================

CreateAxisUnit.py: Part of the CcpNmr Analysis program

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
from memops.gui.BooleanPulldownMenu import BooleanPulldownMenu
from memops.gui.Label import Label
from memops.gui.Entry import Entry
from memops.gui.MessageReporter import showError

from ccpnmr.analysis.popups.BasePopup import BasePopup

from ccpnmr.api import Analysis

class CreateAxisUnitPopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    BasePopup.__init__(self, parent=parent,
                       title='Create axis unit', modal=True, **kw)

  def body(self, master):

    master.grid_columnconfigure(1, weight=1)

    row = 0
    label = Label(master, text='Axis unit name: ', grid=(row, 0))
    tipText = 'Short textual name for the unit, e.g. "KPa" or "ms"'
    self.nameEntry = Entry(master, width=10, grid=(row, 1), tipText=tipText)

    row += 1
    label = Label(master, text='Unit is backwards: ', grid=(row, 0))
    tipText = 'Whether the axis values decrease left to right & bottom to top. For example "ppm" does, but most units do not'
    self.backwardsMenu = BooleanPulldownMenu(master, grid=(row, 1), tipText=tipText)

    row += 1
    tipTexts = ['Make a new unit specification using the stated options and close this popup']
    texts = [ 'Create' ]
    commands = [ self.ok ]
    buttons = UtilityButtonList(master, texts=texts, commands=commands,
                                closeText='Cancel', helpUrl=self.help_url, 
                                grid=(row, 0), gridSpan=(1,2), tipTexts=tipTexts)

    master.grid_rowconfigure(row, weight=1)

  def apply(self):

    unit = self.nameEntry.get()
    if (not unit):
      showError('No unit', 'Need to enter unit')
      return False

    units = [ axisUnit.unit for axisUnit in self.analysisProject.axisUnits ]
    if (unit in units):
      showError('Repeated unit', 'Unit already used')
      return False

    isBackwards = self.backwardsMenu.getSelected()

    self.analysisProject.newAxisUnit(unit=unit, isBackwards=isBackwards)

    return True

