
"""
======================COPYRIGHT/LICENSE START==========================

CreatePanelType.py: Part of the CcpNmr Analysis program

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

from memops.general import Implementation

from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Label import Label
from memops.gui.Entry import Entry
from memops.gui.MessageReporter import showError
from memops.gui.PulldownList import PulldownList
from memops.gui.Separator import Separator

from ccpnmr.analysis.popups.BasePopup import BasePopup

from ccpnmr.api import Analysis

class CreatePanelTypePopup(BasePopup):

  def __init__(self, parent, *args, **kw):

    self.axisType = None

    BasePopup.__init__(self, parent=parent,
                       title='Create panel type', modal=True, **kw)

  def body(self, master):

    master.grid_columnconfigure(1, weight=1)

    row = 0
    label = Label(master, text='Panel name: ', grid=(row,0))
    tipText = 'Short text name for the new axis panel, e.g. "N2"'
    self.name_entry = Entry(master, width=15, grid=(row,1), tipText=tipText)

    row += 1 
    label = Label(master, text='Axis type:', grid=(row,0))
    tipText = 'The type of axis (isotope, time, sampled etc.) represented by panel type'
    self.types_list = PulldownList(master, grid=(row,1), tipText=tipText)

    row += 1
    tipTexts = ['Create a new panel type object with the selected options & close the popup']
    texts = [ 'Create' ]
    commands = [ self.ok ]
    buttons = UtilityButtonList(master, texts=texts, commands=commands, doClone=False,
                                closeText='Cancel', helpUrl=self.help_url, grid=(row,0),
                                gridSpan=(1,2), tipTexts=tipTexts)

    master.grid_rowconfigure(row, weight=1)

    self.administerNotifiers(self.registerNotify)
    self.update()

  def administerNotifiers(self, notifyFunc):

    for func in ('__init__', 'delete', 'setName'):
      notifyFunc(self.update, 'ccpnmr.Analysis.AxisType', func)

  def destroy(self):

    self.administerNotifiers(self.unregisterNotify)

    BasePopup.destroy(self)

  def update(self, *extra):

    axisType = self.axisType
    axisTypes = self.parent.getAxisTypes()
    names = [x.name for x in axisTypes]
    if axisTypes:
      if axisType not in axisTypes:
        self.axisType = axisType = axisTypes[0]
      index = axisTypes.index(axisType)
    else:
      index = 0
      self.axisType = None

    self.types_list.setup(names, axisTypes, index)

  def apply(self):

    name = self.name_entry.get()
    if not name:
      showError('No name', 'Need to enter name', parent=self)
      return False

    names = [ panelType.name for panelType in self.analysisProject.panelTypes ]
    if name in names:
      showError('Repeated name', 'Name already used', parent=self)
      return False

    axisType = self.types_list.getObject()

    if not axisType:
      showError('No axis type', 'Need to create axis type', parent=self)
      return False
    
    self.analysisProject.newPanelType(name=name, axisType=axisType)

    return True

