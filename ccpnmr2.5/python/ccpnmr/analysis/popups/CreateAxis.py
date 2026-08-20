
"""
======================COPYRIGHT/LICENSE START==========================

CreateAxis.py: Part of the CcpNmr Analysis program

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
from memops.gui.Label import Label
from memops.gui.PulldownList import PulldownList

from ccpnmr.analysis.popups.BasePopup import BasePopup

class CreateAxisPopup(BasePopup):

  def __init__(self, parent, **kw):

    self.axisType = None

    BasePopup.__init__(self, parent=parent, title='Create window axis',
                       modal=True, transient=False, **kw)

  def body(self, master):

    row = 0
    label = Label(master, text='Axis type: ', grid=(row, 0))
    tipText = 'Selects what type of measurement is displayed along the window axis'
    self.type_list = PulldownList(master, tipText=tipText, grid=(row, 1))

    row += 1
    tipTexts = ['Make an axis of the selected type in the window & close this popup']
    texts = [ 'Create' ]
    commands = [ self.ok ]
    buttons = UtilityButtonList(master, texts=texts, doClone=False, grid=(row, 0),
                                commands=commands, helpUrl=self.help_url,
                                gridSpan=(1,2), tipTexts=tipTexts)

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
    axisTypes = self.parent.parent.getAxisTypes()
    names = [x.name for x in axisTypes]
    if axisTypes:
      if axisType not in axisTypes:
        self.axisType = axisType = axisTypes[0]
      index = axisTypes.index(axisType)     
    else:
      index = 0
      self.axisType = None

    self.type_list.setup(names, axisTypes, index)

  def apply(self):

    self.axisType = self.type_list.getObject()

    return True
