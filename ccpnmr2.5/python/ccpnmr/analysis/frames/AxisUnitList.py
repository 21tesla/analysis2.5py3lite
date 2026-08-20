
"""
======================COPYRIGHT/LICENSE START==========================

AxisUnitList.py: Part of the CcpNmr Analysis program

Copyright (C) 2005 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

This file contains reserved and/or proprietary information
belonging to the author and/or organisation holding the copyright.
It may not be used, distributed, modified, transmitted, stored,
or in any way accessed, except by members or employees of the CCPN,
and by these people only until 31 December 2005 and in accordance with
the guidelines of the CCPN.
 
A copy of this license can be found in ../../../license/CCPN.license.

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
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.ScrolledListbox import ScrolledListbox

# do not need 'setUnit' since unit is AxisUnit key
notify_funcs = ('__init__', 'delete')

class AxisUnitList(PulldownMenu):
 
  def __init__(self, parent, getAxisUnits, *args, **kw):
 
    self.getAxisUnits = getAxisUnits
 
    PulldownMenu.__init__(self, parent, *args, **kw)
 
    for func in notify_funcs:
      Implementation.registerNotify(self.setAxisUnits, 'ccpnmr.Analysis.AxisUnit', func)

  def destroy(self):

    for func in notify_funcs:
      Implementation.unregisterNotify(self.setAxisUnits, 'ccpnmr.Analysis.AxisUnit', func)

    PulldownMenu.destroy(self)
 
  def setAxisUnits(self, *unit):
 
    axisUnits = self.getAxisUnits()
    names = [ axisUnit.unit for axisUnit in axisUnits ]
    self.replace(names, self.selected_index)

class AxisUnitScrolledList(ScrolledListbox):
 
  def __init__(self, parent, getAxisUnits):
 
    self.getAxisUnits = getAxisUnits
 
    ScrolledListbox.__init__(self, parent, xscroll=False, selectmode=Tkinter.MULTIPLE)
 
    for func in notify_funcs:
      Implementation.registerNotify(self.setAxisUnits, 'ccpnmr.Analysis.AxisUnit', func)

  def destroy(self):

    for func in notify_funcs:
      Implementation.unregisterNotify(self.setAxisUnits, 'ccpnmr.Analysis.AxisUnit', func)

    ScrolledListbox.destroy(self)
 
  def setAxisUnits(self, *unit):
 
    axisUnits = self.getAxisUnits()
    names = [ axisUnit.unit for axisUnit in axisUnits ]
    items = self.getSelectedItems()
    self.setItems(names)
    self.setSelectedItems(items)
