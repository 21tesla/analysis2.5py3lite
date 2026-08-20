
"""
======================COPYRIGHT/LICENSE START==========================

ColorSchemeList.py: Part of the CcpNmr Analysis program

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
from memops.general import Implementation
from memops.gui.PulldownMenu import PulldownMenu
from memops.gui.Color import hexRepr

# do not need for 'setName' since name is ColorScheme key
notify_funcs = ('__init__', 'delete')

class ColorSchemeList(PulldownMenu):
 
  def __init__(self, parent, getColorSchemes, extra_label = '', *args, **kw):
 
    self.getColorSchemes = getColorSchemes
    self.extra_label = extra_label
 
    PulldownMenu.__init__(self, parent, *args, **kw)
 
    for func in notify_funcs:
      Implementation.registerNotify(self.setColorSchemes, 'ccpnmr.Analysis.ColorScheme', func)
 
    self.setColorSchemes()

  def destroy(self):
 
    for func in notify_funcs:
      Implementation.unregisterNotify(self.setColorSchemes, 'ccpnmr.Analysis.ColorScheme', func)

    PulldownMenu.destroy(self)

  def setColorSchemes(self, *color):
 
    schemes = self.getColorSchemes()
    names = [ scheme.name for scheme in schemes ]
    #names.sort()
    if (self.extra_label):
      names = names + [self.extra_label]
    #print 'setColorSchemes', names, self.selected_index
    
    hexColors = []
    for scheme in schemes:
      schemeColors = []
      for c in scheme.colors:
        schemeColors.append(hexRepr(c.r,c.g,c.b))
 
      hexColors.append(schemeColors)
    
    self.replace(names, self.selected_index, colors=hexColors)
