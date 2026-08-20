"""
======================COPYRIGHT/LICENSE START==========================

PeakDrawUtil.py: Part of the CcpNmr Analysis program

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
# do not change below lightly!

peak_draw_methods = [ 'uniform in pixels', 'uniform in ppm',
                      'scaled globally', 'scaled by peaks in peaklist',
                      'line width', 'box width' ]

# note: the below has to be consistent with C code
# matching C code is in c/ccpnmr/analysis/win_peak_list.h:
DRAW_UNIFORM_METHOD  = 0
DRAW_GLOBAL_METHOD   = 1
DRAW_PEAKLIST_METHOD = 2
DRAW_LINE_WIDTH_METHOD = 3

c_peak_draw_methods = {}
c_peak_draw_methods[peak_draw_methods[0]] = DRAW_UNIFORM_METHOD
c_peak_draw_methods[peak_draw_methods[1]] = DRAW_UNIFORM_METHOD
c_peak_draw_methods[peak_draw_methods[2]] = DRAW_GLOBAL_METHOD
c_peak_draw_methods[peak_draw_methods[3]] = DRAW_PEAKLIST_METHOD
c_peak_draw_methods[peak_draw_methods[4]] = DRAW_LINE_WIDTH_METHOD
c_peak_draw_methods[peak_draw_methods[5]] = DRAW_UNIFORM_METHOD

def getCPeakDrawMethod(project):

  method = project.currentAnalysisProject.peakDrawMethod

  return c_peak_draw_methods[method]

