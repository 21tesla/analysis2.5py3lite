"""
======================COPYRIGHT/LICENSE START==========================

PrintTicks.py: Utility code for CCPN code generation framework

Copyright (C) 2005  (CCPN Project)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and automated
software development. Bioinformatics 21, 1678-1684.

===========================REFERENCE END===============================

"""
from memops.universal.Ticks import Ticks

Inside = 'Inside'
Outside = 'Outside'
Top = 'Top'
Bottom = 'Bottom'
Left = 'Left'
Right = 'Right'
tick_locations = [ Inside, Outside ]
tick_placements = [ Top, Bottom, Left, Right ]
##tick_size = 25 # pts

no_tick_text = 'No tick label'

class PrintTicks(Ticks):

  def __init__(self, printHandler, plot_size, tick_location, tick_placement, deltaMajor=None, deltaMinor=None, numberDecimals=None, font=None, tickMajor=None, tickMinor=None):

    #print 'PrintTicks', plot_size, tick_location, tick_placement
    if not font:
      font = 'Helvetica 10'

    if font == no_tick_text:
      font = None
      isTextDrawn = False
      font_size = 10
    else:
      isTextDrawn = True
      font = font.split(' ')[:2]
      font[1] = int(font[1])
      font_size = font[1]

    tick_size = 3 * font_size

    outputHandler = self.outputHandler = printHandler.outputHandler
    outputHandler.save()

    if tick_placement in (Bottom, Top):
      tx = printHandler.x * plot_size[0] / printHandler.w
      if tick_placement == Bottom:
        if tick_location == Inside:
          ty = 0
        else:
          ty = - tick_size
      else:
        if tick_location == Inside:
          ty = plot_size[1] - tick_size
        else:
          ty = plot_size[1]
    else:
      ty = printHandler.y * plot_size[1] / printHandler.h
      if tick_placement == Left:
        if tick_location == Inside:
          tx = 0
        else:
          tx = - tick_size
      else:
        if tick_location == Inside:
          tx = plot_size[0] - tick_size
        else:
          tx = plot_size[0]

    outputHandler.translate(tx, ty)

    if (tick_placement in (Bottom, Top)):
      isHorizontal = True
      length = printHandler.xSize
      width = float(printHandler.xSize * plot_size[0]) / printHandler.w
      height = tick_size
      region = printHandler.xRegion
      #outputHandler.rectangleClip(0, 0, width, tick_size)
    else:
      isHorizontal = False
      length = printHandler.ySize
      width = tick_size
      height = float(printHandler.ySize * plot_size[1]) / printHandler.h
      region = printHandler.yRegion
      #outputHandler.rectangleClip(0, 0, tick_size, height)

    if (tick_placement in (Bottom, Left)):
      if (tick_location == Inside):
        haveTicksOnZeroSide = True
      else:
        haveTicksOnZeroSide = False
    else:
      if (tick_location == Inside):
        haveTicksOnZeroSide = False
      else:
        haveTicksOnZeroSide = True

    Ticks.__init__(self, width=width, height=height, isHorizontal=isHorizontal,
                   font=font, doRuleLines=False, haveYInversion=False,
                   haveTicksOnZeroSide=haveTicksOnZeroSide, isTextDrawn=isTextDrawn,
                   deltaMajor=deltaMajor, deltaMinor=deltaMinor,
                   numberDecimals=numberDecimals, tickMajor=tickMajor, tickMinor=tickMinor)
    if font:
      outputHandler.setFont(*font)
    self.setMaxMajorMinor(length)
    self.setRegion(region) # does actual drawing
    outputHandler.restore()

  def drawLine(self, x0, y0, x1, y1):

    self.outputHandler.drawLine(x0, y0, x1, y1)

  def drawText(self, text, x, y, a, b):

    self.outputHandler.drawText(text, x, y, a, b)

