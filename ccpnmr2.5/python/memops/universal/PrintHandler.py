"""
======================COPYRIGHT/LICENSE START==========================

PrintHandler.py: Utility code for CCPN code generation framework

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

import memops.universal.Output as Output
import memops.universal.Pdf as Pdf
import memops.universal.PostScript as PostScript

# C imports
try:
  import memops.c.PdfHandler as PdfHandler
except:
  PdfHandler = None
try:
  import memops.c.PsHandler as PsHandler
except:
  PsHandler = None

# useful wrapper class which makes it easier to do printing
# these involve an outputHandler (e.g. a PostScript handler)
# and potentially a cHandler (for doing drawing in C world)

class PrintHandler:

  # xRegion and yRegion are regions in display units (backwards if units backwards)
  # xSize and ySize are the sizes of these regions in pixels
  # x, y are offset for lower left corner in pixels
  # outputHandler should have attributes width, height, which is total size in pixels
  def __init__(self, outputHandler, xRegion, yRegion, xSize = None, ySize = None, x = 0, y = 0):

    if (hasattr(outputHandler, 'cHandler')):
      cHandler = outputHandler.cHandler
    else:
      cHandler = None

    if (xSize is None):
      xSize = outputHandler.width

    if (ySize is None):
      ySize = outputHandler.height

    self.outputHandler = outputHandler
    self.xRegion = xRegion
    self.yRegion = yRegion
    self.xSize = xSize
    self.ySize = ySize
    self.x = x
    self.y = y
    self.w = outputHandler.width
    self.h = outputHandler.height
    self.cHandler = cHandler

  # functions below are required interfaces which usual Print handlers don't have

  def makeCurrent(self):

    return True

  def startFront(self):

    pass

  def startBack(self):

    pass

  def expose(self, x, y, w, h):

    pass

  def swapBuffers(self):

    pass

  # mapRanges gets called such that spectrum point (a0, b0)
  # gets mapped to window point (x0, y0), where latter in range (0, 1),
  # and similarly for (a1, b1) and (x1, y1)
  # conversion below so that have correct mapping to corner
  # of plot on paper
  def mapRanges(self, x0, y0, x1, y1, a0, b0, a1, b1):

    (x0, y0, x1, y1, a0, b0, a1, b1) = [ float(z) for z in (x0, y0, x1, y1, a0, b0, a1, b1) ]
    #print 'mapRanges0', x0, y0, x1, y1, a0, b0, a1, b1

    s = self.xSize
    u = self.outputHandler.width
    x = self.x
    r = (a1 - a0) / (s*(x1 - x0))
    c0 = a0 - r * (x + x0*s)
    c1 = a0 + r * (u - x - x0*s)
    #print 'mapRanges1', s, u, x, r

    t = self.ySize
    v = self.outputHandler.height
    y = self.y
    r = (b1 - b0) / (t*(y1 - y0))
    d0 = b0 - r * (y + y0*t)
    d1 = b0 + r * (v - y - y0*t)

    #print 'mapRanges2', t, v, y, r
    #print 'mapRanges3', c0, d0, c1, d1
    self.outputHandler.newRange(c0, d0, c1, d1)
    self.outputHandler.clipRange(a0, b0, a1, b1)
    if (self.cHandler):
      self.cHandler.newRange(c0, d0, c1, d1)
      self.cHandler.clipRange(a0, b0, a1, b1)

  def setBackground(self, color):

    pass

  # pass all other requests onto outputHandler
  def __getattr__(self, name):

    return getattr(self.outputHandler, name)

def determinePlotSize(pixel_width, pixel_height, unit_width, unit_height,
                scaling_choice, output_scaling, w_scaling, h_scaling,
                paper_size, paper_orientation):

  if scaling_choice == Output.scaling_choices[0]:
    (ww, hh, u) = paper_size
    if (paper_orientation != 'Portrait'):
      (ww, hh) = (hh, ww)

    s = output_scaling * Output.pts_conversion[u]
    ww = s * ww
    hh = s * hh

    rw = pixel_width / ww
    rh = pixel_height / hh

    if (rw < rh):
      plot_size = (pixel_width / rh, hh)
    else:
      plot_size = (ww, pixel_height / rw)
  else:
    if scaling_choice in (Output.scaling_choices[1], Output.scaling_choices[3]):
      u = 'cm'
    else:
      u = 'in'
    
    s = Output.pts_conversion[u]
    ww = s * unit_width
    hh = s * unit_height
    if scaling_choice in (Output.scaling_choices[1], Output.scaling_choices[2]):
      ww *= w_scaling
      hh *= h_scaling
    else:
      ww /= w_scaling
      hh /= h_scaling
    plot_size = (ww, hh)

  return plot_size

# You only need getOutputHandler() if need to print multiple canvases
# Otherwise use getPrintHandler()

# width and height in pixels
def getOutputHandler(file_name, pixel_width, pixel_height, unit_width=1.0, unit_height=1.0,
                    scaling_choice=Output.scaling_choices[0], output_scaling=0.8,
                    w_scaling=1.0, h_scaling=1.0, paper_size=None,
                    paper_orientation='Portrait', output_style='Color',
                    output_format='PostScript', border_text=None,
                    x_axis_label='', y_axis_label='',
                    axis_label_font='Times-Roman', axis_label_size='12',
                    axis_label_offset=0, linewidth=Output.default_linewidth,
                    fonts=None, do_outline_box=True):

  if not paper_size:
    paper_size = Output.paper_sizes['A4']

  plot_size = determinePlotSize(pixel_width, pixel_height, unit_width, unit_height,
                scaling_choice, output_scaling, w_scaling, h_scaling,
                paper_size, paper_orientation)
  plot_units = 'pt'

  # TBD: other formats
  if (output_format in PostScript.format_choices):

    outputHandler = PostScript.PostScript(file_name=file_name, plot_size=plot_size,
                      plot_units=plot_units, paper_size=paper_size,
                      paper_orientation=paper_orientation, output_style=output_style,
                      output_format=output_format, border_text=border_text,
                      x_axis_label=x_axis_label, y_axis_label=y_axis_label,
                      axis_label_font=axis_label_font, axis_label_size=axis_label_size,
                      axis_label_offset=axis_label_offset, linewidth=linewidth,
                      fonts=fonts, do_outline_box=do_outline_box)
    if (PsHandler):
      # cHandler = PsHandler.PsHandler(outputHandler.stream, plot_size[0], plot_size[1],
      #                 output_style)
      cHandler = PsHandler.PsHandler(file_name, plot_size[0], plot_size[1],
                                     output_style)
    else:
      cHandler = None

  elif (output_format in Pdf.format_choices):

    outputHandler = Pdf.Pdf(file_name=file_name, plot_size=plot_size,
                      plot_units=plot_units, paper_size=paper_size,
                      paper_orientation=paper_orientation, output_style=output_style,
                      border_text=border_text,
                      x_axis_label=x_axis_label, y_axis_label=y_axis_label,
                      axis_label_font=axis_label_font, axis_label_size=axis_label_size,
                      axis_label_offset=axis_label_offset, linewidth=linewidth,
                      fonts=fonts, do_outline_box=do_outline_box)
    if (PdfHandler):
      # cHandler = PdfHandler.PdfHandler(outputHandler.stream, plot_size[0], plot_size[1],
      #                 output_style)
      cHandler = PdfHandler.PdfHandler(file_name, plot_size[0], plot_size[1],
                                       output_style)
    else:
      cHandler = None

  else:

    raise Exception('output format "%s" not supported yet' % output_format)

  # TBD: not sure the first, second and fourth are correct for new scaling_choices
  outputHandler.width = pixel_width
  outputHandler.height = pixel_height
  outputHandler.cHandler = cHandler
  outputHandler.output_scaling = output_scaling

  return outputHandler

# Use getPrintHandler rather then getOutputHandler if only need one region printed out

# xRegion and yRegion are regions in display units (backwards if units backwards)
# width and height in pixels
def getPrintHandler(file_name, xRegion, yRegion, width, height,
                    output_scaling = None, paper_size = None,
                    paper_orientation = 'Portrait', output_style = 'Color',
                    output_format = 'PostScript', border_text = None,
                    fonts = None, do_outline_box = True):

  outputHandler = getOutputHandler(file_name, width, height,
                    output_scaling=output_scaling,
                    paper_size=paper_size,
                    paper_orientation=paper_orientation,
                    output_style=output_style,
                    output_format=output_format,
                    border_text=border_text,
                    fonts=fonts,
                    do_outline_box=do_outline_box)

  printHandler = PrintHandler(outputHandler, xRegion, yRegion)

  return printHandler
