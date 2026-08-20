LICENSE = """
======================COPYRIGHT/LICENSE START==========================

Output.py: Utility code for CCPN code generation framework

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

other_paper_type = 'Other'
paper_types = [ 'A4', 'A3', 'A', 'B', other_paper_type ]

# mapping from data model PaperFormat enumeration
paper_type_dict = {
  'A3': 'A3',
  'A4': 'A4',
  'A5': 'A4',
  'B4': 'A4',
  'Letter': 'A',
  'Legal': 'A',
  'Tabloid': 'B',
}

# mapping to data model PaperFormat enumeration
paper_type_inverse_dict = {
  'A3': 'A3',
  'A4': 'A4',
  'A': 'Letter',
  'B': 'Tabloid',
}

paper_units = [ 'cm', 'in', 'pt' ]
paper_sizes = { \
  'A4' : (21.0, 29.7, 'cm'),
  'A3' : (29.7, 42.0, 'cm'),
  'A' : (8.5, 11.0, 'in'),
  'B' : (11.0, 17.0, 'in'),
  }
paper_orientations = [ 'Portrait', 'Landscape' ]
style_choices = [ 'Color', 'Black and white' ]

cms_per_in = 2.54
pts_per_in = 72.0
pts_conversion = { \
  'cm': pts_per_in / cms_per_in,
  'in': pts_per_in,
  'pt': 1.0,
  }

default_linewidth = 0.1 # pts

scaling_choices = [ 'Percentage', 'Cms / unit', 'Inches / unit',
                    'Units / cm', 'Units / inch' ]

font_names = ['Arial', 'Courier', 'Helvetica', 'Lucida', 'Times', 'Verdana']
font_sizes = range(6, 16)

class Output:

  def __init__(self, file_name='', plot_size=None, plot_units='cm',
               paper_size=None, paper_orientation='Portrait',
               output_style='Color', border_text=None,
               x_axis_label='', y_axis_label='',
               axis_label_font='Times-Roman', axis_label_size=12,
               axis_label_offset=0, linewidth=default_linewidth,
               fonts=None, do_outline_box=True, file_mode='w'):

    if (not paper_size):
      paper_size = paper_sizes['A4']

    if not fonts:
      fonts = ['Times-Roman']

    if axis_label_font not in fonts:
      fonts.append(axis_label_font)

    (w, h, u) = paper_size
    s = pts_conversion[u]
    w = s * w
    h = s * h
    if (paper_orientation == 'Landscape'):
      (w, h) = (h, w)

    if (plot_size):
      (ww, hh) = plot_size
      s = pts_conversion[plot_units]
      plot_size = (s*ww, s*hh)
    else:
      plot_size = (w, h)

    self.current_range = (0, 0, plot_size[0], plot_size[1])

    # dictionary with key = location for decoration, value = (text, font, size)
    # location can be one of 'nw', 'n', 'ne', 'w', 'e', 'sw', 's', 'e'
    if (not border_text):
      border_text = {}

    self.plot_size = plot_size
    self.paper_size = paper_size
    self.paper_orientation = paper_orientation
    self.output_style = output_style
    self.border_text = border_text
    self.x_axis_label = x_axis_label
    self.y_axis_label = y_axis_label
    self.axis_label_font = axis_label_font
    self.axis_label_size = axis_label_size
    self.axis_label_offset = axis_label_offset
    self.linewidth = linewidth
    self.fonts = fonts
    self.do_outline_box = do_outline_box

    self.paper_size = (w, h)

    self.have_graphics_save = False
    self.file_name = file_name
    self.stream = None
    self._file_mode = 'a'

    if (file_name):
      self.open(file_name, file_mode)

  def open(self, file_name, file_mode = 'w'):

    if self.stream:
      self.close()

    self.file_name = file_name

    self.ax = self.ay = 1.0

    self.file_name = file_name
    with open(file_name, 'w'):
      pass

    self.stream = open(file_name, file_mode)

    self.outputHeader()

    (w, h) = self.paper_size

    # rotate and translate if landscape
    if (self.paper_orientation == 'Landscape'):
      self.rotate(-90)
      self.translate(-w, 0)

    # set color or black and white
    if (self.output_style == 'Black and white'):
      self.use_dashed_lines = True
      self.setBlackWhite()
    else:
      self.use_dashed_lines = False
      self.setColorGray()

    self.setBlack()
    # self.resetLineWidth()
    self.setLineWidth(1.0)

    for location in self.border_text.keys():
      (text, font, size) = self.border_text[location]
      self.setFont(font, size)
      self.doBorder(text, location)

    (ww, hh) = self.plot_size
    x = 0.5 * (w - ww)
    y = 0.5 * (h - hh)

    # translate so that bottom left corner set up
    self.translate(x, y)

    offset = self.axis_label_offset
    if self.x_axis_label:
      self.setFont(self.axis_label_font, self.axis_label_size)
      self.drawText(self.x_axis_label, ww/2, 0.0, 0.5, offset+2.0)

    if self.y_axis_label:
      if not self.x_axis_label:
        self.setFont(self.axis_label_font, self.axis_label_size)
      self.rotate(90)
      self.drawText(self.y_axis_label, hh/2, 0.0, 0.5, -(offset+1.5))
      self.rotate(-90)

    if (self.do_outline_box):
      self.drawRectangle((0, 0, ww, hh))

    self._midStreamClose()

  def _midStreamOpen(self):
    self.stream = open(self.file_name, self._file_mode)

  def _midStreamClose(self):
    if self.stream:
      self.stream.flush()
      self.stream.close()
      self.stream = None

  def close(self):

    # if (not self.stream):
    #   return

    if (self.have_graphics_save):
      self.restore()
      self.have_graphics_save = False

    self.outputTrailer()

    if self.stream:
      self.stream.close()
    self.stream = None
    self.file_name = ''

  # following functions can be overridden in subclass if desired

  def doBorder(self, text, location):

    (w, h) = self.paper_size
    d = 30.0 # points offset in from paper's edge

    if ('w' in location):
      x = d
      a = 0.0
    elif ('e' in location):
      x = w - d
      a = 1.0
    else: 
      x = 0.5 * w
      a = 0.5
 
    if ('s' in location):
      y = d
      b = 0.0
    elif ('n' in location):
      y = h - d
      b = 1.0
    else: 
      y = 0.5 * h
      b = 0.5
 
    #print 'doBorder', w, h, x, a, y, b, text, location
    self.drawText(text, x, y, a, b)

  def newRange(self, x0, y0, x1, y1, clip = False):

    if (self.have_graphics_save):
      self.restore()

    self.save()
    self.have_graphics_save = True

    self.current_range = (x0, y0, x1, y1)

    (w, h) = self.plot_size

    if (clip):
      self.rectangleClip(0, 0, w, h)

    self.ax = w / float(x1-x0)
    bx = - self.ax * x0
    self.ay = h / float(y1-y0)
    by = - self.ay * y0

    self.translate(bx, by)

  def clipRange(self, x0, y0, x1, y1):

    self.rectangleClip(x0, y0, x1-x0, y1-y0)

  def drawRectangle(self, rectangle):

    (x0, y0, width, height) = rectangle
    x1 = x0 + width
    y1 = y0 + height
    self.drawLine(x0, y0, x1, y0)
    self.drawLine(x1, y0, x1, y1)
    self.drawLine(x1, y1, x0, y1)
    self.drawLine(x0, y1, x0, y0)

  def setBlack(self):

    self.setColor((0, 0, 0))

  def setWhite(self):

    self.setColor((1, 1, 1))

  def setLineStyle(self, line_style):

    if (self.use_dashed_lines):
      if (line_style == 'dashed'):
        self.setNonDashed()
      else:
        self.setDashed()

  def resetLineWidth(self):

    self.setLineWidth(self.linewidth)

  # following functions should (in general) be overridden (implemented) in subclass

  def outputHeader(self):

    pass

  def outputTrailer(self):

    pass

  def setBlackWhite(self):

    pass

  def setColorGray(self):

    pass

  def rectangleClip(self, x0, y0, x1, y1):

    pass

  def rotate(self, angle):

    pass

  def translate(self, x, y):

    pass

  def drawLine(self, x0, y0, x1, y1):

    pass

  def drawDashLine(self, x0, y0, x1, y1, dashLength, gapLength):

    pass

  def drawText(self, text, x, y, a, b):

    pass
 
  # r, g, b should be between 0 and 1
  def setColor(self, color):

    pass

  def setFont(self, name, size):

    pass

  def setNonDashed(self):

    pass

  def setDashed(self):

    pass

  def setLineWidth(self, linewidth):

    pass

  def save(self):

    pass

  def restore(self):

    pass

if (__name__ == '__main__'):

  p = PostScript(file_name='test.ps', plot_size=(12.0, 20.0), do_outline_box=True)
  p.newRange(-100.0, -100.0, 100.0, 100.0, True)
  p.drawLine(-100, -100, 300, 150)
  p.setFont('Times-Roman', 14)
  p.drawText('hello world', -100, -100, 0.0, 0.0)
  p.close()
