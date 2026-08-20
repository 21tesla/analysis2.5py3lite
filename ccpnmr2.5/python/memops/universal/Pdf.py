"""
======================COPYRIGHT/LICENSE START==========================

Pdf.py: Utility code for CCPN code generation framework

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
# Note:
#
# Standard Type 1 PDF fonts:
#
# Courier (includes -Bold, -Oblique, -BoldOblique)
# Helvetica (includes -Bold, -Oblique, -BoldOblique)
# Times-Roman (includes Times-Bold, Times-Italic, Times-BoldItalic)
# Symbol
# ZapfDingbats
#

import math
import os

import memops.universal.Output as Output

default_linewidth = Output.default_linewidth

format_choices = [ 'PDF' ]
format_suffixes = [ 'pdf' ]

class Pdf(Output.Output):

  def __init__(self, file_name='', plot_size=None, plot_units='cm',
               paper_size=None, paper_orientation='Portrait',
               output_style='Color', border_text=None,
               x_axis_label='', y_axis_label='',
               axis_label_font='Times-Roman', axis_label_size=12,
               axis_label_offset=0, linewidth=default_linewidth,
               fonts=None, do_outline_box=True):

    self.output_format = 'PDF' # not really needed

    self.nbytes = 0
    self.objs = []

    # TEMP
    ###fonts = ['Times-Roman']

    Output.Output.__init__(self, file_name=file_name, plot_size=plot_size,
                           plot_units=plot_units, paper_size=paper_size,
                           paper_orientation=paper_orientation,
                           output_style=output_style, border_text=border_text,
                           x_axis_label=x_axis_label, y_axis_label=y_axis_label,
                           axis_label_font=axis_label_font, axis_label_size=axis_label_size,
                           axis_label_offset=axis_label_offset, linewidth=linewidth,
                           fonts=fonts, do_outline_box=do_outline_box,
                           file_mode='wb')

    self.setFont('Times-Roman', 12) # arbitrary, should not be needed

  def outputData(self, data):
    """Write output to the stream
    Only keep open whilst writing to ensure no clash with other threads
    """
    self._midStreamOpen()
    self.stream.write(data)
    self._midStreamClose()

    self.nbytes = self.nbytes + len(data)

  def outputRef(self, obj = None):

    if (obj):
      byte = obj[1]
      gen = 0
      typ = 'n'
    else:
      byte = 0
      gen = 65535
      typ = 'f'
    self.outputData('%010d %05d %s \n' % (byte, gen, typ))

  def outputObject(self, n, s):

    if (s[0] == '\n'):
      s = s[1:]

    obj = (n, self.nbytes)
    self.objs.append(obj)

    self.outputData('%d 0 obj\n' % n)
    self.outputData(s)
    self.outputData('endobj\n')

  def outputHeader(self):

    self.outputData('%PDF-1.4\n')
    self.outputData('%' + chr(16*14+2) + chr(16*14+3) + chr(16*12+15) + chr(16*13+3) + '\n')

    # object 1: catalog
    # n = object number
    self.outputObject(1, '''
<<
/Type /Catalog
/Pages 2 0 R
>>
''')
 
    # object 2: pages
    self.outputObject(2, '''
<<
/Type /Pages
/Kids [ 3 0 R ]
/Count 1
>>
''')

    # object 3: page
    nobjs = len(self.fonts)
    # cannot use (w, h) for MediaBox
    # because that only works for Portrait, not Landscape
    (w, h) = self.paper_size
    (w, h) = (min(w, h), max(w, h))
    self.outputObject(3, '''
<<
/Type /Page
/Parent 2 0 R
/Resources 4 0 R
/Contents 5 0 R
/MediaBox [ 0 0 %2.1f %2.1f ]
>>
''' % (w, h))

    # object 4: resources
    s = '/Font <<'
    m = 7 # first font object
    for font in self.fonts:
      #s = s + '/Font << /%s %d 0 R >>\n' % (font, m)
      s = s + ' /%s %d 0 R' % (font, m)
      m = m + 1
    s = s + ' >>\n'
    self.outputObject(4, '''
<<
/ProcSet [ /PDF /Text ]
%s>>
''' % s)

    # font objects (objects 7 through (7 + len(self.fonts)-1))
    m = 7 # first font object
    for font in self.fonts:
      self.outputObject(m, '''
<<
/Type /Font
/Subtype /Type1
/BaseFont /%s
>>
''' % font)
      m = m + 1

    # object 5: contents (i.e. the drawing)
    obj = (5, self.nbytes)
    self.objs.append(obj)
    self.outputData('5 0 obj\n<< /Length 6 0 R >>\nstream\n')
    self.drawingByte = self.nbytes

  def outputTrailer(self):

    # draw stuff in C world so use file size to figure out correct value for self.nbytes
    # self.stream.flush() # make sure all written
    self.nbytes = os.stat(self.file_name)[6]
    length = self.nbytes - self.drawingByte

    self.outputData('endstream\nendobj\n')
   
    # object 6: content length
    self.outputObject(6, '''
%d
''' % length)

    nbytes = self.nbytes
    self.outputData('xref\n')
    self.outputData('0 %d\n' % (1+len(self.objs)))
 
    self.outputRef()
    self.objs.sort()
    for obj in self.objs:
      self.outputRef(obj)
 
    self.outputData('trailer\n')
    self.outputData('<<\n')
    self.outputData('/Size %d\n' % len(self.objs))
    self.outputData('/Root 1 0 R\n')
    # TBD: should do /ID as well
    self.outputData('>>\n')
    self.outputData('startxref\n')
    self.outputData('%d\n' % nbytes)
    self.outputData('%%EOF\n')

  def setBlackWhite(self):

    self.black_white = True

  def setColorGray(self):

    self.black_white = False

  def rectangleClip(self, x0, y0, x1, y1):

    (x0, y0, x1, y1) = [ float(t) for t in (x0, y0, x1, y1) ]

    self.outputData('%3.2f %3.2f m\n' % (x0, y0))
    self.outputData('%3.2f %3.2f l\n' % (x1, y0))
    self.outputData('%3.2f %3.2f l\n' % (x1, y1))
    self.outputData('%3.2f %3.2f l\n' % (x0, y1))
    self.outputData('%3.2f %3.2f l\n' % (x0, y0))
    self.outputData('W n\n')

  def rotate(self, angle):

    angle = math.pi * angle / 180.0
    c = math.cos(angle)
    s = math.sin(angle)

    self.outputData('%3.2f %3.2f %3.2f %3.2f 0 0 cm\n' % (c, s, -s, c))

  def translate(self, x, y):

    self.outputData('1 0 0 1 %3.2f %3.2f cm\n' % (float(x), float(y)))

  def drawLine(self, x0, y0, x1, y1):

    self.outputData('%3.2f %3.2f m\n' % (self.ax*x0, self.ay*y0))
    self.outputData('%3.2f %3.2f l\n' % (self.ax*x1, self.ay*y1))
    self.outputData('S\n')

  NCIRCLE_PTS = 20
  DEGREES_PER_PT = 360 / NCIRCLE_PTS
  RADIANS_PER_PT = DEGREES_PER_PT * math.pi / 180

  def fillCircle(self, x, y, r):

    xx = x + r
    yy = y
    self.outputData('%3.2f %3.2f m\n' % (self.ax*xx, self.ay*yy))
    for i in range(self.NCIRCLE_PTS):
      angle = (i+1) * self.RADIANS_PER_PT
      xx = x + r * math.cos(angle)
      yy = y + r * math.sin(angle)
      self.outputData('%3.2f %3.2f l\n' % (self.ax*xx, self.ay*yy))
    self.outputData('f\n')

  def drawCircle(self, x, y, r):

    xx = x + r
    yy = y
    self.outputData('%3.2f %3.2f m\n' % (self.ax*xx, self.ay*yy))
    for i in range(self.NCIRCLE_PTS):
      angle = (i+1) * self.RADIANS_PER_PT
      xx = x + r * math.cos(angle)
      yy = y + r * math.sin(angle)
      self.outputData('%3.2f %3.2f l\n' % (self.ax*xx, self.ay*yy))
    self.outputData('S\n')

  def drawDashLine(self, x0, y0, x1, y1, dashLength, gapLength):

    self.outputData('[ %d %d ] 0 d\n' % (dashLength, gapLength))
    self.drawLine(x0, y0, x1, y1)
    self.setNonDashed()

  # fudge factors needed to estimate bounding box of text
  # so that can try to correctly place text for a, b non-zero
  xfraction = 0.45
  yfraction = 0.9

  def drawText(self, text, x, y, a, b):

    x = self.ax * x
    y = self.ay * y
    fontsize = self.font[1]
    if a:
      x = x - a * self.xfraction * fontsize * len(text)
    if b:
      y = y - b * self.yfraction * fontsize

    # escape special PS characters
    text = text.replace('\\', '\\\\')
    text = text.replace('(', '\\(')
    text = text.replace(')', '\\)')

    self.outputData('BT\n')
    self.outputData('/%s %d Tf\n' % self.font)
    self.outputData('%3.2f %3.2f Td\n' % (x, y))
    self.outputData('(%s) Tj\n' % text)
    self.outputData('ET\n')
 
  # r, g, b should be between 0 and 1
  def setColor(self, color):

    (r, g, b) = color
    if (self.black_white):
      if (r + g + b < 2.99):
        r = g = b = 0.0
      else:
        r = g = b = 1.0

    self.outputData('%3.2f %3.2f %3.2f RG\n' % (r, g, b))
    self.outputData('%3.2f %3.2f %3.2f rg\n' % (r, g, b))

  def setFont(self, name, size):

    if name == 'Times':
      name = 'Times-Roman'

    self.font = (name, size)

  def setNonDashed(self):

    self.outputData('[] 0 d\n')

  def setDashed(self):

    self.outputData('[ 1 2 ] 0 d\n')

  def setLineWidth(self, linewidth):

    self.outputData('%3.2f w\n' % float(linewidth))

  def save(self):

    self.outputData('q\n')

  def restore(self):

    self.outputData('Q\n')

if (__name__ == '__main__'):

  p = Pdf(file_name='test.pdf', plot_size=(12.0, 20.0), do_outline_box=True)
  p.newRange(-100.0, -100.0, 100.0, 100.0, True)
  p.setColor((1.0, 0.0, 0.0))
  p.drawLine(-200, -200, 300, 150)
  p.setFont('Times-Roman', 14)
  p.setColor((0.0, 0.0, 1.0))
  p.drawText('hello world', -10, -10, 0.0, 0.0)
  p.setColor((0.0, 1.0, 0.0))
  p.drawText('hello world', -10, -10, 1.0, 1.0)
  p.close()
