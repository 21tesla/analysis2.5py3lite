"""
======================COPYRIGHT/LICENSE START==========================

PostScript.py: Utility code for CCPN code generation framework

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
import math

import memops.universal.Output as Output

default_linewidth = Output.default_linewidth

format_choices = [ 'PostScript', 'Encapsulated PostScript' ]
format_suffixes = [ 'ps', 'eps' ]

ps_header = \
        '%!PS-Adobe-3.0\n'
eps_header = \
        '%!PS-Adobe-3.0 EPSF-3.0\n'
ps_comments = \
        '%%%%Creator: CCPN (C) 2003 CCPN\n' \
        '%%%%BoundingBox: 0 0 %2.1f %2.1f\n' \
        '%%%%DocumentNeededResources: font %s\n' \
        '%%%%Pages: 1\n' \
        '%%%%EndComments\n'
ps_setup = \
        '%%%%BeginSetup\n' \
        '1 setlinecap 1 setlinejoin %2.1f setlinewidth\n' \
        '0 setgray [ ] 0 setdash newpath\n' \
        '%%%%EndSetup\n' % default_linewidth
ps_begin = \
        '%%Page: 1 1\n'
ps_end = \
        'PS_end\n' \
        'CCPNSave restore\n' \
        'showpage\n' \
        '%%Trailer\n'
ps_eof = \
        '%%EOF\n'
ps_prolog = \
        '%%BeginProlog\n' \
        '/CCPNSave save def\n' \
        '128 dict begin\n' \
        '/PS_graphics_save { gsave } bind def\n' \
        '/PS_graphics_restore { grestore } bind def\n' \
        '/PS_save { save } bind def\n' \
        '/PS_restore { restore } bind def\n' \
        '/PS_begin { begin } bind def\n' \
        '/PS_end { end } bind def\n' \
        '/PS_show_page { showpage } bind def\n' \
        '/PS_new_path { newpath } bind def\n' \
        '/PS_draw_line { newpath 4 2 roll moveto lineto stroke } bind def\n' \
        '/PS_fill_circle { newpath 0 360 arc fill } bind def\n' \
        '/PS_draw_circle { newpath 0 360 arc stroke } bind def\n' \
        '/PS_fill_triangle { newpath 6 2 roll moveto 4 2 roll lineto lineto closepath fill } bind def \n' \
        '/PS_set_black { 0 setgray } bind def\n' \
        '/PS_set_white { 1 setgray } bind def\n' \
        '/PS_translate { translate } bind def\n' \
        '/PS_rotate { rotate } bind def\n' \
        '/PS_scale { scale } bind def\n' \
        '/PS_nondashed { [ ] 0 setdash } bind def\n' \
        '/PS_dashed { [ 1 2 ] 0 setdash } bind def\n' \
        '/PS_setdash { setdash } bind def\n' \
        '/PS_set_linewidth { setlinewidth } bind def\n' \
        '/NP { newpath } bind def\n' \
        '/M { moveto } bind def\n' \
        '/L { lineto } bind def\n' \
        '/S { stroke } bind def\n' \
        '/CP { closepath } bind def\n' \
        \
        '/PS_black_white {\n' \
        '  /PS_set_color {\n' \
        '    add add 2.99 lt { PS_set_black } { PS_set_white } ifelse\n' \
        '  } bind def\n' \
        '} bind def\n' \
        \
        '/PS_color_gray {\n' \
        '  /PS_set_color { setrgbcolor } bind def\n' \
        '} bind def\n' \
        \
        '/PS_text {\n' \
        '  newpath 0 0 moveto dup false charpath pathbbox\n' \
        '  3 2 roll sub 5 add 3 1 roll sub neg 5 add\n' \
        '  5 4 roll mul exch 4 3 roll mul\n' \
        '  4 3 roll sub neg exch 4 3 roll sub neg exch moveto show\n' \
        '} bind def\n' \
        \
        '/do_rectangle {\n' \
        '  newpath 4 2 roll moveto 1 index 0 rlineto\n' \
        '  0 exch rlineto neg 0 rlineto closepath\n' \
        '} bind def\n' \
        \
        '/PS_draw_rectangle {\n' \
        '  /rectstroke where { pop rectstroke }\n' \
        '  { do_rectangle stroke } ifelse\n' \
        '} bind def\n' \
        \
        '/PS_rectangle_clip {\n' \
        '  /rectclip where { pop rectclip }\n' \
        '  { do_rectangle clip newpath } ifelse\n' \
        '} bind def\n' \
        \
        '/PS_select_font {\n' \
        '  /selectfont where { pop selectfont }\n' \
        '  { exch findfont exch scalefont setfont } ifelse\n' \
        '} bind def\n' \
        '%%EndProlog\n'
        # rectstroke, rectclip, selectfont are level 2 operators

class PostScript(Output.Output):

  def __init__(self, file_name='', plot_size=None, plot_units='cm',
               paper_size=None, paper_orientation='Portrait',
               output_style='Color', border_text=None,
               x_axis_label='', y_axis_label='',
               axis_label_font='Times-Roman', axis_label_size=12,
               axis_label_offset=0, linewidth=default_linewidth,
               fonts=None, do_outline_box=True, output_format='PostScript'):

    self.output_format = output_format
    Output.Output.__init__(self, file_name=file_name, plot_size=plot_size,
                           plot_units=plot_units, paper_size=paper_size,
                           paper_orientation=paper_orientation,
                           output_style=output_style, border_text=border_text,
                           x_axis_label=x_axis_label, y_axis_label=y_axis_label,
                           axis_label_font=axis_label_font, axis_label_size=axis_label_size,
                           axis_label_offset=axis_label_offset, linewidth=linewidth,
                           fonts=fonts, do_outline_box=do_outline_box,
                           file_mode='w')

  def outputData(self, data):
    """Write output to the stream
    Only keep open whilst writing to ensure no clash with other threads
    """
    self._midStreamOpen()
    self.stream.write(data)
    self._midStreamClose()

  def outputHeader(self):

    outputData = self.outputData

    if (self.output_format == 'PostScript'):
      header = ps_header
    else:
      header = eps_header
    outputData(header)

    (w, h) = self.paper_size
    font = ' '.join(self.fonts)
    outputData(ps_comments % (w, h, font))

    outputData(ps_prolog)
    outputData(ps_setup)
    outputData(ps_begin)

  def outputTrailer(self):

    self.outputData(ps_end)
    self.outputData(ps_eof)

  def setBlackWhite(self):

    self.outputData('PS_black_white\n')

  def setColorGray(self):

    self.outputData('PS_color_gray\n')

  def rectangleClip(self, x0, y0, x1, y1):

    (x0, y0, x1, y1) = [ float(t) for t in (x0, y0, x1, y1) ]

    self.outputData('%3.2f %3.2f %3.2f %3.2f PS_rectangle_clip\n' % (x0, y0, x1, y1))

  def rotate(self, angle):

    self.outputData('%3.2f PS_rotate\n' % float(angle))

  def translate(self, x, y):

    self.outputData('%3.2f %3.2f PS_translate\n' % (float(x), float(y)))

  def drawLine(self, x0, y0, x1, y1):

    self.outputData('%3.2f %3.2f %3.2f %3.2f PS_draw_line\n' % \
                      (self.ax*x0, self.ay*y0, self.ax*x1, self.ay*y1))

  def fillCircle(self, x, y, r):

    s = math.sqrt(self.ax * self.ay)

    self.outputData('%3.2f %3.2f %3.2f PS_fill_circle\n' % \
                      (self.ax*x, self.ay*y, s*r))

  def drawCircle(self, x, y, r):

    s = math.sqrt(self.ax * self.ay)

    self.outputData('%3.2f %3.2f %3.2f PS_draw_circle\n' % \
                      (self.ax*x, self.ay*y, s*r))

  def drawDashLine(self, x0, y0, x1, y1, dashLength, gapLength):

    self.outputData('[%d %d] 0 PS_setdash\n' % (dashLength, gapLength))
    self.drawLine(x0, y0, x1, y1)
    self.setNonDashed()

  def drawText(self, text, x, y, a, b):

    # escape special PS characters
    text = text.replace('\\', '\\\\')
    text = text.replace('(', '\\(')
    text = text.replace(')', '\\)')

    self.outputData('%3.2f %3.2f %3.2f %3.2f (%s) PS_text\n' % \
                      (self.ax*x, self.ay*y, a, b, text))
 
  # r, g, b should be between 0 and 1
  def setColor(self, color):

    self.outputData('%3.2f %3.2f %3.2f PS_set_color\n' % color)

  def setFont(self, name, size):

    if name == 'Times':
      name = 'Times-Roman'

    self.outputData('/%s %d PS_select_font\n' % (name, size))

  def setNonDashed(self):

    self.outputData('PS_nondashed\n')

  def setDashed(self):

    self.outputData('PS_dashed\n')

  def setLineWidth(self, linewidth):

    self.outputData('%3.2f PS_set_linewidth\n' % linewidth)

  def save(self):

    self.outputData('PS_graphics_save\n')

  def restore(self):

    self.outputData('PS_graphics_restore\n')

  def drawRectangle(self, rectangle):

    (x, y, width, height) = rectangle
    x *= self.ax
    y *= self.ay
    width *= self.ax
    height *= self.ay
    self.outputData('%3.2f %3.2f %3.2f %3.2f PS_draw_rectangle\n' % (x, y, width, height))

if (__name__ == '__main__'):

  p = PostScript(file_name='test.ps', plot_size=(12.0, 20.0), do_outline_box=True)
  p.newRange(-100.0, -100.0, 100.0, 100.0, True)
  p.drawLine(-100, -100, 300, 150)
  p.setFont('Times-Roman', 14)
  p.drawText('hello world', -100, -100, 0.0, 0.0)
  p.close()
