
"""
======================COPYRIGHT/LICENSE START==========================

PrintFrame.py: <write function here>

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

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
- PDBe website (http://www.ebi.ac.uk/pdbe/)

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
import os
import time

import Tkinter

import memops.universal.Output as Output
import memops.universal.Pdf as Pdf
import memops.universal.PostScript as PostScript
import memops.universal.PrintHandler as PrintHandler
import memops.universal.PrintTicks as PrintTicks

from memops.gui.Button import Button
from memops.gui.CheckButtons import CheckButtons
from memops.gui.Entry import Entry
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FontList import FontList
from memops.gui.Frame import Frame
from memops.gui.FloatEntry import FloatEntry
from memops.gui.IntEntry import IntEntry
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError, showInfo, showYesNo
from memops.gui.PulldownList import PulldownList
from memops.gui.Scale import Scale

paper_types = Output.paper_types
paper_type_dict = Output.paper_type_dict
paper_type_inverse_dict = Output.paper_type_inverse_dict
paper_units = Output.paper_units
paper_sizes = Output.paper_sizes
paper_orientations = Output.paper_orientations
style_choices = Output.style_choices
format_choices = PostScript.format_choices + Pdf.format_choices
format_suffixes = PostScript.format_suffixes + Pdf.format_suffixes
border_decorations = [ 'Time and Date', 'File Name' ]
tick_locations = PrintTicks.tick_locations
tick_placements = PrintTicks.tick_placements

format_options = [format_suffix.upper() for format_suffix in format_suffixes]

format_suffixes = ['.' + format_suffix for format_suffix in format_suffixes]
file_suffixes = {}
for n in range(len(format_choices)):
  file_suffixes[format_choices[n]] = format_suffixes[n]

scaling_choices = Output.scaling_choices
spacing_choices = [ 'Auto', 'Manual' ]
tick_length_choices = [ 'Auto', 'Manual' ]

""" Below are option keys in data model:
  'FileName',
  'Title',
  'PaperSize',
  'OtherHeight',
  'OtherWidth',
  'OtherUnit',
  'Orientation',
  'InColor',
  'OutputFormat',
  'TickOutside',
  'TickPlacement',
  'ShowsDateTime',
  'ShowsFileName',
  'Scaling',
"""

# convert from nsew to top/botton/left/right
def getTickPlacement1(placement):

  result = []

  if 'n' in placement:
    result.append(PrintTicks.Top)

  if 's' in placement:
    result.append(PrintTicks.Bottom)

  if 'e' in placement:
    result.append(PrintTicks.Right)

  if 'w' in placement:
    result.append(PrintTicks.Left)

  return result

# convert from top/botton/left/right to nsew
def getTickPlacement2(placement):

  result = ''

  if PrintTicks.Top in placement:
    result = result + 'n'

  if PrintTicks.Bottom in placement:
    result = result + 's'

  if PrintTicks.Right in placement:
    result = result + 'e'

  if PrintTicks.Left in placement:
    result = result + 'w'

  return result

class PrintFrame(Frame):

  def __init__(self, parent, getOption = None, setOption = None,
               haveTicks = False, doOutlineBox = True, *args, **kw):

    self.getOption = getOption
    self.setOption = setOption
    self.haveTicks = haveTicks
    self.doOutlineBox = doOutlineBox

    Frame.__init__(self, parent=parent, *args, **kw)

    self.file_select_popup = None

    self.getOptionValues()

    try:
      size_index = paper_types.index(self.paper_type)
    except:
      size_index = 0

    try:
      other_unit_index = paper_units.index(self.other_unit)
    except:
      other_unit_index = 0

    try:
      orientation_index = paper_orientations.index(self.paper_orientation)
    except:
      orientation_index = 0

    try:
      style_index = style_choices.index(self.output_style)
    except:
      style_index = 0

    try:
      format_index = format_choices.index(self.output_format)
    except:
      format_index = 0

    if haveTicks:
      try:
        tick_location_index = tick_locations.index(self.tick_location)
      except:
        tick_location_index = 0

    self.grid_columnconfigure(1, weight=1)

    row = 0
    button = Button(self, text='File:', command=self.findFile,
                    tipText='Select location to save print file')
    button.grid(row=row, column=0, sticky='e')
    self.file_entry = Entry(self, width=40, text=self.file_name,
                tipText='Location where file is saved on disk')
    self.file_entry.grid(row=row, column=1, sticky='ew')

    row += 1
    label = Label(self, text='Title:')
    label.grid(row=row, column=0, sticky='e')
    self.title_entry = Entry(self, width=40, text=self.title,
                    tipText='Title of the printout, displayed at top')
    self.title_entry.grid(row=row, column=1, sticky='ew')

    row += 1
    label = Label(self, text='X axis label:')
    label.grid(row=row, column=0, sticky='e')
    self.x_axis_entry = Entry(self, width=40, text=self.x_axis_label,
                    tipText='X axis label for the printout')
    self.x_axis_entry.grid(row=row, column=1, sticky='ew')

    row += 1
    label = Label(self, text='Y axis label:')
    label.grid(row=row, column=0, sticky='e')
    self.y_axis_entry = Entry(self, width=40, text=self.y_axis_label,
                    tipText='Y axis label for the printout')
    self.y_axis_entry.grid(row=row, column=1, sticky='ew')

    row += 1
    frame = Frame(self)
    frame.grid(row=row, column=0, columnspan=2, sticky='ew')
    frame.grid_columnconfigure(4, weight=1)
    
    label = Label(frame, text='Paper size:')
    label.grid(row=0, column=0, sticky='e')
    entries = []
    for t in paper_types:
      if t == Output.other_paper_type:
        entry = t
      else:
        (w, h, u) = paper_sizes[t]
        entry = t + ' (%2.1f %s x %2.1f %s)' % (w, u, h, u)
      entries.append(entry)
    self.size_menu = PulldownList(frame, callback=self.changedSize,
                                  texts=entries, index=size_index,
                                  tipText='The paper size for the printout')
    self.size_menu.grid(row=0, column=1, sticky='w')

    self.other_frame = Frame(frame)

    self.other_frame.grid_columnconfigure(0, weight=1)
    self.other_entry = FloatEntry(self.other_frame, text=self.other_size,
                                  isArray=True,
                                  tipText='The size of the Other paper in both dimensions; this requires two values, space or comma separated')
    self.other_entry.grid(row=0, column=0, sticky='ew')
    self.other_unit_menu= PulldownList(self.other_frame, texts=paper_units,
                                       index=other_unit_index,
                                       tipText='The unit for the Other paper size')
    self.other_unit_menu.grid(row=0, column=1, sticky='ew')

    row += 1
    frame = Frame(self)
    frame.grid(row=row, column=0, columnspan=4, sticky='ew')
    frame.grid_columnconfigure(1, weight=1)
    frame.grid_columnconfigure(3, weight=1)
    frame.grid_columnconfigure(5, weight=1)

    label = Label(frame, text='Orientation:')
    label.grid(row=0, column=0, sticky='e')
    self.orientation_menu = PulldownList(frame, texts=paper_orientations,
                                         index=orientation_index,
                                         tipText='Whether the paper should be set in Portrait or Landscape mode')
    self.orientation_menu.grid(row=0, column=1, sticky='w')

    label = Label(frame, text='  Style:')
    label.grid(row=0, column=2, sticky='e')
    self.style_menu = PulldownList(frame, texts=style_choices,
                                   index=style_index,
                                   tipText='Whether the printout should be in colour or black and white')
    self.style_menu.grid(row=0, column=3, sticky='w')

    label = Label(frame, text='  Format:')
    label.grid(row=0, column=4, sticky='e')
    self.format_menu = PulldownList(frame, callback=self.changedFormat,
                                    texts=format_choices, index=format_index,
                                    tipText='Whether to save as PS, EPS or PDF')

    self.format_menu.grid(row=0, column=5, sticky='w')

    if haveTicks:

      row += 1
      frame = Frame(self)
      frame.grid(row=row, column=0, columnspan=4, sticky='ew')
      frame.grid_columnconfigure(1, weight=1)
      frame.grid_columnconfigure(3, weight=1)
      
      label = Label(frame, text='Tick Location:')
      label.grid(row=0, column=0, sticky='e')
      self.tick_menu = PulldownList(frame, texts=tick_locations,
                                    index=tick_location_index,
                                    tipText='Whether the tick marks appear on the inside or outside of the frame')
      self.tick_menu.grid(row=0, column=1, sticky='w')

      label = Label(frame, text='  Tick Placement:')
      label.grid(row=0, column=2, sticky='e')
      if self.tick_placement is None:
        selected = None
      else:
        selected = [(x in self.tick_placement) for x in tick_placements]
      self.tick_buttons = CheckButtons(frame, entries=tick_placements,
                                       selected=selected,
                                       tipTexts=('Whether the tick marks appear on the top and/or bottom and/or left and/or right',))
      self.tick_buttons.grid(row=0, column=3, sticky='w')

      row += 1
      frame = Frame(self)
      frame.grid(row=row, column=0, columnspan=4, sticky='ew')
      frame.grid_columnconfigure(1, weight=1)
      frame.grid_columnconfigure(3, weight=1)
      
      label = Label(frame, text='Tick Font:')
      label.grid(row=0, column=0, sticky='e')
      self.tick_font_list = FontList(frame, mode='Print',
                            selected=self.tick_font,
                            extraTexts=[PrintTicks.no_tick_text],
                            tipText='The font used for the tick mark labels')
      self.tick_font_list.grid(row=0, column=1, sticky='w')

      label = Label(frame, text='Tick Spacing:')
      label.grid(row=0, column=2, sticky='e')
      # TBD: put preferred choice in data model
      self.spacing_menu = PulldownList(frame, texts=spacing_choices, index=0,
                                       callback=self.changedSpacing,
                                       tipText='Whether the program should automatically calculate the major/minor tick spacings and how many decimal places are used for the ticks, or whether the these are specified manually')
      self.spacing_menu.grid(row=0, column=3, sticky='w')

      ff = self.spacing_frame = Frame(frame)
      ff.grid_columnconfigure(1, weight=1)
      ff.grid_columnconfigure(2, weight=1)
      label = Label(ff, text='Tick Spacing')
      label.grid(row=0, column=0, sticky='w')
      label = Label(ff, text='Major')
      label.grid(row=0, column=1, sticky='ew')
      label = Label(ff, text='Minor')
      label.grid(row=0, column=2, sticky='ew')
      label = Label(ff, text='Decimals')
      label.grid(row=0, column=3, sticky='ew')
      label = Label(ff, text='X:')
      label.grid(row=1, column=0, sticky='w')
      self.x_major_entry = FloatEntry(ff, tipText='The spacing in display units of the major tick marks in the X dimension')
      self.x_major_entry.grid(row=1, column=1, sticky='ew')
      self.x_minor_entry = FloatEntry(ff, tipText='The spacing in display units of the minor tick marks in the X dimension (not printed if left blank)')
      self.x_minor_entry.grid(row=1, column=2, sticky='ew')
      self.x_decimal_entry = IntEntry(ff, tipText='The number of decimal places for the tick numbers in the X dimension')
      self.x_decimal_entry.grid(row=1, column=3, sticky='ew')
      label = Label(ff, text='Y:')
      label.grid(row=2, column=0, sticky='w')
      self.y_major_entry = FloatEntry(ff, tipText='The spacing in display units of the major tick marks in the Y dimension')
      self.y_major_entry.grid(row=2, column=1, sticky='ew')
      self.y_minor_entry = FloatEntry(ff, tipText='The spacing in display units of the minor tick marks in the Y dimension (not printed if left blank)')
      self.y_minor_entry.grid(row=2, column=2, sticky='ew')
      self.y_decimal_entry = IntEntry(ff, tipText='The number of decimal places for the tick numbers in the Y dimension')
      self.y_decimal_entry.grid(row=2, column=3, sticky='ew')

      row += 1
      frame = Frame(self)
      frame.grid(row=row, column=0, columnspan=4, sticky='ew')
      frame.grid_columnconfigure(1, weight=1)
      
      label = Label(frame, text='Tick Length:')
      label.grid(row=0, column=0, sticky='e')
      # TBD: put preferred choice in data model
      self.tick_length_menu = PulldownList(frame, texts=tick_length_choices, index=0,
                                       callback=self.changedLength,
                                       tipText='Whether the program should automatically calculate the major/minor tick lengths, or whether the these are specified manually')
      self.tick_length_menu.grid(row=0, column=1, sticky='w')

      ff = self.length_frame = Frame(frame)
      ff.grid_columnconfigure(1, weight=1)
      label = Label(ff, text='  Major length:')
      label.grid(row=0, column=0, sticky='w')
      self.length_major_entry = FloatEntry(ff, tipText='The length in points of the major tick marks')
      self.length_major_entry.grid(row=0, column=1, sticky='w')
      label = Label(ff, text='Minor length:')
      label.grid(row=0, column=2, sticky='w')
      self.length_minor_entry = FloatEntry(ff, tipText='The length in points of the minor tick marks')
      self.length_minor_entry.grid(row=0, column=3, sticky='w')

    row += 1
    frame = Frame(self)
    frame.grid(row=row, column=0, columnspan=4, sticky='ew')
    frame.grid_columnconfigure(3, weight=1)
    frame.grid_columnconfigure(4, weight=1)

    label = Label(frame, text='Scaling:')
    label.grid(row=0, column=0, sticky='e')
    # TBD: put preferred choice in data model
    self.scaling_menu = PulldownList(frame, texts=scaling_choices, index=0,
                                     callback=self.changedScaling,
                                     tipText='Whether the plot should be scaled as a percentage of the maximum size that would fit on the paper, or instead should be specified by the number of cms or inches per unit')
    self.scaling_menu.grid(row=0, column=1, sticky='ew')

    self.scaling_scale = Scale(frame, orient=Tkinter.HORIZONTAL, value=self.scaling, tipText='The percentage of the maximum size that would fit on the paper that the plot is scaled by')
    self.scaling_scale.grid(row=0, column=2, columnspan=3, sticky='ew')

    self.x_scaling_label = Label(frame, text='X:')
    self.x_scaling_entry = FloatEntry(frame, tipText='The scaling that should be used in the X dimension as cms or inches per unit')
    self.y_scaling_label = Label(frame, text='Y:')
    self.y_scaling_entry = FloatEntry(frame, tipText='The scaling that should be used in the Y dimension as cms or inches per unit')

    row += 1
    frame = Frame(self)
    frame.grid(row=row, column=0, columnspan=4, sticky='w')
    frame.grid_columnconfigure(2, weight=1)
  
    label = Label(frame, text='Include:')
    label.grid(row=0, column=0, sticky='e')
    tipTexts = ('Whether the time and date should be included in the printout',
                'Whether the file name should be included in the printout')
    if self.border_decoration is None:
      selected = None
    else:
      selected = [(x in self.border_decoration) for x in border_decorations]
    self.border_buttons = CheckButtons(frame, entries=border_decorations,
                                       selected=selected,
                                       tipTexts=tipTexts)
    self.border_buttons.grid(row=0, column=1, sticky='w')

    label = Label(frame, text='    Using Font:')
    label.grid(row=0, column=2, sticky='e')
    self.border_font_list = FontList(frame, mode='Print',
                          selected=self.border_font,
                          tipText='The font used for the border texts')
    self.border_font_list.grid(row=0, column=3, sticky='w')

    row += 1
    label = Label(self, text='Line width:')
    label.grid(row=row, column=0, sticky='w')
    self.linewidth_entry = FloatEntry(self, width=10,
                                  text=self.linewidth,
                                  tipText='Line width for drawing')
    self.linewidth_entry.grid(row=row, column=1, sticky='w')
 
  def destroy(self):

    self.setOptionValues()

    if self.file_select_popup:
      self.file_select_popup.destroy()

    Frame.destroy(self)

  def getOptionValues(self):

    getOption = self.getOption
    if getOption:

      file_name = getOption('FileName', defaultValue='')
      title = getOption('Title', defaultValue='')
      x_axis_label = getOption('XAxisLabel', defaultValue='')
      y_axis_label = getOption('YAxisLabel', defaultValue='')
      paper_type = getOption('PaperSize', defaultValue=paper_types[0])
      paper_type = paper_type_dict.get(paper_type, paper_types[0])
      other_height = getOption('OtherHeight', defaultValue=10)
      other_width = getOption('OtherWidth', defaultValue=10)
      other_size = [other_height, other_width]
      other_unit = getOption('OtherUnit', defaultValue=paper_units[0])
      paper_orientation = getOption('Orientation', defaultValue=paper_orientations[0])
      in_color = getOption('InColor', defaultValue=True)
      if in_color:
        output_style = style_choices[0]
      else:
        output_style = style_choices[1]
      format_option = getOption('OutputFormat', defaultValue=format_options[0])
      output_format = format_choices[format_options.index(format_option)]
      if self.haveTicks:
        tick_outside = getOption('TickOutside', defaultValue=tick_locations[0])
        if tick_outside:
          tick_location = tick_locations.index(PrintTicks.Outside)
        else:
          tick_location = tick_locations.index(PrintTicks.Inside)
        tick_placement = getTickPlacement1(getOption('TickPlacement', defaultValue='nsew'))
      dateTime = getOption('ShowsDateTime', defaultValue=True)
      fileName = getOption('ShowsFileName', defaultValue=True)
      border_font = getOption('BorderFont', defaultValue='Helvetica 10')
      border_decoration = []
      if dateTime:
        border_decoration.append(border_decorations[0])
      if fileName:
        border_decoration.append(border_decorations[1])

      if self.haveTicks:
        spacing_choice = getOption('SpacingChoice', defaultValue=spacing_choices[0])
        x_major = getOption('XMajor', defaultValue=1.0)
        x_minor = getOption('XMinor', defaultValue=1.0)
        x_decimal = getOption('XDecimal', defaultValue=3)
        y_major = getOption('YMajor', defaultValue=1.0)
        y_minor = getOption('YMinor', defaultValue=1.0)
        y_decimal = getOption('YDecimal', defaultValue=3)

        tick_length_choice = getOption('TickLengthChoice', defaultValue=tick_length_choices[0])
        tick_major = getOption('TickMajor', defaultValue=10)
        tick_minor = getOption('TickMinor', defaultValue=5)

      scaling_choice = getOption('ScalingChoice', defaultValue=scaling_choices[0])
      scaling = getOption('Scaling', defaultValue=0.7)
      scaling = int(round(100.0 * scaling))
      x_scaling = getOption('XScaling', defaultValue=1.0)
      y_scaling = getOption('YScaling', defaultValue=1.0)

      if self.haveTicks:
        tick_font = getOption('TickFont', defaultValue='Helvetica 10')

      linewidth = getOption('LineWidth', defaultValue=Output.default_linewidth)

    else:
      file_name = ''
      title = ''
      x_axis_label = ''
      y_axis_label = ''
      paper_type = paper_types[0]
      other_unit = paper_units[0]
      other_size = ''
      paper_orientation = paper_orientations[0]
      output_style = style_choices[0]
      output_format = format_choices[0]
      if self.haveTicks:
        tick_location = tick_locations[0]
        tick_placement = tick_placements
      border_decoration = border_decorations
      border_font = 'Helvetica 10'
      if self.haveTicks:
        spacing_choice = spacing_choices[0]
        x_major = 1.0
        x_minor = 1.0
        x_decimal = 3
        y_major = 1.0
        y_minor = 1.0
        y_decimal = 3
        tick_length_choice = tick_length_choices[0]
        tick_major = 10
        tick_minor = 5
      scaling_choice = scaling_choices[0]
      scaling = 70
      x_scaling = 1.0
      y_scaling = 1.0
      if self.haveTicks:
        tick_font = 'Helvetica 10'
      linewidth = Output.default_linewidth

    if not self.haveTicks:
      tick_location = None
      tick_placement = None
      spacing_choice = spacing_choices[0]
      x_major = 1.0
      x_minor = 1.0
      x_decimal = 3
      y_major = 1.0
      y_minor = 1.0
      y_decimal = 3
      tick_font = 'Helvetica 10'
      tick_length_choice = tick_length_choices[0]
      tick_major = 10
      tick_minor = 5

    self.file_name = file_name
    self.title = title
    self.x_axis_label = x_axis_label
    self.y_axis_label = y_axis_label
    self.paper_type = paper_type
    self.other_unit = other_unit
    self.other_size = other_size
    self.paper_orientation = paper_orientation
    self.output_style = output_style
    self.output_format = output_format
    self.tick_location = tick_location
    self.tick_placement = tick_placement
    self.border_decoration = border_decoration
    self.border_font = border_font

    self.spacing_choice = spacing_choice
    self.x_major = x_major
    self.x_minor = x_minor
    self.x_decimal = x_decimal
    self.y_major = y_major
    self.y_minor = y_minor
    self.y_decimal = y_decimal

    self.scaling_choice = scaling_choices[0]
    self.scaling = scaling
    self.x_scaling = x_scaling
    self.y_scaling = y_scaling

    self.tick_font = tick_font
    self.linewidth = linewidth

    self.tick_length_choice = tick_length_choice
    self.tick_major = tick_major
    self.tick_minor = tick_minor

  def setOptionValues(self):

    if not hasattr(self, 'file_entry'):
      # it looks like on destroy can have function called but file_entry deleted already
      return

    self.file_name = file_name = self.file_entry.get()
    self.title = title = self.title_entry.get()
    self.x_axis_label = x_axis_label = self.x_axis_entry.get()
    self.y_axis_label = y_axis_label = self.y_axis_entry.get()

    n = self.size_menu.getSelectedIndex()
    self.paper_type = paper_type = paper_types[n]
    if paper_type == Output.other_paper_type:
      other_size = self.other_entry.get()
      other_unit = self.other_unit_menu.getText()
    else:
      other_size = None
      other_unit = None
    self.other_size = other_size
    self.other_unit = other_unit

    self.paper_orientation = paper_orientation = self.orientation_menu.getText()
    self.output_style = output_style = self.style_menu.getText()
    self.output_format = output_format = self.format_menu.getText()

    if self.haveTicks:
      tick_location = self.tick_menu.getText()
      tick_placement = self.tick_buttons.getSelected()
    else:
      tick_location = tick_placement = None
    self.tick_location = tick_location
    self.tick_placement = tick_placement

    self.border_decoration = border_decoration = self.border_buttons.getSelected()
    self.border_font = border_font = self.border_font_list.getText()

    if self.haveTicks:
      self.spacing_choice = spacing_choice = self.spacing_menu.getText()
      if spacing_choice != spacing_choices[0]:
        self.x_major = self.x_major_entry.get()
        self.x_minor = self.x_minor_entry.get()
        self.x_decimal = self.x_decimal_entry.get()
        self.y_major = self.y_major_entry.get()
        self.y_minor = self.y_minor_entry.get()
        self.y_decimal = self.y_decimal_entry.get()

      self.tick_length_choice = tick_length_choice = self.tick_length_menu.getText()
      if tick_length_choice != tick_length_choices[0]:
        self.tick_major = self.length_major_entry.get()
        self.tick_minor = self.length_minor_entry.get()

    self.scaling_choice = scaling_choice = self.scaling_menu.getText()
    if self.scaling_choice == scaling_choices[0]:
      scaling = self.scaling_scale.get()
      self.scaling = int(round(scaling))
    else:
      self.x_scaling = self.x_scaling_entry.get()
      self.y_scaling = self.y_scaling_entry.get()

    if self.haveTicks:
      self.tick_font = self.tick_font_list.getText()

    self.linewidth = self.linewidth_entry.get()

    setOption = self.setOption
    if setOption:
      setOption('FileName', value=file_name)
      setOption('Title', value=title)
      setOption('XAxisLabel', value=x_axis_label)
      setOption('YAxisLabel', value=y_axis_label)
      if paper_type == Output.other_paper_type:
        setOption('OtherHeight', value=other_size[0])
        setOption('OtherWidth', value=other_size[1])
        setOption('OtherUnit', value=other_unit)
      else:
        paper_type = paper_type_inverse_dict[paper_type]
        setOption('PaperSize', value=paper_type)
      setOption('Orientation', value=paper_orientation)
      in_color = (output_style == style_choices[0])
      setOption('InColor', value=in_color)
      output_format = format_options[format_choices.index(output_format)]
      setOption('OutputFormat', value=output_format)
      if self.haveTicks:
        tick_outside = (tick_location == PrintTicks.Outside)
        setOption('TickOutside', value=tick_outside)
        tick_placement = getTickPlacement2(tick_placement)
        setOption('TickPlacement', value=tick_placement)
      dateTime = (border_decorations[0] in border_decoration)
      fileName = (border_decorations[1] in border_decoration)
      setOption('ShowsDateTime', value=dateTime)
      setOption('ShowsFileName', value=fileName)
      setOption('BorderFont', value=border_font)

      if self.haveTicks:
        setOption('SpacingChoice', value=spacing_choice)
        if spacing_choice != spacing_choices[0]:
          setOption('XMajor', self.x_major)
          setOption('XMinor', self.x_minor)
          setOption('XDecimal', self.x_decimal)
          setOption('YMajor', self.y_major)
          setOption('YMinor', self.y_minor)
          setOption('YDecimal', self.y_decimal)

        setOption('TickLengthChoice', value=tick_length_choice)
        if tick_length_choice != tick_length_choices[0]:
          setOption('TickMajor', self.tick_major)
          setOption('TickMinor', self.tick_minor)

      setOption('ScalingChoice', value=scaling_choice)
      if scaling_choice == scaling_choices[0]:
        setOption('Scaling', value=0.01*self.scaling)
      else:
        setOption('XScaling', value=self.x_scaling)
        setOption('YScaling', value=self.y_scaling)

      if self.haveTicks:
        setOption('TickFont', self.tick_font)

      setOption('LineWidth', self.linewidth)

  def findFile(self):

    if self.file_select_popup:
      self.file_select_popup.open()
    else:
      file_types = [ FileType('All', ['*']),
                     FileType('PostScript', ['*.ps', '*.eps']),
                     FileType('PDF', ['*.pdf', '*.ai']) ]
      self.file_select_popup = FileSelectPopup(self, file_types=file_types)

    file = self.file_select_popup.getFile()
    if file:
      self.file_entry.set(file)

  def changedSize(self, entry):

    if entry == Output.other_paper_type:
      self.other_frame.grid(row=0, column=2, columnspan=2, sticky='w')
    else:
      self.other_frame.grid_forget()

  def changedFormat(self, entry):

    file_suffix = file_suffixes.get(entry)
    if not file_suffix:
      return

    file_name = self.file_entry.get()
    if not file_name:
      return

    for suffix in format_suffixes:
      if file_name.endswith(suffix):
        if suffix != file_suffix:
          n = len(suffix)
          file_name = file_name[:-n] + file_suffix
          self.file_entry.set(file_name)
        break
    else:
      file_name = file_name + file_suffix
      self.file_entry.set(file_name)

  def changedScaling(self, choice):

    if choice == scaling_choices[0]:
      self.scaling_scale.grid(row=0, column=2, columnspan=3, sticky='ew')
      self.x_scaling_label.grid_forget()
      self.x_scaling_entry.grid_forget()
      self.y_scaling_label.grid_forget()
      self.y_scaling_entry.grid_forget()
    else:
      self.scaling_scale.grid_forget()
      self.x_scaling_label.grid(row=0, column=2, sticky='w')
      self.x_scaling_entry.grid(row=0, column=3, columnspan=2, sticky='ew')
      self.y_scaling_label.grid(row=1, column=2, sticky='w')
      self.y_scaling_entry.grid(row=1, column=3, columnspan=2, sticky='ew')

    self.setOptionValues()

  def changedSpacing(self, choice):

    if choice == spacing_choices[0]:
      self.spacing_frame.grid_forget()
    else:
      self.spacing_frame.grid(row=1, column=1, columnspan=3, sticky='ew')

    self.setOptionValues()

  def changedLength(self, choice):

    if choice == tick_length_choices[0]:
      self.length_frame.grid_forget()
    else:
      self.length_frame.grid(row=1, column=0, columnspan=4, sticky='ew')

    self.setOptionValues()

  # width and height are of plot, in pixels
  def getOutputHandler(self, pixel_width, pixel_height, unit_width=1.0, unit_height=1.0, fonts=None):

    if not fonts:
      fonts = []
    else:
      fonts = list(fonts)

    for n in range(len(fonts)):
      if fonts[n] == 'Times':
        fonts[n] = 'Times-Roman'

    self.setOptionValues()

    if not self.file_name:
      showError('No file', 'No file specified', parent=self)
      return None

    x_scaling = y_scaling = 1.0
    if self.scaling_choice != scaling_choices[0]:
      try:
        x_scaling = float(self.x_scaling)
      except:
        showError('Bad X Scaling', 'Specified X Scaling must be floating point', parent=self)
        return None
      try:
        y_scaling = float(self.y_scaling)
      except:
        showError('Bad Y Scaling', 'Specified Y Scaling must be floating point', parent=self)
        return None

    if os.path.exists(self.file_name):
      if not showYesNo('File exists', 'File "%s" exists, overwrite?' % self.file_name, parent=self):
        return None

    if self.paper_type == Output.other_paper_type:
      paper_size = self.other_size + [ self.other_unit ]
    else:
      paper_size = paper_sizes[self.paper_type]

    output_scaling = self.scaling / 100.0

    border_font = self.border_font
    (font, size) = border_font.split()
    size = int(size)
    border_text = {}
    for decoration in self.border_decoration:
      if decoration == 'Time and Date':
        location = 'se'
        text = time.ctime(time.time())
      elif decoration == 'File Name':
        location = 'sw'
        text = self.file_name
      else:
        continue # should not be here
      border_text[location] = (text, font, size)
    if self.title:
      location = 'n'
      border_text[location] = (self.title, font, size+6)

    if font not in fonts:
      fonts.append(font)

    if self.haveTicks and self.tick_location == PrintTicks.Outside:
      axis_label_offset = 2
    else:
      axis_label_offset = 0
     
    outputHandler = PrintHandler.getOutputHandler(self.file_name,
                              pixel_width, pixel_height,
                              unit_width, unit_height,
                              scaling_choice=self.scaling_choice,
                              output_scaling=output_scaling,
                              w_scaling=x_scaling,
                              h_scaling=y_scaling,
                              paper_size=paper_size,
                              paper_orientation=self.paper_orientation,
                              output_style=self.output_style,
                              output_format=self.output_format,
                              border_text=border_text,
                              x_axis_label=self.x_axis_label,
                              y_axis_label=self.y_axis_label,
                              axis_label_font=font,
                              axis_label_size=size,
                              axis_label_offset=axis_label_offset,
                              fonts=fonts,
                              linewidth=self.linewidth,
                              do_outline_box=self.doOutlineBox)

    return outputHandler

  def getAspectRatio(self):

    self.setOptionValues()

    if self.paper_type == Output.other_paper_type:
      paper_size = self.other_size
    else:
      paper_size = paper_sizes[self.paper_type]

    r = paper_size[1] / paper_size[0]
    if self.paper_orientation == 'Landscape':
      r = 1.0 / r

    return r
