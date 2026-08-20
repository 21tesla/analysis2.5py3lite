
"""
======================COPYRIGHT/LICENSE START==========================

PrintPopup.py: <write function here>

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
import Tkinter


from memops.gui.MessageReporter import showError
from memops.gui.PrintFrame import PrintFrame
from memops.gui.Util import createDismissHelpButtonList

from memops.editor.BasePopup import BasePopup

class PrintPopup(BasePopup):

  """
  **Print Window to Output File**

  The purpose of this dialog is to allow the saving of the window
  to a file, in one of the following formats: PostScript (PS),
  Encapsulated PostScript (EPS) or Portable Document Format (PDF).

  At its simplest to print out a window you just need to specify the
  File name, and then click "Save".  But it is likely you will at the
  very least want to change some of the other settings.

  You can specify a Title and a label for the X axis and/or Y axis
  (although these might not make sense dependent on the context).

  You can also specify the Paper size (default A4), the
  Orientation of the paper (default Portrait), whether the printout
  is Color or Black and white (the Style, default Color), and what
  the Format is (PS, EPS or PDF, default PS).

  How the window fits into the paper is determined by the Scaling
  option.  The Percentage option just means that the window is
  scaled by that amount relative to the biggest size it could be and
  still fit on the paper.  The remaining options are if you want to
  specify the cms or inches per unit (normally ppm), or the inverse
  of these.  In this case it could be the case that the window
  actually exceeds the size of the paper.

  You can also include the Time and Date and/or the File Name in the
  printout, and you can specify the font used for this.  The same
  font is used for the Title and X and Y axis labels, except that
  the Title font is 6 pts bigger.
"""

  # width and height are desired plot size in pixels
  def __init__(self, parent, width, height, getOption=None, setOption=None, *args, **kw):
 
    self.outputHandler = None
    self.setPlotSize(width, height)
    self.getOption = getOption
    self.setOption = setOption

    if (not kw.has_key('popup_name')):
      kw['popup_name'] = 'print'

    BasePopup.__init__(self, parent=parent, *args, **kw)

  def body(self, master):

    master.grid_rowconfigure(0, weight=1)
    master.grid_columnconfigure(0, weight=1)

    row = 0

    self.printFrame = PrintFrame(master, getOption=self.getOption, setOption=self.setOption)
    self.printFrame.grid(row=row, column=0, sticky=Tkinter.EW)

    row = row + 1
    texts = [ 'Save' ]
    commands = [ self.ok ]
    buttons = createDismissHelpButtonList(master, texts=texts, commands=commands,
                                          dismiss_text='Cancel')
    buttons.grid(row=row, column=0, sticky=Tkinter.EW)

  # width and height are desired plot size in pixels
  def setPlotSize(self, width, height):

    self.width = width
    self.height = height
 
  def apply(self):

    try:
      self.outputHandler = self.printFrame.getOutputHandler(self.width, self.height)
    except IOError, e:
      showError('IO Error', str(e), parent=self)
      return False
    except:
      showError('Unknown Error', 'unknown error', parent=self)
      return False

    return True
