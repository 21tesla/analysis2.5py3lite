
"""
======================COPYRIGHT/LICENSE START==========================

Frame.py: <write function here>

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

from memops.gui.Base import Base

class Frame(Tkinter.Frame, Base):

  def __init__(self, parent, grid=None, gridSpan=(1,1), expandGrid=None,
               sticky='nsew', docKey=None, tipText=None, createToolTip=False,
               *args, **kw):

    if not kw.get('bg'):
      kw['bg'] = 'grey90'
      
    Tkinter.Frame.__init__(self, parent, *args, **kw)
    Base.__init__(self, docKey=docKey, tipText=tipText, createToolTip=createToolTip)

    self.parent = parent
    
    if grid is not None:
      row, col = grid
      rowSpan, colSpan = gridSpan
      self.grid(row=row, column=col, rowspan=rowSpan,
                columnspan=colSpan, sticky=sticky)

    if expandGrid:
      if type(expandGrid[0]) is type(1):
        row, col = expandGrid
        self.expandGrid(row, col)
      
      else:
        for row, col in expandGrid:
          self.expandGrid(row, col)

  def expandGrid(self, row, col):
  
    if row is not None:
      self.grid_rowconfigure(row, weight=1)
    
    if col is not None:
      self.grid_columnconfigure(col, weight=1 )
    

  def open(self):

    if (hasattr(self.parent, 'open')): # intended use: parent is a BasePopup
      self.parent.open()

  def close(self):

    if (hasattr(self.parent, 'close')): # intended use: parent is a BasePopup
      self.parent.close()

