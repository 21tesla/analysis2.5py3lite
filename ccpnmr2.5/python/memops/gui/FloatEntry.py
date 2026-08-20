
"""
======================COPYRIGHT/LICENSE START==========================

FloatEntry.py: <write function here>

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
import re

from memops.universal.Util import formatFloat

from memops.gui.Entry import Entry

class FloatEntry(Entry):

  def __init__(self, parent, cleanFloat=False, formatPlaces = 5, *args, **kw):

    kw['valToStr'] = lambda x: formatFloat(x, places=formatPlaces)

    if cleanFloat:
      kw['strToVal'] = self.cleanFloat
    else:
      kw['strToVal'] = float

    apply(Entry.__init__, (self, parent) + args, kw)

  def cleanFloat(self, val):
  
    return float(re.sub(r'[^-|0-9|\.|,]','',val) or 0) 

if __name__ == '__main__':

  import Tkinter
  from memops.gui.Button import Button

  def func():

    v = entry1.get()
    print type(v), v
    v = entry2.get()
    print type(v), v

  root = Tkinter.Tk()
 
  entry1 = FloatEntry(root)
  entry1.grid()
  entry2 = FloatEntry(root, isArray=True)
  entry2.grid()
  button = Button(root, text='hit me', command=func)
  button.grid()

  root.mainloop()

