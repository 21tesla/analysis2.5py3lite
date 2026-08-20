"""
======================COPYRIGHT/LICENSE START==========================

Region2D.py: Utility code for CCPN code generation framework

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
from memops.universal.Region1D import Region1D

class Region2D:

  def __init__(self, x_region = None, y_region = None):

    if (x_region is None):
      x_region = Region1D()

    if (y_region is None):
      y_region = Region1D()

    self.x_region = x_region
    self.y_region = y_region

  def copy(self):

    return Region2D(self.x_region, self.y_region)

  def get(self):

    return (self.x_region[0], self.y_region[0], self.x_region[1], self.y_region[1])

  def __getitem__(self, key):

    try:
      n = int(key)
    except:
      raise TypeError, 'sequence index must be integer'

    if (n == 0):
      return self.x_region[0]
    elif (n == 1):
      return self.y_region[0]
    elif (n == 2):
      return self.x_region[1]
    elif (n == 3):
      return self.y_region[1]

    raise IndexError, 'list index out of range'

