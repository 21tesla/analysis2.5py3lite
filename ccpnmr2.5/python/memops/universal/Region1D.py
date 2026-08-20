"""
======================COPYRIGHT/LICENSE START==========================

Region1D.py: Utility code for CCPN code generation framework

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
class Region1D:

  def __init__(self, t0 = 0.0, t1 = 1.0):

    self.set(t0, t1)

  def copy(self):

    return Region1D(self.t0, self.t1)

  def get(self):

    return (self.t0, self.t1)

  def set(self, t0, t1):

    self.t0 = t0
    self.t1 = t1

  def fitInMaxRegion(self, t0, t1, max_region):

    (r0, r1) = max_region
    #print 'fitInMaxRegion1', t0, t1, r0, r1

    if (r0 < r1):
      #t0 = max(t0, r0)
      #t1 = min(t1, r1)
      # if bounds reached try and translate if possible
      if (t0 < r0):
        t1 = min(t1+r0-t0, r1)
        t0 = r0
      elif (t1 > r1):
        t0 = max(t0+r1-t1, r0)
        t1 = r1
    else:
      #t0 = min(t0, r0)
      #t1 = max(t1, r1)
      # if bounds reached try and translate if possible
      if (t0 > r0):
        t1 = max(t1+r0-t0, r1)
        t0 = r0
      elif (t1 < r1):
        t0 = min(t0+r1-t1, r0)
        t1 = r1

    #print 'fitInMaxRegion2', t0, t1, r0, r1
    return (t0, t1)

  def translateFraction(self, f, max_region = None):

    if (f == 0.0):
      return

    d = f * (self.t1 - self.t0)
    t0 = self.t0 + d
    t1 = self.t1 + d

    if (max_region):
      (r0, r1) = max_region
      if (r0 < r1):
        if (t0 < r0):
          t1 = t1 - t0 + r0
          t0 = r0
        elif (t1 > r1):
          t0 = t0 - t1 + r1
          t1 = r1
      else:
        if (t0 > r0):
          t1 = t1 - t0 + r0
          t0 = r0
        elif (t1 < r1):
          t0 = t0 - t1 + r1
          t1 = r1
      # below extra safety
      (t0, t1) = self.fitInMaxRegion(t0, t1, max_region)

    self.t0 = t0
    self.t1 = t1

  def zoom(self, scale, max_region = None, max_extent = None):

    u = 0.5 * (self.t0 + self.t1)
    v = 0.5 * scale * (- self.t0 + self.t1)
    #print 'zoom1', v, scale, self.t0, self.t1
    if (max_extent):
      if (v > 0):
        v = min(v, 0.5*max_extent)
      else:
        v = max(v, -0.5*max_extent)
      #print 'zoom2', v, max_extent

    t0 = u - v
    t1 = u + v

    if (max_region):
      (t0, t1) = self.fitInMaxRegion(t0, t1, max_region)

    self.t0 = t0
    self.t1 = t1

  def resizeFraction(self, f, resize_type, aspect_ratio = 1.0,
                     other_region = None, other_max_region = None,
                     other_max_extent = None):

    #print 'resizeFraction', f, resize_type, self.t0, self.t1

    if (f == 0.0):
      return

    d = f * (self.t1 - self.t0)
    if (resize_type == 0):
      t0 = self.t0 + d
      if (self.t0 < self.t1):
        if (t0 >= self.t1):
          return
      else:
        if (t0 <= self.t1):
          return
      self.t0 = t0
    else:
      t1 = self.t1 + d
      if (self.t0 < self.t1):
        if (t1 <= self.t0):
          return
      else:
        if (t1 >= self.t0):
          return
      self.t1 = t1

    if (other_region):
      (s0, s1) = other_region
      scale = abs((self.t1 - self.t0) * aspect_ratio / (s1 - s0))
      other_region.zoom(scale, other_max_region, other_max_extent)

  def size(self):

    return self.t1 - self.t0

  def __getitem__(self, key):

    try:
      n = int(key)
    except:
      raise TypeError, 'sequence index must be integer'

    if (n == 0):
      return self.t0
    elif (n == 1):
      return self.t1

    raise IndexError, 'list index out of range'

  def __setitem__(self, key, value):

    try:
      n = int(key)
    except:
      raise TypeError, 'sequence index must be integer'

    if (n == 0):
      self.t0 = float(value)
    elif (n == 1):
      self.t1 = float(value)
    else:
      raise IndexError, 'list index out of range'

  def __repr__(self):

    return '(' + str(self.t0) + ', ' + str(self.t1) + ')'
