
"""
======================COPYRIGHT/LICENSE START==========================

RegionScrollbar.py: <write function here>

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

from memops.universal.Region1D import Region1D
from memops.universal.Util import formatDecimals

from memops.gui.Scrollbar import Scrollbar

class RegionScrollbar(Scrollbar):

  def __init__(self, parent, world_region = None, view_region = None,
               callback = None, text_decimals = 2, min_thickness = None,
               *args, **kw):
 
    if (world_region is None):
      world_region = Region1D()
 
    if (view_region is None):
      view_region = Region1D()

    (w0, w1) = world_region
    assert w0 != w1, 'w0 = w1 = %s' % w0
 
    (v0, v1) = view_region
    assert v0 != v1, 'v0 = v1 = %s' % v0
 
    if(w0 < w1):
      assert w0 <= v0 and v0 < v1 and v1 <= w1, 'w0 = %s, w1 = %s, v0 = %s, v1 = %s' % (w0, w1, v0, v1)
    else:
      assert w0 >= v0 and v0 > v1 and v1 >= w1, 'w0 = %s, w1 = %s, v0 = %s, v1 = %s' % (w0, w1, v0, v1)

    self.world_region = world_region
    self.view_region = view_region
    self.myCallback = callback
    self.text_decimals = text_decimals

    kw['text_func'] = self.determineText
    Scrollbar.__init__(self, parent, *args, **kw)

    self.setMinThickness(min_thickness)
    self.setViewRegion(v0, v1, do_callback=True)
 
  def setMinThickness(self, min_thickness):

    if (min_thickness):
      (w0, w1) = self.world_region
      self.min_thickness = min_thickness / abs(w1 - w0)
    else:
      self.min_thickness = None

  def getViewRegion(self):

    return self.view_region

  def setViewRegion(self, v0, v1, do_callback = False):

    (w0, w1) = self.world_region
    if (w0 < w1):
      if (v0 < w0):
        v1 = w0 + v1 - v0
        v0 = w0
        assert v1 <= w1
      elif (v1 > w1):
        v0 = w1 + v0 - v1
        v1 = w1
        assert v0 >= w0
    else:
      if (v0 > w0):
        v1 = w0 + v1 - v0
        v0 = w0
        assert v1 >= w1
      elif (v1 < w1):
        v0 = w1 + v0 - v1
        v1 = w1
        assert v0 <= w0
 
    lo = float(v0-w0) / (w1-w0)
    hi = float(v1-w0) / (w1-w0)

    self.set(lo, hi, do_callback=do_callback)

  def setWorldRegion(self, world_region):

    (v0, v1) = self.view_region
    (w0, w1) = world_region
    (ww0, ww1) = self.world_region
    if (w0 == ww0 and w1 == ww1):
      return
    if (w0 < w1):
      w0 = min(w0, v0)
      w1 = max(w1, v1)
    else:
      w0 = max(w0, v0)
      w1 = min(w1, v1)
    self.world_region = Region1D(w0, w1)
    self.setViewRegion(v0, v1, do_callback=True)

  def set(self, lo, hi, do_callback = False):

    #print 'set1', lo, hi
    Scrollbar.set(self, lo, hi, do_callback=False)
    #print 'set2', self.lo, self.hi
    (w0, w1) = self.world_region
    v0 = w0 + self.lo*(w1-w0)
    v1 = w0 + self.hi*(w1-w0)
    #print 'set3', v0, v1
    self.view_region.set(v0, v1)
    # updateText repeats what Scrollbar.set does but needs to because region updated after set
    self.updateText()
    if (do_callback):
      self.doCallback(self.lo, self.hi)

  def doCallback(self, lo, hi):

    if (self.myCallback):
      self.myCallback(self.view_region)
    
  def determineText(self):

    (v0, v1) = self.view_region
    t0 = formatDecimals(v0, self.text_decimals)
    t1 = formatDecimals(v1, self.text_decimals)

    return (t0, t1)

if __name__ == '__main__':
 
  def callback(view_region):
    print 'callback', view_region

  root = Tkinter.Tk()
 
  wr = Region1D(-0.5, 1.5)
  vr = Region1D(0.2, 0.5)
 
  r = RegionScrollbar(root, world_region=wr, view_region=vr, callback=callback,
                    allow_resize=True, show_text=True)
  r.pack(fill=Tkinter.BOTH, expand=Tkinter.YES)

  def update_world_region():

    (w0, w1) = r.world_region
    w0 = w0 - 1
    w1 = w1 + 1
    print 'update_world_region1', w0, w1
    wr = Region1D(w0, w1)
    print 'update_world_region2', wr
    r.setWorldRegion(wr)

  b = Tkinter.Button(root, text='update world', command=update_world_region)
  b.pack()
 
  root.mainloop()
