"""
======================COPYRIGHT/LICENSE START==========================

Converters.py: Data compatibility handling

Copyright (C) 2007 Rasmus Fogh (CCPN project)
 
=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../../license/LGPL.license.
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

======================COPYRIGHT/LICENSE END============================

To obtain more information about this code:

- CCPN website (http://www.ccpn.ac.uk)

- contact Rasmus Fogh (ccpn@bioc.cam.ac.uk)

=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following reference:

===========================REFERENCE START=============================
Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. Ionides and
Ernest D. Laue (2005). A framework for scientific data modeling and 
automated software development. Bioinformatics 21, 1678-1684.
===========================REFERENCE END===============================
 
"""

# needed  to make it work in Python 2.5
from __future__ import absolute_import

import gc
import time

def doMajorUpgradeToCurrent(oldPath, oldTag, oldName=None, newName=None, 
                            newDir=None, doSave=True):
  """Wrapper to do conversion
  """
  try:
    gc.disable()
    xtime = time.time()
    print('\n+++ Start major upgrade ...')
    
    from memops.general.Constants import currentModelVersion
    curVersion = str(currentModelVersion)
 
    # get old data
    from memops.format.compatibility.part2.Converters2 import getOldData
    oldRoot, oldVersionStr, oldTags = getOldData(oldPath, oldTag, curVersion, 
                                                 oldName=oldName)
    print('+++ Done load old data - time %s' % (time.time() - xtime))

    # make new data
    from memops.format.compatibility.part1.Converters1 import majorUpgradeToCurrent
    result = majorUpgradeToCurrent(oldRoot, oldVersionStr,
                                   newName=newName, newDir=newDir, 
                                   doSave=doSave, oldTags=oldTags)
  finally:
    gc.enable()
    print('+++ Done major upgrade, Total time %s, nObjs %s\n'
          % (time.time() - xtime, len(gc.get_objects())))
  
  return result
