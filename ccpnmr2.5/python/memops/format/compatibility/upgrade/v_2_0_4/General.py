"""
======================COPYRIGHT/LICENSE START==========================

General.py: Data compatibility handling

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
from memops.format.compatibility.upgrade.v_2_0_5 import General as General205

# functions to convert fullKeys before application
fullKeyConverters = {}
fullKeyConverters.update(General205.fullKeyConverters)

# guids of elements that should be treated as old
elemsTreatedAsOld = set((
)).union(General205.elemsTreatedAsOld)

# pairs of element guids that should be treated as matching, e.g. whe n
# a single element must match with several elements in subclasses
elementPairings = [
] 

def extraMapChanges(globalMapping):
  """ Extra map changes specific for a given step
  """
  General205.extraMapChanges(globalMapping)
  
  for guid in (
               # set transferType attributes to delay, to allow resetting
               "www.ccpn.ac.uk_Fogh_2006-08-16-18:22:58_00022", 
               "www.ccpn.ac.uk_Fogh_2006-08-16-18:20:07_00001",
               # set links to ExpPrototype package to delay, to allow resetting
               "www.ccpn.ac.uk_Fogh_2006-08-16-18:20:06_00008",
               "www.ccpn.ac.uk_Fogh_2006-08-16-18:23:00_00002",
               "www.ccpn.ac.uk_Fogh_2006-08-16-18:20:05_00025",
               ):
    globalMapping['mapsByGuid'][guid]['proc'] = 'delay'
    
  for guid in ('www.ccpn.ac.uk_Fogh_2006-08-17-15:11:12_00001', #Macro.path
               'www.ccpn.ac.uk_Fogh_2006-08-16-18:23:11_00011', # FixedResonance.name
               ):
    dd = globalMapping['mapsByGuid'][guid]
    if 'proc' in dd:
      del dd['proc']  # should not be 'proc':'delay' after all.
