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
from memops.format.compatibility.upgrade.v_2_0_b3 import General as General20b3

# functions to convert fullKeys before application
fullKeyConverters = {}
fullKeyConverters.update(General20b3.fullKeyConverters)

# guids of elements that should be treated as old
elemsTreatedAsOld = set((
 'www.ccpn.ac.uk_Fogh_2008-09-24-15:20:52_00002', # spectrumWindowPanes
 'www.ccpn.ac.uk_Fogh_2008-09-24-15:20:52_00008',
 'www.ccpn.ac.uk_Fogh_2008-09-24-15:20:52_00006',
 'www.ccpn.ac.uk_Fogh_2008-09-24-15:20:52_00004'
)).union(General20b3.elemsTreatedAsOld)

# pairs of element guids that should be treated as matching, e.g. whe n
# a single element must match with several elements in subclasses

# these are links from MemopsRoot to individual TopObjects that all have 
# changed their guids in a model change.
elementPairings = []

def extraMapChanges(globalMapping):
  """ Extra map changes specific for a given step
  """
  
  #General20b3.extraMapChanges(globalMapping)
  
  dd = globalMapping.get('ANAL')
  if dd is not None:
    ANALmapChanges(globalMapping)
  
  # Set links Nmr -> nmrExpPrototype to delayed
  for guid in ('www.ccpn.ac.uk_Fogh_2006-08-16-18:20:06_00008',
   'www.ccpn.ac.uk_Fogh_2006-08-16-18:23:00_00002',
   'www.ccpn.ac.uk_Fogh_2006-08-16-18:20:05_00025'
  ):
    globalMapping['mapsByGuid'][guid]['proc'] = 'delay'
  
  General20b3.extraMapChanges(globalMapping)

def ANALmapChanges(globalMapping):
  """ changes for ANAL (ccpnmr.Analysis) package
  """
  
  # make changes for SpectrumWindow -> SpectrumWindowPane split
  # attributes that belong to SpectrumWindow are copied to
  # SpectrumWindowPane and set to 'delay'
  # they are attached to the actual SpectrumWindowPane, and
  # copied up to the proper location in SpectrumWindow by postprocessing
  loadMaps = globalMapping['loadMaps']
  SpWinMap = loadMaps['ANAL.SpectrumWindow']
  SpWinContent = SpWinMap['content']
  SpPaneMap = loadMaps['ANAL.SpectrumWindowPane']
  SpPaneContent = SpPaneMap['content']
  
  # copy maps from Window to WindowPane and set proc to delay
  ll = SpWinMap['headerAttrs'] + SpWinMap['simpleAttrs'] + SpWinMap['cplxAttrs']
  ll.remove('name')
  ll.remove('access')
  ll.remove('applicationData')
  
  for tag in ll:
    dd = SpWinContent[tag]
    dd['proc'] = 'delay'
    SpPaneContent[tag] = dd
  
  # set proc to skip for link into SpectrumWindow
  dd = loadMaps['ANAL.SpectrumWindowGroup']
  dd['content']['spectrumWindows']['proc'] = 'skip'
  
  # set custom constructor
  SpWinMap['constructor'] = constructWindowAndPane
  
  # copy WinPane attributes to Win
  aDict = {
   'headerAttrs':['aspectRatio'],
   'simpleAttrs':['sliceRange'],
   'cplxAttrs':['spectrumWindowViews', 'slicePanels', 'axisPanels']
  }
  for ss in ('headerAttrs', 'simpleAttrs', 'cplxAttrs'):
    ll = aDict[ss]
    SpWinMap[ss].extend(ll)
    for tag in ll:
      SpWinContent[tag] = SpPaneContent[tag]
      loadMaps['ANAL.SpectrumWindow.%s' % tag] = (
        loadMaps['ANAL.SpectrumWindowPane.%s' % tag])
  SpWinMap['children'].extend(aDict['cplxAttrs'])
  
  
  # disable child link setting to SpectrumWindow
  dd = globalMapping['loadMaps']['ANAL.AnalysisProject']
  dd['children'].remove('spectrumWindows')
  #
  return


def constructWindowAndPane(analysisProject):
  """ Special constructor.
  Creates both SpectrumWindow and SpectrumWindowPane
  Called for the slot fo the former but returns the latter
  """
  spectrumWindow = analysisProject.newSpectrumWindow()
  spectrumWindowPane = spectrumWindow.newSpectrumWindowPane(serial=1)
  return spectrumWindowPane
