"""
======================COPYRIGHT/LICENSE START==========================

Minor.py: Data compatibility handling

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

# NB this file will only be used as part of Minor upgrades

from memops.general.Implementation import ApiError

def correctData(topObj, delayDataDict, toNewObjDict, mapping=None):
  """ update topObj object tree using information in delayDataDict
  May be used either to postprocess a file load (minor upgrade)
  or as part of an in-memory data transfer (major upgrade)
  
  topObj is the MemopsRoot in the new tree
  toNewObjDict is _ID:newObj for minor 
    and oldObj/oldObjId:newObj for major upgrades
    
    NB this function does nothing without further additions
  """
  
  
  emptyDict = {}
  emptyList = []
  doGet = delayDataDict.get
  
  pName = topObj.packageName
  
  if pName == 'ccpnmr.Analysis':
    fixAnalysis(topObj, delayDataDict, toNewObjDict, mapping)
  #
  elif pName == 'ccp.nmr.Nmr':
    # Fix Nmr
    fixNmr(topObj, delayDataDict)
  # 
  from memops.format.compatibility.upgrade.v_2_0_b3 import Minor as Minor20b3
  Minor20b3.correctData(topObj, delayDataDict, toNewObjDict, mapping)


def fixNmr(topObj, delayDataDict):
  """ remap defunct RefExperiment
  """
  emptyDict = {}
  emptyList = []
  doGet = delayDataDict.get
  
  memopsRoot = topObj.parent
  topObjByGuid = memopsRoot.__dict__.get('topObjects')
  
  for xpr in doGet(topObj, emptyDict).get('experiments', emptyList):
    expDict = doGet(xpr, emptyDict)
    
    # fix NmrExpPrototype mapping for defunct types
    setNmrExpPrototypeLink(xpr, 'refExperiment', topObjByGuid, delayDataDict,
                           remapPrototypeLink)
    for xd in expDict.get('expDims', emptyList):
      setNmrExpPrototypeLink(xd, 'refExpDim', topObjByGuid, delayDataDict,
                             remapPrototypeLink)
      for xdr in doGet(xd, emptyDict).get('expDimRefs', emptyList):
        setNmrExpPrototypeLink(xdr, 'refExpDimRef', topObjByGuid, delayDataDict,
                               remapPrototypeLink)


def remapPrototypeLink(keyList):
  
  from memops.format.compatibility.upgrade.v_2_0_4 import Minor as Minor204
  
  guid = keyList[0]
  
  # map to different keys
  if guid == 'cam_wb104_2008-01-15-16-06-39_00031':
    # 'H[N[CA[HA]]] removed RefExperiment - change to correct version
    if keyList[1] == 4:
      keyList[1] = 8
  
  Minor204.remapPrototypeLink(keyList)


from memops.format.compatibility.upgrade.v_2_0_a1.Minor import setNmrExpPrototypeLink
    
    
def fixAnalysis(topObj, delayDataDict, toNewObjDict, mapping):
  """ The spectrumWindows slot actually contains spectrumWindowPanes
  With the attributes of the SpectrunmWindow put in as 'delay'
  The fixing function puts the delayed attributes in the right place
  rearranges the delayDataDict to be read correctly when setting child 
  links later, and sets child links from SpectrumWIndow down
  """
  
  simpleAttrs = ['isCanvasLabelShown', 'isCanvasMidpointShown', 'isIconified', 
                 'isXSliceDrawn', 'isXTickShown', 'isYSliceDrawn', 
                 'isYTickShown', 'serial', 'stripAxis', 'useMultiplePeakLists', 
                 'useOverrideRegion']
  
  emptyDict = {}
  emptyList = []
  doGet = delayDataDict.get
  
  dd = doGet(topObj, emptyDict)
  
  spectrumWindows = dd.get('spectrumWindows', emptyList) 
  for swp in spectrumWindows:   # really SpectrumWindowPanes
    
    # get next dictionary
    swpdd = doGet(swp, emptyDict)
    
    # get SpectrumWindow
    sw = swp.spectrumWindow
    swdict = sw.__dict__
    
    # fill in delayDataDict
    delayDataDict[sw] = {'spectrumWindowPanes':[swp]}
    
    # fix SpectrumWindow.name
    swdict['name'] = swp.name
    
    # set simple hicard==1 attrs:
    for tag in simpleAttrs:
      vals = swpdd.get(tag)
      if vals:
        setattr(sw, tag, vals[0])
    
    # set location
    tag = 'location'
    vals = swpdd.get(tag)
    if vals:
      setattr(sw, tag, vals)
    
    # set spectrumWindowGroups crosslink
    tag = 'spectrumWindowGroups'
    vals = swpdd.get(tag)
    if vals:
      setattr(sw, tag, [toNewObjDict.get(x) for x in vals])
    
    # set children, recursively
    from memops.xml.Implementation import linkChildData
    linkChildData(delayDataDict, sw, mapping, linkTopToParent=2)
