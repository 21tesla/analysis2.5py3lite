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

transferTypeMap = {'NOESY':'through-space',
                   'DipDip':'through-space',
                   'SpinDiff':'through-space',
                   'CP':'onebond',
                   'TOCSY':'relayed',
                   'Jnonlocal':'through-space',
                  }

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
  
  
  #
  if pName == 'ccp.nmr.Nmr':
  
    # Fix Nmr
    fixNmr(topObj, delayDataDict)

  elif pName == 'ccp.nmr.NmrExpPrototype':
    # Fix Nmr
    fixNmrExpPrototype(topObj, delayDataDict)

  # 
  from memops.format.compatibility.upgrade.v_2_0_5 import Minor as Minor205
  Minor205.correctData(topObj, delayDataDict, toNewObjDict, mapping)

def fixNmr(topObj, delayDataDict):
  """
  """
  
  # modify transfer type
  emptyDict = {}
  emptyList = []
  doGet = delayDataDict.get
  
  memopsRoot = topObj.parent
  topObjByGuid = memopsRoot.__dict__.get('topObjects')
  
  
  for xpr in doGet(topObj, emptyDict).get('experiments', emptyList):
    expDict = doGet(xpr, emptyDict)
    for expTransfer in expDict.get('expTransfers', emptyList):
      transferType = doGet(expTransfer).get('transferType')
      if transferType:
        transferType = transferType[0]
      if transferType in transferTypeMap:
        transferType = transferTypeMap[transferType]
      setattr(expTransfer, 'transferType', transferType)

  # remap defunct RefExperiment
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
  
  from memops.format.compatibility.upgrade.v_2_0_5 import Minor as Minor205
  
  guid = keyList[0]
  
  # map to different keys
  if guid in ('cam_wb104_2008-01-15-16-06-39_00003',
              "ccpn_rhf22_2009-04-16-18-01-38-822_00001"):
    # COSY variant - remap to [1]
    keyList[0] = "cam_wb104_2008-01-15-16-06-39_00001"
  
  elif guid == "ccpn_rhf22_2009-04-03-14-39-03-775_00001":
    # Remap COCA to CACO
    keyList[0] = "ccpn_rhf22_2009-04-03-14-37-40-370_00001"
    keyList[1] = 2
  
  elif guid == "cam_wb104_2008-01-15-16-06-40_00059":
    # C_C NOESY duplicate
    keyList[0] = "ccpn_rhf22_2009-04-02-16-22-29-921_00001"
  
  elif guid == "ccpn_rhf22_2009-04-17-17-47-01-401_00001":
    # Dept=135
    keyList[0] = "ccpn_rhf22_2009-04-09-17-35-30-819_00001"
  
  elif guid in ("cam_wb104_2008-01-15-16-06-40_00055",
                "cam_wb104_2008-01-15-16-06-40_00056",
                "cam_wb104_2008-01-15-16-06-40_00057"):
    # Superfluous relaxation variant
    if keyList[1] == 2:
      keyList[1] = 1
  
  elif guid == "ccpn_rhf22_2009-05-27-15-54-34-222_00001":
    # Superfluous duplicate
    keyList[0] = "ccpn_rhf22_2009-04-15-17-09-17-471_00001"
    refExpMap = {1:1, 40:38, 41:39, 42:2, 43:5, 44:6}
    keyList[1] = refExpMap[keyList[1]]
  
  elif guid == "ccpn_rhf22_2009-04-14-15-10-00-409_00001":
    # Change 2D filtered NOESY from exp 208 to exp 9
    if keyList[1] == 2:
      keyList[0] = "cam_wb104_2008-01-15-16-06-39_00009"
      keyList[1] = 1
    elif keyList[1] == 4:
      keyList[0] = "cam_wb104_2008-01-15-16-06-39_00009"
      keyList[1] = 3
  
  elif guid == "ccpn_rhf22_2009-04-14-15-18-59-112_00001":
    # Change 2D filtered NOESY from exp 210 to exp 74
    if keyList[1] == 3:
      keyList[0] = "cam_wb104_2008-01-15-16-06-40_00040"
      keyList[1] = 1
  
  elif guid == "Expts_vicky_2010-12-15-15-46-16-259_00001":
    # map exp 279 to exp 175
    keyList[0] =  "ccpn_rhf22_2009-04-09-17-59-59-074_00001"
    if keyList[1] == 3:
      keyList[1] == 4
    elif keyList[1] == 4:
      keyList[1] == 5
  
  elif guid == "Expts_vicky_2010-12-15-15-44-46-557_00001":
    # map exp 278 to exp 174
    keyList[0] =  "ccpn_rhf22_2009-04-09-17-57-16-877_00001"
    if keyList[1] == 3:
      keyList[1] == 4
    elif keyList[1] == 4:
      keyList[1] == 5
  
  elif guid == "Expts_vicky_2010-12-15-15-18-35-679_00001":
    # map exp 275 to exp 148
    keyList[0] =  "ccpn_rhf22_2009-04-03-14-27-15-011_00001"
    if keyList[1] == 3:
      keyList[1] == 4
    elif keyList[1] == 4:
      keyList[1] == 5
  
  elif guid == "Expts_vicky_2010-12-15-15-13-54-881_00001":
    # map exp 274 to exp 102
    keyList[0] =  "cam_wb104_2008-01-15-16-06-40_00063"
  
  Minor205.remapPrototypeLink(keyList)
  
from memops.format.compatibility.upgrade.v_2_0_a1.Minor import setNmrExpPrototypeLink


def fixNmrExpPrototype(topObj, delayDataDict):
  """ remap defunct RefExperiment
  """
  emptyDict = {}
  emptyList = []
  doGet = delayDataDict.get
  
  memopsRoot = topObj.parent
  topObjByGuid = memopsRoot.__dict__.get('topObjects')
  
  for xgr in doGet(topObj, emptyDict).get('expGraphs', emptyList):
    xgrDict = doGet(xgr, emptyDict)
    for expTransfer in xgrDict.get('expTransfers', emptyList):
      transferType = doGet(expTransfer).get('transferType')
      if transferType:
        transferType = transferType[0]
      if transferType in transferTypeMap:
        transferType = transferTypeMap[transferType]
      setattr(expTransfer, 'transferType',transferType)
