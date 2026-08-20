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
  
  if pName == 'ccp.molecule.Validation':
    # Fix Validation
    fixValidation(topObj, delayDataDict)
  
  elif pName == 'memops.Implementation':
    
    # reset names in packageLocators
    for pl in doGet(topObj, emptyDict).get('packageLocators', emptyList):
      if pl.targetName == 'ccp.molecule.StructureValidation':
        pl.__dict__['targetName'] = 'ccp.molecule.Validation'
    
    topObjByGuid = topObj.__dict__.get('topObjects')
    validationStores = [x for x in topObjByGuid.values() 
                        if x.packageShortName == 'VALD']
    for validationStore in validationStores:
      fullLoadValidationStore(validationStore)
  
  # 
  from memops.format.compatibility.upgrade.v_2_0_a3 import Minor as Minor20a3
  Minor20a3.correctData(topObj, delayDataDict, toNewObjDict, mapping)
    
      
fixingValidation = False

def fixValidation(topObj, delayDataDict):
  """ Trigger validation load from correct location
  Fix ValidationStore-Software link
  NB this triggers full load before partial load ha finished
  """
  
  global fixingValidation
  
  if fixingValidation:
    # end of second pass. Do no more
    fixingValidation = False
    
  else:
    # end fo first pass. 
    fixingValidation = True
    
    # fix software link
    root = topObj.root
    methodStore = root.currentMethodStore
    if methodStore is None:
      methodStore = root.newMethodStore(name='auto')
    software = methodStore.findFirstSoftware(name='unknown')
    if software is None:
      software = methodStore.newSoftware(name='unknown', version='unknown')
    #
    topObj.software = software
    
    # reload 
    fullLoadValidationStore(topObj)

def fullLoadValidationStore(topObj):
  """hard load ValidationStore from old location
  """ 
  
  from memops.format.xml import Util as xmlUtil
  from memops.universal import Io as uniIo
  from memops.format.xml import XmlIO
  root = topObj.memopsRoot
  locator = (root.findFirstPackageLocator(targetName='ccp.molecule.Validation')
             or root.findFirstPackageLocator(targetName='any'))
  repository = locator.findFirstRepository()
  #repository = topObj.activeRepositories[0]
  fileLocation = repository.getFileLocation('ccp.molecule.StructureValidation')
  filePath = uniIo.joinPath(fileLocation, xmlUtil.getTopObjectFile(topObj))
  XmlIO.loadFromStream(open(filePath), topObject=topObj, 
                       topObjId=topObj.guid, partialLoad=False)
