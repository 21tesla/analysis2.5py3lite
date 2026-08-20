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
  
  if pName == 'ccp.molecule.MolStructure':
    # Fix MolStructure
    fixMolStructure(topObj, delayDataDict, toNewObjDict)
  
  #
  elif pName == 'ccp.nmr.Nmr':
  
    # Fix Nmr
    fixNmr(topObj, delayDataDict)
  
  from memops.format.compatibility.upgrade.v_2_1_0 import Minor as Minor210
  Minor210.correctData(topObj, delayDataDict, toNewObjDict, mapping)


def fixMolStructure(topObj, delayDataDict, toNewObjDict):
  
  emptyDict = {}
  emptyList = []
  doGet = delayDataDict.get
    
  # make new Atoms if necessary, and make list of Atoms
  allAtoms = []
  for chain in doGet(topObj, emptyDict).get('coordChains', emptyList):
    for res in doGet(chain, emptyDict).get('residues', emptyList):
      for atom in doGet(res, emptyDict).get('atoms', emptyList):
        
        allAtoms.append(atom)
        
        coords= doGet(atom, emptyDict).get('coords', emptyList)
        if coords:
          
          # group coords by their altLocCode
          altLocs = {}
          for coord in coords:
            altLocCode = doGet(coord, emptyDict).get('altLocationCode', 
                                                     emptyList)
            if altLocCode:
              altLocCode = altLocCode[0]
              ll = altLocs.get(altLocCode)
              if ll:
                ll.append(coord)
              else:
                altLocs[altLocCode] = [coord]
          
          del coords[:] # We want them removed before objects get created
          
          # 
          if altLocs:
            altLocCodes = list(sorted(altLocs.keys()))
            atom.altLocationCode = altLocCodes[0]
            
            for altLocCode in altLocCodes[1:]:
              newAtom = res.newAtom(name=atom.name, altLocationCode=altLocCode,
                                    access = atom.access, 
                                    applicationData = atom.applicationData)
              allAtoms.append(newAtom)
              ll = altLocs[altLocCode]
              for coord in ll:
                coord.__dict__['atom'] = newAtom
  
  # set up system - atoms
  topObj.orderedAtoms = allAtoms
  for ii, atom in enumerate(allAtoms):
    atom.index = ii
  
  # set up system - models
  models = doGet(topObj, emptyDict).get('models', emptyList)
  ll = [(x.serial,x) for x in models]
  ll.sort()
  models = [x[1] for x in ll]
  for ii,model in enumerate(models):
    model.index = ii
  nModels = len(models)
  
  # set up system - data matrices
  # NB must set entry in parent __dict__ as this is done while reading
  nAtoms = len(allAtoms)
  xx = topObj.newDataMatrix(name='bFactors', shape=(nModels,nAtoms))
  topObj.__dict__['dataMatrices']['bFactors'] = xx
  xx = topObj.newDataMatrix(name='occupancies', shape=(nModels,nAtoms))
  topObj.__dict__['dataMatrices']['occupancies'] = xx
  xx = topObj.newDataMatrix(name='coordinates', shape=(nModels,nAtoms,3))
  topObj.__dict__['dataMatrices']['coordinates'] = xx
  
  # set data
  for model in models:
    
    occupancies = [1.0] * nAtoms
    bFactors = [0.0] * nAtoms
    coordinates = bFactors * 3
    
    setBFactors = False
    setOccupancies = False
    coordIds = delayDataDict[model].get('coords', emptyList)
    
    #coords = [toNewObjDict.get(x) for x in coordIds]
    
    nfound = 0
    for coordId in coordIds:
      coord = toNewObjDict.pop(coordId)    # Want them gone before final check
      coordDict = delayDataDict.pop(coord) # Want them gone before final check
      index = coord.atom.index
      
      ll = coordDict.get('occupancy')
      if ll:
        setOccupancies = True
        occupancies[index] = ll[0]
      
      ll = coordDict.get('bFactor')
      if ll:
        setBFactors = True
        bFactors[index] = ll[0]
      
      for ii,tag in enumerate(('x','y','z')):
        ll = coordDict.get(tag)
        if ll:
          nfound += 1
          coordinates[3*index + ii] = ll[0]
    
    model.setSubmatrixData('coordinates', coordinates)
    if setOccupancies:
      model.setSubmatrixData('occupancies', occupancies)
    if setBFactors:
      model.setSubmatrixData('bFactors', bFactors)
  #
  topObj.purge()



def fixNmr(topObj, delayDataDict):
  """
  """
  
  emptyDict = {}
  emptyList = []
  doGet = delayDataDict.get
  
  memopsRoot = topObj.parent
  topObjByGuid = memopsRoot.__dict__.get('topObjects')

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
  
  guid = keyList[0]
  
  # map to different keys
  if guid == "ccpn_rhf22_2009-04-02-16-03-22-040_00001":
    # Remap {CA|Cca}CONH 143 to H{CA|Cca}CONH 59
    keyList[0] = "cam_wb104_2008-01-15-16-06-40_00025"
    refExpMap = {1:38, 2:5, 3:39, 8:8, 24:24, }
    keyList[1] = refExpMap[keyList[1]]
  
  elif guid == "ccpn_rhf22_2009-04-02-15-29-19-145_00001":
    # Remap {CA|Cca}NH (142) to H{CA|Cca}NH (70)
    keyList[0] = "cam_wb104_2008-01-15-16-06-40_00036"
    # RefExpMap is identity, just keep them
    keyList[1] = refExpMap[keyList[1]]
  
  elif guid == "ccpn_rhf22_2009-04-03-14-39-42-994_00001":
    # Remap 152 CA[N] to 243 HCA[N]
    keyList[0] = "ccpn_rhf22_2009-04-24-16-19-16-187_00001"
    refExpMap = {1:2}
    keyList[1] = refExpMap[keyList[1]]
  
  elif guid == "ccpn_rhf22_2009-04-02-16-10-54-048_00001":
    # Remap 144 C_cCONH.relayed to 281 HNCO_C.relayed
    # NB refExp 15 (projection) is not mapped. It is assumed it was never used.
    keyList[0] = "Expts_vicky_2010-12-15-16-02-18-326_00001"
    refExpMap = {1:7, 3:8, 6:9, 17:10}
    keyList[1] = refExpMap[keyList[1]]
  
  
from memops.format.compatibility.upgrade.v_2_0_a1.Minor import setNmrExpPrototypeLink
