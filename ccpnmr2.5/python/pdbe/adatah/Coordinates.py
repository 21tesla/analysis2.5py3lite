#
# TODO Should PDB stuff go in here as well?. In separate class?
#

#
# General functions
#

def readListCoordinateFiles(coordFilePath,fileEnd = '.pdb'):
    
  fileList = []
  fileEndLen = len(fileEnd)

  files = os.listdir(coordFilePath)
  
  for fileName in files:
    if fileName[-fileEndLen:] == fileEnd:
      fileList.append(os.path.join(coordFilePath,fileName))
      
  fileList.sort()
  
  return fileList

#
# Generic coordinate handling class - to be used in conjunction with DataHandler
#

import os

from ccpnmr.format.general.Constants import ccpNmr_kw, originalFormat_kw

class CoordinateHandler:
  
  """
  For use with dataHandler classes, can also be used as standalone if self. info set correctly
  """
  
  # WAS FORMERLY: readCnsCoordinates, readPseudoPdbCoordinates
  
  def readCoordinateFiles(self,formatName,coordFileDir, fileEnd = '.pdb', forceNamingSystemName = None):

    #
    # Read the list of coordinate files
    #

    fileList = readListCoordinateFiles(coordFileDir, fileEnd = fileEnd)

    if not fileList:
      raise self.DataHandlerError("No %s coordinate files available!" % formatName)
      return

    keywds = {}

    if self.presets.has_key('readCoordinates'):

      scriptPresets = self.presets['readCoordinates']

      if scriptPresets.has_key('keywds'):

        keywds = scriptPresets['keywds']
        
    if forceNamingSystemName:
      keywds['forceNamingSystemName'] = forceNamingSystemName

    self.formatObjectDict[formatName].readCoordinates(fileList, molSystem = self.molSystem, strucGen = self.strucGen, linkAtoms = False, minimalPrompts = 1, **keywds)
  
    print "  Read %d files from %s in %s format..." % (len(fileList),coordFileDir,formatName)

  def checkCoordinateAtomConsistency(self):

    from ccp.general.Util import findAllSysNamesByChemAtomOrSet

    """
    TODO: this script also has to check atoms that were never included because they are not recognised - set this as
    application data in DataFormat.py!
    """

    models = self.structureEnsemble.models
    missingAtomsInfo = []

    origFormat = None
    appData = self.structureEnsemble.findFirstModel().findFirstApplicationData(application = ccpNmr_kw, keyword = originalFormat_kw)
    if appData:
      origFormat = appData.value

    origNamingSystemName = self.structureEnsemble.atomNamingSystem

    for coordChain in self.structureEnsemble.sortedCoordChains():

      origChainCode = ""
      if origFormat:
        appData = coordChain.findFirstApplicationData(application = origFormat, keyword = 'originalCode')
        if appData:
          origChainCode = appData.value

      for coordResidue in coordChain.sortedResidues():

        origSeqCode = 0
        if origFormat:
          appData = coordResidue.findFirstApplicationData(application = origFormat, keyword = 'originalSeqCode')
          if appData:
            origSeqCode = appData.value

        for coordAtom in coordResidue.sortedAtoms():

          coordModels = set()
          for coord in coordAtom.sortedCoords():
            coordModels.add(coord.model)

          missingModels = models.difference(coordModels)

          if missingModels:

            origAtomName = ""
            if origFormat:
              appData = coordAtom.findFirstApplicationData(application = origFormat, keyword = 'originalName')
              if appData:
                origAtomName = appData.value
              else:
                chemAtomSysNames = findAllSysNamesByChemAtomOrSet(coordResidue.residue.chemCompVar.chemComp,[coordAtom.atom.chemAtom],origNamingSystemName)
                origAtomName = chemAtomSysNames[0].sysName

            modelSerials = [model.serial for model in missingModels]
            modelSerials.sort()

            missingAtomInfo = (coordChain.chain.code,coordResidue.residue.seqId,coordAtom.atom.name,
                               origChainCode,origSeqCode,origAtomName,modelSerials)

            missingAtomsInfo.append(missingAtomInfo)

    return missingAtomsInfo


  def getCoordinateAtomConsistencyWarnings(self):

    warnings = []

    missingAtomsInfo = self.checkCoordinateAtomConsistency()

    for missingAtomInfo in missingAtomsInfo:

      (chainCode,seqId,atomName,origChainCode,origSeqCode,origAtomName,modelSerials) = missingAtomInfo

      message = "Error: Atom '%s.%d.%s' (originally '%s.%d.%s') missing in model(s) %s" % (chainCode,
                                                                                           seqId,
                                                                                           atomName,
                                                                                           origChainCode,
                                                                                           origSeqCode,
                                                                                           origAtomName,
                                                                                           ','.join([str(modelSerial) for modelSerial in modelSerials]))
      warnings.append(message)

    return warnings
    
  def getCoordinateMappingErrors(self):
  
    """
    Finds information about chains/residues/atoms that could not be mapped when importing with the FormatConverter
    """
    
    appData = self.structureEnsemble.findFirstApplicationData(application=ccpNmr_kw,keyword='mappingErrors')

    if appData:
      mappingErrors = eval(appData.value)
      
      # Extra check
      hasMappingErrors = False
      for mapKey in mappingErrors.keys():
        if mappingErrors[mapKey]:
          hasMappingErrors = True
          break
      
      if not hasMappingErrors:
        mappingErrors = {}
      
    else:
      mappingErrors = {}
      
    return mappingErrors
