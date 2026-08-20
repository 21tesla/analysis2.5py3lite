from pdbe.nmrStar.IO.Constants import importVersionSep, mappingImportLocation
from memops.general.Constants import currentModelVersion                                                                   

def reverseDict(origDict):

  newDict = {}
  for origKey in origDict.keys():
    newDict[origDict[origKey]] = origKey

  return newDict

def getCcpnVersion():

  #
  # Get the CCPN version
  #
    
  return str(currentModelVersion)

def getLatestNmrStarVersion():
  
  #
  # Get the latest NMR-STAR version (should be more sophisticated)
  #
  
  version = '3.1'
  
  return version
  
def getNmrStarDict(version):
  
  #
  # Try to get a version specific NMR-STAR dictionary, use default if none
  # available
  #

  importVersion = version.replace('.',importVersionSep)
  
  try:
    importName = "nmrStarDict__%s__" % importVersion
    nmrStarDict = __import__("pdbe.nmrStar.IO.%s" % importName,{},{},importName)
  except:
    #print "  Warning: no NMR-STAR version %s dictionary available, using default." % version
    nmrStarDict = __import__("pdbe.nmrStar.IO.nmrStarDict",{},{},"nmrStarDict")
  
  return nmrStarDict

#
# Ccpn to NMR-STAR import defs
#

def getCcpn2NmrStarImportNames(ccpnVersion,nmrStarVersion):

  importCcpnVersion = ccpnVersion.replace('.',importVersionSep)
  importNmrStarVersion = nmrStarVersion.replace('.',importVersionSep)

  conversionClassFormat = "Ccpn%s_To_NmrStar%s"
  
  specificImport = conversionClassFormat % ('__%s_' % importCcpnVersion,'__%s__' % importNmrStarVersion)
  defaultImport = conversionClassFormat % ('','')
  
  return (specificImport,defaultImport)

def getCcpn2NmrStar(ccpnVersion,nmrStarVersion,exportClass = None):

  (specificImport,defaultImport) = getCcpn2NmrStarImportNames(ccpnVersion,nmrStarVersion)
  
  ConversionClass = getConversionClass(specificImport,defaultImport)

  # Warning here: if using default import, could clash with nmrStarVersion (but then problems anyway...)
  conversionClass = ConversionClass(exportClass, ccpnVersion = ccpnVersion, nmrStarVersion = nmrStarVersion)

  return conversionClass
  
def getCcpn2NmrStarConstants(ccpnVersion,nmrStarVersion):

  (specificImport,defaultImport) = getCcpn2NmrStarImportNames(ccpnVersion,nmrStarVersion)
  
  return getConstants_ByVersion(specificImport,defaultImport) 

#
# NMR-STAR to Ccpn import defs
#

def getNmrStar2CcpnImportNames(nmrStarVersion,ccpnVersion):

  importNmrStarVersion = nmrStarVersion.replace('.',importVersionSep)
  importCcpnVersion = ccpnVersion.replace('.',importVersionSep)

  conversionClassFormat = "NmrStar%s_To_Ccpn%s"
  
  specificImport = conversionClassFormat % ('__%s_' % importNmrStarVersion,'__%s__' % importCcpnVersion)
  defaultImport = conversionClassFormat % ('','')
  
  return (specificImport,defaultImport)

def getNmrStar2Ccpn(nmrStarVersion,ccpnVersion,importClass = None):
  
  (specificImport,defaultImport) = getNmrStar2CcpnImportNames(nmrStarVersion,ccpnVersion)

  ConversionClass = getConversionClass(specificImport,defaultImport)

  # Warning here: if using default import, could clash with nmrStarVersion (but then problems anyway...)
  conversionClass = ConversionClass(importClass, ccpnVersion = ccpnVersion, nmrStarVersion = nmrStarVersion)

  return conversionClass

def getNmrStar2CcpnConstants(nmrStarVersion,ccpnVersion):
  
  # Special case: get the ccpn to nmr-star dictionary and reverse
  (specificImport,defaultImport) = getCcpn2NmrStarImportNames(ccpnVersion,nmrStarVersion)

  ccpn2NmrStarConstants = getConstants_ByVersion(specificImport,defaultImport) 
  
  nmrStar2CcpnDict = {}
  
  for constantsKey in ccpn2NmrStarConstants:
    nmrStar2CcpnDict[constantsKey] = reverseDict(ccpn2NmrStarConstants[constantsKey])
  
  return nmrStar2CcpnDict

#
# General import stuff
#

def getConstants_ByVersion(specificImport,defaultImport):

  #
  # Try to get a version specific constants dictionary, use default if none
  # available
  #
  
  constantsDictNameFormat = "constants_%s"
  
  try:
    convModule = __import__("pdbe.nmrStar.IO.Constants_ByVersion",{},{},constantsDictNameFormat % specificImport)
    conversionDict = getattr(convModule,constantsDictNameFormat % specificImport)
    
  except:
    print "  Warning: No %s constants dictionary available, using %s default." % (specificImport,defaultImport)
    convModule = __import__("pdbe.nmrStar.IO.Constants_ByVersion",{},{},constantsDictNameFormat % defaultImport)
    conversionDict = getattr(convModule,constantsDictNameFormat % defaultImport)  
  
  return conversionDict

def getConversionClass(specificImport,defaultImport):

  #
  # Try to get a version specific conversion class, use default if none
  # available
  #  
  
  try:
    convModule = __import__("%s.%s" % (mappingImportLocation,defaultImport),{},{},specificImport)
    ConversionClass = getattr(convModule,specificImport)
    
  except:
    print "  Warning: No %s conversion class available, using %s default." % (specificImport,defaultImport)
    convModule = __import__("%s.%s" % (mappingImportLocation,defaultImport),{},{},defaultImport)
    ConversionClass = getattr(convModule,defaultImport)  
  
  return ConversionClass

